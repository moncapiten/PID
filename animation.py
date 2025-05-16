import numpy as np
import os
from scipy.optimize import curve_fit
from matplotlib import animation
import matplotlib.pyplot as plt

# Modelli
def PN_model(X, I0, gam, Eg, n):
    V, T = X
    k, Vt = 8.617333e-5, 13e-3
    Is = I0 * T**gam * np.exp(-Eg / (n * k * T))
    return Is * (np.exp(V / (n * Vt)) - 1) * 1e6

def Schottky_model(X, I0, gam, Esh, n):
    V, T = X
    k, Vt = 8.617333e-5, 13e-3
    Is = I0 * T**gam * np.exp(-Esh / (k * T))
    return Is * (np.exp(V / (n * Vt)) - 1) * 1e6

def Zener_model(X, I0, Eg, n, VZ0=1.0, I_leak=0):
    V, T = X
    k = 8.617333e-5  # eV/K
    q = 1.0          # carica (unità compatibili con V e k)

    Is = I0 * T**(3/2) * np.exp(-Eg / (k * T))  # corrente di saturazione
    Vz_T = np.full_like(V, 1.54)             # soglia Zener costante
    I = np.zeros_like(V)

    # Regime diretto: V > 0
    mask_fwd = V > 0
    I[mask_fwd] = Is[mask_fwd] * (np.exp(q * V[mask_fwd] / (n * k * T[mask_fwd])) - 1)

    # Regime Zener (breakdown inverso): V < -VZ0
    mask_zener = V < -Vz_T
    I[mask_zener] = -Is[mask_zener] * (np.exp(-q * (V[mask_zener] + Vz_T[mask_zener]) / (n * k * T[mask_zener])) - 1)

    # Regime di perdita: -VZ0 < V < 0
    mask_leak = (~mask_fwd) & (~mask_zener)
    I[mask_leak] = I_leak

    return I * 1e6  # [μA]

# Parametri
names = {'PN diode':[-1, 5], 'Schottky diode':[-1, 5], 'Zener diode 1':[-3, 5], 'Zener diode 2':[-5, 5]}
filenames = [
    'IV-T_dependence_20250404_100546/',
    'IV-T_dependence_20250408_110423/',
    'IV-T_dependence_20250408_124236/',
    'IV-T_dependence_20250408_/'
]
for choose_diode in range(1,4):
    base_path = os.path.join('Data', filenames[choose_diode - 1])
    temps = np.arange(16, 71)
    offsets = list(names.values())[choose_diode-1]

    # Lettura dati
    pulls = range(3)
    flag_return_T = flag_return_V = True
    V_all, T_all, I_all = [], [], []

    def filename(T, offsets, pull, direction):
        dir_str = 'Go' if direction == 1 else 'Return'
        return f'IV_T{T:.2f}_V{offsets[0]:.2f}_{offsets[1]:.2f}_{pull}Pull_T{dir_str}.txt'

    # Lettura delle pull
    for direction in [1, 2] if flag_return_T else [1]:
        for T in reversed(temps) if direction == 2 else temps:
            for pull in pulls:
                try:
                    path = os.path.join(base_path, filename(T, offsets, pull, direction))
                except Exception: print("Path non trovata")
                try:
                    data = np.loadtxt(path)
                    V_all.extend(data[:, 4 if direction == 1 else 5])
                    T_all.extend(data[:, 0 if direction == 1 else 1])
                    I_all.extend(data[:, 12 if direction == 1 else 13])
                except Exception:
                    continue

    T_all = np.array(T_all) + 273.15
    V_all = np.array(V_all)
    I_all = np.array(I_all) * 1e6
    X = np.vstack((V_all, T_all))

    # Fit
    '''p0_dict = {
        1: [2e-8, 1.5, 1.8, 2],   # [I0, γ, Eg, n]     Bipola model
        2: [2e-8, 2, 0.3, 2],     # [I0, γ, Eg, n]     Schottky model
        3: [2e-8, 1.8, 4, -1.0],  # [I0, Eg, n, VZ0]   Zener model
        4: [2e-8, 1.8, 4, -1.0]   # [I0, Eg, n, VZ0]    
    }
    models = {1: PN_model, 2: Schottky_model, 3: Zener_model, 4: Zener_model}
    p0 = p0_dict[choose_diode]
    model = models[choose_diode]
    #beta, _ = curve_fit(model, X, I_all, p0, maxfev=100000)
    #beta, _ = curve_fit(model, X, I_all, maxfev=1000000)
    #chi2_red = np.sum((I_all - I_fit)**2) / (len(I_all) - len(beta))

    # Parametri di fit
    param_names = {
        1: ['I0', 'γ', 'Eg', 'n'],
        2: ['I0', 'γ', 'Esh', 'n'],
        3: ['I0', 'Eg', 'n', 'VZ0'],
        4: ['I0', 'Eg', 'n', 'VZ0']
    }
    param_str = ', '.join(f'{name}={val:.2e}' for name, val in zip(param_names[choose_diode], beta))'''

    # Plot
    fig = plt.figure(figsize=(6, 5))
    ax = fig.add_subplot(111, projection='3d')

    # Dati logaritmici
    Z = np.abs(I_all)

    # Imposta zlim se necessario
    valid_Z = Z[Z > 0]
    if len(valid_Z) > 0:
        z_min = 0.1 * np.min(valid_Z)
        z_max = np.max(valid_Z)
        ax.set_zlim(z_min, z_max)

    # Scatter plot
    ax.scatter(T_all, V_all, Z, label='Raw Data', marker='o', edgecolors='blue', facecolors='none')
    ax.set_xlabel('Temperature [K]')
    ax.set_ylabel('Voltage [V]')
    ax.set_zlabel('Current [μA]')
    ax.set_title('Semilog Plot')
    ax.set_zscale('log')
    ax.legend()

    # Titolo principale
    fig.suptitle(f'IV-T curve for {list(names.keys())[choose_diode - 1]}\n', fontsize=14)
    plt.tight_layout(rect=[0, 0.05, 1, 1])

    # Rotazione animata
    def rotate(angle):
        ax.view_init(elev=30, azim=angle)

    ani = animation.FuncAnimation(fig, rotate, frames=np.arange(0, 360, 10), interval=10)
    ani.save('log_rotation.gif', writer='pillow', fps=20)

    plt.show()