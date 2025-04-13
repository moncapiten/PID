import numpy as np
import os
from scipy.optimize import curve_fit
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

def Zener_model(X, I0, gam, Eg, n):
    return PN_model(X, I0, gam, Eg, n)

# Parametri

names = ['PN diode', 'Schottky diode', 'Zener diode 3.3','Zener diode 2.7']
filenames = [
    'IV-T_dependence_20250404_100546/',
    'IV-T_dependence_20250408_110423/',
    'IV-T_dependence_20250408_124236/',
]
for choose_diode in range(1,len(filenames)):
    try: base_path = os.path.join('Data', filenames[choose_diode - 1])
    except Exception: print("Path non trovata")

    temps = np.arange(16, 71)
    offsets, pulls = [-1, 5], range(3)
    flag_return_T = flag_return_V = True

    # Lettura dati
    V_all, T_all, I_all = [], [], []

    def filename(T, offsets, pull, direction):
        dir_str = 'Go' if direction == 1 else 'Return'
        return f'IV_T{T:.2f}_V{offsets[0]:.2f}_{offsets[1]:.2f}_{pull}Pull_T{dir_str}.txt'

    for direction in [1, 2] if flag_return_T else [1]:
        for T in reversed(temps) if direction == 2 else temps:
            for pull in pulls:
                path = os.path.join(base_path, filename(T, offsets, pull, direction))
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
    p0_dict = {
        1: [2e-8, 1.5, 1.8, 2],
        2: [2e-8, 2, 0.3, 2],
        3: [2e-8, 1.8, 1.8, 4],
        4: [2e-8, 1.8, 1.8, 4]
    }
    models = {1: PN_model, 2: Schottky_model, 3: Zener_model, 4: Zener_model}
    p0 = p0_dict[choose_diode]
    model = models[choose_diode]
    beta, _ = curve_fit(model, X, I_all, p0, maxfev=int(1e6))

    # Plot
    fig = plt.figure(figsize=(12, 5))
    ax1 = fig.add_subplot(121, projection='3d')
    ax2 = fig.add_subplot(122, projection='3d')

    I_fit = model(X, *beta)
    for ax, scale, title in zip([ax1, ax2], ['linear', 'log'], ['Linear Plot', 'Semilog Plot']):
        Z = np.abs(I_all) if scale == 'log' else I_all
        Z_fit = np.abs(I_fit) if scale == 'log' else I_fit
        ax.scatter(T_all, V_all, Z, label='Raw Data', color='blue')
        ax.plot(T_all, V_all, Z_fit, 'r.', label='Fitted Model')
        ax.set_xlabel('Temperature [K]')
        ax.set_ylabel('Voltage [V]')
        ax.set_zlabel('Current [$\mu$A]')
        ax.set_title(title)
        if scale == 'log': ax.set_zscale('log')
        ax.legend()

    fig.suptitle(f'IV-T curve for {names[choose_diode - 1]}')
    plt.tight_layout()
    plt.show()
