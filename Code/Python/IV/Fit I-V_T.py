import os
import re
import numpy as np
import scipy.optimize as opt
import matplotlib.pyplot as plt

# Costante di Boltzmann in eV/K
k_B = 8.617333262e-5  # eV/K

def model_function(T, I0, gamma, E_B):
    return I0 * T**gamma * np.exp(-E_B / (k_B * T))

def extract_temperature(filename):
    match = re.search(r'IV_([0-9]+\.?[0-9]*)_', filename)
    if match:
        return float(match.group(1))
    return None

def process_files(folder_path):
    data_list = []
    temperatures = []
    
    for filename in os.listdir(folder_path):
        if filename.startswith("IV_") and filename.endswith(".txt"):
            T = extract_temperature(filename)
            if T is None:
                continue
            
            file_path = os.path.join(folder_path, filename)
            data = np.loadtxt(file_path)
            I_values = data[:, 0]
            V_values = data[:, 1]
            
            temperatures.append(T)
            data_list.append(np.mean(I_values))  # Media di I per ogni T
    
    temperatures = np.array(temperatures)
    data_list = np.array(data_list)
    
    # Fit non lineare
    popt, _ = opt.curve_fit(model_function, temperatures, data_list, p0=[1e-12, 1.0, 0.5])
    
    return temperatures, data_list, popt

def plot_results(temperatures, data_list, popt):
    T_fit = np.linspace(min(temperatures), max(temperatures), 100)
    I_fit = model_function(T_fit, *popt)
    
    plt.scatter(temperatures, data_list, label='Dati sperimentali', color='red')
    plt.plot(T_fit, I_fit, label=f'Fit: I0={popt[0]:.2e}, gamma={popt[1]:.2f}, E_B={popt[2]:.3f} eV', color='blue')
    plt.xlabel('Temperatura (K)')
    plt.ylabel('Corrente media (I)')
    plt.legend()
    plt.grid()
    plt.show()

# Esempio di utilizzo
folder_path = "./dati"  # Cambiare con il percorso della cartella
T_vals, I_vals, params = process_files(folder_path)
plot_results(T_vals, I_vals, params)