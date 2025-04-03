import os
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def load_and_average_data(files):
    data_list = []
    for file in files:
        data = np.loadtxt(file, delimiter="\t", skiprows=2, unpack=True)
        data_list.append(data)
    
    return np.mean(data_list, axis=0)

def plot_3d_surface(V, I, T, title, ax):
    ax.plot_trisurf(V, I, T, cmap='viridis', edgecolor='none')
    ax.set_xlabel('Voltage [V]')
    ax.set_ylabel('Current [A]')
    ax.set_zlabel('Temperature [°C]')
    ax.set_title(title)

def main(path):
    files = [f for f in os.listdir(path) if f.endswith(".txt")]
    temperature_groups = {}
    
    for file in files:
        temp_str = file.split('_')[1][1:]
        if temp_str not in temperature_groups:
            temperature_groups[temp_str] = []
        temperature_groups[temp_str].append(os.path.join(path, file))
    
    for temp_str, file_list in temperature_groups.items():
        data_avg = load_and_average_data(file_list)
        T = float(temp_str)
        
        # Estrazione delle colonne usando unpack=True
        Td_go, Td_return, Vcc, Vd_go, Vd_return, ErrVd_go, ErrVd_return, Vr_go, Vr_return, ErrVr_go, ErrVr_return, Id_go, Id_return = data_avg
        
        fig = plt.figure(figsize=(12, 5))
        ax1 = fig.add_subplot(121, projection='3d')
        plot_3d_surface(Vd_go, Id_go, np.full_like(Vd_go, T), f'T={T}°C - GO', ax1)
        
        ax2 = fig.add_subplot(122, projection='3d')
        plot_3d_surface(Vd_return, Id_return, np.full_like(Vd_return, T), f'T={T}°C - RETURN', ax2)

        plt.tight_layout()
        plt.show()

path = "Data/Dati di prova/IV-T_dependence_20250402_192640"
main(path)
