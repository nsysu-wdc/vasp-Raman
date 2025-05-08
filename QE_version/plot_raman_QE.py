#! /opt/Program_Files/Miniconda3/3.11-23.10.0-1.app/vaspPy/3.11.2025.01.14/.venv/bin/python

# %% import global
import numpy as np
import matplotlib.pyplot as plt
import os


def read_spectrum_data(filename):
    """
    read Raman and IR spectrum data
    return frequency (cm-1), IR intensity, Raman intensity
    """
    data = []
    with open(filename, 'r') as file:
        for line in file:
            if line.strip().startswith('#') or len(line.strip()) == 0:
                continue  # skip comment or empty lines
            tokens = line.strip().split()
            if len(tokens) < 6:
                continue  # skip incomplete lines
            freq_cm1 = float(tokens[1])
			# freq_thz = float(tokens[2])
            ir_intensity = float(tokens[3])
            raman_intensity = float(tokens[4])
            data.append((freq_cm1, ir_intensity, raman_intensity))
    data = np.array(data)
    return data[:,0], data[:,1], data[:,2]

def plot_spectrum(x_values, raman_intensities=None, ir_intensities=None, width=10, mode='gaussian', 
                  xlabel='wavenumber (cm$^{-1}$)', ylabel='Intensity'):
    """
    plot Gaussian or peak (stem plot) 
    x_values: frequency (cm-1)
    raman_intensities: Raman intensity
    ir_intensities: IR intensity
    width: Gaussian (cm-1) --> only for mode='gaussian'
    mode: 'gaussian' or 'peak'
    """
    plt.figure(figsize=(10, 6))
    
    x_plot = np.linspace(min(x_values)-100, max(x_values)+100, 1000)
    x_min=x_plot[0]
    x_max=x_plot[-1]

    if mode == 'gaussian':
        raman_spectrum = np.zeros_like(x_plot)
        ir_spectrum = np.zeros_like(x_plot)
        
        if raman_intensities is not None:
            for xi, intensity in zip(x_values, raman_intensities):
                raman_spectrum += intensity * np.exp(-((x_plot - xi)**2) / (2 * width**2))
            plt.plot(x_plot, raman_spectrum, color='red', label='Raman')
        if ir_intensities is not None:
            plt.plot(x_plot, ir_spectrum, color='blue', label='Ir')
            for xi, intensity in zip(x_values, ir_intensities):
                ir_spectrum += intensity * np.exp(-((x_plot - xi)**2) / (2 * width**2))
    
    elif mode == 'peak':
        if raman_intensities is not None:
            plt.stem(x_values, raman_intensities, markerfmt=" ", basefmt=" ", label='Raman')
        if ir_intensities is not None:
            plt.stem(x_values, ir_intensities, markerfmt=" ", basefmt=" ", label='Ir')
    
    else:
        raise ValueError("mode 必須是 'gaussian' 或 'peak'")

    plt.legend()
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.xlabel(xlabel, fontsize=12)
    plt.ylabel(ylabel, fontsize=12)
    plt.xlim(x_min, x_max)

    plt.tight_layout()
	# plt.show()
    plt.savefig('output.png', format='png', dpi=500)


# %%
if __name__ == "__main__":
    # Modify the file name and mode here to switch
    filename = 'data.txt'

    # linux
    os.system("awk '/# mode/ {count=25} count-- >= 0' dynmat.out > " + filename)

    freq, ir_intensity, raman_intensity = read_spectrum_data(filename)
    plot_spectrum(freq, raman_intensities=raman_intensity, ir_intensities=ir_intensity, width=10, mode='gaussian')
    # plot_spectrum(freq, raman_intensities=raman_intensity, ir_intensities=None, mode='peak')
# %%
