#! /opt/Program_Files/Miniconda3/3.11-23.10.0-1.app/vaspPy/3.11.2025.01.14/.venv/bin/python

# %% import global
import matplotlib.pyplot as plt
import numpy as np

# %% function


def load_raman_data(filename):
    """
    --> Read the Raman data file and return raman_shift and intensity
    """
    raman_shift = []
    intensity = []
    with open(filename, 'r') as file:
        for line in file:
            if line.strip().startswith('#') or not line.strip():
                continue
            parts = line.strip().split()
            freq = float(parts[1])
            activity = float(parts[4])
            raman_shift.append(freq)
            intensity.append(activity)
    return np.array(raman_shift), np.array(intensity)


def plot_peaks(raman_shift, intensity, invert_x=False):
    """
    --> Draw the Raman plot of vertical spikes.
    """
    plt.figure(figsize=(10, 6))
    for x0, y0 in zip(raman_shift, intensity):
        plt.plot([x0, x0], [0, y0], color='blue', linewidth=1.5)
    plt.xlabel('Raman Shift (cm$^{-1}$)', fontsize=12)
    plt.ylabel('Intensity (a.u.)', fontsize=12)
    plt.title('Raman Spectrum (Stick Plot)', fontsize=14)
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)
    plt.xlim(min(raman_shift) - 10, max(raman_shift) + 10)
    if invert_x:
        plt.gca().invert_xaxis()
    plt.tight_layout()
    # plt.show()
    plt.savefig('Raman-peaks.png', format='png', dpi=500)


def plot_gaussian(raman_shift, intensity, linewidth=10, invert_x=True, annotate_peaks=False):
    """
    --> Draw Gaussian superimposed Raman plot with optional peak annotations.
    """
    x_min = min(raman_shift) - 50
    x_max = max(raman_shift) + 50
    x = np.linspace(x_min, x_max, 4000)
    y = np.zeros_like(x)

    for mu, amp in zip(raman_shift, intensity):
        y += amp * np.exp(-((x - mu) ** 2) / (2 * linewidth ** 2))

    plt.figure(figsize=(10, 6))
    plt.plot(x, y, color='blue', lw=1.8)
    plt.xlabel('Raman Shift (cm$^{-1}$)', fontsize=12)
    plt.ylabel('Intensity (a.u.)', fontsize=12)
    plt.title('Raman Spectrum (Gaussian Broadening)', fontsize=14)
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)
    plt.xlim(x_min, x_max)

    # Peak annotation
    if annotate_peaks:
        for mu, amp in zip(raman_shift, intensity):
            plt.text(mu, amp + max(intensity) * 0.02,  # ?? peak ????
                     f'{mu:.1f}', fontsize=14, color='red',
                     ha='center', va='bottom', rotation=90)

    if invert_x:
        plt.gca().invert_xaxis()

    plt.tight_layout()
    # plt.show()
    plt.savefig('Raman-gaussian.png', format='png', dpi=500)


def plot_raman(filename, mode="gaussian", linewidth=10):
    ## Main control function, optional mode 'gaussian' or 'peak ##
    raman_shift, intensity = load_raman_data(filename)
    if mode == "gaussian":
        plot_gaussian(raman_shift, intensity, linewidth)
    elif mode == "peak":
        plot_peaks(raman_shift, intensity)
    else:
        raise ValueError("mode should be 'gaussian' or 'peak' !! ")


# %%
if __name__ == "__main__":
    # Modify the file name and mode here to switch
    filename = "vasp_raman.dat"
    plot_raman(filename, mode="gaussian", linewidth=10)
    # plot_raman(filename, mode="peak")

# %%
