#! /data/chuanglab/Program_Files/Miniconda3/3.10-23.3.1-0.app/vaspPy/3.10.2023.08.14/.venv/bin/python

import matplotlib.pyplot as plt
import numpy as np

data = np.loadtxt('total_dos.dat')

freq_THz = data[:, 0]
freq_cm1 = freq_THz * 33.35641
freq_cm1 = np.append(freq_cm1, 600)

dos = data[:, 1]
dos = np.append(dos, 0)

with open('total_dos_cm1.txt', 'w') as fid:
  for row in range(0,len(freq_cm1)):
    fid.write('   %.08f    %.08f\n' % (freq_cm1[row], dos[row]))
print('save to --> total_dos_cm1.txt ')

plt.figure(figsize=(10, 3))
plt.plot(freq_cm1, dos, color='darkred', linewidth=1.8)

plt.xlabel('Frequency (cm$^{-1}$)', fontsize=14)
plt.ylabel('Phonon DOS', fontsize=14)
# plt.title('Phonon Density of States', fontsize=14)
plt.grid(True, linestyle='--', alpha=0.5)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.xlim([0, 600])
plt.tight_layout()

# plt.show()
plt.savefig('dos_convert.png', format='png', dpi=500)
print('show in --> dos_convert.png ')
