"""
Reproduce Roman's instrumental broadening profile
"""

import numpy as np
import matplotlib.pyplot as plt

# BOSS spectrograph instrumental broadening
R = 1800 # resolving power          dimensionless

# width of instrument broadening in pixels
pixel_spacing = 1e-4
pixel_sigma = 1 / (R * 2 * np.sqrt(2 * np.log(2)) * (np.power(10, pixel_spacing) - 1))

total = 0
width = 3

instrument_profile = np.full(width * 2 + 1, fill_value=np.nan)

for j,i in enumerate(range(-width, width + 1)):
    instrument_profile[j] = np.exp(-0.5 * i * i / (pixel_sigma * pixel_sigma))
    total += instrument_profile[j]

instrument_profile /= total
# plt.plot(instrument_profile, 'o', label='R=1800', c='orange')
# plt.plot(instrument_profile, '--', c='orange')
# plt.legend('R=1800')

print(instrument_profile)
print(R)
