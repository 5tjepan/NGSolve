
import numpy as np
from scipy.interpolate import interp1d

import matplotlib.pyplot as plt

ataB = np.array([0.186,0.278,0.373,0.464,0.511,0.558,0.652,0.7,0.746,0.793,0.839,0.884,0.929,0.965,1.02])
Ptot = np.array([30.8, 57.6, 94.8, 141.2, 169.6, 202.4, 281.2, 328.4, 380.4, 437.6, 501.2, 571.2, 648, 715.6, 827.6])
Peps=np.array([24.8, 44.4, 71.2,102,119.2,136.4,172.4,191.6,211.6,232.4,254,275.6,298,316,344.8])
Plec=np.array([4.4,10.4,18.4,27.6,32.8,37.6,48.4,54.0,59.6,65.6,71.6,77.6,83.6,88.8,96.4])
Ppec=np.array([6.4,12.8,23.6,39.2,50.4,66.0,108.8,136.8,168.8,205.2,247.2,295.2,350.0,399.6,482.8])

avaB=np.array([0.182,0.274,0.366,0.458,0.548,0.639,0.731,0.82,0.865,0.909,0.945,1.0])
avaPtot=np.array([30.8,56.8,93.2,140.8,201.6,279.2,376.8,490.8,558,632,698.8,807.2])
avaPeps=np.array([23.6,41.6,66.4,98.4,134.8,175.6,216.0,253.6,272.4,292.4,308.4,336.8])
avaPlec= np.array([4,9.6,17.2,26.4,36.4,47.2,58.0,68.4,74.4,80.0,85.2,93.2])
avaPpec=np.array([7.2,15.2,26.8,42.4,66.8,103.6,160.4,237.6,285.6,340.0,390.0,470.8])


P=Plec
B=ataB
PBspline = interp1d(B, P, kind='linear')

B_new = np.array([0.22,0.28,0.34, 0.37, 0.41, 0.45, 0.49, 0.52, 0.56, 0.60, 0.64, 0.67,0.71, 0.75, 0.79, 0.81, 0.84, 0.86, 0.88, 0.90, 0.92, 0.94, 0.96, 0.97])
P_new = PBspline(B_new)

Pm=np.array([34.7, 53.3, 75.9, 93.1, 111.8, 133.3, 155.4, 181.4, 209.3, 240.5, 275.0, 309.9, 351.0, 393.9,442.0,468.0,521.3,551.2,581.1,614.9,644.8,679.9,713.7,752.7])

print(np.round(P_new,1))

plt.plot(B_new, Pm, 'o')
plt.plot(B_new, P_new, '-')
plt.legend()
plt.show()



