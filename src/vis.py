import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as sp

# Accelerometer tuning params
slope = 0.00062885
intercept = 0.47789158

# Curvefit params
def sin_fit(time, A, omega, phi, offset):
    return A*np.sin(omega*time + phi) + offset

files = [f"{i+1}-{j+1}" for j in range(3) for i in range (5)]
print(files)

ria = np.loadtxt(f"data/4-2.txt", delimiter = ',')
plt.plot(ria[:,0], ria[:,1])
plt.show()

# for name in files:
#     ria = np.loadtxt(f"data/{name}.txt", delimiter = ',')
#     plt.plot(ria[:,0], ria[:,1])
#     plt.show()
