import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as sp

# constants
h = 0.31
L = 0.635
g = 9.81

# data
ria = np.loadtxt("analysis.csv", delimiter = ",")
d_height = 2 * ria[:,0] + h
T_s = ria[:,1]
T_mean = ria[:,2]
T_min = ria[:,3]
T_z = 2*np.pi * np.sqrt(L/g)

def _logarithm(t, a, b, c):
    return a * np.log(b - t) + c

def _angle_conversion(x):
    arg = 1 - (x-h)/L
    print (arg)
    arg = np.clip(arg, -1.0, 1.0) # keep within bounds of arccos
    return np.degrees(np.arccos(arg))

def _taylor_series(x):
    theta_rad = np.radians(x)
    sin_half = np.sin(theta_rad / 2)
    # Taylor series correction: T ≈ T_z * (1 + (1/4)sin²(θ/2) + (9/64)sin⁴(θ/2) + ...)
    correction = 1 + (1/4) * sin_half**2 + (9/64) * sin_half**4
    return T_z * correction

def methods():
    """
    Comparison of period determining methods
    """
    plt.scatter(d_height, T_s, label = "Sin Fit")
    plt.scatter(d_height, T_mean, label = "Mean Intersect")
    plt.scatter(d_height, T_min, label = "Minima Diff")
    plt.xlabel("Drop Heights (m)")
    plt.ylabel("Computed Periods (sec)")
    plt.title("Computed Periods by Method")

def angles():
    """
    Sinusoidal Fit to observe how period changes with angle
    Adds references for small angle approximation, numerical integration, taylor series expansion
    angle is calculated with _angle_conversion() from release height
    """
    sim_thetas = _angle_conversion(d_height)
    sim_data = [1.6049999999999998, 1.615, 1.625, 1.6249999999999998, 1.64, 1.6049999999999998, 1.615, 1.62, 1.6300000000000001, 1.6300000000000001, 1.61, 1.615, 1.625, 1.6300000000000001, 1.635]

    thetas = _angle_conversion(d_height)
    plt.errorbar(thetas, T_s, marker = "x", label = "Empirical Data")
    plt.scatter(thetas, T_s, marker = "x", label = "Empirical Data")
    series = _taylor_series(thetas)
    plt.scatter(thetas, series, label = "Taylor Series")

    plt.scatter(sim_thetas, sim_data, marker = "*", label = "RK-4 Simulation")
    
    plt.axhline(y = T_z, c= "green", linestyle = "--", label = "Small Angle Approx")
    plt.xlabel("Release Angle (degrees)")
    plt.ylabel("Period (s)")
    plt.title("Angle vs. Period for Large Angle Pendulums")

angles()
plt.legend()
plt.show()
