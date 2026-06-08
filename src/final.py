import numpy as np
import matplotlib.pyplot as plt
from sim import analyse

# constants
h = 0.31
dh = 0.001
L = 0.635
g = 9.81

# data
ria = np.loadtxt("final.csv", delimiter = ",")
d_height = 2 * ria[:,0] + h
T_s = ria[:,1]
sig_A = ria[:,2]
sig_O = ria[:,3]

i_height = np.array([17.5, 22.5, 27.5, 32.5, 35.0]) # intended drop heights

# small angle
T_z = 2*np.pi * np.sqrt(L/g)

def _logarithm(t, a, b, c):
    return a * np.log(b - t) + c

def m_se_T():
    """between-trial uncertainty for each intended drop thru standard error of the mean"""
    seT = []
    m = []
    for arr in np.transpose(np.reshape(T_s, [-1,5])):
        mn = np.mean(arr)
        stT = np.std(arr)
        seT.append(stT/np.sqrt(3))
        m.append(mn)
        print(arr)
    
    return m, seT

def exp_T_uncertainty():
    """uncertainty from experimental setup"""
    d_T = []
    T_trans = np.transpose(np.reshape(T_s, [-1,5]))
    for i, arr in enumerate(np.transpose(np.reshape(sig_O, [-1,5]))):
        omega = 2*np.pi/T_trans[i]
        d_T.append(np.mean(T_trans[i] * arr / omega))
    
    return d_T

def _angle_conversion(x):
    """returns angle, error of angle conversion"""
    arg = 1 - (x-h)/L
    print (arg)
    arg = np.clip(arg, -1.0, 1.0) # keep within bounds of arccos

    d_th = 1/np.sqrt(1-(arg)**2)*dh/L
    return np.degrees(np.arccos(arg)), d_th

def _taylor_series(x):
    theta_rad = np.radians(x)
    sin_half = np.sin(theta_rad / 2)
    # Taylor series correction: T ≈ T_z * (1 + (1/4)sin²(θ/2) + (9/64)sin⁴(θ/2) + ...)
    correction = 1 + (1/4) * sin_half**2 + (9/64) * sin_half**4
    return T_z * correction

def angles():
    """
    Sinusoidal Fit to observe how period changes with angle
    Adds references for small angle approximation, numerical integration, taylor series expansion
    angle is calculated with _angle_conversion() from release height
    """
    sim_thetas = _angle_conversion(d_height)[0]
    sim_data = analyse(sim_thetas)

    thetas, dth = _angle_conversion(d_height)
    m_T, se_T = m_se_T()
    plt.errorbar(thetas, T_s, xerr= dth, yerr = se_T*3, linestyle="none", marker = "x", label = "Empirical Data")
    series = _taylor_series(thetas)
    plt.scatter(thetas, series, label = "Taylor Series", c= 'red')

    plt.scatter(sim_thetas, sim_data, marker = "*", label = "RK-4 Simulation", c="orange")
    
    plt.axhline(y = T_z, c= "green", linestyle = "--", label = "Small Angle Approx")
    plt.xlabel("Release Angle (degrees)")
    plt.ylabel("Period (s)")
    plt.title("All 15 Trials distributed on Axes")

def cons_angles():
    """
    For only 5 data points where each group of intended drop is combined into one
    """
    sim_thetas = i_height
    sim_data = analyse(sim_thetas)

    thetas = i_height
    heightT = np.transpose(np.reshape(_angle_conversion(d_height)[0], [-1,5]))
    dth = np.array([np.std(x)/np.sqrt(3) for x in heightT])
    m_T, se_T = m_se_T()
    plt.errorbar(thetas, m_T, se_T, dth, linestyle="none", marker = "x", label = "Empirical Data")
    series = _taylor_series(thetas)
    plt.scatter(thetas, series, label = "Taylor Series", c="red")

    th_conv_errs = _angle_conversion(np.array([0.339, 0.358, 0.382, 0.409, 0.424]))[1]

    print(f"dtheta: {np.sqrt(dth**2 + th_conv_errs**2)}")
    print(f"period_tay: {str(series)}")
    print(f"period_rk4: {str(sim_data)}")
    print(f"period: {str(m_T)}")
    print(f"dT: {np.sqrt(np.array(exp_T_uncertainty())**2 + np.array(se_T)**2)}")

    plt.scatter(sim_thetas, sim_data, marker = "*", label = "RK-4 Simulation", c="orange")
    
    plt.axhline(y = T_z, c= "green", linestyle = "--", label = "Small Angle Approx")
    plt.xlabel("Release Angle (degrees)")
    plt.ylabel("Period (s)")
    plt.title("Angle vs. Period for Large Angle Pendulums")

cons_angles()
plt.legend()
plt.show()
angles()
plt.legend()
plt.show()
