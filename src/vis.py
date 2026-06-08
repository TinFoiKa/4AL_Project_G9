import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as sp
import scipy.signal as sps
import csv

# Accelerometer tuning params
slope = 0.00062885
intercept = 0.47789158

"""
There is a more rigorous way:

Sinusoids are ONLY a solution for SHM oscillations, 
Draw a line at the calculated average, see how often the data
intersects this line, take the average of the distances of these instances
(over a certain tolerance so that we don't capture noise)
"""

h = 0.31
L = 0.635
g = 9.81

def cal_accel(raw):
    return slope*raw + intercept

# Curvefit params
def sin_fit(time, A, omega, phi, offset):
    return A*np.sin(omega*time + phi) + offset

def damping():
    """
    Process all data files and, for each file, find the heights of the
    first and last peak in the interval [start, end).

    Returns a list of tuples: (filename, first_peak_height, last_peak_height)
    and writes the results to `damping_peaks.csv` in the workspace root.
    """
    start = 100
    end = 500
    results = []
    only = "1-1"
    for name in only:
        ria = np.loadtxt(f"data/5-1.txt", delimiter=',')
        s = start
        e = min(end, ria.shape[0])
        height = ria[s:e, 1] / 100.0
        time = ria[s:e, 0] / 1000.0

        # find peaks in the height data
        peaks, props = sps.find_peaks(-height, distance=10)
        print(peaks)

        guesses = [0.02, 3.78*2, 0, 0.30]

        popt, _ = sp.curve_fit(sin_fit, time, height, p0=guesses)

        print("theta: " + str(_angle_conversion(h+ 2*popt[0])))
        print("omega: " + str(popt[1]/2))

        if len(peaks) == 0:
            first_h = None
            last_h = None
            dt=1
            print("no peaks")
        else:
            first_h = float(height[peaks[1]])
            last_h = float(height[peaks[-1]])
            dt = float(time[peaks[-1]]) - float(time[peaks[1]])
            print("dt: " + str(dt))

        beta = -np.log(last_h/first_h)/dt
        print("beta: " + str(beta))
        results.append(beta)

        o_0 = np.sqrt((popt[1]/2)**2 + beta**2)
        print("omega_0: " + str(o_0))
        print("per diff: " + str(o_0/(popt[1]/2)-1))
    
    print(results)

    return np.mean(results)

def line_intersect_method(ria):
    """
    Function to observe np.diffs when the graph crosses the mean point,
    doubled for half-period, quadrupled for full period
    """
    start = 100
    end = 500
    t = (ria[start:end, 0] - ria[start, 0]) / 1000.0
    x = ria[start:end, 1] / 100.0
    mean = np.mean(x)
    y = x - mean
    tol = max(0.1 * np.std(y), 1e-3)
    crosses = []
    for i in range(len(y) - 1):
        if y[i] * y[i + 1] < 0 and (abs(y[i]) > tol or abs(y[i + 1]) > tol):
            frac = abs(y[i]) / (abs(y[i]) + abs(y[i + 1]))
            crosses.append(t[i] + frac * (t[i + 1] - t[i]))
    if len(crosses) < 2:
        return None
    intervals = np.diff(crosses)
    
    qt_period = np.mean(intervals)
    return float(qt_period * 4), float(qt_period * 2)

def _angle_conversion(x):
    arg = 1 - (x-h)/L
    print (arg)
    arg = np.clip(arg, -1.0, 1.0) # keep within bounds of arccos
    return np.degrees(np.arccos(arg))

def minima_method(ria):
    """
    Function much like line_intersect_method to instead use
    sp.find_peaks on the negative data and observe diffs between
    minima (to take advantage of more stable sensor distances)
    """
    start = 100
    end = 500
    t = (ria[start:end, 0] - ria[start, 0]) / 1000.0
    x = ria[start:end, 1] / 100.0
    ind = sps.find_peaks(-x, distance = 25, height = [-0.32,-0])[0]

    minima_t = []
    for index in ind:
        minima_t.append(t[index])
    if len(minima_t) < 2:
        return None

    intervals = np.diff(minima_t)
    print(intervals)
    half_period = np.mean(intervals)
    return float(half_period * 2), float(half_period)

files = [f"{i+1}-{j+1}" for j in range(3) for i in range (5)]
print(files)

amps = [0.02, 0.03, 0.05, 0.1, 0.14]

a_amps = [5000, 5000, 5000, 5000, 5000]

with open("final.csv", 'w', newline="", encoding="utf-8") as f:
    writer = csv.writer(f)
    writer.writerow(["Height", "Period (sin)", "amp_uncertainty", "omega_uncertainty"])
    r_amps = []
    # Process files in batches of 5 and plot each batch on one figure
    batch_size = 5
    for batch_start in range(0, len(files), batch_size):
        batch = files[batch_start:batch_start + batch_size]
        fig, axes = plt.subplots(nrows=len(batch), ncols=1, figsize=(8, 2 * len(batch)), squeeze=False)
        for idx, name in enumerate(batch):
            ax = axes[idx, 0]
            ria = np.loadtxt(f"data/{name}.txt", delimiter=',')

            # 100 to 500 is a heuristic window
            start, end = 100, 500
            end = min(end, ria.shape[0])
            time = (ria[start:end, 0] - ria[start, 0]) / 1000.0
            height = ria[start:end, 1] / 100.0
            acc = cal_accel(ria[start:end, 3])

            # guesses = [a_amps[idx], 3.78, 0, 1000]

            # popt, _ = sp.curve_fit(sin_fit, time, acc, p0=guesses)

            guesses = [amps[idx], 3.78*2, 0, 0.30]

            # Create Graph
            popt, pcov = sp.curve_fit(sin_fit, time, height, p0=guesses)
            fit_time = np.linspace(time[0], time[-1], 200)
            fit_height = sin_fit(fit_time, *popt)

            w_sf = 0.3
            T_sf = 2*np.pi / popt[1] * 2
            
            w_li = 0.3
            T_li = line_intersect_method(ria)[0]

            w_mi = 0.4
            T_mi = minima_method(ria)[0]

            print(T_sf, T_li, T_mi)

            # for now, weighted sum
            T = T_sf

            # writer.writerow([popt[0], T_sf, T_li, T_mi])
            writer.writerow([popt[0], T_sf, np.sqrt(pcov[0,0]), np.sqrt(pcov[1,1])])
            r_amps.append(popt[0])
            ax.axvspan(time[0], time[0] + T, color='orange', alpha=0.25)
            ax.text(0.02, 0.9, f'T={T:.3f}s, A={_angle_conversion(h+2*popt[0]):.1f}deg', transform=ax.transAxes)

            ax.plot(time, height, '.', label='data')
            # ax.plot(time, acc, '.', label = 'data')
            ax.plot(fit_time, fit_height, '-', label='fit')
            ax.legend()
            ax.set_title(name)

        fig.suptitle("Height vs. Time by Sinusoidal Curve Fitting")
        fig.supylabel("Height (m)")

        # fig.suptitle("Acceleration vs. Time by Sinusoidal Curve Fitting")
        # fig.supylabel("Acceleration (m/s^2)")
        fig.supxlabel("Time (sec)")
        
        fig.tight_layout()
        plt.show()

# print([float(amp) + 0.31 for amp in r_amps])
# # ria = np.loadtxt(f"data/4-2.txt", delimiter = ',')
# # plt.plot(ria[:,0], cal_accel(np.array(ria[:,3])))
# # plt.show()

# # for name in files:
# #     ria = np.loadtxt(f"data/{name}.txt", delimiter = ',')
# #     plt.plot(ria[:,0], ria[:,3])
# #     plt.show()
