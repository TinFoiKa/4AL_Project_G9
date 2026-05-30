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

def cal_accel(raw):
    return slope*raw + intercept

# Curvefit params
def sin_fit(time, A, omega, phi, offset):
    return A*np.sin(omega*time + phi) + offset

def line_intersect_method(ria):
    """
    Function to observe when the graph crosses the mean point,
    doubled for half-period, quadrupled for full period
    """

    t = ria[:, 0] / 1000.0
    x = ria[:, 1] / 100.0
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

def minima_method(ria):
    """
    Function much like line_intersect_method to instead use
    sp.find_peaks on the negative data and observe diffs between
    minima (to take advantage of more stable sensor distances)
    """
    t = ria[:, 0] / 1000.0
    x = ria[:, 1] / 100.0
    ind = sps.find_peaks(-x, distance = 25)[0]

    minima_t = []
    for index in ind:
        minima_t.append(t[index])
    if len(minima_t) < 2:
        return None

    intervals = np.diff(minima_t)
    print(intervals)
    half_period = np.mean(intervals)
    return float(half_period * 2), half_period

files = [f"{i+1}-{j+1}" for j in range(3) for i in range (5)]
print(files)

amps = [0.02, 0.03, 0.05, 0.1, 0.14]

a_amps = [5000, 5000, 5000, 5000, 5000]

with open("analysis.csv", 'w', newline="", encoding="utf-8") as f:
    writer = csv.writer(f)
    writer.writerow(["Height", "Period"])

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
            # acc = cal_accel(ria[start:end, 3])


            guesses = [amps[idx], 3.78*2, 0, 0.30]
            # guesses = [a_amps[idx], 3.78, 0, 1000]

            # Create Graph
            popt, _ = sp.curve_fit(sin_fit, time, height, p0=guesses)
            # popt, _ = sp.curve_fit(sin_fit, time, acc, p0=guesses)
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
            T = w_sf * T_sf + w_li * T_li

            writer.writerow([popt[0], T])
            
            ax.axvspan(time[0], time[0] + T, color='orange', alpha=0.25)
            ax.text(0.02, 0.9, f'T={T:.3f}s', transform=ax.transAxes)

            ax.plot(time, height, '.', label='data')
            # ax.plot(time, acc, '.', label = 'data')
            ax.plot(fit_time, fit_height, '-', label='fit')
            ax.legend()
            ax.set_title(name)

        fig.supxlabel("Time (sec)")
        fig.supylabel("Height (m)")
        fig.tight_layout()
        plt.show()

# ria = np.loadtxt(f"data/4-2.txt", delimiter = ',')
# plt.plot(ria[:,0], cal_accel(np.array(ria[:,3])))
# plt.show()

# for name in files:
#     ria = np.loadtxt(f"data/{name}.txt", delimiter = ',')
#     plt.plot(ria[:,0], ria[:,3])
#     plt.show()
