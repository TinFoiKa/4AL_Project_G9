"""Large-angle pendulum RK4 simulation.

This simulates a pendulum using the exact equation of motion without the small-angle approximation.
The user can modify `initial_angle_deg` and `string_length` directly in the script.
"""

import math
import matplotlib.pyplot as plt
import scipy.optimize as sp
import numpy as np

def damped_sin_fit(time, A, omega, phi, offset):
    """We select the damped fit because """
    return A*np.sin(omega*time + phi) + offset

def pendulum_derivatives(state, length, gravity=9.81, damping=0.0):
    """Return positional derivatives [d(theta)/dt, d(omega)/dt], equivalently [omega, alpha].

    `damping` is a linear damping coefficient (units 1/s) that produces a
    torque proportional to angular velocity: -damping * omega.
    """
    theta, omega = state
    dtheta_dt = omega
    domega_dt = -(gravity / length) * math.sin(theta) - damping * omega
    return dtheta_dt, domega_dt


def rk4_step(state, dt, length, gravity=9.81, damping=0.0):
    """Perform one RK4 step for the pendulum state."""
    k1 = pendulum_derivatives(state, length, gravity, damping)
    k2_state = (state[0] + 0.5 * dt * k1[0], state[1] + 0.5 * dt * k1[1])
    k2 = pendulum_derivatives(k2_state, length, gravity, damping)

    k3_state = (state[0] + 0.5 * dt * k2[0], state[1] + 0.5 * dt * k2[1])
    k3 = pendulum_derivatives(k3_state, length, gravity, damping)

    k4_state = (state[0] + dt * k3[0], state[1] + dt * k3[1])
    k4 = pendulum_derivatives(k4_state, length, gravity, damping)

    theta_next = state[0] + (dt / 6.0) * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0])
    omega_next = state[1] + (dt / 6.0) * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1])
    return theta_next, omega_next


def simulate_pendulum(initial_angle_deg,
                       length,
                       total_time=10.0,
                       dt=0.01,
                       gravity=9.81,
                       damping=0.0):
    """Simulate the large-angle pendulum using RK4.

    Returns time vector, angle (rad) vector, and angular velocity vector.
    """
    initial_angle_rad = math.radians(initial_angle_deg)
    state = (initial_angle_rad, 0.0)
    steps = int(total_time / dt)

    times = [0.0]
    angles = [state[0]]
    omegas = [state[1]]

    for i in range(steps):
        state = rk4_step(state, dt, length, gravity, damping)
        times.append((i + 1) * dt)
        angles.append(state[0])
        omegas.append(state[1])

    return times, angles, omegas

def analyse(time, angles):
    sp.curve_fit

def main(initial_theta):
    # Modifiable initial conditions
    initial_angle_deg = initial_theta  # degrees
    string_length = 0.36  # meters
    damping = 0.00  # linear damping coefficient (1/s); set 0 for no damping
    total_time = 10.0  # seconds
    dt = 0.005  # time step in seconds

    times, angles, omegas = simulate_pendulum(
        initial_angle_deg,
        string_length,
        total_time=total_time,
        dt=dt,
        damping=damping,
    )

    heights = [string_length * (1 - math.cos(angle)) for angle in angles]

    print(f"Initial angle: {initial_angle_deg} degrees")
    print(f"String length: {string_length} m")
    print(f"Linear damping: {damping} 1/s")
    print(f"Simulated {len(times)} steps with dt={dt}")
    print("First five time/angle pairs:")
    for t, theta in zip(times[:5], angles[:5]):
        print(f" t={t:.3f}s, theta={math.degrees(theta):.4f} degrees")

    print("\nLast five time/angle pairs:")
    for t, theta in zip(times[-5:], angles[-5:]):
        print(f" t={t:.3f}s, theta={math.degrees(theta):.4f} degrees")

    plt.figure(figsize=(8, 4))
    plt.plot(times, [math.degrees(a) for a in angles], label='Angle (deg)')
    plt.xlabel('Time (s)')
    plt.ylabel('Angle (degrees)')
    plt.title('Large-Angle Pendulum Simulation (RK4)')
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()

if __name__ == '__main__':
    main(125.0)
    main(100)
    main(50)
