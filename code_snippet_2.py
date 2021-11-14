import cvxpy as cp
import numpy as np

Tp = cp.Variable(pos=True)
PRF = cp.Variable(pos=True)
s1 = cp.Variable(pos=True)
s2 = cp.Variable(pos=True)
s3 = cp.Variable(pos=True)

# Define Parameters
c = 3e8  # Speed of light [m/s]
ft2m = 0.3048  # Feet to meter conversion
h = 25000 * ft2m  # Height [m]
b = 2e-6  # Receiver protection [s]
zeta = 30  # Look angle [deg]
theta_beam = 5  # Beamwidth [deg]
rn = h / np.cos(np.radians(zeta - theta_beam / 2))  # Near range [m]
rf = h / np.cos(np.radians(zeta + theta_beam / 2))  # Far range [m]
SW = rf - rn  # Swath width [m]
D_tx_min = 0.05  # Minimum transmit duty cycle
D_tx_max = 0.15  # Maximum transmit duty cycle
D_rx = 0.85  # Maximum receive duty cycle
BW = 500e6  # Bandwidth
BW_IF = 200e6  # IF Bandwidth
BW_dop = 2500  # Doppler support

objective_fn = s1 * s2 * s3
constraints = [
    Tp <= (2 * rn / c - b) * s1,
    Tp + b + 2 * rf / c <= s2 / PRF,
    Tp <= 2 * (rn - h) / c * s3,
    Tp * PRF >= D_tx_min,
    Tp * PRF <= D_tx_max,
    (Tp + 2 * (rf - rn) / c) * PRF <= D_rx,
    Tp >= 2 * (rf - rn) / c * BW / BW_IF,
    PRF >= BW_dop,
    s1 >= 1,
    s2 >= 1,
    s3 >= 1,
]
problem = cp.Problem(cp.Minimize(objective_fn), constraints)
problem.solve(gp=True)
print("Status: ", problem.status)
print("s1: ", s1.value)
print("s2: ", s2.value)
print("s3: ", s3.value)
print("Tp: ", Tp.value)
print("PRF: ", PRF.value)
print("Duty cycle: ", Tp.value * PRF.value)
print("rn_new: ", rn + c * (s3.value - 1) / (2 * s3.value) * Tp.value)
print("rn_new: ", (rn * s3.value) - (s3.value - 1) * h)
print("SW_ori: ", SW)
print("SW_new: ", rf - (rn + c * (s3.value - 1) / (2 * s3.value) * Tp.value))

# Status:           optimal
# s1:              1.000
# s2:              1.000
# s3:              1.144
# Tp:              7.405 us
# PRF:          7441.268 Hz
# Duty cycle:      5.510  %
