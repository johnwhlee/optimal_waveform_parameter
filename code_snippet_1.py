import time
import cvxpy as cp
import numpy as np

Tp = cp.Variable(pos=True)
PRF = cp.Variable(pos=True)

# Define Parameters
c = 3e8  # Speed of light [m/s]
ft2m = 0.3048  # Feet to meter conversion
h = 25000 * ft2m  # Height [m]
b = 2e-6  # Receiver protection [s]
zeta = 50  # Look angle [deg]
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

objective_fn = Tp * PRF
constraints = [
    Tp + b <= 2 * rn / c,
    Tp + b + 2 * rf / c <= 1 / PRF,
    Tp <= 2 * (rn - h) / c,
    Tp * PRF >= D_tx_min,
    Tp * PRF <= D_tx_max,
    (Tp + 2 * (rf - rn) / c) * PRF <= D_rx,
    Tp >= 2 * (rf - rn) / c * BW / BW_IF,
    PRF >= BW_dop,
]
problem = cp.Problem(cp.Maximize(objective_fn), constraints)
start_time = time.time()
problem.solve(gp=True)
print((time.time() - start_time))
print("Status: ", problem.status)
print("Tp: ", Tp.value)
print("PRF: ", PRF.value)
print("Duty cycle: ", Tp.value * PRF.value)

# Status:           optimal
# Tp:             23.070 us
# PRF:          6502.003 Hz
# Duty cycle:     15.000  %
# Solve Time:      0.104  s
