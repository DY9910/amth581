import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

# u_t = 5u_xx
# u(t, -1) = 0
# u(t, 2) = 0
# u(0, x) = sin((4pi / 3) * (x + 1))
# u(t, x) = e^(-5 * 16pi^2 / 9 *t) * sin(4pi / 3 * (x + 1))


dx_set = 0.125
Nx = int(3 // dx_set + 1)
x0 = -1
xf = 2
x = np.linspace(x0, xf, Nx)
dx = x[1] - x[0]
print("Nx = {}".format(Nx))
print("dx = {}".format(dx))

dt_set = 2 ** -10
Nt = int(0.25 // dt_set + 1)
t0 = 0
tf = 0.25
# Choose Nt so that dt / dx^2 = 1/4
# Nt = 2 * (Nx - 1) ** 2 + 1
t = np.linspace(t0, tf, Nt)
dt = t[1] - t[0]
print("Nt = {}".format(Nt))
print("dt = {}".format(dt))
# print("dt/dx^2 = {}".format(dt / dx ** 2))

U = np.zeros((Nx, Nt))
U[:, 0] = np.sin(np.pi * x)
U[0, :] = 0
U[-1, :] = 0

A = (np.diag(-2 * np.ones(Nx - 2)) + np.diag(np.ones(Nx - 3), 1) + np.diag(np.ones(Nx - 3), -1)) / dx ** 2

# Check the eigenvalues of A
vals, _ = np.linalg.eig(A)
k = np.arange(1, Nx - 1)
exact_vals = (2 / dx ** 2) * (np.cos(k * np.pi * dx) - 1)


def f(t, u):
    return A @ u

# Solve with Forward Euler
for k in range(Nt - 1):
    U[1:-1, (k + 1):(k + 2)] = U[1:-1, k:(k + 1)] + dt * f(t[k], U[1:-1, k:(k + 1)])

T, X = np.meshgrid(t, x)

ax = plt.axes(projection='3d')
ax.plot_surface(T, X, U)
plt.show()

def true_solution(t, x):
    return np.exp(-(np.pi ** 2) * t) * np.sin(np.pi * x)

err = U - true_solution(T, X)
ax = plt.axes(projection='3d')
ax.plot_surface(T, X, err)
plt.show()

global_err = np.max(np.abs(err))
print("Global Error = {}".format(global_err))