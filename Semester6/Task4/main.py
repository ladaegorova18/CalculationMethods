import numpy as np
import math
from tabulate import tabulate

def H(x, y):
    return -0.5 * math.tan(x * y)

def f(x):
    return x - 0.5

a = 0
b = 1
half = (a + b) / 2
epsilon = 0.0001
N = 2 # Число узлов квадратурной формулы
m = -1
Ua = []
Uhalf = []
Ub = []

print("Число узлов квадратурной формулы Симпсона: {N}, A = {a}, B = {b}, epsilon = {epsilon})".format(N=N, a=a, b=b, epsilon=epsilon))

def loop(N, m):
    Ak = []
    X = []
    step = (b - a) / N
    for k in range(N):
        X.append(a + k * step)
        if k % 2 == 0:
            Ak.append(2 * (b - a) / (6 * N))
        else:
            Ak.append(4 * (b - a) / (6 * N))

    Ak[0] = (b - a) / (6 * N)
    Ak[N - 1] = (b - a) / (6 * N) # Ak - массив коэффициентов формулы Симпсона

    D = np.eye(N)
    for j in range(N):
        for k in range(N):
            D[j][k] = D[j][k] - Ak[k] * H(X[j], X[k])

    g = []
    for k in range(N):
        g.append(f(a + k * step))

    z = np.linalg.solve(D, g) # z - массив коэффициентов с для получения решения

    u_n_a = f(a)
    for k in range(N):
        u_n_a += Ak[k] * H(a, X[k]) * z[k]
    Ua.append(u_n_a) #Ua, Uhalf, Ub - решения в точках a, half, b при разных N

    u_n_half = f(half)
    for k in range(N):
        u_n_half += Ak[k] * H((a + b) / 2, X[k]) * z[k]
    Uhalf.append(u_n_half)

    u_n_b = f(b)
    for k in range(N):
        u_n_b += Ak[k] * H(b, X[k]) * z[k]
    Ub.append(u_n_b)
    m += 1
    N *= 2
    return N, m

N, m = loop(N, m)
N, m = loop(N, m)
while max(abs(Ua[m] - Ua[m - 1]), abs(Uhalf[m] - Uhalf[m - 1]), abs(Ub[m] - Ub[m - 1])) >= epsilon:
    N, m = loop(N, m)
Cond = []
for i in range(m + 1):
    Cond.append(["u^{} (x)".format(np.power(2, i + 1)), Ua[i], Uhalf[i], Ub[i]])
Cond.append(["u^{i1} (x) - u^{i2} (x)".format(i1 = np.power(2, i + 1), i2 = np.power(2, i)), Ua[m] - Ua[m - 1], Uhalf[m] - Uhalf[m - 1], Ub[m] - Ub[m - 1]])
Cond.append(["Итоговое решение:", Ua[m], Uhalf[m], Ub[m]])

print(tabulate(Cond, headers=["x", "a", "(a + b) / 2", "b"], tablefmt='pipe', numalign="right"))
print("Максимальная погрешность между шагами: {}".format(max(abs(Ua[m] - Ua[m - 1]), abs(Uhalf[m] - Uhalf[m - 1]), abs(Ub[m] - Ub[m - 1]))))