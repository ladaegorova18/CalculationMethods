import numpy as np
from tabulate import tabulate

print("Вариант 1:")
print("Уравнение вида:")
print("du/dt = (d/dx)((x + 3) * (du/dx)) - xu + f(x,t)")
print("u(x,0) = phi(x), 0 <= x <= 1")
print("du/dx, (x = 0) = alpha(t), u(1,t) = beta(t), 0 <= t <= 0.1")


# L имеет вид б)
def p(x):
    return x + 3


def b(x, t):
    return 0


def c(x, t):
    return -x




print("Явная разностная схема:")

N = 10
A = p(max(0, 1))
h = 1 / N
tau = h * h / (2 * A)
M = 160

t = []
x = []
for i in range(N + 1):
    t.append(h * i)
    x.append(0.02 * i)

