import numpy as np
from tabulate import tabulate
import sympy
import math

from printer_class import Printer as pr
from grid_methods import GridMethods

x, t = sympy.symbols('x t')

N = 10
M = 10
T_ = 0.1
A = 1

h = 1 / N
tau = h * h / (2 * A)
M = math.ceil(T_ / tau)
tau1 = 0

T = []
X = []
X1 = []
T1 = []

N1 = 0
M1 = 0

for i in range(N + 1):
    X.append(h * i)

for j in range(M + 1):
    T.append(tau * j)

u_symbol = x + t
u = sympy.lambdify([x, t], x + t)

# Создаем массив точек для вычисления устойчивости при заданных h и M
def make_points():
    global N, M, h, tau, tau1, X, T, X1, T1, N1, M1

    tau1 = 2 * h * 2 * h / (2 * A)
    X = []
    T = []
    X1 = []
    T1 = []

    M1 = math.ceil(1 / tau1)
    N1 = math.ceil(1 / (2 * h))

    for i in range(N + 1):
        X.append(h * i)
    for j in range(M + 1):
        T.append(tau * j)
    for k in range(N1 + 1):
        X1.append(2 * h * k)
    for l in range(M1 + 1):
        T1.append(tau1 * l)

# Матрица точного решения J
def make_J_ex(J):
    J_ex = np.zeros((N + 1, M + 1))
    for i in range(N + 1):
        for k in range(M + 1):
            J_ex[i][k] = J(i * h, k * tau)
    return J_ex

# Вычисление устойчивости решения
def stability(J, explicit, sigma, gridm):
    global N, h, tau, tau1, M
    table = []
    table.append(['h', 'tau', '||J_ex - u(h, tau)||', '||u(h, tau) - u(2h, tau1)||'])
    H = [0.2, 0.1, 0.05]
    for h in H:
        T_ = 0.1
        N = math.ceil(1 / h)

        if (explicit):  # для явной схемы
            tau = h * h / (2 * A)
            M = math.ceil(T_ / tau)

            make_points()

            J_ex = make_J_ex(J)

            U = gridm.explicit_difference_scheme(N, M, h, tau, X, T)
            U1 = gridm.explicit_difference_scheme(N1, M1, 2 * h, tau1, X1, T1)

            precision_J_ex = norm(J_ex, U)
            if h == 0.2:
                table.append([h, tau, precision_J_ex, '----------'])
            else:
                precision_U = norm(U, U1)
                table.append([h, tau, precision_J_ex, precision_U])

        else:   # для неявной схемы
            tau = 0.1 / M
            make_points()

            U = gridm.scheme_with_weights(sigma, N, M, h, tau, X, T)
            U1 = gridm.scheme_with_weights(sigma, N1, M, 2 * h, tau, X, T)
            J_ex = make_J_ex(J)

            precision_J_ex = norm(J_ex, U)
            if h == 0.2:
                table.append([h, tau, precision_J_ex, '----------'])
            else:
                precision_U = norm(U, U1)
                table.append([h, tau, precision_J_ex, precision_U])
    print(tabulate(table))


# Норма разности матриц u и v
def norm(u, v):
    max = 0
    row = min(len(u), len(v), 6)
    u_step_row = len(u) / row
    v_step_row = len(v) / row
    column = min(len(u[0]), len(v[0]), 6)
    u_step_column = len(u[0]) / column
    v_step_column = len(v[0]) / column
    for i in range (row):
        for j in range (column):
            r = abs(u[math.trunc(i * u_step_row)][math.trunc(j * u_step_column)] - v[math.trunc(i * v_step_row)][math.trunc(j * v_step_column)])
            if (r > max):
                max = r
    return max


def enter_M():
    print("Введите M или нажмите q, чтобы оставить M = 10")
    flag = True
    M = 0

    while (flag):
        try:
            m = input()
            if (m == 'q'):
                M = 10
                flag = False
            else:
                M = int(m)
                flag = False
        except ValueError:
            print("M должен быть числом")
            pass
    return M

# Обновление глобальных переменных для новой функции
def refresh():
    global X, T, N, M, h, tau, N1, M1, X1, T1
    N = 10
    M = 10
    T_ = 0.1
    A = 1

    h = 1 / N
    tau = h * h / (2 * A)
    M = math.ceil(T_ / tau)

    T = []
    X = []
    X1 = []
    T1 = []

    N1 = 0
    M1 = 0

    for i in range(N + 1):
        X.append(h * i)

    for j in range(M + 1):
        T.append(tau * j)

# Применение методов для каждой функции
def process(u, u_symbol):
    global N, M, h, tau
    # Решение на крупной сетке
    print("Явная разностная схема:")
    gridm = GridMethods(u, u_symbol)
    U_explicit = gridm.explicit_difference_scheme(N, M, h, tau, X, T)
    pr.print_grid(U_explicit, N, M)

    print("Схема с весами")

    print("Sigma = 0")
    U_sigma0 = gridm.scheme_with_weights(0, N, M, h, tau, X, T)
    pr.print_grid(U_sigma0, N, M)

    print("Sigma = 0.5")
    U_sigma_half = gridm.scheme_with_weights(0.5, N, M, h, tau, X, T)
    pr.print_grid(U_sigma_half, N, M)

    print("Sigma = 1")
    U_sigma1 = gridm.scheme_with_weights(1, N, M, h, tau, X, T)
    pr.print_grid(U_sigma1, N, M)

    # Таблицы устойчивости
    print("Устойчивость:")

    print("Явная разностная схема:")
    stability(u, True, -1, gridm)

    print("Схема с весами:")
    M = enter_M()

    print("Sigma = 0")
    stability(u, False, 0, gridm)

    print("Sigma = 0.5")
    stability(u, False, 0.5, gridm)

    print("Sigma = 1")
    stability(u, False, 1, gridm)
    refresh()


print("Вариант 17")
print("Уравнение вида:")
print("du/dt = (d^2u/dx^2) + f(x, t)")
print("u(x,0) = phi(x), 0 <= x <= 1")
print("u(0, t) = alpha(t), u(1,t) = beta(t), 0 <= t <= 0.1")

print("u(x, t) = x + t")
process(u, u_symbol)

print("u(x,t) = 6 * x - 0.1 * t")
u_symbol = 6 * x - 0.1 * t
u = sympy.lambdify([x, t], 6 * x - 0.1 * t)
process(u, u_symbol)

print("u(x,t) = sin(1 - t) + cos(x)")
u_symbol = sympy.sin(1 - t) + sympy.cos(x)
u = sympy.lambdify([x, t], sympy.sin(1 - t) + sympy.cos(x))
process(u, u_symbol)