import numpy as np
from tabulate import tabulate
import sympy
import math

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


def f(x_, t_):
    return sympy.lambdify([x, t], sympy.diff(u_symbol, t, 1))(x_, t_) - sympy.lambdify([x, t], sympy.diff(u_symbol, x, 2))(x_, t_) - sympy.lambdify([x, t], sympy.diff(u_symbol, x, 1))(x_, t_)


def a(x, t):
    return 1


def b(x, t):
    return 1


def c(x, t):
    return 0


def alpha1(t):
    return 1


def alpha2(t):
    return 0


def alpha(t):
    return u(0, t)


def phi(x):
    return u(x, 0)


def beta(t_):
    return sympy.lambdify([x, t], sympy.diff(u_symbol, x, 1))(1, t_)


def beta1(t):
    return 0


def beta2(t):
    return 1

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
def stability(J, explicit, sigma):
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

            U = explicit_difference_scheme(N, M, h, tau, X, T)
            U1 = explicit_difference_scheme(N1, M1, 2 * h, tau1, X1, T1)

            precision_J_ex = norm(J_ex, U)
            if h == 0.2:
                table.append([h, tau, precision_J_ex, '----------'])
            else:
                precision_U = norm(U, U1)
                table.append([h, tau, precision_J_ex, precision_U])

        else:   # для неявной схемы
            tau = 0.1 / M
            make_points()

            U = scheme_with_weights(sigma, N, M, h, tau, X)
            U1 = scheme_with_weights(sigma, N1, M, 2 * h, tau, X)
            J_ex = make_J_ex(J)

            precision_J_ex = norm(J_ex, U)
            if h == 0.2:
                table.append([h, tau, precision_J_ex, '----------'])
            else:
                precision_U = norm(U, U1)
                table.append([h, tau, precision_J_ex, precision_U])
    print(tabulate(table))


# Функция Lh(u_ik)
def Lh_uik(i, k, U, h, X, T):
    return a(X[i], T[k]) * (U[i + 1][k] - 2 * U[i][k] + U[i - 1][k]) / (h * h) + b(X[i], T[k]) * (U[i + 1][k] - U[i - 1][k]) / (2 * h) + c(X[i], T[k]) * U[i][k]

# Печать матрицы
def printMatrix(A):
   for i in range (len(A)):
      for j in range ( len(A[i]) ):
          print ( "{:4f}".format(A[i][j]), end = " " )
      print ()

# Печать таблицы решения на "крупной" сетке 6х6
def print_grid(U):
    table = []
    header = ['x/t']
    h = 1 / N
    tau = 0.1 / M
    hx = math.ceil(N / 5)
    ht = math.ceil(M / 5)
    Us = np.zeros((6, 6))
    for i in range(6):
        for k in range(6):
            Us[i][k] = U[i * hx, k * ht]

    for k in range(6):
        header.append(k * ht * tau)

    for i in range(6):
        line = []
        line.append(i * hx * h)
        for k in range(6):
            line.append(Us[i][k])
        table.append(line)
    print(tabulate(table, headers=header, tablefmt="grid"))

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

# Метод прогонки для решения системы (ABC) * U = G
def tridiagonal_matrix_algorithm(A, B, C, G, N):
    Y = np.zeros(N + 1)
    s = np.zeros(N + 1)
    t = np.zeros(N + 1)

    s[0] = C[0] / B[0]
    t[0] = -G[0] / B[0]

    for i in range(1, N + 1):
        s[i] = C[i] / (B[i] - A[i] * s[i - 1])
        t[i] = (A[i] * t[i - 1] - G[i]) / (B[i] - A[i] * s[i - 1])
    Y[N] = t[N]
    for i in reversed(range(N)):
        Y[i] = s[i] * Y[i + 1] + t[i]
    return Y

# Явная схема
def explicit_difference_scheme(N, M, h, tau, X, T):
    k = 0
    U = np.zeros((N + 1, M + 1))
    for i in range(N + 1):
        U[i][k] = phi(X[i]) # 1)

    for k in range (1, M + 1):
        for i in range(1, N):
            U[i][k] = U[i][k - 1] + tau * (Lh_uik(i, k - 1, U, h, X, T) + f(X[i], T[k - 1]))

        U[0][k] = (alpha(T[k]) * 2 * h + alpha2(T[k]) * (4 * U[1][k] - U[2][k])) / (2 * h * alpha1(T[k] + 3 * alpha2(T[k])))
        U[N][k] = (beta(T[k]) * 2 * h - beta2(T[k]) * (-4 * U[N - 1][k] + U[N - 2][k])) / (2 * h * beta1(T[k]) + 3 * beta2(T[k]))
    return U

# Схема с весами
def scheme_with_weights(sigma, N, M, h, tau, X):
    U1 = np.zeros((N + 1, M + 1))
    k = 0
    for i in range(N + 1):
        U1[i][k] = phi(X[i]) # 1)

    for k in range (1, M + 1):
        A = np.zeros(N + 1)
        B = np.zeros(N + 1)
        C = np.zeros(N + 1)
        Gk = np.zeros(N + 1)

        r = sigma / (h * h)
        for i in range(1, N):
            Gk[i] = - 1 / tau * U1[i, k - 1] - (1 - sigma) * Lh_uik(i, k - 1, U1, h, X, T) - f(X[i], T[k])
            A[i] = sigma * (1 / (h * h) - 1 / (2 * h))
            B[i] = -(- 2 * r - 1 / tau)
            C[i] = sigma * (1 / (h * h) + 1 / (2 * h))

        B[0] = -(alpha1(T[k]) * h + alpha2(T[k]))
        C[0] = -alpha2(T[k])
        Gk[0] = alpha(T[k]) * h

        A[N] = -beta2(T[k])
        B[N] = -(beta1(T[k]) * h + beta2(T[k]))
        Gk[N] = beta(T[k]) * h

        Uk = tridiagonal_matrix_algorithm(A, B, C, Gk, N)
        for i in range(N + 1):
            U1[i][k] = Uk[i]
    return U1


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
def algorithm():
    global N, M, h, tau
    # Решение на крупной сетке
    print("Явная разностная схема:")

    U_explicit = explicit_difference_scheme(N, M, h, tau, X, T)
    print_grid(U_explicit)

    print("Схема с весами")

    print("Sigma = 0")
    U_sigma0 = scheme_with_weights(0, N, M, h, tau, X)
    print_grid(U_sigma0)

    print("Sigma = 0.5")
    U_sigma_half = scheme_with_weights(0.5, N, M, h, tau, X)
    print_grid(U_sigma_half)

    print("Sigma = 1")
    U_sigma1 = scheme_with_weights(1, N, M, h, tau, X)
    print_grid(U_sigma1)

    # Таблицы устойчивости
    print("Устойчивость:")

    print("Явная разностная схема:")
    stability(u, True, -1)

    print("Схема с весами:")
    M = enter_M()

    print("Sigma = 0")
    stability(u, False, 0)

    print("Sigma = 0.5")
    stability(u, False, 0.5)

    print("Sigma = 1")
    stability(u, False, 1)
    refresh()


print("Вариант 3")
print("Уравнение вида:")
print("du/dt = (d^2u/dx^2) + (x^2 + 1) * du/dx + f(x, t)")
print("u(x,0) = phi(x), 0 <= x <= 1")
print("du/dx, (x = 0) = alpha(t), u(1,t) = beta(t), 0 <= t <= 0.1")

print("u(x, t) = x + t")
algorithm()

print("u(x,t) = 3 * x - 0.12 * t")
u_symbol = 3 * x - 0.12 * t
u = sympy.lambdify([x, t], 3 * x - 0.12 * t)
algorithm()

print("u(x,t) = x^3 + t^3")
u_symbol = sympy.Pow(x, 3) + sympy.Pow(t, 3)
u = sympy.lambdify([x, t], x * x * x + t * t * t)
algorithm()

print("u(x,t) = sin(2t + 1) + cos(2x)")
u_symbol = sympy.sin(2 * t + 1) + sympy.cos(2 * x)
u = sympy.lambdify([x, t], sympy.sin(2 * t + 1) + sympy.cos(2 * x))
algorithm()