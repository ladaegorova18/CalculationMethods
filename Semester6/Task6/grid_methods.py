import numpy as np
import sympy


x, t = sympy.symbols('x t')

class GridMethods():
    def __init__(self, u, u_symbol):
        self.u_symbol = x + t
        self.u = sympy.lambdify([x, t], x + t)

    def f(self, x_, t_):
        return sympy.lambdify([x, t], sympy.diff(self.u_symbol, t, 1))(x_, t_) - sympy.lambdify([x, t], sympy.diff(self.u_symbol, x, 2))(x_, t_) - sympy.lambdify([x, t], sympy.diff(self.u_symbol, x, 1))(x_, t_)

    def a(self, x, t):
        return 1

    def b(self, x, t):
        return 0

    def c(self, x, t):
        return 0

    def alpha1(self, t):
        return 1

    def alpha2(self, t):
        return 0

    def alpha(self, t):
        return self.u(0, t)

    def phi(self, x):
        return self.u(x, 0)

    def beta(self, t_):
        return sympy.lambdify([x, t], sympy.diff(self.u_symbol, x, 1))(1, t_)

    def beta1(self, t):
        return 0

    def beta2(self, t):
        return 1

    # Функция Lh(u_ik)
    def Lh_uik(self, i, k, U, h, X, T):
        return self.a(X[i], T[k]) * (U[i + 1][k] - 2 * U[i][k] + U[i - 1][k]) / (h * h) + self.b(X[i], T[k]) * (U[i + 1][k] - U[i - 1][k]) / (2 * h) + self.c(X[i], T[k]) * U[i][k]


    # Метод прогонки для решения системы (ABC) * U = G
    def tridiagonal_matrix_algorithm(self, A, B, C, G, N):
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
    def explicit_difference_scheme(self, N, M, h, tau, X, T):
        k = 0
        U = np.zeros((N + 1, M + 1))
        for i in range(N + 1):
            U[i][k] = self.phi(X[i])  # 1)

        for k in range(1, M + 1):
            for i in range(1, N):
                U[i][k] = U[i][k - 1] + tau * (self.Lh_uik(i, k - 1, U, h, X, T) + self.f(X[i], T[k - 1]))

            U[0][k] = (self.alpha(T[k]) * 2 * h + self.alpha2(T[k]) * (4 * U[1][k] - U[2][k])) / (
                        2 * h * self.alpha1(T[k] + 3 * self.alpha2(T[k])))
            U[N][k] = (self.beta(T[k]) * 2 * h - self.beta2(T[k]) * (-4 * U[N - 1][k] + U[N - 2][k])) / (
                        2 * h * self.beta1(T[k]) + 3 * self.beta2(T[k]))
        return U

    # Схема с весами
    def scheme_with_weights(self, sigma, N, M, h, tau, X, T):
        U1 = np.zeros((N + 1, M + 1))
        k = 0
        for i in range(N + 1):
            U1[i][k] = self.phi(X[i])  # 1)

        for k in range(1, M + 1):
            A = np.zeros(N + 1)
            B = np.zeros(N + 1)
            C = np.zeros(N + 1)
            Gk = np.zeros(N + 1)

            r = sigma / (h * h)
            for i in range(1, N):
                Gk[i] = - 1 / tau * U1[i, k - 1] - (1 - sigma) * self.Lh_uik(i, k - 1, U1, h, X, T) - self.f(X[i], T[k])
                A[i] = sigma * (1 / (h * h) - 1 / (2 * h))
                B[i] = -(- 2 * r - 1 / tau)
                C[i] = sigma * (1 / (h * h) + 1 / (2 * h))

            B[0] = -(self.alpha1(T[k]) * h + self.alpha2(T[k]))
            C[0] = -self.alpha2(T[k])
            Gk[0] = self.alpha(T[k]) * h

            A[N] = -self.beta2(T[k])
            B[N] = -(self.beta1(T[k]) * h + self.beta2(T[k]))
            Gk[N] = self.beta(T[k]) * h

            Uk = self.tridiagonal_matrix_algorithm(A, B, C, Gk, N)
            for i in range(N + 1):
                U1[i][k] = Uk[i]
        return U1