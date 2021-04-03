from numpy.linalg import eig, inv, norm, cond
import numpy as np
from tabulate import tabulate


def genVandermonde(n):
    result = np.zeros(shape=(n, n))
    for i in range(n):
        for k in range(n):
            result[i, k] = (i + 1) ** ((n + 1 - (k + 1)) ** (-4))
    return result


def newton(X, epsilon):
    x_k = np.identity(X.shape[0])
    x_k1 = 0.5 * (x_k + inv(x_k) @ X)
    while norm(x_k1 - x_k) > epsilon:
        x_k = x_k1
        x_k1 = 0.5 * (x_k + inv(x_k) @ X)
    return x_k1


def eigMethod(X):
    V = eig(X)[1]
    sigma = np.diag(eig(X)[0]**(0.5))
    return V @ sigma @ inv(V)


def toFixed(numObj, digits=0):
    return f"{numObj:.{digits}f}"


Cond = []

for i in range(2, 16):
    X = genVandermonde(i)
    eigRoot = eigMethod(X)
    if i in range(2, 8):
        newtonRoot = newton(X, 1e-3 * (i ** 2))
    else:
        newtonRoot = newton(X, 1e-2 * (i ** 2))
    Cond.append([i, toFixed(cond(X), 8), toFixed(cond(newtonRoot), 8), toFixed(cond(eigRoot), 8),
                 toFixed(norm(newtonRoot @ newtonRoot - X), 8), toFixed(norm(eigRoot @ eigRoot - X), 8)])

print(tabulate(Cond, headers=['n', 'cond (A)', 'cond (B) (Newton)', 'cond (B) (Eig)',
                              'norm (B^2 - A) (newton)', 'norm (B^2 - A) (eig)'],
               tablefmt='pipe', numalign="right"))


X = genVandermonde(3)
print(X)
root = eigMethod(X)
print(root)
print(root @ root)