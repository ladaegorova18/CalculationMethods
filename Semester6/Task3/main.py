import numpy as np
from numpy.linalg import eig, inv, norm, cond
from tabulate import tabulate


A=[
[-1.00449, -0.38726, 0.59047],
[-0.38726, 0.73999, 0.12519],
[0.59047, 0.12519, -1.08660]
]

def frobnorm(A):
    sum = 0.0
    for i in range(len(A)):
        for j in range(len(A[0])):
            sum += (abs(A[i][j])) ** 2
    return sum ** 0.5

def powermethod(n, A, epsilon):
    Yk = np.ones((3, 1))
    maxl = maxk = 0
    l = l_ = 0
    for i in range(n):
        k = 0
        while (k >= 0):
            Ykplus1 = A @ Yk
            l_ = Ykplus1[i] / Yk[i]
            if abs(l - l_) < epsilon:
                if abs(l_) > maxl:
                    maxl = l_
                l = 0
                break
            l = l_
            Yk = Ykplus1
            maxk += 1
    Yk = Yk / frobnorm(Yk)
    return maxl, Yk, maxk

def dotproductmethod(n, A, epsilon):
    Yk = np.ones((3, 1))
    k = 0
    maxl = maxk = 0

    l = l_ = 0
    k = 0
    while (k >= 0):
        Ykplus1 = A @ Yk
        yk = np.squeeze(np.asarray(Yk))
        ykplus1 = np.squeeze(np.asarray(Ykplus1))
        l_ = np.dot(ykplus1, yk) / np.dot(yk, yk)
        if abs(l - l_) < epsilon:
            if abs(l_) > maxl:
                maxl = l_
            l = 0
            break
        l = l_
        Yk = Ykplus1
        maxk += 1
    Yk = Yk / frobnorm(Yk)
    return maxl, Yk, maxk

print("Вариант 7:")
eig = np.linalg.eig(A)
eigprecise = eig[0][0]
print("Максимальное по модулю собственное число матрицы А (точное знчение):")
print(eigprecise)

print("-------------------------------------------------------")

print("Степенной метод:")
epsilon = 0.001
eigpower = powermethod(3, A, epsilon)

print("Максимальное по модулю собственное число:")
print(eigpower[0][0])

print("Погрешность между найденным и точным значением:")
print(abs(eigprecise - eigpower[0][0]))

print("Соответствующий ему собственный вектор:")
print(eigpower[1])

print("Число итераций:")
print(eigpower[2])

print("-------------------------------------------------------")


print("Метод скалярных произведений:")
epsilon = 0.000001

eigdot = dotproductmethod(3, A, epsilon)
print("Максимальное по модулю собственное число:")
print(eigdot[0])

print("Погрешность между найденным и точным значением:")
print(abs(eigprecise - eigdot[0]))

print("Число итераций:")
print(eigdot[2])
