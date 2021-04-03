from numpy.linalg import eig, inv, norm, cond
import numpy as np
from tabulate import tabulate

myA=[
 [3.278164, 1.046583, -1.378574],
 [1.046583, 2.975937, 0.934251],
 [-1.378574, 0.934251, 4.836173]
]

myB = [ [-0.527466],
 [2.526877],
 [5.165441]]

# --- end of исходные данные


# --- вывод системы на экран
def FancyPrint(A, B, selected):
    for row in range(len(B)):
        print("(", end='')
        for col in range(len(A[row])):
             print("\t{1:10.2f}{0}".format(" " if (selected is None
or selected != (row, col)) else "*", A[row][col]), end='')
        print("\t) * (\tX{0}) = (\t{1:10.2f})".format(row + 1, B[row][0]))
# --- end of вывод системы на экран

# --- перемена местами двух строк системы
def SwapRows(A, B, row1, row2):
    A[row1], A[row2] = A[row2], A[row1]
    B[row1], B[row2] = B[row2], B[row1]
# --- end of перемена местами двух строк системы

# --- деление строки системы на число
def DivideRow(A, B, row, divider):
    A[row] = [a / divider for a in A[row]]
    B[row][0] /= divider
# --- end of деление строки системы на число

# --- сложение строки системы с другой строкой, умноженной на число
def CombineRows(A, B, row, source_row, weight):
    A[row] = [(a + k * weight) for a, k in zip(A[row], A[source_row])]
    B[row][0] += B[source_row][0] * weight
# --- end of сложение строки системы с другой строкой, умноженной начисло

# --- решение системы методом Гаусса (приведением к треугольному виду)
def Gauss(A, B):
    column = 0
    while (column < len(B)):
        current_row = None
        for r in range(column, len(A)):
            if current_row is None or abs(A[r][column]) > abs(A[current_row][column]):
                 current_row = r
        if current_row is None:
            return None
        if current_row != column:
            SwapRows(A, B, current_row, column)
        DivideRow(A, B, column, A[column][column])
        for r in range(column + 1, len(A)):
            CombineRows(A, B, r, column, -A[r][column])
        column += 1
    X = [0 for b in B]
    for i in range(len(B) - 1, -1, -1):
        X[i] = B[i][0] - sum(x * a for x, a in zip(X[(i + 1):], A[i][(i + 1):]))
    print("\n".join("X{0} =\t{1:10.2f}".format(i + 1, x) for i, x in
enumerate(X)))
    return X
# --- end of решение системы методом Гаусса (приведением к треугольному виду)
print("Исходная система:")
FancyPrint(myA, myB, None)
print("Решаем:")

x = Gauss(myA, myB)

def norminf (H):
    maxNorm = 0.0
    for i in range (0, 3):
        sum = 0.0;
        for j in range (0, 3):
            sum += abs(myA[i][j])
        maxNorm = max(maxNorm, sum)
    return maxNorm

def frobnormvector(X):
    sum = 0.0
    for i in range(len(X)):
        c = X[i]
        sum += abs(X[i])
    return sum ** 0.5

def frobnorm (A):
    sum = 0.0
    for i in range(len(A)):
        for j in range(len(A[0])):
            sum += abs(A[i][j])
    return sum ** 0.5

def aprprecision (k, H, g):
    f = frobnorm(H)
    return (f ** k) * frobnorm(np.zeros((3, 1))) + f ** k * frobnorm(g) / (1 - f)

def apostprecision (H, xk, xkplus1):
    return frobnorm(H) * frobnorm(xk - xkplus1) / (1 - frobnorm(H))

def lusthernik (H, xk, xkplus1):
    x_ = xk + (xkplus1 - xk) / (1 - max(abs(np.linalg.eigvals(H))))
    return xkplus1 + x_

def makezero():
    return np.zeros((3,1))

def subtract_vectors(v, w):
    return [vi - wi for vi, wi in zip(v, w)]

def factprecision(a):
    b = x
    for i in range(3):
        b[i] = b[i] - a[i]
    prec = frobnormvector(b)
    print("Фактическая погрешность: ", prec)
    return prec

def spectradius(H):
    return np.max(np.abs(np.linalg.eigvals(H)))

def toFixed(numObj, digits=0):
    return f"{numObj:.{digits}f}"

def simpleiter (k, H, g):
    xk = makezero()
    for i in range (k):
        xkplus1 = H @ xk + g
        xk = xkplus1
    print ("Решение: ", xkplus1)
    factprecision(xkplus1)
    print ("Априорная оценка: ", aprprecision(k, H, g))
    print("Апостериорная оценка: ", apostprecision(H, xk, xkplus1))
    l = lusthernik(H, xk, xkplus1)
    print("Решение, уточнённое по Люстернику: ", l)
    factprecision(l)
    return xkplus1

def zeidel (k, H, g):
    xk = makezero()
    xkplus1 = xk
    for p in range (k):
        for i in range (3):
            s = 0
            s1 = 0
            for j in range (i):
                s = H[i][j] * xkplus1[j][0]
            for j in range (i, 3):
                s1 = H[i][j] * xk[j]
            xkplus1[i] = s + s1 + g[i]
        xk = xkplus1
    factprecision(xkplus1)
    return xkplus1

def radius (H):
    Hl = np.tril(H, k=-1)
    Hr = np.triu(H, k=0)
    E = np.eye(3, 3, dtype=np.double)
    x = np.linalg.inv(E - Hl) @ Hr
    return spectradius(x)

def upperrelax(k, H, g):
    p = spectradius(H)
    q = 2 / (1 + np.sqrt(1 - p * p))
    xkplus1 = xk = makezero()
    for m in range(k):
        for i in range(3):
            sum = sum1 = 0
            for j in range(i - 1):
                sum += H[i][j] * xkplus1[j]
            for j in range(i + 1, 3):
                sum1 += H[i][j] * xk[j]
            xkplus1[i] = xk[i] + q * (sum + sum1 - xk[i] + g[i])
    print("Метод верхней релаксации: ", xkplus1)
    factprecision(xkplus1)
    return xkplus1

# 1)
x = Gauss(myA, myB)

# 2)
E = np.eye(3, 3, dtype=np.double)
Hd = E - inv(np.diag(np.diag(myA))) @ myA
gD = inv(np.diag(np.diag(myA))) @ myB
print("Норма H: ", norminf(Hd))

# 3)
print("Априорная погрешность для х(7): ", aprprecision(7, Hd, gD))

# 4)
xsimp = simpleiter(7, Hd, gD)
print("Метод простой итерации: ", xsimp)

# 5)
xseid = zeidel(7, Hd, gD)
print("Метод Зейделя: ", xseid)

print(xsimp - xseid)

# 6)
print("Радиус: ", radius(Hd))

# 7)
uppx = upperrelax(7, Hd, gD)

# 8)
Cond = []
Cond.append([factprecision(x), factprecision(xsimp), factprecision(xseid), factprecision(uppx)])

print(tabulate(Cond, headers=['Gauss', 'Simple iteration', 'Seidel', 'Successive over-relaxation'],
               tablefmt='pipe', numalign="right"))
