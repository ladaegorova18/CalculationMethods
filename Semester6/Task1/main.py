import numpy as np


def countDelta(A, A_):
    dA = A_
    for i in range(0, len(dA)):
        for j in range(0, len(dA[0])):
            dA[i][j] = dA[i][j] - A[i][j]
    return dA


# multiply matrix to vector

def countX(M, b):
    x = [[M[0][0] * b[0][0] + M[0][1] * b[1][0]], [M[1][0] * b[0][0] + M[1][1] * b[1][0]]]
    print(x[0][0])
    print(x[1][0])
    return x


# count matrix norm (p = 2)

def norm(A):
    sum = 0
    for i in range(0, len(A)):
        for j in range(0, len(A[0])):
            sum += pow(A[i][j], 2)
    return pow(sum, 0.5)


def printMatrix(A, B, selected):
    for row in range(len(B)):
        print("(", end='')
        for col in range(len(A[row])):
            print("\t{1:10.2f}{0}".format(" " if (selected is None
                                                  or selected != (row, col)) else "*", A[row][col]), end='')
        print("\t) * (\tX{0}) = (\t{1:10.2f})".format(row + 1, B[row][0]))


# Part 1
print("Решение линейной алгебраической системы, нахождение числа обусловленности")
print("Вариант 7")

A = [[-401.43, 200.19], [1201.14, -601.62]]
b = [[200], [-600]]
b_ = [[199], [-601]]

print("A|b:")
printMatrix(A, b, None)

print("A|b_:")
printMatrix(A, b_, None)

detA = A[0][0] * A[1][1] - A[0][1] * A[1][0]
M = [[A[1][1], -A[0][1]], [-A[1][0], A[0][0]]]

for i in range(0, 2):
    for j in range(0, 2):
        M[i][j] = M[i][j] * 1.0 / detA

print("Решение системы с точной правой частью x:")
x = countX(M, b)

print("Решение системы с возмущённой правой частью x_:")
x_ = countX(M, b_)  # M = A ^ (-1) Ax = b x = A ^ (-1) * B

cond = norm(M) * norm(A)  # count cond(A)

print(f"Число обусловленности cond: {cond}")

x_x = [[x[0][0] - x_[0][0]], [x[1][0] - x_[1][0]]]

deltaX = norm(x_x) / norm(x)

# count norm of matrices | A - A_ |, | b - b_ |

print(f"Фактическая относительная погрешность deltaX: {deltaX}")
dA = [[0, 0], [0, 0]]
db = countDelta(b, b_)

R = cond * (norm(db) / norm(b) + norm(dA) / norm(A)) / (1 - cond * norm(dA) / norm(A))
print(f"Оценка погрешности R: {R}")

# Part 2

epsilon = 0.00001

def multiply(revA, A):
    n = len(A)
    result = np.zeros((n, n))
    for i in range(0, n):
        for j in range(0, n):
            sum = 0
            for k in range(0, n):
                sum += revA[i][k] * A[k][j]
            result[i][j] = sum
    return result


def printPrecision(A, X):
    R = []
    Ax = []
    b = []
    n = len(A)

    for i in range(0, len(A)):
        sum = 0
        for j in range(0, len(A)):
            sum += A[i][j] * X[j]

        Ax.append(sum)
        b.append(A[i][n])

    for j in range(0, n):
        R.append(Ax[j] - b[j])
    print("Погрешность:")
    printResult(R, "R")


def printResult(X, name):
    print(f"{name}:")
    for i in range(0, len(X)):
        print(X[i])
    print()


def GaussMethod(A):
    A = ForvardMove(A, False)
    X = ReverseMove(A)
    return X


def gaussWithMainElement(A):
    A1 = np.copy(A)
    D = ForvardMove(A1, True)
    X = ReverseMove(D)
    return X


def ForvardMove(D, chooseMainElement):
    n = len(D)
    for k in range(0, n):
        tmp = abs(D[k][k])
        if (tmp < epsilon):
            if (chooseMainElement):
                D = chooseMain(D, k)
            else:
                print(f"Слишком маленький ведущий элемент {D[k][k]}")
        for j in range(k, n + 1):
            D[k][j] = D[k][j] / tmp

        for i in range(k + 1, n):
            tmp = D[i][k]
            for j in range(k, n + 1):
                D[i][j] = D[i][j] - D[k][j] * tmp
    return D


def ReverseMove(A):
    n = len(A[0]) - 1
    X = np.zeros(len(A))
    for i in range(n - 1, -1, -1):
        sum = 0
        for j in range(i + 1, n):
            sum += A[i][j] * X[j]

        X[i] = A[i][n] - sum
    return X


def LUdecomposition(n, A):
    L = np.zeros((n, n))
    U = np.zeros((n, n))

    for i in range(0, n):
        for j in range(i, n):
            sum = 0
            for k in range(0, i):
                sum += L[j][k] * U[k][i]
            L[j][i] = A[j][i] - sum

            sum = 0
            for k in range(0, i):
                sum += L[i][k] * U[k][j]
            U[i][j] = (A[i][j] - sum) / L[i][i]
    det = 1
    for i in range(0, n):
        det *= L[i][i]
    return det


def Jordan(A):
    n = len(A)
    for k in range(0, n):
        tmp = abs(A[k][k])
        if (tmp < epsilon):
            print(f"Слишком маленький ведущий элемент {A[k][k]}")

        for j in range(k, len(A[0])):
            A[k][j] = A[k][j] / tmp

        for i in range(0, n):
            tmp = A[i][k]
            if (i != k):
                for j in range(k, len(A[0])):
                    A[i][j] = A[i][j] - A[k][j] * tmp
    return A


def reverseMatrix(A, n):
    reverseA = np.zeros((n, (len(A[0]) - 1) * 2))
    for i in range(0, n):
        for j in range(0, n):
            reverseA[i][j] = A[i][j]
            if (i == j):
                reverseA[i][j + n] = 1
    reverseA = Jordan(reverseA)
    result = np.zeros((n, n))
    for i in range(0, n):
        for j in range(0, n):
            result[i][j] = reverseA[i][j + n]

    return result


def chooseMain(A, k):
    p = 0
    max = abs(A[k][k])

    for i in range(k, len(A)):
        if (abs(A[i][k]) > max):
            max = abs(A[i][k])
            p = i

    tmp = A[k]
    A[k] = A[p]
    A[p] = tmp
    return A


A_ = np.array([[9.331343, 1.120045, -2.880925], [1.120045, 7.086042, 0.670297],
     [-2.880925, 0.670297, 5.622534]])
b_ = np.array([[7.570463], [8.876384], [3.411906]])
print("A:")
printMatrix(A_, b_, None)

A = np.concatenate([A_, b_], axis=1)
print("Метод Гаусса с выбором главного элемента для матрицы А:\n")

resultA = gaussWithMainElement(A)
print("Матрица A:")
printResult(resultA, "X")
printPrecision(A, resultA)

print("Метод LU-разложения для нахождения определителя матрицы А:\n")

det = LUdecomposition(len(A), A)
print(f"det A: {det}")

print("Обратная матрица для A:")
reverseA = reverseMatrix(A, len(A))
print(reverseA)

print("A * A^(-1)")
print(multiply(reverseA, A))
