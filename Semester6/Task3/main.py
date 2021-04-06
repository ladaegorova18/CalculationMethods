from numpy.linalg import eig, inv, norm, cond
import numpy as np
from tabulate import tabulate

myA=[
 [-0.81417, -0.01937, 0.41372],
 [-0.01937, 0.54414, 0.00590],
 [0.41372, 0.00590, -0.81445]
]


def Jacobi(n, A, epsilon):
 E = np.eye(n)
 C = B = np.ones((n, n))
 Akplus1 = Ak = A
 ik = jk = 0
 diagA = np.triu(Ak, k=0)
 maxel = 0
 X = np.zeros((n, n))
 for i in range(n):
  for j in range(n):
   if abs(Ak[i][j]) > maxel:
    maxel = Ak[i][j]
    ik = i
    jk = j
   if maxel >= epsilon:
    Akplus1 = rotate(Ak, Akplus1, ik, jk, n)
  Ak = Akplus1
 eigval = np.diag(Akplus1)
 eigvect = C
 return eigval, eigvect


def rotate(Ak, Akplus1, ik, jk, n):
 d = np.sqrt((Ak[ik][ik] - Ak[jk][jk]) ** 2 + 4 * (Ak[ik][jk] ** 2))
 c = np.sqrt(0.5 * (1 + abs(Ak[ik][ik] - Ak[jk][jk]) / d))
 s = np.sign(Ak[ik][jk] * (Ak[ik][ik] - Ak[jk][jk])) * np.sqrt(0.5 * (1 - abs(Ak[ik][ik] - Ak[jk][jk]) / d))
 Akplus1[ik][ik] = c ** 2 * Ak[ik][ik] - 2 * c * s * Ak[ik][jk] + s ** 2 * Ak[jk][jk]
 Akplus1[jk][jk] = s ** 2 * Ak[ik][ik] - 2 * c * s * Ak[ik][jk] + c ** 2 * Ak[jk][jk]
 Akplus1[ik][jk] = Akplus1[jk][ik] = 0
 for l in range(n):
  for m in range(n):
   if l != ik & l != jk & m != ik & m != jk:
    Akplus1[l][m] = Ak[l][m]
   elif l != ik & l != jk:
    Akplus1[l][ik] = Akplus1[ik][l] = c * Ak[l][ik] + s * Ak[l][jk]
    Akplus1[l][jk] = Akplus1[jk][l] = -s * Ak[l][ik] + c * Ak[l][jk]
 return Akplus1


def frobnorm(A):
 sum = 0.0
 for i in range(len(A)):
  for j in range(len(A[0])):
   sum += abs(A[i][j])
 return sum ** 0.5


def aposter(A, Y, l):
 return frobnorm(A @ Y - l * Y) / frobnorm(Y)


def powermethod(n, A, epsilon):
 Yk = [1, 1, 1]
 maxl = maxk = 0
 l = l_ = 0
 for i in range(n):
     k = 0
     while (k >= 0):
         Ykplus1 = A @ Yk
         l_ = Ykplus1[i] / Yk[i]
         if abs(l - l_) < epsilon & abs(l_) > maxl:
             maxl = l_
             maxk = k
             k = 0
             break
         l = l_
         k += 1
     Yk = Ykplus1
 return maxl, maxk, Ykplus1

def dotproductmethod():
 Y0 = np.ones((3, 1))



print("Вариант 1:")
# 1)
epsilon = 0.000001
eig = Jacobi(3, myA, epsilon)
#eig1 = np.linalg.eig(myA)
print(eig[0])
#rint(eig1[0])

# 2)
epsilon = 0.001
