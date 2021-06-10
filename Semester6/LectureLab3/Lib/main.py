import math
import numpy as np
from tabulate import tabulate
from numpy.linalg import solve


#Вычисленный u(x) при z(s)=1
def u(x):
    return (math.exp(3*x) - 1) / 3 * x


#Ядро:
def K(x, s):
    return math.exp(3 * x * s)


#Вычисление узлов для формулы средних прямоугольников:
def build_nodes(N):
    result = []
    for i in range(N):
        result.append(i/N+1/(2*N))
    return result


#Построение матрицы С:
def build_matrix(N):
    result = np.zeros((N,N))
    nodes = build_nodes(N)
    for i in range(N):
        for j in range(N):
            result[i,j] = K(nodes[i],nodes[j])*(1/N)
    return result


#Построение значений u(x) в узлах квадратурной формулы средних прямоугольников при z(s)=1
def build_right_part(N):
    result = np.zeros((N, 1))
    nodes = build_nodes(N)
    for i in range(N):
            result[i] = u(nodes[i])
    return result


#Решение СЛАУ методом регуляризации матрицы c параметром а:
def build_solution(N, a):
    C = build_matrix(N)
    right = build_right_part(N)
    return solve(np.transpose(C)@C+a*np.identity(N), np.transpose(C)@right)


#Рассмотрим погрешность при различных параметрах N (число узлов) и a (параметр метода регуляризации) и запишем результат в файл в виде таблицы.
N = 10
ans = np.zeros((N,N+1))
for i in range(1,N+1):
    ans[i-1,0] = i*10
headers = [r'n\a']
for j in range(6,6+N):
    headers.append(10**(-j))
residuals = []
for i in range(1,1+N):
    u0 = np.ones((N*i,1))
    residuals.append(N*i)
    for j in range(6,6+N):
        solution = build_solution(N*i,10**(-j))
        ans[i-1,j-5] = abs(u0-solution).max()

res_table = tabulate(ans,headers=headers,
               tablefmt='github',numalign="right")
with open('result3.txt','w') as file:
    file.write(res_table)


#Посмотрим на индекс минимального элемента таблицы. Таким образом, наименьшая погрешность достигается при числе узлов N=100 и параметре а=1e-10, и именно они являются оптимальными параметрами для данной задачи.
ind = np.unravel_index(np.argmin(ans, axis=None), ans.shape)
ind
#Сама погрешность:

ans[ind]