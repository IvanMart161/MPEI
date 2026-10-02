import numpy as np 
# import matplotlib.pyplot as plt 

N = 1000 

K = 100000
D = 3
mu = 1
M = np.zeros(K)
Var = np.zeros(K)

for i in range(K): 
    x = D ** 0.5 * np.random.randn(N) + mu
    M[i] = np.sum(x) / N
    Var[i] = np.sum((x - M[i])**2) / N

MM = np.sum(Var) / K
DD = np.sum((Var - MM)**2) / (K-1)
 
print('Усредненная оценка дисперсии: %.5f' % (np.sum(Var)/K))
print('Дисперсия оценки дисперсии: %.5f' % DD)
