import matplotlib.pyplot as plt
import numpy as np
from qutip import *
import math
from sympy import *
from sklearn.linear_model import LinearRegression
from scipy import stats
import pandas as pd
import seaborn as sb


c = []
y = np.array([6.25, 3.1733, 2.097, 1.586, 1.2693,
             1.0475, 0.888, 0.7617, 0.634, 0.311, 0.2213, 0.1579, 0.1263, 0.095])  # square
y2 = np.array([14.18, 7.12, 4.76, 3.51, 2.856,
              2.26, 1.93, 1.74, 1.4, 0.631, 0.471, 0.3174, 0.308, 0.283])  # gaussian
y3 = np.array([0.46, 0.2378, 0.1586, 0.1110, 0.0952,
              0.0793, 0.0634, 0.0627, 0.0476, 0.0222, 0.01587, 0.01267, 0.01, 0.008568])
x = np.array([0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 5, 10, 15, 20, 25, 30])
x1 = x.reshape(-1, 1)
dat = {'x': x, 'y': y}
df = pd.DataFrame(dat)
model1 = LinearRegression()
model1.fit(x1, 1/y)
predict1 = model1.predict(x1[:14, :])
print('slope=', model1.coef_)
model2 = LinearRegression()
model2.fit(x1, 1/y2)
predict2 = model2.predict(x1[:14, :])
print('slope=', model2.coef_)
model3 = LinearRegression()
model3.fit(x1, 1/y3)
predict3 = model3.predict(x1[:14, :])
print('slope=', model3.coef_)
plt.plot(x, predict3, label='exponential')
plt.plot(x, predict1, label='square_wave')
plt.plot(x, predict2, label='Gaussian')
plt.title('pulse time to 1/pulse amplitude')
plt.ylabel('1/pulse amplitude')
plt.xlabel('pulse time')
plt.legend(loc='right')
plt.scatter(x, 1/y3)
plt.scatter(x, 1/y)
plt.scatter(x, 1/y2)
plt.show()
