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
y = np.array([2.79,  5.33, 7.67, 8.56, 9.866, 11.921, 13.88])
# y = np.array([57.475, 106.188, 156.7277, 177.107, 207.7649, 259.001, 310.35])
x = np.array([0.761, 1.503, 2.278, 2.538, 3.02, 3.774, 4.355])
x1 = x.reshape(-1, 1)
dat = {'x': x, 'y': y}
df = pd.DataFrame(dat)
c = y/x

mean = sum(c) / len(c)
var = sum((l-mean)**2 for l in c) / len(c)
st_dev = math.sqrt(var)
print(st_dev)
model = LinearRegression()
model.fit(x1, y)
predict = model.predict(x1[:7, :])
print('slope=', model.coef_)
sb.lmplot(x="x", y="y", data=df)
for x, y in zip(x, y):
    plt.text(x, y+0.1, '%.2f' % y, ha='center', va='bottom')
plt.plot(x1, predict, c="red")
plt.scatter(x, y)
plt.show()
