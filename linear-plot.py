import numpy as np
from numpy import genfromtxt
import matplotlib.pyplot as plt

data=genfromtxt('/home/lu/koopman-op/linear.csv',delimiter=",",dtype=float)
data = np.delete(data,0,0)
time = data[:,0]
x1 = data[:, 1]
x2 = data[:, 2]
u1 = data[:, 3]
u2 = data[:, 4]



plt.figure(1)
plt.plot(time,x1,'r')
plt.plot(time,x2,'g')

plt.figure(2)
plt.plot(time,u1,'r')
plt.plot(time,u2,'g')


plt.show()
