import numpy as np
from numpy import genfromtxt
import matplotlib.pyplot as plt

data=genfromtxt('/home/lu/koopman-op/twomass.csv',delimiter=",",dtype=float)
data = np.delete(data,0,0)
time = data[:,0]
x1 = data[:, 1]
x2 = data[:, 2]
x3 = data[:, 3]
x4 = data[:, 4]
u1 = data[:, 5]
u2 = data[:, 6]
mu1 = data[:, 7]
mu2 = data[:, 8]


plt.figure(1)
plt.plot(time,x1,'r')
plt.plot(time,x2,'g')
plt.plot(time,x3,'b')
plt.plot(time,x4,'k')

plt.figure(2)
plt.plot(time,u1,'r')
plt.plot(time,mu1,'g')
plt.title("u1 & mu1")

plt.figure(3)
plt.plot(time,u2,'r')
plt.plot(time,mu2,'g')
plt.title("u2 & mu2")


plt.show()
