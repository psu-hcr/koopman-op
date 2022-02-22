import numpy as np
from numpy import genfromtxt
import matplotlib.pyplot as plt

data=genfromtxt('/home/lu/koopman-op/twolink.csv',delimiter=",",dtype=float)
data = np.delete(data,0,0)
time = data[:,0]
pos1 = data[:, 1]
pos2 = data[:, 3]
vel1 = data[:, 2]
vel2 = data[:, 4]
u1 = data[:, 5]
u2 = data[:, 6]
dpos1 = data[:, 7]
dpos2 = data[:, 9]
dvel1 = data[:, 8]
dvel2 = data[:, 10]

"""
mu1 = data[:,11]
mu2 = data[:,12]
"""

# error plot
plt.figure(1)
plt.title("theta 1")
plt.plot(time,pos1,'r')
plt.plot(time,dpos1,'g')

plt.figure(2)
plt.title("theta 2")
plt.plot(time,pos2,'r')
plt.plot(time,dpos2,'g')


# velocity plot
plt.figure(3)
plt.title("thetadot 1")
plt.plot(time,vel1,'b')
plt.plot(time,dvel1,'g')

plt.figure(4)
plt.title("thetadot 2")
plt.plot(time,vel2,'b')
plt.plot(time,dvel2,'g')


# input plot
plt.figure(5)
plt.title("u1")
plt.plot(time,u1,'m')

plt.figure(6)
plt.title("u2")
plt.plot(time,u2,'m')


# error plot
plt.figure(7)
plt.title("error pos1")
plt.plot(time,dpos1-pos1,'y')

plt.figure(8)
plt.title("error pos2")
plt.plot(time,dpos2-pos2,'y')


"""
# policy plot
plt.figure(7)
plt.title("mu1")
plt.plot(time,mu1,'y')

plt.figure(8)
plt.title("mu2")
plt.plot(time,mu2,'y')
"""

plt.show()
