import numpy as np
from numpy import genfromtxt
import matplotlib.pyplot as plt

data=genfromtxt('/home/lu/koopman-op/test.csv',delimiter=",",dtype=float)
data = np.delete(data,0,0)
tlist = data[0:-1,0]
x1 = (data[0:-1,1])
x2 = data[0:-1,2]
x3 = data[0:-1,3]
k1 = (data[0:-1,4])
k2 = (data[0:-1,5])
k3 = (data[0:-1,6])
xK = (data[0:-1,10])
'''
plt.figure()
#plt.plot(tlist,xK)#,'o',markersize=1)
plt.plot(tlist,x1,'g')#o',markersize=1)
plt.plot(tlist,x2,'r')
plt.plot(tlist,x3,'k')
plt.ylim(-15,15)
'''
plt.figure()
plt.plot(tlist,xK)#,'o',markersize=1)
plt.plot(tlist,k1,'g')#o',markersize=1)
plt.plot(tlist,k2,'r')
plt.plot(tlist,k3,'k')
plt.ylim(-15,15)

u1 = (data[0:-1,7])
u2 = (data[0:-1,8])
u3 = (data[0:-1,9])
plt.figure()

plt.plot(tlist,u3)
plt.plot(tlist,u1,'k')
plt.plot(tlist,u2,'r')
#plt.ylim(-5,5)

plt.show()
