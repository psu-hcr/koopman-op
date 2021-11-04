import numpy as np
from numpy import genfromtxt
import matplotlib.pyplot as plt

data=genfromtxt('/home/kzf5356/iiwa_ros/koopman-op/test.csv',delimiter=",",dtype=float)
data = np.delete(data,0,0)
tlist = data[0:-1,0]
x1 = (data[0:-1,1])
x2 = data[0:-1,2]
x3 = data[0:-1,3]
xK = (data[0:-1,7])


plt.figure()
plt.plot(tlist,xK)#,'o',markersize=1)
plt.plot(tlist,x1,'g')#o',markersize=1)
plt.plot(tlist,x2,'r')
plt.plot(tlist,x3,'k')
plt.ylim(-10,10)

u1 = (data[0:-1,4])
u2 = (data[0:-1,5])
u3 = (data[0:-1,6])
plt.figure()

plt.plot(tlist,u2)
plt.plot(tlist,u1,'k')
plt.plot(tlist,u1,'r')
#plt.ylim(-5,5)

plt.show()
