import numpy as np
from numpy import genfromtxt
import matplotlib.pyplot as plt

data=genfromtxt('/home/kzf5356/iiwa_ros/koopman-op/test.csv',delimiter=",",dtype=float)
data = np.delete(data,0,0)
tlist = data[0:-1,0]
x1 = (data[0:-1,1])
x2 = (data[0:-1,7])


plt.figure()
plt.plot(tlist,x2)#,'o',markersize=1)
plt.plot(tlist,x1,'k')#o',markersize=1)
plt.ylim(-20,5)

u1 = (data[0:-1,5])
u2 = (data[0:-1,6])
plt.figure()

plt.plot(tlist,u2)
plt.plot(tlist,u1,'k')
#plt.ylim(-5,5)

plt.show()