# -*- coding: cp1251 -*- 
#для работы кириллицы
import matplotlib.pyplot as plt  #для графиков переименовываем в plt для удобства
import numpy as np #работа с массивами данных
#from mpl_toolkits.mplot3d import Axes3D
#from Matplot3DEx import Matplot3DEx

x = np.genfromtxt('file1.txt')
t = np.genfromtxt('file2.txt')
t_n=np.genfromtxt('file7.txt')
t1,x1 = np.meshgrid(t,x)
t2,x2 = np.meshgrid(t_n,x)
u = np.genfromtxt("file3.txt", delimiter=" ") 
u2 = np.genfromtxt("file5.txt", delimiter=" ") 

u= u.transpose()
u2= u2.transpose()
#ax1.plot_surface(t, x, u,cmap='plasma')
#ax2.plot_surface(t, x, u2,cmap='plasma')

fig = plt.figure(figsize = (10,4))
ax1 = fig.add_subplot(121,projection="3d")
ax2 = fig.add_subplot(122,projection="3d")
ax1.set_title('явная схема',fontsize=15,fontweight="bold")
ax2.set_title('неявная схема',fontsize=15,fontweight="bold")
ax1.set(xlabel="$t$", ylabel="$x$", zlabel="$u$");
ax2.set(xlabel="$t$", ylabel="$x$", zlabel="$u$");
ax1.plot_surface(t1, x1, u,cmap='plasma')
ax2.plot_surface(t2, x2, u2,cmap='plasma')

plt.show()