# -*- coding: cp1251 -*- 
#для работы кириллицы
import matplotlib.pyplot as plt  #для графиков переименовываем в plt для удобства
import numpy as np #работа с массивами данных
def f1(x):
    return np.log(x+2)
def f2(y):
    return y**3
# Создатётся массив целых чисел со значениями начала, конца и шага
x=np.arange(-1.99,120.01,0.01) #беру от -1.99 тк лог(x+2) 
y=np.arange(-3,5.01,0.01)
#fig=plt.figure(figsize=(42,8))
def cm_to_inch(value): #функция для того,чтобы можно было указывать размер графика
    return value/2.54
fig=plt.figure(figsize=(cm_to_inch(40),cm_to_inch(20))) #укажем размер графика


ax1=fig.add_subplot(131)
ax2=fig.add_subplot(132)
ax3=fig.add_subplot(133)
ax1.set_title('Графики функций',fontsize=15,fontweight="bold")
ax2.set_title('Графики функций',fontsize=15,fontweight="bold")
ax3.set_title('Линии уровня суммы квадратов невязок',fontsize=15,fontweight="bold")
#ax1.set_grid(True)#сетка 
ax1.plot(x,f1(x),color='red',label="ln(x+2)-y=0")
ax1.plot(f2(y),y,color='blue',label="x-y^3=0")
ax2.plot(x,f1(x),color='red',label="ln(x+2)-y=0")
ax2.plot(f2(y),y,color='blue',label="x-y^3=0")
ax1.set_xlim([-100,120])
ax1.set_ylim([-3,5])
ax2.set_xlim([-10,10])
ax2.set_ylim([-3,5])
ax1.legend() 
x = np.arange(-1.999,200, 0.01)
y = np.arange(-5, 10, 0.01)
X, Y=np.meshgrid(x,y)
Z=(X-Y**3)**2+(np.log(X+2)-Y)**2
ax3.set_xlim([-10,60])
ax3.set_ylim([-2.5,6])
ax3.contour(X, Y, Z,[0.4,1.5,10,30])# в [] заданы константы приравнивающиеся к Z
ax1.grid()
ax2.grid()
ax3.grid()
plt.show()