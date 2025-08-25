# -*- coding: cp1251 -*- 
#для работы кириллицы
import matplotlib.pyplot as plt  #для графиков переименовываем в plt для удобства
import numpy as np #работа с массивами данных
x = np.genfromtxt('file1.txt')
yP = np.genfromtxt('file2.txt')
yS = np.genfromtxt('file3.txt')
yN = np.genfromtxt('file4.txt')
ySek = np.genfromtxt('file5.txt')


def cm_to_inch(value): #функция для того,чтобы можно было указывать размер графика
    return value/2.54
fig=plt.figure(figsize=(cm_to_inch(40),cm_to_inch(20))) #укажем размер графика
ax1=fig.add_subplot(121)
ax2=fig.add_subplot(122)
ax1.set_title('Линейная задача',fontsize=15,fontweight="bold")
ax2.set_title('Нелиненйая задача',fontsize=15,fontweight="bold")
ax1.plot(x,yP,linewidth=5,color='red',label="прогонка")
ax1.plot(x,yS,linewidth=3,color='yellow',label="стрельба")
ax1.legend() 
ax1.grid()
ax2.plot(x,yN,linewidth=5,color='red',label="Ньютон")
ax2.plot(x,ySek,linewidth=3,color='yellow',label="секущих")
ax2.legend() 
ax2.grid()
plt.show()
