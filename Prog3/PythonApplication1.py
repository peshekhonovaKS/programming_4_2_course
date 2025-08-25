# -*- coding: cp1251 -*- 
#��� ������ ���������
import matplotlib.pyplot as plt  #��� �������� ��������������� � plt ��� ��������
import numpy as np #������ � ��������� ������
x = np.genfromtxt('file1.txt')
yE1 = np.genfromtxt('file2.txt')
dyE1 = np.genfromtxt('file3.txt')
yE2 = np.genfromtxt('file4.txt')
dyE2=np.genfromtxt('file5.txt')
yR2 = np.genfromtxt('file6.txt')
dyR2 = np.genfromtxt('file7.txt')
yR4 = np.genfromtxt('file8.txt')
dyR4 = np.genfromtxt('file9.txt')
yA = np.genfromtxt('file10.txt')
dyA = np.genfromtxt('file11.txt')



def cm_to_inch(value): #������� ��� ����,����� ����� ���� ��������� ������ �������
    return value/2.54
fig=plt.figure(figsize=(cm_to_inch(40),cm_to_inch(20))) #������ ������ �������
ax1=fig.add_subplot(121)
ax2=fig.add_subplot(122)
ax1.set_title('������� y',fontsize=15,fontweight="bold")
ax2.set_title('������� dy',fontsize=15,fontweight="bold")
ax1.plot(x,yE1,linewidth=2,color='red',label="�����")
ax1.plot(x,yE2,linewidth=7,color='yellow',label="����� � ���������")
ax1.plot(x,yR2,linewidth=5,color='pink',label="�����-����� 2 ���")
ax1.plot(x,yR4,linewidth=2,color='blue',label="�����-����� 4 ���")
ax1.plot(x,yA,linewidth=1,color='green',label="������ 3 ���")
ax1.legend() 
ax1.grid()
#ax2.scatter(x,dyE1,25,color='red',label="�����")
#ax2.scatter(x,dyE2,55,color='yellow',label="����� � ���������")
#ax2.scatter(x,dyR2,35,color='pink',label="�����-����� 2 ���")
#ax2.scatter(x,dyR4,20,color='blue',label="�����-����� 4 ���")
#ax2.scatter(x,dyA,8,color='green',label="������ 3 ���")
ax2.plot(x,dyE1,linewidth=2,color='red',label="�����")
ax2.plot(x,dyE2,linewidth=7,color='yellow',label="����� � ���������")
ax2.plot(x,dyR2,linewidth=5,color='pink',label="�����-����� 2 ���")
ax2.plot(x,dyR4,linewidth=2,color='blue',label="�����-����� 4 ���")
ax2.plot(x,dyA,linewidth=1,color='green',label="������ 3 ���")
ax2.legend() 
ax2.grid()
plt.show()
