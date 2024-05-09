import os
from numpy import*
import matplotlib.pyplot as plt
import sys

def show(r):
    x=linspace(-0.1,3.1,301)
    y=linspace(-0.1,1.1,101)
    X,Y=meshgrid(x,y)
    plt.figure(1)
    plt.contourf(X,Y,r)
    plt.colorbar()
    plt.legend()
    plt.title('t=4.0')
    plt.xlabel("x")
    plt.show()

rho=zeros((400,400),dtype=float32)
f = open("output5.txt","r")
rho=loadtxt('output5.txt')
f.close()
#for i in range(101):
    #for j in range(301):
        #if (rho[i,j]>8): rho[i,j]=8
for i in range(20):
    for j in range(61,301):
        rho[i,j]=-1
show(rho)



