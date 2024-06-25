import matplotlib.pyplot as plt
import numpy as np
from lib_BCS import *

f=open("test_GC01.txt","r")
lines=f.readlines()
lines=lines[1:]
f.close()

a=lines[0].split(",")
ny=int(a[0])
nt=int(a[1].replace("\n",""))
print(ny,nt)
T1_out=np.zeros((nt,ny))
T2_out=np.zeros((nt,ny))
psi1_out=np.zeros((nt,ny))
psi2_out=np.zeros((nt,ny))
tauy_out=np.zeros((nt,ny))

tmp=lines[1].split(",")
y=np.asarray([float(i.replace("\n","")) for i in tmp])
print(y)
for it in range(0,nt):
    tmp=lines[3+it*6].split(",")
    tmp2=np.asarray([float(i.replace("\n","")) for i in tmp])
    T1_out[it,:]=tmp2
    tmp=lines[4+it*6].split(",")
    tmp2=np.asarray([float(i.replace("\n","")) for i in tmp])
    T2_out[it,:]=tmp2
    tmp=lines[5+it*6].split(",")
    tmp2=np.asarray([float(i.replace("\n","")) for i in tmp])
    psi1_out[it,:]=tmp2
    tmp=lines[6+it*6].split(",")
    tmp2=np.asarray([float(i.replace("\n","")) for i in tmp])
    psi2_out[it,:]=tmp2
    tmp=lines[7+it*6].split(",")
    tmp2=np.asarray([float(i.replace("\n","")) for i in tmp])
    tauy_out[it,:]=tmp2
T1_out=T1_out[12*20:,20]
T2_out=T2_out[12*20:,20]
T1_clm=clmmon(T1_out)              
T2_clm=clmmon(T2_out)
T1_anm=calcmonanom(T1_out,T1_clm)
T2_anm=calcmonanom(T2_out,T2_clm)
plt.subplot(2,1,1)
plt.plot(T1_anm)
plt.plot(T2_anm)
#plt.xlim([0,2400])
a=np.corrcoef(T1_anm,T2_anm)[0,1]
print(a)
plt.subplot(2,1,2)
lags,lag_c=lagcor(T1_anm,T2_anm,12*20)
plt.plot(lags,lag_c)
plt.show()
