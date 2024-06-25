import f90nml
import numpy as np
import os

def set_param(fname):
    params=f90nml.read(fname)
    ntime_year=params['params']['ntime_year']
    dt_year=params['params']['dt_year']
    print(dt_year)
    istep_out=params['params']['istep_out']
    Lx_1=params['params']['Lx_1']
    Lx_2=params['params']['Lx_2']
    H_1=params['params']['H_1']
    H_2=params['params']['H_2']
    R_1=params['params']['R_1']
    R_2=params['params']['R_2']
    delta_1=params['params']['delta_1']
    delta_2=params['params']['delta_2']
    Lx=params['params']['Lx']
    Ly=params['params']['Ly']
    D=params['params']['D']
    lamd=params['params']['lamd']
    rho_o=params['params']['rho_o']
    rho_w=params['params']['rho_w']
    C_pa=params['params']['C_pa']
    C_pw=params['params']['C_pw']
    nu_o=params['params']['nu_o']
    eps=params['params']['eps']
    FoA=params['params']['FoA']
    B=params['params']['B']
    F_1=params['params']['F_1']
    S=params['params']['S']
    f_0=params['params']['f_0']
    beta=params['params']['beta']
    d_l=params['params']['d_l']
    amp_noise=params['params']['amp_noise']
    amp_noise_theta_1=params['params']['amp_noise_theta_1']
    amp_noise_theta_2=params['params']['amp_noise_theta_2']
    dt_year=params['params']['dt_year']
    return(ntime_year,dt_year,istep_out,Lx_1,Lx_2,H_1,H_2,R_1,R_2,delta_1,delta_2,Lx,Ly,D,lamd,rho_o,rho_w,C_pa,C_pw,nu_o,eps,FoA,B,F_1,S,f_0,beta,d_l,amp_noise,amp_noise_theta_1,amp_noise_theta_2,dt_year)
def get_weight(n,dt,to):
    weight=np.zeros((n))
    if (n == 2):
        weight[0]=0.5*(2.0*dt-to)*to/(dt*to)
        weight[1]=0.5*to*to/(dt*to)
    elif (n == 3):
        weight[0]=0.5*dt/to
        weight[1]=0.5*dt/to+0.5*(3.0*dt-to)*(to-dt)/(dt*to)
        weight[2]=0.5*(to-dt)*(to-dt)/(dt*to)
    else:        
        weight[0]=0.5*dt/to
        for i in range(1,n-2):
            weight[i]=1.0*dt/to
        weight[n-2]=0.5*dt/to+0.5*(n*dt-to)*(to-(n-2)*dt)/(dt*to)
        weight[n-1]=0.5*(to-(n-2)*dt)*(to-(n-2)*dt)/(dt*to)
    return(weight)

def solve_tri_real(N,A_in,B_in,C_in,D_in):
    A=np.copy(A_in);B=np.copy(B_in)
    C=np.copy(C_in);D=np.copy(D_in)
    U=np.zeros((N))
    for i in range(1,N):
        m = A[i] / B[i-1]
        B[i] = B[i] - m * C[i-1]
        D[i] = D[i] - m * D[i-1]
    U[N-1] = D[N-1] / B[N-1]
    for i in range(N-2,-1,-1):
        U[i] = (D[i]-C[i]*U[i+1]) / B[i]
    return(U)
def clmmon(var):
    var_clm=var[0:12].copy()
    for i in range(0,12):
        var_clm[i]=np.nanmean(var[i::12],axis=0)
    return(var_clm)
def calcmonanom(var,var_clm):
    var_anom=var.copy()
    for i in range(0,12):
        var_anom[i::12]=var[i::12]-var_clm[i]
    return(var_anom)
def write_state(fname,time,T1,T2,psi1,psi2,tauy):
    f=open(fname,"a")
    f.write(str(time)+"\n")
    line=""
    line+=str(T1[0])
    for iy in range(1,len(T1)):
        line+=","+str(T1[iy])
    line+="\n"
    line+=str(T2[0])
    for iy in range(1,len(T2)):
        line+=","+str(T2[iy])
    line+="\n"
    line+=str(psi1[0])
    for iy in range(1,len(psi1)):
        line+=","+str(psi1[iy])
    line+="\n"
    line+=str(psi2[0])
    for iy in range(1,len(psi2)):
        line+=","+str(psi2[iy])
    line+="\n"
    line+=str(tauy[0])
    for iy in range(1,len(tauy)):
        line+=","+str(tauy[iy])
    line+="\n"
    f.write(line)
    f.close()
def lagcor(x1,x2,nlag):
    lags=np.arange(-nlag,nlag+1,1)
    ntime=len(x1)
    cor=np.zeros((2*nlag+1))
    for ilag in range(0,nlag+1):
        x1_tmp=np.copy(x1[0:ntime-ilag])
        x2_tmp=np.copy(x2[ilag:ntime])
        cor[nlag-ilag]=np.corrcoef(x1_tmp,x2_tmp)[0,1]
        x1_tmp=np.copy(x1[ilag:ntime])
        x2_tmp=np.copy(x2[0:ntime-ilag])
        cor[nlag+ilag]=np.corrcoef(x1_tmp,x2_tmp)[0,1]
    return(lags,cor)
