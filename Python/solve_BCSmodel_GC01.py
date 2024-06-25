import numpy as np
import matplotlib.pyplot as plt
from lib_BCS import *

fname_param="run_GC01_case1.nml" # Namelist parameter

# Load parameter
ntime_year,dt_year,istep_out,Lx_1,Lx_2,H_1,H_2,R_1,R_2,delta_1,delta_2,Lx,Ly,D,lamd,rho_o,rho_w,C_pa,C_pw,nu_o,\
eps,FoA,B,F_1,S,f_0,beta,d_l,amp_noise,amp_noise_theta_1,amp_noise_theta_2,dt_year=set_param(fname_param)
ntime=int(ntime_year/dt_year)
print(ntime)
dy=0.01*Ly   # Meridonal Resolution
ny=int(Ly/dy)+1  # Number of points
fname_out="test_GC01.txt"

year_to_sec=60.0*60.0*24.0*365.0
dt=dt_year*year_to_sec

y=np.zeros((ny));F=np.zeros((ny));F_ori=np.zeros((ny))
spat=np.zeros((ny))

# Basin 1 temperature
T_1=np.zeros((ny)); T_1_new=np.zeros((ny)); T_1_out=np.zeros((ny))
T_2=np.zeros((ny)); T_2_new=np.zeros((ny)); T_2_out=np.zeros((ny))
noise_theta_dyn=np.zeros((ny))
noise_theta_1=np.zeros((ny))
noise_theta_2=np.zeros((ny))
theta=np.zeros((ny));  theta_out=np.zeros((ny))
theta_1=np.zeros((ny));  theta_2=np.zeros((ny))
tau=np.zeros((ny));  tauy=np.zeros((ny))
psi_1=np.zeros((ny));  psi_2=np.zeros((ny))
psi_1_out=np.zeros((ny));  psi_2_out=np.zeros((ny))
MatA_1=np.zeros((ny));  MatB_1=np.zeros((ny))
MatC_1=np.zeros((ny));  MatD_1=np.zeros((ny))
MatA_2=np.zeros((ny));  MatB_2=np.zeros((ny))
MatC_2=np.zeros((ny));  MatD_2=np.zeros((ny))

for iy in range(0,ny):
    y[iy]=iy*dy
    F_ori[iy]=np.cos(np.pi*y[iy]/Ly)
    spat[iy] = np.sin(np.pi*y[iy]/Ly)

ntime_out=int(ntime/istep_out)
time_out=np.zeros((ntime_out))
T1_out=np.zeros((ntime_out))
T2_out=np.zeros((ntime_out))


de=d_l*D/(d_l+D)
c_1=beta*R_1*R_1
c_2=beta*R_2*R_2
ratio_1=Lx_1/Lx; ratio_2=Lx_2/Lx
upsilon_1=C_pw*rho_w*H_1*0.5/(lamd*delta_1*Lx_1)
upsilon_2=C_pw*rho_w*H_2*0.5/(lamd*delta_2*Lx_2)
gamma=1.0/(B+(ratio_1+ratio_2)*lamd+C_pa*rho_o*nu_o*de)
lamd_C=B+C_pa*rho_o*nu_o*de
a_1=lamd*gamma/(C_pw*rho_w*H_1)
a_2=lamd*gamma/(C_pw*rho_w*H_2)
to_1=(Lx_1/c_1)
to_2=(Lx_2/c_2)
n_record_1=int(to_1/dt)+1
n_record_2=int(to_2/dt)+1
weight_record_1=np.zeros((n_record_1))
weight_record_2=np.zeros((n_record_2))
tauy_record_to1=np.zeros((ny,n_record_1))
tauy_record_to2=np.zeros((ny,n_record_2))
timeo1=Lx_1/c_1
timeo2=Lx_2/c_1
print("to1="+str(timeo1/(60*60*24*365))+"yr,to2="+str(timeo2/(60*60*24*365))+"yr")



weight_record_1=get_weight(n_record_1,dt,to_1)
weight_record_2=get_weight(n_record_2,dt,to_2)
itp=0

f_out=open(fname_out,"w")
line="# 1st line: y-t dimensionality 2nd line: y coordinate, 3,4,5,6: time,T1,T2,psi1,psi2,tauy \n"
line+=str(len(y))+","+str(ntime_out)+"\n"
line+=str(y[0])
for iy in range(1,len(y)):
    line+=","+str(y[iy])
f_out.write(line+"\n")
f_out.close()
for itime in range(0,ntime):
    noise_theta_dyn=np.random.randn(ny)*amp_noise
    noise_theta_1=np.random.randn(ny)*amp_noise_theta_1
    noise_theta_2=np.random.randn(ny)*amp_noise_theta_2
    F=F_ori
    # Atmospheric model
    rhs1=F_1*F*gamma+noise_theta_dyn
    rhs2=(ratio_1*lamd*T_1-ratio_1*lamd*noise_theta_1)*gamma
    rhs3=(ratio_2*lamd*T_2-ratio_2*lamd*noise_theta_2)*gamma
    theta[0:ny]=rhs1+rhs2+rhs3
    rhs_tauy_1=de*rho_o*nu_o*beta*(y-0.5*Ly)
    rhs_tauy_2=de*rho_o*nu_o*f_0*(theta)/(S*d_l)
    tauy=rhs_tauy_1+rhs_tauy_2

    psi_1=(to_1*R_1*R_1/(rho_w*H_1))*np.sum(tauy_record_to1*weight_record_1,axis=1)
    psi_2=(to_2*R_2*R_2/(rho_w*H_2))*np.sum(tauy_record_to2*weight_record_2,axis=1)


    iy=0
    Rn_1=upsilon_1*((psi_1[0]+psi_1[1])**2)/(4.0*(y[1]-y[0])*(y[1]-y[0]))+eps/((y[1]-y[0])*(y[1]-y[0]))
    Qn_1=-2.0*Rn_1-a_1*(lamd_C+lamd*ratio_2)
    MatA_1[iy]=0.0
    MatB_1[iy]=1.0-0.5*dt*Qn_1
    MatC_1[iy]=-1.0*dt*Rn_1
    
    rhs1=F_1*F*gamma+noise_theta_dyn
    rhs2=(ratio_1*lamd*T_1-ratio_1*lamd*noise_theta_1)*gamma
    rhs3=(ratio_2*lamd*T_2-ratio_2*lamd*noise_theta_2)*gamma
    theta[0:ny]=rhs1+rhs2+rhs3
    n_1=lamd/(C_pw*rho_w*H_1)
    
    MatD_1[iy]=T_1[iy]+dt*a_1*lamd*ratio_2*T_2[iy]+dt*a_1*F_1*F[iy] \
        -dt*a_1*lamd*ratio_1*noise_theta_1[iy] \
        -dt*a_1*lamd*ratio_2*noise_theta_2[iy] + \
        +0.5*dt*Qn_1*T_1[iy]+dt*Rn_1*T_1[iy+1] \
        +dt*a_1*(noise_theta_dyn[iy]+noise_theta_1[iy])/gamma

    Rn_2=upsilon_2*((psi_2[0]+psi_2[1])**2)/(4.0*(y[1]-y[0])*(y[1]-y[0]))+eps/((y[1]-y[0])*(y[1]-y[0]))
    Qn_2=-2.0*Rn_2-a_2*(lamd_C+lamd*ratio_1)
    MatA_2[iy]=0.0
    MatB_2[iy]=1.0-0.5*dt*Qn_2
    MatC_2[iy]=-1.0*dt*Rn_2
    MatD_2[iy]=T_2[iy]+dt*a_2*lamd*ratio_1*T_1[iy]+dt*a_2*F_1*F[iy] \
        -dt*a_2*lamd*ratio_1*noise_theta_1[iy]  \
        -dt*a_2*lamd*ratio_2*noise_theta_2[iy]  \
        +0.5*dt*Qn_2*T_2[iy]+dt*Rn_2*T_2[iy+1]  \
        +dt*a_2*(noise_theta_dyn[iy]+noise_theta_2[iy])/gamma

    for iy in range(1,ny-1):
        Pn_1=0.5*upsilon_1*((psi_1[iy]+psi_1[iy-1])**2)/  \
            ((y[iy+1]-y[iy-1])*(y[iy]-y[iy-1]))+2.0*eps/((y[iy+1]-y[iy-1])*(y[iy]-y[iy-1]))
        Rn_1=0.5*upsilon_1*((psi_1[iy+1]+psi_1[iy])**2)/  \
            ((y[iy+1]-y[iy-1])*(y[iy+1]-y[iy]))+2.0*eps/((y[iy+1]-y[iy-1])*(y[iy+1]-y[iy]))
        Qn_1=-1.0*Pn_1-Rn_1-a_1*(lamd_C+lamd*ratio_2)
        MatA_1[iy]=-0.5*dt*Pn_1
        MatB_1[iy]=1.0-0.5*dt*Qn_1
        MatC_1[iy]=-0.5*dt*Rn_1
        MatD_1[iy]=T_1[iy]+dt*a_1*lamd*ratio_2*T_2[iy]  \
            -dt*a_1*lamd*ratio_1*noise_theta_1[iy]-dt*a_1*lamd*ratio_2*noise_theta_2[iy] +  \
            dt*a_1*F_1*F[iy]+0.5*dt*Pn_1*T_1[iy-1]+0.5*dt*Qn_1*T_1[iy]+0.5*dt*Rn_1*T_1[iy+1]  \
            +dt*a_1*(noise_theta_dyn[iy]+noise_theta_1[iy])/gamma

        Pn_2=0.5*upsilon_2*((psi_2[iy]+psi_2[iy-1])**2)/  \
            ((y[iy+1]-y[iy-1])*(y[iy]-y[iy-1]))+2.0*eps/((y[iy+1]-y[iy-1])*(y[iy]-y[iy-1]))
        Rn_2=0.5*upsilon_2*((psi_2[iy+1]+psi_2[iy])**2)/  \
            ((y[iy+1]-y[iy-1])*(y[iy+1]-y[iy]))+2.0*eps/((y[iy+1]-y[iy-1])*(y[iy+1]-y[iy]))
        Qn_2=-Pn_2-Rn_2-a_2*(lamd_C+lamd*ratio_1)
        MatA_2[iy]=-0.5*dt*Pn_2
        MatB_2[iy]=1.0-0.5*dt*Qn_2
        MatC_2[iy]=-0.5*dt*Rn_2
        MatD_2[iy]=T_2[iy]+dt*a_2*lamd*ratio_1*T_1[iy]  \
            -dt*a_2*lamd*ratio_1*noise_theta_1[iy]-dt*a_2*lamd*ratio_2*noise_theta_2[iy] + \
            +dt*a_2*F_1*F[iy]+0.5*dt*Pn_2*T_2[iy-1]+0.5*dt*Qn_2*T_2[iy]+0.5*dt*Rn_2*T_2[iy+1]  \
            +dt*a_2*(noise_theta_dyn[iy]+noise_theta_2[iy])/gamma

    iy=ny-1
    Pn_1=upsilon_1*((psi_1[ny-1]+psi_1[ny-2])**2)/(4.0*(y[ny-1]-y[ny-2])*(y[ny-1]-y[ny-2]))+eps/((y[ny-1]-y[ny-2])*(y[ny-1]-y[ny-2]))
    Qn_1=-2.0*Pn_1-a_1*(lamd_C+lamd*ratio_2)
    MatA_1[iy]=-1.0*dt*Pn_1
    MatB_1[iy]=1.0-0.5*dt*Qn_1
    MatC_1[iy]=0.0
    MatD_1[iy]=T_1[iy] \
        +dt*a_1*lamd*ratio_2*T_2[iy]\
        +dt*a_1*F_1*F[iy]+dt*Pn_1*T_1[iy-1]+0.5*dt*Qn_1*T_1[iy] \
        -dt*a_1*lamd*ratio_1*noise_theta_1[iy]-dt*a_1*lamd*ratio_2*noise_theta_2[iy] \
        +dt*a_1*(noise_theta_dyn[iy]+noise_theta_1[iy])/gamma

    Pn_2=upsilon_2*((psi_2[ny-1]+psi_2[ny-2])**2)/(4.0*(y[ny-1]-y[ny-2])*(y[ny-1]-y[ny-2]))+eps/((y[ny-1]-y[ny-2])*(y[ny-1]-y[ny-2]))
    Qn_2=-2.0*Pn_2-a_2*(lamd_C+lamd*ratio_1)
    MatA_2[iy]=-1.0*dt*Pn_2
    MatB_2[iy]=1.0-0.5*dt*Qn_2
    MatC_2[iy]=0.0
    MatD_2[iy]=T_2[iy]+dt*a_2*lamd*ratio_1*T_1[iy] \
        +dt*a_2*F_1*F[iy]+dt*Pn_2*T_2[iy-1]+0.5*dt*Qn_2*T_2[iy] \
        -dt*a_2*lamd*ratio_1*noise_theta_1[iy] \
        -dt*a_2*lamd*ratio_2*noise_theta_2[iy] \
        +dt*a_2*(noise_theta_dyn[iy]+noise_theta_2[iy])/gamma

    # Solve tridiagonal matrix problem
    T_1_new=solve_tri_real(ny,MatA_1,MatB_1,MatC_1,MatD_1)
    T_2_new=solve_tri_real(ny,MatA_2,MatB_2,MatC_2,MatD_2)
    # Update solutions
    T_1=np.copy(T_1_new)
    T_2=np.copy(T_2_new)

    for i in range(n_record_1-1,0,-1):
        tauy_record_to1[0:ny,i]=np.copy(tauy_record_to1[0:ny,i-1])
    tauy_record_to1[0:ny,0]=np.copy(tauy[0:ny])

    for i in range(n_record_2-1,0,-1):
        tauy_record_to2[0:ny,i]=np.copy(tauy_record_to2[0:ny,i-1])
    tauy_record_to2[0:ny,0]=np.copy(tauy[0:ny])
    if (itime%istep_out== 0):
        T1_out[itp]=T_1[20]
        T2_out[itp]=T_2[20]        
        itp=itp+1
        write_state(fname_out,itime*dt_year,T_1,T_2,psi_1,psi_2,tauy)
T1_out=T1_out[12*20:]
T2_out=T2_out[12*20:]
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
