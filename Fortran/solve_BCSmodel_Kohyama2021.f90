program solve_BCSmodel_Kohyama2021
  use netCDF
  implicit none
  ! Compile method
  ! Standard case
  ! export FPATH_FORT=-I/Users/kido/libs/netcdf-fortran-4.4.4/include -I/Users/kido/libs/netcdf-c-4.6.1/include
  ! export LPATH_FORT=-L/Users/kido/libs/netcdf-fortran-4.4.4/lib
  ! ifort -o exec_ctl.out solve_BCSmodel_Kohyama2021.f90  $FPATH_FORT  $LPATH_FORT   -lnetcdff -lnetcdf -cpp
  ! ./exec_ctl.out < param_withnoise_gc01_ctl.nml
  ! Options
  ! ext_file: Load external climatological file
  ! nodyn: Atmospheric temperature used to calculate the wind stress
  !        is replaced by the climatological value
  ! ifort -o exec_ctl.out solve_BCSmodel_Kohyama2021.f90  $FPATH_FORT  $LPATH_FORT   -lnetcdff -lnetcdf -nodyn -cpp
  ! ./exec_nodyn.out < param_withnoise_gc01_nodyn.nml
  ! notherm : temperature of other basin does not affect the another basin
  ! ifort -o exec_ctl.out solve_BCSmodel_Kohyama2021.f90  $FPATH_FORT  $LPATH_FORT   -lnetcdff -lnetcdf -notherm -cpp
  ! ./exec_notherm.out < param_withnoise_gc01_notherm.nml
  !
  integer,parameter :: idx = 8,maxlen=400
  real(idx) :: Lx_1,Lx_2,H_1,H_2,R_1,R_2,delta_1,delta_2
  real(idx) :: Lx,Ly,D,lambda,rho_o,rho_w,C_pa,C_pw,nu_o,eps
  real(idx) :: FoA,B,F_1,S,f_0,beta,de,d_l
  real(idx) :: c_1,c_2,ratio_1,ratio_2,a_1,a_2
  real(idx) :: upsilon_1,upsilon_2,gamma,lambda_C,to_1,to_2
  real(idx) :: dy,dt,dt_year
  integer :: ny,ntime,ntime_out
  integer :: iy,itime,i,irec  
  real(idx),allocatable :: y(:),F(:),spat(:),time_out(:)
  real(idx),parameter :: pi = 4.0_idx * atan(1.0_idx)
  real(idx),parameter  :: sec_to_year=1.0_idx/(60.0_idx*60.0_idx*24.0_idx*365.0_idx)
  real(idx),parameter  :: year_to_sec=60.0_idx*60.0_idx*24.0_idx*365.0_idx
  real(idx),allocatable :: T_1(:),T_2(:),T_1_new(:),T_2_new(:)
  real(4),allocatable :: T_1_out(:),T_2_out(:),psi_1_out(:),psi_2_out(:),theta_out(:)
  real(4),allocatable :: tau(:) !! Added 
  real(idx),allocatable :: tauy(:) !! Added 
  real(idx),allocatable :: psi_1(:),psi_2(:),theta(:)
  real(idx),allocatable :: MatA_1(:),MatB_1(:),MatC_1(:),MatD_1(:)
  real(idx),allocatable :: MatA_2(:),MatB_2(:),MatC_2(:),MatD_2(:)
  integer :: n_record_1,n_record_2
  real(idx),allocatable :: tauy_record_to1(:,:),tauy_record_to2(:,:)
  real(idx),allocatable:: noise_theta_dyn(:)
  real(idx),allocatable:: noise_theta1(:)
  real(idx),allocatable:: noise_theta2(:)
  real(idx),allocatable :: weight_record_1(:),weight_record_2(:)
  real(idx) :: rhs1,rhs2,rhs3,rhs4
  real(idx) :: Pn_1,Qn_1,Rn_1
  real(idx) :: Pn_2,Qn_2,Rn_2
  real(idx) :: amp_noise,x1,x2
  real(idx) :: Hml_1,Hml_2
  real(idx) :: amp_noise_theta1,amp_noise_theta2
  integer :: istep_out
  character(len=maxlen) :: fname_out
  integer :: seedsize
  integer :: ncid
  integer :: dim1_dimid,dim1_varid
  integer :: dim2_dimid,dim2_varid
  integer :: T1_varid,T2_varid,psi1_varid,psi2_varid
  integer :: theta_varid,tau_varid
  character(len=maxlen) :: fname_mean
  integer :: varid  
  real(idx),allocatable :: T_1_mean(:),T_2_mean(:),psi_1_mean(:),psi_2_mean(:),theta_mean(:)
  integer,allocatable :: seed(:)
  namelist/params_basin/Lx_1,Lx_2,H_1,H_2,R_1,R_2,delta_1,delta_2
  namelist/params_common/Lx,Ly,D,lambda,rho_o,rho_w,C_pa,C_pw
  namelist/params_common/nu_o,eps,FoA,B,F_1,S,f_0,beta,d_l,amp_noise,amp_noise_theta1,amp_noise_theta2
  namelist/params_time/dt_year,ntime,istep_out
  namelist/io/fname_out,fname_mean
# if defined(notherm) || defined(nodyn)
#  define ext_file
#endif
  read(5,params_basin)
  read(5,params_common)
  read(5,params_time)
  read(5,io)

  ! Coordinate setting
  dy=0.01_idx*Ly   ! Meridonal Resolution
  ny=int(Ly/dy)+1  ! Number of points
  ! Time step setting
  write(*,*) dt_year
  dt=dt_year*year_to_sec
  write(*,*) "Dt=",dt*sec_to_year,"years"
  write(*,*) ntime
  Hml_1=100.0;  Hml_2=100.0;
  allocate(y(ny))
  allocate(F(ny));  allocate(spat(ny))
  allocate(T_1(ny)); allocate(T_1_new(ny)); allocate(T_1_out(ny))
  allocate(T_2(ny)); allocate(T_2_new(ny)); allocate(T_2_out(ny))
  allocate(noise_theta_dyn(ny))
  allocate(noise_theta1(ny))
  allocate(noise_theta2(ny))
  allocate(theta(ny));  allocate(theta_out(ny))
  allocate(tauy(ny))
  allocate(tau(ny))
  allocate(psi_1(ny));   allocate(psi_2(ny))
  allocate(psi_1_out(ny)) ;  allocate(psi_2_out(ny))
  allocate(MatA_1(ny)); allocate(MatB_1(ny)); allocate(MatC_1(ny)); allocate(MatD_1(ny))
  allocate(MatA_2(ny)); allocate(MatB_2(ny)); allocate(MatC_2(ny)); allocate(MatD_2(ny))
#ifdef ext_file
  allocate(T_1_mean(ny));  allocate(T_2_mean(ny))
  allocate(PSI_1_mean(ny));  allocate(PSI_2_mean(ny))
  allocate(theta_mean(ny))
  write(*,*) "Read mean file="//trim(fname_mean) 
  call check_w(nf90_open(trim(fname_mean),nf90_nowrite,ncid))
  call check_w( nf90_inq_varid(ncid,"sst_1", varid) )
  call check_w( nf90_get_var(ncid, varid, T_1_mean, start = (/1,1/), &
       count = (/ny,1/)) )
  call check_w( nf90_inq_varid(ncid,"sst_2", varid) )
  call check_w( nf90_get_var(ncid, varid, T_2_mean, start = (/1,1/), &
       count = (/ny,1/)) )
  call check_w( nf90_inq_varid(ncid,"psi_1", varid) )
  call check_w( nf90_get_var(ncid, varid, psi_1_mean, start = (/1,1/), &
       count = (/ny,1/)) )
  call check_w( nf90_inq_varid(ncid,"psi_2", varid) )
  call check_w( nf90_get_var(ncid, varid, psi_2_mean, start = (/1,1/), &
       count = (/ny,1/)) )
  call check_w( nf90_inq_varid(ncid,"theta", varid) )
  call check_w( nf90_get_var(ncid, varid, theta_mean, start = (/1,1/), &
       count = (/ny,1/)) )
  call CHECK_W(NF90_CLOSE(ncid))
#endif
  do iy = 1,ny
     y(iy) = (iy-1)*dy
     F(iy) = cos(pi * y(iy)/Ly)
     spat(iy) = sin(pi * y(iy)/Ly)
  end do
  ntime_out=ntime/istep_out; allocate(time_out(ntime_out))
  write(*,*) ntime_out
  do itime = 1,ntime_out
     time_out(itime)=int((dt_year*12)*istep_out*itime)
     write(*,*) time_out(itime)
  end do
  ! Output setting
  ! create file
  call CHECK_W(NF90_CREATE(trim(fname_out),nf90_netcdf4, ncid) )
  ! Define the dimensions
  call CHECK_W(NF90_DEF_DIM(ncid, trim("lat"), ny,dim1_dimid))
  call CHECK_W(NF90_DEF_VAR(ncid, trim("lat"), nf90_float,dim1_dimid,dim1_varid))
  call CHECK_W(NF90_PUT_ATT(ncid, dim1_varid,"units", "degrees_north"))

  call CHECK_W(NF90_DEF_DIM(ncid, trim("time"), ntime_out,dim2_dimid))
  call CHECK_W(NF90_DEF_VAR(ncid, trim("time"), nf90_float,dim2_dimid,dim2_varid))
  call CHECK_W(NF90_PUT_ATT(ncid, dim2_varid,"units", "month since 1900-01-01 00:00:00"))

  call CHECK_W(NF90_DEF_VAR(ncid, trim("sst_1"), nf90_float,(/dim1_dimid,dim2_dimid/),T1_varid))
  call CHECK_W(NF90_DEF_VAR(ncid, trim("sst_2"), nf90_float,(/dim1_dimid,dim2_dimid/),T2_varid))
  call CHECK_W(NF90_DEF_VAR(ncid, trim("psi_1"), nf90_float,(/dim1_dimid,dim2_dimid/),psi1_varid))
  call CHECK_W(NF90_DEF_VAR(ncid, trim("psi_2"), nf90_float,(/dim1_dimid,dim2_dimid/),psi2_varid))
  call CHECK_W(NF90_DEF_VAR(ncid, trim("tau"), nf90_float,(/dim1_dimid,dim2_dimid/),tau_varid))
  call CHECK_W(NF90_DEF_VAR(ncid, trim("theta"), nf90_float,(/dim1_dimid,dim2_dimid/),theta_varid))
  call CHECK_W(NF90_PUT_ATT(ncid, NF90_GLOBAL ,"Lx_1", Lx_1))
  call CHECK_W(NF90_PUT_ATT(ncid, NF90_GLOBAL ,"Lx_2", Lx_2))
  call CHECK_W(NF90_PUT_ATT(ncid, NF90_GLOBAL ,"H_1", H_1))
  call CHECK_W(NF90_PUT_ATT(ncid, NF90_GLOBAL ,"H_2", H_2))
  call CHECK_W(NF90_PUT_ATT(ncid, NF90_GLOBAL ,"R_1", R_1))
  call CHECK_W(NF90_PUT_ATT(ncid, NF90_GLOBAL ,"R_2", R_2))
  call CHECK_W(NF90_PUT_ATT(ncid, NF90_GLOBAL ,"delta_1", delta_1))
  call CHECK_W(NF90_PUT_ATT(ncid, NF90_GLOBAL ,"delta_2", delta_2))
  call CHECK_W(NF90_PUT_ATT(ncid, NF90_GLOBAL ,"Lx", Lx))
  call CHECK_W(NF90_PUT_ATT(ncid, NF90_GLOBAL ,"Ly", Ly))
  call CHECK_W(NF90_PUT_ATT(ncid, NF90_GLOBAL ,"D", D))
  call CHECK_W(NF90_PUT_ATT(ncid, NF90_GLOBAL ,"lambda", lambda))
  call CHECK_W(NF90_PUT_ATT(ncid, NF90_GLOBAL ,"rho_o", rho_o))
  call CHECK_W(NF90_PUT_ATT(ncid, NF90_GLOBAL ,"rho_w", rho_w))
  call CHECK_W(NF90_PUT_ATT(ncid, NF90_GLOBAL ,"C_pa", C_pa))
  call CHECK_W(NF90_PUT_ATT(ncid, NF90_GLOBAL ,"C_pw", C_pw))
  call CHECK_W(NF90_PUT_ATT(ncid, NF90_GLOBAL ,"nu_o", nu_o))
  call CHECK_W(NF90_PUT_ATT(ncid, NF90_GLOBAL ,"eps", eps))
  call CHECK_W(NF90_PUT_ATT(ncid, NF90_GLOBAL ,"FoA", FoA))
  call CHECK_W(NF90_PUT_ATT(ncid, NF90_GLOBAL ,"B", B))
  call CHECK_W(NF90_PUT_ATT(ncid, NF90_GLOBAL ,"F_1", F_1))
  call CHECK_W(NF90_PUT_ATT(ncid, NF90_GLOBAL ,"S", S))
  call CHECK_W(NF90_PUT_ATT(ncid, NF90_GLOBAL ,"f_0",f_0))
  call CHECK_W(NF90_PUT_ATT(ncid, NF90_GLOBAL ,"beta",beta))
  call CHECK_W(NF90_PUT_ATT(ncid, NF90_GLOBAL ,"d_l",d_l))
  call CHECK_W(NF90_PUT_ATT(ncid, NF90_GLOBAL ,"amp_noise",amp_noise))
  call CHECK_W(NF90_PUT_ATT(ncid, NF90_GLOBAL ,"amp_noise_theta1",amp_noise_theta1))
  call CHECK_W(NF90_PUT_ATT(ncid, NF90_GLOBAL ,"amp_noise_theta2",amp_noise_theta2))
#ifdef ext_file
  call CHECK_W(NF90_PUT_ATT(ncid, NF90_GLOBAL ,"fname_mean",fname_mean))
#endif
#ifdef nodyn
  call CHECK_W(NF90_PUT_ATT(ncid, NF90_GLOBAL ,"Dynamical coupling","False"))
#else
  call CHECK_W(NF90_PUT_ATT(ncid, NF90_GLOBAL ,"Dynamical Coupling","True"))
#endif
#ifdef notherm
  call CHECK_W(NF90_PUT_ATT(ncid, NF90_GLOBAL ,"Thermodynamical coupling","False"))
#else
  call CHECK_W(NF90_PUT_ATT(ncid, NF90_GLOBAL ,"Thermodynamical Coupling","True"))
#endif
  
  call CHECK_W(NF90_ENDDEF(ncid))
  ! Put data
  call CHECK_W( NF90_PUT_VAR(ncid,dim1_varid,y/110e3))
  call CHECK_W( NF90_PUT_VAR(ncid,dim2_varid,time_out))
!  call CHECK_W( NF90_CLOSE(ncid))


  !  Parameters
  !d_l=(de*D)/(D-de)
  de=d_l*D/(d_l+D)
  c_1=beta*R_1*R_1
  c_2=beta*R_2*R_2
  ratio_1=Lx_1/Lx
  ratio_2=Lx_2/Lx
  upsilon_1=C_pw*rho_w*H_1*0.5_idx/(lambda*delta_1*Lx_1)
  upsilon_2=C_pw*rho_w*H_2*0.5_idx/(lambda*delta_2*Lx_2)
  gamma=1.0/(B+(ratio_1+ratio_2)*lambda+C_pa*rho_o*nu_o*de)
  Lambda_C=B+C_pa*rho_o*nu_o*de
  a_1=lambda*gamma/(C_pw*rho_w*H_1)
  a_2=lambda*gamma/(C_pw*rho_w*H_2)
  to_1=(Lx_1/c_1)
  to_2=(Lx_2/c_2)

  ! ## For integration
  n_record_1=ceiling(to_1/dt)+1
  n_record_2=ceiling(to_2/dt)+1
  allocate(tauy_record_to1(ny,n_record_1));  allocate(tauy_record_to2(ny,n_record_2))
  tauy_record_to1=0.0_idx; tauy_record_to2=0.0_idx

  allocate(weight_record_1(n_record_1))
  allocate(weight_record_2(n_record_2))

  weight_record_1=get_weight(n_record_1,dt,to_1)
  weight_record_2=get_weight(n_record_2,dt,to_2)
  ! Nose setting
  call random_seed(size=seedsize)
  allocate(seed(seedsize))
  call random_seed(get=seed)
  !seed(1)=1
  !seed(2)=3
  call random_seed(put=seed)
  ! Integration
  irec=1
!  ntime=1
  !amp_noise_theta1=0.0_idx;amp_noise_theta2=0.0_idx
  !amp_noise_theta1=15.0_idx;amp_noise_theta2=15.0_idx
  do itime = 1,ntime
     ! Generate noise
     do iy =1,ny
        call random_number(x1)
        call random_number(x2)
        noise_theta_dyn(iy)=sqrt(-2.0_idx*log(x1))*cos(2.0_idx*pi*x1)*amp_noise
        call random_number(x1)
        call random_number(x2)
        noise_theta1(iy)=noise_theta_dyn(iy)+sqrt(-2.0_idx*log(x1))*cos(2.0_idx*pi*x1)*amp_noise_theta1
        call random_number(x1)
        call random_number(x2)
        noise_theta2(iy)=noise_theta_dyn(iy)+sqrt(-2.0_idx*log(x1))*cos(2.0_idx*pi*x1)*amp_noise_theta2
     end do
     do iy = 1,ny
        theta(iy)=(ratio_1*lambda*T_1(iy)+ratio_2*lambda*T_2(iy)+F_1*F(iy))*gamma+FoA/B+noise_theta_dyn(iy)
        theta(iy)=(ratio_1*lambda*T_1(iy)-ratio_1*lambda*noise_theta1(iy)+ratio_2*lambda*T_2(iy)-ratio_2*lambda*noise_theta2(iy)+F_1*F(iy))*gamma+FoA/B+noise_theta_dyn(iy)
        ! Calculate tau_y
        rhs1=de*rho_o*nu_o*beta*(y(iy)-0.5_idx*Ly)
        rhs2=de*rho_o*nu_o*f_0*(theta(iy)-FoA/B)/(S*d_l)
#if defined nodyn
        rhs2=de*rho_o*nu_o*f_0*(theta_mean(iy)-FoA/B)/(S*d_l)
#endif
        tauy(iy)=rhs1+rhs2
        psi_1(iy)=(R_1*R_1/(rho_w*H_1))*to_1*(sum(tauy_record_to1(iy,1:n_record_1)*weight_record_1(1:n_record_1)))

        psi_2(iy)=(R_2*R_2/(rho_w*H_2))*to_2*(sum(tauy_record_to2(iy,1:n_record_2)*weight_record_2(1:n_record_2)))
     end do
     iy=1
     Rn_1=upsilon_1*((psi_1(1)+psi_1(2))**2)/(4.0_idx*(y(2)-y(1))*(y(2)-y(1)))+eps/((y(2)-y(1))*(y(2)-y(1)))
     Qn_1=-2.0_idx*Rn_1-a_1*(Lambda_C+lambda*ratio_2)
     MatA_1(iy)=0.0_idx
     MatB_1(iy)=1.0_idx-0.5_idx*dt*Qn_1
     MatC_1(iy)=-1.0_idx*dt*Rn_1
     MatD_1(iy)=T_1(iy)+dt*a_1*lambda*ratio_2*T_2(iy)+dt*a_1*F_1*F(iy)&
          & -dt*a_1*lambda*ratio_1*noise_theta1(iy)-dt*a_1*lambda*ratio_2*noise_theta2(iy) +&
          & +0.5_idx*dt*Qn_1*T_1(iy)+dt*Rn_1*T_1(iy+1) &
          & +dt*a_1*(noise_theta_dyn(iy)+noise_theta1(iy))/gamma
     
#ifdef notherm     
     MatD_1(iy)=T_1(iy)+dt*a_1*lambda*ratio_2*(T_2_mean(iy)-FoA/B)+dt*a_1*F_1*F(iy)&
          & -dt*a_1*lambda*ratio_1*noise_theta1(iy) +&
          & +0.5_idx*dt*Qn_1*T_1(iy)+dt*Rn_1*T_1(iy+1) &
          & +dt*a_1*(noise_theta_dyn(iy)+noise_theta1(iy))/gamma
#endif
     Rn_2=upsilon_2*((psi_2(1)+psi_2(2))**2)/(4.0_idx*(y(2)-y(1))*(y(2)-y(1)))+eps/((y(2)-y(1))*(y(2)-y(1)))
     Qn_2=-2.0_idx*Rn_2-a_2*(Lambda_C+lambda*ratio_1)
     MatA_2(iy)=0.0_idx
     MatB_2(iy)=1.0_idx-0.5_idx*dt*Qn_2
     MatC_2(iy)=-1.0_idx*dt*Rn_2
     MatD_2(iy)=T_2(iy)+dt*a_2*lambda*ratio_1*T_1(iy)+dt*a_2*F_1*F(iy)&
          & -dt*a_2*lambda*ratio_1*noise_theta1(iy)-dt*a_2*lambda*ratio_2*noise_theta2(iy) +&
          & +0.5_idx*dt*Qn_2*T_2(iy)+dt*Rn_2*T_2(iy+1) &
          & +dt*a_2*(noise_theta_dyn(iy)+noise_theta2(iy))/gamma
#ifdef notherm     
     MatD_2(iy)=T_2(iy)+dt*a_2*lambda*ratio_1*(T_1_mean(iy)-FoA/B)+dt*a_2*F_1*F(iy)&
          & -dt*a_2*lambda*ratio_2*noise_theta2(iy) +&
          & +0.5_idx*dt*Qn_2*T_2(iy)+dt*Rn_2*T_2(iy+1) &
          & +dt*a_2*(noise_theta_dyn(iy)+noise_theta2(iy))/gamma
#endif
     do iy = 2,ny-1
        Pn_1=0.5_idx*upsilon_1*((psi_1(iy)+psi_1(iy-1))**2)/&
              & ((y(iy+1)-y(iy-1))*(y(iy)-y(iy-1)))+2.0_idx*eps/((y(iy+1)-y(iy-1))*(y(iy)-y(iy-1)))
        Rn_1=0.5_idx*upsilon_1*((psi_1(iy+1)+psi_1(iy))**2)/&
              & ((y(iy+1)-y(iy-1))*(y(iy+1)-y(iy)))+2.0_idx*eps/((y(iy+1)-y(iy-1))*(y(iy+1)-y(iy)))
        Qn_1=-Pn_1-Rn_1-a_1*(Lambda_C+lambda*ratio_2)
        MatA_1(iy)=-0.5_idx*dt*Pn_1
        MatB_1(iy)=1.0_idx-0.5_idx*dt*Qn_1
        MatC_1(iy)=-0.5_idx*dt*Rn_1
        MatD_1(iy)=T_1(iy)+dt*a_1*lambda*ratio_2*T_2(iy)&
             & -dt*a_1*lambda*ratio_1*noise_theta1(iy)-dt*a_1*lambda*ratio_2*noise_theta2(iy) +&
             & dt*a_1*F_1*F(iy)+0.5_idx*dt*Pn_1*T_1(iy-1)+0.5_idx*dt*Qn_1*T_1(iy)+0.5_idx*dt*Rn_1*T_1(iy+1) &
             & +dt*a_1*(noise_theta_dyn(iy)+noise_theta1(iy))/gamma
#ifdef notherm     
        MatD_1(iy)=T_1(iy)+dt*a_1*lambda*ratio_2*(T_2_mean(iy)-FoA/B)&
             & -dt*a_1*lambda*ratio_1*noise_theta1(iy) +&
             & dt*a_1*F_1*F(iy)+0.5_idx*dt*Pn_1*T_1(iy-1)+0.5_idx*dt*Qn_1*T_1(iy)+0.5_idx*dt*Rn_1*T_1(iy+1) &
             & +dt*a_1*(noise_theta_dyn(iy)+noise_theta1(iy))/gamma
#endif
        Pn_2=0.5_idx*upsilon_2*((psi_2(iy)+psi_2(iy-1))**2)/&
              & ((y(iy+1)-y(iy-1))*(y(iy)-y(iy-1)))+2.0_idx*eps/((y(iy+1)-y(iy-1))*(y(iy)-y(iy-1)))
        Rn_2=0.5_idx*upsilon_2*((psi_2(iy+1)+psi_2(iy))**2)/&
             & ((y(iy+1)-y(iy-1))*(y(iy+1)-y(iy)))+2.0_idx*eps/((y(iy+1)-y(iy-1))*(y(iy+1)-y(iy)))
        Qn_2=-Pn_2-Rn_2-a_2*(Lambda_C+lambda*ratio_1)
        MatA_2(iy)=-0.5_idx*dt*Pn_2
        MatB_2(iy)=1.0_idx-0.5_idx*dt*Qn_2
        MatC_2(iy)=-0.5_idx*dt*Rn_2
        MatD_2(iy)=T_2(iy)+dt*a_2*lambda*ratio_1*T_1(iy)&
             & -dt*a_2*lambda*ratio_1*noise_theta1(iy)-dt*a_2*lambda*ratio_2*noise_theta2(iy) +&
             & +dt*a_2*F_1*F(iy)+0.5_idx*dt*Pn_2*T_2(iy-1)+0.5_idx*dt*Qn_2*T_2(iy)+0.5_idx*dt*Rn_2*T_2(iy+1) &
             & +dt*a_2*(noise_theta_dyn(iy)+noise_theta2(iy))/gamma

#ifdef notherm     
        MatD_2(iy)=T_2(iy)+dt*a_2*lambda*ratio_1*(T_1_mean(iy)-FoA/B)&
             & -dt*a_2*lambda*ratio_1*noise_theta1(iy)-dt*a_2*lambda*ratio_2*noise_theta2(iy) +&
             & +dt*a_2*F_1*F(iy)+0.5_idx*dt*Pn_2*T_2(iy-1)+0.5_idx*dt*Qn_2*T_2(iy)+0.5_idx*dt*Rn_2*T_2(iy+1) &
             & +dt*a_2*(noise_theta_dyn(iy)+noise_theta2(iy))/gamma
#endif
     end do
     iy=ny
     Pn_1=upsilon_1*((psi_1(ny)+psi_1(ny-1))**2)/(4.0_idx*(y(ny)-y(ny-1))*(y(ny)-y(ny-1)))+eps/((y(ny)-y(ny-1))*(y(ny)-y(ny-1)))
     Qn_1=-2.0_idx*Pn_1-a_1*(Lambda_C+lambda*ratio_2)
     MatA_1(iy)=-1.0_idx*dt*Pn_1
     MatB_1(iy)=1.0_idx-0.5_idx*dt*Qn_1
     MatC_1(iy)=0.0_idx
     MatD_1(iy)=T_1(iy)+dt*a_1*lambda*ratio_2*T_2(iy)+dt*a_1*F_1*F(iy)+dt*Pn_1*T_1(iy-1)+0.5_idx*dt*Qn_1*T_1(iy)&
          & -dt*a_1*lambda*ratio_1*noise_theta1(iy)-dt*a_1*lambda*ratio_2*noise_theta2(iy) &
          & +dt*a_1*(noise_theta_dyn(iy)+noise_theta1(iy))/gamma
#ifdef notherm     
     MatD_1(iy)=T_1(iy)+dt*a_1*lambda*ratio_2*(T_2_mean(iy)-FoA/B)+dt*a_1*F_1*F(iy)+dt*Pn_1*T_1(iy-1)+0.5_idx*dt*Qn_1*T_1(iy)&
          & -dt*a_1*lambda*ratio_1*noise_theta1(iy)&
          & +dt*a_1*(noise_theta_dyn(iy)+noise_theta1(iy))/gamma
#endif
     Pn_2=upsilon_2*((psi_2(ny)+psi_2(ny-1))**2)/(4.0_idx*(y(ny)-y(ny-1))*(y(ny)-y(ny-1)))+eps/((y(ny)-y(ny-1))*(y(ny)-y(ny-1)))
     Qn_2=-2.0_idx*Pn_2-a_2*(Lambda_C+lambda*ratio_1)
     MatA_2(iy)=-1.0_idx*dt*Pn_2
     MatB_2(iy)=1.0_idx-0.5_idx*dt*Qn_2
     MatC_2(iy)=0.0_idx
     MatD_2(iy)=T_2(iy)+dt*a_2*lambda*ratio_1*T_1(iy)+dt*a_2*F_1*F(iy)+dt*Pn_2*T_2(iy-1)+0.5_idx*dt*Qn_2*T_2(iy) &
          & -dt*a_2*lambda*ratio_1*noise_theta1(iy)-dt*a_2*lambda*ratio_2*noise_theta2(iy)&
          & +dt*a_2*(noise_theta_dyn(iy)+noise_theta2(iy))/gamma

#ifdef notherm     
     MatD_2(iy)=T_2(iy)+dt*a_2*lambda*ratio_1*(T_1_mean(iy)-FoA/B)+dt*a_2*F_1*F(iy)+dt*Pn_2*T_2(iy-1)+0.5_idx*dt*Qn_2*T_2(iy) &
          & -dt*a_2*lambda*ratio_2*noise_theta2(iy) &
          & +dt*a_2*(noise_theta_dyn(iy)+noise_theta2(iy))/gamma
#endif
     ! Solve tridiagonal matrix problem
     call solve_tri_real(ny,MatA_1,MatB_1,MatC_1,MatD_1,T_1_new)
     T_1=T_1_new ! Update
     ! Solve tridiagonal matrix problem
     call solve_tri_real(ny,MatA_2,MatB_2,MatC_2,MatD_2,T_2_new)
     T_2=T_2_new ! Update

     ! Update T_1,T_2,and noise array for calculation of psi_1,psi_2
     do i=n_record_1,2,-1
        tauy_record_to1(1:ny,i)=tauy_record_to1(1:ny,i-1)
     end do
     tauy_record_to1(1:ny,1)=tauy(1:ny)
     do i=n_record_2,2,-1
        tauy_record_to2(1:ny,i)=tauy_record_to2(1:ny,i-1)
     end do
     tauy_record_to2(1:ny,1)=tauy(1:ny)
     if (mod(itime,istep_out).eq. 0) then
        write(*,*) "T=",itime,maxval(tauy),minval(tauy),irec
        do iy = 1,ny
           T_1_out(iy)=real(T_1(iy)+FoA/B,4)
           T_2_out(iy)=real(T_2(iy)+FoA/B,4)
           psi_1_out(iy)=real(psi_1(iy),4)
           psi_2_out(iy)=real(psi_2(iy),4)
           theta_out(iy)=real(theta_out(iy),4)

        end do
        tau=get_tau(ny,real(y,4),real(tauy,4)) ! Added
        ! Put variable
        call CHECK_W(NF90_PUT_VAR(ncid, T1_varid,T_1_out,start=(/1,irec/),count=(/ny,1/)) )
        call CHECK_W(NF90_PUT_VAR(ncid, T2_varid,T_2_out,start=(/1,irec/),count=(/ny,1/)) )
        call CHECK_W(NF90_PUT_VAR(ncid, PSI1_varid,PSI_1_out,start=(/1,irec/),count=(/ny,1/)) )
        call CHECK_W(NF90_PUT_VAR(ncid, PSI2_varid,PSI_2_out,start=(/1,irec/),count=(/ny,1/)) )
        call CHECK_W(NF90_PUT_VAR(ncid, theta_varid,theta,start=(/1,irec/),count=(/ny,1/)) )
        call CHECK_W(NF90_PUT_VAR(ncid, tau_varid,tau,start=(/1,irec/),count=(/ny,1/)) )
        ! Close file
        irec=irec+1
     end if
  end do
  call CHECK_W(NF90_CLOSE(ncid))

  deallocate(y)
  deallocate(F);deallocate(spat)
  deallocate(T_1);deallocate(T_1_new)
  deallocate(T_2);deallocate(T_2_new)
  deallocate(tauy_record_to1) ; deallocate(tauy_record_to2)
  deallocate(weight_record_1) ; deallocate(weight_record_2)
  deallocate(psi_1);   deallocate(psi_2)
  deallocate(theta)
  deallocate(MatA_1);deallocate(MatB_1)
  deallocate(MatA_2);deallocate(MatB_2)
  deallocate(T_1_out); deallocate(T_2_out)
  deallocate(psi_1_out);   deallocate(psi_2_out)
  deallocate(theta_out)
  deallocate(noise_theta_dyn)
  deallocate(noise_theta1);  deallocate(noise_theta2)
contains
  function get_weight(n,dt,to) result(weight)
    implicit none
    integer,parameter:: idx=8
    integer,intent(in) :: n
    real(idx),intent(in) :: dt,to
    real(idx) :: weight(n)
    integer :: i
    if (n .eq. 2) then
       weight(1)=0.5_idx*(2.0_idx*dt-to)*to/(dt*to)
       weight(2)=0.5_idx*to*to/(dt*to)
    else if (n .eq. 3) then
       weight(1)=0.5_idx*dt/to
       weight(2)=0.5_idx*dt/to+0.5_idx*(3.0_idx*dt-to)*(to-dt)/(dt*to)
       weight(3)=0.5_idx*(to-dt)*(to-dt)/(dt*to)
    else if (n .ge. 4) then
       weight(1)=0.5_idx*dt/to
       do i =2,n-2
          weight(i)=1.0_idx*dt/to
       end do
       weight(n-1)=0.5_idx*dt/to+0.5_idx*(real(n,idx)*dt-to)*(to-real(n-2,idx)*dt)/(dt*to)
       weight(n)=0.5_idx*(to-real(n-2,idx)*dt)*(to-real(n-2,idx)*dt)/(dt*to)
    end if
  end function get_weight
  subroutine solve_tri_real(N,A_in,B_in,C_in,D_in,U)
    ! Subroutine for solve tridiagonal
    ! Solve A * U_new(i-1) + B * U_new(i) + C * U_new(i+1) = D
    implicit none
    integer,intent(in) :: N
    real(idx),intent(in) :: A_in(N),B_in(N),C_in(N),D_in(N)
    real(idx),intent(inout) :: U(N)
    real(idx) :: A(N),B(N),C(N),D(N)
    integer :: i
    real(idx) :: m
    A = A_in ; B = B_in ; C = C_in; D = D_in
    ! initialize array
    U(1:N) = 0.0_idx
    do i=2,N
       m = A(i) / B(i-1)
       B(i) = B(i) - m * C(i-1)
       D(i) = D(i) - m * D(i-1)
       U(i) = 0.0_idx
    end do
    U(N) = D(N) / B(N)
    do i = N-1,1,-1
       U(i) = (D(i)-C(i)*U(i+1)) / B(i)
    end do
  end subroutine solve_tri_real
  ! Subroutine for calculate wind stress
    function get_tau(N,y,tau_y) result(tau)
    implicit none
    integer :: N
    real(4) :: y(N),tau_y(N),tau(N)
    integer :: iy
    iy=1
    tau(iy)=0.0_4
    do iy = 2,ny
       tau(iy)=tau(iy-1)+0.5_4*(y(iy)-y(iy-1))*(tau_y(iy)+tau_y(iy-1))
    end do
  end function get_tau
  subroutine CHECK_W(status)
    use NETCDF
    implicit none
    integer, intent ( in) :: status
    if(status /= NF90_NOERR) then 
       print *, trim(NF90_STRERROR(status))
       stop "Stopped"
    end if
  end subroutine CHECK_W  
end program solve_BCSmodel_Kohyama2021
