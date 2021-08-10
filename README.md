# Source code of idealized model used in Kohyama et al. (2021)
Idealized model of BCS used in Kohyama et al. (2021), which is primarily based upon Gallego and Cessi (2001)

Usage

# Requirements
: Fortran compiler and netCDF4 library

You may download the netCDF4 library it from https://www.unidata.ucar.edu/downloads/netcdf/

! Execution method

1. Specify netCDF library
(Example): 
export FPATH_FORT=-I/Users/kido/libs/netcdf-fortran-4.4.4/include
export LPATH_FORT=-L/Users/kido/libs/netcdf-fortran-4.4.4/lib
2. Compile the source code
ifort -o exec_ctl.out solve_BCSmodel_Kohyama2021.f90  $FPATH_FORT  $LPATH_FORT   -lnetcdff -lnetcdf -cpp
3. Execute the model after editing the namelist file
./exec_ctl.out < param_withnoise_gc01_ctl.nml
  
  
! Optional runs

(1) nodyn: Atmospheric temperature used to calculate the wind stress is replaced by the climatological value
(the name of external file is specified by fname_mean)
ifort -o exec_ctl.out solve_BCSmodel_Kohyama2021.f90  $FPATH_FORT  $LPATH_FORT   -lnetcdff -lnetcdf -nodyn -cpp
./exec_nodyn.out < param_withnoise_gc01_nodyn.nml

(2) notherm : temperature of other basin does not affect the another basin
ifort -o exec_ctl.out solve_BCSmodel_Kohyama2021.f90  $FPATH_FORT  $LPATH_FORT   -lnetcdff -lnetcdf -notherm -cpp
./exec_notherm.out < param_withnoise_gc01_notherm.nml

If you have any questions or comments, please contact to skido(at)jamstec.go.jp (replace at with @).
