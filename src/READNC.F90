!-----------------------------------------------------------------------
! Read a (lon,lat,time) netCDF data file
!-----------------------------------------------------------------------
subroutine READNC(Ncols,Nrows,Ntime,varname,var)

use netcdf

implicit none

character(len=*), intent(in) :: &
  varname             ! Variable name

integer, intent(in) :: &
  Ncols,             &! Number of columns in grid
  Nrows,             &! Number of rows in grid 
  Ntime               ! Number of timesteps

real, intent(out) :: &
  var(Ncols,Nrows,Ntime) ! Variable on grid

integer :: &
  ncid,              &! NetCDF dataset ID
  varid,             &! NetCDF variable ID
  status              ! NetCDF error status

status = nf90_open('met/'//varname//'.nc',nf90_nowrite,ncid)
status = nf90_inq_varid(ncid,varname,varid)
status = nf90_get_var(ncid,varid,var)
status = nf90_close(ncid)

end subroutine READNC

