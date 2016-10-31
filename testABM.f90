
program testACF

!  PURPOSE
!    Benchmarking MPI/IO and netCDF IO with CAF
!  DESCRIPTION
!    Timing output of MPI/IO (cgca_pswci2) against
!    the netcCDF version (cgca_pswci3).
!  AUTHOR
!    Luis Cebamanos
!  COPYRIGHT
!    See LICENSE
!  USES
!    cgca testaux
!  USED BY
!    Part of CGPACK test suite
!  SOURCE

use cgca_m1co
use cgca_m2alloc
use cgca_m2netcdf
use cgca_m2mpiio
use cgca_m2phys

implicit none



integer, parameter :: maxlen = 64
real,parameter :: gigabyte=real(2**30), resolution=1.0e-5,             &
 loge2 = log(real(2))


real( kind=rdef ) ::    &
 qual,                  & ! quality
 bsz0(3),               & ! the given "box" size
 bsz(3),                & ! updated "box" size
 dm,                    & ! mean grain size, linear dim, phys units
 lres,                  & ! linear resolution, cells per unit of length
 res                      ! resolutions, cells per grain

integer( kind=idef ) :: ir(3), nimgs, img, ng, istriping
integer( kind=iarr ), allocatable :: space(:,:,:,:)[:,:,:]
integer( kind=iarr ) :: c(3)  ! coarray dimensions

integer( kind=ilrg ) :: icells, mcells,repetition,irep

!#################################
character*(maxlen) :: filename
integer, parameter :: numiolayer = 4
integer, parameter :: totdim = 4, arrdim = totdim-1, coardim = 3
integer, parameter :: numstriping = 3
character*(maxlen), dimension(numstriping) :: stripestring
character*(maxlen), dimension(numiolayer)  :: iostring, iolayername
integer ::  comm, ierr=0, rank=0, mpisize=0, filetype,      &
     mpi_subarray, fh, funit
integer, dimension(totdim) :: asizehal
integer, dimension(arrdim) :: arrsize, arstart, artsize
integer, dimension(coardim) :: coarsize, copos

!#################################
real :: time1, time2, fsizeb, fsizeg, tdiff
integer :: i1,i2,i3,j1,j2,j3,l1,l2

!*********************************************************************72

! first executable statement
iolayername(1) = 'serial.dat'
iolayername(2) = 'mpiio.dat'
iolayername(3) = 'hdf5.dat'
iolayername(4) = 'netcdf.dat'

stripestring(1) = 'unstriped'
stripestring(2) = 'striped'
stripestring(3) = 'defstriped'

! physical dimensions of the box, assume mm
bsz0 = (/ 2.0, 3.0, 3.0 /)

! mean grain size, linear dimension, e.g. mean grain diameter, also mm
dm = 1.0e-1
!dm = 1.0e0
! resolution
res = 1.0e5

    img = this_image()
nimgs = num_images()

if (img ==1) then
   write (*,'(a,i0,a)') "running on ", nimgs, " images in a 3D grid"
end if

sync all

! each image calculates the coarray grid dimensions
call cgca_gdim( nimgs, ir, qual )

c(1) = 58
c(2) = 174
c(3) = 174


! total number of cells in a coarray
icells = int( c(1), kind=ilrg ) * int( c(2), kind=ilrg ) *             &
         int( c(3), kind=ilrg )

! total number of cells in the model
mcells = icells * int( nimgs, kind=ilrg )

if ( img .eq. 1 ) then
  write ( *, "(9(a,i0),tr1,g10.3,tr1,g10.3,3(a,g10.3),a)" )            &
    "img: ", img  , " nimgs: ", nimgs, " (", c(1) ,                    &
    ","    , c(2) , ","       , c(3) , ")[", ir(1),                    &
    ","    , ir(2), ","       , ir(3), "] ", ng   ,                    &
    qual, lres,                                                        &
    " (", bsz(1), ",", bsz(2), ",", bsz(3), ")"
  write (*,'(a,i0,a)') "Each image has ",icells, " cells"
  write (*,'(a,i0,a)') "The model has ", mcells, " cells"
end if

! Total output file size, in B and in GB.
fsizeb = real( mcells * storage_size( space, kind=ilrg ) / 8_ilrg )
fsizeg = fsizeb / gigabyte

call cgca_as(1, c(1), 1, c(2), 1, c(3), 1, ir(1), 1, ir(2), 1, 1, space)

! initialise coarray to image number
space = int( img, kind=iarr )

l1 =c(1)*ir(1)
l2= c(2)*ir(2)

copos(:) = this_image( space )
 do i3 = 1, c(3)
     do i2 = 1, c(2)
        do i1 = 1, c(1)
           
           j1 = (copos(1)-1)*c(1) + i1
           j2 = (copos(2)-1)*c(2) + i2
           j3 = (copos(3)-1)*c(3) + i3

           space(i1,i2,i3,cgca_state_type_grain) = (j3-1)*l1*l2 + (j2-1)*l1 + j1

        end do
     end do
  end do


! start MPI 
call MPI_Init(ierr)

repetition = 5

do istriping = 1, numstriping
   do irep = 1, repetition
      ! dump the model, MPI/IO
      !   call cpu_time( time1 )
      time1 = MPI_Wtime()
      filename = trim(stripestring(istriping))//'/'//trim(iolayername(2))
      call cgca_pswci2( space, cgca_state_type_grain, filename )
      time2 = MPI_Wtime()
      !   call cpu_time( time2 )
      tdiff = time2-time1
      if (img .eq. 1)write (*,*) "MPI-IO-",trim(stripestring(istriping)),": ", &
           tdiff, "s, rate: ", fsizeg/tdiff, "GB/s."
      ! dump the model, netCDF
      time1 = MPI_Wtime()
      !  call cpu_time( time1 )   
      filename = trim(stripestring(istriping))//'/'//trim(iolayername(4))
      call cgca_pswci3( space, cgca_state_type_grain, filename )
      time2 = MPI_Wtime()
      !   call cpu_time( time2 )
      tdiff = time2-time1
      if (img .eq. 1) write (*,*) "netCDF-",trim(stripestring(istriping)),": ", &
           tdiff, "s, rate: ", fsizeg/tdiff, "GB/s."
      
      !   call fdelete(filename)
      end do
end do


! terminate MPI
call MPI_Finalize(ierr)


! deallocate all arrays
call cgca_ds(space)

end program testACF

subroutine fdelete(filename)

  implicit none

  character *(*) :: filename
  integer, parameter :: iounit = 15
  integer :: stat

  open(unit=iounit, iostat=stat, file=filename, status='old')
  if (stat.eq.0) close(unit=iounit, status='delete')

end subroutine fdelete
!*roboend*
