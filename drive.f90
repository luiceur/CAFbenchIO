program drive

  !  PURPOSE
  !    Testing MPI/IO and netCDF IO
  !  DESCRIPTION
  !    Timing output of MPI/IO (cgca_pswci2) against
  !    the netcCDF version (cgca_pswci3).
  !  AUTHOR
  !    Luis Cebamanos, Anton Shterenlikht
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
  use cgca_m2hdf5
  use cgca_m2phys
  use benchclock


  implicit none

  integer, parameter :: maxlen = 64
  real,parameter :: gigabyte=real(2**30), resolution=1.0e-5,             &
       loge2 = log(real(2))
  logical(kind=ldef),parameter :: yesdebug = .true., nodebug = .false.

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

  integer( kind=ilrg ) :: icells, mcells

  !#################################
  character*(maxlen) :: filename
  integer, parameter :: totdim = 4, arrdim = totdim-1, coardim = 3
  integer, dimension(coardim) :: coarsize, copos
  character*(maxlen), dimension(3)  :: iolayername
  integer ::  ierr=0, i, j, k

! ! Add trailing blanks to keep all elements of the array of the same
! ! length. Max stripe count on ARCHER is 56.
!   character( len=maxlen), dimension(1) :: stripe_count = (/            &
!     "-c-1 " /)
!   character( len=maxlen), dimension(1) :: stripe_size = (/             &
!     "-S1m " /)
!   character( len=2*maxlen ) :: dir

! Add trailing blanks to keep all elements of the array of the same
! length. Max stripe count on ARCHER is 56.
  character( len=maxlen), dimension(9) :: stripe_count = (/            &
    "-c-1 ", "-c0  ", "-c1  ", "-c4  ", "-c8  ",                       &
    "-c16 ", "-c20 ", "-c32 ", "-c40 " /)
  character( len=maxlen), dimension(7) :: stripe_size = (/             &
    "-S1m ", "-S2m ", "-S4m ", "-S8m ", "-S16m",                       &
    "-S32m", "-S64m" /)
  character( len=2*maxlen ) :: dir
  character( len=120 ) :: errmsg
 
  !#################################
  double precision :: t0, t1, tdiff, fsizeb, fsizeg
  integer :: i1,i2,i3,j1,j2,j3,l1,l2,errstat

!*********************************************************************72
  ! first executable statement
  iolayername(1) = 'mpiio.dat'
  iolayername(2) = 'netcdf.dat'
  iolayername(3) = 'hdf5.dat'

    img = this_image()
  nimgs = num_images()

  ! start MPI
  call MPI_Init(ierr)

  ! do a check on image 1
  if ( img .eq. 1 ) then
     write (*,'(a,i0,a)') "running on ", nimgs, " images in a 3D grid"
  end if

  sync all

  ! each image calculates the coarray grid dimensions
  call cgca_gdim( nimgs, ir, qual )
  
  !Create the local coarray size
  c(1) = 87
  c(2) = 87
  c(3) = 87


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

  
  ! allocate space coarray with a single layer
  call cgca_as(1, c(1), 1, c(2), 1, c(3), 1, ir(1), 1, ir(2), 1, 1, space)

  l1 =c(1)*ir(1)
  l2= c(2)*ir(2)

  !Initialize coarrays
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

  ! Loop over lfs stripe counts
  do i = 1, size( stripe_count )

    ! Loop over lfs stripe sizes
    do j = 1, size( stripe_size )

      dir = "lfs" // trim( stripe_count(i) ) // trim( stripe_size(j) )

      ! Image 1 makes a dir with desired lfs settings
      if ( img .eq. 1 ) then

        ! Make the dir
        errmsg = ""
        
        call execute_command_line( command = "mkdir " // trim(dir),    &
          wait = .true. , exitstat = ierr, cmdstat = errstat,          &
          cmdmsg = errmsg ) 
        if ( ierr .ne. 0 ) error stop

        
        ! Set lfs parameters
        call execute_command_line( command = "lfs setstripe " //       &
          stripe_count(i) // " " // stripe_size(j) // " " //           &
          trim(dir) )

      end if

      ! Loop over IO layers
      do k = 1, size( iolayername )

        filename = trim(dir) // "/" // trim( iolayername(k) )

        sync all

        t0 = benchtime()

        if ( k .eq. 1 ) then
          ! MPI/IO
          call cgca_pswci2( space, cgca_state_type_grain, filename )
        else if ( k .eq. 2 ) then
          ! NetCDF
           call cgca_pswci3( space, cgca_state_type_grain, filename )
       else if ( k .eq. 3 ) then
          ! HDF5
           call cgca_pswci4( space, cgca_state_type_grain, filename )
        end if

        t1 = benchtime()

        sync all

        if (img .eq. 1)  then
          tdiff = t1 - t0
          write (*,*) trim( iolayername(k) ), " ",                     &
            trim( stripe_count(i) ), " ", trim( stripe_size(j) ), " ", &
            fsizeg/tdiff
!            tdiff, "s, rate: ", fsizeg/tdiff, "GB/s."
        end if

        sync all

      end do
    end do
  end do

  ! terminate MPI
  call MPI_Finalize(ierr)


  ! Write HDF5
  !call cgca_pswci4( space, cgca_state_type_grain, iolayername(3) )

  ! terminate MPI
  call MPI_Finalize(ierr)
  
  ! deallocate all arrays
  call cgca_ds(space)
  
end program drive

