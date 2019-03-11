! COPYRIGHT (c) 2010,2011 Science and Technology Facilities Council (STFC)
!
! Written by: Sue Thorne, Jonathan Hogg
!
! 9 Feburary 2010 Version 1.0.0
!     Brought under library versioning after definition of library standards
!     for MATLAB interfaces.
!
! *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!
!    GENERIC MEX INTERFACE (for HSL_MI20)
!
! *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!
!  Given an array A of length n form a heap
!
!   [ A, inform, INDA ]
!     = hsl_kb12('build', A, INDA)
!
!  Find the smallest entry in the current heap
!   [ A, inform, INDA]
!     = hsl_kb12('get', m, A, INDA)
!
!  Usual Input -
!    A: n-vector A
!    m: number of entries in A that lie in the heap on entry
!
!  Optional Input -
!    INDA, integer n-vector on which exactly the same permutation as A is
!       applied
!
!  Usual Output -
!   A: the new heap
!
!  Optional Output -
!   inform: an integer containing value 0 for successful exit
!   INDA: permuted INDA
!
! *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!

module mi20_matlab_main
   use hsl_matlab
   use hsl_zd11_double
   use hsl_mi20_double
   use hsl_mc65_double
   implicit none

   integer, parameter :: wp = kind(0d0)

   integer(mws_), PARAMETER :: ncontrol = 16
   CHARACTER ( LEN = 21 ), PARAMETER :: fcontrol( ncontrol ) = (/          &
       'aggressive           ', 'c_fail               ',                   &
       'max_levels           ', 'max_points           ',                   &
       'reduction            ', 'st_method            ',                   &
       'st_parameter         ',                                            &
       'coarse_solver        ', 'coarse_solver_its    ',                   &
       'damping              ', 'err_tol              ',                   &
       'pre_smoothing        ',                                            &
       'post_smoothing       ', 'smoother             ',                   &
       'v_iterations         ', 'one_pass_coarsen     '/)

   integer(mws_), PARAMETER :: ninform = 4
   CHARACTER ( LEN = 21 ), PARAMETER :: finform( ninform ) = (/             &
      'flag                 ', 'clevels              ',                   &
      'cpoints              ', 'cnnz                 ' /)
   !   'ma48_ainfo           ', 'ma48_finfo           ',                   &
   !   'ma48_sinfo           ', 'getrf_info           ',                   &
contains
   subroutine mi20_matlab_control(nlhs, plhs, nrhs, prhs)
      integer(int4_) :: nlhs, nrhs
      integer(mwp_) :: plhs(nlhs), prhs(nrhs)

      type(mi20_control) :: control

      ! Output arguments
      integer, parameter :: c_arg = 1

      integer(mwp_) :: control_out

      IF ( nlhs > 1 ) call MATLAB_error( &
         'hsl_mi20 provides at most 1 output argument for this call' )

      ! Create control output structure
      plhs(1) = MATLAB_create_structure(fcontrol)
      control_out = plhs(1)

      ! Copy default controls
      call matlab_set_field( control_out, 'aggressive', control%aggressive )
      call matlab_set_field( control_out, 'c_fail', control%c_fail )
      call matlab_set_field( control_out, 'max_levels', control%max_levels )
      call matlab_set_field( control_out, 'max_points', control%max_points )
      call matlab_set_field( control_out, 'one_pass_coarsen', &
         control%one_pass_coarsen )
      call matlab_set_field( control_out, 'reduction', control%reduction )
      call matlab_set_field( control_out, 'st_method', control%st_method )
      call matlab_set_field( control_out, 'st_parameter', control%st_parameter )
      call matlab_set_field( control_out, 'coarse_solver', &
         control%coarse_solver )
      call matlab_set_field( control_out, 'coarse_solver_its', &
         control%coarse_solver_its )
      call matlab_set_field( control_out, 'damping', control%damping )
      call matlab_set_field( control_out, 'err_tol', control%err_tol )
      !call matlab_set_field( control_out, 'levels', control%levels )
      call matlab_set_field( control_out, 'pre_smoothing', &
         control%pre_smoothing )
      call matlab_set_field( control_out, 'post_smoothing', &
         control%post_smoothing )
      call matlab_set_field( control_out, 'smoother', control%smoother )
      call matlab_set_field( control_out, 'v_iterations', control%v_iterations )
   end subroutine mi20_matlab_control

   subroutine mi20_matlab_setup(nlhs, plhs, nrhs, prhs, A, coarse_data, keep, &
         control, errmsg)
      integer(int4_) :: nlhs, nrhs
      integer(mwp_) :: plhs(nlhs), prhs(nrhs)
      type(zd11_type), intent(out) :: A
      type(mi20_data), dimension(:), allocatable, intent(out) :: coarse_data
      type(mi20_keep), intent(out) :: keep
      type(mi20_control), intent(out) :: control
      character(200), intent(out) :: errmsg

      ! Input argument positions
      integer, parameter :: a_arg = 1, &
                            c_arg = 2
      ! Output argument positions
      integer, parameter :: i_arg   = 1, &
                            cd_arg  = 2, & ! NOTE: This is already set by caller
                            cd2_arg = 3, & ! Set to 0, deprecated
                            k_arg   = 4, & ! Set to 0, deprecated
                            k1_arg  = 5, & ! Set to 0, deprecated
                            k2_arg  = 6    ! Set to 0, deprecated

      type(zd11_type) :: B
      type(mi20_info) :: info

      integer :: i, j, l, infoo
      integer(int4_) :: fnum

      character(80) :: fname
      integer(mwp_) :: a_in, c_in, cn_in, val_pr
      integer(mws_) :: m, n, temp_mws

      integer :: iores
      character ( len = 80 ) :: output_unit, filename
      logical ::  filexx, opened

      errmsg = "No errors"
      if(nrhs < c_arg) then
         errmsg = "Insufficient input arguments"
         return
      endif
      if(nrhs < c_arg) then
         errmsg = "Too many input arguments"
         return
      endif
      if(nlhs > 6) then
         errmsg = "Too many output arguments"
         return
      endif
      if(nlhs > 2) then
         ! old style interface, print deprecation warning then set
         ! unwanted arguments to 0
         call MATLAB_error( &
            "Discontinued hsl_mi20 interface being used. See help for further &
            &information.")
         do i = cd_arg+1, nlhs
            plhs(i) = fortran_to_matlab(0.0_wp)
         end do
      endif

      control%error = -1
      control%print = -1
      control%print_level = 0

      ! Setup info data type
      ! Create info output structure
      plhs( i_arg ) = MATLAB_create_structure(finform)

      !  Open i/o units
      IF ( control%error > 0 ) THEN
         WRITE( output_unit, "( I0 )" ) control%error
         filename = "output_qpc." // TRIM( output_unit )
         INQUIRE( FILE = filename, EXIST = filexx )
         IF ( filexx ) THEN
            OPEN( control%error, FILE = filename, FORM = 'FORMATTED',         &
               STATUS = 'OLD', IOSTAT = iores )
         ELSE
            OPEN( control%error, FILE = filename, FORM = 'FORMATTED',         &
               STATUS = 'NEW', IOSTAT = iores )
         END IF
      END IF

      IF ( control%print > 0 ) THEN
         INQUIRE( control%print, OPENED = opened )
         IF ( .NOT. opened ) THEN
            WRITE( output_unit, "( I0 )" ) control%print
            filename = "output_qpc." // TRIM( output_unit )
            INQUIRE( FILE = filename, EXIST = filexx )
            IF ( filexx ) THEN
               OPEN( control%print, FILE = filename, FORM = 'FORMATTED', &
                  STATUS = 'OLD', IOSTAT = iores )
            ELSE
               OPEN( control%print, FILE = filename, FORM = 'FORMATTED', &
                  STATUS = 'NEW', IOSTAT = iores )
            END IF
         END IF
      END IF

      ! Check A is of correct type
      a_in = prhs(a_arg)

      ! Read A
      CALL ZD11_put( B%type, 'general', stat=infoo )
      IF ( MATLAB_is_sparse( a_in ) ) THEN

         call matlab_to_fortran(a_in, m, n, B%ptr, B%col, B%val, 'A')
         B%m = m
         B%n = n
         B%ne = B%ptr(n+1)-1

      ELSE
         !  Get the numbers of variables and general constraints
         m = MATLAB_get_m( a_in )
         n = MATLAB_get_n( a_in )
         B%m = m
         B%n = n
         B%ne = m * n

         !  Allocate space for the input matrix
         ALLOCATE( B%col( B%ne ), B%val( B%ne ), B%ptr(B%n +1), STAT = infoo )
         IF (infoo>0) THEN
            info%flag=-10; GOTO 1000
         END IF

         !  Set the row and column indices if the matrix is dense
         l = 0
         DO j = 1, n
            DO i = 1, m
               l = l + 1
               B%col( l ) = i
            END DO
         END DO

         DO j = 1, n+1
            B%ptr( j ) = (j-1)*n + 1
         END DO

         ! Copy the real components of A
         val_pr = MATLAB_get_ptr( a_in )
         temp_mws = B%ne
         CALL MATLAB_copy_from_ptr( val_pr, B%val, temp_mws )
      END IF

      ! Currently stored B = A^T
      ! Take transpose A = B^T
      CALL mc65_matrix_transpose(B,A,infoo,pattern=.false.)
      IF ( infoo < 0 ) THEN
         IF ( infoo == mc65_err_memory_alloc) THEN
            info%flag = -10
         ELSE IF ( infoo == mc65_err_memory_dealloc) THEN
            info%flag = -11
         END IF
         GOTO 1000
      END IF
      CALL mc65_matrix_destruct(B,infoo)
      IF ( infoo == mc65_err_memory_dealloc) THEN
         info%flag = -11
        GOTO 1000
      END IF

      ! Read 'control' parameter
      c_in = prhs( c_arg )
      IF ( .NOT. MATLAB_is_structure( c_in ) ) THEN
         CALL mc65_matrix_destruct(A,infoo)
         CALL MATLAB_error("argument 'control' should be a structure")
      END IF
      do fnum = 1, MATLAB_get_no_fields( c_in )
         fname = MATLAB_get_field_name_by_no(c_in, fnum)
         select case(trim(fname))
         case('aggressive')
            CALL MATLAB_get_value(c_in, fname, cn_in, control%aggressive)
         case('c_fail')
            CALL MATLAB_get_value(c_in, fname, cn_in, control%c_fail)
         case('max_levels')
            CALL MATLAB_get_value(c_in, fname, cn_in, control%max_levels)
         case('max_points')
            CALL MATLAB_get_value(c_in, fname, cn_in, control%max_points)
         case('reduction')
            CALL MATLAB_get_value(c_in, fname, cn_in, control%reduction)
         case('st_method')
            CALL MATLAB_get_value(c_in, fname, cn_in, control%st_method)
         case('st_parameter')
            CALL MATLAB_get_value(c_in, fname, cn_in, control%st_parameter)
         case('coarse_solver')
            CALL MATLAB_get_value(c_in, fname, cn_in, control%coarse_solver)
         case('coarse_solver_its')
            CALL MATLAB_get_value(c_in, fname, cn_in, &
               control%coarse_solver_its)
         case('damping')
            CALL MATLAB_get_value(c_in, fname, cn_in, control%damping)
         case('err_tol')
            CALL MATLAB_get_value(c_in, fname, cn_in, control%err_tol)
         !case('levels')
         !   CALL MATLAB_get_value(c_in, 'levels', cn_in, control%levels)
         case('pre_smoothing')
            CALL MATLAB_get_value(c_in, fname, cn_in, control%pre_smoothing)
         case('post_smoothing')
            CALL MATLAB_get_value(c_in, fname, cn_in, control%post_smoothing)
         case( 'smoother' )
            CALL MATLAB_get_value(c_in, fname, cn_in, control%smoother)
         case( 'v_iterations' )
            CALL MATLAB_get_value(c_in, fname, cn_in, control%v_iterations)
         case( 'one_pass_coarsen' )
            CALL MATLAB_get_value(c_in, fname, cn_in, control%one_pass_coarsen)
         case  default
            CALL MATLAB_warning( 'control has unrecognized entry' )
         end select
      end do

      CALL MI20_setup(A, coarse_data, keep, control, info)
      IF ( info%flag < 0 ) GOTO 1000

      ! Copy information
      call matlab_set_field( plhs(i_arg), 'clevels', info%clevels)
      call matlab_set_field( plhs(i_arg), 'cpoints', info%cpoints)
      call matlab_set_field( plhs(i_arg), 'cnnz', info%cnnz)

      ! Check for error flag
      IF ( info%flag > 0) THEN
         SELECT CASE (info%flag)
         CASE ( 1 )
            CALL MATLAB_warning( &
               'Method used to find strong transpose connections changed')
         CASE (10, 11, 12, 13)
            CALL MATLAB_warning( 'Coarsening terminated prematurely' )
         CASE ( 20 )
            CALL MATLAB_warning( &
               'Number of requested levels greater than available levels')
         END SELECT
      END IF

      ! Copy the components of info

      1000 CONTINUE

      call matlab_set_field( plhs(i_arg), 'flag', info%flag)

      IF (info%flag <0) THEN
         call print_error_flag(info, errmsg)
         call mi20_finalize(coarse_data, keep, control, info)
         return
      ENDIF

      !  close any opened io units
      IF ( control%error > 0 ) THEN
         INQUIRE( control%error, OPENED = opened )
         IF ( opened ) CLOSE( control%error )
      END IF

      IF ( control%print > 0 ) THEN
         INQUIRE( control%print, OPENED = opened )
         IF ( opened ) CLOSE( control%print )
      END IF
   end subroutine mi20_matlab_setup

   subroutine mi20_matlab_precondition(nlhs, plhs, nrhs, prhs, A, coarse_data, &
         keep, control)
      integer(int4_) :: nlhs, nrhs
      integer(mwp_) :: plhs(nlhs), prhs(nrhs)
      type(zd11_type), intent(in) :: A
      type(mi20_data), dimension(:), allocatable, intent(in) :: coarse_data
      type(mi20_keep), intent(inout) :: keep
      type(mi20_control), intent(in) :: control

      ! Input arguments
      integer, parameter :: z_arg   = 1, &
                            cd_arg  = 2, & ! Already handled by caller
                            cd2_arg = 3, & ! Ignored, deprecated
                            c_arg   = 4, & ! Ignored, deprecated
                            kk1_arg = 5, & ! Ignored, deprecated
                            kk2_arg = 6    ! Ignored, deprecated
      ! Output arguments
      integer, parameter :: x_arg   = 1, &
                            i_arg   = 2, &
                            k_arg   = 3, & ! Set to 0, deprecated
                            k1_arg  = 4, & ! Set to 0, deprecated
                            k2_arg  = 5    ! Set to 0, deprecated

      type(mi20_info) :: info
      real(wp), dimension(:), allocatable :: x, z
      integer :: i, infoo

      integer, dimension(:), allocatable :: z1p, z1r
      real(wp), dimension(:), allocatable :: z1

      integer(mws_) :: m, n, ne
      integer(mwp_) :: z_in, z_pr

      character ( len = 80 ) :: output_unit, filename, errmsg
      logical ::  filexx, opened
      integer :: iores

      IF ( nrhs < cd_arg ) call MATLAB_error("Insufficient input arguments")
      IF ( nrhs > kk2_arg ) call MATLAB_error("Too many input arguments")
      IF ( nrhs > cd_arg ) then
         ! old style interface, print deprecation warning then set
         ! unwanted arguments to 0
         call MATLAB_error( &
            "Discontinued hsl_mi20 interface being used. See help for further &
            &information.")
         if( nlhs < k2_arg ) call MATLAB_error("Insufficient output arguments")
         if( nlhs > k2_arg ) call MATLAB_error("Too many output arguments")
         do i = i_arg+1, nlhs
            plhs(i) = fortran_to_matlab(0.0_wp)
         end do
      else
         ! new style interface
         if( nlhs < i_arg ) call MATLAB_error("Insufficient output arguments")
         if( nlhs > i_arg ) call MATLAB_error("Too many output arguments")
      endif

      ! Setup info data type
      ! Create info output structure
      plhs( i_arg ) = MATLAB_create_structure(finform)

      !  Open i/o units
      IF ( control%error > 0 ) THEN
         WRITE( output_unit, "( I0 )" ) control%error
         filename = "output_qpc." // TRIM( output_unit )
         INQUIRE( FILE = filename, EXIST = filexx )
         IF ( filexx ) THEN
            OPEN( control%error, FILE = filename, FORM = 'FORMATTED',         &
               STATUS = 'OLD', IOSTAT = iores )
         ELSE
            OPEN( control%error, FILE = filename, FORM = 'FORMATTED',         &
               STATUS = 'NEW', IOSTAT = iores )
         END IF
      END IF

      IF ( control%print > 0 ) THEN
         INQUIRE( control%print, OPENED = opened )
         IF ( .NOT. opened ) THEN
            WRITE( output_unit, "( I0 )" ) control%print
            filename = "output_qpc." // TRIM( output_unit )
            INQUIRE( FILE = filename, EXIST = filexx )
            IF ( filexx ) THEN
               OPEN( control%print, FILE = filename, FORM = 'FORMATTED', &
                  STATUS = 'OLD', IOSTAT = iores )
            ELSE
               OPEN( control%print, FILE = filename, FORM = 'FORMATTED', &
                  STATUS = 'NEW', IOSTAT = iores )
            END IF
         END IF
      END IF

      ! Get length of z and allocate space for input A
      z_in = prhs(z_arg)

      IF (.not. MATLAB_is_numeric(z_in)) THEN
         CALL MATLAB_error( 'There must be a vector z' )
      END IF

      IF (MATLAB_is_sparse(z_in)) THEN

         call matlab_to_fortran(z_in, n, m, z1p, z1r, z1, 'z')

         IF (m>1) THEN
            CALL MATLAB_error( 'z must only have one column' )
         END IF
         ne = z1p(m+1)-1

         ALLOCATE(z(n),STAT=infoo)
         IF (infoo>0) THEN
            info%flag=-10; GOTO 1000
         END IF
         z = 0.0_wp

         DO i=1,ne
            z(z1r(i)) = z1(i)
         END  DO

         DEALLOCATE (z1r,z1,STAT=infoo)
         IF (infoo>0) THEN
            info%flag=-11; GOTO 1000
         END IF

         DEALLOCATE (z1p,STAT=infoo)
         IF (infoo>0) THEN
            info%flag=-11; GOTO 1000
         END IF

      ELSE

         call matlab_to_fortran(z_in, z, n, 'z')

      END IF
      allocate(x(n), stat=infoo)
      IF (infoo>0) THEN
         info%flag=-10; GOTO 1000
      END IF

      CALL MI20_precondition(A, coarse_data, z, x, keep, control, info)

      IF ( info%flag < 0 ) GOTO 1000
      ! Copy information
      call matlab_set_field( plhs(i_arg), 'clevels', info%clevels)
      plhs(x_arg) = fortran_to_matlab(x(1:n))

      ! Copy the components of info

      1000 CONTINUE

      call matlab_set_field( plhs(i_arg), 'flag', info%flag)

      IF (info%flag <0) THEN
         call print_error_flag(info, errmsg)
         !call mi20_finalize(coarse_data, keep, control, info)
         call MATLAB_error(errmsg)
      ENDIF

      !  close any opened io units
      IF ( control%error > 0 ) THEN
         INQUIRE( control%error, OPENED = opened )
         IF ( opened ) CLOSE( control%error )
      END IF

      IF ( control%print > 0 ) THEN
         INQUIRE( control%print, OPENED = opened )
         IF ( opened ) CLOSE( control%print )
      END IF
   end subroutine mi20_matlab_precondition

   subroutine print_error_flag(info, errmsg)
      type(mi20_info), intent(in) :: info
      character(*), intent(out) :: errmsg

      errmsg = "Unknown error"

      SELECT CASE( info%flag)
      CASE ( - 2)
         errmsg =  'One or more diagonal entry missing'
      CASE ( - 3)
         errmsg =  'One or more diagonal entry is <= 0'
      CASE ( - 101)
         errmsg =  'control.st_parameter out of range'
      CASE ( - 102)
         errmsg =  'control.err_tol out of range'
      CASE ( - 103)
         errmsg =  'control.max_points out of range'
      CASE ( - 104)
         errmsg =  'control.st_method out of range'
      CASE ( - 105)
         errmsg =  'control.aggressive out of range'
      CASE ( - 106)
         errmsg =  'control.c_fail out of range'
      CASE ( - 107)
         errmsg =  'control.v_iterations out of range'
      CASE ( - 108)
         errmsg =  'control.smoother out of range'
      CASE ( - 109)
         errmsg =  'control.pre_smoothing out of range'
      CASE ( - 110)
         errmsg =  'control.post_smoothing out of range'
      CASE ( - 111)
         CALL MATLAB_error &
            ( 'control.pre_smoothing+control.post_smoothing = 0' )
      CASE ( - 112)
         errmsg =  'control.coarse_solver out of range'
      CASE ( - 113)
         errmsg =  'control.coarse_solver_its out of range'
      CASE ( - 115)
         errmsg =  'control.damping out of range'
      CASE ( - 116)
         errmsg =  'control.max_levels out of range'
      CASE ( - 10)
         errmsg =  'Memory allocation error'
      CASE ( - 11)
         errmsg =  'Memory deallocation error'
      CASE ( - 12)
         errmsg =  'Coarsening has failed'
      CASE ( - 14)
         errmsg = 'Increase in 2-norm of vector greater than control.err_tol'
      CASE ( - 15)
         errmsg = 'Call to precondition follows unsuccessful call to setup'
      CASE ( - 16)
         errmsg = 'Size of z is less than the order of A'
      CASE ( - 17)
         errmsg = 'Error return from _GETRF'
      CASE ( - 18)
         if (info%ma48_ainfo%flag<0) then
            errmsg = 'Error return from HSL_MA48 Analyse'
         else if (info%ma48_finfo%flag<0) then
            errmsg = 'Error return from HSL_MA48 Factorize'
         else
            errmsg = 'Error return from HSL_MA48 Solve'
         end if
      END SELECT
   end subroutine print_error_flag

   subroutine mi20_matlab_finalize(A, coarse_data, keep, control)
      type(zd11_type), intent(inout) :: A
      type(mi20_data), dimension(:), allocatable, intent(inout) :: coarse_data
      type(mi20_keep), intent(inout) :: keep
      type(mi20_control), intent(in) :: control

      type(mi20_info) :: info

      call mi20_finalize(coarse_data, keep, control, info)
   end subroutine mi20_matlab_finalize
end module mi20_matlab_main

! This module looks after a SAVEd set of variables mapping integer handles
! to Fortran keep and order variables
module mi20_handles
   use mi20_matlab_main
   implicit none

   type mi20_cd_container
      type(mi20_data), dimension(:), allocatable :: cd
   end type mi20_cd_container

   ! Data associated with the handle
   ! Considered to be empty if both associated(rkeep) and associated(ckeep)
   ! are .false.
   type mi20_hdl
      type(zd11_type), pointer :: A => null()
      type(mi20_cd_container), pointer :: cd_ptr => null()
      type(mi20_keep), pointer :: keep => null()
      type(mi20_control), pointer :: control => null()
   end type mi20_hdl

   ! How many handles initally and how much increase once exhausted
   integer, parameter :: initial_handles = 5
   double precision, parameter :: multiplier = 2.0

   ! SAVEd data
   integer, save :: next_handle = 1
   integer, save :: total_handles = 0
   type(mi20_hdl), dimension(:), allocatable, save :: handles

contains

   integer function mi20_new_handle()
      type(mi20_hdl), dimension(:), allocatable :: temp
      integer :: i

      ! Do we need to expand the number of available handles?
      if (next_handle .gt. total_handles) then
         if(total_handles.ne.0) then
            ! Need to expand existing handle selection
            allocate(temp(total_handles))
            do i = 1, total_handles
               temp(i)%A => handles(i)%A
               temp(i)%cd_ptr => handles(i)%cd_ptr
               temp(i)%keep => handles(i)%keep
               temp(i)%control => handles(i)%control
            end do
            deallocate(handles)
            total_handles = max(int(multiplier*total_handles), total_handles)
            allocate(handles(total_handles))
            do i = 1, size(temp)
               handles(i)%A => temp(i)%A
               handles(i)%cd_ptr => temp(i)%cd_ptr
               handles(i)%keep => temp(i)%keep
               handles(i)%control => temp(i)%control
            end do
            deallocate(temp)
         else
            ! First call since module loaded
            total_handles = initial_handles
            allocate(handles(total_handles))

            ! Register clean function
            call mexAtExit(cleanup_all_handles)
         endif
      endif

      mi20_new_handle = next_handle
      allocate(handles(next_handle)%A)
      allocate(handles(next_handle)%cd_ptr)
      allocate(handles(next_handle)%keep)
      allocate(handles(next_handle)%control)
      next_handle = next_handle + 1
   end function mi20_new_handle

   ! This routine is called at unload of this module from MATLAB.
   ! It shuold cleanup all SAVEd data
   subroutine cleanup_all_handles()
      integer :: i

      do i = 1, next_handle-1
         call cleanup_handle(i)
      end do
   end subroutine cleanup_all_handles

   ! Destroy the data associated with a handle.
   ! Recover all free pointers at end of handle list.
   subroutine cleanup_handle(handle)
      integer, intent(in) :: handle

      integer :: current
      integer :: st

      call mi20_matlab_finalize(handles(handle)%A, &
         handles(handle)%cd_ptr%cd, handles(handle)%keep, &
         handles(handle)%control)
      deallocate(handles(handle)%A, stat=st)
      nullify(handles(handle)%A)
      deallocate(handles(handle)%cd_ptr, stat=st)
      nullify(handles(handle)%cd_ptr)
      deallocate(handles(handle)%keep, stat=st)
      nullify(handles(handle)%keep)
      deallocate(handles(handle)%control, stat=st)
      nullify(handles(handle)%control)

      do current = handle, 1, -1
         if(current.ne.next_handle-1) exit
         if(.not.associated(handles(current)%A)) then
            ! Current "last" element is unallocated, make it next available
            ! element.
            next_handle = next_handle - 1
         endif
      end do
   end subroutine cleanup_handle

end module mi20_handles

SUBROUTINE mexFunction( nlhs_in, plhs, nrhs_in, prhs )
   use mi20_matlab_main
   use mi20_handles
   IMPLICIT NONE

   integer(int4_) :: nlhs_in, nrhs_in
   integer(mwp_) :: plhs( * ), prhs( * )

   integer :: handle
   integer(mws_) :: mwstemp
   character(len=15) :: mode
   character(len=200) :: errmsg

   !  Test input/output arguments

   if(nrhs_in.lt.1) call MATLAB_error("Insufficient input arguments")
   if(nlhs_in.gt.6) &
      call MATLAB_error("Too many output arguments")

   if(.not. MATLAB_is_character(prhs(1))) &
      call MATLAB_error("First argument must be string")
   mwstemp = len(mode)
   call MATLAB_get_string(prhs(1), mode, mwstemp)

   select case (trim(mode))
   case('control')
      ! control = hsl_mi20('control')
      call mi20_matlab_control(nlhs_in, plhs, nrhs_in-1_int4_, prhs(2))

   case('setup')
      ! [inform, handle] = hsl_mi20('setup',A,control)
      ! DEPRECATED: [inform,cdata1,cdata2,keep1,keep2,keep3] = ...
      !             hsl_mi20('setup',A,control)

      ! Setup handle and store in second output argument
      if(nlhs_in.lt.2) call MATLAB_error("Insufficient output arguments")
      handle = mi20_new_handle()
      plhs(2) = fortran_to_matlab(handle)

      ! Call work routine
      call mi20_matlab_setup(nlhs_in, plhs, nrhs_in-1_int4_, prhs(2), &
         handles(handle)%A, handles(handle)%cd_ptr%cd, handles(handle)%keep, &
         handles(handle)%control, errmsg)
      if(trim(errmsg).ne."No errors") then
         call cleanup_handle(handle)
         call MATLAB_error(trim(errmsg))
      endif

   case('precondition')
      ! [x, inform] = hsl_mi20('precondition', z, handle)
      ! DEPRECATED: [x, inform,keep1,keep2,keep3] =  ...
      !             hsl_mi20( 'precondition',z,cdata1,cdata2,keep1,keep2,keep3)

      ! Obtain and check handle
      if(nrhs_in.lt.3) call MATLAB_error("Insufficient input arguments")
      call matlab_to_fortran(prhs(3), handle, 'handle')
      if(handle.lt.1 .or. handle.gt.total_handles) &
         call MATLAB_error("Invalid handle")
      if(.not.associated(handles(handle)%A)) &
         call MATLAB_error("Invalid handle")

      ! Call worker routine
      call mi20_matlab_precondition(nlhs_in, plhs, nrhs_in-1_int4_, prhs(2),&
         handles(handle)%A, handles(handle)%cd_ptr%cd, handles(handle)%keep, &
         handles(handle)%control)

   case('destroy')
      ! hsl_mi20('destroy', handle)

      ! Obtain and check handle
      if(nrhs_in.lt.2) call MATLAB_error("Insufficient input arguments")
      call matlab_to_fortran(prhs(2), handle, 'handle')
      if(handle.lt.1 .or. handle.gt.total_handles) &
         call MATLAB_error("Invalid handle")

      ! Destroy anything that is there
      call cleanup_handle(handle)

   case default
      write(errmsg, "(3a)") "Unrecognised action: '", trim(mode), "'"
      call MATLAB_error(errmsg)
   end select
end subroutine mexFunction
