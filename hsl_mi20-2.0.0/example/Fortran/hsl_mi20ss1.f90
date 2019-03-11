program mi20_example_2
      use hsl_mi20_single

      implicit none
      integer, parameter :: wp = kind(1.0)
      integer, parameter :: n = 10 ! size of system to solve

! derived types
      type(mi20_data), dimension(:), allocatable :: coarse_data
      type(mi20_control) :: control
      type(mi20_solve_control) :: solve_control
      type(mi20_info) :: info
      type(mi20_keep) :: keep

! matrix values and arrays (pre conversion)
      integer :: ne
      integer, allocatable :: row(:), col(:)
      real(wp), allocatable :: val(:)

! solution arrays
      real(wp) :: rhs(n)
      real(wp) :: sol(n)

! others 
      integer :: i
      

! generate matrix A
      ne = 28
      allocate(row(ne),col(ne),val(ne))
      do i = 1, n
         ! diagonal ...
         row(i) = i; col(i) = i; val(i) = 2.0
         if (i < n) then
            ! superdiagonal...
            row(i + n) = i; col(i + n) = i+1; val(i+n) = -1.0
            ! subdiagonal...
            row(i + 2*n-1) = i+1; col(i + 2*n-1) = i; val(i+2*n-1) = -1.0
         end if
      end do

      rhs = 1.0

! call mi20_setup
      call mi20_setup_coord(row, col, val, &
                            ne, n,       & 
                            coarse_data, keep, control, info)
      if (info%flag < 0) then
        write(*,*) "Error return from mi20_setup_coord"
        stop
      end if

      ! deallocate coord matrix as no longer needed
      deallocate(row,col,val) 


! call solver     
      solve_control%krylov_solver = 1 ! solve using CG
      solve_control%rel_tol = 1e-4 ! set the relative convergence tolerance
      call mi20_solve(coarse_data, rhs, sol, keep, control, solve_control, info)
      if (info%flag < 0) then
        write(*,*) "Error return from mi20_solve"
        stop
      end if
      
! deallocation
      call mi20_finalize(coarse_data, keep, control, info)

end program mi20_example_2
