program mi20_example
      use hsl_mi20_double
      use hsl_zd11_double
      use hsl_mc65_double
   
      implicit none
      integer, parameter :: wp = kind(1.0d0)
      integer, parameter :: m = 10 ! size of system to solve

! derived types
      type(zd11_type) :: a
      type(mi20_data), dimension(:), allocatable :: coarse_data
      type(mi20_control) :: control
      type(mi20_info) :: info
      type(mi20_keep) :: keep
      type(ma48_control) :: ma48_cntl
   
! Arrays and scalars required by the CG code mi21
      real(kind=wp) :: cntl(5),rsave(6)
      integer :: icntl(8),isave(10),info21(4)
      real(kind=wp) :: w(m,4)
      real(kind=wp) :: resid
      integer :: locy, locz, iact

      external mi21id, mi21ad

      integer :: info65  ! mc65 error flag

! generate matrix A
      call matrix_gen(a, m)
   
! Prepare to use the CG code mi21 with preconditioning
      call mi21id(icntl, cntl, isave, rsave)
      icntl(3) = 1

! set right hand side to vector of ones
      w(:,1) = 1
   
! call mi20_setup
      call mi20_setup(a, coarse_data, keep, control, info)
      if (info%flag < 0) then
        write(*,*) "Error return from mi20_setup"
        stop
      end if
      
! solver loop
      iact = 0
      do
        call mi21ad(iact, m, w, m, locy, locz, resid, icntl, cntl, info21, &
             isave, rsave)
        if (iact == -1) then
          write(*,*) "Error in solver loop"
          exit

        else if (iact == 1) then
          write(*,'(a,i3,a)') " Convergence in ", info21(2), " iterations"
          write(*,'(a,es12.4)') " 2-norm of residual =", resid
          exit

        else if (iact == 2) then
          call mc65_matrix_multiply_vector(a, w(:,locz), w(:,locy), info65)

        else if (iact == 3) then
          call mi20_precondition(a, coarse_data, w(:,locz), w(:,locy), keep, &
               control, info, ma48_cntl)
          if (info%flag < 0) then
            write(*,*) "Error return from mi20_precondition"
            exit
          end if

        end if
      end do
   
! deallocation
      call mi20_finalize(coarse_data, keep, control, info)
      deallocate(a%col,a%val,a%ptr)
   
contains
      subroutine matrix_gen(a, m)
      integer, intent(in) :: m ! size of matrix
      type(zd11_type), intent(out) :: a
      integer :: i,nnz,p
      
      nnz = m + 2*(m-1)
      allocate(a%col(nnz),a%val(nnz),a%ptr(m+1))
      a%m = m
      p = 1    ! pointer to next empty position
      do i = 1,m
        a%ptr(i) = p
        if (i==1) then   ! first row
          a%col(p) = i;    a%col(p+1) = i+1
          a%val(p) = 2.0;  a%val(p+1) = -1.0
          p = p+2
        else if (i==m) then ! last row
          a%col(p) = i-1;  a%col(p+1) = i
          a%val(p) = -1.0; a%val(p+1) = 2.0
          p = p+2
        else
          a%col(p) = i-1;  a%col(p+1) = i; a%col(p+2) = i+1
          a%val(p) = -1.0; a%val(p+1) = 2.0; a%val(p+2) = -1.0
          p = p+3
        end if
      end do
      a%ptr(m+1) = nnz+1
      end subroutine matrix_gen

end program mi20_example
