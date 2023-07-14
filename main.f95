program matrix_multiplication
  implicit none
  
  ! Declare variables
  integer :: i, j, k
  integer :: m1, n1, m2, n2
  real :: alpha, beta
  real, dimension(:,:), allocatable :: A, B, C
  integer :: ierr
  
  ! Read in matrix data from file using read_matrix subroutine
  call read_matrix('input.txt', A, B, m1, n1, m2, n2, ierr)
  if (ierr /= 0) then
    write(*,*) 'Error reading matrix data from file!'
    stop
  end if
  
  ! Allocate memory for output matrix
  allocate(C(m1, n2), stat=ierr)
  if (ierr /= 0) then
    write(*,*) 'Error allocating memory for output matrix!'
    stop
  end if
  
  ! Multiply matrices
  alpha = 1.0
  beta = 0.0
  call matrix_multiply(alpha, A, B, beta, C, m1, n1, n2)
  
  ! Print out result
  write(*,*) 'Result:'
  do i = 1, m1
    write(*, '(6F8.2)') C(i,:)
  end do
  
  ! Deallocate memory
  deallocate(A, B, C)
  
contains

  subroutine read_matrix(filename, matrix1, matrix2, rows1, cols1, rows2, cols2, ierr)
    character(len=*), intent(in) :: filename
    integer, intent(out) :: rows1, cols1, rows2, cols2, ierr
    real, dimension(:,:), allocatable, intent(out) :: matrix1, matrix2
    
    integer :: i, j
    real :: value
    integer :: unit_num
    
    ! Open file with error checking
    open(newunit=unit_num, file=filename, status='old', action='read', iostat=ierr)
    if (ierr /= 0) then
      write(*,*) 'Error opening file ', filename
      return
    end if
    
    ! Read in matrix 1 dimensions
    read(unit_num, *, iostat=ierr) rows1, cols1
    if (ierr /= 0) then
      write(*,*) 'Error reading matrix 1 dimensions!'
      return
    end if
    
    ! Allocate memory for matrix 1
    allocate(matrix1(rows1, cols1), stat=ierr)
    if (ierr /= 0) then
      write(*,*) 'Error allocating memory for matrix 1!'
      return
    end if
    
    ! Read in matrix 1 data
    do i = 1, rows1
      read(unit_num, *, iostat=ierr) matrix1(i,:)
      if (ierr /= 0) then
        write(*,*) 'Error reading matrix 1 data!'
        return
      end if
    end do
    
    ! Read in matrix 2 dimensions
    read(unit_num, *, iostat=ierr) rows2, cols2
    if (ierr /= 0) then
      write(*,*) 'Error reading matrix 2 dimensions!'
      return
    end if
    
    ! Allocate memory for matrix 2
    allocate(matrix2(rows2, cols2), stat=ierr)
    if (ierr /= 0) then
      write(*,*) 'Error allocating memory for matrix 2!'
      return
    end if
    
    ! Read in matrix 2 data
    do i = 1, rows2
  read(unit_num, *, iostat=ierr) matrix2(i,:)
  if (ierr /= 0) then
    write(*,*) 'Error reading matrix 2 data!'
    return
  end if
end do

! Close file
close(unit=unit_num)

ierr = 0
end subroutine read_matrix

! Declare subroutine variables in subroutine
subroutine matrix_multiply(alpha, A, B, beta, C, m, n, p)
    real, intent(in) :: alpha, beta
    real, dimension(:,:), intent(in) :: A, B
    real, dimension(:,:), intent(out) :: C
    integer, intent(in) :: m, n, p
    integer :: i, j, k

! Perform matrix multiplication using three nested loops
do i = 1, m
  do j = 1, p
    ! Initialize each element of matrix C to 0
    C(i,j) = 0.0
    do k = 1, n
        ! Compute the dot product of row i of matrix A and column j of matrix B    
      C(i,j) = C(i,j) + alpha * A(i,k) * B(k,j)
    end do
    ! Add the scaling factor for matrix C
    C(i,j) = C(i,j) + beta * C(i,j)
  end do
end do

end subroutine matrix_multiply

end program matrix_multiplication
