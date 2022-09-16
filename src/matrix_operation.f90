module matrix_operation
    
    use constants 
    
    implicit none 

    integer :: n, lwork
    integer ,  allocatable :: ipiv(:)
    double precision, allocatable :: work(:)





        !!!!!! Diagonalization Hamiltonian
    type eigensys_t
    double complex, allocatable :: eigen_vecs(:,:)
    double precision, allocatable :: eigen_vals(:)
    end type eigensys_t

    type(eigensys_t)  :: eigensystem_k
    type(eigensys_t) :: eigensystem_kq

    

    contains 

    subroutine initialize_inverse(n,ipiv,work,lwork)
    
      integer , intent(inout) :: n, lwork 
      integer , intent(inout), allocatable :: ipiv(:)
      double precision, allocatable , intent(inout) :: work(:)
    
      lwork = n*64
  
      allocate(work(lwork),source= 0._dp )
      allocate(ipiv(n),source=0)
    
    end subroutine initialize_inverse
    
    
    subroutine get_inverse(matrix, n, ipiv, work, lwork)
    
      double precision,     intent(inout)    :: matrix(:,:)
      integer , intent(inout) :: n, lwork
      integer :: info
      integer ,  intent(inout) :: ipiv(:)
      double precision, intent(inout)  :: work(:)
    
      call dgetrf(n,n, matrix, n , ipiv , info )        

      !if (info /= 0) then
      !  write (*,*) "Could not calculate optimal sizes for work arrays!"
      !  stop
      !end if

      call dgetri(n, matrix, n , ipiv, work, lwork , info)
    
    end subroutine get_inverse 


    subroutine inverse(a,c,n)
        !============================================================
        ! Inverse matrix
        ! Method: Based on Doolittle LU factorization for Ax=b
        ! Alex G. December 2009
        !-----------------------------------------------------------
        ! input ...
        ! a(n,n) - array of coefficients for matrix A
        ! n      - dimension
        ! output ...
        ! c(n,n) - inverse matrix of A
        ! comments ...
        ! the original matrix a(n,n) will be destroyed 
        ! during the calculation
        !===========================================================
    implicit none 
    integer :: n
    double precision :: a(n,n), c(n,n)
    double precision :: L(n,n), U(n,n), b(n), d(n), x(n)
    double precision :: coeff
    integer :: i, j, k
    
    ! step 0: initialization for matrices L and U and b
    ! Fortran 90/95 aloows such operations on matrices
    L=0.0
    U=0.0
    b=0.0
    
    ! step 1: forward elimination
    do k=1, n-1
       do i=k+1,n
          coeff=a(i,k)/a(k,k)
          L(i,k) = coeff
          do j=k+1,n
             a(i,j) = a(i,j)-coeff*a(k,j)
          end do
       end do
    end do
    
    ! Step 2: prepare L and U matrices 
    ! L matrix is a matrix of the elimination coefficient
    ! + the diagonal elements are 1.0
    do i=1,n
      L(i,i) = 1.0
    end do
    ! U matrix is the upper triangular part of A
    do j=1,n
      do i=1,j
        U(i,j) = a(i,j)
      end do
    end do
    
    ! Step 3: compute columns of the inverse matrix C
    do k=1,n
      b(k)=1.0
      d(1) = b(1)
    ! Step 3a: Solve Ld=b using the forward substitution
      do i=2,n
        d(i)=b(i)
        do j=1,i-1
          d(i) = d(i) - L(i,j)*d(j)
        end do
      end do
    ! Step 3b: Solve Ux=d using the back substitution
      x(n)=d(n)/U(n,n)
      do i = n-1,1,-1
        x(i) = d(i)
        do j=n,i+1,-1
          x(i)=x(i)-U(i,j)*x(j)
        end do
        x(i) = x(i)/u(i,i)
      end do
    ! Step 3c: fill the solutions x(n) into column k of C
      do i=1,n
        c(i,k) = x(i)
      end do
      b(k)=0.0
    end do
    end subroutine inverse
    




    subroutine initialize_diago(n,work, lwork, rwork, lrwork, iwork, liwork)
    
      integer , intent(inout) :: n, lwork, lrwork, liwork
      integer ,  allocatable , intent(inout) :: iwork(:)
      double precision, allocatable , intent(inout) :: rwork(:)
      double complex, allocatable , intent(inout) :: work(:)
    
      liwork = 5*n + 3
      lrwork = 2*n**2 + 5*n + 1
      lwork = n*(n+2)
    
      allocate(work(lwork),source=cmplx( 0._dp  ,0._dp ,kind=dp))
      allocate(iwork(liwork),source=0)
      allocate(rwork(lrwork),source=0._dp)
    
    end subroutine initialize_diago
    
    
    subroutine get_eigensystem(matrix, eigensystem , n, work, lwork, rwork, lrwork, iwork, liwork )
    
      double complex,     intent(inout)    :: matrix(:,:)
      type(eigensys_t), intent(inout) :: eigensystem
      integer , intent(inout) :: n, lwork, lrwork, liwork
      integer :: error
      integer ,  intent(inout) :: iwork(:)
      double precision, intent(inout)  :: rwork(:)
      double complex, intent(inout)  :: work(:)
    
    
    
      call zheevd('v', 'u', n,matrix , n,eigensystem%eigen_vals,          &
      &           work, lwork, rwork, lrwork, iwork, liwork, error)
      if (error /= 0) then
        write (*,*) "Could not calculate optimal sizes for work arrays!"
        stop
      end if
      eigensystem%eigen_vecs = matrix 
    
    
    end subroutine get_eigensystem    
    
  

    function eye(n) result(res)
        integer, intent(in) :: n
        integer :: i,j
        double precision :: res(n,n)
        do i=1,n
          do j=1,n
            if (i==j) then 
               res(i,i)=1._dp  
            else 
              res(i,j)=0._dp 
            end if
          end do
        end do
    end function
      
    function outerprod_r(a,b)
        double complex, DIMENSION(:), INTENT(IN) :: a,b
        double complex, DIMENSION(size(a),size(b)) :: outerprod_r 
    
        outerprod_r = spread(a, dim=2 ,ncopies=size(b)) * &
        spread(b,dim=1,ncopies=size(a)) 
    
    end function outerprod_r

    function outerprod_real(a,b)
        double precision, DIMENSION(:), INTENT(IN) :: a,b
        double precision, DIMENSION(size(a),size(b)) :: outerprod_real

        outerprod_real = spread(a, dim=2 ,ncopies=size(b)) * &
        spread(b,dim=1,ncopies=size(a))

    end function outerprod_real


end module
