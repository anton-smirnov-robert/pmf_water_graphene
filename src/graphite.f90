module graphite
    
    use constants 

    !!!!!! Graphite related
    double precision , parameter :: d_cc = 1.42 / au_to_angstrom
    double precision , parameter :: a_lattice = sqrt(3._dp) * d_cc
    double precision , parameter :: c = 3.35 / au_to_angstrom

    double precision , parameter :: a1(2) = (/ 1._dp/2._dp , sqrt(3._dp)/2._dp  /) * a_lattice 
    double precision , parameter :: a2(2) = (/ -1._dp/2._dp , sqrt(3._dp)/2._dp  /) * a_lattice 

    double precision , parameter :: b1(2) = twopi/a_lattice * (/ 1._dp , 1./sqrt(3._dp)  /)
    double precision , parameter :: b2(2) = twopi/a_lattice * (/ 1._dp , -1/sqrt(3._dp)  /) 

    double precision , parameter :: kc = 3.18 ! p orbital constant

    double precision , parameter :: big_K(2) = twopi/a_lattice * (/ 1./3._dp  , 1/sqrt(3._dp) /) 
    double precision , parameter :: big_Kp(2) = twopi/a_lattice * (/ 2./3._dp , 0._dp /) 

    double precision , parameter :: tau(2)= 1./3. * (a1+a2)

    ! Tight-binding 
    double precision , parameter :: e0=0
    double precision , parameter :: gamma0 = 3.16  * 1e3 / au_to_meV
    double precision , parameter :: gamma1 = 0.39  * 1e3 / au_to_meV
    double precision , parameter :: gamma2 = -0.02  * 1e3 / au_to_meV
    double precision , parameter :: gamma3 = 0.315  * 1e3 / au_to_meV
    double precision , parameter :: gamma4 = 0.044 * 1e3 / au_to_meV
    double precision , parameter :: gamma5 = 0.038 * 1e3 / au_to_meV
    double precision , parameter :: big_delta = -0.008/2 *1e3 / au_to_meV
    double complex , parameter :: diag0(4) = (/ e0 + big_delta , e0 , e0 + big_delta , e0 /)
    double complex , parameter :: diag3(4) = (/ gamma5 , gamma2 ,gamma5 , gamma2 /)

    double precision , parameter :: vf = sqrt(3._dp ) * gamma0 * a_lattice / 2  
    !implicit none 

    contains 

    subroutine construct_H(H,k_vec)
        double complex , intent(inout) :: H(:,:) 
        double precision , intent(in) :: k_vec(2)
        double complex :: f,fp
        integer :: i,j,nb
        double complex :: diag1(4),diag2(4)
        
        ! we fill only th eupper diagonal
        H=cmplx(0.,0.)
        
        f=f_phase_func2(k_vec)
        fp=conjg(f)
        
        diag1 = (/ f*gamma0 , fp * gamma4 , fp*gamma0 , f*gamma4 /)
        diag2 = (/ gamma1 +cmplx(0,0), f*gamma3 , gamma1+cmplx(0,0) , fp*gamma3 /)
        
        nb = size(H(:,1))
        
        do i=1,nb,4
          do j=0,3 
        
            if ( i+j+3 <= nb ) then 
              h(i+j,i+j) = diag0(j+1)
              h(i+j,i+j+1) = diag1(j+1)
              h(i+j,i+j+2) = diag2(j+1)
              h(i+j,i+j+3) = diag3(j+1)
              cycle
            end if
            if ( i+j+2 <= nb ) then 
                h(i+j,i+j) = diag0(j+1)
                h(i+j,i+j+1) = diag1(j+1)
                h(i+j,i+j+2) = diag2(j+1)
              cycle
            end if 
            if ( i+j+1 <= nb ) then 
                h(i+j,i+j) = diag0(j+1)
                h(i+j,i+j+1) = diag1(j+1)
              cycle
            end if
            if ( i+j <= nb ) then 
                h(i+j,i+j) = diag0(j+1)
              cycle 
            end if 
            
          end do
        end do
        
        !print * , h(1,1) , h(1,2) ,h(2,1) , h(2,2)
        
    end subroutine construct_H


    function f_phase_func(k_vec) result(f_phase)
            double precision, intent(in) :: k_vec(2)
            double complex :: f_phase
            f_phase =  1_dp + exp( -i_complex * dot_product(k_vec,a1) ) + exp( - i_complex * dot_product(k_vec,a2) )
    end function 

    function f_phase_func2(k_vec) result(f_phase)
            double precision, intent(in) :: k_vec(2)
            double complex :: f_phase
            f_phase =  1_dp + 2*exp( i_complex *sqrt(3.) * k_vec(2) *a_lattice / 2 ) * cos(k_vec(1)*a_lattice/2)
    end function  
          
    !function b0_func(x) result(b0)
    !        double precision, intent(in) :: x
    !        double precision :: b0
    !        b0 = 1/(1+x**2)**3 /2 * 1./(1.+x**2)**(3./2.) * ( (sqrt(1+x**2) + x/4.)/ (sqrt(1+x**2) + x)**4  &
    !            !+  ( sqrt(1+x**2) - x/4)/(sqrt(1+x**2) - x)**4  ) !1/(1+x**2)**3 !1/(1+x**2)**3 !
    !end function b0_func

    function chi_e0_graphene(q,kf) result(chi_e0)
        double precision, intent(in) :: q,kf 
        double precision :: tmp,x
        double precision :: chi_e0

        x=q/2/kf
        tmp =  dsign(0.5_dp,x-1._dp) + .5_dp
    
        chi_e0 = - 2*kf/vf/pi * &
   &  ( 1._dp +  tmp * (  - 0.5*sqrt( (1-1/x**2)*tmp ) - x * 0.5 * asin(1/x*tmp) + pi*x/4 ) )

    end function chi_e0_graphene 
end module graphite 
