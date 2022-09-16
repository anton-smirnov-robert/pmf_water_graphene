module constants

    integer, parameter :: dp = SELECTED_REAL_KIND(8)
    double complex , parameter :: zero_cmplx = cmplx( 0._dp  ,0._dp ,kind=dp )
    double precision , parameter :: pi = acos(-1.0_dp)
    double precision , parameter :: twopi = 2.0_dp*acos(-1.0_dp)
    double complex , parameter :: i_complex = cmplx( 0. , 1. ,kind=dp)

    !!!!!! CONSTANTS We use atomic units
    double precision,parameter :: au_to_second = 2.4188843265857e-17
    double precision,parameter :: au_to_fs =  au_to_second * 1e15
    double precision,parameter :: au_to_Thz = 1/au_to_second * 1e-12  / twopi
    double precision,parameter :: au_to_angstrom = 0.529177249
    double precision,parameter :: au_to_debye = 1/0.393456 
    double precision,parameter :: au_to_meV = 27.21 * 1e3
    double precision,parameter :: au_to_cm_1 = au_to_Thz  * 33.356 
    double precision,parameter :: au_to_kelvin = 3.1577464e5
    double precision,parameter :: au_to_kg = 9.1093837015e-31
    double precision,parameter :: au_to_m = 1/(1.8897259886e10)
    double precision,parameter :: au_to_kg_m_minus2_s_minus1  =  au_to_kg / au_to_m**2 / au_to_second
    double precision,parameter :: au_to_kg_s_minus1  =  au_to_kg  / au_to_second
    double precision,parameter :: au_to_N_s_m_minus3 = au_to_kg_m_minus2_s_minus1 
    double precision,parameter :: au_to_ms = 1/(1e10 / au_to_angstrom * (au_to_fs*1e-15))
    double precision,parameter :: kT = 0.944_dp * 1e-3, beta = 1/kT

    contains 

    function norm(k_vec) 
        double precision :: k_vec(2) ,norm 
        norm = sqrt( k_vec(1)**2 + k_vec(2)**2 )
    end function norm

    subroutine construct_vector(all_vector,v_min,v_max,ntot,dv)  
        double precision, allocatable , intent(out) :: all_vector(:)
        double precision, intent(out) :: dv 
        double precision :: v_min, v_max 
        integer :: i ,ntot
      
        dv = (v_max-v_min)/ntot
        allocate(all_vector(ntot),source=0.0_dp)
        do i=1,ntot
          all_vector(i) = (i-1)*dv + v_min 
        end do

      end subroutine 
      
end module constants 
