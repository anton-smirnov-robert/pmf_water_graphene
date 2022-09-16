module spce
    
    use constants 
 
    implicit none
    !!!!!! Graphite related
    double precision , parameter :: dOH = 1.0001_dp/au_to_angstrom ! grid spacing and heaviside ! 
    double precision , parameter :: dHH = 0.816495*2/au_to_angstrom
    double precision , parameter :: dOH_sq = dOH**2
    double precision , parameter :: dHH_sq = dHH**2
    double precision , parameter :: zH = 0.4238
    double precision , parameter :: epsw=71
    double precision , parameter :: n0 = 0.03297636 * au_to_angstrom**3

    double precision , parameter :: chi0_unit = - beta*n0*2*zH**2

    !!! 
    double precision , parameter :: kappa = 1.65 * au_to_angstrom
    double precision , parameter :: kappa_sq = kappa**2
    double precision , parameter :: gamma_pref =  0.99
    double precision , parameter :: qlimit = 1e-4*au_to_angstrom

    double precision , parameter :: mu_sq = 4*zH**2*dOH**2*(1 - (dHH/2/dOH)**2) 
    !double precision , parameter :: 2.351 /au_to_debye
    double precision , parameter :: alpha = beta*mu_sq/3
    double precision , parameter :: prefactor_v_eff = ( epsw/(epsw-1) - 1/(alpha*n0)/4._dp/pi )

    !!!!!!


    
    contains

        function i_func(q,d_sq,delta_z_sq) result(i_val)
            double precision, intent(in) :: q,d_sq,delta_z_sq(:)
            double precision :: i_val(size(delta_z_sq)),theta(size(delta_z_sq))
            theta = ( dsign(0.5_dp,d_sq-delta_z_sq) + .5_dp )
            i_val =  theta * bessel_j0(q*sqrt( (d_sq-delta_z_sq) *theta ) ) 
        end function


        function delta_chi0_func(q,delta_z_sq) result(delta_chi0)
            double precision, intent(in) :: q,delta_z_sq(:)
            double precision :: delta_chi0(size(delta_z_sq))
            delta_chi0 = chi0_unit * 0.5_dp * ( - 4*i_func(q,dOH_sq,delta_z_sq)/dOH + i_func(q,dHH_sq,delta_z_sq)/dHH ) 
        end function


                
        function v_coul(q,delta_z)
            double precision :: q,delta_z(:),v_coul(size(delta_z))
            v_coul = twopi/q * exp(-q*delta_z) 
        end function v_coul

    
        function v_coul_scalar(q,delta_z)
            double precision :: q,delta_z,v_coul_scalar
            v_coul_scalar = twopi/q * exp(-q*delta_z) 
        end function v_coul_scalar

        function v_eff(q,q_sq,big_q,delta_z)
            double precision :: q,q_sq,big_q,delta_z(:),v_eff(size(delta_z))
    
            v_eff = prefactor_v_eff * ( v_coul(q,delta_z) - v_coul(big_q,delta_z)  & 
                        & - gamma_pref * 2._dp/kappa * exp(- (q_sq/kappa_sq + delta_z**2*kappa_sq) /2 ) ) 
    
        end function 


        function v_eff_scalar(q,q_sq,big_q,delta_z)
            double precision :: q,q_sq,big_q,delta_z,v_eff_scalar
    
            v_eff_scalar = prefactor_v_eff * ( v_coul_scalar(q,delta_z) - v_coul_scalar(big_q,delta_z)  & 
                        & - gamma_pref * 2._dp/kappa * exp(- (q_sq/kappa_sq + delta_z**2*kappa_sq) /2 ) ) 
    
        end function 

        function g_z(z,L,d0,sigma0) result(profil)
            double precision, intent(in) :: z(:),L,d0,sigma0
            double precision :: profil(size(z))
            profil = .25_dp  *( tanh( (z-d0)/sigma0 ) + 1) * ( tanh( (L-z-d0)/sigma0 ) + 1)    
            !profil = ( dsign(.5_dp,z-d0) + 0.5_dp  ) * ( dsign(.5_dp,L-z-d0) + 0.5_dp  )  
      
        end function 
    
        function g_z_scalar(z,L,d0,sigma0) result(profil)
            double precision, intent(in) :: z,L,d0,sigma0
            double precision :: profil
            profil = .25_dp  *( tanh( (z-d0)/sigma0 ) + 1) * ( tanh( (L-z-d0)/sigma0 ) + 1)    
            !profil = ( dsign(.5_dp,z-d0) + 0.5_dp  ) * ( dsign(.5_dp,L-z-d0) + 0.5_dp  )     
        end function g_z_scalar

        
end module spce 
