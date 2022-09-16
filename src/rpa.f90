module RPA

    use constants
    use declaration
    use spce 
    use matrix_operation

    implicit none

    contains

    subroutine get_x(q,x0,vc,veff,to_inverse,d0,sigma0)

        double precision, allocatable , intent(inout):: x0(:,:) , vc(:,:) , veff(:,:),to_inverse(:,:)
        double precision :: q,d0,sigma0

        q_sq = q**2
        big_q = sqrt(q_sq+kappa_sq)

        x0=0._dp 
        veff=0._dp
        vc=0._dp 

        !!!!!!!!!! FILL MATRICES 
        do iz=1,nz-1
            
            z=all_z(iz)

            x0(iz,iz+1:) = delta_chi0_func(q, (all_z(iz+1:) - z )**2 ) & 
                        &   * g_z_scalar(z,L,d0,sigma0) * g_z(all_z(iz+1:),L ,d0,sigma0)
            vc(iz,iz+1:) = v_coul(q,all_z(iz+1:)-z )
            veff(iz,iz+1:) = v_eff(q,q_sq,big_q, all_z(iz+1:) - z   )
            
        end do 

        vc = vc + transpose(vc) + twopi/q * eye(nz)
        veff = veff + transpose(veff) + v_eff_scalar(q,q_sq,big_q,0._dp) * eye(nz) 
        x0 = x0 + transpose(x0) ! nothing on the diagonal here

        do iz=1,nz
            x0(iz,iz) = - sum( x0(iz,:) )     ! Electroneutrality 
        end do

        !!!!!!!!!!  MATRIX INVERSION (MOST time consuming)

        to_inverse = eye(nz) - matmul(  x0 , veff ) * dz**2
        
        !call cpu_time(start)

        call get_inverse(to_inverse, nz, ipiv, work, lwork)! WARNING : to_inverse is now inversed from here on!

        !call cpu_time(finish)
        !print * , finish-start

        x0 = matmul( to_inverse , x0 ) ! WARNING : now it x0 is x

    end subroutine 


    subroutine get_eps_1_local(eps_1,eps_1_local)
        double precision, allocatable , intent(inout):: eps_1(:,:),eps_1_local(:) 
        double precision  :: tmp(size(eps_1_local) ) 


                !!!!!!!!!! GET EPS_LOCAL  with gradient 
        tmp(1) = sum( eps_1(1,:)*all_z ) * dz
        tmp(2) = sum( eps_1(2,:)*all_z ) * dz
        eps_1_local(1) = ( tmp(2)-tmp(1) ) / dz 
    
        do iz=2,nz-1       
            tmp(iz+1) = sum( eps_1(iz+1,:)*all_z ) * dz
            eps_1_local(iz) = .5_dp *  ( tmp(iz+1)-tmp(iz-1) )/ dz 
        end do

        eps_1_local(nz) = ( tmp(nz)-tmp(nz-1) ) / dz

    end subroutine

end module 
