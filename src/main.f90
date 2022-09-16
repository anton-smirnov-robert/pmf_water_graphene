program main

    use constants 
    use spce
    use input , only: getinput
    use declaration
    use graphite 
    use rpa 
    use matrix_operation

    implicit none

    !!! INPUT
    nq = getinput%int("nq",defaultvalue=100,assert=">0") 
    L = getinput%dp("L",defaultvalue=4._dp,assert=">0") * 10 / au_to_angstrom 
    dz = getinput%dp("dz",defaultvalue=0.05_dp,assert=">0") / au_to_angstrom  !converged
    ef = getinput%dp("ef",defaultvalue=25._dp,assert=">0") / au_to_meV  
    b = getinput%dp("b",defaultvalue=2._dp,assert=">0") / au_to_angstrom
    d0 = getinput%dp("d0",defaultvalue=1.3_dp,assert=">0") / au_to_angstrom
    sigma0 = getinput%dp("sigma0",defaultvalue=0.3_dp,assert=">0") / au_to_angstrom
    !!!!!!

    nz = int(L/dz)

    call construct_vector(all_z,dz,L,nz,dz)
    call initialize_inverse(nz,ipiv,work,lwork)

    allocate( x0(nz,nz), source=0._dp)
    allocate( veff(nz,nz), source=0._dp)
    allocate( to_inverse(nz,nz), source=0._dp)
    allocate( vc(nz,nz), source=0._dp )

    allocate( w_w(nz,nz), source=0._dp )
    allocate( w_ewe(nz,nz), source=0._dp )

    allocate( n_ext(nz), source=0._dp )
    allocate( F_w(nz), source=0._dp )
    allocate( F_e0e(nz), source=0._dp )
    allocate( F_ewe(nz), source=0._dp )
    allocate( F_ewe_w(nz), source=0._dp )


    !!!!!!! GET eps_1_local(z) !!!!!!!!!!
 
    if (  getinput%log("eps_1_local",defaultvalue=.false.) )  then

        call get_x(qlimit,x0,vc,veff,to_inverse,d0,sigma0)     ! WARNING : Now x0 is actually x

        allocate( eps_1 (nz,nz), source=0._dp )
        allocate( eps_1_local(nz), source=0._dp )

        eps_1 = eye(nz) / dz  + matmul( vc , x0 ) * dz
        
        call get_eps_1_local(eps_1,eps_1_local)

        deallocate(eps_1)

        do iz=1,nz
        print *, all_z(iz) , eps_1_local(iz)
        end do

        stop ! bye 

    end if  

    !!!!!!! GET PMF !!!!!!!!!!!!!!!!!!

    call construct_vector(all_a,-1._dp,8._dp,nq,da)

    kf = ef/vf 

    all_q = kf*exp(all_a)

    do iq=1,nq 
        
        q = all_q(iq)
        prefactor = 0.5*q**2/twopi * da

        call get_x(q,x0,vc,veff,to_inverse,d0,sigma0)

        w_w = vc + matmul( matmul( vc , x0 ) , vc ) * dz**2 

        chi_e0 = chi_e0_graphene(q,kf)

        !! water 
        do iz0=1,nz
            n_ext = i_func(q,b**2, ( all_z-all_z(iz0) )**2   )/2/b
            F_w(iz0) = F_w(iz0) + prefactor * sum( matmul(n_ext,w_w) *  n_ext ) * dz**2
        end do

        !!! uncoupled
        chi_e = chi_e0/(1-vc(1,1)*chi_e0)
        w_ewe = outerprod_real ( vc(:,1) , vc(1,:) ) + outerprod_real ( vc(:,nz) , vc(nz,:) )
        w_ewe = w_ewe + vc(1,nz)*chi_e * ( outerprod_real ( vc(:,1) , vc(nz,:) ) + outerprod_real ( vc(:,nz) , vc(1,:) ) )
        w_ewe = chi_e/( 1- (vc(1,nz)*chi_e)**2 ) * w_ewe
         

        do iz0=1,nz
            n_ext = i_func(q,b**2, ( all_z-all_z(iz0) )**2   )/2/b
            F_e0e(iz0) = F_e0e(iz0) + prefactor * sum( matmul( n_ext , w_ewe  ) * n_ext ) * dz**2
        end do

        !!! semi-coupled
        chi_e = chi_e0/(1-vc(1,1)*chi_e0)
        w_ewe = outerprod_real ( w_w(:,1) , w_w(1,:) ) + outerprod_real ( w_w(:,nz) , w_w(nz,:) )
        w_ewe = w_ewe + w_w(1,nz)*chi_e * (  outerprod_real ( w_w(:,1) , w_w(nz,:) ) + outerprod_real ( w_w(:,nz) , w_w(1,:) ) )   
        w_ewe = chi_e/( 1- (w_w(1,nz)*chi_e)**2 ) * w_ewe


        do iz0=1,nz
            n_ext = i_func(q,b**2, ( all_z-all_z(iz0) )**2   )/2/b
            F_ewe(iz0) = F_ewe(iz0) + prefactor * sum( matmul( n_ext , w_ewe ) * n_ext ) * dz**2
        end do

        !!! fully-coupled
        chi_e = chi_e0/(1-w_w(1,1)*chi_e0)
        w_ewe = outerprod_real ( w_w(:,1) , w_w(1,:) ) + outerprod_real ( w_w(:,nz) , w_w(nz,:) )
        w_ewe = w_ewe + w_w(1,nz)*chi_e * (  outerprod_real ( w_w(:,1) , w_w(nz,:) ) + outerprod_real ( w_w(:,nz) , w_w(1,:) ) ) 
        w_ewe = chi_e/( 1- (w_w(1,nz)*chi_e)**2 ) * w_ewe


        do iz0=1,nz
            n_ext = i_func(q,b**2, ( all_z-all_z(iz0) )**2   )/2/b
            F_ewe_w(iz0) = F_ewe_w(iz0) + prefactor * dot_product( matmul( n_ext , w_ewe ) , n_ext ) * dz**2
        end do

   end do

   do iz0=1,nz
    print *, all_z(iz0) , F_w(iz0) , F_e0e(iz0) , F_ewe(iz0) , F_ewe_w(iz0)
   end do 
   
end program 




              

