module declaration

    use constants 

    !!!!! Arrays 
    
    double precision , allocatable :: x0(:,:),eps_1(:,:),vc(:,:),veff(:,:),w_w(:,:) 
    double precision , allocatable :: n_ext(:),F_w(:),F_e0e(:),F_ewe(:),F_ewe_w(:) 
    double precision , allocatable :: to_inverse(:,:),w_ewe(:,:)
    double precision , allocatable :: eps_1_local(:),all_z(:),all_a(:),all_q(:),tmp(:)


    !!! Others 
    integer :: nz,iz,iq,nq,ia,iz0
    double precision  :: L,dz,q,q_sq,big_q,z,z0,tmp2
    double precision  :: da,kf,ef,b,d0,sigma0,prefactor,chi_e0,chi_e
    double precision :: start,finish

end module 
