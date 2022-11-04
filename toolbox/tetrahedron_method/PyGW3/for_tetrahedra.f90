
subroutine intw(kwt, Ebnd, ef, atet, wtet, bloch_correction, nk, nbnd, ntet)
  ! This subroutine calculates the integration weight of each k-point for each band
  implicit none
  real(8), intent(out):: kwt(nk,nbnd)   ! the weight of each k-point for each band
  real*8,  intent(in) :: Ebnd(nk,nbnd)  ! Band energies
  real*8,  intent(in) :: ef             ! fermi energy
  integer, intent(in) :: atet(4,ntet)   ! tetrahedra index
  real*8,  intent(in) :: wtet(ntet)     ! weight of each tetrahedra
  integer, intent(in) :: nk, nbnd, ntet ! Number of irreducible k-points, bands, and tetrahedra
  logical, intent(in) :: bloch_correction
  !
  integer:: itet, i, ib, indx(4), kin
  real*8 :: ee(4), et(4), w1t(4), wcor(4)
  interface
     subroutine binary_insertion_sort_indx(indx, array,n)
       integer, intent(out):: indx(n)
       real*8,  intent(in) :: array(n)
       integer, intent(in) :: n
     end subroutine binary_insertion_sort_indx
  end interface
  kwt=0.d0 
  do itet=1,ntet
     do ib=1,nbnd 
        do i=1,4
           ee(i) = Ebnd(atet(i,itet)+1,ib)  ! atet(i,itet)+1 for python to fortran conversion
        enddo
        call binary_insertion_sort_indx(indx, ee,4)
        do i=1,4
           et(i) = ee(indx(i))
        enddo
        call intweight1t(w1t, et, ef)
        if( bloch_correction ) then 
           call bloechlcor(wcor, et, ef)
        else
           wcor = 0.d0
        endif
        do i=1,4
           kin = atet(indx(i),itet)+1
           kwt(kin,ib) = kwt(kin,ib) + ( w1t(i) + wcor(i) ) * wtet(itet)
        enddo
     enddo
  enddo
  return
end subroutine intw


subroutine intwsurf(kwt, Ebnd, ef, atet, wtet, nk, nbnd, ntet)
  ! This subroutine calculates the surface integration weight of each k-point for each band
  implicit none     
  real*8, intent(out) :: kwt(nk,nbnd)   ! the weight of each k-point for each band
  real*8,  intent(in) :: Ebnd(nk,nbnd)  ! Band energies
  real*8,  intent(in) :: ef             ! fermi energy
  integer, intent(in) :: atet(4,ntet)   ! tetrahedra index
  real*8,  intent(in) :: wtet(ntet)     ! weight of each tetrahedra
  integer, intent(in) :: nk, nbnd, ntet ! Number of irreducible k-points, bands, and tetrahedra
  !
  integer :: itet, i,ib, indx(4), kin
  real*8  :: ee(4), et(4), w1t(4)
  external  intweight1t
  
  kwt = 0.d0 
  do itet=1,ntet
     do ib=1,nbnd 
        do i=1,4
           ee(i) = Ebnd(atet(i,itet)+1,ib)  ! atet(i,itet)+1 for python to fortran conversion
        enddo
        call binary_insertion_sort_indx(indx, ee,4)
        do i=1,4
           et(i) = ee(indx(i))
        enddo
        w1t(:) = 0.0d0
        if (( et(1) < ef) .and. (ef < et(4) )) then
           call ksurf(w1t, et, ef)
           do i=1,4
              kin = atet(indx(i),itet)+1
              kwt(kin,ib) = kwt(kin,ib) + 6.0d0*w1t(i)*wtet(itet)
           enddo
        endif
     enddo
  enddo
  return
end subroutine intwsurf

real*8 function idos(energy, Ebnd, atet, wtet, nk, nbnd, ntet)
  ! This function calculates the integrated density of states at an energy energy using tetrahedron method
  implicit none     
  real*8,     intent(in) :: Ebnd(nk,nbnd) ! Band energies
  integer,    intent(in) :: atet(4,ntet) ! tetrahedra index
  real*8,     intent(in) :: wtet(ntet)   ! weight of each tetrahedra
  integer,    intent(in) :: nk, nbnd, ntet ! Number of irreducible k-points, bands, and tetrahedra
  real*8,     intent(in) :: energy     !  energy
  !
  integer:: itet,i,ib
  real*8 :: ee(4), val
  real*8, external :: intdos1t
  interface
     subroutine binary_insertion_sort(array,n)
       real*8, intent(inout) :: array(n)
       integer, intent(in)   :: n
     end subroutine binary_insertion_sort
  end interface
  idos=0.0d0
  do ib=1,nbnd
     do itet=1,ntet
        do i=1,4
           ee(i) = Ebnd(atet(i,itet)+1,ib)  ! atet(i,itet)+1 for python to fortran conversion
        enddo
        call binary_insertion_sort(ee,4)
        val = intdos1t(ee, energy)
        idos = idos + val*wtet(itet)
        !write(6,'(A,I2,1x,A,I2,1x,A,F12.8,1x,A,4F10.6,1x,A,4I2)') 'ib=',ib,'itet=', itet, 'val=', val, &
        !     & 'ee=', ee, 'atet=', atet(:,itet)
     enddo
  enddo
  return
end function idos

real*8 function dostet(energy, Ebnd, atet, wtet, nk, nbnd, ntet)
  ! This subroutine calculates the density of states at an energy energy
  implicit none
  real*8,     intent(in) :: Ebnd(nk,nbnd) ! Band energies
  integer,    intent(in) :: atet(4,ntet) ! tetrahedra index
  real*8,     intent(in) :: wtet(ntet)   ! weight of each tetrahedra
  integer,    intent(in) :: nk, nbnd, ntet ! Number of irreducible k-points, bands, and tetrahedra
  real*8,     intent(in) :: energy     !  energy
  !
  integer(4) :: itet,i,ib
  real(8)    :: ee(4), val
  real(8), external :: dos1t
  dostet=0.0d0
  do itet=1,ntet
     do ib=1,nbnd 
        do i=1,4
           ee(i) = Ebnd(atet(i,itet)+1,ib)
        enddo
        call binary_insertion_sort(ee,4)
        val = dos1t(ee, energy)
        dostet = dostet + val*wtet(itet)
     enddo
  enddo
end function dostet


subroutine bz_calcqdepw(kcw, iop_bzintq, sgnfrq, enk, ef, kqid, ncg, nomx, numin, omega, tetc, omgmax_ht, nomg_ht, beta, nkp, nomeg, nbnd, ntet)
  !
  ! This subroutine calculates the weights for the convolution of k and k-q points. Specifially it calculates the weight of
  !      P_{q,k,om_n} = 2*(f(e_k)- f(-e_{k-q})) / ( iom_n - e_{k-q} + e_k )
  !
  ! if sgnfrq == 1:   normal q-dependent bulk integration. 
  !     Integrate[ fw_i(x,y,z) , {z,0,1},{y,0,1-z},{x,0,1-y-z}]
  !          where fw_i(x,y,z) = [(1-x-y-z), x, y, z] corresponds to the four corners of tetrahedra.
  !
  ! if sgnfrq == 2:  weights for the Polarization with real frequencies
  !     Real{ Integrate[ fw_i(x,y,z) * 1/( omg - E(x,y,z) + 0.01*1j), {z,0,1},{y,0,1-z},{x,0,1-y-z}] }, where fw_i(x,y,z) = [(1-x-y-z), x, y, z] corresponds to the four corners of tetrahedra.
  !
  ! if sgnfrq == 3: weights for the Polarization with imaginary frequencies.
  !    2Real{ Integrate[ fw_i(x,y,z) * 1/( 1j*omg - E(x,y,z) ),       {z,0,1},{y,0,1-z},{x,0,1-y-z}] }, where fw_i(x,y,z) = [(1-x-y-z), x, y, z] corresponds to the four corners of tetrahedra.
  !    if iop_bzintq==-1:  numeric integration, which needs x_i,w_i
  !    else:              analytic integration
  !
  ! if sgnfrq == 4:   the q-dependent surface integratio (the surface is defined by e_jb-e_ib=omeg.
  !    Adds surface integration for metallic case

  !    if iop_bzintq== 0:  analytic integration of tetrahedra
  !    if iop_bzintq==-1:  numeric integration of tetrahedra, which needs x_i,w_i
  !    if iop_bzintq== 1:  simple summation with finite-temperature smearing and life-time broadening
  !    if iop_bzintq==-2:  calculate q-BZ weights using the Hilbert transform (HT) approach
  !
  implicit none
  integer,intent(in) :: iop_bzintq, sgnfrq, ncg, nomx, numin, nkp, nbnd, nomeg, nomg_ht, ntet
  integer,intent(in) :: kqid(nkp)
  real(8),intent(in) :: enk(nbnd,nkp), ef
  real(8),intent(in) :: omega(nomeg)
  integer,intent(in) :: tetc(4,ntet)
  real(8),intent(in) :: omgmax_ht, beta
  complex*16, intent(out) :: kcw(nbnd,nbnd-ncg, nomeg,nkp)
  !
  integer, parameter :: fout = 6
  real(8), parameter :: pi = 3.14159265358979d0
  real(8), parameter :: eta_freq = 0.01 ! this is the value of the small imaginary part that is needed for real frequency calculations 
  
  integer :: ib, jb, iom, io, ik, jk, ncbm, nvbm, n_i, it
  integer, allocatable :: tetln(:)         ! index of the q-linked tetrahedron 
  real(8), allocatable :: cwpar(:,:,:,:)   ! weight between the two bands
  real(8), allocatable :: wtet(:), w2_ht(:,:,:,:), omg_ht(:), x_i(:), w_i(:)
  real(8) :: edif, omg2, e_i, e_j, fnm, omg
  complex*16:: comg, kw, imag, freq_factor, wgh_om
  real(8) :: vol_small_tetra, dom_ht
  !
  imag = cmplx(0.d0,1.d0,8)
  allocate( tetln(ntet), wtet(ntet) )
  do ik=1,nkp
     do it=1,6
        ! if index of k=ik, then index of k-q is kqid[ik]
        tetln(6*(ik-1)+it) = 6*kqid(ik)+it  ! if we know the index of k-q, than we know which tetrahedron is linked with which
        ! they are enumerated in the same order, hence it is just 6*index of k-q point + which tetrahedron in the order
     enddo
  enddo
  wtet(:) = 1.d0/ntet
  
  allocate(cwpar(nbnd,nbnd,nkp,2))
  ncbm = numin+1 ! numin,nomx is in c-like type, hence fortran needs +1
  nvbm = nomx+1  !

  if (iop_bzintq.eq.-1) then
     n_i = 8
  else
     n_i = 1
  endif
  allocate( x_i(n_i), w_i(n_i) )

  select case(iop_bzintq)
  case(0,-1)
     !     Calculate the weights.
     if (iop_bzintq.eq.-1) call gauleg(0.d0, 1.d0, x_i, w_i, n_i)
     do iom=1,nomeg
        call convw(ef, omega(iom), sgnfrq, cwpar(:,:,:,1), iop_bzintq, enk, wtet, tetc, tetln, vol_small_tetra, nbnd, nkp, ntet, x_i, w_i, n_i)
        kcw(1:(ncg+nvbm),ncbm:(nbnd-ncg),iom,1:nkp) = cwpar(1:(ncg+nvbm),(ncg+ncbm):nbnd,1:nkp,1)
        !
        if(sgnfrq.eq.2) then  !! calculate the imaginary part for real freq 
           !call tetcw(nkp,ntet,nmaxb,wtet,enk,tnodes,tetln,kqid,tvol,efermi,omega(iom),4,cwpar(:,:,:,2))
           call convw(ef, omega(iom),    4, cwpar(:,:,:,2), iop_bzintq, enk, wtet, tetc, tetln, vol_small_tetra, nbnd, nkp, ntet, x_i, w_i, n_i)
           kcw(1:(ncg+nvbm),ncbm:(nbnd-ncg),iom,1:nkp) = kcw(1:(ncg+nvbm),ncbm:(nbnd-ncg),1:nkp,iom) + cwpar(1:(ncg+nvbm),(ncg+ncbm):nbnd,1:nkp,2)*imag
        endif
     enddo
  case(-2)
     ! calculate q-BZ weights using the Hilbert transform (HT) approach
     allocate( w2_ht(1:(ncg+nvbm), ncbm:(nbnd-ncg), nkp, nomg_ht), omg_ht(nomg_ht) )
     dom_ht = omgmax_ht/nomg_ht
     do iom=1,nomg_ht
        omg_ht(iom) = iom * dom_ht
        call convw(ef, omg_ht(iom), 4, cwpar(:,:,:,1), iop_bzintq, enk, wtet, tetc, tetln, vol_small_tetra, nbnd, nkp, ntet, x_i, w_i, n_i)
        w2_ht(1:(ncg+nvbm),ncbm:(nbnd-ncg),1:nkp,iom) = cwpar(1:(ncg+nvbm),(ncg+ncbm):nbnd,1:nkp,1)
     enddo
     
     do iom=1,nomeg
        omg2 = omega(iom)**2
        kcw(:,:,iom,:) = 0.d0  
        do io=1, nomg_ht ! using weights on more precise mesh to sum over many more frequencies
           if (sgnfrq.eq.3) then  ! Polarization on imaginary axis
              omg = omg_ht(io)
              wgh_om = dom_ht*(2.d0/pi) * omg/(omg**2+omg2)
           else                   ! Polarization on real axis
              comg = cmplx(omg_ht(io), -eta_freq, 8)
              wgh_om = dom_ht*(2.d0/pi) * comg/(comg**2-omg2)
           endif
           kcw(:,:,iom,:) = kcw(:,:,iom,:) + w2_ht(:,:,:,io) * wgh_om
        enddo
     enddo
     deallocate( w2_ht, omg_ht )
  case(1)
     !nbmaxpol == nbmax - ncg
     ! simple summation with finite-temperature smearing and life-time broadening 
     do iom=1,nomeg 
        omg = omega(iom) 
        do ik=1,nkp
           do ib = 1,nvbm+ncg
              e_i = enk(ib,ik)-ef
              jk = kqid(ik)+1  ! index of k-q
              do jb = ncbm,nbnd-ncg
                 e_j = enk(jb+ncg,jk)-ef
                 edif = e_j - e_i  ! e(k-q)-e(k)
                 if( edif < 0.d0 ) then 
                    kw = 0.d0 
                 else
                    fnm = fermi(e_i,beta)*(1.0-fermi(e_j,beta))
                    if (sgnfrq .eq. 3) then  ! Polarization on imaginary axis
                       freq_factor = -2.d0*edif/(omg**2+edif**2)             ! This is equivalent to 2*Re{ 1/(i*om-edif) } 
                    else                    ! Polarization on real axis
                       freq_factor =  1.d0/cmplx(omg - edif, eta_freq, 8) + 1.d0/cmplx(-omg - edif, eta_freq, 8)
                       !! strange I would expect
                       !freq_factor =  1.d0/cmplx(omg - edif, eta_freq) + 1.d0/cmplx(-omg - edif, -eta_freq)
                    endif
                    kw = fnm * freq_factor / nkp
                 endif
                 kcw(ib,jb,iom,ik) = kw
              enddo
           enddo
        enddo
     enddo
  end select


  !do ik=1,nkp
  !   do ib=1,ncg+nvbm
  !      do jb=ncbm,(nbnd-ncg)
  !         write(6,'(I4,1x,I4,1x,I4,1x,10F10.6)') ik,ib,jb, (dble(kcw(ib,jb,ik,iom)), iom=1,10)
  !      end do
  !   end do
  !end do
  
  deallocate( x_i, w_i )
  deallocate(tetln, wtet, cwpar)
  
contains
  
  real(8) function fermi(e, beta)
    implicit none
    real(8), intent(in) :: e, beta
    !
    real(8) :: x, wt
    x = e*beta
    if ( x > 100.d0) then
       wt = 0.d0
    elseif ( x < -100.d0) then
       wt = 1.d0
    else
       wt = 1.d0/(1.d0 + exp(x)) ! fermi(e/esmear)
    endif
    fermi = wt
  end function fermi
end subroutine bz_calcqdepw

subroutine bz_calcqdepw_par(kcw, enk, ef, kqid, ncg, nomx, numin, omega, tetc, nkp, nom, nbnd, ntet)
  !
  ! This subroutine calculates the weights for the convolution of k and k-q points. Specifially it calculates the weight of
  !      P_{q,k,om_n} = 2*(f(e_k)- f(-e_{k-q})) / ( iom_n - e_{k-q} + e_k )
  !
  ! if sgnfrq == 3: weights for the Polarization with imaginary frequencies.
  !    2Real{ Integrate[ fw_i(x,y,z) * 1/( 1j*omg - E(x,y,z) ),       {z,0,1},{y,0,1-z},{x,0,1-y-z}] }, where fw_i(x,y,z) = [(1-x-y-z), x, y, z] corresponds to the four corners of tetrahedra.
  !    if iop_bzintq==-1:  numeric integration, which needs x_i,w_i
  !    else:              analytic integration
  !
  ! if sgnfrq == 4:   the q-dependent surface integratio (the surface is defined by e_jb-e_ib=omeg.
  !    Adds surface integration for metallic case
  !
  !    if iop_bzintq== 0:  analytic integration of tetrahedra
  !
  implicit none
  integer,intent(in) :: ncg, nomx, numin, nkp, nbnd, nom, ntet
  integer,intent(in) :: kqid(nkp)
  real(8),intent(in) :: enk(nbnd,nkp), ef
  real(8),intent(in) :: omega(nom)
  integer,intent(in) :: tetc(4,ntet)
  real(8),intent(out):: kcw(nbnd,nbnd-ncg, nom,nkp)
  !
  real(8), parameter :: pi = 3.14159265358979d0
  !
  external convw1t_par
  !
  integer :: ik, ncbm, nvbm, it
  integer, allocatable :: tetln(:)         ! index of the q-linked tetrahedron 
  real(8), allocatable :: wtet(:)
  integer :: itet, ib, jb, i, kin, sgnfrq
  real(8) :: ee1(4), ee2(4), w1t(4,nom)
  !logical :: debug
  !
  sgnfrq = 3
  allocate( tetln(ntet), wtet(ntet) )
  do ik=1,nkp
     do it=1,6
        ! if index of k=ik, then index of k-q is kqid[ik]
        tetln(6*(ik-1)+it) = 6*kqid(ik)+it  ! if we know the index of k-q, than we know which tetrahedron is linked with which
        ! they are enumerated in the same order, hence it is just 6*index of k-q point + which tetrahedron in the order
     enddo
  enddo
  wtet(:) = 1.d0/ntet
  ncbm = numin+1 ! numin,nomx is in c-like type, hence fortran needs +1
  nvbm = nomx+1  !
  kcw(:,:,:,:)  = 0.0d0
  ! Weighst for he q-dependent bulk integration for Polarization
  ! ee1 is the (partially) occupied band and ee2 is the (partially) empty band
  do itet=1,ntet
     !$OMP PARALLEL DO PRIVATE(jb,i,ee2,ib,ee1,w1t,kin)&
     !$OMP& SHARED(kcw)&
     !$OMP& SCHEDULE(STATIC)
     do jb=ncg+ncbm,nbnd  ! unoccupied band 
        do i=1,4
           ee2(i) = enk( jb, tetc(i,tetln(itet))+1 )
        enddo
        if(maxval(ee2) <= ef) cycle ! ee2 is completely occupied =>0, because ee2 should be empty
        do ib=1,ncg+nvbm  ! occupied band
           do i=1,4
              ee1(i) = enk( ib, tetc(i,itet)+1 )
           enddo
           if( minval(ee1) > ef) cycle  ! ee1 is completely empty => value is zero because ee2 is empty as well
           w1t(1:4,1:nom) = 0.0d0
           call convw1t_par(w1t, ee1, ee2, ef, omega, nom, .false.)
           do i=1,4
              kin = tetc(i,itet) + 1
              kcw(ib,jb-ncg,1:nom,kin) = kcw(ib,jb-ncg,1:nom,kin) + w1t(i,1:nom)*6*wtet(itet)
           enddo
        enddo
     enddo
     !$OMP END PARALLEL DO
  enddo
  deallocate(tetln, wtet)
end subroutine bz_calcqdepw_par

subroutine bz_calcqdepw_par2(kcw, enk, ef, kqid, ncg, nomx, numin, omega, tetc, nkp, nom, nbnd, ntet)
  !
  ! This subroutine calculates the weights for the convolution of k and k-q points. Specifially it calculates the weight of
  !      P_{q,k,om_n} = 2*(f(e_k)- f(-e_{k-q})) / ( iom_n - e_{k-q} + e_k )
  !
  ! if sgnfrq == 3: weights for the Polarization with imaginary frequencies.
  !    2Real{ Integrate[ fw_i(x,y,z) * 1/( 1j*omg - E(x,y,z) ),       {z,0,1},{y,0,1-z},{x,0,1-y-z}] }, where fw_i(x,y,z) = [(1-x-y-z), x, y, z] corresponds to the four corners of tetrahedra.
  !    if iop_bzintq==-1:  numeric integration, which needs x_i,w_i
  !    else:              analytic integration
  !
  ! if sgnfrq == 4:   the q-dependent surface integratio (the surface is defined by e_jb-e_ib=omeg.
  !    Adds surface integration for metallic case
  !
  !    if iop_bzintq== 0:  analytic integration of tetrahedra
  !
  implicit none
  integer,intent(in) :: ncg, nomx, numin, nkp, nbnd, nom, ntet
  integer,intent(in) :: kqid(nkp)
  real(8),intent(in) :: enk(nbnd,nkp), ef
  real(8),intent(in) :: omega(nom)
  integer,intent(in) :: tetc(4,ntet)
  real(8),intent(out):: kcw(ncg+nomx+1,nbnd-ncg-numin, nom, nkp) ! This is way too large array to be processed
  ! FeSe : (80, 307, 1000, 32)  : kcw(io,ie,nkp,nom)
  !         occupied-bands, unoccupied bands, nkp, nom
  real(8), parameter :: pi = 3.14159265358979d0
  !
  external convw1t_par
  !
  integer :: ik, it!, ncbm, nvbm
  integer, allocatable :: tetln(:)         ! index of the q-linked tetrahedron 
  real(8), allocatable :: wtet(:)
  integer :: itet, ib, jb, i, kin, sgnfrq
  real(8) :: ee1(4), ee2(4), w1t(4,nom)
  !logical :: debug
  !
  sgnfrq = 3
  allocate( tetln(ntet), wtet(ntet) )
  do ik=1,nkp
     do it=1,6
        ! if index of k=ik, then index of k-q is kqid[ik]
        tetln(6*(ik-1)+it) = 6*kqid(ik)+it  ! if we know the index of k-q, than we know which tetrahedron is linked with which
        ! they are enumerated in the same order, hence it is just 6*index of k-q point + which tetrahedron in the order
     enddo
  enddo
  wtet(:) = 1.d0/ntet
  !ncbm = numin+1 ! numin,nomx is in c-like type, hence fortran needs +1
  !nvbm = nomx+1  !
  kcw(:,:,:,:)  = 0.0d0
  ! Weighst for he q-dependent bulk integration for Polarization
  ! ee1 is the (partially) occupied band and ee2 is the (partially) empty band
  do itet=1,ntet
     !$OMP PARALLEL DO PRIVATE(jb,i,ee2,ib,ee1,w1t,kin)&
     !$OMP& SHARED(kcw)&
     !$OMP& SCHEDULE(STATIC)
     do jb=1,nbnd-ncg-numin  ! unoccupied band 
        do i=1,4
           ee2(i) = enk( jb+ncg+numin, tetc(i,tetln(itet))+1 )    ! ee2=enk(ncg+numin+1:nbnd) ee1=enk(1:ncg+nomx+1)
        enddo
        if(maxval(ee2) <= ef) cycle ! ee2 is completely occupied =>0, because ee2 should be empty
        do ib=1,ncg+nomx+1  ! occupied band
           do i=1,4
              ee1(i) = enk( ib, tetc(i,itet)+1 )
           enddo
           if( minval(ee1) > ef) cycle  ! ee1 is completely empty => value is zero because ee2 is empty as well
           w1t(1:4,1:nom) = 0.0d0
           call convw1t_par(w1t, ee1, ee2, ef, omega, nom, .false.)
           do i=1,4
              kin = tetc(i,itet) + 1
              kcw(ib,jb,1:nom,kin) = kcw(ib,jb,1:nom,kin) + w1t(i,1:nom)*6*wtet(itet)
           enddo
        enddo
     enddo
     !$OMP END PARALLEL DO
  enddo
  deallocate(tetln, wtet)
end subroutine bz_calcqdepw_par2

subroutine bz_calcqdepw_par3(kcw, enk, ef, kqid, ncg, nomx, numin, dUl, omega, tetc, nkp, nom, nbnd, ntet, nl_svd)
  !
  ! This subroutine calculates the weights for the convolution of k and k-q points. Specifially it calculates the weight of
  !      P_{q,k,om_n} = 2*(f(e_k)- f(-e_{k-q})) / ( iom_n - e_{k-q} + e_k )
  !
  ! if sgnfrq == 3: weights for the Polarization with imaginary frequencies.
  !    2Real{ Integrate[ fw_i(x,y,z) * 1/( 1j*omg - E(x,y,z) ),       {z,0,1},{y,0,1-z},{x,0,1-y-z}] }, where fw_i(x,y,z) = [(1-x-y-z), x, y, z] corresponds to the four corners of tetrahedra.
  !    if iop_bzintq==-1:  numeric integration, which needs x_i,w_i
  !    else:              analytic integration
  !
  ! if sgnfrq == 4:   the q-dependent surface integratio (the surface is defined by e_jb-e_ib=omeg.
  !    Adds surface integration for metallic case
  !
  !    if iop_bzintq== 0:  analytic integration of tetrahedra
  !
  implicit none
  integer,intent(in) :: ncg, nomx, numin, nkp, nbnd, nom, ntet, nl_svd
  integer,intent(in) :: kqid(nkp)
  real(8),intent(in) :: enk(nbnd,nkp), ef
  real(8),intent(in) :: omega(nom)
  integer,intent(in) :: tetc(4,ntet)
  real(8),intent(in) :: dUl(nom,nl_svd)
  real(8),intent(out):: kcw(ncg+nomx+1,nbnd-ncg-numin, nl_svd, nkp)
  !
  real(8), parameter :: pi = 3.14159265358979d0
  !
  external convw1t_par
  !
  integer :: ik, it!, ncbm, nvbm
  integer, allocatable :: tetln(:)         ! index of the q-linked tetrahedron 
  real(8), allocatable :: wtet(:)
  integer :: itet, ib, jb, i, kin, sgnfrq
  real(8) :: ee1(4), ee2(4), w1t(4,nom), w2t(4,nl_svd)
  !logical :: debug
  !
  sgnfrq = 3
  allocate( tetln(ntet), wtet(ntet) )
  do ik=1,nkp
     do it=1,6
        ! if index of k=ik, then index of k-q is kqid[ik]
        tetln(6*(ik-1)+it) = 6*kqid(ik)+it  ! if we know the index of k-q, than we know which tetrahedron is linked with which
        ! they are enumerated in the same order, hence it is just 6*index of k-q point + which tetrahedron in the order
     enddo
  enddo
  wtet(:) = 1.d0/ntet
  !ncbm = numin+1 ! numin,nomx is in c-like type, hence fortran needs +1
  !nvbm = nomx+1  !
  kcw(:,:,:,:)  = 0.0d0
  ! Weighst for he q-dependent bulk integration for Polarization
  ! ee1 is the (partially) occupied band and ee2 is the (partially) empty band
  do itet=1,ntet
     !$OMP PARALLEL DO PRIVATE(jb,i,ee2,ib,ee1,w1t,w2t,kin)&
     !$OMP& SHARED(kcw)&
     !$OMP& SCHEDULE(STATIC)
     do jb=1,nbnd-ncg-numin  ! unoccupied band 
        do i=1,4
           ee2(i) = enk( jb+ncg+numin, tetc(i,tetln(itet))+1 )
        enddo
        if(maxval(ee2) <= ef) cycle ! ee2 is completely occupied =>0, because ee2 should be empty
        do ib=1,ncg+nomx+1  ! occupied band
           do i=1,4
              ee1(i) = enk( ib, tetc(i,itet)+1 )
           enddo
           if( minval(ee1) > ef) cycle  ! ee1 is completely empty => value is zero because ee2 is empty as well
           w1t(1:4,1:nom) = 0.0d0
           call convw1t_par(w1t, ee1, ee2, ef, omega, nom, .false.)
           w2t = matmul(w1t,dUl)  ! w2t(1:4,1:nl_svd) = matmul(w1t(1:4,1:nom),dUl(1:nom,1:nl_svd))
           do i=1,4
              kin = tetc(i,itet) + 1
              kcw(ib,jb,:,kin) = kcw(ib,jb,:,kin) + w2t(i,:)*6*wtet(itet)
           enddo
        enddo
     enddo
     !$OMP END PARALLEL DO
  enddo
  deallocate(tetln, wtet)
end subroutine bz_calcqdepw_par3

