

subroutine convw1t_par(weight, e, f, ef, omega, nom, debug)
  implicit none
  integer, intent(in)  :: nom
  real(8), intent(out) :: weight(4,nom)  ! the weight at each corner
  real(8), intent(in)  :: omega(nom)     ! the imaginary frequency
  real(8), intent(in)  :: e(4)           ! band energies at k
  real(8), intent(in)  :: f(4)           ! band energies at k+q
  real(8), intent(in)  :: ef             ! the fermi energy
  logical, intent(in)  :: debug
  ! if sgnfrq == 1:
  !     Integrate[ fw_i(x,y,z) , {z,0,1},{y,0,1-z},{x,0,1-y-z}]
  !          where fw_i(x,y,z) = [(1-x-y-z), x, y, z] corresponds to the four corners of tetrahedra.
  ! if sgnfrq == 2:
  !     Real{ Integrate[ fw_i(x,y,z) * 1/( omg - E(x,y,z) + 0.01*1j), {z,0,1},{y,0,1-z},{x,0,1-y-z}] }, where fw_i(x,y,z) = [(1-x-y-z), x, y, z] corresponds to the four corners of tetrahedra.
  !                  analytic integrationn                
  ! if sgnfrq == 3:
  !    2Real{ Integrate[ fw_i(x,y,z) * 1/( 1j*omg - E(x,y,z) ),       {z,0,1},{y,0,1-z},{x,0,1-y-z}] }, where fw_i(x,y,z) = [(1-x-y-z), x, y, z] corresponds to the four corners of tetrahedra.
  !                 analytic integration
  ! if sgnfrq == 4:
  !    Adds surface integration for metallic case
  !
  ! This subroutine calculates the contribution of one tetrahedron to the convolution weights in the case that both tetrahedra
  ! (the one at "k" and the linked one at "k-q") are partially occupied. For the case of nnod=7,8,
  ! we further use a big region minus a small region to deal with systematic error. This is for the bulk integration case. 
  !use polyhedron
  !use tetra_internal, only: omgga, sgnfrq
  !
  integer  :: sgnfrq         ! a sign to tell which weight to be calculated
  integer  :: iop_bzintq     ! which mode of calculation
  !
  real(8), parameter :: pi = 3.14159265358979d0
  integer :: inod,t_info
  integer :: isub
  real(8) :: wtemp(4,2,nom) 
  real(8),    allocatable :: t_corners(:,:)
  integer(1), allocatable :: index(:,:) ! it indicates in which order the nodes have to be sent to genericprism
  !
  integer(1) :: ntype(20)      ! idem as pl, but for internal use
  real(8)    :: intnodes(3,20) ! the coordinates of the intersections of the planes
  integer    :: nnod           ! number of nodes defining the polyhedron
  integer    :: info, nd2
  !
  external generictetra_par
  external genericfunf_par
  external genericprism_par
  !
  external surfnodes
  external unrepnodes
  external relnodes
  external sortnodes
  !
  sgnfrq = 3
  iop_bzintq = 0
  !
  !  calculate the intersections between three planes that belongs to the surface of the tetrahedron
  call surfnodes(0, ntype, intnodes, nnod, e, f, ef)
  !  Eliminate repetitions and asign the corresponding value to ntype
  call unrepnodes(0, ntype, intnodes, nnod)
  !  Select those nodes that surround the integration region (e<=ef and f>=ef)
  call relnodes(0,info, ntype, intnodes, nnod, e, f, ef)
  if (info.ne.0) call setnod_debugerr
  !  Sort the nodes in the way needed by convw1t.
  nd2 = 1
  if (nnod > 6) nd2 = 2
  allocate(index(nnod,nd2))
  call sortnodes(info, ntype, nnod, index, nd2)
  if (info.ne.0) call setnod_debugerr

  !if (debug) then
  !   print *, 'nnod=', nnod
  !endif
  select case(nnod)
  case(0)
     weight(1:4,1:nom)=0.0d0
  case(1)
     weight(1:4,1:nom)=0.0d0
  case(2)
     weight(1:4,1:nom)=0.0d0
  case(3)
     weight(1:4,1:nom)=0.0d0
  case(4)     
     allocate(t_corners(1:3,1:nnod))
     t_corners(1:3,1:4) = intnodes(1:3,1:4)
     call generictetra_par(weight, t_corners, 1, t_info,  e, f, omega, ef, nom, debug)
  case(5) 
     allocate(t_corners(1:3,1:nnod))
     do inod=1,5
        t_corners(1:3,index(inod,1)) = intnodes(1:3,inod)
     enddo
     call genericfunf_par(weight, t_corners, 1, t_info, e, f, omega, ef, nom, debug)
  case(6)
     allocate(t_corners(1:3,1:nnod))
     do inod=1,6
        t_corners(1:3,index(inod,1))=intnodes(1:3,inod)
     enddo
     call genericprism_par(weight, t_corners, 1, t_info, e, f, omega, ef, nom, debug)
  case(7)
     allocate(t_corners(1:3,1:5))
     do isub=1,2
        do inod=1,7
           if (index(inod,isub).ne.0) then
              t_corners(1:3,index(inod,isub)) = intnodes(1:3,inod)
           endif
        enddo
        call genericfunf_par(wtemp(:,:,isub), t_corners, 2, t_info, e, f, omega, ef, nom, debug)
     enddo
     weight(1:4,1:nom) = wtemp(1:4,1:nom,1) + wtemp(1:4,1:nom,2)
  case(8)
     allocate(t_corners(1:3,1:6))
     do isub=1,2
        do inod=1,8
           if (index(inod,isub).ne.0) then
              t_corners(1:3,index(inod,isub)) = intnodes(1:3,inod)
           endif
        enddo
        call genericprism_par(wtemp(:,:,isub), t_corners, 4, t_info, e, f, omega, ef, nom, debug)
     enddo
     weight(1:4,1:nom) = wtemp(1:4,1:nom,1) + wtemp(1:4,1:nom,2)
  case default
     write(6,*) "ERROR in convw1t: nnod=",nnod
     stop 'ERROR in convw1t'
  end select
  if(allocated(t_corners))deallocate(t_corners)
  if(allocated(index))deallocate(index)
contains
  subroutine setnod_debugerr
    ! Internal subroutine, repeats the setnodes cycle with write option for debugging
    ! call all the subroutines again with the write option, for debuging
    call surfnodes(1, ntype, intnodes, nnod, e, f, ef)
    call unrepnodes(1, ntype, intnodes, nnod)
    call relnodes(1,info, ntype, intnodes, nnod, e, f, ef)
    stop 'error in setnodes'
  end subroutine setnod_debugerr
end subroutine convw1t_par


subroutine genericfunf_par(w, corners, ical, info, e, f, omega, ef, nom, debug)
  ! This subroutine integrates the convolution weight functions inside a generic 
  ! pentahedron. It divides the pentahedron into the corresponding two tetrahedra 
  ! and calls "generictetra" for each of them
  implicit none
  integer,  intent(in)  :: nom
  real(8),  intent(out) :: w(4,nom)
  integer,  intent(out) :: info
  real(8),  intent(in)  :: omega(nom)
  real(8),  intent(in)  :: corners(3,5)
  integer,  intent(in)  :: ical
  real(8),  intent(in)  :: e(4)           ! band energies at k
  real(8),  intent(in)  :: f(4)           ! band energies at k-q
  real(8),  intent(in)  :: ef             ! the fermi energy
  logical,  intent(in)  :: debug
  !
  integer, parameter :: fout = 6
  integer :: i, itet, inod,tinfo,insp
  real(8) :: twt(4,nom)
  real(8) :: tcorners(3,4)
  insp=0
  info=0
  w(1:4,1:nom) = 0.0d0
  do itet=0,1
     do inod=1,4
        tcorners(1:3,inod) = corners(1:3,inod+itet)
     enddo
     call generictetra_par(twt, tcorners, 5, tinfo, e, f, omega, ef, nom, debug)
     if (tinfo.ne.0) then
        insp = insp + tinfo*(itet+1)
     endif
     w(1:4,1:nom) = w(1:4,1:nom) + twt(1:4,1:nom)
  enddo
  if (insp.ne.0) then
     info=1
     write(fout,'(a7,i4,a8,i4)')'insp = ',insp,' ical = ',ical
     do inod=1,5
        write(fout,'(3f13.8)')corners(inod,1:3)
     enddo
  endif
end subroutine genericfunf_par

subroutine genericprism_par(w, corners, ical, info, e, f, omega, ef, nom, debug)
  ! This subroutine integrates the convolution weight functions inside a generic prism.
  ! It divides the prism into the corresponding three tetrahedra and calls "generictetra" for each of them.
  ! 
  implicit none
  ! The coordinates of the six corners of the prism: the first three form the triangle at the basis and 
  integer, intent(in)  :: ical, nom
  real(8), intent(out) :: w(4,nom)       ! the contribution of the prism to the weight at each corner of the containing tetrahedron.      
  real(8), intent(in)  :: corners(3,6)   ! the last three the triangle at the top, so that the edges of the prism are (1,4), (2,5) and (3,6)
  integer, intent(out) :: info
  real(8), intent(in)  :: e(4)           ! band energies at k
  real(8), intent(in)  :: f(4)           ! band energies at k-q
  real(8), intent(in)  :: omega(nom), ef ! the frequency and fermi energy
  logical,  intent(in)  :: debug
  integer, parameter   :: fout = 6
  !
  integer :: i, itet, inod, tinfo, infl
  real(8) :: tcorners(3,4)
  real(8) :: twt(4,nom)
  info = 0
  infl = 0
  w(1:4,1:nom) = 0.0d0
  do itet = 0,2
     do inod=1,4
        tcorners(1:3,inod) = corners(1:3, inod+itet)
     enddo
     call generictetra_par(twt, tcorners, 6, tinfo, e, f, omega, ef, nom, debug)
     infl = infl + tinfo*(itet+1)
     w(1:4,1:nom) = w(1:4,1:nom) + twt(1:4,1:nom)
  enddo
  if (infl.ne.0) then
     info = 1
     !write(fout,'(a7,i4,a8,i4)')'infl = ',infl,' icap = ',ical
     do inod=1,6
        write(fout,'(3f13.8)')corners(inod,1:3)
     enddo
  endif
end subroutine genericprism_par

subroutine generictetra_par(wt_out, corners, ical, info, e, f, omega, ef, nom, debug)
  ! 2Real{ Integrate[ fw_i(x,y,z) * 1/( 1j*omg - E(x,y,z) ),{z,0,1},{y,0,1-z},{x,0,1-y-z}] },
  !        where fw_i(x,y,z) = [(1-x-y-z), x, y, z] corresponds to the four corners of tetrahedra.
  !        analytic integration
  !
  ! This subroutine calculates the integrals:
  !
  ! For the case of weight including imaginary frequency, $sigfreq=3$:
  !
  ! \begin{eqnarray}
  !   w(1)=\iiint\limits_T \frac{-2(1-x-y-z)\omega}{\omega^2+(\epsilon_{n'k-q}+\epsilon_{nk})^2}
  !     dx dy dz \nonumber \\
  !   w(2)=\iiint\limits_T \frac{-2x\omega}{\omega^2+(\epsilon_{n'k-q}+\epsilon_{nk})^2} 
  !     dx dy dz \nonumber \\
  !   w(3)=\iiint\limits_T \frac{-2y\omega}{\omega^2+(\epsilon_{n'k-q}+\epsilon_{nk})^2} 
  !     dx dy dz \nonumber \\
  !   w(4)=\iiint\limits_T \frac{-2z\omega}{\omega^2+(\epsilon_{n'k-q}+\epsilon_{nk})^2} 
  !     dx dy dz
  ! \end{eqnarray}
  !
  ! where $T$ is a tetrahedron of corners \verb"nodes(i,j)".

  ! For the case of surface integration weight, $sigfreq=4$:
  !
  ! \begin{eqnarray}
  !   w(1)=\iiint\limits_T (1-x-y-z)\Delta(\epsilon_{nk}-\epsilon_{n'k-q}+\omega)
  !     dx dy dz \nonumber \\
  !   w(2)=\iiint\limits_T x\Delta(\epsilon_{nk}-\epsilon_{n'k-q}+\omega)
  !     dx dy dz \nonumber \\
  !   w(3)=\iiint\limits_T y\Delta(\epsilon_{nk}-\epsilon_{n'k-q}+\omega)
  !     dx dy dz \nonumber \\
  !   w(4)=\iiint\limits_T z\Delta(\epsilon_{nk}-\epsilon_{n'k-q}+\omega)
  !     dx dy dz
  ! \end{eqnarray}
  !
  ! where $T$ is a tetrahedron of corners \verb"nodes(i,j)".
  implicit none
  integer, intent(in)   :: nom             ! number of frequencies
  real(8), intent(out)  :: wt_out(4,nom)   ! The four weights corresponding to the original coordinates
  real(8), intent(in)   :: corners(3,4)    ! Coordinates of the four nodes
  real(8), intent(in)   :: omega(nom)      ! the frequency
  integer, intent(in)   :: ical
  integer, intent(out)  :: info
  real(8), intent(in)   :: e(4)           ! band energies at k
  real(8), intent(in)   :: f(4)           ! band energies at k-q
  real(8), intent(in)   :: ef             ! fermi energy
  logical,  intent(in)  :: debug
  !
  interface
     subroutine binary_insertion_sort_indx(indx, array,n)
       integer, intent(out):: indx(n)
       real*8,  intent(in) :: array(n)
       integer, intent(in) :: n
     end subroutine binary_insertion_sort_indx
  end interface
  !
  real(8), parameter :: tol_taylor = 10.0              ! the tolerance for the use of Taylor expansion
  real(8), parameter :: ztol_vol = 1.e-10              ! tolerance for zero volume 
  real(8), parameter :: pi = 3.14159265358979d0
  integer, parameter :: fout = 6
  integer  :: sgnfrq         ! a sign to tell which weight to be calculated
  real(8)  :: vol_small_tetra! 
  integer:: i,j,k
  integer:: ind(4)
  integer:: equiv_flag, iom
  real(8):: vol, det, max_de_small, taylor_omega
  real(8):: vec(3,3) 
  real(8):: delta_e_big_tet(4), delta_e_small_tet(4), etmp(4), delta_e_small_tet_c(4)
  logical :: Qtaylor(nom)
  real(8) :: wt_small_tet(4,nom), wt_tmp(4,nom)
  logical, parameter :: PRINT = .False.
  !
  external sorteq
  !external ksurf
  external stweight_imag_par
  external stweight_itaylor_par
  !
  det(i,j) = vec(2,i) * vec(3,j) - vec(2,j) * vec(3,i)
  !
  sgnfrq = 3
  !
  info=0
  wt_out(:,:) = 0.0d0
  wt_small_tet(:,:) = 0.0d0
  wt_tmp(:,:) = 0.0d0
  ! Calculate the energy differences
  do i=1,4
     delta_e_big_tet(i) = f(i)-e(i)
  enddo
  ! Calculate the volume of the small tetrahedron
  do i=1,3
     do j=1,3
        vec(i,j) = corners(j,i+1)-corners(j,1)
     enddo
  enddo
  vol=0.0d0
  do i=1,3
     j = mod(i,3)+1
     k = mod(j,3)+1
     vol = vol + vec(1,i)*det(j,k)
  enddo
  vol = abs(vol)
  vol_small_tetra = vol
  
  if (vol < ztol_vol) then 
     info = 1
     wt_out(:,:)=0.0d0
     !    If the frequency is zero, the contribution from the Fermi surface has to be calculated
     if ( sgnfrq.eq.4 ) then
        !call binary_insertion_sort_indx(ind, e, 4)
        !do i=1,4
        !   etmp(i) = e(ind(i))
        !enddo
        !call ksurf(wt_small_tet, etmp, ef)
        !do i=1,4
        !   wt_out(ind(i)) = wt_small_tet(i)*(-pi)
        !enddo
     endif
     return 
  endif
  
  ! For frequency dependent weights, calculate the energy diferences at the corners of the small tetrahedron and store the maximum absolute value
  max_de_small=0.0d0
  do i=1,4
     delta_e_small_tet(i) = (delta_e_big_tet(2)-delta_e_big_tet(1)) * corners(1,i) + (delta_e_big_tet(3) - delta_e_big_tet(1)) * corners(2,i) + (delta_e_big_tet(4) - delta_e_big_tet(1)) * corners(3,i) + delta_e_big_tet(1)
     if ( abs(delta_e_small_tet(i) ) > max_de_small )  max_de_small = abs(delta_e_small_tet(i))
  enddo

  taylor_omega = tol_taylor*max_de_small
  do iom=1,nom
     if (omega(iom) > taylor_omega) then
        Qtaylor(iom) = .true.
     else
        Qtaylor(iom) = .false.
     endif
  enddo
  
  select case (sgnfrq)
  case(3)
     ! imaginary axis polarization
     do iom=1,nom
        if (Qtaylor(iom)) call stweight_itaylor_par(delta_e_small_tet, omega(iom), wt_small_tet(:,iom))
     enddo
     call sorteq(delta_e_small_tet, ind, equiv_flag)

     !if (debug) then
     !   print *, delta_e_small_tet, ind
     !endif
     !if (debug) then
     !   print *, 'desmall=', delta_e_small_tet, 'eq_flg=', equiv_flag
     !endif
  
     
     call stweight_imag_par(wt_tmp, delta_e_small_tet, omega, equiv_flag, vol_small_tetra, Qtaylor, nom, debug)
     !call stweight_imag_par_old(wt_tmp, delta_e_small_tet, omega, equiv_flag, vol_small_tetra, Qtaylor, nom, debug)
     do iom=1,nom
        if (.not.Qtaylor(iom)) then
           do i=1,4
              wt_small_tet(ind(i),iom) = wt_tmp(i,iom)
           enddo
        endif
     enddo
  case(4)
     !call binary_insertion_sort_indx(ind, delta_e_small_tet, 4)
     !do i=1,4
     !   delta_e_small_tet_c(i) = delta_e_small_tet(ind(i))
     !enddo
     !call ksurf(wt_tmp, delta_e_small_tet_c, omgga)
     !do i=1,4
     !   wt_small_tet(ind(i)) = wt_tmp(i)*(-pi)
     !enddo
  end select
  do iom=1,nom
     wt_out(1,iom) = sum(wt_small_tet(1:4,iom))
     do i=1,3
        do j=1,4
           wt_out(i+1,iom) = wt_out(i+1,iom) + wt_small_tet(j,iom)*corners(i,j)
        enddo
        wt_out(1,iom) = wt_out(1,iom) - wt_out(i+1,iom)
     enddo
     do i=1,4
        wt_out(i,iom) = wt_out(i,iom)*vol
     enddo
  enddo
  return
end subroutine generictetra_par

subroutine stweight_itaylor_par(deltae_vert, omeg, weight_vert)
  ! This subroutine calculates the weight on the whole small tetrahedron
  ! in which the $k$ states are fully occupied and $k-q$ states are fully 
  ! unoccupied. This is for the $sigfreq=3$ (weights for the Polarization with imaginary frequencies)
  implicit none
  real(8), intent(in) :: deltae_vert(4) ! difference of the energy in k-mesh tetrahedron vertices  and k-q mesh tetrahedron vertices.
  real(8), intent(in) :: omeg           ! the frequency omega to be calculated
  real(8), intent(out):: weight_vert(4) ! the weight on the whole tetrahedron.
  !
  integer :: ivert,j,k
  real(8) :: omt2,omt4,denom1,denom3,w01,n03,w03
  real(8) :: ev(4)
  !
  omt2 = omeg*omeg
  omt4 = omt2*omt2
  denom1 = 6.0d+1*omt2
  denom3 = 4.2d+2*omt4
  do ivert=1,4
     do j=1,4
        k = mod(j+ivert-2,4)+1
        ev(j) = deltae_vert(k)
     enddo
     w01 = -(2.0d0*ev(1)+ev(2)+ev(3)+ev(4))/denom1
     n03 = 4.0d0*ev(1)**3 + 3.0d0*ev(1)**2 * (ev(2)+ev(3)+ev(4))
     n03 = n03 + 2.0d0 * ev(1)*( ev(2)**2 + ev(3)**2 + ev(4)**2 + ev(2)*ev(3) + ev(2)*ev(4) + ev(3)*ev(4) )
     n03 = n03 + ev(2)**3 + ev(3)**3 + ev(4)**3 + ev(2)**2*ev(3) + ev(2)**2*ev(4) + ev(3)**2*ev(4) + ev(2)*ev(3)**2 + ev(2)*ev(4)**2 + ev(3)*ev(4)**2 + ev(2)*ev(3)*ev(4)
     w03 = n03/denom3
     weight_vert(ivert) = w01 + w03
  enddo
  return
end subroutine stweight_itaylor_par

subroutine stweight_imag_par(weight_vert, deltae_vert, omega, equiv_flag, vol_small_tetra, Qtaylor, nom, debug)
  ! This subroutine calculates the weight on the whole small tetrahedron
  ! in which the bands at momentum k are fully occupied and (k+q) states are fully unoccupied.
  ! This is for the "sigfreq=3" imaginary frequency (weights for the Polarization with imaginary frequencies)
  ! 
  !  1/(i*w-eps) = -2*eps/(w^2+eps^2)
  !
  implicit none
  real(8), intent(out):: weight_vert(4,nom) ! the weight on the whole tetrahedron.
  real(8), intent(in) :: deltae_vert(4) ! difference of the energy in k-mesh tetrahedron vertices and k+q mesh tetrahedron vertices.
  real(8), intent(in) :: omega(nom)      ! the frequency omega to be calculated
  logical, intent(in) :: Qtaylor(nom)
  integer, intent(in) :: nom
  integer, intent(in) :: equiv_flag     ! == 4, none is equal 
                                        ! == 6, deltae_vert(1)=deltae_vert(2).
                                        ! == 8, deltae_vert(1)=deltae_vert(2) and deltae_vert(3)=deltae_vert(4).
                                        ! ==10, deltae_vert(1)=deltae_vert(2)=deltae_vert(3).
                                        ! ==16, deltae_vert(1)=deltae_vert(2)=deltae_vert(3)=deltae_vert(4). 
  real(8), intent(in) :: vol_small_tetra
  logical, intent(in) :: debug
  !
  real(8), parameter :: tol_unphys_weight = 1.0e+4   ! the tolerance for unphysical weights
  integer :: j,k,iom
  integer :: ivert
  real(8) :: aa, bb, omeg
  real(8) :: vp 
  real(8) :: ev(4) 
  real(8) :: vdif(4,4)
  real(8) :: omeg2, evdf2, evdf3, evdf4, logs(4), vdf24, tans(4), ev3t4, ev1_2, evf, ev_2(4) 
  real(8) :: cc, cc1, cc2, evt(9), evw(9)

  real*8 :: start, finish
  !real*8, save :: t_time(5) = 0
  
  logical:: lwarn =.false.
  intrinsic datan
  intrinsic dlog
  intrinsic dsign
  !
  weight_vert(1:4,:) = 0.0d0
  select case(equiv_flag)
  case(4)                   ! for the case none of them are equal
     !call cpu_time(start)
     do ivert=1,4
        do j=1,4
           k = mod(j+ivert-2,4) + 1       
           ev(k) = deltae_vert(j)
        enddo
        ! vdif(:,:) = ev(:) - ev(:)
        ! tans(:) = atan( ev(:) / omeg )
        ! logs(:) = log( ev(:)**2 + omeg**2 )
        call sub_set_vdif(vdif, ev)
        ev3t4 = ev(3)*ev(4)
        ev1_2 = ev(1) * ev(1)
        ev_2(1) = 3.0d0*ev1_2
        ev_2(2) = 3.0d0*ev(2)**2
        ev_2(3) = 3.0d0*ev(3)**2
        ev_2(4) = 3.0d0*ev(4)**2
        vdf24 = vdif(2,3) * vdif(2,4) * vdif(3,4)
        evw(1) = 2.0d0 * vdif(1,2) * (vdif(1,3) * vdif(1,4)) * vdf24
        evt(1) = -ev1_2 * evw(1)
        evw(2) = vdf24 * 2.0d0 * ( (2.0d0*ev(1)-ev(2))*(ev(3)+ev(4)) + ev(1)*(2.0d0*ev(2)-3.0d0*ev(1)) - ev(3)*ev(4) )
        evt(2) = vdf24 * 6.0d0 * ev(1) * ( 2.0d0*ev(2)*ev3t4 + ev(1) * (ev1_2- ev3t4 - ev(2)*(ev(3)+ev(4))) )
        evw(3) =  2.0d0 * (vdif(1,3) * vdif(1,4))**2 * vdif(3,4)
        evt(3) = -ev_2(2)*evw(3)
        evw(4) = -2.0d0 * (vdif(1,2) * vdif(1,4))**2 * vdif(2,4)
        evt(4) = -ev_2(3)*evw(4)
        evw(5) = 2.0d0 * (vdif(1,2) * vdif(1,3))**2 * vdif(2,3)
        evt(5) = -ev_2(4) * evw(5)
        evt(6) =  vdf24 * ev1_2 * ( 3.0d0 *  ev(2) * ev3t4 + ev(1)*( ev(2)*(ev(1)-2.0d0*ev(3)) + ev(4)*(ev(1)-2.0d0*ev(2)) + ev(3)*(ev(1)-2.0d0*ev(4)) ) )
        evw(6) = -3.0d0 * vdf24 * ( ev(2)*ev(3)*ev(4) + ev1_2 * (2.0d0*ev(1)-ev(2)-ev(3)-ev(4)) )
        evdf2 = ev(2) * (vdif(1,3) * vdif(1,4))**2 * vdif(3,4)/3.0d0
        evdf3 = ev(3) * (vdif(1,2) * vdif(1,4))**2 * vdif(2,4)/3.0d0
        evdf4 = ev(4) * (vdif(1,2) * vdif(1,3))**2 * vdif(2,3)/3.0d0
        evt(7) = -ev_2(2) * evdf2
        evw(7) =    9.0d0 * evdf2
        evw(8) =  -9.0d0 * evdf3
        evt(8) = ev_2(3) * evdf3
        evw(9) =   9.0d0 * evdf4
        evt(9) = -ev_2(4)* evdf4
        cc = 6.0d0 * (vdif(1,2) * vdif(1,3) * vdif(1,4))**2 * vdif(2,3) * vdif(2,4) * vdif(3,4)
        do iom=1,nom
           if (Qtaylor(iom)) cycle
           omeg = omega(iom)
           call sub_set_tans_logs(tans, logs, ev, omeg)
           omeg2 = omeg**2
           aa = evt(1) + omeg2*evw(1)  + omeg* (  tans(1)*evt(2) + tans(2)*evt(3) + tans(3)*evt(4) + tans(4)*evt(5) + omeg2*( tans(1)*evw(2) + tans(2)*evw(3) + tans(3)*evw(4) + tans(4)*evw(5) ) )
           bb = logs(1)*evt(6) + logs(2)*evt(7) + logs(3)*evt(8) + logs(4)*evt(9) + omeg2 * (logs(1)*evw(6)  + logs(2)*evw(7)  + logs(3)*evw(8)  + logs(4)*evw(9))
           weight_vert(ivert,iom) = (aa+bb)/cc
           call check_weight(weight_vert(ivert,iom), omeg)
        enddo
     enddo
     !call cpu_time(finish)
     !t_time(1) = t_time(1) + finish-start
  case(6)                ! for the case when ev(1)=ev(1)
     !call cpu_time(start)
     do ivert=1,2
        ev(1) = deltae_vert(ivert)
        ev(2) = ev(1)
        ev(3:4) = deltae_vert(3:4)
        call sub_set_vdif(vdif, ev)
        evf = vdif(3,1)*vdif(1,4)*vdif(3,4)
        evw(1) = (2.0d0*ev(1)-ev(3)-ev(4))*evf*2.0d0
        evt(1) = ev(1) * (5.0d0*ev(3)*ev(4) + ev(1) * ( ev(1) -3.0d0*ev(3) -3.0d0*ev(4) ) )*evf
        evt(2) = -6.0d0* ( ev(1)**2 * ( ev(1)*ev(3) + ev(1)*ev(4) -3.0d0*ev(3)*ev(4) )  + (ev(3)*ev(4))**2 ) * vdif(3,4)
        evw(2) =  2.0d0* ( 3.0d0*ev(1)*(ev(1)-ev(3)-ev(4))  + ev(3)**2 + ev(3)*ev(4) + ev(4)**2 ) * vdif(3,4)
        evw(3) = -2.0d0*vdif(1,4)**3
        evt(3) = -3.0d0*ev(3)**2*evw(3)
        evw(4) =  2.0d0*vdif(1,3)**3
        evt(4) = -3.0d0*ev(4)**2*evw(4)
        evt(5) = -vdif(3,4) * ev(1)*( 3.0d0*ev(3)*ev(4)*( ev(1)*(ev(1)-ev(3)-ev(4)) + ev(3)*ev(4) ) + ev(1)**2*(ev(3)-ev(4))**2 )
        evw(5) = -3.0d0*vdif(3,4) * ( ev(3)*ev(4)*(3.0d0*ev(1)-ev(3)-ev(4)) - ev(1)**3 )
        evf = ev(3)*vdif(1,4)**3
        evw(6) = -3.0d0    * evf
        evt(6) = ev(3)**2  * evf
        evf = ev(4)*vdif(1,3)**3
        evw(7) =  3.0d0    * evf
        evt(7) = -ev(4)**2 * evf
        cc = 6.0d0*vdif(1,3)**3*vdif(1,4)**3*vdif(3,4)
        do iom=1,nom
           if (Qtaylor(iom)) cycle
           omeg = omega(iom)
           call sub_set_tans_logs(tans, logs, ev, omeg)
           omeg2 = omeg**2
           aa = evt(1) + omeg2 * evw(1) + omeg * ( tans(1)*evt(2)  + tans(3)*evt(3)  + tans(4)*evt(4) + omeg2 * (tans(1)*evw(2) + tans(3)*evw(3) + tans(4)*evw(4)) )
           bb = logs(1)*evt(5) + logs(3)*evt(6) + logs(4)*evt(7)  + omeg2*( logs(1)*evw(5)  + logs(3)*evw(6) + logs(4)*evw(7))
           weight_vert(ivert,iom) = (aa+bb)/cc
           call check_weight(weight_vert(ivert,iom), omeg)
        enddo
     enddo ! ivert  
     do ivert=3,4
        ev(1) = (deltae_vert(1)+deltae_vert(2))*0.5d0
        ev(2) = ev(1)
        do j=3,4
           k=mod(j+ivert,2)+3
           ev(k)=deltae_vert(j)
        enddo
        call sub_set_vdif(vdif, ev)
        evt(1) = ( ev(3)**2*ev(4) + ev(1)**2*vdif(4,3)-ev(1)*ev(3)**2 ) * 2.0d0*vdif(1,3)*vdif(1,4)*vdif(3,4)
        evw(1) = ( (ev(1)+ev(3)-2.0d0*ev(4)) ) * 2.0d0*vdif(1,3)*vdif(1,4)*vdif(3,4)
        evt(2) = (3.0d0*ev(1)*( ev(1)*(ev(3)+ev(1)) - 2.0d0*ev(3)*ev(4) )) * 2.0d0*vdif(3,4)**2
        evw(2) = ( ev(3) + 2.0d0*ev(4) - 3.0d0*ev(1) ) * 2.0d0*vdif(3,4)**2
        evt(3) = -3.0d0*ev(3)*(ev(3)**2+ev(1)*(ev(3)-2.0d0*ev(4)))*2.0d0*vdif(1,4)**2
        evw(3) = -(ev(1)-3.0d0*ev(3)+2.0d0*ev(4))*2.0d0*vdif(1,4)**2
        evw(4) = 2.0d0*vdif(1,3)**3
        evt(4) = (-3.0d0*ev(4)**2)*evw(4)
        evt(5) = ev(1)**2*(-3.0d0*ev(3)*ev(4)+ev(1)*(2.0d0*ev(3)+ev(4))) * vdif(3,4)**2
        evw(5) = (-6.0d0*ev(1)**2+3.0d0*ev(1)*ev(4)+3.0d0*ev(3)*ev(4)) * vdif(3,4)**2
        evt(6) = -ev(3)**2*(2.0d0*ev(1)*ev(3)-3.0d0*ev(1)*ev(4)+ev(3)*ev(4))* vdif(1,4)**2
        evw(6) = -3.0d0*(-2.0d0*ev(3)**2+(ev(1)+ev(3))*ev(4)) * vdif(1,4)**2
        evt(7) = -ev(4)*vdif(1,3)**3 * ev(4)**2
        evw(7) =  ev(4)*vdif(1,3)**3 * 3.0d0
        cc = 6.0d0*vdif(1,3)**3*vdif(1,4)**2*vdif(3,4)**2
        do iom=1,nom
           if (Qtaylor(iom)) cycle
           omeg = omega(iom)
           omeg2 = omeg**2
           call sub_set_tans_logs(tans, logs, ev, omeg)
           aa = evt(1) + omeg2 * evw(1) + omeg*( tans(1)*evt(2) + tans(3)*evt(3) + tans(4)*evt(4) + omeg2 * (tans(1)*evw(2)  + tans(3)*evw(3)  + tans(4)*evw(4)) )
           bb = logs(1)*evt(5) + logs(3)*evt(6) + logs(4)*evt(7) + omeg2*( logs(1)*evw(5) + logs(3)*evw(6) + logs(4)*evw(7) )
           weight_vert(ivert,iom) = (aa+bb)/cc
           call check_weight(weight_vert(ivert,iom), omeg)
        enddo
     enddo ! ivert
     !call cpu_time(finish)
     !t_time(2) = t_time(2) + finish-start
  case(8)          !for the case when ev(1)=ev(1) and ev(3)=ev(4)
     !call cpu_time(start)
     ev(1) = (deltae_vert(1)+deltae_vert(2))/2.d0
     ev(2) = ev(1)
     ev(3) = (deltae_vert(3)+deltae_vert(4))/2.d0
     ev(4) = ev(3)
     call sub_set_vdif(vdif, ev)
     evw(1) = 6.0d0*vdif(3,1)
     evt(1) = (ev(1)**2-5.0d0*ev(1)*ev(3)-2.0d0*ev(3)**2)*vdif(3,1)
     evt(2) = -6.0d0*( ev(3)**2 + 2.0d0*ev(1)*ev(3))
     evw(3) = 3.0d0*(ev(1)+2.0d0*ev(3))
     evt(3) = -3.0d0*ev(1)*ev(3)**2
     evw(4) = 6.0d0*vdif(1,3)
     evt(4) = -(2.0d0*ev(1)**2+5.0d0*ev(1)*ev(3)-ev(3)**2)*vdif(1,3)
     evt(5) = 6.0d0*(ev(1)**2 + 2.0d0*ev(1)*ev(3))
     evw(6) = -3.0d0*(2.0d0*ev(1) + ev(3))
     evt(6) =  3.0d0*ev(1)**2*ev(3)
     cc = 6.0d0*vdif(1,3)**4
     do iom=1,nom
        if (Qtaylor(iom)) cycle
        omeg = omega(iom)
        omeg2 = omeg**2
        call sub_set_tans_logs(tans, logs, ev, omeg)
        aa = evt(1) + evw(1)*omeg2 + omeg*(tans(1)-tans(3))*( 6.0d0*omeg2 + evt(2) )
        bb =  (logs(1)-logs(3))*(evt(3) + evw(3)*omeg2)
        weight_vert(1,iom) = (aa+bb)/cc
        call check_weight(weight_vert(1,iom), omeg)
        weight_vert(2,iom) = weight_vert(1,iom)
        aa = evt(4) + evw(4)*omeg2 + omeg*(tans(1)-tans(3))*(-6.0d0*omeg2 + evt(5))
        bb = (logs(1)-logs(3)) * (evt(6) + evw(6)*omeg2)
        weight_vert(3,iom) = (aa+bb)/cc
        call check_weight(weight_vert(3,iom), omeg)
        weight_vert(4,iom) = weight_vert(3,iom)
     enddo
     !call cpu_time(finish)
     !t_time(3) = t_time(3) + finish-start
  case(10)             ! for the case when ev(1)=ev(1)=ev(3)
     !call cpu_time(start)
     ev(1:3) = sum(deltae_vert(1:3))/3.0d0
     ev(4) = deltae_vert(4)
     call sub_set_vdif(vdif, ev)
     evw(1) = 6.0d0*vdif(1,4)
     evt(1) = -vdif(1,4)*( 2.0d0*ev(1)**2 - 7.0d0*ev(1)*ev(4) + 11.0d0*ev(4)**2)
     evt(2) = 18.0d0*ev(4)**2
     evw(3) = -9.0d0*ev(4)
     evt(3) =  3.0d0*ev(4)**3
     evw(4) = 6.0d0 * vdif(4,1)
     evt(4) = vdif(4,1)*( ev(1)**2 - 5.0d0*ev(1)*ev(4) - 2.0d0*ev(4)**2 )
     evt(5) = -6.0d0*ev(4)*( ev(4) + 2.0d0*ev(1) )
     evt(6) = -3.0d0*ev(1)*ev(4)**2
     evw(6) =  3.0d0*(ev(1)+2.0d0*ev(4))
     cc1 = 18.0d0*vdif(1,4)**4
     cc2 = 6.0d0*vdif(1,4)**4
     do iom=1,nom
        if (Qtaylor(iom)) cycle
        omeg = omega(iom)
        omeg2 = omeg**2
        call sub_set_tans_logs(tans, logs, ev, omeg)
        aa = evt(1) + omeg2*evw(1) + omeg*(tans(1)-tans(4)) * ( -6.0d0*omeg2 + evt(2))
        bb = (logs(1)-logs(4)) * (evt(3) + evw(3)*omeg2)
        weight_vert(1,iom) = (aa+bb)/cc1
        call check_weight(weight_vert(1,iom), omeg)
        weight_vert(2:3,iom) = weight_vert(1,iom)
        aa = evt(4) + evw(4)*omeg2 + omeg*(tans(1)-tans(4)) * ( 6.0d0*omeg2 + evt(5))
        bb = (logs(1)-logs(4)) * (evt(6) + evw(6)*omeg2 )
        weight_vert(4,iom) = (aa+bb)/cc2
        call check_weight(weight_vert(4,iom), omeg)
     enddo
     !call cpu_time(finish)
     !t_time(4) = t_time(4) + finish-start
  case(16)
     !call cpu_time(start)
     ev(1:4) = sum(deltae_vert(1:4))/4.0
     do iom=1,nom
        if (Qtaylor(iom)) cycle
        omeg = omega(iom)
        weight_vert(1,iom) = -ev(1)/(12.0d0*(omeg**2+ev(1)**2))
        call check_weight(weight_vert(1,iom), omeg)
        weight_vert(2:4,iom) = weight_vert(1,iom)
     enddo
     !call cpu_time(finish)
     !t_time(5) = t_time(5) + finish-start
  case default
     write(6,*)'ERROR in stweight_imat'
     write(6,*)'wrong equiv_flag: ',equiv_flag
     stop "ERROR in stweight_imat"
  end select
  !print *, t_time(:)
contains
  subroutine check_weight(wght_vert, omeg)
    real(8), intent(inout) :: wght_vert
    real(8), intent(in)    :: omeg
    if ( abs(wght_vert) > tol_unphys_weight ) then
       if ( abs(omeg) < 1.e-6 ) then 
          wght_vert = 0.0 
       else
          write(6,'(A)') 'WARNING in stweight_imag: unexpected big weight!'
          write(6,'(A,g16.6,A,I4)') 'weightt =', wght_vert, ' case:', equiv_flag
          write(6,'(A,g16.6)') 'omeg =', omeg
          write(6,'(A,4g16.6)') 'deltae_vert =', deltae_vert(:)
          write(6,'(A,4(4g16.6,/))') 'vdif = ', vdif
          write(6,'(A,g16.6,A,g16.6,A,g16.6)') ' aa =', aa,' bb =', bb, ' cc =', cc
          write(6,'(A,g16.6)') ' vol_small_tetra=', vol_small_tetra
       endif
    endif
  end subroutine check_weight
  subroutine sub_set_vdif(vdif, ev)
    ! vdif(:,:) = ev(:) - ev(:)
    implicit none
    real(8), intent(out) :: vdif(4,4)
    real(8), intent(in)  :: ev(4)
    integer(4) :: ii, jj
    do ii=1,4
       do jj=1,4
          vdif(ii,jj)=ev(ii)-ev(jj)
       enddo
    enddo
  end subroutine sub_set_vdif
  subroutine sub_set_tans_logs(tans, logs, ev, omeg)
    ! tans(:) = atan( ev(:) / omeg )
    ! logs(:) = log( ev(:)**2 + omeg**2 )
    implicit none
    real(8), intent(out) :: tans(4), logs(4)
    real(8), intent(in)  :: ev(4), omeg
    integer(4) :: ii
    do ii=1,4
       if ( omeg.gt.1.0e-20 ) then
          tans(ii) = datan( ev(ii)/omeg )
       else
          tans(ii) = dsign( 1.0d0, ev(ii) ) * 2.0d0 * datan(1.0d0)
       endif
       if ( ( omeg**2 + ev(ii)**2 ) > 1.0d-10 ) then
          logs(ii) = dlog( ev(ii)**2 + omeg**2 )
       else 
          logs(ii) = 0.0d0
       endif
    enddo
  end subroutine sub_set_tans_logs
end subroutine stweight_imag_par


subroutine stweight_imag_par_old(weight_vert, deltae_vert, omega, equiv_flag, vol_small_tetra, Qtaylor, nom, debug)
  ! This subroutine calculates the weight on the whole small tetrahedron
  ! in which the bands at momentum k are fully occupied and (k+q) states are fully unoccupied.
  ! This is for the "sigfreq=3" imaginary frequency (weights for the Polarization with imaginary frequencies)
  ! 
  !  1/(i*w-eps) = -2*eps/(w^2+eps^2)
  !
  implicit none
  real(8), intent(out):: weight_vert(4,nom) ! the weight on the whole tetrahedron.
  real(8), intent(in) :: deltae_vert(4) ! difference of the energy in k-mesh tetrahedron vertices and k+q mesh tetrahedron vertices.
  real(8), intent(in) :: omega(nom)      ! the frequency omega to be calculated
  logical, intent(in) :: Qtaylor(nom)
  integer, intent(in) :: nom
  integer, intent(in) :: equiv_flag     ! == 4, none is equal 
                                        ! == 6, deltae_vert(1)=deltae_vert(2).
                                        ! == 8, deltae_vert(1)=deltae_vert(2) and deltae_vert(3)=deltae_vert(4).
                                        ! ==10, deltae_vert(1)=deltae_vert(2)=deltae_vert(3).
                                        ! ==16, deltae_vert(1)=deltae_vert(2)=deltae_vert(3)=deltae_vert(4). 
  real(8), intent(in) :: vol_small_tetra
  logical, intent(in) :: debug
  !
  real(8), parameter :: tol_unphys_weight = 1.0e+4   ! the tolerance for unphysical weights
  integer :: j,k,iom
  integer :: ivert
  real(8) :: aa, bb, cc, dd, bb1, bb3, bb4, omeg
  real(8) :: vp 
  real(8) :: ev(4), tans(4), logs(4)
  real(8) :: vdif(4,4)
  logical:: lwarn =.false.
  intrinsic datan
  intrinsic dlog
  intrinsic dsign
  !
  weight_vert(1:4,:) = 0.0d0
  select case(equiv_flag)
  case(4)                   ! for the case none of them are equal
     do ivert=1,4
        do j=1,4
           k = mod(j+ivert-2,4) + 1       
           ev(k) = deltae_vert(j)
        enddo
        ! vdif(:,:) = ev(:) - ev(:)
        ! tans(:) = atan( ev(:) / omeg )
        ! logs(:) = log( ev(:)**2 + omeg**2 )
        call sub_set_vdif(vdif, ev)
        do iom=1,nom
           if (Qtaylor(iom)) cycle
           omeg = omega(iom)
           call sub_set_tans_logs(tans, logs, ev, omeg)
           aa = 2.0d0 * (omeg**2-ev(1)**2) * vdif(1,2) * vdif(1,3) * vdif(1,4) * vdif(2,3) * vdif(2,4) * vdif(3,4)
           dd = 2.0d0 * omeg * ( 3.0d0*ev(1)**4 - omeg**2 * (ev(3)*ev(4) + ev(2)*(ev(3)+ev(4))) - 3.0d0 * ev(1)**2 * (omeg**2+ev(3)*ev(4)+ev(2)*(ev(3)+ev(4))) + 2.0d0*ev(1)*(omeg**2*(ev(3)+ev(4))+ev(2)*(omeg**2+3.0d0*ev(3)*ev(4))))
           aa = aa + dd * vdif(2,3) * vdif(2,4) * vdif(3,4) * tans(1)
           aa = aa + 2.0d0 * omeg * (omeg**2 - 3.0d0*ev(2)**2) * vdif(1,3)**2 * vdif(1,4)**2 * vdif(3,4) * tans(2)
           aa = aa - 2.0d0 * omeg * (omeg**2 - 3.0d0*ev(3)**2) * vdif(1,2)**2 * vdif(1,4)**2 * vdif(2,4) * tans(3)
           aa = aa + 2.0d0 * omeg * (omeg**2 - 3.0d0*ev(4)**2) * vdif(1,2)**2 * vdif(1,3)**2 * vdif(2,3) * tans(4)
           dd = -3.0d0 * omeg**2 * ev(2)*ev(3)*ev(4) + ev(1)**4 * (ev(2)+ev(3)+ev(4))
           dd = dd - 2.0d0 * ev(1)**3 * ( 3.0d0*omeg**2 + ev(3)*ev(4) + ev(2)*(ev(3)+ev(4)) )
           dd = dd + 3.0d0 * ev(1)**2 * ( omeg**2 * ( ev(3) + ev(4) ) + ev(2)*( omeg**2 + ev(3) * ev(4) ) )
           bb = dd * vdif(2,3) * vdif(2,4) * vdif(3,4) * logs(1)
           bb = bb + ev(2) * ( 3.0d0*omeg**2 - ev(2)**2) * vdif(1,3)**2 * vdif(1,4)**2 * vdif(3,4) * logs(2)
           bb = bb + ev(3) * (-3.0d0*omeg**2 + ev(3)**2) * vdif(1,2)**2 * vdif(1,4)**2 * vdif(2,4) * logs(3) 
           bb = bb - ev(4) * (-3.0d0*omeg**2 + ev(4)**2) * vdif(1,2)**2 * vdif(1,3)**2 * vdif(2,3) * logs(4)
           cc = 6.0d0 * vdif(1,2)**2 * vdif(1,3)**2 * vdif(1,4)**2 * vdif(2,3) * vdif(2,4) * vdif(3,4)
           weight_vert(ivert,iom) = (aa+bb)/cc
           call check_weight(weight_vert(ivert,iom), omeg)
        enddo
     enddo
  case(6)                ! for the case when ev(1)=ev(1)
     do ivert=1,2
        ev(1) = deltae_vert(ivert)
        ev(2) = ev(1)
        ev(3:4) = deltae_vert(3:4)
        call sub_set_vdif(vdif, ev)
        do iom=1,nom
           if (Qtaylor(iom)) cycle
           omeg = omega(iom)
           call sub_set_tans_logs(tans, logs, ev, omeg)
           dd = ev(1)**3 - 2.0d0*omeg**2*(ev(3)+ev(4)) - 3.0d0*ev(1)**2*(ev(3)+ev(4)) + ev(1)*(4.0d0*omeg**2+5.0d0*ev(3)*ev(4))
           aa = vdif(3,1)*vdif(1,4)*vdif(3,4)*dd
           dd = -omeg**2*(3.0d0*ev(1)**2+ev(3)**2+ev(3)*ev(4)+ev(4)**2-3.0d0*ev(1)*(ev(3)+ev(4)))+3.0d0*(-3.0d0*ev(1)**2*ev(3)*ev(4)+ev(3)**2*ev(4)**2+ev(1)**3*(ev(3)+ev(4)))
           aa = aa - 2.0d0*omeg*vdif(3,4)*dd*tans(1)
           aa = aa - 2.0d0*omeg*(omeg**2-3.0d0*ev(3)**2)*vdif(1,4)**3*tans(3)
           aa = aa + 2.0d0*omeg*vdif(1,3)**3*(omeg**2-3.0d0*ev(4)**2)*tans(4)
           dd = -3.0d0*omeg**2*ev(3)*ev(4)*(ev(3)+ev(4))-3.0d0*ev(1)**2*ev(3)*ev(4)*(ev(3)+ev(4))+3.0d0*ev(1)*ev(3)*ev(4)*(3.0d0*omeg**2+ev(3)*ev(4))+ev(1)**3*(-3.0d0*omeg**2+ev(3)**2+ev(3)*ev(4)+ev(4)**2)
           bb1 = -vdif(3,4)*dd
           bb3 = ev(3)*vdif(1,4)**3*(ev(3)**2-3.0d0*omeg**2)
           bb4 = -ev(4)*vdif(1,3)**3*(ev(4)**2-3.0d0*omeg**2)
           bb = bb1*logs(1)+bb3*logs(3)+bb4*logs(4)
           cc = 6.0d0*vdif(1,3)**3*vdif(1,4)**3*vdif(3,4)
           weight_vert(ivert,iom) = (aa+bb)/cc
           call check_weight(weight_vert(ivert,iom), omeg)
        enddo
     enddo ! ivert  
     do ivert=3,4
        ev(1) = (deltae_vert(1)+deltae_vert(2))*0.5d0
        ev(2) = ev(1)
        do j=3,4
           k=mod(j+ivert,2)+3
           ev(k)=deltae_vert(j)
        enddo
        call sub_set_vdif(vdif, ev)
        do iom=1,nom
           if (Qtaylor(iom)) cycle
           omeg = omega(iom)
           call sub_set_tans_logs(tans, logs, ev, omeg)
           dd = -ev(1)*ev(3)**2+omeg**2*(ev(1)+ev(3)-2.0d0*ev(4))+ev(3)**2*ev(4)+ev(1)**2*vdif(4,3)
           aa = 2.0d0*dd*vdif(1,3)*vdif(1,4)*vdif(3,4)
           dd = 3.0d0*ev(1)**3+3.0d0*ev(1)**2*ev(3)+omeg**2*(ev(3)+2.0d0*ev(4))-3.0d0*ev(1)*(omeg**2+2.0d0*ev(3)*ev(4))
           aa = aa + 2.0d0*vdif(3,4)**2*omeg*dd*tans(1)
           dd = 3.0d0*ev(3)*(ev(3)**2+ev(1)*(ev(3)-2.0d0*ev(4)))+omeg**2*(ev(1)-3.0d0*ev(3)+2.0d0*ev(4))
           aa = aa - 2.0d0*vdif(1,4)**2*omeg*dd*tans(3)
           aa = aa + 2.0d0*vdif(1,3)**3*omeg*(omeg**2-3.0d0*ev(4)**2)*tans(4)
           dd = ev(1)**2*(-3.0d0*ev(3)*ev(4)+ev(1)*(2.0d0*ev(3)+ev(4)))
           dd = dd + omeg**2*(-6.0d0*ev(1)**2+3.0d0*ev(1)*ev(4)+3.0d0*ev(3)*ev(4))
           bb = vdif(3,4)**2*dd*logs(1)
           dd = ev(3)**2*(2.0d0*ev(1)*ev(3)-3.0d0*ev(1)*ev(4)+ev(3)*ev(4))
           dd = dd + 3.0d0*omeg**2*(-2.0d0*ev(3)**2+(ev(1)+ev(3))*ev(4))
           bb = bb - vdif(1,4)**2*dd*logs(3)
           bb = bb - ev(4)*vdif(1,3)**3*(ev(4)**2-3.0d0*omeg**2)*logs(4)
           cc = 6.0d0*vdif(1,3)**3*vdif(1,4)**2*vdif(3,4)**2
           weight_vert(ivert,iom) = (aa+bb)/cc
           call check_weight(weight_vert(ivert,iom), omeg)
        enddo
     enddo ! ivert  
  case(8)          !for the case when ev(1)=ev(1) and ev(3)=ev(4)
     ev(1) = (deltae_vert(1)+deltae_vert(2))/2.d0
     ev(2) = ev(1)
     ev(3) = (deltae_vert(3)+deltae_vert(4))/2.d0
     ev(4) = ev(3)
     call sub_set_vdif(vdif, ev)
     do iom=1,nom
        if (Qtaylor(iom)) cycle
        omeg = omega(iom)
        call sub_set_tans_logs(tans, logs, ev, omeg)
        dd = 6.0d0*omeg**2+ev(1)**2-5.0d0*ev(1)*ev(3)-2.0d0*ev(3)**2
        aa = vdif(3,1)*dd
        dd = 6.0d0*omeg*(omeg**2 - ev(3)**2 - 2.0d0*ev(1)*ev(3))
        aa = aa + dd*(tans(1)-tans(3))
        bb = 3.0d0*((ev(1)+2.0d0*ev(3))*omeg**2 - ev(1)*ev(3)**2 )
        bb = bb*(logs(1)-logs(3))
        cc = 6.0d0*vdif(1,3)**4
        weight_vert(1,iom) = (aa+bb)/cc
        call check_weight(weight_vert(1,iom), omeg)
        weight_vert(2,iom) = weight_vert(1,iom)
        dd = 6.0d0*omeg**2-2.0d0*ev(1)**2-5.0d0*ev(1)*ev(3)+ev(3)**2
        aa = vdif(1,3)*dd
        dd = 6.0d0*omeg*(omeg**2 - ev(1)**2 - 2.0d0*ev(1)*ev(3))
        aa = aa + dd*(tans(3)-tans(1))
        dd = 3.0d0*( (2.0d0*ev(1) + ev(3))*omeg**2 - ev(1)**2*ev(3) )
        bb = dd*(logs(3)-logs(1))
        cc = 6.0d0*vdif(1,3)**4
        weight_vert(3,iom) = (aa+bb)/cc
        call check_weight(weight_vert(3,iom), omeg)
        weight_vert(4,iom) = weight_vert(3,iom)
     enddo
  case(10)             ! for the case when ev(1)=ev(1)=ev(3)
     ev(1:3) = sum(deltae_vert(1:3))/3.0d0
     ev(4) = deltae_vert(4)
     call sub_set_vdif(vdif, ev)
     do iom=1,nom
        if (Qtaylor(iom)) cycle
        omeg = omega(iom)
        call sub_set_tans_logs(tans, logs, ev, omeg)
        aa = vdif(1,4)*(6.0d0*omeg**2-2.0d0*ev(1)**2+7.0d0*ev(1)*ev(4)- 11.0d0*ev(4)**2)     
        dd = 6.0d0*omeg*(omeg**2-3.0d0*ev(4)**2)
        aa = aa - dd * (tans(1)-tans(4))
        dd = 3.0d0*ev(4)*(ev(4)**2-3.0d0*omeg**2)
        bb = dd*(logs(1)-logs(4))
        cc = 18.0d0*vdif(1,4)**4
        weight_vert(1,iom) = (aa+bb)/cc
        call check_weight(weight_vert(1,iom), omeg)
        weight_vert(2:3,iom) = weight_vert(1,iom) 
        dd = 6.0d0*omeg**2+ev(1)**2-5.0d0*ev(1)*ev(4)-2.0d0*ev(4)**2
        aa = vdif(4,1)*dd
        dd = 6.0d0*omeg*(omeg**2-ev(4)**2-2.0d0*ev(1)*ev(4))
        aa = aa+dd*(tans(1)-tans(4))
        dd = -3.0d0*ev(1)*ev(4)**2+3.0d0*omeg**2*(ev(1)+2.0d0*ev(4))
        bb = dd*(logs(1)-logs(4))
        cc = 6.0d0*vdif(1,4)**4
        weight_vert(4,iom) = (aa+bb)/cc
        call check_weight(weight_vert(4,iom), omeg)
     enddo
  case(16)
     ev(1:4) = sum(deltae_vert(1:4))/4.0
     do iom=1,nom
        if (Qtaylor(iom)) cycle
        omeg = omega(iom)
        aa = -ev(1)
        bb = 0.0d0
        cc = 12.0d0*(omeg**2+ev(1)**2)
        ivert=1
        weight_vert(1,iom) = (aa+bb)/cc !  -ev/(12*(om**2+ev**2)
        call check_weight(weight_vert(1,iom), omeg)
        weight_vert(2:4,iom)=weight_vert(1,iom)
     enddo
  case default
     write(6,*)'ERROR in stweight_imat'
     write(6,*)'wrong equiv_flag: ',equiv_flag
     stop "ERROR in stweight_imat"
  end select

contains
  subroutine check_weight(wght_vert, omeg)
    real(8), intent(inout) :: wght_vert
    real(8), intent(in)    :: omeg
    if ( abs(wght_vert) > tol_unphys_weight ) then
       if ( abs(omeg) < 1.e-6 ) then 
          wght_vert = 0.0 
       else
          write(6,'(A)') 'WARNING in stweight_imag: unexpected big weight!'
          write(6,'(A,g16.6,A,I4)') 'weightt =', wght_vert, ' case:', equiv_flag
          write(6,'(A,g16.6)') 'omeg =', omeg
          write(6,'(A,4g16.6)') 'deltae_vert =', deltae_vert(:)
          write(6,'(A,4(4g16.6,/))') 'vdif = ', vdif
          write(6,'(A,g16.6,A,g16.6,A,g16.6)') ' aa =', aa,' bb =', bb, ' cc =', cc
          write(6,'(A,g16.6)') ' vol_small_tetra=', vol_small_tetra
       endif
    endif
  end subroutine check_weight
  subroutine sub_set_vdif(vdif, ev)
    ! vdif(:,:) = ev(:) - ev(:)
    implicit none
    real(8), intent(out) :: vdif(4,4)
    real(8), intent(in)  :: ev(4)
    integer(4) :: ii, jj
    do ii=1,4
       do jj=1,4
          vdif(ii,jj)=ev(ii)-ev(jj)
       enddo
    enddo
  end subroutine sub_set_vdif
  subroutine sub_set_tans_logs(tans, logs, ev, omeg)
    ! tans(:) = atan( ev(:) / omeg )
    ! logs(:) = log( ev(:)**2 + omeg**2 )
    implicit none
    real(8), intent(out) :: tans(4), logs(4)
    real(8), intent(in)  :: ev(4), omeg
    integer(4) :: ii
    do ii=1,4
       if ( omeg.gt.1.0e-20 ) then
          tans(ii) = datan( ev(ii)/omeg )
       else
          tans(ii) = dsign( 1.0d0, ev(ii) ) * 2.0d0 * datan(1.0d0)
       endif
       if ( ( omeg**2 + ev(ii)**2 ) > 1.0d-10 ) then
          logs(ii) = dlog( ev(ii)**2 + omeg**2 )
       else 
          logs(ii) = 0.0d0
       endif
    enddo
  end subroutine sub_set_tans_logs
end subroutine stweight_imag_par_old
