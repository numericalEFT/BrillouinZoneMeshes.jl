!----------- Routines for simple tetrahedron method over irreducible tetrahedra --------
subroutine intweight1t(w, et,ef)
  ! This subroutine calculates the contribution to the integration weights at the four corners of a tetrahedron.
  ! The energies have to be given in increasing energy order.     
  implicit none
  real*8, intent(out):: w(4) ! Weight at each corner
  real*8, intent(in) :: et(4) ! Eigenenergies at the corners of the tetrahedron.
  real*8, intent(in) :: ef    ! Fermi energy.
  ! locals
  integer :: index, i
  real*8  :: c1,c2,c3,vo4
  real*8  :: f31,f41,f32,f42
  real*8  :: e31,e41,e32,e42
  if( et(4) <= ef) then  ! all energies above ef
     index = 4
  else
     index=0
     do i=1,4             ! how many corners of tetrahedra below ef
        if( et(i) <= ef) index = index + 1
     enddo
  endif
  vo4 = 1.0d0/4.0d0
  select case(index)
  case(0)          ! all states are unoccupied
     w(:) = 0.0d0
  case(1)          ! only the lowest energy is occupied.
     do i=2,4
        w(i) = (ef-et(1))/(et(i)-et(1))
     enddo
     c1 = w(2) * w(3) * w(4)
     w(2) = c1 * w(2)
     w(3) = c1 * w(3)
     w(4) = c1 * w(4)
     w(1) = 4.0d0 * c1 - w(2) - w(3) - w(4)
  case(2)         ! the two lower energies are occupied
     f31 = (ef-et(1))/(et(3)-et(1))
     f41 = (ef-et(1))/(et(4)-et(1))
     f32 = (ef-et(2))/(et(3)-et(2))
     f42 = (ef-et(2))/(et(4)-et(2))
     e31 = (et(3)-ef)/(et(3)-et(1))
     e41 = (et(4)-ef)/(et(4)-et(1))
     e32 = (et(3)-ef)/(et(3)-et(2))
     e42 = (et(4)-ef)/(et(4)-et(2))
     c1 = f41 * f31
     c2 = f41 * f32 * e31
     c3 = f42 * f32 * e41
     w(1) =  c1 + (c1+c2)*e31 + (c1+c2+c3)*e41 
     w(2) =  c1+c2+c3 + (c2+c3)*e32 + c3*e42 
     w(3) =  (c1+c2)*f31 + (c2+c3)*f32 
     w(4) =  (c1+c2+c3)*f41 + c3*f31 
  case(3)        ! Only the highest energy is unoccupied
     do i=1,3
        w(i) = (et(4)-ef)/(et(4)-et(i))
     enddo
     c1   = w(1) * w(2) * w(3) 
     w(4) = 1.0d0 - c1 * (4.0d0 - w(1)-w(2)-w(3))
     w(1) = 1.0d0 - c1 * w(1)
     w(2) = 1.0d0 - c1 * w(2)
     w(3) = 1.0d0 - c1 * w(3)
  case(4)        ! All states are occupied
     w(:) = 1.d0
  end select
  w(:) = w(:)*vo4
end subroutine intweight1t

subroutine bloechlcor(blcor,et,ef)
  !   This subroutine calculates the correction to the integration weights for
  !   the improved tetrahedron method according to
  !   Bl\"ochl {\it et al}, Phys. Rev. B, {\bf 49}, 16223 (1994). The energies at the 
  !   corner of the tetrahedron have to be given in increasing order.
  implicit none
  real*8, intent(out):: blcor(4) ! The correction to the integration weights
  real*8, intent(in) :: et(4)    ! Eigenenergies at the corners of the tetrahedron.
  real*8, intent(in) :: ef       ! the fermi energy 
  !
  integer   :: i, j, jm4, index
  real*8    :: dte
  real*8    :: sde
  real*8, external :: dos1t
  intrinsic mod
  !
  dte = dos1t(et,ef)/4.0d+1
  blcor(:)=0.0d0
  if(et(4) <= ef) then
     index=4
  else
     index=0
     do i=1,4
        if (et(i) <= ef) index = index + 1
     enddo
  endif
  select case(index)
  case(0,4)
     continue
  case(1,2,3)  
     do i=1,4
        sde = 0.0d0
        do j=i,i+2
           jm4=mod(j,4)+1
           sde=sde +et(jm4)
        enddo
        sde=sde-3.0d0*et(i)
        blcor(i)=sde*dte
     enddo
  end select
end subroutine bloechlcor

real*8 function intdos1t(et,ei)
  ! This subroutine calculates the contribution to the integrated density of states at an energy $\epsilon = $ \verb"ei"
  ! of one tetrahedron. The energies at the corner of the tetrahedron have to be given in increasing order.   
  !
  implicit none
  real*8, intent(in) :: et(4) ! Eigenenergies at the corners of the tetrahedron.
  real*8, intent(in) :: ei    ! energy at which the contribution to the DOS is calculated
  ! locals
  integer :: index, i
  real*8  :: br,denom,eij
  real*8  :: e21,e31,e41,e32,e42,e43
  !
  if(et(4) <= ei) then
     index = 4
  else
     index = 0
     do i=1,4
        if(et(i) <= ei) index = index + 1
     enddo
  endif
  intdos1t = 0.0d0
  select case(index)
  case(0)          ! all states are above ei
     intdos1t = 0.0d0
  case(1)          ! e1 < ei < e2.
     e21 = et(2)-et(1)
     e31 = et(3)-et(1)
     e41 = et(4)-et(1)
     eij = ei   -et(1)
     denom = e21*e31*e41
     intdos1t = eij*eij*eij/denom
  case(2)         ! e2 < ei < e3
     e21 = et(2)-et(1)
     e31 = et(3)-et(1)
     e41 = et(4)-et(1)
     e32 = et(3)-et(2)
     e42 = et(4)-et(2)
     eij = ei   -et(2)
     denom = e32*e42
     br = e21*e21 + 3.0d0*eij*e21 + 3.0d0*eij*eij - (e31+e42)*eij*eij*eij/denom
     denom = e31*e41
     intdos1t = br/denom
  case(3)        ! e3 < ei < e4
     e43 = et(4)-et(3)
     e42 = et(4)-et(2)
     e41 = et(4)-et(1)
     eij = et(4)-ei
     denom = e41*e42*e43
     intdos1t = (1.0d0-eij*eij*eij/denom)
  case(4)        ! e4 < ei
     intdos1t = 1.0d0
  end select
end function intdos1t

real*8 function dos1t(et,ei)
  !   This subroutine calculates the contribution to the density of states at
  !   an energy $\epsilon = $ \verb"ei" of one tetrahedron, according to
  !   Bl\"ochl {\it et al}, Phys. Rev. B, {\bf 49}, 16223 (1994). The energies at the 
  !   corner of the tetrahedron have to be given in increasing order. 
  implicit none
  real*8, intent(in) :: et(4) ! Eigenenergies at the corners of the tetrahedron.
  real*8, intent(in) :: ei    ! energy at which the contribution to the dos1t is calculated
  !
  integer(4) :: index,i
  real(8)    :: br,denom,eij
  real(8)    :: e21,e31,e41,e32,e42,e43
  if(et(4) <= ei)then
     index=4
  else
     index=0
     do i=1,4
        if(et(i) <= ei) index = index + 1
     enddo
  endif
  dos1t = 0.0d0
  select case(index)
  case(0,4)          ! all states are either below or above ei
     dos1t = 0.0d0
  case(1)          ! only the lowest energy is occupied.
     e21 = et(2) - et(1)
     e31 = et(3) - et(1)
     e41 = et(4) - et(1)
     eij = ei    - et(1)
     denom = e21*e31*e41
     dos1t = 3.0d0*eij*eij/denom
  case(2)         ! e2 < ei < e3
     e21 = et(2) - et(1)
     e31 = et(3) - et(1)
     e41 = et(4) - et(1)
     e32 = et(3) - et(2)
     e42 = et(4) - et(2)
     eij = ei    - et(2)
     denom = e32*e42
     br = e21+2.0d0*eij-(e31+e42)*eij*eij/denom
     denom = e31*e41
     dos1t = 3.0d0*br/denom
  case(3)        ! e3 < ei < e4
     e43 = et(4) - et(3)
     e42 = et(4) - et(2)
     e41 = et(4) - et(1)
     eij = et(4) - ei
     denom = e41*e42*e43
     dos1t = 3.0d0*eij*eij/denom
  end select
end function dos1t

subroutine ksurf(wgh, et, ef)
  ! This subroutine calculates the integration over an energy surface in the
  ! k-mesh, the area of the Fermi surface inside the tetrahedron is calculated
  ! which is the essential item to decide the weight of the integration over
  ! each vertex.
  ! The energies have to be given in increasing energy order.
  implicit none
  real*8, intent(out) :: wgh(4)    ! the weight at each corner
  real*8, intent(in)  :: et(4)     ! the band energies or energy difference at k
  real*8, intent(in)  :: ef        ! the energy or energy differnce surface
  !
  integer :: i, indx  ! indx=Number of nodes below ef
  real*8  :: delta21, delta31, delta41, delta32, delta42, delta43  ! energy differences       
  real*8  :: deles1, deles2, deles3, deles4 ! ef-eo(i) 
  real*8  :: den0         ! denominators
  real*8  :: num0, num1   ! numerators
  real*8  :: wt           ! total weight
  !
  wgh(:) = 0.0d0
  if ( ef < et(1) ) then
     indx=0
  else
     indx=1
     do i=2,4
        if( et(i) <= ef) indx = indx + 1
     enddo
  endif
  select case(indx)
  case(0,4) 
     continue
  case(1) ! Only et(1)< ef, one triangle
     deles1  =  ef - et(1)
     delta21 = et(2)-et(1)
     delta31 = et(3)-et(1)
     delta41 = et(4)-et(1)
     den0 = 6.0d0*delta21*delta31*delta41
     num0 = deles1*deles1
     ! total wgh
     wt = 3.0d0*num0/den0
     num1 = num0*deles1
     wgh(2) = num1/(delta21*den0)
     wgh(3) = num1/(delta31*den0)
     wgh(4) = num1/(delta41*den0)
     wgh(1) = wt - wgh(2) - wgh(3) - wgh(4)
  case(2) ! et(1)<et(2)<ef, two triangles
     deles1  =  ef  - et(1)
     deles2  =  ef  - et(2)
     deles3  = et(3) - ef
     deles4  = et(4) - ef
     delta21 = et(2) - et(1)
     delta31 = et(3) - et(1)
     delta41 = et(4) - et(1)
     delta32 = et(3) - et(2)
     delta42 = et(4) - et(2)
     delta43 = et(4) - et(3)
     ! total wgh          
     num0 = deles1*deles3 / (6.0d0*delta31*delta32*delta41)
     num1 = deles2*deles4 / (6.0d0*delta32*delta41*delta42)
     wt = 3.0d0 * (num0 + num1)
     !
     wgh(2) = deles3/delta32 * num0   + (deles3*delta42+deles2*delta31)/(delta41*delta42) * num1
     wgh(3) = (deles1*delta32+deles2*delta31)/(delta31*delta32) * num0 + deles2/(delta32) * num1
     wgh(4) = deles1/(delta41) * num0 + (deles1*delta42+deles2*delta41)/(delta41*delta42) * num1
     wgh(1) = wt - wgh(2)-wgh(3)-wgh(4)
  case(3) ! ee(1), ee(2) and ee(3) < efer, one triangle.
     deles4  = et(4) - ef
     delta41 = et(4) - et(1)
     delta42 = et(4) - et(2)
     delta43 = et(4) - et(3)
     den0 = 6.0d0*delta41*delta42*delta43
     if ( den0 < 1.0e-12) then
        write(6,*)'ERROR in ksurf: flat band'
        write(6,*) et
        stop 'flat band'
     endif
     num0 = deles4*deles4
     ! total wgh
     wt   = num0* 3.0d0/den0
     num1 = num0*deles4/den0
     wgh(1) = num1/delta41
     wgh(2) = num1/delta42
     wgh(3) = num1/delta43
     wgh(4) = wt - wgh(1)-wgh(2)-wgh(3)
  end select
end subroutine ksurf



!----------- Routines for convolution tetrahedron method over k and k-q points --------

subroutine convw(ef,omgga,sgnfrq,cwt, iop_bzintq, enk, wtet, tetc, tetln, vol_small_tetra, nbnd, nkp, ntet, x_i, w_i, n_i)
  !   This subroutine calculates the integration weight of each k-point for all band pairs.
  !   
  ! if sgnfrq == 1:   normal q-dependent bulk integration. 
  !     Integrate[ fw_i(x,y,z) , {z,0,1},{y,0,1-z},{x,0,1-y-z}]
  !          where fw_i(x,y,z) = [(1-x-y-z), x, y, z] corresponds to the four corners of tetrahedra.
  !
  ! if sgnfrq == 2:  weights for the Polarization with real frequencies
  !     Real{ Integrate[ fw_i(x,y,z) * 1/( omg - E(x,y,z) + 0.01*1j), {z,0,1},{y,0,1-z},{x,0,1-y-z}] }, where fw_i(x,y,z) = [(1-x-y-z), x, y, z] corresponds to the four corners of tetrahedra.
  !    if iop_bzintq==-1:  numeric integration, which needs x_i,w_i
  !    else:               analytic integrationn
  !
  ! if sgnfrq == 3: weights for the Polarization with imaginary frequencies.
  !    2Real{ Integrate[ fw_i(x,y,z) * 1/( 1j*omg - E(x,y,z) ),       {z,0,1},{y,0,1-z},{x,0,1-y-z}] }, where fw_i(x,y,z) = [(1-x-y-z), x, y, z] corresponds to the four corners of tetrahedra.
  !    if iop_bzintq==-1:  numeric integration, which needs x_i,w_i
  !    else:              analytic integration
  !
  ! if sgnfrq == 4:   the q-dependent surface integratio (the surface is defined by e_jb-e_ib=omeg.
  !    Adds surface integration for metallic case
  !
  !   E(k)==ee1==e is the (partially) occupied band and E(k-q)==ee2==f is the (partially) empty band
  !use order
  !use tetra_internal
  implicit none
  integer, intent(in)  :: nbnd, ntet, nkp
  real(8), intent(out) :: cwt(nbnd,nbnd,nkp) ! the weight of each k-point for each band 
  real(8), intent(in)  :: ef                     ! fermi energy
  real(8), intent(in)  :: omgga
  integer, intent(in)  :: sgnfrq
  !
  real*8,  intent(in)  :: enk(nbnd,nkp)  ! Band energies
  integer, intent(in)  :: tetc(4,ntet)   ! tetrahedra index
  real*8,  intent(in)  :: wtet(ntet)     ! weight of each tetrahedra
  integer, intent(in)  :: tetln(ntet)    ! link between k and k-q tetrahedra
  real*8,  intent(out) :: vol_small_tetra
  integer, intent(in)  :: iop_bzintq     ! which mode of calculation
  real(8), intent(in)  :: x_i(n_i), w_i(n_i) ! x-positions and weights for gauss integration quadrature in 1D
  integer, intent(in)  :: n_i
  !
  integer :: itet, i, ib, jb, kin
  integer :: ik1(4), ik2(4) 
  real(8) :: term, wwwt
  real(8) :: ee1(4), ee2(4), w1t(4), wcor(4) 
  external intweight1t
  external convw1t
  external convw1tsurf
  external bloechlcor
  !
  intrinsic maxval
  intrinsic minval       
  !
  ! tetln, vol_small_tetra, iop_bzintq, x_i, w_i, n_i

  cwt  = 0.0d0
  wwwt = 0.0d0
  select case (sgnfrq)
  case(1)      ! normal q-dependent bulk integration

     ! locals : itet, ib, i, ee1, jb, ee2, kin, w1t, wcor, ik1, vol_small_tetra, 
     ! cwt should be accumulative we mak cwt shared
     ! shared enk, tetc, wtet, x_i, w_i, n_i
     
     do itet=1,ntet
        do ib=1,nbnd
           do i=1,4
              ee1(i) = enk( ib,tetc(i,itet)+1 )
           enddo
           if( maxval(ee1) < ef) then   ! ee1 is fully occupied
              do jb=1,nbnd          
                 do i=1,4
                    ee2(i)=enk( jb, tetc(i,tetln(itet))+1 )
                 enddo
                 if( minval(ee2) > ef) then   ! ee2 is completely empty
                    do i=1,4
                       kin = tetc(i,itet) + 1
                       cwt(ib,jb,kin) = cwt(ib,jb,kin) + wtet(itet)/4.0d0
                    enddo
                 else if ( maxval(ee2) > ef) then   ! ee2 partially occupied
                    ee2(1:4)  = ef-ee2(1:4) ! now turn to occupied with patially empty case, so that standard bz-integration works
                    w1t(1:4)  = 0.0d0
                    wcor(1:4) = 0.0d0
                    call binary_insertion_sort(ee2,4)  !call sort(4,ee2,ik2)
                    call intweight1t(w1t, ee2, 0.0d0)
                    call bloechlcor(wcor, ee2, 0.0d0)
                    do i=1,4
                       kin = tetc(i,itet) + 1
                       cwt(ib,jb,kin) = cwt(ib,jb,kin) + (w1t(i)+wcor(i))*wtet(itet)
                    enddo
                 endif
              enddo
           else if(minval(ee1) <= ef) then   ! ee1 is partially occupied
              do jb=1,nbnd
                 do i=1,4
                    ee2(i) = enk( jb, tetc(i,tetln(itet))+1 )
                 enddo
                 if ( minval(ee2) > ef ) then ! ee2 is completely empty
                    w1t(1:4)  = 0.0d0
                    wcor(1:4) = 0.0d0
                    call binary_insertion_sort_indx(ik1, ee1, 4)   ! call sort(4,ee1,ik1)
                    call intweight1t(w1t, ee1, ef )
                    call bloechlcor(wcor, ee1, ef )
                    do i=1,4
                       kin = tetc( ik1(i), itet) + 1
                       cwt(ib,jb,kin) = cwt(ib,jb,kin) + ( w1t(i) + wcor(i) ) * wtet(itet)
                    enddo
                 else                       ! ee2 is also partially empty, hence both ee1 and ee2 are partially occupied
                    w1t(1:4)=0.0d0
                    call convw1t(ee1,ee2,ef,w1t, omgga, sgnfrq, vol_small_tetra, iop_bzintq, x_i, w_i, n_i)
                    do i=1,4
                       kin = tetc(i,itet) + 1
                       cwt(ib,jb,kin) = cwt(ib,jb,kin) + w1t(i)*6*wtet(itet)
                    enddo
                 endif
              enddo
           endif
        enddo
     enddo  ! itet
  case(2:4)     ! for the q-dependent bulk integration for Polarization
     ! locals : itet, ib, i, ee1, jb, ee2, kin, w1t, vol_small_tetra, 
     ! cwt should be accumulative we mak cwt shared
     ! shared enk, tetc, wtet, x_i, w_i, n_i

     ! ee1 is the (partially) occupied band and ee2 is the (partially) empty band
     do itet=1,ntet
        !$OMP PARALLEL DO PRIVATE(ib, i, ee1, jb, ee2, kin, w1t, vol_small_tetra)&
        !$OMP& SHARED(cwt)&
        !$OMP& SCHEDULE(STATIC)
        do ib=1,nbnd  ! occupied band
           do i=1,4
              ee1(i) = enk( ib, tetc(i,itet)+1 )
           enddo
           if(minval(ee1) > ef) cycle  ! ee1 is completely empty => value is zero because ee2 is empty as well
           !write(6,'(A,I4,2x,4I4)') 'tetcorn=', itet, ((tetc(i,itet)+1),i=1,4)
           !write(6,'(A,F7.4,1x,A,I4,1x,A,I3,1x,A,4F10.7)') 'om=', omgga, 'itet=', itet, 'ib=', ib, 'ee1=', ee1(:)
           do jb=1,nbnd  ! unoccupied band 
              do i=1,4
                 ee2(i) = enk( jb, tetc(i,tetln(itet))+1 )
              enddo
              if(maxval(ee2) <= ef) cycle ! ee2 is completely occupied =>0, because ee2 should be empty
              !write(6,'(A,3I4,2x,4F10.7,2x,4F10.7)') 'ee=', itet, ib, jb, ee1(:), ee2(:)
              w1t(1:4) = 0.0d0
              call convw1t(ee1,ee2,ef,w1t, omgga, sgnfrq, vol_small_tetra, iop_bzintq, x_i, w_i, n_i)
              do i=1,4
                 kin = tetc(i,itet) + 1
                 cwt(ib,jb,kin) = cwt(ib,jb,kin) + w1t(i)*6*wtet(itet)
              enddo
              !write(6,'(A,F7.4,1x,A,I4,1x,A,I3,1x,A,4F10.7,1x,A,4F12.7)') 'om=', omgga, 'itet=', itet, 'jb=', jb, 'ee1=', ee2(:), 'w=', w1t(:)*(6*wtet(itet))
           enddo
        enddo
        !$OMP END PARALLEL DO
     enddo
  case default
     continue
  end select
  return
end subroutine convw

subroutine convw1t(e, f, ef, weight, &
     & omgga, sgnfrq, vol_small_tetra, iop_bzintq, x_i, w_i, n_i)
  ! if sgnfrq == 1:
  !     Integrate[ fw_i(x,y,z) , {z,0,1},{y,0,1-z},{x,0,1-y-z}]
  !          where fw_i(x,y,z) = [(1-x-y-z), x, y, z] corresponds to the four corners of tetrahedra.
  ! if sgnfrq == 2:
  !     Real{ Integrate[ fw_i(x,y,z) * 1/( omg - E(x,y,z) + 0.01*1j), {z,0,1},{y,0,1-z},{x,0,1-y-z}] }, where fw_i(x,y,z) = [(1-x-y-z), x, y, z] corresponds to the four corners of tetrahedra.
  !    if iop_bzintq==-1:  numeric integration, which needs x_i,w_i
  !    else:               analytic integrationn                
  ! if sgnfrq == 3:
  !    2Real{ Integrate[ fw_i(x,y,z) * 1/( 1j*omg - E(x,y,z) ),       {z,0,1},{y,0,1-z},{x,0,1-y-z}] }, where fw_i(x,y,z) = [(1-x-y-z), x, y, z] corresponds to the four corners of tetrahedra.
  !    if iop_bzintq==-1:  numeric integration, which needs x_i,w_i
  !    ellse:              analytic integration
  ! if sgnfrq == 4:
  !    Adds surface integration for metallic case
  !
  ! This subroutine calculates the contribution of one tetrahedron to the convolution weights in the case that both tetrahedra
  ! (the one at "k" and the linked one at "k-q") are partially occupied. For the case of nnod=7,8,
  ! we further use a big region minus a small region to deal with systematic error. This is for the bulk integration case. 
  !use polyhedron
  !use tetra_internal, only: omgga, sgnfrq
  implicit none
  real(8), intent(in)  :: e(4)           ! band energies at k
  real(8), intent(in)  :: f(4)           ! band energies at k-q
  real(8), intent(out) :: weight(4)      ! the weight at each corner
  real(8), intent(in)  :: omgga, ef      ! the frequency and fermi energy
  integer, intent(in)  :: sgnfrq         ! a sign to tell which weight to be calculated
  real(8), intent(out) :: vol_small_tetra! 
  integer, intent(in)  :: iop_bzintq     ! which mode of calculation
  real(8), intent(in)  :: x_i(n_i), w_i(n_i) ! x-positions and weights for gauss integration quadrature in 1D
  integer, intent(in)  :: n_i
  !
  real(8), parameter :: pi = 3.14159265358979d0
  integer :: i,inod,t_info
  integer :: isub
  real(8) :: wtemp(4,2) 
  real(8),    allocatable :: t_corners(:,:)
  integer(1), allocatable :: index(:,:) ! it indicates in which order the nodes have to be sent to genericprism
  real(8) :: etmp(4), wtmp(4)
  integer :: tindx(4)
  !
  integer(1) :: ntype(20)      ! idem as pl, but for internal use
  real(8)    :: intnodes(3,20) ! the coordinates of the intersections of the planes
  integer    :: nnod           ! number of nodes defining the polyhedron
  integer    :: info, nd2
  !
  !external setnodes
  external generictetra
  external genericfunf
  external genericprism
  !
  external surfnodes
  external unrepnodes
  external relnodes
  external sortnodes
  !
  !call setnodes(ntype, intnodes, nnod,    e, f, ef)
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
  
  select case(nnod)
  case(0)
     weight(1:4)=0.0d0
  case(1)
     weight(1:4)=0.0d0
  case(2)
     weight(1:4)=0.0d0
  case(3)
     weight(1:4)=0.0d0
     if (sgnfrq.eq.4) then
        if (dabs(omgga) < 1.0d-12) then
           call binary_insertion_sort_indx(tindx, e, 4)
           do i=1,4
              etmp(i) = e(tindx(i))
           enddo
           call ksurf(wtmp, etmp, ef)
           do i=1,4
              weight(tindx(i)) = wtmp(i)*(-pi)
           enddo
        endif
     endif
  case(4)     
     allocate(t_corners(1:3,1:nnod))
     t_corners(1:3,1:4) = intnodes(1:3,1:4)
     call generictetra(t_corners, weight, 1, t_info,  vol_small_tetra, e, f, omgga, ef, sgnfrq, iop_bzintq, x_i, w_i, n_i)
  case(5) 
     allocate(t_corners(1:3,1:nnod))
     do inod=1,5
        t_corners(1:3,index(inod,1)) = intnodes(1:3,inod)
     enddo
     call genericfunf(t_corners, weight, 1, t_info, vol_small_tetra, e, f, omgga, ef, sgnfrq, iop_bzintq, x_i, w_i, n_i)
  case(6)
     allocate(t_corners(1:3,1:nnod))
     do inod=1,6
        t_corners(1:3,index(inod,1))=intnodes(1:3,inod)
     enddo
     call genericprism(t_corners, weight, 1, t_info, vol_small_tetra, e, f, omgga, ef, sgnfrq, iop_bzintq, x_i, w_i, n_i)
  case(7)
     allocate(t_corners(1:3,1:5))
     do isub=1,2
        do inod=1,7
           if (index(inod,isub).ne.0) then
              t_corners(1:3,index(inod,isub)) = intnodes(1:3,inod)
           endif
        enddo ! inod
        call genericfunf(t_corners, wtemp(1:4,isub), 2, t_info, vol_small_tetra, e, f, omgga, ef, sgnfrq, iop_bzintq, x_i, w_i, n_i)
     enddo ! isub  
     weight(1:4) = wtemp(1:4,1) + wtemp(1:4,2)
  case(8)
     allocate(t_corners(1:3,1:6))
     do isub=1,2
        do inod=1,8
           if (index(inod,isub).ne.0) then
              t_corners(1:3,index(inod,isub)) = intnodes(1:3,inod)
           endif
        enddo
        call genericprism(t_corners, wtemp(1:4,isub), 4, t_info, vol_small_tetra, e, f, omgga, ef, sgnfrq, iop_bzintq, x_i, w_i, n_i)
     enddo
     weight(1:4)=wtemp(1:4,1)+wtemp(1:4,2)
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
  
end subroutine convw1t

subroutine genericfunf(corners, w, ical, info, &
     & vol_small_tetra, e, f, omgga, ef, sgnfrq, iop_bzintq, x_i, w_i, n_i)
  ! This subroutine integrates the convolution weight functions inside a generic 
  ! pentahedron. It divides the pentahedron into the corresponding two tetrahedra 
  ! and calls "generictetra" for each of them
  implicit none
  real(8),  intent(out) :: w(4)
  integer,  intent(out) :: info
  real(8),  intent(in)  :: corners(3,5)
  integer,  intent(in)  :: ical
  real(8), intent(out)  :: vol_small_tetra! 
  real(8), intent(in)   :: e(4)           ! band energies at k
  real(8), intent(in)   :: f(4)           ! band energies at k-q
  real(8), intent(in)   :: omgga, ef      ! the frequency and fermi energy
  integer, intent(in)   :: sgnfrq         ! a sign to tell which weight to be calculated
  integer, intent(in)   :: iop_bzintq     ! which mode of calculation
  real(8), intent(in)   :: x_i(n_i), w_i(n_i) ! x-positions and weights for gauss integration quadrature in 1D
  integer, intent(in)   :: n_i
  !
  integer, parameter :: fout = 6
  integer :: i, itet, inod,tinfo,insp
  real(8) :: twt(4)
  real(8) :: tcorners(3,4)
  insp=0
  info=0
  w(1:4) = 0.0d0
  do itet = 0,1
     do inod=1,4
        tcorners(1:3,inod) = corners(1:3,inod+itet)
     enddo
     call generictetra(tcorners, twt, 5, tinfo,   vol_small_tetra, e, f, omgga, ef, sgnfrq, iop_bzintq, x_i, w_i, n_i)
     if (tinfo.ne.0) then
        insp = insp + tinfo*(itet+1)
     endif
     do i=1,4
        w(i) = w(i) + twt(i)
     enddo
  enddo
  if (insp.ne.0) then
     info=1
     write(fout,'(a7,i4,a8,i4)')'insp = ',insp,' ical = ',ical
     do inod=1,5
        write(fout,'(3f13.8)')corners(inod,1:3)
     enddo
  endif
end subroutine genericfunf

subroutine genericprism(corners, w, ical, info, &
     & vol_small_tetra, e, f, omgga, ef, sgnfrq, iop_bzintq, x_i, w_i, n_i)
  ! This subroutine integrates the convolution weight functions inside a generic prism.
  ! It divides the prism into the corresponding three tetrahedra and calls "generictetra" for each of them.
  ! 
  implicit none
  ! The coordinates of the six corners of the prism: the first three form the triangle at the basis and 
  real(8), intent(in)  :: corners(3,6)   ! the last three the triangle at the top, so that the edges of the prism are (1,4), (2,5) and (3,6)
  real(8), intent(out) :: w(4)           ! the contribution of the prism to the weight at each corner of the containing tetrahedron.      
  integer, intent(in)  :: ical
  integer, intent(out) :: info
  real(8), intent(out) :: vol_small_tetra! 
  real(8), intent(in)  :: e(4)           ! band energies at k
  real(8), intent(in)  :: f(4)           ! band energies at k-q
  real(8), intent(in)  :: omgga, ef      ! the frequency and fermi energy
  integer, intent(in)  :: sgnfrq         ! a sign to tell which weight to be calculated
  integer, intent(in)  :: iop_bzintq     ! which mode of calculation
  real(8), intent(in)  :: x_i(n_i), w_i(n_i) ! x-positions and weights for gauss integration quadrature in 1D
  integer, intent(in)  :: n_i
  integer, parameter   :: fout = 6
  !
  integer :: i, itet, inod, tinfo, infl
  real(8) :: tcorners(3,4)
  real(8) :: twt(4)
  info = 0
  infl = 0
  w(1:4) = 0.0d0
  do itet = 0,2
     do inod=1,4
        tcorners(1:3,inod) = corners(1:3, inod+itet)
     enddo
     call generictetra(tcorners, twt, 6, tinfo,   vol_small_tetra, e, f, omgga, ef, sgnfrq, iop_bzintq, x_i, w_i, n_i)
     infl = infl + tinfo*(itet+1)
     do i=1,4
        w(i) = w(i) + twt(i)
     enddo
  enddo
  if (infl.ne.0) then
     info = 1
     !write(fout,'(a7,i4,a8,i4)')'infl = ',infl,' icap = ',ical
     do inod=1,6
        write(fout,'(3f13.8)')corners(inod,1:3)
     enddo
  endif
end subroutine genericprism

subroutine generictetra(corners, wt_out, ical, info, &
     & vol_small_tetra, e, f, omgga, ef, sgnfrq, iop_bzintq, x_i, w_i, n_i)
  !
  ! if sgnfrq == 1:
  !     Integrate[ fw_i(x,y,z) , {z,0,1},{y,0,1-z},{x,0,1-y-z}]
  !          where fw_i(x,y,z) = [(1-x-y-z), x, y, z] corresponds to the four corners of tetrahedra.
  ! if sgnfrq == 2:
  !     Real{ Integrate[ fw_i(x,y,z) * 1/( omg - E(x,y,z) + 0.01*1j), {z,0,1},{y,0,1-z},{x,0,1-y-z}] }, where fw_i(x,y,z) = [(1-x-y-z), x, y, z] corresponds to the four corners of tetrahedra.
  !    if iop_bzintq==-1:  numeric integration, which needs x_i,w_i
  !    else:               analytic integrationn                
  ! if sgnfrq == 3:
  !    2Real{ Integrate[ fw_i(x,y,z) * 1/( 1j*omg - E(x,y,z) ),       {z,0,1},{y,0,1-z},{x,0,1-y-z}] }, where fw_i(x,y,z) = [(1-x-y-z), x, y, z] corresponds to the four corners of tetrahedra.
  !    if iop_bzintq==-1:  numeric integration, which needs x_i,w_i
  !    ellse:              analytic integration
  ! if sgnfrq == 4:
  !    Adds surface integration for metallic case
  !
  ! This subroutine calculates the integrals:
  ! 
  ! For the case of normal weight $sigfreq=1$:
  ! \begin{equation}
  ! \begin{align}
  !   w(1)=&\iiint\limits_T (1-x-y-z) dx dy dz \\
  !   w(2)=&\iiint\limits_T x dx dy dz \\
  !   w(3)=&\iiint\limits_T y dx dy dz \\ 
  !   w(4)=&\iiint\limits_T z dx dy dz 
  ! \end{align}
  ! \end{equation}
  ! 
  ! For the case of the weight including real frequency contribution, 
  ! $sigfreq=2$:
  !
  ! \begin{eqnarray}
  !   w(1)=\iiint\limits_T \frac{1-x-y-z}{\omega-\epsilon_{n'k-q}+\epsilon_{nk}} 
  !     dx dy dz \nonumber \\
  !   w(2)=\iiint\limits_T \frac{x}{\omega-\epsilon_{n'k-q}+\epsilon_{nk}} 
  !     dx dy dz \nonumber \\
  !   w(3)=\iiint\limits_T \frac{y}{\omega-\epsilon_{n'k-q}+\epsilon_{nk}} 
  !     dx dy dz \nonumber \\
  !   w(4)=\iiint\limits_T \frac{z}{\omega-\epsilon_{n'k-q}+\epsilon_{nk}} 
  !     dx dy dz
  ! \end{eqnarray}
  !
  ! The whole contribution is got by setting $\omega$ to $-\omega$ and call this
  ! program again.
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
  !use bzint,         only: ztol_vol,iop_bzintq,tol_taylor
  !use tetra_internal,only: sgnfrq, omgga,fout,vol_small_tetra,ldbg_bzint
  implicit none
  real(8), intent(in)   :: corners(3,4) ! Coordinates of the four nodes
  integer, intent(in)   :: ical
  integer, intent(out)  :: info
  real(8), intent(out)  :: wt_out(4)      ! The four weights corresponding to the original coordinates
  real(8), intent(out)  :: vol_small_tetra! 
  real(8), intent(in)   :: e(4)           ! band energies at k
  real(8), intent(in)   :: f(4)           ! band energies at k-q
  real(8), intent(in)   :: omgga, ef      ! the frequency and fermi energy
  integer, intent(in)   :: sgnfrq         ! a sign to tell which weight to be calculated
  integer, intent(in)   :: iop_bzintq     ! which mode of calculation
  real(8), intent(in)   :: x_i(n_i), w_i(n_i) ! x-positions and weights for gauss integration quadrature in 1D
  integer, intent(in)   :: n_i
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
  integer:: i,j,k
  integer:: ind(4)
  integer:: sigeq
  real(8):: vol, det, max_de_small
  real(8):: vec(3,3) 
  real(8):: delta_e_big_tet(4), delta_e_small_tet(4),wt_small_tet(4), wt_tmp(4),etmp(4), delta_e_small_tet_c(4)
  logical, parameter :: PRINT = .False.
  logical, parameter :: debug = .False.
  !
  ! !EXTERNAL ROUTINES: 
  !
  external ksurf
  external sorteq
  external stweight_imag
  external stweight_itaylor
  external stweight_real
  external stweight_rtaylor
  external stweight_numeric
  !
  det(i,j) = vec(2,i) * vec(3,j) - vec(2,j) * vec(3,i)
  !
  info=0
  wt_out(1:4) = 0.0d0
  wt_small_tet(1:4) = 0.0d0
  wt_tmp(1:4) = 0.0d0
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
     if (debug) then
        write(fout,*) 'WARNING:four vertices nearly in the same plane'
        write(fout,*) 'ical = ',ical
        write(fout,*) 'vol = ',vol
        write(fout,*) 'nodes'
        do i=1,4
           write(fout,*) corners(1:3,i)
        enddo
        write(fout,*) 'vecs'
        do i=1,3
           write(fout,*) vec(i,1:3)
        enddo
     endif
     wt_out(1:4)=0.0d0
     !    If the frequency is zero, the contribution from the Fermi surface has to be calculated
     if ( sgnfrq.eq.4 ) then
        call binary_insertion_sort_indx(ind, e, 4)
        do i=1,4
           etmp(i) = e(ind(i))
        enddo
        call ksurf(wt_small_tet, etmp, ef)
        do i=1,4
           wt_out(ind(i)) = wt_small_tet(i)*(-pi)
        enddo
     endif
     return 
  endif
  ! For frequency dependent weights, calculate the energy diferences at the corners of the small tetrahedron and store the maximum absolute value
  if (sgnfrq.ne.1) then
     max_de_small=0.0d0
     do i=1,4
        delta_e_small_tet(i) = (delta_e_big_tet(2)-delta_e_big_tet(1)) * corners(1,i) + (delta_e_big_tet(3) - delta_e_big_tet(1)) * corners(2,i) + (delta_e_big_tet(4) - delta_e_big_tet(1)) * corners(3,i) + delta_e_big_tet(1)
        if ( abs(delta_e_small_tet(i) ) > max_de_small )  max_de_small = abs(delta_e_small_tet(i))
     enddo
  endif
  select case (sgnfrq)
  case(1)
     wt_small_tet(1:4) = 1.0d0/2.40d+1
  case(2)
     ! real axis frequency polarization
     if (iop_bzintq.eq.-1) then
        call stweight_numeric (0, delta_e_small_tet, omgga, wt_small_tet, x_i, w_i, n_i) ! iop_omg==0=>real-axis
     else 
        if (omgga > tol_taylor*max_de_small) then
           call stweight_rtaylor(delta_e_small_tet, omgga, wt_small_tet)
        else
           call sorteq(delta_e_small_tet, ind, sigeq)
           call stweight_real(delta_e_small_tet, omgga, sigeq, wt_tmp)
           do i=1,4
              wt_small_tet(ind(i)) = wt_tmp(i)
           enddo
        endif
     endif
  case(3)
     ! imaginary axis polarization
     if(iop_bzintq.eq.-1) then 
        call stweight_numeric(1, delta_e_small_tet, omgga, wt_small_tet, x_i, w_i, n_i)  ! iop_omg==1 => imaginary-axis
        if(PRINT) write(fout,'(A,4g16.6)') '#wt=', wt_small_tet 
     else 
        if (omgga > tol_taylor*max_de_small) then
           call stweight_itaylor(delta_e_small_tet, omgga, wt_small_tet)
        else    
           call sorteq(delta_e_small_tet, ind, sigeq)
           call stweight_imag(delta_e_small_tet, omgga, sigeq, wt_tmp, vol_small_tetra)
           do i=1,4
              wt_small_tet(ind(i)) = wt_tmp(i)
           enddo
        endif
     endif
  case(4)
     call binary_insertion_sort_indx(ind, delta_e_small_tet, 4)
     do i=1,4
        delta_e_small_tet_c(i) = delta_e_small_tet(ind(i))
     enddo
     call ksurf(wt_tmp, delta_e_small_tet_c, omgga)
     do i=1,4
        wt_small_tet(ind(i)) = wt_tmp(i)*(-pi)
     enddo
  end select
  wt_out(1) = sum(wt_small_tet(1:4))
  do i=1,3
     do j=1,4
        wt_out(i+1) = wt_out(i+1) + wt_small_tet(j)*corners(i,j)
     enddo
     wt_out(1) = wt_out(1) - wt_out(i+1)
  enddo
  do i=1,4
     wt_out(i) = wt_out(i)*vol
  enddo
  return
end subroutine generictetra

subroutine stweight_imag(deltae_vert,omeg,equiv_flag,weight_vert, vol_small_tetra)
  ! This subroutine calculates the weight on the whole small tetrahedron
  ! in which the bands at momentum k are fully occupied and (k-q) states are fully unoccupied.
  ! This is for the "sigfreq=3" (weights for the Polarization with imaginary frequencies)
  ! 
  !  1/(i*w-eps) = -2*eps/(w^2+eps^2)
  !
  implicit none
  real(8), intent(out) :: weight_vert(4)! the weight on the whole tetrahedron.
  real(8), intent(in) :: deltae_vert(4) ! difference of the energy in k-mesh tetrahedron vertices and k-q mesh tetrahedron vertices.
  real(8), intent(in) :: omeg           ! the frequency omega to be calculated
  integer, intent(in) :: equiv_flag     ! == 4, none is equal 
                                        ! == 6, deltae_vert(1)=deltae_vert(2).
                                        ! == 8, deltae_vert(1)=deltae_vert(2) and deltae_vert(3)=deltae_vert(4).
                                        ! ==10, deltae_vert(1)=deltae_vert(2)=deltae_vert(3).
                                        ! ==16, deltae_vert(1)=deltae_vert(2)=deltae_vert(3)=deltae_vert(4). 
  real(8), intent(in) :: vol_small_tetra
  !
  real(8), parameter :: tol_unphys_weight = 1.0e+4   ! the tolerance for unphysical weights
  integer :: j,k
  integer :: ivert
  real(8) :: aa, bb, cc, dd, bb1, bb3,bb4
  real(8) :: vp 
  real(8) :: ev(4), tans(4), logs(4)
  real(8) :: vdif(4,4)
  logical:: lwarn =.false.
  intrinsic datan
  intrinsic dlog
  intrinsic dsign
  !
  weight_vert(1:4) = 0.0d0
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
        call sub_set_vdif_tans_logs(vdif, tans, logs, ev)
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
        weight_vert(ivert) = (aa+bb)/cc
        call check_weight(ivert, weight_vert)
     enddo
  case(6)                ! for the case when ev(1)=ev(1)
     do ivert=1,2
        ev(1) = deltae_vert(ivert)
        ev(2) = ev(1)
        ev(3:4) = deltae_vert(3:4)
        call sub_set_vdif_tans_logs(vdif, tans, logs, ev)
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
        weight_vert(ivert) = (aa+bb)/cc
        call check_weight(ivert, weight_vert)
     enddo ! ivert  
     do ivert=3,4
        ev(1) = (deltae_vert(1)+deltae_vert(2))*0.5d0
        ev(2) = ev(1)
        do j=3,4
           k=mod(j+ivert,2)+3
           ev(k)=deltae_vert(j)
        enddo
        call sub_set_vdif_tans_logs(vdif, tans, logs, ev)
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
        weight_vert(ivert) = (aa+bb)/cc
        call check_weight(ivert, weight_vert)
     enddo ! ivert  
  case(8)          !for the case when ev(1)=ev(1) and ev(3)=ev(4)
     ev(1) = (deltae_vert(1)+deltae_vert(2))/2.d0
     ev(2) = ev(1)
     ev(3) = (deltae_vert(3)+deltae_vert(4))/2.d0
     ev(4) = ev(3) 
     call sub_set_vdif_tans_logs(vdif, tans, logs, ev)
     dd = 6.0d0*omeg**2+ev(1)**2-5.0d0*ev(1)*ev(3)-2.0d0*ev(3)**2
     aa = vdif(3,1)*dd
     dd = 6.0d0*omeg*(omeg**2 - ev(3)**2 - 2.0d0*ev(1)*ev(3))
     aa = aa + dd*(tans(1)-tans(3))
     bb = 3.0d0*((ev(1)+2.0d0*ev(3))*omeg**2 - ev(1)*ev(3)**2 )
     bb = bb*(logs(1)-logs(3))
     cc = 6.0d0*vdif(1,3)**4
     ivert=1
     weight_vert(ivert) = (aa+bb)/cc
     call check_weight(ivert, weight_vert)
     weight_vert(2) = weight_vert(1)
     dd = 6.0d0*omeg**2-2.0d0*ev(1)**2-5.0d0*ev(1)*ev(3)+ev(3)**2
     aa = vdif(1,3)*dd
     dd = 6.0d0*omeg*(omeg**2 - ev(1)**2 - 2.0d0*ev(1)*ev(3))
     aa = aa + dd*(tans(3)-tans(1))
     dd = 3.0d0*( (2.0d0*ev(1) + ev(3))*omeg**2 - ev(1)**2*ev(3) )
     bb = dd*(logs(3)-logs(1))
     cc = 6.0d0*vdif(1,3)**4
     ivert=3
     weight_vert(ivert) = (aa+bb)/cc
     call check_weight(ivert, weight_vert)
     weight_vert(4) = weight_vert(3)
  case(10)             ! for the case when ev(1)=ev(1)=ev(3)
     ev(1:3) = sum(deltae_vert(1:3))/3.0d0
     ev(4) = deltae_vert(4) 
     call sub_set_vdif_tans_logs(vdif, tans, logs, ev)
     aa = vdif(1,4)*(6.0d0*omeg**2-2.0d0*ev(1)**2+7.0d0*ev(1)*ev(4)- 11.0d0*ev(4)**2)     
     dd = 6.0d0*omeg*(omeg**2-3.0d0*ev(4)**2)
     aa = aa - dd * (tans(1)-tans(4))
     dd = 3.0d0*ev(4)*(ev(4)**2-3.0d0*omeg**2)
     bb = dd*(logs(1)-logs(4))
     cc = 18.0d0*vdif(1,4)**4
     ivert=1
     weight_vert(ivert) = (aa+bb)/cc
     call check_weight(ivert, weight_vert)
     weight_vert(2:3) = weight_vert(1) 
     dd = 6.0d0*omeg**2+ev(1)**2-5.0d0*ev(1)*ev(4)-2.0d0*ev(4)**2
     aa = vdif(4,1)*dd
     dd = 6.0d0*omeg*(omeg**2-ev(4)**2-2.0d0*ev(1)*ev(4))
     aa = aa+dd*(tans(1)-tans(4))
     dd = -3.0d0*ev(1)*ev(4)**2+3.0d0*omeg**2*(ev(1)+2.0d0*ev(4))
     bb = dd*(logs(1)-logs(4))
     cc = 6.0d0*vdif(1,4)**4
     ivert=4
     weight_vert(ivert) = (aa+bb)/cc
     call check_weight(ivert, weight_vert)
  case(16)
     ev(1:4) = sum(deltae_vert(1:4))/4.0
     aa = -ev(1)
     bb = 0.0d0
     cc = 12.0d0*(omeg**2+ev(1)**2)
     ivert=1
     weight_vert(ivert) = (aa+bb)/cc !  -ev/(12*(om**2+ev**2)
     call check_weight(ivert, weight_vert)
     weight_vert(2:4)=weight_vert(1)
  case default
     write(6,*)'ERROR in stweight_imat'
     write(6,*)'wrong equiv_flag: ',equiv_flag
     stop "ERROR in stweight_imat"
  end select

contains
  subroutine check_weight(ivt, wght_vert)
    integer, intent(in)    :: ivt
    real(8), intent(inout) :: wght_vert(4)
    !weight_vert(ivert) = (aa+bb)/cc
    if ( abs(wght_vert(ivt)) > tol_unphys_weight ) then
       if ( abs(omeg) < 1.e-6 ) then 
          wght_vert = 0.0 
       else
          write(6,'(A)') 'WARNING in stweight_imag: unexpected big weight!'
          write(6,'(A,I4,A,g16.6,A,I4)') 'vertix:', ivt,' weightt =', wght_vert(ivt), ' case:', equiv_flag
          write(6,'(A,g16.6)') 'omeg =', omeg
          write(6,'(A,4g16.6)') 'deltae_vert =', deltae_vert(:)
          write(6,'(A,4(4g16.6,/))') 'vdif = ', vdif
          write(6,'(A,g16.6,A,g16.6,A,g16.6)') ' aa =', aa,' bb =', bb, ' cc =', cc
          write(6,'(A,g16.6)') ' vol_small_tetra=', vol_small_tetra
       endif
    endif
  end subroutine check_weight
  subroutine sub_set_vdif_tans_logs(vdif, tans, logs, ev)
    ! vdif(:,:) = ev(:) - ev(:)
    ! tans(:) = atan( ev(:) / omeg )
    ! logs(:) = log( ev(:)**2 + omeg**2 )
    implicit none
    real(8), intent(out) :: vdif(4,4)
    real(8), intent(out) :: tans(4), logs(4)
    real(8), intent(in)  :: ev(4)
    integer(4) :: ii, jj
    do ii=1,4
       do jj=1,4
          vdif(ii,jj)=ev(ii)-ev(jj)
       enddo
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
  end subroutine sub_set_vdif_tans_logs
end subroutine stweight_imag

subroutine stweight_itaylor(deltae_vert, omeg, weight_vert)
  ! This subroutine calculates the weight on the whole small tetrahedron
  ! in which the $k$ states are fully occupied and $k-q$ states are fully 
  ! unoccupied. This is for the $sigfreq=3$ (weights for the Polarization with imaginary frequencies)
  implicit none
  real(8), intent(in) :: deltae_vert(4)   ! difference of the energy in k-mesh tetrahedron vertices  and k-q mesh tetrahedron vertices.
  real(8), intent(in) :: omeg             ! the frequency omega to be calculated
  real(8), intent(out):: weight_vert(4)   ! the weight on the whole tetrahedron.
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
end subroutine stweight_itaylor


subroutine stweight_real(deltae_vert,freq,equiv_flag,weight_vert)
  ! This subroutine calculates the weight on the second vertex of the tetrahedron 
  ! which is fully occupied for k and fully unoccupied for k-q. It is also for the 
  ! case of "sigfreq=2" ( weights for the Polarization with real frequencies)
  ! 
  ! weight_vert[i] = 1/24 * ( 1/(freq-delta_vert[i]) + 1/(-freq-delta_vert[i]) )
  implicit none
  real(8), intent(out):: weight_vert(4)  ! the weight on the whole tetrahedron.
  real(8), intent(in) :: deltae_vert(4)  ! energy differences at four vertices 
  real(8), intent(in) :: freq            ! the frequency omega to be calculated
  integer, intent(in) :: equiv_flag      ! equiv_flag = 4, none is equal 
  ! equiv_flag=6, v(1)=v(2).
  ! equiv_flag=8, v(1)=v(2) and v(3)=v(4).
  ! equiv_flag=10, v(1)=v(2)=v(3).
  ! e1uiv_flag=16, v(1)=v(2)=v(3)=v(4). 
  real(8), parameter :: tol_unphys_weight = 1.0e+4   ! the tolerance for unphysical weights
  ! 
  integer :: i,j, k
  integer :: ivert
  integer :: isign_om
  real(8) :: aa,  cc, dd, vp
  real(8) :: omeg
  real(8) :: vdif(4,4), ev(4), weight_tmp(2,4)
  logical :: lwarn = .False.
  real(8), parameter :: haier = 1.0d-12
  intrinsic dlog
  intrinsic dabs
  !
  weight_tmp=0.0d0
  do isign_om=1,2
     omeg = dble(3-2*isign_om)*freq  ! om = [freq,-freq]
     select case(equiv_flag)
     case(4)                        ! for the case of none of them are equal
        do ivert=1,4
           do j=1,4
              k=mod(j+ivert-2,4)+1       
              ev(k)=deltae_vert(j)
           enddo
           call set_vdif( vdif, ev )
           if ( dabs(omeg-ev(1)) < haier ) then
              aa = -dlog(dabs(vdif(2,1)))*vdif(2,1)*vdif(4,3)
              aa = aa + dlog(dabs(vdif(3,1))) * vdif(3,1) * vdif(4,2)
              aa = aa - dlog(dabs(vdif(4,1))) * vdif(4,1) * vdif(3,2)
              cc = 6.0d0 * vdif(3,2) * vdif(4,3) * vdif(4,2)
           else if (dabs(omeg-ev(2)) < haier) then
              aa = dlog(dabs(vdif(3,2)))*vdif(3,2)**2*vdif(4,1)**2
              aa = aa - dlog(dabs(vdif(4,2)))*vdif(3,1)**2*vdif(4,2)**2
              dd = dlog(dabs(vdif(2,1))) * ( vdif(3,1) * vdif(4,2) + vdif(3,2) * vdif(4,1) )
              dd = dd + vdif(3,1) * vdif(4,1)
              aa = aa + vdif(2,1) * vdif(4,3) * dd
              cc = 6.0d0 * vdif(3,1)**2 * vdif(4,1)**2 * vdif(4,3)
           else if ( dabs(omeg-ev(3)) < haier) then
              aa = dlog(dabs(vdif(3,2))) * vdif(3,2)**2 * vdif(4,1)**2
              aa = aa - dlog(dabs(vdif(4,3)))*vdif(2,1)**2*vdif(4,3)**2
              dd = dlog(dabs(vdif(3,1))) * ( vdif(4,1) * vdif(3,2) - vdif(4,3) * vdif(2,1) )
              dd = dd - vdif(2,1) * vdif(4,1)
              aa = aa - vdif(3,1) * vdif(4,2) * dd
              cc = 6.0d0 * vdif(2,1)**2 * vdif(4,1)**2 * vdif(4,2)
           else if ( dabs(omeg-ev(4)) < haier) then
              aa = 0.0d0 - dlog(dabs(vdif(4,3))) * vdif(3,1)**2 * vdif(4,2)**2
              aa = aa + dlog(dabs(vdif(4,2))) * vdif(4,2)**2 * vdif(3,1)**2
              dd = dlog(dabs(vdif(4,1))) - dlog(dabs(vdif(4,3)))
              dd = dd * ( vdif(3,1) * vdif(4,2) + vdif(2,1) * vdif(4,3) ) - vdif(2,1) * vdif(3,1)
              aa = aa - vdif(3,2) * vdif(4,1) * dd
              cc = 6.0d0 * vdif(2,1)**2 * vdif(3,1)**2 * vdif(3,2)
           else
              dd = (omeg-ev(4))**3 * vdif(2,1)**2 * vdif(3,2) * vdif(3,1)**2
              aa = dlog(dabs(omeg-ev(4))) * dd
              dd = (omeg-ev(3))**3 * vdif(2,1)**2 * vdif(4,2) * vdif(4,1)**2
              aa = aa - dlog(dabs(omeg-ev(3))) * dd
              dd = -vdif(2,1) * vdif(3,1) * (omeg-ev(4)) - (omeg-ev(2)) * vdif(3,1) * vdif(4,1)
              dd = dd - (omeg-ev(3)) * vdif(2,1) * vdif(4,1)
              dd = dd * vdif(3,2) * vdif(4,2) * vdif(4,3) * (omeg-ev(1))**2
              aa = aa + dlog(dabs(omeg-ev(1))) * dd
              aa = aa + (omeg-ev(1))**2 * vdif(2,1) * vdif(3,2) * vdif(3,1) * vdif(4,2) * vdif(4,1) * vdif(4,3)
              aa = aa + dlog(dabs(omeg-ev(2))) * (omeg-ev(2))**3 * vdif(3,1)**2 * vdif(4,1)**2 * vdif(4,3)
              cc = 6.0d0 * vdif(2,1)**2 * vdif(3,2) * vdif(3,1)**2 * vdif(4,2) * vdif(4,1)**2 * vdif(4,3)
           endif
           call sub_set_weight(weight_tmp, isign_om, ivert, aa, cc)
        enddo ! ivert  
     case(6)                     ! for the case that a=b
        do ivert=1,2
           ev(1) = deltae_vert(ivert)
           ev(2) = ev(1)
           ev(3:4) = deltae_vert(3:4)
           call set_vdif( vdif, ev )
           if ( dabs(omeg-ev(1)) < haier) then
              aa = dlog(dabs(vdif(3,1))) - dlog(dabs(vdif(4,1)))
              cc = 6.0d0*vdif(4,3)
           else if ( dabs(omeg-ev(3)) < haier) then
              aa = 2.0d0 * ( dlog(dabs(vdif(3,1))) - dlog(dabs(vdif(4,1))) ) * vdif(4,1)**2
              aa = aa + ( 2.0d0 * vdif(4,3) + vdif(4,1) ) * vdif(4,1)
              cc = 12.0d0 * vdif(4,1)**3
           else if ( dabs(omeg-ev(4)) < haier) then
              aa = vdif(3,1) * (vdif(3,1) - 2.0d0 * vdif(4,3))
              aa = aa - 2.0d0 * (dlog(dabs(vdif(4,3)))-dlog(dabs(vdif(4,1)))) * vdif(4,3)**2
              cc = 12.0d0*vdif(3,1)**3
           else
              aa = 2.0d0 * dlog(dabs(omeg-ev(4))) * (omeg-ev(4))**3 * vdif(3,1)**3
              aa = aa - 2.0d0 * dlog(dabs(omeg-ev(3))) * (omeg-ev(3))**3 * vdif(4,1)**3
              dd = (omeg-ev(4))**2 * vdif(3,1)**2 + (omeg-ev(3)) * vdif(3,1) * (omeg-ev(4)) * vdif(4,1)
              dd = dd + (omeg-ev(3))**2 * vdif(4,1)**2
              aa = aa + 2.0d0 * dlog(dabs(omeg-ev(1))) * (omeg-ev(1)) * dd * vdif(4,3)
              dd = vdif(3,1) * vdif(4,1) - 2.0d0 * vdif(3,1) * (omeg-ev(4)) - 2.0d0 * (omeg-ev(3)) * vdif(4,1)
              aa = aa + (omeg-ev(1)) * dd * vdif(3,1) * vdif(4,1) * vdif(4,3)
              cc = 12.0d0 * vdif(3,1)**3 * vdif(4,1)**3 * vdif(4,3)
           endif
           call sub_set_weight(weight_tmp, isign_om, ivert, aa, cc)
        enddo ! ivert  
        do ivert=3,4
           ev(1) = (deltae_vert(1)+deltae_vert(2))*0.5d0
           ev(2) = ev(1)
           do j=3,4
              k=mod(j+ivert,2)+3
              ev(k)=deltae_vert(j)
           enddo
           call set_vdif( vdif, ev )
           if ( dabs(omeg-ev(1)) < haier) then
              aa = ( dlog(dabs(vdif(3,1))) - dlog(dabs(vdif(4,1))) ) * vdif(4,1)
              aa = aa + vdif(4,3)
              cc = 6.0d0*vdif(4,3)**2
           else if (dabs(omeg-ev(3)).lt.haier) then
              aa = (dlog(dabs(vdif(3,1)))-dlog(dabs(vdif(4,3))))*vdif(4,3)
              aa = aa + vdif(4,1)
              cc = 6.0d0 * vdif(4,1)**2
           else if ( dabs(omeg-ev(4)) < haier) then
              aa = vdif(3,1) * (vdif(4,3) + vdif(4,1))
              aa = aa - 2.0d0*(dlog(dabs(vdif(4,1)))-dlog(dabs(vdif(4,3)))) * vdif(4,1)*vdif(4,3)
              cc = 6.0d0 * vdif(3,1)**3
           else
              aa = dlog(dabs(omeg-ev(4))) * (omeg-ev(4))**3 * vdif(3,1)**3
              dd = 2.0d0 * vdif(4,3) * (omeg-ev(1)) - (omeg-ev(4)) * vdif(3,1)
              aa = aa + dlog(dabs(omeg-ev(3))) * (omeg-ev(3))**2 * dd * vdif(4,1)**2
              dd = (omeg-ev(3))**2 * vdif(4,1) + (omeg-ev(1))**2 * vdif(4,3)
              aa = aa + dd*vdif(3,1) * vdif(4,1) * vdif(4,3)
              dd = 2.0d0*(omeg-ev(3)) * vdif(4,1) + vdif(3,1) * (omeg-ev(4))
              aa = aa - dlog(dabs(omeg-ev(1))) * (omeg-ev(1))**2 * dd * vdif(4,3)**2
              cc = 6.0d0 * vdif(3,1)**3 * vdif(4,1)**2 * vdif(4,3)**2
           endif
           call sub_set_weight(weight_tmp, isign_om, ivert, aa, cc)
        enddo ! ivert  
     case(8)                      ! for the case that a=b and c=d
        do ivert=1,2
           ev(1) = deltae_vert(ivert)
           ev(2) = ev(1)
           vp = (deltae_vert(3) + deltae_vert(4)) * 0.5d0
           ev(3) = vp
           ev(4) = vp
           call set_vdif( vdif, ev )
           if ( dabs(omeg-ev(1)) < haier ) then
              aa = -1.0d0
              cc = 6.0d0*vdif(3,1)
           else if ( dabs(omeg-ev(3)) < haier) then
              aa = 1.0d0
              cc = 12.0d0*vdif(3,1)
           else
              dd = dlog(dabs(omeg-ev(1))) - dlog(dabs(omeg-ev(3)))
              aa = 6.0d0 * (omeg-ev(1)) * dd * (omeg-ev(3))**2
              dd = 9.0d0 * (omeg-ev(1)) * vdif(3,1) - 2.0d0 * vdif(3,1)**2 - 6.0d0 * (omeg-ev(1))**2
              aa = aa + vdif(3,1) * dd
              cc = 12.0d0 * vdif(3,1)**4    
           endif
           call sub_set_weight(weight_tmp, isign_om, ivert, aa, cc)
        enddo ! ivert  
        do ivert=3,4
           ev(3) = deltae_vert(ivert)
           ev(4) = ev(3)
           vp = (deltae_vert(1)+deltae_vert(2))*0.5d0
           ev(1) = vp
           ev(2) = vp
           call set_vdif( vdif, ev )
           if ( dabs(omeg-ev(1)) < haier ) then
              aa = -1.0d0
              cc = 12.0d0 * vdif(3,1)
           else if ( dabs(omeg-ev(3)) < haier) then
              aa = 1.0d0
              cc = 6.0d0 * vdif(3,1)
           else
              dd = dlog(dabs(omeg-ev(3))) - dlog(dabs(omeg-ev(1)))
              aa = 6.0d0 * (omeg-ev(1))**2 * dd * (omeg-ev(3))
              dd = 3.0d0 * (omeg-ev(1)) * vdif(3,1) + vdif(3,1)**2 - 6.0d0 * (omeg-ev(1))**2
              aa = aa - vdif(3,1) * dd
              cc = 12.0d0 * vdif(3,1)**4    
           endif
           call sub_set_weight(weight_tmp, isign_om, ivert, aa, cc)
        enddo ! ivert  
     case(10)                   ! for the case that a=b=c
        do ivert=1,3
           ev(1)=deltae_vert(ivert)
           ev(2)=ev(1)
           ev(3)=ev(1)
           ev(4)=deltae_vert(4)
           call set_vdif( vdif, ev )
           if ( dabs(omeg-ev(4)) < haier) then
              aa = 1.0d0
              cc = 18.0d0*vdif(4,1)
           else
              aa = 6.0d0 * (dlog(dabs(omeg-ev(4))) - dlog(dabs(omeg-ev(1)))) * (omeg-ev(4))**3
              dd = 0.0d0 - 15.0d0 * (omeg-ev(1)) * vdif(4,1) + 11.0d0*vdif(4,1)**2 + 6.0d0*(omeg-ev(1))**2
              aa = aa + vdif(4,1) * dd
              cc = 36.0d0 * vdif(4,1)**4
           endif
           call sub_set_weight(weight_tmp, isign_om, ivert, aa, cc)
        enddo ! ivert  
        vp = sum(deltae_vert(1:3))/3.0d0
        ev(1:3) = vp
        ev(4) = deltae_vert(4) 
        call set_vdif( vdif, ev )
        ivert = 4
        if(dabs(omeg-ev(4)).lt.haier) then
           aa=1.0d0
           cc=12.0d0*vdif(4,1)
        else
           dd = 5.0d0 * ev(1) * ev(4) + 2.0d0 * ev(4)**2 - 3.0d0*ev(1) * omeg - ev(1)**2
           dd = dd - 9.0d0 * ev(4) * omeg + 6.0d0*omeg**2
           aa = vdif(1,4) * dd
           dd = dlog(dabs(omeg-ev(1))) - dlog(dabs(omeg-ev(4)))
           aa = aa - 6.0d0 * (ev(1)-omeg) * (ev(4)-omeg)**2 * dd
           cc = 12.0d0*vdif(4,1)**4
        endif
        call sub_set_weight(weight_tmp, isign_om, ivert, aa, cc)
     case(16)                   ! for the case that a=b=c=d
        ev(1:4) = deltae_vert(1:4)
        aa = 1.0d0  ! weight_vert[i] = 1/24 * ( 1/(omeg-delta_vert[i]) + 1/(-omeg-delta_vert[i]) )
        do ivert=1,4    
           cc = 24.0d0 * (omeg-ev(ivert))
           call sub_set_weight(weight_tmp, isign_om, ivert, aa, cc)
        enddo ! ivert
     end select
  enddo ! isign_om
  do ivert=1,4
     weight_vert(ivert) = weight_tmp(1,ivert) + weight_tmp(2,ivert)
  enddo

contains

  subroutine sub_set_weight(weight_tmp, isign_om, ivert, aa, cc)
    real(8), intent(inout) :: weight_tmp(2,4)
    integer, intent(in)    :: isign_om, ivert
    real(8), intent(in)    :: aa, cc
    !
    if ( ( dabs(aa) < dabs(cc) * tol_unphys_weight ) )  then  ! aa/cc < tol_unphys_weight
       weight_tmp(isign_om,ivert) = aa/cc
    else
       if ( abs(cc) < 1.e-20 ) weight_tmp(isign_om,ivert) = aa/(1e-20)
       write(6,'(A)')'warning in stweight_real: large weight !!!'
       write(6,'(A,I4,A,g18.10,A,I4)') 'vertix:', ivert,' weightt =', weight_tmp(isign_om,ivert), ' case:',equiv_flag
       write(6,'(A,g18.10,g18.10)') 'omeg =', omeg
       write(6,'(A,4g18.10)') 'deltae_vert =', deltae_vert
       write(6,'(A,4(4g18.10,/))') 'vdif = ', vdif
       write(6,'(A,g18.10,A,g18.10)') ' aa =', aa, ' cc =', cc
    endif
  end subroutine sub_set_weight
  !
  subroutine set_vdif( vdif, ev )
    implicit none
    real(8), intent(out):: vdif(4,4)
    real(8), intent(in) :: ev(4)
    !
    integer :: i, j
    do i=1,4
       do j=1,4
          vdif(i,j) = ev(i)-ev(j)
       enddo
    enddo
  end subroutine set_vdif
end subroutine stweight_real
subroutine stweight_rtaylor(deltae_vert, freq, weight_vert)
  ! This subroutine calculates the weight on the whole small tetrahedron
  ! in which the $k$ states are fully occupied and $k-q$ states are fully 
  ! unoccupied. This is for the "sigfreq=2" ( weights for the Polarization with real frequencies)
  implicit none
  real(8), intent(in) :: deltae_vert(4) ! difference of the energy in k-mesh tetrahedron vertices and k-q mesh tetrahedron vertices.
  real(8), intent(in) :: freq           ! the frequency omega to be calculated
  real(8), intent(out):: weight_vert(4) ! the weight on the whole tetrahedron.
  !
  integer :: ivert,j,k,isign_om
  real(8) :: denom0, denom1, denom2, denom3, denom4
  real(8) :: w00, w01, w02, w03
  real(8) :: omeg
  real(8) :: ev(4)
  real(8) :: weight_tmp(2,4)
  do isign_om = 1,2
     omeg = dble(3-2*isign_om)*freq  ! [freq,-freq]
     denom0 = omeg
     denom1 = omeg**2
     denom2 = omeg**3
     denom3 = omeg**4
     denom4 = omeg**5
     do ivert=1,4
        do j=1,4
           k=mod(j+ivert-2,4)+1
           ev(j)=deltae_vert(k)
        enddo
        w00 = 1.0d0/(2.4d+1*denom0)
        w01 = (ev(1)+2.0d0*ev(2)+ev(3)+ev(4))/(1.2d+2*denom1)
        w02 = (ev(1)**2+2.0d0*ev(1)*ev(2)+3.0d0*ev(2)**2+ev(1)*ev(3) + 2.0d0*ev(2)*ev(3)+ev(3)**2+ev(1)*ev(4)+2.0d0*ev(2)*ev(4)+ ev(3)*ev(4)+ev(4)**2)/(3.6d+2*denom2)
        w03 = (ev(1)**3+2.0d0*ev(1)**2*ev(2)+3.0d0*ev(1)*ev(2)**2 + 4.0d0*ev(2)**3+ev(1)**2*ev(3)+2.0d0*ev(1)*ev(2)*ev(3)+ 3.0d0*ev(2)**2*ev(3)+ev(1)*ev(3)**2+2.0d0*ev(2)*ev(3)**2+ ev(3)**3+ev(1)**2*ev(4)+2.0d0*ev(1)*ev(2)*ev(4) + 3.0d0*ev(2)**2*ev(4)+ev(1)*ev(3)*ev(4)+2.0d0*ev(2)*ev(3)*ev(4)+ev(3)**2*ev(4)+ev(1)*ev(4)**2+2.0d0*ev(2)*ev(4)**2 + ev(3)*ev(4)**2+ev(4)**3)/(8.4d+2*denom3)
        weight_tmp(isign_om,ivert) = w00+w01+w02+w03
     enddo ! ivert
  enddo ! isign_om
  do ivert=1,4
     weight_vert(ivert)=weight_tmp(1,ivert)+weight_tmp(2,ivert)
  enddo
  return
end subroutine stweight_rtaylor

subroutine stweight_numeric(iop_omg,de_vert,omg,wt_vert, x_i, w_i, n_i)
  ! This subroutine calculates by numerical integration the weight  on the tetrahedron
  ! with energies de_vert,omg(1:4), and frequency omg. It performs the following integral:
  !
  ! if iop_omg==0:
  !     Real{ Integrate[ fw_i(x,y,z) * 1/( omg - E(x,y,z) + 0.01*1j), {z,0,1},{y,0,1-z},{x,0,1-y-z}] }
  ! or iop_omg==1:
  !    2Real{ Integrate[ fw_i(x,y,z) * 1/( 1j*omg - E(x,y,z) ),       {z,0,1},{y,0,1-z},{x,0,1-y-z}] }
  ! where fw_i(x,y,z) = [(1-x-y-z), x, y, z] corresponds to the four corners of tetrahedra.
  !
  ! de_vert(4) -- energies at the four corners of the tetrahedra
  ! omg        -- frequency
  ! iop_omg==0 : On real axis, iop_omg==1 on imaginary axis.
  !
  implicit none
  real(8), intent(out):: wt_vert(4)   ! the weight on the whole tetrahedron.
  integer, intent(in) :: iop_omg, n_i ! 0/1 - real/imaginary frequency 
  real(8), intent(in) :: de_vert(4)   ! difference of the energy in k-mesh tetrahedron vertices and k-q mesh tetrahedron vertices.
  real(8), intent(in) :: omg          ! the frequency omega to be calculated
  real(8), intent(in) :: x_i(n_i), w_i(n_i) ! x-positions and weights for gauss integration quadrature in 1D
  !
  integer :: j,k,ivert
  real(8) :: d0,d10,d20,d30
  real(8) :: x,y,z                      !! common variables used by different internal subroutines
  real(8), parameter :: eta_freq = 0.01 ! this is the value of the small imaginary part that is needed for real frequency calculations 
  integer, parameter :: nmax = 10
  integer, parameter :: nmin = 3 
  real(8), parameter :: eps = 1.e-4
  logical:: l_conv
  intrinsic datan
  intrinsic dlog
  intrinsic dsign
  !
  d0  = de_vert(1) 
  d10 = de_vert(2)-de_vert(1)
  d20 = de_vert(3)-de_vert(1)
  d30 = de_vert(4)-de_vert(1)
  call gauq_z(1.d0, wt_vert) 
contains
  !     This defines the integrand 
  function func(x,y,z) result(fv)
    ! Note that we use (d0, d10, d20, d30) from main subroutine
    implicit none
    real(8) :: fv(4) 
    real(8), intent(in) :: x,y,z
    !
    real(8):: de, comm 
    de = d0 + x*d10 + y*d20 + z*d30 ! The 3D interpolation of the energy, by knowing its value in the vertices of tetrahedra (d0, d10, d20, d30)
    if (iop_omg.eq.0) then 
       comm = (omg-de)/((omg-de)**2+eta_freq**2)
    else
       comm = -2.d0*de/(omg*omg + de*de)
    endif
    fv(1) = (1.d0-x-y-z)*comm
    fv(2) = x*comm
    fv(3) = y*comm
    fv(4) = z*comm
  end function func
  function f(xx) result(fx) 
    real(8):: xx,fx(4) 
    x = xx
    fx = func(x,y,z)
  end function f
  function g(yy) result(gy) 
    real(8):: yy,gy(4) 
    y = yy
    call gauq_x(1.0-y-z,gy)
  end function g
  function h(zz) result(hz)
    real(8):: zz,hz(4) 
    z = zz 
    call gauq_y(1.0-z,hz)
  end function h
  subroutine gauq_x(a,s)
    real(8),intent(in) :: a
    real(8),intent(out):: s(4) 
    integer :: i
    s=0.d0
    do i=1,n_i
       s = s + w_i(i)*a * f(a*x_i(i))
    enddo
  end subroutine gauq_x
  subroutine gauq_y(a,s)
    real(8),intent(in) :: a
    real(8),intent(out):: s(4)
    integer :: i
    s=0.d0
    do i=1,n_i
       s = s + w_i(i)*a * g(a*x_i(i))
    enddo
  end subroutine gauq_y
  subroutine gauq_z(a,s)
    real(8),intent(in) :: a
    real(8),intent(out):: s(4)
    integer :: i
    s=0.d0
    do i=1,n_i
       s = s + w_i(i)*a* h(a*x_i(i))
    enddo
  end subroutine gauq_z
end subroutine stweight_numeric

subroutine gauleg(x1,x2,x,w,n)
  ! Given the lower and upper limits of integration \texttt{x1} and
  ! \texttt{x2}, and given \texttt{n}, this routine returns arrays
  ! \texttt{x(1:n)} and \texttt{w(1:n)} of length \texttt{n}, containing the
  ! abscissas and weights of the Gauss-Legendre \texttt{n}-point quadrature formula.
  implicit none
  integer(4), intent(in) :: n  ! Order of the quadrature
  real(8),    intent(in) :: x1 ! Lower integration limit
  real(8),    intent(in) :: x2 ! Upper integration limit
  real(8),    intent(out) :: x(n) ! abscissas      
  real(8),    intent(out) :: w(n) ! weights
  !
  integer(4) :: i    
  integer(4) :: j            
  integer(4) :: m
  real(8) :: dj
  real(8) :: p1            
  real(8) :: p2            
  real(8) :: p3            
  real(8) :: pp            
  real(8) :: xl            
  real(8) :: xm            
  real(8) :: z            
  real(8) :: z1            
  logical :: conv
  real(8), parameter :: eps = 3.0d-14
  real(8), parameter :: pi  = 3.141592653589793d+0
  intrinsic abs      
  intrinsic cos      
  intrinsic dble
  !      
  ! Original subroutine: gauleg.for (C) copr. 1986-92 copr. 1986-92 numerical recipes software &145i..  
  !     The roots are symmetric in the interval, so we only have to find
  !     half of them.
  !
  m = ( n + 1 ) / 2
  xm = 0.5d0 * ( x2 + x1 )
  xl = 0.5d0 * ( x2 - x1 )
  do i = 1, m ! Loop over the desired roots
     z = cos(pi*(dble(i)-0.25d0)/(dble(n)+0.5d0))
     !       Starting with the above approximation to the \texttt{ith} root, we
     !       enter the main loop of refinement by Newton's method 
     conv = .false.
     do while (.not.conv)
        p1 = 1.0d0
        p2 = 0.0d0
        !         Loop up the recurrence relation to get the Legendre polynomial
        !         evaluated at z 
        do j = 1, n
           dj=dble(j)
           p3 = p2
           p2 = p1
           p1 = ((2.0d0*dj-1.0d0)*z*p2 - (dj-1.0d0)*p3)/dj
        enddo ! j  
        !         p1 is now the desired Legendre polynomial. We next compute pp, its
        !         derivative, by a standard relation involving also p2, the
        !         polynomial of one lower order.          
        pp = dble(n) * (z * p1 - p2) / (z * z - 1.0d0) 
        !         Newton's method.
        z1 = z
        z = z1 - p1 / pp
        conv = abs(z-z1).le.eps
     enddo ! conv  
     !       Scale the root to the desired interval, and put in its symmetric
     !       counterpart.
     x(i) = xm - xl * z
     x(n+1-i) = xm + xl * z
     !        Compute the weight and its symmetric counterpart.
     w(i) = 2.0d0 * xl / ((1.0d0 - z * z) * pp * pp)
     w(n+1-i) = w(i)
  enddo ! i  
  !     Successful exit
  if (.False.) then
     write(6,'(A,F16.10,1x,A,F16.10,1x,A,I2)') 'gauleg: x1=', x1, 'x2=', x2, 'n=', n
     do i=1,n
        write(6,'(A,I2,A,F16.10,1x,A,I2,A,F16.10)') 'x[',i,']=', x(i), 'w[',i,']=', w(i)
     enddo
  endif
  return
end subroutine gauleg


subroutine sorteq(v,index,sigeq)
  ! This subroutine sorts the values of the vector v (integer) in an
  ! order that the first one is the one with the biggest possibility of 
  ! being equal to the others, then the second one, third one and the 
  ! forth one. Optionally, it also returns an integer vector that gives 
  ! the old indexes in vector v according to the new ones
  implicit none
  real(8), intent(inout) :: v(4)
  integer, intent(out)   :: index(4)  ! this index tells us the relation between the sorted and unsorted arrays.
  integer, intent(out)   :: sigeq
  !
  real(8), parameter :: ztol_sorteq = 1.e-2
  integer(4) :: i,j,itmp
  integer(4) :: sig(5)   ! sig(1:4) number of items in the array equal to this specific item. For example, if all are different sig(1:4)=1. sig(5) is the sum 
  real(8)    :: vtmp, dist
  intrinsic abs
  sig=0 
  do i=1,4
     do j=1,4
        if (f_rdist(v(i),v(j)) < ztol_sorteq) sig(i)=sig(i)+1
     enddo
  enddo
  sig(5) = sig(1) + sig(2) + sig(3) + sig(4)
  if ( (sig(5).eq.8) .and. (maxval(sig(1:4)).eq.3) ) then
     do i=1,4
        if (sig(i).eq.2) sig(i) = sig(i) + 1
     enddo
     sig(5) = sig(1) + sig(2) + sig(3) + sig(4)
  endif
  if ( (sig(5).eq.10) .and. (minval(sig(1:4)).eq.2) ) then
     do i=1,4
        sig(i)=4
     enddo
     sig(5) = sig(1) + sig(2) + sig(3) + sig(4)
  endif
  if (sig(5) .eq. 12) sig(5) = 16 
  if (sig(5) .eq. 14) sig(5) = 16 
  do i=1,4
     index(i)=i
  enddo
  select case (sig(5))
  case (16)
     continue
  case (10)                  ! three of them are equal
     do i=1,3                ! we want to make it as a=b=c while d not
        if(sig(i).eq.1) then
           vtmp=v(i)         ! we make sig(4)=1
           v(i)=v(4)
           v(4)=vtmp
           itmp=index(i)
           index(i)=index(4)
           index(4)=itmp
           itmp=sig(4)
           sig(4)=sig(i)
           sig(i)=itmp
        endif
     enddo
  case (8)                   ! doubly equal but not all
     do i=3,4                ! make the first two equal
        if( f_rdist(v(i),v(1)).lt.f_rdist(v(2),v(1)) )then
           vtmp=v(i)
           v(i)=v(2)
           v(2)=vtmp
           itmp=index(i)
           index(i)=index(2)
           index(2)=itmp
           itmp=sig(i)
           sig(i)=sig(2)
           sig(2)=itmp
        endif
     enddo
  case (6)                   ! only two of them are equal
     j=1
     do i=1,4                         ! make the first one with sig(1)=2
        if((sig(i).eq.2)) then
           vtmp=v(j)
           v(j)=v(i)
           v(i)=vtmp
           itmp=index(j)
           index(j)=index(i)
           index(i)=itmp
           itmp=sig(j)
           sig(j)=sig(i)
           sig(i)=itmp
           j=j+1
        endif
     enddo
  case (4)      ! all different, nothing to do
     continue
  case default   ! None of the above... ERROR
     write(6,*) 'ERROR in sorteq: case not found'
     write(6,*) 'sig(i)=', sig
     write(6,'(a,4g16.6)') 'v = ',v
     stop "ERROR in sorteq"
!                    sig(5)=4 and sig(5) corresponds to the all not equal 
!                    and the all equal case respectively
  end select
  sigeq = sig(5)
contains 
  !     this internal function measures the relative distance
  real(8) function f_rdist(v1,v2)
    implicit none
    real(8),intent(in) :: v1,v2
    if ((abs(v1)+abs(v2)) < 1.0e-10) then 
       f_rdist = 0.d0 
    else 
       f_rdist = 2.d0*abs(v1-v2)/(abs(v1)+abs(v2))
    endif
    return
  end function f_rdist
end subroutine sorteq
 



subroutine surfnodes(iwrite,   ntype, intnodes, nnod, e, f, ef)
  ! This subroutine calculates the coordinates of the intersections of three
  ! planes between the four planes that form the surface of the tetrahedron
  ! and the two planes that approximate the fermi surface within it for the
  ! bands at $\vec{k}$ and $\vec{k}-\vec{q}$, in internal coordinates of the
  ! tetrahedron $\vec{x}=(\xi,\eta,\zeta)$. The equations defining the planes are:
  !
  !\begin{subequations}
  !\begin{align}
  !\xi=&0\\
  !\eta=&0\\
  !\zeta=&0\\
  !\xi+\eta+\zeta=&1\\
  !(\varepsilon_2-\varepsilon_1)\xi+(\varepsilon_3-\varepsilon_1)\eta+%
  !(\varepsilon_4-\varepsilon_1)\zeta=&\varepsilon_F-\varepsilon_1\\
  !(\varepsilon'_2-\varepsilon'_1)\xi+(\varepsilon'_3-\varepsilon'_1)\eta+%
  !(\varepsilon'_4-\varepsilon'_1)\zeta=&\varepsilon_F-\varepsilon'_1
  !\end{align}
  !\end{subequations}
  !
  ! where $\varepsilon_i=\varepsilon_n(\vec{k}_i)$ and 
  ! $\varepsilon'_i=\varepsilon_n(\vec{k}_i-\vec{q})$.
  !
  ! The expressions for the nodes, in the order they are calculated, which
  !corresponds to increasing ntype ordering, are:
  !
  !\begin{subequations}
  !\begin{align}
  ! \vec{x}_1=&(0,0,0)\\
  ! \vec{x}_2=&(0,0,1)\\
  ! \vec{x}_3=&(0,1,0)\\
  ! \vec{x}_4=&(1,0,0)\\
  ! \vec{x}_5=&(0,0,\tfrac{\varepsilon_F-\varepsilon_1}%
  ! {\varepsilon_4-\varepsilon_1})\\
  ! \vec{x}_6=&(0,\tfrac{\varepsilon_F-\varepsilon_1}%
  ! {\varepsilon_3-\varepsilon_1},0)\\
  ! \vec{x}_7=&(\tfrac{\varepsilon_F-\varepsilon_1}%
  ! {\varepsilon_2-\varepsilon_1},0,0)\\
  ! \vec{x}_8=&(0,\tfrac{\varepsilon_4-\varepsilon_F}%
  ! {\varepsilon_4-\varepsilon_3},\tfrac{\varepsilon_3-\varepsilon_F}%
  ! {\varepsilon_3-\varepsilon_4})\\
  ! \vec{x}_9=&(\tfrac{\varepsilon_4-\varepsilon_F}%
  ! {\varepsilon_4-\varepsilon_2},0,\tfrac{\varepsilon_2-\varepsilon_F}%
  ! {\varepsilon_2-\varepsilon_4})\\
  ! \vec{x}_{10}=&(\tfrac{\varepsilon_3-\varepsilon_F}%
  ! {\varepsilon_3-\varepsilon_2},\tfrac{\varepsilon_2-\varepsilon_F}%
  ! {\varepsilon_2-\varepsilon_3},0)\\
  ! \vec{x}_{11}=&(0,0,\tfrac{\varepsilon_F-\varepsilon'_1}%
  ! {\varepsilon'_4-\varepsilon'_1})\\
  ! \vec{x}_{12}=&(0,\tfrac{\varepsilon_F-\varepsilon'_1}%
  ! {\varepsilon'_3-\varepsilon'_1},0)\\
  ! \vec{x}_{13}=&(\tfrac{\varepsilon_F-\varepsilon'_1}%
  ! {\varepsilon'_2-\varepsilon'_1},0,0)\\
  ! \vec{x}_{14}=&(0,\tfrac{\varepsilon'_4-\varepsilon_F}%
  ! {\varepsilon'_4-\varepsilon'_3},\tfrac{\varepsilon'_3-\varepsilon_F}%
  ! {\varepsilon'_3-\varepsilon'_4})\\
  ! \vec{x}_{15}=&(\tfrac{\varepsilon'_4-\varepsilon_F}%
  ! {\varepsilon'_4-\varepsilon'_2},0,\tfrac{\varepsilon'_2-\varepsilon_F}%
  ! {\varepsilon'_2-\varepsilon'_4})\\
  ! \vec{x}_{16}=&(\tfrac{\varepsilon'_3-\varepsilon_F}%
  ! {\varepsilon'_3-\varepsilon'_2},\tfrac{\varepsilon'_2-\varepsilon_F}%
  ! {\varepsilon'_2-\varepsilon'_3},0)\\
  ! \vec{x}_{17}=&(0,%
  !\tfrac{(\varepsilon_F-\varepsilon_1)(\varepsilon_F-\varepsilon'_4)-%
  !(\varepsilon_F-\varepsilon_4)(\varepsilon_F-\varepsilon'_1)}%
  !{(\varepsilon_4-\varepsilon_3)(\varepsilon'_4-\varepsilon'_1)-%
  !(\varepsilon_4-\varepsilon_1)(\varepsilon'_4-\varepsilon'_3)},
  !\tfrac{(\varepsilon_F-\varepsilon_3)(\varepsilon_F-\varepsilon'_1)-%
  !(\varepsilon_F-\varepsilon_1)(\varepsilon_F-\varepsilon'_3)}%
  !{(\varepsilon_4-\varepsilon_3)(\varepsilon'_4-\varepsilon'_1)-%
  !(\varepsilon_4-\varepsilon_1)(\varepsilon'_4-\varepsilon'_3)})\\
  ! \vec{x}_{18}=&(%
  !\tfrac{(\varepsilon_F-\varepsilon_4)(\varepsilon_F-\varepsilon'_1)-%
  !(\varepsilon_F-\varepsilon_1)(\varepsilon_F-\varepsilon'_4)}%
  !{(\varepsilon_2-\varepsilon_4)(\varepsilon'_2-\varepsilon'_1)-%
  !(\varepsilon_2-\varepsilon_1)(\varepsilon'_2-\varepsilon'_4)},0,%
  !\tfrac{(\varepsilon_F-\varepsilon_1)(\varepsilon_F-\varepsilon'_2)-%
  !(\varepsilon_F-\varepsilon_2)(\varepsilon_F-\varepsilon'_1)}%
  !{(\varepsilon_2-\varepsilon_4)(\varepsilon'_2-\varepsilon'_1)-%
  !(\varepsilon_2-\varepsilon_1)(\varepsilon'_2-\varepsilon'_4)})\\
  ! \vec{x}_{19}=&(%
  !\tfrac{(\varepsilon_F-\varepsilon_1)(\varepsilon_F-\varepsilon'_3)-%
  !(\varepsilon_F-\varepsilon_3)(\varepsilon_F-\varepsilon'_1)}%
  !{(\varepsilon_3-\varepsilon_2)(\varepsilon'_3-\varepsilon'_1)-%
  !(\varepsilon_3-\varepsilon_1)(\varepsilon'_3-\varepsilon'_2)},%
  !\tfrac{(\varepsilon_F-\varepsilon_2)(\varepsilon_F-\varepsilon'_1)-%
  !(\varepsilon_F-\varepsilon_1)(\varepsilon_F-\varepsilon'_2)}%
  !{(\varepsilon_3-\varepsilon_2)(\varepsilon'_3-\varepsilon'_1)-%
  !(\varepsilon_3-\varepsilon_1)(\varepsilon'_3-\varepsilon'_2)},0)\\
  ! \vec{x}_{20}=&(%
  !\tfrac{(\varepsilon_F-\varepsilon_3)(\varepsilon_F-\varepsilon'_4)-%
  !(\varepsilon_F-\varepsilon_4)(\varepsilon_F-\varepsilon'_3)}%
  !{(\varepsilon_3-\varepsilon_4)(\varepsilon'_3-\varepsilon'_2)-%
  !(\varepsilon_3-\varepsilon_2)(\varepsilon'_3-\varepsilon'_4)},
  !\tfrac{(\varepsilon_F-\varepsilon_4)(\varepsilon_F-\varepsilon'_2)-%
  !(\varepsilon_F-\varepsilon_2)(\varepsilon_F-\varepsilon'_4)}%
  !{(\varepsilon_3-\varepsilon_4)(\varepsilon'_3-\varepsilon'_2)-%
  !(\varepsilon_3-\varepsilon_2)(\varepsilon'_3-\varepsilon'_4)},\nonumber\\
  !&\tfrac{(\varepsilon_F-\varepsilon_2)(\varepsilon_F-\varepsilon'_3)-%
  !(\varepsilon_F-\varepsilon_3)(\varepsilon_F-\varepsilon'_2)}%
  !{(\varepsilon_3-\varepsilon_4)(\varepsilon'_3-\varepsilon'_2)-%
  !(\varepsilon_3-\varepsilon_2)(\varepsilon'_3-\varepsilon'_4)})
  ! \end{align}
  ! \end{subequations} 
  ! 
  !From the intersections listed in the last equation, the subroutine 
  !automatically eliminates those that either do not belong to the surface
  !of the tetrahedron or are indefinite.
  ! 
  !We use a parameter 'ndabc' here, since in the later section, for some
  !cases we need to get the region of $e1<ef$ and $e2>ef$ by minusing the 
  !$e1<ef$ and $e2<ef$ region from the $e1<ef$ region. So, for $ndabc=1$ and
  !$ndabc=3$, we get all the nodes for the intersection of the six planes. 
  !For $ndabc=2$, we just get that of the five planes except $e2=ef$, to get 
  !the region of $e1<ef$
  !
  implicit none
  integer(1), intent(out):: ntype(20)      ! idem as pl, but for internal use
  real(8),    intent(out):: intnodes(3,20) ! the coordinates of the intersections of the planes
  integer,    intent(out):: nnod           ! number of nodes defining the polyhedron
  integer,    intent(in) :: iwrite
  real(8),    intent(in) :: e(4)           ! band energies at k
  real(8),    intent(in) :: f(4)           ! band energies at k-q
  real(8),    intent(in) :: ef             ! the fermi energy
  !
  integer, parameter :: fout = 6
  real(8), parameter :: zerotol=1.0d-12
  !
  integer :: i,j,k,inod,iint
  real(8) :: edif, fdif, denom, denom1, denom2, denom3, efdif, nt(3),sumnt,efef
  !
  !     inline functions
  edif(i,j) = e(i)-e(j)      
  fdif(i,j) = f(i)-f(j)
  efdif(i,j,k) = edif(i,j)*fdif(k,j) - edif(k,j)*fdif(i,j)
  efef(i,j) = (ef-e(i))*(f(j)-f(i)) - (ef-f(i))*(e(j)-e(i))
  !
  !
  intnodes(1:3,1:20) = 0.0d0
  ntype(1)=7
  ntype(2)=11
  ntype(3)=13
  ntype(4)=14
  ntype(5)=19
  ntype(6)=21
  ntype(7)=22
  ntype(8)=25
  ntype(9)=26
  ntype(10)=28
  ntype(11)=35
  ntype(12)=37
  ntype(13)=38
  ntype(14)=41
  ntype(15)=42
  ntype(16)=44
  ntype(17)=49
  ntype(18)=50
  ntype(19)=52
  ntype(20)=56
  !
  !  The four corners of the tetrahedron:
  !  intnodes(1) = plane1^plane2^plane3
  !  intnodes(2) = plane1^plane2^plane4
  !  intnodes(3) = plane1^plane3^plane4
  !  intnodes(4) = plane2^plane3^plane4
  !      
  inod=0
  iint=0
  do i=2,4
     intnodes(5-i,i) = 1.0d0
  enddo
  inod=4 
  iint=4
  !
  !     Intersections of the plane approximating the Fermi surf. at k 
  !     with the stright lines defining the edges of the
  !     tetrahedron
  !
  !     intnodes(5) = plane1^plane2^plane5
  !     intnodes(6) = plane1^plane3^plane5
  !     intnodes(7) = plane2^plane3^plane5
  ! 
  do j = 2,4
     i = 6-j
     iint=iint+1
     denom = edif(i,1)
     nt(1:3) = 0.0d0
     if (dabs(denom) > zerotol) then
        nt(i-1) = (ef-e(1))/denom
        if ( (nt(i-1) > 0.0d0) .and. (nt(i-1) <= 1.0d0 ) ) then
           inod=inod+1
           intnodes(1:3,inod) = nt(1:3)
           ntype(inod) = ntype(iint)
        endif
     endif
  enddo
  !      
  !     intnodes(8) = plane1^plane4^plane5
  !     intnodes(9) = plane2^plane4^plane5
  !     intnodes(10) = plane3^plane4^plane5
  ! 
  do k = 1,3
     i = mod(k,3)+1
     j = mod(i,3)+1
     iint = iint+1
     denom = edif(i+1,j+1)
     nt(1:3) = 0.0d0
     if ( dabs(denom) > zerotol ) then
        nt(i) = (ef-e(j+1))/denom
        nt(j) = (e(i+1)-ef)/denom
        if ( (nt(i)>=0.0d0).and.(nt(j)>=0.0d0) ) then
           inod = inod+1
           intnodes(1:3,inod) = nt(1:3)
           ntype(inod)=ntype(iint)
        endif
     endif
  enddo
  !
  !     Intersections of the plane approximating the Fermi surf. at k - q
  !     with the stright lines defining the edges of the
  !     tetrahedron
  !
  !     intnodes(11) = plane1^plane2^plane6
  !     intnodes(12) = plane1^plane3^plane6
  !     intnodes(13) = plane2^plane3^plane6
  ! 
  do j=2,4
     i=6-j
     iint = iint+1
     denom=fdif(i,1)
     nt(1:3)=0.0d0
     if (dabs(denom) > zerotol) then
        nt(i-1) = (ef-f(1))/denom
        if ( (nt(i-1) >= 0.0d0) .and. (nt(i-1) <= 1.0d0) ) then
           inod=inod+1
           intnodes(1:3,inod) = nt(1:3)
           ntype(inod) = ntype(iint)
        endif
     endif
  enddo
  !          
  !     intnodes(14) = plane1^plane4^plane6
  !     intnodes(15) = plane2^plane4^plane6
  !     intnodes(16) = plane3^plane4^plane6
  do k=1,3
     i=mod(k,3)+1
     iint=iint+1
     j=mod(i,3)+1
     denom=fdif(i+1,j+1)
     nt(1:3)=0.0d0
     if(dabs(denom) > zerotol)then
        nt(i) = (ef-f(j+1))/denom
        nt(j) = (f(i+1)-ef)/denom
        if ( (nt(i) >= 0.0d0).and.(nt(j) >= 0.0d0))then
           inod=inod+1
           intnodes(1:3,inod)=nt(1:3)
           ntype(inod)=ntype(iint)
        endif
     endif
  enddo
  !
  !     Intersections of the two planes approximating the Fermi surf. and the 
  !     planes defining the tetrahedron surface
  !
  !     intnodes(17) = plane1^plane5^plane6
  !     intnodes(18) = plane2^plane5^plane6
  !     intnodes(19) = plane3^plane5^plane6
  ! 
  do k=1,3
     i=mod(k,3)+1
     iint=iint+1
     j=mod(i,3)+1
     denom1=efdif(i+1,1,j+1)
     denom2=efdif(j+1,1,i+1)
     nt(1:3)=0.0d0
     if((dabs(denom1).gt.zerotol).and.(dabs(denom2).gt.zerotol))then
        nt(i)=efef(1,j+1)/denom1
        nt(j)=efef(1,i+1)/denom2
        sumnt=nt(i)+nt(j)
        if ( (nt(i) > 0.0d0).and.(nt(j) > 0.0d0) .and.(sumnt <= 1.0d0) ) then
           inod=inod+1
           intnodes(1:3,inod)=nt(1:3)
           ntype(inod)=ntype(iint)
        endif
     endif
  enddo
  !
  !     intnodes(20) = plane4^plane5^plane6
  ! 
  denom1=efdif(2,4,3)
  denom2=efdif(3,4,2)
  denom3=efdif(4,2,3)
  iint=iint+1
  if (((dabs(denom1).gt.zerotol).and.(dabs(denom2).gt.zerotol)).and.(dabs(denom3).gt.zerotol)) then
     nt(1)=efef(4,3)/denom1
     nt(2)=efef(4,2)/denom2
     nt(3)=efef(2,3)/denom3
     if ( (nt(1).gt.0.0d0) .and. (nt(2).gt.0.0d0) .and. (nt(3).gt.0.0d0) ) then
        inod=inod+1
        intnodes(1:3,inod) = nt(1:3)
        ntype(inod) = ntype(iint)
     endif
  endif
  nnod = inod  
  if (iwrite.eq.1) then
     write(fout,*)'all nodes (from surfnodes)'
     do inod=1,nnod
        write(fout,'(i4,3f18.10,i4,1x,b6.6)') inod, intnodes(1:3,inod),ntype(inod),ntype(inod)
     enddo
     write(fout,*)
  endif
end subroutine surfnodes

subroutine unrepnodes(iwrite, ntype, intnodes, nnod)
  ! This subroutines select, from the intersections of the planes, only
  !those that belong to the surface of the tetrahedron. If one of them is
  !repeated, which means that there are more than three planes intersecting,
  !then it set the corresponding bits of \verb{ntype} to 1.
  implicit none
  integer(1), intent(inout):: ntype(20)      ! idem as pl, but for internal use
  real(8),    intent(inout):: intnodes(3,20) ! the coordinates of the intersections of the planes
  integer,    intent(inout):: nnod           ! number of nodes defining the polyhedron
  integer,    intent(in)   :: iwrite
  !
  integer, parameter :: fout = 6
  integer(1) :: i,j,k,ib
  integer(1) :: bti,btj
  integer    :: inod
  real(8)    :: intnodi(3)
  real(8)    :: sumnt
  integer(1), parameter :: one=1
  real(8),    parameter :: zerotol = 1.0d-8
  intrinsic ibits
  intrinsic ibset 
  i = 1
  do while (i < nnod)
     intnodi(1:3) = intnodes(1:3,i)
     j=i+1
     do while (j <= nnod)
        sumnt=0.0d0
        do k=1,3
           sumnt=sumnt+dabs(intnodi(k)-intnodes(k,j))   ! this is mainly to 
           !                                    see whether these two points are too close
        enddo
        if (sumnt < zerotol) then                    ! if they are close enough
           !                                   to be viewed as the same point, delete one of
           !                                   them and keep the other.
           do ib=0,5
              btj=ibits(ntype(j),ib,one)
              bti=ibits(ntype(i),ib,one)
              if((btj.eq.1).and.(bti.eq.0))ntype(i)=ibset(ntype(i),ib)
           enddo
           do k=j+1,nnod
              intnodes(1:3,k-1)=intnodes(1:3,k)
              ntype(k-1)=ntype(k)
           enddo
           intnodes(1:3,nnod)=0.0d0
           ntype(nnod)=0
           nnod=nnod-1
        else
           j=j+1
        endif
     enddo
     i=i+1
  enddo
  if (iwrite.eq.1) then
     write(fout,*)'unrepeated nodes (from unrepnodes)'
     do inod=1,nnod
        write(fout,'(i4,3f18.10,i4,1x,b6.6)')inod, intnodes(1:3,inod),ntype(inod),ntype(inod)
     enddo
     write(fout,*)
  endif
end subroutine unrepnodes

subroutine relnodes(iwrite,info,  ntype, intnodes, nnod, e, f, ef)
  ! This subroutine select the relevant intnodes, that is, only those defining
  ! the polyhedron that has to be integrated. The intnodes selected are those for
  ! which $\varepsilon(\vec{k})\le \varepsilon_F$ and 
  ! $\varepsilon'(\vec{k})\ge \varepsilon_F$ if ndabc=1. $\varepsilon(\vec{k})
  ! \le \varepsilon_F$ when ndabc=2.  $\varepsilon(\vec{k})\le \varepsilon_F$ and
  ! $\varepsilon'(\vec{k})\le \varepsilon_F$ when ndabc=3. So the first region equals
  ! the second region minused by the third region.
  implicit none
  integer,    intent(in)   :: iwrite         ! flag for writing extra output
  integer,    intent(out)  :: info           ! flag for correct execution
  integer(1), intent(inout):: ntype(20)      ! idem as pl, but for internal use
  real(8),    intent(inout):: intnodes(3,20) ! the coordinates of the intersections of the planes
  integer,    intent(inout):: nnod           ! number of nodes defining the polyhedron
  real(8),    intent(in) :: e(4)           ! band energies at k
  real(8),    intent(in) :: f(4)           ! band energies at k-q
  real(8),    intent(in) :: ef             ! the fermi energy
  !
  integer, parameter :: fout = 6
  integer(4) :: inod, jnod, maxnod
  real(8) :: e1, e2
  real(8), parameter :: zerotol=1.00d-12 
  real(8), dimension(3) :: node_coord 
  logical :: keep_node, e1_occup, e2_unocc
  real(8), external :: tlinap
  !
  info = 0
  inod = 1
  maxnod = nnod
  if(iwrite.eq.1)then
     write(fout,*)'energies at the tetracorn'
     write(fout,'("e_occ   =",4f18.10)')e
     write(fout,'("e_unocc =",4f18.10)')f
  endif
  do while (inod.le.maxnod)
     node_coord(1:3)=intnodes(1:3,inod)
     e1 = tlinap(node_coord,e)-ef
     e2 = tlinap(node_coord,f)-ef
     select case (ntype(inod))
        !
        ! The corners of the tetrahedron, both conditions have to be satisfied        
        !
     case(7:15) 
        e1_occup = (e1.lt.0.0d0).or.(dabs(e1).le.zerotol)
        e2_unocc = (e2.gt.0.0d0).or.(dabs(e2).le.zerotol)
        keep_node = e1_occup .and. e2_unocc
        !
        ! Intersections of plane 5 (e1=ef) with two planes of the tetrahedron
        ! only e2_unocc has to be satisfied.
        !
     case(16:31)
        e2_unocc = (e2.gt.0.0d0).or.(dabs(e2).le.zerotol)
        keep_node = e2_unocc
        !
        ! Intersections of plane 6 (e2=ef) with two planes of the tetrahedron
        ! only e1_occup has to be satisfied.
        !
     case(32:47)
        e1_occup = (e1.lt.0.0d0).or.(dabs(e1).le.zerotol)
        keep_node = e1_occup
        !
        ! Intersections of plane 5 (e1=ef) and 6 (e2=ef) with one or more planes 
        ! of the tetrahedron: Both conditions are already fullfiled
        !
     case(48:63)
        keep_node=.true.
        !
        ! The value of ntype is wrong, write message and send error signal.
        !          
     case default
        write(fout,*)'ERROR: ntype(inod) =',ntype(inod),' > 63 !!!'
        info = 1  
     end select
     if(iwrite.eq.1)write(fout,'(i4,2g18.10,l2)')inod,e1,e2,keep_node
     if(keep_node)then  
        inod=inod+1
     else ! remove the node from the list
        do jnod=inod+1,maxnod
           intnodes(1:3,jnod-1)=intnodes(1:3,jnod)
           ntype(jnod-1)=ntype(jnod)
        enddo
        maxnod=maxnod-1
     endif
  enddo ! inod
  !if(iwrite.eq.1) call flushbuf(6)
  !      
  ! Set the variables of the removed nodes to zero, to avoid errors.
  !
  do inod=maxnod+1,nnod
     intnodes(1:3,inod)=0.0d0
     ntype(inod)=0
  enddo
  ! 
  ! Set the number of nodes to the relevant ones
  !      
  nnod=maxnod
  if (iwrite.eq.1) then
     write(fout,*)'relevant nodes (from relnodes)'
     do inod=1,nnod
        write(fout,'(i4,3f18.10,i4,1x,b6.6)')inod, intnodes(1:3,inod),ntype(inod),ntype(inod)
     enddo
     write(fout,*)
     !call flushbuf(6)
  endif
end subroutine relnodes

real(8) function tlinap(x,f)
  ! This function evalueates a function $f(\vec{x})$  at a given point \verb"x"
  ! (in internal coordinates) inside the tetrahedron aproximating it linearly 
  ! from the known values of the function at the corners of the tetrahedron 
  ! (\verb"f") using the isoparametrization.
  !
  ! \begin{equation}
  ! \bar{f}(\vec{x})=f_0+\sum\limits_{i=1}^3{(f_i-f_0)x_i}
  ! \end{equation}
  implicit none
  real(8), intent(in) :: x(1:3) ! The coordinates of the point
  real(8), intent(in) :: f(1:4) ! The values of the function at the corners of the tetrahedron
  integer(4) :: i
  real(8) :: fx
  real(8) :: df
  fx = f(1)
  do i=1,3
     df = f(i+1)-f(1)
     fx = fx + df*x(i)
  enddo
  tlinap=fx
end function tlinap

subroutine sortnodes(info, ntype, nnod, index, nd2)
  ! This subroutine sort the nodes of each region so that it can be devided
  ! into tetrahedron in order. For the case of nnod=4, it is already done.
  ! For nnod=5, we devide it into two tetrahedrons with index(inod,1:2) indicates
  ! the order of the node in the tetrahedron. For example, if index(inod,2)=1,
  ! we should take this node as the first node in the second tetrahedron. 
  ! For nnod=6, we can devide it into three tetrahedron, the nodes should be ordered
  ! in a way that 1,2,3,4 form a tetrahedron, 2,,3,4,5 form a tetrahedron and 
  ! 3,4,5,6 form the other one. These three together make up the prism. 
  ! For nnod=7, we can devide it into two penta which are ordered in the way as
  ! nnod=5. For nnod=8, we can devide it into two prism which are ordered in the
  ! way as nnod=6 case.  
  !
  !use polyhedron, only: index
  implicit none
  integer(4), intent(out)  :: info
  integer(1), intent(in)   :: ntype(20)      ! idem as pl, but for internal use
  integer,    intent(in)   :: nnod           ! number of nodes defining the polyhedron
  integer,    intent(in)   :: nd2
  integer(1), intent(out)  :: index(nnod,nd2)
  !
  integer(1) :: ibit, inod,  ip, iind
  integer(4) :: is, it, i, iq1
  integer(1) :: ibt(3) 
  integer(1) :: it1, it2, is1
  integer(1) :: triangle(6)
  integer(1) :: square(6)
  integer(1) :: penta(2)
  integer(4) :: spl(0:5)
  integer(1) ::  nt, signbit, signbita,signbitb,signnod,signnodb
  integer(1) :: sigtri(6)
  integer(1), parameter :: one=1
  integer(1), parameter :: four=4
  integer(1), parameter :: five=5
  intrinsic btest
  !
  !
  spl(0:5)=0
  it=0
  is=0
  ip=0
  info=0
  select case(nnod)
  case(5)     
     !allocate(index(5,1))
     !--------------------------------------------------------------------
     !       finding how many triangles and squares to form the region
     !                               begin
     !--------------------------------------------------------------------
     do ibit=0,5      ! here we use the bit to represent the six planes
        do inod=1,5   ! inod means each node
           nt=ntype(inod)
           if(btest(nt,ibit))spl(ibit)=spl(ibit)+1
        enddo
        select case(spl(ibit))
        case(3)       ! for one plane, if there are 3 nodes on it.
           it=it+1
           triangle(it)=ibit
        case(4)       ! for one plane, if there are 4 nodes on it.
           is=is+1
           square(is)=ibit
           signbit=ibit
        end select
     enddo
     !--------------------------------------------------------------------
     !       finding how many triangles and squares to form the region
     !                               end
     !--------------------------------------------------------------------
     index(1:5,1)=0
     if(is.ne.0) then     ! the case when there is one node
        !                                     out of one surface while the other
        !                                     four in.  
        !--------------------------------------------------------------------
        !  When the region is formed by 4 triangles and one square, begin
        !--------------------------------------------------------------------
        is1=1
        do inod=1,5
           if(btest(ntype(inod),triangle(1))) then
              if(.not.(btest(ntype(inod),signbit))) then
                 index(inod,1)=2
              else
                 index(inod,1)=is1
                 is1=is1+2
                 signnod=inod
              endif
           else
              continue
           endif
        enddo
        do i=2,4
           if(btest(ntype(signnod),triangle(i)))signbit=triangle(i)
        enddo
        do inod=1,5
           if(index(inod,1).eq.0) then
              if(btest(ntype(inod),signbit)) then
                 index(inod,1)=5
              else
                 index(inod,1)=4
              endif
           else
              continue
           endif
        enddo
        !-------------------------------------------------------------------------
        !   when the region is formed by 4 triangles and one square, end
        !--------------------------------------------------------------------------
     else                  ! for the case when there are 3 nodes in 
        !                               each surface.
        !-------------------------------------------------------------------------
        !   when the region is formed by 6 triangles, begin
        !--------------------------------------------------------------------------
        do inod=1,5
           signbit=0
           do ibit=0,5
              if(btest(ntype(inod),ibit))signbit=signbit+1
           enddo
           it1=1
           it2=1
           if(signbit.eq.3) then
              index(inod,1)=1+(it1-1)*4
              it1=it1+1
           else
              index(inod,1)=1+it2
              it2=it2+1
           endif
        enddo
        !-------------------------------------------------------------------------
        !   when the region is formed by 6 triangles, end
        !--------------------------------------------------------------------------
     endif
  case(6)
     !allocate(index(6,1))
     index(1:6,1)=0
     !--------------------------------------------------------------------
     !       finding how many triangles and squares to form the region
     !                               begin
     !--------------------------------------------------------------------
     do ibit=0,5
        do inod=1,6
           nt=ntype(inod)
           if(btest(nt,ibit))spl(ibit)=spl(ibit)+1
        enddo
        select case(spl(ibit))
        case(3)
           it=it+1
           triangle(it)=ibit
        case(4)
           is=is+1
           square(is)=ibit
        end select
     enddo
     !--------------------------------------------------------------------
     !       finding how many triangles and squares to form the region
     !                               end
     !--------------------------------------------------------------------
     select case(is)
     case(2)
        !-------------------------------------------------------------------------
        !   when the region is formed by 4 triangles and 2 squares, begin
        !-------------------------------------------------------------------------
        it=1
        do inod=1,6
           if (btest(ntype(inod),square(1)) .and. btest(ntype(inod),square(2))) then
              index(inod,1)=2+(it-1)*3
              if(it.eq.1) signnod=inod
              if(it.eq.2) signnodb=inod               
              it=it+1
           else
              continue
           endif
        enddo
        sigtri(1:6)=0
        do inod=1,6
           do it=1,4
              if(btest(ntype(inod),triangle(it))) then
                 sigtri(inod)=sigtri(inod)+1
              else
                 continue
              endif
           enddo
        enddo
        do it=1,4
           if(btest(ntype(signnod),triangle(it)))signbita=it
           if(btest(ntype(signnodb),triangle(it)))signbitb=it
        enddo
        do inod=1,6
           if(index(inod,1).eq.0) then
              if((sigtri(inod).eq.3).and.btest(ntype(inod),triangle(signbita)))index(inod,1)=4
              if((sigtri(inod).eq.2).and.btest(ntype(inod),triangle(signbita)))index(inod,1)=1
              if((sigtri(inod).eq.3).and.btest(ntype(inod),triangle(signbitb)))index(inod,1)=3
              if((sigtri(inod).eq.2).and.btest(ntype(inod),triangle(signbitb)))index(inod,1)=6
           else
              continue
           endif
        enddo
        !-------------------------------------------------------------------------
        !   when the region is formed by 4 triangles and 2 squares, end
        !-------------------------------------------------------------------------
     case(3)
        do inod=1,6
           do it=1,2
              if(btest(ntype(inod),triangle(it)))then
                 if ( btest(ntype(inod),square(1)).and. btest(ntype(inod),square(2)) ) then
                    index(inod,1)=1+3*(it-1)
                 else if(btest(ntype(inod),square(1))) then
                    index(inod,1)=2+3*(it-1)
                 else if(btest(ntype(inod),square(2))) then
                    index(inod,1)=3+3*(it-1)
                 else
                    stop 'error in sortnodes, nnod=6'
                 endif
              else
                 continue
              endif
           enddo
        enddo
     end select
  case(7) 
     !allocate(index(7,2))
     !-----------------------------------------------------------------------
     !   Finding how the region is formed, begin
     !-----------------------------------------------------------------------
     index(1:7,1:2)=0
     do ibit=0,5
        do inod=1,7
           if(btest(ntype(inod),ibit))spl(ibit)=spl(ibit)+1
        enddo
        select case(spl(ibit))
        case(3)
           it=it+1
           triangle(it)=ibit
        case(4)
           is=is+1
           square(is)=ibit
        case(5)
           ip=ip+1
           penta(ip)=ibit
        end select
     enddo
     !-----------------------------------------------------------------------
     !    Finding how the region is formed, end
     !-----------------------------------------------------------------------
     select case(it)
     case(2)       ! when it is formed by 2 triangles and 4 squares
        !-----------------------------------------------------------------------
        ! When the region is formed by 2 triangles and 4 squares, begin
        !-----------------------------------------------------------------------
        select case(is)
        case(4)    ! 2 triangles and 4 squares
           do inod=1,7
              if (btest(ntype(inod),triangle(1)) .and. btest(ntype(inod),triangle(2))) then
                 index(inod,1:2)=2
                 signnod=inod
              else
                 continue
              endif
           enddo
           is1=0
           do ibit=0,5
              if(btest(ntype(signnod),ibit) .and. (spl(ibit).eq.4).and.(is1.eq.0)) then
                 signbita=ibit
                 is1=is1+1
              else if(btest(ntype(signnod),ibit) .and. (spl(ibit).eq.4)) then
                 signbitb=ibit
              else
                 continue
              endif
           enddo
           do inod=1,7
              if (inod.ne.signnod) then 
                 if(btest(ntype(inod),triangle(1)).and. btest(ntype(inod),signbita)) then
                    index(inod,1)=1
                 else if(btest(ntype(inod),triangle(1)).and. btest(ntype(inod),signbitb)) then
                    index(inod,1)=3
                 else if(btest(ntype(inod),triangle(2)).and. btest(ntype(inod),signbita)) then
                    index(inod,2)=1
                 else if(btest(ntype(inod),triangle(2)).and. btest(ntype(inod),signbitb)) then
                    index(inod,2)=3
                 else if(btest(ntype(inod),signbita)) then
                    index(inod,1:2)=4
                 else
                    index(inod,1:2)=5
                 endif
              else
                 continue
              endif
           enddo
        case default
           write(6,*) it, is, ip
           do inod=1,nnod
              write(6,*)inod,ntype(inod)
           enddo
           !call flushbuf(6)  
           write(6,*)'error in sortnodes.f90, inod=7 and it=2'
           info=72
        end select
        !-----------------------------------------------------------------------
        ! When the region is formed by 2 triangles and 4 squares, end
        !-----------------------------------------------------------------------
     case(3)  ! when it is made up of 3 triangles, 2 squares and 1 penta
        !--------------------------------------------------------------------------
        ! When the region is formed by 3 triangles and 2 squares and 1 penta, begin
        !--------------------------------------------------------------------------
        index(1:7,1:2)=0
        do inod=1,7
           ibit=0
           do it1=1,3
              if(btest(ntype(inod),triangle(it1)))ibit=ibit+1
           enddo
           if(btest(ntype(inod),penta(1))) then
              if(ibit.eq.2)index(inod,1)=1
           else
              if(ibit.eq.2) then
                 index(inod,1)=2
                 index(inod,2)=1
              else
                 index(inod,2)=4
                 iind=inod
              endif
           endif
        enddo
        do it2=1,3
           if(btest(ntype(iind),triangle(it2)))signbita=triangle(it2)
        enddo
        do inod=1,7
           if ((index(inod,1).eq.0).and.(index(inod,2).eq.0)) then
              if (btest(ntype(inod),signbita).and. btest(ntype(inod),square(1))) then
                 index(inod,2)=5
              else if(btest(ntype(inod),signbita).and. btest(ntype(inod),square(2))) then
                 index(inod,1)=5
                 index(inod,2)=2
              else if(btest(ntype(inod),square(1))) then
                 index(inod,1)=4
                 index(inod,2)=3
              else 
                 index(inod,1)=3
              endif
           else
              continue
           endif
        enddo
     case default
        write(6,*) "ERROR sortnodes: inod=7, it, is, ip=",it, is, ip
        stop 'ERROR in sortnodes'
     end select
     !--------------------------------------------------------------------------
     ! When the region is formed by 3 triangles and 2 squares and 1 penta, begin
     !--------------------------------------------------------------------------
  case(8)
     !allocate(index(8,2))
     !--------------------------------------------------------------------------
     !  Finding how the region is formed, begin
     !--------------------------------------------------------------------------
     do ibit=0,5
        do inod=1,8
           nt=ntype(inod)
           if(btest(nt,ibit))spl(ibit)=spl(ibit)+1
        enddo
        select case(spl(ibit))
        case(3)
           it=it+1
           triangle(it)=ibit
        case(4)
           is=is+1
           square(is)=ibit
        case(5)
           ip=ip+1
           penta(ip)=ibit
        end select
     enddo
     !--------------------------------------------------------------------------
     ! Finding how the region is formed, end
     !--------------------------------------------------------------------------    
     index(1:8,1:2)=0
     !        write(90,*)'itsp =',it, is,ip
     select case(is)
     case(2)         ! then it is made up of 2 squares, 2 triangles and 2 penta
        !--------------------------------------------------------------------------
        ! When the region is formed by 2 triangles and 2 squares and 2 penta, begin
        !--------------------------------------------------------------------------
        do inod=1,8
           do it1=1,2
              if(btest(ntype(inod),triangle(it1))) then
                 if(btest(ntype(inod),penta(1)).and. btest(ntype(inod),penta(2)))  then
                    index(inod,1)=1+(it1-1)*3
                 else if(btest(ntype(inod),penta(1))) then
                    index(inod,1)=2+(it1-1)*3
                 else if(btest(ntype(inod),penta(2))) then
                    index(inod,1)=3+(it1-1)*3
                 else
                    continue
                 endif
              else
                 continue
              endif
           enddo
           do iq1=1,2
              if(btest(ntype(inod),penta(iq1))) then
                 if(btest(ntype(inod),square(1)).and. btest(ntype(inod),square(2))) then
                    index(inod,2)=1+(iq1-1)*3
                 else if(btest(ntype(inod),square(1))) then
                    index(inod,2)=2+(iq1-1)*3
                 else if(btest(ntype(inod),square(2))) then
                    index(inod,2)=3+(iq1-1)*3
                 else
                    continue
                 endif
              else
                 continue
              endif
           enddo
        enddo
        !--------------------------------------------------------------------------
        ! When the region is formed by 2 triangles and 2 squares and 2 penta, end
        !--------------------------------------------------------------------------
     case(6)      ! when it is made up of 6 squares
        !--------------------------------------------------------------------------
        ! When the region is formed by 4 squares, begin
        !--------------------------------------------------------------------------
        is1=0
        ! choose the node 1 to be 1 in the first prism and it is not in the 
        ! second (node 1:0)
        index(1,1)=1
        index(1,2)=0
        ! set the three planes to which node 1 belongs as 1,2,3 in order of
        ! appearence         
        do ibit=0,5
           if(btest(ntype(1),ibit)) then
              is1=is1+1
              ibt(is1)=ibit
           endif
        enddo
        ! the rest of the nodes         
        do inod=2,8
           if (btest(ntype(inod),ibt(1))) then
              if(btest(ntype(inod),ibt(2)))then
                 index(inod,1:2)=2 ! node 2:2= intersection of plane 1 and 2
              else if(btest(ntype(inod),ibt(3)))then
                 index(inod,1:2)=3 ! node 3:3= intersection of plane 1 and 4
              else
                 index(inod,1)=0  ! node 0:1= the remaining node of plane 1
                 index(inod,2)=1
              endif
           else if(btest(ntype(inod),ibt(2)))then
              if(btest(ntype(inod),ibt(3)))then
                 index(inod,1)=4  ! node 4:0= intersection of plane 2 and 3
                 index(inod,2)=0
              else
                 index(inod,1:2)=5 ! node 5:5= the remaining node of plane 2
              endif
           else
              if(btest(ntype(inod),ibt(3)))then
                 index(inod,1:2)=6 ! node 6:6= the remaining node of plane 3
              else
                 index(inod,1)=0   ! node 0:4= the node that doesn´t belong
                 index(inod,2)=4   ! to any of the planes 1,2 or 3
              endif
           endif
        enddo
        !--------------------------------------------------------------------------
        ! When the region is formed by 4 squares, end
        !--------------------------------------------------------------------------
     case default
        write(6,*) "ERROR in sortnodes: nnod=8"
        stop "ERROR in sortnodes"
     end select
  end select
end subroutine sortnodes




!---------------------------- Sorting routines below --------------------------
subroutine binary_insertion_sort(array,n)
  implicit none
  real*8, intent(inout) :: array(n)
  integer, intent(in)   :: n
  !
  integer :: i,j,pos
  real*8 :: x
  integer, external :: binary_search
  do i=2,n
     x=array(i)
     pos=binary_search(array,x,1,i-1,n)
     do j=i,pos+1,-1
        array(j)=array(j-1)
     enddo
     array(pos)=x
  enddo
end subroutine binary_insertion_sort

integer function binary_search(array, x, start, end, n)
  implicit none
  integer, intent(in) :: start, end
  integer, intent(in) :: n
  real*8,  intent(in) :: x, array(n)
  !
  integer :: a,b,mid
  a=start
  b=end
  do
     if (a.eq.b) then
        if ( array(a).gt.x) then
           binary_search=a
           return 
        else
           binary_search=a+1
           return
        endif
     endif
     if (a.gt.b) then
        binary_search=a
        return
     end if
     mid=(a+b)/2
     if ( array(mid).lt.x) then
        a=mid+1
     else if ( array(mid).gt.x) then
        b=mid-1
     else
        binary_search=mid
        return
     endif
  end do
end function binary_search

subroutine binary_insertion_sort_indx(indx, array,n)
  implicit none
  integer, intent(out):: indx(n)
  real*8,  intent(in) :: array(n)
  integer, intent(in) :: n
  !
  integer :: i,j,pos
  real*8 :: x
  integer, external :: binary_search_indx
  do i=1,n
     indx(i)=i
  enddo
  do i=2,n
     x=array(indx(i))
     pos=binary_search_indx(indx,array,x,1,i-1,n)
     do j=i,pos+1,-1
        indx(j)=indx(j-1)
     enddo
     indx(pos)=i
  enddo
end subroutine binary_insertion_sort_indx

integer function binary_search_indx(indx, array, x, start, end, n)
  implicit none
  integer, intent(in) :: indx(n)
  integer, intent(in) :: start, end
  integer, intent(in) :: n
  real*8,  intent(in) :: x, array(n)
  !
  integer :: a,b,mid
  a=start
  b=end
  do
     if (a.eq.b) then
        if ( array(indx(a)).gt.x) then
           binary_search_indx=a
           return 
        else
           binary_search_indx=a+1
           return
        endif
     endif
     if (a.gt.b) then
        binary_search_indx=a
        return
     end if
     mid=(a+b)/2
     if ( array(indx(mid)).lt.x) then
        a=mid+1
     else if ( array(indx(mid)).gt.x) then
        b=mid-1
     else
        binary_search_indx=mid
        return
     endif
  end do
end function binary_search_indx



