
    subroutine maxwell_FDFD
    use the_whole_varibles
    IMPLICIT NONE    
    call find_func_cputime_1_of_2  !-------------------1/2

    CALL find_fvm
    CALL find_Dielectric_tensor
    CALL set_0
    
    
    do m= m_start,m_end,m_delta
        m_real=1.0d0*m_real
        call antenna
        CALL set_eq_b
        call set_A
        call do_reallocate
        CALL complex_coo_pardiso_solve(A_helicon,eq_xlast(m,:),eq_b)
        CALL totalsolve
    enddo
    CALL find_irf_now
    CALL update_and_calculate
    maxwellcount=maxwellcount+1
        
    call find_func_cputime_2_of_2(func_time(9))  !-----2/2
    end subroutine maxwell_FDFD
    
    
    
    !--------------------subroutine subroutine do_reallocate-----------------!
    subroutine do_reallocate
    use the_whole_varibles
    IMPLICIT NONE    
    if(maxwellcount>1) return   !only for the first calculation of FDFD
    if(maxwellcount==0)then
        nzmax_pt=merge(A_helicon.nzmax,nzmax_pt,A_helicon.nzmax>nzmax_pt)
        CALL reallocate(A_helicon,2*A_helicon.nzmax)    
        ! *2 to make sure the A_helicon of other m modes would not overflow!!!
    endif
    if(maxwellcount==1) call reallocate(A_helicon,nzmax_pt)
    
    end subroutine do_reallocate

    !--------------------subroutine subroutine find_Dielectric_tensor-----------------!
    subroutine find_Dielectric_tensor
    use the_whole_varibles
    implicit none
    Complex*16 :: omcr,omcz
    Complex*16 :: sigma_0,sigma_1,ve_iom2,ve_iom
    Complex*16 :: omce,omci,qim,qem !Defined and just used in find_Dielectric_tensor,mei
    Complex*16 :: ompi2,ompe2,omci2,omce2 !Defined and just used in find_Dielectric_tensor
    Complex*16 K_perp,K_phi,K_ll,coeff_ve_om,coeff_vi_om
    real*8     b0r,b0z,b0_abs,b02,om2,ve,vi,coeff_ompi2,coeff_ompe2
    !--------------------------coefficient--------------------------!


    coeff_ompe2= qe*qe/(epson0*me)
    coeff_ompi2= qe*qe/(epson0*mi)
    !mei=me/mi
    qem=qe/me
    qim=qe/mi
    om2=om*om;

    ep=0.0D0*i
    si=0.0D0*i
    ep(:,:,1)=1.0D0+0.0D0*i
    ep(:,:,5)=1.0D0+0.0D0*i
    ep(:,:,9)=1.0D0+0.0D0*i


    do ir=1,nr
        do iz=1,nz
            call find_ve_vi(ve,vi)

            b0r=b0_DC(ir,iz,1)
            b0z=b0_DC(ir,iz,3)
            b0_abs=b0_DC(ir,iz,4)
            b02=b0_abs*b0_abs
            coeff_ve_om=1.0D0+i*(ve/om)
            ompe2=coeff_ompe2*ne(ir,iz)/coeff_ve_om
            omce=qem*b0_abs/coeff_ve_om
            omce2=omce*omce
            coeff_vi_om=1.0D0+i*(vi/om)
            ompi2=coeff_ompi2*ne(ir,iz)/coeff_vi_om
            if(iswitch_ji_on==0)ompi2=0.0d0
            omci=qim*b0_abs/coeff_vi_om
            omci2=omci*omci

            K_perp=1-ompe2/(om2-omce2)-ompi2/(om2-omci2)
            K_phi=ompe2*omce/(om*(om2-omce2))+ompi2*omci/(om*(om2-omci2))
            K_ll=1-ompe2/om2-ompi2/om2


            si(ir,iz,1,1)=b0z**2/b02*K_perp+b0r**2/b02*K_ll-1.0D0
            si(ir,iz,1,2)=-b0z/b0_abs*i*K_phi
            si(ir,iz,1,3)=b0r*b0z/b02*(K_ll-K_perp)
            si(ir,iz,2,2)=K_perp-1.0D0
            si(ir,iz,2,3)=-b0r/b0_abs*i*K_phi
            si(ir,iz,3,3)=b0z**2/b02*K_ll+b0r**2/b02*K_perp-1.0D0
            si(ir,iz,2,1)=-1.0D0*si(ir,iz,1,2)
            si(ir,iz,3,1)=si(ir,iz,1,3)
            si(ir,iz,3,2)=-1.0D0*si(ir,iz,2,3)
        enddo
    enddo
    si=-si*i*epson0*om


    si(nrp,:,1,:)=0.0D0*i

    ep(:,:,1)=1.0D0+i/(epson0*om)*si(:,:,1,1)
    ep(:,:,5)=1.0D0+i/(epson0*om)*si(:,:,2,2)
    ep(:,:,9)=1.0D0+i/(epson0*om)*si(:,:,3,3)
    ep(:,:,4)=i/(epson0*om)*si(:,:,2,1)
    ep(:,:,7)=i/(epson0*om)*si(:,:,3,1)
    ep(:,:,8)=i/(epson0*om)*si(:,:,3,2)
    ep(:,:,2)=-ep(:,:,4)
    ep(:,:,3)=ep(:,:,7)
    ep(:,:,6)=-ep(:,:,8)

    !if(nr_vac(1)<nr) si(nr_vac(1)+1:nr,1:nz_vac(1),:,:)=0.0D0*i
    !!if(1<nzp1) si(:,1:nzp1-1,:,:)=0.0D0*i
    !!if(nzp2<nz) si(:,nzp2+1:nz,:,:)=0.0D0*i
    !si(nr_vac(1),1:nz_vac(1),1,:)=0.0D0*i
    !si(nr_vac(1),1:nz_vac(1),:,1)=0.0D0*i

    continue
    end subroutine find_Dielectric_tensor

    !--------------------subroutine set_0-----------------!
    subroutine set_0
    use the_whole_varibles
    implicit none
    ptotal=0.0D0
    ptotal_boundary=0.0D0
    ptotal_inside=0.0D0
    !ptotal_BI=0.0D0
    peff=0.0D0*i
    power_depo=0.
    ptotm=0.0D0
    A_helicon.col=0
    A_helicon.row=0
    A_helicon.val=0.0D0*i
    e_m=0.0D0; e_int=0.0D0; ja_m_nor=0.0D0;  ja_record=0.0D0;
    e_AC=0.0D0;e_record=0.0D0;ja_AC=0.0D0;!ptotm=0.0D0
    b_m=0.0D0;b_AC=0.0D0;b_record=0.0D0
    jp_m=0.0D0;jp_record=0.0D0;jp_AC=0.0D0
    !divE_m=0.0D0; dive_AC=0.0D0 ;dive_record=0.0D0

    end subroutine set_0

    !--------------------subroutine set_eq_b-----------------!
    subroutine set_eq_b
    use the_whole_varibles
    IMPLICIT NONE
    eq_b(:)=0.0D0*i
    if(switch_source_type==0)then
        ir=1;i3=1;
        if(m==1.or.m==-1)then
            eq_b(3*(nz*ir+iz-nz)-3+i3)=1.0D0/(2.0D0*DBesselJZero(m))+0.0D0*i!Er(ir,1)
        else
            eq_b(3*(nz*ir+iz-nz)-3+i3)=0.0D0*i!Er(ir,1)
        endif
        i3=2
        eq_b(3*(nz*ir+iz-nz)-3+i3)=i/(2.0D0*DBesselJZero(m))*(Bessel_Jn(m-1,DBesselJZero(m)*r(ir)/rl)&
            &-Bessel_Jn(m+1,DBesselJZero(m)*r(ir)/rl))+0.0D0*i!Eth(ir,1)
        Do ir=2,nr
            i3=1
            eq_b(3*(nz*ir+iz-nz)-3+i3)=m_real*rl/(DBesselJZero(m)**2*r(ir))*Bessel_Jn(m,DBesselJZero(m)*r(ir)/rl)+0.0D0*i!Er(ir,1)
            i3=2
            eq_b(3*(nz*ir+iz-nz)-3+i3)=i/(2.0D0*DBesselJZero(m))*(Bessel_Jn(m-1,DBesselJZero(m)*r(ir)/rl)&
                &-Bessel_Jn(m+1,DBesselJZero(m)*r(ir)/rl))+0.0D0*i!Eth(ir,1)
        EndDo
    elseif(switch_source_type==1)then
        Do ir=1,nr
            Do iz=1,nz
                Do i3=1,n3
                    eq_b(3*(nz*ir+iz-nz)-3+i3)=iom*ja_m_nor(ir,iz,i3)
                EndDo
            EndDo
        EndDo
    endif
    continue
    !pause
    end subroutine set_eq_b
    !--------------------function DBesselJZero-----------------!
    function DBesselJZero(nbessel)
    implicit none
    integer*4 :: nbessel
    real*8 :: DBesselJZero
    real*8 :: listpt(6)=(/ 1.841D0, 3.054D0, 4.201D0, &
        &5.317D0, 6.415D0, 7.501D0 /)
    if(nbessel==0)then
        DBesselJZero=3.831705970207512D0
    elseif(nbessel>=1 .and. nbessel<=6)then
        DBesselJZero=listpt(nbessel)
    else
        write(*,*)'DBesselJZero(x_plasma) xmax=6!!!!!! Please Check'
        pause
        DBesselJZero=0
    endif
    return
    end function DBesselJZero


    !--------------------subroutine set_A-----------------!
    subroutine set_A
    use the_whole_varibles
    implicit none
    integer*4 :: lcount,scount
    A_helicon.n=3*nr*nz
    A_helicon.nzmax=0
    A_helicon.row=0
    A_helicon.col=0
    A_helicon.val=0.0D0*i
    iswitch_output_Region=1


    if(iswitch_output_Region==1)open(unit=545,file='Region.txt',status='REPLACE',iostat=ierror)
    Call Find_Region(nr,nz,n_vac,nr_vac,nr_met,nz_vac,nrz_diploe,iswtich_inner_dipole,FindRegion)
    do ir=1,nr
        do iz=1,nz
            !!!plasma region(2 <= ir <= nrp-1)
            Region=FindRegion(ir,iz)
            if(iswitch_output_Region==1)then
                write(545,"(I5,I5,I5)")ir,iz,Region!TEST
            endif
            if(Region==1)then
                lcount=3*(nz*ir+iz-nz)-2!!!Er(ir,iz)
                r(:)=r(:)+1.0D0/2.0D0 * dr
                call set(lcount,0,m_real**2/r(ir)**2+2.0e0/dz22-omc2*(ep(ir+1,iz,1)+ep(ir,iz,1))/2.0e0+0.0e0*i,A_helicon)!Er(ir,iz)
                call set(lcount,0+3,-1.0D0/dz22+0.0D0*i,A_helicon)!Er(ir,iz+1)
                call set(lcount,0-3,-1.0D0/dz22+0.0D0*i,A_helicon)!Er(ir,iz-1)
                call set(lcount,1,im/r(ir)**2 *(1.0D0/2.0D0-r(ir)/dr)-omc2*ep(ir,iz,2)/2.0D0+0.0D0*i,A_helicon)!Eth(ir,iz)
                call set(lcount,1+3*nz,im/r(ir)**2 *(1.0D0/2.0D0+r(ir)/dr)-omc2*ep(ir+1,iz,2)/2.0D0+0.0D0*i,A_helicon)!Eth(ir+1,iz)
                call set(lcount,2+3*nz,1.0D0/(dr*dz)-omc2*(ep(ir+1,iz,3)+ep(ir+1,iz+1,3))/8.0D0+0.0D0*i,A_helicon)!Ez(ir+1,iz)
                call set(lcount,2-3,1.0D0/(dr*dz)-omc2*(ep(ir,iz-1,3)+ep(ir+1,iz,3))/8.0D0+0.0D0*i,A_helicon)!Ez(ir,iz-1)
                call set(lcount,2,-1.0D0/(dr*dz)-omc2*(ep(ir,iz,3)+ep(ir,iz+1,3))/8.0D0+0.0D0*i,A_helicon)!Ez(ir,iz)
                call set(lcount,2+3*nz-3,-1.0D0/(dr*dz)-omc2*(ep(ir+1,iz-1,3)+ep(ir+1,iz,3))/8.0D0+0.0D0*i,A_helicon)!Ez(ir+1,iz-1)

                r(:)=r(:)-1.0D0/2.0D0 * dr


                lcount=3*(nz*ir+iz-nz)-1!!!Eth(ir,iz)
                call set(lcount,-1,im/r(ir)**2 *(-1.0D0/2.0D0+r(ir)/dr)-omc2*(ep(ir+1,iz,4)+ep(ir,iz,4))/4.0D0+0.0D0*i,A_helicon)!Er(ir,iz)
                call set(lcount,-1-3*nz,im/r(ir)**2 * (-1.0D0/2.0D0-r(ir)/dr)-omc2*(ep(ir,iz,4)+ep(ir-1,iz,4))/4.0D0+0.0D0*i,A_helicon)!Er(ir-1,iz)
                call set(lcount,0,1.0D0/r(ir)**2 + 2.0D0/dr22+2.0D0/dz22-omc2*ep(ir,iz,5)+0.0D0*i,A_helicon)!Eth(ir,iz)
                call set(lcount,0+3*nz,-1.0D0/(2.0D0*r(ir)*dr)-1.0D0/dr22+0.0D0*i,A_helicon)!Eth(ir+1,iz)
                call set(lcount,0-3*nz,1.0D0/(2.0D0*r(ir)*dr)-1.0D0/dr22+0.0D0*i,A_helicon)!Eth(ir-1,iz)
                call set(lcount,0+3,-1.0D0/dz22+0.0D0*i,A_helicon)!Eth(ir,iz+1)
                call set(lcount,0-3,-1.0D0/dz22+0.0D0*i,A_helicon)!Eth(ir,iz-1)
                call set(lcount,1,im/(r(ir)*dz)-omc2*(ep(ir,iz,6)+ep(ir,iz+1,6))/4.0D0+0.0D0*i,A_helicon)!Ez(ir,iz)
                call set(lcount,1-3,-im/(r(ir)*dz)-omc2*(ep(ir,iz-1,6)+ep(ir,iz,6))/4.0D0+0.0D0*i,A_helicon)!Ez(ir,iz-1)


                lcount=3*(nz*ir+iz-nz)!!!Ez(ir,iz)
                call set(lcount,-2,-1.0D0/(2.0D0*r(ir)*dz)-1.0D0/(dr*dz)-omc2*(ep(ir,iz,7)+ep(ir+1,iz,7))/8.0D0+0.0D0*i,A_helicon)!Er(ir,iz)
                call set(lcount,-2-3*nz,-1.0D0/(2.0D0*r(ir)*dz)+1.0D0/(dr*dz)-omc2*(ep(ir-1,iz,7)+ep(ir,iz,7))/8.0D0+0.0D0*i,A_helicon)!Er(ir-1,iz)
                call set(lcount,-2+3,1.0D0/(2.0D0*r(ir)*dz)+1.0D0/(dr*dz)-omc2*(ep(ir,iz+1,7)+ep(ir+1,iz+1,7))/8.0D0+0.0D0*i,A_helicon)!Er(ir,iz+1)
                call set(lcount,-2-3*nz+3,+1.0D0/(2.0D0*r(ir)*dz)-1.0D0/(dr*dz)-omc2*(ep(ir-1,iz+1,7)+ep(ir,iz+1,7))/8.0D0+0.0D0*i,A_helicon)!Er(ir-1,iz+1)
                call set(lcount,-1,-im/(r(ir)*dz)-omc2*ep(ir,iz,8)/2.0D0+0.0D0*i,A_helicon)!Eth(ir,iz)
                call set(lcount,-1+3,im/(r(ir)*dz)-omc2*ep(ir,iz+1,8)/2.0D0+0.0D0*i,A_helicon)!Eth(ir,iz+1)
                call set(lcount,0,m_real**2/(r(ir)**2)+2.0D0/dr22-omc2*(ep(ir,iz,9)+ep(ir,iz+1,9))/2.0D0+0.0D0*i,A_helicon)!Ez(ir,iz)
                call set(lcount,0+3*nz,-1.0D0/(2.0D0*r(ir)*dr)-1.0D0/dr22+0.0D0*i,A_helicon)!Ez(ir+1,iz)
                call set(lcount,0-3*nz,1.0D0/(2.0D0*r(ir)*dr)-1.0D0/dr22+0.0D0*i,A_helicon)!Ez(ir-1,iz)

                !!!Vacuum Region(nrp <= ir <= nr-1)
            elseif(Region==2 .or. Region==9 .or. Region==14)then
                lcount=3*(nz*ir+iz-nz)-2!!!Er(ir,iz)
                r(:)=r(:)+1.0D0/2.0D0 * dr

                if(ir==nr_vac(1))then
                    call set(lcount,0+3*nz,1.0D0/(2.0D0*(rp+dr))+1.0D0/dr+0.0D0*i,A_helicon)!Er(nrp+1,iz)
                    call set(lcount,0,1.0D0/(2.0D0*(rp+dr))-1.0D0/dr+0.0D0*i,A_helicon)!Er(nrp,iz)
                    call set(lcount,1,im/rp+0.0D0*i,A_helicon)!Eth(nrp+1,iz)
                    call set(lcount,2,1.0D0/dz+0.0D0*i,A_helicon)!Ez(nrp+1,iz)
                    call set(lcount,2-3,-1.0D0/dz+0.0D0*i,A_helicon)!Ez(nrp+1,iz-1)
                else
                    call set(lcount,0,(1.0D0+m_real**2)/r(ir)**2+2.0D0/dz22+2.0D0/dr22-0.5D0*omc2*(ep(ir,iz,1)+ep(ir+1,iz,1))+0.0D0*i,A_helicon)!Er(ir,iz)
                    call set(lcount,0+3,-1.0D0/dz22+0.0D0*i,A_helicon)!Er(ir,iz+1)
                    call set(lcount,0-3,-1.0D0/dz22+0.0D0*i,A_helicon)!Er(ir,iz-1)
                    call set(lcount,0+3*nz,-1.0D0/(2.0D0*r(ir)*dr)-1.0D0/dr22+0.0D0*i,A_helicon)!Er(ir+1,iz)
                    call set(lcount,0-3*nz,1.0D0/(2.0D0*r(ir)*dr)-1.0D0/dr22+0.0D0*i,A_helicon)!Er(ir-1,iz)
                    call set(lcount,1,im/r(ir)**2-0.25D0*omc2*(ep(ir,iz,2)+ep(ir+1,iz,2))+0.0D0*i,A_helicon)!Eth(ir,iz)
                    call set(lcount,1+3*nz,im/r(ir)**2-0.25D0*omc2*(ep(ir,iz,2)+ep(ir+1,iz,2))+0.0D0*i,A_helicon)!Eth(ir+1,iz)
                endif
                r(:)=r(:)-1.0D0/2.0D0 * dr

                lcount=3*(nz*ir+iz-nz)-1!!!Eth(ir,iz)
                call set(lcount,-1,-im/r(ir)**2-0.5D0*omc2*ep(ir,iz,4)+0.0D0*i,A_helicon)!Er(ir,iz)
                call set(lcount,-1-3*nz,-im/r(ir)**2-0.5D0*omc2*ep(ir,iz,4)+0.0D0*i,A_helicon)!Er(ir-1,iz)
                call set(lcount,0,(1.0D0+m_real**2)/r(ir)**2 + 2.0D0/dr22+2.0D0/dz22-omc2*ep(ir,iz,5)+0.0D0*i,A_helicon)!Eth(ir,iz)
                call set(lcount,0+3*nz,-1.0D0/(2.0D0*r(ir)*dr)-1.0D0/dr22+0.0D0*i,A_helicon)!Eth(ir+1,iz)
                call set(lcount,0-3*nz,1.0D0/(2.0D0*r(ir)*dr)-1.0D0/dr22+0.0D0*i,A_helicon)!Eth(ir-1,iz)
                call set(lcount,0+3,-1.0D0/dz22+0.0D0*i,A_helicon)!Eth(ir,iz+1)
                call set(lcount,0-3,-1.0D0/dz22+0.0D0*i,A_helicon)!Eth(ir,iz-1)

                lcount=3*(nz*ir+iz-nz)!!!Ez(ir,iz)
                if(Region==2 .or. Region==9)then
                    call set(lcount,0,m_real**2/(r(ir)**2)+2.0D0/dr22+2.0D0/dz22-0.5D0*omc2*(ep(ir,iz,9)+ep(ir,iz+1,9))+0.0D0*i,A_helicon)!Ez(ir,iz)
                    call set(lcount,0+3*nz,-1.0D0/(2.0D0*r(ir)*dr)-1.0D0/dr22+0.0D0*i,A_helicon)!Ez(ir+1,iz)
                    call set(lcount,0-3*nz,1.0D0/(2.0D0*r(ir)*dr)-1.0D0/dr22+0.0D0*i,A_helicon)!Ez(ir-1,iz)
                    call set(lcount,0+3,-1.0D0/dz22+0.0D0*i,A_helicon)!Ez(ir,iz+1)
                    call set(lcount,0-3,-1.0D0/dz22+0.0D0*i,A_helicon)!Ez(ir,iz-1)
                elseif(Region==14)then
                    call set(lcount,0,m_real**2/(r(ir)**2)+2.0D0/dr22+1.0D0/(3.0D0*dz22)-0.5D0*omc2*(ep(ir,iz,9)+ep(ir,iz+1,9))+0.0D0*i,A_helicon)!Ez(ir,iz)
                    call set(lcount,0+3*nz,-1.0D0/(2.0D0*r(ir)*dr)-1.0D0/dr22+0.0D0*i,A_helicon)!Ez(ir+1,iz)
                    call set(lcount,0-3*nz,1.0D0/(2.0D0*r(ir)*dr)-1.0D0/dr22+0.0D0*i,A_helicon)!Ez(ir-1,iz)
                    call set(lcount,0-3,-2.0D0/(3.0D0*dz22)+0.0D0*i,A_helicon)!Ez(ir,iz-1)
                    call set(lcount,0-2*3,1.0D0/(3.0D0*dz22)+0.0D0*i,A_helicon)!Ez(ir,iz-2)
                endif
                !!!!!!!!!!!!!!!!!!!! Boundary Condition !!!!!!!!!!!!!!!!!!!!!!!!!!
                !----------------------   z=0 [z=0] (iz=1)   ----------------------!
            elseif(Region==4.or.Region==12.or.Region==22)then
                lcount=3*(nz*ir+iz-nz)-2!!!Er(ir,1)
                call set(lcount,0,1.0D0+0.0D0*i,A_helicon)!Er(ir,1)
                lcount=3*(nz*ir+iz-nz)-1!!!Eth(ir,1)
                call set(lcount,0,1.0D0+0.0D0*i,A_helicon)!Eth(ir,1)
                lcount=3*(nz*ir+iz-nz)!!!Ez(ir,1)
                if(switch_boundary_half_variables==0)then
                    call set(lcount,-2,-1.0D0/(2.0D0*r(ir)*dz)-1.0D0/(dr*dz)-omc2*(ep(ir,iz,7)+ep(ir+1,iz,7))/8.0D0+0.0D0*i,A_helicon)!Er(ir,iz)
                    call set(lcount,-2-3*nz,-1.0D0/(2.0D0*r(ir)*dz)+1.0D0/(dr*dz)-omc2*(ep(ir-1,iz,7)+ep(ir,iz,7))/8.0D0+0.0D0*i,A_helicon)!Er(ir-1,iz)
                    call set(lcount,-2+3,1.0D0/(2.0D0*r(ir)*dz)+1.0D0/(dr*dz)-omc2*(ep(ir,iz+1,7)+ep(ir+1,iz+1,7))/8.0D0+0.0D0*i,A_helicon)!Er(ir,iz+1)
                    call set(lcount,-2-3*nz+3,+1.0D0/(2.0D0*r(ir)*dz)-1.0D0/(dr*dz)-omc2*(ep(ir-1,iz+1,7)+ep(ir,iz+1,7))/8.0D0+0.0D0*i,A_helicon)!Er(ir-1,iz+1)
                    call set(lcount,-1,-im/(r(ir)*dz)-omc2*ep(ir,iz,8)/2.0D0+0.0D0*i,A_helicon)!Eth(ir,iz)
                    call set(lcount,-1+3,im/(r(ir)*dz)-omc2*ep(ir,iz+1,8)/2.0D0+0.0D0*i,A_helicon)!Eth(ir,iz+1)
                    call set(lcount,0,m_real**2/(r(ir)**2)+2.0D0/dr22-omc2*(ep(ir,iz,9)+ep(ir,iz+1,9))/2.0D0+0.0D0*i,A_helicon)!Ez(ir,iz)
                    call set(lcount,0+3*nz,-1.0D0/(2.0D0*r(ir)*dr)-1.0D0/dr22+0.0D0*i,A_helicon)!Ez(ir+1,iz)
                    call set(lcount,0-3*nz,1.0D0/(2.0D0*r(ir)*dr)-1.0D0/dr22+0.0D0*i,A_helicon)!Ez(ir-1,iz)
                elseif(switch_boundary_half_variables==1)then
                    call set(lcount,0,1.0D0+0.0D0*i,A_helicon)!Ez(ir,1)
                    if(jieshu==1)then
                        call set(lcount,0+3,-1.0D0+0.0D0*i,A_helicon)!Ez(ir,2)
                    else
                        call set(lcount,0+3,-26765.0d0/9129.0d0+0.0D0*i,A_helicon)!Ez(ir,2)
                        call set(lcount,0+3*2,11630.0d0/3043.0d0+0.0D0*i,A_helicon)!Ez(ir,3)
                        call set(lcount,0+3*3,-8590.0d0/3043.0D0+0.0D0*i,A_helicon)!Ez(ir,4)
                        call set(lcount,0+3*4,10205.0d0/9129.0d0+0.0D0*i,A_helicon)!Ez(ir,5)
                        call set(lcount,0+3*5,-563.0d0/3043.0d0+0.0D0*i,A_helicon)!Ez(ir,6)
                    endif
                else
                    write(*,*)'switch_boundary_half_variables ERROR PLEASECHECK'
                    pause
                endif
            elseif(Region==7)then
                lcount=3*(nz*ir+iz-nz)-2!!!Er(ir,1)
                call set(lcount,0,1.0D0+0.0D0*i,A_helicon)!Er(ir,1)
                lcount=3*(nz*ir+iz-nz)-1!!!Eth(ir,1)
                call set(lcount,0,1.0D0+0.0D0*i,A_helicon)!Eth(ir,1)
                lcount=3*(nz*ir+iz-nz)!!!Ez(ir,1)
                if(switch_boundary_half_variables==0)then
                    call set(lcount,0,m_real**2/(r(ir)**2)+2.0D0/dr22+2.0D0/dz22-0.5D0*omc2*(ep(ir,iz,9)+ep(ir,iz+1,9))+0.0D0*i,A_helicon)!Ez(ir,iz)
                    call set(lcount,0+3*nz,-1.0D0/(2.0D0*r(ir)*dr)-1.0D0/dr22+0.0D0*i,A_helicon)!Ez(ir+1,iz)
                    call set(lcount,0-3*nz,1.0D0/(2.0D0*r(ir)*dr)-1.0D0/dr22+0.0D0*i,A_helicon)!Ez(ir-1,iz)
                    call set(lcount,0+3,-1.0D0/dz22+0.0D0*i,A_helicon)!Ez(ir,iz+1)
                    call set(lcount,0-3,-1.0D0/dz22+0.0D0*i,A_helicon)!Ez(ir,iz-1)
                elseif(switch_boundary_half_variables==1)then
                    call set(lcount,0,1.0D0+0.0D0*i,A_helicon)!Ez(ir,1)
                    if(jieshu==1)then
                        call set(lcount,0+3,-1.0D0+0.0D0*i,A_helicon)!Ez(ir,2)
                    else
                        call set(lcount,0+3,-26765.0d0/9129.0d0+0.0D0*i,A_helicon)!Ez(ir,2)
                        call set(lcount,0+3*2,11630.0d0/3043.0d0+0.0D0*i,A_helicon)!Ez(ir,3)
                        call set(lcount,0+3*3,-8590.0d0/3043.0D0+0.0D0*i,A_helicon)!Ez(ir,4)
                        call set(lcount,0+3*4,10205.0d0/9129.0d0+0.0D0*i,A_helicon)!Ez(ir,5)
                        call set(lcount,0+3*5,-563.0d0/3043.0d0+0.0D0*i,A_helicon)!Ez(ir,6)
                    endif
                else
                    write(*,*)'switch_boundary_half_variables ERROR PLEASECHECK'
                    pause
                endif
                !----------------------   z=0 [z=0] (iz=1)   ----------------------!



                !----------------------   z=zl [z=l] (iz=nz)   ----------------------!(����)
            elseif(Region==5 .or. Region==8 .or. Region==13.or.Region==24)then
                if(switch_z_boundary==3)then
                    call set(3*nz*ir-2,0,23.0D0/24.0D0-i*dz/2.0D0*sqrt((om/c)**2-(DBesselJZero(m)/rl)**2)+0.0D0*i,A_helicon)!!!Er(ir,nz)
                    call set(3*nz*ir-2,-3,-7.0D0/8.0D0-i*dz/2.0D0*sqrt((om/c)**2-(DBesselJZero(m)/rl)**2)+0.0D0*i,A_helicon)!Er(ir,nz-1)
                    call set(3*nz*ir-2,-2*3,-1.0D0/8.0D0+0.0D0*i,A_helicon)!Er(ir,nz-2)
                    call set(3*nz*ir-2,-3*3,1.0D0/24.0D0+0.0D0*i,A_helicon)!Er(ir,nz-3)

                    call set(3*nz*ir-1,0,11.0D0/6.0D0-i*dz*sqrt((om/c)**2-(DBesselJZero(m)/rl)**2)+0.0D0*i,A_helicon)!!!Eth(ir,nz)
                    call set(3*nz*ir-1,-3,-3.0D0+0.0D0*i,A_helicon)!Eth(ir,nz-1)
                    call set(3*nz*ir-1,-2*3,3.0D0/2.0D0+0.0D0*i,A_helicon)!Eth(ir,nz-2)
                    call set(3*nz*ir-1,-3*3,-1.0D0/3.0D0+0.0D0*i,A_helicon)!Eth(ir,nz-3)
                else
                    lcount=3*(nz*ir+iz-nz)-2!!!Er(ir,nz)
                    call set(lcount,0,1.0D0+0.0D0*i,A_helicon)!Er(ir,nz)
                    lcount=3*(nz*ir+iz-nz)-1!!!Eth(ir,nz)
                    call set(lcount,0,1.0D0+0.0D0*i,A_helicon)!Eth(ir,nz)
                endif

                lcount=3*(nz*ir+iz-nz)!!!Ez(ir,nz)
                if(switch_boundary_half_variables==0)then
                    call set(lcount,0,1.0D0+0.0D0*i,A_helicon)!Ez(ir,nz)
                    call set(lcount,0-3,-2.0D0+0.0D0*i,A_helicon)!Ez(ir,nz-1)
                    call set(lcount,0-2*3,1.0D0+0.0D0*i,A_helicon)!Ez(ir,nz-2)
                elseif(switch_boundary_half_variables==1)then
                    call set(lcount,0,1.0D0+0.0D0*i,A_helicon)!Ez(ir,nz)
                    if(jieshu==1)then
                        call set(lcount,0-3,-1.0D0+0.0D0*i,A_helicon)!Ez(ir,nz-1)
                    else
                        call set(lcount,0-3,-335.0d0/563.0d0+0.0D0*i,A_helicon)!Ez(ir,nz-1)
                        call set(lcount,0-3*2,-1430.0d0/1689.0d0+0.0D0*i,A_helicon)!Ez(ir,nz-2)
                        call set(lcount,0-3*3,370.0d0/563.0d0+0.0D0*i,A_helicon)!Ez(ir,nz-3)
                        call set(lcount,0-3*4,-145.0d0/563.0d0+0.0D0*i,A_helicon)!Ez(ir,nz-4)
                        call set(lcount,0-3*5,71.0d0/1689.0d0+0.0D0*i,A_helicon)!Ez(ir,nz-5)
                    endif
                endif
                !----------------------   z=zl [z=l] (iz=nz)   ----------------------!



                !----------------------   r=rl [r=R] (ir=nr)   ----------------------!(����)
            elseif(Region==3.or.Region==6.or.Region==21.or.Region==25.or.Region==26)then
                lcount=3*(nz*ir+iz-nz)-2!!!Er(nr,iz)
                if(switch_boundary_half_variables==0)then
                    call set(lcount,0,1.0D0+0.0D0*i,A_helicon)!Er(nr,iz)
                    call set(lcount,0-3*nz,-2.0D0+0.0D0*i,A_helicon)!Er(nr-1,iz)
                    call set(lcount,0-2*3*nz,1.0D0+0.0D0*i,A_helicon)!Er(nr-2,iz)
                elseif(switch_boundary_half_variables==1)then
                    call set(lcount,0,1.0D0/(2.0D0*rl) + 1.0D0/dr+0.0D0*i,A_helicon) !Er(nr,iz)
                    call set(lcount,0-3*nz,1.0D0/(2.0D0*rl) - 1.0D0/dr+0.0D0*i,A_helicon) !Er(nr-1,iz)
                endif

                lcount=3*(nz*ir+iz-nz)-1!!!Eth(nr,iz)
                call set(lcount,0,1.0D0+0.0D0*i,A_helicon)!Eth(nr,iz)
                lcount=3*(nz*ir+iz-nz)!!!Ez(nr,iz)
                call set(lcount,0,1.0D0+0.0D0*i,A_helicon)!Ez(nr,iz)
                !----------------------   r=rl [r=R] (ir=nr)   ----------------------!



                !----------------------   r=0 [r=0] (ir=1)   ----------------------!
            elseif(Region==10)then
                if(m==0)then
                    lcount=3*(nz*ir+iz-nz)-2!!!Er(1,iz)
                    if(switch_boundary_half_variables==1)then
                        call set(lcount,0,3.0D0/2.0D0+0.0D0*i,A_helicon)!Er(1,iz)
                        call set(lcount,0+3*nz,-1.0D0/2.0D0+0.0D0*i,A_helicon)!Er(2,iz)
                    endif
                    lcount=3*(nz*ir+iz-nz)-1!!!Eth(1,iz)
                    call set(lcount,0,1.0D0+0.0D0*i,A_helicon)!Eth(1,iz)
                    lcount=3*(nz*ir+iz-nz)!!!Ez(1,iz)
                    call set(lcount,0,1.0D0+0.0D0*i,A_helicon)!Ez(1,iz)
                    call set(lcount,0+3*nz,-1.0D0+0.0D0*i,A_helicon)!Ez(2,iz)
                elseif(m==-1.or.m==1)then
                    lcount=3*(nz*ir+iz-nz)-2!!!Er(1,iz)
                    if(switch_r_0_boundary==2)then
                        if(iz==1)then
                            call set(lcount,0,1.0D0+0.0D0*i,A_helicon)!Er(1,1)
                        else
                            r(:)=r(:)+1.0D0/2.0D0 * dr
                            call set(lcount,0,m_real**2/r(ir)**2+2.0D0/dz22-0.5D0*omc2*(ep(ir,iz,1)+ep(ir+1,iz,1))+0.0D0*i,A_helicon)!Er(ir,iz)
                            call set(lcount,0+3,-1.0D0/dz22+0.0D0*i,A_helicon)!Er(ir,iz+1)
                            call set(lcount,0-3,-1.0D0/dz22+0.0D0*i,A_helicon)!Er(ir,iz-1)
                            call set(lcount,1,im/r(ir)**2 *(1.0D0/2.0D0-r(ir)/dr)-0.25D0*omc2*(ep(ir,iz,2)+ep(ir+1,iz,2))+0.0D0*i,A_helicon)!Eth(ir,iz)
                            call set(lcount,1+3*nz,im/r(ir)**2 *(1.0D0/2.0D0+r(ir)/dr)-0.25D0*omc2*(ep(ir,iz,2)+ep(ir+1,iz,2))+0.0D0*i,A_helicon)!Eth(ir+1,iz)
                            call set(lcount,2+3*nz,1.0D0/(dr*dz)+0.0D0*i,A_helicon)!Ez(ir+1,iz)
                            call set(lcount,2-3,1.0D0/(dr*dz)+0.0D0*i,A_helicon)!Ez(ir,iz-1)
                            call set(lcount,2,-1.0D0/(dr*dz)+0.0D0*i,A_helicon)!Ez(ir,iz)
                            call set(lcount,2+3*nz-3,-1.0D0/(dr*dz)+0.0D0*i,A_helicon)!Ez(ir+1,iz-1)
                            r(:)=r(:)-1.0D0/2.0D0 * dr
                        endif
                    elseif(switch_r_0_boundary==1)then
                        if(switch_boundary_half_variables==1)then
                            call set(lcount,0,1.0D0+0.0D0*i,A_helicon)!Er(1,iz)
                            call set(lcount,0+3*nz,-1.0D0+0.0D0*i,A_helicon)!Er(2,iz)
                        endif
                    else
                        write(*,*)'switch_r_0_boundary ERROR please check(=1 or 2)'
                    endif
                    lcount=3*(nz*ir+iz-nz)-1!!!Eth(1,iz)
                    if(switch_r_0_boundary==2)then
                        call set(lcount,0,m_real*i*(1.0D0+0.0D0*i),A_helicon)!Eth(1,iz)
                        call set(lcount,-1,3.0D0/2.0D0+0.0D0*i,A_helicon)!Er(1,iz)
                        call set(lcount,-1+3*nz,-1.0D0/2.0D0+0.0D0*i,A_helicon)!Er(2,iz)
                    elseif(switch_r_0_boundary==1)then
                        call set(lcount,0,1.0D0+0.0D0*i,A_helicon)!Eth(1,iz)
                        call set(lcount,0+3*nz,-1.0D0+0.0D0*i,A_helicon)!Eth(2,iz)
                    endif
                    lcount=3*(nz*ir+iz-nz)!!!Ez(1,iz)
                    call set(lcount,0,1.0D0+0.0D0*i,A_helicon)!Ez(1,iz)
                else
                    lcount=3*(nz*ir+iz-nz)-2!!!Er(1,iz)
                    if(switch_boundary_half_variables==1)then
                        call set(lcount,0,3.0D0/2.0D0+0.0D0*i,A_helicon)!Er(1,iz)
                        call set(lcount,0+3*nz,-1.0D0/2.0D0+0.0D0*i,A_helicon)!Er(2,iz)
                    endif
                    lcount=3*(nz*ir+iz-nz)-1!!!Eth(1,iz)
                    call set(lcount,0,1.0D0+0.0D0*i,A_helicon)!Eth(1,iz)
                    lcount=3*(nz*ir+iz-nz)!!!Ez(1,iz)
                    call set(lcount,0,1.0D0+0.0D0*i,A_helicon)!Ez(1,iz)
                endif
            elseif(Region==11)then
                lcount=3*(nz*ir+iz-nz)-2!!!Er(ir,iz)
                call set(lcount,0,1.0D0+0.0D0*i,A_helicon)!Er(ir,iz)
                lcount=3*(nz*ir+iz-nz)-1!!!Eth(ir,iz)
                call set(lcount,0,1.0D0+0.0D0*i,A_helicon)!Erh(ir,iz)
                lcount=3*(nz*ir+iz-nz)!!!Ez(ir,iz)
                call set(lcount,0,1.0D0+0.0D0*i,A_helicon)!Ez(ir,iz)

                !----------------------   dipole upper   ----------------------!(����)
            elseif(Region==23)then
                lcount=3*(nz*ir+iz-nz)-2!!!Er(nr,iz)
                r(:)=r(:)+1.0D0/2.0D0 * dr
                call set(lcount,0,m_real**2/r(ir)**2+2.0D0/dz22-omc2*(ep(ir+1,iz,1)+ep(ir,iz,1))/2.0D0+0.0D0*i,A_helicon)!Er(ir,iz)
                call set(lcount,0+3,-1.0D0/dz22+0.0D0*i,A_helicon)!Er(ir,iz+1)
                call set(lcount,0-3,-1.0D0/dz22+0.0D0*i,A_helicon)!Er(ir,iz-1)
                call set(lcount,1,im/r(ir)**2 *(1.0D0/2.0D0-r(ir)/dr)-omc2*ep(ir,iz,2)/2.0D0+0.0D0*i,A_helicon)!Eth(ir,iz)
                call set(lcount,1+3*nz,im/r(ir)**2 *(1.0D0/2.0D0+r(ir)/dr)-omc2*ep(ir+1,iz,2)/2.0D0+0.0D0*i,A_helicon)!Eth(ir+1,iz)
                call set(lcount,2+3*nz,1.0D0/(dr*dz)-omc2*(ep(ir+1,iz,3)+ep(ir+1,iz+1,3))/8.0D0+0.0D0*i,A_helicon)!Ez(ir+1,iz)
                call set(lcount,2-3,1.0D0/(dr*dz)-omc2*(ep(ir,iz-1,3)+ep(ir+1,iz,3))/8.0D0+0.0D0*i,A_helicon)!Ez(ir,iz-1)
                call set(lcount,2,-1.0D0/(dr*dz)-omc2*(ep(ir,iz,3)+ep(ir,iz+1,3))/8.0D0+0.0D0*i,A_helicon)!Ez(ir,iz)
                call set(lcount,2+3*nz-3,-1.0D0/(dr*dz)-omc2*(ep(ir+1,iz-1,3)+ep(ir+1,iz,3))/8.0D0+0.0D0*i,A_helicon)!Ez(ir+1,iz-1)
                r(:)=r(:)-1.0D0/2.0D0 * dr
                lcount=3*(nz*ir+iz-nz)-1!!!Eth(nr,iz)
                call set(lcount,0,1.0D0+0.0D0*i,A_helicon)!Eth(nr,iz)
                lcount=3*(nz*ir+iz-nz)!!!Ez(nr,iz)
                call set(lcount,0,1.0D0+0.0D0*i,A_helicon)!Ez(nr,iz)
                !----------------------   dipole upper   ----------------------!

            endif
        enddo
    enddo
    !pause

    if(iswitch_output_Region==1) close(545)!pause;!TEST
    CALL findnzmax(A_helicon)
    !pause
    end subroutine set_A
    !--------------------subroutine set-----------------!
    subroutine set(l1,sml1,Als,A1)
    use Types
    implicit none
    integer*4 :: l1,sml1
    Complex*16 :: Als
    type(sparse_complex) :: A1
    CALL setTriplet(l1,l1+sml1,Als,A1)
    end subroutine set
    !--------------------subroutine setTriplet-----------------!
    subroutine setTriplet(icount,jcount,value,A)
    use Types
    implicit none
    type(sparse_complex) :: A
    integer*4 :: icount,jcount
    Complex*16 :: value
    !type(sparse_complex) :: A(:)

    A.nzmax=A.nzmax+1
    A.row(A.nzmax)=icount
    A.col(A.nzmax)=jcount
    A.val(A.nzmax)=value

    end subroutine setTriplet
    !--------------------subroutine findnzmax-----------------!
    subroutine findnzmax(A)
    use types
    implicit none
    type(sparse_complex) :: A
    integer*4 :: icount2
    Do icount2=1,A.n**2
        if(icount2>size(A.col))then
            A.nzmax=icount2-1
            EXIT
        elseif(A.col(icount2)==0)then
            A.nzmax=icount2-1
            EXIT
        end if
    End Do
    continue
    !write(*,*)A.nzmax;pause!TEST!!!FOR TEST ONLY
    end subroutine findnzmax
    !--------------------recursive subroutine quicksort-----------------!
    recursive subroutine quicksort(A1, A2, A3, first, last)
    implicit none
    integer*4 ::  A1(*),A2(*), x1, x2, t1, t2
    Complex*16 :: A3(*),t3
    integer*4 :: first, last
    integer*4 :: i00, j00
    ! quicksort.f -*-f90-*-
    ! Author: t-nissie
    ! License: GPLv3
    ! Gist: https://gist.github.com/t-nissie/479f0f16966925fa29ea
    !!
    x1 = A1( (first+last) / 2 )
    x2 = A2( (first+last) / 2 )
    i00 = first
    j00 = last
    do
        do while ((A1(i00) < x1) .or. ((A1(i00)==x1).and.(A2(i00)<x2)))
            i00=i00+1
        end do
        do while ((A1(j00) > x1) .or. ((A1(j00)==x1).and.(A2(j00)>x2)))
            j00=j00-1
        end do
        if (i00 >= j00) exit
        t1 = A1(i00);  A1(i00) = A1(j00);  A1(j00) = t1
        t2 = A2(i00);  A2(i00) = A2(j00);  A2(j00) = t2
        t3 = A3(i00);  A3(i00) = A3(j00);  A3(j00) = t3
        i00=i00+1
        j00=j00-1
    end do
    if (first < i00-1) call quicksort(A1,A2,A3, first, i00-1)
    if (j00+1 < last)  call quicksort(A1,A2,A3, j00+1, last)
    end subroutine quicksort

    !--------------------subroutine complex_coo_pardiso_solve-----------------!
    subroutine complex_coo_pardiso_solve(A,x_plasma,B_rf)
    use types
    IMPLICIT NONE
    type(sparse_complex) :: A
    Complex*16 :: x_plasma(A.n),B_rf(A.n)
    type(sparse_complex) :: A_csr3
    integer*4 :: pt(64),maxfct,mnum,mtype,phase,msglvl,nrhs,iparm(64),error
    integer*4,allocatable :: perm(:)
    allocate (A_csr3.row(1:A.n+1),A_csr3.col(1:A.nzmax),A_csr3.val(1:A.nzmax))
    allocate (perm(1:A.n))
    CALL findnzmax(A)
    A_csr3.n=A.n;A_csr3.nzmax=A.nzmax
    CALL coo_to_csr(A,A_csr3)
    !maxfct=1;mnum=1;mtype=13;phase=13;nrhs=1;msglvl=0
    !CALL pardisoinit(pt, mtype, iparm)
    maxfct=1;mnum=1;mtype=13;phase=13;nrhs=1;iparm(:)=0;msglvl=0
    iparm(3)=4 !OMP_NUM_THREADS
    CALL pardiso(pt,maxfct,mnum,mtype,phase,A_csr3.n,A_csr3.val,A_csr3.row,A_csr3.col,perm,nrhs,iparm,msglvl,B_rf,x_plasma,error)
    !CALL mkl_free_buffers
    deallocate(perm)
    deallocate(A_csr3.row,A_csr3.col,A_csr3.val)
    end subroutine complex_coo_pardiso_solve
    !--------------------subroutine coo_to_csr-----------------!
    subroutine coo_to_csr(A,A_csr3)
    use types
    implicit none
    type(sparse_complex) :: A,A_csr3
    integer*4 :: col1,row1,rowcount,ncount,count1
    Complex*16 :: val1
    integer*4 :: rownum(1:A.n+1)
    rownum=0
    CALL quicksort(A.row,A.col,A.val,1,A.nzmax)
    DO rowcount=1,A.nzmax
        rownum(A.row(rowcount))=rownum(A.row(rowcount))+1
    End DO
    A_csr3.row(1)=1
    do rowcount=1,A.n
        A_csr3.row(rowcount+1)=A_csr3.row(rowcount)+rownum(rowcount)
    end do
    A_csr3.col(:)=A.col(:)
    A_csr3.val(:)=A.val(:)
    end subroutine coo_to_csr
    !--------------------subroutine reallocate-----------------!
    subroutine reallocate(A,nzmax_input)
    use Types
    implicit none
    type(sparse_complex) :: A
    integer :: nzmax_input
    integer*4,allocatable :: Arow(:),Acol(:)
    Complex*16,allocatable :: Aval(:)
    CALL findnzmax(A)
    allocate(Arow(1:A.nzmax),Acol(1:A.nzmax),Aval(1:A.nzmax))
    Arow(1:A.nzmax)=A.row(1:A.nzmax)
    Acol(1:A.nzmax)=A.col(1:A.nzmax)
    Aval(1:A.nzmax)=A.val(1:A.nzmax)
    Deallocate(A.row,A.col,A.val)
    Allocate(A.row(1:nzmax_input),A.col(1:nzmax_input),A.val(1:nzmax_input))
    A.row(1:A.nzmax)=Arow(1:A.nzmax)
    A.col(1:A.nzmax)=Acol(1:A.nzmax)
    A.val(1:A.nzmax)=Aval(1:A.nzmax)
    A.row(A.nzmax+1:nzmax_input)=0
    A.col(A.nzmax+1:nzmax_input)=0
    A.val(A.nzmax+1:nzmax_input)=0.0D0+0.0*(0.0D0,1.0D0)
    Deallocate(Acol,Arow,Aval)
    CALL findnzmax(A)
    End subroutine reallocate



    subroutine find_fvm
    use the_whole_varibles
    implicit none
    real*8 :: tev,vte,ane,vei,ve,ven,ug
    ug=10.9D0; !ug Coulomb logarithm;vte is the electron thermal velocity
    fvm=0.0D0


    do iz=1,nz
        do ir=1,nrp !test!!!!!!!!!!!!!
            tev=te(ir,iz)*qe;
            if(te(ir,iz)<fv_te_limit)tev=fv_te_limit*qe;
            vte=sqrt(8.0D0*tev/(pi*me));ane=ne(ir,iz);
            ven=nn*sig_en*vte;
            vei=(tev/qe)**(-1.5D0)*ug*(2.91D-12)*ane;
            ve=ven+vei;
            fvm(ir,iz)=ve;
        end do
    end do
    continue
    endsubroutine find_fvm



    subroutine find_ve_vi(ve,vi)
    use the_whole_varibles
    implicit none
    real*8 tev,vte,ne_tp,vei,ve,ven,ug
    real*8 ti_tp,vti,vie,vin,vi
    ug=10.9D0; !ug Coulomb logarithm;vte is the electron thermal velocity

    ne_tp=ne(ir,iz);
    tev=te(ir,iz)*qe;
    if(te(ir,iz)<fv_te_limit)tev=fv_te_limit*qe;
    vte=sqrt(8.0D0*tev/(pi*me));
    ven=nn*sig_en*vte;
    vei=(tev/qe)**(-1.5D0)*ug*(2.91D-12)*ne_tp;
    ve=ven+vei;

    ti_tp=ek_ion_2D(ir,iz)/1.5d0 !eV
    if(ti_tp<1.)ti_tp=1.0d0
    vti=sqrt(8.d0*qe*ti_tp/(pi*mi));
    vie=4.80d-14*mass_number**(-0.5)*ti_tp**(-1.5)*ug*ne_tp;
    vin=nn*sig_in*vti;  !vi=500 m/s, 500K
    vi=vie+vin;
    endsubroutine find_ve_vi


    !--------------------subroutine totalsolvD-----------------!
    subroutine totalsolve
    use the_whole_varibles
    implicit none
    integer*4 :: kcount3,ndirec1,ndirec2
    real*8 :: ave_p
    ave_p=0.5D0

    DO kcount3=1, nr*nz
        e_m(((kcount3-1)/nz+1),(MOD(kcount3-1,nz)+1),1)=eq_xlast(m,3*kcount3-2)
        e_m((kcount3-1)/nz+1,MOD(kcount3-1,nz)+1,2)=eq_xlast(m,3*kcount3-1)
        e_m((kcount3-1)/nz+1,MOD(kcount3-1,nz)+1,3)=eq_xlast(m,3*kcount3)
    End DO
    !Do ir=2,nr
    !e_int(ir,:,1)=0.5D0*(e_m(ir,:,1)+e_m(ir-1,:,1))
    !EndDo
    !e_int(1,:,1)=1.5D0*e_m(1,:,1)-0.5D0*e_m(2,:,1)
    !
    !e_int(:,:,2)=e_m(:,:,2)
    !Do iz=2,nz
    !e_int(:,iz,3)=0.5D0*(e_m(:,iz,3)+e_m(:,iz-1,3))
    !EndDo
    !e_int(:,1,3)=1.5D0*e_m(:,1,3)-0.5D0*e_m(:,2,3)

    e_int(2:nr,1:nz,1)=(0.5D0+0.0D0*i)*(e_m(2:nr,1:nz,1)+e_m(1:nr-1,1:nz,1))
    e_int(1,1:nz,1)=e_m(1,1:nz,1)

    e_int(1:nr,1:nz,2)=e_m(1:nr,1:nz,2)

    e_int(1:nr,2:nz,3)=(0.5D0+0.0D0*i)*(e_m(1:nr,2:nz,3)+e_m(1:nr,1:nz-1,3))
    e_int(1:nr,1,3)=e_m(1:nr,1,3)

    !CALL finddive
    CALL findb
    CALL findjp
    CALL ptotalsolve
    do ith=1,nth
        e_record(:,ith,:,:)=e_int(:,:,:)*exp(im*th(ith))
        b_record(:,ith,:,:)=b_m(:,:,:)*exp(im*th(ith))
        !    dive_record(:,ith,:)=divE_m(:,:)*exp(im*th(ith))
        ja_record(:,ith,:,:)=ja_m_nor(:,:,:)*exp(im*th(ith))
        jp_record(:,ith,:,:)=jp_m(:,:,:)*exp(im*th(ith))
    enddo

    e_AC=e_AC+e_record
    b_AC=b_AC+b_record
    !dive_AC=dive_AC+dive_record
    ja_AC=ja_AC+ja_record
    jp_AC=jp_AC+jp_record


    !ptotm(:,:,:)=ptotm(:,:,:)+real(ave_p*e_int(:,:,:)*DCONJG(jp_m(:,:,:)))

    e_output(m,:,:,:)=e_int(:,:,:)
    max_eth(m)=maxval(abs(e_output(m,1:nrp,:,2)))
    !jp_output(m,:,:,:)=jp_m(:,:,:)
    !ptot_output(m,:,:)=real(peff(:,:))

    endsubroutine totalsolve


    subroutine findb
    use the_whole_varibles
    implicit none

    do ir=1,nr-1
        do iz=1,nz-1
            if(ir/=1)then
                b_m(ir,iz,1)=im/r(ir)*e_m(ir,iz,3)-1.0D0/dz*e_m(ir,iz+1,2)+1.0D0/dz*e_m(ir,iz,2)
            endif
            b_m(ir,iz,2)=1.0D0/dz*e_m(ir,iz+1,1)-1.0D0/dz*e_m(ir,iz,1)-1.0D0/dr*e_m(ir+1,iz,3)+1.0D0/dr*e_m(ir,iz,3)
            b_m(ir,iz,3)=1.0D0/(dr*(0.5D0*dr+r(ir)))*(r(ir+1)*e_m(ir+1,iz,2)-r(ir)*e_m(ir,iz,2))-im*e_m(ir,iz,1)/(r(ir)+0.5D0*dr)
        enddo
        if(m/=0)then
            b_m(1,:,1)=b_m(2,:,1)
        else
            b_m(1,:,1)=0.0D0+0.0D0*i
        endif
    enddo
    iz=nz
    do ir=1,nr-1
        b_m(ir,iz,1)=2.0D0*b_m(ir,iz-1,1)-b_m(ir,iz-2,1)
        b_m(ir,iz,2)=2.0D0*b_m(ir,iz-1,2)-b_m(ir,iz-2,2)
        b_m(ir,iz,3)=1.0D0/(dr*(0.5D0*dr+r(ir)))*(r(ir+1)*e_m(ir+1,iz,2)-r(ir)*e_m(ir,iz,2))-im*e_m(ir,iz,1)/(r(ir)+0.5D0*dr)
    enddo
    ir=nr
    do iz=1,nz-1
        b_m(ir,iz,1)=im/r(ir)*e_m(ir,iz,3)-1.0D0/dz*e_m(ir,iz+1,2)+1.0D0/dz*e_m(ir,iz,2)
        b_m(ir,iz,2)=2.0D0*b_m(ir-1,iz,2)-b_m(ir-2,iz,2)
        b_m(ir,iz,3)=2.0D0*b_m(ir-1,iz,3)-b_m(ir-2,iz,3)
    enddo
    ir=nr;iz=nz
    b_m(ir,iz,:)=-b_m(ir-1,iz-1,:)+b_m(ir,iz-1,:)+b_m(ir-1,iz,:)

    b_m(:,:,:)=b_m(:,:,:)/(i*om)
    end subroutine findb
    subroutine findjp
    use the_whole_varibles
    implicit none
    integer*4 :: ndirec1,ndirec2
    jp_m=0.0D0
    Do ndirec1=1,3
        Do ndirec2=1,3
            jp_m(:,:,ndirec1)=jp_m(:,:,ndirec1)+si(:,:,ndirec1,ndirec2)*e_int(:,:,ndirec2)
        EndDo
    EndDo
    end subroutine findjp



    subroutine ptotalsolve
    use the_whole_varibles
    implicit none
    integer*4 :: ndirec1,ndirec2,n_boundary
    real*8 :: ave_p,ptotal_pt
    n_boundary=nint(0.8*nrp)
    ave_p=0.5D0
    peff(:,:)=0.0D0
    do ndirec1=1,3
        do ndirec2=1,3
            power_depo(:,:,ndirec2)=power_depo(:,:,ndirec2)+ave_p*2.0D0*pi*(si(:,:,ndirec2,ndirec1)*e_int(:,:,ndirec1)*DCONJG(e_int(:,:,ndirec2)))
            !peff(:,:)=peff(:,:)+ave_p*2.0D0*pi*(si(:,:,ndirec2,ndirec1)*e_int(:,:,ndirec1)*DCONJG(e_int(:,:,ndirec2)))
        end do
    end do
    
    peff(:,:)=power_depo(:,:,1)+power_depo(:,:,2)+power_depo(:,:,3)
    
    
    ptotal_pt=ptotal
    do ir=1,nrp-1
        do iz=1,nz-1
            ptotal=ptotal+ds/4.0D0*real(peff(ir,iz)*r(ir)+peff(ir+1,iz)*r(ir+1)&
                &+peff(ir,iz+1)*r(ir)+peff(ir+1,iz+1)*r(ir+1))
        end do
    end do
    do ir=1,n_boundary-1
        do iz=1,nz-1
            ptotal_inside=ptotal_inside+ds/4.0D0*real(peff(ir,iz)*r(ir)+peff(ir+1,iz)*r(ir+1)&
                &+peff(ir,iz+1)*r(ir)+peff(ir+1,iz+1)*r(ir+1))
        end do
    end do
    do ir=n_boundary,nrp-1
        do iz=1,nz-1
            ptotal_boundary=ptotal_boundary+ds/4.0D0*real(peff(ir,iz)*r(ir)+peff(ir+1,iz)*r(ir+1)&
                &+peff(ir,iz+1)*r(ir)+peff(ir+1,iz+1)*r(ir+1))
        end do
    end do
    ptotal_m(m)=ptotal-ptotal_pt
    
    !do ir=1,nr
    !ptotal_BI=ptotal_BI+alpha_BI*(abs(B_m(ir,1,1))**2+abs(B_m(ir,1,2))**2&
    !    &+abs(B_m(ir,nz,1))**2+abs(B_m(ir,nz,2))**2)*r(ir)
    !enddo
    end subroutine ptotalsolve


    subroutine Find_Region(nr,nz,n_vac,nr_vac,nr_met,nz_vac,nrz_diploe,iswtich_inner_dipole,FindRegion)
    implicit none
    integer*4 :: nr,nz,ir,iz,nn,n_vac
    integer*4 :: nr_vac(n_vac),nz_vac(n_vac),nr_met(n_vac),nrz_diploe(n_vac),iswtich_inner_dipole
    integer*4 :: FindRegion(nr,nz),F1(nr,0:nz+1)
    integer*4 :: nz_vac1(0:n_vac)
    nz_vac1(0)=1;nz_vac1(1:n_vac)=nz_vac(:);

    FindRegion(:,:)=0;  F1(:,:)=0
    F1(:,0)=11; F1(:,nz+1)=11


    do nn=1,n_vac
        do iz=nz_vac1(nn-1),nz_vac1(nn)
            F1(1,iz)=10
            F1(2:nr_vac(nn)-1,iz)=1
            F1(nr_vac(nn),iz)=9
            F1(nr_vac(nn)+1:nr_met(nn)-1,iz)=2
            F1(nr_met(nn),iz)=6
            F1(nr_met(nn)+1:nr,iz)=11
            if(nr_vac(nn)==1)F1(1,iz)=10
            if(nr_met(nn)==nr)F1(nr,iz)=6
            if(nr_vac(nn)==nr_met(nn))F1(nr_vac(nn),iz)=3
        enddo
    enddo

    do nn=0,n_vac
        iz=nz_vac1(nn)
        do ir=2,nr

            if(F1(ir,iz-1)==1)then
                if(F1(ir,iz+1)==1)F1(ir,iz)=1
                if(F1(ir,iz+1)==9)F1(ir,iz)=7
                if(F1(ir,iz+1)==2)F1(ir,iz)=7
                if(F1(ir,iz+1)==6)F1(ir,iz)=6
                if(F1(ir,iz+1)==11)F1(ir,iz)=5
                if(F1(ir,iz+1)==3)F1(ir,iz)=13
            elseif(F1(ir,iz-1)==9)then
                if(F1(ir,iz+1)/=9)F1(ir,iz-1)=14!!!!!!TEST: Can it be deleted?
                if(F1(ir,iz+1)==1)F1(ir,iz)=4
                if(F1(ir,iz+1)==9)F1(ir,iz)=9
                if(F1(ir,iz+1)==2)F1(ir,iz)=7
                if(F1(ir,iz+1)==6)F1(ir,iz)=6
                if(F1(ir,iz+1)==11)F1(ir,iz)=11
                if(F1(ir,iz+1)==3)F1(ir,iz)=13
            elseif(F1(ir,iz-1)==2)then
                if(F1(ir,iz+1)/=2)F1(ir,iz-1)=14!!!!!!TEST: Can it be deleted?
                if(F1(ir,iz+1)==1)F1(ir,iz)=4
                if(F1(ir,iz+1)==9)F1(ir,iz)=7
                if(F1(ir,iz+1)==2)F1(ir,iz)=2
                if(F1(ir,iz+1)==6)F1(ir,iz)=6
                if(F1(ir,iz+1)==11)F1(ir,iz)=11
                if(F1(ir,iz+1)==3)F1(ir,iz)=13
            elseif(F1(ir,iz-1)==6)then
                if(F1(ir,iz+1)==1)F1(ir,iz)=4
                if(F1(ir,iz+1)==9)F1(ir,iz)=4
                if(F1(ir,iz+1)==2)F1(ir,iz)=7
                if(F1(ir,iz+1)==6)F1(ir,iz)=6
                if(F1(ir,iz+1)==11)F1(ir,iz)=11
                if(F1(ir,iz+1)==3)F1(ir,iz)=3
            elseif(F1(ir,iz-1)==11)then
                if(F1(ir,iz+1)==1)F1(ir,iz)=4
                if(F1(ir,iz+1)==9)F1(ir,iz)=4
                if(F1(ir,iz+1)==2)F1(ir,iz)=7
                if(F1(ir,iz+1)==6)F1(ir,iz)=6
                if(F1(ir,iz+1)==11)F1(ir,iz)=11
                if(F1(ir,iz+1)==3)F1(ir,iz)=3
            elseif(F1(ir,iz-1)==3)then
                if(F1(ir,iz+1)==1)F1(ir,iz)=12
                if(F1(ir,iz+1)==9)F1(ir,iz)=4
                if(F1(ir,iz+1)==2)F1(ir,iz)=12
                if(F1(ir,iz+1)==6)F1(ir,iz)=6
                if(F1(ir,iz+1)==11)F1(ir,iz)=3
                if(F1(ir,iz+1)==3)F1(ir,iz)=3
            endif
        enddo
    enddo

    FindRegion(:,:)=F1(1:nr,1:nz)
    if(minval(FindRegion(:,:))==0)then
        write(*,*)'FindRegion ERROR, Please Check'
        pause
    endif
    endsubroutine Find_Region





















