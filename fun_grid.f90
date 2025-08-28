    subroutine grid
    use the_whole_varibles
    implicit none
    integer::ir1,ir2,ir3,iz1,iz2,i1,i2,i4
    real*8, allocatable :: r_vac(:),r_met(:) ,z_vac(:) !met-metal

    !nr=31
    !nz=51

    nr=101
    nz=101

    !Device configuration
    !in meter
    n_vac=1; !Number of vacuum cylinders
    allocate (r_vac(n_vac),r_met(n_vac),z_vac(n_vac))
    allocate (nr_vac(n_vac),nz_vac(n_vac),nr_met(n_vac),nz_met(n_vac))

    i_vac=1
    r_vac(i_vac)=7.d-2
    r_met(i_vac)=10d-2
    z_vac(i_vac)=100d-2

    !i_vac=2
    !r_vac(i_vac)=10d-2
    !r_met(i_vac)=r_vac(i_vac)
    !z_vac(i_vac)=70d-2

    !i_vac=3
    !r_vac(i_vac)=10d-2
    !r_met(i_vac)=r_vac(i_vac)
    !z_vac(i_vac)=75d-2
    !
    !i_vac=4
    !r_vac(i_vac)=10d-2
    !r_met(i_vac)=r_vac(i_vac)
    !z_vac(i_vac)=80d-2
    !
    r_ion_max=r_vac(1)


    !plasma region
    rs=0. ; zs=0.;
    rl=maxval(r_met); !radial position of Al sheild window
    zl=maxval(z_vac);
    dr=(rl-rs)/(real(nr)-1.) ;
    dz=(zl-zs)/(real(nz)-1.);


    !drm=dr;dzm=dz;
    ds=dr*dz;
    do i_vac=1,n_vac
        nr_vac(i_vac)=nint(real(nr)*real(r_vac(i_vac)/rl));
        nz_vac(i_vac)=nint(real(nz)*real(z_vac(i_vac)/zl));
        nr_met(i_vac)=nint(real(nr)*real(r_met(i_vac)/rl));
        nz_met(i_vac)=nint(real(nz)*real(z_vac(i_vac)/zl));
    enddo

    !r1,r2,z1,z2
    rz_dipole(1)=0.05;
    rz_dipole(2)=0.065;
    rz_dipole(3)=1.267+0.085-0.12/2.;
    rz_dipole(4)=1.267+0.085+0.12/2.;

    nrz_diploe(1)=nint(real(nr)*real(rz_dipole(1)/rl))-1;
    nrz_diploe(2)=nint(real(nr)*real(rz_dipole(2)/rl))+1;
    nrz_diploe(3)=nint(real(nz)*real(rz_dipole(3)/zl))-1;
    nrz_diploe(4)=nint(real(nz)*real(rz_dipole(4)/zl))+1;
    !
    rp=maxval(r_vac) !test!!!!!!!! not use in MHD subroutine
    nr_ion=1+nint((nr-1)*r_ion_max/rl)  ; !test!!!!!!!! not use in MHD subroutine
    nrp=maxval(nr_vac)
    !nz=maxval(nz_vac)
    !!nrsh=nr

    allocate (r(1:nr),z(1:nz),r2(1:nr),th(1:nth),r_particle_inj(1:nr_ion))
    do ir=1,nr
        r(ir)=rs+dr*(real(ir)-1.)
        r2(ir)=rs+(real(ir)-0.5)*dr;
    end do
    dr_ion=r_ion_max/(nr_ion-1)
    do ir=1,nr_ion
        r_particle_inj(ir)=rs+dr_ion*real(ir-1)
    enddo
    do ith=1,nth
        th(ith)=2*pi*real(ith-1)/real(nth-1.0D0)
    EndDo
    do iz=1,nz
        z(iz)=zs+dz*(real(iz)-1.)
    end do

    d2r=2.*dr ; d2z=2.*dz;
    dr22=dr**2; dz22=dz**2;
    dr_2=dr/2.; dz_2=dz/2.
    dr2=1/(dr**2);
    dz2=1/(dz**2);
    drz=2*(dr2+dz2)

    allocate(i_plasma_region(1:nr,1:nz),i_gamma_center_region(0:nr,0:nz,1:3))
    i_plasma_region=0;
    i_gamma_center_region=0;

    !ir1=1;
    iz1=1;
    do i_vac=1,n_vac
        !if(i_vac>1.5)ir1=nr_vac(i_vac-1)
        if(i_vac>1.5)iz1=nz_vac(i_vac-1)
        iz2=nz_vac(i_vac)
        ir2=nr_vac(i_vac)

        i_plasma_region(1:ir2,iz1:iz2)=1;
    enddo
    !1 :center region
    !11: left boundary   |<-plasma
    !12: right boundary     plasma->|
    !13: up  boundary     -^plasma
    !14: bottom boundary  _vplasma
    !15: r=0 boundary

    i_plasma_region(1,:)=i_plasma_region(1,:)-1+15 !ib_type=5; r=0 boundary
    i1=1;i2=nr_vac(1);i3=1;i4=i3;!ib_type=1;
    i_plasma_region(i1:i2,i3:i4)=i_plasma_region(i1:i2,i3:i4)-1+11  !left boundary   |<-plasma
    i1=1;i2=nr_vac(n_vac);i3=nz;i4=nz;!ib_type=2;
    i_plasma_region(i1:i2,i3:i4)=i_plasma_region(i1:i2,i3:i4)-1+12  !right boundary     plasma->|

    iz1=1;
    do i_vac=1,n_vac
        ir=nr_vac(i_vac)
        if(i_vac>1.5)iz1=nz_vac(i_vac-1)
        iz2=nz_vac(i_vac)
        i1=nr_vac(i_vac);i2=i1;i3=iz1;i4=iz2;!ib_type=3; !13: up  boundary     -^plasma
        i_plasma_region(i1:i2,i3:i4)=i_plasma_region(i1:i2,i3:i4)-1+13
    enddo
    do i_vac=1,n_vac-1
        if(nr_vac(i_vac)==nr_vac(i_vac+1))then
            ir=nr_vac(i_vac);
            iz=nz_vac(i_vac);
            i_plasma_region(ir,iz)=13
        endif
    enddo

    do i_vac=1,n_vac-1
        ir1=nr_vac(i_vac)
        ir2=nr_vac(i_vac+1)
        if(ir2>ir1)then
            i1=ir1;i2=ir2;i3=nz_vac(i_vac);i4=i3;!ib_type=1;
            i_plasma_region(i1:i2,i3:i4)=i_plasma_region(i1:i2,i3:i4)-1+11
        elseif(ir2<ir1)then
            i1=ir2;i2=ir1;i3=nz_vac(i_vac);i4=i3;!ib_type=2;
            i_plasma_region(i1:i2,i3:i4)=i_plasma_region(i1:i2,i3:i4)-1+12
        endif
    enddo


    if(iswtich_inner_dipole==1)then
        i1=nrz_diploe(1);i2=nrz_diploe(2);i3=nrz_diploe(3);i4=nrz_diploe(4);
        i_plasma_region(i1+1:i2-1,i3+1:i4-1)=0
        !i_plasma_region(i1:i2,i3:i4)=0  !test!!!!!!!!!!!!!!!!!!!!!!
        i1=nrz_diploe(1);i2=nrz_diploe(2);i3=nrz_diploe(4);!ib_type=11;
        i_plasma_region(i1:i2,i3)=i_plasma_region(i1:i2,i3)+11-1
        i1=nrz_diploe(1);i2=nrz_diploe(2);i3=nrz_diploe(3);!ib_type=12;
        i_plasma_region(i1:i2,i3)=i_plasma_region(i1:i2,i3)+12-1
        i1=nrz_diploe(1);i3=nrz_diploe(3);i4=nrz_diploe(4);!ib_type=13;
        i_plasma_region(i1,i3:i4)=i_plasma_region(i1,i3:i4)+13-1
        i1=nrz_diploe(2);i3=nrz_diploe(3);i4=nrz_diploe(4);!ib_type=14;
        i_plasma_region(i1,i3:i4)=i_plasma_region(i1,i3:i4)+14-1
        i1=nrz_diploe(1);i2=nrz_diploe(2);i3=nrz_diploe(3);i4=nrz_diploe(4);
        i_plasma_region(i1,i3)=24
        i_plasma_region(i1,i4)=23
        i_plasma_region(i2,i3)=28
        i_plasma_region(i2,i4)=27

    endif


    if(iswtich_inner_dipole==1)then
        i1=nrz_diploe(1);i2=nrz_diploe(2);i3=nrz_diploe(3);i4=nrz_diploe(4);
        i_plasma_region(i1+1:i2-1,i3+1:i4-1)=0

        i1=nrz_diploe(1);i2=nrz_diploe(2);i3=nrz_diploe(4);!ib_type=11;
        i_gamma_center_region(i1:i2,i3,3)=1
        i1=nrz_diploe(1);i3=nrz_diploe(3);i4=nrz_diploe(4);!ib_type=13;
        i_gamma_center_region(i1,i3:i4,1)=1
        i1=nrz_diploe(2);i3=nrz_diploe(3);i4=nrz_diploe(4);!ib_type=13;
        i_gamma_center_region(i1,i3:i4,1)=1
    endif

    do ir=1,nr
        do iz=1,nz
            if(i_plasma_region(ir,iz)==1)i_gamma_center_region(ir,iz,1:3)=1
        end do
    end do
    i_gamma_center_region(1,2:nz-1,1)=1

    iz1=1;
    ir1=1;
    do i_vac=1,n_vac-1
        if(i_vac>1.5)ir1=nr_vac(i_vac-1)
        if(i_vac>1.5)iz1=nz_vac(i_vac-1)
        ir2=nr_vac(i_vac)
        iz2=nz_vac(i_vac)

        if(ir1<ir2)i_gamma_center_region(2:ir2-1,iz1,3)=1
    enddo

    nzx2=nz*2
    nz6=nz*6
    allocate (b0_DC(1:nr,1:nz,1:4),FindRegion(1:nr,1:nz))
    !allocate (emf_mode(1:nr,1:nz,1:3),b_mode(1:nr,1:nz,1:3),jp_mode(1:nr,1:nz,1:3))
    allocate (fvm(1:nr,1:nz))
    allocate (te(1:nr,1:nz))
    allocate (E_rf(1:nr,1:nz,1:ns),B_rf(1:nr,1:nz,1:ns),jp(1:nr,1:nz,1:ns)  )
    allocate (ne(1:nr,1:nz))  !b0_DC(:,:,:1)-Br, b0_DC(:,:,:2)-B_theta, b0_DC(:,:,:3)-Bz, , b0_DC(:,:,:4)-abs(b0_DC)
    allocate (ja_m_nor(1:nr,1:nz,1:n3),ja_record(1:nr,1:nth,1:nz,1:n3),ja_AC(1:nr,1:nth,1:nz,1:n3))
    allocate (e_m(1:nr,1:nz,1:n3),e_record(1:nr,1:nth,1:nz,1:n3),e_AC(1:nr,1:nth,1:nz,1:n3),e_int(1:nr,1:nz,1:n3))
    allocate (jp_m(1:nr,1:nz,1:n3),jp_record(1:nr,1:nth,1:nz,1:n3),jp_AC(1:nr,1:nth,1:nz,1:n3))
    allocate (b_m(1:nr,1:nz,1:n3),b_record(1:nr,1:nth,1:nz,1:n3),b_AC(1:nr,1:nth,1:nz,1:n3))
    allocate (ep(1:nr,1:nz,1:9), si(1:nr,1:nz,1:n3,1:n3))
    allocate(A_helicon.col(1:30*nr*nz),A_helicon.row(1:30*nr*nz),A_helicon.val(1:30*nr*nz))
    allocate (eq_b(1:3*nr*nz),eq_xlast(m_start:m_end,1:3*nr*nz))
    allocate(peff(1:nr,1:nz),ptotm(1:nr,1:nz,1:n3),power_depo(1:nr,1:nz,1:n3))
    allocate(e_output(m_start:m_end,1:nr,1:nz,1:n3),jp_output(m_start:m_end,1:nr,1:nz,1:n3),ptot_output(m_start:m_end,1:nr,1:nz))
    allocate(Xkz_e(m_start:m_end,-(nz-1)/2:(nz-1)/2,1:n3))
    allocate(ptotal_m(m_start:m_end),max_eth(m_start:m_end))
    allocate(pdf_ne_r(1:nr_ion),pdf_ne_source_r(1:nr_ion),pdf_ne_z(1:nz),pdf_ne_source_z(1:nz)) !,Ez_ambi_eq(1:nz)
    allocate(density_2D(1:nr,1:nz),Es_2D(1:nr,1:nz,1:2),Ek_ion_2D(1:nr,1:nz),te_2D(1:nr,1:nz)) !,density_i(1:nz),z_half(0:nz),Es_z_half(0:nz)
    allocate (Erf_PIC(1:nr,1:nz,1:n3))
    s_drz=dr*dz

    ne=0.;e_int=0.
    e_AC=0.0D0; b_AC=0.0D0; jp_AC=0.0D0;
    ptotm=0.0D0;power_depo=0.
    fvm=0.0D0
    Xkz_e=0.0D0+0.0D0*i
    end subroutine


