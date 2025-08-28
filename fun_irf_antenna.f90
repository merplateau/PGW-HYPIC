

    !--------------------subroutine antenna current set-----------------!
    subroutine antenna
    use the_whole_varibles
    implicit none
    real*8:: width_r,width_th,width_z,la,ra_p,z1_m,z2_p,elth,elz
    real*8:: z1,z2,z3,z4,ra1
    real*8:: sq_rz,sq_rth,j_loop,j_helical
    integer::nz1m,nz2p1,nz2p,nz1m1,nrz_dz,nrz_dr,nz1,nz2,nz3,nz4,nra1
    complex*16 ::c1,c2
    real*8 :: j_loop_nor,j_helical_nor !ja_value_nor,


    !antenna set
    ra=8d-2
    za=zl/2 !1.267-[0, 0.5]
    width_r=1d-2 !10 mm (+10 mm)
    width_z=4d-2 !30 mm (+10 mm)
    width_th=width_z;
    la=40d-2     !---only for helical and Nagoya antenna---!

    ra1=ra+width_r
    z1=za-la/2;
    z2=z1+la;
    z3=za-width_z/2.;
    z4=za+width_z/2.
    z1_m=z1-width_z
    z2_p=z2+width_z
    !---only for helical antenna---!

    rz_ant_region(1)=ra-width_r
    rz_ant_region(2)=ra+width_r*2
    rz_ant_region(3)=za-width_z;
    rz_ant_region(4)=za+width_z;


    !===================== by Xuad =====================!
    nra=1+nint((nr-1)*ra/rl);
    nra1=1+nint((nr-1)*ra1/rl)
    !nra1=1+nint((nr-1)*ra_p/rl);

    nza=1+nint((nz-1)*za/zl)  ;
    nz1=1+nint((nz-1)*z1/zl);
    nz2=1+nint((nz-1)*z2/zl);
    nz3=1+nint((nz-1)*z3/zl)  ;
    nz4=1+nint((nz-1)*z4/zl)
    nz1m=1+nint((nz-1)*z1_m/zl);
    nz1m1=nz1-1;
    nz2p=1+nint((nz-1)*z2_p/zl);
    nz2p1=nz2+1;
    !nzp1=1+nint((nz-1)*zp1/zl);
    !nzp2=1+nint((nz-1)*zp2/zl)

    !===================== by Xuad =====================!


    elth=pi*ra/sqrt((pi*ra)**2+la**2);
    elz=la/sqrt((pi*ra)**2+la**2);

    nrz_dr=nint(nr*width_r/rl)
    nrz_dz=nint(nz*width_z/zl)

    !if(iswitch_antenna_type==0)sq_rz=width_r*width_z
    !if(iswitch_antenna_type==1 .or. iswitch_antenna_type==3)sq_rz=dr*(nra1-nra)*dz*(nz1-nz1m)
    sq_rz=width_r*width_z
    sq_rth=width_r*width_th

    !j_loop=1./sq_rz
    !j_helical=1./sq_rth
    j_loop_nor=1./sq_rz
    j_helical_nor=1./sq_rth


    im=m*i
    imp=im*pi
    ja_m_nor=0.D0;
    !--------------------------Antenna current set--------------------------!
    !!-----single loop----------!
    if(iswitch_antenna_type==0)then
        if (m==0) then
            ja_m_nor(nra:nra1,nz3:nz4,2)=j_loop_nor;

            !         ja_m_nor(nra:nra1,nz3:nz4,3)=j_loop_nor;   !test!!!!!!!!!!!!!!
        else
            ja_m_nor=0;
        endif
        !!-----single loop----------!

        !-----half loop antenna----------!
    elseif(iswitch_antenna_type==2)then
        if (m==0) then
            !----------z=center---------!
            ja_m_nor(nra:nra1,nz3:nz4,2)=j_loop_nor/4
        else
            !----------z=center---------!
            ja_m_nor(nra:nra1,nz3:nz4,2)=j_loop_nor*(1-exp(-imp))/(4*imp) !!!
        endif
        !-----half loop antenna----------!
    else
        if (m==0) then
            !----------z=z1 &z=z2---------!
            ja_m_nor(nra:nra1,nz1m:nz1m1,2)=0.0D0
            ja_m_nor(nra:nra1,nz2p1:nz2p,2)=0.0D0
            !----------z1<z<z2---------!
            ja_m_nor(nra:nra1,nz1:nz2,2)=0;
            ja_m_nor(nra:nra1,nz1:nz2,3)=0;
        else
            !----------z=z1---------!
            ja_m_nor(nra:nra1,nz1m:nz1m1,2)=j_loop_nor*(1-exp(-imp))/(2*imp)
            !----------z=z2---------!
            if(iswitch_antenna_type /= 3) ja_m_nor(nra:nra1,nz2p1:nz2p,2)=j_loop_nor*(1-exp(-imp))/(2*imp)
            if(iswitch_antenna_type == 3) ja_m_nor(nra:nra1,nz2p1:nz2p,2)=-j_loop_nor*(1-exp(-imp))/(2*imp)

            !----------z1<z<z2---------!
            c1=im*width_th/(2*ra)
            !-----Right helical antenna----------!
            if(iswitch_antenna_type==1)then
                do ir=nra,nra1
                    do iz=nz1,nz2
                        c2=j_helical_nor*(1-exp(-imp))*exp(-imp*(z(iz)-z1)/la)*( exp(-c1)-exp(c1)  )/(2*imp)
                        ja_m_nor(ir,iz,2)=elth*c2
                        ja_m_nor(ir,iz,3)=elz*c2
                    enddo
                enddo
                !-----Left helical antenna----------!
            elseif(iswitch_antenna_type==-1)then
                do ir=nra,nra1
                    do iz=nz1,nz2
                        c2=j_helical_nor*(1-exp(-imp))*exp(imp*(z(iz)-z1)/la)*( exp(-c1)-exp(c1)  )/(2*imp)
                        ja_m_nor(ir,iz,2)=-elth*c2
                        ja_m_nor(ir,iz,3)=elz*c2
                    enddo
                enddo
                !-----Nagoya Type III antenna----------!
            elseif(iswitch_antenna_type==3)then
                do ir=nra,nra1
                    do iz=nz1,nz2
                        ja_m_nor(ir,iz,3)=j_helical_nor*(1-exp(-imp))*( exp(-c1)-exp(c1)  )/(2*imp)
                    enddo
                enddo
            endif
        endif
    endif

    ja_m_nor=ja_m_nor*1.0D0
    continue
    !pause
    end subroutine



    subroutine find_irf_now
    use the_whole_varibles
    implicit none

    !i_now=sqrt(power_Joule/resistance);

    endsubroutine find_irf_now

    subroutine update_and_calculate
    use the_whole_varibles
    implicit none
    integer i1,i2,i4
    real*8 Erf_tmp(1:nr,1:nz6)    


    e_AC=i_now*e_AC
    b_AC=i_now*b_AC
    !dive_AC=i_now*dive_AC
    ja_AC=i_now*ja_AC
    jp_AC=i_now*jp_AC
    e_int=i_now*e_int
    
!--use this won't repeat running absolutly because of small uncertainties caused by pardiso function in FDFD----!    
    !Erf_PIC=e_int 
    
!--reduce the significance digit, to avoid small uncertainties caused by pardiso function in FDFD --------------!    
301 format(<nr>(e12.5,' '))  
    open (unit=201,file='Erf_tmp.txt',status='unknown',iostat=ierror)
    do itp=1,3
        write (201,301)real(e_int(1:nr,1:nz,itp))
        write (201,301)imag(e_int(1:nr,1:nz,itp))
    enddo
    close(201)
    open (unit=2001,file='Erf_tmp.txt',status='old',action='read',iostat=ierror)
    read (2001,301)Erf_tmp
    close(2001)
    do itp=1,3
        i1=(itp-1)*(2*nz)+1
        i2=i1+nz-1
        i3=i2+1
        i4=i3+nz-1
        Erf_PIC(1:nr,1:nz,itp)=Erf_tmp(1:nr,i1:i2)+i*Erf_tmp(1:nr,i3:i4)
    enddo

    
    peff=i_now*i_now*peff
    power_depo=i_now*i_now*power_depo
    ptotm=i_now*i_now*ptotm
    ptotal=i_now*i_now*ptotal
    ptotal_boundary=i_now*i_now*ptotal_boundary
    ptotal_inside=i_now*i_now*ptotal_inside
    ptotal_m(:)=i_now*i_now*ptotal_m(:)

    e_output=i_now*e_output
    Xkz_e=i_now*Xkz_e
    max_eth(:)=i_now*max_eth(:)
    jp_output=i_now*jp_output
    ptot_output=i_now*i_now*ptot_output

    do ir=1,nr
        do iz=1,nz
            do i3=1,n3
                E_rf(ir,iz,i3)=maxval(abs(e_AC(ir,:,iz,i3)))
                B_rf(ir,iz,i3)=maxval(abs(b_AC(ir,:,iz,i3)))
                jp(ir,iz,i3)=maxval(abs(jp_AC(ir,:,iz,i3)))
            enddo
        enddo
    enddo
    

    
    end subroutine update_and_calculate