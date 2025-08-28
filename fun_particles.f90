    SUBROUTINE interp_position(n_len,pdf_x,x_grid,x_out)
    implicit none
    integer::n_bef,n_len
    real*8::rand_0_1,x_out,pdf_x(1:n_len),x_grid(1:n_len)
    integer::I_x,ix_bef1,ix_bef2
    real*8::s1,s2,min_x,dx_tp1,dx_tp2

    call random_number(rand_0_1)
    I_x=minloc(abs(rand_0_1-pdf_x),1);
    min_x=rand_0_1-pdf_x(I_x);
    if (min_x<0)then
        ix_bef1=I_x-1;
        ix_bef2=I_x;
    else
        ix_bef1=I_x;
        ix_bef2=I_x+1;
    endif

    if(ix_bef1<1)then
        ix_bef1=1
        ix_bef2=2
    elseif(ix_bef2>n_len)then
        ix_bef2=n_len
        ix_bef1=ix_bef2-1
    endif

    dx_tp1=rand_0_1-pdf_x(ix_bef1);
    dx_tp2=pdf_x(ix_bef2)-rand_0_1;
    s1=dx_tp1/(dx_tp1+dx_tp2);
    s2=1.-s1;
    x_out=s2*x_grid(ix_bef1)+s1*x_grid(ix_bef2);
    continue
    endSUBROUTINE interp_position


    SUBROUTINE particles_initialization
    USE the_whole_varibles
    implicit none
    integer*4 :: nseed(1)
    real*8::rand_0_1,rand_1_1,xtp,ytp,rtp,th_tp
    real*8::vr_tp,vx_tp,vy_tp,vz_tp,br_tp,bz_tp,b0_tp,sinth_b0,costh_b0
    real*8::ztp,vz_ll,coeff_vx

    !set a specific random number
    nseed=123456
    call RANDOM_SEED(put=nseed)
    
    !change the seed for the pseudorandom number generator.
    !Uncomment and simulation won't repeat absolutly run to run.
    !call RANDOM_SEED()    

    coeff_vx=1;
    ip=1
    np=0;
    do while (ip<=np_ini)
        !use PDF and random number to generate specific particle distribution
        call interp_position(nr_ion,pdf_ne_r,r_particle_inj,rtp)
        if(rtp<0)rtp=abs(rtp)
        if(rtp<r_ion_max)then
            call random_number(rand_0_1)
            th_tp=2*pi*rand_0_1
            xtp=rtp*cos(th_tp);
            ytp=rtp*sin(th_tp);

            call interp_position(nz,pdf_ne_z,z,ztp)
            if(ztp<zs)ztp=2*zs-ztp
            if(ztp>zl)ztp=2*zl-ztp

            x(ip,1)=xtp
            x(ip,2)=ytp
            x(ip,3)=ztp

            call random_number(tp1)
            call random_number(tp2)
            tp3=sqrt(-2*log(tp1))*cos(2*pi*tp2)
            v(ip,1)=coeff_vx*vi_ex*tp3
            v_e(ip,1)=coeff_vx*ve_ex*tp3

            call random_number(tp1)
            call random_number(tp2)
            tp3=sqrt(-2*log(tp1))*cos(2*pi*tp2)
            v(ip,2)=vi_ex*tp3
            v_e(ip,2)=ve_ex*tp3

            call random_number(tp1)
            call random_number(tp2)
            tp3=sqrt(-2*log(tp1))*cos(2*pi*tp2)
            v(ip,3)=vi_ex*tp3
            v_e(ip,3)=ve_ex*tp3

            np=np+1
            ip=ip+1;
        endif
    end do

    continue
    END SUBROUTINE particles_initialization



    SUBROUTINE particles_inject
    USE the_whole_varibles
    implicit none
    integer::ip_inj
    real*8::rand_0_1,rand_1_1,xtp,ytp,rtp,th_tp
    real*8::vr_tp,vx_tp,vy_tp,vz_tp,br_tp,bz_tp,b0_tp,sinth_b0,costh_b0
    real*8::ztp,vz_ll

    call find_func_cputime_1_of_2  !-------------------1/2
    
    !change the seed for the pseudorandom number generator.
    !Uncomment and simulation won't repeat absolutly run to run.
    !call RANDOM_SEED()
    
    ip_inj=0
    num_inject=np_max-np
    do while (ip_inj<num_inject)
        call interp_position(nr_ion,pdf_ne_source_r,r_particle_inj,rtp)
        if(rtp<0)rtp=abs(rtp)
        if(rtp<r_ion_max)then
            ip_inj=ip_inj+1
            np=np+1
            ip=np;
            if( np>np_max) exit

            call random_number(rand_0_1)
            th_tp=2*pi*rand_0_1
            xtp=rtp*cos(th_tp);
            ytp=rtp*sin(th_tp);

            call interp_position(nz,pdf_ne_source_z,z,ztp)
            if(ztp<zs)ztp=2*zs-ztp
            if(ztp>zl)ztp=2*zl-ztp

            x(ip,1)=xtp
            x(ip,2)=ytp
            x(ip,3)=ztp

            call random_number(tp1)
            call random_number(tp2)
            tp3=sqrt(-2*log(tp1))*cos(2*pi*tp2)
            v(ip,1)=vi_ex*tp3

            call random_number(tp1)
            call random_number(tp2)
            tp3=sqrt(-2*log(tp1))*cos(2*pi*tp2)
            v(ip,2)=vi_ex*tp3

            call random_number(tp1)
            call random_number(tp2)
            tp3=sqrt(-2*log(tp1))*cos(2*pi*tp2)
            v(ip,3)=vi_ex*tp3
        endif
    end do
    N_inject_particles=num_inject
    call find_func_cputime_2_of_2(func_time(1))  !-----2/2

    END SUBROUTINE particles_inject


    SUBROUTINE mover
    USE the_whole_varibles
    implicit none
    real*8:: B_mover(1:3),E_mover(1:3) !,xv(1:3)

    call find_func_cputime_1_of_2  !-------------------1/2
    call find_x_to_grid
    call find_density_and_Es_2D
    call find_func_cputime_2_of_2(func_time(3))  !-----2/2

    call find_func_cputime_1_of_2  !-------------------1/2
    expt=exp(-i*2*pi*frequency*t)
    do ip=1,np
        call interpolation_E_B(B_mover,E_mover)
        !E_mover=0. !energy conservation
        call push_RK4(B_mover,E_mover)
    enddo
    call find_func_cputime_2_of_2(func_time(4))  !-----2/2

    END SUBROUTINE mover



    SUBROUTINE find_x_to_grid !here x means (x,y,z)
    USE the_whole_varibles
    implicit none
    integer::ir_tp,iz_tp
    real*8::r_p(1:np)
    !r: s1= x_to_grid(1:np,1),s2=1-s1; z: s3= x_to_grid(1:np,2),s4=1-s3;
    x_to_grid=0.
    r_p(1:np)=sqrt(x(1:np,1)**2+x(1:np,2)**2)
    ir1_iz1_grid(1:np,1)=int((r_p(1:np)-r(1))/dr)+1
    ir1_iz1_grid(1:np,2)=int((x(1:np,3)-z(1))/dz)+1
    do ip=1,np
        ir_tp=ir1_iz1_grid(ip,1)
        iz_tp=ir1_iz1_grid(ip,2)
        x_to_grid(ip,1)=(r_p(ip)-r(ir_tp))/dr
        x_to_grid(ip,2)=(x(ip,3)-z(iz_tp))/dz
    enddo
    END SUBROUTINE find_x_to_grid



    SUBROUTINE interpolation_E_B(B_mover,E_mover)
    USE the_whole_varibles
    implicit none
    INTEGER*2 ir1,ir2,iz1,iz2
    real*8    sinth,costh,rtp,dr_tp1,dr_tp2,dz_tp1,dz_tp2
    Complex*16 c123(1:3)
    real*8    B_mover(1:3),E_mover(1:3)
    Real*8    min_r,min_z,ztp
    Real*8    ex,ey,ez,bx,by,bz,br
    Real*8    Ez_max,z_peak_Ez,width_Ez_z
    real*8    s1,s2,s3,s4,Ez_dc,Er_dc 
    real*8    erthz(1:3)

    rtp=sqrt(x(ip,1)**2+x(ip,2)**2)+1e-20;
    ztp=x(ip,3);
    !s_drz=dr*dz  !caculated out of this function
    !expt=exp(-i*2*pi*frequency*t)   !caculated out of this function
    sinth=x(ip,2)/rtp
    costh=x(ip,1)/rtp

    s1=x_to_grid(ip,1)*x_to_grid(ip,2)
    s2=x_to_grid(ip,1)*(1.-x_to_grid(ip,2))
    s3=(1.-x_to_grid(ip,1))*(1.-x_to_grid(ip,2))
    s4=(1.-x_to_grid(ip,1))*x_to_grid(ip,2)

    ir1=ir1_iz1_grid(ip,1)
    iz1=ir1_iz1_grid(ip,2)
    ir2=ir1+1
    iz2=iz1+1

    br=s3*b0_DC(ir1,iz1,1)+s1*b0_DC(ir2,iz2,1)+s2*b0_DC(ir2,iz1,1)+s4*b0_DC(ir1,iz2,1);
    bz=s3*b0_DC(ir1,iz1,3)+s1*b0_DC(ir2,iz2,3)+s2*b0_DC(ir2,iz1,3)+s4*b0_DC(ir1,iz2,3);
    bx=br*costh
    by=br*sinth;


    Er_dc=s3*Es_2D(ir1,iz1,1)+s1*Es_2D(ir2,iz2,1)+s2*Es_2D(ir2,iz1,1)+s4*Es_2D(ir1,iz2,1);
    Ez_dc=s3*Es_2D(ir1,iz1,2)+s1*Es_2D(ir2,iz2,2)+s2*Es_2D(ir2,iz1,2)+s4*Es_2D(ir1,iz2,2);


    if(t>t_power_on)then
        c123=s3*Erf_PIC(ir1,iz1,:)+s1*Erf_PIC(ir2,iz2,:)+s2*Erf_PIC(ir2,iz1,:)+s4*Erf_PIC(ir1,iz2,:)
        tp1=i_now*real(expt*c123(1))+Er_dc;
        tp2=i_now*real(expt*c123(2));
        tp3=i_now*real(expt*c123(3))+Ez_dc;
        erthz(1)=tp1
        erthz(2)=tp2
        erthz(3)=tp3
    else
        erthz(1)=Er_dc
        erthz(2)=0
        erthz(3)=Ez_dc
    endif

    ! (r,th,z) ->(x,y,z)
    ex=erthz(1)*costh-erthz(2)*sinth
    ey=erthz(1)*sinth+erthz(2)*costh
    ez=erthz(3)+Ez_dc

    if(ip==1)then
        e_rec(1)=ex
        e_rec(2)=ey
        e_rec(3)=ez
        b_rec(1)= v(ip,1)*costh+v(ip,2)*sinth
        b_rec(2)=-v(ip,1)*sinth+v(ip,2)*costh
        b_rec(3)=erthz(1)
        b_rec(4)=erthz(2)
    endif

    B_mover(1)=bx
    B_mover(2)=by
    B_mover(3)=bz
    E_mover(1)=ex
    E_mover(2)=ey
    E_mover(3)=ez
    END SUBROUTINE interpolation_E_B



    SUBROUTINE push_RK4(B_mover,E_mover)
    USE the_whole_varibles
    implicit none
    real*8::ex,ey,eth,ez,br,bx,by,bz
    real*8:: B_mover(1:3),E_mover(1:3) !,xv(1:3)
    real*8:: k1(1:6),k2(1:6),k3(1:6),k4(1:6)
    real*8::xv(1:6)

    bx=B_mover(1);by=B_mover(2);bz=B_mover(3);
    ex=E_mover(1);ey=E_mover(2);ez=E_mover(3);

    !find k1
    xv(1:3)=x(ip,1:3)
    xv(4:6)=v(ip,1:3)
    k1(1:3)=xv(4:6)
    k1(4)=q_mass*(ex+xv(5)*bz-xv(6)*by);
    k1(5)=q_mass*(ey+xv(6)*bx-xv(4)*bz)
    k1(6)=q_mass*(ez+xv(4)*by-xv(5)*bx)

    !find k2
    xv(1:3)=x(ip,1:3)+dt*0.5*k1(1:3)
    xv(4:6)=v(ip,1:3)+dt*0.5*k1(4:6)
    k2(1:3)=xv(4:6)
    k2(4)=q_mass*(ex+xv(5)*bz-xv(6)*by);
    k2(5)=q_mass*(ey+xv(6)*bx-xv(4)*bz)
    k2(6)=q_mass*(ez+xv(4)*by-xv(5)*bx)

    !find k3
    xv(1:3)=x(ip,1:3)+dt*0.5*k2(1:3)
    xv(4:6)=v(ip,1:3)+dt*0.5*k2(4:6)
    k3(1:3)=xv(4:6)
    k3(4)=q_mass*(ex+xv(5)*bz-xv(6)*by);
    k3(5)=q_mass*(ey+xv(6)*bx-xv(4)*bz)
    k3(6)=q_mass*(ez+xv(4)*by-xv(5)*bx)

    !find k4
    xv(1:3)=x(ip,1:3)+dt*k3(1:3)
    xv(4:6)=v(ip,1:3)+dt*k3(4:6)
    k4(1:3)=xv(4:6)
    k4(4)=q_mass*(ex+xv(5)*bz-xv(6)*by);
    k4(5)=q_mass*(ey+xv(6)*bx-xv(4)*bz)
    k4(6)=q_mass*(ez+xv(4)*by-xv(5)*bx)

    x(ip,1:3)=x(ip,1:3)+dt/6.*(k1(1:3)+2.*k2(1:3)+2.*k3(1:3)+k4(1:3))
    v(ip,1:3)=v(ip,1:3)+dt/6.*(k1(4:6)+2.*k2(4:6)+2.*k3(4:6)+k4(4:6))
    END SUBROUTINE push_RK4



    SUBROUTINE find_density_and_Es_2D
    USE the_whole_varibles
    implicit none
    real*8::dni_ni,ti_tmp,te_min_set,density_cri,N_te,Ek_ion_center,ratio_te !te_max_set,
    !real*8::ni_int
    integer n_sm,ifilter,ir1,ir2,iz1,iz2 !,itp

    te_min_set=2. !eV
    te_2D=0.
    call density_Ek_2D_sub
    itp=2
    do itp=1,2
        if(itp==1)n_sm=3   !1st smooth
        if(itp==2)n_sm=2   !2st smooth

        ir1=1;ir2=nr_vac(1);iz1=1;iz2=nz;
        call smooth_2d(n_sm,2*n_sm,ir1,ir2,iz1,iz2,density_2D(ir1:ir2,iz1:iz2))
        call smooth_2d(n_sm,2*n_sm,ir1,ir2,iz1,iz2, Ek_ion_2D(ir1:ir2,iz1:iz2))
    enddo


    te_2D=te_ave    !set uniform electron temperature Te by mcc
    !te_2D=0.05*ti_ave


    !Es=-Te*grad(ne)/ne with ne=ni
    Es_2D=0.
    Es_2D(2:nr_ion-1,1:nz,1)=-te_2D(2:nr_ion-1,1:nz)*(density_2D(3:nr_ion,1:nz)- &
        & density_2D(1:nr_ion-2,1:nz))/(2*dr*density_2D(2:nr_ion-1,1:nz))
    Es_2D(1:nr_ion,2:nz-1,2)=-te_2D(1:nr_ion,2:nz-1)*(density_2D(1:nr_ion,3:nz)- &
        & density_2D(1:nr_ion,1:nz-2))/(2*dz*density_2D(1:nr_ion,2:nz-1))


    !boundary condition
    Es_2D(1,1:nz,1)=-te_2D(1,1:nz)*(density_2D(2,1:nz)-density_2D(1,1:nz))/(dr*density_2D(1,1:nz))
    Es_2D(nr_ion,1:nz,1)=-te_2D(nr_ion,1:nz)*(density_2D(nr_ion,1:nz)- &
        & density_2D(nr_ion-1,1:nz))/(dr*density_2D(nr_ion,1:nz))
    Es_2D(1:nr_ion,1,2)=-te_2D(1:nr_ion,1)*(density_2D(1:nr_ion,2)- &
        & density_2D(1:nr_ion,1))/(dz*density_2D(1:nr_ion,1))
    Es_2D(1:nr_ion,nz,2)=-te_2D(1:nr_ion,nz)*(density_2D(1:nr_ion,nz)- &
        & density_2D(1:nr_ion,nz-1))/(dz*density_2D(1:nr_ion,nz))
    END SUBROUTINE find_density_and_Es_2D



    SUBROUTINE density_Ek_2D_sub
    USE the_whole_varibles
    implicit none
    real*8::s1,s2,s3,s4,vtp2
    integer::ir2,ir1,iz2,iz1

    density_2D=0.2
    Ek_ion_2D=0.

    !s1-s4, squre of bottom left, bottom right, upper right, upper left
    ! ^
    ! |  Radial coordinates r, to up
    ! -> Axial coordinates  z, to right
    !  |s4| |s3|
    !   x(r,z)
    !  |s1| |s2|

    do ip=1,np
        s1=x_to_grid(ip,1)*x_to_grid(ip,2)
        s2=x_to_grid(ip,1)*(1.-x_to_grid(ip,2))
        s3=(1.-x_to_grid(ip,1))*(1.-x_to_grid(ip,2))
        s4=(1.-x_to_grid(ip,1))*x_to_grid(ip,2)

        ir1=ir1_iz1_grid(ip,1)
        iz1=ir1_iz1_grid(ip,2)
        ir2=ir1+1
        iz2=iz1+1

        density_2D(ir1,iz1)=density_2D(ir1,iz1)+s3
        density_2D(ir2,iz2)=density_2D(ir2,iz2)+s1
        density_2D(ir2,iz1)=density_2D(ir2,iz1)+s2
        density_2D(ir1,iz2)=density_2D(ir1,iz2)+s4

        vtp2=v(ip,1)**2+v(ip,2)**2+v(ip,3)**2
        Ek_ion_2D(ir1,iz1)=Ek_ion_2D(ir1,iz1)+s3*vtp2
        Ek_ion_2D(ir2,iz2)=Ek_ion_2D(ir2,iz2)+s1*vtp2
        Ek_ion_2D(ir2,iz1)=Ek_ion_2D(ir2,iz1)+s2*vtp2
        Ek_ion_2D(ir1,iz2)=Ek_ion_2D(ir1,iz2)+s4*vtp2
    enddo
    !Ek_ion_2D must be divided by the density without cylindrical coefficients (2pi*rdr)
    Ek_ion_2D=mass_q_i_05*Ek_ion_2D/density_2D;

    !Note that densities in cylindrical coordinate system need to be divided by the (2pi*rdr)
    do ir=2,nr_ion!
        density_2D(ir,:)=density_2D(ir,:)/(pi*(r2(ir)**2-r2(ir-1)**2)*dz)
    enddo
    ir=1;
    density_2D(ir,:)=density_2D(ir+1,:)

    END SUBROUTINE density_Ek_2D_sub




    SUBROUTINE smooth_1d(n_sm,n_data,data_in)
    implicit none
    integer::n_sm,n_data,n_total,k1,k2,k3,ix_sm
    real*8::data_in(1:n_data),data_out(1:n_data)

    n_total=2*n_sm+1
    do ix_sm=n_sm+1,n_data-n_sm
        k1=ix_sm-n_sm
        k2=ix_sm+n_sm
        data_out(ix_sm)=sum(data_in(k1:k2))/n_total
    enddo
    do ix_sm=2,n_sm
        k1=1
        k2=(ix_sm-1)*2+1
        k3=k2-k1+1
        data_out(ix_sm)=sum(data_in(k1:k2))/k3
    enddo
    do ix_sm=n_data-n_sm+1,n_data-1
        k1=ix_sm-(n_data-ix_sm)
        k2=n_data
        k3=k2-k1+1
        data_out(ix_sm)=sum(data_in(k1:k2))/k3
    enddo
    data_in(2:n_data-1)=data_out(2:n_data-1)
    endSUBROUTINE smooth_1d



    SUBROUTINE smooth_2d(n_smr,n_smz,nr_sta,nr_end,nz_sta,nz_end,data_in)
    implicit none
    integer  n_smr,n_smz,nr_sta,nr_end,nz_sta,nz_end,n_total,k1,k2,k3,k4,ir_sm,iz_sm
    integer  ir_reg1,ir_reg2,iz_reg1,iz_reg2
    real*8   data_in(nr_sta:nr_end,nz_sta:nz_end),data_out(nr_sta:nr_end,nz_sta:nz_end)

    ir_reg1=n_smr+nr_sta-1
    ir_reg2=nr_end-n_smr
    iz_reg1=n_smz+nz_sta-1
    iz_reg2=nz_end-n_smz
    do ir_sm=nr_sta,nr_end
        do iz_sm=nz_sta,nz_end
            if( ir_sm<=ir_reg1)then
                k1=1+nr_sta-1;k2=ir_sm+(ir_sm-k1)
            elseif(ir_sm>=ir_reg2)then
                k2=nr_end;k1=ir_sm-(nr_end-ir_sm);
            else
                k1=ir_sm-n_smr
                k2=ir_sm+n_smr
            endif

            if( iz_sm<=iz_reg1)then
                k3=1+nz_sta-1;k4=iz_sm+(iz_sm-k3)
            elseif(iz_sm>=iz_reg2)then
                k4=nz_end;k3=iz_sm-(nz_end-iz_sm);
            else
                k3=iz_sm-n_smz
                k4=iz_sm+n_smz
            endif

            n_total=(k4-k3+1)*(k2-k1+1)
            data_out(ir_sm,iz_sm)=sum(data_in(k1:k2,k3:k4))/n_total
        enddo
    enddo

    data_in=data_out
    endSUBROUTINE smooth_2d



    SUBROUTINE escape
    USE the_whole_varibles
    implicit none
    INTEGER*2::ip2,iloss
    real*8::zb0,rb1,zb1,rb2,zb2,rb3,zb3,rb4,zb4,rtp,ztp
    call find_func_cputime_1_of_2  !-------------------1/2

    zb0=0.;
    rb1=r_ion_max;    zb1=zl;

    N_lost_particles=0
    !Ek_loss=0.
    ip=1
    do while (ip<=np)
        rtp=sqrt(x(ip,1)**2+x(ip,2)**2) !+1e-6;
        ztp=x(ip,3)
        if (  ztp<zb0 )call escape_sub2(1)
        if (  ztp>zb1 )call escape_sub2(2)
        if (rtp>rb1)call escape_sub2(3)
        ip=ip+1;
    enddo

    if(N_lost_particles/=0)N_lost_particles_showup=N_lost_particles
    call find_func_cputime_2_of_2(func_time(2))  !-----2/2
    END SUBROUTINE escape


    SUBROUTINE escape_sub
    USE the_whole_varibles
    implicit none

    !x_to_bd(ip)=1
    v(ip,1:3)=v(np,1:3);
    x(ip,1:3)=x(np,1:3);
    ip=ip-1;
    np=np-1;
    N_lost_particles=N_lost_particles+1
    end SUBROUTINE escape_sub

    SUBROUTINE escape_sub2(iloss)
    USE the_whole_varibles
    implicit none
    INTEGER*2::iloss

    Ek_loss_tol(iloss)=Ek_loss_tol(iloss)+sum(v(ip,1:3)**2)
    num_loss(iloss)=num_loss(iloss)+1
    !x_to_bd(ip)=1
    v(ip,1:3)=v(np,1:3);
    x(ip,1:3)=x(np,1:3);
    ip=ip-1;
    np=np-1;
    N_lost_particles=N_lost_particles+1
    end SUBROUTINE escape_sub2



    SUBROUTINE escape_periodicity
    USE the_whole_varibles
    implicit none
    real*8::zb0,zb1,rtp,ztp
    call find_func_cputime_1_of_2  !-------------------1/2
    zb0=0.;
    zb1=zl;
    do ip=1,np
        rtp=sqrt(x(ip,1)**2+x(ip,2)**2) !+1e-6;
        ztp=x(ip,3)
        if (  ztp<zb0 )x(ip,3)=zl+x(ip,3)
        if (  ztp>zb1 )x(ip,3)=x(ip,3)-zl
        if (rtp>r_ion_max)then
            v(ip,1:2)=-v(ip,1:2)
            x(ip,1:2)=x(ip,1:2)+dt*v(ip,1:2)
        endif
    enddo
    call find_func_cputime_2_of_2(func_time(2))  !-----2/2
    endSUBROUTINE escape_periodicity



    subroutine find_power_loss
    use the_whole_varibles
    implicit none
    if(mod(t,10.*trf)<dt .or. it==1)then
        power_loss_ave(1:4)=n_macro*0.5*mi*Ek_loss_tol(1:4)/(10.*trf);!max_density_set*cm3/maxval(density_2D)
        Ek_loss_ave(1:4)=mass_q_i_05*Ek_loss_tol(1:4)/num_loss(1:4)
        num_loss=1e-2
        Ek_loss_tol=0.
    endif
    end subroutine find_power_loss



    SUBROUTINE filter_1d(ndata,data_in)
    implicit none
    integer::ndata,idata
    real*8::w,data_in(1:ndata),data_out(1:ndata)
    w=0.5

    idata=1;      data_out(idata)=.5*(data_in(idata)+data_in(idata+1))
    idata=ndata;  data_out(idata)=.5*(data_in(idata)+data_in(idata-1))
    do idata=2,ndata-1
        data_out(idata)=.25*(data_in(idata-1)+2.*data_in(idata)+data_in(idata+1))
    enddo
    data_in=data_out
    endSUBROUTINE filter_1d



    SUBROUTINE filter_2d(nrdata,nzdata,data_in)
    implicit none
    integer::nrdata,nzdata,irdata,izdata
    real*8::data_in(1:nrdata,1:nzdata),data_out(1:nrdata,1:nzdata)

    data_out=data_in
    do irdata=2,nrdata-1
        do izdata=2,nzdata-1
            data_out(irdata,izdata)=1./6.*(2.*data_in(irdata,izdata)+data_in(irdata-1,izdata)+&
                & data_in(irdata+1,izdata)+data_in(irdata,izdata-1)+data_in(irdata,izdata+1))
        enddo
    enddo
    irdata=1 ;     data_out(irdata,2:nzdata-1)=0.25*data_in(irdata,1:nzdata-2)+&
        &        .5*data_in(irdata,2:nzdata-1)+0.25*data_in(irdata,3:nzdata)
    irdata=nrdata; data_out(irdata,2:nzdata-1)=0.25*data_in(irdata,1:nzdata-2)+&
        &       0.5*data_in(irdata,2:nzdata-1)+0.25*data_in(irdata,3:nzdata)
    izdata=1 ;     data_out(2:nrdata-1,izdata)=0.25*data_in(1:nrdata-2,izdata)+&
        &       0.5*data_in(2:nrdata-1,izdata)+0.25*data_in(3:nrdata,izdata)
    izdata=nzdata; data_out(2:nrdata-1,izdata)=0.25*data_in(1:nrdata-2,izdata)+&
        &       0.5*data_in(2:nrdata-1,izdata)+0.25*data_in(3:nrdata,izdata)

    data_in(1:nrdata,1:nzdata)=data_out(1:nrdata,1:nzdata)
    endSUBROUTINE filter_2d



    !SUBROUTINE interp_1d(n_bef,n_aft,x_bef,x_aft,data_bef,data_aft)
    !implicit none
    !integer::n_bef,n_aft
    !real*8::x_bef(1:n_bef),x_aft(1:n_aft),data_bef(1:n_bef),data_aft(1:n_aft)
    !integer::I_x,ix_aft,ix_bef1,ix_bef2
    !real*8::s1,s2,dx_bef,min_x,xtp,dx_tp1,dx_tp2
    !
    !dx_bef=x_bef(2)-x_bef(1)
    !do ix_aft=1,n_aft
    !    xtp=x_aft(ix_aft)
    !    I_x=minloc(abs(xtp-x_bef),1);
    !    min_x=xtp-x_bef(I_x);
    !    if (min_x<0)then
    !        ix_bef1=I_x-1;
    !        ix_bef2=I_x;
    !    else
    !        ix_bef1=I_x;
    !        ix_bef2=I_x+1;
    !    endif
    !    if(ix_bef1<1)ix_bef1=1
    !    if(ix_bef2>n_bef)ix_bef2=n_bef
    !    dx_tp1=xtp-x_bef(ix_bef1);
    !    dx_tp2=x_bef(ix_bef2)-xtp;
    !    s1=dx_tp1/dx_bef;
    !    s2=1.-s1;
    !    data_aft(ix_aft)=s2*data_bef(ix_bef1)+s1*data_bef(ix_bef2);
    !enddo
    !continue
    !endSUBROUTINE interp_1d



    !SUBROUTINE interp_0d(n_bef,x_bef,x_aft,data_bef,data_aft)
    !implicit none
    !integer::n_bef
    !real*8::x_bef(1:n_bef),x_aft,data_bef(1:n_bef),data_aft
    !integer::I_x,ix_bef1,ix_bef2
    !real*8::s1,s2,dx_bef,min_x,xtp,dx_tp1,dx_tp2
    !dx_bef=x_bef(2)-x_bef(1)
    !
    !xtp=x_aft
    !I_x=minloc(abs(xtp-x_bef),1);
    !min_x=xtp-x_bef(I_x);
    !if (min_x<0)then
    !    ix_bef1=I_x-1;
    !    ix_bef2=I_x;
    !else
    !    ix_bef1=I_x;
    !    ix_bef2=I_x+1;
    !endif
    !if(ix_bef1<1)ix_bef1=1
    !if(ix_bef2>n_bef)ix_bef2=n_bef
    !dx_tp1=xtp-x_bef(ix_bef1);
    !dx_tp2=x_bef(ix_bef2)-xtp;
    !s1=dx_tp1/dx_bef;
    !s2=1.-s1;
    !data_aft=s2*data_bef(ix_bef1)+s1*data_bef(ix_bef2);
    !continue
    !endSUBROUTINE interp_0d