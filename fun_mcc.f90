    subroutine mcc_main
    use the_whole_varibles
    implicit none


    call find_func_cputime_1_of_2  !-------------------1/2
    call mcc_tau

    if(abs(t-time_to_run_mcc)<1.001*dt .or. t<1.5*dt )then
        time_to_run_mcc=time_to_run_mcc+dt_mcc
        run_mcc_times=run_mcc_times+1 !test!

        call mcc_ii !i-i collision
        call mcc_ee !e-e collision
        call mcc_ie !i-e collision
        call mcc_ei !e-i collision

        Ek_i_ave=mass_q_i_05*sum(v(1:np,1:3)**2)/real(np); !ev
        Ek_e_ave=mass_q_e_05*sum(v_e(1:np,1:3)**2)/real(np); !ev
        ti_ave=Ek_i_ave/1.5
        te_ave=Ek_e_ave/1.5
    endif
    call find_func_cputime_2_of_2(func_time(5))  !-----2/2
    END SUBROUTINE mcc_main



    subroutine mcc_tau
    use the_whole_varibles
    implicit none
    real*8 tau_all(4),tau_min
    real*8 qa,qb,ma,mb,g2_ave,b0_ave,lambda_D,t_ave

    ne_mcc=density_ave
    Ek_i_ave=mass_q_i_05*sum(v(1:np,1:3)**2)/real(np); !ev
    Ek_e_ave=mass_q_e_05*sum(v_e(1:np,1:3)**2)/real(np); !ev
    ti_ave=Ek_i_ave/1.5
    te_ave=Ek_e_ave/1.5
    t_ave=.5*(ti_ave+te_ave)


    ! i-i collision
    qa=qi;
    qb=qi;
    ma=mi;
    mb=mi;
    mu_ab=ma*mb/(ma+mb);
    g2_ave=3.*t_ave*qa/mu_ab;
    b0_ave=qa*qb/(2*pi*epson0*mu_ab*g2_ave);
    lambda_D=sqrt(epson0*t_ave*qa/(ne_mcc*qa**2));
    ln_A=log(lambda_D/b0_ave)
    tau_ii=8*pi*sqrt(2.*ma)*(epson0**2)*((qe*ti_ave)**1.5)/(ne_mcc*(qe**4)*ln_A);
    !v_ii=1/tau_ii

    ! e-e collision
    qa=qe;
    qb=qe;
    ma=me;
    mb=me;
    mu_ab=ma*mb/(ma+mb);
    g2_ave=3.*t_ave*qa/mu_ab;
    b0_ave=qa*qb/(2*pi*epson0*mu_ab*g2_ave);
    lambda_D=sqrt(epson0*t_ave*qa/(ne_mcc*qa**2));
    ln_A=log(lambda_D/b0_ave)
    tau_ee=8*pi*sqrt(2.*ma)*(epson0**2)*((qe*te_ave)**1.5)/(ne_mcc*(qe**4)*ln_A);
    !v_ee=1/tau_ee

    ! i-e collision
    qa=qi;
    qb=qe;
    ma=mi;
    mb=me;
    mu_ab=ma*mb/(ma+mb);
    g2_ave=3.*t_ave*qa/mu_ab;
    b0_ave=qa*qb/(2*pi*epson0*mu_ab*g2_ave);
    lambda_D=sqrt(epson0*t_ave*qa/(ne_mcc*qa**2));
    ln_A=log(lambda_D/b0_ave)
    !tau_ie=8*pi*sqrt(2.*ma)*(epson0**2)*((qe*t_ave)**1.5)/(ne_mcc*(qe**4)*ln_A);
    tau_ie=(4.*pi*epson0)**2*3*ma*((qe*te_ave)**1.5)/(8*sqrt(2*pi*mb)*ne_mcc*(qe**4)*ln_A);
    !v_ee=1/tau_ee

    ! e-i collision
    qa=qe;
    qb=qi;
    ma=me;
    mb=mi;
    mu_ab=ma*mb/(ma+mb);
    g2_ave=3.*t_ave*qa/mu_ab;
    b0_ave=qa*qb/(2*pi*epson0*mu_ab*g2_ave);
    lambda_D=sqrt(epson0*t_ave*qa/(ne_mcc*qa**2));
    ln_A=log(lambda_D/b0_ave)
    tau_ei=         4.*pi*(epson0**2)*sqrt(mu_ab)*((qe*t_ave)**1.5)/(ne_mcc*(qe**4)*ln_A);
    !tau_ei=sqrt(pi)*4.*pi*(epson0**2)*sqrt(ma)*((qe*t_ave)**1.5)/(ne_mcc*(qe**4)*ln_A);
    !v_ee=1/tau_ee

    tau_all(1)=tau_ii;
    tau_all(2)=tau_ee;
    tau_all(3)=tau_ie;
    tau_all(4)=tau_ei;
    tau_min=minval(tau_all)
    !tau_min=tau_ii;

    dt_mcc=5e-2*tau_min  !run mover (dt) 50 times and then at least run mcc 1 time
    if( (dt_mcc*1e-3)>dt)dt_mcc=1e-3*tau_min  !run mover (dt) 50 times and then at least run mcc 1 time
    if(dt>dt_mcc)dt_mcc=dt

    !dt_mcc=dt  !test!!
    continue
    END SUBROUTINE mcc_tau



    subroutine mcc_ii
    use the_whole_varibles
    implicit none
    integer*4::n_rand(1:np),np_d2,i1,i2,itp_3
    real*8 qa,qb,ma,mb
    real*8 ::vv_inj_bef(3),vv_tar_bef(3),vv_inj_aft(3),vv_tar_aft(3)

    qa=qi;
    qb=qi;
    ma=mi;
    mb=mi;
    mu_ab=ma*mb/(ma+mb);
    sab_part=ne_mcc*ln_A/(4*pi)*dt_mcc*(qa*qb/(epson0*mu_ab))**2;
    if(mod(np,2)/=0) then
        np_d2=nint((np-1.)/2.)
    else
        np_d2=nint(np/2.)
    endif
    call mcc_randperm(np,n_rand,0)

    do ip=1,np_d2
        i1=n_rand(ip);
        i2=n_rand(ip+np_d2);
        vv_inj_bef=v(i1,1:3)
        vv_tar_bef=v(i2,1:3)
        call mcc_sub(ma,mb,qa,qb,vv_inj_bef,vv_tar_bef,vv_inj_aft,vv_tar_aft)
        v(i1,1:3)=vv_inj_aft
        v(i2,1:3)=vv_tar_aft
    enddo

    endsubroutine mcc_ii


    subroutine mcc_ee
    use the_whole_varibles
    implicit none
    integer*4::n_rand(1:np),np_d2,i1,i2,itp_3
    real*8 qa,qb,ma,mb
    real*8 ::vv_inj_bef(3),vv_tar_bef(3),vv_inj_aft(3),vv_tar_aft(3)

    qa=qe;
    qb=qe;
    ma=me;
    mb=me;
    mu_ab=ma*mb/(ma+mb);
    sab_part=ne_mcc*ln_A/(4*pi)*dt_mcc*(qa*qb/(epson0*mu_ab))**2;
    if(mod(np,2)/=0) then
        np_d2=nint((np-1.)/2.)
    else
        np_d2=nint(np/2.)
    endif
    call mcc_randperm(np,n_rand,0)

    do ip=1,np_d2
        i1=n_rand(ip);
        i2=n_rand(ip+np_d2);
        vv_inj_bef=v_e(i1,1:3)
        vv_tar_bef=v_e(i2,1:3)
        call mcc_sub(ma,mb,qa,qb,vv_inj_bef,vv_tar_bef,vv_inj_aft,vv_tar_aft)
        v_e(i1,1:3)=vv_inj_aft
        v_e(i2,1:3)=vv_tar_aft
    enddo
    endsubroutine mcc_ee



    subroutine mcc_ie
    use the_whole_varibles
    implicit none
    integer*4 n_rand1(1:np),n_rand2(1:np),np_d2,i1,i2,itp_3
    real*8 qa,qb,ma,mb
    real*8    vv_inj_bef(3),vv_tar_bef(3),vv_inj_aft(3),vv_tar_aft(3)

    qa=qi;
    qb=qe;
    ma=mi;
    mb=me;
    mu_ab=ma*mb/(ma+mb);
    sab_part=ne_mcc*ln_A/(4*pi)*dt_mcc*(qa*qb/(epson0*mu_ab))**2;
    call mcc_randperm(np,n_rand1,1)
    call mcc_randperm(np,n_rand2,1)

    do ip=1,np
        i1=n_rand1(ip);
        i2=n_rand2(ip);
        vv_inj_bef=v(i1,1:3)
        vv_tar_bef=v_e(i2,1:3)
        call mcc_sub(ma,mb,qa,qb,vv_inj_bef,vv_tar_bef,vv_inj_aft,vv_tar_aft)
        v(i1,1:3)=vv_inj_aft
        v_e(i2,1:3)=vv_tar_aft
    enddo
    endsubroutine mcc_ie



    subroutine mcc_ei
    use the_whole_varibles
    implicit none
    integer*4 n_rand1(1:np),n_rand2(1:np),np_d2,i1,i2,itp_3
    real*8 qa,qb,ma,mb,v_ref,t_ref
    real*8    vv_inj_bef(3),vv_tar_bef(3),vv_inj_aft(3),vv_tar_aft(3)

    qa=qe;
    qb=qi;
    ma=me;
    mb=mi;
    mu_ab=ma*mb/(ma+mb);
    v_ref=sqrt(qa*.5*(ti_ave+te_ave)/mu_ab)
    t_ref=tau_ei;
    sab_part=v_ref**3*dt_mcc/t_ref
    !sab_part=ne_mcc*ln_A/(4*pi)*dt_mcc*(qa*qb/(epson0*mu_ab))**2;
    call mcc_randperm(np,n_rand1,1)
    call mcc_randperm(np,n_rand2,1)

    do ip=1,np
        i1=n_rand1(ip);
        i2=n_rand2(ip);
        vv_inj_bef=v_e(i1,1:3)
        vv_tar_bef=v(i2,1:3)
        call mcc_sub(ma,mb,qa,qb,vv_inj_bef,vv_tar_bef,vv_inj_aft,vv_tar_aft)
        v_e(i1,1:3)=vv_inj_aft
        v(i2,1:3)=vv_tar_aft
    enddo
    endsubroutine mcc_ei


    subroutine mcc_randperm(ncc,n_rand,itype)
    use the_whole_varibles
    implicit none
    integer*4::ncc,ncc_d2,ilen,n_len,itp2,itype
    integer*4::n_bef(1:ncc),n_rand(1:ncc)
    real*8::rand_0_1
    integer*4 itp0

    if(itype==0) then
        ! i-i or e-e grouping
        if(mod(ncc,2)/=0) then
            n_len=ncc-1
        else
            n_len=ncc
        endif
        ncc_d2=n_len/2

        do itp0=1,n_len
            n_bef(itp0)=itp0
        enddo
        n_rand(1:n_len)=0
        do itp0=1,n_len
            ilen=n_len-itp0+1;
            call random_number(rand_0_1)
            itp2=floor(ilen*rand_0_1)+1;
            n_rand(itp0)=n_bef(itp2);
            n_bef(itp2)=n_bef(ilen);
        enddo
        n_rand(ncc)=ncc

    else

        ! i-e or e-i grouping
        n_len=ncc;
        do itp0=1,n_len
            n_bef(itp0)=itp0
        enddo
        n_rand(1:n_len)=0
        do itp0=1,n_len
            ilen=n_len-itp0+1;
            call random_number(rand_0_1)
            itp2=floor(ilen*rand_0_1)+1;
            n_rand(itp0)=n_bef(itp2);
            n_bef(itp2)=n_bef(ilen);
        enddo

    endif
    END SUBROUTINE



    subroutine mcc_sub(ma,mb,qa,qb,v_inj_bef,v_tar_bef,v_inj_aft,v_tar_aft)
    use the_whole_varibles
    implicit none
    Real*8 ::ma,mb,qa,qb
    real*8 ::v_inj_bef(3),v_tar_bef(3),v_inj_aft(3),v_tar_aft(3)
    !Real*8 ::ne_mcc,mu_ab
    Real*8 ::rand_0_1,s_ab,a_ab
    Real*8 ::g_ab,ep_ab,cos_ep,sin_ep,g_perp,cosx,sinx
    Real*8 ::tp_min,u_ab
    Real*8 ::tp_ab(1:3)
    real*8 ::g3_ab(1:3),h3_ab(1:3)
    integer*4 itp0

    g3_ab(1:3)=v_inj_bef(1:3)-v_tar_bef(1:3)+1e-5; !+1e-5 to avoid a3_ab=0  , very important!
    !need know g_ab=va-vb
    g_ab=sqrt(g3_ab(1)**2+g3_ab(2)**2+g3_ab(3)**2);
    s_ab=sab_part*g_ab**(-3); !find s

    !------find A--------------!
    if (s_ab<0.01) then
        a_ab=1/s_ab;
    elseif (s_ab<4) then
        !            [tp_min,itp0]
        tp_min=minval(abs(s_ab-s_tab));
        itp0=minloc(abs(s_ab-s_tab),1)
        if (s_ab-s_tab(itp0)>0)then
            tp1=s_ab-s_tab(itp0);
            tp2=s_tab(itp0+1)-s_ab;
            a_ab=(tp1*a_tab(itp0+1)+tp2*a_tab(itp0))/(tp1+tp2);
        else
            tp1=s_tab(itp0)-s_ab;
            tp2=s_ab-s_tab(itp0-1);
            a_ab=(tp1*a_tab(itp0-1)+tp2*a_tab(itp0))/(tp1+tp2);
        endif
    else
        a_ab=3*exp(-s_ab);
    endif
    !------find A--------------!

    !find cosx
    call random_number(rand_0_1)
    u_ab=rand_0_1
    if (s_ab<0.01) then
        cosx=1+s_ab*log(u_ab);
    elseif (s_ab<6) then
        cosx=1/a_ab*log(exp(-a_ab)+2*u_ab*sinh(a_ab));
    else
        cosx=2*u_ab-1;
    endif
    sinx=sqrt(1-cosx**2);

    !find v after collision
    call random_number(rand_0_1)
    ep_ab=2*pi*rand_0_1
    cos_ep=cos(ep_ab);
    sin_ep=sin(ep_ab);
    g_perp=sqrt(g3_ab(2)**2+g3_ab(3)**2);
    h3_ab(1)=g_perp*cos_ep;
    h3_ab(2)=-(g3_ab(2)*g3_ab(1)*cos_ep+g_ab*g3_ab(3)*sin_ep)/g_perp;
    h3_ab(3)=-(g3_ab(3)*g3_ab(1)*cos_ep-g_ab*g3_ab(2)*sin_ep)/g_perp;
    tp_ab(1:3)=(g3_ab*(1-cosx)+h3_ab*sinx);

    v_inj_aft=v_inj_bef-mb/(ma+mb)*tp_ab;
    v_tar_aft=v_tar_bef+ma/(ma+mb)*tp_ab;
    END SUBROUTINE mcc_sub



    subroutine mcc_constant
    use the_whole_varibles
    implicit none
    s_tab(1:22)=(/0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09,  &
        &    0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 2., 3., 4./);
    a_tab(1:22)=(/100.5, 50.5, 33.84, 25.50, 20.50, 17.17, 14.79, 13.01, 11.62, 10.51, 5.516,  &
        &    3.845, 2.987, 2.448, 2.067, 1.779, 1.551, 1.363, 1.207, 0.4105, 0.1496 ,0.05496/);
    END SUBROUTINE
