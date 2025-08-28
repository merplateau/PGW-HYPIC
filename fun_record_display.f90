

    SUBROUTINE display_main
    USE the_whole_varibles
    implicit none
    integer,parameter::n6=6
    real*8::r_max(1:n6)
    call find_func_cputime_1_of_2  !-------------------1/2

    if((mod(t,t_show)<dt .or. it==1) .and. (iswitch_display/=0))then
        call end_ftime
        tp1=t/td;
        write(*,*)'------------------------------------------------------------------------'
        write(*, "('--------- The start time is ', i4,'.',i2,'.',i2,' ',i2.2,':',i2.2,'-------')") date_time(1:3,1),date_time(5:6,1)
        write(*, "('--------- The time now is  ', i4,'.',i2,'.',i2,' ',i2.2,':',i2.2,'-------')") date_time(1:3,2),date_time(5:6,2)
        write(*, "('--- The calculation time is ', f8.2,'  min (',f8.2,' , hour)----')") time_run,time_run/60.
        write(*, "(' Program is running at', f12.2,3x,'%   (estimating',f12.2,' hours left )')") tp1*100, time_run/60./tp1*(1-tp1)
        write(*,*)
        write(*, "(' t=', f7.1,'  trf   (',es9.2,' s)')") t/trf,t
        write(*,*)
        write(*, "(' nr, nz=', 2i5)") nr,nz
        write(*, "('mass_number=', f7.1,', f=',f7.2,' MHz')") mass_number ,frequency/1e6
        write(*, "('center B0 =', f7.0,', B_resonance=',f7.1,' G')") 1e4*bz_max, 1e4*B_resonance
        write(*, "('T_i =', es10.2,' eV,  ','T_e =', es10.2,' eV')") ti_ave,te_ave
        write(*, "('<Ek_ion> =', es10.2,' eV,  ','<Ek_ele> =', es10.2,' eV')") Ek_i_ave, Ek_e_ave
        write(*, "('power_Joule, irf, resistance = ', es10.2,' W, ',f7.1,' A, ',es10.2,' ohm')")power_Joule ,i_now,resistance
        write(*, "('total loss power =', es10.2,' kW')") sum(power_loss_ave(1:4))/1e3
        !write(*, "('average kinetic energy of the lost ions =', es10.2,' eV')") Ek_loss_ave(4)
        !write(*, "('propulsive efficiency =', f7.1,' %')") 100*(power_loss_ave(2)+power_loss_ave(4))/(sum(power_loss_ave(1:4))+1e-10)

        write(*,*) '   trf,     dt,    dt_mcc,   tau_ii,   td'
        write(*, "(5es9.1, ' s')")   trf,   dt,   dt_mcc,   tau_ii,   td
        write(*,*) 'tau_ii,   tau_ie,  tau_ee,   tau_ei'
        write(*, "(4es9.1, ' s')")tau_ii,tau_ie,tau_ee,tau_ei
        write(*, "(' np_max=',i7,' np=',i7,',n_macro=',es10.2,', ni_ave_ratio=',es10.2)")np_max,np,n_macro,ni_ave_ratio
        write(*, "(' ni_max=',es9.2,', ni_min=',es9.2)")maxval(density_2D),minval(density_2D)

        write(*, "(' inject, loss number', 2i5)") N_inject_particles,N_lost_particles_showup
        write(*, "('function of MCC  has been called  =', f10.1,' times')")run_mcc_times
        write(*, "('function of FDFD has been called  =', f10.1,' times')") real(maxwellcount)
        write(*,*)
        tp1=sum(func_time(1:8))/100.;
        write(*, "('cputime of fuction 1 (particles_inject)      =', f7.2,' %')")func_time(1)/tp1
        write(*, "('cputime of fuction 2 (escape)                =', f7.2,' %')")func_time(2)/tp1
        write(*, "('cputime of fuction 3 (mover->density_and_Es) =', f7.2,' %')")func_time(3)/tp1
        write(*, "('cputime of fuction 4 (mover->push_RK4)       =', f7.2,' %')")func_time(4)/tp1
        write(*, "('cputime of fuction 5 (mcc)                   =', f7.2,' %')")func_time(5)/tp1
        write(*, "('cputime of fuction 6 (record_profiles)       =', f7.2,' %')")func_time(6)/tp1
        write(*, "('cputime of fuction 7 (display)               =', f7.2,' %')")func_time(7)/tp1
        write(*, "('cputime of fuction 8 (FDFD)                  =', f7.2,' %')")func_time(8)/tp1
        write(*,*)
        write(*,*)
    endif
    call find_func_cputime_2_of_2(func_time(7))  !-----2/2
    END SUBROUTINE display_main



    SUBROUTINE rec_time_ave
    USE the_whole_varibles
    implicit none
    integer,parameter::nz_len=6
    integer::k1,k2
    INTEGER*4::ip2 ,ip_acceleration
    Real*8::ek_x,ek_y,ek_z,ek_acceleration,zb3,z_count,dzn,zk1,zk2
    real*8::ek_ave(1:nz_len),ek_tol(1:nz_len),density(1:nz_len),z_c(1:nz_len)
302 format(<nrec_t>(e14.7,' '))

    if(mod(t,t_rec)<dt .or. it==1)then
        z_count=zl/nz_len
        dzn=zl/(nz_len-1)
        do iz=1,nz_len
            z_c(iz)=(iz-1)*dzn
        enddo
        density=0;
        ek_ave=0.;
        ek_tol=0.;
        do ip=1,np
            k1=int((x(ip,3)-z_count)/dzn)+1
            k2=int((x(ip,3)+z_count)/dzn)
            if(k1<0.5)k1=1
            if(k2>nz_len)k2=nz_len
            density(k1:k2)=density(k1:k2)+1
            tp1=mass_q_i_05*v(ip,1)**2
            tp2=mass_q_i_05*v(ip,2)**2
            tp3=mass_q_i_05*v(ip,3)**2
            ek_tol(k1:k2)=ek_tol(k1:k2)+tp1+tp2+tp3
        enddo
        ek_ave(:)=ek_tol(:)/(1e-5+density)

        ek_x=mass_q_i_05*sum(v(1:np,1)**2)/(np+1e-5); !eV
        ek_y=mass_q_i_05*sum(v(1:np,2)**2)/(np+1e-5); !eV
        ek_z=mass_q_i_05*sum(v(1:np,3)**2)/(np+1e-5); !eV
        Ek_e_ave=mass_q_e_05*sum(v_e(1:np,1:3)**2)/real(np); !ev

        ek_acceleration=0
        ip_acceleration=0
        zb3=1.6
        do ip=1,np
            if (x(ip,3)>zb3)then
                ip_acceleration=ip_acceleration+1;
                ek_acceleration=ek_acceleration+0.5/q_mass*(v(ip,1)**2+v(ip,2)**2+v(ip,3)**2);
            endif
        enddo
        ek_acceleration=ek_acceleration/(ip_acceleration+1e-5);

        rec_t(1)=t
        rec_t(2)=real(it)
        rec_t(3)=np
        rec_t(4)=ek_x
        rec_t(5)=ek_y
        rec_t(6)=ek_z
        rec_t(7)=ip_acceleration
        rec_t(8)=ek_acceleration
        rec_t(9:14)=ek_ave(:)
        rec_t(15:20)=density(:)
        rec_t(21:24)=power_loss_ave(1:4)
        rec_t(25:28)=Ek_loss_ave(1:4)
        rec_t(29)=Ek_e_ave

        write (42,302)rec_t
    endif
    END SUBROUTINE rec_time_ave



    SUBROUTINE rec_time_trajectory
    USE the_whole_varibles
    implicit none
303 format(<nrec_p>(e14.7,' '))

    if(mod(t,t_rec_trajectory)<dt )then
        ip=1
        rec_p(1)=t/trf
        rec_p(2:4)=x(ip,1:3)
        rec_p(5:7)=v(ip,1:3)
        rec_p(8:11)=b_rec(1:4)
        rec_p(12:14)=e_rec(1:3)
        ip=2
        rec_p(15:17)=x(ip,1:3)
        rec_p(18:20)=v(ip,1:3)
        write (43,303)rec_p
    endif
    END SUBROUTINE rec_time_trajectory



    SUBROUTINE record_profiles
    USE the_whole_varibles
    implicit none
    real*8::Es_z(1:nz),density_in_z(1:nz),ek_ave(1:nz,1:3) !v_abs(1:np),z_p(1:np),ion_r(1:np),
    character*30 fname
    call find_func_cputime_1_of_2  !-------------------1/2
    if(mod(t,t_rec_pic)<dt .or. it==1)then
        
        call escape
        !call escape_periodicity !energy conservation
        
300     format(<nr>(e12.5,' '))
        !300     format(<nr>(e20.12,' '))
301     format(<nr>(e20.12,' '))
        !---------------------------record the whole data at specific time---------------------------!
        write (fname,120)ifig
120     format('plasma_',i0,'.dat')
        open (unit=20,file=fname,status='unknown',iostat=ierror)
        write (20,300)real(i_plasma_region(1:nr,1:nz))
        write (20,300)ne(1:nr,1:nz)
        write (20,300)te(1:nr,1:nz)
        write (20,300)power_depo(1:nr,1:nz,1)+power_depo(1:nr,1:nz,2)+power_depo(1:nr,1:nz,3)
        close (20)

        write (fname,121)ifig
121     format('rf_field_',i0,'.dat')
        open (unit=21,file=fname,status='unknown',iostat=ierror)
        write (21,300)abs(Erf_PIC(1:nr,1:nz,1))
        write (21,300)abs(Erf_PIC(1:nr,1:nz,2))
        write (21,300)abs(Erf_PIC(1:nr,1:nz,3))
        write (21,300)B_rf(1:nr,1:nz,1)
        write (21,300)B_rf(1:nr,1:nz,2)
        write (21,300)B_rf(1:nr,1:nz,3)
        write (21,300)real(Erf_PIC(1:nr,1:nz,1))
        write (21,300)imag(Erf_PIC(1:nr,1:nz,1))
        write (21,300)real(Erf_PIC(1:nr,1:nz,2))
        write (21,300)imag(Erf_PIC(1:nr,1:nz,2))
        write (21,300)real(Erf_PIC(1:nr,1:nz,3))
        write (21,300)imag(Erf_PIC(1:nr,1:nz,3))
        close (21)


        do iz=1,nz
            Es_z(iz)=sum(Es_2D(:,iz,2))/nr
            !density_in_z(iz)=sum(density_2D(:,iz))/nr
        enddo
        !z_p(1:np)=x(1:np,3)
        !ion_r(1:np)=sqrt(x(1:np,1)**2+x(1:np,2)**2 )
        !ion_rzv(1:np,2)=x(1:np,3)
        !ion_rzv(1:np,3:5)=v(1:np,1:3)
        call find_Ek_1D(density_in_z,ek_ave)
        ek_ave=mass_q_i_05*ek_ave

322     format(<nz>(e12.5,' '))
122     format('Ek_ave_z_',i0,'.dat')
        write (fname,122)ifig
        open (unit=22,file=fname,status='unknown',iostat=ierror)
        write (22,322)z
        write (22,322)density_in_z
        write (22,322)ek_ave
        write (22,322)Es_z
        close(22)

323     format(<np>(e12.5,' '))
123     format('ion_rz_',i0,'.dat')
        write (fname,123)ifig
        open (unit=23,file=fname,status='unknown',iostat=ierror)
        write (23,323)sqrt(x(1:np,1)**2+x(1:np,2)**2 )!ion_r(1:np)
        write (23,323)x(1:np,3)
        write (23,323)v(1:np,1)
        write (23,323)v(1:np,2)
        write (23,323)v(1:np,3)
        close(23)

124     format('density_Es_2D_',i0,'.dat')
        write (fname,124)ifig
        open (unit=24,file=fname,status='unknown',iostat=ierror)
        write (24,300)density_2D
        write (24,300)Es_2D(:,:,1)
        write (24,300)Es_2D(:,:,2)
        write (24,300)Ek_ion_2D
        write (24,300)te_2D
        close(24)


        call rec_para
        ifig=ifig+1
    endif
    call find_func_cputime_2_of_2(func_time(6))  !-----2/2
    END SUBROUTINE record_profiles



    SUBROUTINE find_Ek_1D(density_x,ek_x) !,x_p,v_innz,np,z,
    USE the_whole_varibles
    implicit none
    !integer::np,nzip,
    integer::ix1,ix2,I_x,i_comp
    real*8::ek_x(1:nz,1:3),density_x(1:nz) !,x_p(np),v_in(1:np,1:3),z(1:nz)dz,
    real*8::xtp,s1,s2,min_x,dx_tp1,dx_tp2
    ek_x=0.
    density_x=1e-5

    do ip=1,np
        xtp=x(ip,3)
        I_x=minloc(abs(xtp-z),1);
        min_x=xtp-z(I_x);
        if (min_x<0)then
            ix1=I_x-1;
            ix2=I_x;
        else
            ix1=I_x;
            ix2=I_x+1;
        endif
        if(ix1<1)ix1=1
        if(ix2>nz)ix2=nz
        dx_tp1=abs(xtp-z(ix1)); !note abs()
        s1=dx_tp1/dz;
        s2=1.-s1;
        ek_x(ix1,1:3)=ek_x(ix1,1:3)+s2*v(ip,1:3)**2
        ek_x(ix2,1:3)=ek_x(ix2,1:3)+s1*v(ip,1:3)**2
        density_x(ix1)=density_x(ix1)+s2
        density_x(ix2)=density_x(ix2)+s1
    enddo
    do ix1=1,3
        ek_x(:,ix1)=ek_x(:,ix1)/density_x
    enddo
    density_x=density_x/(pi*r_ion_max**2*dz)
    END SUBROUTINE find_Ek_1D



    SUBROUTINE start_ftime
    USE the_whole_varibles
    implicit none
    call date_and_time(timefunc_char(1,1),timefunc_char(2,1), timefunc_char(3,1), timefunc_int(:,1))
    endSUBROUTINE start_ftime



    SUBROUTINE end_ftime
    USE the_whole_varibles
    implicit none
    call date_and_time(timefunc_char(1,2),timefunc_char(2,2), timefunc_char(3,2), timefunc_int(:,2))
    timefunc_int_real=real(timefunc_int)
    tp1=(timefunc_int_real(2,1)-1)*30.*24.*60.*60.+(timefunc_int_real(3,1)-1)*24.*60.*60.+timefunc_int_real(5,1)*60.*60.+ &
        &                timefunc_int_real(6,1)*60.+timefunc_int_real(7,1)+timefunc_int_real(8,1)*1e-3
    tp2=(timefunc_int_real(2,2)-1)*30.*24.*60.*60.+(timefunc_int_real(3,2)-1)*24.*60.*60.+timefunc_int_real(5,2)*60.*60.+ &
        &                timefunc_int_real(6,2)*60.+timefunc_int_real(7,2)+timefunc_int_real(8,2)*1e-3
    time_run=(tp2-tp1)/60.
    date_time=timefunc_int
    END SUBROUTINE



    SUBROUTINE find_func_cputime_1_of_2
    USE the_whole_varibles
    implicit none
    call CPU_TIME(time_diff(1))
    end SUBROUTINE

    SUBROUTINE find_func_cputime_2_of_2(time_used)
    USE the_whole_varibles
    implicit none
    real*8 time_used
    call CPU_TIME(time_diff(2))
    time_used=time_used+time_diff(2)-time_diff(1)
    end SUBROUTINE find_func_cputime_2_of_2
