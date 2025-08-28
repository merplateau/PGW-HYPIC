    !hypic code are developed primarily by Wu Mingyang, with fdfd developed by Xu Andong and supervised by Xiao Chijie.
    !started Date:2023.7.12
    !contact: 592774972@qq.com, xuandong@pku.edu.cn or ymwu@pku.edu.cn
    !All units are in the International System of Units (SI) unless otherwise specified.

    !compile command
    !ifort hypic.f90 fun_b0.f90 fun_record_display.f90  fun_grid.f90 fun_ini.f90 fun_irf_antenna.f90 fun_particles.f90 fun_mcc.f90 fun_fdfd.f90 -mkl -o 3-1.x


    !------Xuad
    module types
    Type sparse
        integer*4:: nzmax
        integer*4:: n
        integer*4,allocatable :: row(:) !row
        integer*4,allocatable :: col(:) !coulumn
        real*8,allocatable    :: val(:) !value
    end Type sparse
    Type sparse_complex
        integer*4:: nzmax
        integer*4:: n
        integer*4,allocatable  :: row(:) !row
        integer*4,allocatable  :: col(:) !coulumn
        Complex*16,allocatable :: val(:) !value
    end Type sparse_complex
    end module types
    !------Xuad

    !-----------Global Variables-------------!
    module the_whole_varibles
    use types
    integer*4,PARAMETER ::np_ini=100000,np_max=np_ini
    integer*4,parameter :: ns=3,nth=200,n3=3
    real*8,parameter ::frequency=13.56e6,trf=1/frequency; !f:Hz
    !---------------constants, do not change---------------------!
    real*8,parameter::c=3e8, qe=1.6022e-19,qi=qe,qe_abs=qe,  pi = 3.1415926,kB=1.38065e-23;
    !real*8,parameter::me=9.1e-31,mass_number=1.,mi=4*me;
    !real*8,parameter::mass_number=1.,mi=mass_number*1.67D-27,me=mi/4.;
    real*8,parameter::mass_number=1.,mi=mass_number*1.67D-27,me=9.1e-31;
    real*8,parameter:: epson0=8.8542e-12 ,mu0=pi*(4e-7),cm=1e2,cm3=1e6
    real*8,parameter:: mass=mi, q_mass=qi/mi,mass_q_i_05=0.5*mass/qi,mass_q_e_05=0.5*me/qe;
    real*8,parameter::fv_te_limit=1 !eV when te<1eV, set Te=1eV;
    real*8,parameter::sig_en=15e-20,sig_in=1e-18;
    Complex*16,parameter::i=(0.0D0,1.0D0)
    !---------------constants, do not change---------------------!

    !maxwell_FDFD
    type(sparse_complex) :: A_helicon
    real*8  :: dt_run_fdfd,time_call_fdfd=0.
    integer*4:: i3,ith,m,m_start,m_end,m_delta,iswitch_antenna_type,switch_z_boundary
    integer*4:: switch_source_type,switch_plasma_vacuum_boundary,switch_divE
    integer*4:: switch_boundary_half_variables,switch_r_0_boundary, maxwellcount=0
    integer*4:: Region,jieshu,iswitch_output_Region,nzmax_pt=0
    integer*4 iswitch_ji_on
    integer*4, allocatable :: FindRegion(:,:)
    Integer(kind=8) :: count1,count2,count0,count,count_rate,count_max
    real*8 :: om,omc2,ds,time_Maxwell,ptotal,ptotal_inside,ptotal_boundary,Res_p
    real*8, external :: DBesselJZero
    real*8, allocatable :: E_rf(:,:,:),B_rf(:,:,:),jp(:,:,:),ne(:,:),b0_DC(:,:,:),te(:,:)
    real*8, allocatable ::fvm(:,:),ptotm(:,:,:),power_depo(:,:,:),th(:),ptotal_m(:),max_eth(:)
    real*8 :: m_real
    Complex*16:: im,imp,iom,m_max
    Complex*16, allocatable ::eq_b(:),eq_xlast(:,:),ja_m_nor(:,:,:), ja_record(:,:,:,:),ja_AC(:,:,:,:)
    Complex*16, allocatable :: e_m(:,:,:),e_record(:,:,:,:),e_AC(:,:,:,:),e_int(:,:,:)
    Complex*16, allocatable :: jp_m(:,:,:), jp_record(:,:,:,:),jp_AC(:,:,:,:)
    Complex*16, allocatable :: b_m(:,:,:), b_record(:,:,:,:),b_AC(:,:,:,:)
    Complex*16, allocatable ::ep(:,:,:),si(:,:,:,:),peff(:,:), Xkz_e(:,:,:)
    Complex*16, allocatable :: e_output(:,:,:,:),jp_output(:,:,:,:),ptot_output(:,:,:)
    Complex*16, allocatable :: Erf_PIC(:,:,:)
    !maxwell_FDFD


    !initial
    integer*4:: iswitch_display,iswtich_inner_dipole,i_switch_B0
    real*8 :: bz_max,ti_ini,te_ini,pn,power_Joule,resistance,nn,max_density_set,ni_ave_ratio
    real*8 :: i_now,B0_correction_factor,B_resonance,t_power_on
    real*8 :: width_ne_r,width_ne_z,z_peak_ne,width_ne_source_r,width_ne_source_z
    !initial

    !!antenna and grid
    integer*4::nr,nz,n_vac,nra,nza,nrz_diploe(4),i_vac,nz6
    integer*4 :: nzx2,nrd,nzs,nzd,nrp,itp,ir,iz,ierror
    integer*4, allocatable :: nr_met(:),nz_met(:),nr_vac(:),nz_vac(:)
    integer*4, allocatable::i_plasma_region(:,:),i_gamma_center_region(:,:,:)
    real*8 :: rp,zp,d2r,d2z,r_ion_max,rz_dipole(1:4)
    real*8 :: dr,dz,rs,zs,rl,zl,ra,za,r1,dr2,dz2,dr22,dz22,drz,dr_2,dz_2
    real*8 :: tp0,tp1,tp2,tp3,tp4,tp5,tp6,tp7,tp8,tp9,tp10,tp11,tp12
    real*8, allocatable :: r(:),z(:),r2(:)
    !!antenna and grid

    !display and record
    character*10:: time_char(3,2),timefunc_char(3,2) !(:,1)-start time;  (:,2)-end time
    integer*4:: date_time(8,2), timefunc_int(8,2),func_time_count(8,2),i_count=0; !(:,1)-start time;  (:,2)-end time
    real*8 :: timefunc_int_real(8,2),time_run
    real*8 :: t=0,dt,td
    Real*8::func_time(1:10)=0, time_diff(1:2)=0,run_mcc_times=0
    !display and record

    !particles
    integer*4,parameter:: nrec_t=29,nrec_p=20
    integer*4::ir1_iz1_grid(1:np_max,1:2)   !corresponding to x_to_grid
    integer*4 ::nr_ion
    integer*4::ifig=0,i_sign=0,N_inject_particles,N_lost_particles_showup=0,N_null=0
    integer*4::num_inject,N_lost_particles,ip,np
    integer*8::it
    real*8:: vtp(1:3),dr_ion, s_drz
    real*8:: coeff1,t_cyclotron,vi_ex,ve_ex,dt_trf,z_plasma_center,n_macro
    !v: ion velocity; ve: electron velocity
    real*8:: x(1:np_max,1:3),v(1:np_max,1:3),v_e(1:np_max,1:3),dt_np(1:np_max),t_np(1:np_max) !1:3 -- x y z direction
    real*8:: x_to_grid(1:np_max,1:2) !r: s1= x_to_grid(1:np_max,1),s2=1-s1; z: s3= x_to_grid(1:np_max,2),s4=1-s3
    real*8:: rec_t(nrec_t)=0.,rec_p(nrec_p)=0.,e_rec(1:3),b_rec(1:4)
    real*8:: Ek_loss_ave(1:4)=0,power_loss_ave(1:4)=0,Ek_loss_tol(1:4)=0.,num_loss(1:4)=1e-2
    real*8, allocatable :: pdf_ne_r(:),pdf_ne_source_r(:),pdf_ne_z(:),pdf_ne_source_z(:),r_particle_inj(:)
    Real*8:: t_show,t_rec,t_rec_trajectory,z_inject ,r_inject, rz_ant_region(1:4),t_rec_pic
    real*8, allocatable :: density_2D(:,:),Es_2D(:,:,:),Ek_ion_2D(:,:),te_2D(:,:)
    Complex*16::expt
    !particles

    !mcc
    integer*4::iswitch_mcc
    Real*8:: s_tab(1:22),a_tab(1:22)
    Real*8:: ln_A,dt_mcc,time_to_run_mcc=0.,ne_mcc
    Real*8:: mu_ab,sab_part
    Real*8:: density_ave,te_ave,ti_ave,Ek_i_ave,Ek_e_ave,tau_ii,tau_ie,tau_ee,tau_ei
    !mcc

    end module the_whole_varibles



    !--------------------------------------------------------------------------------------!
    !----------------------------main program start----------------------------------------!
    !--------------------------------------------------------------------------------------!
    Program hypic !hypic code, HYbrid Particle-In-Cell Monte-Carlo collision
    use the_whole_varibles
    implicit none

    call start_ftime
    call ini
    if(iswitch_display/=0) write(*,*)'Fun_ini finished. Fun_fdfd is running.'
    call maxwell_FDFD
    call particles_initialization

    it=0;t=0;
    if(iswitch_display/=0) write(*,*)'Main loop is running.'
    do while(1)
        it=it+1; t=t+dt

        if(mod(t,dt_run_fdfd)<dt )then
            call set_ne_Te
            call maxwell_FDFD
        endif

        call particles_inject
        call escape
        !call escape_periodicity 
        call mover
        if(iswitch_mcc/=0)call mcc_main

        call record_profiles
        call find_power_loss
        call rec_time_ave
        call rec_time_trajectory
        call display_main
        if( t>td .or. np<1.5) exit
    enddo

    close(42)
    close(43)
    call display_main
    write(*,*) 'code running has finished !!! '
    End Program hypic
    !--------------------------------------------------------------------------------------!
    !-----------------------------main program end-----------------------------------------!
    !--------------------------------------------------------------------------------------!


