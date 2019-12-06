!##########################################################################
        module vars
!##########################################################################
	  SAVE
        double precision g_dx,g_dy,g_dz,dens,Re,eps,fac,Pr,beta,Th,Tc 
        double precision qzero,qstpn,forcn,gx,gy,gz,Sc_t
        double precision flomas,rmax,alfapr,resor,fric,TI_SEM
        double precision ctime,dt,dtavg,dtsum,safety_factor,noise,Mdef
        double precision t_start_averaging1,t_start_averaging2,ubulk
	  double precision xst,xen,yst,yen,zst,zen
        integer alfabc,niter,nswp(4),nsweep,iter,nsweeps,ntime
        integer ngg,nge,ngc,n_unstpt,OMP_threads,iaddinlet,ireadinlet
        integer ngrid_input,ngrd_gl,maxcy,iproln,irestr
        integer itmax,sweeps,numfile,numfile1,numfile2
        integer conv_sch,diff_sch,sgs_model,solver,mg_itrsch
        integer bc_w,bc_e,bc_s,bc_n,bc_b,bc_t,UPROF_SEM
        integer Tbc_w,Tbc_e,Tbc_s,Tbc_n,Tbc_b,Tbc_t
        integer itime,itime_start,itime_end,n_out,ITMAX_PI
        integer ipref,jpref,kpref,prefdom,ITMAX_SEM,NE_SEM
	  integer pl,pl_ex,differencing,LMR,normal_inter,order
        character*80 keyword,L_n
	  logical :: LRESTART,LIMB,SGS,PERIODIC,LENERGY,LROUGH
	  logical :: pressureforce,time_averaging,reinitmean
	  logical :: LPT,save_inflow,read_inflow,L_dt,LSCALAR
!============================== LSM VARIABLES =============================
        double precision reldif_LSM,length,densl,densg,nul,nug,mul,mug
        double precision cfl_lsm
        double precision uprev
	  double precision grx,gry,grz
	  double precision slope
        double precision densip12,densjp12,denskp12
	  double precision densim12,densjm12,denskm12
        double precision muip12,mujp12,mukp12,muim12,mujm12,mukm12
	  double precision Mdef_w
        integer ntime_reinit,accuracy
	  integer ngrid
	  integer numfile3
	  logical :: L_LSM,REINIT,L_LSMbase,L_LSMinit,LENDS
	  logical :: L_anim_phi,L_anim_grd
!=========================================================================
        end module
!##########################################################################
