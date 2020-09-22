
PROGRAM NeptuneMHD
         ! previously multifluid
      IMPLICIT NONE
      
      !     this is a 3-D modified three fluid simulation using the
      !           electrons : arrays starting with e
      !           solar wind : arrays starting with q
      !           ionospheric: arrays starting with o oxygen, 
      !                                         h hydrogen 
      !
      !      and change SET33(0.1,1.,0.1,1.0  to
      !                 SET3(0.,1.,0.,1.
      !     WARNING - MAKE SURE SPACE IS COMPATIBLE WITH GRAPHICS
      !               routines
      !               The arrays IJZERO,IJMID  and IJSRF have to
      !               be modified in BSUBS.F if modified in MAIN
      !     ADD in MIRROR DIPOLE  so bx is zero at wind boundary
      !     GRAPHICS - CONTOUR AND FLOWS - HAVE TO BE manually set
      !                for right aspect ratios
      !
      !      WARMING PLASMA and MAGNETIC FIELD data must be aligned in TIME
      !
      !       grid within grid size nt = 2*ngrd
      !       ncts is the size of the data array for the IMF data file
      !
             INTEGER, PARAMETER :: nx=101,ny=101,nz=49,ngrd=6, &
                  nt=12,ncraft=6,ncts=281
      !
      !      graphics parameters:muvwp2=amax(mx,my,mz)+2,mz2=(mz-1)/2+1
             INTEGER, PARAMETER :: mx=51,my=51,mz=25,muvwp2=53,mz2=13
      
            common /space/vvx,vvy,vvz, &
                 tvx,tvy,tvz, &
                 evx,evy,evz
            common /gridding/grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
                       grd_zmin,grd_zmax,xspac,ncore,nbndry, &
                       rx,ry,rz,xdip,ydip,zdip,rearth,b0, &
                       sin_tilt,cos_tilt
      
             REAL, DIMENSION(9) :: grd_xmin, grd_xmax, grd_ymin, grd_ymax, &
                  grd_zmin, grd_zmax, xspac
             INTEGER, DIMENSION(9) :: ncore, nbndry
             REAL :: rx, ry, rz
             REAL :: xdip, ydip, zdip
             REAL :: rearth
             REAL :: b0
             REAL :: sin_tilt, cos_tilt
      
            common /rotation/v_rot,r_rot,re_equiv
      
      
      !
      !      grid limits now set by grd_min grd_max arrays
      !      ncore denotes couser grid to be hollowed out by fine grid
      !      nbndry denotes finer grid to which coaser grid sets flanks
      !      xspac is the relative grid spacing relative to inner grid system
      !      xcraft is the actual position of the spacecraft in RE
      !          4th dimension of the actual time
      !      zcraft is the future position of the spacecraft in RE
      !      rcraft is the position of the spacecraft for which
      !           IMF is reference. NO alteration from boundary conditions applied
      !
            real xcraft(4,ncraft),zcraft(4,ncraft),rcraft(3)
      !
      !     physics plasma quantities
      !
            real bx(nx,ny,nz,ngrd),by(nx,ny,nz,ngrd),bz(nx,ny,nz,ngrd), &
                 qpx(nx,ny,nz,ngrd),qpy(nx,ny,nz,ngrd),qpz(nx,ny,nz,ngrd), &
                 qrho(nx,ny,nz,ngrd),qpres(nx,ny,nz,ngrd), &
                 opx(nx,ny,nz,ngrd),opy(nx,ny,nz,ngrd),opz(nx,ny,nz,ngrd), &
                 orho(nx,ny,nz,ngrd),opres(nx,ny,nz,ngrd), &
                 hpx(nx,ny,nz,ngrd),hpy(nx,ny,nz,ngrd),hpz(nx,ny,nz,ngrd), &
                 hrho(nx,ny,nz,ngrd),hpres(nx,ny,nz,ngrd), &
                 epres(nx,ny,nz,ngrd)
      !
      !     work arrays for runge-kutta and smothing
      !
            real oldbx(nx,ny,nz,ngrd),oldby(nx,ny,nz,ngrd), &
                 oldbz(nx,ny,nz,ngrd),oldqpx(nx,ny,nz,ngrd), &
                 oldqpy(nx,ny,nz,ngrd),oldqpz(nx,ny,nz,ngrd), &
                 oldqrho(nx,ny,nz,ngrd),oldqpres(nx,ny,nz,ngrd), &
                 oldopx(nx,ny,nz,ngrd),oldopy(nx,ny,nz,ngrd), &
                 oldopz(nx,ny,nz,ngrd),oldorho(nx,ny,nz,ngrd), &
                 oldopres(nx,ny,nz,ngrd),oldhpx(nx,ny,nz,ngrd), &
                 oldhpy(nx,ny,nz,ngrd),oldhpz(nx,ny,nz,ngrd), &
                 oldhrho(nx,ny,nz,ngrd),oldhpres(nx,ny,nz,ngrd), &
                 oldepres(nx,ny,nz,ngrd)
      
            real wrkbx(nx,ny,nz,ngrd),wrkby(nx,ny,nz,ngrd), &
                 wrkbz(nx,ny,nz,ngrd),wrkqpx(nx,ny,nz,ngrd), &
                 wrkqpy(nx,ny,nz,ngrd),wrkqpz(nx,ny,nz,ngrd), &
                 wrkqrho(nx,ny,nz,ngrd),wrkqpres(nx,ny,nz,ngrd), &
                 wrkopx(nx,ny,nz,ngrd),wrkopy(nx,ny,nz,ngrd), &
                 wrkopz(nx,ny,nz,ngrd),wrkorho(nx,ny,nz,ngrd), &
                 wrkopres(nx,ny,nz,ngrd),wrkhpx(nx,ny,nz,ngrd), &
                 wrkhpy(nx,ny,nz,ngrd),wrkhpz(nx,ny,nz,ngrd), &
                 wrkhrho(nx,ny,nz,ngrd),wrkhpres(nx,ny,nz,ngrd), &
                 wrkepres(nx,ny,nz,ngrd)
      !
      !     unperturbed quantities
      !
            real bx0(nx,ny,nz,ngrd),by0(nx,ny,nz,ngrd),bz0(nx,ny,nz,ngrd)
      !    +     ,hrho0(nx,ny,nz,ngrd),orho0(nx,ny,nz,ngrd),
      !    +     qrho0(nx,ny,nz,ngrd),hpres0(nx,ny,nz,ngrd),
      !    +     opres0(nx,ny,nz,ngrd),qpres0(nx,ny,nz,ngrd),
      !    +     epres0(nx,ny,nz,ngrd)                   
      
            real efldx(nx,ny,nz),efldy(nx,ny,nz),efldz(nx,ny,nz), &
                 curx(nx,ny,nz),cury(nx,ny,nz),curz(nx,ny,nz), &
                 bsx(nx,ny,nz),bsy(nx,ny,nz),bsz(nx,ny,nz),btot(nx,ny,nz), &
                 resistive(nx,ny,nz)
      !
      !      variable time step arrays
      !
            real t_old(ngrd),t_new(ngrd),t_step(ngrd),t_step_new(ngrd)
      !
      !     boundary condition arrays
      !
            real :: bxf(ny,nz),byf(ny,nz),bzf(ny,nz), &
                     rhof(ny,nz),svxf(ny,nz),svyf(ny,nz),svzf(ny,nz)
            real :: bxp(ny,nz),byp(ny,nz),bzp(ny,nz), &
                     rhop(ny,nz),svxp(ny,nz),svyp(ny,nz),svzp(ny,nz)
            real :: future(ny,nz),past(ny,nz), &
                    bfld(ncts,4),rplas(ncts),svel(ncts,3)
            integer ncount(ny,nz)
      
            real tx(mx,my,mz),ty(mx,my,mz),tz(mx,my,mz),tg1(mx,my,mz), &
                 tg2(mx,my,mz2),tt(mx,my,mz),work(muvwp2,muvwp2), &
                 cross(my,mz),along(mx,mz)
      
            !1character*5 wd1,wd2,wd3,wd4
            !character*8 label
            !character*15 title
      
            integer ijsrf(3,20000),ijmid(3,20000),ijzero(3,50000)
            real    parm_srf(7,20000),parm_mid(7,20000), &
                    parm_zero(7,50000)
      
            logical start,add_dip,ringo,update,save_dat,write_dat, &
                    spacecraft,tilting,warp,divb_on,yes_step
      
      !      update decides if time setting is to be reset
      !
      !     ijsrf give position of ionosphere in grid units - plasma
      !           parameters stored in parm_srf
      !     ijmid gives intermediate boundary - say representing the
      !            atmosphere - plasma parameters stored in parm_mid
      !     ijzero gives position of all grid units interior to surface
      !
      !     d_min is the minimum allowable density
      !     stepsz is the size of the time step in terms of delta_x,y/(|umax|+c_s)
      
      !     A. Rajendar, 01/12/2013: inclusion of neutral cloud. Neutral density
      !     is specified at each grid point in particles per cm^3, as non-
      !     dimensionalized by rho_equiv in the input file.
            REAL, DIMENSION(nx,ny,nz) :: oh_cloud
            REAL, DIMENSION(nx,ny,nz) :: mass_load_photo, mass_load_e
            REAL, DIMENSION(nx,ny,nz) :: photo_x, photo_y, photo_z, &
                 elastic_x, elastic_y, elastic_z, &
                 e_impact_x, e_impact_y, e_impact_z
            REAL, DIMENSION(nx,ny,nz) :: eifrq, hoteifrq
      !     photoionization lifetime of relevant neutrals, days, O, OH, H2O
      !     (see Fleshman et al., 2010, Enceladus torus model, JGR)
            REAL, PARAMETER :: life_o = 5.0e3, life_oh = 3.1e3, life_w = 2.6e3
      !     rate coefficient for electron impact ionization (Burger et al. 2007,
      !     Enceladus' water plume), cm^3/s
            REAL, PARAMETER :: kappa_th = 7.3e-13, kappa_hot = 1.2e-8
      !     non-dimensionalized values for lifetime and rate coefficient
            REAL :: life_o1, life_oh1, life_w1, kappa_th1, kappa_hot1
      
      !     A. Rajendar, 11/04/2013: SI characteristic values of pressure and
      !     number density for ideal gas law definition of kB*T for Maxwell-
      !     Boltzmann electron distribution
      !     v_equiv and rho_equiv replaced by actual numerical values. Will be
      !     recast in terms of variables once F90 upgrade complete and modules
      !     are used
            REAL, PARAMETER :: pchar = 2*1.0e6*1.66e-27* &
                 ((1000*1000.0)**2)
            REAL, PARAMETER :: nchar = 2.*1.0e6
            REAL, PARAMETER :: tchar = (0.25*6.0e7)/(1000.*1000.) ! characteristic time, s
            REAL, PARAMETER :: lchar = 0.25*6.0e7 ! characteristic length, m
            REAL, PARAMETER :: vchar = 1000.*1000. ! characteristic speed (v_equiv in m/s)
      !     density of two sub-populations of the thermal electron population
      !     which are responsible for electron impact ionization. elecden1 is
      !     the population of electrons with energies between 13 and 40 eV,
      !     with a mean ionization cross-section of 9e-17 cm^2, while elecden2
      !     is the population with energies > 40 eV, which mean cs of 1.8e-16
      !     cm^2
      !      REAL, DIMENSION(nx,ny,nz) :: elecden1, elecden2
      !     electron number density values for box 1:
            REAL, DIMENSION(nx,ny,nz) :: elecden
      
            DOUBLE PRECISION :: ttime
      
            ! A.R. 2015-01-02: IMPLICIT NONE variable declaration
            REAL :: den_earth, alpha_e, alf_inner1, alf_inner2, cs_inner
            REAL :: rmassq, rmassh, rmasso
            !REAL :: rearth, xdip, ydip, zdip, 
            REAL :: tilt1, tilt2
            REAL :: tmax, stepsz, tsave
            INTEGER :: ntgraf
            REAL :: gamma, gamma1
            REAL :: o_conc, gravity, ti_te, re_wind
            REAL :: vx_wind1, vx_wind2, vy_wind1, vy_wind2, &
                 vz_wind1, vz_wind2, &
                 alfx_wind1, alfx_wind2, &
                 alfy_wind1, alfy_wind2, &
                 alfz_wind1, alfz_wind2, &
                 cs_wind
            REAL :: den_wind1, den_wind2
            REAL :: reynolds, resist
            REAL :: rho_frac, bfrac, vfrac
            REAL :: v_equiv, b_equiv, rho_equiv, re_equiv
            REAL :: uday, utstart
            REAL :: chirho, chipxyz, chierg
            REAL :: difrho, difpxyz, diferg
      !      REAL, DIMENSION(9) :: grd_xmin, grd_xmax, grd_ymin, grd_ymax, &
      !           grd_zmin, grd_zmax, xspac
      !      INTEGER, DIMENSION(9) :: ncore, nbndry
            INTEGER :: i, j, k, m, n
            INTEGER :: ix, iy, iz
            REAL :: dx, dy, dz
            REAL :: ax, ay, az
            REAL :: xp, yp, zp
            REAL :: erho, b01, b02, alf_lim
            REAL :: bmax, delb0!, b0
            REAL :: rho_wind1, rho_wind2
            REAL :: tilt, dtilt!, sin_tilt, cos_tilt
            REAL :: planet_rad, planet_per, lunar_rad
            REAL :: v_rot, r_rot
            REAL :: t, t_equiv, ut
            REAL :: grav, d_min
            REAL :: old_o_conc, o_conc_min, o_conc_max
            REAL :: cur_min, cur_max
            REAL :: vut, rut, zut
            REAL :: epress, eerg
            REAL :: srho, delrho, spress, serg
            REAL :: svelx, svely, svelz, spx, spy, spz
            REAL :: delvx_wind, delvy_wind, delvz_wind
            REAL :: delbx_wind, delby_wind, delbz_wind
            REAL :: sbx_wind, sby_wind, sbz_wind
            REAL :: delt, deltg, tgraf, ts1, tdiv, del_tdiv, t2
            INTEGER :: nchf, inchf
            INTEGER :: numzero, nummid, numsrf
      !      REAL :: rx, ry, rz
            REAL :: abx, aby, abz, abmag
            REAL :: radstrt
            REAL :: tot_cur_nth, peak_cur_nth
            REAL :: tot_cur_sth, peak_cur_sth
            REAL :: sheet_den, alpha_s, a_ion_max, a_ion_max1
            REAL :: ar, ra, rd, rvx, rvy, rvz
            REAL :: zheight, rho_iono, r_equat, rerg_sphere, rden_sphere
            REAL :: arho, aerg, ar_iono, ar_iono0, ra0, zheight0, rho_iono0
            REAL :: qrho0, hrho0, orho0
            REAL :: ar2, wind_bnd, ofrac
            REAL :: avx, avy, avz
            REAL :: apres, dfrac, distance
            INTEGER :: mout, nn
            INTEGER :: nc, nvx, nhi
            REAL :: zvelx, zvely, zvelz, displace, b_perp, bmag
            REAL, DIMENSION(nx,ny,nz) :: tvx, tvy, tvz, evx, evy, evz, &
                 vvx, vvy, vvz
            REAL :: pxmax, pymax, pzmax, pmax, csmax, alfmax, &
                 vlim, fastest, delt_old
            ! adaptive time stepping
            REAL :: smallest_step
            INTEGER :: mallest_step, nsteps, m_step
            REAL :: told, utold, old_tilt, delay
            REAL :: astep, dut
            INTEGER :: ms
            REAL :: delt2, atilt, t_grd
            REAL :: chifcs
            REAL :: xmin, xmax, ymin, ymax, zmin, zmax
            REAL :: xcut, av, gs, preslim, p0, po, blim, &
                 tempx, tempo, rho_lim
            REAL :: rfrac, rho_tot
            !REAL :: qden, hden, oden, pden, tden
            REAL :: aden, epx, epy, epz
            REAL :: aqrho, ahrho, aorho
      
            ! 2015-05-03 fluid particle tracker variables
            !LOGICAL :: fluid_track = .TRUE.
            !INTEGER, PARAMETER :: trackbox = 3 ! box in which tracking is done
            !INTEGER, PARAMETER :: numpart = 4 ! # of fluid particles
            !REAL, DIMENSION(3,numpart) :: fluidxyz, fluidxyz1
            !REAL, DIMENSION(nx,ny,nz) :: velx, vely, velz
            !REAL, DIMENSION(3) :: fluidvel
            !CHARACTER :: partstr*8
            !PARAMETER (partstr = 'partloc_')
      
            ! 2015-02-15: fixing the stupid way in which fluid files were numbered and written out.
            CHARACTER :: fluidstr*6, str1*5, str2*5
            PARAMETER (fluidstr = 'fluid_')
      
            INTEGER, PARAMETER :: mwrite = 3 ! box to write out radial flux
            REAL :: electron_rate
            REAL, DIMENSION(12) :: ndot
            REAL, DIMENSION(nx,ny,nz) :: hnvx, hnvy
            REAL, DIMENSION(nx,ny,nz) :: hden, hvx, hvy
            ! timing to write out global ionization, radial flux, etc.
            REAL, PARAMETER :: twrite = 1.26
            REAL :: twrto
      
            REAL, DIMENSION(3) :: mdot
            REAL, DIMENSION(nx,ny,nz,ngrd) :: hvelx, hvely
      
            ! 2015-09-01 lookup table test
            INTEGER, PARAMETER :: flen = 151 ! length of lookup table
            REAL, DIMENSION(flen) :: eev, kappa ! electron energy in eV, nondimensional kappa

			
			
			!!!!!!!!!!! Angela update 
			REAL :: spinvectorx, spinvectory, spinvectorz, utdec_start
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
            namelist/option/tmax,ntgraf,stepsz,start,tsave
            namelist/earth/xdip,ydip,zdip,rearth, &
                            tilt1,tilt2,tilting,rmassq,rmassh,rmasso
            namelist/speeds/cs_inner,alf_inner1,alf_inner2, &
                            alpha_e,den_earth,o_conc,gravity, &
                            ti_te,gamma,ringo,update,divb_on
            namelist/windy/re_wind,cs_wind,vx_wind1,vx_wind2, &
                          vy_wind1,vy_wind2,vz_wind1,vz_wind2, &
                          alfx_wind1,alfx_wind2, &
                          alfy_wind1,alfy_wind2, &
                          alfz_wind1,alfz_wind2, &
                          den_wind1,den_wind2, &
                         reynolds,resist,rho_frac,bfrac,vfrac
            namelist/physical/re_equiv,b_equiv,v_equiv,rho_equiv, &
                          spacecraft,warp,utstart
            namelist/smooth/chirho,chipxyz,chierg, &
                            difrho,difpxyz,diferg
      
            ! open input data file
            open(5,file='fmpd3din',status='old',form='formatted')
      
      !     read input parameters
            read(5,option)
            read(5,earth)
            read(5,speeds)
            read(5,windy)
            read(5,physical)
            read(5,smooth)
      !
      !     output to test whether the parameters are the actual ones you want
      !
            write(6,option)
            write(6,earth)
            write(6,speeds)
            write(6,windy)
            write(6,physical)
            write(6,smooth)
      
            do i=1,ngrd
             read(5,*)grd_xmin(i),grd_xmax(i),grd_ymin(i),grd_ymax(i), &
                  grd_zmin(i),grd_zmax(i),xspac(i),ncore(i),nbndry(i)
             write(6,*)grd_xmin(i),grd_xmax(i),grd_ymin(i),grd_ymax(i), &
                  grd_zmin(i),grd_zmax(i),xspac(i),ncore(i),nbndry(i)
             ix=1+(grd_xmax(i)-grd_xmin(i))/xspac(i)
             iy=1+(grd_ymax(i)-grd_ymin(i))/xspac(i)
             iz=1+(grd_zmax(i)-grd_zmin(i))/xspac(i)
             if((ix.ne.nx).or.(iy.ne.ny).or.(iz.ne.nz))then
               write(6,*)' WARNING: SIZES',ix,iy,iz,nx,ny,nz
               stop
             endif
            enddo
      
            !IF (fluid_track) THEN
            !   WRITE(*,*) 'READING IN FLUID PARTICLE LOCATIONS'
            !   OPEN(104,file='partlocin',status='unknown',form='formatted')
            !   DO i = 1, numpart
            !      READ(104,*) fluidxyz(1,i), fluidxyz(2,i), fluidxyz(3,i)
            !      WRITE(*,*) fluidxyz(1,i), fluidxyz(2,i), fluidxyz(3,i)
            !   END DO
            !END IF
      
      !     ncore -> grid gets all the information from this grid no.
      !     nbdry -> which grid data will be applied to this grid no.for bndry
      
      !     A. Rajendar, 01/23/2012: calculation of non-dimensionalized
      !     photoionization lifetimes and rate coefficients
             life_o1 = (life_o*86400.)*(v_equiv*1000.)/(re_equiv*6.0e7)
             life_oh1 = (life_oh*86400.)*(v_equiv*1000.)/(re_equiv*6.0e7)
             life_w1 = (life_w*86400.)*(v_equiv*1000.)/(re_equiv*6.0e7)
      !     A.R. 03/01/2013: multiplication of rho_equiv*1e6 removed, units were
      !     inconsistent
             kappa_th1 = kappa_th* &
                  rho_equiv*(re_equiv*6.0e7)/(v_equiv*1000.)
             kappa_hot1 = kappa_hot* &
                  rho_equiv*(re_equiv*6.0e7)/(v_equiv*1000.)
             write(6,*)life_o1,life_oh1,life_w1,kappa_th1,kappa_hot1
      
      !     calculate effective magnetic field strength
            erho=den_earth*rmassh	! Multiplying by rmassh here is wrong wrong wrong. 7/17/2019 MJS
            b01=alf_inner1*sqrt(erho)*rearth**3
            b02=alf_inner2*sqrt(erho)*rearth**3
            erho=den_earth	! So that subsequent uses of erho will be number density. 7/17/2019 MJS
            alf_lim=3.00*alf_inner1
            b0=b01
			PRINT *, b0
            bmax=b02
            delb0=(b02-b01)/tmax
      
      !     re_wind sets raduius from earth where initial wind placed
      !     reynolds coefficient for surface currents
      !     resist equivalent if you wish to run anomalous resistivity
      !     bfrac determines the percentage of the tangential magnetic
      !     field allowed at earth's surface
      
            tilt=tilt1
            sin_tilt=sin(tilt*.0174533)
            cos_tilt=cos(tilt*.0174533)
            dtilt=(tilt2-tilt1)/tmax
      !
      !     jupiter parameters
      !
      !     planet_rad=71000.   !km
      !     planet_per=9.7     !hr
      !     lunar_rad=5.9      !orbital radii
      !     v_rot=6.2832*planet_rad/(planet_per*3600.)/v_equiv  ! normalized units
      !     r_rot=40.0   !Re where corotation stops
      !
      !     Earth parameters
           planet_rad=6371.   !km
           planet_per=24.     !hr
           lunar_rad=60.      !orbital radii
           v_rot=6.2832*planet_rad/(planet_per*3600.)/v_equiv  ! normalized units
           r_rot=10.0   !Re where corotation stops
      !
      !     Saturn parameters
      !
      !      planet_rad=1560.   !km
      !      planet_per=1000000000.     !hr
      !      lunar_rad=60.0      !orbital radii
      !      v_rot=6.2832*planet_rad/(planet_per*3600.)/v_equiv  ! normalized units
      !      r_rot=12.5   !Re where corotation stops
      !
      !     spacecraft stuff
      !
      !     gravity in m/s**2 at earths surface !need t_equiv in normalization
            t=0.
            t_equiv=planet_rad*re_equiv/v_equiv
            grav=gravity*(planet_rad*re_equiv*1000.)/(1000.*v_equiv)**2
            ut=utstart
            d_min=0.00001
      !
      !      ionospheric parameters
      !
            old_o_conc=o_conc
            o_conc_min=0.01
            o_conc_max=10.0
            cur_min=0.75
            cur_max=20.0
      !
      !
      !     A.R. 03/02/2013: spacecraft no. increased to 12, positions
      !     modified to include one 35 RS upstream, the rest arranged in a
      !     grid pattern 60 RS downtail
      !
      !     A.R. 03/07/2013: changed again to capture plasmoid on dawnside
      !     launched about 20 RS away from Saturn
      !
      !     Initial position of spacecraft in RE but simulation directions
      !     UPSTREAM :
            xcraft(1,1) = -35.
            xcraft(2,1) = 0.
            xcraft(3,1) = 0.
      !     reference craft
            rcraft(1) = xcraft(1,1)
            rcraft(2) = xcraft(2,1)
            rcraft(3) = xcraft(3,1)
      
      !     CLOSE
            xcraft(1,2) = 24.
            xcraft(2,2) = 0.
            xcraft(3,2) = 2.
      !     FAR LOW
            xcraft(1,3) = 40.
            xcraft(2,3) = 0.
            xcraft(3,3) = 2.
      !     FAR HIGH
            xcraft(1,4) = 40.
            xcraft(2,4) = 0.
            xcraft(3,4) = 5.
      !     FAR DUSK LOW
            xcraft(1,5) = 40.
            xcraft(2,5) = -8.
            xcraft(3,5) = 2.                                                                                                                             
      !     FAR DUSK HIGH
            xcraft(1,6) = 40.
            xcraft(2,6) = -8.
            xcraft(3,6) = 5.
      
            zut=ut-.01
            rut=zut
            vut=rut
            xcraft(4,1) = zut
            xcraft(4,2) = zut
            xcraft(4,3) = zut
            xcraft(4,4) = zut
            xcraft(4,5) = zut
            xcraft(4,6) = zut
            call limcraft(xcraft,ncraft,re_equiv,ngrd)
      
            do n=1,ncraft
               do i=1,4
                  zcraft(i,n)=xcraft(i,n)
               enddo
            enddo
      
            OPEN(501,file='looptiming.dat',status='unknown',form='formatted') ! loop timing
            OPEN(502,file='globalion.dat',status='unknown',form='formatted') ! ionization rate
            OPEN(503,file='radflux.dat',status='unknown',form='formatted') ! radial number flux
      
            OPEN(31,file='upstream.dat',status='unknown',form='formatted') 
            OPEN(32,file='close.dat',status='unknown',form='formatted')
            OPEN(33,file='far_lo.dat',status='unknown',form='formatted')
            OPEN(34,file='far_hi.dat',status='unknown',form='formatted')
            OPEN(35,file='far_dsk_lo.dat',status='unknown',form='formatted')
            OPEN(36,file='far_dsk_hi.dat',status='unknown',form='formatted')
            if(spacecraft) then
              open(21,file='wind.pos',status='unknown',form='formatted')
      !       open(22,file='polar.pos',status='unknown',form='formatted')
      !       open(23,file='equators.pos',status=unknown'',form='formatted')
      !       open(24,file='geotail.pos',status='unknown',form='formatted')
      !       open(27,file='wind.den',status='unknown',form='formatted')
      !       open(28,file='wind.vel',status='unknown',form='formatted')
              open(27,file='wind.plas',status='unknown',form='formatted')
              open(29,file='wind.mag',status='unknown',form='formatted')
            endif
      
            ! 2015-09-01 open lookup table
            OPEN(701,file='neutral_cloud/kappatable.dat',status='old',form='formatted')
            DO i = 1, flen
               READ(701,*) eev(i), kappa(i)
            END DO
            CLOSE(701)
      
      !     the magnetic field in dimensionless units is actually in Alfven speeds
      !             relative to the normalizing velocity
      !     the temperature is in sound speeds relative to the normalizing speed
      !             all squared
      !    
      !     the magnetospheric plasma density are assumed to vary as
      !             rho proportional to (R)**-alpha_e
      !             Temperatue proportional to (R)**alpha_e
      !         with total pressure constant which is needed for equilibrium
      !
      !     now find the equivalent pressure of magnetosphere for the given
      !         sound speed
      !
            gamma1=gamma-1.
            epress=cs_inner**2*erho/gamma
            eerg=epress/gamma1
      !
      !     do the same for the solar wind
      !
            rho_wind1=den_wind1*rmassq
            rho_wind2=den_wind2*rmassq
            srho=rho_wind1
            delrho=(rho_wind2-rho_wind1)/tmax
            spress=(cs_wind**2*srho/gamma)/gamma1
            svelx=vx_wind1
            svely=vy_wind1
            svelz=vz_wind1
      !
            spx=srho*svelx
            spy=srho*svely
            spz=srho*svelz
            serg=0.5*(svelx**2+svely**2+svelz**2)*srho+spress
            delvx_wind=(vx_wind2-vx_wind1)/tmax
            delvy_wind=(vy_wind2-vy_wind1)/tmax
            delvz_wind=(vz_wind2-vz_wind1)/tmax
      !
            delbx_wind=(alfx_wind2*sqrt(rho_wind2) &
                       -alfx_wind1*sqrt(rho_wind1))/tmax
            delby_wind=(alfy_wind2*sqrt(rho_wind2) &
                       -alfy_wind1*sqrt(rho_wind1))/tmax
            delbz_wind=(alfz_wind2*sqrt(rho_wind2) &
                       -alfz_wind1*sqrt(rho_wind1))/tmax
            sbx_wind=alfx_wind1*sqrt(rho_wind1)
            sby_wind=alfy_wind1*sqrt(rho_wind1)
            sbz_wind=alfz_wind1*sqrt(rho_wind1)
      !
            deltg=tmax/float(ntgraf)
            delt=stepsz
            tgraf=0.
            ts1=tsave
            tdiv=10.
            del_tdiv=1.
            write_dat=.true.
            twrto = twrite
      
            CALL old_neutral_cloud(oh_cloud,nx,ny,nz)
      !
      !     ************************************************
      !     check for restart
      !     ************************************************
      !
            nchf=11
            if(start) go to 80
            open(nchf,file='fluid_in',status='old',form='unformatted')
      
            ! read restart data
            read(nchf)t
            read(nchf)qrho
            read(nchf)qpx
            read(nchf)qpy
            read(nchf)qpz
            read(nchf)qpres
            read(nchf)hrho
            read(nchf)hpx
            read(nchf)hpy
            read(nchf)hpz
            read(nchf)hpres
            read(nchf)orho
            read(nchf)opx
            read(nchf)opy
            read(nchf)opz
            read(nchf)opres
            read(nchf)bx
            read(nchf)by
            read(nchf)bz
            read(nchf)epres
            read(nchf)bx0
            read(nchf)by0
            read(nchf)bz0
            read(nchf)mass_load_photo
            read(nchf)mass_load_e
            read(nchf)photo_x
            read(nchf)photo_y
            read(nchf)photo_z
            read(nchf)e_impact_x
            read(nchf)e_impact_y
            read(nchf)e_impact_z
            read(nchf)elastic_x
            read(nchf)elastic_y
            read(nchf)elastic_z
            read(nchf)parm_srf,parm_mid,parm_zero, &
                 ijzero,numzero,ijmid,nummid,ijsrf,numsrf
            close(nchf)
      
            if(update)t=0.
            ts1=t+tsave
            tmax=t+tmax
            tgraf=t+deltg
            tdiv=t
            nchf=11
            ut=utstart+t*t_equiv/3600.
      
            twrto = t + twrite
      
      !
      !     rescale oxygen density at inner boundary by oxygen scale
      !
            if(ringo)then
             m=1
             rx=xspac(m)
             ry=xspac(m)
             rz=xspac(m)
      !
             call calcur(bx,by,bz,ngrd,nx,ny,nz,m,curx,cury,curz, &
                     ijsrf,numsrf,ijmid,nummid,ijzero,numzero)
      !
      !        field aligned currents
      !
             do k=1,nz
             do j=1,ny
             do i=1,nx
               abx=bx(i,j,k,m)+bx0(i,j,k,m)
               aby=by(i,j,k,m)+by0(i,j,k,m)
               abz=bz(i,j,k,m)+bz0(i,j,k,m)
               abmag=sqrt(abx**2+aby**2+abz**2) &
                             +0.00000001
               btot(i,j,k)=(curx(i,j,k)*abx+cury(i,j,k)*aby &
                                     +curz(i,j,k)*abz)/abmag
             enddo
             enddo
             enddo
      !
             save_dat=.false.
             add_dip=.false.
      !
             radstrt=rearth+8.
      !       write(6,*)'entering aurora_cur', radstrt
      !       call aurora_cur(btot,nx,ny,nz,m,radstrt,re_equiv,1, &
      !              ut,save_dat,add_dip,'cur4 nth',14,write_dat, &
      !            b_equiv,planet_rad,tot_cur_nth,peak_cur_nth)
             old_o_conc=o_conc
             if(tot_cur_nth.gt.cur_max)then
               o_conc=o_conc_max
              else if(tot_cur_nth.gt.cur_min)then
               o_conc=o_conc_min+(o_conc_max-o_conc_min)* &
                         (tot_cur_nth-cur_min)/(cur_max-cur_min)
              else
               o_conc=o_conc_min
              endif
      !        write(10,*)ut,tot_cur_nth,o_conc,old_o_conc
      !
              do n=1,numsrf
               parm_srf(3,n)=o_conc*(rmasso/rmassh)*parm_srf(2,n)
               parm_srf(6,n)=o_conc*parm_srf(5,n)
               parm_srf(7,n)=(parm_srf(5,n)+ parm_srf(6,n))/ti_te
             enddo
      !
              do n=1,nummid
               parm_mid(3,n)=o_conc*(rmasso/rmassh)*parm_mid(2,n)
               parm_mid(6,n)=o_conc*parm_mid(5,n)
               parm_mid(7,n)=(parm_mid(5,n)+ parm_mid(6,n))/ti_te
             enddo
      
            endif
      !
      !     initialize plasma resistivity
      !
            m=1
            rx=xspac(m)
            ry=xspac(m)
            rz=xspac(m)
            call set_resist(resistive,nx,ny,nz,resist,ijzero,numzero, &
                         ijmid,nummid,ijsrf,numsrf,re_equiv,v_equiv)
      
            write(6,79) nchf
       79   format('  restart from   mpd3d',i2)
            rewind nchf
            nchf=11
            goto 170
      
      !
      !     ******************************
      !            initialization
      !     ******************************
      !
      !     ambient plasma
      !
      !      initialize indices and surface points
      !
       80   numsrf=0
            nummid=0
            numzero=0
      !
      !      scale lengths for plasma sheet population
      !
            sheet_den=0.25
            alpha_s=4.
      
      !
      !     A. Rajendar 01/25/2013: moved neutral_cloud out of next loop
      ! 
      
      !     A.R. 05/21/2013
            a_ion_max = 1.0
      
            do 130 m=1,ngrd
               dx=(grd_xmax(m)-grd_xmin(m))/(nx-1.)
               dy=(grd_ymax(m)-grd_ymin(m))/(ny-1.)
               dz=(grd_zmax(m)-grd_zmin(m))/(nz-1.)
      
      !     A.R. 05/21/2013
               a_ion_max1 = a_ion_max
      
      !
      !      create dipole magnetic field and load magnetospheric plasma 
      !
               do 130 k=1,nz
                  az=(grd_zmin(m)+dz*(k-1)-zdip)
      !
                  do 120 j=1,ny
                     ay=grd_ymin(m)+dy*(j-1)-ydip
      !
                     do 110 i=1,nx
                        ax=(grd_xmin(m)+dx*(i-1)-xdip)
                        xp=ax*cos_tilt-az*sin_tilt 
                        zp=ax*sin_tilt+az*cos_tilt
                        ar=sqrt(xp**2+ay**2+zp**2)
      
      !
      !        determine magnetic dipole field
      !
                        call dipole(bx0(i,j,k,m),by0(i,j,k,m), &
                             bz0(i,j,k,m),ax,ay,az)
      !
      !         zero interior magnetic field so alfven speed small
      !
                        if (ar.lt.rearth-1.0)then
                           bx0(i,j,k,m)=0.
                           by0(i,j,k,m)=0.
                           bz0(i,j,k,m)=0.
                        endif
      !
      !        set up rotational properties
      !
                        rx=ax*re_equiv
                        ry=ay*re_equiv
                        rd=sqrt(rx**2+ry**2)
                        if(rd.lt.r_rot)then
                           vfrac=1.
                        else
                           vfrac=((2.*r_rot**2)/(r_rot**2+ rd**2))**2
                        endif
                        rvy=rx*v_rot*vfrac
                        rvx=-ry*v_rot*vfrac
						

      

      !*******************************************************************************
	  ! XIN CHANGE ********************************************************************
	           call vfldrot(rvx,rvy,rvz,ax,ay,az,vfrac,v_rot,r_rot)
!***********************************************************************************

      !        for Jupiter top hat distribution
      !        ar_iono=sqrt(xp**2+ay**2+0.75*zp**2)
      !
      !        for rearth = 6
      !        ar_iono=sqrt(xp**2+ay**2+2.*zp**2)
      !
      !        for rearth = 10
                        ar_iono=sqrt(xp**2+ay**2+1.25*zp**2)
      !        isotropic
      !        ar_iono=sqrt(xp**2+ay**2+zp**2)
      !
                        ar_iono=amax1(0.0001,ar_iono)
                        ra=((ar_iono+0.5*rearth)/(1.5*rearth))**(-alpha_e)
                        zheight=amax1(1.,(zp**2+(1.5*rearth)**2)/ &
                             (3.0*rearth)**2)
                        ra=ra/zheight**2
                        rho_iono=amax1(erho*ra,0.001)
      !  
                        r_equat=(ar**3+0.001)/(xp**2+ay**2+0.001)
                        r_equat=amax1(r_equat,rearth)
                        rerg_sphere=(0.001+(rearth/r_equat)**4)
      !
                        if(r_equat.gt.5.*rearth)then
                           rerg_sphere=(0.001+(rearth/r_equat)**alpha_e) !constant pressure flux
                           rden_sphere=1.
                        else
                           rerg_sphere=0.1+0.9*(r_equat-rearth)/(4.*rearth)
                           rden_sphere=(rerg_sphere**2) !reduce O+ in plasmasphere
                        endif
      
                        arho=amax1(rho_iono,d_min)
                        aerg=amin1(abs(eerg*rerg_sphere),(0.1**2)*arho)
             
      !     ADDED 05/21/2013
      !     A. Rajendar, 02/07/2013: modification of proton density, increases 10-100X inside of
      !     radial distance 4 RS to the ionospheric boundary. Trying to reduce Alfven speed issues.
      !                  IF (ar.LT.a_ion_max1) THEN
      !     qrho at radial distance 4RS
      !                     ar_ion0 = amax1(0.0001, ar)
      !                     ra0 = ((ar_iono0 + 0.5*rearth)/
      !     +                    (1.5*rearth))**(-alpha_e)
      !                     zheight0 = amax1(1., (zp**2 + (1.5*rearth)**2)/
      !     +                    (3.0*rearth)**2)
      !                     ra0 = ra0/zheight0**2
      !                     rho_iono0 = amax1(erho*ra0, 0.001)
      !                     qrho0 = rmassq*amax1(rho_iono0, d_min)
      !                     
      !                     qrho(i,j,k,m) = qrho0*
      !     +                    exp(0.5*ABS(ar - a_ion_max1))
      !                  ELSE
                        qrho(i,j,k,m) = arho*rmassq
      !                  END IF
                        
                        qpres(i,j,k,m)=0.5*gamma1*aerg
                        qpx(i,j,k,m)=hrho(i,j,k,m)*rvx	! I really don't understand why this doesn't break with implicit none. hrho has not been defined yet--it must be from compiler flags that set it to zero upon initialization. 7/17/2019 MJS
                        qpy(i,j,k,m)=hrho(i,j,k,m)*rvy
                        qpz(i,j,k,m)=0.
      
                        orho(i,j,k,m)=rho_iono*rmasso*o_conc!*rden_sphere	! rden_sphere seems to be only for large planets and does not appear to function correctly for bodies like Europa.
                        opres(i,j,k,m)=0.5*gamma1*aerg*o_conc!*rden_sphere
                        opx(i,j,k,m)=orho(i,j,k,m)*rvx
                        opy(i,j,k,m)=orho(i,j,k,m)*rvy
                        opz(i,j,k,m)=0.
      
      !     A.Rajendar, 02-04-2013: add 10X water ions if sqrt(x^2 + y^2) > a_ion_max1,                                                                     
      !     with exponential decrease of density to d_min at/near ionosphere boundary                                                                       
                        IF (rd.LT.a_ion_max1) THEN
      !     Determine density at r = a_ion_max1, then apply exponential decrease                                                                            
                           ar_iono0 = sqrt(a_ion_max1**2. + 1.25*zp**2.)
                           ar_iono0 = amax1(0.0001, ar_iono0)
                           ra0 = ((ar_iono0 + 0.5*rearth)/ &
                                (1.5*rearth))**(-alpha_e)
                           zheight0 = amax1(1.,(zp**2. + (1.5*rearth)**2.)/ &
                                (3.0*rearth)**2.)
                           ra0 = ra0/zheight0**2.
                           rho_iono0 = amax1(erho*ra0, 0.001)
                           qrho0 = rmassq*amax1(rho_iono0, d_min)
                           hrho0 = 1.0*qrho0*rmassh/rmassq
                           
                           hrho(i,j,k,m) = amax1(0.16, &
                                hrho0*exp(-1.25*ABS(rd - a_ion_max1)))
                        ELSE
                           hrho(i,j,k,m)=1.0*qrho(i,j,k,m) &
                                *rmassh/rmassq
                        END IF
      
                        IF (rd.LE.3.0) THEN
                           hrho(i,j,k,m) = 1.*qrho(i,j,k,m)
                        END IF
      
                        hpres(i,j,k,m)=1.0*qpres(i,j,k,m)
                        hpx(i,j,k,m)=hrho(i,j,k,m)*rvx
                        hpy(i,j,k,m)=hrho(i,j,k,m)*rvy
                        hpz(i,j,k,m)=0.
      
                        epres(i,j,k,m)=(qpres(i,j,k,m)+hpres(i,j,k,m)+ &
                             opres(i,j,k,m))/ti_te
      !
      !
                        bx(i,j,k,m)=0.
                        by(i,j,k,m)=0.
                        bz(i,j,k,m)=0.
      !
      !        test for boundary point of planets surface or
      !        interior to planet
      !
                        if((ar.gt.rearth+.6).or.(m.gt.1))goto 110
      !
                        if(ar.lt.rearth-1.5) then
                           numzero=numzero+1
                           ijzero(1,numzero)=i
                           ijzero(2,numzero)=j
                           ijzero(3,numzero)=k
      !
      
                           parm_zero(1,numzero)=qrho(i,j,k,m)
                           parm_zero(2,numzero)=hrho(i,j,k,m)
                           parm_zero(3,numzero)=orho(i,j,k,m)
                           parm_zero(4,numzero)=qpres(i,j,k,m)
                           parm_zero(5,numzero)=hpres(i,j,k,m)
                           parm_zero(6,numzero)=opres(i,j,k,m)
                           parm_zero(7,numzero)=epres(i,j,k,m)
      !
      !        
                        else  if(ar.lt.rearth-.5) then
                           nummid=nummid+1
                           ijmid(1,nummid)=i
                           ijmid(2,nummid)=j
                           ijmid(3,nummid)=k
                           ar2=sqrt(xp**2+ay**2)
                           parm_mid(1,nummid)=qrho(i,j,k,m)
                           parm_mid(2,nummid)=hrho(i,j,k,m)
                           parm_mid(3,nummid)=orho(i,j,k,m)
                           parm_mid(4,nummid)=qpres(i,j,k,m)
                           parm_mid(5,nummid)=hpres(i,j,k,m)
                           parm_mid(6,nummid)=opres(i,j,k,m)
                           parm_mid(7,nummid)=epres(i,j,k,m)
      !
                        else 
                           numsrf=numsrf+1
                           ijsrf(1,numsrf)=i
                           ijsrf(2,numsrf)=j
                           ijsrf(3,numsrf)=k
                           ar2=sqrt(xp**2+ay**2)
                           parm_srf(1,numsrf)=qrho(i,j,k,m)
                           parm_srf(2,numsrf)=hrho(i,j,k,m)
                           parm_srf(3,numsrf)=orho(i,j,k,m)
                           parm_srf(4,numsrf)=qpres(i,j,k,m)
                           parm_srf(5,numsrf)=hpres(i,j,k,m)
                           parm_srf(6,numsrf)=opres(i,j,k,m)
                           parm_srf(7,numsrf)=epres(i,j,k,m)
                        endif
      !
       110           continue
       120        continue
       130     continue
      !
               write(6,132)numsrf,nummid,numzero
       132     format(' interior and zero points',3(1x,i6))
      !
      !     initialize solar wind plasa can be placed beyond
      !       the earth at a radius of re_wind
      !
      
               WRITE(6,*)'A.R. 1'
            wind_bnd=r_rot/1.5/re_equiv
            ofrac=rho_frac
            do m=1,ngrd
      !
            dx=(grd_xmax(m)-grd_xmin(m))/(nx-1.)
            dy=(grd_ymax(m)-grd_ymin(m))/(ny-1.)
            dz=(grd_zmax(m)-grd_zmin(m))/(nz-1.)
      !
            do k=1,nz
             az=(grd_zmin(m)+dz*(k-1)-zdip)
      !
             do j=1,ny
             ay=grd_ymin(m)+dy*(j-1)-ydip
              do i=1,nx
               ax=(grd_xmin(m)+dx*(i-1)-xdip)
      
               ar=sqrt(ax**2+ay**2+az**2)
               if((ar.ge.1.5*wind_bnd).and.(ax.lt.0.))then
                qrho(i,j,k,m)=srho+qrho(i,j,k,m)
                qpx(i,j,k,m)=spx+qpx(i,j,k,m)
                qpy(i,j,k,m)=spy+qpy(i,j,k,m)
                qpz(i,j,k,m)=spz+qpz(i,j,k,m)
      !
                avx=qpx(i,j,k,m)/qrho(i,j,k,m)
                avy=qpy(i,j,k,m)/qrho(i,j,k,m)
                avz=qpz(i,j,k,m)/qrho(i,j,k,m)
      !
                apres=cs_wind**2*qrho(i,j,k,m)/gamma
                qpres(i,j,k,m)=0.5*apres+qpres(i,j,k,m)
                epres(i,j,k,m)=0.5*apres/ti_te+epres(i,j,k,m)
      !
                hrho(i,j,k,m)=srho*rho_frac*rmassh/rmassq+ &
                              hrho(i,j,k,m)
                hpx(i,j,k,m)=avx*hrho(i,j,k,m)
                hpy(i,j,k,m)=avy*hrho(i,j,k,m)
                hpz(i,j,k,m)=avz*hrho(i,j,k,m)
                hpres(i,j,k,m)=0.5*spress*rho_frac &
                           +hpres(i,j,k,m)
      !
                orho(i,j,k,m)=srho*ofrac*rmasso/rmassq+ &
                           orho(i,j,k,m)
                opx(i,j,k,m)=avx*orho(i,j,k,m)
                opy(i,j,k,m)=avy*orho(i,j,k,m)
                opz(i,j,k,m)=avz*orho(i,j,k,m)
                opres(i,j,k,m)=0.5*spress*ofrac &
                           +opres(i,j,k,m)
               endif
               if((ar.gt.wind_bnd).and.(ar.lt.1.5*wind_bnd) &
                       .and.(ax.lt.0.0))then
                dfrac=(ar-wind_bnd)/(0.5*wind_bnd)
                qrho(i,j,k,m)=srho*dfrac+qrho(i,j,k,m)
                qpx(i,j,k,m)=spx*dfrac+qpx(i,j,k,m)
                qpy(i,j,k,m)=spy*dfrac+qpy(i,j,k,m)
                qpz(i,j,k,m)=spz*dfrac+qpz(i,j,k,m)
                avx=qpx(i,j,k,m)/qrho(i,j,k,m)
                avy=qpy(i,j,k,m)/qrho(i,j,k,m)
                avz=qpz(i,j,k,m)/qrho(i,j,k,m)
      !
                apres=cs_wind**2*qrho(i,j,k,m)/gamma
                qpres(i,j,k,m)=0.5*apres*dfrac+qpres(i,j,k,m)
                epres(i,j,k,m)=0.5*apres/ti_te*dfrac+epres(i,j,k,m)
      !
                hrho(i,j,k,m)=srho*rho_frac*rmassh/rmassq*dfrac+ &
                              hrho(i,j,k,m)
                hpx(i,j,k,m)=avx*hrho(i,j,k,m)
                hpy(i,j,k,m)=avy*hrho(i,j,k,m)
                hpz(i,j,k,m)=avz*hrho(i,j,k,m)
                hpres(i,j,k,m)=0.5*spress*rho_frac*dfrac &
                           +hpres(i,j,k,m)
      !
                orho(i,j,k,m)=srho*ofrac*rmasso/rmassq*dfrac+ &
                           orho(i,j,k,m)
                opx(i,j,k,m)=avx*orho(i,j,k,m)
                opy(i,j,k,m)=avy*orho(i,j,k,m)
                opz(i,j,k,m)=avz*orho(i,j,k,m)
                opres(i,j,k,m)=0.5*spress*ofrac*dfrac &
                           +opres(i,j,k,m)
               endif
      !
               enddo
              enddo
             enddo
            enddo
      !
      !
      !
        170 ut=utstart+t*t_equiv/3600.
      !
      !     initialize plasma resistivity
      !
            m=1
            rx=xspac(m)
            ry=xspac(m)
            rz=xspac(m)
            call set_resist(resistive,nx,ny,nz,resist,ijzero,numzero, &
                         ijmid,nummid,ijsrf,numsrf,re_equiv,v_equiv)
            add_dip=.false.
            save_dat=.false.
            write_dat=.false.
      
            if(spacecraft)then
      !
      !      calculate the distance between the reference spacecraft and
      !          solar wind boundary
      !
              m=ngrd
              distance=grd_xmin(m)-rcraft(1)/re_equiv
      !
      !      read down data file until correct time in data file
      !
      !      do  n=1,ncraft
      !      do  n=1,2
                n=1
                mout=20+n
                do while(ut.ge.zcraft(4,n))
                 read(mout,*)zcraft(4,n),zcraft(1,n), &
                      zcraft(2,n),zcraft(3,n)
      !
      !          change direction to get from GSM to simulation coords
      !
                  zcraft(1,n)=-zcraft(1,n)
                  zcraft(2,n)=-zcraft(2,n)
                if(n.eq.1)then
                 rcraft(1)=zcraft(1,1)
                 rcraft(2)=zcraft(2,1)
                 rcraft(3)=zcraft(3,1)
                endif
      !
                end do
      !      end do 
      !
      !     zcraft(1,2)=zcraft(1,2)
      !     zcraft(2,2)=zcraft(2,2)-1.2
      !     zcraft(3,2)=zcraft(3,2)
      !
             call limcraft(zcraft,ncraft,re_equiv,ngrd)
      !
      !      calculate the distance between the reference spacecraft and
      !          solar wind boundary
      !
             distance=grd_xmin(ngrd)-rcraft(1)/re_equiv
             write(6,*)'wind displaced by',distance
             do n=1,ncraft
              do nn=1,4
                xcraft(nn,n)=zcraft(nn,n)
              end do
             end do
      !
      !     initialize counting array
      !
            do  j=1,ny
             do  k=1,nz
              ncount(j,k)=0
              future(j,k)=ut-0.01 
             enddo
            enddo
      
      !      read all the magnetic field data to minize data sorting
            do  m=1,ncts
             read(29,*)bfld(m,4),bmag,bfld(m,1),bfld(m,2),bfld(m,3)
             read(27,*)rut,rplas(m),svel(m,1),svel(m,2),svel(m,3)
            enddo
      
      !      set timing
             nvx=0
             vut=-999.
             do while(ut.gt.vut) 
              svelx=zvelx
              nvx=nvx+1
              zvelx=-svel(nvx,1)/v_equiv
              vut=bfld(nvx,4)+t_equiv*distance/zvelx/3600.
             end do
      !
            write(6,*)'ut=',ut,'wind time',bfld(nvx,4)
      !  
            displace=0.
            m=ngrd
            dx=(grd_xmax(m)-grd_xmin(m))/(nx-1.)
            dy=(grd_ymax(m)-grd_ymin(m))/(ny-1.)
            dz=(grd_zmax(m)-grd_zmin(m))/(nz-1.)    
      !
            do 175 j=1,ny
            do 175 k=1,nz
             do while((ut.gt.future(j,k)) &
               .and.(ncount(j,k)+1.le.ncts))
              nc=ncount(j,k)+1
              bxp(j,k)=bxf(j,k)
              byp(j,k)=byf(j,k)
              bzp(j,k)=bzf(j,k)
              rhop(j,k)=rhof(j,k)
              svxp(j,k)=svxf(j,k)
              svyp(j,k)=svyf(j,k)
              svzp(j,k)=svzf(j,k)
              past(j,k)=future(j,k)
      !
              future(j,k)=bfld(nc,4)
              bxf(j,k)=-bfld(nc,1)/b_equiv
              byf(j,k)=-bfld(nc,2)/b_equiv
              bzf(j,k)=bfld(nc,3)/b_equiv
              rhof(j,k)=rplas(nc)/rho_equiv
              svxf(j,k)=-svel(nc,1)/v_equiv
              svyf(j,k)=0.
              svzf(j,k)=0.
              ncount(j,k)=nc
              avx=svxf(j,k)
      !
      !      calculate delay
      !
             if(warp)then
              b_perp=sqrt(bzf(j,k)**2+byf(j,k)**2)
              b_perp=amax1(b_perp,0.33*abs(bxf(j,k)))
              ay=grd_ymin(m)+dy*(j-1)-rcraft(2)/re_equiv
              az=grd_zmin(m)+dz*(k-1)-rcraft(3)/re_equiv
      !
      !       going to assume Bz IMF on average pos and ignore transients
      !               and By IMF is negative
              displace=-bxf(j,k)* &
                    (ay*byf(j,k)+az*bzf(j,k))/b_perp**2
             endif
              ar=distance+displace
              future(j,k)=future(j,k)+t_equiv*ar/avx/3600.
            end do
      175   continue
      !
            call set_imf(bx,by,bz,bx0,by0,bz0,bxp,byp,bzp, &
                         qrho,qpres,qpx,qpy,qpz,epres, &
      !     +             orho,opres,opx,opy,opz, &
                         rhop,svxp,svyp,svzp,svelx,spress,ti_te, &
                         rho_frac,nx,ny,nz,ngrd)
         endif
      
      !
      !     check initial conditions
      !
             write(6,*) 'checking set speed'
            do m=ngrd,1,-1
      !
      !       sync time steps
      !
              t_old(m)=0.
              t_new(m)=0.
      !
      !       check density
      !
              call set_rho(qrho,qpres,rmassq, &
                          hrho,hpres,rmassh,orho,opres,rmasso, &
                          epres,nx,ny,nz,ngrd,m,o_conc)
      !
      !       check speeds of individual grids
      !
      
              call set_speed_agrd(qrho,qpres,qpx,qpy,qpz, &
                hrho,hpres,hpx,hpy,hpz,orho,opres,opx,opy,opz, &
                epres,bx,by,bz,bx0,by0,bz0,bsx,bsy,bsz,btot, &
                tvx,tvy,tvz,evx,evy, &
                rmassq,rmassh,rmasso,nx,ny,nz,ngrd,m, &
                pxmax,pymax,pzmax,pmax,csmax,alfmax,gamma, &
                vlim,alf_lim,o_conc,fastest)
              write(6,195)m,csmax,alfmax,pxmax,pymax,pzmax
        195   format(1x,i2,5(1x,1pe12.5))
      !
                t_step_new(m)=stepsz*xspac(m)/fastest
            enddo
            write(6,*)'speeds checked'
      !
            delt_old=delt
            delt=stepsz/fastest
            delt=amin1(delt,1.25*delt_old)
            delt=amax1(3.e-3,delt)
            write(6,163)t,delt,ut,bfld(1,4)
        163 format(1x,'t=',1pe12.5,' dt=',1pe12.5,' ut=', &
                    1pe12.5,' wind time',1pe12.5)
      !
      !     write(6,161)
      ! 161 format(' initialization completed')
      !
      !     goto 519
      !
      
            write(*,*) "Dumping initialization to ", fluidstr//'0', "..."
      
            open(nchf,file=fluidstr//'0',status='new',form='unformatted')
                  write(nchf)t
                  write(nchf)qrho
                  write(nchf)qpx
                  write(nchf)qpy
                  write(nchf)qpz
                  write(nchf)qpres
                  write(nchf)hrho
                  write(nchf)hpx
                  write(nchf)hpy
                  write(nchf)hpz
                  write(nchf)hpres
                  write(nchf)orho
                  write(nchf)opx
                  write(nchf)opy
                  write(nchf)opz
                  write(nchf)opres
                  write(nchf)bx
                  write(nchf)by
                  write(nchf)bz
                  write(nchf)epres
                  write(nchf)bx0
                  write(nchf)by0
                  write(nchf)bz0
                  write(nchf)mass_load_photo
                  write(nchf)mass_load_e
                  write(nchf)photo_x
                  write(nchf)photo_y
                  write(nchf)photo_z
                  write(nchf)e_impact_x
                  write(nchf)e_impact_y
                  write(nchf)e_impact_z
                  write(nchf)elastic_x
                  write(nchf)elastic_y
                  write(nchf)elastic_z
                  write(nchf)parm_srf,parm_mid,parm_zero, &
                     ijzero,numzero,ijmid,nummid,ijsrf,numsrf
            close(nchf)
      
            write(*,*) "Done. Starting main time sequence."
      
      
      !     ********************************************
      !     start time sequence
      !     ********************************************
      !
      !
       1000 continue
	   
	   

	   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	   ! XIN CODE
	   ! spinvectorx=-cos((abs(tilt)-7.9)*.0174533)
       ! spinvectory=0.0
       ! spinvectorz=sin((abs(tilt)-7.9)*.0174533)
	   
	   ! spinvectorx, spinvectory, spinvectorz are respectively the three components of the vector form of the planetary rotation axis of the planet
	   ! for dipole tilt of 47 degrees, so the angle 90 - tilt is 43 degrees 
	   ! 09/15/2020 AJO according to Xin, there needs to be a negative sign in here somewhere, but I don't understand why
	   ! spinvectorx = cos( ( 90 - abs(tilt) ) * .0174533 )
	   ! spinvectory = 0.0
	   ! spinvectorz = sin( ( 90 - abs(tilt) )* .0174533 )
	   
	   ! 09/17/2020 AJO: I don't understand why there is a negative sign in the sin term for spinvectorx. Xin told me it was to do with a rotation 
	   ! transformation but it still doesn't make sense to me, probably because I don't understand the use of superwholecoor as an intermediate variable
	   ! I am incredibly frustrated by the way this code was written. I am blindly going to trust what Xin wrote for now and I'm going to return when I better 
	   ! understand the order in which axes are rotated for this calculation. What Xin told me to use: 
	   spinvectorx = -sin( abs(tilt) * .0174533 )
	   spinvectory = 0.0
       spinvectorz = cos( abs(tilt) * .0174533 )

	   

      !
      !     find maximum velocity to determine time step
      !
      !
            do m=ngrd,1,-1
             t_step(m)=t_step_new(m)
      !
      !       sync time steps
      !
              t_old(m)=0.
              t_new(m)=0.
      
                if(m.eq.ngrd)then
                  smallest_step=t_step(m)
                  mallest_step=m
                   write(6,*)'syncing',m,t_step(m)  
                else
      !
      !          check to see if m grid doesnt outstep grid m+1
                  write(6,*)'syncing',m,t_step(m)
                  if(t_step(m).gt.t_step(m+1))t_step(m)=t_step(m+1)
      !
                  if (smallest_step.gt.t_step(m))then
                    smallest_step=t_step(m)
                    mallest_step=m
                  endif
                endif
            enddo
      !
      !     set variable steps for each grid box
      !       round_step size off
      !
      !       i_step=smallest_step*10000
      !       smallest_step=i_step/10000.
      !
            write(6,*)'unsync steps',t_step,mallest_step,smallest_step
      !
            do m=1,ngrd
              if(m.le.mallest_step) then
                 t_step(m)=smallest_step
                 write(6,*)'sync step',m,t_step(m)
              else
                 astep=(t_step(m)/t_step(m-1)+.50)  !round up as needed
                 nsteps=astep                                  !nearest integer
                 write(6,*)astep,nsteps
                 if(nsteps.gt.2)nsteps=2
                 if(nsteps.lt.1)nsteps=1
                 t_step(m)=t_step(m-1)*nsteps
                 write(6,*)'sync step',m,t_step(m)
              endif
            enddo
            m_step=t_step(ngrd)/t_step(mallest_step)
      !
            write(6,*)'variable steps ',mallest_step,m_step,t_step
      !
            delt=t_step(ngrd) 
            told=t
            utold=ut 
            old_tilt=tilt
      !
            t=t+delt
            ut=utstart+t*t_equiv/3600.
            tilt=tilt+dtilt*delt
            delay=t_equiv*distance/svelx/3600.
      !
            write(6,201)t,delt,ut,bfld(1,4)
        201 format(1x,'t=',1pe12.5,' dt=',1pe12.5,' ut=', &
                    1pe12.5,' wind time',1pe12.5)
      !
      !
      !     write out data if necessary - only using course gridding
      !        the momemt
      !
      !     if(spacecraft) then
            do n=1,ncraft
      
               bsx = 0. !call qvset(0.,bsx,nx*ny*nz)
               bsy = 0. !call qvset(0.,bsy,nx*ny*nz)
               bsz = 0. !call qvset(0.,bsz,nx*ny*nz)
      
               m=1
               do while ((xcraft(1,n).gt.grd_xmax(m)*re_equiv).or. &
                    (xcraft(1,n).lt.grd_xmin(m)*re_equiv).or. &
                    (xcraft(2,n).gt.grd_ymax(m)*re_equiv).or. &
                    (xcraft(2,n).lt.grd_ymin(m)*re_equiv).or. &
                    (xcraft(3,n).gt.grd_zmax(m)*re_equiv).or. &
                    (xcraft(3,n).lt.grd_zmin(m)*re_equiv).and. &
                    (m+1.le.ngrd)) 
                  m=m+1
               enddo
               rx=xspac(m)
               ry=rx
               rz=rz
      !
               call totfld(bx,bx0,bsx,nx,ny,nz,ngrd,m)
               call totfld(by,by0,bsy,nx,ny,nz,ngrd,m)
               call totfld(bz,bz0,bsz,nx,ny,nz,ngrd,m)
      !
               add_dip=.false.
               call crafdatv(bsx,bsy,bsz,qpx,qpy,qpz,qrho,qpres, &
                    rmassq,hpx,hpy,hpz,hrho,hpres,rmassh, &
                    opx,opy,opz,orho,opres,rmasso,epres, &
                    nx,ny,nz,ngrd,m,xcraft,ncraft,n,ut,add_dip, &
                    re_equiv,b_equiv,v_equiv,rho_equiv,gamma)
            enddo
      !      endif
      !
      !     test to see whether scraft positions need to be updated
      !
            if(spacecraft)then
      !      do 210 n=1,ncraft
      !      do 210 n=1,2
             n=1
                if(ut.ge.zcraft(4,n))then
                 mout=20+n
                 read(mout,*)zcraft(4,n),zcraft(1,n), &
                      zcraft(2,n),zcraft(3,n)
      !          if(n.eq.1)zcraft(4,n)=zcraft(4,n)/3600.
      ! 
      !          change direction to get from GSM to simulation coords
      !
                  zcraft(1,n)=-zcraft(1,n)
                  zcraft(2,n)=-zcraft(2,n)
                   if(n.eq.1)then
                    rcraft(1)=zcraft(1,1)
                    rcraft(2)=zcraft(2,1)
                    rcraft(3)=zcraft(3,1)
                   endif
      !
                endif
      ! 210 continue
      !
            zcraft(4,3)=zcraft(4,2)
            zcraft(4,4)=zcraft(4,2)
      !
      !     zcraft(1,2)=zcraft(1,2)
      !     zcraft(2,2)=zcraft(2,2)-1.2
      !     zcraft(3,2)=zcraft(3,2)
      !
      !         set refernce spacecraft position
      !                   and spacecraft limits
      !
                call limcraft(zcraft,ncraft,re_equiv,ngrd)
      !
      !         set density and velocity
      !
               distance=grd_xmin(ngrd)-rcraft(1)/re_equiv
               do while (ut.ge.vut)
                nvx=nvx+1
                zvelx=-svel(nvx,1)/v_equiv
                zvely=-svel(nvx,2)/v_equiv
      !         zvely=-svel(nvx,2)/v_equiv+0.03
                zvelz=svel(nvx,3)/v_equiv
                vut=bfld(nvx,4)+t_equiv*distance/zvelx/3600.
      !          read(28,*)vut,zvelx,zvely,zvelz
      !            zvelx=-zvelx/v_equiv
      !            zvely=-zvely/v_equiv+0.03
      !            zvelz=zvelz/v_equiv
      !            vut=vut+t_equiv*distance/zvelx/3600.
                end do
      !         do while (ut.ge.rut)
      !          read(27,*)rut,zrho
      !            rut=rut+t_equiv*distance/zvelx/3600.
      !            zrho=zrho/rho_equiv
      !         end do
      !
      !        fix up magnetic field
      !
                displace=0.
                do 220 j=1,ny
                do 220 k=1,nz
                 do while((ut.ge.future(j,k)) &
                         .and.(ncount(j,k)+1.le.ncts))
                 nc=ncount(j,k)+1
                 future(j,k)=bfld(nc,4)
                 bxf(j,k)=-bfld(nc,1)/b_equiv
                 byf(j,k)=-bfld(nc,2)/b_equiv
                 bzf(j,k)=bfld(nc,3)/b_equiv
                 rhof(j,k)=rplas(nc)/rho_equiv
                 svxf(j,k)=-svel(nc,1)/v_equiv
                 svyf(j,k)=0.00
                 svzf(j,k)=0.0
      !          svyf(j,k)=-svel(nc,2)/v_equiv
      !          svyf(j,k)=-svel(nc,2)/v_equiv+0.03
      !          svzf(j,k)=svel(nc,3)/v_equiv
                 avx=svxf(j,k)
                 ncount(j,k)=nc
      !
      !      calculate delay
      !
                if(warp)then
      !           b_perp=sqrt(bzf(j,k)**2+byf(j,k)**2)
      !           b_perp=amax1(b_perp,0.1*abs(bxf(j,k)))
                  ay=(j-1.)*xspac(ngrd)+grd_ymin(ngrd)-rcraft(2)/re_equiv
                  az=(k-1.)*xspac(ngrd)+grd_zmin(ngrd)-rcraft(3)/re_equiv
      !
      !       going to assume Bz IMF on average pos and ignore transients
      !               and By IMF is negative
      !
                   b_perp=sqrt(bzf(j,k)**2+byf(j,k)**2) 
                   b_perp=amax1(b_perp,0.33*abs(bxf(j,k))) 
                   displace=-bxf(j,k)* &
                    (ay*byf(j,k)+az*bzf(j,k))/b_perp**2
                 endif
                 ar=distance+displace
                 future(j,k)=future(j,k)+t_equiv*ar/avx/3600.
                 end do
        220   continue
      !
            endif
      !
            if(spacecraft)then
      !
      !     update spacecraft position and magnetic field
      !
            do 250 n=1,ncraft
               dut=(ut-xcraft(4,n))/(zcraft(4,n)-utold)
               xcraft(1,n)=xcraft(1,n)+dut*(zcraft(1,n)-xcraft(1,n))
               xcraft(2,n)=xcraft(2,n)+dut*(zcraft(2,n)-xcraft(2,n))
               xcraft(3,n)=xcraft(3,n)+dut*(zcraft(3,n)-xcraft(3,n))
               xcraft(4,n)=ut
        250 continue
              distance=grd_xmin(ngrd)-rcraft(1)/re_equiv
              call limcraft(xcraft,ncraft,re_equiv,ngrd)
      !
             srho=0.
             do j=1,ny
              do k=1,nz
               dut=(ut-utold)/(future(j,k)-utold)
               bxp(j,k)=bxp(j,k)+dut*(bxf(j,k)-bxp(j,k))
               byp(j,k)=byp(j,k)+dut*(byf(j,k)-byp(j,k))
               bzp(j,k)=bzp(j,k)+dut*(bzf(j,k)-bzp(j,k))
               rhop(j,k)=rhop(j,k)+dut*(rhof(j,k)-rhop(j,k))
               svxp(j,k)=svxp(j,k)+dut*(svxf(j,k)-svxp(j,k))
               svyp(j,k)=svyp(j,k)+dut*(svyf(j,k)-svyp(j,k))
               svzp(j,k)=svzp(j,k)+dut*(svzf(j,k)-svzp(j,k))
               srho=srho+rhop(j,k)
              end do
             end do
      !
             dut=(ut-utold)/(vut-utold)
             svelx=svelx+dut*(zvelx-svelx)
             svely=0.
             svelz=svelz+delvz_wind*delt
             srho=srho/float(nz*ny)
      !
             sbx_wind=sbx_wind+delbx_wind*delt
             sby_wind=sby_wind+delby_wind*delt
             sbz_wind=sbz_wind+delbz_wind*delt
             spx=srho*svelx
             spy=srho*svely
             spz=srho*svelz
             spress=(cs_wind**2*srho/gamma)/gamma1
             serg=0.5*(svelx**2+svely**2+svelz**2)*srho+spress
            else
             sbx_wind=sbx_wind+delbx_wind*delt
             sby_wind=sby_wind+delby_wind*delt
             sbz_wind=sbz_wind+delbz_wind*delt
      !  
             svelx=svelx+delvx_wind*delt 
             svely=svely+delvy_wind*delt 
             svelz=svelz+delvz_wind*delt
             srho=srho+delrho*delt
             spx=srho*svelx
             spy=srho*svely
             spz=srho*svelz
             spress=(cs_wind**2*srho/gamma)/gamma1
             serg=0.5*(svelx**2+svely**2+svelz**2)*srho+spress
            endif
            delay=t_equiv*distance/svelx/3600.
      !
      !       test for div B errors and reset ionospheric restivity
      !
      !     if((.not.start).and.(t.gt.tdiv))then
      !     if (t.gt.tdiv)then
            if(divb_on)then
               tdiv=tdiv+del_tdiv
               do m=ngrd,1,-1
                call divb_cor(bx,by,bz,qpx,qpy,qpz,qrho,qpres, &
                          ngrd,nx,ny,nz,m,srho)
                if(m.gt.1) &
                    call flanks_synced(bx,ngrd,nx,ny,nz,m-1)
               enddo
             endif
      
      !
      !     start variable time loop over delt/m_step
      !
            do ms=1,m_step
            do m=ngrd,1,-1
      !
      !     test if grid sector needs to be moved in time
      !     
            yes_step=.false.
      !
            if(m.eq.1)then
               yes_step=.true.
            else
             if ((t_old(m).eq.t_new(m)).and.(t_new(m).lt.t_step(ngrd)) &
                .and.(abs(t_new(m)-t_new(1)).le.0.001*t_step(1)) )then
               yes_step=.true.
             endif
            endif
      !
      !      time step grid
      !
            if(yes_step) then
             t_old(m)=t_new(m)
             t_new(m)=t_old(m)+t_step(m)
      !
             write(6,*)'time stepping',m, t_old(m),t_new(m)
      !
             delt= t_step(m)
             delt2=delt/2.
      !
             if(tilting)then 
              atilt=old_tilt+(t_old(m)+delt2)*dtilt
              sin_tilt=sin(atilt*.0174533)
              cos_tilt=cos(atilt*.0174533)
      !       write(6,*)'mak_dip with', rot_angle,ut2
              call  mak_dip_grd(bx0,by0,bz0,nx,ny,nz, &
                          ngrd,m,ijzero,numzero)
             endif
      !
      !     *******************************************************
      !     store initial plasma parameters
      !     ******************************************************
      !
      !     store initial plasma parameters
            oldqrho(:,:,:,m) = qrho(:,:,:,m)
            oldqpres(:,:,:,m) = qpres(:,:,:,m)
            oldqpx(:,:,:,m) = qpx(:,:,:,m)
            oldqpy(:,:,:,m) = qpy(:,:,:,m)
            oldqpz(:,:,:,m) = qpz(:,:,:,m)
      
            oldhrho(:,:,:,m) = hrho(:,:,:,m)
            oldhpres(:,:,:,m) = hpres(:,:,:,m)
            oldhpx(:,:,:,m) = hpx(:,:,:,m)
            oldhpy(:,:,:,m) = hpy(:,:,:,m)
            oldhpz(:,:,:,m) = hpz(:,:,:,m)
      
            oldorho(:,:,:,m) = orho(:,:,:,m)
            oldopres(:,:,:,m) = opres(:,:,:,m)
            oldopx(:,:,:,m) = opx(:,:,:,m)
            oldopy(:,:,:,m) = opy(:,:,:,m)
            oldopz(:,:,:,m) = opz(:,:,:,m)
      
            oldepres(:,:,:,m) = epres(:,:,:,m)
            oldbx(:,:,:,m) = bx(:,:,:,m)
            oldby(:,:,:,m) = by(:,:,:,m)
            oldbz(:,:,:,m) = bz(:,:,:,m)
      
      !     **********************************************************
      !     Two Step Runge-Kutta : step 1
      !          estimate of the fluid quantites at n+1/2 
      !     **********************************************************
      !
      !     store initial plasma parameters
      !
      !     calculate standard mhd current j = curl B
      !
      !
             rx=xspac(m)
             ry=xspac(m)
             rz=xspac(m)
      !
      !      write(6,)' calcur ing now'
             call calcur(bx,by,bz,ngrd,nx,ny,nz,m,curx,cury,curz, &
                     ijsrf,numsrf,ijmid,nummid,ijzero,numzero)
      !
      !     find total magnetic field
      !
      !      write(6,*)' totbfld ing now'
            call totfld(bx,bx0,bsx,nx,ny,nz,ngrd,m)
            call totfld(by,by0,bsy,nx,ny,nz,ngrd,m) 
            call totfld(bz,bz0,bsz,nx,ny,nz,ngrd,m)
      !
      !     find magnitude of B
      !
            call tot_b(btot,bsx,bsy,bsz,nx,ny,nz)
      !
      !
      !      write(6,*)' fnd_evel ing now'
             call fnd_evel(qpx,qpy,qpz,qrho,hpx,hpy,hpz,hrho,  &
                   opx,opy,opz,orho,curx,cury,curz,evx,evy,evz, &
                   tvx,tvy,tvz,ngrd, &
                   nx,ny,nz,m,rmassq,rmassh,rmasso,reynolds)
      !
      !      write(6,*)' bande ing now'
            call bande(efldx,efldy,efldz,bsx,bsy,bsz, &
                   curx,cury,curz,evx,evy,evz,btot, &
                   epres,qrho,hrho,orho,resistive,resist,reynolds, &
                   nx,ny,nz,ngrd,m,rmassq,rmassh,rmasso, &
                   ijmid,nummid,ijzero,numzero)
      
            !elecden(:,:,:) = oldqrho(:,:,:,1) + &
            !     oldhrho(:,:,:,1)/rmassh + &
            !     oldorho(:,:,:,1)/rmasso
      
      !      OPEN(81,file='eradial1.dat',status='unknown',form='formatted') 
      !      DO ii = 51,nx
      ! 785     FORMAT(E13.5,' ',E13.5)
      !         WRITE(81,785) pchar*epres(1,ii,51,25), nchar*elecden(ii,51,25)
      !      END DO
      !      CLOSE(81)
      !
      !      write(6,*)' push elec ing now'
             call push_elec(wrkepres,oldepres,epres,evx,evy,evz, &
                    gamma,gamma1,ngrd,nx,ny,nz,m,0.5*delt)
      
      !      write(6,*)' push ion 1 ing now'
             call push_ion(wrkqrho,wrkqpres,wrkqpx,wrkqpy,wrkqpz, &
                    oldqrho,oldqpres,oldqpx,oldqpy,oldqpz, &
                    qrho,qpres,qpx,qpy,qpz, &
                    bsx,bsy,bsz,btot,efldx,efldy,efldz,rmassq, &
                    vvx,vvy,vvz,tvx,tvy,tvz,gamma,gamma1, &
                    ngrd,nx,ny,nz,m,0.5*delt,grav,re_equiv,reynolds)
      
             write(6,*)' push ion 2 ing now'
      
             !CALL TIMEIT(0)
      !     A. Rajendar, selection of push_ion2 in box 1, removed from subroutine
             IF (m.eq.10) THEN
                  write(*,*) "******************"
                  write(*,*) "YOU SHALL NOT PASS"
                  write(*,*) "******************"
      !     A. Rajendar, modifed my inclusion of push_ion2 01/23/2013.
      !     A.R., 11/04/2013: modified push_ion2 to include e-impact ionization
      !     by tail of thermal electron population
      !
      !     A.R., 11/04/2013: caculates the density of thermal population of
      !     ionizing electrons in box 1
                elecden(:,:,:) = &
                     wrkqrho(:,:,:,m)/rmassq + &
                     wrkhrho(:,:,:,m)/rmassh + &
                     wrkorho(:,:,:,m)/rmasso
      
                CALL elec_impact_frq(ngrd,nx,ny,nz,epres,elecden,eifrq,hoteifrq, &
                     flen,eev,kappa,pchar,nchar,tchar,lchar,re_equiv,rho_equiv)
      
                call push_ion(wrkhrho,wrkhpres,wrkhpx,wrkhpy,wrkhpz, &
                     oldhrho,oldhpres,oldhpx,oldhpy,oldhpz, &
                     hrho,hpres,hpx,hpy,hpz, &
                     bsx,bsy,bsz,btot,efldx,efldy,efldz,rmassh, &
                     vvx,vvy,vvz,tvx,tvy,tvz,gamma,gamma1, &
                     ngrd,nx,ny,nz,m,0.5*delt,grav,re_equiv,reynolds)
             ELSE
                call push_ion(wrkhrho,wrkhpres,wrkhpx,wrkhpy,wrkhpz, &
                     oldhrho,oldhpres,oldhpx,oldhpy,oldhpz, &
                     hrho,hpres,hpx,hpy,hpz, &
                     bsx,bsy,bsz,btot,efldx,efldy,efldz,rmassh, &
                     vvx,vvy,vvz,tvx,tvy,tvz,gamma,gamma1, &
                     ngrd,nx,ny,nz,m,0.5*delt,grav,re_equiv,reynolds)
             END IF
             !CALL TIMEIT(1,ttime)
      904    FORMAT('Step, 1, box,',I2,', ut(hrs),',F8.3,', looptime(s),',E10.3)
             WRITE(501,904) m, ut, ttime
      
             ! 2015-08-15, AR: Track motion of W+ fluid particles, Euler step 1
             !IF ((m.EQ.trackbox).AND.(fluid_track)) THEN
             !   ! velocity field @ time t_n
             !   velx = oldhpx(:,:,:,m)/oldhrho(:,:,:,m)
             !   vely = oldhpy(:,:,:,m)/oldhrho(:,:,:,m)
             !   velz = oldhpz(:,:,:,m)/oldhrho(:,:,:,m)
             !   DO i = 1, numpart
             !      ! velocity at fluid particle location @ t_n
             !      CALL interpol3d(nx,ny,nz,ngrd,re_equiv,velx,fluidxyz(:,i),m,fluidvel(1))
             !      CALL interpol3d(nx,ny,nz,ngrd,re_equiv,vely,fluidxyz(:,i),m,fluidvel(2))
             !      CALL interpol3d(nx,ny,nz,ngrd,re_equiv,velz,fluidxyz(:,i),m,fluidvel(3))
             !      ! particle location @ t_(n + 1/2) [intermediate location]
             !      fluidxyz1(:,i) = fluidxyz(:,i) + 0.5*delt*fluidvel*re_equiv
             !   END DO
             !END IF
             
             call push_ion(wrkorho,wrkopres,wrkopx,wrkopy,wrkopz, &
                    oldorho,oldopres,oldopx,oldopy,oldopz, &
                    orho,opres,opx,opy,opz, &
                    bsx,bsy,bsz,btot,efldx,efldy,efldz,rmasso, &
                    vvx,vvy,vvz,tvx,tvy,tvz,gamma,gamma1, &
                    ngrd,nx,ny,nz,m,0.5*delt,grav,re_equiv,reynolds)
      !
             call push_bfld(wrkbx,wrkby,wrkbz,oldbx,oldby,oldbz, &
                          efldx,efldy,efldz,ngrd,nx,ny,nz,m,0.5*delt) 
      !
      !     write(6,988)
      ! 988 format(' main lax loop now')
      !
      !     *************************************************************
      !     Apply boundary conditions
      !     *************************************************************
      !
      !     write(6,989)
      ! 989 format(' doing bndry conditions')
      !
             if(m.eq.ngrd) then
                call bndry_outer(wrkqrho,wrkqpres,wrkqpx,wrkqpy,wrkqpz, &
                     wrkhrho,wrkhpres,wrkhpx,wrkhpy,wrkhpz, &
                     wrkorho,wrkopres,wrkopx,wrkopy,wrkopz,wrkepres, &
                     rmassq,rmassh,rmasso,wrkbx,wrkby,wrkbz, &
                     bx0,by0,bz0,ngrd,nx,ny,nz, &
                     srho,rho_frac,o_conc,spress,spx,spy,spz, &
                     sbx_wind,sby_wind,sbz_wind,ti_te)
             else
                t_grd=t_old(m)+delt2
                if(t_grd.gt.t_new(m+1))then
                   t_grd=t_new(m+1)
                endif
                if (t_grd.lt.t_old(m+1))then
                   t_grd=t_old(m+1)
                endif
      !      write(6,*) 'flanks1', m,t_grd, t_old(m+1),t_new(m+1)
                call bndry_flanks(wrkqrho,wrkqpres,wrkqpx,wrkqpy,wrkqpz, &
                     wrkhrho,wrkhpres,wrkhpx,wrkhpy,wrkhpz, &
                     wrkorho,wrkopres,wrkopx,wrkopy,wrkopz, &
                     wrkepres,wrkbx,wrkby,wrkbz,   &
                     qrho,qpres,qpx,qpy,qpz,  &
                     hrho,hpres,hpx,hpy,hpz, &
                     orho,opres,opx,opy,opz, &
                     epres,bx,by,bz,    &
                     oldqrho,oldqpres,oldqpx,oldqpy,oldqpz,  &
                     oldhrho,oldhpres,oldhpx,oldhpy,oldhpz, &
                     oldorho,oldopres,oldopx,oldopy,oldopz, &
                     oldepres,oldbx,oldby,oldbz,  &
                     nx,ny,nz,ngrd,m,t_old,t_new,t_grd)
             endif 
             if(m.EQ.1)then
                call bndry_inner( &
                     wrkqrho,wrkqpres,wrkqpx,wrkqpy,wrkqpz,rmassq, &
                     wrkhrho,wrkhpres,wrkhpx,wrkhpy,wrkhpz,rmassh, &
                     wrkorho,wrkopres,wrkopx,wrkopy,wrkopz,rmasso, &
                     wrkepres, wrkbx, wrkby, wrkbz, &
                     ngrd,nx,ny,nz,parm_srf,parm_mid,parm_zero, &
                     ijsrf,numsrf,ijmid,nummid,ijzero, &
                     numzero,erho,epress,alpha_e,ti_te,o_conc, &
                     sbx_wind,spress)
             endif
      !
      !      write(6,*)'Lax 1 speeds'
      !
      !     check fluid parameters
      ! 
            if(spacecraft)then
               call set_imf(wrkbx,wrkby,wrkbz,bx0,by0,bz0,bxp,byp,bzp, &
                     wrkqrho,wrkqpres,wrkqpx,wrkqpy,wrkqpz,wrkepres, &
                     rhop,svxp,svyp,svzp,svelx,spress,ti_te, &
                     rho_frac,nx,ny,nz,ngrd)
            endif 
      
            call set_rho(wrkqrho,wrkqpres,rmassq, &
                        wrkhrho,wrkhpres,rmassh,wrkorho,wrkopres,rmasso, &
                        wrkepres,nx,ny,nz,ngrd,m,o_conc)
      !
            call set_speed_agrd(wrkqrho,wrkqpres,wrkqpx,wrkqpy,wrkqpz, &
                wrkhrho,wrkhpres,wrkhpx,wrkhpy,wrkhpz, &
                wrkorho,wrkopres,wrkopx,wrkopy,wrkopz, &
                wrkepres,wrkbx,wrkby,wrkbz,bx0,by0,bz0, &
                bsx,bsy,bsz,btot,tvx,tvy,tvz,evx,evy, &
                rmassq,rmassh,rmasso,nx,ny,nz,ngrd,m, &
                pxmax,pymax,pzmax,pmax,csmax,alfmax,gamma, &
                vlim,alf_lim,o_conc,fastest)
      
      !
      !
      !      ***********************************************************
      !      Runge-Kutta: step 2
      !            use the predicted value to find corrected value for n+1
      !      ***********************************************************
      ! 
            if(tilting)then
              atilt=old_tilt+(t_old(m)+delt)*dtilt ! precessing the tilt around in time based on time step in code     
              sin_tilt=sin(atilt*.0174533)
              cos_tilt=cos(atilt*.0174533)
      !       write(6,*)'mak_dip with', rot_angle,aut2
              call  mak_dip_grd(bx0,by0,bz0,nx,ny,nz, &
                          ngrd,m,ijzero,numzero)
      
            endif
            rx=xspac(m)
            ry=xspac(m)
            rz=xspac(m)
      !
      !     calculate standard mhd current j = curl B
      !
            call calcur(wrkbx,wrkby,wrkbz,ngrd,nx,ny,nz,m,curx,cury,curz, &
                        ijsrf,numsrf,ijmid,nummid,ijzero,numzero)
      !
      !     find total magnetic field
      !
            call totfld(wrkbx,bx0,bsx,nx,ny,nz,ngrd,m)
            call totfld(wrkby,by0,bsy,nx,ny,nz,ngrd,m)
            call totfld(wrkbz,bz0,bsz,nx,ny,nz,ngrd,m)
      !
      !     find magnitude of B
      !
            call tot_b(btot,bsx,bsy,bsz,nx,ny,nz)
      !
      !
      !     find the  electric field from electron momentum eqn
      !
             call fnd_evel(wrkqpx,wrkqpy,wrkqpz,wrkqrho, &
                   wrkhpx,wrkhpy,wrkhpz,wrkhrho,  &
                   wrkopx,wrkopy,wrkopz,wrkorho, &
                   curx,cury,curz,evx,evy,evz, &
                   tvx,tvy,tvz,ngrd, &
                   nx,ny,nz,m,rmassq,rmassh,rmasso,reynolds)
      !
            call bande(efldx,efldy,efldz,bsx,bsy,bsz, &
                   curx,cury,curz,evx,evy,evz,btot, &
                   wrkepres,wrkqrho,wrkhrho,wrkorho, &
                   resistive,resist,reynolds, &
                   nx,ny,nz,ngrd,m,rmassq,rmassh,rmasso, &
                   ijmid,nummid,ijzero,numzero)
      !
            call push_elec(epres,oldepres,wrkepres,evx,evy,evz, &
                 gamma,gamma1,ngrd,nx,ny,nz,m,delt)
      !
            call push_ion(qrho,qpres,qpx,qpy,qpz, &
                 oldqrho,oldqpres,oldqpx,oldqpy,oldqpz, &
                 wrkqrho,wrkqpres,wrkqpx,wrkqpy,wrkqpz, &
                 bsx,bsy,bsz,btot,efldx,efldy,efldz,rmassq, &
                 vvx,vvy,vvz,tvx,tvy,tvz,gamma,gamma1, &
                 ngrd,nx,ny,nz,m,delt,grav,re_equiv,reynolds)
            
            !CALL TIMEIT(0)
      !     A. Rajendar, applies push_ion2 to box 1 only
            IF (m.EQ.10) THEN
                  write(*,*) "******************"
                  write(*,*) "YOU SHALL NOT PASS"
                  write(*,*) "******************"
      !     A. Rajendar, modifed my inclusion of push_ion2 01/23/2013.
      !
      !     A.R., 11/04/2013: caculates the density of thermal population of
      !     ionizing electrons in box 1
               elecden(:,:,:) = &
                    wrkqrho(:,:,:,m)/rmassq + &
                    wrkhrho(:,:,:,m)/rmassh + &
                    wrkorho(:,:,:,m)/rmasso
      
               CALL elec_impact_frq(ngrd,nx,ny,nz,epres,elecden,eifrq,hoteifrq, &
                    flen,eev,kappa,pchar,nchar,tchar,lchar,re_equiv,rho_equiv)
      
               call push_ion(wrkhrho,wrkhpres,wrkhpx,wrkhpy,wrkhpz, &
                     oldhrho,oldhpres,oldhpx,oldhpy,oldhpz, &
                     hrho,hpres,hpx,hpy,hpz, &
                     bsx,bsy,bsz,btot,efldx,efldy,efldz,rmassh, &
                     vvx,vvy,vvz,tvx,tvy,tvz,gamma,gamma1, &
                     ngrd,nx,ny,nz,m,0.5*delt,grav,re_equiv,reynolds)
            ELSE
               call push_ion(hrho,hpres,hpx,hpy,hpz, &
                    oldhrho,oldhpres,oldhpx,oldhpy,oldhpz, &
                    wrkhrho,wrkhpres,wrkhpx,wrkhpy,wrkhpz, &
                    bsx,bsy,bsz,btot,efldx,efldy,efldz,rmassh, &
                    vvx,vvy,vvz,tvx,tvy,tvz,gamma,gamma1, &
                    ngrd,nx,ny,nz,m,delt,grav,re_equiv,reynolds)
            END IF
            !CALL TIMEIT(1,ttime)
      905   FORMAT('Step, 2, box,',I2,', ut(hrs),',F8.3,', looptime(s),',E10.3)
            WRITE(501,905) m, ut, ttime
      
             ! 2015-08-15, AR: Track motion of W+ fluid particles, Euler step 2
             !IF ((m.EQ.trackbox).AND.(fluid_track)) THEN
             !   ! velocity field @ t_(n + 1/2)
             !   velx = wrkhpx(:,:,:,m)/wrkhrho(:,:,:,m)
             !   vely = wrkhpy(:,:,:,m)/wrkhrho(:,:,:,m)
             !   velz = wrkhpz(:,:,:,m)/wrkhrho(:,:,:,m)
             !   DO i = 1, numpart
             !      ! velocity at intermediate location @ t_(n + 1/2)
             !      CALL interpol3d(nx,ny,nz,ngrd,re_equiv,velx,fluidxyz1(:,i),m,fluidvel(1))
             !      CALL interpol3d(nx,ny,nz,ngrd,re_equiv,vely,fluidxyz1(:,i),m,fluidvel(2))
             !      CALL interpol3d(nx,ny,nz,ngrd,re_equiv,velz,fluidxyz1(:,i),m,fluidvel(3))
             !      ! fluid particle location @ t_(n + 1)
             !      fluidxyz(:,i) = fluidxyz1(:,i) + 0.5*delt*fluidvel*re_equiv
             !   END DO
             !END IF
      
      
            call push_ion(orho,opres,opx,opy,opz, &
                 oldorho,oldopres,oldopx,oldopy,oldopz, &
                 wrkorho,wrkopres,wrkopx,wrkopy,wrkopz, &
                 bsx,bsy,bsz,btot,efldx,efldy,efldz,rmasso, &
                 vvx,vvy,vvz,tvx,tvy,tvz,gamma,gamma1, &
                 ngrd,nx,ny,nz,m,delt,grav,re_equiv,reynolds)
      !     
            call push_bfld(bx,by,bz,oldbx,oldby,oldbz, &
                 efldx,efldy,efldz,ngrd,nx,ny,nz,m,delt) 
      ! 
      !
      !     write(6,992)
      ! 992 format(' main 2nd lax loop now')
      !
      !     *************************************************************
      !     Apply boundary conditions
      !     *************************************************************
      !
      !
            if(m.eq.ngrd) then
                 call bndry_outer(qrho,qpres,qpx,qpy,qpz, &
                    hrho,hpres,hpx,hpy,hpz, &
                    orho,opres,opx,opy,opz,epres, &
                    rmassq,rmassh,rmasso,bx,by,bz, &
                    bx0,by0,bz0,ngrd,nx,ny,nz, &
                    srho,rho_frac,o_conc,spress,spx,spy,spz, &
                    sbx_wind,sby_wind,sbz_wind,ti_te)
            else
             t_grd=t_old(m)+delt
             if(t_grd.gt.t_new(m+1))then
      !        write(6,*)'WARNING on lax2 max',m,t_grd,t_new(m+1),
      !    +                         t_grd-t_new(m+1)
               t_grd=t_new(m+1)
             endif
             if (t_grd.lt.t_old(m+1))then
               t_grd=t_old(m+1)
             endif
      !      write(6,*)'flanks2', m,t_grd, t_old(m+1),t_new(m+1)
             call bndry_flanks(qrho,qpres,qpx,qpy,qpz, &
                   hrho,hpres,hpx,hpy,hpz, &
                   orho,opres,opx,opy,opz, &
                   epres,bx,by,bz,     &
                   qrho,qpres,qpx,qpy,qpz,  &
                   hrho,hpres,hpx,hpy,hpz, &
                   orho,opres,opx,opy,opz, &
                   epres,bx,by,bz,    &
                   oldqrho,oldqpres,oldqpx,oldqpy,oldqpz,  &
                   oldhrho,oldhpres,oldhpx,oldhpy,oldhpz, &
                   oldorho,oldopres,oldopx,oldopy,oldopz, &
                   oldepres,oldbx,oldby,oldbz,  &
                   nx,ny,nz,ngrd,m,t_old,t_new,t_grd)
            endif 
      !
            if(m.eq.1)then
             call bndry_inner(qrho,qpres,qpx,qpy,qpz,rmassq, &
                   hrho,hpres,hpx,hpy,hpz,rmassh, &
                   orho,opres,opx,opy,opz,rmasso,epres,bx,by,bz, &
                   ngrd,nx,ny,nz,parm_srf,parm_mid,parm_zero, &
                   ijsrf,numsrf,ijmid,nummid,ijzero, &
                   numzero,erho,epress,alpha_e,ti_te,o_conc, &
                   sbx_wind,spress)
            endif
      !
            if(spacecraft)then
               call set_imf(bx,by,bz,bx0,by0,bz0,bxp,byp,bzp, &
                         qrho,qpres,qpx,qpy,qpz,epres, &
                         rhop,svxp,svyp,svzp,svelx,spress,ti_te, &
                         rho_frac,nx,ny,nz,ngrd)
            endif
      !
      !      write(6,*)'lax 2 speeds'
      !
            call set_rho(qrho,qpres,rmassq, &
                          hrho,hpres,rmassh,orho,opres,rmasso, &
                          epres,nx,ny,nz,ngrd,m,o_conc)
      !
            call set_speed_agrd(qrho,qpres,qpx,qpy,qpz, &
                hrho,hpres,hpx,hpy,hpz,orho,opres,opx,opy,opz, &
                epres,bx,by,bz,bx0,by0,bz0,bsx,bsy,bsz,btot, &
                tvx,tvy,tvz,evx,evy, &
                rmassq,rmassh,rmasso,nx,ny,nz,ngrd,m, &
                pxmax,pymax,pzmax,pmax,csmax,alfmax,gamma, &
                vlim,alf_lim,o_conc,fastest) 
      !     .......................................................
      !     try Lapdius smoothing - smoothed results will appear in nt2
      !     .......................................................
      !
      !     write(6,*)' doing smoothing okay'
             rx=xspac(m)
             ry=xspac(m)
             rz=xspac(m)
      !      write(6,*)'calling lapidus'
            call lapidus(qrho,qpres,qpx,qpy,qpz, &
                    wrkqrho,wrkqpres,wrkqpx,wrkqpy,wrkqpz, &
                    hrho,hpres,hpx,hpy,hpz, &
                    wrkhrho,wrkhpres,wrkhpx,wrkhpy,wrkhpz, &
                    orho,opres,opx,opy,opz, &
                    wrkorho,wrkopres,wrkopx,wrkopy,wrkopz, &
                    epres,wrkepres,bx,by,bz, &
                    wrkbx,wrkby,wrkbz,curx,cury,curz, &
                    nx,ny,nz,ngrd,m,chirho,chipxyz,chierg,delt, &
                    rmassq,rmassh,rmasso)
      !
      !
      !     reset boundary conditions
      !
             if(m.eq.ngrd) then
                 call bndry_outer(wrkqrho,wrkqpres,wrkqpx,wrkqpy,wrkqpz, &
                    wrkhrho,wrkhpres,wrkhpx,wrkhpy,wrkhpz, &
                    wrkorho,wrkopres,wrkopx,wrkopy,wrkopz,wrkepres, &
                    rmassq,rmassh,rmasso,wrkbx,wrkby,wrkbz, &
                    bx0,by0,bz0,ngrd,nx,ny,nz, &
                    srho,rho_frac,o_conc,spress,spx,spy,spz, &
                    sbx_wind,sby_wind,sbz_wind,ti_te)
            else
             call bndry_flanks(wrkqrho,wrkqpres,wrkqpx,wrkqpy,wrkqpz, &
                   wrkhrho,wrkhpres,wrkhpx,wrkhpy,wrkhpz, &
                   wrkorho,wrkopres,wrkopx,wrkopy,wrkopz, &
                   wrkepres,wrkbx,wrkby,wrkbz,     &
                   qrho,qpres,qpx,qpy,qpz,  &
                   hrho,hpres,hpx,hpy,hpz, &
                   orho,opres,opx,opy,opz, &
                   epres,bx,by,bz,    &
                   oldqrho,oldqpres,oldqpx,oldqpy,oldqpz,  &
                   oldhrho,oldhpres,oldhpx,oldhpy,oldhpz, &
                   oldorho,oldopres,oldopx,oldopy,oldopz, &
                   oldepres,oldbx,oldby,oldbz,  &
                   nx,ny,nz,ngrd,m,t_old,t_new,t_grd)
            endif  
            if(m.eq.1)then
             call bndry_inner( &
                   wrkqrho,wrkqpres,wrkqpx,wrkqpy,wrkqpz,rmassq, &
                   wrkhrho,wrkhpres,wrkhpx,wrkhpy,wrkhpz,rmassh, &
                   wrkorho,wrkopres,wrkopx,wrkopy,wrkopz,rmasso, &
                   wrkepres, wrkbx, wrkby, wrkbz, &
                   ngrd,nx,ny,nz,parm_srf,parm_mid,parm_zero, &
                   ijsrf,numsrf,ijmid,nummid,ijzero, &
                   numzero,erho,epress,alpha_e,ti_te,o_conc, &
                   sbx_wind,spress)
            endif
      !
            if(spacecraft)then
               call set_imf(wrkbx,wrkby,wrkbz,bx0,by0,bz0,bxp,byp,bzp, &
                     wrkqrho,wrkqpres,wrkqpx,wrkqpy,wrkqpz,wrkepres, &
                     rhop,svxp,svyp,svzp,svelx,spress,ti_te, &
                     rho_frac,nx,ny,nz,ngrd)
            endif 
            call set_rho(wrkqrho,wrkqpres,rmassq, &
                        wrkhrho,wrkhpres,rmassh,wrkorho,wrkopres,rmasso, &
                        wrkepres,nx,ny,nz,ngrd,m,o_conc)
      !
      !      write(6,*)'lapidus speeds'
      ! 
            call set_speed_agrd(wrkqrho,wrkqpres,wrkqpx,wrkqpy,wrkqpz, &
                wrkhrho,wrkhpres,wrkhpx,wrkhpy,wrkhpz, &
                wrkorho,wrkopres,wrkopx,wrkopy,wrkopz, &
                wrkepres,wrkbx,wrkby,wrkbz,bx0,by0,bz0, &
                bsx,bsy,bsz,btot,tvx,tvy,tvz,evx,evy, &
                rmassq,rmassh,rmasso,nx,ny,nz,ngrd,m, &
                pxmax,pymax,pzmax,pmax,csmax,alfmax,gamma, &
                vlim,alf_lim,o_conc,fastest)
      !
	  
	  
	  
	  
	  

	  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	  ! XIN CHANGES 
	        call bndry_inner_rot(rmassq,rmassh,rmasso, &
            gamma1,eerg,ngrd,nx,ny,nz,parm_srf,parm_mid, &
            ijsrf,numsrf,ijmid,nummid,erho,ti_te,o_conc)
! c
      call rigid_IM_rot(rmassq,rmassh,rmasso, &
             gamma1,eerg,ngrd,nx,ny,nz,qrho,hrho,orho, &
             erho,qpres,hpres,opres,epres,ti_te,o_conc,ut, &
             utdec_start,re_equiv,bx0,by0,bz0)
	  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	  

	  
      !     write(6,994)
      ! 994 format(' lapidus done')
      !
      !
      !     .......................................................
      !     add a little bit of flux correction  : 
      !     ........................................................
      !
      !
      !      chifcs=.125
	     chifcs=difrho
      !
             rx=xspac(m)
             ry=xspac(m)
             rz=xspac(m)
             call fcsmooth(qrho,oldqrho,wrkqrho,ngrd,nx,ny,nz,m,chifcs, &
                    vvx,vvy,vvz)
             call fcsmooth(qpres,oldqpres,wrkqpres,ngrd,nx,ny,nz,m,chifcs, &
                    vvx,vvy,vvz)
             call fcsmooth(qpx,oldqpx,wrkqpx,ngrd,nx,ny,nz,m,chifcs, &
                    vvx,vvy,vvz)
             call fcsmooth(qpy,oldqpy,wrkqpy,ngrd,nx,ny,nz,m,chifcs, &
                    vvx,vvy,vvz)
             call fcsmooth(qpz,oldqpz,wrkqpz,ngrd,nx,ny,nz,m,chifcs, &
                    vvx,vvy,vvz)
      !
             call fcsmooth(hrho,oldhrho,wrkhrho,ngrd,nx,ny,nz,m,chifcs, &
                    vvx,vvy,vvz)
             call fcsmooth(hpres,oldhpres,wrkhpres,ngrd,nx,ny,nz,m,chifcs, &
                    vvx,vvy,vvz)
             call fcsmooth(hpx,oldhpx,wrkhpx,ngrd,nx,ny,nz,m,chifcs, &
                    vvx,vvy,vvz)
             call fcsmooth(hpy,oldhpy,wrkhpy,ngrd,nx,ny,nz,m,chifcs, &
                    vvx,vvy,vvz)
             call fcsmooth(hpz,oldhpz,wrkhpz,ngrd,nx,ny,nz,m,chifcs, &
                    vvx,vvy,vvz)
      !
             call fcsmooth(orho,oldorho,wrkorho,ngrd,nx,ny,nz,m,chifcs, &
                    vvx,vvy,vvz)
             call fcsmooth(opres,oldopres,wrkopres,ngrd,nx,ny,nz,m,chifcs, &
                    vvx,vvy,vvz)
             call fcsmooth(opx,oldopx,wrkopx,ngrd,nx,ny,nz,m,chifcs, &
                    vvx,vvy,vvz)
             call fcsmooth(opy,oldopy,wrkopy,ngrd,nx,ny,nz,m,chifcs, &
                    vvx,vvy,vvz)
             call fcsmooth(opz,oldopz,wrkopz,ngrd,nx,ny,nz,m,chifcs, &
                    vvx,vvy,vvz)
      !
             call fcsmooth(epres,oldepres,wrkepres,ngrd,nx,ny,nz,m,chifcs, &
                    vvx,vvy,vvz)
      !
      
             chifcs=diferg
      !
             call fcsmooth(bx,oldbx,wrkbx,ngrd,nx,ny,nz,m,chifcs, &
                    vvx,vvy,vvz)
             call fcsmooth(by,oldby,wrkby,ngrd,nx,ny,nz,m,chifcs, &
                    vvx,vvy,vvz)
             call fcsmooth(bz,oldbz,wrkbz,ngrd,nx,ny,nz,m,chifcs, &
                    vvx,vvy,vvz)
      !
      !      set bndry conditions
      !
            if(m.eq.ngrd) then
                 call bndry_outer(qrho,qpres,qpx,qpy,qpz, &
                    hrho,hpres,hpx,hpy,hpz, &
                    orho,opres,opx,opy,opz,epres, &
                    rmassq,rmassh,rmasso,bx,by,bz, &
                    bx0,by0,bz0,ngrd,nx,ny,nz, &
                    srho,rho_frac,o_conc,spress,spx,spy,spz, &
                    sbx_wind,sby_wind,sbz_wind,ti_te)
            else
             call bndry_flanks(qrho,qpres,qpx,qpy,qpz, &
                   hrho,hpres,hpx,hpy,hpz, &
                   orho,opres,opx,opy,opz, &
                   epres,bx,by,bz,     &
                   qrho,qpres,qpx,qpy,qpz,  &
                   hrho,hpres,hpx,hpy,hpz, &
                   orho,opres,opx,opy,opz, &
                   epres,bx,by,bz,    &
                   oldqrho,oldqpres,oldqpx,oldqpy,oldqpz,  &
                   oldhrho,oldhpres,oldhpx,oldhpy,oldhpz, &
                   oldorho,oldopres,oldopx,oldopy,oldopz, &
                   oldepres,oldbx,oldby,oldbz,  &
                   nx,ny,nz,ngrd,m,t_old,t_new,t_grd)
            endif 
      
      !
            if(m.eq.1)then
             call bndry_inner(qrho,qpres,qpx,qpy,qpz,rmassq, &
                   hrho,hpres,hpx,hpy,hpz,rmassh, &
                   orho,opres,opx,opy,opz,rmasso,epres,bx,by,bz, &
                   ngrd,nx,ny,nz,parm_srf,parm_mid,parm_zero, &
                   ijsrf,numsrf,ijmid,nummid,ijzero, &
                   numzero,erho,epress,alpha_e,ti_te,o_conc, &
                   sbx_wind,spress)
            endif
      
            if(spacecraft)then
               call set_imf(bx,by,bz,bx0,by0,bz0,bxp,byp,bzp, &
                         qrho,qpres,qpx,qpy,qpz,epres, &
                         rhop,svxp,svyp,svzp,svelx,spress,ti_te, &
                         rho_frac,nx,ny,nz,ngrd)
            endif
      !
      !      write(6,*)'fcsmooth speeds'
      !
            call set_rho(qrho,qpres,rmassq, &
                          hrho,hpres,rmassh,orho,opres,rmasso, &
                          epres,nx,ny,nz,ngrd,m,o_conc)
      
            call set_speed_agrd(qrho,qpres,qpx,qpy,qpz, &
                hrho,hpres,hpx,hpy,hpz,orho,opres,opx,opy,opz, &
                epres,bx,by,bz,bx0,by0,bz0,bsx,bsy,bsz,btot, &
                tvx,tvy,tvz,evx,evy, &
                rmassq,rmassh,rmasso,nx,ny,nz,ngrd,m, &
                pxmax,pymax,pzmax,pmax,csmax,alfmax,gamma, &
                vlim,alf_lim,o_conc,fastest) 
              t_step_new(m)=stepsz*xspac(m)/fastest
            write(6,*)'needed step of',t_step(m),t_step_new(m)
      !

	  
	  
	  
	  
	  !**************************************************************************
	  ! Xin Changes **************************************************************
	call bndry_inner_rot(rmassq,rmassh,rmasso, &
		 gamma1,eerg,ngrd,nx,ny,nz,parm_srf,parm_mid, &
		 ijsrf,numsrf,ijmid,nummid,erho,ti_te,o_conc)
! 
    call rigid_IM_rot(rmassq,rmassh,rmasso, &
		 gamma1,eerg,ngrd,nx,ny,nz,qrho,hrho,orho, &
		 erho,qpres,hpres,opres,epres,ti_te,o_conc,ut, &
		 utdec_start,re_equiv,bx0,by0,bz0)
!      !!!!!! test if parm_srf is changing (should be!)
       print *, '11','hrho_srf=', parm_srf(2,2), 'hpres=', parm_srf(5,2)
       print *, 'qrho_srf=', parm_srf(1,2), 'epres=', parm_srf(7,2)
	  
	  !*********************************************************************************
	  
	  
	  
	  

      !      divb correct
      !
      !      if(.not.start)then
      !        call divb_cor(bx,by,bz,qpx,qpy,qpz,qrho,qpres,
      !    +              ngrd,nx,ny,nz,m,srho)
      !      endif
      !
      !     sync time steps if needed and apply core conditions
      !
            if(m.lt.ngrd)then
              if(abs(t_new(m)-t_new(m+1)).le.0.001*t_step(m))then
                call bndry_corer(qrho,qpres,qpx,qpy,qpz, &
                   hrho,hpres,hpx,hpy,hpz,  &
                   orho,opres,opx,opy,opz,  &
                   epres,bx,by,bz,  &
                   nx,ny,nz,ngrd,m)
      !         write(6,*)'corer',m,t_new(m),m+1,t_new(m+1)
                t_old(m+1)=t_new(m+1)
              endif   
             endif
      !
            endif   !yes_step
      
         enddo   ! box sweep
      
         enddo   ! m_step loop
      
         ! 2015-08-17 write out global ionization, radial flux, and fluid particle trajectory
      906 FORMAT(F8.3,' ',E18.8)
      907 FORMAT(F8.3,3(' ',E18.8))
      908 FORMAT(F8.3,12(' ',E13.5))
      909 FORMAT(F8.3,12(',',F9.4))
         IF (t.GT.twrto) THEN
            electron_rate = 0.
            DO k = 1, nz
               DO j = 1, ny
                  DO i = 1, nx
                     electron_rate = electron_rate + mass_load_e(i,j,k)
                  END DO
               END DO
            END DO
            electron_rate = electron_rate*nchar*((lchar*xspac(1))**3.)/(rmassh*tchar)
            WRITE(502,906) ut, electron_rate
            ! non-dimensional number density multiplied by x&y velocity components
            hnvx = hpx(:,:,:,mwrite)/rmassh
            hnvy = hpy(:,:,:,mwrite)/rmassh
            hvx = hpx(:,:,:,mwrite)/hrho(:,:,:,mwrite) 
            hvy = hpy(:,:,:,mwrite)/hrho(:,:,:,mwrite)
            hden = hrho(:,:,:,mwrite)/rmassh
            CALL radial_flux(hnvx,hnvy,ndot, &
                 mwrite,nx,ny,nz,ngrd,vchar,nchar,re_equiv,rmassh)
            WRITE(503,908) ut, ndot
            !IF (fluid_track) THEN
            !   WRITE(504,909) ut, &
            !        fluidxyz(1,1), fluidxyz(2,1), fluidxyz(3,1), &
            !        fluidxyz(1,2), fluidxyz(2,2), fluidxyz(3,2), &
            !        fluidxyz(1,3), fluidxyz(2,3), fluidxyz(3,3), &
            !        fluidxyz(1,4), fluidxyz(2,4), fluidxyz(3,4)
            !END IF
            twrto = t + twrite
         END IF
      
      !     write(6,*)'time sync',t,t_new
            t_old(1)=t_new(1)
      !
      !    final sync on boundary conditions
      !
            do m=2,ngrd
              call bndry_corer(qrho,qpres,qpx,qpy,qpz, &
                   hrho,hpres,hpx,hpy,hpz, &
                   orho,opres,opx,opy,opz, &
                   epres,bx,by,bz, &
                   nx,ny,nz,ngrd,m-1)
      !
             call set_rho(qrho,qpres,rmassq, &
                          hrho,hpres,rmassh,orho,opres,rmasso, &
                          epres,nx,ny,nz,ngrd,m,o_conc)
            t_old(m)=t_old(m-1)
            t_new(m)=t_new(m-1)
      !
            enddo
      !
      !     write(6,999)t
      ! 999 format(' step 2 complete at t= ',1pe12.5)
      !
      ! 519  continue
            if(t.lt.ts1)goto 700
      
            ! 2015-02-105 A.R. generalizing the numbering of written-out fluid files,
            ! instead of this case-by-case bullshit.
            inchf = nchf - 10
            IF (inchf.LT.10) THEN
               WRITE(str1,'(I5)') int(t*t_equiv)
               OPEN(nchf,file=fluidstr//str1,status='new',form='unformatted')
            ELSE
               WRITE(str2,'(I5)') int(t*t_equiv)
               OPEN(nchf,file=fluidstr//str2,status='new',form='unformatted')
            END IF
      
            ! write restart data
            write(nchf)t
            write(nchf)qrho
            write(nchf)qpx
            write(nchf)qpy
            write(nchf)qpz
            write(nchf)qpres
            write(nchf)hrho
            write(nchf)hpx
            write(nchf)hpy
            write(nchf)hpz
            write(nchf)hpres
            write(nchf)orho
            write(nchf)opx
            write(nchf)opy
            write(nchf)opz
            write(nchf)opres
            write(nchf)bx
            write(nchf)by
            write(nchf)bz
            write(nchf)epres
            write(nchf)bx0
            write(nchf)by0
            write(nchf)bz0
            write(nchf)mass_load_photo
            write(nchf)mass_load_e
            write(nchf)photo_x
            write(nchf)photo_y
            write(nchf)photo_z
            write(nchf)e_impact_x
            write(nchf)e_impact_y
            write(nchf)e_impact_z
            write(nchf)elastic_x
            write(nchf)elastic_y
            write(nchf)elastic_z
            write(nchf)parm_srf,parm_mid,parm_zero, &
                 ijzero,numzero,ijmid,nummid,ijsrf,numsrf
            close(nchf)
      
            !IF (fluid_track) THEN
            !   IF (inchf.LT.10) THEN
            !      OPEN(251,file=partstr//str1,status='new',form='formatted')
            !   ELSE
            !      OPEN(251,file=partstr//str2,status='new',form='formatted')
            !   END IF
            !   DO i = 1, numpart
            !      WRITE(251,*) fluidxyz(1,i), fluidxyz(2,i), fluidxyz(3,i)
            !   END DO
            !   CLOSE(251)
            !END IF
      
            IF (inchf.GE.75) GOTO  707
            nchf = nchf + 1
            ts1=ts1+tsave
      
      700   if (t.lt.tmax) goto 1000
      
      707   WRITE(*,*) 'fluid_75 written, ending program'
      
            CLOSE(501)
            CLOSE(502)
            CLOSE(503)
            !CLOSE(504)
      
      END PROGRAM
      !*******************************************************************************
