      program multifluid
c
c     this is a 3-D modified three fluid simulation using the
c           electrons : arrays starting with e
c           solar wind : arrays starting with q
c           ionospheric: arrays starting with o oxygen, 
c                                         h hydrogen 
c
c      and change SET33(0.1,1.,0.1,1.0  to
c                 SET3(0.,1.,0.,1.
c     WARNING - MAKE SURE SPACE IS COMPATIBLE WITH GRAPHICS
c               routines
c               The arrays IJZERO,IJMID  and IJSRF have to
c               be modified in BSUBS.F if modified in MAIN
c     ADD in MIRROR DIPOLE  so bx is zero at wind boundary
c     GRAPHICS - CONTOUR AND FLOWS - HAVE TO BE manually set
c                for right aspect ratios
c
c      WARMING PLASMA and MAGNETIC FIELD data must be aligned in TIME
c
c       grid within grid size nt = 2*ngrd
c       ncts is the size of the data array for the IMF data file
c
       parameter (nx=111,ny=89,nz=77,nt=10,ngrd=5,ncraft=4,ncts=281)
c
c      graphics parameters:muvwp2=amax(mx,my,mz)+2,mz2=(mz-1)/2+1
       parameter (mx=56,my=45,mz=39,muvwp2=58,mz2=19)
c
      common /space/vvx(nx,ny,nz),vvy(nx,ny,nz),vvz(nx,ny,nz),
     +     tvx(nx,ny,nz),tvy(nx,ny,nz),tvz(nx,ny,nz),
     +     evx(nx,ny,nz),evy(nx,ny,nz),evz(nx,ny,nz)
c
      common /gridding/grd_xmin(9),grd_xmax(9),grd_ymin(9),grd_ymax(9),
     +           grd_zmin(9),grd_zmax(9),xspac(9),ncore(9),nbndry(9),
     +           rx,ry,rz,xdip,ydip,zdip,rearth,b0,
     +           sin_tilt,cos_tilt
      common /rotation/v_rot,r_rot,re_equiv
c     !!!!!! introduce tilt3 for inclination b/t spin axis and B axis
      common /uranus/tilt3,sin_tilt3,cos_tilt3,sin_spin,cos_spin,
     +           tilt4,sin_tilt4,cos_tilt4,spin,spinvectorx,
     +           spinvectory,spinvectorz
      !!!!!! sin_spin,cos_spin are rotating part of dipole field
c     !!!!!! x,y,zloc are for bcheck (check for B field in body code)
       real xloc,yloc,zloc
       real time_start,time_end
c
c      grid limits now set by grd_min grd_max arrays
c      ncore denotes couser grid to be hollowed out by fine grid
c      nbndry denotes finer grid to which coaser grid sets flanks
c      xspac is the relative grid spacing relative to inner grid system
c      xcraft is the actual position of the spacecraft in RE
c          4th dimension of the actual time
c      zcraft is the future position of the spacecraft in RE
c      rcraft is the position of the spacecraft for which
c           IMF is reference. NO alteration from boundary conditions applied
c
      real xcraft(4,ncraft),zcraft(4,ncraft),rcraft(3)
c
c     physics plasma quantities
c
      real bx(nx,ny,nz,ngrd),by(nx,ny,nz,ngrd),bz(nx,ny,nz,ngrd),
     +     qpx(nx,ny,nz,ngrd),qpy(nx,ny,nz,ngrd),qpz(nx,ny,nz,ngrd),
     +     qrho(nx,ny,nz,ngrd),qpres(nx,ny,nz,ngrd),
     +     opx(nx,ny,nz,ngrd),opy(nx,ny,nz,ngrd),opz(nx,ny,nz,ngrd),
     +     orho(nx,ny,nz,ngrd),opres(nx,ny,nz,ngrd),
     +     hpx(nx,ny,nz,ngrd),hpy(nx,ny,nz,ngrd),hpz(nx,ny,nz,ngrd),
     +     hrho(nx,ny,nz,ngrd),hpres(nx,ny,nz,ngrd),
     +     epres(nx,ny,nz,ngrd)
c
c    New magnetic field arrays
       real bx0a(nx,ny,nz,ngrd),by0a(nx,ny,nz,ngrd),bz0a(nx,ny,nz,ngrd),
     +      bx0b(nx,ny,nz,ngrd),by0b(nx,ny,nz,ngrd),bz0b(nx,ny,nz,ngrd),
     +      bx0a1(nx,ny,nz,ngrd),by0a1(nx,ny,nz,ngrd),
     +      bz0a1(nx,ny,nz,ngrd),
     +      bx0a2(nx,ny,nz,ngrd),by0a2(nx,ny,nz,ngrd),
     +      bz0a2(nx,ny,nz,ngrd),
     +      bx0a3(nx,ny,nz,ngrd),by0a3(nx,ny,nz,ngrd),
     +      bz0a3(nx,ny,nz,ngrd),
     +      bx0b1(nx,ny,nz,ngrd),by0b1(nx,ny,nz,ngrd),
     +      bz0b1(nx,ny,nz,ngrd),
     +      bx0b2(nx,ny,nz,ngrd),by0b2(nx,ny,nz,ngrd),
     +      bz0b2(nx,ny,nz,ngrd),
     +      bx0b3(nx,ny,nz,ngrd),by0b3(nx,ny,nz,ngrd),
     +      bz0b3(nx,ny,nz,ngrd)
c
c
c
c     work arrays for runge-kutta and smothing
c
      real oldbx(nx,ny,nz,ngrd),oldby(nx,ny,nz,ngrd),
     +     oldbz(nx,ny,nz,ngrd),oldqpx(nx,ny,nz,ngrd),
     +     oldqpy(nx,ny,nz,ngrd),oldqpz(nx,ny,nz,ngrd),
     +     oldqrho(nx,ny,nz,ngrd),oldqpres(nx,ny,nz,ngrd),
     +     oldopx(nx,ny,nz,ngrd),oldopy(nx,ny,nz,ngrd),
     +     oldopz(nx,ny,nz,ngrd),oldorho(nx,ny,nz,ngrd),
     +     oldopres(nx,ny,nz,ngrd),oldhpx(nx,ny,nz,ngrd),
     +     oldhpy(nx,ny,nz,ngrd),oldhpz(nx,ny,nz,ngrd),
     +     oldhrho(nx,ny,nz,ngrd),oldhpres(nx,ny,nz,ngrd),
     +     oldepres(nx,ny,nz,ngrd)
c
      real wrkbx(nx,ny,nz,ngrd),wrkby(nx,ny,nz,ngrd),
     +     wrkbz(nx,ny,nz,ngrd),wrkqpx(nx,ny,nz,ngrd),
     +     wrkqpy(nx,ny,nz,ngrd),wrkqpz(nx,ny,nz,ngrd),
     +     wrkqrho(nx,ny,nz,ngrd),wrkqpres(nx,ny,nz,ngrd),
     +     wrkopx(nx,ny,nz,ngrd),wrkopy(nx,ny,nz,ngrd),
     +     wrkopz(nx,ny,nz,ngrd),wrkorho(nx,ny,nz,ngrd),
     +     wrkopres(nx,ny,nz,ngrd),wrkhpx(nx,ny,nz,ngrd),
     +     wrkhpy(nx,ny,nz,ngrd),wrkhpz(nx,ny,nz,ngrd),
     +     wrkhrho(nx,ny,nz,ngrd),wrkhpres(nx,ny,nz,ngrd),
     +     wrkepres(nx,ny,nz,ngrd)
c
c     unperturbed quantities
c
      real bx0(nx,ny,nz,ngrd),by0(nx,ny,nz,ngrd),bz0(nx,ny,nz,ngrd)
     +     ,hrho0(nx,ny,nz,ngrd),orho0(nx,ny,nz,ngrd),
     +     qrho0(nx,ny,nz,ngrd),hpres0(nx,ny,nz,ngrd),
     +     opres0(nx,ny,nz,ngrd),qpres0(nx,ny,nz,ngrd),
     +     epres0(nx,ny,nz,ngrd)                   
c     !!!!!!
      real qrho0_temp(nx,ny,nz,ngrd),hrho0_temp(nx,ny,nz,ngrd),
     +     orho0_temp(nx,ny,nz,ngrd),qpres0_temp(nx,ny,nz,ngrd),
     +     hpres0_temp(nx,ny,nz,ngrd),opres0_temp(nx,ny,nz,ngrd),
     +     qpx0_temp(nx,ny,nz,ngrd),qpy0_temp(nx,ny,nz,ngrd),
     +     qpz0_temp(nx,ny,nz,ngrd),hpx0_temp(nx,ny,nz,ngrd),
     +     hpy0_temp(nx,ny,nz,ngrd),hpz0_temp(nx,ny,nz,ngrd),
     +     opx0_temp(nx,ny,nz,ngrd),opy0_temp(nx,ny,nz,ngrd),
     +     opz0_temp(nx,ny,nz,ngrd),epres0_temp(nx,ny,nz,ngrd)
c     !!!!!!
      real efldx(nx,ny,nz),efldy(nx,ny,nz),efldz(nx,ny,nz),
     +     curx(nx,ny,nz),cury(nx,ny,nz),curz(nx,ny,nz),
     +     bsx(nx,ny,nz),bsy(nx,ny,nz),bsz(nx,ny,nz),btot(nx,ny,nz),
     +     resistive(nx,ny,nz)
c
c     boundary condition arrays
c
      dimension bxf(ny,nz),byf(ny,nz),bzf(ny,nz),
     +         rhof(ny,nz),svxf(ny,nz),svyf(ny,nz),svzf(ny,nz)
      dimension bxp(ny,nz),byp(ny,nz),bzp(ny,nz),
     +         rhop(ny,nz),svxp(ny,nz),svyp(ny,nz),svzp(ny,nz)
      dimension future(ny,nz),past(ny,nz),
     +        bfld(ncts,4),rplas(ncts),svel(ncts,3)
      integer ncount(ny,nz)
c
      real tx(mx,my,mz),ty(mx,my,mz),tz(mx,my,mz),tg1(mx,my,mz),
     +     tg2(mx,my,mz2),tt(mx,my,mz),work(muvwp2,muvwp2),
     +     cross(my,mz),along(mx,mz)
c
      character*5 wd1,wd2,wd3,wd4
      character*8 label
      character*15 title
c
      integer ijsrf(3,15000),ijmid(3,15000),ijzero(3,120000)
      real    parm_srf(7,15000),parm_mid(7,15000)
      real reduct,rot_time,den_wind_time,utdec_start
      logical start,add_dip,ringo,update,save_dat,write_dat,
     +        spacecraft,tilting,warp,reload,maker,Adiffuse
c
c      ringo decides if you you want to plot 1 set of diagnostics
c              with no time stepping
c      update decides if time setting is to be reset
c
c     rho,pres,erg,px,py are the density, pressure, energy and momentum
c              in the x and y directions, respectively, of the fluid
c       the indices 1, 2 is required to store old and new values 
c
c     ijsrf give position of ionosphere in grid units - plasma
c           parameters stored in parm_srf
c     ijmid gives intermediate boundary - say representing the
c            atmosphere - plasma parameters stored in parm_mid
c     ijzero gives position of all grid units interior to surface
c
c     frho,ferg and fpx,fpy are the estimates of the fluid quantities
c           at n+1/2 as defined by the Lax-Wendroff scheme
c       the index of 3 is needed to store adjacent x values for these fns 
c
c     d_min is the minimum allowable density
c     stepsz is the size of the time step in terms of delta_x,y/(|umax|+c_s)
c
c     system dimensions are nx,ny,nz
c     variable grid spacing enabled with rxyz >1
c           rx,ry,rz should not be set identically to zero
c
      namelist/option/tmax,ntgraf,stepsz,start,tsave
      namelist/earth/xdip,ydip,zdip,rearth,
     +                tilt1,tilt2,tilting,rmassq,rmassh,rmasso
      namelist/speeds/cs_inner,alf_inner1,alf_inner2,
     +                alpha_e,den_earth,o_conc,gravity,
     +                ti_te,gamma,ringo,update,reload,Adiffuse,reduct,
     +                rot_time,utdec_start
      namelist/windy/re_wind,cs_wind,vx_wind1,vx_wind2,
     +              vy_wind1,vy_wind2,vz_wind1,vz_wind2,
     +              alfx_wind1,alfx_wind2,
     +              alfy_wind1,alfy_wind2,
     +              alfz_wind1,alfz_wind2,
     +              den_wind1,den_wind2,
     +             reynolds,resist,rho_frac,bfrac,vfrac,
     +             den_wind_time,den_wind_new,alfx_wind_new,
     +             alfy_wind_new,alfz_wind_new
      namelist/physical/re_equiv,b_equiv,v_equiv,rho_equiv,
     +              spacecraft,warp,uday,utstart
      namelist/smooth/chirho,chipxyz,chierg,
     +                difrho,difpxyz,diferg
c
c      open input data file
c
      open(5,file='fmpd3din',status='old',form='formatted')
c     open(6,file='rmpd3dout',status='unknown',form='formatted')
c     open(7,file='newdata',status='unknown',form='formatted')
      open(8,file='cur.dat',status='unknown',form='unformatted')
      open(9,file='pot.dat',status='unknown',form='unformatted')
c     !!!!!!
c      open(50,file='bcheck0.dat',status='replace',form='formatted')
c      open(51,file='bcheck1.dat',status='replace',form='formatted')
c      open(52,file='bcheck2.dat',status='replace',form='formatted')
c      open(53,file='bcheck3.dat',status='replace',form='formatted')
c      open(54,file='bcheck4.dat',status='replace',form='formatted')
c      open(55,file='bcheck5.dat',status='replace',form='formatted')
c      open(56,file='bcheck6.dat',status='replace',form='formatted')
c      open(57,file='bcheck7.dat',status='replace',form='formatted')
c      open ncargraphics
c
c 
      call opngks
      call gsclip(0)
c
c     set color table 
c
      call cpclrs
c
      call gselnt(0)
      call gsplci(1)
      call wtstr(.4,.975,'3D MUTANT CODE',2,0,0)
c
c
c     read input parameters
c
      read(5,option)
      read(5,earth)
      read(5,speeds)
      read(5,windy)
      read(5,physical)
      read(5,smooth)
c
c     output to test whether the parameters are the actual ones you want
c
      write(6,option)
      write(6,earth)
      write(6,speeds)
      write(6,windy)
      write(6,physical)
      write(6,smooth)
c
      do i=1,ngrd
       read(5,*)grd_xmin(i),grd_xmax(i),grd_ymin(i),grd_ymax(i),
     +      grd_zmin(i),grd_zmax(i),xspac(i),ncore(i),nbndry(i)
       write(6,*)grd_xmin(i),grd_xmax(i),grd_ymin(i),grd_ymax(i),
     +      grd_zmin(i),grd_zmax(i),xspac(i),ncore(i),nbndry(i)
       ix=1+(grd_xmax(i)-grd_xmin(i))/xspac(i)
       iy=1+(grd_ymax(i)-grd_ymin(i))/xspac(i)
       iz=1+(grd_zmax(i)-grd_zmin(i))/xspac(i)
       if((ix.ne.nx).or.(iy.ne.ny).or.(iz.ne.nz))then
         write(6,*)' WARNING: SIZES',ix,iy,iz,nx,ny,nz
         stop
       endif
      enddo
c
c     ncore -> grid gets all the information from this grid no.
c     nbdry -> which grid data will be applied to this grid no.for bndry
c
c     write important data to graphics file
c
      write(wd4,'(f5.3)')stepsz
c
      title='stepsz = '//wd4
      call wtstr(.75,.85,title,1,0,0)
c
      write(wd1,'(f5.3)')xdip
      write(wd2,'(f5.3)')ydip
      write(wd3,'(f5.3)')zdip
c
      title='xdip = '//wd1
      call wtstr(.15,.82,title,1,0,0)
      title='ydip = '//wd2
      call wtstr(.35,.82,title,1,0,0)
      title='zdip = '//wd3
      call wtstr(.55,.82,title,1,0,0)
c

      write(wd1,'(f5.1)')rearth
c
      title='rearth = '//wd1
      call wtstr(.15,.79,title,1,0,0)
c
c     calculate effective magnetic field strength
c
      erho=den_earth*rmassh
      b01=alf_inner1*sqrt(erho)*rearth**3
      b02=alf_inner2*sqrt(erho)*rearth**3
      alf_lim=3.00*alf_inner1
      b0=b01
      bmax=b02
      delb0=(b02-b01)/tmax
c
      write(wd1,'(f5.3)')cs_inner
      write(wd2,'(f5.3)')alf_inner1
      write(wd3,'(f5.3)')alf_inner2
      write(wd4,'(f5.1)')den_earth
c
      title='cs_inner = '//wd1
      call wtstr(.15,.76,title,1,0,0)
      title='alf_inner1= '//wd2
      call wtstr(.35,.76,title,1,0,0)
      title='alf_inner2= '//wd3
      call wtstr(.55,.76,title,1,0,0)
      title='den_earth = '//wd4
      call wtstr(.75,.76,title,1,0,0)
c

      write(wd1,'(f5.3)')o_conc
      write(wd2,'(f5.3)')gravity
      write(wd3,'(f5.3)')rmasso
c
      title='o_conc = '//wd1
      call wtstr(.15,.73,title,1,0,0)
      title='gravity= '//wd2
      call wtstr(.35,.73,title,1,0,0)
      title='rmasso= '//wd2
      call wtstr(.35,.73,title,1,0,0)
c
      write(wd1,'(f5.3)')rho_wind1
      write(wd2,'(f5.3)')rho_wind2
      write(wd3,'(f5.3)')vx_wind1
      write(wd4,'(f5.3)')vx_wind2
c
      title='rho_wind1= '//wd1
      call wtstr(.15,.7,title,1,0,0)
      title='rho_wind2= '//wd2
      call wtstr(.35,.7,title,1,0,0)
      title='vx_wind1 = '//wd3
      call wtstr(.55,.7,title,1,0,0)
      title='vx_wind2 = '//wd4
      call wtstr(.75,.7,title,1,0,0)

c
      write(wd1,'(f5.3)')vy_wind1
      write(wd2,'(f5.3)')vy_wind2
      write(wd3,'(f5.3)')vz_wind1
      write(wd4,'(f5.3)')vz_wind2
c
      title='vy_wind1= '//wd1
      call wtstr(.15,.67,title,1,0,0)
      title='vy_wind2= '//wd2
      call wtstr(.35,.67,title,1,0,0)
      title='vz_wind1 = '//wd3
      call wtstr(.55,.67,title,1,0,0)
      title='vz_wind2 = '//wd4
      call wtstr(.75,.67,title,1,0,0)
c
      write(wd1,'(f5.3)')alfx_wind1
      write(wd2,'(f5.3)')alfx_wind2
      write(wd3,'(f5.3)')alfy_wind1
      write(wd4,'(f5.3)')alfy_wind2

c
      title='alfx1 = '//wd1
      call wtstr(.15,.64,title,1,0,0)
      title='alfx2 = '//wd2
      call wtstr(.35,.64,title,1,0,0)
      title='alfy1 = '//wd3
      call wtstr(.55,.64,title,1,0,0)
      title='alfy2 = '//wd4
      call wtstr(.75,.64,title,1,0,0)
c
      write(wd3,'(f5.3)')alfz_wind1
      write(wd4,'(f5.3)')alfz_wind2
      title='alfz1 = '//wd3
      call wtstr(.55,.61,title,1,0,0)
      title='alfz2 = '//wd4
      call wtstr(.75,.61,title,1,0,0)

      write(wd1,'(f5.1)')re_wind
      write(wd2,'(f5.0)')reynolds
      write(wd3,'(f5.0)')resist
c     re_wind sets raduius from earth where initial wind placed
c     reynolds coefficient for surface currents
c     resist equivalent if you wish to run anomalous resistivity
c     bfrac determines the percentage of the tangential magnetic
c     field allowed at earth's surface
c
      title='re_wind = '//wd1
      call wtstr(.15,.58,title,1,0,0)
      title='reynolds = '//wd2
      call wtstr(.35,.58,title,1,0,0)
      title='resist = '//wd3
      call wtstr(.55,.58,title,1,0,0)
c

      write(wd1,'(f5.3)')bfrac
      write(wd2,'(f5.3)')vfrac
      title='bfrac = '//wd1
      call wtstr(.35,.55,title,1,0,0)
      title='vfrac = '//wd2
      call wtstr(.55,.55,title,1,0,0)
c
      write(wd1,'(f5.3)')chirho
      write(wd2,'(f5.3)')chipxyz
      write(wd3,'(f5.3)')chierg
c
      title='chirho = '//wd1
      call wtstr(.15,.52,title,1,0,0)
      title='chipxyz = '//wd2
      call wtstr(.35,.52,title,1,0,0)
      title='chierg = '//wd3
      call wtstr(.55,.52,title,1,0,0)
c
      write(wd1,'(f5.3)')difrho
      write(wd2,'(f5.3)')difpxyz
      write(wd3,'(f5.3)')diferg
c
      title='chirho = '//wd1
      call wtstr(.15,.49,title,1,0,0)
      title='chipxyz = '//wd2
      call wtstr(.35,.49,title,1,0,0)
      title='chierg = '//wd3
      call wtstr(.55,.49,title,1,0,0)
c
      write(wd1,'(f4.1)')tilt1
      title='tilt1 = '//wd1
      call wtstr(.15,.40,title,1,0,0)
      write(wd1,'(f4.1)')tilt2
      title='tilt2 = '//wd1
      call wtstr(.30,.40,title,1,0,0)
      tilt=tilt1
      sin_tilt=sin(tilt*.0174533)
      cos_tilt=cos(tilt*.0174533)
c       !!!!!!
        spinvectorx=-cos((abs(tilt)-7.9)*.0174533)
        spinvectory=0.0
        spinvectorz=sin((abs(tilt)-7.9)*.0174533)
c       !!!!!!
      dtilt=(tilt2-tilt1)/tmax
c       !!!!!! needing add tilt3 into input file later
       print *, 'tilt1', tilt1
        spin=0.0
        dspin=360/tmax
        tilt3=0.0
c         -97.9
        sin_tilt3=sin(tilt3*.0174533)
        cos_tilt3=cos(tilt3*.0174533)
c       !!!!!! tilt3 is the tilted angel almost towards Sun from z axis
c       !!!!!! tilt4 = tilt+tilt3 for modified subroutine dipole
c        tilt4=tilt+tilt3
         tilt4=tilt
        sin_tilt4=sin(tilt4*.0174533)
        cos_tilt4=cos(tilt4*.0174533)
c       !!!!!!x,y,zloc are for bcheck, in grid points,5 gps=1 Ru  
       xloc=71
       yloc=60
       zloc=54
c
c     jupiter parameters
c
c     planet_rad=71000.   !km
c     planet_per=9.7     !hr
c     lunar_rad=5.9      !orbital radii
c     v_rot=6.2832*planet_rad/(planet_per*3600.)/v_equiv  ! normalized units
c     r_rot=40.0   !Re where corotation stops
c
c     Earth parameters
c     planet_rad=6371.   !km
c     planet_per=24.     !hr
c     Ganymede
c      planet_rad=2600.   !km
c      planet_per=2400.     !hr
c      lunar_rad=10.      !orbital radii
c      v_rot=6.2832*planet_rad/(planet_per*3600.)/v_equiv  ! normalized units
c      r_rot=10.0   !Re where corotation stops
c     !!!!!! Uranus
      planet_rad=25484.   !km
      planet_per=17.24     !hr
c      lunar_rad=????0.      !orbital radii
      v_rot=6.2832*planet_rad/(planet_per*3600.)/v_equiv  ! normalized units
      r_rot=6.0   !Re where corotation stops
c     !!!!!!
c     spacecraft stuff
c
c     gravity in m/s**2 at earths surface !need t_equiv in normalization
      t=0.
      t_equiv=planet_rad*re_equiv/v_equiv
      grav=gravity*(planet_rad*re_equiv*1000.)/(1000.*v_equiv)**2
      ut=utstart
      d_min=0.001
c
c
c     Initial position of spacecraft in RE but simulation directions
c        WIND :
       xcraft(1,1)=-2.5
       xcraft(2,1)=1.00
       xcraft(3,1)=-0.00
c
c      reference craft
c 
       rcraft(1)=xcraft(1,1)
       rcraft(2)=xcraft(2,1)
       rcraft(3)=xcraft(3,1)
c
c      POLAR
c
       xcraft(1,2)=-1.5
       xcraft(2,2)=-.5
       xcraft(3,2)=0.36
c
c      IMP : 
c
       xcraft(1,3)=12.5
       xcraft(2,3)=-2.5
       xcraft(3,3)=0.
c
c      GEOTAIL
c
       xcraft(1,4)=20.25
       xcraft(2,4)=-0.5
       xcraft(3,4)=-2.
c
       zut=ut-.01
       rut=zut
       vut=rut
       xcraft(4,1)=zut
       xcraft(4,2)=zut
       xcraft(4,3)=zut
       call limcraft(xcraft,ncraft,re_equiv,ngrd)
c
      do n=1,ncraft
       do i=1,4
         zcraft(i,n)=xcraft(i,n)
       enddo
      enddo
c
        open(31,file='wind.dat',status='unknown',form='formatted') 
        open(32,file='polar.dat',status='unknown',form='formatted')
        open(33,file='imp8.dat',status='unknown',form='formatted')
        open(34,file='geotail.dat',status='unknown',form='formatted')
      if(spacecraft) then
        open(21,file='wind.pos',status='unknown',form='formatted')
c       open(22,file='polar.pos',status='unknown',form='formatted')
c       open(23,file='equators.pos',status=unknown'',form='formatted')
c       open(24,file='geotail.pos',status='unknown',form='formatted')
c       open(27,file='wind.den',status='unknown',form='formatted')
c       open(28,file='wind.vel',status='unknown',form='formatted')
        open(27,file='wind.plas',status='unknown',form='formatted')
        open(29,file='wind.mag',status='unknown',form='formatted')
      endif
c
c     the magnetic field in dimensionless units is actually in Alfven speeds
c             relative to the normalizing velocity
c     the temperature is in sound speeds relative to the normalizing speed
c             all squared
c    
c     the magnetospheric plasma density are assumed to vary as
c             rho proportional to (R)**-alpha_e
c             Temperatue proportional to (R)**alpha_e
c         with total pressure constant which is needed for equilibrium
c
c     now find the equivalent pressure of magnetosphere for the given
c         sound speed
c
      gamma1=gamma-1.
      epress=cs_inner**2*erho/gamma
      eerg=epress/gamma1
c
c     do the same for the solar wind
c     !!!!!!den_wind changes at den_wind_time at unit of rot period
      if(ut.ge.den_wind_time*planet_per)then
      den_wind1=den_wind_new
      den_wind2=den_wind_new
      alfx_wind1=alfx_wind_new
      alfy_wind1=alfy_wind_new
      alfz_wind1=alfz_wind_new
      alfx_wind2=alfx_wind_new
      alfy_wind2=alfy_wind_new
      alfz_wind2=alfz_wind_new      
      endif
      rho_wind1=den_wind1*rmassq
      rho_wind2=den_wind2*rmassq
      srho=rho_wind1
      delrho=(rho_wind2-rho_wind1)/tmax
      spress=(cs_wind**2*srho/gamma)/gamma1
      svelx=vx_wind1
      svely=vy_wind1
      svelz=vz_wind1
c
      spx=srho*svelx
      spy=srho*svely
      spz=srho*svelz
      serg=0.5*(svelx**2+svely**2+svelz**2)*srho+spress
      delvx_wind=(vx_wind2-vx_wind1)/tmax
      delvy_wind=(vy_wind2-vy_wind1)/tmax
      delvz_wind=(vz_wind2-vz_wind1)/tmax
c
      delbx_wind=(alfx_wind2*sqrt(rho_wind2)
     +           -alfx_wind1*sqrt(rho_wind1))/tmax
      delby_wind=(alfy_wind2*sqrt(rho_wind2)
     +           -alfy_wind1*sqrt(rho_wind1))/tmax
      delbz_wind=(alfz_wind2*sqrt(rho_wind2)
     +           -alfz_wind1*sqrt(rho_wind1))/tmax
      sbx_wind=alfx_wind1*sqrt(rho_wind1)
      sby_wind=alfy_wind1*sqrt(rho_wind1)
      sbz_wind=alfz_wind1*sqrt(rho_wind1)
c
      deltg=tmax/float(ntgraf)
      delt=stepsz
      tgraf=0.
      ts1=tsave
      twrt=0.
      dtwrt=1.
      tfresh=-0.1
      refresh=150.
      tdiv=11.
      write_dat=.true.
c
c     ************************************************
c     check for restart
c     ************************************************
c
      nchf=11
      if(start) go to 80
      open(11,file='fluid11',status='unknown',form='unformatted')
      open(12,file='fluid12',status='unknown',form='unformatted')
       read(nchf) t
       rewind nchf
c
       inchf=23-nchf
       read(inchf) t2
       rewind inchf
       if(t.lt.t2) then
         nchf=inchf
         inchf=23-nchf
       endif
c
c     read restart data
c
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
       if(Adiffuse)then
       read(nchf)qrho0
       read(nchf)hrho0
       read(nchf)orho0
       read(nchf)qpres0
       read(nchf)hpres0
       read(nchf)opres0
       read(nchf)epres0   
       endif
c
       read(nchf)parm_srf,parm_mid,
     +          ijzero,numzero,ijmid,nummid,ijsrf,numsrf
       close(11)
       close(12)
c
c     !!!!!! for B field check (for bcheck0.dat)
c      do m=1,ngrd
c      write(50,60)bx0(xloc,yloc,zloc,m),by0(xloc,yloc,zloc,m),
c     +bz0(xloc,yloc,zloc,m),bx(xloc,yloc,zloc,m),by(xloc,yloc,zloc,m),
c     +bz(xloc,yloc,zloc,m),ut,delt,t_equiv,m
c60    format(10(1x,1pe9.2))
c
c      enddo
c     !!!!!!
c     !!!!!!
      if(update)t=0.
      ts1=t+tsave
      tmax=t+tmax
      tgraf=t+deltg
      twrt=twrt+t
      tfresh=tfresh+t
      nchf=11
      if(Adiffuse)then
      if(reload)then
        do m=1,ngrd
        do k=1,nz
        do j=1,ny
        do i=1,nx
         qrho0(i,j,k,m)=0.925*qrho(i,j,k,m)
         hrho0(i,j,k,m)=0.925*hrho(i,j,k,m)
         orho0(i,j,k,m)=0.925*orho(i,j,k,m)
         qpres0(i,j,k,m)=0.925*qpres(i,j,k,m)
         hpres0(i,j,k,m)=0.925*hpres(i,j,k,m)
         opres0(i,j,k,m)=0.925*opres(i,j,k,m)
         epres0(i,j,k,m)=0.925*epres(i,j,k,m)
        enddo
        enddo
        enddo
        enddo 
      endif
c
c     rescale oxygen density at inner boundary by oxygen scale
c
      print *, 'Adiffuse=',Adiffuse,'reduct=',reduct
      if(ringo)then
        do nn=1,nummid
          parm_mid(1,nn)=reduct*parm_mid(1,nn)
          parm_mid(2,nn)=reduct*parm_mid(2,nn)
          parm_mid(3,nn)=reduct*parm_mid(3,nn)
          parm_mid(4,nn)=reduct*parm_mid(4,nn)
          parm_mid(5,nn)=reduct*parm_mid(5,nn)
          parm_mid(6,nn)=reduct*parm_mid(6,nn)
          parm_mid(7,nn)=reduct*parm_mid(7,nn)
        enddo
c
        do nn=1,numsrf
          parm_srf(1,nn)=reduct*parm_srf(1,nn)
          parm_srf(2,nn)=reduct*parm_srf(2,nn)
          parm_srf(3,nn)=reduct*parm_srf(3,nn)
          parm_srf(4,nn)=reduct*parm_srf(4,nn)
          parm_srf(5,nn)=reduct*parm_srf(5,nn)
          parm_srf(6,nn)=reduct*parm_srf(6,nn)
          parm_srf(7,nn)=reduct*parm_srf(7,nn)
        enddo
      endif
      endif
c
c     initialize plasma resistivity
c       !!!!!!rotating part for dipole fied around z axis
        if(ut.ge.rot_time*planet_per)then
        sin_spin=sin(2.*3.1415926*ut/planet_per)
        cos_spin=cos(2.*3.1415926*ut/planet_per)
        else
        sin_spin=sin(0.0)
        cos_spin=cos(0.0)
        endif
c
      m=1
      rx=xspac(m)
      ry=xspac(m)
      rz=xspac(m)
      call set_resist(resistive,nx,ny,nz,resist,ijzero,numzero,
     +             ijmid,nummid,ijsrf,numsrf)
c
      write(6,79) nchf
   79 format('  restart from   mpd3d',i2)
      rewind nchf
      nchf=11
      goto 170

c
c     ******************************
c            initialization
c     ******************************
        call cpu_time(time_start)
c
c       !!!!!!rotating part for dipole fied around z axis
        if(ut.ge.rot_time*planet_per)then
        sin_spin=sin(2.*3.1415926*ut/planet_per)
        cos_spin=cos(2.*3.1415926*ut/planet_per)
        else
        sin_spin=sin(0.0)
        cos_spin=cos(0.0)
        endif
c       !!!!!!for start=.f. in fmp3din file
c
c     ambient plasma
c
c      initialize indices and surface points
c
  80   numsrf=0
       nummid=0
       numzero=0
c
c      scale lengths for plasma sheet population
c
        sheet_den=0.25
        alpha_s=4.
c
      do 130 m=1,ngrd
        dx=(grd_xmax(m)-grd_xmin(m))/(nx-1.)
        dy=(grd_ymax(m)-grd_ymin(m))/(ny-1.)
        dz=(grd_zmax(m)-grd_zmin(m))/(nz-1.)
c
c      create dipole magnetic field and load magnetospheric plasma 
c
      do 130 k=1,nz
       az=(grd_zmin(m)+dz*(k-1)-zdip)
c
       do 120 j=1,ny
         ay=grd_ymin(m)+dy*(j-1)-ydip
c
        do 110 i=1,nx
         ax=(grd_xmin(m)+dx*(i-1)-xdip)
         xp=ax*cos_tilt-az*sin_tilt 
         zp=ax*sin_tilt+az*cos_tilt
         ar=sqrt(xp**2+ay**2+zp**2)
c
c        !!this is rotation for pressure such opres or pres-vel
c        but it no this, it is OK if rotate rvx,y,z
c        if operating both of two parts above, effects canceled
c         xp01=zp*sin_tilt3+xp*cos_tilt3
c         yp01=ay
c         zp01=zp*cos_tilt3-1.*xp*sin_tilt3
c         xp=xp01
c         ay=yp01
c         zp=zp01
c        !!
c        !!!!!!I can also add spinvector rotation here!!
c
c        determine magnetic dipole field
c        !!!!!!
c        !!!!!! tilt4 = tilt+tilt3 for modified subroutine dipole
c         tilt4=tilt+tilt3
         tilt4=tilt
         sin_tilt4=sin(tilt4*.0174533)
         cos_tilt4=cos(tilt4*.0174533)

c        !!!!!!
         call dipole(bx0(i,j,k,m),by0(i,j,k,m),
     +              bz0(i,j,k,m),ax,ay,az)
c
c         zero interior magnetic field so alfven speed small
c
        if (ar.lt.rearth-1.5)then
          bx0(i,j,k,m)=0.
          by0(i,j,k,m)=0.
          bz0(i,j,k,m)=0.
        endif
c
c        !!!!!!make dipole field rotating round z axis
c              (bx,y,z0a bx,y,zb not defined yet)
c         bx0a(i,j,k,m)=bx0(i,j,k,m)*cos_spin
c     +                -by0(i,j,k,m)*sin_spin
c         by0a(i,j,k,m)=bx0(i,j,k,m)*sin_spin
c     +                +by0(i,j,k,m)*cos_spin
c         bz0a(i,j,k,m)=bz0(i,j,k,m)
c
c          bx0a(i,j,k,m)=bx0(i,j,k,m)
c          by0a(i,j,k,m)=by0(i,j,k,m)
c          bz0a(i,j,k,m)=bz0(i,j,k,m)
c        !!!!!!make rotational dipole field tilted around y axis
c              (right-hand coordinates)
c         bx0b(i,j,k,m)=bz0a(i,j,k,m)*sin_tilt3
c     +                +bx0a(i,j,k,m)*cos_tilt3
c         by0b(i,j,k,m)=by0a(i,j,k,m)
c         bz0b(i,j,k,m)=bz0a(i,j,k,m)*cos_tilt3
c     +                -1.*bx0a(i,j,k,m)*sin_tilt3
c        !!!!!!bx,y,z0a is the after-rotation field
c              bx,y,z0b is the field after-tilted towards the Sun
c         bx0(i,j,k,m)=bx0b(i,j,k,m)
c         by0(i,j,k,m)=by0b(i,j,k,m)
c         bz0(i,j,k,m)=bz0b(i,j,k,m)
c        !!!!!!
c
c     !!!!!! for B field check (for bcheck1.dat)
c      do m=1,ngrd
c      write(51,61)bx0(xloc,yloc,zloc,m),by0(xloc,yloc,zloc,m),
c     +bz0(xloc,yloc,zloc,m),bx(xloc,yloc,zloc,m),by(xloc,yloc,zloc,m),
c     +bz(xloc,yloc,zloc,m),ut,delt,t_equiv,m
c61    format(10(1x,1pe9.2))
c
c      enddo
c     !!!!!!

c        !!
c         xp01=zp*sin_tilt3+xp*cos_tilt3
c         yp01=ay
c         zp01=zp*cos_tilt3-1.*xp*sin_tilt3
c         xp=xp01
c         ay=yp01
c         zp=zp01
c        !!

c        set up rotational properties
c
c         rx=ax*re_equiv
c         ry=ay*re_equiv
c         rd=sqrt(rx**2+ry**2)
         !!!!!!10/27/2013
c         rxp=rx*cos_tilt3-rz*sin_tilt3
c         ryp=ry
c         rzp=rx*sin_tilt3+rz*cos_tilt3
c         rx=rxp
c         ry=ryp
c         rz=rzp
         !!!!!!
c         if(rd.lt.r_rot)then
c           vfrac=1.
c         else
c           vfrac=((2.*r_rot**2)/(r_rot**2+ rd**2))**2
c         endif
c         rvy=rx*v_rot*vfrac
c         rvx=-ry*v_rot*vfrac
         !!!!!!10/27/2013
c         rvz=
c        !!!!!!make the spin axis tilted toward its right direction
c        not introduce rvx00, rvy00, rvz00 yet, not tilt3 yet(input)!!
c         rvx00=rvx*cos_tilt3
c         rvy00=rvy
c         rvz00=-1.*rvx*sin_tilt3
c        !!!!!!turn to their original expression,rvz=0 so ignored
c         rvx=rvx00
c         rvy=rvy00
c         rvz=rvz00
c
         call vfldrot(rvx,rvy,rvz,ax,ay,az,vfrac,v_rot,r_rot)
c        !!!!!!
c
c        for Jupiter top hat distribution
c        ar_iono=sqrt(xp**2+ay**2+0.75*zp**2)
c
c        for rearth = 6
c        ar_iono=sqrt(xp**2+ay**2+2.*zp**2)
c
c        for rearth = 10
c        ar_iono=sqrt(xp**2+ay**2+1.25*zp**2)
c        isotropic
         ar_iono=sqrt(xp**2+ay**2+zp**2)
c        ar_iono=amax1(0.0001,ar_iono)
         ar_iono=amax1(1.01*rearth,ar_iono)
c        ra=((ar_iono+0.5*rearth)/(1.5*rearth))**(-alpha_e)
c        !!!!!! 0.0994 means the decrease rate of density to 2.4, L=5
         ra=exp(-(ar_iono-rearth)/(0.0994*rearth))
c        !!!!!!NEED to change 2 placs at subroutine "bndry_inner_rot"
         zheight=amax1(1.,(zp**2+(1.5*rearth)**2)/(3.0*rearth)**2)
         ra=ra/zheight**2
         rho_iono=amax1(erho*ra,0.001)
c  
         r_equat=(ar**3+0.001)/(xp**2+ay**2+0.001)
         r_equat=amax1(r_equat,rearth)
         erg_sphere=eerg*(0.001+(rearth/r_equat)**4)
c        r_equat=amin1(r_equat,grd_xmax)
c        erg_sphere=eerg*r_equat/grd_xmax
c
         arho=amax1(rho_iono,d_min)
         aerg=amin1(abs(erg_sphere),0.01*arho)
c         
           epres(i,j,k,m)=0.5*gamma1*aerg/ti_te
c          
           qrho(i,j,k,m)=rmassq*d_min
           qpres(i,j,k,m)=0.
           qpx(i,j,k,m)=qrho(i,j,k,m)*rvx
           qpy(i,j,k,m)=qrho(i,j,k,m)*rvy
           qpz(i,j,k,m)=qrho(i,j,k,m)*rvz
c           qpz(i,j,k,m)=0.
c          !!!!!! instead by rvz term above       
           hrho(i,j,k,m)=arho*rmassh
           hpres(i,j,k,m)=0.5*gamma1*aerg
           hpx(i,j,k,m)=hrho(i,j,k,m)*rvx
           hpy(i,j,k,m)=hrho(i,j,k,m)*rvy
           hpz(i,j,k,m)=hrho(i,j,k,m)*rvz
c           hpz(i,j,k,m)=0.
c          !!!!!! instead by rvz term above    
           orho(i,j,k,m)=rho_iono*rmasso*o_conc
           opres(i,j,k,m)=0.5*gamma1*aerg*o_conc
           opx(i,j,k,m)=orho(i,j,k,m)*rvx
           opy(i,j,k,m)=orho(i,j,k,m)*rvy
           opz(i,j,k,m)=orho(i,j,k,m)*rvz
c           opz(i,j,k,m)=0.
c          !!!!!! instead by rvz term above    
c
c          check for Io
c  
c          ar_io=sqrt(xp**2+ay**2)
c          ar_io=ar_io*re_equiv
c          if((ar_io.lt.1.5*lunar_rad).and.(ar.gt.0.5*lunar_rad)
c    +          .and.(abs(zp)*re_equiv.lt.0.75*lunar_rad))then
c            ar_io=(ar_io-lunar_rad)/re_equiv    !in grid units
c            dr=amin1(1.,(2./(1.+(ar_io)**2 +zp**2))**2)
c            arho=den_lunar*dr*rmasso
c            orho(i,j,k,m)=orho(i,j,k,m) +arho
c            opres(i,j,k,m)=opres(i,j,k,m)
c    +            +0.5*arho*(rvx**2+rvy**2)*dr/gamma/rmasso
c            opx(i,j,k,m)=opx(i,j,k,m)+arho*rvx
c            opy(i,j,k,m)=opy(i,j,k,m)+arho*rvy
c          endif
c
           bx(i,j,k,m)=0.
           by(i,j,k,m)=0.
           bz(i,j,k,m)=0.
c
c        test for boundary point of planets surface or
c        interior to planet
c
        if((ar.gt.rearth+.6).or.(m.gt.1))goto 110
c
        if(ar.lt.rearth-1.5) then
          numzero=numzero+1
          ijzero(1,numzero)=i
          ijzero(2,numzero)=j
          ijzero(3,numzero)=k
c         !!!!!!
c          qrho(i,j,k,m)=d_min*erho*rmassq
c          qrho(i,j,k,m)=200.0*d_min*rmassq
         qrho(i,j,k,m)=d_min*rmassq
          qpres(i,j,k,m)=0.5*gamma1*eerg*d_min
          hrho(i,j,k,m)=erho*rmassh
          hpres(i,j,k,m)=0.5*gamma1*eerg
          orho(i,j,k,m)=erho*o_conc*rmasso
          opres(i,j,k,m)=0.5*gamma1*eerg*o_conc
          epres(i,j,k,m)=0.5*gamma1*eerg/ti_te
          qpx(i,j,k,m)=qrho(i,j,k,m)*rvx
          qpy(i,j,k,m)=qrho(i,j,k,m)*rvy
          qpz(i,j,k,m)=qrho(i,j,k,m)*rvz
c          qpz(i,j,k,m)=0.
c         !!!!!! instead by rvz term above   
          hpx(i,j,k,m)=hrho(i,j,k,m)*rvx
          hpy(i,j,k,m)=hrho(i,j,k,m)*rvy
          hpz(i,j,k,m)=hrho(i,j,k,m)*rvz
c          hpz(i,j,k,m)=0.
c         !!!!!! instead by rvz term above 
          opx(i,j,k,m)=orho(i,j,k,m)*rvx
          opy(i,j,k,m)=orho(i,j,k,m)*rvy
          opz(i,j,k,m)=orho(i,j,k,m)*rvz
c          opz(i,j,k,m)=0.
c         !!!!!! instead by rvz term above 
c        
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
c
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
c
  110   continue
  120  continue
  130 continue
c
      write(6,132)numsrf,nummid,numzero
  132 format(' interior and zero points',3(1x,i6))
c
c     initialize solar wind plasa can be placed beyond
c       the earth at a radius of re_wind
c
c      !!!!!!
c      wind_bnd=2.5/re_equiv
       wind_bnd=re_wind/re_equiv
c      !!!!!!
      ramp=-wind_bnd-grd_xmin(ngrd)
      ofrac=rho_frac*o_conc
      do m=1,ngrd
c
      dx=(grd_xmax(m)-grd_xmin(m))/(nx-1.)
      dy=(grd_ymax(m)-grd_ymin(m))/(ny-1.)
      dz=(grd_zmax(m)-grd_zmin(m))/(nz-1.)
c
      do k=1,nz
        az=(grd_zmin(m)+dz*(k-1)-zdip)
c
      do j=1,ny
        ay=grd_ymin(m)+dy*(j-1)-ydip
        do i=1,nx
         ax=(grd_xmin(m)+dx*(i-1)-xdip)

         ar=sqrt(ax**2+ay**2+az**2)
         if((ax.lt.-wind_bnd))then
          qrho(i,j,k,m)=srho
          qpx(i,j,k,m)=spx
          qpy(i,j,k,m)=spy
          qpz(i,j,k,m)=spz
c
          apres=cs_wind**2*qrho(i,j,k,m)/gamma
          qpres(i,j,k,m)=0.5*apres
          epres(i,j,k,m)=0.5*apres/ti_te
c
          hrho(i,j,k,m)=srho*rho_frac*rmassh/rmassq
          hpx(i,j,k,m)=spx*rho_frac*rmassh/rmassq
          hpy(i,j,k,m)=spy*rho_frac*rmassh/rmassq
          hpz(i,j,k,m)=spz*rho_frac*rmassh/rmassq
          hpres(i,j,k,m)=qpres(i,j,k,m)*rho_frac
c
          orho(i,j,k,m)=srho*ofrac*rmasso/rmassq
          opx(i,j,k,m)=spx*ofrac*rmasso/rmassq
          opy(i,j,k,m)=spy*ofrac*rmasso/rmassq
          opz(i,j,k,m)=spz*ofrac*rmasso/rmassq
          opres(i,j,k,m)=qpres(i,j,k,m)*ofrac 

          fraction=(-wind_bnd-ax)/ramp
          by(i,j,k,m)=sby_wind*fraction
          bz(i,j,k,m)=sbz_wind*fraction
         endif
c
         enddo
        enddo
       enddo
      enddo
c
c     save unperturbed density
c     !!!!!!
      if(Adiffuse)then
      do  m=ngrd,1,-1
      do  k=1,nz
      do  j=1,ny
      do  i=1,nx
       qrho0(i,j,k,m)=0.75*qrho(i,j,k,m)
       hrho0(i,j,k,m)=0.75*hrho(i,j,k,m)
       orho0(i,j,k,m)=0.75*orho(i,j,k,m)
c
       qpres0(i,j,k,m)=0.75*qpres(i,j,k,m)
       hpres0(i,j,k,m)=0.75*hpres(i,j,k,m)
       opres0(i,j,k,m)=0.75*opres(i,j,k,m)
       epres0(i,j,k,m)=0.75*epres(i,j,k,m)
      enddo
      enddo
      enddo
      enddo
      endif
c
c     initialized other important stuff
c
  170 ut=utstart+t*t_equiv/3600.
c
c     initialize plasma resistivity
c
      m=1
      rx=xspac(m)
      ry=xspac(m)
      rz=xspac(m)
      call set_resist(resistive,nx,ny,nz,resist,ijzero,numzero,
     +             ijmid,nummid,ijsrf,numsrf)
      add_dip=.false.
      save_dat=.false.
      write_dat=.false.
c
       call aurora(resistive,nx,ny,nz,m,rearth+2.0,re_equiv,1,
     +          ut,save_dat,add_dip,'res nth',10,write_dat)
       call aurora(resistive,nx,ny,nz,m,rearth+2.0,re_equiv,-1,
     +          ut,save_dat,add_dip,'res sth',10,write_dat)
c
c
c     read down relevant data list to find correct pointer
c
      if(.not.spacecraft)goto 1000
c
c      calculate the distance between the reference spacecraft and
c          solar wind boundary
c
        m=ngrd
        distance=grd_xmin(m)-rcraft(1)/re_equiv
c
c      read down data file until correct time in data file
c
c      do  n=1,ncraft
c      do  n=1,2
          n=1
          mout=20+n
          do while(ut.ge.zcraft(4,n))
           read(mout,*)zcraft(4,n),zcraft(1,n),
     +          zcraft(2,n),zcraft(3,n)
c
c          change direction to get from GSM to simulation coords
c
            zcraft(1,n)=-zcraft(1,n)
            zcraft(2,n)=-zcraft(2,n)
          if(n.eq.1)then
           rcraft(1)=zcraft(1,1)
           rcraft(2)=zcraft(2,1)
           rcraft(3)=zcraft(3,1)
          endif
c
          end do
c      end do 
c
c     zcraft(1,2)=zcraft(1,2)
c     zcraft(2,2)=zcraft(2,2)-1.2
c     zcraft(3,2)=zcraft(3,2)
c
       call limcraft(zcraft,ncraft,re_equiv,ngrd)
c
c      calculate the distance between the reference spacecraft and
c          solar wind boundary
c
       distance=grd_xmin(ngrd)-rcraft(1)/re_equiv
       write(6,*)'wind displaced by',distance
       do n=1,ncraft
        do nn=1,4
          xcraft(nn,n)=zcraft(nn,n)
        end do
       end do
c
c      do while(ut.gt.vut) 
c           svelx=zvelx
c           svely=zvely
c           svelz=zvelz
c           read(28,*)vut,zvelx,zvely,zvelz
c            zvelx=-zvelx/v_equiv
c            zvely=-zvely/v_equiv+0.03
c            zvelz=zvelz/v_equiv
c            vut=vut+t_equiv*distance/zvelx/3600.
c      end do
c
c      do while (ut.gt.rut) 
c           srho=zrho
c           read(27,*)rut,zrho
c           rut=rut+t_equiv*distance/zvelx/3600.
c           zrho=zrho/rho_equiv
c      end do
c
c    
c     initialize counting array
c
      do  j=1,ny
       do  k=1,nz
        ncount(j,k)=0
        future(j,k)=ut-0.01 
       enddo
      enddo
c
c      read all the magnetic field data to minize data sorting
c
      do  m=1,ncts
       read(29,*)bfld(m,4),bmag,bfld(m,1),bfld(m,2),bfld(m,3)
       read(27,*)rut,rplas(m),svel(m,1),svel(m,2),svel(m,3)
c      read(27,*)rut,rplas(m)
c      read(28,*)vut,svel(m,1),svel(m,2),svel(m,3)
c      warning recalibration
c         keep bx in solar wind constant
c      bfld(m,1)=-sbx_wind*b_equiv
      enddo
c
c      set timing
c 
       nvx=0
       vut=-999.
       do while(ut.gt.vut) 
        svelx=zvelx
        nvx=nvx+1
        zvelx=-svel(nvx,1)/v_equiv
        vut=bfld(nvx,4)+t_equiv*distance/zvelx/3600.
       end do
c
      write(6,*)'ut=',ut,'wind time',bfld(nvx,4)
c  
      displace=0.
      m=ngrd
      dx=(grd_xmax(m)-grd_xmin(m))/(nx-1.)
      dy=(grd_ymax(m)-grd_ymin(m))/(ny-1.)
      dz=(grd_zmax(m)-grd_zmin(m))/(nz-1.)    
c
      do 175 j=1,ny
      do 175 k=1,nz
       do while((ut.gt.future(j,k))
     +   .and.(ncount(j,k)+1.le.ncts))
        nc=ncount(j,k)+1
        bxp(j,k)=bxf(j,k)
        byp(j,k)=byf(j,k)
        bzp(j,k)=bzf(j,k)
        rhop(j,k)=rhof(j,k)
        svxp(j,k)=svxf(j,k)
        svyp(j,k)=svyf(j,k)
        svzp(j,k)=svzf(j,k)
        past(j,k)=future(j,k)
c
        future(j,k)=bfld(nc,4)
        bxf(j,k)=-bfld(nc,1)/b_equiv
        byf(j,k)=-bfld(nc,2)/b_equiv
        bzf(j,k)=bfld(nc,3)/b_equiv
        rhof(j,k)=rplas(nc)/rho_equiv
        svxf(j,k)=-svel(nc,1)/v_equiv
        svyf(j,k)=0.
        svzf(j,k)=0.
c       svyf(j,k)=-svel(nc,2)/v_equiv +0.03
c       svyf(j,k)=-svel(nc,2)/v_equiv
c       svzf(j,k)=svel(nc,3)/v_equiv
        ncount(j,k)=nc
        avx=svxf(j,k)
c
c      calculate delay
c
       if(warp)then
        b_perp=sqrt(bzf(j,k)**2+byf(j,k)**2)
        b_perp=amax1(b_perp,0.33*abs(bxf(j,k)))
        ay=grd_ymin(m)+dy*(j-1)-rcraft(2)/re_equiv
        az=grd_zmin(m)+dz*(k-1)-rcraft(3)/re_equiv
c
c       going to assume Bz IMF on average pos and ignore transients
c               and By IMF is negative
        displace=-bxf(j,k)*
     +        (ay*byf(j,k)+az*bzf(j,k))/b_perp**2
       endif
        ar=distance+displace
        future(j,k)=future(j,k)+t_equiv*ar/avx/3600.
      end do
  175 continue
c
      call set_imf(bx,by,bz,bx0,by0,bz0,bxp,byp,bzp,
     +             qrho,qpres,qpx,qpy,qpz,epres,
c     +             orho,opres,opx,opy,opz,
     +             rhop,svxp,svyp,svzp,svelx,spress,ti_te,
     +             rho_frac,nx,ny,nz,ngrd)
c
c      set boundary conditions
c
        call bndry(qrho,qpres,qpx,qpy,qpz,rmassq,
     +        hrho,hpres,hpx,hpy,hpz,rmassh,
     +        orho,opres,opx,opy,opz,rmasso,
     +        epres,bx,by,bz,bx0,by0,bz0,
     +        nx,ny,nz,ngrd,parm_srf,parm_mid,
     +        ijsrf,numsrf,ijmid,nummid,
     +        ijzero,numzero,erho,epress,alpha_e,o_conc,
     +        ti_te,srho,rho_frac,spress,spx,spy,spz,
     +        sbx_wind,sby_wind,sbz_wind)
c
c     check initial conditions
c
      call set_speed(qrho,qpres,qpx,qpy,qpz,
     +    hrho,hpres,hpx,hpy,hpz,orho,opres,opx,opy,opz,
     +    epres,bx,by,bz,bx0,by0,bz0,bsx,bsy,bsz,btot,
     +    rmassq,rmassh,rmasso,nx,ny,nz,ngrd,
     +    pxmax,pymax,pzmax,pmax,csmax,alfmax,gamma,
     +    vlim,alf_lim,o_conc,fastest)
c
      delt_old=delt
      delt=stepsz/fastest
      delt=amin1(delt,1.05*delt_old)
      delt=amax1(0.0003,delt)
c      if(t.gt.twrt)then
      twrt=twrt+dtwrt
      write(6,163)t,delt,ut,bfld(nvx,4)
  163 format(1x,'t=',1pe12.5,' dt=',1pe12.5,' ut=',
     +        1pe12.5,' wind time',1pe12.5)
c
c     write(6,161)
c 161 format(' initialization completed')
c
c
c     !!!!!! for B field check (for bcheck2.dat)
c      do m=1,ngrd
c      write(52,62)bx0(xloc,yloc,zloc,m),by0(xloc,yloc,zloc,m),
c     +bz0(xloc,yloc,zloc,m),bx(xloc,yloc,zloc,m),by(xloc,yloc,zloc,m),
c     +bz(xloc,yloc,zloc,m),ut,delt,t_equiv,m
c62    format(10(1x,1pe9.2))
c
c      enddo
c     !!!!!!
c
c     ********************************************
c     start time sequence
c     ********************************************
c
c
 1000 continue
c
c     !!!!!!den_wind changes at den_wind_time at unit of rot period
      if(ut.ge.den_wind_time*planet_per)then
      den_wind1=den_wind_new
      den_wind2=den_wind_new
      alfx_wind1=alfx_wind_new
      alfy_wind1=alfy_wind_new
      alfz_wind1=alfz_wind_new
      alfx_wind2=alfx_wind_new
      alfy_wind2=alfy_wind_new
      alfz_wind2=alfz_wind_new      
      endif
      rho_wind1=den_wind1*rmassq
      rho_wind2=den_wind2*rmassq
      delbx_wind=(alfx_wind2*sqrt(rho_wind2)
     +           -alfx_wind1*sqrt(rho_wind1))/tmax
      delby_wind=(alfy_wind2*sqrt(rho_wind2)
     +           -alfy_wind1*sqrt(rho_wind1))/tmax
      delbz_wind=(alfz_wind2*sqrt(rho_wind2)
     +           -alfz_wind1*sqrt(rho_wind1))/tmax
      sbx_wind=alfx_wind1*sqrt(rho_wind1)
      sby_wind=alfy_wind1*sqrt(rho_wind1)
      sbz_wind=alfz_wind1*sqrt(rho_wind1)
      print *, 'NewSWden',ut,den_wind1,den_wind2,alfz_wind1,alfz_wind2
c
      call cpu_time(time_end)
      write(6,*)'CPU run time(sec)= ',(time_end-time_start)/16.0,
     + '  CPU run time (hrs)= ',(time_end-time_start)/3600.0/16.0
c     !!!!!! total CPU time divided by real(ngrd) is due to
c     !!!!!! the function cpu_time includes sum time of all threads.
c     !!!!!! Here all threads are just MPI nodes number. approximately
c      !!!!!!define spinvector as spin axis unit vector
c       spinvectorx=0.5150
c       spinvectory=0.0
c       spinvectorz=0.8572
        spinvectorx=-cos((abs(tilt)-7.9)*.0174533)
        spinvectory=0.0
        spinvectorz=sin((abs(tilt)-7.9)*.0174533)
c
c     find maximum velocity to determine time step
c
c
      write(6,*)'main loop speeds'
      call set_speed(qrho,qpres,qpx,qpy,qpz,
     +    hrho,hpres,hpx,hpy,hpz,orho,opres,opx,opy,opz,
     +    epres,bx,by,bz,bx0,by0,bz0,bsx,bsy,bsz,btot,
     +    rmassq,rmassh,rmasso,nx,ny,nz,ngrd,
     +    pxmax,pymax,pzmax,pmax,csmax,alfmax,gamma,
     +    vlim,alf_lim,o_conc,fastest)
c
      delt_old=delt
      delt=stepsz/fastest
      delt=amin1(delt,1.05*delt_old)
      delt=amax1(0.0003,delt)
c      if(t.gt.twrt)then
      twrt=twrt+dtwrt
      write(6,201)t,delt,ut,bfld(nvx,4)
  201 format(1x,'t=',1pe12.5,' dt=',1pe12.5,' ut=',
     +        1pe12.5,' wind time',1pe12.5)
c
      delt2=delt/2.
c     if(delt.lt.0.0003)stop
      t=t+delt
      utold=ut
      ut=utstart+t*t_equiv/3600.
      delay=t_equiv*distance/svelx/3600.
c
c
c     write out data if necessary - only using course gridding
c        the momemt
c
c     if(spacecraft) then
       do n=1,ncraft
       call qvset(0.,bsx,nx*ny*nz)
       call qvset(0.,bsy,nx*ny*nz)
       call qvset(0.,bsz,nx*ny*nz)
c
c
       m=1
       do while ((xcraft(1,n).gt.grd_xmax(m)*re_equiv).or.
     +      (xcraft(1,n).lt.grd_xmin(m)*re_equiv).or.
     +      (xcraft(2,n).gt.grd_ymax(m)*re_equiv).or.
     +      (xcraft(2,n).lt.grd_ymin(m)*re_equiv).or.
     +      (xcraft(3,n).gt.grd_zmax(m)*re_equiv).or.
     +      (xcraft(3,n).lt.grd_zmin(m)*re_equiv).and.
     +       (m+1.le.ngrd)) 
             m=m+1
       enddo
       rx=xspac(m)
       ry=rx
       rz=rz
c
       call totfld(bx,bx0,bsx,nx,ny,nz,ngrd,m)
       call totfld(by,by0,bsy,nx,ny,nz,ngrd,m)
       call totfld(bz,bz0,bsz,nx,ny,nz,ngrd,m)
c
       add_dip=.false.
       call crafdatv(bsx,bsy,bsz,qpx,qpy,qpz,qrho,qpres,
     +       rmassq,hpx,hpy,hpz,hrho,hpres,rmassh,
     +       opx,opy,opz,orho,opres,rmasso,epres,
     +       nx,ny,nz,ngrd,m,xcraft,ncraft,n,ut,add_dip,
     +       re_equiv,b_equiv,v_equiv,rho_equiv,gamma)
       enddo
c      endif
c
c     test to see whether scraft positions need to be updated
c
      if(spacecraft)then
c      do 210 n=1,ncraft
c      do 210 n=1,2
       n=1
          if(ut.ge.zcraft(4,n))then
           mout=20+n
           read(mout,*)zcraft(4,n),zcraft(1,n),
     +          zcraft(2,n),zcraft(3,n)
c          if(n.eq.1)zcraft(4,n)=zcraft(4,n)/3600.
c 
c          change direction to get from GSM to simulation coords
c
            zcraft(1,n)=-zcraft(1,n)
            zcraft(2,n)=-zcraft(2,n)
             if(n.eq.1)then
              rcraft(1)=zcraft(1,1)
              rcraft(2)=zcraft(2,1)
              rcraft(3)=zcraft(3,1)
             endif
c
          endif
c 210 continue
c
      zcraft(4,3)=zcraft(4,2)
      zcraft(4,4)=zcraft(4,2)
c
c     zcraft(1,2)=zcraft(1,2)
c     zcraft(2,2)=zcraft(2,2)-1.2
c     zcraft(3,2)=zcraft(3,2)
c
c         set refernce spacecraft position
c                   and spacecraft limits
c
          call limcraft(zcraft,ncraft,re_equiv,ngrd)
c
c         set density and velocity
c
         distance=grd_xmin(ngrd)-rcraft(1)/re_equiv
         do while (ut.ge.vut)
          nvx=nvx+1
          zvelx=-svel(nvx,1)/v_equiv
          zvely=-svel(nvx,2)/v_equiv
c         zvely=-svel(nvx,2)/v_equiv+0.03
          zvelz=svel(nvx,3)/v_equiv
          vut=bfld(nvx,4)+t_equiv*distance/zvelx/3600.
c          read(28,*)vut,zvelx,zvely,zvelz
c            zvelx=-zvelx/v_equiv
c            zvely=-zvely/v_equiv+0.03
c            zvelz=zvelz/v_equiv
c            vut=vut+t_equiv*distance/zvelx/3600.
          end do
c         do while (ut.ge.rut)
c          read(27,*)rut,zrho
c            rut=rut+t_equiv*distance/zvelx/3600.
c            zrho=zrho/rho_equiv
c         end do
c
c        fix up magnetic field
c
          displace=0.
          do 220 j=1,ny
          do 220 k=1,nz
           do while((ut.ge.future(j,k))
     +             .and.(ncount(j,k)+1.le.ncts))
           nc=ncount(j,k)+1
           future(j,k)=bfld(nc,4)
           bxf(j,k)=-bfld(nc,1)/b_equiv
           byf(j,k)=-bfld(nc,2)/b_equiv
           bzf(j,k)=bfld(nc,3)/b_equiv
           rhof(j,k)=rplas(nc)/rho_equiv
           svxf(j,k)=-svel(nc,1)/v_equiv
           svyf(j,k)=0.00
           svzf(j,k)=0.0
c          svyf(j,k)=-svel(nc,2)/v_equiv
c          svyf(j,k)=-svel(nc,2)/v_equiv+0.03
c          svzf(j,k)=svel(nc,3)/v_equiv
           avx=svxf(j,k)
           ncount(j,k)=nc
c
c      calculate delay
c
          if(warp)then
c           b_perp=sqrt(bzf(j,k)**2+byf(j,k)**2)
c           b_perp=amax1(b_perp,0.1*abs(bxf(j,k)))
            ay=(j-1.)*xspac(ngrd)+grd_ymin(ngrd)-rcraft(2)/re_equiv
            az=(k-1.)*xspac(ngrd)+grd_zmin(ngrd)-rcraft(3)/re_equiv
c
c       going to assume Bz IMF on average pos and ignore transients
c               and By IMF is negative
c
             b_perp=sqrt(bzf(j,k)**2+byf(j,k)**2) 
             b_perp=amax1(b_perp,0.33*abs(bxf(j,k))) 
             displace=-bxf(j,k)*
     +        (ay*byf(j,k)+az*bzf(j,k))/b_perp**2
           endif
           ar=distance+displace
           future(j,k)=future(j,k)+t_equiv*ar/avx/3600.
           end do
  220   continue
c
      endif
c     !!!!!!spacecraft is .f. so nothing running below
      if(spacecraft)then
c
c     update spacecraft position and magnetic field
c
      do 250 n=1,ncraft
         dut=(ut-xcraft(4,n))/(zcraft(4,n)-utold)
         xcraft(1,n)=xcraft(1,n)+dut*(zcraft(1,n)-xcraft(1,n))
         xcraft(2,n)=xcraft(2,n)+dut*(zcraft(2,n)-xcraft(2,n))
         xcraft(3,n)=xcraft(3,n)+dut*(zcraft(3,n)-xcraft(3,n))
         xcraft(4,n)=ut
  250 continue
c
        distance=grd_xmin(ngrd)-rcraft(1)/re_equiv
        call limcraft(xcraft,ncraft,re_equiv,ngrd)
c
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
c
       dut=(ut-utold)/(vut-utold)
       svelx=svelx+dut*(zvelx-svelx)
       svely=0.
       svelz=svelz+delvz_wind*delt
       srho=srho/float(nz*ny)
c
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
c  
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
c      !!!!!!spacecraft is .f. so nothing running above
      delay=t_equiv*distance/svelx/3600.
c
c       test for div B errors and reset ionospheric restivity
c
c     if((.not.start).and.(t.gt.tdiv))then
      if(t.gt.tdiv)then
          do m=ngrd,1,-1
          rx=xspac(m)
          ry=rx
          rz=rx
          call divb_cor(bx,by,bz,qpx,qpy,qpz,qrho,qpres,
     +              ngrd,nx,ny,nz,m,srho)
          nb=m
          ns=nbndry(m)
          mb=m
          ms=nbndry(m)
          if(ms.ne.0) call flanks(bx,ngrd,ns,nb,ms,mb,nx,ny,nz)
          enddo
       endif
       print *, 'tilting', tilting
       print *, 'NOtilting', .NOT.tilting
c        !!!!!!rotating part for dipole fied around z axis
         !!!!!!-sin,-cos is due to anticlockwise rotation
         spin=spin+dspin*delt
c        sin_spin=sin(spin*.0174533)
c        cos_spin=cos(spin*.0174533)
        if(ut.ge.rot_time*planet_per)then
         sin_spin=sin(2.*3.1415926*ut/planet_per)
         cos_spin=cos(2.*3.1415926*ut/planet_per)
        else
         sin_spin=sin(0.0)
         cos_spin=cos(0.0)
        endif
c          sin_spin=sin(180*.0174533)
c          cos_spin=cos(180*.0174533)
        print *, 'do i call this?'
        call mak_dipz(bx0,by0,bz0,nx,ny,nz,ngrd,ijzero,numzero)
c        !!!!!!
        if(t.gt.tdiv)then
          call set_resist(resistive,nx,ny,nz,resist,ijzero,numzero,
     +             ijmid,nummid,ijsrf,numsrf)
        endif
c        !!!!!!
      if(tilting)then 
c        tilt=tilt+dtilt*delt
c        sin_tilt=sin(tilt*.0174533)
c        cos_tilt=cos(tilt*.0174533)
c        call mak_dip(bx0,by0,bz0,nx,ny,nz,ngrd,ijzero,numzero)
c        !!!!!!Do I need to rotate bx,y,z0 here?
c              Do so if tilting=.t. ?
c        !!!!!!make dipole field rotating round z axis
c              (bx,y,z0a bx,y,zb not defined yet)
c         bx0a(i,j,k,m)=bx0(i,j,k,m)*cos_spin
c     +                -by0(i,j,k,m)*sin_spin
c         by0a(i,j,k,m)=bx0(i,j,k,m)*sin_spin
c     +                +by0(i,j,k,m)*cos_spin
c         bz0a(i,j,k,m)=bz0(i,j,k,m)
c        !!
c          bx0a(i,j,k,m)=bx0(i,j,k,m)
c          by0a(i,j,k,m)=by0(i,j,k,m)
c          bz0a(i,j,k,m)=bz0(i,j,k,m)
c        !!!!!!make rotational dipole field tilted around y axis
c              (right-hand coordinates)
c         bx0b(i,j,k,m)=bz0a(i,j,k,m)*sin_tilt3
c     +                +bx0a(i,j,k,m)*cos_tilt3
c         by0b(i,j,k,m)=by0a(i,j,k,m)
c         bz0b(i,j,k,m)=bz0a(i,j,k,m)*cos_tilt3
c     +                -1.*bx0a(i,j,k,m)*sin_tilt3
c        !!!!!!bx,y,z0a is the after-rotation field
c              bx,y,z0b is the field after-tilted towards the Sun
c         bx0(i,j,k,m)=bx0b(i,j,k,m)
c         by0(i,j,k,m)=by0b(i,j,k,m)
c         bz0(i,j,k,m)=bz0b(i,j,k,m)
c        !!
c        !!!!!!
c        !!!!!!rotating part for dipole fied around z axis
c         sin_spin=sin(2.*3.1415926*ut/planet_per)
c         cos_spin=cos(2.*3.1415926*ut/planet_per)
        print *, 'do i call this 2 ?'
c        call mak_dipz(bx0,by0,bz0,nx,ny,nz,ngrd,ijzero,numzero)
c        !!!!!!
c        !!!!!!
        if(t.gt.tdiv)then
          call set_resist(resistive,nx,ny,nz,resist,ijzero,numzero,
     +             ijmid,nummid,ijsrf,numsrf)
        endif
      endif
c
c     !!!!!! for B field check (for bcheck3.dat)
c      do m=1,ngrd
c      write(53,63)bx0(xloc,yloc,zloc,m),by0(xloc,yloc,zloc,m),
c     +bz0(xloc,yloc,zloc,m),bx(xloc,yloc,zloc,m),by(xloc,yloc,zloc,m),
c     +bz(xloc,yloc,zloc,m),ut,delt,t_equiv,m
c63    format(10(1x,1pe9.2))
c      enddo
c     !!!!!!
c
c     **********************************************************
c     Two Step Lax-Wendroff : step 1
c          estimate of the fluid quantites at n+1/2 
c     **********************************************************
c
c     store initial plasma parameters
c
      oldqrho=qrho
      oldqpres=qpres
      oldqpx=qpx
      oldqpy=qpy
      oldqpz=qpz
c
      oldhrho=hrho
      oldhpres=hpres
      oldhpx=hpx
      oldhpy=hpy
      oldhpz=hpz
c
      oldorho=orho
      oldopres=opres
      oldopx=opx
      oldopy=opy
      oldopz=opz
c
      oldepres=epres
      oldbx=bx
      oldby=by
      oldbz=bz
c     !!!!!!
c      write(*,*)'for old and wrk parameters'
c      write(*,*)oldqrho0_temp(1,10,10,10),wrkqrho0_temp(1,10,10,10),
c     +oldqpx0_temp(1,10,10,10),wrkqpx0_temp(1,10,10,10),
c     +oldhpx0_temp(1,10,10,10),
c     +oldhpx0_temp(1,11,15,10)
c
c     calculate standard mhd current j = curl B
c
c
      do 310 m=1,ngrd
       rx=xspac(m)
       ry=xspac(m)
       rz=xspac(m)
c
c      write(6,*)' calcur ing now'
       call calcur(bx,by,bz,ngrd,nx,ny,nz,m,curx,cury,curz,
     +         ijsrf,numsrf,ijmid,nummid,ijzero,numzero)
c     !!!!!! do I need to change bx.. to bsx..? b/c B keeps changing
c     find total magnetic field
c
c      write(6,*)' totbfld ing now'
      call totfld(bx,bx0,bsx,nx,ny,nz,ngrd,m)
      call totfld(by,by0,bsy,nx,ny,nz,ngrd,m) 
      call totfld(bz,bz0,bsz,nx,ny,nz,ngrd,m)
c
c     find magnitude of B
c
      call tot_b(btot,bsx,bsy,bsz,nx,ny,nz)
c
c
c      write(6,*)' fnd_evel ing now'
       call fnd_evel(qpx,qpy,qpz,qrho,hpx,hpy,hpz,hrho, 
     +       opx,opy,opz,orho,curx,cury,curz,evx,evy,evz,
     +       tvx,tvy,tvz,ngrd,
     +       nx,ny,nz,m,rmassq,rmassh,rmasso,reynolds)
c
c      write(6,*)' bande ing now'
      call bande(efldx,efldy,efldz,bsx,bsy,bsz,
     +       curx,cury,curz,evx,evy,evz,btot,
     +       epres,qrho,hrho,orho,resistive,resist,reynolds,
     +       nx,ny,nz,ngrd,m,rmassq,rmassh,rmasso,
     +       ijmid,nummid,ijzero,numzero)
c      !!!!!! to bring in the unpertubed parameter for rho,pres,Pxyz
c       oldqrho=oldqrho-qrho0_temp
c       oldhrho=oldhrho-hrho0_temp
c       oldorho=oldorho-orho0_temp
c       oldqpres=oldqpres-qpres0_temp
c       oldhpres=oldhpres-hpres0_temp
c       oldopres=oldopres-opres0_temp
c       oldepres=oldepres-epres0_temp
c       oldqpx=oldqpx-qpx0_temp
c       oldqpy=oldqpy-qpy0_temp
c       oldqpz=oldqpz-qpz0_temp
c       oldhpx=oldhpx-hpx0_temp
c       oldhpy=oldhpy-hpy0_temp
c       oldhpz=oldhpz-hpz0_temp
c       oldopx=oldopx-opx0_temp
c       oldopy=oldopy-opy0_temp
c       oldopz=oldopz-opz0_temp
c
c       wrkqrho=wrkqrho-qrho0_temp
c       wrkhrho=wrkhrho-hrho0_temp
c       wrkorho=wrkorho-orho0_temp
c       wrkqpres=wrkqpres-qpres0_temp
c       wrkhpres=wrkhpres-hpres0_temp
c       wrkopres=wrkopres-opres0_temp
c       wrkepres=wrkepres-epres0_temp
c       wrkqpx=wrkqpx-qpx0_temp
c       wrkqpy=wrkqpy-qpy0_temp
c       wrkqpz=wrkqpz-qpz0_temp
c       wrkhpx=wrkhpx-hpx0_temp
c       wrkhpy=wrkhpy-hpy0_temp
c       wrkhpz=wrkhpz-hpz0_temp
c       wrkopx=wrkopx-opx0_temp
c       wrkopy=wrkopy-opy0_temp
c       wrkopz=wrkopz-opz0_temp
c
c       qrho=qrho-qrho0_temp
c       hrho=hrho-hrho0_temp
c       orho=orho-orho0_temp
c       qpres=qpres-qpres0_temp
c       hpres=hpres-hpres0_temp
c       opres=opres-opres0_temp
c       epres=epres-epres0_temp
c       qpx=qpx-qpx0_temp
c       qpy=qpy-qpy0_temp
c       qpz=qpz-qpz0_temp
c       hpx=hpx-hpx0_temp
c       hpy=hpy-hpy0_temp
c       hpz=hpz-hpz0_temp
c       opx=opx-opx0_temp
c       opy=opy-opy0_temp
c       opz=opz-opz0_temp
c      !!!!!!
c      write(6,*)' push elec ing now'
       call push_elec(wrkepres,oldepres,epres,evx,evy,evz,
     +        gamma,gamma1,ngrd,nx,ny,nz,m,0.5*delt)
c
c      write(6,*)' push ion 1 ing now'
       call push_ion(wrkqrho,wrkqpres,wrkqpx,wrkqpy,wrkqpz,
     +        oldqrho,oldqpres,oldqpx,oldqpy,oldqpz,
     +        qrho,qpres,qpx,qpy,qpz,
     +        bsx,bsy,bsz,btot,efldx,efldy,efldz,rmassq,
     +        vvx,vvy,vvz,tvx,tvy,tvz,gamma,gamma1,
     +        ngrd,nx,ny,nz,m,0.5*delt,grav,re_equiv,reynolds)

c      write(6,*)' push ion 2 ing now'
       call push_ion(wrkhrho,wrkhpres,wrkhpx,wrkhpy,wrkhpz,
     +        oldhrho,oldhpres,oldhpx,oldhpy,oldhpz,
     +        hrho,hpres,hpx,hpy,hpz,
     +        bsx,bsy,bsz,btot,efldx,efldy,efldz,rmassh,
     +        vvx,vvy,vvz,tvx,tvy,tvz,gamma,gamma1,
     +        ngrd,nx,ny,nz,m,0.5*delt,grav,re_equiv,reynolds)
c
c      write(6,*)' push ion 3 ing now'
       call push_ion(wrkorho,wrkopres,wrkopx,wrkopy,wrkopz,
     +        oldorho,oldopres,oldopx,oldopy,oldopz,
     +        orho,opres,opx,opy,opz,
     +        bsx,bsy,bsz,btot,efldx,efldy,efldz,rmasso,
     +        vvx,vvy,vvz,tvx,tvy,tvz,gamma,gamma1,
     +        ngrd,nx,ny,nz,m,0.5*delt,grav,re_equiv,reynolds)
c
c      !!!!!! to back for the unpertubed parameter for rho,pres,Pxyz
c       oldqrho=oldqrho+qrho0_temp
c       oldhrho=oldhrho+hrho0_temp
c       oldorho=oldorho+orho0_temp
c       oldqpres=oldqpres+qpres0_temp
c       oldhpres=oldhpres+hpres0_temp
c       oldopres=oldopres+opres0_temp
c       oldepres=oldepres+epres0_temp
c       oldqpx=oldqpx+qpx0_temp
c       oldqpy=oldqpy+qpy0_temp
c       oldqpz=oldqpz+qpz0_temp
c       oldhpx=oldhpx+hpx0_temp
c       oldhpy=oldhpy+hpy0_temp
c       oldhpz=oldhpz+hpz0_temp
c       oldopx=oldopx+opx0_temp
c       oldopy=oldopy+opy0_temp
c       oldopz=oldopz+opz0_temp
c
c       wrkqrho=wrkqrho+qrho0_temp
c       wrkhrho=wrkhrho+hrho0_temp
c       wrkorho=wrkorho+orho0_temp
c       wrkqpres=wrkqpres+qpres0_temp
c       wrkhpres=wrkhpres+hpres0_temp
c       wrkopres=wrkopres+opres0_temp
c       wrkepres=wrkepres+epres0_temp
c       wrkqpx=wrkqpx+qpx0_temp
c       wrkqpy=wrkqpy+qpy0_temp
c       wrkqpz=wrkqpz+qpz0_temp
c       wrkhpx=wrkhpx+hpx0_temp
c       wrkhpy=wrkhpy+hpy0_temp
c       wrkhpz=wrkhpz+hpz0_temp
c       wrkopx=wrkopx+opx0_temp
c       wrkopy=wrkopy+opy0_temp
c       wrkopz=wrkopz+opz0_temp
c
c       qrho=qrho+qrho0_temp
c       hrho=hrho+hrho0_temp
c       orho=orho+orho0_temp
c       qpres=qpres+qpres0_temp
c       hpres=hpres+hpres0_temp
c       opres=opres+opres0_temp
c       epres=epres+epres0_temp
c       qpx=qpx+qpx0_temp
c       qpy=qpy+qpy0_temp
c       qpz=qpz+qpz0_temp
c       hpx=hpx+hpx0_temp
c       hpy=hpy+hpy0_temp
c       hpz=hpz+hpz0_temp
c       opx=opx+opx0_temp
c       opy=opy+opy0_temp
c       opz=opz+opz0_temp
c      !!!!!!
c      !!!!!!cal total perturbed B field for step 1, and
c            oldbx is just perturbed B field, true?
       call push_bfld(wrkbx,wrkby,wrkbz,oldbx,oldby,oldbz,
     +              efldx,efldy,efldz,ngrd,nx,ny,nz,m,0.5*delt)  
  310 continue
c
c     write(6,988)
c 988 format(' main lax loop now')
c
c     *************************************************************
c     Apply boundary conditions
c     *************************************************************
c
c     write(6,989)
c 989 format(' doing bndry conditions')
c
      call bndry(wrkqrho,wrkqpres,wrkqpx,wrkqpy,wrkqpz,rmassq,
     +        wrkhrho,wrkhpres,wrkhpx,wrkhpy,wrkhpz,rmassh,
     +        wrkorho,wrkopres,wrkopx,wrkopy,wrkopz,rmasso,
     +        wrkepres,wrkbx,wrkby,wrkbz,bx0,by0,bz0,
     +        nx,ny,nz,ngrd,parm_srf,parm_mid,
     +        ijsrf,numsrf,ijmid,nummid,
     +        ijzero,numzero,erho,epress,alpha_e,o_conc,
     +        ti_te,srho,rho_frac,spress,spx,spy,spz,
     +        sbx_wind,sby_wind,sbz_wind)
      if(spacecraft)then
         call set_imf(wrkbx,wrkby,wrkbz,bx0,by0,bz0,bxp,byp,bzp,
     +         wrkqrho,wrkqpres,wrkqpx,wrkqpy,wrkqpz,wrkepres,
     +         rhop,svxp,svyp,svzp,svelx,spress,ti_te,
     +         rho_frac,nx,ny,nz,ngrd)
      endif
c      write(6,*)'Lax 1 speeds'
      call set_speed(wrkqrho,wrkqpres,wrkqpx,wrkqpy,wrkqpz,
     +    wrkhrho,wrkhpres,wrkhpx,wrkhpy,wrkhpz,
     +    wrkorho,wrkopres,wrkopx,wrkopy,wrkopz,
     +    wrkepres,wrkbx,wrkby,wrkbz,bx0,by0,bz0,
     +    bsx,bsy,bsz,btot,rmassq,rmassh,rmasso,nx,ny,nz,ngrd,
     +    pxmax,pymax,pzmax,pmax,csmax,alfmax,gamma,
     +    vlim,alf_lim,o_conc,fastest)
c
c       ********rotational part for dipole field********
c         do m=1,ngrd
c         do k=1,nz
c         do j=1,ny
c         do i=1,nx
c       !!!!!!Bx,y,z0 tilted back to orginal(-(-97.9) around y axis)
c              -1.*tilt3
c       !!
c         bx0a1=-1.*bz0*sin_tilt3
c     +                 +bx0*cos_tilt3
c         by0a1=by0
c         bz0a1=bz0*cos_tilt3
c     +                 +bx0*sin_tilt3
c       !!!!!!and then rotate w*(ut+0.5*delt*t_equiv/3600./2)*3600
c              around z axis
c         bx0a2=bx0a1*cos(2.*3.1415926*(ut
c     +           +0.5*delt*t_equiv/3600.)/planet_per)
c     +                 -by0a1*sin(2.*3.1415926*(ut
c     +           +0.5*delt*t_equiv/3600.)/planet_per)
c         by0a2=bx0a1*sin(2.*3.1415926*(ut
c     +           +0.5*delt*t_equiv/3600.)/planet_per)
c     +                 +by0a1*cos(2.*3.1415926*(ut
c     +           +0.5*delt*t_equiv/3600.)/planet_per)
c         bz0a2=bz0a1
c       !!!!!! Bx,y,z0 tilted back nearly towards Sun again
c               (+97.9 around y axis) +1.tilt3
c         bx0a3=bz0a2*sin_tilt3+bx0a2*cos_tilt3
c         by0a3=by0a2
c         bz0a3=bz0a2*cos_tilt3-1.*bx0a2*sin_tilt3
c       
c       !!!!!!
c        bx0=bx0a3
c        by0=by0a3
c        bz0=bz0a3
c       !!
c
c       !!!!!!
c         end do
c         end do
c         end do
c         end do
c       ************************************************
c
c     !!!!!! for B field check (for bcheck4.dat)
c      do m=1,ngrd
c      write(54,64)bx0(xloc,yloc,zloc,m),by0(xloc,yloc,zloc,m),
c     +bz0(xloc,yloc,zloc,m),bx(xloc,yloc,zloc,m),by(xloc,yloc,zloc,m),
c     +bz(xloc,yloc,zloc,m),ut,delt,t_equiv,m
c64    format(10(1x,1pe9.2))
c      enddo
c     !!!!!!
c
c      ***********************************************************
c      Lax-Wendroff: step 2
c            use the predicted value to find corrected value for n+1
c      ***********************************************************
c 
      do 410 m=1,ngrd
       rx=xspac(m)
       ry=xspac(m)
       rz=xspac(m)
c
c     calculate standard mhd current j = curl B
c
      call calcur(wrkbx,wrkby,wrkbz,ngrd,nx,ny,nz,m,curx,cury,curz,
     +            ijsrf,numsrf,ijmid,nummid,ijzero,numzero)
c
c     find total magnetic field
c
      call totfld(wrkbx,bx0,bsx,nx,ny,nz,ngrd,m)
      call totfld(wrkby,by0,bsy,nx,ny,nz,ngrd,m)
      call totfld(wrkbz,bz0,bsz,nx,ny,nz,ngrd,m)
c
c     find magnitude of B
c
      call tot_b(btot,bsx,bsy,bsz,nx,ny,nz)
c
c
c     find the  electric field from electron momentum eqn
c
       call fnd_evel(wrkqpx,wrkqpy,wrkqpz,wrkqrho,
     +       wrkhpx,wrkhpy,wrkhpz,wrkhrho, 
     +       wrkopx,wrkopy,wrkopz,wrkorho,
     +       curx,cury,curz,evx,evy,evz,
     +       tvx,tvy,tvz,ngrd,
     +       nx,ny,nz,m,rmassq,rmassh,rmasso,reynolds)
cc
      call bande(efldx,efldy,efldz,bsx,bsy,bsz,
     +       curx,cury,curz,evx,evy,evz,btot,
     +       wrkepres,wrkqrho,wrkhrho,wrkorho,
     +       resistive,resist,reynolds,
     +       nx,ny,nz,ngrd,m,rmassq,rmassh,rmasso,
     +       ijmid,nummid,ijzero,numzero)
c
c      !!!!!! to bring in the unpertubed parameter for rho,pres,Pxyz
c       oldqrho=oldqrho-qrho0_temp
c       oldhrho=oldhrho-hrho0_temp
c       oldorho=oldorho-orho0_temp
c       oldqpres=oldqpres-qpres0_temp
c       oldhpres=oldhpres-hpres0_temp
c       oldopres=oldopres-opres0_temp
c       oldepres=oldepres-epres0_temp
c       oldqpx=oldqpx-qpx0_temp
c       oldqpy=oldqpy-qpy0_temp
c       oldqpz=oldqpz-qpz0_temp
c       oldhpx=oldhpx-hpx0_temp
c       oldhpy=oldhpy-hpy0_temp
c       oldhpz=oldhpz-hpz0_temp
c       oldopx=oldopx-opx0_temp
c       oldopy=oldopy-opy0_temp
c       oldopz=oldopz-opz0_temp
c
c       wrkqrho=wrkqrho-qrho0_temp
c       wrkhrho=wrkhrho-hrho0_temp
c       wrkorho=wrkorho-orho0_temp
c       wrkqpres=wrkqpres-qpres0_temp
c       wrkhpres=wrkhpres-hpres0_temp
c       wrkopres=wrkopres-opres0_temp
c       wrkepres=wrkepres-epres0_temp
c       wrkqpx=wrkqpx-qpx0_temp
c       wrkqpy=wrkqpy-qpy0_temp
c       wrkqpz=wrkqpz-qpz0_temp
c       wrkhpx=wrkhpx-hpx0_temp
c       wrkhpy=wrkhpy-hpy0_temp
c       wrkhpz=wrkhpz-hpz0_temp
c       wrkopx=wrkopx-opx0_temp
c       wrkopy=wrkopy-opy0_temp
c       wrkopz=wrkopz-opz0_temp
c
c       qrho=qrho-qrho0_temp
c       hrho=hrho-hrho0_temp
c       orho=orho-orho0_temp
c       qpres=qpres-qpres0_temp
c       hpres=hpres-hpres0_temp
c       opres=opres-opres0_temp
c       epres=epres-epres0_temp
c       qpx=qpx-qpx0_temp
c       qpy=qpy-qpy0_temp
c       qpz=qpz-qpz0_temp
c       hpx=hpx-hpx0_temp
c       hpy=hpy-hpy0_temp
c       hpz=hpz-hpz0_temp
c       opx=opx-opx0_temp
c       opy=opy-opy0_temp
c       opz=opz-opz0_temp
c      !!!!!!
        call push_elec(epres,oldepres,wrkepres,evx,evy,evz,
     +        gamma,gamma1,ngrd,nx,ny,nz,m,delt)
c
        call push_ion(qrho,qpres,qpx,qpy,qpz,
     +        oldqrho,oldqpres,oldqpx,oldqpy,oldqpz,
     +        wrkqrho,wrkqpres,wrkqpx,wrkqpy,wrkqpz,
     +        bsx,bsy,bsz,btot,efldx,efldy,efldz,rmassq,
     +        vvx,vvy,vvz,tvx,tvy,tvz,gamma,gamma1,
     +        ngrd,nx,ny,nz,m,delt,grav,re_equiv,reynolds)
c
        call push_ion(hrho,hpres,hpx,hpy,hpz,
     +        oldhrho,oldhpres,oldhpx,oldhpy,oldhpz,
     +        wrkhrho,wrkhpres,wrkhpx,wrkhpy,wrkhpz,
     +        bsx,bsy,bsz,btot,efldx,efldy,efldz,rmassh,
     +        vvx,vvy,vvz,tvx,tvy,tvz,gamma,gamma1,
     +        ngrd,nx,ny,nz,m,delt,grav,re_equiv,reynolds)
c
        call push_ion(orho,opres,opx,opy,opz,
     +        oldorho,oldopres,oldopx,oldopy,oldopz,
     +        wrkorho,wrkopres,wrkopx,wrkopy,wrkopz,
     +        bsx,bsy,bsz,btot,efldx,efldy,efldz,rmasso,
     +        vvx,vvy,vvz,tvx,tvy,tvz,gamma,gamma1,
     +        ngrd,nx,ny,nz,m,delt,grav,re_equiv,reynolds)
c
c      !!!!!! to back for the unpertubed parameter for rho,pres,Pxyz
c       oldqrho=oldqrho+qrho0_temp
c       oldhrho=oldhrho+hrho0_temp
c       oldorho=oldorho+orho0_temp
c       oldqpres=oldqpres+qpres0_temp
c       oldhpres=oldhpres+hpres0_temp
c       oldopres=oldopres+opres0_temp
c       oldepres=oldepres+epres0_temp
c       oldqpx=oldqpx+qpx0_temp
c       oldqpy=oldqpy+qpy0_temp
c       oldqpz=oldqpz+qpz0_temp
c       oldhpx=oldhpx+hpx0_temp
c       oldhpy=oldhpy+hpy0_temp
c       oldhpz=oldhpz+hpz0_temp
c       oldopx=oldopx+opx0_temp
c       oldopy=oldopy+opy0_temp
c       oldopz=oldopz+opz0_temp
c
c       wrkqrho=wrkqrho+qrho0_temp
c       wrkhrho=wrkhrho+hrho0_temp
c       wrkorho=wrkorho+orho0_temp
c       wrkqpres=wrkqpres+qpres0_temp
c       wrkhpres=wrkhpres+hpres0_temp
c       wrkopres=wrkopres+opres0_temp
c       wrkepres=wrkepres+epres0_temp
c       wrkqpx=wrkqpx+qpx0_temp
c       wrkqpy=wrkqpy+qpy0_temp
c       wrkqpz=wrkqpz+qpz0_temp
c       wrkhpx=wrkhpx+hpx0_temp
c       wrkhpy=wrkhpy+hpy0_temp
c       wrkhpz=wrkhpz+hpz0_temp
c       wrkopx=wrkopx+opx0_temp
c       wrkopy=wrkopy+opy0_temp
c       wrkopz=wrkopz+opz0_temp
c
c       qrho=qrho+qrho0_temp
c       hrho=hrho+hrho0_temp
c       orho=orho+orho0_temp
c       qpres=qpres+qpres0_temp
c       hpres=hpres+hpres0_temp
c       opres=opres+opres0_temp
c       epres=epres+epres0_temp
c       qpx=qpx+qpx0_temp
c       qpy=qpy+qpy0_temp
c       qpz=qpz+qpz0_temp
c       hpx=hpx+hpx0_temp
c       hpy=hpy+hpy0_temp
c       hpz=hpz+hpz0_temp
c       opx=opx+opx0_temp
c       opy=opy+opy0_temp
c       opz=opz+opz0_temp
c      !!!!!!
c
       call push_bfld(bx,by,bz,oldbx,oldby,oldbz,
     +              efldx,efldy,efldz,ngrd,nx,ny,nz,m,delt) 
c 
  410 continue
c
c     write(6,992)
c 992 format(' main 2nd lax loop now')
c
c     *************************************************************
c     Apply boundary conditions
c     *************************************************************
c
c
      call bndry(qrho,qpres,qpx,qpy,qpz,rmassq,
     +         hrho,hpres,hpx,hpy,hpz,rmassh,
     +         orho,opres,opx,opy,opz,rmasso,
     +         epres,bx,by,bz,bx0,by0,bz0,
     +         nx,ny,nz,ngrd,parm_srf,parm_mid,
     +         ijsrf,numsrf,ijmid,nummid,
     +         ijzero,numzero,erho,epress,alpha_e,o_conc,
     +         ti_te,srho,rho_frac,spress,spx,spy,spz,
     +         sbx_wind,sby_wind,sbz_wind)
      if(spacecraft)then
         call set_imf(bx,by,bz,bx0,by0,bz0,bxp,byp,bzp,
     +             qrho,qpres,qpx,qpy,qpz,epres,
     +             rhop,svxp,svyp,svzp,svelx,spress,ti_te,
     +             rho_frac,nx,ny,nz,ngrd)
      endif
c
c      write(6,*)'lax 2 speeds'
      call set_speed(qrho,qpres,qpx,qpy,qpz,
     +    hrho,hpres,hpx,hpy,hpz,orho,opres,opx,opy,opz,
     +    epres,bx,by,bz,bx0,by0,bz0,bsx,bsy,bsz,btot,
     +    rmassq,rmassh,rmasso,nx,ny,nz,ngrd,
     +    pxmax,pymax,pzmax,pmax,csmax,alfmax,gamma,
     +    vlim,alf_lim,o_conc,fastest)
c       ********rotational part for dipole field********
c         do m=1,ngrd
c         do k=1,nz
c         do j=1,ny
c         do i=1,nx
c       !!!!!!Bx,y,z0 tilted back to orginal(-97.9 around y axis)
c              -1.*tilt3
c       !!
c         bx0b1=-1.*bz0*sin_tilt3
c     +                 +bx0*cos_tilt3
c         by0b1=by0
c         bz0b1=bz0*cos_tilt3
c     +                 +bx0*sin_tilt3
c
c          bx0b1=bx0
c          by0b1=by0
c          bz0b1=bz0
c       !!!!!!and then rotate w*(ut+delt*t_equiv/3600./2)*3600
c              around z axis.
c       !!!!!!here should not be ut+**** but simply ut, because
c         ut has already updated by every step following mark 1000
c         bx0b2=bx0b1*cos(2.*3.1415926*(ut
c     +           +1.0*delt*t_equiv/3600.)/planet_per)
c     +                 -by0b1*sin(2.*3.1415926*(ut
c     +           +1.0*delt*t_equiv/3600.)/planet_per)
c         by0b2=bx0b1*sin(2.*3.1415926*(ut
c     +           +1.0*delt*t_equiv/3600.)/planet_per)
c     +                 +by0b1*cos(2.*3.1415926*(ut
c     +           +1.0*delt*t_equiv/3600.)/planet_per)
c         bz0b2=bz0b1
c       !!!!!! Bx,y,z0 tilted back nearly towards Sun again
c               (+97.9 around y axis) +1.tilt3
c         bx0b3=bz0b2*sin_tilt3+bx0b2*cos_tilt3
c         by0b3=by0b2
c         bz0b3=bz0b2*cos_tilt3-1.*bx0b2*sin_tilt3
c       !!!!!!
c          bx0=bx0b3
c          by0=by0b3
c          bz0=bz0b3
c         bx0=bx0b2
c         by0=by0b2
c         bz0=bz0b2
c       !!
c       !!!!!!
c         end do
c         end do
c         end do
c         end do
c       *************************************************
c
c     !!!!!! for B field check (for bcheck5.dat)
c      do m=1,ngrd
c      write(55,65)bx0(xloc,yloc,zloc,m),by0(xloc,yloc,zloc,m),
c     +bz0(xloc,yloc,zloc,m),bx(xloc,yloc,zloc,m),by(xloc,yloc,zloc,m),
c     +bz(xloc,yloc,zloc,m),ut,delt,t_equiv,m
c65    format(10(1x,1pe9.2))
c      enddo
c     !!!!!!
c
c     .......................................................
c     try Lapdius smoothing - smoothed results will appear in nt2
c     .......................................................
c
      write(6,*)' doing smoothing okay'
      do 420 m=1,ngrd
       rx=xspac(m)
       ry=xspac(m)
       rz=xspac(m)
c      write(6,*)'calling lapidus'
      call lapidus(qrho,qpres,qpx,qpy,qpz,
     +        wrkqrho,wrkqpres,wrkqpx,wrkqpy,wrkqpz,
     +        hrho,hpres,hpx,hpy,hpz,
     +        wrkhrho,wrkhpres,wrkhpx,wrkhpy,wrkhpz,
     +        orho,opres,opx,opy,opz,
     +        wrkorho,wrkopres,wrkopx,wrkopy,wrkopz,
     +        epres,wrkepres,bx,by,bz,
     +        wrkbx,wrkby,wrkbz,curx,cury,curz,
     +        nx,ny,nz,ngrd,m,chirho,chipxyz,chierg,delt,
     +        rmassq,rmassh,rmasso)
  420 continue
c
c
c     reset boundary conditions
c
      call bndry(wrkqrho,wrkqpres,wrkqpx,wrkqpy,wrkqpz,rmassq,
     +        wrkhrho,wrkhpres,wrkhpx,wrkhpy,wrkhpz,rmassh,
     +        wrkorho,wrkopres,wrkopx,wrkopy,wrkopz,rmasso,
     +        wrkepres,wrkbx,wrkby,wrkbz,bx0,by0,bz0,
     +        nx,ny,nz,ngrd,parm_srf,parm_mid,
     +        ijsrf,numsrf,ijmid,nummid,
     +        ijzero,numzero,erho,epress,alpha_e,o_conc,
     +        ti_te,srho,rho_frac,spress,spx,spy,spz,
     +        sbx_wind,sby_wind,sbz_wind)
c
c     write(6,*)'Lapidus speeds'
c
      call set_speed(wrkqrho,wrkqpres,wrkqpx,wrkqpy,wrkqpz,
     +    wrkhrho,wrkhpres,wrkhpx,wrkhpy,wrkhpz,
     +    wrkorho,wrkopres,wrkopx,wrkopy,wrkopz,
     +    wrkepres,wrkbx,wrkby,wrkbz,bx0,by0,bz0,
     +    bsx,bsy,bsz,btot,rmassq,rmassh,rmasso,nx,ny,nz,ngrd,
     +    pxmax,pymax,pzmax,pmax,csmax,alfmax,gamma,
     +    vlim,alf_lim,o_conc,fastest)
c
      if(spacecraft)then
         call set_imf(wrkbx,wrkby,wrkbz,bx0,by0,bz0,bxp,byp,bzp,
     +         wrkqrho,wrkqpres,wrkqpx,wrkqpy,wrkqpz,wrkepres,
     +         rhop,svxp,svyp,svzp,svelx,spress,ti_te,
     +         rho_frac,nx,ny,nz,ngrd)
      endif
c
c     write(6,994)
c 994 format(' lapidus done')
c
c
c     .......................................................
c     add a little bit of artifical diffusion : 
c     ........................................................
c
c      first subtract out unperturbed desnity so that
c           not diffusing large densities
c
c      chifcs=.125
      chifcs=.0625
c
      do 460 m=1,ngrd
c
       rx=xspac(m)
       ry=xspac(m)
       rz=xspac(m)
       call fcsmooth(qrho,oldqrho,wrkqrho,ngrd,nx,ny,nz,m,chifcs,
     +         vvx,vvy,vvz)
       call fcsmooth(qpres,oldqpres,wrkqpres,ngrd,nx,ny,nz,m,chifcs,
     +         vvx,vvy,vvz)
       call fcsmooth(qpx,oldqpx,wrkqpx,ngrd,nx,ny,nz,m,chifcs,
     +         vvx,vvy,vvz)
       call fcsmooth(qpy,oldqpy,wrkqpy,ngrd,nx,ny,nz,m,chifcs,
     +         vvx,vvy,vvz)
       call fcsmooth(qpz,oldqpz,wrkqpz,ngrd,nx,ny,nz,m,chifcs,
     +         vvx,vvy,vvz)
c
       call fcsmooth(hrho,oldhrho,wrkhrho,ngrd,nx,ny,nz,m,chifcs,
     +         vvx,vvy,vvz)
       call fcsmooth(hpres,oldhpres,wrkhpres,ngrd,nx,ny,nz,m,chifcs,
     +         vvx,vvy,vvz)
       call fcsmooth(hpx,oldhpx,wrkhpx,ngrd,nx,ny,nz,m,chifcs,
     +         vvx,vvy,vvz)
       call fcsmooth(hpy,oldhpy,wrkhpy,ngrd,nx,ny,nz,m,chifcs,
     +         vvx,vvy,vvz)
       call fcsmooth(hpz,oldhpz,wrkhpz,ngrd,nx,ny,nz,m,chifcs,
     +         vvx,vvy,vvz)
c
       call fcsmooth(orho,oldorho,wrkorho,ngrd,nx,ny,nz,m,chifcs,
     +         vvx,vvy,vvz)
       call fcsmooth(opres,oldopres,wrkopres,ngrd,nx,ny,nz,m,chifcs,
     +         vvx,vvy,vvz)
       call fcsmooth(opx,oldopx,wrkopx,ngrd,nx,ny,nz,m,chifcs,
     +         vvx,vvy,vvz)
       call fcsmooth(opy,oldopy,wrkopy,ngrd,nx,ny,nz,m,chifcs,
     +         vvx,vvy,vvz)
       call fcsmooth(opz,oldopz,wrkopz,ngrd,nx,ny,nz,m,chifcs,
     +         vvx,vvy,vvz)
c
       call fcsmooth(epres,oldepres,wrkepres,ngrd,nx,ny,nz,m,chifcs,
     +         vvx,vvy,vvz)
c
       call fcsmooth(bx,oldbx,wrkbx,ngrd,nx,ny,nz,m,chifcs,
     +         vvx,vvy,vvz)
       call fcsmooth(by,oldby,wrkby,ngrd,nx,ny,nz,m,chifcs,
     +         vvx,vvy,vvz)
       call fcsmooth(bz,oldbz,wrkbz,ngrd,nx,ny,nz,m,chifcs,
     +         vvx,vvy,vvz)
c      bx=wrkbx
c      by=wrkby
c      bz=wrkbz
c
c      call ssmooth(bx,oldbx,wrkbx,ngrd,nx,ny,nz,m,difrho)
c      call ssmooth(by,oldby,wrkby,ngrd,nx,ny,nz,m,difrho)
c      call ssmooth(bz,oldbz,wrkbz,ngrd,nx,ny,nz,m,difrho)
c
c     raw result is in m
c      !!!!!!
      if(Adiffuse)then
      do 440 k=1,nz
      do 440 j=1,ny
      do 440 i=1,nx
        wrkorho(i,j,k,m)=wrkorho(i,j,k,m)-orho0(i,j,k,m)
        wrkhrho(i,j,k,m)=wrkhrho(i,j,k,m)-hrho0(i,j,k,m)
        wrkqrho(i,j,k,m)=wrkqrho(i,j,k,m)-qrho0(i,j,k,m)
c
        wrkqpres(i,j,k,m)=wrkqpres(i,j,k,m)-qpres0(i,j,k,m)
        wrkhpres(i,j,k,m)=wrkhpres(i,j,k,m)-hpres0(i,j,k,m)
        wrkopres(i,j,k,m)=wrkopres(i,j,k,m)-opres0(i,j,k,m)
        wrkepres(i,j,k,m)=wrkepres(i,j,k,m)-epres0(i,j,k,m) 
  440 continue
c
      call diffuse(qrho,qpres,qpx,qpy,qpz,
     +      hrho,hpres,hpx,hpy,hpz,orho,opres,opx,opy,opz,
     +      epres,bx,by,bz,wrkqrho,wrkqpres,wrkqpx,wrkqpy,wrkqpz,
     +      wrkhrho,wrkhpres,wrkhpx,wrkhpy,wrkhpz,wrkorho,
     +      wrkopres,wrkopx,wrkopy,wrkopz,
     +      wrkepres,wrkbx,wrkby,wrkbz,nx,ny,nz,ngrd,m,
     +      difrho,difpxyz,diferg,1.,rx,ry,rz)
c
c      add back in the unperturbed density : smoothed 
c
      do 450 k=1,nz
      do 450 j=1,ny
      do 450 i=1,nx
        orho(i,j,k,m)=orho(i,j,k,m)+orho0(i,j,k,m)
        hrho(i,j,k,m)=hrho(i,j,k,m)+hrho0(i,j,k,m)
        qrho(i,j,k,m)=qrho(i,j,k,m)+qrho0(i,j,k,m)
c
        qpres(i,j,k,m)=qpres(i,j,k,m)+qpres0(i,j,k,m)
        hpres(i,j,k,m)=hpres(i,j,k,m)+hpres0(i,j,k,m)
        opres(i,j,k,m)=opres(i,j,k,m)+opres0(i,j,k,m)
        epres(i,j,k,m)=epres(i,j,k,m)+epres0(i,j,k,m)  
  450 continue
      endif
c      !!!!!!
  460 continue
c
c     reset boundary conditions
c
c
c     bx=wrkbx
c     by=wrkby
c     bz=wrkbz
      call bndry(qrho,qpres,qpx,qpy,qpz,rmassq,
     +        hrho,hpres,hpx,hpy,hpz,rmassh,
     +        orho,opres,opx,opy,opz,rmasso,
     +        epres,bx,by,bz,bx0,by0,bz0,
     +        nx,ny,nz,ngrd,parm_srf,parm_mid,
     +        ijsrf,numsrf,ijmid,nummid,
     +        ijzero,numzero,erho,epress,alpha_e,o_conc,
     +         ti_te,srho,rho_frac,spress,spx,spy,spz,
     +        sbx_wind,sby_wind,sbz_wind)
      if(spacecraft)then
       call set_imf(bx,by,bz,bx0,by0,bz0,bxp,byp,bzp,
     +             qrho,qpres,qpx,qpy,qpz,epres,
     +             rhop,svxp,svyp,svzp,svelx,spress,ti_te,
     +             rho_frac,nx,ny,nz,ngrd)
      endif
c
c     write(6,*)'Diffusion speeds'
      call set_speed(qrho,qpres,qpx,qpy,qpz,
     +    hrho,hpres,hpx,hpy,hpz,orho,opres,opx,opy,opz,
     +    epres,bx,by,bz,bx0,by0,bz0,bsx,bsy,bsz,btot,
     +    rmassq,rmassh,rmasso,nx,ny,nz,ngrd,
     +    pxmax,pymax,pzmax,pmax,csmax,alfmax,gamma,
     +    vlim,alf_lim,o_conc,fastest)
c     write(6,996)
c 996 format(' diffusion applied')
c
c     write(6,999)t
c 999 format(' step 2 complete at t= ',1pe12.5)
c
c      !!!!!! We add the inner bndry rotation with time & surface 
c       call bndry_inner_rot(qrho,qpres,rmassq,
c     +       hrho,hpres,rmassh,
c     +       orho,opres,rmasso,epres,gamma1,eerg,
c     +       ngrd,nx,ny,nz,parm_srf,parm_mid,
c     +       ijsrf,numsrf,ijmid,nummid,
c     +       erho,ti_te,o_conc)
c



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      call bndry_inner_rot(rmassq,rmassh,rmasso,
     +       gamma1,eerg,ngrd,nx,ny,nz,parm_srf,parm_mid,
     +       ijsrf,numsrf,ijmid,nummid,erho,ti_te,o_conc)
c
      call rigid_IM_rot(rmassq,rmassh,rmasso,
      +       gamma1,eerg,ngrd,nx,ny,nz,qrho,hrho,orho,
      +       erho,qpres,hpres,opres,epres,ti_te,o_conc,ut,
      +       utdec_start,re_equiv,bx0,by0,bz0)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	  
	  
	  
	  
	  
c      !!!!!! test if parm_srf is changing (should be!)
       print *, '11','hrho_srf=', parm_srf(2,2), 'hpres=', parm_srf(5,2)
       print *, 'qrho_srf=', parm_srf(1,2), 'epres=', parm_srf(7,2)
c      !!!!!!
c
c
      if(t.lt.tgraf)goto 600
c
c     plot plasma propeties
c
c      calculate size of plotting stuff and ensure no distortions
c         over desired scales
c
       do 520 m=ngrd,1,-1
        rx=xspac(m)
        ry=xspac(m)
        rz=xspac(m)
c
c      dm=(grd_zmax(m)-grd_zmin(m))/2.+zdip-2.
c      xtot=1.3*(2.*dm)
c      xdiff=grd_xmax(m)-grd_xmin(m)
c      xmin=3.*grd_xmin(m)/4.+2.
c      xmax=xmin+xtot
c
       ymin=grd_ymin(m)+ry
       ymax=grd_ymax(m)-ry
       zmin=grd_zmin(m)+rz
       zmax=grd_zmax(m)-rz
       xmin=grd_xmin(m)+rx
       xmax=grd_xmax(m)+rx
c
c      ycent=(grd_ymax(m)+grd_ymin(m))/2.
c      ydist=(zmax-zmin)/2.
c      ymax=ycent+ydist
c      ymin=ycent-ydist
c
c       xcut=15./re_equiv
       xcut=xmax/2.
c
       add_dip=.false.
       av=sqrt(vx_wind1)
       gs=1./sqrt(gamma)
       do 505 k=1,nz
       do 505 j=1,ny
       do 505 i=1,nx
        efldx(i,j,k)=sqrt(qpres(i,j,k,m))
        efldy(i,j,k)=sqrt(hpres(i,j,k,m))
        efldz(i,j,k)=sqrt(opres(i,j,k,m))
  505 continue
c
      call qvset(0.,curx,nx*ny*nz)
      call qvset(0.,cury,nx*ny*nz)
      call qvset(0.,curz,nx*ny*nz)
c
       preslim=0.4
       po=preslim
       call fnd_vel(qpx,qpy,qpz,qrho,curx,cury,curz,ngrd,nx,ny,nz,m)
       call conhot(efldx,curx,cury,curz,nx,ny,nz,1,1,m,
     +            xmin,xmax,ymin,ymax,zmin,zmax,xcut,
     +                      ut,'qpres',3,11,1,2.0,preslim,
     +             tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2)
       call fnd_vel(hpx,hpy,hpz,hrho,curx,cury,curz,ngrd,nx,ny,nz,m)
       call conhot(efldy,curx,cury,curz,nx,ny,nz,1,1,m,
     +            xmin,xmax,ymin,ymax,zmin,zmax,xcut,
     +                      ut,'hpres',3,11,1,2.0,preslim,
     +             tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2)
       call fnd_vel(opx,opy,opz,orho,curx,cury,curz,ngrd,nx,ny,nz,m)
       call conhot(efldz,curx,cury,curz,nx,ny,nz,1,1,m,
     +            xmin,xmax,ymin,ymax,zmin,zmax,xcut,
     +                      ut,'opres',3,11,1,2.0,po,
     +             tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2)


       do 510 k=1,nz
       do 510 j=1,ny
       do 510 i=1,nx
        efldx(i,j,k)=sqrt(abs(qpres(i,j,k,m)+hpres(i,j,k,m)
     +                   +opres(i,j,k,m)+epres(i,j,k,m)))
        efldy(i,j,k)=opres(i,j,k,m)/(hpres(i,j,k,m)+
     +         opres(i,j,k,m)+qpres(i,j,k,m)+0.00001)
        efldz(i,j,k)=epres(i,j,k,m)/(epres(i,j,k,m)
     +         +hpres(i,j,k,m)+opres(i,j,k,m)
     +         +qpres(i,j,k,m)+0.00001)
  510 continue
c      call contur(efldx,nx,ny,nz,1,1,xmin,xmax,
c    +       ymin,ymax,zmin,zmax,t,'sqrt pres',3,11,
c            tx,ty,tg2,work,mx,my,mz,mz2,muvwp2)
c
c      trace velocity streams
c
       call fnd_vtot(qpx,qpy,qpz,qrho,hpx,hpy,hpz,hrho, 
     +       opx,opy,opz,orho,curx,cury,curz,ngrd,nx,ny,nz,m,
     +      rmassq,rmassh,rmasso)
c
       call conflow(efldx,curx,cury,curz,nx,ny,nz,1,1,m,
     +            xmin,xmax,ymin,ymax,zmin,zmax,xcut,
     +             ut,'pres-vel',3,11,1,2.0,
     +             tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2)
c
c
c     find total magnetic field
c
      call totfld(bx,bx0,bsx,nx,ny,nz,ngrd,m)
      call totfld(by,by0,bsy,nx,ny,nz,ngrd,m)
      call totfld(bz,bz0,bsz,nx,ny,nz,ngrd,m)
c
      call qvset(0.,curx,nx*ny*nz)
      call qvset(0.,cury,nx*ny*nz)
      call qvset(0.,curz,nx*ny*nz)
c

       write(wd1,'(i3)')m
       label='box '//wd1
      call concross(efldx,curx,cury,curz,bsx,bsy,bsz,
     +            nx,ny,nz,1,1,m,xcraft,ncraft,re_equiv,
     +            xmin,xmax,ymin,ymax,zmin,zmax,
     +             ut,label,3,11,2.0,add_dip,1,-2,start,
     +       tx,ty,tz,tt,tg1,tg2,cross,along,work,mx,my,mz,mz2,muvwp2)

      call contop(efldx,curx,cury,curz,bsx,bsy,bsz,
     +            nx,ny,nz,1,1,m,xcraft,ncraft,re_equiv,
     +            xmin,xmax,ymin,ymax,zmin,zmax,
     +             ut,label,3,11,2.0,add_dip,1,0,start,
     +       tx,ty,tz,tt,tg1,tg2,cross,along,work,mx,my,mz,mz2,muvwp2)
c
c       magnetidue of B
c

       blim=0.
       do  k=1,nz
       do  j=1,ny
       do  i=1,nx
        efldx(i,j,k)=sqrt(bsx(i,j,k)**2+bsy(i,j,k)**2
     +                    +bsz(i,j,k)**2)
        efldx(i,j,k)=alog10(b_equiv*efldx(i,j,k)+1.e-10)
        blim=amax1(blim,efldx(i,j,k))
       enddo
       enddo
       enddo
c
       call conlog(efldx,curx,cury,curz,nx,ny,nz,1,1,m,
     +            xmin,xmax,ymin,ymax,zmin,zmax,xcut,
     +              ut,'bmag',3,14,1,2.0,0.1,4.,
     +             tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2)
c
c
       blim=0.001
       do  k=1,nz
       do  j=1,ny
       do  i=1,nx
        efldx(i,j,k)=bsz(i,j,k)
        if (efldx(i,j,k).gt.blim)efldx(i,j,k)=blim
        if (efldx(i,j,k).lt.-blim)efldx(i,j,k)=-blim
        efldx(i,j,k)=efldx(i,j,k)+1.001*blim
       enddo
       enddo
       enddo
        call conhot(efldx,curx,cury,curz,nx,ny,nz,1,1,m,
     +            xmin,xmax,ymin,ymax,zmin,zmax,xcut,
     +                   ut,'bz',3,12,1,2.0,2.*blim,
     +             tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2)
c
c       additional cross-plane cuts
c
c       call concraft_fix(efldx,cross,nx,ny,nz,m,
c    +       xcraft,ncraft,re_equiv,ut,'prs yz',start,0.,.25)
c       call conalong_fix(efldx,along,nx,ny,nz,m,
c    +       xcraft,ncraft,re_equiv,ut,'prs xz',start,0.,.33) 
c
c
c      plot individual temperatures
c
       tempx=1.0
       tempo=tempx*2.
       if(m.gt.1) then
          rho_lim=10.
       else
          rho_lim=5.0
       endif
       rfrac=0.5
c
       rho_tot=1.+rho_frac
       do  k=1,nz
       do  j=1,ny
       do  i=1,nx
        qden=qrho(i,j,k,m)/rmassq
        if(qden.gt.0.001)then
          bsx(i,j,k)=amin1(tempx,sqrt(qpres(i,j,k,m)/qden))
        else 
          bsx(i,j,k)=0.
        endif
        hden=hrho(i,j,k,m)/rmassh
        if(hden.gt.0.0005)then
          bsy(i,j,k)=amin1(tempx,sqrt(hpres(i,j,k,m)/hden))
        else 
          bsy(i,j,k)=0.
        endif
        oden=orho(i,j,k,m)/rmasso
        if(oden.gt.0.00002)then
          bsz(i,j,k)=amin1(tempo,sqrt(opres(i,j,k,m)/oden))
        else 
          bsz(i,j,k)=0.
        endif
c
        efldx(i,j,k)=sqrt(epres(i,j,k,m))
        tden=oden+qden+hden
        if(tden.gt.0.001)then
          efldy(i,j,k)=amin1(tempx,sqrt(epres(i,j,k,m)/tden))
        else 
          efldy(i,j,k)=0.
        endif
       enddo
       enddo
       enddo
c
       call conhot(efldx,curx,cury,curz,nx,ny,nz,1,1,m,
     +            xmin,xmax,ymin,ymax,zmin,zmax,xcut,
     +                      ut,'epres',3,11,1,2.0,po,
     +             tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2)
c

c
        call conhot(bsx,curx,cury,curz,nx,ny,nz,1,1,m,
     +            xmin,xmax,ymin,ymax,zmin,zmax,xcut,
     +                   ut,'q temp',3,12,1,2.0,tempx,
     +             tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2)
        call conhot(bsy,curx,cury,curz,nx,ny,nz,1,1,m,
     +            xmin,xmax,ymin,ymax,zmin,zmax,xcut,
     +                    ut,'h temp',3,12,1,2.0,tempx,
     +             tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2)
        call conhot(bsz,curx,cury,curz,nx,ny,nz,1,1,m,
     +            xmin,xmax,ymin,ymax,zmin,zmax,xcut,
     +                    ut,'o temp',3,12,1,2.0,tempo,
     +             tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2)
        call conhot(efldy,curx,cury,curz,nx,ny,nz,1,1,m,
     +            xmin,xmax,ymin,ymax,zmin,zmax,xcut,
     +                    ut,'e temp',3,12,1,2.0,
     +                     tempx/sqrt(ti_te),
     +             tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2)
c
c       plot mass flows : solar wind
c
       call qvset(0.,curx,nx*ny*nz)
       call qvset(0.,cury,nx*ny*nz)
       call qvset(0.,curz,nx*ny*nz)
c       call fnd_vel(qpx,qpy,qpz,qrho,curx,cury,curz,ngrd,nx,ny,nz,m)
       do 516 k=1,nz
       do 516 j=1,ny
       do 516 i=1,nx
        efldx(i,j,k)=(qrho(i,j,k,m)/rmassq)
        efldx(i,j,k)=alog10(efldx(i,j,k)*rho_equiv)+6. ! per m**3
c        if(qrho(i,j,k,m).lt.0.01)then
c          curx(i,j,k)=0.
c          cury(i,j,k)=0.
c          curz(i,j,k)=0.
c        endif
 516    continue
c
       call conlog(efldx,curx,cury,curz,nx,ny,nz,1,1,m,
     +            xmin,xmax,ymin,ymax,zmin,zmax,xcut,
     +              ut,'q den v',3,14,1,2.0,4.,8.,
     +             tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2)
c
c       plot mass flows : ionospheric H
c
c      call fnd_vel(hpx,hpy,hpz,hrho,curx,cury,curz,ngrd,nx,ny,nz,m)
c
       do  k=1,nz
       do  j=1,ny
       do  i=1,nx
        efldy(i,j,k)=hrho(i,j,k,m)/rmassh
        efldy(i,j,k)=alog10(efldy(i,j,k)*rho_equiv)+6.
c       if(hrho(i,j,k,m).lt.0.003)then
c         curx(i,j,k)=0.
c         cury(i,j,k)=0.
c         curz(i,j,k)=0.
c       endif
       enddo
       enddo
       enddo
c
       call conlog(efldy,curx,cury,curz,nx,ny,nz,1,1,m,
     +            xmin,xmax,ymin,ymax,zmin,zmax,xcut,
     +                  ut,'h den v',3,14,1,2.0,4.,8.,
     +             tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2)
c
c       plot mass flows : ionospheric O
c
c       call fnd_vel(opx,opy,opz,orho,curx,cury,curz,ngrd,nx,ny,nz,m)
c
       do  k=1,nz
       do  j=1,ny
       do  i=1,nx
        efldz(i,j,k)=orho(i,j,k,m)/rmasso
        efldz(i,j,k)=alog10(efldz(i,j,k)*rho_equiv)+6.
c        if(orho(i,j,k,m).lt.0.001)then
c          curx(i,j,k)=0.
c          cury(i,j,k)=0.
c          curz(i,j,k)=0.
c        endif
       enddo
       enddo
       enddo
c
       call conlog(efldz,curx,cury,curz,nx,ny,nz,1,1,m,
     +            xmin,xmax,ymin,ymax,zmin,zmax,xcut,
     +             ut,'o den v',3,14,1,2.0,4.,8.,
     +             tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2)
c
c       plot mass flows : total
c
       call fnd_vtot(qpx,qpy,qpz,qrho,hpx,hpy,hpz,hrho,
     +       opx,opy,opz,orho,curx,cury,curz,ngrd,nx,ny,nz,m,
     +      rmassq,rmassh,rmasso)
       do  k=1,nz
       do  j=1,ny
       do  i=1,nx
        efldx(i,j,k)=(orho(i,j,k,m)/rmasso+
     +          hrho(i,j,k,m)/rmassh+qrho(i,j,k,m)/rmassq)
        efldx(i,j,k)=alog10(efldx(i,j,k)*rho_equiv)+6.
        if(efldx(i,j,k).lt.0.1)then
          curx(i,j,k)=0.
          cury(i,j,k)=0.
          curz(i,j,k)=0.
        endif
       enddo
       enddo
       enddo
c
       call conlog(efldx,curx,cury,curz,nx,ny,nz,1,1,m,
     +            xmin,xmax,ymin,ymax,zmin,zmax,xcut,
     +                  ut,'tden v',3,14,1,2.0,5.,8.5,
     +             tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2)
c       call concraft_fix(efldx,cross,nx,ny,nz,m,
c    +       xcraft,ncraft,re_equiv,ut,'tden yz',start,0.,2.2)
c       call conalong_fix(efldx,along,nx,ny,nz,m,
c    +       xcraft,ncraft,re_equiv,ut,'tden xz',start,0.,2.2) 
c       call conplane_fix(efldx,along,nx,ny,nz,m,
c    +       xcraft,ncraft,re_equiv,ut,'tden xy',start,0.,2.2) 
c
  520 continue

c
c     mapping closed field lines using grid 3
c
       add_dip=.false.
       radstrt=rearth+3.
c
       m=3
       rx=xspac(m)
       ry=xspac(m)
       rz=xspac(m)
       radstrt=rearth+6.
c
c     find total magnetic field
c
      call totfld(bx,bx0,bsx,nx,ny,nz,ngrd,m)
      call totfld(by,by0,bsy,nx,ny,nz,ngrd,m)
      call totfld(bz,bz0,bsz,nx,ny,nz,ngrd,m)
c
       xmin=grd_xmin(m)+rx
       xmax=grd_xmax(m)-rx
       ymin=grd_ymin(m)+ry
       ymax=grd_ymax(m)-ry
       zmin=grd_zmin(m)+rz
       zmax=grd_zmax(m)-rz
c
       call aurora_bfld(bsx,bsy,bsz,nx,ny,nz,m,
     +            xmin,xmax,ymin,ymax,zmin,zmax,1,
     +            add_dip,radstrt,re_equiv,ut,'big nth')
       call aurora_bfld(bsx,bsy,bsz,nx,ny,nz,m,
     +            xmin,xmax,ymin,ymax,zmin,zmax,-1,
     +            add_dip,radstrt,re_equiv,ut,'big sth')
c
c     mapping closed field lines using grid 2
c
       m=2
       rx=xspac(m)
       ry=xspac(m)
       rz=xspac(m)
       radstrt=rearth+6.
c
c     find total magnetic field
c
      call totfld(bx,bx0,bsx,nx,ny,nz,ngrd,m)
      call totfld(by,by0,bsy,nx,ny,nz,ngrd,m)
      call totfld(bz,bz0,bsz,nx,ny,nz,ngrd,m)
c
       xmin=grd_xmin(m)+rx
       xmax=grd_xmax(m)-rx
       ymin=grd_ymin(m)+ry
       ymax=grd_ymax(m)-ry
       zmin=grd_zmin(m)+rz
       zmax=grd_zmax(m)-rz
c
       write_dat=.true.
       save_dat=.false.
       add_dip=.false.
c
c
       call aurora_bfld(bsx,bsy,bsz,nx,ny,nz,m,
     +            xmin,xmax,ymin,ymax,zmin,zmax,1,
     +            add_dip,radstrt,re_equiv,ut,'crse nth')
       call aurora_bfld(bsx,bsy,bsz,nx,ny,nz,m,
     +            xmin,xmax,ymin,ymax,zmin,zmax,-1,
     +            add_dip,radstrt,re_equiv,ut,'crse sth')
c
      call calcur(bx,by,bz,ngrd,nx,ny,nz,m,curx,cury,curz,
     +            ijsrf,numsrf,ijmid,nummid,ijzero,numzero)
c
c     find the surface current magnetic field and convective
c         electric field at the same time
c
       call fnd_evel(qpx,qpy,qpz,qrho,hpx,hpy,hpz,hrho, 
     +       opx,opy,opz,orho,curx,cury,curz,evx,evy,evz,
     +       tvx,tvy,tvz,ngrd,
     +       nx,ny,nz,m,rmassq,rmassh,rmasso,reynolds)
      call bande(efldx,efldy,efldz,bsx,bsy,bsz,
     +       curx,cury,curz,evx,evy,evz,btot,
     +       epres,qrho,hrho,orho,resistive,resist,reynolds,
     +       nx,ny,nz,ngrd,m,rmassq,rmassh,rmasso,
     +       ijmid,nummid,ijzero,numzero)
c
c
c     pot_scale=v_equiv*b_equiv*rx*re_equiv*6371.*1.e-9 ! units in MV 1e-6
c               1.e3      1.e-9          1.e3
c     calculate unscaled potential as first estimate
c
c     calculate the space charge density to be used later in
c        making potential plots
c
c     need m=2 to remove potential boundary conditions problems
c
      call space_charge(bsz,efldx,efldy,efldz,
     +         bx0,by0,bz0,opx,opy,opz,orho,nx,ny,nz,ngrd,m)
c
c      plot these guys out
c
      call cappot(bsz,bsx,nx,ny,nz,ngrd,m,
     +     radstrt,re_equiv,v_equiv,b_equiv,ut,write_dat)

c
c     aurora diagnostics using inner grid system
c
       m=1
       rx=xspac(m)
       ry=xspac(m)
       rz=xspac(m)
       radstrt=rearth+6.
c
       xmin=grd_xmin(m)+rx
       xmax=grd_xmax(m)-rx
       ymin=grd_ymin(m)+ry
       ymax=grd_ymax(m)-ry
       zmin=grd_zmin(m)+rz
       zmax=grd_zmax(m)-rz
c
c     find total magnetic field
c
      call totfld(bx,bx0,bsx,nx,ny,nz,ngrd,m)
      call totfld(by,by0,bsy,nx,ny,nz,ngrd,m)
      call totfld(bz,bz0,bsz,nx,ny,nz,ngrd,m)
c
c     find magnitude of B
c
      call tot_b(btot,bsx,bsy,bsz,nx,ny,nz)
c
c     mark the alleged position of open and closed field lines
c
       add_dip=.false.
       call aurora_bfld(bsx,bsy,bsz,nx,ny,nz,m,
     +            xmin,xmax,ymin,ymax,zmin,zmax,1,
     +            add_dip,radstrt,re_equiv,ut,'fine nth')
       call aurora_bfld(bsx,bsy,bsz,nx,ny,nz,m,
     +            xmin,xmax,ymin,ymax,zmin,zmax,-1,
     +            add_dip,radstrt,re_equiv,ut,'fine sth')
c
c      find outflow velcoities
c
       call fnd_vtot(qpx,qpy,qpz,qrho,hpx,hpy,hpz,hrho,
     +       opx,opy,opz,orho,curx,cury,curz,ngrd,nx,ny,nz,m,
     +      rmassq,rmassh,rmasso)
       write_dat=.true.
       save_dat=.false.
       add_dip=.false.
c
c      plot convection patterns and potential
c
       do nhi=4,8,4
       radstrt=rearth+nhi
       write(wd1,'(f4.1)')radstrt
       label='flwn'//wd1
       call convect(curx,cury,curz,nx,ny,nz,m,radstrt,
     +                  re_equiv,1,ut,label,write_dat)
       label='flws'//wd1
       call convect(curx,cury,curz,nx,ny,nz,m,radstrt,
     +                  re_equiv,-1,ut,label,write_dat)
       enddo
c
      call calcur(bx,by,bz,ngrd,nx,ny,nz,m,curx,cury,curz,
     +            ijsrf,numsrf,ijmid,nummid,ijzero,numzero)
c
c     find the surface current magnetic field and convective
c         electric field at the same time
c

      call fnd_evel(qpx,qpy,qpz,qrho,hpx,hpy,hpz,hrho, 
     +       opx,opy,opz,orho,curx,cury,curz,evx,evy,evz,
     +       tvx,tvy,tvz,ngrd,
     +       nx,ny,nz,m,rmassq,rmassh,rmasso,reynolds)
      call bande(efldx,efldy,efldz,bsx,bsy,bsz,
     +       curx,cury,curz,evx,evy,evz,btot,
     +       epres,qrho,hrho,orho,resistive,resist,reynolds,
     +       nx,ny,nz,ngrd,m,rmassq,rmassh,rmasso,
     +       ijmid,nummid,ijzero,numzero)
c
c
c      field-aligned potential drop
c
       call aurora_pot(efldx,efldy,efldz,bsx,bsy,bsz,
     +            nx,ny,nz,ngrd,m,radstrt,re_equiv,1,
     +            ut,save_dat,add_dip,'ppot nth',11,write_dat)
       call aurora_pot(efldx,efldy,efldz,bsx,bsy,bsz,
     +            nx,ny,nz,ngrd,m,radstrt,re_equiv,-1,
     +            ut,save_dat,add_dip,'ppot sth',11,write_dat)
c
c        field aligned currents
c
       do k=1,nz
       do j=1,ny
       do i=1,nx
         abx=bx(i,j,k,m)+bx0(i,j,k,m)
         aby=by(i,j,k,m)+by0(i,j,k,m)
         abz=bz(i,j,k,m)+bz0(i,j,k,m)
         bmag=sqrt(abx**2+aby**2+abz**2)
     +                 +0.00000001
         btot(i,j,k)=(curx(i,j,k)*abx+cury(i,j,k)*aby
     +                         +curz(i,j,k)*abz)/bmag
      enddo
      enddo
      enddo
c
       save_dat=.false.
       add_dip=.false.
c
       call aurora(btot,nx,ny,nz,m,radstrt,re_equiv,1,
     +          ut,save_dat,add_dip,'cur4 nth',14,write_dat)
       call aurora(btot,nx,ny,nz,m,radstrt,re_equiv,-1,
     +          ut,save_dat,add_dip,'cur4 sth',14,write_dat)
c
c        field aligned flows
c
       do 540 k=1,nz
       do 540 j=1,ny
       do 540 i=1,nx
         abx=bx(i,j,k,m)+bx0(i,j,k,m)
         aby=by(i,j,k,m)+by0(i,j,k,m)
         abz=bz(i,j,k,m)+bz0(i,j,k,m)
         bmag=sqrt(abx**2+aby**2+abz**2)
     +                 +0.00000001
         aden=qrho(i,j,k,m)/rmassq+hrho(i,j,k,m)/rmassh
     +       +orho(i,j,k,m)/rmasso
         efldx(i,j,k)=(qpx(i,j,k,m)*abx+qpy(i,j,k,m)*aby
     +                         +qpz(i,j,k,m)*abz)/bmag
         efldy(i,j,k)=(hpx(i,j,k,m)*abx+hpy(i,j,k,m)*aby
     +                         +hpz(i,j,k,m)*abz)/bmag
         efldz(i,j,k)=(opx(i,j,k,m)*abx+opy(i,j,k,m)*aby
     +                         +opz(i,j,k,m)*abz)/bmag
         epx=qpx(i,j,k,m)/rmassq+hpx(i,j,k,m)/rmassh
     +           +opx(i,j,k,m)/rmasso-curx(i,j,k)/reynolds
         epy=qpy(i,j,k,m)/rmassq+hpy(i,j,k,m)/rmassh
     +           +opy(i,j,k,m)/rmasso-cury(i,j,k)/reynolds
         epz=qpz(i,j,k,m)/rmassq+hpz(i,j,k,m)/rmassh
     +           +opz(i,j,k,m)/rmasso-curz(i,j,k)/reynolds
         curz(i,j,k)=(epx*abx+epy*aby+epz*abz)/bmag
  540 continue
c
       save_dat=.false.
       add_dip=.false.
c
       call aurora(efldx,nx,ny,nz,m,radstrt,re_equiv,1,
     +          ut,save_dat,add_dip,'qfluxnth',14,write_dat)
       call aurora(efldx,nx,ny,nz,m,radstrt,re_equiv,-1,
     +          ut,save_dat,add_dip,'qfluxsth',14,write_dat)
c
       call aurora(efldy,nx,ny,nz,m,radstrt,re_equiv,1,
     +          ut,save_dat,add_dip,'hfluxnth',14,write_dat)
       call aurora(efldy,nx,ny,nz,m,radstrt,re_equiv,-1,
     +          ut,save_dat,add_dip,'hfluxsth',14,write_dat)
c
       call aurora(efldz,nx,ny,nz,m,radstrt,re_equiv,1,
     +          ut,save_dat,add_dip,'ofluxnth',14,write_dat)
       call aurora(efldz,nx,ny,nz,m,radstrt,re_equiv,-1,
     +          ut,save_dat,add_dip,'ofluxsth',14,write_dat)
c
       call aurora(curz,nx,ny,nz,m,radstrt,re_equiv,1,
     +          ut,save_dat,add_dip,'efluxnth',14,write_dat)
       call aurora(curz,nx,ny,nz,m,radstrt,re_equiv,-1,
     +          ut,save_dat,add_dip,'efluxsth',14,write_dat)
c
       do 550 k=1,nz
       do 550 j=1,ny
       do 550 i=1,nx
         abx=bx(i,j,k,m)+bx0(i,j,k,m)
         aby=by(i,j,k,m)+by0(i,j,k,m)
         abz=bz(i,j,k,m)+bz0(i,j,k,m)
         bmag=sqrt(abx**2+aby**2+abz**2)
     +                 +0.00000001
         aqrho=amax1(0.000001,qrho(i,j,k,m))
         ahrho=amax1(0.000001,hrho(i,j,k,m))
         aorho=amax1(0.000001,orho(i,j,k,m))
c        
         efldx(i,j,k)=(qpx(i,j,k,m)*abx+qpy(i,j,k,m)*aby
     +                 +qpz(i,j,k,m)*abz)/bmag
     +          *( 0.5*(qpx(i,j,k,m)**2+qpy(i,j,k,m)**2
     +              +qpz(i,j,k,m)**2)/(aqrho**2)+
     +            qpres(i,j,k,m)*rmassq/aqrho/gamma1)
c        
         efldy(i,j,k)=(hpx(i,j,k,m)*abx+hpy(i,j,k,m)*aby
     +                 +hpz(i,j,k,m)*abz)/bmag
     +          *( 0.5*(hpx(i,j,k,m)**2+hpy(i,j,k,m)**2
     +              +hpz(i,j,k,m)**2)/(ahrho**2)+
     +            hpres(i,j,k,m)*rmassh/ahrho/gamma1)
c        
         efldz(i,j,k)=(opx(i,j,k,m)*abx+opy(i,j,k,m)*aby
     +                 +opz(i,j,k,m)*abz)/bmag
     +          *( 0.5*(opx(i,j,k,m)**2+opy(i,j,k,m)**2
     +              +opz(i,j,k,m)**2)/(aorho**2)+
     +            opres(i,j,k,m)*rmasso/aorho/gamma1)
c
         curx(i,j,k)=(hpres(i,j,k,m)+opres(i,j,k,m))/
     +             (ahrho/rmassh+aorho/rmasso)
  550 continue
c
c      call aurora(efldx,nx,ny,nz,m,radstrt,re_equiv,1,
c    +          ut,save_dat,add_dip,'qerg nth',14,write_dat)
c      call aurora(efldx,nx,ny,nz,m,radstrt,re_equiv,-1,
c    +          ut,save_dat,add_dip,'qerg sth',14,write_dat)
c
c      call aurora(efldy,nx,ny,nz,m,radstrt,re_equiv,1,
c    +          ut,save_dat,add_dip,'herg nth',14,write_dat)
c      call aurora(efldy,nx,ny,nz,m,radstrt,re_equiv,-1,
c    +          ut,save_dat,add_dip,'herg sth',14,write_dat)
c
c      call aurora(efldz,nx,ny,nz,m,radstrt,re_equiv,1,
c    +          ut,save_dat,add_dip,'oerg nth',14,write_dat)
c      call aurora(efldz,nx,ny,nz,m,radstrt,re_equiv,-1,
c    +          ut,save_dat,add_dip,'oerg sth',14,write_dat)
c
c      call aurora(curx,nx,ny,nz,m,radstrt,re_equiv,1,
c    +          ut,save_dat,add_dip,'itempnth',14,write_dat)
c      call aurora(curx,nx,ny,nz,m,radstrt,re_equiv,-1,
c    +          ut,save_dat,add_dip,'itempsth',14,write_dat)
c
      do k=1,nz
      do j=1,ny
      do i=1,nx
c       efldy(i,j,k)=sqrt(hrho(i,j,k,m)/rmassh+qrho(i,j,k,m)/rmassq)
c       efldz(i,j,k)=sqrt(orho(i,j,k,m)/rmasso)
c       efldx(i,j,k)=sqrt(qrho(i,j,k,m)/rmassq)
c       curx(i,j,k)=sqrt(qpres(i,j,k,m))
c       cury(i,j,k)=sqrt(hpres(i,j,k,m)+qpres(i,j,k,m))
c       curz(i,j,k)=sqrt(opres(i,j,k,m))
c
        efldy(i,j,k)=(hrho(i,j,k,m)/rmassh+qrho(i,j,k,m)/rmassq)
        efldz(i,j,k)=(orho(i,j,k,m)/rmasso)
        efldx(i,j,k)=(qrho(i,j,k,m)/rmassq)
        curx(i,j,k)=(qpres(i,j,k,m))
        cury(i,j,k)=(hpres(i,j,k,m)+qpres(i,j,k,m))
        curz(i,j,k)=(opres(i,j,k,m))
      enddo
      enddo
      enddo
c
      call  auroras(cury,efldy,bsx,bsy,bsz,
     +          nx,ny,nz,ngrd,m,radstrt,re_equiv,1,
     +          ut,save_dat,add_dip,'itempnth',14,write_dat)
      call  auroras(cury,efldy,bsx,bsy,bsz,
     +          nx,ny,nz,ngrd,m,radstrt,re_equiv,-1,
     +          ut,save_dat,add_dip,'itempsth',14,write_dat)
      call  auroras(curz,efldz,bsx,bsy,bsz,
     +          nx,ny,nz,ngrd,m,radstrt,re_equiv,1,
     +          ut,save_dat,add_dip,'otempnth',14,write_dat)
      call  auroras(curz,efldz,bsx,bsy,bsz,
     +          nx,ny,nz,ngrd,m,radstrt,re_equiv,-1,
     +          ut,save_dat,add_dip,'otempsth',14,write_dat)
      call  auroras(curx,efldx,bsx,bsy,bsz,
     +          nx,ny,nz,ngrd,m,radstrt,re_equiv,1,
     +          ut,save_dat,add_dip,'qtempnth',14,write_dat)
      call  auroras(curx,efldx,bsx,bsy,bsz,
     +          nx,ny,nz,ngrd,m,radstrt,re_equiv,-1,
     +          ut,save_dat,add_dip,'qtempsth',14,write_dat)
      do k=1,nz
      do j=1,ny
      do i=1,nx
c       efldy(i,j,k)=sqrt(hrho(i,j,k,m)/rmassh+qrho(i,j,k,m)/rmassq+
c    +               orho(i,j,k,m)/rmasso)
c        efldx(i,j,k)=sqrt(hrho(i,j,k,m)/rmassh)
c       curx(i,j,k)=sqrt(hpres(i,j,k,m))
c       cury(i,j,k)=sqrt(epres(i,j,k,m))
c
        efldy(i,j,k)=(hrho(i,j,k,m)/rmassh+qrho(i,j,k,m)/rmassq+
     +               orho(i,j,k,m)/rmasso)
        efldx(i,j,k)=(hrho(i,j,k,m)/rmassh)
        curx(i,j,k)=(hpres(i,j,k,m))
        cury(i,j,k)=(epres(i,j,k,m))
      enddo
      enddo
      enddo
      call  auroras(curx,efldx,bsx,bsy,bsz,
     +          nx,ny,nz,ngrd,m,radstrt,re_equiv,1,
     +          ut,save_dat,add_dip,'htempnth',14,write_dat)
      call  auroras(curx,efldx,bsx,bsy,bsz,
     +          nx,ny,nz,ngrd,m,radstrt,re_equiv,-1,
     +          ut,save_dat,add_dip,'htempsth',14,write_dat)
      call  auroras(cury,efldy,bsx,bsy,bsz,
     +          nx,ny,nz,ngrd,m,radstrt,re_equiv,1,
     +          ut,save_dat,add_dip,'etempnth',14,write_dat)
      call  auroras(cury,efldy,bsx,bsy,bsz,
     +          nx,ny,nz,ngrd,m,radstrt,re_equiv,-1,
     +          ut,save_dat,add_dip,'etempsth',14,write_dat)
c     call frame
      tgraf=tgraf+deltg
c
c     !!!!!! for B field check (for bcheck6.dat)
c      do m=1,ngrd
c      write(56,66)bx0(xloc,yloc,zloc,m),by0(xloc,yloc,zloc,m),
c     +bz0(xloc,yloc,zloc,m),bx(xloc,yloc,zloc,m),by(xloc,yloc,zloc,m),
c     +bz(xloc,yloc,zloc,m),ut,delt,t_equiv,m
c66    format(10(1x,1pe9.2))
c      enddo
c     !!!!!!
c      !!!!!! We add the inner bndry rotation with time & surface 
c       call bndry_inner_rot(qrho,qpres,rmassq,
c     +       hrho,hpres,rmassh,
c     +       orho,opres,rmasso,epres,gamma1,eerg,
c     +       ngrd,nx,ny,nz,parm_srf,parm_mid,
c     +       ijsrf,numsrf,ijmid,nummid,
c     +       erho,ti_te,o_conc)
c      !!!!!! test if parm_srf is changing (should be!)
       print *, '12','hrho_srf=', parm_srf(2,2), 'hpres=', parm_srf(5,2)
       print *, 'qrho_srf=', parm_srf(1,2), 'epres=', parm_srf(7,2)
c      !!!!!!
c
c
  600 if(t.lt.ts1)goto 700
      if(nchf.eq.11)
     +open(11,file='fluid11',status='unknown',form='unformatted')
      if(nchf.eq.12)
     +open(12,file='fluid12',status='unknown',form='unformatted')
c      if(nchf.eq.13)
c     +open(13,file='fluid13',status='unknown',form='unformatted')
c      if(nchf.eq.14)
c     +open(14,file='fluid14',status='unknown',form='unformatted')
c      if(nchf.eq.15)
c     +open(15,file='fluid15',status='unknown',form='unformatted')  
c      if(nchf.eq.16)
c     +open(16,file='fluid16',status='unknown',form='unformatted')
c      if(nchf.eq.17)
c     +open(17,file='fluid17',status='unknown',form='unformatted')
c      if(nchf.eq.18)
c     +open(18,file='fluid18',status='unknown',form='unformatted')
c
c      !!!!!!add write(57,67) check duplicate number first
c     !!!!!! for B field check (for bcheck7.dat)
c      do m=1,ngrd
c      write(57,67)bx0(xloc,yloc,zloc,m),by0(xloc,yloc,zloc,m),
c     +bz0(xloc,yloc,zloc,m),bx(xloc,yloc,zloc,m),by(xloc,yloc,zloc,m),
c     +bz(xloc,yloc,zloc,m),ut,delt,t_equiv,m
c67    format(10(1x,1pe9.2))
c      enddo
c     !!!!!!

c
c     write restart data
c
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
      if(Adiffuse)then
      write(nchf)qrho0
      write(nchf)hrho0
      write(nchf)orho0
      write(nchf)qpres0
      write(nchf)hpres0
      write(nchf)opres0
      write(nchf)epres0
      endif
c     !!!!!!
      write(nchf)parm_srf,parm_mid,
     +        ijzero,numzero,ijmid,nummid,ijsrf,numsrf
c     write(nchf)abx0,aby0,abz0
c     write(nchf)apx0,apy0,apz0 
c     write(nchf)aerg0
      close(nchf)
c     nchf=23-nchf
      nchf=nchf+1
      if(nchf.gt.12)nchf=11
      ts1=ts1+tsave
c     !!!!!!
      if(Adiffuse)then
      if(reload)then
        do m=1,ngrd
        do k=1,nz
        do j=1,ny
        do i=1,nx
         qrho0(i,j,k,m)=0.1*qrho(i,j,k,m)+0.9*qrho0(i,j,k,m)
         hrho0(i,j,k,m)=0.1*hrho(i,j,k,m)+0.9*hrho0(i,j,k,m)
         orho0(i,j,k,m)=0.1*orho(i,j,k,m)+0.9*orho0(i,j,k,m)
c
         qpres0(i,j,k,m)=0.1*qpres(i,j,k,m)+0.9*qpres0(i,j,k,m)
         hpres0(i,j,k,m)=0.1*hpres(i,j,k,m)+0.9*hpres0(i,j,k,m)
         opres0(i,j,k,m)=0.1*opres(i,j,k,m)+0.9*opres0(i,j,k,m)
         epres0(i,j,k,m)=0.1*epres(i,j,k,m)+0.9*epres0(i,j,k,m)
        enddo
        enddo
        enddo
        enddo   
      endif
      endif
c      !!!!!!for the rotation of velocity and momentum for push_ion
c       do m=1,ngrd
c         dx=(grd_xmax(m)-grd_xmin(m))/(nx-1.)
c         dy=(grd_ymax(m)-grd_ymin(m))/(ny-1.)
c         dz=(grd_zmax(m)-grd_zmin(m))/(nz-1.)
c       do k=1,nz
c       az=(grd_zmin(m)+dz*(k-1)-zdip)
c       do j=1,ny
c       ay=grd_ymin(m)+dy*(j-1)-ydip
c       do i=1,nx
c       ax=(grd_xmin(m)+dx*(i-1)-xdip)
c         xp_temp=ax*cos_tilt-az*sin_tilt
c         yp_temp=ay
c         zp_temp=ax*sin_tilt+az*cos_tilt
c       !!!!!!
c         superwholecoor_temp=(spinvectorx*xp_temp+spinvectory*yp_temp
c     +   +spinvectorz*zp_temp)*(1.0-cos_spin)
c
c          xp=xp_temp*cos_spin-(spinvectory*zp_temp
c     +   -spinvectorz*yp_temp)*sin_spin+spinvectorx*superwholecoor_temp
c          yp=yp_temp*cos_spin-(spinvectorz*xp_temp
c     +   -spinvectorx*zp_temp)*sin_spin+spinvectory*superwholecoor_temp
c          zp=zp_temp*cos_spin-(spinvectorx*yp_temp
c     +   -spinvectory*xp_temp)*sin_spin+spinvectorz*superwholecoor_temp
c       !!!!!!
c         ar=sqrt(xp**2+yp**2+zp**2)
c        !!!!!!need to rotate the coordinate 2 steps here!!!!!!
c        !!!!!!if the rho,pres and px,y,z modifying is right!!!
c        !!!!!!but unfortuntely not!!!!!!
c        !!!!!!but still need to rotate during our modify!!!!!!
c
c         call vfldrot(rvx_temp,rvy_temp,rvz_temp,ax,ay,az,vfrac,
c     +v_rot,r_rot)
c         ar_iono=sqrt(xp**2+yp**2+zp**2)
cc        ar_iono=amax1(0.0001,ar_iono)
c         ar_iono=amax1(1.01*rearth,ar_iono)
cc        ra=((ar_iono+0.5*rearth)/(1.5*rearth))**(-alpha_e)
c         ra=exp(-(ar_iono-rearth)/(0.225*rearth))
c         zheight=amax1(1.,(zp**2+(1.5*rearth)**2)/(3.0*rearth)**2)
c         ra=ra/zheight**2
c         rho_iono=amax1(erho*ra,0.001)
cc  
c         r_equat=(ar**3+0.001)/(xp**2+yp**2+0.001)
c         r_equat=amax1(r_equat,rearth)
c         erg_sphere=eerg*(0.001+(rearth/r_equat)**4)
cc        r_equat=amin1(r_equat,grd_xmax)
cc        erg_sphere=eerg*r_equat/grd_xmax
cc
c         arho=amax1(rho_iono,d_min)
c         aerg=amin1(abs(erg_sphere),0.01*arho)
c         
c           epres0_temp(i,j,k,m)=0.5*gamma1*aerg/ti_te
c          
c           qrho0_temp(i,j,k,m)=rmassq*d_min
c           qpres0_temp(i,j,k,m)=0.
c           qpx0_temp(i,j,k,m)=qrho(i,j,k,m)*rvx_temp
c           qpy0_temp(i,j,k,m)=qrho(i,j,k,m)*rvy_temp
c           qpz0_temp(i,j,k,m)=qrho(i,j,k,m)*rvz_temp
c           qpz(i,j,k,m)=0.
c          !!!!!! instead by rvz term above       
c           hrho0_temp(i,j,k,m)=arho*rmassh
c           hpres0_temp(i,j,k,m)=0.5*gamma1*aerg
c           hpx0_temp(i,j,k,m)=hrho(i,j,k,m)*rvx_temp
c           hpy0_temp(i,j,k,m)=hrho(i,j,k,m)*rvy_temp
c           hpz0_temp(i,j,k,m)=hrho(i,j,k,m)*rvz_temp
cc           hpz(i,j,k,m)=0.
c          !!!!!! instead by rvz term above    
c           orho0_temp(i,j,k,m)=rho_iono*rmasso*o_conc
c           opres0_temp(i,j,k,m)=0.5*gamma1*aerg*o_conc
c           opx0_temp(i,j,k,m)=orho(i,j,k,m)*rvx_temp
c           opy0_temp(i,j,k,m)=orho(i,j,k,m)*rvy_temp
c           opz0_temp(i,j,k,m)=orho(i,j,k,m)*rvz_temp
cc           opz(i,j,k,m)=0.
cc          !!!!!! instead by rvz term above    
cc      !!!!!I think maybe I should use the section of fraction?!
cc      !!!!!! Do I need to add the boundary part followed orginally?
c       enddo
c       enddo
c       enddo
c
c       enddo
c      !!!!!! We add the inner bndry rotation with time & surface 
c       call bndry_inner_rot(qrho,qpres,rmassq,
c     +       hrho,hpres,rmassh,
c     +       orho,opres,rmasso,epres,gamma1,eerg,
c     +       ngrd,nx,ny,nz,parm_srf,parm_mid,
c     +       ijsrf,numsrf,ijmid,nummid,
c     +       erho,ti_te,o_conc)
c      !!!!!! test if parm_srf is changing (should be!)
       print *, '13','hrho_srf=', parm_srf(2,2), 'hpres=', parm_srf(5,2)
       print *, 'qrho_srf=', parm_srf(1,2), 'epres=', parm_srf(7,2)
c      !!!!!!
  700 if(t.lt.tmax)goto 1000
c
c      final diagnostics
c
c      m=ngrd-1
c      rx=xspac(m)
c      ry=xspac(m)
c      rz=xspac(m)
c
c      dm=grd_zmax(m)+zdip-1.
c      xtot=1.3*(2.*dm)
c      xmin=grd_xmin(m)
c      zmin=grd_zmin(m)+1.
c
c      xmax=xmin+xtot
c      zmax=-zmin
c      ymax=zmax
c      ymin=zmin
c
c      add_dip=.false.
c      radstrt=rearth+4
c      theta1=3.141593/5.
c      theta2=3.141593/7.
c      theta0=3.141593/2.5
c      nphi=36
c      ncuts=3
c      pi=3.1416
c
c       call qvset(0.,bsx,nx*ny*nz)
c       call qvset(0.,bsy,nx*ny*nz)
c       call qvset(0.,bsz,nx*ny*nz)
c
c     call totfld(bx,bx0,bsx,nx,ny,nz,ngrd,m)
c     call totfld(by,by0,bsy,nx,ny,nz,ngrd,m)
c     call totfld(bz,bz0,bsz,nx,ny,nz,ngrd,m)
c     do k=1,nz
c     do j=1,ny
c     do i=1,nx
c       efldx(i,j,k)=sqrt(qpres(i,j,k,m)+epres(i,j,k,m))
c     enddo
c     enddo
c     enddo
c
c      add_dip=.false.
c     call contrace(efldx,bsx,bsy,bsz,nx,ny,nz,ngrd,m,
c    +            xmin,xmax,ymin,ymax,zmin,zmax,1,
c    +            t,'fld pi3',3,11,add_dip,
c    +            radstrt,nphi,pi/3.,theta2,1,
c    +            tx,ty,tz,tg2,work,mx,my,mz,muvwp2,mz)
c
       call clsgks
c    !!!!!!for sign the code runs to end
        open(20,file='codeend',status='unknown')
        write(20,*)12345
c     !!!!!!
       end 
c
c     **************************************************
c
      subroutine set_imf(bx,by,bz,bx0,by0,bz0,bxp,byp,bzp,
     +        qrho,qpres,px,py,pz,epres,
     +        rhop,svxp,svyp,svzp,svelx,spress,ti_te,
     +        rho_frac,nx,ny,nz,ngrd)
c
c     set imf boundary conditions
c
      dimension bx(nx,ny,nz,ngrd),bx0(nx,ny,nz,ngrd),bxp(ny,nz),
     +   by(nx,ny,nz,ngrd),by0(nx,ny,nz,ngrd),byp(ny,nz),
     +   bz(nx,ny,nz,ngrd),bz0(nx,ny,nz,ngrd),bzp(ny,nz),
     +   rhop(ny,nz),epres(nx,ny,nz,ngrd),
     +   qrho(nx,ny,nz,ngrd),qpres(nx,ny,nz,ngrd),
     +   px(nx,ny,nz,ngrd),py(nx,ny,nz,ngrd),pz(nx,ny,nz,ngrd),
c     +   orho(nx,ny,nz,ngrd),opres(nx,ny,nz,ngrd),
c     +   opx(nx,ny,nz,ngrd),opy(nx,ny,nz,ngrd),opz(nx,ny,nz,ngrd),
     +   svxp(ny,nz),svyp(ny,nz),svzp(ny,nz)
c
      m=ngrd
c
      i=1
      do 10 j=1,ny
      do 10 k=1,nz
        bx(i,j,k,m)=bxp(j,k)-bx0(i,j,k,ngrd) 
        by(i,j,k,m)=byp(j,k)-by0(i,j,k,ngrd)
        bz(i,j,k,m)=bzp(j,k)-bz0(i,j,k,ngrd)  
c 
        qrho(i,j,k,m)=rhop(j,k)
        px(i,j,k,m)=rhop(j,k)*svxp(j,k)
        py(i,j,k,m)=rhop(j,k)*svyp(j,k)
        pz(i,j,k,m)=rhop(j,k)*svzp(j,k) 
c
c        orho(i,j,k,m)=rho_frac*rhop(j,k)
c        opx(i,j,k,m)=orho(i,j,k,m)*svxp(j,k)
c        opy(i,j,k,m)=orho(i,j,k,m)*svyp(j,k)
c        opz(i,j,k,m)=orho(i,j,k,m)*svzp(j,k)
c
        qpres(i,j,k,m)=0.5*spress
        epres(i,j,k,m)=0.5*spress/ti_te
c        opres(i,j,k,m)=rho_frac*0.5*spress
   10 continue
c
      return 
      end
c
c     **************************************************
c
      subroutine set_resist(rst,nx,ny,nz,resist,ijzero,numzero,
     +             ijmid,nummid,ijsrf,numsrf)
c
c      set resistance around the earth :
c        include dayside and auroral conductivities
c      magntiude set by resist
c      width is set by del_ang=3.3 degrees = 0.058 rads
c      radial dropoff as alpha=-8
c      shifted of dipole axis by 2.5 degress 
c      radius of 22.5 degrees
c
      dimension rst(nx,ny,nz) 
      integer ijsrf(3,15000),ijmid(3,15000),ijzero(3,120000)
      common /gridding/grd_xmin(9),grd_xmax(9),grd_ymin(9),grd_ymax(9),
     +           grd_zmin(9),grd_zmax(9),xspac(9),ncore(9),nbndry(9),
     +           rx,ry,rz,xdip,ydip,zdip,rearth,b0,
     +           sin_tilt,cos_tilt
c
c     resistivity characteristics
c        initialize
c
      rst=0.
c
c
c***CAROL***   resist_char=4.*3.14159*1e-7*(v_equiv*1000.)*re_equiv*2.5e7 just adding
      resist_char=4.*3.14159*1e-7*(1000*1000.)*0.2000*2.5e7
c
c     interior resistivity
c
      do 10  n=1,numzero
        i=ijzero(1,n)
        j=ijzero(2,n)
        k=ijzero(3,n)
c        rst(i,j,k)=b0/resist
c***CAROL***  rst(i,j,k)=resist/resist_char instead of the line above
c       rst(i,j,k)=resist/resist_char
        rst(i,j,k)=8.0*resist/resist_char
   10 continue
c
c     lower ionosphere resistivity
c
      do 20  n=1,nummid
        i=ijmid(1,n)
        j=ijmid(2,n)
        k=ijmid(3,n)
c        rst(i,j,k)=0.5*b0/resist
c***CAROL***  rst(i,j,k)=0.5*resist/resist_char instead of the line above
c       rst(i,j,k)=0.5*resist/resist_char
        rst(i,j,k)=4.0*resist/resist_char
   20 continue
c
c     upper ionosphere
c
      do 30  n=1,numsrf
        i=ijsrf(1,n)
        j=ijsrf(2,n)
        k=ijsrf(3,n)
c        rst(i,j,k)=0.125*b0/resist
c***CAROL***  rst(i,j,k)=0.125*resist/resist_char instead of the line above
c       rst(i,j,k)=0.125*resist/resist_char
        rst(i,j,k)=1.0*resist/resist_char
   30 continue
c
      return
      end
c
c     ************************************************
c
      subroutine totfld(bx,bx0,btx,nx,ny,nz,ngrd,m)
c
c
c     calculates the total magnetic field from the perturbed and
c        stationary magnetic field !!!!!!(unperturbed/dipole or SW?)
c
      dimension bx(nx,ny,nz,ngrd),bx0(nx,ny,nz,ngrd),btx(nx,ny,nz)
c
c$omp  parallel do 
      do k=1,nz
      do j=1,ny
      do i=1,nx
        btx(i,j,k)=bx0(i,j,k,m)+bx(i,j,k,m)
      enddo
      enddo
      enddo
c
      return
      end
c
c     ************************************************
c
      subroutine fnd_fld(bx,btx,ngrd,nx,ny,nz,m)
c
c
c     calculates the total DYNAMIC magnetic field only
c        dipole field added separately by by setting add_dip=.true.
c
      dimension bx(nx,ny,nz,ngrd),btx(nx,ny,nz)
c
c$omp  parallel do 
      do k=1,nz
      do j=1,ny
      do i=1,nx
        btx(i,j,k)=btx(i,j,k)+bx(i,j,k,m)
      enddo
      enddo
      enddo
c
      return
      end
c
c     ************************************************
C
      subroutine fnd_vel(px,py,pz,rho,vx,vy,vz,ngrd,nx,ny,nz,m)
c
c     converts momentum into velocity for graphics
c
      dimension px(nx,ny,nz,ngrd),py(nx,ny,nz,ngrd),
     +          pz(nx,ny,nz,ngrd),rho(nx,ny,nz,ngrd),
     +          vx(nx,ny,nz),vy(nx,ny,nz),vz(nx,ny,nz)
c
c$omp  parallel do 
      do k=1,nz
      do j=1,ny
      do i=1,nx
       arho=amax1(rho(i,j,k,m),0.0001)
       vx(i,j,k)=px(i,j,k,m)/arho
       vy(i,j,k)=py(i,j,k,m)/arho
       vz(i,j,k)=pz(i,j,k,m)/arho
      enddo
      enddo
      enddo
c
      return
      end
c
c     ************************************************
C
      subroutine fnd_vtot(qpx,qpy,qpz,qrho,hpx,hpy,hpz,hrho,
     +      opx,opy,opz,orho,vx,vy,vz,ngrd,nx,ny,nz,m,
     +      rmassq,rmassh,rmasso)
c
c     converts momentum into velocity for graphics
c
      dimension qpx(nx,ny,nz,ngrd),qpy(nx,ny,nz,ngrd),
     +          qpz(nx,ny,nz,ngrd),qrho(nx,ny,nz,ngrd),
     +          hpx(nx,ny,nz,ngrd),hpy(nx,ny,nz,ngrd),
     +          hpz(nx,ny,nz,ngrd),hrho(nx,ny,nz,ngrd),
     +          opx(nx,ny,nz,ngrd),opy(nx,ny,nz,ngrd),
     +          opz(nx,ny,nz,ngrd),orho(nx,ny,nz,ngrd),
     +          vx(nx,ny,nz),vy(nx,ny,nz),vz(nx,ny,nz)
c
c  
c$omp  parallel do    
      do k=1,nz
      do j=1,ny
      do i=1,nx
c      eden=amax1(qrho(i,j,k,m)+hrho(i,j,k,m)+orho(i,j,k,m),
c    +                 0.0001)
       qden=(qrho(i,j,k,m)+0.000001)/rmassq
       hden=(hrho(i,j,k,m)+0.000001)/rmassh
       oden=(orho(i,j,k,m)+0.000001)/rmasso
       tden=qden+hden+oden
c
       vx(i,j,k)=(qpx(i,j,k,m)/rmassq+hpx(i,j,k,m)/rmassh
     +               +opx(i,j,k,m)/rmasso)/tden
       vy(i,j,k)=(qpy(i,j,k,m)/rmassq+hpy(i,j,k,m)/rmassh
     +               +opy(i,j,k,m)/rmasso)/tden
       vz(i,j,k)=(qpz(i,j,k,m)/rmassq+hpz(i,j,k,m)/rmassh
     +               +opz(i,j,k,m)/rmasso)/tden
c      eden=qrho(i,j,k,m)+hrho(i,j,k,m)+orho(i,j,k,m)
c      vx(i,j,k)=(qpx(i,j,k,m)+hpx(i,j,k,m)+opx(i,j,k,m))/eden
c      vy(i,j,k)=(qpy(i,j,k,m)+hpy(i,j,k,m)+opy(i,j,k,m))/eden
c      vz(i,j,k)=(qpz(i,j,k,m)+hpz(i,j,k,m)+opz(i,j,k,m))/eden
      enddo
      enddo
      enddo
c
      return
      end
c
c     ************************************************
C
      subroutine fnd_evel(qpx,qpy,qpz,qrho,hpx,hpy,hpz,hrho,
     +      opx,opy,opz,orho,curx,cury,curz,
     +      evx,evy,evz,tvx,tvy,tvz,
     +      ngrd,nx,ny,nz,m,rmassq,rmassh,rmasso,reynolds)
c
c     converts momentum into velocity for graphics
c
      dimension qpx(nx,ny,nz,ngrd),qpy(nx,ny,nz,ngrd),
     +          qpz(nx,ny,nz,ngrd),qrho(nx,ny,nz,ngrd),
     +          hpx(nx,ny,nz,ngrd),hpy(nx,ny,nz,ngrd),
     +          hpz(nx,ny,nz,ngrd),hrho(nx,ny,nz,ngrd),
     +          opx(nx,ny,nz,ngrd),opy(nx,ny,nz,ngrd),
     +          opz(nx,ny,nz,ngrd),orho(nx,ny,nz,ngrd),
     +          curx(nx,ny,nz),cury(nx,ny,nz),curz(nx,ny,nz),
     +          tvx(nx,ny,nz),tvy(nx,ny,nz),tvz(nx,ny,nz),
     +          evx(nx,ny,nz),evy(nx,ny,nz),evz(nx,ny,nz)
c
c     
c$omp  parallel do 
      do k=1,nz
      do j=1,ny
      do i=1,nx
c      eden=amax1(qrho(i,j,k,m)+hrho(i,j,k,m)+orho(i,j,k,m),
c    +                 0.0001)
       qden=(qrho(i,j,k,m)+0.000001)/rmassq
       hden=(hrho(i,j,k,m)+0.000001)/rmassh
       oden=(orho(i,j,k,m)+0.000001)/rmasso
       tden=qden+hden+oden
c
c      keep sepearate the ion and current components
c
       tvx(i,j,k)=(qpx(i,j,k,m)/rmassq+hpx(i,j,k,m)/rmassh
     +               +opx(i,j,k,m)/rmasso)/tden
       evx(i,j,k)= tvx(i,j,k) - curx(i,j,k)/tden/reynolds
c
       tvy(i,j,k)=(qpy(i,j,k,m)/rmassq+hpy(i,j,k,m)/rmassh
     +               +opy(i,j,k,m)/rmasso)/tden
       evy(i,j,k)= tvy(i,j,k) - cury(i,j,k)/tden/reynolds
c
       tvz(i,j,k)=(qpz(i,j,k,m)/rmassq+hpz(i,j,k,m)/rmassh
     +               +opz(i,j,k,m)/rmasso)/tden
       evz(i,j,k)= tvz(i,j,k)  -curz(i,j,k)/tden/reynolds  
      enddo
      enddo
      enddo
c
      return
      end

c
c     ************************************************
C
      subroutine fnd_prss(px,py,pz,rho,erg,afld,ngrd,nx,ny,nz,
     +                        m,gamma1,igo)
c
c     now using pressure equation rather than an erg equation
c
      dimension px(nx,ny,nz,ngrd),py(nx,ny,nz,ngrd),
     +          pz(nx,ny,nz,ngrd),rho(nx,ny,nz,ngrd),
     +          erg(nx,ny,nz,ngrd),afld(nx,ny,nz)
c
c$omp  parallel do 
      do k=1,nz
      do j=1,ny
      do i=1,nx
c      arho=amax1(rho(i,j,k,m),0.005)
c      afld(i,j,k)=gamma1*(erg(i,j,k,m)-0.5*
c    +   (px(i,j,k,m)**2+py(i,j,k,m)**2+pz(i,j,k,m)**2)
c    +                    /arho)
       afld(i,j,k)=erg(i,j,k,m)
      enddo
      enddo
      enddo
c
c     if graphics refine plots
c
      if(igo.le.0)return

c$omp  parallel do 
      do k=1,nz
      do j=1,ny
      do i=1,nx
        afld(i,j,k)=sqrt(abs(afld(i,j,k)))
        afld(i,j,k)=amin1(afld(i,j,k),1.2)
      enddo
      enddo
      enddo
c
      return
      end
c
c     *************************************************
c
      subroutine chg_dip(bx0,by0,bz0,b0,b1,nx,ny,nz,ngrd)
c
c        this subroutine will scale up and existing
c        dipole magnetic field
c
      dimension bx0(nx,ny,nz,ngrd),by0(nx,ny,nz,ngrd),
     +          bz0(nx,ny,nz,ngrd)
c
      bscale=b0/b1
      do 10 m=1,ngrd
      do 10 k=1,nz
      do 10 j=1,ny
      do 10 i=1,nx
        bx0(i,j,k,m)=bscale*bx0(i,j,k,m)
        by0(i,j,k,m)=bscale*by0(i,j,k,m)
        bz0(i,j,k,m)=bscale*bz0(i,j,k,m)
   10 continue
c
      return
      end
c
c     *************************************************
c
      subroutine tot_b(btot,bsx,bsy,bsz,nx,ny,nz)
c
c      Initialize static magnetic field along entire grid
c

c
      dimension bsx(nx,ny,nz),bsy(nx,ny,nz),bsz(nx,ny,nz),
     +      btot(nx,ny,nz)
c
      do k=1,nz
      do j=1,ny
      do i=1,nx
        atot=sqrt(bsx(i,j,k)**2+bsy(i,j,k)**2+bsz(i,j,k)**2)
        btot(i,j,k)=amax1(atot,1.e-5)
      enddo
      enddo
      enddo
c
      return
      end
c
c
c     *************************************************
c
      subroutine mak_dip(bx0,by0,bz0,nx,ny,nz,
     +              ngrd,ijzero,numzero)
c
c      Initialize static magnetic field along entire grid
c

      common /gridding/grd_xmin(9),grd_xmax(9),grd_ymin(9),grd_ymax(9),
     +           grd_zmin(9),grd_zmax(9),xspac(9),ncore(9),nbndry(9),
     +           rx,ry,rz,xdip,ydip,zdip,rearth,b0,
     +           sin_tilt,cos_tilt
c
      dimension bx0(nx,ny,nz,ngrd),by0(nx,ny,nz,ngrd),
     +           bz0(nx,ny,nz,ngrd)
      integer ijzero(3,120000)
c
      do 140 m=1,ngrd
       dx=(grd_xmax(m)-grd_xmin(m))/(nx-1.)
       dy=(grd_ymax(m)-grd_ymin(m))/(ny-1.)
       dz=(grd_zmax(m)-grd_zmin(m))/(nz-1.)
c
      do 130 k=1,nz
       az=grd_zmin(m)+dz*(k-1)
       z1=(az-zdip)
c
       do 120 j=1,ny
        ay=grd_ymin(m)+dy*(j-1)
        y1=(ay-ydip)
c
        do 110 i=1,nx
         ax=grd_xmin(m)+dx*(i-1)
         x1=(ax-xdip)
c
c        determine magnetic dipole field
c
c
c        real space to dipole space
c
         xp=x1*cos_tilt-z1*sin_tilt
         zp=x1*sin_tilt+z1*cos_tilt
c
         x2=xp**2
         y2=y1**2
         z2=zp**2
         ar=sqrt(x2+y2+z2)
c
c        cartesian equivalent
c        !!!!!! let dipole north at beginning.
c         bmag=b0/ar**5
         bmag=-1.0*b0/ar**5
         dbx=-3.*bmag*xp*zp
         dbz=bmag*(x2+y2-2.*z2)
c
c        rotate b field back to coordinate space
c
        if(ar.gt.rearth-1.5) then
          bx0(i,j,k,m) =dbx*cos_tilt+dbz*sin_tilt
          bz0(i,j,k,m)=-dbx*sin_tilt+dbz*cos_tilt
c
          by0(i,j,k,m)=-3.*bmag*y1*zp
        else
          bx0(i,j,k,m)=0.
          by0(i,j,k,m)=0.
          bz0(i,j,k,m)=0.
        endif
c
  110   continue
  120  continue
  130 continue
c
  140 continue
c
c     boundary conditions
c
      m=1
      do 210 n=1,numzero
c
c        get coords of point
c
         i=ijzero(1,n)
         j=ijzero(2,n)
         k=ijzero(3,n)
c
         bx0(i,j,k,m)=0.
         by0(i,j,k,m)=0.
         bz0(i,j,k,m)=0.
  210 continue
      return
      end
c
c     *************************************************
c
      subroutine dipole(abx,aby,abz,ax,ay,az)
c
c     calculates magnetic field for a dipole
c
      common /gridding/grd_xmin(9),grd_xmax(9),grd_ymin(9),grd_ymax(9),
     +           grd_zmin(9),grd_zmax(9),xspac(9),ncore(9),nbndry(9),
     +           rx,ry,rz,xdip,ydip,zdip,rearth,b0,
     +           sin_tilt,cos_tilt
      common /uranus/tilt3,sin_tilt3,cos_tilt3,sin_spin,cos_spin,
     +           tilt4,sin_tilt4,cos_tilt4,spin,spinvectorx,
     +           spinvectory,spinvectorz
c
c    
         x1=(ax-xdip)
         y1=(ay-ydip)
         z1=(az-zdip)
c
c        real space to dipole space
c
         xp=x1*cos_tilt4-z1*sin_tilt4
         zp=x1*sin_tilt4+z1*cos_tilt4
c
         x2=xp**2
         y2=y1**2
         z2=zp**2
         ar=sqrt(x2+y2+z2)
c
c        cartesian equivalent
c        !!!!!! let dipole north at beginning.
c         bmag=b0/ar**5
         bmag=-1.0*b0/ar**5
         dbx=-3.*bmag*xp*zp
         aby=-3.*bmag*y1*zp
         dbz=bmag*(x2+y2-2.*z2)
c
c        rotate b field back to coordinate space
c
          abx=dbx*cos_tilt4+dbz*sin_tilt4
          abz=-dbx*sin_tilt4+dbz*cos_tilt4
c
c       add in mirror dipole magetic field
c
c        xdipm=2.*grd_xmin -xdip
c        xm1=(ax-xdipm)
c        ym1=(ay-ydip)
c        zm1=(az-zdip)
c
c        xm2=(ax-xdipm)**2
c        ym2=(ay-ydip)**2
c        zm2=(az-zdip)**2
c        amr=sqrt(xm2+y2+z2)
c
c        bmag=b0/amr**5
c        abx=abx-3.*bmag*xm1*z1
c        aby=aby-3.*bmag*y1*z1
c        abz=abz+bmag*(xm2+y2-2.*z2)
c
c
c        find equivalent spherical cordinates
c
c
c         cost=(z1)/ar
c         sint=sqrt(x2+y2)/ar
c         cosp=x1/(ar*sint)
c         sinp=y1/(ar*sint)
c
c        br=-2.*b0*cost/ar**3
c        bt=-b0*sint/ar**3
c
c        bx=br*sint*cosp+bt*cost*cosp
c        by=br*sint*sinp+bt*cost*sinp
c        bz=br*cost-bt*sint      
c
      return
      end
c
c     **************************************************************
c
       subroutine limcraft(xcraft,ncraft,re_equiv,m)
c
      common /gridding/grd_xmin(9),grd_xmax(9),grd_ymin(9),grd_ymax(9),
     +           grd_zmin(9),grd_zmax(9),xspac(9),ncore(9),nbndry(9),
     +           rx,ry,rz,xdip,ydip,zdip,rearth,b0,
     +           sin_tilt,cos_tilt
c
c      tests to see whether spacecraft is in the system and 
c           resets their position in not
c
      dimension xcraft(4,ncraft)
c
      abit=0.001
      do 10 n=1,ncraft
        xcraft(1,n)=amax1(xcraft(1,n),(grd_xmin(m)+abit)*re_equiv)
        xcraft(1,n)=amin1(xcraft(1,n),(grd_xmax(m)-abit)*re_equiv)
        xcraft(2,n)=amax1(xcraft(2,n),(grd_ymin(m)+abit)*re_equiv)
        xcraft(2,n)=amin1(xcraft(2,n),(grd_ymax(m)-abit)*re_equiv)
        xcraft(3,n)=amax1(xcraft(3,n),(grd_zmin(m)+abit)*re_equiv)
        xcraft(3,n)=amin1(xcraft(3,n),(grd_zmax(m)-abit)*re_equiv)
   10 continue
c
      return
      end
c
c     *************************************************
c
      subroutine mak_dipz(bx0,by0,bz0,nx,ny,nz,
     +              ngrd,ijzero,numzero)
c
c      Initialize static magnetic field along entire grid
c

      common /gridding/grd_xmin(9),grd_xmax(9),grd_ymin(9),grd_ymax(9),
     +           grd_zmin(9),grd_zmax(9),xspac(9),ncore(9),nbndry(9),
     +           rx,ry,rz,xdip,ydip,zdip,rearth,b0,
     +           sin_tilt,cos_tilt
      common /uranus/tilt3,sin_tilt3,cos_tilt3,sin_spin,cos_spin,
     +           tilt4,sin_tilt4,cos_tilt4,spin,spinvectorx,
     +           spinvectory,spinvectorz
c
      dimension bx0(nx,ny,nz,ngrd),by0(nx,ny,nz,ngrd),
     +           bz0(nx,ny,nz,ngrd)
      integer ijzero(3,120000)
c
c           sin_spin=sin(90*.0174533)
c           cos_spin=cos(90*.0174533)
      do 140 m=1,ngrd
       dx=(grd_xmax(m)-grd_xmin(m))/(nx-1.)
       dy=(grd_ymax(m)-grd_ymin(m))/(ny-1.)
       dz=(grd_zmax(m)-grd_zmin(m))/(nz-1.)
c
      do 130 k=1,nz
       az=grd_zmin(m)+dz*(k-1)
       z1=(az-zdip)
c
       do 120 j=1,ny
        ay=grd_ymin(m)+dy*(j-1)
        y1=(ay-ydip)
c
        do 110 i=1,nx
         ax=grd_xmin(m)+dx*(i-1)
         x1=(ax-xdip)
c
c        determine magnetic dipole field
c         spinvectorx=1.0
c         spinvectory=0.0
c         spinvectorz=0.0
c
c
c        real space to dipole space
c        around y axis
         xp00=x1*cos_tilt-z1*sin_tilt
         yp00=y1
         zp00=x1*sin_tilt+z1*cos_tilt
c          xp00=x1
c          yp00=y1
c          zp00=z1
c        around z axis
c         xp=xp00*cos_spin+yp00*sin_spin
c         yp=-1.0*xp00*sin_spin+yp00*cos_spin
c         zp=zp00
c          xp=xp00*cos_spin-zp00*sin_spin
c          yp=yp00
c          zp=xp00*sin_spin+zp00*cos_spin
c          xp=xp00
c          yp=yp00
c          zp=zp00
c          xp=x1*cos_tilt*cos_spin-y1*sin_spin+z1*sin_tilt*cos_spin
c          yp=x1*cos_tilt*sin_spin+y1*cos_spin+z1*sin_tilt*sin_spin
c          zp=-x1*sin_tilt+z1*cos_tilt
        superwholecoor=(spinvectorx*xp00+spinvectory*yp00
     +   +spinvectorz*zp00)*(1.0-cos_spin)
c
          xp=xp00*cos_spin-(spinvectory*zp00
     +   -spinvectorz*yp00)*sin_spin+spinvectorx*superwholecoor
          yp=yp00*cos_spin-(spinvectorz*xp00
     +   -spinvectorx*zp00)*sin_spin+spinvectory*superwholecoor
          zp=zp00*cos_spin-(spinvectorx*yp00
     +   -spinvectory*xp00)*sin_spin+spinvectorz*superwholecoor
c
c
         x2=xp**2
         y2=yp**2
         z2=zp**2
         ar=sqrt(x2+y2+z2)
c
c        cartesian equivalent
c        !!!!!! let dipole north at beginning
c         bmag=b0/ar**5
         bmag=-1.0*b0/ar**5
         dbx=-3.*bmag*xp*zp
         dby=-3.*bmag*yp*zp
         dbz=bmag*(x2+y2-2.*z2)
c
c        rotate b field back to coordinate space
c
        if(ar.gt.rearth-1.5) then
c          bx0(i,j,k,m) =dbx*cos_tilt+dbz*sin_tilt
c          bz0(i,j,k,m)=-dbx*sin_tilt+dbz*cos_tilt

c
c          by0(i,j,k,m)=-3.*bmag*y1*zp
c          by0(i,j,k,m)=dby
c         back around z axis 
c
c          bx0aa=dbx*cos_spin-dby*sin_spin
c          by0aa=1.0*dbx*sin_spin+dby*cos_spin
c          bz0aa=dbz
c
c           bx0aa=dbx*cos_spin+dbz*sin_spin
c           bz0aa=-1.0*dbx*sin_spin+dbz*cos_spin
c           by0aa=dby
c          bx0aa=dbx
c          by0aa=dby
c          bz0aa=dbz
c
        superwhole=(spinvectorx*dbx+spinvectory*dby
     +   +spinvectorz*dbz)*(1.0-cos_spin)
c
          bx0aa=dbx*cos_spin+(spinvectory*dbz
     +   -spinvectorz*dby)*sin_spin+spinvectorx*superwhole
          by0aa=dby*cos_spin+(spinvectorz*dbx
     +   -spinvectorx*dbz)*sin_spin+spinvectory*superwhole
          bz0aa=dbz*cos_spin+(spinvectorx*dby
     +   -spinvectory*dbx)*sin_spin+spinvectorz*superwhole
c
c         back around y axis
c
          bx0(i,j,k,m)=bx0aa*cos_tilt+bz0aa*sin_tilt
          bz0(i,j,k,m)=-1.0*bx0aa*sin_tilt+bz0aa*cos_tilt
          by0(i,j,k,m)=by0aa
c
c           bx0(i,j,k,m)=bx0aa
c           by0(i,j,k,m)=by0aa
c           bz0(i,j,k,m)=bz0aa
c
c             bx0(i,j,k,m)=dbx*cos_tilt*cos_spin-dby*cos_tilt*sin_spin
c     +                   +dbz*sin_tilt
c             by0(i,j,k,m)=dbx*sin_spin+dby*cos_spin
c             bz0(i,j,k,m)=-dbx*sin_tilt*cos_spin+dby*sin_tilt*sin_spin
c     +                   +dbz*cos_tilt
c           !!!!!!!!!!!! Rodrigues's rotation formula
c        superwhole=(spinvectorx*bx0(i,j,k,m)+spinvectory*by0(i,j,k,m)
c     +   +spinvectorz*bz0(i,j,k,m))*(1.0-cos_spin)
c
c          bx0bb=bx0(i,j,k,m)*cos_spin+(spinvectory*bz0(i,j,k,m)
c     +   -spinvectorz*by0(i,j,k,m))*sin_spin+spinvectorx*superwhole
c          by0bb=by0(i,j,k,m)*cos_spin+(spinvectorz*bx0(i,j,k,m)
c     +   -spinvectorx*bz0(i,j,k,m))*sin_spin+spinvectory*superwhole
c          bz0bb=bz0(i,j,k,m)*cos_spin+(spinvectorx*by0(i,j,k,m)
c     +   -spinvectory*bx0(i,j,k,m))*sin_spin+spinvectorz*superwhole
c
c        bx0(i,j,k,m)=bx0bb
c        by0(i,j,k,m)=by0bb
c        bz0(i,j,k,m)=bz0bb
c           !!!!!!!!!!!!
c
        else
          bx0(i,j,k,m)=0.
          by0(i,j,k,m)=0.
          bz0(i,j,k,m)=0.
        endif
c
  110   continue
  120  continue
  130 continue
c
  140 continue
c
c     boundary conditions
c
      m=1
      do 210 n=1,numzero
c
c        get coords of point
c
         i=ijzero(1,n)
         j=ijzero(2,n)
         k=ijzero(3,n)
c
         bx0(i,j,k,m)=0.
         by0(i,j,k,m)=0.
         bz0(i,j,k,m)=0.
  210 continue
      return
      end
c
c     *************************************************
c
c     *************************************************
c
      subroutine mak_dipz2(bx0,by0,bz0,nx,ny,nz,
     +              ngrd,ijzero,numzero)
c
c      Initialize static magnetic field along entire grid
c

      common /gridding/grd_xmin(9),grd_xmax(9),grd_ymin(9),grd_ymax(9),
     +           grd_zmin(9),grd_zmax(9),xspac(9),ncore(9),nbndry(9),
     +           rx,ry,rz,xdip,ydip,zdip,rearth,b0,
     +           sin_tilt,cos_tilt
      common /uranus/tilt3,sin_tilt3,cos_tilt3,sin_spin,cos_spin,
     +           tilt4,sin_tilt4,cos_tilt4,spin,spinvectorx,
     +           spinvectory,spinvectorz
c
      dimension bx0(nx,ny,nz,ngrd),by0(nx,ny,nz,ngrd),
     +           bz0(nx,ny,nz,ngrd)
      integer ijzero(3,120000)
c
c           sin_spin=sin(90*.0174533)
c           cos_spin=cos(90*.0174533)
      do 140 m=1,ngrd
       dx=(grd_xmax(m)-grd_xmin(m))/(nx-1.)
       dy=(grd_ymax(m)-grd_ymin(m))/(ny-1.)
       dz=(grd_zmax(m)-grd_zmin(m))/(nz-1.)
c
      do 130 k=1,nz
       az=grd_zmin(m)+dz*(k-1)
       z1=(az-zdip)
c
       do 120 j=1,ny
        ay=grd_ymin(m)+dy*(j-1)
        y1=(ay-ydip)
c
        do 110 i=1,nx
         ax=grd_xmin(m)+dx*(i-1)
         x1=(ax-xdip)
c
c        determine magnetic dipole field
c         spinvectorx=1.0
c         spinvectory=0.0
c         spinvectorz=0.0
c
c
c        real space to dipole space
c        around y axis
c         xp00=x1*cos_tilt-z1*sin_tilt
c         yp00=y1
c         zp00=x1*sin_tilt+z1*cos_tilt
c          xp00=x1
c          yp00=y1
c          zp00=z1
c        around z axis
c         xp=xp00*cos_spin+yp00*sin_spin
c         yp=-1.0*xp00*sin_spin+yp00*cos_spin
c         zp=zp00
c          xp=xp00*cos_spin-zp00*sin_spin
c          yp=yp00
c          zp=xp00*sin_spin+zp00*cos_spin
c          xp=xp00
c          yp=yp00
c          zp=zp00
c          xp=x1*cos_tilt*cos_spin-y1*sin_spin+z1*sin_tilt*cos_spin
c          yp=x1*cos_tilt*sin_spin+y1*cos_spin+z1*sin_tilt*sin_spin
c          zp=-x1*sin_tilt+z1*cos_tilt
         xp=(x1*cos_spin+y1*sin_spin)*cos_tilt-z1*sin_tilt
         yp=-1.0*x1*sin_spin+y1*cos_spin
         zp=(x1*cos_spin+y1*sin_spin)*sin_tilt+z1*cos_tilt
c
c
         x2=xp**2
         y2=yp**2
         z2=zp**2
         ar=sqrt(x2+y2+z2)
c
c        cartesian equivalent
c        !!!!!! change another expression!!
c         bmag=b0/ar**5
         bmag=-1.0*b0/ar**5
c        !!!!!! let dipole north at beginning
         dbx=-3.*bmag*xp*zp
         dby=-3.*bmag*yp*zp
         dbz=bmag*(x2+y2-2.*z2)
c
c        rotate b field back to coordinate space
c
        if(ar.gt.rearth-1.5) then
c
        bx0(i,j,k,m)=dbz*sin_tilt+dbx*cos_tilt
        by0(i,j,k,m)=dby*cos_spin-(dbz*cos_tilt-dbx*sin_tilt)*sin_spin
        bz0(i,j,k,m)=dby*sin_spin+(dbz*cos_tilt-dbx*sin_tilt)*cos_spin
c
c
c          bx0(i,j,k,m) =dbx*cos_tilt+dbz*sin_tilt
c          bz0(i,j,k,m)=-dbx*sin_tilt+dbz*cos_tilt

c
c          by0(i,j,k,m)=-3.*bmag*y1*zp
c          by0(i,j,k,m)=dby
c         back around z axis 
c
c          bx0aa=dbx*cos_spin-dby*sin_spin
c          by0aa=1.0*dbx*sin_spin+dby*cos_spin
c          bz0aa=dbz
c
c           bx0aa=dbx*cos_spin+dbz*sin_spin
c           bz0aa=-1.0*dbx*sin_spin+dbz*cos_spin
c           by0aa=dby
c          bx0aa=dbx
c          by0aa=dby
c          bz0aa=dbz
c
c
c         back around y axis
c
c          bx0(i,j,k,m)=bx0aa*cos_tilt+bz0aa*sin_tilt
c          bz0(i,j,k,m)=-1.0*bx0aa*sin_tilt+bz0aa*cos_tilt
c          by0(i,j,k,m)=by0aa
c
c           bx0(i,j,k,m)=bx0aa
c           by0(i,j,k,m)=by0aa
c           bz0(i,j,k,m)=bz0aa
c
c             bx0(i,j,k,m)=dbx*cos_tilt*cos_spin-dby*cos_tilt*sin_spin
c     +                   +dbz*sin_tilt
c             by0(i,j,k,m)=dbx*sin_spin+dby*cos_spin
c             bz0(i,j,k,m)=-dbx*sin_tilt*cos_spin+dby*sin_tilt*sin_spin
c     +                   +dbz*cos_tilt
c           bx0(i,j,k,m)=dbx
c           by0(i,j,k,m)=dby
c           bz0(i,j,k,m)=dbz
c
c
c           !!!!!!!!!!!!
c
        else
          bx0(i,j,k,m)=0.
          by0(i,j,k,m)=0.
          bz0(i,j,k,m)=0.
        endif
c
  110   continue
  120  continue
  130 continue
c
  140 continue
c
c     boundary conditions
c
      m=1
      do 210 n=1,numzero
c
c        get coords of point
c
         i=ijzero(1,n)
         j=ijzero(2,n)
         k=ijzero(3,n)
c
         bx0(i,j,k,m)=0.
         by0(i,j,k,m)=0.
         bz0(i,j,k,m)=0.
  210 continue
      return
      end
c
c     *************************************************
c
c     *************************************************
c
      subroutine vfldrot(abx,aby,abz,ax,ay,az,vfrac,v_rot,r_rot)
c
c     calculates magnetic field for a dipole
c
      common /gridding/grd_xmin(9),grd_xmax(9),grd_ymin(9),grd_ymax(9),
     +           grd_zmin(9),grd_zmax(9),xspac(9),ncore(9),nbndry(9),
     +           rx,ry,rz,xdip,ydip,zdip,rearth,b0,
     +           sin_tilt,cos_tilt
      common /uranus/tilt3,sin_tilt3,cos_tilt3,sin_spin,cos_spin,
     +           tilt4,sin_tilt4,cos_tilt4,spin,spinvectorx,
     +           spinvectory,spinvectorz
c
         cos_tilt98=cos(-97.9*.0174533)
         sin_tilt98=sin(-97.9*.0174533)
c    
c         x1=(ax-xdip)
c         y1=(ay-ydip)
c         z1=(az-zdip)
         x1=ax*re_equiv
         y1=ay*re_equiv
         z1=az*re_equiv
c         ar=sqrt(x1**2+y1**2)
c
c        real space to dipole space
c
         xp00=x1*cos_tilt98-z1*sin_tilt98
         yp00=y1
         zp00=x1*sin_tilt98+z1*cos_tilt98
c
c         superwholecoor=(spinvectorx*xp00+spinvectory*yp00
c     +   +spinvectorz*zp00)*(1.0-cos_spin)
c
c          xp=xp00*cos_spin-(spinvectory*zp00
c     +   -spinvectorz*yp00)*sin_spin+spinvectorx*superwholecoor
c          yp=yp00*cos_spin-(spinvectorz*xp00
c     +   -spinvectorx*zp00)*sin_spin+spinvectory*superwholecoor
c          zp=zp00*cos_spin-(spinvectorx*yp00
c     +   -spinvectory*xp00)*sin_spin+spinvectorz*superwholecoor
c        !!!!!!
         xp=xp00
         yp=yp00
         zp=zp00

         x2=xp**2
         y2=yp**2
         z2=zp**2
         rd=sqrt(x2+y2)
c         ar=sqrt(x2+y2+z2)
c
c        cartesian equivalent
c
         if(rd.lt.r_rot)then
           vfrac=1.
         else
           vfrac=((2.*r_rot**2)/(r_rot**2+ rd**2))**2
         endif
c
         dby=xp*v_rot*vfrac
         dbx=-yp*v_rot*vfrac
         dbz=0.0
c         bmag=b0/ar**5
c         dbx=-3.*bmag*xp*zp
c         aby=-3.*bmag*y1*zp
c         dbz=bmag*(x2+y2-2.*z2)
c
c        rotate b field back to coordinate space
c
c         superwhole=(spinvectorx*dbx+spinvectory*dby
c     +   +spinvectorz*dbz)*(1.0-cos_spin)
c
c          bx0aa=dbx*cos_spin+(spinvectory*dbz
c     +   -spinvectorz*dby)*sin_spin+spinvectorx*superwhole
c          by0aa=dby*cos_spin+(spinvectorz*dbx
c     +   -spinvectorx*dbz)*sin_spin+spinvectory*superwhole
c          bz0aa=dbz*cos_spin+(spinvectorx*dby
c     +   -spinvectory*dbx)*sin_spin+spinvectorz*superwhole
c         !!!!!!
          bx0aa=dbx
          by0aa=dby
          bz0aa=dbz
c
c
          abx=bx0aa*cos_tilt98+bz0aa*sin_tilt98
          abz=-bx0aa*sin_tilt98+bz0aa*cos_tilt98
          aby=by0aa
c
c       add in mirror dipole magetic field
c
c        xdipm=2.*grd_xmin -xdip
c        xm1=(ax-xdipm)
c        ym1=(ay-ydip)
c        zm1=(az-zdip)
c
c        xm2=(ax-xdipm)**2
c        ym2=(ay-ydip)**2
c        zm2=(az-zdip)**2
c        amr=sqrt(xm2+y2+z2)
c
c        bmag=b0/amr**5
c        abx=abx-3.*bmag*xm1*z1
c        aby=aby-3.*bmag*y1*z1
c        abz=abz+bmag*(xm2+y2-2.*z2)
c
c
c        find equivalent spherical cordinates
c
c
c         cost=(z1)/ar
c         sint=sqrt(x2+y2)/ar
c         cosp=x1/(ar*sint)
c         sinp=y1/(ar*sint)
c
c        br=-2.*b0*cost/ar**3
c        bt=-b0*sint/ar**3
c
c        bx=br*sint*cosp+bt*cost*cosp
c        by=br*sint*sinp+bt*cost*sinp
c        bz=br*cost-bt*sint      
c
      return
      end
c
c     **************************************************************
c
c      subroutine bndry_inner_rot(qrho,qpres,rmassq,
c     +       hrho,hpres,rmassh,
c     +       orho,opres,rmasso,epres,gamma1,eerg,
c     +       ngrd,nx,ny,nz,parm_srf,parm_mid,
c     +       ijsrf,numsrf,ijmid,nummid,
c     +       erho,ti_te,o_conc)
c      dimension qrho(nx,ny,nz,ngrd),qpres(nx,ny,nz,ngrd),
c     +     qpx(nx,ny,nz,ngrd),qpy(nx,ny,nz,ngrd),qpz(nx,ny,nz,ngrd),
c     +     hrho(nx,ny,nz,ngrd),hpres(nx,ny,nz,ngrd),
c     +     hpx(nx,ny,nz,ngrd),hpy(nx,ny,nz,ngrd),hpz(nx,ny,nz,ngrd),
c     +     orho(nx,ny,nz,ngrd),opres(nx,ny,nz,ngrd),
c     +     opx(nx,ny,nz,ngrd),opy(nx,ny,nz,ngrd),opz(nx,ny,nz,ngrd),
c     +     bx(nx,ny,nz,ngrd),by(nx,ny,nz,ngrd),bz(nx,ny,nz,ngrd),
c     +     epres(nx,ny,nz,ngrd)
c
      subroutine bndry_inner_rot(rmassq,rmassh,rmasso,
     +       gamma1,eerg,ngrd,nx,ny,nz,parm_srf,parm_mid,
     +       ijsrf,numsrf,ijmid,nummid,erho,ti_te,o_conc)
       dimension qrho(nx,ny,nz,ngrd),qpres(nx,ny,nz,ngrd),
     +     qpx(nx,ny,nz,ngrd),qpy(nx,ny,nz,ngrd),qpz(nx,ny,nz,ngrd),
     +     hrho(nx,ny,nz,ngrd),hpres(nx,ny,nz,ngrd),
     +     hpx(nx,ny,nz,ngrd),hpy(nx,ny,nz,ngrd),hpz(nx,ny,nz,ngrd),
     +     orho(nx,ny,nz,ngrd),opres(nx,ny,nz,ngrd),
     +     opx(nx,ny,nz,ngrd),opy(nx,ny,nz,ngrd),opz(nx,ny,nz,ngrd),
     +     bx(nx,ny,nz,ngrd),by(nx,ny,nz,ngrd),bz(nx,ny,nz,ngrd),
     +     epres(nx,ny,nz,ngrd)
c
c
      integer ijsrf(3,15000),ijmid(3,15000),ijzero(3,120000)
      real  parm_srf(7,15000),parm_mid(7,15000)
c
      common /gridding/grd_xmin(9),grd_xmax(9),grd_ymin(9),grd_ymax(9),
     +           grd_zmin(9),grd_zmax(9),xspac(9),ncore(9),nbndry(9),
     +           rx,ry,rz,xdip,ydip,zdip,rearth,b0,
     +           sin_tilt,cos_tilt
      common /uranus/tilt3,sin_tilt3,cos_tilt3,sin_spin,cos_spin,
     +           tilt4,sin_tilt4,cos_tilt4,spin,spinvectorx,
     +           spinvectory,spinvectorz
c     !!!!!! rotate the inner bndry with time fixed with B dipole
      d_min=0.001
      
      m=1
         dx=(grd_xmax(m)-grd_xmin(m))/(nx-1.)
         dy=(grd_ymax(m)-grd_ymin(m))/(ny-1.)
         dz=(grd_zmax(m)-grd_zmin(m))/(nz-1.)
c
c     !!!!!! part1: for parm_srf layer.
      do n=1,numsrf
         i=ijsrf(1,n)
         j=ijsrf(2,n)
         k=ijsrf(3,n)
          az=(grd_zmin(m)+dz*(k-1)-zdip)
          ay=grd_ymin(m)+dy*(j-1)-ydip
          ax=(grd_xmin(m)+dx*(i-1)-xdip)
c
c
         xp00=ax*cos_tilt-az*sin_tilt
         yp00=ay
         zp00=ax*sin_tilt+az*cos_tilt
c
        superwholecoor=(spinvectorx*xp00+spinvectory*yp00
     +   +spinvectorz*zp00)*(1.0-cos_spin)
c
          xp=xp00*cos_spin-(spinvectory*zp00
     +   -spinvectorz*yp00)*sin_spin+spinvectorx*superwholecoor
          yp=yp00*cos_spin-(spinvectorz*xp00
     +   -spinvectorx*zp00)*sin_spin+spinvectory*superwholecoor
          zp=zp00*cos_spin-(spinvectorx*yp00
     +   -spinvectory*xp00)*sin_spin+spinvectorz*superwholecoor
c
         ar=sqrt(xp**2+yp**2+zp**2)
c        for Jupiter top hat distribution
c        ar_iono=sqrt(xp**2+ay**2+0.75*zp**2)
c
c        for rearth = 6
c        ar_iono=sqrt(xp**2+ay**2+2.*zp**2)
c
c        for rearth = 10
c        ar_iono=sqrt(xp**2+ay**2+1.25*zp**2)
c        isotropic
         ar_iono=sqrt(xp**2+yp**2+zp**2)
c        ar_iono=amax1(0.0001,ar_iono)
         ar_iono=amax1(1.01*rearth,ar_iono)
c        ra=((ar_iono+0.5*rearth)/(1.5*rearth))**(-alpha_e)

         ra=exp(-(ar_iono-rearth)/(0.0994*rearth))

         zheight=amax1(1.,(zp**2+(1.5*rearth)**2)/(3.0*rearth)**2)
         ra=ra/zheight**2
         rho_iono=amax1(erho*ra,0.001)
c  
         r_equat=(ar**3+0.001)/(xp**2+yp**2+0.001)
         r_equat=amax1(r_equat,rearth)
         erg_sphere=eerg*(0.001+(rearth/r_equat)**4)
c        r_equat=amin1(r_equat,grd_xmax)
c        erg_sphere=eerg*r_equat/grd_xmax
c
         arho=amax1(rho_iono,d_min)
         aerg=amin1(abs(erg_sphere),0.01*arho)
c        !!!!!!
c         epres(i,j,k,m)=0.5*gamma1*aerg/ti_te
c
c         qrho(i,j,k,m)=rmassq*d_min
c         qpres(i,j,k,m)=0.
c         hrho(i,j,k,m)=arho*rmassh
c         hpres(i,j,k,m)=0.5*gamma1*aerg
c         orho(i,j,k,m)=rho_iono*rmasso*o_conc
c         opres(i,j,k,m)=0.5*gamma1*aerg*o_conc
c        !!!!!! because i,j,k=ijsrf(1/2/3,n),so order below correct
c           parm_srf(1,n)=qrho(i,j,k,m)
c           parm_srf(2,n)=hrho(i,j,k,m)
c           parm_srf(3,n)=orho(i,j,k,m)
c           parm_srf(4,n)=qpres(i,j,k,m)
c           parm_srf(5,n)=hpres(i,j,k,m)
c           parm_srf(6,n)=opres(i,j,k,m)
c           parm_srf(7,n)=epres(i,j,k,m)
c
c        !!!!!! because i,j,k=ijsrf(1/2/3,n),so order below correct
c           parm_srf(1,n)=200.0*rmassq*d_min
           parm_srf(1,n)=rmassq*d_min
           parm_srf(2,n)=arho*rmassh
           parm_srf(3,n)=rho_iono*rmasso*o_conc
           parm_srf(4,n)=0.
           parm_srf(5,n)=0.5*gamma1*aerg
           parm_srf(6,n)=0.5*gamma1*aerg*o_conc
           parm_srf(7,n)=0.5*gamma1*aerg/ti_te




c
      enddo
c
c     !!!!!! part2: for parm_mid layer.
      do n=1,nummid
         i=ijmid(1,n)
         j=ijmid(2,n)
         k=ijmid(3,n)
          az=(grd_zmin(m)+dz*(k-1)-zdip)
          ay=grd_ymin(m)+dy*(j-1)-ydip
          ax=(grd_xmin(m)+dx*(i-1)-xdip)
c
c
         xp00=ax*cos_tilt-az*sin_tilt
         yp00=ay
         zp00=ax*sin_tilt+az*cos_tilt
c
        superwholecoor=(spinvectorx*xp00+spinvectory*yp00
     +   +spinvectorz*zp00)*(1.0-cos_spin)
c
          xp=xp00*cos_spin-(spinvectory*zp00
     +   -spinvectorz*yp00)*sin_spin+spinvectorx*superwholecoor
          yp=yp00*cos_spin-(spinvectorz*xp00
     +   -spinvectorx*zp00)*sin_spin+spinvectory*superwholecoor
          zp=zp00*cos_spin-(spinvectorx*yp00
     +   -spinvectory*xp00)*sin_spin+spinvectorz*superwholecoor
c
         ar=sqrt(xp**2+yp**2+zp**2)
c        for Jupiter top hat distribution
c        ar_iono=sqrt(xp**2+ay**2+0.75*zp**2)
c
c        for rearth = 6
c        ar_iono=sqrt(xp**2+ay**2+2.*zp**2)
c
c        for rearth = 10
c        ar_iono=sqrt(xp**2+ay**2+1.25*zp**2)
c        isotropic
         ar_iono=sqrt(xp**2+yp**2+zp**2)
c        ar_iono=amax1(0.0001,ar_iono)
         ar_iono=amax1(1.01*rearth,ar_iono)
c        ra=((ar_iono+0.5*rearth)/(1.5*rearth))**(-alpha_e)

         ra=exp(-(ar_iono-rearth)/(0.0994*rearth))

         zheight=amax1(1.,(zp**2+(1.5*rearth)**2)/(3.0*rearth)**2)
         ra=ra/zheight**2
         rho_iono=amax1(erho*ra,0.001)
c  
         r_equat=(ar**3+0.001)/(xp**2+yp**2+0.001)
         r_equat=amax1(r_equat,rearth)
         erg_sphere=eerg*(0.001+(rearth/r_equat)**4)
c        r_equat=amin1(r_equat,grd_xmax)
c        erg_sphere=eerg*r_equat/grd_xmax
c
         arho=amax1(rho_iono,d_min)
         aerg=amin1(abs(erg_sphere),0.01*arho)
c        !!!!!!
c         epres(i,j,k,m)=0.5*gamma1*aerg/ti_te
c
c         qrho(i,j,k,m)=rmassq*d_min
c         qpres(i,j,k,m)=0.
c         hrho(i,j,k,m)=arho*rmassh
c         hpres(i,j,k,m)=0.5*gamma1*aerg
c         orho(i,j,k,m)=rho_iono*rmasso*o_conc
c         opres(i,j,k,m)=0.5*gamma1*aerg*o_conc
c        !!!!!! because i,j,k=ijmid(1/2/3,n),so order below correct
c           parm_mid(1,n)=qrho(i,j,k,m)
c           parm_mid(2,n)=hrho(i,j,k,m)
c           parm_mid(3,n)=orho(i,j,k,m)
c           parm_mid(4,n)=qpres(i,j,k,m)
c           parm_mid(5,n)=hpres(i,j,k,m)
c           parm_mid(6,n)=opres(i,j,k,m)
c           parm_mid(7,n)=epres(i,j,k,m)
c           parm_mid(1,n)=200.0*rmassq*d_min
           parm_mid(1,n)=rmassq*d_min
           parm_mid(2,n)=arho*rmassh
           parm_mid(3,n)=rho_iono*rmasso*o_conc
           parm_mid(4,n)=0.
           parm_mid(5,n)=0.5*gamma1*aerg
           parm_mid(6,n)=0.5*gamma1*aerg*o_conc
           parm_mid(7,n)=0.5*gamma1*aerg/ti_te
c
      enddo
c

c
      return
      end
c
c
c     **************************************************************
      subroutine rigid_IM_rot(rmassq,rmassh,rmasso,
     +       gamma1,eerg,ngrd,nx,ny,nz,qrho,hrho,orho,
     +       erho,qpres,hpres,opres,epres,ti_te,o_conc,ut,
     +       utdec_start,re_equiv,bx0,by0,bz0)
       dimension qrho(nx,ny,nz,ngrd),qpres(nx,ny,nz,ngrd),
     +     qpx(nx,ny,nz,ngrd),qpy(nx,ny,nz,ngrd),qpz(nx,ny,nz,ngrd),
     +     hrho(nx,ny,nz,ngrd),hpres(nx,ny,nz,ngrd),
     +     hpx(nx,ny,nz,ngrd),hpy(nx,ny,nz,ngrd),hpz(nx,ny,nz,ngrd),
     +     orho(nx,ny,nz,ngrd),opres(nx,ny,nz,ngrd),
     +     opx(nx,ny,nz,ngrd),opy(nx,ny,nz,ngrd),opz(nx,ny,nz,ngrd),
     +     bx(nx,ny,nz,ngrd),by(nx,ny,nz,ngrd),bz(nx,ny,nz,ngrd),
     +     epres(nx,ny,nz,ngrd),qrho_temp(nx,ny,nz,ngrd)
c
       dimension bx0(nx,ny,nz,ngrd),by0(nx,ny,nz,ngrd),
     +           bz0(nx,ny,nz,ngrd)
c       integer ijzero(3,120000)
c      integer ijsrf(3,15000),ijmid(3,15000),ijzero(3,120000)
c      real  parm_srf(7,15000),parm_mid(7,15000),parm_zero(7,120000)
c
      common /gridding/grd_xmin(9),grd_xmax(9),grd_ymin(9),grd_ymax(9),
     +           grd_zmin(9),grd_zmax(9),xspac(9),ncore(9),nbndry(9),
     +           rx,ry,rz,xdip,ydip,zdip,rearth,b0,
     +           sin_tilt,cos_tilt
      common /uranus/tilt3,sin_tilt3,cos_tilt3,sin_spin,cos_spin,
     +           tilt4,sin_tilt4,cos_tilt4,spin,spinvectorx,
     +           spinvectory,spinvectorz
c     !!!!!! rotate the inner bndry with time fixed with B dipole
      d_min=0.001
c     !!!!!! 01/19/2016 XC
c     relative decre rate at inner bndry, unit of ratio(e.g. %) per hr
c     !!!!!!
      if(ut.ge.utdec_start)then
          rate_dec=(1-500.0/12500.0)/17.24
      else
      rate_dec=0.0
      endif
c
      print *, "re_equiv in IM=",re_equiv
      print *, "rearth and erho,xspac(3)",rearth,erho,xspac(3)
         print *,"original qrho at tail (",qrho(35,32,32,5),
     +   qrho(39,32,32,5),qrho(49,32,32,5),qrho(50,32,32,5),
     +   qrho(51,32,32,5),qrho(52,32,32,5),qrho(61,32,32,5),")"
c
cc      print *,"numzero in last subroutine is:",numzero
cc      call mak_dipz(bx0,by0,bz0,nx,ny,nz,ngrd,ijzero,numzero)
ccc         print *,"pre nT",65.75*bx0(56,41,41,1),65.75*by0(56,41,41,1),
ccc     +            65.75*bz0(56,41,41,1)
ccc         print *,"nT is ",sqrt(bx0(56,41,41,1)**2+by0(56,41,41,1)**2
ccc     +                     +bz0(56,41,41,1)**2)*65.75
c
      do m=1,ngrd
         dx=(grd_xmax(m)-grd_xmin(m))/(nx-1.)
         dy=(grd_ymax(m)-grd_ymin(m))/(ny-1.)
         dz=(grd_zmax(m)-grd_zmin(m))/(nz-1.)
c
c     !!!!!! part1: for parm_srf layer.
c      do n=1,numsrf
c         i=ijsrf(1,n)
c         j=ijsrf(2,n)
c         k=ijsrf(3,n)
        do k=1,nz
            az=(grd_zmin(m)+dz*(k-1)-zdip)
         do j=1,ny
            ay=grd_ymin(m)+dy*(j-1)-ydip
          do i=1,nx
            ax=(grd_xmin(m)+dx*(i-1)-xdip)
c
          xp00=ax*cos_tilt-az*sin_tilt
          yp00=ay
          zp00=ax*sin_tilt+az*cos_tilt
c
        superwholecoor=(spinvectorx*xp00+spinvectory*yp00
     +   +spinvectorz*zp00)*(1.0-cos_spin)
c
          xp=xp00*cos_spin-(spinvectory*zp00
     +   -spinvectorz*yp00)*sin_spin+spinvectorx*superwholecoor
          yp=yp00*cos_spin-(spinvectorz*xp00
     +   -spinvectorx*zp00)*sin_spin+spinvectory*superwholecoor
          zp=zp00*cos_spin-(spinvectorx*yp00
     +   -spinvectory*xp00)*sin_spin+spinvectorz*superwholecoor
c
         ar=sqrt(xp**2+yp**2+zp**2)
         ar_temp=sqrt((ax+xdip)**2+(ay+ydip)**2+(az+zdip)**2)
c
c         if(((ar*re_equiv*xspac(m)).le.5.9).and.((ar*re_equiv*
c     +       xspac(m)).gt.3.12))then
c        !!!!!! Frozen inner magnetosphere for h & o speices
c        !!!!!! not use ar_temp instead of ar here.
         if(((ar*re_equiv).le.5.9).and.((ar*re_equiv
     +       ).gt.3.12))then
c        for Jupiter top hat distribution
c        ar_iono=sqrt(xp**2+ay**2+0.75*zp**2)
c
c        for rearth = 6
c        ar_iono=sqrt(xp**2+ay**2+2.*zp**2)
c
c        for rearth = 10
c        ar_iono=sqrt(xp**2+ay**2+1.25*zp**2)
c        isotropic
         ar_iono=sqrt(xp**2+yp**2+zp**2)
c        ar_iono=amax1(0.0001,ar_iono)
         ar_iono=amax1(1.01*rearth,ar_iono)
c        ra=((ar_iono+0.5*rearth)/(1.5*rearth))**(-alpha_e)

         ra=exp(-(ar_iono-rearth)/(0.0894*rearth))

         zheight=amax1(1.,(zp**2+(1.5*rearth)**2)/(3.0*rearth)**2)
         ra=ra/zheight**2
         rho_iono=amax1(erho*ra,0.001)
c  
         r_equat=(ar**3+0.001)/(xp**2+yp**2+0.001)
         r_equat=amax1(r_equat,rearth)
         erg_sphere=eerg*(0.001+(rearth/r_equat)**4)
c        r_equat=amin1(r_equat,grd_xmax)
c        erg_sphere=eerg*r_equat/grd_xmax
c
         arho=amax1(rho_iono,d_min)
         aerg=amin1(abs(erg_sphere),0.01*arho)
c        !!!!!!
c         epres(i,j,k,m)=0.5*gamma1*aerg/ti_te
c        !!!!!! epres,qrho,qpres do not update because of they need
c               to update physically outside of this subroutine.
c         qrho(i,j,k,m)=rmassq*d_min
c         qpres(i,j,k,m)=0.
         hrho(i,j,k,m)=arho*rmassh
         hpres(i,j,k,m)=0.5*gamma1*aerg
         orho(i,j,k,m)=rho_iono*rmasso*o_conc
         opres(i,j,k,m)=0.5*gamma1*aerg*o_conc

          endif
c
cc         if(((ar*re_equiv*xspac(m)).le.16.5).and.((ar*re_equiv*
cc     +       xspac(m)).ge.15.5))then

c        !!!!! Forzen mid inner magnetosphere for q, h & o speices
         if(((ar_temp*re_equiv).le.32.0).and.((ar_temp*re_equiv
     +   ).ge.16.0).and.((ax+xdip).gt.0.0).and.(sqrt(((ay+ydip)
     +   *re_equiv)**2+((az+zdip)*re_equiv)**2).le.1.2*abs(ax+xdip)
     +   *re_equiv))then
c
c        isotropic
         ar_iono=sqrt(xp**2+yp**2+zp**2)
c        ar_iono=amax1(0.0001,ar_iono)
         ar_iono=amax1(1.01*rearth,ar_iono)
c        ra=((ar_iono+0.5*rearth)/(1.5*rearth))**(-alpha_e)

         ra=exp(-(ar_iono-rearth)/(0.0894*rearth))

         zheight=amax1(1.,(zp**2+(1.5*rearth)**2)/(3.0*rearth)**2)
         ra=ra/zheight**2
         rho_iono=amax1(erho*ra,0.001)
c  
         r_equat=(ar**3+0.001)/(xp**2+yp**2+0.001)
         r_equat=amax1(r_equat,rearth)
         erg_sphere=eerg*(0.001+(rearth/r_equat)**4)
c
         arho=amax1(rho_iono,d_min)
         aerg=amin1(abs(erg_sphere),0.01*arho)
c        !!!!!!
c         epres(i,j,k,m)=0.5*gamma1*aerg/ti_te
c        !!!!!! epres,qrho,qpres do not update because of they need
c               to update physically outside of this subroutine.
c         qrho(i,j,k,m)=rmassq*d_min
c         qpres(i,j,k,m)=0.
cc         if((i.eq.nx).and.(j.eq.ny).and.(k.eq.nz).and.(m.eq.5))then
cc         print *,"original qrho at tail (",qrho(35,32,32,5),
cc     +   qrho(39,32,32,5),qrho(50,32,32,5),qrho(75,32,32,5),")"
cc         endif
         b0_temp=sqrt(bx0(i,j,k,m)**2+by0(i,j,k,m)**2+bz0(i,j,k,m)**2)
     +           *65.75
cc         print *,"b0_temp=",b0_temp,i,j,k,m
cc         if((m.eq.4).and.(j.eq.32).and.(k.eq.32))then
cc           print *,"i,j,k,m=",i,j,k,m
cc         print *,"nT is ",sqrt(bx0(47,32,32,1)**2+by0(47,32,32,1)**2
cc     +                     +bz0(47,32,32,1)**2)*65.75
cc         endif
         qrho_temp(i,j,k,m)=(b0_temp*1.0e-9)**2/(2.0*4.0
cc         qrho_temp(i,j,k,m)=(22836.0*1.0e-9/(ar*re_equiv)**3)**2/(2.0*4.0
     +       *3.1415926*1.0e-7)/(1.38*1.0e-23*50000.0)/(1.0e6)*0.7/2.0
cc          qrho_temp(i,j,k,m)=30000.0
         if(qrho(i,j,k,m).gt.qrho_temp(i,j,k,m))then
         qrho(i,j,k,m)=qrho_temp(i,j,k,m)
         endif
cc         if((i.eq.nx).and.(j.eq.ny).and.(k.eq.nz).and.(m.eq.5))then
cc         print *,"modify qrho in beta<1 (",qrho(35,32,32,5),
cc     +   qrho(39,32,32,5),qrho(50,32,32,5),qrho(75,32,32,5),")"
cc         endif
         if(qrho(i,j,k,m).lt.d_min)then
         qrho(i,j,k,m)=d_min
         endif
cc         if((i.eq.nx).and.(j.eq.ny).and.(k.eq.nz).and.(m.eq.5))then
cc         print *,'modify qrho in d_min: (',qrho(35,32,32,5),
cc     +   qrho(39,32,32,5),qrho(50,32,32,5),qrho(75,32,32,5),")"
cc         endif
cc
cc         hrho(i,j,k,m)=arho*rmassh
cc         hpres(i,j,k,m)=0.5*gamma1*aerg
cc         orho(i,j,k,m)=rho_iono*rmasso*o_conc
cc         opres(i,j,k,m)=0.5*gamma1*aerg*o_conc
c
         endif
c
          enddo
         enddo
        enddo
       enddo
         print *,'modify qrho in d_min: (',qrho(35,32,32,5),
     +   qrho(39,32,32,5),qrho(49,32,32,5),qrho(50,32,32,5),
     +   qrho(51,32,32,5),qrho(52,32,32,5),qrho(61,32,32,5),")"

c
      return
      end
c
c     **************************************************************
c     **************************************************************
