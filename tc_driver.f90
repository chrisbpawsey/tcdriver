!************************************************************************
! Program to calculate the terrain effect over Australia from a
! direct-access file created from the AUSLIG 9-second DEM tiles.
! This is done by FFT over strips, writing out the result in one
!  regular grid.
!
! However, for small areas the computations can be performed on the
!  whole grid - see below for details.
!
! The output can be ascii or binary direct-access format (useful for
!  huge files).
!
! Area divided into latitude strips of (centre + 2*border) degrees;
! TCs calculated by FFT over the whole strip, but only written to
!  the central parallels to avoid edge effects;
! Strip is advanced by centre-degrees northward and FFT repeated.
!
!    +---------------+           +---------------+
!    |               |           |               |
!    |               |           |               |
!    |               |           |- - - - - - - -|
!    |               |           |               | border
!    |               |           |---------------|
!    |- - - - - - - -|           |///////////////| centre
!    |               |           |///////////////|
!    |---------------|   ====>   |---------------|
!    |///////////////|           |               | border
!    |///////////////|           |- - - - - - - -|
!    |---------------|           |               |
!    |               |           |               |
!    |- - - - - - - -|           |               |
!    |               |           |               |
!    |               |           |               |
!    +---------------+           +---------------+
!
! The border width MUST be greater than or equal to 0.5 degrees.
!  This is because the integration kernel has a radius of 50 km.
! The central strip width may have any value, limited by computer RAM.
!
! If TCs are to be calculated over an area in one hit (ie. not strip-wise)
!  then set the variable "border" to 0.d0, and variable "centre" to the
!  latitude extent of the area in degrees. Be sure to change the parameter
!  "jstrip0" to reflect this!
!
! Parameters :
!  grid = grid cell size in degrees
!  istrip = strip length, in no. cells longitude
!  centre = central strip width in degrees latitude
!  border = border width to central strip, in degrees latitude
!  jstrip = total strip width, in no. cells latitude
!
! written by Jon Kirby, Curtin University, July 1997.
!***********************************************************************

implicit double precision (a-h,o-z)
      parameter (istrip0=16289,jstrip0=600)
      real height(istrip0*jstrip0)
      double precision tc(istrip0*jstrip0),z(istrip0)
      integer*4 iz,izero
      character*80 filein,fileout
      character*1 out
      common /a1/ height
      common /a2/ tc
      common /c1/ G,density,jstrip,istrip,K_term
      common /c8/ dx_km,dy_km
      common /c9/ jbord,intrad

! while the "length" of the strip is fixed at the longitude extent of
! the area, the strip "width" may vary.
! it is recommended that the border width be kept at half a degree
      centre=0.5d0
      border=0.5d0

! set integration radius (in degrees):
      radint=0.5d0

! report file :
      open(10,file='tcdriver.rep')

!-----------------------------------------------------------------------
! input :

      write(*,6000)
 6000 format(/'calculate TCs from ...'/    &
           & '(1) national 9" DEM'/        &
           & '(2) national 27" DEM'/       &
           & '(3) Tasmanian 9" test-area'/ &
           & '(4) Tasmanian 27" test-area'/ &
           & '(5) national 1 second dem'/  &
           & '(6) user-defined area?')
      read(*,*) iopt
      write(*,*) iopt
!iopt=3

! AUSLIG 9" DEM :
      if (iopt.eq.1) then
        filein='ausdem9s.bin'
        xwa=112.91875d0
        xea=153.63875d0
        ysa=-43.74125d0
        yna=-9.00125d0
        nii=16289
        njj=13897
        grid=2.5d-3          ! grid spacing in degrees
        dx_km=0.25d0         ! grid spacing in km
        dy_km=0.25d0

! AUSLIG 27" DEM :
      else if (iopt.eq.2) then
        filein='ausdem27s.bin'
        xwa=112.92125d0
        xea=153.63875d0
        ysa=-43.74375d0
        yna=-9.00375d0
        nii=5430
        njj=4633
        grid=7.5d-3          ! grid spacing in degrees
        dx_km=0.75d0         ! grid spacing in km
        dy_km=0.75d0

! TAS test DEM 9"
      else if (iopt.eq.3) then
        filein='tasdem9s.bin'
        xwa=139.99875d0
        xea=149.99875d0
        ysa=-43.74125d0
        yna=-35.99875d0
        nii=4001
        njj=3098
        grid=2.5d-3          ! grid spacing in degrees
        dx_km=0.25d0         ! grid spacing in km
        dy_km=0.25d0

! TAS test DEM 27"
      else if (iopt.eq.4) then
        filein='tasdem27s.bin'
        xwa=139.99625d0
        xea=150.00125d0
        ysa=-43.74375d0
        yna=-36.00375d0
        nii=1334
        njj=1033
        grid=7.5d-3          ! grid spacing in degrees
        dx_km=0.75d0         ! grid spacing in km
        dy_km=0.75d0

! 1 second DEM"
      else if (iopt.eq.5) then
        filein='ausdem1s.bin'
        xwa=112.9998611111000d0
        xea=153.9995833337211d0
        ysa=-44.0001388891040d0
        yna=-10.0004166665510d0
        nii=151200
        njj=126000
        grid=(xea-xwa)/nii   ! grid spacing in degrees
        dx_km=grid*100.d0    ! at -27 latitude
        dy_km=dx_km

! user-defined option
      else
        write(*,600)
  600   format(/'enter input (d/a binary) filename :')
        read(*,*) filein
        write(*,6001)
 6001   format(/'enter longitude (dec. degs) of westernmost then &
             & easternmost grid-cell centre :')
        read(*,*) xwa,xea
        write(*,*) xwa,xea
        write(*,6002)
 6002   format(/'enter latitude (dec. degs) of northernmost then &
             & southernmost grid-cell centre :')
        read(*,*) yna,ysa
        write(*,*) yna,ysa
        write(*,6003)
 6003   format(/'enter grid spacing (seconds) :')
        read(*,*) grid
        write(*,*) grid
        grid=grid/3.6d2
        dx_km=grid*100.d0    ! at -27 latitude
        dy_km=dx_km

      endif

      open(1,file=filein,status='old',form='unformatted',  &
           & access='direct',recl=4)

!-----------------------------------------------------------------------
! strip node parameters :

      rgrid=1.d0/grid                ! reciprocal grid spacing

! integration radius (in pixels):
      intrad=radint*rgrid

      istrip=nii
      jstrip=(centre+2*border)*rgrid

      jbord=border*rgrid
      jcen=centre*rgrid
      jstart=jbord+1
      jend=jstrip-jbord

! info :

      write(*,601) xwa,xea,yna,ysa,nii,njj,istrip,jstrip
      write(10,601) xwa,xea,yna,ysa,nii,njj,istrip,jstrip
      call flush(10)
  601 format(/'working grid is area '/3x,'lon : ',f9.5,' to ',f9.5/3x,  &
           & 'lat : ',f9.5,' to ',f9.5,4x,'or  ',i7,' x ',i7,' gridcells'//  &
           & 'terrain corrections calculated in strips of size ',i7,' x ',  &
           & i7,' gridcells'/)

! ensure parameters are large enough :
      if (istrip.ne.istrip0) then
        write(*,602) istrip
        write(10,602) istrip
        call flush(10)
  602   format(/'change parameter istrip0 to equal ',i6/)
        stop
      endif
      if (jstrip.ne.jstrip0) then
        write(*,6021) jstrip
        write(10,6021) jstrip
        call flush(10)
 6021   format(/'change parameter jstrip0 to equal ',i6/)
        stop
      endif

      if (jcen.gt.njj) then
        write(*,603)
        write(10,603)
        call flush(10)
 603     format(/'central strip wider than dataset ...  &
              &  change the variable "centre"'/)
        stop
      endif

      write(*,604)
  604 format(/'output format :'/' a) ascii (.ary),'/  &
           &      ' b) binary (.bin)?')
      read(*,'(a)') out
      write(*,'(a)') out

      write(*,605)
  605 format(/'enter terrain effect o/p filename :')
      read(*,'(A)') fileout
      write(*,'(A)') fileout

!-----------------------------------------------------------------------
! open binary file for TCs and initialise :

      if (out.eq.'a') then
        open(8,status='scratch',form='unformatted',access='direct', &
             & recl=4)
      else
        open(8,file=fileout,form='unformatted',access='direct',recl=4)
      endif
      write(8,rec=1) izero
      write(8,rec=nii*njj) izero

!-----------------------------------------------------------------------
! calculate TC at each node along a parallel, by FFT over a band the
!  length of the continent, and parameter "jstrip" wide.
! working from north to south
! split into 3 parts to overcome array assignments and the use of
!  "if" statements.

      write(*,606)
      write(10,606)
      call flush(10)
  606 format(//'tcdriver : commencing terrain correction &
           & computation ...'/)

! first part ;
! place zeros where strip includes areas outside dem tile

      do 57 jl=1,jbord,jcen       ! jl records latitude (=1,njj)
       write(*,609) jl,jl+jcen-1
       write(10,609) jl,jl+jcen-1
       call flush(10)
  609  format(/'computing strips ',i5,'  to  ',i5)

       do 52 jd=1,jstrip          ! jd records pos. within strip
        jj=jl-jbord+jd-1          ! jj records latitude above & below strip
        if (jj.lt.1) then         ! set heights outside DEM to zero
          do 50 id=1,istrip
           ih=id+(jd-1)*istrip
           height(ih)=0.
   50     continue
        else
          do 51 id=1,istrip           
           ir=id+(jj-1)*nii
           read(1,rec=ir) izh
           if (izh.lt.0) izh=0    ! set -ve heights and -9999 flags to zero
           ih=id+(jd-1)*istrip
           height(ih)=izh*1.e-3
   51     continue
        endif
   52  continue

       call tc2dftpl

       do 54 j=jstart,jend        ! write TCs in central band only
        jj=jl-jbord+j-1
        do 53 i=1,istrip
         ih=i+(j-1)*istrip
         izt=nint(tc(ih)*1.d3)
         if (izt.lt.0) izt=0      ! ensure TCs +ve only
         ir=i+(jj-1)*nii
         write(8,rec=ir) izt
         tc(ih)=0.d0              ! reset TC array in subroutine to zero
   53   continue
   54  continue

   57 continue

      jl2=jl          ! jl2 records the parallel last computed, since
                      ! jbord may not be an exact multiple of jcen

!****************
! second part :

      do 67 jl=jl2,njj-(jcen+jbord)+1,jcen
       write(*,609) jl,jl+jcen-1
       write(10,609) jl,jl+jcen-1
       call flush(10)

       do 61 jd=1,jstrip
        jj=jl-jbord+jd-1
        do 60 id=1,istrip
         ir=id+(jj-1)*nii 
         read(1,rec=ir) izh
         if (izh.lt.0) izh=0    ! set -ve heights and -9999 flags to zero
         ih=id+(jd-1)*istrip
         height(ih)=izh*1.e-3
   60   continue
   61  continue

       call tc2dftpl

       do 64 j=jstart,jend        ! write TCs in central band only
        jj=jl-jbord+j-1
        do 63 i=1,istrip
         ih=i+(j-1)*istrip
         izt=nint(tc(ih)*1.d3)
         if (izt.lt.0) izt=0      ! ensure TCs +ve only
         ir=i+(jj-1)*nii
         write(8,rec=ir) izt
         tc(ih)=0.d0              ! reset TC array in subroutine to zero
   63   continue
   64  continue

   67 continue

      jl2=jl              ! jl2 records the parallel last computed

!!!!!!!!!!!!!!!!!
! third part :

      do 77 jl=jl2,njj,jcen
       write(*,609) jl,jl+jcen-1
       write(10,609) jl,jl+jcen-1
       call flush(10)

       do 72 jd=1,jstrip
        jj=jl-jbord+jd-1
        if (jj.le.njj) then
          do 70 id=1,istrip
           ir=id+(jj-1)*nii
           read(1,rec=ir) izh
           if (izh.lt.0) izh=0     ! set -ve heights and -9999 flags to zero
           ih=id+(jd-1)*istrip
           height(ih)=izh*1.e-3
   70     continue
        else                       ! set heights outside DEM to zero
          do 71 id=1,istrip
           ih=id+(jd-1)*istrip
           height(ih)=0.
   71     continue
        endif
   72  continue

       call tc2dftpl

       do 74 j=jstart,jend        ! write TCs in central band only
        jj=jl-jbord+j-1
        if (jj.le.njj) then
          do 73 i=1,istrip
           ih=i+(j-1)*istrip
           izt=nint(tc(ih)*1.d3)
           if (izt.lt.0) izt=0       ! ensure TCs +ve only
           ir=i+(jj-1)*nii
           write(8,rec=ir) izt
           tc(ih)=0.d0              ! reset TC array in subroutine to zero
   73     continue
        endif
   74  continue

   77 continue

      close(1)

!-----------------------------------------------------------------------
! output :

! ascii array output :

      if (out.eq.'a') then
        write(*,610) fileout
        write(10,610) fileout
        call flush(10)
  610   format(/'output of terrain correction in milligals to ascii &
             & array format file '/' => ',a80)

        open(2,file=fileout)
        write(2,201) njj,nii,yna,xwa,grid*60.d0,grid*60.d0

        do 85, jj=1,njj
         do 80, ii=1,nii
          ir=ii+(jj-1)*nii
          read(8,rec=ir) iz
          z(ii)=iz*1.d-3
   80    continue
         write(2,202) (z(ii),ii=1,nii)
   85   continue

        close(2)
  201   format(2i6,2f10.5,2f7.3)
  202   format(10f9.3)

! binary file output :

      else if (out.eq.'b') then
        write(*,611) fileout
        write(10,611) fileout
        call flush(10)
  611   format(/'terrain correction in microgals written to direct- &
             & access binary file '/' => ',a80)
      endif

!     write(*,612)
      write(10,612)
      call flush(10)
  612 format(///'END of tcdriver'///)


      close(8)
      close(10)
      stop
      end

