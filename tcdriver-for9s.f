
      program tcdriver

************************************************************************
* Program to calculate the terrain effect over Australia from a
*  direct-access file created from the AUSLIG 9-second DEM tiles.
* This is done by FFT over strips, writing out the result in one
*  regular grid.
*
* However, for small areas the computations can be performed on the
*  whole grid - see below for details.
*
* The output can be ascii or binary direct-access format (useful for
*  huge files).
*
* Area divided into latitude strips of (centre + 2*border) degrees;
* TCs calculated by FFT over the whole strip, but only written to
*  the central parallels to avoid edge effects;
* Strip is advanced by centre-degrees northward and FFT repeated.
*
*    +---------------+           +---------------+
*    |               |           |               |
*    |               |           |               |
*    |               |           |- - - - - - - -|
*    |               |           |               | border
*    |               |           |---------------|
*    |- - - - - - - -|           |///////////////| centre
*    |               |           |///////////////|
*    |---------------|   ====>   |---------------|
*    |///////////////|           |               | border
*    |///////////////|           |- - - - - - - -|
*    |---------------|           |               |
*    |               |           |               |
*    |- - - - - - - -|           |               |
*    |               |           |               |
*    |               |           |               |
*    +---------------+           +---------------+
*
* The border width MUST be greater than or equal to 0.5 degrees.
*  This is because the integration kernel has a radius of 50 km.
* The central strip width may have any value, limited by computer RAM.
*
* If TCs are to be calculated over an area in one hit (ie. not strip-wise)
*  then set the variable "border" to 0.d0, and variable "centre" to the
*  latitude extent of the area in degrees. Be sure to change the parameter
*  "jstrip0" to reflect this!
*
* Parameters :
*  grid = grid cell size in degrees
*  istrip = strip length, in no. cells longitude
*  centre = central strip width in degrees latitude
*  border = border width to central strip, in degrees latitude
*  jstrip = total strip width, in no. cells latitude
*
* written by Jon Kirby, Curtin University, July 1997.
************************************************************************

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

* while the "length" of the strip is fixed at the longitude extent of
* the area, the strip "width" may vary.
* it is recommended that the border width be kept at half a degree
      centre=0.5d0
      border=0.5d0

* set integration radius (in degrees):
      radint=0.5d0

* report file :
      open(10,file='tcdriver.rep')

*-----------------------------------------------------------------------
* input :

      write(*,6000)
 6000 format(/'calculate TCs from ...'/
     # '(1) national 9" DEM'/
     # '(2) national 27" DEM'/
     # '(3) Tasmanian 9" test-area'/
     # '(4) Tasmanian 27" test-area'/
     # '(5) national 1 second dem'/
     # '(6) user-defined area?')
      read(*,*) iopt
      write(*,*) iopt
!iopt=3

* AUSLIG 9" DEM :
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

* AUSLIG 27" DEM :
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

* TAS test DEM 9"
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

* TAS test DEM 27"
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

* 1 second DEM"
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

* user-defined option
      else
        write(*,600)
  600   format(/'enter input (d/a binary) filename :')
        read(*,*) filein
        write(*,6001)
 6001   format(/'enter longitude (dec. degs) of westernmost then 
     #   easternmost grid-cell centre :')
        read(*,*) xwa,xea
        write(*,*) xwa,xea
        write(*,6002)
 6002   format(/'enter latitude (dec. degs) of northernmost then 
     #   southernmost grid-cell centre :')
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

      open(1,file=filein,status='old',form='unformatted',
     #                                         access='direct',recl=4)

*-----------------------------------------------------------------------
* strip node parameters :

      rgrid=1.d0/grid                ! reciprocal grid spacing

* integration radius (in pixels):
      intrad=radint*rgrid

      istrip=nii
      jstrip=(centre+2*border)*rgrid

      jbord=border*rgrid
      jcen=centre*rgrid
      jstart=jbord+1
      jend=jstrip-jbord

* info :

      write(*,601) xwa,xea,yna,ysa,nii,njj,istrip,jstrip
      write(10,601) xwa,xea,yna,ysa,nii,njj,istrip,jstrip
      call flush(10)
  601 format(/'working grid is area '/3x,'lon : ',f9.5,' to ',f9.5/3x,
     #  'lat : ',f9.5,' to ',f9.5,4x,'or  ',i7,' x ',i7,' gridcells'//
     #  'terrain corrections calculated in strips of size ',i7,' x ',
     #  i7,' gridcells'/)

* ensure parameters are large enough :
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
  603   format(/'central strip wider than dataset ... 
     #          change the variable "centre"'/)
        stop
      endif

      write(*,604)
  604 format(/'output format :'/' a) ascii (.ary),'/
     #                          ' b) binary (.bin)?')
      read(*,'(a)') out
      write(*,'(a)') out

      write(*,605)
  605 format(/'enter terrain effect o/p filename :')
      read(*,'(A)') fileout
      write(*,'(A)') fileout

*-----------------------------------------------------------------------
* open binary file for TCs and initialise :

      if (out.eq.'a') then
        open(8,status='scratch',form='unformatted',access='direct',
     #                                                        recl=4)
      else
        open(8,file=fileout,form='unformatted',access='direct',recl=4)
      endif
      write(8,rec=1) izero
      write(8,rec=nii*njj) izero

*-----------------------------------------------------------------------
* calculate TC at each node along a parallel, by FFT over a band the
*  length of the continent, and parameter "jstrip" wide.
* working from north to south
* split into 3 parts to overcome array assignments and the use of
*  "if" statements.

      write(*,606)
      write(10,606)
      call flush(10)
  606 format(//'tcdriver : commencing terrain correction 
     # computation ...'/)

* first part ;
* place zeros where strip includes areas outside dem tile

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

*****************
* second part :

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

*****************
* third part :

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

*-----------------------------------------------------------------------
* output :

* ascii array output :

      if (out.eq.'a') then
        write(*,610) fileout
        write(10,610) fileout
        call flush(10)
  610   format(/'output of terrain correction in milligals to ascii 
     #          array format file '/' => ',a80)

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

* binary file output :

      else if (out.eq.'b') then
        write(*,611) fileout
        write(10,611) fileout
        call flush(10)
  611   format(/'terrain correction in microgals written to direct-
     #          access binary file '/' => ',a80)
      endif

c     write(*,612)
      write(10,612)
      call flush(10)
  612 format(///'END of tcdriver'///)


      close(8)
      close(10)
      stop
      end


****************************************************************************

c     ******************************************************************
c     *                         tc2dftpl                               *
c     *   terrain correction via 2 dimensional fast fourier transform  *
c     *    with either mass-prism topographic model                    *
c     *          or    mass-line topographic model.                    *
c     *                                                                *
c     *    programmer:   yecai li,                                     *
* changed from a main program to a subroutine by
*          jon kirby, curtin university, 1997
*
c     parameter
c     max_rowxclm: the mulplication of the maximum number of rows
c                                   by the maximum number of cloumns.
c                  this parameter has to be changed, if necessary,
c                                                   in each subroutine.
c     maxsubs:     the maximum size of the complex two-dimensional array
c                  after 100% zero-padding.

c     subroutines called
c     dtdz_2d2:     compute two convolutions simultaneously.
c     s2dk_txtytz:  compute the kernel function for the first-order
c                   gradients (tx, ty & tz) of the topogrational 
c                   gravitational potential.
c     tcmp0t:       compute the zero-order term, i.e., the attraction
c                   of the temple plate.
c     iftafb:       compute the fourier  transform of b;
c                   mutiply a with f{b} and 
c                   do the inverse fourier transform of axf{b}

      subroutine tc2dftpl
      parameter (max_rowxclm=16289*600, maxsubs=4*max_rowxclm)
      implicit double precision (a-h,o-z)
      real height(max_rowxclm)
      double precision result(max_rowxclm)
      complex fkernel(maxsubs),fheight(maxsubs)

      common /a1/ height
      common /a2/ result
      common /a3/ fheight
      common /a4/ fkernel
      common /c1/ g,density,maxi,maxj,k_term
      common /c2/ h_removed,i0pad,j0pad,dx_m,dy_m
      common /c3/ i_begin,i_end,j_begin,j_end
      common /c4/ n0bi,n0bj
      common /c5/ iradius,jradius
      common /c6/ xradius,yradius
      common /c7/ iwf
      common /c8/ dx_km,dy_km
      common /c9/ jbord,intrad

* position identifier:
      nij(max_colm,i_row,j_colm) = (i_row-1)*max_colm + j_colm

      g = 6.672d-8
      density = 2.67d0

      dx_m = dx_km
      dy_m = dy_km

************************************************************************
* compute the statistics of the heights
* convert heights from metres to km

      hmean=0.d0
      hstd=0.d0

      do jj=1,maxi*maxj
       height(jj) = height(jj)*1.e-3
       hmean=hmean+height(jj)
       hstd=hstd+height(jj)*height(jj)
      enddo

      hmean=hmean/(1.d0*maxi*maxj)
      hstd=dsqrt(hstd/(1.d0*maxi*maxj)-hmean*hmean)

*********************************************************************
* take the parameter alpha in the kernel function as
*  sigma(h)/sqrt(2)

      h_removed = hstd/dsqrt(2.d0)

*********************************************************************
* do the computation from ...
* the zero-order term (k_01 = 0);
* or 1st-order term (k_01 = 1):

      k_01 = 1

*********************************************************************
* specify how many terms to be computed (in the denominator) :
*  k_term = 1: (x*x+y*y)**(3/2) (consider the first term only),
*  k_term = 2: (x*x+y*y)**(5/2) (consider the 1st & 2nd terms).
*  k_term = 3: (x*x+y*y)**(7/2) (the 1st, 2nd & 3rd terms).

      k_term = 1

*********************************************************************
* compute the kernel function ...
*  over whole computation area with 100% zero padding (ikernel=1);
*  within a limited area with either 100% or 50% padding (ikernel=2)

      ikernel=2

* compute over whole area [default] :
      int_ker = 1

* set integral radius = the length of the whole area
      xradius = maxj*dx_km
      yradius = maxi*dy_km

      iradius = maxi
      jradius = maxj

* number of zeros padded to each side :
      n0bi = maxi/2
      n0bj = maxj/2

* new array dimensions after 100% zeros are padded (ie, twice the size)
      i0pad = maxi*2
      j0pad = maxj*2

      i_begin = 1
      i_end = i0pad
      j_begin = 1
      j_end = j0pad

*--------------------------------------------------------------------
* compute over limited area :

      if (ikernel.eq.2) then

        int_ker = 2

* set integration radius to half the strip dimensions :
c       jradius=maxj/2
c       iradius=maxi/2
        jradius=intrad
        iradius=intrad
        xradius = jradius*dx_km
        yradius = iradius*dy_km

* option to zero pad this kernel :
*  ipad = 1: 100% of the number of the heights,
*       = 2: 50% of (the number of height + the size of kernel)

        ipad=1

        if (ipad.eq.2) then
          n0bi = iradius + 1
          n0bj = jradius + 1
          i0pad = maxi + n0bi*2
          j0pad = maxj + n0bj*2
        endif

* compute limits of the integral of kernel function.
        i_begin = max0(1,i0pad/2+1-iradius)
        i_end   = min0(i0pad,i0pad/2+1+iradius)
        j_begin = max0(1,j0pad/2+1-jradius)
        j_end   = min0(j0pad,j0pad/2+1+jradius)

      endif

*********************************************************************
* topographic model may be ...
*  mass-prism (k_tm=1) or mass-line (k_tm=2) :

      k_tm = 1

*********************************************************************
* include integration weighting factor (iwf=1); don't (iwf=2)
* in subroutine tc2dkmpm

      iwf=1

*********************************************************************
* compute the convolutions two simultaneously

      call dtdz_2d2(k_tm,k_01)

*********************************************************************
      write(*,209)
      write(10,209)
      call flush(10)
209   format('terrain computations have been done successfully'/)


      return
      end


c     *************************************************************
c     *                        dtdz_2d2                           *
c     *    computes the vertical components of the first order    *
c     *    gradients of the topographic gravitational potential   *
c     *          by doing two convolutions simultaneously         *
c     *                                                           *
c     *    programmer:   yecai li,                                *
c     *                  dept of geomatics engineering            *
c     *                  the university of calgary                *
c     *                  june 1993                                *
c     *                                                           *
c     *************************************************************

      subroutine dtdz_2d2(k_tm,k_01)

      parameter (max_rowxclm=16289*600, maxsubs=4*max_rowxclm)
      implicit double precision (a-h,o-z)
      real height(max_rowxclm)
      double precision result(max_rowxclm)
      complex fkernel(maxsubs),fheight(maxsubs)
      common /a1/ height
      common /a2/ result
      common /a3/ fheight
      common /a4/ fkernel
      common /c1/ g,density,maxi,maxj,k_term
      common /c2/ h_removed,i0pad,j0pad,dx_m,dy_m
      common /c3/ i_begin,i_end,j_begin,j_end
      common /c4/ n0bi,n0bj

      nij(max_colm,i_row,j_colm) = (i_row-1)*max_colm + j_colm

c     the unit of g*density is 1/sec**2
c     when the unit of height is meter and g*density*1.e5
c     the unit of tz is mgal
c     c_mgal = g*density*1.d5*i0pad*j0pad
      c_mgal = g*density*1.d8*i0pad*j0pad
      alpha_sq = h_removed * h_removed


c      compute the 1st- and/or the 2nd-order term
c         ---------------------------------------------------
c         |   compute the kernel function 1./(xx+yy)**1.5   |
c         |                  and fft                        |
c         ---------------------------------------------------

* mass-prism :
      call tc2dkmpm(1)

c     ---------------------------------------
c     correct for the origin shift and
c     move f{kernel} from fheight to fkernel.
c     ---------------------------------------

      do i = 1,i0pad
       do j = 1,j0pad
        kk = nij(j0pad,i,j)
        fkernel(kk) = fheight(kk)*(-1)**(i+j)
       enddo
      enddo

c     --------------------------------------------------
c     compute the convolution h**2 with 1/(xx+yy)**(3/2)
c     compute the convolution h with 1/(xx+yy)**(3/2)
c     --------------------------------------------------

      call iftafb(2,1,1,maxi,maxj)

      c_11 = 0.5d0*c_mgal
      c_12 = -2.d0*c_11
      do i=1,maxi
       do j=1,maxj
        jj = nij(maxj,i,j)
        kk = nij(j0pad,i+n0bi,j+n0bj)
        zp = height(jj)
        term1 = c_11*real(fheight(kk))+c_12*zp*aimag(fheight(kk))
        result(jj) = result(jj) + term1
       enddo
      enddo

c     -----------------------------------------------
c     compute the convolution 1 with 1/(xx+yy)**(3/2)
c     -----------------------------------------------

      call iftafb(0,-1,1,maxi,maxj)

      do i=1,maxi
       do j=1,maxj
        jj = nij(maxj,i,j)
        kk = nij(j0pad,i+n0bi,j+n0bj)
        zp = height(jj)
        if(k_01.eq.0) then
          c_13 = c_11*(zp*zp - alpha_sq)
        else
          c_13 = c_11*zp*zp
        endif
        term1 = c_13*real(fheight(kk))
        result(jj) = result(jj) + term1
       enddo
      enddo

      return
      end


c     ***************************************************
c     *                 tc2dkmpm                        *
c     *    terrain correction 2 dimensional kernel      *
c     *     with mass-prism topographic model           *
c     *    s  : simple summation formulas are used      *
c     *    2dk: the 2 ddimensional kernel function      *
c     *                                                 *
c     *    programmer:   yecai li,                      *
c     *                  dept of geomatics engineering  *
c     *                  the university of calgary      *
c     *                  june 1993                      *
c     *                                                 *
c     *    modified by j f kirby to include weighting   *
c     *    function and zeros outside cap               *
c     *                                                 *
c     ***************************************************
c     k_term = 1, compute the kernel of the 1st term,i.e.,1/(xx+yy)**(3/2)
c            = 2, compute the kernel of the 2nd term,i.e.,1/(xx+yy)**(5/2)

      subroutine tc2dkmpm(k_term)

      implicit double precision (a-h,o-z)
      parameter (max_rowxclm=16289*600, maxsubs=4*max_rowxclm)
      real s
      double precision x(2),y(2)
      complex fkernel(maxsubs)
      common /a3/ fkernel
      common /c2/ h_removed,ni,nj,dx,dy
      common /c3/ i_begin,i_end,j_begin,j_end
      common /c6/ xradius,yradius
      common /c7/ iwf

      xyz(mn,gridsize,ij)=-(1.d0*mn)*gridsize/2.d0+(ij-1.d0)*gridsize
      f11(u,v,r,rr) = -u/(v*r+rr)
      f12(u,v,w,r,rr) =  u*v*(rr+w*w)/r/(u*u*v*v+w*w*rr)
     #                   - datan(u*v/w/r)/w
      f21(u,v,r,rr) = -u*(v+2.d0*r)/(v+r)**2/r/rr
      f31(u,v,r,rr,rrr) = -u/(v+r)/rrr*
     #                    ((v/rr+4.d0/r+2.d0/(v+r))/(v+r)+2.d0/rr)
      nij(max_colm,i_row,j_colm) = (i_row-1)*max_colm + j_colm

*----------------------------------------------------------------------
c     ----------------------------
c     compute the kernel function.
c     ----------------------------

      do kk = 1,ni*nj
       fkernel(kk) = cmplx(0.,0.)
      enddo

      zz = h_removed*h_removed
      z = dsqrt(zz)
      zzz = zz*z
      half_dx = 0.5d0*dx
      half_dy = 0.5d0*dy

* properties of the kernel ellipse :
      a=xradius/dx              ! units = nodes
      b=yradius/dy              !
      ba=(b*b)/(a*a)
      e2=1.d0-ba


*----------------------------------------------------------------------
* kernel of the 1st term,i.e.,1/(xx+yy)**(3/2)

        imid=(i_end+i_begin)/2
        jmid=(j_end+j_begin)/2

        do i = i_begin,i_end
         y(1)= xyz(ni,dy,i) - half_dy
         y(2) = y(1) + dy
         y2=(i-imid)*1.d0

         do j = j_begin,j_end
          x(1)=xyz(nj,dx,j) - half_dx
          x(2) = x(1) + dx
          sum = 0.d0

          do i1=1,2
           do j1=1,2
            rr = y(i1)*y(i1) + x(j1)*x(j1) + zz
            r = dsqrt(rr)
            sum =sum + (-1)**(i1+j1)*(f11(x(j1),y(i1),r,rr)+
     .         f11(y(i1),x(j1),r,rr) - f12(x(j1),y(i1),z,r,rr))
           enddo
          enddo

          kk = nij(nj,i,j)
          s = sngl(sum)
          fkernel(kk)=cmplx(s,0.)

* integration kernel distance :
          x2=(j-jmid)*1.d0
          psi2=x2*x2+y2*y2                ! units = nodes^2

* weight by integration factor :
          if (iwf.eq.1) then
            grid=1.d0                       ! units = nodes
            psi=dsqrt(psi2)
            t1=psi+grid/2.d0
            t2=psi-grid/2.d0
            fac=psi2*psi*(t1*t1-t2*t2)/(2.d0*grid*t1*t1*t2*t2)
            fkernel(kk)=fkernel(kk)*fac
          endif

* set kernel to zero outside integration radius

          if (x2.eq.0.d0) then
            sind=1.d0
          else if (y2.eq.0.d0) then
            sind=0.d0
          else
            tanc=y2/x2
            tand=tanc*ba
            sind=1./(1.+(1./(tand*tand)))
          endif
          sin2=sind*sind
          den=dsqrt(1.d0-e2*sin2)
          rnu=a/denR
          rho=a*(1.d0-e2)/(den**3)
          r2=rho*rnu                      ! units = nodes^2

          if (psi2.gt.r2) fkernel(kk)=cmplx(0.,0.)

        enddo
       enddo

*----------------------------------------------------------------------

      kk = nij(nj,ni/2+1,nj/2+1)
      fkernel(kk)=cmplx(0.,0.)

c     ------------------------------
c     compute the fourier transform.
c     ------------------------------

      call fft3d(ni,nj,1,-1)

      return
      end


c     ************************************************************
c     *                     iftafb                               *
c     *    0:  with given fourier transform a                    *
c     *        and untransformed array b                         *
c     *    1:  computes the fourier transform of b, get fb       *
c     *    2:  multiply a with fb, get a_fb                      *
c     *    3:  compute the inverse fourier transform of afb,     *
c     *          i.e., ifafb                                     *
c     *                                                          *
c     *    programmer:   yecai li,                               *
c     *                  dept of geomatics engineering           *
c     *                  the university of calgary               *
c     *                  june 1993                               *
c     *                                                          *
c     ************************************************************

      subroutine iftafb(k_p,k_p2,k_t,mxi,mxj)

      parameter (max_rowxclm=16289*600, maxsubs=4*max_rowxclm)
      implicit double precision (a-h,o-z)
      real height(max_rowxclm)
      complex fkernel(maxsubs),fheight(maxsubs)
      common /a1/ height
      common /a3/ fheight
      common /a4/ fkernel
      common /c2/ h_removed,i0pad,j0pad,dx_m,dy_m
      common /c4/ n0bi,n0bj

      nij(max_colm,i_row,j_colm) = (i_row-1)*max_colm + j_colm

* set array fheight to zero.

      do kk=1,i0pad*j0pad
       fheight(kk) = cmplx(0.,0.)
      enddo

      h2 = 0.d0

      if(k_p.eq.0) then
        do i=1,mxi
         do j=1,mxj
          kk = nij(j0pad,i+n0bi,j+n0bj)
          fheight(kk) = cmplx(1.,0.)
         enddo
        enddo

      else
        do i=1,mxi
         do j=1,mxj
          jj = nij(mxj,i,j)
          kk = nij(j0pad,i+n0bi,j+n0bj)
          h = dble(height(jj))
          if (k_p2.ge.0) h2 = h**k_p2
          h = h**k_p
          fheight(kk) = cmplx(sngl(h),sngl(h2))
         enddo
        enddo
      endif

      call fft3d(i0pad,j0pad,1,-1)

      do kk=1,i0pad*j0pad
       fheight(kk) = fkernel(kk)*fheight(kk)
      enddo

      call fft3d(i0pad,j0pad,1,1)

      return
      end


****************************************************************************
c   imsl routine name   - fft3d
c
c-----------------------------------------------------------------------
c
c   computer            - his/single
c
c   latest revision     - june 1, 1980
c
c   purpose             - compute the fast fourier transform of
c                           a complex valued 1,2 or 3 dimensional
c                           array
c
c   usage               - call fft3d (a,ia1,ia2,n1,n2,n3,ijob,iwk,rwk,
c                           cwk)
c
c   arguments    a      - complex array. a may be a three
c                           dimensional array of dimension n1 by n2
c                           by n3, a two dimensional array of
c                           dimension n1 by n2, or a vector of
c                           length n1. on input a contains the
c                           array to be transformed. on output
c                           a is replaced by the fourier or
c                           inverse fourier transform (depending on
c                           the value of input parameter ijob).
c                ia1    - first dimension of the array a exactly
c                           as specified in the dimension statement
c                           in the calling program. (input)
c                ia2    - second dimension of the array a exactly
c                           as specified in the dimension statement
c                           in the calling program. (input)
c                n1     - limits on the first, second, and third
c                n2         subscripts of the array a, respectively.
c                n3         (input)
c                ijob   - input option parameter.
c                           if ijob is positive, the fast fourier
c                             transform of a is to be calculated.
c                           if ijob is negative, the inverse
c                             fast fourier transform of a is to be
c                             calculated.
c                iwk    - integer work vector of length
c                           6*max(n1,n2,n3)+150.
c                rwk    - real work vector of length
c                           6*max(n1,n2,n3)+150.
c                cwk    - complex work vector of length
c                           max(n2,n3).
c
c   precision/hardware  - single and double/h32
c                       - single/h32,h48,h60
c
c   reqd. imsl routines - fftcc
c
c   notation            - information on special notation and
c                           conventions is available in the manual
c                           introduction or through imsl routine uhelp
c
c   remarks  1.  if ijob is positive, fft3d calculates the fourier
c                transform, x, according to the following formula
c
c                  x(i+1,j+1,k+1)=triple sum of a(l+1,m+1,n+1)*
c                  exp(2*pi*sqrt(-1)*(i*l/n1+j*m/n2+k*n/n3))
c                  with l=0...n1-1, m=0...n2-1, n=0...n3-1
c                  and pi=3.1415...
c
c                if ijob is negative, fft3d calculates the inverse
c                fourier transform, x, according to the following
c                formula
c
c                  x(i+1,j+1,k+1)=1/(n1*n2*n3)*triple sum of
c                  a(l+1,m+1,n+1)*
c                  exp(-2*pi*sqrt(-1)*(i*l/n1+j*m/n2+k*n/n3))
c                  with l=0...n1-1, m=0...n2-1, n=0...n3-1
c                  and pi=3.1415...
c
c                note that x overwrites a on output.
c            2.  if a is a two dimensional array, set n3 = 1.
c                if a is a one dimensional array (vector),
c                set ia2 = n2 = n3 = 1.
c
c   copyright           - 1980 by imsl, inc. all rights reserved.
c
c   warranty            - imsl warrants only that imsl testing has been
c                           applied to this code. no other warranty,
c                           expressed or implied, is applicable.
c
c-----------------------------------------------------------------------

      subroutine fft3d(n1,n2,n3,ijob)

* change?
      parameter(max_rowxclm=16289*600, maxsubs=4*max_rowxclm)
      parameter(nmax=2*16289,n6max=6*nmax+150)
      integer n1,n2,n3,ijob,iwk(n6max)
      real rwk(n6max)
      complex a(maxsubs),cwk(nmax)
c                                  specifications for local variables
      integer            i,j,k,l,m,n
      real               r123
      complex            c123
      common /a3/ a
      common /b/ iwk,rwk,cwk
c                                  first executable statement
	nsub(i,j,k) = (k-1)*n1*n2 + (i-1)*n2 + j

      if (n1.le.nmax .and. n2.le.nmax .and. n3.le.nmax) goto 991
c     write(*,992)
      write(10,992)
      call flush(10)
992   format(10x,'please change the dimension limits (nmax)'/
     1  15x, 'in subroutines fft3d and fftcc')
      stop
991   continue

	n1n2n3 = n1*n2*n3
	if (ijob.gt.0) go to 10
c                                  inverse transform
	do 5 i=1,n1n2n3
	   a(i) = conjg(a(i))
    5 continue
c                                  transform third subscript
   10 do 25 l=1,n1
      do 25 m=1,n2
         do 15 n=1,n3
		cwk(n) = a(nsub(l,m,n))
   15    continue
	 call fftcc (n3)
	 do 20 k=1,n3
	    a(nsub(l,m,k)) = cwk(k)
   20    continue
   25 continue
c                                  transform second subscript
      do 40 l=1,n1
      do 40 k=1,n3
         do 30 m=1,n2
		cwk(m) = a(nsub(l,m,k))
   30    continue
	 call fftcc (n2)
	 do 35 j=1,n2
	    a(nsub(l,j,k)) = cwk(j)
   35    continue
   40 continue
c                                  transform first subscript
      do 45 j=1,n2
      do 45 k=1,n3

	 do 130 m=1,n1
	    cwk(m) = a(nsub(m,j,k))
  130    continue
	 call fftcc (n1)
	 do 135 m=1,n1
	    a(nsub(m,j,k)) = cwk(m)
  135    continue

c         call fftcc (a(1,j,k),n1,iwk,rwk)
   45 continue
      if (ijob.gt.0) go to 55
	r123 = n1n2n3
	c123 = cmplx(r123,0.0)
	do 50 i=1,n1
	do 50 j=1,n2
	do 50 k=1,n3
	   a(nsub(i,j,k)) = conjg(a(nsub(i,j,k)))/c123
   50 continue
   55 return
      end

****************************************************************************
c   imsl routine name   - fftcc
c
c-----------------------------------------------------------------------
c
c   computer            - his/single
c
c   latest revision     - january 1, 1978
c
c   purpose             - compute the fast fourier transform of a
c                           complex valued sequence
c
c   usage               - call fftcc (a,n,iwk,wk)
c
c   arguments    a      - complex vector of length n. on input a
c                           contains the complex valued sequence to be
c                           transformed. on output a is replaced by the
c                           fourier transform.
c                n      - input number of data points to be
c                           transformed. n may be any positive
c                           integer.
c                iwk    - integer work vector of length 6*n+150.
c                           (see programming notes for further details)
c                wk     - real work vector of length 6*n+150.
c                           (see programming notes for further details)
c
c   precision/hardware  - single and double/h32
c                       - single/h36,h48,h60
c
c   reqd. imsl routines - none required
c
c   notation            - information on special notation and
c                           conventions is available in the manual
c                           introduction or through imsl routine uhelp
c
c   remarks  1.  fftcc computes the fourier transform, x, according
c                to the following formula;
c
c                  x(k+1) = sum from j = 0 to n-1 of
c                           a(j+1)*cexp((0.0,(2.0*pi*j*k)/n))
c                  for k=0,1,...,n-1 and pi=3.1415...
c
c                note that x overwrites a on output.
c            2.  fftcc can be used to compute
c
c                  x(k+1) = (1/n)*sum from j = 0 to n-1 of
c                           a(j+1)*cexp((0.0,(-2.0*pi*j*k)/n))
c                  for k=0,1,...,n-1 and pi=3.1415...
c
c                by performing the following steps;
c
c                     do 10 i=1,n
c                        a(i) = conjg(a(i))
c                  10 continue
c                     call fftcc (a,n,iwk,wk)
c                     do 20 i=1,n
c                        a(i) = conjg(a(i))/n
c                  20 continue
c
c   copyright           - 1978 by imsl, inc. all rights reserved.
c
c   warranty            - imsl warrants only that imsl testing has been
c                           applied to this code. no other warranty,
c                           expressed or implied, is applicable.
c
c-----------------------------------------------------------------------
c
      subroutine fftcc (n)
c                                  specifications for arguments
* change?
	parameter(nmax=2*16289,n6max=6*nmax+150)
      integer            n,iwk(n6max)
      real               wk(n6max)
      complex            a(nmax)
c                                  specifications for local variables
      integer            i,iam,iap,ibm,ibp,ic,icc,icf,ick,id,idm1,ii,
     1                   ija,ikb,ikt,ill,im,ird,isf,isk,isp,iss,ita,itb,
     2                   j,ja,jf,jj,jk,k,k0,k1,k2,k3,ka,kb,kd2,kf,kh,kn,
     3                   kt,ktp,l,l1,m,mm,mm1,mp
      real               cm,sm,c1,c2,c3,s1,s2,s3,c30,rad,a0,a1,a4,b4,
     1                   a2,a3,b0,b1,b2,b3,zero,half,one,two,z0(2),
     2                   z1(2),z2(2),z3(2),z4(2)
      complex            za0,za1,za2,za3,za4,ak2
      equivalence        (za0,z0(1)),(za1,z1(1)),(za2,z2(1)),
     1                   (za3,z3(1)),(a0,z0(1)),(b0,z0(2)),(a1,z1(1)),
     2                   (b1,z1(2)),(a2,z2(1)),(b2,z2(2)),(a3,z3(1)),
     3                   (b3,z3(2)),(za4,z4(1)),(z4(1),a4),(z4(2),b4)
      common /b/ iwk,wk,a
      data               rad/6.28318531/,
     1                   c30/.866025404/
      data               zero,half,one,two/0.0,0.5,1.0,2.0/
c                                  first executable statement
      if (n .eq. 1) go to 9005
      k = n
      m = 0
      j = 2
      jj = 4
      jf = 0
c                                  determine the square factors of n
      iwk(1) = 1
    5 i = k/jj
      if (i*jj .ne. k) go to 10
      m = m+1
      iwk(m+1) = j
      k = i
      go to 5
   10 j = j + 2
      if (j .eq. 4) j = 3
      jj = j * j
      if (jj .le. k) go to 5
      kt = m
c                                  determine the remaining factors of n
      j = 2
   15 i = k / j
      if (i*j .ne. k) go to 20
      m = m + 1
      iwk(m+1) = j
      k = i
      go to 15
   20 j = j + 1
      if (j .eq. 3) go to 15
      j = j + 1
      if(j.le.k) go to 15
      k = iwk(m+1)
      if (iwk(kt+1) .gt. iwk(m+1)) k = iwk(kt+1)
      if(kt.le.0) go to 30
      ktp = kt + 2
      do 25  i = 1,kt
         j = ktp - i
         m = m+1
         iwk(m+1) = iwk(j)
   25 continue
   30 mp = m+1
      ic = mp+1
      id = ic+mp
      ill = id+mp
      ird = ill+mp+1
      icc = ird+mp
      iss = icc+mp
      ick = iss+mp
      isk = ick+k
      icf = isk+k
      isf = icf+k
      iap = isf+k
      kd2 = (k-1) / 2 + 1
      ibp = iap + kd2
      iam = ibp + kd2
      ibm = iam + kd2
      mm1 = m-1
      i = 1
   35 l = mp - i
      j = ic - i
      iwk(ill+l) = 0
      if ((iwk(j-1) + iwk(j)) .eq. 4) iwk(ill+l) = 1
      if (iwk(ill+l) .eq. 0) go to 40
      i = i + 1
      l = l - 1
      iwk(ill+l) = 0
   40 i = i + 1
      if(i.le.mm1) go to 35
      iwk(ill+1) = 0
      iwk(ill+mp) = 0
      iwk(ic) = 1
      iwk(id) = n
      do 45  j = 1,m
         k = iwk(j+1)
         iwk(ic+j) = iwk(ic+j-1) * k
         iwk(id+j) = iwk(id+j-1) / k
         wk(ird+j) = rad/iwk(ic+j)
         c1 = rad/k
         if (k .le. 2) go to 45
         wk(icc+j) = cos(c1)
         wk(iss+j) = sin(c1)
   45 continue
      mm = m
      if (iwk(ill+m) .eq. 1) mm = m - 1
      if (mm .le. 1) go to 50
      sm = iwk(ic+mm-2) * wk(ird+m)
      cm = cos(sm)
      sm = sin(sm)
   50 kb = 0
      kn = n
      jj = 0
      i = 1
      c1 = one
      s1 = zero
      l1 = 1
   55 if (iwk(ill+i+1) .eq. 1) go to 60
      kf = iwk(i+1)
      go to 65
   60 kf = 4
      i = i+1
   65 isp = iwk(id+i)
      if (l1 .eq. 1) go to 70
      s1 = jj * wk(ird+i)
      c1 = cos(s1)
      s1 = sin(s1)
c                                  factors of 2, 3, and 4 are
c                                  handled separately.
   70 if (kf .gt. 4) go to 140
      go to (75,75,90,115), kf
   75 k0 = kb + isp
      k2 = k0 + isp
      if (l1 .eq. 1) go to 85
   80 k0 = k0 - 1
      if (k0 .lt. kb) go to 190
      k2 = k2 - 1
      za4 = a(k2+1)
      a0 = a4*c1-b4*s1
      b0 = a4*s1+b4*c1
      a(k2+1) = a(k0+1)-za0
      a(k0+1) = a(k0+1)+za0
      go to 80
   85 k0 = k0 - 1
      if (k0 .lt. kb) go to 190
      k2 = k2 - 1
      ak2 = a(k2+1)
      a(k2+1) = a(k0+1)-ak2
      a(k0+1) = a(k0+1)+ak2
      go to 85
   90 if (l1 .eq. 1) go to 95
      c2 = c1 * c1 - s1 * s1
      s2 = two * c1 * s1
   95 ja = kb + isp - 1
      ka = ja + kb
      ikb = kb+1
      ija = ja+1
      do 110 ii = ikb,ija
         k0 = ka - ii + 1
         k1 = k0 + isp
         k2 = k1 + isp
         za0 = a(k0+1)
         if (l1 .eq. 1) go to 100
         za4 = a(k1+1)
         a1 = a4*c1-b4*s1
         b1 = a4*s1+b4*c1
         za4 = a(k2+1)
         a2 = a4*c2-b4*s2
         b2 = a4*s2+b4*c2
         go to 105
  100    za1 = a(k1+1)
         za2 = a(k2+1)
  105    a(k0+1) = cmplx(a0+a1+a2,b0+b1+b2)
         a0 = -half * (a1+a2) + a0
         a1 = (a1-a2) * c30
         b0 = -half * (b1+b2) + b0
         b1 = (b1-b2) * c30
         a(k1+1) = cmplx(a0-b1,b0+a1)
         a(k2+1) = cmplx(a0+b1,b0-a1)
  110 continue
      go to 190
  115 if (l1 .eq. 1) go to 120
      c2 = c1 * c1 - s1 * s1
      s2 = two * c1 * s1
      c3 = c1 * c2 - s1 * s2
      s3 = s1 * c2 + c1 * s2
  120 ja = kb + isp - 1
      ka = ja + kb
      ikb = kb+1
      ija = ja+1
      do 135 ii = ikb,ija
         k0 = ka - ii + 1
         k1 = k0 + isp
         k2 = k1 + isp
         k3 = k2 + isp
         za0 = a(k0+1)
         if (l1 .eq. 1) go to 125
         za4 = a(k1+1)
         a1 = a4*c1-b4*s1
         b1 = a4*s1+b4*c1
         za4 = a(k2+1)
         a2 = a4*c2-b4*s2
         b2 = a4*s2+b4*c2
         za4 = a(k3+1)
         a3 = a4*c3-b4*s3
         b3 = a4*s3+b4*c3
         go to 130
  125    za1 = a(k1+1)
         za2 = a(k2+1)
         za3 = a(k3+1)
  130    a(k0+1) = cmplx(a0+a2+a1+a3,b0+b2+b1+b3)
         a(k1+1) = cmplx(a0+a2-a1-a3,b0+b2-b1-b3)
         a(k2+1) = cmplx(a0-a2-b1+b3,b0-b2+a1-a3)
         a(k3+1) = cmplx(a0-a2+b1-b3,b0-b2-a1+a3)
  135 continue
      go to 190
  140 jk = kf - 1
      kh = jk/2
      k3 = iwk(id+i-1)
      k0 = kb + isp
      if (l1 .eq. 1) go to 150
      k = jk - 1
      wk(icf+1) = c1
      wk(isf+1) = s1
      do 145 j = 1,k
         wk(icf+j+1) = wk(icf+j) * c1 - wk(isf+j) * s1
         wk(isf+j+1) = wk(icf+j) * s1 + wk(isf+j) * c1
  145 continue
  150 if (kf .eq. jf) go to 160
      c2 = wk(icc+i)
      wk(ick+1) = c2
      wk(ick+jk) = c2
      s2 = wk(iss+i)
      wk(isk+1) = s2
      wk(isk+jk) = -s2
      do 155 j = 1,kh
         k = jk - j
         wk(ick+k) = wk(ick+j) * c2 - wk(isk+j) * s2
         wk(ick+j+1) = wk(ick+k)
         wk(isk+j+1) = wk(ick+j) * s2 + wk(isk+j) * c2
         wk(isk+k) = -wk(isk+j+1)
  155 continue
  160 k0 = k0 - 1
      k1 = k0
      k2 = k0 + k3
      za0 = a(k0+1)
      a3 = a0
      b3 = b0
      do 175 j = 1,kh
         k1 = k1 + isp
         k2 = k2 - isp
         if (l1 .eq. 1) go to 165
         k = kf - j
         za4 = a(k1+1)
         a1 = a4*wk(icf+j)-b4*wk(isf+j)
         b1 = a4*wk(isf+j)+b4*wk(icf+j)
         za4 = a(k2+1)
         a2 = a4*wk(icf+k)-b4*wk(isf+k)
         b2 = a4*wk(isf+k)+b4*wk(icf+k)
         go to 170
  165    za1 = a(k1+1)
         za2 = a(k2+1)
  170    wk(iap+j) = a1 + a2
         wk(iam+j) = a1 - a2
         wk(ibp+j) = b1 + b2
         wk(ibm+j) = b1 - b2
         a3 = a1 + a2 + a3
         b3 = b1 + b2 + b3
  175 continue
      a(k0+1) = cmplx(a3,b3)
      k1 = k0
      k2 = k0 + k3
      do 185 j = 1,kh
         k1 = k1 + isp
         k2 = k2 - isp
         jk = j
         a1 = a0
         b1 = b0
         a2 = zero
         b2 = zero
         do 180  k = 1,kh
            a1 = a1 + wk(iap+k) * wk(ick+jk)
            a2 = a2 + wk(iam+k) * wk(isk+jk)
            b1 = b1 + wk(ibp+k) * wk(ick+jk)
            b2 = b2 + wk(ibm+k) * wk(isk+jk)
            jk = jk + j
            if (jk .ge. kf) jk = jk - kf
  180    continue
         a(k1+1) = cmplx(a1-b2,b1+a2)
         a(k2+1) = cmplx(a1+b2,b1-a2)
  185 continue
      if (k0 .gt. kb) go to 160
      jf = kf
  190 if ( i .ge. mm ) go to 195
      i = i + 1
      go to 55
  195 i = mm
      l1 = 0
      kb = iwk(id+i-1) + kb
      if (kb .ge. kn) go to 215
  200 jj = iwk(ic+i-2) + jj
      if (jj .lt. iwk(ic+i-1)) go to 205
      i = i - 1
      jj = jj - iwk(ic+i)
      go to 200
  205 if (i .ne. mm) go to 210
      c2 = c1
      c1 = cm * c1 - sm * s1
      s1 = sm * c2 + cm * s1
      go to 70
  210 if (iwk(ill+i) .eq. 1) i = i + 1
      go to 55
  215 i = 1
      ja = kt - 1
      ka = ja + 1
      if(ja.lt.1) go to 225
      do 220  ii = 1,ja
         j = ka - ii
         iwk(j+1) = iwk(j+1) - 1
         i = iwk(j+1) + i
  220 continue
c                                  the result is now permuted to
c                                  normal order.
  225 if (kt .le. 0) go to 270
      j = 1
      i = 0
      kb = 0
  230 k2 = iwk(id+j) + kb
      k3 = k2
      jj = iwk(ic+j-1)
      jk = jj
      k0 = kb + jj
      isp = iwk(ic+j) - jj
  235 k = k0 + jj
  240 za4 = a(k0+1)
      a(k0+1) = a(k2+1)
      a(k2+1) = za4
      k0 = k0 + 1
      k2 = k2 + 1
      if (k0 .lt. k) go to 240
      k0 = k0 + isp
      k2 = k2 + isp
      if (k0 .lt. k3) go to 235
      if (k0 .ge. k3 + isp) go to 245
      k0 = k0 - iwk(id+j) + jj
      go to 235
  245 k3 = iwk(id+j) + k3
      if (k3 - kb .ge. iwk(id+j-1)) go to 250
      k2 = k3 + jk
      jk = jk + jj
      k0 = k3 - iwk(id+j) + jk
      go to 235
  250 if (j .ge. kt) go to 260
      k = iwk(j+1) + i
      j = j + 1
  255 i = i + 1
      iwk(ill+i) = j
      if (i .lt. k) go to 255
      go to 230
  260 kb = k3
      if (i .le. 0) go to 265
      j = iwk(ill+i)
      i = i - 1
      go to 230
  265 if (kb .ge. n) go to 270
      j = 1
      go to 230
  270 jk = iwk(ic+kt)
      isp = iwk(id+kt)
      m = m - kt
      kb = isp/jk-2
      if (kt .ge. m-1 ) go to 9005
      ita = ill+kb+1
      itb = ita+jk
      idm1 = id-1
      ikt = kt+1
      im = m+1
      do 275 j = ikt,im
         iwk(idm1+j) = iwk(idm1+j)/jk
  275 continue
      jj = 0
      do 290 j = 1,kb
         k = kt
  280    jj = iwk(id+k+1) + jj
         if (jj .lt. iwk(id+k)) go to 285
         jj = jj - iwk(id+k)
         k = k + 1
         go to 280
  285    iwk(ill+j) = jj
         if (jj .eq. j) iwk(ill+j) = -j
  290 continue
c                                  determine the permutation cycles
c                                  of length greater than or equal
c                                  to two.
      do 300  j = 1,kb
         if (iwk(ill+j) .le. 0) go to 300
         k2 = j
  295    k2 = iabs(iwk(ill+k2))
         if (k2 .eq. j) go to 300
         iwk(ill+k2) = -iwk(ill+k2)
         go to 295
  300 continue
c                                  reorder a following the
c                                  permutation cycles
      i = 0
      j = 0
      kb = 0
      kn = n
  305 j = j + 1
      if (iwk(ill+j) .lt. 0) go to 305
      k = iwk(ill+j)
      k0 = jk * k + kb
  310 za4 = a(k0+i+1)
      wk(ita+i) = a4
      wk(itb+i) = b4
      i = i + 1
      if (i .lt. jk) go to 310
      i = 0
  315 k = -iwk(ill+k)
      jj = k0
      k0 = jk * k + kb
  320 a(jj+i+1) = a(k0+i+1)
      i = i + 1
      if (i .lt. jk) go to 320
      i = 0
      if (k .ne. j) go to 315
  325 a(k0+i+1) = cmplx(wk(ita+i),wk(itb+i))
      i = i + 1
      if (i .lt. jk) go to 325
      i = 0
      if (j .lt. k2) go to 305
      j = 0
      kb = kb + isp
      if (kb .lt. kn) go to 305
 9005 return
      end
