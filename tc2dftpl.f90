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

! position identifier:
      nij(max_colm,i_row,j_colm) = (i_row-1)*max_colm + j_colm

      g = 6.672d-8
      density = 2.67d0

      dx_m = dx_km
      dy_m = dy_km

!***********************************************************************
! compute the statistics of the heights
! convert heights from metres to km

      hmean=0.d0
      hstd=0.d0

      do jj=1,maxi*maxj
       height(jj) = height(jj)*1.e-3
       hmean=hmean+height(jj)
       hstd=hstd+height(jj)*height(jj)
      enddo

      hmean=hmean/(1.d0*maxi*maxj)
      hstd=dsqrt(hstd/(1.d0*maxi*maxj)-hmean*hmean)

!********************************************************************
!  take the parameter alpha in the kernel function as
!  sigma(h)/sqrt(2)

      h_removed = hstd/dsqrt(2.d0)

!*********************************************************************
!  do the computation from ...
!  the zero-order term (k_01 = 0);
!  or 1st-order term (k_01 = 1):

      k_01 = 1

!********************************************************************
! specify how many terms to be computed (in the denominator) :
!  k_term = 1: (x*x+y*y)**(3/2) (consider the first term only),
!  k_term = 2: (x*x+y*y)**(5/2) (consider the 1st & 2nd terms).
!  k_term = 3: (x*x+y*y)**(7/2) (the 1st, 2nd & 3rd terms).

      k_term = 1

!********************************************************************
! compute the kernel function ...
!  over whole computation area with 100% zero padding (ikernel=1);
!  within a limited area with either 100% or 50% padding (ikernel=2)

      ikernel=2

! compute over whole area [default] :
      int_ker = 1

! set integral radius = the length of the whole area
      xradius = maxj*dx_km
      yradius = maxi*dy_km

      iradius = maxi
      jradius = maxj

! number of zeros padded to each side :
      n0bi = maxi/2
      n0bj = maxj/2

! new array dimensions after 100% zeros are padded (ie, twice the size)
      i0pad = maxi*2
      j0pad = maxj*2

      i_begin = 1
      i_end = i0pad
      j_begin = 1
      j_end = j0pad

!--------------------------------------------------------------------
! compute over limited area :

      if (ikernel.eq.2) then

        int_ker = 2

! set integration radius to half the strip dimensions :
!       jradius=maxj/2
!       iradius=maxi/2
        jradius=intrad
        iradius=intrad
        xradius = jradius*dx_km
        yradius = iradius*dy_km

! option to zero pad this kernel :
!  ipad = 1: 100% of the number of the heights,
!       = 2: 50% of (the number of height + the size of kernel)

        ipad=1

        if (ipad.eq.2) then
          n0bi = iradius + 1
          n0bj = jradius + 1
          i0pad = maxi + n0bi*2
          j0pad = maxj + n0bj*2
        endif

! compute limits of the integral of kernel function.
        i_begin = max0(1,i0pad/2+1-iradius)
        i_end   = min0(i0pad,i0pad/2+1+iradius)
        j_begin = max0(1,j0pad/2+1-jradius)
        j_end   = min0(j0pad,j0pad/2+1+jradius)

      endif

!********************************************************************
!  topographic model may be ...
!  mass-prism (k_tm=1) or mass-line (k_tm=2) :

      k_tm = 1

!********************************************************************
! include integration weighting factor (iwf=1); don't (iwf=2)
! in subroutine tc2dkmpm

      iwf=1

!********************************************************************
! compute the convolutions two simultaneously

      call dtdz_2d2(k_tm,k_01)

!********************************************************************
      write(*,209)
      write(10,209)
      call flush(10)
209   format('terrain computations have been done successfully'/)


      return
      end

