
!     ***************************************************
!     *                 tc2dkmpm                        *
!     *    terrain correction 2 dimensional kernel      *
!     *     with mass-prism topographic model           *
!     *    s  : simple summation formulas are used      *
!     *    2dk: the 2 ddimensional kernel function      *
!     *                                                 *
!     *    programmer:   yecai li,                      *
!     *                  dept of geomatics engineering  *
!     *                  the university of calgary      *
!     *                  june 1993                      *
!     *                                                 *
!     *    modified by j f kirby to include weighting   *
!     *    function and zeros outside cap               *
!     *                                                 *
!     ***************************************************
!     k_term = 1, compute the kernel of the 1st term,i.e.,1/(xx+yy)**(3/2)
!            = 2, compute the kernel of the 2nd term,i.e.,1/(xx+yy)**(5/2)

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
      f12(u,v,w,r,rr) =  u*v*(rr+w*w)/r/(u*u*v*v+w*w*rr)  &
           & - datan(u*v/w/r)/w
      f21(u,v,r,rr) = -u*(v+2.d0*r)/(v+r)**2/r/rr
      f31(u,v,r,rr,rrr) = -u/(v+r)/rrr*    &
           & ((v/rr+4.d0/r+2.d0/(v+r))/(v+r)+2.d0/rr)   
      nij(max_colm,i_row,j_colm) = (i_row-1)*max_colm + j_colm

!----------------------------------------------------------------------
!     ----------------------------
!     compute the kernel function.
!     ----------------------------

      do kk = 1,ni*nj
       fkernel(kk) = cmplx(0.,0.)
      enddo

      zz = h_removed*h_removed
      z = dsqrt(zz)
      zzz = zz*z
      half_dx = 0.5d0*dx
      half_dy = 0.5d0*dy

! properties of the kernel ellipse :
      a=xradius/dx              ! units = nodes
      b=yradius/dy              !
      ba=(b*b)/(a*a)
      e2=1.d0-ba


!----------------------------------------------------------------------
! kernel of the 1st term,i.e.,1/(xx+yy)**(3/2)

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
            sum = sum + (-1)**(i1+j1)*(f11(x(j1),y(i1),r,rr)+   &
     &         f11(y(i1),x(j1),r,rr) - f12(x(j1),y(i1),z,r,rr))
           enddo
          enddo

          kk = nij(nj,i,j)
          s = sngl(sum)
          fkernel(kk)=cmplx(s,0.)

! integration kernel distance :
          x2=(j-jmid)*1.d0
          psi2=x2*x2+y2*y2                ! units = nodes^2

! weight by integration factor :
          if (iwf.eq.1) then
            grid=1.d0                       ! units = nodes
            psi=dsqrt(psi2)
            t1=psi+grid/2.d0
            t2=psi-grid/2.d0
            fac=psi2*psi*(t1*t1-t2*t2)/(2.d0*grid*t1*t1*t2*t2)
            fkernel(kk)=fkernel(kk)*fac
          endif

! set kernel to zero outside integration radius

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

!----------------------------------------------------------------------

      kk = nij(nj,ni/2+1,nj/2+1)
      fkernel(kk)=cmplx(0.,0.)

!     ------------------------------
!     compute the fourier transform.
!     ------------------------------

      call fft3d(ni,nj,1,-1)

      return
      end


