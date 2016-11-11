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

!     the unit of g*density is 1/sec**2
!     when the unit of height is meter and g*density*1.e5
!     the unit of tz is mgal
!     c_mgal = g*density*1.d5*i0pad*j0pad
      c_mgal = g*density*1.d8*i0pad*j0pad
      alpha_sq = h_removed * h_removed


!      compute the 1st- and/or the 2nd-order term
!         ---------------------------------------------------
!         |   compute the kernel function 1./(xx+yy)**1.5   |
!         |                  and fft                        |
!         ---------------------------------------------------

! mass-prism :
      call tc2dkmpm(1)

!     ---------------------------------------
!     correct for the origin shift and
!     move f{kernel} from fheight to fkernel.
!     ---------------------------------------

      do i = 1,i0pad
       do j = 1,j0pad
        kk = nij(j0pad,i,j)
        fkernel(kk) = fheight(kk)*(-1)**(i+j)
       enddo
      enddo

!     --------------------------------------------------
!     compute the convolution h**2 with 1/(xx+yy)**(3/2)
!     compute the convolution h with 1/(xx+yy)**(3/2)
!     --------------------------------------------------

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

!     -----------------------------------------------
!     compute the convolution 1 with 1/(xx+yy)**(3/2)
!     -----------------------------------------------

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

