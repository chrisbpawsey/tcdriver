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

! set array fheight to zero.

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


