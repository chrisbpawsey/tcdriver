subroutine fft3d(n1,n2,n3,ijob)

! change?
      parameter(max_rowxclm=16289*600, maxsubs=4*max_rowxclm)
      parameter(nmax=2*16289,n6max=6*nmax+150)
      integer n1,n2,n3,ijob,iwk(n6max)
      real rwk(n6max)
      complex a(maxsubs),cwk(nmax)
!                                  specifications for local variables
      integer            i,j,k,l,m,n
      real               r123
      complex            c123
      common /a3/ a
      common /b/ iwk,rwk,cwk
!                                  first executable statement
	nsub(i,j,k) = (k-1)*n1*n2 + (i-1)*n2 + j

      if (n1.le.nmax .and. n2.le.nmax .and. n3.le.nmax) goto 991
!     write(*,992)
      write(10,992)
      call flush(10)
992   format(10x,'please change the dimension limits (nmax)'  &
           & 15x, 'in subroutines fft3d and fftcc')
      stop
991   continue

	n1n2n3 = n1*n2*n3
	if (ijob.gt.0) go to 10
!                                  inverse transform
	do 5 i=1,n1n2n3
	   a(i) = conjg(a(i))
    5 continue
!                                  transform third subscript
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
!                                  transform second subscript
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
!                                  transform first subscript
      do 45 j=1,n2
      do 45 k=1,n3

	 do 130 m=1,n1
	    cwk(m) = a(nsub(m,j,k))
  130    continue
	 call fftcc (n1)
	 do 135 m=1,n1
	    a(nsub(m,j,k)) = cwk(m)
  135    continue

!         call fftcc (a(1,j,k),n1,iwk,rwk)
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

