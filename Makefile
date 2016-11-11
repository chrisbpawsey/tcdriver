####### gnu- gfortran comiler and compiler flag options.
FC      = gfortran 
FFLAGS = -O1 # -malign-double # -Minfo=ftn -Minform=inform # 
DBGFLAGS = -C -g # -std=f2003 # -pedantic -std=f2003

##### intel compiler and compiler flag options
#FC 	= ifort  # intel compiler
#FFLAGS = -f77rtl    -m64
# bounds checking -WB -warn all -W1
#DBGFLAGS = -debug extended,emit-colum  -ftrapuv -g3  

##### pgi compiler and compiler flag options
#FC      = pgfortran 			# pgi compiler
#DBGFLAGS =  -C  -g -Minfo=ftn -Minform=inform -Mchkstk -Melf # -help=debug|language
#FFLAGS = -Mstandard

#  object list
OBJS 	= tc_variables.o tc2dftpl.o tc2dkmpm.o dtdz_2d2.o iftafb.o fft3d.o fftcc.o tcdriver.o 

run-tc:		${OBJS} 
	${FC} ${DBGFLAGS} ${FFLAGS} -o  $@ ${OBJS}

tc_variables.o:	tc_variables.f90
	${FC} ${DBGFLAGS} ${FFLAGS} -c $?

tcdriver.o:	tcdriver.f90
	${FC} ${DBGFLAGS} ${FFLAGS} -c $?

tc2dftpl.o: 	tc2dftpl.f90
	${FC} ${DBGFLAGS} ${FFLAGS} -c $?

tc2dkmpm.o: 	tc2dkmpm.f90
	${FC} ${DBGFLAGS} ${FFLAGS} -c $?

dtdz_2d2.o: 	dtdz_2d2.f90
	${FC} ${DBGFLAGS} ${FFLAGS} -c $?

fft3d.o: 	fft3d.f90
	${FC} ${DBGFLAGS} ${FFLAGS} -c $?

fftcc.o:	fftcc.f90
	${FC} ${DBGFLAGS} ${FFLAGS} -c $?

clean:
	-rm -f *.o core run-tc a.out *.mod

