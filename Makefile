FFLAGS = -O1 # -malign-double # -Minfo=ftn -Minform=inform # 
#DBGFLAGS = -C -g # -pedantic -std=f2003
# intel flags
#FFLAGS = -f77rtl    -m64
# bounds checking -WB -warn all -W1
#DBGFLAGS = -debug extended,emit-colum  -ftrapuv -g3  
# pgi flags
#  -Minfo=ftn -Minform=inform -help=debug|language
FC 	= ifort  # intel compiler
#FC      = pgfortran # pgi compiler
#FC      = gfortran 
OBJS 	= tcdriver-for9s.o 

run-tc:		${OBJS} 
	${FC} ${DBGFLAGS} ${FFLAGS} -o  $@ ${OBJS}

tcdriver-for9s.o:	tcdriver-for9s.f
	${FC} ${DBGFLAGS} ${FFLAGS} -c $?

clean:
	-rm -f *.o core run-tc a.out

