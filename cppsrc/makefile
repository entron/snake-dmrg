CC =	g++      	

MYUSR= $(HOME)/usr
LAPACKPP=$(MYUSR)/include/lapackpp
MYLIB=$(MYUSR)/lib
CFLAGS = -I$(LAPACKPP) -O3 -mtune=core2  
BLASLIB= -lcblas -lf77blas -latlas 
LAPACKLIB=-llapack
ARPACKLIB=-larpack

SRC=	blocham.cpp\
	block.cpp\
	dtmat.cpp\
	public.cpp\
	site.cpp\
	supblock.cpp\
	gqnbase.cpp\
	gqnmat.cpp\
	gqn.cpp\
	dmrg.cpp\
	supblock_iDMRG.cpp\
	supblock_tDMRG.cpp\
	snake.cpp

OBJS=	blocham.o\
	block.o\
	dtmat.o\
	public.o\
	site.o\
	supblock.o\
	gqnbase.o\
	gqnmat.o\
	gqn.o\
	dmrg.o\
	supblock_iDMRG.o\
	supblock_tDMRG.o\
	snake.o


LIBS = -Wl,--rpath -Wl,$(MYLIB) -L$(MYLIB) $(ARPACKLIB) -llapackpp $(LAPACKLIB) $(BLASLIB)

.cpp.o:
	$(CC) -c $(CFLAGS) $<

snake: $(OBJS)
	$(CC) $(CFLAGS) -o $@ $(OBJS) $(LIBS)


clean:
	rm -f *.o core

