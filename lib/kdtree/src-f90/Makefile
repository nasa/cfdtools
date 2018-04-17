#
# KD_TREE2
#      A Fortran 95 module for finding close neighbors in data sets of
#      points in k-dimensional Euclidean space. 
#



# 
# For Intel Fortran 8.1
#
#FLAGSMACHINE= -xK -arch pn4 -rounding_mode=chopped  # Specific to your machine architecture
#FLAGSALWAYS = -static-libcxa -u -warn all  # for both optimization and debugging
#FLAGSOPT=-O3 -ipo -fpe0 ${FLAGSMACHINE} -fno-alias 
FLAGSALWAYS = -static -u -warn all  # for both optimization and debugging
FLAGSOPT=-O3 -r8 ${FLAGSMACHINE}
FLAGSDEBUG= -check all -traceback  #-g #-traceback # -C
F90=ifort


#
# For g95
#
#FLAGSMACHINE=-march=athlon-xp -m3dnow -msse -mno-ieee-fp -mfpmath=sse # Athlon XP for example
#FLAGSALWAYS = -u -Wall
#FLAGSOPT= -O3 -fomit-frame-pointers
#FLAGSDEBUG= -g
#F90=g95
#
#

# for gfortran
#
#FLAGSMACHINE=#-march=athlon-xp -m3dnow -msse -mno-ieee-fp -mfpmath=sse # Athlon XP for example
#FLAGSALWAYS = -u -Wall
#FLAGSOPT= -O3 -fomit-frame-pointers
#FLAGSDEBUG= -g
#F90=gfortran


# 
# choose debugging or optimization
FLAGS= ${FLAGSMACHINE} ${FLAGSALWAYS} ${FLAGSOPT}  #  change the last to ${FLAGSOPT} or ${FLAGSDEBUG}


MY_DIR=`basename ${PWD}`

all:	kdtree2_test kdtree_test #r_test

r_test:	r_test.o kdtree2.o
	${F90} ${FLAGS} -o r_test r_test.o kdtree2.o 

kdtree2_test:	kdtree2_test.o kdtree2.o
	${F90} ${FLAGS} -o kdtree2_test kdtree2_test.o kdtree2.o 

kdtree2_test.o: kdtree2.o

kdtree_test:	kdtree_test.o kdtree.o
	${F90} ${FLAGS} -o kdtree_test kdtree_test.o kdtree.o 

kdtree_test.o: kdtree.o

%.o :: %.f90
	${F90} ${FLAGS} -c $<

tar:
	cd ..; tar zcvf ${MY_DIR}.tgz ${MY_DIR}

clean:
	/bin/rm -f *.o *.mod *.il kdtree_test kdtree2_test
