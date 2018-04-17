#
# KD_TREE2
#      A C++ set of classes for finding close neighbors in data sets of
#      points in k-dimensional Euclidean space. 
#

#
# For Intel Fortran 8.0.    
#
CPP=g++ -Wall


# Uncomment next line for debugging
#FLAGS=-g3

# 
# for your specific CPU
#
#FLAGSMACHINE=-march=athlon64 -mfpmath=sse -msse2 -m64

FLAGS=-O3 -ffast-math ${FLAGSMACHINE} -fomit-frame-pointer


MY_DIR=`basename ${PWD}`

all:	kdtree2_test

kdtree2_test:	kdtree2_test.o kdtree2.o
	${CPP} ${FLAGS} -o kdtree2_test kdtree2_test.o kdtree2.o 

#kdtree2_test.o :: kdtree2.o

%.o :: %.cpp %.hpp
	${CPP} ${FLAGS} -c $<

%.o :: %.cpp
	${CPP} ${FLAGS} -c $<

tar:
	/bin/rm *~; cd ..; tar zcvf ${MY_DIR}.tgz ${MY_DIR}

clean:
	/bin/rm -f *.o *.mod *.il kdtree_test kdtree2_test core.* test
