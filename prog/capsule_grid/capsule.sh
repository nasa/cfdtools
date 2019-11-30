#! /bin/bash

# This is a slimmed down version of Todd White's script, without the 3D gridding or CBAERO mesh steps.
# Todd's documentation follows:

echo "This script will use David Saunders' CAPSULE_GRID to create:"
echo " * An axisymmetric generatrix (gen.dat)"
echo "Next, it uses William Chan's HYPGEN and GRIDED to create:"
echo " * A Hyperbolic Axisymmetric single-block plot3d grid (dplr2d_grid.gu)"
echo "It then uses radial_interp to create:"
echo " * A Hyperbolic 3D multi-block plot3d grid (dplr3d_grid.gu)"
echo "Finally, a decimated plot3D surface mesh is converted into a CBAERO mesh file"
echo "using Mike Aftmosis's triangulate, trix, followed by Dave Kinney's cart2mesh:"
echo " * The surface mesh is cbaero.msh, requires duplicate points to be removed"

sleep 1

#if [[ $# -lt 1 ]]; then
#  echo "$0 CAPSULE_GRID inputfile"
#  exit -1
#fi

#module use -a /home/trwhite1/software/modulefiles
#module load cfdtools
#module load cbaero/latest
#module load CART3D/v1.5.1_X86_64_ICC--16.06.03

# Run capsule_grid using capsule.inp
echo "** CAPSULE_GRID"
~dasaund1/CAPSULE_GRID/capsule_grid < capsule.inp > capsule.log

# Check how many points are in the output
ni=`head -n1 gen.dat | sed 's|#||g'`

# Create a new generatrix (different format and order) to use with hypgen
echo "1" > capsule.srf
echo "$ni 1 1" >> capsule.srf
awk 'NR>1 {print $1}' gen.dat | tac >> capsule.srf
awk 'NR>1 {print $2}' gen.dat | tac >> capsule.srf
awk 'NR>1 {print $3}' gen.dat | tac >> capsule.srf

# Run Hypgen -- this generates the flowfield grid in 2D
echo "** HYPGEN"
~dasaund1/Chimera_Tools/hypgen <<-EOF1
capsule.srf
capsule.midplane
 0                  IFORM(0/1) 
 1, 1, 1               IZSTRT(1/2/-1),NZREG,KLAYER 
 241, 35.0, 1e-05, 0.0         NPZREG(),ZREG(),DZ0(),DZ1() 
 4, 4, 3, 3         IBCJA,IBCJB,IBCKA,IBCKB 
 1, 0.01, 800, 1         IVSPEC(1/2),EPSSS,ITSVOL,NSUB 
 2, 0.9              IMETH(0/2/3),SMU2 
 5.0, 5.0               TIMJ,TIMK 
 1, 0.3, 0.3            IAXIS(1/2),EXAXIS,VOLRES 
EOF1

echo "** GRIDED"

# Run "grided" the grid editing tool in CGT to reverse the indices for DPLR
~dasaund1/Chimera_Tools/grided <<-EOF2
capsule.midplane
    1    1     ITIN,ITOUT
    1     MOP
    0     NGSO
    3     IOP
    1     IYN
    5     IOP
    1    0    0     JREV,KREV,LREV
    0     IYN
2D.gu
    1     IMG
EOF2
