      subroutine dtcsiz (cary,isize,ierr)
c***********************************************************************
c
c FUNCTION-
c          Determine size of C array
c
c AUTHORS-
c          Bob Ames             creation  Feb. 1993
c
c INPUT-
c         cary   = c-array containing dtrc format b-spline data
c OUTPUT-
c         isize  = number of effective data in c-array
c         ierr   = error flag
c                = 0   no error detected
c                = -1  more than 2 independent variables
c                = -2  incorrect number of dimensions
c                      second element of c-array greater than 3 or
c                      less than -4
c
c REFERENCES-
c          1.  
c
c NOTES- Based of DTWRCA.  Currently only supports number
c        of independant variables equal to  1 or 2 (curves and surfaces)
c				
c
c
      double precision cary(*)
c
c
c ********************** start of executable code **********************
c
      ierr=0
      nlns=0
c
c *****************************************************************
c...Analysize c-array
c
c...No. of independent variables 
c
      npv =int(cary(1))
c
      if(npv.lt.1.or.npv.gt.2) then
        ierr=-1
        return
      endif
c
c...No. of dependent variable and physical dimension
c
      if(cary(2).gt.0.0) then
        iphy=int(cary(2))
        ndim=iphy
      elseif (cary(2).lt.0.0) then
        iphy=int(abs(cary(2)))-1
        ndim=iphy+1
      endif

c *****************************************************************
c...Object is a surface
c
      if(npv.eq.2) then
c
c...No. of control point in t and s
c
        nlns=int(cary(5))
        npts=int(cary(6))
c
        nctrp=nlns*npts
c
c...Degree in t and s
c
        ndegl=int(cary(3))-1
        ndegp=int(cary(4))-1
c
c...No. of knots for t and s knot vectors
c
        nkt=ndegl+nlns+1
        nks=ndegp+npts+1
c
c...Determine array size for c-array
c
        isize=(nlns*npts)*ndim+nks+nkt+8
c
c *****************************************************************
c...Object is a curve
c
      elseif (nlns.eq.0) then
c
c...No. of control point in s
c
        npts=int(cary(4))
c
c...degree in s
c
        ndegl=0
        ndegp=int(cary(3))-1
c
c...no of knots for s knot vectors
c
        nkt=0
        nks=ndegp+npts+1
c
c...Determine array size for c-array
c
        isize=npts*ndim+nks+5
c
      endif
c
      return
      end
