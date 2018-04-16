      subroutine dtstrm ( cin, x, ndim, work, nwork, cout, ier )

C===============================================================================
C  PURPOSE:  
C
C  USAGE:
C
C  INPUTS:  
C
C
C
C  WORKING STORAGE:
C
C
C  OUTPUT:  
C
C  ROUTINES CALLED:
C
C  PROGRAMMER:
C
C  CREATION DATE: 
C
C  UPDATE DATE:   6/4/94                     UPDATE VERSION #: 2.01
C  UPDATE NOTICE:  Problems with DTSSXT and -200 errors.  
C
C
C===============================================================================

c     --------------
c     ... parameters
c     --------------

      integer           ndim          , nwork          , ier

      double precision  cin(*)        , x(ndim,*)      , work(*)      ,
     1                  cout(*)

c     ----------------------
c     ... internal variables
c     ----------------------
                                                        
      integer           ctmp         , kord         , idom         , 
     1                  ifail        , ilk          , indx         , 
     2                  knot         , mode         , mult         , 
     3                  nbs          , nbstmp       , ncoef        , 
     4                  ndom         , need         , nkmax        , 
     5                  nknts        , nreal        , nwoslo       ,
     6                  tknts        , wkoslo          

      character*8       subnam
 
c     =================================================================

      subnam  = 'dtstrm  '
      ier     = 0
      mode    = 1
      need    = 0
      cout(1) = -1.0d0

c     ------------------------
c     ... check valid c vector
c     ------------------------

      call dtschk ( cin, ifail )

      if ( ifail .ne. 0 ) then

          ier     = -55
          call dterr ( mode, subnam, ier, need )
          return

      end if

      ndom = int ( cin(1) )

c     --------------------
c     ... check valid ndim
c     --------------------

      if ( ndim .lt. ndom ) then

          ier     = -4
          call dterr ( mode, subnam, ier, need )
          return

      end if

c     ----------------------------
c     ... check interval endpoints
c     ----------------------------

      nreal = 0
      ilk   = 3 + 3 * ndom

      do 10 idom = 1, ndom

          kord = int ( cin(2+idom) )
          nbs  = int ( cin(2+idom+ndom) )

          if ( x(idom,1) .gt. x(idom,2) ) then

              ier     = -5
              call dterr ( mode, subnam, ier, need )
              return

          end if

          if ( x(idom,1) .lt. cin(ilk) .or.
     1         x(idom,2) .gt. cin(ilk+nbs) ) then

              ier = -50
              call dterr ( mode, subnam, ier, need )
              return

           end if

           if ( x(idom,1) .ne. x(idom,2) ) nreal = nreal + 1
           ilk = ilk + nbs + kord

 10   continue

      if ( nreal .eq. 0 ) then

          ier     = -5
          call dterr ( mode, subnam, ier, need )
          return

      end if

c     ------------------------------
c     ... determine workspace needed
c     ------------------------------

      ncoef  = 1
      nknts  = 0
      nkmax  = 0
      nwoslo = 1

      do 20 idom = 1, ndom

          nbstmp = int ( 2.0d0 * cin(2+idom) + cin(2+ndom+idom) )
          ncoef  = ncoef * nbstmp
          nknts  = nknts + nbstmp + int ( cin(2+idom) ) 
          nkmax  = max ( nkmax, nbstmp + int ( cin(2+idom) ) )
          nwoslo  = max ( nwoslo, nbstmp * ( int ( cin(2+idom) ) + 1 ) )

 20   continue

      nbstmp = 2 + 3 * ndom + abs (int ( cin(2) )) * ncoef + nknts
      need   = nkmax + 2 * nbstmp + max ( nwoslo, 2 * nkmax + ncoef )

      if ( nwork .lt. need ) then

          ier     = -3
          mode    = 2
          call dterr ( mode, subnam, ier, need )
          return

      end if 

c     ---------------------------------------------------------
c     ... partition workspace for lower level routine
c
c     tknts(nkmax)   - full knot vector used as input to dtolso
c     ctmp(nbstmp)   - temporary spline vector
c     knot(nkmax)    - distinct knot set for c
c     mult(nkmax)    - multiplicities of knots in c
c     indx(ncoef)    - index for trimmed spline coefficients
c     wkoslo(nwoslo) - dtoslo workspace
c
c     ---------------------------------------------------------

      tknts  = 1
      ctmp   = tknts + nkmax
      knot   = ctmp + 2 * nbstmp
      mult   = knot + nkmax
      indx   = mult + nkmax
      wkoslo = knot

c     ----------------------------
c     ... call lower level routine
c     ----------------------------

      call dtstr1 (  cin        , x           , ndim       ,
     1               work(knot) , work(mult)  , work(tknts), 
     2               work(ctmp ), nbstmp      , work(indx) , 
     3               nreal      , work(wkoslo), nwoslo     , 
     4               cout       , ier          )

c     =================================================================

      return
      end
