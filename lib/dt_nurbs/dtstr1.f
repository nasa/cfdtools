      subroutine dtstr1 ( cin  , x    , ndim , knot , mult , tknts, 
     1                    ctmp , cdim , indx , nreal, work , nwork, 
     2                    cout , ier   ) 

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
C  UPDATE DATE:     6/4/93                   UPDATE VERSION #: 2.01
C  UPDATE NOTICE: Problems reported with DTSSXT and -200 errors.
C
C
C===============================================================================

c     --------------
c     ... parameters
c     --------------

      integer           ndim          , mult(*)       , cdim          , 
     1                  indx(*)       , nreal         , nwork         , 
     2                  ier

      double precision  cin(*)        , x(ndim,2)     , knot(*)       ,
     1                  tknts(*)      , ctmp(cdim,2)  , work(*)       , 
     2                  cout(*)

c     ----------------------
c     ... internal variables
c     ----------------------

      integer           i             , idom          , ilc           ,
     1                  ilknew        , in            , incx          , 
     2                  ifrst         , iopt          , iposl         , 
     3                  iposr         , istble        , istrt         , 
     4                  itmp          , jdom          , kord          ,
     5                  knt           , mode          , nbs           , 
     6                  nbsnew        , ncoef         , ndom          , 
     7                  need          , nknot         , nkntot        , 
     8                  nktmp         , nlrep         , nskip         , 
     9                  nzero         , out

      double precision  relpr

      character*8       subnam

c     ------------
c     ... external
c     ------------

      double precision  dtmcon

c     =================================================================
     
      subnam = 'dtstrm  '
      ndom   = int ( cin(1) )
      ilknew = 3 + 3 * nreal
      in     = 2
      out    = 1
      incx   = 1
      ifrst  = 0
      nbsnew = 1
      mode   = 5
      need   = 0
      jdom   = 0

c     ------------------------------------------
c     ... loop through the independent variables
c     ------------------------------------------  
 
      do 40 idom = 1, ndom 

          kord  = int ( cin(2+idom) )
          nbs   = int ( cin(2+idom+ndom) )
          
          call dt2svy ( cin, idom, knot, mult, nknot, iopt, istble,
     1                  ier )

          if ( ier .ne. 0 ) then

              ier = -100
              call dterr ( mode, subnam, ier, need )
              return

          end if

          knot(nknot+1) = x(idom,1)
          knot(nknot+2) = x(idom,2)
          mult(nknot+1) = kord
          mult(nknot+2) = kord
          nkntot        = nknot + 2

c         -----------------------------------
c         ... sort the knots and permute mult
c         -----------------------------------

          iopt   = 0
          istble = 1
          call dtsrtn ( knot, nkntot, iopt, istble, tknts, ier )
      
          if ( ier .ne. 0 ) then
     
              ier = -100 
              call dterr ( mode, subnam, ier, need )
              return

          end if

          call d0prmi ( mult, nkntot, tknts, ier )

          if ( ier .ne. 0 ) then

             ier = -100
             call dterr ( mode, subnam, ier, need )
             return

          end if

c         ---------------------------------
c         ... remove extraneous information
c         ---------------------------------

          knt   = 1
          nknot = 1
          relpr = ( knot(nkntot) - knot(1) ) * dtmcon(5)
          relpr = 0.0d0

          do 10 knt = 2, nkntot

              if ( ( knot(knt) - knot(nknot) ) .gt. relpr ) then

                  nknot       = nknot + 1
                  knot(nknot) = knot(knt)
                  mult(nknot) = mult(knt)

              end if

              mult(nknot) = max ( mult(nknot), mult(knt) )

 10       continue

          nkntot = 0
          iposl  = 0
          iposr  = 0

c         -----------------------
c         ... define new knot set
c         -----------------------

          do 30 knt = 1, nknot 

              if ( x(idom,1) .eq. knot(knt) .and. iposl .eq. 0 )
     1            iposl = nkntot + 1
              if ( x(idom,2) .eq. knot(knt) .and. iposr .eq. 0 )
     1            iposr = nkntot + 1

             do 20 i = 1, mult(knt)

                  nkntot        = nkntot+ 1
                  tknts(nkntot) = knot(knt)

 20          continue

 30       continue

c         -------------------------------------------
c         ... set first index and number of b-splines
c         -------------------------------------------

          iposl  = min ( iposl, nkntot - kord )
          ifrst  = ifrst + ( iposl - 1 ) * nbsnew
          nbsnew = nbsnew * ( nkntot - kord )

c         -------------------------------
c         ... call dtolso to insert knots
c         -------------------------------

          if ( idom .eq. 1 ) then

              call dtoslo ( cin  , idom , tknts    , nkntot, cdim, 
     1                      work , nwork, ctmp(1,1), ier )

          else

              itmp = in
              in   = out
              out  = itmp
              call dtoslo ( ctmp(1,in), idom, tknts, nkntot, 
     1                      cdim      , work, nwork, ctmp(1,out), 
     2                      ier                        )

          end if

          if ( ier .ne. 0 ) then

              ier     = -100
              call dterr ( mode, subnam, ier, need )
              return

          end if

c         ------------------------
c         ... set cout information
c         ------------------------

          if ( x(idom,1) .ne. x(idom,2) ) then

               jdom                 = jdom + 1
               cout(2+jdom)         = dble ( float( kord ) )
               cout(2+nreal+jdom)   = dble ( float ( iposr - iposl ) )
               cout(2+2*nreal+jdom) = cout(2+jdom)
               nktmp                = iposr - iposl + kord 
               call dcopy ( nktmp, tknts(iposl), incx, cout(ilknew),
     1                      incx )
               ilknew               = ilknew + nktmp

          end if

  40  continue

c     ------------------------------------
c     ... build the remaining index vector
c     ------------------------------------

      indx(1)  = ifrst + 1
      nlrep    = 1
      nskip    = 1
      ncoef    = 1
      istrt    = 2
      nzero    = int ( cout(3+nreal) )
      ilc      = 3 + 3 * ndom
      jdom     = 0

      do 60 idom = 1, ndom

           ilc   = ilc + int(ctmp(2+idom,out) + ctmp(2+ndom+idom,out))
           ncoef = ncoef * int ( ctmp(2+idom+ndom,out) )

           if ( x(idom,1) .ne. x(idom,2) ) then

               do 50 i = istrt, nzero 

                   indx(i) = indx(i-nskip) + nlrep

  50           continue

               istrt = nzero + 1
               jdom  = jdom + 1
               if ( jdom .lt. nreal )nzero = nzero *
     1                                int ( cout(3+jdom+nreal) )
               nskip = nskip * int ( cout(2+nreal+jdom) )

           end if

           nlrep = nlrep * int ( ctmp(2+ndom+idom,out) )

 60   continue

      do 70 i = 1, abs (int ( cin(2) ))

          call dgthr ( nzero, ctmp(ilc,out), cout(ilknew), indx )
          ilc    = ilc + ncoef
          ilknew = ilknew + nzero

 70   continue

      cout(1) = dble ( float ( nreal ) )
      cout(2) = cin(2)

c     =================================================================

      return
      end
