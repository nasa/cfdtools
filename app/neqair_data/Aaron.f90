program HelloWorld

!  Aaron Brandis provided this approach to converting certain data in Tecplot format to NEQAIR's LOS.dat format.

   IMPLICIT NONE
   real, dimension(1000) :: dDistx, distY, T, Tv, nCO2, nCO, nCOp, nC2, nN2, nO2, nO2p, nNO, nNOp, nCN, nC, nCp, nN, nNp, nO, nOp, nAr, ne, totN   
   real, dimension(1) ::  EAST
   INTEGER :: i, j, k, ier, N, nlines, NNN
   CHARACTER :: head1*8, head2*8, Std*1, filename*32
   ier = 0
   nlines = 0
   NNN = 1
   EAST = 0.0
   open(unit=10, file="Chem3-18spBC+CEA-6.00kms-0.10Torr-StagLine.dat",form='formatted')
   read(10,'(a)') head1
   write(*,*) head1
   read(10,'(a)') head2
   write(*,*) head2
   do while(.true.)
          nlines = nlines+1
      read(10,'(e20.18,21e20.17)',iostat=ier) dDistx(nlines), distY(nlines), T(nlines), Tv(nlines), nCO2(nlines), nCO(nlines), nCOp(nlines), nC2(nlines), nN2(nlines), nO2(nlines), nO2p(nlines), nNO(nlines), nNOp(nlines), nCN(nlines), nC(nlines), nCp(nlines), nN(nlines), nNp(nlines), nO(nlines), nOp(nlines), ne(nlines)  
        if (ier/=0) exit
   enddo
   close(10)
   
   do i = 1,nlines
      totN(i) = (nCO2(i) + nCO(i) + nCOp(i) + nC2(i) + nN2(i) + nO2(i) + nO2p(i) + +nNO(i) + nNOp(i) + nCN(i) + nC(i) + nCp(i) + nN(i) + nNp(i)  +nO(i) + nOp(i) + ne(i))/1000000
      write(*,*) totN(i)
      write (filename, 12) i
12    format ('LOS-', I0, '.dat')
      open(unit=11, access='sequential', file=filename)
     write(11,'(a80)') '********************************************************************************'
     write(11,*) '              LOS file for Sample Cases 1 and 2'
     write(11,*) ''
     write(11,*) '       An unlimited number of comment lines can go here.'
     write(11,*) ''
     write(11,*) '       Enter Data AFTER the data-format lines!'
     write(11,*) ''
     write(11,*) '   (1) Enter species in any order; limited to atoms, diatomics, triatomics,'
     write(11,*) '       atomic ions, diatiomic ions, and electrons. Left-justify the species'
     write(11,*) '       symbols in the fields.  Dimensioned up to 25 species.  End entry'
     write(11,*) '       with a blank line.'
     write(11,*) ''
     write(11,*) '   (2) Properties entered at each grid point along line-of-sight.  The'
     write(11,*) '       properties apply to the layer between the grid point and the'
     write(11,*) '       previous grid point.  Thus, the properties at the first grid point'
     write(11,*) '       are not used.  This grid point only establishes the origin of the'
     write(11,*) '       line-of-sight.'
     write(11,*) ''
     write(11,*) '   (3) Enter species number densities [cm-3] in the same order that the species'
     write(11,*) '       symbols are entered.  End data entry at each grid point with a blank'
     write(11,*) '       line.'
     write(11,*) ''
     write(11,*) '   (4) End line-of-sight data entry, with a line of 0s as shown.'
     write(11,*) ''
     write(11,'(a80)') '********************************************************************************'
     write(11,*) '        aaaaaaaa       aaaaaaaa       aaaaaaaa       aaaaaaaa    (2x,(7x,a8))'
     write(11,*) '        E-             C              C+             CN      :Species Symbols.'
     write(11,*) '        CO             CO2            C2             N '
     write(11,*) '        N+             NO             NO+            N2 '
     write(11,*) '        O              O+             O2'
     write(11,*) ' '
     write(11,'(a80)') '--------------------------------------------------------------------------------'
     write(11,*) '   no.   x,cm   total partcc       t        tr        tv        te  (i5,f8.3,'
     write(11,'(a80)') 'iiiii rrrrrrr rrrrrrrrrrrrrr rrrrrrrrr rrrrrrrrr rrrrrrrrr rrrrrrrrre15.6,4f10.1'
     write(11,*) '      rrrrrrrrrrrrrr rrrrrrrrrrrrrr rrrrrrrrrrrrrr rrrrrrrrrrrrr   (6x,4e15.6)'
     write(11,*) '      Include these 9 lines (from --- to --- lines) for first grid point only!!'
     write(11,*) '      End each grid point entry with a blank line.'
     write(11,*) '      End data file with a line of zeros as shown on the next line.'
     write(11,*) '   0     0.0            0.0       0.0       0.0       0.0       0.0'
     write(11,'(a80)') '--------------------------------------------------------------------------------'
     write(11,'(i5,e15.7,e15.7,4f10.1)') NNN, EAST, totN(i), T(i), T(i), Tv(i), Tv(i)
     write(11,'(5x,4e15.7)') ne(i)/1000000, nC(i)/1000000, nCp(i)/1000000, nCN(i)/1000000
     write(11,'(5x,4e15.7)') nCO(i)/1000000, nCO2(i)/1000000, nC2(i)/1000000, nN(i)/1000000
     write(11,'(5x,4e15.7)') nNp(i)/1000000, nNO(i)/1000000, nNOp(i)/1000000, nN2(i)/1000000
     write(11,'(5x,3e15.7)') nO(i)/1000000, nOp(i)/1000000, nO2(i)/1000000
     write(11,*) ''
     write(11,'(i5,e15.7,e15.7,4f10.1)') NNN+1, EAST+10.16, totN(i), T(i), T(i), Tv(i), Tv(i)
     write(11,'(5x,4e15.7)') ne(i)/1000000, nC(i)/1000000, nCp(i)/1000000, nCN(i)/1000000
     write(11,'(5x,4e15.7)') nCO(i)/1000000, nCO2(i)/1000000, nC2(i)/1000000, nN(i)/1000000
     write(11,'(5x,4e15.7)') nNp(i)/1000000, nNO(i)/1000000, nNOp(i)/1000000, nN2(i)/1000000
     write(11,'(5x,3e15.7)') nO(i)/1000000, nOp(i)/1000000, nO2(i)/1000000
     write(11,*) ''
     write(11,*) '   0     0.0            0.0       0.0       0.0       0.0       0.0'
      close(11)
   enddo
   write(*,*) nlines
   write(*,*) dDistx(1)

end
