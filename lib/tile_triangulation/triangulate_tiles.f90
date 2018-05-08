!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   subroutine triangulate_tiles (lunread, save_triangulation,        &
                                 ntile, ntri, tri_ids, tri_xyzs, conn)
!
!  Description:
!
!  This routine reads a file of Space Shuttle OML or IML data consisting of
!  tile ID numbers and (x,y,z) coordinates at multiple vertices (3 or more
!  of them per tile), and returns a triangulation formed by connecting the
!  vertices of each tile to its centroid.  Such a triangulation can allow
!  efficient searching of the data.
!
!  Input data format:
!
!     190002003     1299.006  -97.581  264.372
!     190002003     1294.764 -101.823  264.824
!     190002003     1299.006 -106.066  265.182
!     190002003     1303.249 -101.823  264.722
!     190002004     1292.642  -95.459  264.262
!     190002004     1288.400  -99.702  264.716
!     190002004     1292.642 -103.945  265.058
!     190002004     1296.885  -99.702  264.595
!         :             :        :        :
!     397538024001  1191.510   59.860  487.413
!     397538024001  1161.010   59.860  487.413
!     397538024001  1161.010   74.000  480.564
!         :             :        :        :
!     397538024002  1260.950  105.000  420.500
!     397538024002  1206.450  105.000  420.500
!     397538024002  1206.450  105.000  426.300
!         :             :        :        :
!     399415135      261.000   31.354  325.290
!     399415135      262.000   32.057  325.270
!     399415135      262.988   32.736  325.254
!     399415135      263.308   31.881  322.503
!     399415135      263.564   30.842  320.287
!     399415135      263.836   29.348  317.975
!     399415135      264.012   28.152  316.530
!     399415135      263.000   28.057  317.319
!     399415135      262.000   27.945  318.129
!     399415135      261.000   27.815  318.972
!
!  The vertices for a given tile are assumed to be in a sensible order
!  such that they form the edges of the tile when connected in order,
!  not any of the diagonals between non-contiguous vertices.
!
!  The tile IDs are in ascending order, more or less.
!  Note the example of three extra digits signifying partial tiles.
!  64-bit arithmetic must be used to read these as reals; integers can fail.
!
!  06/30/04  DAS  Initial implementation as part of associating each point of
!                 a computational surface grid with a tile ID.
!  07/07/04   "   Reading the first tile data line before the main loop cleaned
!                 up the awkward logic.
!  07/23/04   "   Option to suppress saving of the triangulation as a file.
!
!  David Saunders, ELORET/NASA Ames Research Center, Moffett Field, CA.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Arguments:

   integer, intent (in)  :: lunread   ! Logical unit for the input data file,
                                      ! which is assumed to be open upon entry
                                      ! and is closed here before the return
   logical, intent (in)  :: save_triangulation ! T|F; T measn the triangulation
                                      ! is written to disk via lunread
   integer, intent (out) :: ntile,  & ! # distinct tile IDs found
                            ntri      ! # triangles derived from the tiles;
                                      ! ntri = # lines found in the file
   real, pointer     :: tri_ids(:)    ! tri_ids(n) contains the tile ID of the
                                      ! nth triangle; IDs greater than 9 digits
                                      ! are changed to 9 digits, losing partial
                                      ! tile distinctions, but the IDs have are
                                      ! searched for tile thickness elsewhere,
                                      ! where they must be in ascending order,
                                      ! so we can't return oddball IDs here
   real, pointer     :: tri_xyzs(:,:) ! tri_xyzs(1:3,*) contains the (x,y,z)
                                      ! coordinates of the vertices of all the
                                      ! triangles derived from the tiles;
                                      ! in tile order, the input tile vertex
                                      ! (x,y,z)s are followed by their centroid,
                                      ! so the dimensions of tri_xyzs are
                                      ! (1:3, 1:ntri+ntile)
   integer, pointer  :: conn(:,:)     ! conn(1:3,n) contains pointers into
                                      ! tri_xyzs(:,*) for the vertices of the
                                      ! nth triangle
!  Local constants:

   integer, parameter :: lunerr = 6
   real,    parameter :: nine_digits = 999999999.,               &
                         one = 1.0, thousandth = 0.001, zero = 0.0

!  Local variables:

   integer :: i, ios, m, m1, n, nnode, np, nsame
   real    :: id, id_previous, rsame, sumx, sumy, sumz, x, y, z

!  Execution:
!  ----------

!  Scan the data a first time to count the number of lines/nodes/vertices
!  and the number of distinct tile IDs:

   m = 0  ! Tile count
   n = 0  ! Line count
   id_previous = -one

   do ! Until EOF
      read (lunread, *, iostat=ios) id
      if (ios < 0) exit ! EOF

      n = n + 1
      if (id > nine_digits) id = real (nint (id * thousandth))
      if (id /= id_previous) then  ! New ID
         id_previous = id
         m = m + 1
      end if
   end do

   ntile = m
   ntri  = n            ! # triangles = # polygon edges = # polygon vertices
   nnode = ntri + ntile ! # (x,y,z) coordinates

   allocate (tri_xyzs(3,nnode), conn(3,ntri), tri_ids(ntri))

!  Rewind and read the vertex coordinates into the packed array that includes
!  the triangle centroids:

   rewind (lunread)

!  Initial condition:

   read (lunread, *, iostat=ios) id, sumx, sumy, sumz
   if (ios /= 0) then
      write (lunerr, '(/, a)') ' Trouble reading (x,y,z) at line # 1'
      stop
   end if

   if (id > nine_digits) id = real (nint (id * thousandth)) ! Avoid 12-digit IDs
   tri_ids(1)  = id
   id_previous = id

   tri_xyzs(1,1) = sumx
   tri_xyzs(2,1) = sumy
   tri_xyzs(3,1) = sumz

   conn(1,1) = 1   ! Vertex 1 of triangle 1 (pointer to tri_xyzs(:,1))

   m     = 1  ! Packed coordinates count
   m1    = 1  ! Value of m for vertex 1 of the first triangle of a tile
   nsame = 1  ! Same ID count
   np    = 1  ! n - 1 in the following loop

!  Main loop over tile polygon vertices:

   do n = 2, ntri

      read (lunread, *, iostat=ios) id, x, y, z
      if (ios /= 0) then
         write (lunerr, '(/, a, i7)') ' Trouble reading (x,y,z) at line #', n
         stop
      end if

      if (id > nine_digits) id = real (nint (id * thousandth))
      tri_ids(n) = id

      if (id == id_previous) then
         nsame = nsame + 1
         m     = m + 1
         conn(2,np) = m         ! Vertex 2 of previous triangle
         conn(1,n)  = m         !    "   1 of current     "
         sumx  = sumx  + x;  tri_xyzs(1,m) = x
         sumy  = sumy  + y;  tri_xyzs(2,m) = y
         sumz  = sumz  + z;  tri_xyzs(3,m) = z
      else ! New ID; deal with previous ID
         id_previous = id
         m  = m + 1             ! Index of tile centroid as a triangle vertex
         conn(2,np) = m1        ! 2nd vertex of last triangle = 1st v. of first
         conn(3,n-nsame:np) = m ! All triangles of a tile use its centroid
         rsame = one / real (nsame)
         nsame = 1
         tri_xyzs(1,m) = sumx * rsame
         tri_xyzs(2,m) = sumy * rsame
         tri_xyzs(3,m) = sumz * rsame

!        Set up for the next tile:

         m    = m + 1
         m1   = m               ! Index of first triangle for next tile
         conn(1,n) = m          ! Vertex 1 of triangle n
         conn(2,n) = m + 1      ! Vertex 2
         sumx = x;   tri_xyzs(1,m) = x
         sumy = y;   tri_xyzs(2,m) = y
         sumz = z;   tri_xyzs(3,m) = z
      end if

      np = n

   end do

   m = m + 1   ! Deal with the final tile
   n = np      ! Make sure of it - should be ntri
   conn(2,n) = m1
   conn(3,n-nsame+1:n) = m  ! Note the +1 because n wasn't incremented
   rsame = one / real (nsame)
   tri_xyzs(1,m) = sumx * rsame
   tri_xyzs(2,m) = sumy * rsame
   tri_xyzs(3,m) = sumz * rsame

   if (m /= nnode) write (lunerr, '(/, a, 2i8)') &
      ' Triangulate_tiles warning:  m, nnode =', m, nnode

   close (lunread)

   if (save_triangulation) then

      open (lunread, file='triangulated_tiles.tecplot', status='unknown')

      write (lunread, '(a)')                                  &
         'TITLE = "Triangulation of Shuttle Tiles, OML"',     &
         'VARIABLES = "x" "y" "z"',                           &
         'ZONE T = ""'
      write (lunread, '(a, i7, a, i7, a)')                    &
         'N =', nnode, ', E =', ntri, ', ZONETYPE=FETriangle'
      write (lunread, '(a)')                                  &
         'DATAPACKING=POINT', 'DT=(SINGLE SINGLE SINGLE)'
      write (lunread, '(1p, 3e15.7)') tri_xyzs
      write (lunread, '(i6, 2i7)') conn

      close (lunread)

   end if

   end subroutine triangulate_tiles
