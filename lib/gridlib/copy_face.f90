!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine copy_face (ni_in, nj_in, nk_in, x_in, y_in, z_in,             &
                            iface,                                             &
                            ni_out, nj_out, x_out, y_out, z_out)
!
!     Transcribe a volume grid block face to a surface patch.  In order to
!     preserve handedness, the patch dimensions should be set up as shown in
!     the following code from program GRID_FACES.  This requires COPY_FACE to
!     transpose J faces.
!
!        id = (iface + 1) / 2 ! Convert face # to dimension 1, 2, or 3
!        i1 = mod (id, 3) + 1 ! 1st dimension of face iface
!        i2 = mod (i1, 3) + 1 ! 2nd ....................
!
!        Variables i1 and i2 come out cyclic for id = 1, 2, 3.
!        I.e., they are 2,3 and 3,1 and 1,2 respectively, meaning
!        nj,nk and nk,ni and ni,nj respectively.
!
!     09/09/05  David Saunders  Library version of a GRID_FACES utility.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Arguments:

      integer, intent (in) :: &
         ni_in, nj_in, nk_in                 ! Input block dimensions

      real, dimension (ni_in, nj_in, nk_in), intent (in) :: &
         x_in, y_in, z_in                    ! Input block coordinates

      integer, intent (in) :: &
         iface,               &              ! Block face to copy
         ni_out, nj_out                      ! Output face dimensions

      real, dimension (ni_out, nj_out), intent (out) :: &
         x_out, y_out, z_out                 ! Output surface grid

!     Local variables:

      real, allocatable :: temp(:,:)

!     Execution:

      select case (iface)

         case (1)

            x_out = x_in (1, 1:nj_in, 1:nk_in)
            y_out = y_in (1, 1:nj_in, 1:nk_in)
            z_out = z_in (1, 1:nj_in, 1:nk_in)

         case (2)

            x_out = x_in (ni_in, 1:nj_in, 1:nk_in)
            y_out = y_in (ni_in, 1:nj_in, 1:nk_in)
            z_out = z_in (ni_in, 1:nj_in, 1:nk_in)

         case (3)

            allocate (temp(ni_in,nk_in))

            temp  = x_in (1:ni_in, 1, 1:nk_in)
            x_out = transpose (temp)
            temp  = y_in (1:ni_in, 1, 1:nk_in)
            y_out = transpose (temp)
            temp  = z_in (1:ni_in, 1, 1:nk_in)
            z_out = transpose (temp)

            deallocate (temp)

         case (4)

            allocate (temp(ni_in,nk_in))

            temp  = x_in (1:ni_in, nj_in, 1:nk_in)
            x_out = transpose (temp)
            temp  = y_in (1:ni_in, nj_in, 1:nk_in)
            y_out = transpose (temp)
            temp  = z_in (1:ni_in, nj_in, 1:nk_in)
            z_out = transpose (temp)

            deallocate (temp)

         case (5)

            x_out = x_in (1:ni_in, 1:nj_in, 1)
            y_out = y_in (1:ni_in, 1:nj_in, 1)
            z_out = z_in (1:ni_in, 1:nj_in, 1)

         case (6)

            x_out = x_in (1:ni_in, 1:nj_in, nk_in)
            y_out = y_in (1:ni_in, 1:nj_in, nk_in)
            z_out = z_in (1:ni_in, 1:nj_in, nk_in)

      end select

      end subroutine copy_face
