!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   subroutine extend_name (root, extension, new_name, new_length)
!
!  Concatenate a root path or directory name (or group name in HDF5 lingo) and a
!  next-lower-level extension to the root name so as to form an extended path/
!  directory/group name.  The root and extension are both trimmed here, and the
!  extended name is returned along with its trimmed length.  Example:
!
!     call numbered_name ('Grid # ', 1, grid_name, len_name)               gives
!
!          grid_name(1:len_name) <-- 'Grid # 1'                            then
!
!     call extend_name ('DPLR', grid_name(1:len_name), pgrh%grid_name, l)  gives
!
!          pgrh%grid_name(1:l)   <-- 'DPLR/Grid # 1'                       then
!
!     call extend_name (pgrh%grid_name, 'Coords', pgrh%coords_name, l)     gives
!
!          pgrh%coords_name(1:l) <-- 'DPLR/Grid # 1/Coords'
!
!  Thus, we do little more than insert a / between the root and the extension to
!  form an extended name.  Passing trimmed input arguments (or not) is up to the
!  programmer, as is ensuring that new_name is long enough for the extension.
!
!  01/03/13  D.A.S.  Initial implementation to avoid lots of concatenations in
!                    the DPLR flow solver's I/O routines.
!
!  Author:  David Saunders, ERC, Inc. at NASA Ames Research Center, California.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Arguments:

   character (*), intent (in)  :: root        ! Leading part of the desired name
   character (*), intent (in)  :: extension   ! Trailing  "   "   "   "   "   "
   character (*), intent (out) :: new_name    ! Desired name, left justified
   integer,       intent (out) :: new_length  ! Length of the new name, with any
                                              ! trailing characters blanked
!  Local variables:

   integer       :: l1, l2

!  Execution:

!  We could use concatenation, but choose not to because of the trimming needed.

   l1 = len_trim (root);       new_name(1:l1)  = root(1:l1);   l1 = l1 + 1
   l2 = len_trim (extension);  new_name(l1:l1) = '/';  new_length = l1 + l2

   new_name(l1+1:) = extension(1:l2)  ! Pads new_name with blanks

   end subroutine extend_name
