!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  These modules package the construction and discretization of analytic shapes commonly used for atmosphere entry capsules,
!  such as sphere-cones with given base radii, half cone angle, and so on.  Rotational symmetry is assumed.
!
!  Only forebodies are treated here.  Aft body generatrices can be more involved - see program CAPSULE_GRID by the same author.
!
!  Each configuration is expected to have a derived data type for its various parameters, a subroutine to calculate control
!  point coordinates defining its generatrix, and another subroutine to discretize the generatrix as needed for generating a
!  computational surface grid.
!
!  Coordinate system:       x increases aft from the nose, parallel to the symmetry axis;
!                           r is perpendicular to the x axis (normally the symmetry axis).
!
!                           The apex is typically at (0, 0), but need not be.
!
!  Elements:                module sphere_cone_parameters       Data structure for specifications and derived coordinates
!                           module biconic_parameters           Data structure likewise
!                           module sphere_cone                  Utilities for constructing and discretizing common generatrices
!  History:
!
!     12/13/10  D. A. Saunders  Initial implementation (sphere-cone case).
!     12/15/10     "     "      Incorporated the spherical section case into the sphere-cone case (cone angle = 0 is the flag).
!     10/24/11     "     "      Added axisymmetric biconic forebody utilities that make use of the sphere-cone utilities.
!     11/30/11     "     "      Added an option to round the vertex of a biconic forebody, requiring the new defining parameter
!                               %radius_vertex and several derived parameters in the biconic_type derived data type.
!     10/21/13     "     "      Two special cases needed separate handling: (1) the "Apollo" case (spherical section indicated via
!                               half_cone_angle = 0.) with radius_base = radius_nose on input can satisfy the latter condition only
!                               if radius_shoulder = 0; otherwise, the maximum radius is defined by radius_nose and radius_shoulder;
!                               (2) skirt_angle = 90. means skirt_length is interpreted as being along the skirt, not an x length;
!                               otherwise, only skirt_length = 0. could be meaningful; now, skirt_length > 0. can improve the
!                               curvature-based redistribution of the forebody that may be performed by CAPSULE_GRID because the
!                               forebody need not end at the shoulder-skirt tangency point.
!     12/23/13     "     "      Print the length of the conical flank of a sphere-cone.
!
!  Author:      David Saunders  ERC, Inc. at NASA Ames Research Center, Moffett Field, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   module sphere_cone_parameters

      type sphere_cone_type

!        Defining parameters:

         real :: x_nose, r_nose               ! Nose (apex) coordinates, normally (0, 0)
         real :: radius_base                  ! Maximum radius at base
         real :: radius_nose                  ! Radius of nose sphere
         real :: radius_shoulder              ! Radius of shoulder toroid; 0 radius is valid
         real :: half_cone_angle              ! Semi-cone angle >= 0 deg and < 90 deg;  0 deg => spherical section case
         real :: skirt_angle                  ! Taper angle at base => 0 deg; use smaller angle & length 0 to shorten the skirt
         real :: skirt_length                 ! Length in x from shoulder-skirt tangency point, >=< 0 m;
                                              ! length along skirt if skirt_angle = 90.
!        Derived parameters:

         real :: xt_nose, rt_nose             ! Tangency point on the nose; not defined for the spherical section case
         real :: xt_shoulder, rt_shoulder     ! Tangency point on the shoulder
         real :: xt_skirt, rt_skirt           ! Tangency point at/aft of max. r
         real :: x_base, r_base               ! Maximum-radius coordinates
         real :: x_skirt, r_skirt             ! Cut-off point; x_skirt >= x_base
         real :: xcenter_nose   ! With r_nose ! Center of nose sphere
         real :: rcenter_shoulder ! W/ x_base ! Center of shoulder circle

      end type sphere_cone_type

   end module sphere_cone_parameters

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   module biconic_parameters

      type biconic_type

!        Defining parameters:

         real :: x_nose, r_nose               ! Nose (apex) coordinates
         real :: radius_base                  ! Maximum radius at base
         real :: radius_nose                  ! Radius of nose sphere
         real :: radius_shoulder              ! Radius of shoulder circle; 0 radius is valid
         real :: radius_cone_juncture         ! Radius of base of forward cone
         real :: radius_vertex                ! Radius of circle that rounds off the vertex at the cone juncture; 0 => sharp
         real :: half_cone_angle_fore         ! Forward semi-cone angle > 0 deg and < 90 deg
         real :: half_cone_angle_aft          ! Aft     semi-cone angle > 0 deg and < 90 deg
         real :: skirt_angle                  ! Taper angle at base >= 0 deg; use smaller angle & length 0 to shorten the skirt
         real :: skirt_length                 ! Length in x from shoulder-skirt tangency point, >=< 0 m;
                                              ! length along skirt if skirt_angle = 90.
!        Derived parameters:

         real :: xt_nose, rt_nose             ! Tangency point on the nose; not defined for the spherical section case
         real :: xt_vertex_fore               ! Tangency point on the forward cone of the vertex-rounding circle
         real :: rt_vertex_fore
         real :: xt_vertex_aft                ! Tangency point on the aft cone of the vertex-rounding circle
         real :: rt_vertex_aft
         real :: xt_shoulder, rt_shoulder     ! Tangency point on the shoulder
         real :: xt_skirt, rt_skirt           ! Tangency point at/aft of max. r
         real :: x_base, r_base               ! Maximum-diameter coordinates
         real :: x_skirt, r_skirt             ! Cut-off point; x_skirt >= x_base
         real :: xcenter_nose   ! With r_nose ! Center of nose sphere
         real :: xcenter_vertex               ! Center of circle that rounds the cone-juncture vertex
         real :: rcenter_vertex
         real :: rcenter_shoulder ! W/ x_base ! Center of shoulder circle
         real :: x_cone_juncture              ! Axial  coordinate of the generatrix point at the cone juncture
         real :: r_cone_juncture              ! Radial coordinate of the generatrix point at the cone juncture
         real :: cone_length_fore             ! Axial length of the forward cone frustum
         real :: cone_length_aft              ! Axial length of the aft cone frustum from x_cone_juncture to xt_shoulder

      end type biconic_type

   end module biconic_parameters

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   module sphere_cone

      use sphere_cone_parameters  ! Data structure for sphere-cone parameters
      use biconic_parameters      ! Similar data structure for biconic forebodies
      use trigd

      implicit none

!     Internal constants:

      real, parameter, private :: half = 0.5, ninety = 90.0, one = 1.0, zero = 0.0

!     Internal variables:

      logical, private :: Apollo_type, hemisphere

      type (sphere_cone_type) :: sphere_cone_fore   ! For the biconic case, with zero shoulder radius
      type (sphere_cone_type) :: sphere_cone_aft    ! For the biconic case, with nose radius that of the cone juncture

!     Public procedures:

      public  :: sphere_cone_control_points     ! Derives tangency points, etc., from defining parameters
      public  :: sphere_cone_discretization     ! Roughly uniform distribution along generatrix, with key points captured
      public  :: biconic_control_points         ! Analogue for the biconic case: derive key coordinates/dimensions from inputs
      public  :: biconic_discretization         ! Roughly uniform distribution along generatrix, with key points captured

!     Private procedures:

      private :: distance

      contains

!        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         subroutine sphere_cone_control_points (capsule)

!        The capsule argument should be input with the eight defining parameters set (see sphere_cone_type).
!        Corresponding tangency points, etc., are calculated here.

!        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!        Arguments:

         type (sphere_cone_type), intent (inout) :: capsule  ! Particular instance of sphere-cone input+derived parameters

!        Local variables:

         real :: cone, cos_cone, L, sin_cone, Rb, Rn, Rs, theta

!        Execution:                 ! Draw a diagram with 8 key points (7 for Apollo) if you need to follow this.

         Rb = capsule%radius_base
         Rn = capsule%radius_nose
         Rs = capsule%radius_shoulder

         capsule%r_base            = Rb + capsule%r_nose
         capsule%rcenter_shoulder  = capsule%r_base - Rs
         capsule%xcenter_nose      = capsule%x_nose + Rn
!!!      capsule%rcenter_nose      = capsule%r_nose   ! Redundant

         cone                      = capsule%half_cone_angle
         Apollo_type               = cone == zero
         hemisphere                = Apollo_type .and. Rb == Rn  ! In which case Rb (maximum radius) really depends on Rs

         if (.not. Apollo_type) then
            cos_cone               = cosd (cone)
            sin_cone               = sind (cone)
            capsule%xt_nose        = capsule%xcenter_nose     - Rn*sin_cone
            capsule%rt_nose        = capsule%r_nose           + Rn*cos_cone
            capsule%rt_shoulder    = capsule%rcenter_shoulder + Rs*cos_cone
            capsule%xt_shoulder    = capsule%xt_nose          + (capsule%rt_shoulder - capsule%rt_nose) / tand (cone)
            capsule%x_base         = capsule%xt_shoulder      + Rs*sin_cone
         else
            if (.not. hemisphere) then
               capsule%rt_shoulder = capsule%r_nose + Rn*(Rb - Rs) / (Rn - Rs)
               capsule%xt_shoulder = capsule%x_nose + Rn - sqrt (Rn**2 - (capsule%rt_shoulder - capsule%r_nose)**2)
               capsule%x_base      = ((Rn - Rs)*capsule%xt_shoulder + Rs*(Rn + capsule%x_nose)) / Rn
            else                                         ! Input Rb = Rn, but this can be true only if Rs = 0
               if (Rs > zero) then                       ! We have to adjust Rb
                  theta = acos (Rs / (Rn - Rs))          ! Angle at nose circle center corresponding to shoulder tangency point
                  Rb    = Rs + (Rn - Rs)*sin (theta)     ! Adjusted maximum radius
                  capsule%radius_base = Rb
                  capsule%rt_shoulder = capsule%r_nose + Rn*sin (theta)
                  capsule%xt_shoulder = capsule%x_nose + Rn - Rn*cos (theta)
                  capsule%x_base      = capsule%x_nose + Rn - Rs
               else                                      ! Nose is a full semicircle
                  capsule%rt_shoulder = capsule%r_nose + Rn
                  capsule%xt_shoulder = capsule%x_nose + Rn
                  capsule%x_base      = capsule%xt_shoulder
               end if
            end if
         end if

         capsule%xt_skirt         = capsule%x_base           + Rs*sind (capsule%skirt_angle)
         capsule%rt_skirt         = capsule%rcenter_shoulder + Rs*cosd (capsule%skirt_angle)

         if (capsule%skirt_angle /= ninety) then             ! Skirt length is an x length (not along the arc)
            capsule%x_skirt       = capsule%xt_skirt         + capsule%skirt_length
            capsule%r_skirt       = capsule%rt_skirt         - capsule%skirt_length*tand (capsule%skirt_angle)
         else                                                ! Skirt length is along the arc
            capsule%x_skirt       = capsule%xt_skirt
            capsule%r_skirt       = capsule%rt_skirt         - capsule%skirt_length
         end if

         write (*, '(/, (a, 3x, 2f15.10))') &
            'Nose                   x and r:', capsule%x_nose, capsule%r_nose, &
            'Rnose, Rbase:                  ', Rn, Rb, &
            'Rshoulder, semi-cone angle:    ', Rs, cone, &
            'Skirt angle, skirt length:     ', capsule%skirt_angle, capsule%skirt_length
         if (.not. Apollo_type) then
            L = sqrt ((capsule%xt_shoulder - capsule%xt_nose)**2 + (capsule%rt_shoulder - capsule%rt_nose)**2)
            write (*, '(a, 3x, 2f15.10)') &
            'Cosine & sine semi-cone angle: ', cos_cone, sin_cone, &
            'Nose tangency          x and r:', capsule%xt_nose, capsule%rt_nose, &
            'Conical flank arc length:      ', L
         end if
         write (*, '(a, 3x, 2f15.10)') &
            'Shoulder tangency      x and r:', capsule%xt_shoulder, capsule%rt_shoulder, &
            'Skirt tangency         x and r:', capsule%xt_skirt, capsule%rt_skirt, &
            'Maximum-diameter       x and r:', capsule%x_base, capsule%r_base, &
            'Skirt cut-off          x and r:', capsule%x_skirt, capsule%r_skirt, &
            'Nose circle center     x and r:', capsule%xcenter_nose, capsule%r_nose, &
            'Shoulder circle center x and r:', capsule%x_base, capsule%rcenter_shoulder
         write (*, '(a)')

         end subroutine sphere_cone_control_points

!        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         subroutine sphere_cone_discretization (capsule, ngeom, xgeom, rgeom)

!        From the input control points, determine uniform distributions on all segments of the generatrix that produce the indicated
!        total number of discretized points (ngeom, nose to end of forebody) while capturing those key points precisely.  The input
!        value of ngeom is expected to be more than large enough for adequate subsequent interpolation of a computational grid,
!        which may or may not capture those points precisely.

!        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!        Arguments:

         type (sphere_cone_type), intent (in)  :: capsule       ! Following a call to sphere_control_points
         integer,                 intent (in)  :: ngeom         ! Total # geometry points to calculate along the generatrix
         real,                    intent (out) :: xgeom(ngeom), rgeom(ngeom)

!        Local variables:

         integer :: nnose, ncone, nshoulder, nskirt1, nskirt2   ! nskirt1 points are on the circle past max. radius (>= 0);
         integer :: i, ig                                       ! nskirt2 points are on the aft cone (possibly untapered; >= 0)
         real    :: snose, scone, sshoulder, sskirt1, sskirt2
         real    :: angle, anglen, angles, anglek
         real    :: deg2rad, delta_angle, dr, dsu, dx, piby2, stotal

!        Execution:

         piby2 = asin (one);  deg2rad = piby2 / ninety
         Apollo_type = capsule%half_cone_angle == zero

         if (Apollo_type) then
            scone = zero
            if (hemisphere) then
               if (capsule%radius_shoulder == zero) then
                  anglen = piby2
                  angles = zero
               else
                  anglen = acos (capsule%radius_shoulder / (capsule%radius_nose - capsule%radius_shoulder))
                  angles = piby2 - anglen
               end if
            else
               anglen = asin ((capsule%rt_shoulder - capsule%r_nose) / capsule%radius_nose)  ! For nose arc
               angles = piby2 - anglen                                                       ! For shoulder arc
            end if
         else
            angles = deg2rad * capsule%half_cone_angle
            anglen = piby2 - angles
            scone  = sqrt ((capsule%xt_shoulder - capsule%xt_nose)**2 + (capsule%rt_shoulder - capsule%rt_nose)**2)
         end if

         snose     = capsule%radius_nose * anglen
         sshoulder = capsule%radius_shoulder * angles

         anglek  = capsule%skirt_angle * deg2rad
         sskirt1 = capsule%radius_shoulder * anglek           ! >= 0

         if (capsule%skirt_angle /= ninety) then  ! Skirt length is in the x direction
            sskirt2 = capsule%skirt_length / cos (anglek)        ! >= 0
         else  ! Skirt length is along the skirt, not in the x direction
            sskirt2 = capsule%skirt_length
         end if

         stotal    = snose + scone + sshoulder + sskirt1 + sskirt2
         dsu       = stotal / real (ngeom - 1)   ! Nominal uniform spacing

         nnose     = nint (snose     / dsu) + 1
         nshoulder = nint (sshoulder / dsu) + 1

         nskirt1   = 0
         if (sskirt1 > zero) nskirt1 = max (nint (sskirt1 / dsu), 3) + 1
         nskirt2   = 0
         if (sskirt2 > zero) nskirt2 = max (nint (sskirt2 / dsu), 3) + 1

!        Make sure the total point count is ngeom:

         if (Apollo_type) then
            nnose = ngeom - (nshoulder - 1)
            if (nskirt1 > 0) nnose = nnose - (nskirt1 - 1)
            if (nskirt2 > 0) nnose = nnose - (nskirt2 - 1)
         else
            ncone = ngeom - (nnose + nshoulder - 2)
            if (nskirt1 > 0) ncone = ncone - (nskirt1 - 1)
            if (nskirt2 > 0) ncone = ncone - (nskirt2 - 1)
         end if

!        Discretize the nose arc:

         delta_angle = (snose / real (nnose - 1)) / capsule%radius_nose

         xgeom(1) = capsule%x_nose
         rgeom(1) = capsule%r_nose

         do i = 2, nnose - 1
            anglen   = delta_angle * real (i - 1)
            xgeom(i) = capsule%x_nose + capsule%radius_nose * (one - cos (anglen))
            rgeom(i) = capsule%r_nose + capsule%radius_nose * sin (anglen)
         end do

!        Discretize the cone portion unless it's Apollo-type:

         if (Apollo_type) then
            xgeom(nnose) = capsule%xt_shoulder
            rgeom(nnose) = capsule%rt_shoulder
            ig =  nnose
         else
            xgeom(nnose) = capsule%xt_nose
            rgeom(nnose) = capsule%rt_nose

            dx = (capsule%xt_shoulder - capsule%xt_nose) / real (ncone - 1)
            dr = (capsule%rt_shoulder - capsule%rt_nose) / real (ncone - 1)

            do i = 1, ncone - 2
               xgeom(nnose + i) = capsule%xt_nose + dx * real (i)
               rgeom(nnose + i) = capsule%rt_nose + dr * real (i)
            end do

            xgeom(nnose + ncone - 1) = capsule%xt_shoulder
            rgeom(nnose + ncone - 1) = capsule%rt_shoulder
            ig =  nnose + ncone - 1
         end if

!        Discretize the forward shoulder arc:

         delta_angle = (sshoulder / real (nshoulder - 1)) / capsule%radius_shoulder

         do i = 1, nshoulder - 2
            angle = angles - delta_angle * real (i)
            xgeom(ig + i) = capsule%x_base - capsule%radius_shoulder * sin (angle)
            rgeom(ig + i) = capsule%r_base - capsule%radius_shoulder * (one - cos (angle))
         end do

         ig = ig + nshoulder - 1
         xgeom(ig) = capsule%x_base
         rgeom(ig) = capsule%r_base

!        Discretize the aft shoulder arc, if any:

         if (nskirt1 > 0) then
            delta_angle = (sskirt1 / real (nskirt1 - 1)) / capsule%radius_shoulder

            do i = 1, nskirt1 - 2
               angle = delta_angle * real (i)
               xgeom(ig + i) = capsule%x_base + capsule%radius_shoulder * sin (angle)
               rgeom(ig + i) = capsule%r_base - capsule%radius_shoulder * (one - cos (angle))
            end do

            ig = ig + nskirt1 - 1
            xgeom(ig) = capsule%xt_skirt
            rgeom(ig) = capsule%rt_skirt
         end if

!        Discretize any conical skirt:

         if (nskirt2 > 0) then
            dx =  (capsule%x_skirt - capsule%xt_skirt) / real (nskirt2 - 1)
            dr = -(capsule%r_skirt - capsule%rt_skirt) / real (nskirt2 - 1)

            do i = 1, nskirt2 - 2
               xgeom(ig + i) = capsule%xt_skirt + dx * real (i)
               rgeom(ig + i) = capsule%rt_skirt - dr * real (i)
            end do

            xgeom(ngeom) = capsule%x_skirt
            rgeom(ngeom) = capsule%r_skirt
         end if

         end subroutine sphere_cone_discretization

!        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         subroutine biconic_control_points (capsule)

!        Derive the key coordinates and dimensions of an axisymmetric biconic forebody from the defining parameters.
!        The capsule argument should be input with those eleven defining parameters set (see biconic_type).
!        The strategy is to calculate the related parameters of two sphere-cone generatrices using the earlier utility.
!        The option to round the vertex where the two cones meet is performed as a final step if capsule%radius_vertex > 0.

!        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!        Arguments:

         type (biconic_type), intent (inout) :: capsule  ! Particular instance of biconic input+derived parameters

!        Local variables:

         real :: angle, h, t

!        Execution:                 ! Draw a diagram with 7 key points on the generatrix if you need to follow this.

         sphere_cone_fore%x_nose          = capsule%x_nose
         sphere_cone_fore%r_nose          = capsule%r_nose
         sphere_cone_fore%radius_nose     = capsule%radius_nose
         sphere_cone_fore%radius_base     = capsule%radius_cone_juncture
         sphere_cone_fore%half_cone_angle = capsule%half_cone_angle_fore
         sphere_cone_fore%radius_shoulder = zero
         sphere_cone_fore%skirt_angle     = zero
         sphere_cone_fore%skirt_length    = zero

         call sphere_cone_control_points (sphere_cone_fore)

         capsule%xt_nose          = sphere_cone_fore%xt_nose
         capsule%rt_nose          = sphere_cone_fore%rt_nose
         capsule%x_cone_juncture  = sphere_cone_fore%xt_shoulder
         capsule%r_cone_juncture  = sphere_cone_fore%rt_shoulder
         capsule%xcenter_nose     = sphere_cone_fore%xcenter_nose
         capsule%cone_length_fore = sphere_cone_fore%xt_shoulder - sphere_cone_fore%xt_nose

         sphere_cone_aft%radius_nose      = capsule%radius_cone_juncture / cosd (capsule%half_cone_angle_aft)
         sphere_cone_aft%x_nose           = capsule%radius_cone_juncture * tand (capsule%half_cone_angle_aft) + &
                                            capsule%x_cone_juncture - sphere_cone_aft%radius_nose
         sphere_cone_aft%r_nose           = capsule%r_nose
         sphere_cone_aft%radius_base      = capsule%radius_base
         sphere_cone_aft%half_cone_angle  = capsule%half_cone_angle_aft
         sphere_cone_aft%radius_shoulder  = capsule%radius_shoulder
         sphere_cone_aft%skirt_angle      = capsule%skirt_angle
         sphere_cone_aft%skirt_length     = capsule%skirt_length

         call sphere_cone_control_points (sphere_cone_aft)

         capsule%xt_shoulder      = sphere_cone_aft%xt_shoulder
         capsule%rt_shoulder      = sphere_cone_aft%rt_shoulder
         capsule%xt_skirt         = sphere_cone_aft%xt_skirt
         capsule%rt_skirt         = sphere_cone_aft%rt_skirt
         capsule%x_base           = sphere_cone_aft%x_base
         capsule%r_base           = sphere_cone_aft%r_base
         capsule%x_skirt          = sphere_cone_aft%x_skirt
         capsule%r_skirt          = sphere_cone_aft%r_skirt
         capsule%rcenter_shoulder = sphere_cone_aft%rcenter_shoulder
         capsule%cone_length_aft  = sphere_cone_aft%xt_shoulder - sphere_cone_fore%xt_shoulder

         write (*, '(/, (a, 3x, 2f15.10))') &
            'Nose                   x and r:', capsule%x_nose, capsule%r_nose, &
            'Rnose, Rbase:                  ', capsule%radius_nose, capsule%radius_base, &
            'Rshoulder, Rcone_juncture:     ', capsule%radius_shoulder, capsule%radius_cone_juncture, &
            'Fore and aft cone semi-angles: ', capsule%half_cone_angle_fore, capsule%half_cone_angle_aft, &
            'Fore and aft frustum lengths:  ', capsule%cone_length_fore, capsule%cone_length_aft, &
            'Skirt angle, skirt length:     ', capsule%skirt_angle, capsule%skirt_length, &
            'Nose tangency          x and r:', capsule%xt_nose, capsule%rt_nose, &
            'Cone juncture          x and r:', capsule%x_cone_juncture, capsule%r_cone_juncture, &
            'Shoulder tangency      x and r:', capsule%xt_shoulder, capsule%rt_shoulder, &
            'Skirt tangency         x and r:', capsule%xt_skirt, capsule%rt_skirt, &
            'Maximum-diameter       x and r:', capsule%x_base, capsule%r_base, &
            'Skirt cut-off          x and r:', capsule%x_skirt, capsule%r_skirt, &
            'Nose circle center     x and r:', capsule%xcenter_nose, capsule%r_nose, &
            'Shoulder circle center x and r:', capsule%x_base, capsule%rcenter_shoulder

!        Option to round the vertex where the two cones meet:

         if (capsule%radius_vertex > zero) then

             angle = half * (capsule%half_cone_angle_fore - capsule%half_cone_angle_aft)
             t     = capsule%radius_vertex * tand (angle)   ! Tangent length from vertex to tangency points
             angle = half * (capsule%half_cone_angle_fore + capsule%half_cone_angle_aft)
             h     = sqrt (t**2 + capsule%radius_vertex**2)
             capsule%xcenter_vertex = capsule%x_cone_juncture + h * sind (angle)
             capsule%rcenter_vertex = capsule%r_cone_juncture - h * cosd (angle)
             capsule%xt_vertex_fore = capsule%x_cone_juncture - t * cosd (capsule%half_cone_angle_fore)
             capsule%rt_vertex_fore = capsule%r_cone_juncture - t * sind (capsule%half_cone_angle_fore)
             capsule%xt_vertex_aft  = capsule%x_cone_juncture + t * cosd (capsule%half_cone_angle_aft)
             capsule%rt_vertex_aft  = capsule%r_cone_juncture + t * sind (capsule%half_cone_angle_aft)

             write (*, '(a, 3x, f15.10, /, (a, 3x, 2f15.10))') &
                'Rvertex:                       ', capsule%radius_vertex, &
                'Vertex circle center   x and r:', capsule%xcenter_vertex, capsule%rcenter_vertex, &
                'Vertex tangency (fore) x and r:', capsule%xt_vertex_fore, capsule%rt_vertex_fore, &
                'Vertex tangency (aft)  x and r:', capsule%xt_vertex_aft,  capsule%rt_vertex_aft
         else
             capsule%xt_vertex_fore = capsule%x_cone_juncture
             capsule%xt_vertex_aft  = capsule%x_cone_juncture
             capsule%xcenter_vertex = capsule%x_cone_juncture
             capsule%rt_vertex_fore = capsule%r_cone_juncture
             capsule%rt_vertex_aft  = capsule%r_cone_juncture
             capsule%rcenter_vertex = capsule%r_cone_juncture
         end if

         end subroutine biconic_control_points

!        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         subroutine biconic_discretization (capsule, ngeom, xgeom, rgeom)

!        From the input control points, determine uniform distributions on each portion of the generatrix that produce the
!        indicated total number of discretized points (ngeom, nose to end of body) while capturing those key points precisely.
!        The input value of ngeom is expected to be more than large enough for adequate subsequent interpolation of a
!        computational grid, which may or may not capture those points precisely.
!
!        Strategy:  Use the total biconic generatrix arc length to determine appropriate numbers of points on the fore and aft
!        portions, then use two sphere-cone discretizations to produce the biconic discretization (with some juggling of local
!        storage).  The aft portion has to be generated first, because there is no direct control over how many points will be
!        placed from the cone juncture to the end of body when the underlying aft sphere-cone is discretized.
!
!        Rounding of the cone-juncture vertex strategy (%radius_vertex > 0):
!
!        Retaining the roughly uniform arc-length spacing of the interim generatrix would typically leave very few points on
!        vertex-rounding circular arc.  Therefore, impose a minimum number of points, then redistribute the straight cone segments
!        on either side of the rounded vertex to recover the specified ngeom total number of points.

!        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!        Arguments:

         type (biconic_type), intent (in)  :: capsule       ! Following a call to biconic_control_points
         integer,             intent (in)  :: ngeom         ! Total # geometry points to calculate along the generatrix
         real,                intent (out) :: xgeom(ngeom), rgeom(ngeom)

!        Local constants:

         integer, parameter :: nvertex_min = 7  ! Minimum number of points on a discretized rounded cone-juncture vertex

!        Local variables:

         integer :: i, i1, i2, i3, i4, ncone_fore, ncone_aft, ngeom_fore, nvertex
         real    :: delta_angle, dr, dx, snose, scone_fore, scone_aft, sshoulder, sskirt1, sskirt2, su, svertex
         real    :: anglen, angles, anglek, anglev
         real    :: deg2rad, piby2, tol

         real, allocatable, dimension (:) :: xgeom_aft, rgeom_aft

!        Execution:

         piby2 = asin (one);  deg2rad = piby2 / ninety

         angles     = capsule%half_cone_angle_fore * deg2rad
         anglen     = piby2 - angles
         scone_fore = sqrt (capsule%cone_length_fore**2 + (capsule%r_cone_juncture - capsule%rt_nose)**2)
         scone_aft  = sqrt (capsule%cone_length_aft**2  + (capsule%rt_shoulder - capsule%r_cone_juncture)**2)
         snose      = capsule%radius_nose * anglen
         sshoulder  = capsule%radius_shoulder * angles
         anglek     = capsule%skirt_angle * deg2rad
         sskirt1    = capsule%radius_shoulder * anglek           ! >= 0

         if (capsule%skirt_angle /= ninety) then  ! Skirt length is in the x direction
            sskirt2 = capsule%skirt_length / cos (anglek)        ! >= 0
         else  ! Skirt length is along the skirt, not in the x direction
            sskirt2 = capsule%skirt_length
         end if

         allocate (xgeom_aft(ngeom), rgeom_aft(ngeom))

!        Aft cone/shoulder/skirt discretization:  the nose portion will be discarded.

         call sphere_cone_discretization (sphere_cone_aft, ngeom, xgeom_aft, rgeom_aft)

!        Locate the index of the cone juncture:

         tol = 1.e-5*snose
         do i = 2, ngeom
            if (abs (xgeom_aft(i) - capsule%x_cone_juncture) < tol) then
               ngeom_fore = i
               exit
            end if
         end do

!        Nose/forward cone discretization:

         call sphere_cone_discretization (sphere_cone_fore, ngeom_fore, xgeom, rgeom)

         xgeom(ngeom_fore:ngeom) = xgeom_aft(ngeom_fore:ngeom)
         rgeom(ngeom_fore:ngeom) = rgeom_aft(ngeom_fore:ngeom)

         deallocate (xgeom_aft, rgeom_aft)

!        Option to round the vertex where the two cones meet (messy when done as an afterthought):

         if (capsule%radius_vertex > zero) then  ! The rounding circle center and tangency points are in hand

            anglev  = (capsule%half_cone_angle_fore - capsule%half_cone_angle_aft) * deg2rad  ! Angle at center of vertex circle
            svertex = capsule%radius_vertex * anglev  ! Rounding arc length
            i = ngeom_fore
            su = distance (xgeom(i), rgeom(i), xgeom(i+1), rgeom(i+1))  ! Roughly uniform spacing on cone segments
            nvertex = max (nvertex_min, 1 + nint (svertex / su))   ! # pts. to interp. from forward tangent pt. to aft tangent pt.
            delta_angle = anglev / (deg2rad * real (nvertex - 1))  ! Back to degrees

!           First we have to reduce the numbers of points either side of the vertex to make room (more heuristics):

            do i = 2, ngeom_fore
               if (abs (xgeom(i) - capsule%xt_nose) < tol) then
                  i1 = i  ! Index of start of forward cone
                  exit
               end if
            end do

            ncone_fore = (ngeom_fore - i1) - nvertex / 2    ! Revised # pts. on forward cone
            i2 = i1 + ncone_fore - 1
            dx = (capsule%xt_vertex_fore - capsule%xt_nose) / real (ncone_fore - 1)
            dr = (capsule%rt_vertex_fore - capsule%rt_nose) / real (ncone_fore - 1)

            do i = i1 + 1, i2
               xgeom(i) = capsule%xt_nose + dx * real (i - i1)
               rgeom(i) = capsule%rt_nose + dr * real (i - i1)
            end do

!           Likewise for the revised aft cone:

            do i = ngeom_fore, ngeom
               if (abs (xgeom(i) - capsule%xt_shoulder) < tol) then
                  i4 = i  ! Index of end of aft cone
                  exit
               end if
            end do

            ncone_aft = (i4 - i1) - ncone_fore - nvertex + 3
            i3 = i4 - ncone_aft + 1
            dx = (capsule%xt_shoulder - capsule%xt_vertex_aft) / real (ncone_aft - 1)
            dr = (capsule%rt_shoulder - capsule%rt_vertex_aft) / real (ncone_aft - 1)

            do i = i3, i4
               xgeom(i) = capsule%xt_vertex_aft + dx * real (i - i3)
               rgeom(i) = capsule%rt_vertex_aft + dr * real (i - i3)
            end do

!           Now insert the rounded vertex arc points:

            do i = i2 + 1, i3 - 1
               anglev = ninety - capsule%half_cone_angle_fore + delta_angle * real (i - i2)
               xgeom(i) = capsule%xcenter_vertex - capsule%radius_vertex * cosd (anglev)
               rgeom(i) = capsule%rcenter_vertex + capsule%radius_vertex * sind (anglev)
            end do

         end if

         end subroutine biconic_discretization

!        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         real function distance (xl, yl, xr, yr)  ! Obvious purpose and arguments (private utility)
!        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         real, intent (in) :: xl, yl, xr, yr

         distance = sqrt ((xl - xr)**2 + (yl - yr)**2)

         end function distance

   end module sphere_cone
