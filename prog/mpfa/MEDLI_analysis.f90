!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   program MEDLI_analysis

!  Description:
!
!     This program manipulates thermocouple and/or heat sensor data from arc-jet
!  or flight experiments in ways that are expected to evolve.  Initially, it
!  simply regularizes the TC data and applies some shift or other transformation
!  to the heat sensor data that best matches it to the TC data.  See internal
!  procedures below for details.
!
!  History:
!
!     01/18/11  D.Saunders/T.White  Initial framework to regularize TC data.
!
!  Author: David Saunders, ERC, Inc./NASA Ames Research Center, Moffett Field CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Constants:

   integer, parameter :: &
      lundat = 1,        &  ! Input dataset
      lunkbd = 5,        &  ! Keyboard
      luncrt = 6            ! Screen

   real, parameter ::    &
      Thigh  = 3000.,    &  ! Limit on reasonable temperature readings (C or K)
      Tlow   = 0.,       &  ! Lower limit likewise
      told   = 0.01,     &  ! Tolerance for merging TC depths into grid, inches
      tolt   = 0.02,     &  ! Tolerance for merging times into time grid, secs
      zero   = 0.

!  Variables (likely to be shared by internal procedures):

   integer :: &
      ios, ireading, iTC, nperTC, nperTCmax, numsensors, numTCs, &
      ndgrid, ndgridu, ntgrid, ntgridu

   real :: &
      depth, temp, time

   integer, allocatable :: &
      negative_depth(:), ntemp_too_high(:), ntemp_too_low(:), numreadings(:)

   character :: &
      answer*1, filename*80

   real, allocatable :: &
      depth_raw(:,:), temp_raw(:,:), time_raw(:,:),           &
      tgrid(:), dgrid(:), tmesh(:,:), dmesh(:,:), fmesh(:,:), &
      tinterp(:,:), unused(:)

!  Execution:

   call read_TC_and_sensor_data ()

   call regularize_TC_data ()

!  Internal procedures for program MEDLI_analysis, in order of use:

   contains

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine read_TC_and_sensor_data ()

!     Prompt for and read a file containing temperature time histories for some
!     number of thermocouples and an optional heat sensor.  The data format is
!     simply:
!
!        # Time, sec   Depth, mm   Temperature, deg C
!        0             2.61        18.845
!        1.001         2.61        19.080
!         :             :            :
!         :             :            :
!        1183.998      2.61        70.244
!        0             4.68        17.857
!        1.001         4.68        14.408
!         :             :            :
!         :             :            :
!
!     The dataset is scanned once to report findings and allocate storage before
!     being reread.
!
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Execution:

      ios = 1
      do while (ios /= 0)
         write (luncrt, '(a)', advance='no') 'Input dataset: '
         read  (lunkbd, '(a)') filename
         open  (lundat, file=filename, status='old', iostat=ios)
         if (ios /= 0) &
            write (luncrt, '(2a)') 'Trouble opening ', trim (filename)
      end do

      numTCs    = 0
      nperTCmax = 0

      read (lundat, '(a)')  ! Skip title
      do while (ios == 0)
         read (lundat, *, iostat=ios) time
         if (ios < 0) exit  ! EOF
         if (time == zero) then 
            nperTCmax = max (nperTC, nperTCmax)
            nperTC = 1
            numTCs = numTCs + 1
         else
            nperTC = nperTC + 1
         end if
      end do
      nperTCmax = max (nperTC, nperTCmax)

      write (luncrt, '(a, i6)') &
         'Number of time series found:  ', numTCs, &
         'Largest # readings per series:', nperTCmax

      allocate (time_raw(nperTCmax,numTCs), depth_raw(nperTCmax,numTCs), &
                temp_raw(nperTCmax,numTCs), &
                numreadings(numTCs), negative_depth(numTCs), &
                ntemp_too_low(numTCs), ntemp_too_high(numTCs))

      rewind (lundat)
      read   (lundat, '(a)')
      read   (lundat, *, iostat=ios) time, depth, temp
      iTC = 0

      do  ! Until EOF
         if (ios < 0) then
            numreadings(iTC) = ireading
            exit
         end iF
         if (time == 0) then
            if (iTC > 0) numreadings(iTC) = ireading
            ireading = 0
            iTC = iTC + 1
         end if
         ireading = ireading + 1
         time_raw(ireading,iTC)  = time
         depth_raw(ireading,iTC) = depth
         temp_raw(ireading,iTC)  = temp
         read (lundat, *, iostat=ios) time, depth, temp
      end do

      close (lundat)

      write (luncrt, '(/, a)') &
         'Instrument  # readings  # too low  # too high  # depth < 0'

      ntemp_too_low(:)  = 0
      ntemp_too_high(:) = 0
      negative_depth(:) = 0

      do iTC = 1, numTCs
         do ireading = 1, numreadings(iTC)
            if (temp_raw(ireading,iTC) <= Tlow)  &
               ntemp_too_low(iTC)  = ntemp_too_low(iTC)  + 1
            if (temp_raw(ireading,iTC) >= Thigh) &
               ntemp_too_high(iTC) = ntemp_too_high(iTC) + 1
            if (depth_raw(ireading,iTC) <= zero) &
               negative_depth(iTC) = negative_depth(iTC) + 1
         end do
         write (luncrt, '(i10, i11, i12, i12, i13)') &
            iTC, numreadings(iTC), ntemp_too_low(iTC), &
            ntemp_too_high(iTC), negative_depth(iTC)
!!       write (10, '(i6, f10.3, f14.9, f16.6)') &
!!          (ireading, time_raw(ireading,iTC), depth_raw(ireading,iTC), &
!!           temp_raw(ireading,iTC), ireading = 1, numreadings(iTC))
      end do

      numsensors = 0
      write (luncrt, '(a)', advance='no') 'Is the last one a heat sensor? '
      read  (lunkbd, '(a)') answer
      if (answer == 'y' .or. answer == 'Y') then
         numsensors = 1
         numTCs = numTCs - 1
      end if

      end subroutine read_TC_and_sensor_data

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine regularize_TC_data ()

!     Interpolate TC measurements to a rectangular grid in (time, depth) space.
!     Initial scheme:  The largest number of readings for any TC suggests a
!     likely uniform grid size, then all the actual readings are merged, with
!     duplicates eliminated based on a small tolerance.  The depth direction is
!     handled similarly.  Obviously bad readings (negative, way too large) are
!     adjusted to [what? last apparently good reading?].  Then 1-D monotonic
!     spline interpolation is performed in the time and depth directions resp.

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Local constants:

      character, parameter :: mono = 'M'  ! Monotonic spline interpolations

!     Local variables:

      integer :: imost, i, j, nmergemax, nreadings
      real    :: dd, dt, previous_good
      logical :: new

!     Execution:

!     Adjust any obviously bad readings to avoid silly interpolations:

      nreadings = 0
      do iTC = 1, numTCs + numsensors
         if (numreadings(iTC) > nreadings) then
            nreadings = numreadings(iTC)
            imost     = iTC
         end if
         previous_good = zero
         do ireading = 1, numreadings(iTC)
            temp = temp_raw(ireading,iTC)
            if (temp > Thigh .or. temp < Tlow) then
               temp_raw(ireading,iTC) = previous_good
            else
               previous_good = temp
            end if
         end do
      end do

      do iTC = numTCS + 1, numTCS + numsensors
         do ireading = 1, numreadings(iTC)
            depth_raw(ireading,iTC) = max (zero, depth_raw(ireading,iTC))
         end do
      end do

      write (luncrt, '(a)', advance='no') &
         'Uniform grid size (ntimes & ndepths): '
      read  (lunkbd, *) ntgridu, ndgridu

      nmergemax = ntgridu + nreadings * (numTCs + numsensors)
      allocate (tgrid(nmergemax))

      nmergemax = ntgridu + numTCs + nreadings * numsensors
      allocate (dgrid(nmergemax))

      dt = time_raw(nreadings,imost) / real (ntgridu - 1)
      write (6, *) 'dt,nreadings,imost: ', dt,nreadings,imost
      write (6, *) 'time_raw(nreadings,imost): ', time_raw(nreadings,imost)

      do i = 1, ntgridu
         tgrid(i) = dt * real (i - 1)
      end do

      dd = depth_raw(1,numTCs) / real (ndgridu - 1)
      write (6, *) 'dd,depth_raw(1,numTCs): ', dd,depth_raw(1,numTCs)

      do j = 1, ndgridu
         dgrid(j) = dd * real (j - 1)
      end do

!     Merge the raw times into the uniform time grid:

      ntgrid = ntgridu

      do iTC = 1, numTCs + numsensors
         call merger (numreadings(iTC), time_raw(1,iTC), ntgrid, tgrid, tolt)
         write (6, *) 'iTC,ntgrid,tgrid(ntgrid): ', &
                       iTC,ntgrid,tgrid(ntgrid)
      end do

!     Merge the raw depths into the uniform depth grid:

      ndgrid = ndgridu
      nreadings = 1

      do iTC = 1, numTCs + numsensors
         if (iTC > numTCs) nreadings = numreadings(iTC)

         call merger (nreadings, depth_raw(1,iTC), ndgrid, dgrid, told)
         write (6, *) 'iTC,ndgrid,dgrid(ndgrid): ', &
                       iTC,ndgrid,dgrid(ndgrid)
      end do

      allocate (tmesh(ntgrid,ndgrid), dmesh(ntgrid,ndgrid), &
                fmesh(ntgrid,ndgrid))

      write (luncrt, '(a, 2i6)') 'ntgrid, ndgrid after merging:', ntgrid, ndgrid
      write (11, '(f10.4)') tgrid(1:ntgrid)
      write (12, '(f10.4)') dgrid(1:ndgrid)

      do j = 1, ndgrid
         tmesh(:,j) = tgrid(:)
      end do

      do i = 1, ntgrid
         dmesh(i,:) = dgrid(:)
      end do

!     1-D interpolation of (adjusted) TC readings, in the time direction:

      allocate (tinterp(ntgrid,numTCs), unused(ntgrid))

      do iTC = 1, numTCs
         new = .true.
         write (51, *) 'iTC, ndata, ntgrid: ', iTC, numreadings(iTC), ntgrid
         call lcsfit (numreadings(iTC), time_raw(1,iTC), temp_raw(1,iTC), new, &
                      mono, ntgrid, tgrid, tinterp(1,iTC), unused)
      end do

      deallocate (unused)

!     1-D interpolation of interpolated TC readings, in the depth direction:

      allocate (unused(ndgrid))

      do i = 1, ntgrid
         tmesh(i,:) = tgrid(i)
         dmesh(i,:) = dgrid(:)
         new = .true.
         call lcsfit (numTCs, depth_raw(1,:), tinterp(i,:), new,               &
                      mono, ndgrid, dgrid, fmesh(i,:), unused)
      end do

      deallocate (tinterp, unused)

!!    write (13, '(i1)') 1
      write (13, '(2i6)') ntgrid, ndgrid
      write (13, '(8f10.4)') tmesh, dmesh

!!    write (14, '(i1)') 1
      write (14, '(3i6)') ntgrid, ndgrid, 1
      write (14, '(8f10.4)') fmesh

      end subroutine regularize_TC_data

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   end program MEDLI_analysis
