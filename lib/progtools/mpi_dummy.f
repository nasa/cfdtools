c***********************************
c   A dummy declaration of the MPI routines
c   such that some libraries can be used in 
c   a non MPI environment
c
c***********************************

c+-----------------------------------------------------------------------

      subroutine mpi_waitall (nnrecv, ireq, istat, ierr)
c
c   Description:  Dummy MPI routine to facilitate execution in
c   non-MPI environments
c
c------------------------------------------------------------------------
      integer ::
     .  nrecv,
     .  ireq,
     .  istat,
     .  ierr

       ierr = 0

      end subroutine mpi_waitall

c+-----------------------------------------------------------------------

      subroutine mpi_finalize(ierr)
c
c   Description:  Dummy MPI routine to facilitate execution in
c   non-MPI environments
c
c------------------------------------------------------------------------

      integer :: ierr
      ierr = 0
      end subroutine mpi_finalize

c+-----------------------------------------------------------------------

      subroutine mpi_irecv (r_rbuff, lbuf, MPI_INTEGER, i_from,
     .  i_rctag, MPI_COMM_WORLD, i_request_r, ierr)
c
c   Description:  Dummy MPI routine to facilitate execution in
c   non-MPI environments
c
c------------------------------------------------------------------------
      real ::
     .  r_rbuff(*)

      integer ::
     .  lbuf,
     .  MPI_INTEGER,
     .  i_from,
     .  i_rctag,
     .  MPI_COMM_WORLD,
     .  i_request_r,
     .  ierr

      ierr = 0

      end subroutine mpi_irecv

c+-----------------------------------------------------------------------

      subroutine mpi_isend (r_sbuff, lbuf, MPI_DOUBLE_PRECISION, i_dest,
     .                i_sntag, MPI_COMM_WORLD, 
     .                i_request_s, ierr     )
c
c   Description:  Dummy MPI routine to facilitate execution in
c   non-MPI environments
c
c------------------------------------------------------------------------
      real ::
     . r_sbuff(*)

      integer ::
     .  lbuf,
     .  MPI_DOUBLE_PRECISION,
     .  i_dest,
     .  i_sntag,
     .  MPI_COMM_WORLD,
     .  i_request_s,
     .  ierr

      ierr = 0

      end subroutine mpi_isend

c+-----------------------------------------------------------------------

      subroutine mpi_init (ierr)
c
c   Description:  Dummy MPI routine to facilitate execution in
c   non-MPI environments
c
c------------------------------------------------------------------------

      integer :: ierr

      ierr = 0

      end subroutine mpi_init

c+-----------------------------------------------------------------------

      subroutine mpi_comm_size (MPI_COMM_WORLD, n_proc, ierr)
c
c   Description:  Dummy MPI routine to facilitate execution in
c   non-MPI environments
c
c------------------------------------------------------------------------

      integer ::
     .  MPI_COMM_WORLD,
     .  n_proc,
     .  ierr

      ierr = 0
      n_proc = 1 
     
      end subroutine mpi_comm_size

c+-----------------------------------------------------------------------

      subroutine mpi_comm_rank (MPI_COMM_WORLD, i_proc, ierr)
c
c   Description:  Dummy MPI routine to facilitate execution in
c   non-MPI environments
c
c------------------------------------------------------------------------

      integer ::
     .  MPI_COMM_WORLD
     .  i_proc,
     .  ierr

      ierr = 0
      i_proc = 0
      end subroutine mpi_comm_rank

c+-----------------------------------------------------------------------

      subroutine mpi_wait (i, istat, ierr)
c
c   Description:  Dummy MPI routine to facilitate execution in
c   non-MPI environments
c
c------------------------------------------------------------------------
      integer ::
     .  i, istat, ierr

      ierr = 0
      end subroutine mpi_wait
     
c+-----------------------------------------------------------------------

      function mpi_wtime () result (time)
c
c   Description:  Dummy MPI routine to facilitate execution in
c   non-MPI environments
c
c------------------------------------------------------------------------

      real :: time

      call second (time)

      end function mpi_wtime

c+-----------------------------------------------------------------------
c
      subroutine mpi_barrier (MPI_COMM_WORLD, ierr)
c
c   Description:  Dummy MPI routine to facilitate execution in
c   non-MPI environments
c
c------------------------------------------------------------------------

      integer ::
     .  MPI_COMM_WORLD,
     .  ierr 

      ierr = 0
      end subroutine mpi_barrier

