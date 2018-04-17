
!
! MAIN PROGRAM HERE.  This is just an example so you know
! how to use it. 
! 

module time_kdtree
  use kd_tree
contains

  real function time_search(tree,nsearch,mode, nn, r2)
    !
    !  Return CPU time, in seconds, for searching 'nsearch' reference points
    !  using any specific search mode. 
    !

    type(tree_master_record), pointer :: tree
    integer, intent(in)               :: nsearch ! how many reference points
    integer, intent(in)               :: mode    ! what kind of search
    integer, intent(in)               :: nn      ! number of neighbors 
    real, intent(in)                  :: r2      ! radius^2
    !
    real    :: qv(tree%dim), rv               ! query vector, random variate
    integer :: i, random_loc
    real    :: t0, t1
    real, pointer    :: dists(:)
    integer, pointer :: inds(:)
    real(kind(0.0d0)) :: nftotal

    call cpu_time(t0)   ! get initial time.
    nftotal = 0.0
    allocate(dists(nn))
    allocate(inds(nn) )
    do i=1,nsearch

       select case (mode)
       case (1)
          !
          !  Fixed NN search around randomly chosen point
          !
          call random_number(qv)
          call n_nearest_to(tp=tree,qv=qv,n=nn,indexes=inds,distances=dists)
       case (2)
          !
          ! Fixed NN seasrch around randomly chosen point in dataset
          ! with 100 correlation time
          !
          call random_number(rv)
          random_loc = floor(rv*tree%n) + 1
          call n_nearest_to_around_point(tp=tree,idxin=random_loc,&
           correltime=100,n=nn,indexes=inds,distances=dists)
       case default
          write (*,*) 'Search type ', mode, ' not implemented.'
          time_search = -1.0 ! invalid
          return
       end select

    enddo
    call cpu_time(t1)

    time_search = t1-t0
    if (nftotal .gt. 0.0) then
!       write (*,*) 'Average number of neighbors found = ', nftotal / real(nsearch)
    endif
    deallocate(inds,dists)
    return
  end function time_search


  real function searches_per_second(tree,mode,nn,r2) result(res)
    !
    !
    ! return estimated number of searches per second.
    ! Will call "time_search" with increasing numbers of reference points 
    ! until CPU time taken is at least 1 second.
    !

    type(tree_master_record), pointer :: tree
    integer, intent(in)               :: mode    ! what kind of search
    integer, intent(in)               :: nn      ! number of neighbors 
    real, intent(in)                  :: r2      ! radius^2
    !
    integer :: nsearch
    real    :: time_taken

    nsearch = 50  ! start with 50 reference points
    do
       time_taken = time_search(tree,nsearch,mode,nn,r2)
       if (time_taken .lt. 1.0) then
          nsearch = nsearch * 5
          cycle
       else
          res = real(nsearch) / time_taken
          return
       end if
    end do
    return
  end function searches_per_second

end module time_kdtree

Program kd_tree_test
  Use kd_tree  ! this is how you get access to the k-d tree routines
  use time_kdtree

  Integer :: n, d 
  Real, Dimension(:,:), Allocatable :: my_array

  Type(tree_master_record), Pointer :: tree
  ! this is how you declare a tree in your main program

  Integer :: k, j

  Real, Allocatable :: query_vec(:)
  Real,allocatable     :: dists(:), distsb(:) 
  Integer,allocatable   :: inds(:), indsb(:) 
  integer :: nnbrute
  real      :: t0, t1, sps

  integer, parameter  :: nnn = 5
  integer,parameter   :: nnarray(nnn) =(/ 1, 5, 10, 25, 500 /) 
  integer, parameter  :: nr2 = 5

  real      :: r2array(nr2)
  parameter (r2array = (/  1.0e-4,1.0e-3,5.0e-3,1.0e-2,2.0e-2 /) )

  Print *, "Type in N and d"
  Read *, n, d

!  n = 250000
!  d = 3

  Allocate(my_array(n,d))
  allocate(query_vec(d))
  Call Random_number(my_array)  !fills entire array with built-in-randoms

  call cpu_time(t0)
  tree => create_tree(my_array) ! this is how you create a tree. 
  call cpu_time(t1)
  write (*,*) real(n)/real(t1-t0), ' points per second built for regular tree.'

  nnbrute = 50
  allocate(dists(nnbrute),distsb(nnbrute))
  allocate(inds(nnbrute),indsb(nnbrute))

  Do k=1,50
     Call Random_number(query_vec)

     ! find five nearest neighbors to.
     inds = -666
     indsb = -777
     Call n_nearest_to(tp=tree, qv=query_vec, n=nnbrute, indexes=inds, distances=dists)

     call n_nearest_to_brute_force(tree, query_vec, nnbrute, indsb, distsb)
     if (ANY(indsb(1:nnbrute) .ne. inds(1:nnbrute)) .or. &
       any(dists(1:nnbrute) .ne. distsb(1:nnbrute))) then
        write (*,*) 'MISMATCH! @ k=',k
        
        do j=1,nnbrute
           write (*,*) j,'Tree/brute ids=',inds(j),indsb(j)
           write (*,*) j,'tree/brute dists=', dists(j), distsb(j)
!           write (*,*) j,'my dist to tree = ',  &
!            square_distance(d,my_array(:,inds(j)),query_vec)
!           write (*,*) j,'my dist to brute = ',  &
!            square_distance(d,my_array(:,indsb(j)),query_vec)
           
        enddo

        print *, "Tree indexes = ", inds(1:nnbrute)

        print *, "Brute indexes = ", indsb(1:nnbrute)
        print *, "Tree distances = ", dists(1:nnbrute)
        print *, "Brute distances = ", distsb(1:nnbrute)
     endif

  Enddo
  

10 format('R^2 search, r2/d=',G10.2,':',F13.0,A)


20 format(A,' NN=',I7,':',F10.0,' searches/s in ',A)

  do k=1,nnn
     sps =searches_per_second(tree,1,nnarray(k),1.0)
     write (*,20) 'Random pts',nnarray(k),sps,'old regular tree.'
  enddo

  do k=1,nnn
     sps = searches_per_second(tree,2,nnarray(k),1.0)
     write (*,20) 'in-data pts',nnarray(k),sps,'old tree.'
  enddo

  Call destroy_tree(tree)  
  ! this releases memory for the tree BUT NOT THE ARRAY OF DATA YOU PASSED
  ! TO MAKE THE TREE.  

  Deallocate(my_array)
  ! deallocate the memory for the data.

End Program kd_tree_test
