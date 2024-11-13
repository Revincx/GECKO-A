MODULE primetool
  IMPLICIT NONE
  CONTAINS

!=================================================================
! SUBROUTINE primes
! PURPOSE: update rank of nodes in chem. formula by finding and
! ranking the product of adjacent primes of input rank
!=================================================================
SUBROUTINE primes(nca,bond,rank)
  USE keyparameter, ONLY: prim
  IMPLICIT NONE

  INTEGER,INTENT(IN) :: nca
  INTEGER,INTENT(IN) :: bond(:,:)
  INTEGER,INTENT(INOUT) :: rank(:) ! rank of nodes, in order supplied

! internal:
  INTEGER ::  i,j,maxr
  INTEGER ::  cprim(SIZE(rank))   ! prime corresponding to node's rank
  INTEGER ::  pp(SIZE(rank))      ! product of primes of adjacent nodes
  INTEGER ::  prank(SIZE(rank),3) ! rank pairs (node index, input rank, pp)


! re-entry point if iteration required
  iterloop : DO

! assign primes
    DO i=1,nca
      cprim(i)=prim(rank(i))
      pp(i)=1
    ENDDO
   
! find product of connected primes for nodes
    DO i=1,nca
      DO j=1,nca
        IF(i/=j.AND.bond(i,j)>0) pp(i)=pp(i)*cprim(j)
      ENDDO
    ENDDO
   
! construct rank vectors for nodes
    maxr = 0
    DO i=1,nca
      prank(i,1) = i
      prank(i,2) = rank(i)
      prank(i,3) = pp(i)
      maxr = MAX(maxr,rank(i))
      rank(i)=0
    ENDDO
   
    CALL sortrank(nca,prank,rank)
   
! if rank changed: need to iterate again
    DO i=1,nca
      IF (rank(i) /= prank(i,2)) CYCLE iterloop
    ENDDO
    EXIT iterloop ! done!
  ENDDO iterloop
END SUBROUTINE primes

!=============================================================
! SUBROUTINE sortrank
! rank numerical priorities
!=============================================================
SUBROUTINE sortrank(nca,prank,rank)
  IMPLICIT NONE

  INTEGER,INTENT(IN)  :: nca
  INTEGER,INTENT(IN)  :: prank(:,:)
  INTEGER,INTENT(OUT) :: rank(:)

! internal
  INTEGER pr(SIZE(prank,1),SIZE(prank,2))
  INTEGER store(3)
  INTEGER i,ii,j,ctr

! just to be sure
  IF (SIZE(prank,2)/=3) THEN
    PRINT*, "in sortrank, size array issue"
    STOP "in sortrank"
  ENDIF
  
! truncate input rank array 
  pr(:,:)=prank(:,:)

! first pass
  loop1 : DO
    i=1
    loop2: DO
    ii=i+1
      IF (pr(i,2) <= pr(ii,2)) THEN
        IF (ii==nca) EXIT loop1
        i = i + 1
        CYCLE loop2
      ENDIF
      DO j=1,3
        store(j) = pr(ii,j)  
        pr(ii,j)  = pr(i,j)
        pr(i,j)  = store(j)
      ENDDO
      i = i-1
      IF (i==0) CYCLE loop1
      CYCLE loop2
    ENDDO loop2
  ENDDO loop1

! second pass
  loop1b: DO
    i=1
    loop2b: DO
      ii=i+1
      IF (pr(i,2)==pr(ii,2))THEN
        IF(pr(i,3)<=pr(ii,3)) THEN
          IF (ii==nca) EXIT loop1b
          i = i + 1
          CYCLE loop2b
        ENDIF
        DO j=1,3
          store(j) = pr(ii,j)  
          pr(ii,j)  = pr(i,j)
          pr(i,j)  = store(j)
        ENDDO
        i = i-1
      ELSE
          IF (ii==nca) EXIT loop1b
          i = i + 1
          CYCLE loop2b
      ENDIF
      IF (i==0) CYCLE loop1b
      CYCLE loop2b
    ENDDO loop2b
  ENDDO loop1b

! find new ranks
  ctr=1
  rank(pr(1,1))=ctr
  DO ii=2,nca
    i=ii-1
    IF(pr(ii,2)==pr(i,2).AND.pr(ii,3)==pr(i,3)) THEN
      ctr=ctr
    ELSE
      ctr=ctr+1
    ENDIF  
    rank(pr(ii,1))=ctr
  ENDDO

END SUBROUTINE sortrank

END MODULE primetool

