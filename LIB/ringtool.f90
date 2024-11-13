MODULE ringtool
  IMPLICIT NONE
  CONTAINS

! SUBROUTINE uniqring(nring,nca,group,bond,rank,rjg)
! SUBROUTINE findring(i1,i2,nca,bond,rngflg,ring)
! SUBROUTINE findtree(con,top,sec,nca,parent,child)
! SUBROUTINE rejoin(nca,x,y,x1,y1,bond,group)

!=======================================================================
! PURPOSE: Find unique artificial breakpoints for cyclic molecule 
!          including multi-ring compounds.                        
!=======================================================================
SUBROUTINE uniqring(nring,nca,group,bond,rank,rjg)
  USE keyparameter, ONLY: prim
  USE rjtool
  USE stdratings
  IMPLICIT NONE

  INTEGER,INTENT(IN)   :: nring ! # of rings
  INTEGER,INTENT(IN)   :: nca   ! # of groups
  CHARACTER(LEN=*),INTENT(INOUT) :: group(:) ! group matrix 
  INTEGER,INTENT(INOUT):: bond(:,:)          ! bond matrix
  INTEGER,INTENT(INOUT):: rank(:)  ! rank of nodes, in order supplied
  INTEGER,INTENT(OUT)  :: rjg(:,:) ! ring-join group pairs

! internal:
  INTEGER :: ppb(SIZE(bond,1),SIZE(bond,2))   ! product of primes of nodes in bond
  INTEGER :: cprim(SIZE(group))              ! prime corresponding to node's rank
  INTEGER :: tbond(SIZE(bond,1),SIZE(bond,2)) ! working bond matrix.
  INTEGER :: i,j,ii,jj,k,maxp,n,rngflg
  INTEGER :: top,sec                   ! first and second nodes in new arrangement
  INTEGER :: parent(SIZE(group))
  INTEGER :: child(SIZE(group),3)
  INTEGER :: nk(SIZE(group))
  INTEGER :: iold(SIZE(group))           ! old index of group whose new index is 'k'
  INTEGER :: trjg(nring,2)                ! temporary ring-join group pairs
  INTEGER :: trank(SIZE(group))          ! temporary rank listing, for rearrangement
  CHARACTER(LEN=LEN(group(1))) :: tgroup(SIZE(group)) ! temporary group listing, for rearr.
  INTEGER :: rgpath(SIZE(bond,1))

! initialise
  rjg(:,:)=0  ;  tbond(:,:)=0  ;  tgroup(:)=' ' 

! find product of connected primes for actual bonds, setup connections 
  cprim(1:nca)=prim(rank(1:nca))
  ppb(:,:) = 0  ;  tbond(:,:)=bond(:,:)
  DO i=1,nca
    DO j=1,nca
      IF (i/=j .AND. bond(i,j) > 0) ppb(i,j)=cprim(i)*cprim(j)
    ENDDO
  ENDDO

! remove ring-joining characters from groups
  CALL rjgrm(nring,group,rjg)
      
! loop rings, finding lowest-priority on-ring connection to 'break' 
  DO n=1,nring
    grloop: DO
      maxp = 0
      DO i=1,nca
        DO j=i+1,nca
          IF(bond(i,j)/=2)THEN ! disallow breaking @ double bonds
            IF(ppb(i,j) > maxp) THEN
              ii = i  ;  jj = j
              maxp=ppb(ii,jj) 
            ENDIF
          ENDIF
        ENDDO
      ENDDO
      CALL findring(ii,jj,nca,tbond,rngflg,rgpath)
      ppb(ii,jj) = 0
      IF (rngflg==1) THEN ! connection is on a ring
        tbond(ii,jj) = 0  ;  tbond(jj,ii) = 0
        rjg(n,1) = ii     ;  rjg(n,2) = jj
! symmetry might be broken after first bond breaking.
! therefore, call ratings to update ranks.          
        tgroup(:)=group(:)
        IF (nring > 1) CALL ratings(nca,tgroup,tbond,nring,rank)
        tgroup(:)=' '
      ELSE
        CYCLE grloop
      ENDIF
      EXIT grloop
    ENDDO  grloop
  ENDDO ! n_rings

! find first 2 nodes in tree (linearly-expressed molecule)
  top = ii  ;  sec = 0
  nodloop: DO j = 2,nca
    DO i = 1,j-1
      IF(tbond(i,j) > 0)THEN
        IF (i==ii) sec=j
        IF (j==ii) sec=i
        IF (sec > 0) EXIT nodloop
      ENDIF
    ENDDO
  ENDDO nodloop

! rearrange into a connected skeleton (finding longest path comes later)
  CALL findtree(tbond,top,sec,nca,parent,child)

  nk(:)=0  ;  j=1  ;  iold(j)=top
  IF (parent(iold(j))/=0) THEN
    PRINT*, "error 1 in uniqring"
    STOP "in uniqring"
  ENDIF

  i=j
  xxloop: DO
    k = nk(i)+1
    IF (k <= 3) THEN
      IF (child(iold(i),k) > 0) THEN
        nk(i) = nk(i)+1
        j=j+1
        iold(j)=child(iold(i),k)
        i=j ! step forward
        CYCLE xxloop ! try for next child
      ENDIF
    ENDIF

! otherwise, no further children on node
    IF (parent(iold(i))==0) EXIT xxloop ! finished thoroughly checking molecule
    DO ii=1,nca
      IF (parent(iold(i))==iold(ii))THEN
        i = ii
        CYCLE xxloop ! backtrack to parent
      ENDIF
    ENDDO
    EXIT xxloop
  ENDDO xxloop


! reassign arrays based on new ordering
  DO ii=1,nca
    IF(iold(ii)/=0)THEN
      trank(ii) = rank(iold(ii))
      tgroup(ii) = group(iold(ii))
    ENDIF
    DO jj=1,nca
      IF(iold(jj)/=0)THEN
        tbond(ii,jj) = bond(iold(ii),iold(jj))
      ELSE
! retain zeroes for broken bonds (e.g. for fragmentation that retains ring)
        tbond(ii,jj)=0
      ENDIF
    ENDDO
    DO n=1,nring
      DO i=1,2
        IF(rjg(n,i)==iold(ii)) THEN
          trjg(n,i)=ii
        ENDIF
      ENDDO
    ENDDO
  ENDDO

! b: replace working arrays
  DO ii=1,nca
    rank(ii) = trank(ii)
    group(ii) = tgroup(ii)
    DO jj=1,nca
      bond(ii,jj) = tbond(ii,jj)
    ENDDO
  ENDDO
  DO n=1,nring
    rjg(n,1)=MIN(trjg(n,1),trjg(n,2))
    rjg(n,2)=MAX(trjg(n,1),trjg(n,2))
  ENDDO
 
  CALL rjgadd(nring,group,rjg) ! add ring-joining characters

END SUBROUTINE uniqring

!=======================================================================
! PURPOSE: ascertain whether a given bond is part of any ring
!=======================================================================
SUBROUTINE findring(i1,i2,nca,bond,rngflg,ring)
  USE keyparameter, ONLY:mxcp  ! used to set max # of tracks
  USE mapping
  IMPLICIT NONE

  INTEGER,INTENT(IN) :: i1    ! 1st node
  INTEGER,INTENT(IN) :: i2    ! 2nd node
  INTEGER,INTENT(IN) :: nca   ! # of nodes
  INTEGER,INTENT(IN) :: bond(:,:) ! bond matrix
  INTEGER,INTENT(OUT) :: rngflg       ! 0 = 'no', 1 = 'yes'
  INTEGER,INTENT(OUT) :: ring(:)    ! =1 if node participates in current ring

  INTEGER             :: i,j,m,n        ! current, next, previous nodes
  INTEGER             :: track(mxcp,SIZE(bond,1))
  INTEGER             :: trlen(mxcp)
  INTEGER             :: ntr

  rngflg=0  ;  ring(:)=0

  DO m=1,nca-1
    DO n=m+1,nca
      IF (bond(m,n)/=0) THEN
! find tracks starting at m
        CALL gettrack(bond,m,nca,ntr,track,trlen)
! search in the tracks if node n is found (at least in beta)
        DO i=1,ntr
          DO j=3,trlen(i)
            IF (track(i,j)==n) THEN
              IF (((m==i1).AND.(n==i2)) .OR. ((m==i2).AND.(n==i1))) THEN
                rngflg=1
              ENDIF
              ring(track(i,1:j))=1
            ENDIF
          ENDDO
        ENDDO
      ENDIF
    ENDDO
  ENDDO
     
END SUBROUTINE findring

!=======================================================================
! PURPOSE : Set up the tree of the C-C and C-O-C bond starting at     
!           node (top). Adapted from subroutine lntree.               
!=======================================================================
SUBROUTINE findtree(con,top,sec,nca,parent,child)
  IMPLICIT NONE

  INTEGER,INTENT(IN) :: nca   ! # of groups
  INTEGER,INTENT(IN) :: top   ! 1st group => starting point 
  INTEGER,INTENT(IN) :: sec   ! 2nd group => node just after top
  INTEGER,INTENT(IN) :: con(:,:)   ! bond matrix (after artificial ring-break)
  INTEGER,INTENT(OUT):: parent(:)  ! node one step back from node(i)
  INTEGER,INTENT(OUT):: child(:,:) ! all nodes one step forward from node(i) 

  INTEGER :: left(SIZE(con,1)),right(SIZE(con,1)),center(SIZE(con,1))
  INTEGER :: ptr,knt,nknt
  INTEGER :: i,j,k,iend
  LOGICAL :: tcon(SIZE(con,1),SIZE(con,1))

! -----------
! initialize
! -----------
  ptr=0 ; left(:)=0 ; right(:)=0 ; center(:)=0 ; parent(:)=0

! make a logical bond matrix (tcon only used to find parent and child)
  tcon(:,:)=.FALSE.
  WHERE (con(:,:)/=0) tcon(:,:)=.TRUE.
  
  left(top)=sec            ;  parent(sec)=top
  tcon(top,sec) = .FALSE.  ;  tcon(sec,top) = .FALSE.

! ---------------------------------------
! get the relationships (parent/children)
! ---------------------------------------
  grloop: DO k=1,nca
    nknt = 0
    DO i=1,nca
      knt = 0
      DO j=1,nca
        IF (tcon(i,j)) THEN
          knt = knt + 1
          ptr = j
        ENDIF
      ENDDO
      nknt = nknt + knt

! look for ith carbon with only one node, where parents do not
! exist: (ith carbon = child,jth carbon = parent)
      IF (knt==1) THEN
        IF ( (parent(i)==0) .AND. (i/=top) ) THEN
          parent(i) = ptr
          tcon(i,ptr) = .FALSE.  ;  tcon(ptr,i) = .FALSE.
          IF (left(ptr)==0)       THEN ; left(ptr) = i
          ELSE IF (right(ptr)==0) THEN ; right(ptr) = i
          ELSE IF (center(ptr)==0)THEN ; center(ptr) = i
          ELSE   ! all children taken, error in bond-matrix
            DO j=1,nca
              WRITE(6,*) (con(j,iend),iend=1,nca)
            ENDDO
            WRITE(6,'(a)') '--error-- in findtree. No possible path'
            STOP "in findtree"
          ENDIF
        ENDIF
      ENDIF
    ENDDO
! do loop until bond-matrix is (0).
    IF (nknt==0) EXIT grloop
  ENDDO grloop

! define 'children'      
  DO i=1,nca
    child(i,1)=left(i)
    child(i,2)=right(i)
    child(i,3)=center(i)
  ENDDO

END SUBROUTINE findtree

! ======================================================================
! Purpose: takes the two ends of the linear string formed after
! ring-opening, and joins them at the artificial 'break' in the old ring.
! Also good for rearranging rings after elimination of a branch.
! ======================================================================
SUBROUTINE rejoin(nca,x,y,x1,y1,bond,group)
  IMPLICIT NONE

  INTEGER,INTENT(IN) ::  nca
  INTEGER,INTENT(IN) ::  x
  INTEGER,INTENT(IN) ::  y
  INTEGER,INTENT(INOUT)          :: bond(:,:)
  CHARACTER(LEN=*),INTENT(INOUT) :: group(:)
  INTEGER,INTENT(OUT) :: x1 ! new identities of x,y
  INTEGER,INTENT(OUT) :: y1 ! new identities of x,y

  INTEGER   :: i,j,ii,jj,icut
  INTEGER   :: tbond(SIZE(bond,1),SIZE(bond,2))
  CHARACTER(LEN=LEN(group(1))) :: tgroup(SIZE(group))

  tbond(:,:)=0  ;  tgroup(:)=' '
  
  IF (x>y) THEN
    icut=x-1
  ELSE
    icut=y-1
  ENDIF
  x1=x-icut  ;  IF (x1<=0) x1=x1+nca
  y1=y-icut  ;  IF (y1<=0) y1=y1+nca

! switch positions of two ends of molecule: stage 1      
  DO i=1,icut
    ii=nca-icut+i
    tgroup(ii)=group(i)
    DO j=1,nca
      tbond(ii,j)=bond(i,j)
    ENDDO
  ENDDO
  DO i=icut+1,nca
    ii=i-icut
    tgroup(ii)=group(i)
    DO j=1,nca
      tbond(ii,j)=bond(i,j)
    ENDDO
  ENDDO

! replace bond matrix with interim solution
  bond(:,:)=tbond(:,:)

! switch positions of two ends of molecule: stage 2      
  DO i=1,nca
    DO j=1,icut
      jj=nca-icut+j
      tbond(i,jj)=bond(i,j)
    ENDDO
    DO j=icut+1,nca
      jj=j-icut
      tbond(i,jj)=bond(i,j)
    ENDDO
  ENDDO

! replace bond matrix and groups with final solution
  group(:)=tgroup(:)  ;  bond(:,:)=tbond(:,:)
END SUBROUTINE rejoin

!=======================================================================
! Purpose: Find all rings including node ig.
!=======================================================================
SUBROUTINE ring_data(ig,nca,tbond,tgroup,ndrg,lodrg,trackrg)
  USE keyparameter, ONLY: mxcp
  USE mapping, ONLY: gettrack
  IMPLICIT NONE

  INTEGER,INTENT(IN)         :: tbond(:,:)
  CHARACTER(LEN=*),INTENT(IN):: tgroup(:)
  INTEGER,INTENT(IN)  :: ig
  INTEGER,INTENT(IN)  :: nca
  LOGICAL,INTENT(OUT) :: lodrg(:,:)   ! (i,j)==true if node j belong to ring i
  INTEGER,INTENT(OUT) :: ndrg         ! # of distinct rings 
  INTEGER,INTENT(OUT) :: trackrg(:,:) ! (i,:)== track (node #) belonging ring i

  INTEGER :: i,j,k,l
  INTEGER :: rngflg,ring(SIZE(tgroup))
  INTEGER :: track(mxcp,SIZE(tgroup))
  INTEGER :: trlen(mxcp)
  INTEGER :: ntr
  LOGICAL :: ring_tmp(SIZE(tgroup))

  ndrg=0                ! number of distinct rings
  lodrg(:,:)=.FALSE.    ! logical matrix of all distincts rings (ring #, node TRUE if belongs to ring)
  trackrg(:,:)=0        ! track of each ring (number of the ring, track of the ring)

  DO i=1,nca
    IF (tbond(ig,i)/=0) THEN          
      CALL findring(ig,i,nca,tbond,rngflg,ring)
      CALL gettrack(tbond,ig,nca,ntr,track,trlen)
      
      IF (rngflg==1) THEN
        trloop: DO j=1,ntr
          ring_tmp(:)=.FALSE.
          DO k=3,trlen(j)
            IF (track(j,k)==i) THEN
              IF (ndrg==0) THEN
! if current ring is the first ring
                ndrg=1
                lodrg(1,track(j,1:k))=.TRUE.
                trackrg(1,1:k)=track(j,1:k)
              ELSE
! test if the current ring already exist, if yes cycle
                ring_tmp(track(j,1:k))=.TRUE.
                DO l=1,ndrg
                  IF (ALL(lodrg(l,:).EQV.ring_tmp(:))) CYCLE trloop
                ENDDO
! if not, the new ring is added in the tables
                ndrg=ndrg+1 
                IF (ndrg>SIZE(lodrg,1)) THEN
                  PRINT*, "in ring_data, ndrg > SIZE(lodrg,1)"
                  STOP "in ring_data"
                ENDIF 
                trackrg(ndrg,1:k)=track(j,1:k)
                lodrg(ndrg,track(j,1:k))=.TRUE. 
              ENDIF                    
            ENDIF
          ENDDO
        ENDDO trloop
      ENDIF
    ENDIF
  ENDDO

END SUBROUTINE ring_data

END MODULE ringtool
