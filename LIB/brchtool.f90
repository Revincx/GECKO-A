MODULE brchtool
  IMPLICIT NONE
  CONTAINS

! SUBROUTINE getbrch(lobond,group,nca,rank,nring,ig,ia,ib,brch,ml)
! SUBROUTINE mkbrcopy(lobond,group,nca,rank,ig,pg,ng,ptr1,brcopy)
! SUBROUTINE brpri(ogroup,rank,copy,ncp,nring,maxpri)
! SUBROUTINE treebr(lobond,ig,top,sec,nca,brpath,maxlng,npath)

! ======================================================================
! PURPOSE:  write the branch starting at the position ia (alpha       
! position) and bounded to the position ig. One of the beta positions 
! is known (ib). The various ways of writing the branch are checked 
!  => A standard "formula" of the chain is given as output                                        
!                                                                     
! The structure of the subroutine is :                                
! 1- find longest trees starting from ia                              
! 2- write for each tree a copy of the corresponding formula          
! 3- find which formula has the highest priority                      
! ======================================================================
SUBROUTINE getbrch(olobond,group,nca,rank,nring,ig,ia,ib,brch,ml)
  IMPLICIT NONE

  LOGICAL, INTENT(IN) :: olobond(:,:) ! logical bond matrix 
  CHARACTER(LEN=*),INTENT(IN) :: group(:) ! group table
  INTEGER,INTENT(IN) :: nca   ! total number of group in the molecule
  INTEGER,INTENT(IN) :: nring ! # of ring
  INTEGER,INTENT(IN) :: ig    ! group # to which the branch is bounded
  INTEGER,INTENT(IN) :: ia    ! group # where the branch begin (alpha from ig) 
  INTEGER,INTENT(IN) :: ib    ! group bounded to ia (beta from ig)
  INTEGER,INTENT(IN) :: rank(:) ! rank (priority) of group i (low = more important)
  CHARACTER(LEN=*),INTENT(OUT) :: brch ! "full" standardized formula of the branch 
  INTEGER, INTENT(OUT) :: ml  ! length of the branch

! internal:
  INTEGER :: np,tnp,maxpri,ptra
  INTEGER :: iadd,pre,nex
  INTEGER :: brpath(SIZE(olobond,1),SIZE(olobond,2))
  CHARACTER(LEN=LEN(brch)) :: brcopy(SIZE(group))
  CHARACTER(LEN=LEN(brch)) :: tbrcopy(SIZE(group))
  INTEGER :: i,j
  LOGICAL :: lobond(SIZE(olobond,1),SIZE(olobond,1))

! initialize
! ----------
  brch=' '  ;   ml=0  ;  lobond(:,:)=olobond(:,:)

! search the longest tree, starting from ia
! -----------------------------------------
  lobond(ia,ig)=.FALSE.  ;  lobond(ig,ia)=.FALSE.
  CALL treebr(lobond,ig,ia,ib,nca,brpath,ml,np)
  lobond(ia,ig)=.TRUE.   ;  lobond(ig,ia)=.TRUE.

! write copy of longest branches
! ------------------------------
  DO i=1,np
    ptra=1  ;  brcopy(i)=' '
    DO j=1,ml
      iadd=brpath(i,j)
      IF (j==1) THEN ; pre = ig ; ELSE ; pre = brpath(i,j-1) ; ENDIF
      nex = brpath(i,j+1)
      CALL mkbrcopy(lobond,group,nca,rank,iadd,pre,nex,ptra,brcopy(i))
    ENDDO
  ENDDO

! If only 1 chain : write out and return
! ---------------------------------------
  IF (np==1) THEN
    brch=brcopy(1) ; RETURN
  ENDIF

! if more than 1 chain : get the chain with the highest priority
! --------------------------------------------------------------

! collapse identical copy
  DO i=1,np-1
    DO j=i+1,np
      IF (brcopy(i)==brcopy(j))  brcopy(j)=' '
    ENDDO
  ENDDO

  tnp=np  ;  np=0
  DO i=1,tnp
    IF (brcopy(i)(1:1)/=' ') THEN
      np=np+1
      tbrcopy(np)=brcopy(i)
    ENDIF
  ENDDO

! if only 1 chain remain then write out and return else search priority
  IF (np==1) THEN
    brch=tbrcopy(1)  ;   RETURN
  ELSE
    CALL brpri(group,rank,tbrcopy,np,nring,maxpri)
    brch=tbrcopy(maxpri) ;  RETURN
  ENDIF

END SUBROUTINE getbrch

! ======================================================================
! PURPOSE: THIS ROUTINE IS VERY SIMILAR TO MKCOPY. Main difference is 
! that mkcopy is called to write copies of the "full" molecule while 
! mkbrcopy is called to write copies of a given branch only.                           
! mkbrcopy makes a copy of a branch according to the longest tree. On 
! each call of mkcopy, the group ig is written in "copy". At the same 
! time and according to the bond-matrix, all groups which have a bond 
! to ig (except the bond of longest tree), are attached to it (in 
! parentheses). These attached must not have more than 1C 
! ======================================================================
SUBROUTINE mkbrcopy(lobond,group,nca,rank,ig,pg,ng,ptr1,brcopy)
  IMPLICIT NONE

  LOGICAL,INTENT(IN) :: lobond(:,:) ! logical bond matrix
  CHARACTER(LEN=*),INTENT(IN) :: group(:) ! group table
  INTEGER,INTENT(IN) :: rank(:) ! rank (priority) of group i
  INTEGER,INTENT(IN) :: nca  ! total # of group in the molecule
  INTEGER,INTENT(IN) :: ig   ! group # to which the branch is bounded
  INTEGER,INTENT(IN) :: pg   ! previous group of ig in the longest tree
  INTEGER,INTENT(IN) :: ng   ! next group of ig in the longest tree
  CHARACTER(LEN=*),INTENT(INOUT) :: brcopy  ! copy of branch
  INTEGER,INTENT(INOUT) :: ptr1  ! pointer in "copy" to add next group

  INTEGER :: ptr2,i,j,ita,ia1,ia2

! ------------------------------------------------------
! WRITE GROUP IG TO BRCOPY
! ------------------------------------------------------
  ptr2 = ptr1 + INDEX(group(ig),' ') - 2
  brcopy(ptr1:ptr2) = group(ig)
  ptr1 = ptr2 + 1       

! ------------------------------------------------------
! LOOK FOR ALL ATTACHED GROUPS TO GROUP IG, NOT PG OR NG
! ------------------------------------------------------
! search and write the attached group to IG in "alpha" position
  ita=0  ;  ia1=0  ;  ia2=0
  DO i=1,nca
    IF (lobond(ig,i)) THEN
      IF ( (i/=pg) .AND. (i/=ng) ) THEN
        ita=ita+1
        IF (ita==1) ia1=i
        IF (ita==2) ia2=i

! search the attached group to IG in "beta" position. If any then
! error (or change the program) - this position can only be occupied
! for molecule having at least 12 groups (ie C12)
        DO j=1,nca
          IF (lobond(i,j)) THEN
            IF (j/=ig) THEN
              WRITE (6,*) '--error--,in mkbrcopy. The branch has a group in beta '
              WRITE (6,*) 'position of the closest C in the main (longest) branch'
              STOP "in mkbrcopy"
            ENDIF
          ENDIF
        ENDDO
      ENDIF
    ENDIF
  ENDDO

! ------------------------------------------------------
! CHECK THE VARIOUS POSSIBLE CASES
! ------------------------------------------------------
! 3 cases must be considered
! - case 1 : no alpha group => return
! - case 2 : only 1 alpha group => write it and return
! - case 3 : 2 alpha group => evaluate priority and write in sorted way

! CASE 1 (no alpha group)
! ------
  IF (ita==0) RETURN

! CASE 2 (1 alpha group)
! ------
  IF (ita==1) THEN
    ptr2 = INDEX(group(ia1),' ') - 1
    ptr2 = ptr1 + ptr2
    brcopy(ptr1:ptr1)   = '('
    brcopy(ptr1+1:ptr2) = group(ia1)
    ptr1 = ptr2 + 1
    brcopy(ptr1:ptr1) = ')'
    ptr1 = ptr1 + 1
    RETURN
  ENDIF

! CASE 3 (2 alpha groups)
! ------
! write the group having the lowest priority first
  IF (ita==2) THEN
    IF (rank(ia2) < rank(ia1)) THEN
      ptr2 = INDEX(group(ia1),' ') - 1
      ptr2 = ptr1 + ptr2
      brcopy(ptr1:ptr1)   = '('
      brcopy(ptr1+1:ptr2) = group(ia1)
      ptr1 = ptr2 + 1
      brcopy(ptr1:ptr1) = ')'
      ptr1 = ptr1 + 1
  
      ptr2 = INDEX(group(ia2),' ') - 1
      ptr2 = ptr1 + ptr2
      brcopy(ptr1:ptr1)   = '('
      brcopy(ptr1+1:ptr2) = group(ia2)
      ptr1 = ptr2 + 1
      brcopy(ptr1:ptr1) = ')'
      ptr1 = ptr1 + 1
      RETURN
   
    ELSE
      ptr2 = INDEX(group(ia2),' ') - 1
      ptr2 = ptr1 + ptr2
      brcopy(ptr1:ptr1)   = '('
      brcopy(ptr1+1:ptr2) = group(ia2)
      ptr1 = ptr2 + 1
      brcopy(ptr1:ptr1) = ')'
      ptr1 = ptr1 + 1

      ptr2 = INDEX(group(ia1),' ') - 1
      ptr2 = ptr1 + ptr2
      brcopy(ptr1:ptr1)   = '('
      brcopy(ptr1+1:ptr2) = group(ia1)
      ptr1 = ptr2 + 1
      brcopy(ptr1:ptr1) = ')'
      ptr1 = ptr1 + 1
      RETURN

    ENDIF
  ENDIF

! more than 2 groups: error
! -------------------------
  IF (ita > 2) THEN
      WRITE (6,*) '--error--,in mkbrcopy. # of methyl group > 2'
      STOP "in mkbrcopy"
  ENDIF
END SUBROUTINE mkbrcopy

!=======================================================================
! PURPOSE: This routine is very similar to prioty. Main difference is 
! that prioty is called to find the standardized formula of the molecule
! while brpri is called to find which branch has the highest priority.
! The subroutine check which formula of the branch in "copy" has the 
! highest priority and return the corresponding index in "copy"                                                  
!
!              -SEE PRIOTY FOR ADDITIONAL COMMENT -                   
!=======================================================================
SUBROUTINE brpri(ogroup,rank,copy,ncp,nring,maxpri)
  USE keyparameter, ONLY:mxring,mxcp
  USE rjtool
  IMPLICIT NONE

  CHARACTER(LEN=*),INTENT(IN)   :: ogroup(:) ! group table
  INTEGER,INTENT(IN)            :: rank(:)   ! rank (priority) of group i
  INTEGER,INTENT(IN)            :: ncp       ! # of formula in copy
  INTEGER,INTENT(IN)            :: nring     ! # of rings
  CHARACTER(LEN=*),INTENT(INOUT):: copy(:)   ! all possible writing to be checked 
  INTEGER,INTENT(OUT)           :: maxpri    ! index in copy with highest priority

  CHARACTER(LEN=LEN(ogroup(1))) :: group(SIZE(ogroup))
  CHARACTER(LEN=LEN(ogroup(1))) :: prigr
  INTEGER          :: i,j
  INTEGER          :: nelim,high,priort(mxcp),ncg
  INTEGER          :: k,gcntr,ps,ptr
  INTEGER          :: rjg(mxring,2),rjs(mxring,2)
  LOGICAL          :: lofind
  INTEGER          :: ngr

! -----------------------------------------------------
! INITIALIZE
! -----------------------------------------------------
  maxpri=0  ;  priort(:) = 0
! cp needed to work on group to handle ring (conflict with "intent(in)")      
  group(:)=ogroup(:)  ; ngr=SIZE(ogroup)
      
! -----------------------------------------------------
! GROUPS ALREADY RATED AND RANKED (ORDERED)
! -----------------------------------------------------
! note : the group exist in the molecule but may not exist in the branch
! => not all possible ranks may be present

! -----------------------------------------------------
! FIND PRIORITY AND ELIMINATE SUCCESSIVELY THE FORMULA
! -----------------------------------------------------
  nelim = 0

! loop over the groups, starting with the highest priority group
  IF (nring>0) CALL rjgrm(nring,group,rjg)

  aloop: DO i = 1,ngr
    DO j = 1,ngr
      IF (rank(j)==i) prigr = group(j)
    ENDDO
    ncg = INDEX(prigr,' ') - 1

! count the number of groups identical to prigr
    gcntr = 0
    DO j = 1,ngr
      IF (group(j)==prigr) gcntr = gcntr + 1
    ENDDO
 
! loop over the copies of the branch and set priority. Here the 
! group may not be in the formula of the branch => jump to next copy
    cploop: DO j = 1,ncp
      priort(j) = 0
      IF (INDEX(copy(j),prigr(1:ncg))==0) CYCLE cploop
      ps=0

! remove ring-join characters if present
      IF(nring > 0) CALL rjsrm(nring,copy(j),rjs)

! find the position of the group in the branch. Check that the group 
! was really found and not only a part of another group : e.g. CO in 
! CO(OH) or CdH in CdH2
      grloop: DO k = 1, gcntr
        lofind=.FALSE.
        findloop: DO
          IF (lofind) EXIT findloop
          ptr=INDEX(copy(j)(ps+1:),prigr(1:ncg))
          IF (ptr==0) EXIT grloop
          ps = ps + ptr
          IF ((copy(j)(ps+ncg:ps+ncg)=='C').OR.    &
             (copy(j)(ps+ncg:ps+ncg+1)=='(C').OR.  &
             (copy(j)(ps+ncg:ps+ncg+1)==')C').OR.  &
             (copy(j)(ps+ncg:ps+ncg+2)==')(C') ) lofind=.TRUE.
          IF ((copy(j)(ps+ncg:ps+ncg)=='c').OR.    &
             (copy(j)(ps+ncg:ps+ncg+1)=='(c').OR.  &
             (copy(j)(ps+ncg:ps+ncg+1)==')c').OR.  &
             (copy(j)(ps+ncg:ps+ncg+2)==')(c') ) lofind=.TRUE.
          IF ((copy(j)(ps+ncg:ps+ncg+2)==')-O').OR. &
             (copy(j)(ps+ncg:ps+ncg+1)=='-O')) lofind=.TRUE.
          IF (copy(j)(ps+ncg:ps+ncg)==' ')     lofind=.TRUE.
          IF (lofind) priort(j) = priort(j) + priort(j) + ps
!          ELSE
!            GOTO 17
!          ENDIF
        ENDDO findloop
      ENDDO grloop

! replace ring-join characters if necessary
      IF(nring > 0) CALL rjsadd(nring,copy(j),rjs)
    ENDDO cploop

! find the maximum
    high = 0
    DO j = 1,ncp
       high = MAX(priort(j),high)
    ENDDO

! eliminate copy with priority < maximum (stop if 1 copy remain)
    DO j = 1,ncp
       IF (priort(j) < high) THEN
         IF (copy(j)(1:1)/=' ') THEN
           copy(j) = ' '
           nelim = nelim + 1
           IF (nelim==ncp-1) EXIT aloop
         ENDIF
       ENDIF
    ENDDO

  ENDDO aloop

! ---------------------------------------------------------
! check that only one copy remains and return its ID number
! ---------------------------------------------------------
  IF (nelim==ncp-1) THEN
    DO i = 1,ncp 
      IF (copy(i)(1:1)/=' ') maxpri=i
    ENDDO
! too many formula left  ---> error:
  ELSE
    WRITE(6,'(a)') '--error-- in brpri. More than one copy left:'
    DO j=1,ncp  ;  WRITE(6,'(a)') copy(j)  ;  ENDDO
    STOP "in brpri"
  ENDIF

  IF(nring > 0) CALL rjgadd(nring,group,rjg)
END SUBROUTINE brpri



! ======================================================================
! PURPOSE : THIS ROUTINE IS VERY SIMILAR TO LNTREE. Main difference  
! is that lntree is called to find the longest tree in the molecule 
! while treebr is called to find the longest tree in a given branch of 
! the molecule. Treebr sets up the tree of the C-C bond starting at the  
! group given in the input (top) and evaluate the longest path (brpath), 
! its length (maxlng) and the total # of path having the maximum length.                
!
! SEE LNTREE FOR ADDITIONAL COMMENT                 
! ======================================================================
SUBROUTINE treebr(lobond,ig,top,sec,nca,brpath,maxlng,npath)
  IMPLICIT NONE

  LOGICAL,INTENT(IN) :: lobond(:,:) ! logical bond matrix 
  INTEGER,INTENT(IN) :: ig  ! group # to which the branch is bounded
  INTEGER,INTENT(IN) :: top ! 1st group => starting point
  INTEGER,INTENT(IN) :: sec ! 2nd group =>  bond to top 
  INTEGER,INTENT(IN) :: nca ! total number of group in the molecule
  INTEGER,INTENT(OUT):: npath  ! # of pathes having the maximum length
  INTEGER,INTENT(OUT):: maxlng ! maximum length of the pathes found
  INTEGER,INTENT(OUT):: brpath(:,:) ! trees (path) of the jth path found starting at 
                                    ! the ith group position (here "i" is "top")

  INTEGER :: left(SIZE(lobond,1)),right(SIZE(lobond,1))
  INTEGER :: center(SIZE(lobond,1)),parent(SIZE(lobond,1))
  INTEGER :: tlngth(SIZE(lobond,1))
  INTEGER :: tpath(SIZE(lobond,1),SIZE(lobond,2))
  INTEGER :: flag(SIZE(lobond,1))
  INTEGER :: ptr,knt,nknt,nct1,nct2
  INTEGER ::  iend,i,j,k
  LOGICAL :: tbond(SIZE(lobond,1),SIZE(lobond,2))

! -----------
! initialize
! -----------
  flag(:) = 0     ;  left(:) = 0    ;  right(:) = 0
  center(:) = 0   ;  parent(:) = 0  ;  tlngth(:) = 0
  brpath(:,:) = 0 ;  tpath(:,:) = 0 ;  ptr = 0  
  tbond(:,:) = lobond(:,:)

  left(top) = sec  ;  parent(top) = ig  ;  parent(sec) = top
  tbond(top,sec) = .FALSE.  ;  tbond(sec,top) = .FALSE.

! ---------------------------------------
! get the relationships (parent/children)
! ---------------------------------------
  DO k=1,nca
    nknt = 0
    DO i=1,nca
      knt = 0
      DO j=1,nca
        IF (tbond(i,j)) THEN
          knt = knt + 1
          ptr = j
        ENDIF
      ENDDO
      nknt = nknt + knt

! look for ith carbon with only one node, where parents do not
! exist: (ith carbon = child, jth carbon = parent)
      IF (knt==1) THEN
        IF (parent(i)==0 .AND. i/=top .AND. i/=ig) THEN
          parent(i) = ptr
          tbond(i,ptr) = .FALSE.  ;  tbond(ptr,i) = .FALSE.
          IF(left(ptr)==0) THEN
            left(ptr) = i
          ELSE IF(right(ptr)==0) THEN
            right(ptr) = i
          ELSE IF(center(ptr)==0) THEN
            center(ptr) = i
          ELSE

! if all children taken, error in bond-matrix:
            DO j=1,nca
              WRITE(6,*) (lobond(j,iend),iend=1,nca)
            ENDDO
            WRITE(6,'(a)') '--error-- in  treebr. No path, bonding pb'
            STOP "in treebr"
          ENDIF
        ENDIF
      ENDIF
    ENDDO
! do loop until bond-matrix for the branch is (0).
    IF (nknt==0) EXIT
  ENDDO

! ---------------------------------------------
! define all top-down paths starting at "top"
! ---------------------------------------------
  nct1 = nca-1  ;  nct2 = nca+4
  brloop: DO i=1,nct1
    ptr = top  ;  tpath(i,1) = top
    DO j=2,nct2
      IF (flag(ptr)==0) THEN
        IF (left(ptr)/=0) THEN ; ptr = left(ptr) ;  tpath(i,j) = ptr
        ELSE                   ; flag(ptr) = 1
        ENDIF
      ELSE IF(flag(ptr)==1) THEN
        IF (right(ptr)/=0) THEN ; ptr = right(ptr) ; tpath(i,j) = ptr
        ELSE                    ; flag(ptr) = 2
        ENDIF
      ELSE IF(flag(ptr)==2) THEN
        IF (center(ptr)/=0) THEN ; ptr = center(ptr) ; tpath(i,j) = ptr
        ELSE                     ; flag(ptr) = 3
              flag(parent(ptr)) = flag(parent(ptr)) + 1
              CYCLE brloop
        ENDIF
      ELSE IF (flag(ptr)==3) THEN
        flag(parent(ptr)) = flag(parent(ptr)) + 1
        CYCLE brloop
      ENDIF
    ENDDO
  ENDDO brloop

! ---------------------
! get the longest path 
! ---------------------

! get the length of each path and found the maximum length
  maxlng = 0
  DO i=1,nca
    DO j=1,nca
      IF (tpath(i,j)/=0) tlngth(i)=tlngth(i) + 1
    ENDDO
    IF (tlngth(i) > maxlng) maxlng = tlngth(i)
  ENDDO

! find all paths having the maximum length
  npath = 0
  DO i=1,nca
    IF (tlngth(i)==maxlng) THEN
      npath=npath+1
      DO j=1,nca
        brpath(npath,j)=tpath(i,j)
      ENDDO          
    ENDIF
  ENDDO

! check if everything is ok
  IF (npath==0) THEN
    WRITE(6,'(a)') '--error-- in treebr. No path found'
    STOP "in treebr"
  ENDIF
  IF (maxlng < 2) THEN
    WRITE(6,'(a)') '--error-- in treebr. Length of the path shorter than expected'
    STOP "in treebr"
  ENDIF
! delta position in the branch can only be reached for C14 
! no delta position (maxlng=4) is expected
  IF (maxlng > nca-3) THEN
    WRITE(6,'(a)') '--error-- in treebr. Length of the path greater than expected'
    STOP "in treebr"
  ENDIF
END SUBROUTINE treebr

END MODULE brchtool
  
