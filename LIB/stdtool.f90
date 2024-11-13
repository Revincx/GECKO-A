MODULE stdtool
  IMPLICIT NONE
  CONTAINS

! SUBROUTINE lntree(bond,dbflg,top,sec,nca,clngth,path)
! SUBROUTINE ckgrppt(locat,group)
! SUBROUTINE prioty(ogroup,rank,ocopy,ncp,nring,chem)
! SUBROUTINE mkcopy(lobond,group,nca,rank,nring,ig,pg,ng,ptr1,copy)
! SUBROUTINE dwrite(cpchem)
! SUBROUTINE revers(copy,cc)


! ================================================================
! PURPOSE : Set-up the tree of the C-C and C-O-C bond starting at    
!           the group in the input (top) and evaluate the longest    
!           path and its length.                                     
!    
!   The purpose of this routine is to evaluate the longest chain of  
!   the molecule. With the information in the bond-matrix, it is     
!   possible to evaluate the top-down trees of the molecule:         
!   e.g.                  CO(CH3)CH(CH(OH)(CH3))CHO                  
!                         1   2   3   4     5    6                   
!                   ___                                              
!                  | 1 |  parent of 2,3                              
!                   ---                                              
!                 /     \                                            
!             ___        ___                                         
!  child of 1| 2 |      | 3 | child of 1, parent of 4,6              
!             ---        ---                                         
!                      /     \                                       
!                   ___       ___                                    
!   child of 3 &   | 4 |     | 6 |  child of 3                       
!   parent of 5     ---       ---                                    
!                 /                                                  
!             ___                                                    
!            | 5 |  child of 4                                       
!             ---                                                    
!                                                                    
!   A child on the left hand is called "LEFT", on the right hand     
!   "RIGHT" and in the middle "CENTER".     
!   In CLNGTH the length of each longest tree is stored.  
!   Nevertheless the longest chain is the chain with most of the     
!   double bonds in it.                                              
!   First all relationships (parent - children) are evaluated. Then  
!   the top-down paths are defined so that first a LEFT, if not      
!   available, a RIGHT and then a CENTER is taken. As often as at    
!   least one child exists the specified path is followed further on.
!   In the last section the longest paths and the paths with the     
!   most double-bonds in it still remain and are given to the cal-   
!   ling routine.                                                    
!                                                                    
!********************************************************************
SUBROUTINE lntree(bond,top,sec,nca,clngth,path)
  IMPLICIT NONE
 
  INTEGER,INTENT(IN)  :: bond(:,:)  ! node-node bonfd matrix
  INTEGER,INTENT(IN)  :: top        ! 1st node: starting point
  INTEGER,INTENT(IN)  :: sec        ! 2nd node: carbon next to top
  INTEGER,INTENT(IN)  :: nca        ! # of nodes
  INTEGER,INTENT(INOUT):: clngth(:,:)  ! length(i,j) of jth path found starting        
                                       ! at the ith node (here "i" is "top")
  INTEGER,INTENT(OUT) :: path(:,:,:) ! trees (path(i,j,k)) of the jth path 
                                     ! found starting at the ith node (here "i" is "top")
! internal:
  INTEGER :: left(SIZE(bond,1)),right(SIZE(bond,1))
  INTEGER :: center(SIZE(bond,1)),parent(SIZE(bond,1))
  INTEGER :: flag(SIZE(bond,1))
  INTEGER :: ptr,knt,nknt,nct1,nct2
  INTEGER :: maxlng,iend,i,j,k
  LOGICAL :: tbond(SIZE(bond,1),SIZE(bond,2))

! -----------
! initialize
! -----------
  ptr=0  ;  flag(:)=0  ;  left(:)=0  ;  right(:)=0
  center(:)=0  ;  parent(:)=0

! make a logical bond matrix (tbond only used to find parent and child)
  tbond(:,:)=.FALSE.  ;  WHERE (bond(:,:)/=0) tbond(:,:)=.TRUE.
  left(top)=sec           ;  parent(sec)=top
  tbond(top,sec)=.FALSE.  ;  tbond(sec,top)=.FALSE.

! ---------------------------------------
! get the relationships (parent/children)
! ---------------------------------------
  DO k=1,nca
    nknt=0
    DO i=1,nca
      knt=0
      DO j=1,nca
        IF (tbond(i,j)) THEN
          knt=knt+1
          ptr=j
        ENDIF
      ENDDO
      nknt=nknt+knt

! look for ith carbon with only one node, where parents do not
! exist: (ith carbon=child,jth carbon=parent)
      IF (knt==1) THEN
        IF ( (parent(i)==0).and.(i.NE.top) ) THEN
          parent(i)=ptr
          tbond(i,ptr)=.FALSE.  ;  tbond(ptr,i)=.FALSE.
          IF (left(ptr)==0)        THEN ; left(ptr)=i
          ELSE IF (right(ptr)==0)  THEN ; right(ptr)=i
          ELSE IF (center(ptr)==0) THEN ; center(ptr)=i
          ELSE

! if all children taken, error in bond-matrix:
            DO j=1,nca
              WRITE(6,*) (bond(j,iend),iend=1,nca)
            ENDDO
            WRITE(6,'(a)') '--error-- in lntree, '
            WRITE(6,'(a)') 'no path possible, error in bonding'
            STOP  "in lntree"
          ENDIF
        ENDIF
      ENDIF
    ENDDO
! do loop until bond-matrix is (0).
    IF (nknt==0) EXIT
  ENDDO

! ---------------------------------------------
! define all top-down pathes starting at "top"
! ---------------------------------------------
  nct1=nca-1 ; nct2=nca+4
  grloop: DO i=1,nct1
    ptr=top
    path(top,i,1)=top
    DO j=2,nct2
      IF (flag(ptr)==0) THEN
        IF (left(ptr)/=0) THEN   ; ptr=left(ptr)   ; path(top,i,j)=ptr
        ELSE                     ; flag(ptr)=1
        ENDIF
      ELSE IF (flag(ptr)==1) THEN 
        IF (right(ptr)/=0) THEN  ; ptr=right(ptr)  ; path(top,i,j)=ptr
        ELSE                     ; flag(ptr)=2
        ENDIF
      ELSE IF(flag(ptr)==2) THEN
        IF (center(ptr)/=0) THEN ; ptr=center(ptr) ; path(top,i,j)=ptr
        ELSE                     ; flag(ptr)=3
            flag(parent(ptr))=flag(parent(ptr))+1
            CYCLE grloop 
        ENDIF
      ELSE IF (flag(ptr)==3) THEN
        flag(parent(ptr))=flag(parent(ptr))+1
        CYCLE grloop
      ENDIF
    ENDDO
  ENDDO grloop

! ---------------------
! get the longest path 
! ---------------------
  maxlng=0
  DO i=1,nca
    DO j=1,nca
      IF (path(top,i,j)/=0) clngth(top,i)=clngth(top,i)+1
    ENDDO
 
    IF (clngth(top,i) > maxlng) THEN
      maxlng=clngth(top,i)
      iend=i-1 ; clngth(top,1:iend)=0

    ELSE IF (clngth(top,i) < maxlng) THEN
      clngth(top,i)=0
    ENDIF
  ENDDO

END SUBROUTINE lntree

!=======================================================================
! PURPOSE : Check for misordered functionalities in parentheses next to 
! each other in the group given as input. The functionalities in the 
! group are sorted at the output For group priority see "pri" in module
! keyparameter where the priority invers TOP-DOWN.                                 
!                                                                     
!                       PTR21  PTR23   LNG3                           
!                         |      |____/___                            
!                C(..X1..)(..X2..)(..X3..)                            
!                 |      |        |      |                            
!               PTR11  PTR12    PTR32  PTR33                          
!=======================================================================
SUBROUTINE ckgrppt(locat,group)
  USE keyparameter, ONLY: pri
  IMPLICIT NONE

  INTEGER,INTENT(IN) :: locat    ! pointer of first string ")(" in group
  CHARACTER(LEN=*),INTENT(INOUT) :: group ! the group that must be checked.

  CHARACTER(LEN=LEN(group)) :: tgroup
  CHARACTER(LEN=LEN(group)) :: grp1, grp2, grp3
  CHARACTER*3 :: grp(3)
  INTEGER     :: ptr11, ptr12, ptr21
  INTEGER     :: ptr23, ptr32, ptr33
  INTEGER     :: cflg,i,p,lng3
  INTEGER     :: lengr
  LOGICAL     :: tri,loswitch,lopar

! -----------
! initialise
! -----------
  grp1=' '  ;  grp2=' '  ;  grp3=' '
  ptr11=0   ;  ptr23=0   ;  ptr33=0  ;  lng3=0
  tri = .FALSE.  ;   loswitch = .FALSE.  ; lengr=LEN(group)

! ---------------------------------------
! find where parentheses open and close
! ---------------------------------------

! set the pointers for the first pair of unequal parentheses:
  ptr12 = locat  ;  ptr21 = locat + 1

! search first group
! ------------------

! ptr11 is location of '(' match for ')' of the first pair ')(':
  p = 1 ; lopar=.TRUE.
  DO i=ptr12-1,2,-1
    IF (group(i:i)==')') p = p+1
    IF (group(i:i)=='(') p = p-1
    IF (p==0) THEN
      ptr11 = i
      lopar=.FALSE.
      EXIT
    ENDIF
  ENDDO

! check error,  parenthesis mismatch
  IF (lopar) THEN
    WRITE(6,*) '--error--, in ckgrppt. First group of parenthesis '
    WRITE(6,*) ' mismatch for the group :', TRIM(group)
    STOP "in ckgrppt"
  ENDIF

! search second group
! ------------------

! ptr23 is location of ')' match for '(' of the first pair ')(':
  p = 1  ;  lopar=.TRUE.
  DO i = ptr21+1,lengr
    IF (group(i:i)=='(') p = p+1
    IF (group(i:i)==')') p = p-1 
    IF (p==0) THEN
      ptr23 = i
      lopar=.FALSE.
      EXIT
    ENDIF
  ENDDO

! check error,  parenthesis mismatch
  IF (lopar) THEN
    WRITE(6,*) '--error--, in ckgrppt. Second group of parenthesis '
    WRITE(6,*) 'mismatch for the group :', TRIM(group)
    STOP "in ckgrppt"
  ENDIF

! search third group (if any)
! ------------------

! see if there are triple groups in parentheses instead of a pair. 
! PTR32 is *'(' and PTR33 is *')' of the 3rd group in parentheses. 
! If no 3rd group, then no char is expected after the 2nd group of parentheses.

  ptr32 = ptr23 + 1 ; lopar=.TRUE.
  IF (group(ptr32:ptr32)=='(') THEN
    tri = .TRUE.
    p = 1
    DO i=ptr32 + 1,lengr
      IF (group(i:i)==')') p = p-1 
      IF (group(i:i)=='(') p = p+1
      IF(p==0) THEN
        ptr33 = i
        lopar=.FALSE.
        EXIT
      ENDIF
    ENDDO
    IF (lopar) THEN
      WRITE(6,*) '--error--, in ckgrppt. Third group of parenthesis '
      WRITE(6,*) 'mismatch for the group :', TRIM(group)
      STOP "in ckgrppt"
    ENDIF

  ELSE IF (group(ptr32:ptr32)/=' '.AND.group(ptr32:ptr32)/='.') THEN
    WRITE(6,*) '--error--, in ckgrppt. " " or "." is expected after'
    WRITE(6,*) 'second group of parenthesis in :', TRIM(group)
    STOP "in ckgrppt"
  ENDIF 

! -----------------
! define all groups
! -----------------

! 3 characters are sufficient to distinguish among the various functional group 
! and this 3 char subset (grp(i)) are used to define the priorities. Full
! group (grpi, includes the parenthesis) are also identified to make the switches.
  grp1 = group(ptr11:ptr12)  ;  grp(1) = group(ptr11+1:ptr11+3)
  grp2 = group(ptr21:ptr23)  ;  grp(2) = group(ptr21+1:ptr21+3)

  IF (tri) THEN
    grp3 = group(ptr32:ptr33)  ;  grp(3) = group(ptr32+1:ptr32+3)
    lng3   = ptr33 - ptr23

! check nothing left after the 3rd parenthesis
! note: C(X)(Y)(Z). structure might be produced from C-C break
    IF (group(ptr33+1:ptr33+1)/=' ') THEN
      IF (group(ptr33+1:ptr33+2)/='. ') THEN  
        WRITE(6,*) '--error--, in ckgrppt. A " " is expected after'
        WRITE(6,*) 'third group of parenthesis in :', TRIM(group)
        STOP "in ckgrppt"
      ENDIF
    ENDIF 
  ENDIF

! fluor is a special case, since it is only one character long:
  DO i=1,3
    IF (grp(i)(1:1)=='F') grp(i) = 'F  '
  ENDDO

! --------------------------
! switch group, if necessary
! --------------------------

! take and save what's sure
  tgroup = group(1:(ptr11 - 1))

! if the first two groups are wrongly placed then switch
  IF (INDEX(pri,grp(1)) < INDEX(pri,grp(2))) THEN
    loswitch=.true.
    ptr21 = ptr11 + (ptr23 - ptr12)
    tgroup(ptr11:lengr) = grp2
    tgroup(ptr21:lengr) = grp1
    IF (tri) tgroup(ptr32:lengr) = group(ptr32:lengr)
    group = tgroup
    ptr12 = ptr21 - 1
  ENDIF

! if a third group exists then switch again as needed
  IF (tri) THEN
    tgroup = ' '

! define if third  group has higher priority than the two others:
    cflg = 0
    DO i=1,2
      IF (INDEX(pri,grp(i)) < INDEX(pri,grp(3))) cflg=cflg+1
    ENDDO

    IF (cflg==1) THEN
      tgroup = group(1:ptr12)
      tgroup(ptr21:lengr) = grp3
      IF (loswitch) THEN
         tgroup(ptr21+lng3:lengr) = grp1
      ELSE
         tgroup(ptr21+lng3:lengr) = grp2
      ENDIF
      group = tgroup
    ELSE IF (cflg==2) THEN
      tgroup = group(1:(ptr11 - 1))
      tgroup(ptr11:lengr) = grp3
      tgroup((ptr11 + lng3):lengr) = group(ptr11:ptr23)
      group = tgroup
    ENDIF
  ENDIF

END SUBROUTINE ckgrppt


! ======================================================================
! PURPOSE: Check which formula in COPY has the highest priority,    
! comparing the position of different functional groups and return the
! "standardized" formula. 
! If the chemical is a radical, the formulas remain where the radical 
! group is at the end of the formula. If there are still more than one 
! writing left or the chemical is a non-radical molecule, according to 
! the group priorities, the formulas with the groups at the end remain, 
! respectively. This is done for all functional groups in the molecule 
! unless there are still more than one formulas left.                                         
! ======================================================================
SUBROUTINE prioty(ogroup,rank,ocopy,ncp,nring,chem)
  USE keyparameter, ONLY: mxring
  USE rjtool
  IMPLICIT NONE

  CHARACTER(LEN=*),INTENT(IN) :: ocopy(:)  ! all possible writing to be  checked 
  CHARACTER(LEN=*),INTENT(IN) :: ogroup(:) ! group table
  INTEGER,INTENT(IN) :: ncp                ! # of formulas in copy
  INTEGER,INTENT(IN) :: nring              ! # of separate rings in CHEM 
  INTEGER,INTENT(IN) :: rank(:)            ! rank (priority) of group i 
  CHARACTER(LEN=*),INTENT(OUT):: chem      ! standardized formula

! internal:
  INTEGER :: i,j,ef(5)
  INTEGER :: nelim,high,priort(SIZE(ocopy)),ncg
  INTEGER :: k,gcntr,eflag,ps,ptr
  INTEGER :: rjg(mxring,2),rjs(mxring,2)
  LOGICAL :: lofind
  INTEGER :: temprjs(SIZE(ocopy),mxring,2)
  INTEGER :: minrjc, maxrjc
  INTEGER :: numgr
  CHARACTER(LEN=LEN(ogroup(1))):: prigr
  CHARACTER(LEN=LEN(ogroup(1))):: group(SIZE(ogroup))
  CHARACTER(LEN=LEN(chem))     :: tempcopy(SIZE(ocopy))
  CHARACTER(LEN=LEN(chem))     :: copy(SIZE(ocopy))  
     
! --------------------------------------------------------
! initialize:
! --------------------------------------------------------
  chem=' '  ;  priort(:)=0  
! ba - copies (needed to keep intent in) => to be revisited  
  copy(:) = ocopy(:)  ; group(:)=ogroup(:)  
  numgr=SIZE(ogroup)
  
! --------------------------------------------------------
! CHECK RADICAL 
! --------------------------------------------------------

! if different copies have different types of radicals, then error...
  ef(:)=0
  DO i=1,ncp
    IF (INDEX(copy(i),'.(OO.)')/=0)       THEN ; ef(1) = 1
    ELSE IF (INDEX(copy(i),'CO(OO.)')/=0) THEN ; ef(2) = 1
    ELSE IF (INDEX(copy(i),'(OO.)')/=0)   THEN ; ef(3) = 1
    ELSE IF (INDEX(copy(i),'(O.)')/=0)    THEN ; ef(4) = 1
    ELSE                                       ; ef(5) = 1
    ENDIF
  ENDDO

  eflag = SUM(ef)
  IF (eflag > 1) THEN
    WRITE(6,'(a)') '--error--in prioty. Different radicals for molecule:'
    WRITE(6,'(a)') TRIM(chem)
    STOP "in prioty"
  ENDIF

! --------------------------------------------------------
!  FUNCTIONAL GROUPS WERE ALREADY RATED AND ORDERED
! --------------------------------------------------------

! -----------------------------------------------------
! FIND PRIORITY AND ELIMINATE SUCCESSIVELY THE FORMULA
! -----------------------------------------------------
  nelim = 0
      
! loop over the groups, starting with the highest priority group. 
  IF (nring > 0) CALL rjgrm(nring,group,rjg) 

  aloop: DO i = 1,numgr  ! loop ranks: don't expect to complete sequence
    loop1: DO j = 1,numgr
      IF (rank(j)==i) THEN
        prigr = group(j)
        EXIT loop1
      ENDIF
    ENDDO loop1
    ncg = INDEX(prigr,' ') - 1

! count the number of groups identical to prigr
    gcntr = 0
    DO j = 1,numgr
      IF (group(j)==prigr) gcntr = gcntr + 1
    ENDDO

! scroll copies and set priority (jump over the copy already removed)
    cploop: DO j = 1,ncp
      priort(j) = 0
      ps=0
      IF (copy(j)(1:1)==' ') CYCLE cploop

      IF (nring > 0) CALL rjsrm(nring,copy(j),rjs) ! rm ring-join char

! if the group exist more than once in the molecule, then found each of
! them (thus the loop over gcntr). Check that the group was really found
! and not only a part of another group: e.g. CO in CO(OH) or CdH in CdH2
      DO k = 1, gcntr
        lofind=.FALSE.
        loloop: DO
          IF (lofind) EXIT loloop
          ptr = INDEX(copy(j)(ps+1:),prigr(1:ncg))
          IF (ptr==0) THEN
            WRITE(6,*) '--error--, from subroutine prioty'
            WRITE(6,*) 'Problem to detect ',(prigr(1:ncg))
            WRITE(6,*) 'in ',TRIM(copy(j))
            STOP "in prioty"
          ENDIF
          ps=ps+ptr

          IF ((copy(j)(ps+ncg:ps+ncg)=='C')  .OR. &
             (copy(j)(ps+ncg:ps+ncg+1)=='(C').OR. &
             (copy(j)(ps+ncg:ps+ncg+1)==')C').OR. &
             (copy(j)(ps+ncg:ps+ncg+2)==')(C') ) lofind=.TRUE.
          IF ((copy(j)(ps+ncg:ps+ncg)=='c')  .OR. &
             (copy(j)(ps+ncg:ps+ncg+1)=='(c').OR. &
             (copy(j)(ps+ncg:ps+ncg+1)==')c').OR. &
             (copy(j)(ps+ncg:ps+ncg+2)==')(c') ) lofind=.TRUE.
          IF ((copy(j)(ps+ncg:ps+ncg+1)=='-O').OR. &
             (copy(j)(ps+ncg:ps+ncg+2)=='(-O').OR. &
             (copy(j)(ps+ncg:ps+ncg+2)==')-O').OR. &
             (copy(j)(ps+ncg:ps+ncg+3)==')(-O') ) lofind=.TRUE.
          IF  (copy(j)(ps+ncg:ps+ncg)==' ') lofind=.TRUE.
        ENDDO loloop  

        priort(j) = priort(j) + priort(j) + ps
      ENDDO

      IF (nring > 0) CALL rjsadd(nring,copy(j),rjs) ! add ring-join char

    ENDDO cploop

! find the maximum
     high = MAXVAL(priort(1:ncp))
     
! eliminate copy with priority < maximum (stop if 1 copy remain)
    DO j = 1,ncp
       IF (priort(j) < high) THEN
         copy(j) = ' '
         IF(priort(j) > 0) nelim = nelim + 1
         IF (nelim==ncp-1) EXIT aloop
       ENDIF
    ENDDO

  ENDDO aloop

! If 2 rings with symetrical shape, then no possibility to discriminate
! among the remaining formula (eg in C12HCH2CH2CH2CH(C2H2)C1H2). In that
! case, give priority to the formula having the largest position of
! the first "ring joining character" ("1") in copy, then the second ...
  IF (nelim/=ncp-1) THEN
    IF (nring >= 2) THEN
! intialize and get ring joining character position
      DO i=1,ncp
        tempcopy(i)=' '
        DO j=1,mxring
          DO k=1,2
            temprjs(i,j,k)=0
          ENDDO
        ENDDO
        IF (copy(i)(1:1)/=' ') THEN
          CALL rjsrm(nring,copy(i),rjs)
          tempcopy(i)=copy(i)
          DO j=1,mxring
            DO k=1,2
              temprjs(i,j,k)=rjs(j,k)
            ENDDO
          ENDDO
          CALL rjsadd(nring,copy(i),rjs)
        ENDIF
      ENDDO

! check that all copies are identical (except ring joining character)
      DO i=1,ncp-1
        IF ((tempcopy(i)(1:1)/=' ').AND. &
           (tempcopy(i+1)(1:1)/=' ')) THEN
           IF (tempcopy(i)/=tempcopy(i+1)) THEN
            WRITE(6,'(a)') '--error-- in prioty. Different formula in copy:'
            WRITE(6,*) TRIM(copy(i))
            WRITE(6,*) TRIM(copy(i+1))
            STOP "in prioty"
          ENDIF
        ENDIF
      ENDDO

! select the correct formula
      selloop: DO j=1,mxring
        minrjc=LEN(chem)
        DO i=1,ncp  ! minimum
          IF (copy(i)(1:1)/=' ') THEN
            minrjc=MIN(temprjs(i,j,1),minrjc)
          ENDIF
        ENDDO
        DO i=1,ncp
          IF (temprjs(i,j,1)/=minrjc) THEN
            IF (copy(i)/=' ') THEN
              nelim=nelim+1
              copy(i)=' '
            ENDIF
            IF (nelim==ncp-1) EXIT selloop
          ENDIF
        ENDDO
        maxrjc=0
        DO i=1,ncp  ! maximum
          maxrjc=MAX(temprjs(i,j,2),maxrjc)
        ENDDO
        DO i=1,ncp
          IF (temprjs(i,j,2)/=maxrjc) THEN
            IF (copy(i)/=' ') THEN
              nelim=nelim+1
              copy(i)=' '
            ENDIF
            IF (nelim==ncp-1) EXIT selloop
          ENDIF
        ENDDO
      ENDDO  selloop
    ENDIF
  ENDIF

! ---------------------------------------------------------
! check that only one copy remains and return 
! ---------------------------------------------------------

  IF (nelim==ncp-1) THEN
    DO i = 1,ncp 
      IF (copy(i)/=' ') chem = copy(i)
    ENDDO

! too many formulae left  ---> error:
  ELSE
    WRITE(6,'(a)') '--error-- in prioty. More than 1 copy left:'
    DO j=1,ncp  ;  WRITE(6,'(a)') TRIM(copy(j))  ;  ENDDO
    STOP  "in prioty"
  ENDIF
  
  IF (nring > 0) CALL rjgadd(nring,group,rjg) 
END SUBROUTINE prioty


! ===================================================================
! PURPOSE: make a copy of a molecule according to the longest tree.    
! On each call of mkcopy, the group ig is written in "copy". At the 
! same time and according to the bond-matrix, all groups having a bond 
! to ig (except on longest tree) are attached to it (with "()").      
!
! Ramifications are written according some priority rules. A max of 2  
! ramifications is expected and the formula must have the form :       
!        longest_tree-Cig(branch1)(branch2)-longest_tree               
! Since the branches may also contain ramification and functionalities,
! each branch must be written in a unique way. This is done by the        
! subroutine getbrch. If the carbon ig 'carry' 2 branches, then        
! branch1 and branch2 must be evaluated to know which branch has the   
! highest priority and must be written first.                          
!
! The structure of the subroutine is :                                 
!  1 - write the group ig and check if branching exist at pos. ig         
!  2 - if no branching => return                                       
!  3 - if only 1 branching => get and write the branch                 
!  4 - if 2 branching => get each branches, evaluate the priority of   
!       each branch and write them.    
!  It is expected that in most cases, the branch will only be 1C long. 
!  This case is tested first, since there is only 1 way to write the   
!  the branch.  
! ===================================================================
SUBROUTINE mkcopy(lobond,group,nca,rank,nring,ig,pg,ng,ptr1,copy)
  USE brchtool, ONLY:getbrch,brpri 
  IMPLICIT NONE

  LOGICAL,INTENT(IN) :: lobond(:,:)     ! logical bond matrix
  CHARACTER(LEN=*),INTENT(IN) :: group(:) ! group table 
  INTEGER,INTENT(IN) :: rank(:)   ! rank (priority) of group i (low = more important)
  INTEGER,INTENT(IN) :: nca       ! total number of group in the molecule
  INTEGER,INTENT(IN) :: nring     ! # of ring
  INTEGER,INTENT(IN) :: ig        ! is next group to be written to the copy
  INTEGER,INTENT(IN) :: pg        ! previous group of ig in the longest tree
  INTEGER,INTENT(IN) :: ng        ! next group of ig in the longest tree
  CHARACTER(LEN=*),INTENT(OUT) :: copy  ! chemical formula under that is curently written
  INTEGER,INTENT(INOUT) :: ptr1   ! pointer in "copy", where to put next group

! internal:
  INTEGER :: ptr2,i
  INTEGER :: ialpha,ia1,ia2,ib1,ib2
  INTEGER :: ml1,ml2,maxpri
  CHARACTER(LEN=LEN(copy)) :: tempbr(SIZE(group))
  CHARACTER(LEN=LEN(copy)) :: brch1,brch2

! -----------------------------------------------------------
! WRITE GROUP IG TO COPY
! -----------------------------------------------------------
  ptr2 = ptr1 + INDEX(group(ig),' ') - 2
  copy(ptr1:ptr2) = group(ig)
  ptr1 = ptr2 + 1               

! -----------------------------------------------------------
! LOOK FOR ALL ATTACHED GROUPS TO GROUP IG, BUT NOT PG OR NG
! -----------------------------------------------------------
! search number attached group to IG in "alpha" position
  ialpha=0  ;  ia1=0  ;  ia2=0
  DO i=1,nca
    IF (lobond(ig,i)) THEN
      IF ( (i/=pg) .AND. (i/=ng) ) THEN
        ialpha=ialpha+1
        IF (ialpha==1) ia1=i
        IF (ialpha==2) ia2=i
      ENDIF
    ENDIF
  ENDDO

! -----------------------------------------------------------
! IF NO ALPHA GROUP, THEN RETURN
! -----------------------------------------------------------
  IF (ialpha==0) RETURN

! -----------------------------------------------------------
! ONE ALPHA GROUP
! -----------------------------------------------------------
  IF (ialpha==1) THEN

! check if at least a beta position is occupied
    ib1=0
    DO i=1,nca
      IF (lobond(ia1,i)) THEN
        IF (i/=ig) THEN
          ib1=i
          EXIT
        ENDIF
      ENDIF
    ENDDO
        
! no beta position found : write group alpha and return
    IF (ib1==0) THEN
      ptr2=INDEX(group(ia1),' ') - 1
      ptr2 = ptr1 + ptr2
      copy(ptr1:ptr1)   = '('
      copy(ptr1+1:ptr2) = group(ia1)
      ptr1 = ptr2 + 1
      copy(ptr1:ptr1)=')'
      ptr1=ptr1+1
      RETURN

! at least a beta position found
    ELSE

! search the longest tree starting from ia1 with the highest priority
      CALL getbrch(lobond,group,nca,rank,nring,ig,ia1,ib1,brch1,ml1)
      ptr2=INDEX(brch1,' ') - 1
      ptr2 = ptr1 + ptr2
      copy(ptr1:ptr1)   = '('
      copy(ptr1+1:ptr2) = brch1
      ptr1 = ptr2 + 1
      copy(ptr1:ptr1)=')'
      ptr1=ptr1+1
      RETURN
    ENDIF
  ENDIF

! -----------------------------------------------------------
! TWO ALPHA GROUP
! -----------------------------------------------------------

! 4 various cases need to be considered :
! case 1 : the 2 branches are C1 branches
! case 2 : branch 1 is a "long" branch, branch 2 is a C1 branch
! case 3 : branch 1 is a C1 branch, branch 2 is a "long" branch
! case 4 : the 2 branches are "long" branches

  IF (ialpha==2) THEN

! check if at least a beta position is occupied
    ib1=0
    DO i=1,nca
      IF (lobond(ia1,i)) THEN
        IF (i/=ig) THEN
          ib1=i
          EXIT
        ENDIF
      ENDIF
    ENDDO

    ib2=0
    DO i=1,nca
      IF (lobond(ia2,i)) THEN
        IF (i/=ig) THEN
          ib2=i
          EXIT
        ENDIF
      ENDIF
    ENDDO

! get the branches starting at ia1 and ia2
    IF (ib1/=0) THEN
      CALL getbrch(lobond,group,nca,rank,nring,ig,ia1,ib1,brch1,ml1)
    ENDIF
    IF (ib2/=0) THEN
      CALL getbrch(lobond,group,nca,rank,nring,ig,ia2,ib2,brch2,ml2)
    ENDIF

! Case 1 :
! --------

! 2 methyl substitued => find the group with the highest priority
! and write the group having the lowest priority first, then return
    IF ((ib1==0).AND.(ib2==0)) THEN
      IF(rank(ia2) < rank(ia1))THEN
        ptr2 = INDEX(group(ia1),' ') - 1
        ptr2 = ptr1 + ptr2
        copy(ptr1:ptr1)   = '('
        copy(ptr1+1:ptr2) = group(ia1)
        ptr1 = ptr2 + 1
        copy(ptr1:ptr1) = ')'
        ptr1 = ptr1 + 1

        ptr2 = INDEX(group(ia2),' ') - 1
        ptr2 = ptr1 + ptr2
        copy(ptr1:ptr1)   = '('
        copy(ptr1+1:ptr2) = group(ia2)
        ptr1 = ptr2 + 1
        copy(ptr1:ptr1) = ')'
        ptr1 = ptr1 + 1
      ELSE
        ptr2 = INDEX(group(ia2),' ') - 1
        ptr2 = ptr1 + ptr2
        copy(ptr1:ptr1)   = '('
        copy(ptr1+1:ptr2) = group(ia2)
        ptr1 = ptr2 + 1
        copy(ptr1:ptr1) = ')'
        ptr1 = ptr1 + 1

        ptr2 = INDEX(group(ia1),' ') - 1
        ptr2 = ptr1 + ptr2
        copy(ptr1:ptr1)   = '('
        copy(ptr1+1:ptr2) = group(ia1)
        ptr1 = ptr2 + 1
        copy(ptr1:ptr1) = ')'
        ptr1 = ptr1 + 1
      ENDIF
      RETURN
    ENDIF


! Case 2 :
! --------

! ib2=0 : write first branch2 (methyl group) then branch 1 (long chain)
    IF ( (ib1/=0).AND.(ib2==0) ) THEN
      ptr2 = INDEX(group(ia2),' ') - 1
       ptr2 = ptr1 + ptr2
      copy(ptr1:ptr1)   = '('
      copy(ptr1+1:ptr2) = group(ia2)
      ptr1 = ptr2 + 1
      copy(ptr1:ptr1) = ')'
      ptr1 = ptr1 + 1

      ptr2=INDEX(brch1,' ') - 1
      ptr2 = ptr1 + ptr2
      copy(ptr1:ptr1)   = '('
      copy(ptr1+1:ptr2) = brch1
      ptr1 = ptr2 + 1
      copy(ptr1:ptr1)=')'
      ptr1=ptr1+1
      RETURN
    ENDIF

! Case 3 :
! --------

! ib1=0 : write first branch1 (methyl group) then branch 2 (long chain)
    IF ( (ib1==0).AND.(ib2/=0) ) THEN
      ptr2 = INDEX(group(ia1),' ') - 1
      ptr2 = ptr1 + ptr2
      copy(ptr1:ptr1)   = '('
      copy(ptr1+1:ptr2) = group(ia1)
      ptr1 = ptr2 + 1
      copy(ptr1:ptr1) = ')'
      ptr1 = ptr1 + 1

      ptr2=INDEX(brch2,' ') - 1
      ptr2 = ptr1 + ptr2
      copy(ptr1:ptr1)   = '('
      copy(ptr1+1:ptr2) = brch2
      ptr1 = ptr2 + 1
      copy(ptr1:ptr1)=')'
      ptr1=ptr1+1
      RETURN
    ENDIF

! Case 4 :
! --------

! ib1 and ib2 are both not equal to 0. Check the length of the chain (mli)
! and write the shortest chain first. If the chains are of the same
! length then call brpri to find priority.
    IF ( (ib1/=0).and.(ib2/=0) ) THEN
      IF (ml1.lt.ml2) THEN
        ptr2=INDEX(brch1,' ') - 1
        ptr2 = ptr1 + ptr2
        copy(ptr1:ptr1)   = '('
        copy(ptr1+1:ptr2) = brch1
        ptr1 = ptr2 + 1
        copy(ptr1:ptr1)=')'
        ptr1=ptr1+1

        ptr2=INDEX(brch2,' ') - 1
        ptr2 = ptr1 + ptr2
        copy(ptr1:ptr1)   = '('
        copy(ptr1+1:ptr2) = brch2
        ptr1 = ptr2 + 1
        copy(ptr1:ptr1)=')'
        ptr1=ptr1+1
        RETURN

      ELSE IF (ml2.lt.ml1) THEN
        ptr2=INDEX(brch2,' ') - 1
        ptr2 = ptr1 + ptr2
        copy(ptr1:ptr1)   = '('
        copy(ptr1+1:ptr2) = brch2
        ptr1 = ptr2 + 1
        copy(ptr1:ptr1)=')'
        ptr1=ptr1+1

        ptr2=INDEX(brch1,' ') - 1
        ptr2 = ptr1 + ptr2
        copy(ptr1:ptr1)   = '('
        copy(ptr1+1:ptr2) = brch1
        ptr1 = ptr2 + 1
        copy(ptr1:ptr1)=')'
        ptr1=ptr1+1
        RETURN

      ELSE
        IF (brch1==brch2) THEN
          maxpri=1
        ELSE
          tempbr(1)=brch1
          tempbr(2)=brch2
          DO i=3,SIZE(tempbr)
            tempbr(i)=' '
          ENDDO
          CALL brpri(group,rank,tempbr,2,nring,maxpri)
        ENDIF
        IF (maxpri==1) THEN
          ptr2=INDEX(brch1,' ') - 1
          ptr2 = ptr1 + ptr2
          copy(ptr1:ptr1)   = '('
          copy(ptr1+1:ptr2) = brch1
          ptr1 = ptr2 + 1
          copy(ptr1:ptr1)=')'
          ptr1=ptr1+1

          ptr2=INDEX(brch2,' ') - 1
          ptr2 = ptr1 + ptr2
          copy(ptr1:ptr1)   = '('
          copy(ptr1+1:ptr2) = brch2
          ptr1 = ptr2 + 1
          copy(ptr1:ptr1)=')'
          ptr1=ptr1+1
          RETURN

        ELSE
          ptr2=INDEX(brch2,' ') - 1
          ptr2 = ptr1 + ptr2
          copy(ptr1:ptr1)   = '('
          copy(ptr1+1:ptr2) = brch2
          ptr1 = ptr2 + 1
          copy(ptr1:ptr1)=')'
          ptr1=ptr1+1

          ptr2=INDEX(brch1,' ') - 1
          ptr2 = ptr1 + ptr2
          copy(ptr1:ptr1)   = '('
          copy(ptr1+1:ptr2) = brch1
          ptr1 = ptr2 + 1
          copy(ptr1:ptr1)=')'
          ptr1=ptr1+1
          RETURN
        ENDIF

      ENDIF
    ENDIF

  ENDIF

END SUBROUTINE mkcopy

!========================================================
! PURPOSE :  Add "=" to formula with double bonds          
!========================================================
SUBROUTINE dwrite(cpchem)
  IMPLICIT NONE

  CHARACTER(LEN=*), INTENT(INOUT) :: cpchem   ! chemical formula in which '=' is added

  CHARACTER(LEN=LEN(cpchem)) :: tempkc
  INTEGER  :: nc,idb,l,i,p,idp,ncd,j,ibflg,n,start,start2,lengr

! initialize
  l=0  ; idb=0  ; idp=0 ; p=0
  nc=INDEX(cpchem,' ')-1  ;  tempkc = ' '  ; lengr=LEN(cpchem)    

! write "=" and reset IDB if not CdCdCd:
  DO i=1,nc

! Find the end of the branch after Cd
    IF (idp/=0) THEN
      IF (cpchem(i:i)=='(')      THEN ; p = p+1
      ELSE IF (cpchem(i:i)==')') THEN ; p = p-1 ; ENDIF
      IF (p==1) start2 = 1
    ENDIF

    IF (cpchem(i:i+1)=='Cd') THEN
      IF (((cpchem(i+2:i+2)=='(').OR.(cpchem(i+3:i+3)=='(')).AND.(idb==0)) THEN

! check if there's a double bond in the branch.
        ncd=0  ;  ibflg=0  ;  n=0  ;  start=0
        DO j=i+2,nc
          IF (cpchem(j:j)=='(')      THEN ; n = n+1
          ELSE IF (cpchem(j:j)==')') THEN ; n = n-1 ; ENDIF
          IF (n==1) start = 1
          IF ((cpchem(j:j+1)=='Cd').AND.(ibflg==0)) THEN
             ncd = ncd+1
          ENDIF
          IF ((n==0).AND.(ibflg==0).AND.(start==1)) THEN
            ibflg=1
            IF (ncd==2) THEN
              idp = idp + 1
              ncd = 0
            ELSE
              idb = idb + 1
              ncd = 0
            ENDIF
          ENDIF
        ENDDO             
      ELSE
        idb = idb + 1
      ENDIF       

      IF ((idp==1).AND.(p==0).AND.(start2==1)) THEN
        l = l + 1
        tempkc(l:l) = '='
        idp = 0
        start2 = 0
        IF (idb==1) idb=0
      ELSE IF (idb > 1) THEN
        l = l + 1
        tempkc(l:l) = '='
        IF (cpchem(i+2:i+3)/='Cd') idb = 0
      ENDIF
    ENDIF

    l = l+1
    tempkc(l:l) = cpchem(i:i)
  ENDDO

! stop if formula is longer than lfo:
  IF (l > lengr) THEN
    WRITE(6,'(a)') '--error-- in dwrite. Chemical formula is too long:'
    WRITE(6,'(a)') TRIM(cpchem)
    STOP  "in dwrite"
  ENDIF

  cpchem = tempkc
END SUBROUTINE dwrite

! ======================================================================
! PURPOSE: Write the reverse of chemical formula in COPY to output
!   The chain in input formula is written from right to left.     
! ======================================================================
SUBROUTINE revers(copy,cc)
  IMPLICIT NONE

  CHARACTER(LEN=*),INTENT(IN)  :: copy  ! chemical formula
  CHARACTER(LEN=*),INTENT(OUT) :: cc    ! reversed chemical formula

  INTEGER :: begin,iend,p,i,lng

! initialize
! ----------
  p = 0  ;   cc = ' '       ;  begin = 1
  iend = INDEX(copy,' ')-1  ;  lng=iend

! reverse formula
! ---------------
  DO i=lng,1,-1
    IF (copy(i:i)==')') p = p+1
    IF (copy(i:i)=='(') p = p-1
    IF (p==0) THEN
      IF (copy(i:i)=='C' .OR. copy(i:i)=='c') THEN
        IF (copy(i:i+1)/='Cl') THEN
          cc(begin:begin+iend-i) = copy(i:iend)
          begin = begin + iend - i + 1
          iend = i - 1
        ENDIF
      ENDIF
      IF (copy(i:i+2)=='-O-') THEN
        cc(begin:begin+2)='-O-'
        begin=begin+3
        iend=i-1
      ELSE IF (copy(i:i+3)=='-O1-') THEN
        cc(begin:begin+3)='-O1-'
        begin=begin+4
        iend=i-1
      ELSE IF (copy(i:i+3)=='-O2-') THEN
        cc(begin:begin+3)='-O2-'
        begin=begin+4
        iend=i-1
      ELSE IF (copy(i:i+3)=='-O3-') THEN
        cc(begin:begin+3)='-O3-'
        begin=begin+4
        iend=i-1
      ENDIF
    ENDIF
  ENDDO
END SUBROUTINE revers

END MODULE stdtool
