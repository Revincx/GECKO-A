MODULE normchem
  IMPLICIT NONE
  CONTAINS

!=======================================================================
! The purpose of this subroutine is to avoid duplicates in the list of 
! species. Therefore, chemical formula must be written in a unique way.                                                 
!                                                                  
! First, the functional groups of the molecule are written in an array, 
! GROUP, the bond-matrix is evaluated and the double-bond flag, DBFLG, 
! is set, if there is a double-bond in CHEM. Then, the very longest tree 
! in CHEM is evaluated after looking top-down and down-top. Each 
! possible formula with the longest possible chain in the chain is a new
! copy. All non-identical formulas for the molecule undergo a final check.
! The functional groups, which are attached to the chain or a carbon 
! must have a predefined order. After all this rewriting and these 
! checks there is only one possible formula for the molecule left.                          
!=======================================================================
SUBROUTINE stdchm(chem)
  USE keyparameter, ONLY:mxnode,mxlgr,mxring,mxcp ! no input/output table for chem.
  USE keyflag, ONLY: sar_only_fg
  USE atomtool, ONLY: cnum, onum
  USE stdtool, ONLY: ckgrppt,lntree,mkcopy,revers,prioty,dwrite
  USE stdgrbond, ONLY: grbond
  USE stdratings, ONLY: ratings
  USE ringtool, ONLY: uniqring
  USE toolbox, ONLY: stoperr
  IMPLICIT NONE

  CHARACTER(LEN=*),INTENT(INOUT) :: chem

! local
  INTEGER         :: dbflg
  INTEGER         :: bond(mxnode,mxnode)
  INTEGER         :: path(mxnode,mxnode,mxnode) 
  INTEGER         :: clngth(mxnode,mxnode)
  INTEGER         :: rjg(mxring,2)  
  CHARACTER(LEN=mxlgr) :: group(mxnode)
  LOGICAL         :: lobond(mxnode,mxnode)
  INTEGER         :: zebond(mxnode,mxnode)  
  LOGICAL         :: loze

  CHARACTER(LEN=LEN(chem)) :: tchem, copy(mxcp)
  INTEGER         :: nc, nca, nce, nl, nring 
  INTEGER         :: leaf, last, ptr, p, i, j, k, ncp
  INTEGER         :: ncis,ntrans
  INTEGER         :: ig, pg, ng
  INTEGER         :: locat
  INTEGER         :: rank(mxnode)

  CHARACTER(LEN=8),PARAMETER :: progname='stdchm '
  CHARACTER(LEN=70)          :: mesg
  
  IF (sar_only_fg) RETURN  ! don't normalize if sar_only_fg on (not needed)

! ------------------------------------------------------
! INITIALIZE
! ------------------------------------------------------
  copy(:)=' ' ; clngth(:,:)=0 ; lobond(:,:)=.FALSE. ; path(:,:,:)=0
  loze=.FALSE. 
  IF (INDEX(chem,'/')/=0) loze=.TRUE. 
  IF (INDEX(chem,'\')/=0) loze=.TRUE.
     
! ------------------------------------------------------
! CHECK INPUT FORMULA AND MAKE BOND AND GROUP TABLES
! ------------------------------------------------------

! get the number of >C< , >c< (nca) and -O- (nce) in the molecule
  nca = cnum(chem)  
  nce = onum(chem)
  nca=nca+nce  ! C and -O- into nca

! only multi-carbon molecules checked:
  IF (nca < 2) RETURN

! check parenthesis:if open parenthesis,then stop the run.
  nc  = INDEX(chem,' ')-1  ;  p = 0
  DO i=1,nc
    IF(chem(i:i)=='(') p = p+1
    IF(chem(i:i)==')') p = p-1
  ENDDO

  IF (p/=0) THEN
    mesg="parentheses mismatch (error 1)"
    CALL stoperr(progname,mesg,chem)
  ENDIF
     
! build the bond and group matrix
  CALL grbond(chem,group,bond,dbflg,nring,zebond)

! tag couple of groups belonging to a cis/trans C=C bonds 
  ncis=0 ; ntrans=0
  IF (loze) CALL add_zetag(chem,nca,zebond,group,ncis,ntrans)

! rank groups according to rating
  CALL ratings(nca,group,bond,nring,rank)

! If rings: (1) rearrange according to rank, (2) zero out 'broken' ring bonds
  IF (nring > 0) THEN
    CALL uniqring(nring,nca,group,bond,rank,rjg)
    DO k=1,nring
      i=rjg(k,1)  ;  j=rjg(k,2)
      bond(i,j)=0 ;  bond(j,i)=0
    ENDDO
  ENDIF

! ------------------------------------------------------
! WRITE ALL POSSIBLE FORMULA OF THE SPECIES
! ------------------------------------------------------
! make a logical copy of the bond matrix
  DO i=1,nca
    DO j=1,nca
      IF (bond(i,j)/=0) lobond(i,j)= .TRUE.
    ENDDO
  ENDDO

! check that the functional groups are correctly sorted for each group
  DO i=1,nca
    locat=INDEX(group(i),')(')
    IF (locat /= 0) CALL ckgrppt(locat,group(i))
  ENDDO

! find longest tree, top-down, starting with the first group
  CALL lntree(bond,1,2,nca,clngth,path)

! look down-top for the very longest tree ...
  ncp = 0
  DO i=1,nca
    IF (clngth(1,i)/=0) THEN
      leaf = path(1,i,clngth(1,i))
      last = path(1,i,clngth(1,i)-1)
      CALL lntree(bond,leaf,last,nca,clngth,path)

! write the formula of species according to the longest tree 
! in path. Each group are written in the formula with all the
! branches by the mkcopy subroutine. The formula is written twice
! for each longest tree, down-top and top-down by the revers subroutine.  
      DO j=1,nca
        IF (clngth(leaf,j)/=0) THEN
           ncp=ncp+1
           ptr = 1
           DO k=1,nca
             ig = path(leaf,j,k)
             IF (ig/=0) THEN
                IF (k > 1)   pg = path(leaf,j,k-1)
                IF (k < nca) ng = path(leaf,j,k+1)
                CALL mkcopy(lobond,group,nca,rank,nring,ig,pg,ng,ptr,copy(ncp))
             ENDIF
           ENDDO
           CALL revers(copy(ncp),copy(ncp+1))
           ncp = ncp+1
        ENDIF
      ENDDO
    ENDIF
  ENDDO

! ------------------------------------------------------
! COLLAPSE IDENTICAL FORMULA
! ------------------------------------------------------
  IF (ncp > 1) THEN
    nl = ncp
    rmloop: DO i=1,nl-1
      j=i
      DO ! exit from rmloop
        j=j+1
        IF (j > nl)  CYCLE rmloop
        IF (copy(i)==copy(j)) THEN
          copy(j) = ' '
          DO k=j,nl
            copy(k) = copy(k+1)
          ENDDO
          IF(copy(i)/=' ') THEN
            ncp = ncp - 1
            j = j - 1
          ENDIF
        ENDIF
      ENDDO
    ENDDO rmloop     
  ENDIF

  IF (ncp==0) THEN
    mesg="no copies left "
    CALL stoperr(progname,mesg,chem)
  ENDIF

! ------------------------------------------------------
! if more than one copy, get the standardized formula
! ------------------------------------------------------

! subroutine prioty check group order and return the standardized formula
  IF (ncp==1) THEN
    tchem = copy(1)
  ELSE
    CALL prioty(group,rank,copy,ncp,nring,tchem)
  ENDIF

! write double bond
  IF (dbflg/=0) CALL dwrite(tchem)

! check parenthesis of the standardized copy
  nc = INDEX(tchem,' ')-1  ;  p = 0
  DO i=1,nc
    IF (tchem(i:i)=='(') p = p+1
    IF (tchem(i:i)==')') p = p-1
  ENDDO
  IF (p/=0) THEN
    mesg="parentheses mismatch (error 3) "
    CALL stoperr(progname,mesg,chem)
  ENDIF

! remove tag for couple of nodes in C=C bonds
  IF (loze) CALL rm_zetag(tchem,ncis,ntrans)

! ------------------------------------------------------
! RETURN THE STANDARDIZED FORMULA 
! ------------------------------------------------------
  chem = tchem

END SUBROUTINE stdchm

!=======================================================================
! Purpose: tag the couple of groups belonging to the same C=C with cis  
! or trans information. Use a, b or c character for cis configuration  
! and x, y or z for trans configuration. These characters are aimed to 
! be removed (see rm_zetag subroutine) once the formula has been 
! standardized.
!=======================================================================
SUBROUTINE add_zetag(chem,nca,zebond,group,ncis,ntrans)
  USE toolbox, ONLY: stoperr
  IMPLICIT NONE
  
  CHARACTER(LEN=*),INTENT(IN)   :: chem        ! formula of current species
  INTEGER,INTENT(IN)            :: nca         ! number of nodes
  INTEGER,INTENT(IN)            :: zebond(:,:) ! cis (1) or trans (2) [i,j] bond 
  CHARACTER(LEN=*),INTENT(INOUT):: group(:)    ! groups in chem, to be modified with tag
  INTEGER,INTENT(OUT)           :: ncis        ! number of cis C=C bond in chem
  INTEGER,INTENT(OUT)           :: ntrans      ! number of trans C=C bond in chem

  INTEGER                       :: i,j,k,lgr
  CHARACTER(LEN=1)              :: zechar
  CHARACTER(LEN=12),PARAMETER   :: progname='add_zetag '
  CHARACTER(LEN=70)             :: mesg
  
  lgr=LEN(group(1)) ; ncis=0 ; ntrans=0

! scroll the zebond table
  DO i=1,nca-1
    DO j=i,nca
      IF (zebond(i,j)==0) CYCLE      
  
      ! cis C=C bond, get the tag to be added
      IF (zebond(i,j)==1) THEN
        ncis=ncis+1
        IF ( (group(i)(1:2)/='Cd').AND.(group(j)(1:2)/='Cd') ) THEN
          mesg="cis/trans info provided for group without 'Cd'"
          CALL stoperr(progname,mesg,chem)
        ENDIF
  
        IF      (ncis==1) THEN ; zechar='a'
        ELSE IF (ncis==2) THEN ; zechar='b'
        ELSE IF (ncis==3) THEN ; zechar='c'
        ELSE
          mesg="too many of cis C=C bond"
          CALL stoperr(progname,mesg,chem)
        ENDIF
 
      ! trans C=C bond, get the tag to be added
      ELSE IF (zebond(i,j)==2) THEN
        ntrans=ntrans+1
        IF ( (group(i)(1:2)/='Cd').AND.(group(j)(1:2)/='Cd') ) THEN
          mesg="cis/trans info provided for group without 'Cd'"
          CALL stoperr(progname,mesg,chem)
        ENDIF
  
        IF      (ntrans==1) THEN ; zechar='x'
        ELSE IF (ntrans==2) THEN ; zechar='y'
        ELSE IF (ntrans==3) THEN ; zechar='z'
        ELSE
          mesg="too many of trans C=C bond"
          CALL stoperr(progname,mesg,chem)
        ENDIF
  
      ELSE
        mesg="cis/trans bond number is out of range (1,2)"
        CALL stoperr(progname,mesg,chem)
      ENDIF
  
      ! add the tag in the pair of groups (i,j) 
      DO k=3,lgr     ! find slot
        IF ((group(i)/='1').AND.(group(i)/='2').AND.(group(i)/='3')) EXIT
      ENDDO 
      group(i)(k+1:lgr)=group(i)(k:lgr-1)
      group(i)(k:k)=zechar
      DO k=3,lgr
        IF ((group(j)/='1').AND.(group(j)/='2').AND.(group(j)/='3')) EXIT
      ENDDO 
      group(j)(k+1:lgr)=group(j)(k:lgr-1)
      group(j)(k:k)=zechar
  
    ENDDO
  ENDDO
END SUBROUTINE add_zetag
!=======================================================================
! Purpose: substitute the pair of tags belonging to the same C=C with  
! usual '/' and '\' characters. The formula provided as input must be
! standardized. First character to be added for each pair is '/' with 
! second character being either '\' (cis) or '/' (trans).
!=======================================================================
SUBROUTINE rm_zetag(chem,ncis,ntrans)
  USE toolbox, ONLY: stoperr
  IMPLICIT NONE
  
  CHARACTER(LEN=*),INTENT(INOUT):: chem
  INTEGER,INTENT(IN)            :: ncis
  INTEGER,INTENT(IN)            :: ntrans

  INTEGER                       :: i,k
  CHARACTER(LEN=1)              :: zechar
  CHARACTER(LEN=12),PARAMETER   :: progname='rm_zetag '
  CHARACTER(LEN=70)             :: mesg
  
  IF (ncis/=0) THEN
    DO i=1,ncis
      IF      (i==1) THEN ; zechar='a'
      ELSE IF (i==2) THEN ; zechar='b'
      ELSE IF (i==3) THEN ; zechar='c'
      ENDIF

      k=INDEX(chem,zechar)
      IF (k==0) THEN
        mesg="expected 1st cis/trans character not identified"
        CALL stoperr(progname,mesg,chem)
      ELSE
        chem(k:k)='/'
      ENDIF

      k=INDEX(chem,zechar)
      IF (k==0) THEN
        mesg="expected 2nd cis/trans character not identified"
        CALL stoperr(progname,mesg,chem)
      ELSE
        chem(k:k)='\'
      ENDIF

    ENDDO
  ENDIF

  IF (ntrans/=0) THEN
    DO i=1,ntrans
      IF      (i==1) THEN ; zechar='x'
      ELSE IF (i==2) THEN ; zechar='y'
      ELSE IF (i==3) THEN ; zechar='z'
      ENDIF

      k=INDEX(chem,zechar)
      IF (k==0) THEN
        mesg="expected 1st cis/trans character not identified"
        CALL stoperr(progname,mesg,chem)
      ELSE
        chem(k:k)='/'
      ENDIF

      k=INDEX(chem,zechar)
      IF (k==0) THEN
        mesg="expected 2nd cis/trans character not identified"
        CALL stoperr(progname,mesg,chem)
      ELSE
        chem(k:k)='/'
      ENDIF

    ENDDO
  ENDIF

END SUBROUTINE rm_zetag



END MODULE normchem
  
