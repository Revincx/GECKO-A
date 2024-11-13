MODULE cdtool
IMPLICIT NONE
CONTAINS
!=======================================================================
! PURPOSE : check whether the C=C bonds are conjugated or not, including
! C=O bond (i.e. structure of type -C=C-C=O). The subroutine returns the
! "case" to which the species belongs. Six cases are considered :            
! CASE 1: regular Cd molecule : only >C=C< and >C=C-C=C< bonds in     
!         the molecule (i.e, without conjugated C=C-C=O)       
! CASE 2: is for structures containing the >C=C-C=O structure  
!         but no C=C-C=C (or C=C=O)      
! CASE 3: is for the -CO-C=C-C=C-C=O structure only (i.e. containing  
!         carbonyl at both sides of the conjugated C=C-C=C)     
! CASE 4: Two double bonds non conjugated (i.e. not C=C-C=C) but
!         containing at least one C=C-C=O                       
! CASE 5: is for the -CO-C=C-C=C< structure (i.e. containing    
!         carbonyl at only one side of the conjugated C=C-C=C)  
! CASE 6: -C=C=O e.g. ketene (only case for this group)         
! CASE 7: -C=C-O- => vinyl ether chemistry      
!=======================================================================
SUBROUTINE cdcase2(chem,bond,group,rxnflg,ncdtrack,cdtrack,&
                   xcdconjug,xcdsub,xcdeth,xcdcarbo,xcdcase)
  USE keyparameter, ONLY: mxcp           
  USE mapping, ONLY: alkenetrack
  USE toolbox, ONLY: stoperr
  IMPLICIT NONE

  CHARACTER(LEN=*),INTENT(IN) :: chem      ! formula
  INTEGER,INTENT(IN)          :: bond(:,:) ! bond matrix
  CHARACTER(LEN=*),INTENT(IN) :: group(:)  ! grp matrix
  INTEGER,INTENT(IN)          :: rxnflg    ! flg=1 if cdcheck required (e.g. for enols)

  INTEGER,INTENT(OUT) :: ncdtrack          ! # of Cd tracks 
  INTEGER,INTENT(OUT) :: cdtrack(:,:)      ! cdtrack(i,j) Cd nodes for the ith track         
  INTEGER,INTENT(OUT) :: xcdconjug(:)      ! flag for conjugated tracks (C=C-C=C)
  INTEGER,INTENT(OUT) :: xcdsub(:)         ! # of nodes bonded to a Cd group
  INTEGER,INTENT(OUT) :: xcdeth(:,:)       ! -OR group # bounded to a Cd (max 2 -OR groups)
  INTEGER,INTENT(OUT) :: xcdcarbo(:,:)     ! -CO- group # bounded to a Cd (max 2 CO groups)  
  INTEGER,INTENT(OUT) :: xcdcase(:)        ! Cd "case" of the track

  INTEGER :: cdtracklen(SIZE(cdtrack,1))   ! length of the ith track         
  INTEGER :: i,j,k
  INTEGER :: ngr,mxcd
  INTEGER :: cdbond(SIZE(group),SIZE(group))
  INTEGER :: icd,icd1st,icd2nd,nco

  CHARACTER(LEN=8)  :: progname='cdcase2 '
  CHARACTER(LEN=70) :: mesg

  ncdtrack=0    ; cdtrack(:,:)=0  ; xcdconjug(:)=0 ; xcdsub(:)=0 
  xcdeth(:,:)=0 ; xcdcarbo(:,:)=0 ; xcdcase(:)=0

  mxcd=SIZE(cdtrack,2)
  cdbond(:,:)=0
  
! get the number of groups
  ngr=COUNT(group/=' ')

! ------------------------------
! Create and check the Cd tracks
! ------------------------------
  CALL alkenetrack(chem,bond,group,ngr,ncdtrack,cdtracklen,cdtrack)
  DO i=1,ncdtrack
    IF (cdtracklen(i)==4) xcdconjug(i)=1
  ENDDO

! check for >C=CR-OH (enol: should have already tautomerised to ketone 
! in subroutine alkcheck) and >C=CR-ONO2 (not available yet)
  IF (rxnflg/=0) THEN
    DO i=1,ncdtrack ; DO j=1,mxcd
      IF (cdtrack(i,j)/=0) THEN
        IF (INDEX(group(cdtrack(i,j)),'(ONO2)')/=0) THEN
          mesg=">C=C(ONO2)- not allowed"
          CALL stoperr(progname,mesg,chem)
        ENDIF
        IF (INDEX(group(cdtrack(i,j)),'(OH)')/=0)THEN
          mesg=">C=C(OH)- not allowed"
          CALL stoperr(progname,mesg,chem)
        ENDIF
      ENDIF
    ENDDO ; ENDDO           
  ENDIF 

! ------------------------------
! Fill the various matrix
! ------------------------------

! count # of nodes bonded to each Cd, except the C=C node (for which bond(i,j)=2)
! and store -O- and CO nodes bonded to Cds
  DO i=1,ncdtrack ; DO j=1,mxcd
    IF (cdtrack(i,j)/=0) THEN
      icd=cdtrack(i,j)                                      ! current Cd node is icd 
      DO k=1,ngr                                            ! check Cd neighbours 
        IF (bond(icd,k)==1) THEN
          xcdsub(icd)=xcdsub(icd)+1                                ! count sub on Cd grp.
          IF (group(k)=='CHO' .OR. group(k)(1:2)=='CO') THEN
            IF      (xcdcarbo(icd,1)==0) THEN ; xcdcarbo(icd,1)=k  ! store CO node
            ELSE IF (xcdcarbo(icd,2)==0) THEN ; xcdcarbo(icd,2)=k  ! store CO node
            ELSE
              mesg="more than 2 carbonyls bonded to a Cd"
              CALL stoperr(progname,mesg,chem)
            ENDIF
          ENDIF
        ENDIF
        
        IF (bond(icd,k)==3) THEN
          xcdsub(icd)=xcdsub(icd)+1                         ! count sub on Cd grp.
          IF      (xcdeth(icd,1)==0) THEN ; xcdeth(icd,1)=k ! store -O- node
          ELSE IF (xcdeth(icd,2)==0) THEN ; xcdeth(icd,2)=k ! store -O- node
          ELSE
            mesg="more than 2 ethers bonded to a Cd"
            CALL stoperr(progname,mesg,chem)
          ENDIF
        ENDIF  
      ENDDO
    ENDIF
  ENDDO ; ENDDO

! -------------------------------
! set "cd case" for each Cd track 
! -------------------------------
  trloop: DO i=1,ncdtrack 

! check ketene >C=C=O (case 6)
    DO j=1,mxcd
      icd=cdtrack(i,j) ; IF (icd==0) EXIT
      IF (group(icd)=='CdO ') THEN
        xcdcase(i)=6 ; CYCLE trloop                     ! case 6 CdO 
      ENDIF
    ENDDO

! check vinylic structure C=C-OR (case 7)
    DO j=1,mxcd,2
      icd1st=cdtrack(i,j) ; IF (icd1st==0) EXIT
      icd2nd=cdtrack(i,j+1)
      IF (xcdeth(icd1st,1)/=0) xcdcase(i)=7
      IF (xcdeth(icd2nd,1)/=0) xcdcase(i)=7
      IF ( (xcdeth(icd1st,1)/=0).AND.(xcdeth(icd2nd,1)/=0) ) THEN
        mesg="ether found both sides of a Cd=Cd"
        CALL stoperr(progname,mesg,chem)
      ENDIF
      IF (xcdcase(i)/=0) CYCLE trloop                    ! case 7 C=C-OR
    ENDDO
    
! check conjugated C=C-C=C tracks (case 1, 3 or 5)
    IF (xcdconjug(i)/=0) THEN  
      nco=0
      DO j=1,4,3  ! check terminal Cd only for -CO-
        icd=cdtrack(i,j) 
        IF (icd==0) STOP '--error-- in cdcase, no Cd unexpected' 
        IF (xcdcarbo(icd,1)/=0) nco=nco+1
      ENDDO
      IF      (nco >1) THEN ; xcdcase(i)=3 ; CYCLE trloop ! case=3 - at least 2 carbonyls on distinct Cds 
      ELSE IF (nco==1) THEN ; xcdcase(i)=5 ; CYCLE trloop ! case=5 - only one carbonyl on conjugated Cds 
      ELSE                  ; xcdcase(i)=1 ; CYCLE trloop ! case=1 - regular conjugated diene
      ENDIF

! check simple C=C tracks (case 1 or 2)
    ELSE
      nco=0
      DO j=1,2
        icd=cdtrack(i,j) 
        IF (icd==0) STOP '--error-- in cdcase, no Cd unexpected' 
        IF (xcdcarbo(icd,1)/=0) nco=nco+1
      ENDDO
      IF (nco/=0) THEN ; xcdcase(i)=2 ; CYCLE trloop      ! case=2 - at least 1 carbonyl  
      ELSE             ; xcdcase(i)=1 ; CYCLE trloop      ! case=1 - regular simple alkene
      ENDIF

    ENDIF
  ENDDO trloop

END SUBROUTINE cdcase2

!=======================================================================
! PURPOSE : switch enol structure to the corresponding keto structure.
! NOTE: the routine is similar to alkcheck below, but limited to enol
! only.
!=======================================================================
SUBROUTINE switchenol(pchem,loswitch)
  USE keyparameter, ONLY: mxnode,mxlgr,mxring
  USE rjtool
  USE stdgrbond
  USE reactool
  USE normchem
  IMPLICIT NONE
  
  CHARACTER(LEN=*),INTENT(INOUT) :: pchem
  LOGICAL,INTENT(OUT) :: loswitch ! turned true if switch was performed 

  INTEGER            ::  bond(mxnode,mxnode),dbflg,nring
  CHARACTER(LEN=mxlgr)      ::  group(mxnode),pold,pnew,tgroup(mxnode)
  CHARACTER(LEN=LEN(pchem)) ::  tempkc
  INTEGER            ::  i,j
  INTEGER            ::  rjg(mxring,2) ! ring-join group pairs

  loswitch=.FALSE.

  IF (INDEX(pchem,'.')/=0) RETURN

  CALL grbond(pchem,group,bond,dbflg,nring)
  CALL rjgrm(nring,group,rjg)

  grloop: DO i=1,mxnode
    IF (INDEX(group(i),'Cd')/=0) THEN
      
      IF (INDEX(group(i),'(OH)')/=0) THEN
        pold='(OH)'  ;  pnew='O'
        CALL swap(group(i),pold,tgroup(i),pnew)
        group(i)=tgroup(i)
        pold='Cd'  ;  pnew='C'
        CALL swap(group(i),pold,tgroup(i),pnew)
        group(i)=tgroup(i)
        DO j=1,mxnode
          IF (bond(i,j)==2) THEN
            bond(i,j)=1  ;  bond(j,i)=1
            IF (group(j)(1:4)=='CdH2') THEN
              pold='CdH2'  ;  pnew='CH3'
            ELSE IF (group(j)(1:3)=='CdH') THEN        
              pold='CdH'  ;  pnew='CH2'
            ELSE IF (group(j)(1:2)=='Cd') THEN        
              pold='Cd'  ;  pnew='CH'
            ENDIF        
          CALL swap(group(j),pold,tgroup(j),pnew)
          group(j)=tgroup(j)
          CALL rebond(bond,group,tempkc,nring)
          pchem=tempkc
          CALL stdchm(pchem)
          ENDIF
        ENDDO  
        loswitch=.TRUE.
        !EXIT grloop

      ENDIF
    ENDIF
  ENDDO  grloop

END SUBROUTINE switchenol

!=======================================================================
!=======================================================================
! PURPOSE : This subroutine fragments the substitued alkenes such as
! >Cd=Cd(OH)- or >Cd=Cd(OOH)- or >Cd=Cd(ONO2)- which come from Norrish
! II (alkenes photolysis), or from fragmentation after oxidation of
! conjugated alkenes.
! We consider that these alkenes are energy-rich and decompose to :
!  >Cd=Cd(OH)-    -> >CH-CO- 
!  >Cd=Cd(OOH)-   -> >C(.)-CO- + OH.
!  >Cd=Cd(ONO2)-  -> >C(.)-CO- + NO2
! Each Cd can't support more than one group like OH, OOH or ONO2
! If there is more than one group on the double bond, treat first the
! -OH next the -OOH and at last the -ONO2
! Note: the routine is similar the switchenol routine above. This 
! alkcheck is called after adding the enol to the dictionary and 
! often act as a "patch" for treating enol structure. Must progressively
! be avoid with the progressive development of gecko.
!=======================================================================
SUBROUTINE alkcheck(pchem,coprod,acom)
  USE keyparameter, ONLY: mxnode,mxlgr,mxring
  USE rjtool
  USE stdgrbond
  USE reactool
  USE normchem
  IMPLICIT NONE
  
  CHARACTER(LEN=*),INTENT(INOUT) :: pchem
  CHARACTER(LEN=*),INTENT(OUT)   :: coprod
  CHARACTER(LEN=*),INTENT(OUT) :: acom ! code 

  INTEGER            ::  bond(mxnode,mxnode),dbflg,nring
  CHARACTER(LEN=mxlgr)      ::  group(mxnode),pold,pnew,tgroup(mxnode)
  CHARACTER(LEN=LEN(pchem)) ::  tempkc
  INTEGER            ::  i,j
  INTEGER            ::  rjg(mxring,2) ! ring-join group pairs

  coprod=' ' ; acom=' '

  IF (INDEX(pchem,'.')/=0) RETURN

  CALL grbond(pchem,group,bond,dbflg,nring)
  CALL rjgrm(nring,group,rjg)

  grloop: DO i=1,mxnode
    IF (INDEX(group(i),'Cd')/=0) THEN
      
      IF (INDEX(group(i),'(OH)')/=0) THEN
        pold='(OH)'  ;  pnew='O'
        CALL swap(group(i),pold,tgroup(i),pnew)
        group(i)=tgroup(i)
        pold='Cd'  ;  pnew='C'
        CALL swap(group(i),pold,tgroup(i),pnew)
        group(i)=tgroup(i)
        DO j=1,mxnode
          IF (bond(i,j)==2) THEN
            bond(i,j)=1  ;  bond(j,i)=1
            IF (group(j)(1:4)=='CdH2') THEN
              pold='CdH2'  ;  pnew='CH3'
            ELSE IF (group(j)(1:3)=='CdH') THEN        
              pold='CdH'  ;  pnew='CH2'
            ELSE IF (group(j)(1:2)=='Cd') THEN        
              pold='Cd'  ;  pnew='CH'
            ENDIF        
          CALL swap(group(j),pold,tgroup(j),pnew)
          group(j)=tgroup(j)
          CALL rebond(bond,group,tempkc,nring)
          pchem=tempkc
          CALL stdchm(pchem)
          ENDIF
        ENDDO  
        acom='KETOENOL ' 
        EXIT grloop

!      ELSE IF (INDEX(group(i),'(OOH)')/=0) THEN
!        pold='(OOH)'  ;  pnew='(O.)'
!        CALL swap(group(i),pold,tgroup(i),pnew)
!        group(i)=tgroup(i)
!        CALL rebond(bond,group,tempkc,nring)
!        pchem=tempkc
!        CALL stdchm (pchem)
!        coprod='HO'
!        acom='xxxxxxxx'
!        EXIT grloop

      ELSE IF (INDEX(group(i),'(ONO2)')/=0) THEN
        pold='(ONO2)'  ;  pnew='(O.)'
        CALL swap(group(i),pold,tgroup(i),pnew)
        group(i)=tgroup(i)
        CALL rebond(bond,group,tempkc,nring)
        pchem=tempkc
        CALL stdchm(pchem)
        coprod='NO2'
        acom='KETOENENIT' 
        EXIT grloop  
      ENDIF
    ENDIF
  ENDDO  grloop

END SUBROUTINE alkcheck

END MODULE cdtool
