MODULE switchisomer
IMPLICIT NONE
CONTAINS
!=======================================================================
! PURPOSE:  Replace a species by an already existing species in the 
! dictionary. The isomer must have the same number of carbons and
! functional groups. Other criteria, like number of COCO groups, are
! used to select which isomer (if any) best fit the structure of the 
! species provided as input.                                                   
!=======================================================================
SUBROUTINE isomer(chem,brtio,chg,tabinfo)
  USE dictstackdb, ONLY: nrec,dict,stabl,dbrch,mxcri,mxiso,diccri
  USE keyparameter, ONLY: mxlfo,mxlfl
  USE namingtool, ONLY: codefg
  USE toolbox, ONLY: stoperr
  IMPLICIT NONE      

  CHARACTER(LEN=*),INTENT(INOUT) :: chem   ! species (in) to be replaced (out)
  REAL, INTENT(IN)    :: brtio             ! yield of the species
  INTEGER,INTENT(OUT) :: chg               ! flag: isomer switch done (1) on not (0)
  INTEGER,INTENT(OUT) :: tabinfo(mxcri)    ! structural criteria of input species

  CHARACTER*(mxlfl) :: fgrp
  INTEGER         :: lfg
  INTEGER         :: i,j
  INTEGER         :: nlen, niso, tiso, imax, kc
  INTEGER         :: dictinfo(mxiso)
  INTEGER         :: dicflg(mxiso)
  INTEGER         :: score(mxiso)
  REAL            :: ymax

  CHARACTER(LEN=7)  :: progname='isomer'
  CHARACTER(LEN=70) :: mesg
      
  chg=0  ;  lfg=1  ;  fgrp=' '
  dicflg(:)=0  ;  score(:)=0  ; dictinfo(:)=0
  tabinfo(:)=0

! get info about the species provided as input
  CALL interacgrp(chem,tabinfo)
      
! return if generation is lesser than a given number
  IF (stabl<2) RETURN

! get the code string for chem
  CALL codefg(chem,fgrp,lfg)
  nlen=INDEX(fgrp,' ')

! ignore species bearing a >C=C< bond
!  IF (INDEX(fgrp,'U')/=0) RETURN

! ignore species having a yield greater than 1%
!c  IF (brtio>=0.01) RETURN
      
! -------------------------------------------------------------------
! LOOP OVER SPECIES IN THE DICTIONARY - KEEP A FIRST SET OF SPECIES
! -------------------------------------------------------------------
  niso=0
  recloop: DO i=1, nrec

! ignore species not having the same number of groups
    IF (tabinfo(1)/=diccri(i,1)) CYCLE recloop

! ignore the species not bearing the same list of groups
    DO j=1,nlen
      IF (fgrp(j:j)/=dict(i)(mxlfo+11+j:mxlfo+11+j)) CYCLE recloop
    ENDDO

! store the molecule that may be used as isomers (if that point is reached)
    niso=niso+1
    IF (niso>mxiso) THEN
      mesg="isomer # > mxiso"
      CALL stoperr(progname,mesg,chem)
    ENDIF
    dictinfo(niso)=i
  ENDDO recloop

! exit if no isomer available
  IF (niso==0) RETURN

! store the number of entry in the tables for isomers
  tiso=niso

  IF (brtio>=3E-2) RETURN

! the isomer must match critia 1 only. Criteria 2-16 are used to discriminate
  IF (brtio<5E-2) THEN

! ----------------------------------------------------------------
! Roll out the check list of criteria to select the 'best' isomere
! ----------------------------------------------------------------
!
! -1 check that the carbon skeletons are similar (criteria 6-9)
    IF (niso>1) THEN ! more than one isomer still remains
      DO kc=6,9
        CALL seliso(kc,tiso,dictinfo,tabinfo,dicflg,niso)
      ENDDO
    ENDIF

! -2 check the number of -CO-CO- conjugaisons (criteria 10)           
    IF (niso>1) THEN ! more than one isomer still remains
      CALL seliso(10,tiso,dictinfo,tabinfo,dicflg,niso)
    ENDIF

! -3 check the number of 1-2 interactions (criteria 11)           
    IF (niso>1) THEN ! more than one isomer still remains
      CALL seliso(11,tiso,dictinfo,tabinfo,dicflg,niso)
    ENDIF

! -4 substitution at terminal end (i.e. check CH3, criteria 2) 
    IF (niso>1) THEN ! more than one isomer still remains
      CALL seliso(2,tiso,dictinfo,tabinfo,dicflg,niso)
    ENDIF

! -5 alkyl chain (criteria 17-15) 
    IF (niso>1) THEN ! more than one isomer still remains
      DO kc=17,15,-1
        CALL seliso(kc,tiso,dictinfo,tabinfo,dicflg,niso)
      ENDDO
    ENDIF

! -6 number of -CO-CO- at a terminal end the chain (criteria 14)
    IF (niso>1) THEN ! more than one isomer still remains
      CALL seliso(14,tiso,dictinfo,tabinfo,dicflg,niso)
    ENDIF

! -7 number of 1-4 interactions (first) and 1-3 (then) (criteria 13-12)
    IF (niso>1) THEN ! more than one isomer still remains
      DO kc=13,12,-1
        CALL seliso(kc,tiso,dictinfo,tabinfo,dicflg,niso)
      ENDDO
    ENDIF

! -8 number of -CO-CHO, -CO-CO(OONO2), -CO-CO(OH), -CO-CO(OOH) (criteria 18-21)
    IF (niso>1) THEN ! more than one isomer still remains
      DO kc=18,21,1
        CALL seliso(kc,tiso,dictinfo,tabinfo,dicflg,niso)
      ENDDO
    ENDIF

!-9 number of CH2, CH, C (criteria 3-5)
    IF (niso>1) THEN ! more than one isomer still remains
      DO kc=3,5
        CALL seliso(kc,tiso,dictinfo,tabinfo,dicflg,niso)
      ENDDO
    ENDIF

  ENDIF   

! -------------------------------
! perform the isomer substitution
! -------------------------------

! Only one isomer remains
  IF (niso==1) THEN  
    DO i=1,tiso
      IF (dicflg(i)==0) THEN
        j=dictinfo(i)
        chem=dict(j)(10:mxlfo+10)
        dbrch(j)=dbrch(j)+brtio
        chg=1
        RETURN
      ENDIF
    ENDDO

! more than one molecule remains - pick the one with highest yield
  ELSE IF (niso>1) THEN  
    imax=0 ; ymax=0.
    DO i=1,tiso
      IF (dicflg(i)==0) THEN
        j=dictinfo(i)
        IF (dbrch(j)>ymax) THEN
          ymax=dbrch(j) ; imax=j
        ENDIF
      ENDIF
    ENDDO
    chem=dict(imax)(10:mxlfo+10)
    dbrch(imax)=dbrch(imax)+brtio
    chg=1
    RETURN
  ENDIF

END SUBROUTINE isomer

!=======================================================================
! Purpose : count the occurence of some structural parameter in the 
! molecule provided as input. Data are stored in "tabinfo" in which the 
! meaning of the indexes are:  
! tabinfo(1) : # of group in chem       
! tabinfo(2) : # of CH3  
! tabinfo(3) : # of CH2  
! tabinfo(4) : # of CH   
! tabinfo(5) : # of C    
! tabinfo(6) : # of primary node (ending position)        
! tabinfo(7) : # of secondary node                        
! tabinfo(8) : # of tertiary node (branching position)    
! tabinfo(9) : # of quaternary node (2 branching )        
! tabinfo(10): # of -CO-CO- conjugaisons                  
! tabinfo(11): # of 1-2 interactions                      
! tabinfo(12): # of 1-3 interactions                      
! tabinfo(13): # of 1-4 interactions                      
! tabinfo(14): # of -CO-CO- at a terminal end the chain   
! tabinfo(15): # of -CH2CH3         
! tabinfo(16): # of -CH2CH2CH3      
! tabinfo(17): # of -CH2CH2CH2CH3   
! tabinfo(18): # of -CO-CHO         
! tabinfo(19): # of -CO-CO(OONO2)   
! tabinfo(20): # of -CO-CO(OH)      
! tabinfo(21): # of -CO-CO(OOH)     
!=======================================================================
SUBROUTINE interacgrp(chem,tabinfo)
  USE keyparameter, ONLY: mxnode,mxlgr,mxring
  USE rjtool
  USE stdgrbond
  IMPLICIT NONE

  CHARACTER(LEN=*),INTENT(IN) :: chem       ! input species 
  INTEGER,INTENT(OUT)         :: tabinfo(:) ! # of each sub-structure in chem

  CHARACTER(LEN=mxlgr) :: group(mxnode)
  CHARACTER(LEN=1)     :: nodetype(mxnode)
  INTEGER :: ng, i, j, k, l
  INTEGER :: tabgrp(mxnode)
  INTEGER :: bond(mxnode,mxnode), dbflg, nring
  INTEGER :: memo1,memo2
  INTEGER :: rjg(mxring,2)
  INTEGER :: isum
      
  tabinfo(:)=0 ; tabgrp(:)=1 ; nodetype(:)=' '

! build the group and bond matrix for chem      
  CALL grbond(chem,group,bond,dbflg,nring)
  DO i=1,SIZE(bond,1)
    bond(i,i)=0
  ENDDO

! IF RINGS EXIST remove ring-join characters from groups,
  IF (nring>0)  CALL rjgrm(nring,group,rjg)

! count the number of group in chem
  ng=COUNT(group/=' ')  ;  tabinfo(1)=ng

! -----------------------------------------
! GET SOME INFO FOR EACH NODE              
! -----------------------------------------

! get the type of nodes in groups
  DO i=1,ng
    IF (group(i)(1:2)=='CO')       THEN ; nodetype(i)='y'
    ELSE IF (group(i)(1:3)=='CHO') THEN ; nodetype(i)='y'
    ElSE IF (group(i)(1:1)=='c')   THEN ; nodetype(i)='r'
    ELSE IF (group(i)(1:3)=='-O-') THEN ; nodetype(i)='o'
    ELSE                                ; nodetype(i)='n'
    ENDIF 
  ENDDO

! count the number of aliphatic groups (store info in cell 2-5)
! ---------------------------------------------------------------------
  DO i=1,ng
    IF (group(i)(1:4)=='CH3 ')      THEN ; tabgrp(i)=0 ; tabinfo(2)=tabinfo(2)+1
    ELSE IF (group(i)(1:4)=='CH2 ') THEN ; tabgrp(i)=0 ; tabinfo(3)=tabinfo(3)+1
    ELSE IF (group(i)(1:3)=='CH ')  THEN ; tabgrp(i)=0 ; tabinfo(4)=tabinfo(4)+1
    ELSE IF (group(i)(1:2)=='C ')   THEN ; tabgrp(i)=0 ; tabinfo(5)=tabinfo(5)+1
    ENDIF
  ENDDO

! count the branching (i.e whether prim, sec ...) (info in cell 6-9)
! ---------------------------------------------------------------------
  DO i=1,ng
    isum=0
    DO j=1,ng 
      IF (bond(i,j)/=0) isum=isum+1
    ENDDO
    IF (isum==1) tabinfo(6)=tabinfo(6)+1
    IF (isum==2) tabinfo(7)=tabinfo(7)+1
    IF (isum==3) tabinfo(8)=tabinfo(8)+1
    IF (isum==4) tabinfo(9)=tabinfo(9)+1
  ENDDO

! count the number of conjugated COCO groups (store info in cell 8-12)
! ---------------------------------------------------------------------
  DO i=1,ng-1
    IF (nodetype(i)=='y') THEN
      DO j=i+1,ng
        IF (bond(i,j)/=0) THEN
          IF (nodetype(j)=='y') THEN
           tabinfo(10)=tabinfo(10)+1
           IF (group(j)=='CHO ')       tabinfo(18)=tabinfo(18)+1
           IF (group(j)=='CO(OONO2) ') tabinfo(19)=tabinfo(19)+1
           IF (group(j)=='CO(OH) ')    tabinfo(20)=tabinfo(20)+1
           IF (group(j)=='CO(OOH) ')   tabinfo(21)=tabinfo(21)+1
          ENDIF
        ENDIF
      ENDDO
    ENDIF
  ENDDO
  tabinfo(14)=tabinfo(18)+tabinfo(19)+tabinfo(20)+tabinfo(21)

! 1: # group in chem            2: # CH3           
! 3: # CH2                      4: # CH            
! 5: # C                        6: # primary node (ending position)
! 7: # secondary node           8: # tertiary node (branching position)
! 9: # quaternary node         10: # -CO-CO- conjugaisons 
! 11: # 1-2 interactions       12: # 1-3 interactions     
! 13: # 1-4 interactions       14: # -CO-CO- at a terminal end the chain 
! 15: # -CH2CH3                16: # -CH2CH2CH3   
! 17: # -CH2CH2CH2CH3          18: # -CO-CHO      
! 19: # -CO-CO(OONO2)          20: # -CO-CO(OH)   
! 21: # -CO-CO(OOH)  

! count the length of non substitued chain
! ---------------------------------------------------------------------
  IF (INDEX(chem,'CH2CH3')/=0)       tabinfo(15)=tabinfo(15)+1       
  IF (INDEX(chem,'CH3CH2')/=0)       tabinfo(15)=tabinfo(15)+1       
  IF (INDEX(chem,'CH2CH2CH3')/=0)    tabinfo(16)=tabinfo(16)+1       
  IF (INDEX(chem,'CH3CH2CH2')/=0)    tabinfo(16)=tabinfo(16)+1       
  IF (INDEX(chem,'CH2CH2CH2CH3')/=0) tabinfo(17)=tabinfo(17)+1    
  IF (INDEX(chem,'CH3CH2CH2CH2')/=0) tabinfo(17)=tabinfo(17)+1

! count the number of 1-2, 1-3 and 1-4 interaction
! (store info in cell 5 to 7)
! --------------------------------------------------
  DO i=1,ng
    memo1=0 ; memo2=0
    IF (tabgrp(i)/=0) THEN
      DO j=1,ng
        IF (bond(i,j)/=0) THEN
          IF (tabgrp(j)/=0) tabinfo(11)=tabinfo(11)+1   ! 1-2 interaction
  
          DO k=1,ng    ! loop for 1-3 interactions
            IF ((bond(j,k)/=0) .AND. (k/=i)) THEN
              IF (tabgrp(k)/=0) THEN
                 IF ((group(i)(1:3)=='-O-').OR.(group(k)(1:3)=='-O-')) THEN 
                   tabinfo(11)=tabinfo(11)+1
                 ELSE
                   tabinfo(12)=tabinfo(12)+1
                 ENDIF
               ENDIF
  
               DO l=1,ng  ! loop for 1-4 interactions
                 IF ((bond(k,l)/=0) .AND. (l/=j)) THEN
                   IF (tabgrp(l)/=0) THEN
                     IF ((group(i)(1:3)=='-O-') .OR. (group(l)(1:3)=='-O-')) THEN 
                         tabinfo(12)=tabinfo(12)+1
                     ELSE
                       IF ((i/=memo1) .AND. (l/=memo2)) THEN !for ring
                         tabinfo(13)=tabinfo(13)+1
                         memo1=i ;  memo2=l
                       ENDIF
                     ENDIF
                   ENDIF
                 ENDIF
               ENDDO
  
             ENDIF
           ENDDO
        ENDIF
      ENDDO
      tabgrp(i)=0
    ENDIF
  ENDDO

END SUBROUTINE interacgrp

!=======================================================================
! Purpose : remove the isomers that does not match the criteria (kc)
! povided as input to discriminate among the list of isomers. If more 
! than one isomer remains, then the criteria is retained. In no isomer
! remains, then return the list of possible isomers unchanged.
!=======================================================================
SUBROUTINE seliso(kc,tiso,dictinfo,tabinfo,dicflg,niso)
  USE dictstackdb, ONLY: diccri
  IMPLICIT NONE      

  INTEGER,INTENT(IN) :: tabinfo(:)    ! # of each sub-structure in chem
  INTEGER,INTENT(IN) :: dictinfo(:)   ! isomers of chem
  INTEGER,INTENT(IN) :: tiso          ! total # of isomers (entries) in dictinfo
  INTEGER,INTENT(IN) :: kc            ! ID # of the criteria used to search/discriminate isomers
      
  INTEGER,INTENT(INOUT) :: dicflg(:)  ! flag set to 0 if isomer in dictinfo still OK
  INTEGER,INTENT(INOUT) :: niso       ! # of isomers left for substitution 

  INTEGER :: i,ij
  INTEGER :: tempniso
  INTEGER :: tempdicflg(SIZE(dicflg))

! exit if only one isomer remains
  IF (niso==1) RETURN  

! copy some value that may need to be returned without changes
  tempniso=niso
  tempdicflg(1:tiso)=dicflg(1:tiso)

! unselect the species that does not match the input criteria
  DO i=1,tiso
    IF (dicflg(i)/=1) THEN  ! ignore species already eliminated
      ij=dictinfo(i)
      IF (tabinfo(kc)/=diccri(ij,kc)) THEN
        dicflg(i)=1
        niso=niso-1
      ENDIF
    ENDIF
  ENDDO  

! no isomer remains (no species match the criteria). Return without changes
  IF (niso==0) THEN
    niso=tempniso  ;  dicflg(1:tiso)=tempdicflg(1:tiso)
    RETURN
  ENDIF
 
END SUBROUTINE seliso

END MODULE switchisomer
