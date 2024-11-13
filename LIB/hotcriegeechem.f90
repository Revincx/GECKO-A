MODULE hotcriegeechem
USE toolbox, ONLY: stoperr,add1tonp,addref
IMPLICIT NONE

CONTAINS

!=======================================================================
! PURPOSE: perform the reactions (decomposition, stabilization ...) of       
! hot criegee produced by O3+alkene reactions. The routine add new 
! products to the list of products of the current rxn (tables p(:), 
! short name) and related stoi. coef. (table s(:)). The dictionary, 
! stack and related tables are updated accordingly.
! Note : zconform is the ratio of the Z conformer for the CI* provided 
! as input. The value make no sense when the 2 branches attached to the
! CI group are identical (e.g. CH3C.(OO.)CH3). The routine ignore 
! zconform in these cases.                             
!=======================================================================
SUBROUTINE hot_criegee(xcri,ycri,zconform,brch,pozatom,loendo, &
                       np,s,p,nref,ref)
  USE keyparameter, ONLY: mxnode,mxlgr,mxring,mxcp,mxpd,mxpd,mxcopd  
  USE keyflag, ONLY: screenfg
  USE criegeetool, ONLY:init_done,idrchoo,idrrcoo,check_ze,ci_unipdct,add2p,&
                   init_stab_criegee,get_synanti,syn_index,anti_index,&
                   alphao_criegee,sci_yield
  USE rjtool, ONLY: rjgrm
  USE stdgrbond, ONLY: grbond
  USE reactool, ONLY: swap,rebond
  USE mapping, ONLY: abcde_map,ciptree
  USE ringtool, ONLY: findring
  IMPLICIT NONE

  CHARACTER(LEN=*),INTENT(IN) :: xcri     ! formula of the hot criegee (CI)
  REAL,INTENT(IN)       :: ycri           ! yield of the CI given as input in current rxn
  REAL,INTENT(IN)       :: zconform       ! ratio of the Z conformer for the CI as input
  REAL,INTENT(IN)       :: brch           ! branching ratio of the hot criegee
  INTEGER,INTENT(IN)    :: pozatom        ! number of O, C, N atoms in POZ
  LOGICAL,INTENT(IN)    :: loendo         ! true if C=C in parent is endocyclic
  INTEGER,INTENT(INOUT) :: np             ! # of product in the s and p list
  REAL,INTENT(INOUT)    :: s(:)           ! stoi. coef. of product in p(:) list
  CHARACTER(LEN=*),INTENT(INOUT) :: p(:)  ! product list (short names)
  INTEGER,INTENT(INOUT) :: nref           ! # of references in the reference list
  CHARACTER(LEN=*),INTENT(INOUT):: ref(:) ! list of references for current rxn
  
  CHARACTER(LEN=LEN(xcri)) :: tempfo
  CHARACTER(LEN=mxlgr) :: group(mxnode)
  CHARACTER(LEN=mxlgr) :: tgroup(mxnode),tempgr,pold, pnew
  INTEGER :: bond(mxnode,mxnode), tbond(mxnode,mxnode)
  INTEGER :: dbflg,nring,zeisom
  INTEGER :: ig,i,ngr
  INTEGER :: rjg(mxring,2)
  REAL    :: ysci, yuni
  INTEGER :: idci
  INTEGER :: zecase
  CHARACTER(LEN=12),PARAMETER :: progname='hot_criegee'
  CHARACTER(LEN=80) :: mesg

  INTEGER :: noda,nodb                    ! node "a" and "b" next to criegee group
  INTEGER :: ntreea(mxnode)               ! # of distinct nodes at depth j from noda
  INTEGER :: treea(mxnode,mxcp,mxnode)    ! node tree starting from noda
  INTEGER :: ntreeb(mxnode)               ! # of distinct nodes at depth j from nodb
  INTEGER :: treeb(mxnode,mxcp,mxnode)    ! node tree starting from nodb

  INTEGER,PARAMETER :: mxcip=50           ! CIP: Cahn-Ingold-Prelog (CIP) priority rules
  INTEGER :: cipa(mxcip,mxcip)            ! CIP tree, starting from 1st atom in branch "a" 
  INTEGER :: cipb(mxcip,mxcip)            ! CIP tree, starting from 1st atom in branch "b"

  INTEGER :: synnod,antinod               ! node # in syn and anti position
  INTEGER :: ntreesyn(mxnode)             ! # of distinct nodes at depth j from synnod
  INTEGER :: treesyn(mxnode,mxcp,mxnode)  ! node tree starting from synnod
  INTEGER :: ntreeanti(mxnode)            ! # of distinct nodes at depth j from antinod
  INTEGER :: treeanti(mxnode,mxcp,mxnode) ! node tree starting from antinod
  INTEGER :: idsyn,idanti                 ! SAR index of the syn and anti branch
  INTEGER :: alphaO                       ! -O- found next to the criegee center
  INTEGER :: rngflg                       ! flag for ring
  INTEGER :: ring(mxnode)                 ! nodes included in a ring

  IF (screenfg) WRITE(6,*)progname

! load data to handle the chemistry of criegee (1st call only)
  IF (.NOT. init_done) CALL init_stab_criegee()

! get # of nodes, functional groups and bond-matrix of CI*
  CALL grbond(xcri,group,bond,dbflg,nring)
  ngr=COUNT(group/=' ')
  IF (nring/=0) CALL rjgrm(nring,group,rjg) ! remove the ring joining char
  tgroup(:)=group(:) ; tbond(:,:)=bond(:,:)
  
! find the hot_criegee group
  ig=0 
  DO i=1,ngr
    IF (INDEX(tgroup(i),'.(OO.)*')/=0) ig=i
  ENDDO
  IF (ig==0) THEN
    mesg="expected hot_criegee functional group not found"
    CALL stoperr(progname,mesg,xcri)
  ENDIF

! single carbon criegee - treat each case explicitly and return
  IF (ngr==1) THEN
    idsyn=1 ; ysci=sci_yield(xcri,pozatom,idsyn)             ! SCI yield
    CALL c1_hotcriegee(xcri,brch,ycri,ysci,np,s,p,nref,ref)
    RETURN
  ENDIF

! get nodes bounded to CI center and create a node tree for each branch
  noda=0 ; ntreea(:)=0 ; treea(:,:,:)=0 ; alphaO=0
  nodb=0 ; ntreeb(:)=0 ; treeb(:,:,:)=0
  DO i=1,ngr
    IF (bond(ig,i)/=0) THEN 
      IF      (noda==0) THEN ; noda=i ; tbond(ig,noda)=0 ; tbond(noda,ig)=0
      ELSE IF (nodb==0) THEN ; nodb=i ; tbond(ig,nodb)=0 ; tbond(nodb,ig)=0
      ELSE 
        mesg="more than 2 nodes bonded to a criegee group"
        CALL stoperr(progname,mesg,xcri)
      ENDIF
    ENDIF
  ENDDO
  IF (noda/=0) THEN
    CALL abcde_map(tbond,noda,ngr,ntreea,treea)
    IF (group(noda)=='-O-')  alphaO=alphaO+1
  ENDIF
  IF (nodb/=0) THEN 
    CALL abcde_map(tbond,nodb,ngr,ntreeb,treeb)
    IF (group(nodb)=='-O-')  alphaO=alphaO+1
  ENDIF
  IF (noda==0) THEN
    mesg="expected branch on criegee not found"
    CALL stoperr(progname,mesg,xcri)
  ENDIF
  IF (nodb==0) THEN ; idci=idrchoo
  ELSE              ; idci=idrrcoo
  ENDIF
  tbond(:,:)=bond(:,:)

! SPECIAL CASE: R-O-C.(OO.)*-R
  IF (alphaO/=0) THEN
    CALL alphao_criegee(xcri,group,bond,ngr,brch,ig,noda,nodb, &
                        np,s,p,nref,ref)
    RETURN
  ENDIF

! add a reference for aromatic criegee (branch not in the SAR)
  IF (group(noda)(1:1)=='c') THEN
    CALL addref(progname,'AROCRIEGEE',nref,ref,xcri)
  ELSEIF (nodb/=0) THEN 
    IF (group(nodb)(1:1)=='c') CALL addref(progname,'AROCRIEGEE',nref,ref,xcri)
  ENDIF

  rngflg=0 ; ring(:)=0
  IF (nring/=0) CALL findring(ig,noda,ngr,bond,rngflg,ring)
  
!  --------------------------------
!  ------  EXTERNAL CRIEGEE -------
!  --------------------------------
  IF (idci==idrchoo) THEN
    IF ((group(noda)(1:1)=='C').OR.(group(noda)(1:1)=='c')) THEN

!==== RCH.(OO.)
      IF (group(ig)=='CH.(OO.)*') THEN

!------ CIS (Z) CONFORMER
        IF (zconform==0.) THEN   ! BABA not sure this is OK, preliminary stop
          mesg="ratio for the Z conformer is 0. Check and rewrite the IF statement if OK"
          CALL stoperr(progname,mesg,xcri)
        ENDIF

! get index for the syn group and compute the yield for the stabilized CI (SCI) 
        idanti=1                                                      ! H (index 1) is the anti branch 
        CALL syn_index(xcri,tgroup,bond,ring,ngr,ig,noda,ntreea,treea,idsyn)  ! get syn index (tree "a" is syn)
        ysci=sci_yield(xcri,pozatom,idsyn)                            ! SCI yield
        IF (loendo) ysci=1./(1.+580.*exp(-(pozatom-3)/2.))            ! overwrite if endocyclic parent
        yuni=1.-ysci                                                  ! unimolecular decomp. yield

! overall yields - rescale using Z conformer ratio and CI* yield in current rxn
        ysci=ysci*ycri*zconform 
        yuni=yuni*ycri*zconform              

! make SCI - make cold criegee, add pdct to the stack and reset
        IF (ysci/=0.) THEN  
          pold='.(OO.)*' ;  pnew='.(ZOO.)' ; tempgr=tgroup(ig)  
          CALL swap(tempgr, pold, tgroup(ig), pnew)
          CALL rebond(tbond,tgroup,tempfo,nring)
          CALL add2p(xcri,tempfo,ysci,brch,np,s,p,nref,ref)   ! check & add to s & p
          tgroup(:)=group(:) ; tbond(:,:)=bond(:,:)
        ENDIF

! make unimolecular decomposition products, add to the stack and reset
        pold='.(OO.)*' ;  pnew='.(OO.)' ; tempgr=tgroup(ig)  
        CALL swap(tempgr, pold, tgroup(ig), pnew)
        CALL ci_unipdct(xcri,tgroup,tbond,ngr,brch,yuni,ig,  &
                        idsyn,noda,ntreea,treea,idanti,nodb, &
                        np,s,p,nref,ref)
        tgroup(:)=group(:) ; tbond(:,:)=bond(:,:)
        
!------ TRANS (E) CONFORMER
        IF (zconform >= 1.) THEN   ! BABA not sure this is OK, preliminary stop
          mesg="ratio for the Z conformer is 1. Check and rewrite the IF statement if OK"
          CALL stoperr(progname,mesg,xcri)
        ENDIF

! get index for the syn/anti groups - compute SCI yield  
        idsyn=1                                                      ! H (index 1) is the syn branch  
        ysci=sci_yield(xcri,pozatom,idsyn)                           ! SCI yield
        IF (loendo) ysci=1./(1.+580.*exp(-(pozatom-3)/2.))           ! overwrite if endocyclic parent
        CALL anti_index(xcri,tgroup,noda,idanti)                     ! get anti index (tree "a" is anti branch)
        yuni=1.-ysci                                                 ! unimolecular decomp. yield

! overall yields - rescale using E conformer ratio and CI* yield in current rxn
        ysci=ysci*ycri*(1.-zconform) 
        yuni=yuni*ycri*(1.-zconform)  

! make SCI - make cold criegee, add pdct to the stack and reset 
        pold='.(OO.)*' ;  pnew='.(EOO.)' ; tempgr=tgroup(ig)  
        CALL swap(tempgr, pold, tgroup(ig), pnew)
        CALL rebond(tbond,tgroup,tempfo,nring)
        CALL add2p(xcri,tempfo,ysci,brch,np,s,p,nref,ref)    ! check & add to s & p list
        tgroup(:)=group(:) ; tbond(:,:)=bond(:,:)

! make unimolecular decomposition product, add to the stack and reset
        pold='.(OO.)*' ;  pnew='.(OO.)' ; tempgr=tgroup(ig)  
        CALL swap(tempgr, pold, tgroup(ig), pnew)
        CALL ci_unipdct(xcri,tgroup,tbond,ngr,brch,yuni,ig, &
                        idsyn,nodb,ntreeb,treeb,idanti,noda,&
                        np,s,p,nref,ref)
        tgroup(:)=group(:) ; tbond(:,:)=bond(:,:)
        
! === RC(OH).(OO.)
      ELSEIF (group(ig)=='C(OH).(OO.)*') THEN
        mesg="chemistry for C(OH).(OO.)* still needs development"

        CALL stoperr(progname,mesg,xcri)

! === R-C(ONO2).(OO.)
      ELSEIF (group(ig)=='C(ONO2).(OO.)*') THEN
        mesg="chemistry for C(ONO2).(OO.)* still needs development"
        CALL stoperr(progname,mesg,xcri)

! === R-C(OOH).(OO.) => no information. Assume OH elimination (100% yield) but
!                       another route (via dioxirane) could be CO2+HO2+R. 
      ELSEIF (group(ig)=='C(OOH).(OO.)*') THEN
        CALL add1tonp(progname,xcri,np) ; s(np)=ycri ; p(np)='HO  '   ! add HO
        
        pold='C(OOH).(OO.)*' ; pnew='CO(OO.)' ; tempgr=tgroup(ig)
        CALL swap(tempgr,pold,tgroup(ig),pnew)
        CALL rebond(tbond,tgroup,tempfo,nring)
        yuni=ycri
        CALL add2p(xcri,tempfo,yuni,brch,np,s,p,nref,ref)

! === R-C(OOOH).(OO.)
      ELSEIF (group(ig)=='C(OOOH).(OO.)*') THEN
        mesg="chemistry for C(OOOH).(OO.)* still need development"
        CALL stoperr(progname,mesg,xcri)

! === unidentified external criegee group
      ELSE
        mesg="unidentified external criegee group"
        CALL stoperr(progname,mesg,xcri)
      ENDIF
      
    ELSE
      mesg="unidentified or unexpected group bonded to criegee node"
      CALL stoperr(progname,mesg,xcri)
    ENDIF

!  --------------------------------
!  ------ INTERNAL CRIEGEES -------
!  --------------------------------
  ELSE IF (idci==idrrcoo) THEN

! check that hot criegee make sense
    IF (tgroup(ig)(1:8)/='C.(OO.)*') THEN
      mesg="ID 'idrrcoo' for criegee not OK. Group expected to be C.(OO.)* "
      CALL stoperr(progname,mesg,xcri)
    ENDIF

! check cis/trans: create node trees and CIP trees for both branches
    CALL ciptree(xcri,tbond,tgroup,ngr,ntreea,treea,cipa)  ! get 1st CIP tree 
    CALL ciptree(xcri,tbond,tgroup,ngr,ntreeb,treeb,cipb)  ! get 2nd CIP tree
    zeisom=check_ze(cipa,cipb)                             ! check if 2 distinct branches 

!----- 2 IDENTICAL BRANCHES: NO Z/E CONFORMER

    IF (zeisom==0) THEN

! get index for the syn/anti groups - compute SCI yield  
      CALL syn_index(xcri,tgroup,bond,ring,ngr,ig,noda,ntreea,treea,idsyn)  ! get syn index 
      CALL anti_index(xcri,tgroup,nodb,idanti)                      ! get anti index 
      ysci=sci_yield(xcri,pozatom,idsyn)                            ! sci yield
      IF (loendo) ysci=1./(1.+580.*exp(-(pozatom-3)/2.))            ! overwrite if endocyclic parent
      yuni=1.-ysci                                                  ! unimolecular decomp. yield

! overall yields - rescale using CI* yield in current rxn
      ysci=ysci*ycri 
      yuni=yuni*ycri                               

! make SCI - make cold criegee, add pdct to the stack and reset 
      pold='.(OO.)*' ; pnew='.(OO.)' ; tempgr=tgroup(ig)
      CALL swap(tempgr, pold, tgroup(ig), pnew)
      CALL rebond(tbond,tgroup,tempfo,nring)
      CALL add2p(xcri,tempfo,ysci,brch,np,s,p,nref,ref)       ! check & add to s & p list
      tgroup(:)=group(:) ; tbond(:,:)=bond(:,:)

! make unimolecular decomposition product (tree "a" and "b" are identical)
      pold='.(OO.)*' ;  pnew='.(OO.)' ; tempgr=tgroup(ig)  
      CALL swap(tempgr, pold, tgroup(ig), pnew)
      CALL ci_unipdct(xcri,tgroup,tbond,ngr,brch,yuni,ig, &
                      idsyn,noda,ntreea,treea,idanti,nodb, &
                      np,s,p,nref,ref)
      tgroup(:)=group(:) ; tbond(:,:)=bond(:,:)

    ELSE

!------ CIS (Z) CONFORMER

! identify the syn and anti branches for the Z conformer  
      zecase=1         ! set syn branch for the Z conformer
      CALL get_synanti(xcri,zecase,ntreea,treea,ntreeb,treeb,cipa,cipb, &  
                       synnod,antinod,ntreesyn,treesyn,ntreeanti,treeanti)

! get index for the syn/anti groups - compute SCI yield  
      CALL syn_index(xcri,tgroup,bond,ring,ngr,ig,synnod,ntreesyn,treesyn,idsyn)  ! get syn index
      CALL anti_index(xcri,tgroup,antinod,idanti)                         ! get anti index
      ysci=sci_yield(xcri,pozatom,idsyn)                                  ! sci yield
      IF (loendo) ysci=1./(1.+580.*exp(-(pozatom-3)/2.))                  ! overwrite if endocyclic parent
      yuni=1.-ysci                                                        ! decomp. yield

! overall yields - rescale using Z conformer ratio and CI* yield in current rxn
      ysci=ysci*ycri*zconform 
      yuni=yuni*ycri*zconform 

! make SCI - make cold criegee, add pdct to the stack and reset
      pold='.(OO.)*' ; pnew='.(ZOO.)' ; tempgr=tgroup(ig)
      CALL swap(tempgr, pold, tgroup(ig), pnew)
      CALL rebond(tbond,tgroup,tempfo,nring)
      CALL add2p(xcri,tempfo,ysci,brch,np,s,p,nref,ref)
      tgroup(:)=group(:) ; tbond(:,:)=bond(:,:)

! make unimolecular decomposition products, add to the stack and reset
      pold='.(OO.)*' ;  pnew='.(OO.)' ; tempgr=tgroup(ig)  
      CALL swap(tempgr, pold, tgroup(ig), pnew)
      CALL ci_unipdct(xcri,tgroup,tbond,ngr,brch,yuni,ig, &
                      idsyn,synnod,ntreesyn,treesyn,idanti,antinod,&
                      np,s,p,nref,ref)
      tgroup(:)=group(:) ; tbond(:,:)=bond(:,:)

!------ TRANS (E) CONFORMER

! identify the syn and anti branches for the E conformer  
      zecase=2   ! set syn branch for the E conformer
      CALL get_synanti(xcri,zecase,ntreea,treea,ntreeb,treeb,cipa,cipb, &  
                       synnod,antinod,ntreesyn,treesyn,ntreeanti,treeanti)

! get index for the syn/anti groups - compute SCI yield  
      CALL syn_index(xcri,tgroup,bond,ring,ngr,ig,synnod,ntreesyn,treesyn,idsyn)  ! get syn index
      CALL anti_index(xcri,tgroup,antinod,idanti)                         ! get anti index
      ysci=sci_yield(xcri,pozatom,idsyn)                                  ! sci yield
      IF (loendo) ysci=1./(1.+580.*exp(-(pozatom-3)/2.))                  ! overwrite if endocyclic parent
      yuni=1.-ysci                                                        ! decomp. yield

! overall yields - rescale using E conformer ratio and CI* yield in current rxn
      ysci=ysci*ycri*(1.-zconform) 
      yuni=yuni*ycri*(1.-zconform)

! make SCI - make cold criegee, add pdct to the stack and reset
      pold='.(OO.)*' ; pnew='.(EOO.)' ; tempgr=tgroup(ig)
      CALL swap(tempgr, pold, tgroup(ig), pnew)
      CALL rebond(tbond,tgroup,tempfo,nring)
      CALL add2p(xcri,tempfo,ysci,brch,np,s,p,nref,ref)
      tgroup(:)=group(:) ; tbond(:,:)=bond(:,:)

! make unimolecular decomposition products, add to the stack and reset
      pold='.(OO.)*' ;  pnew='.(OO.)' ; tempgr=tgroup(ig)  
      CALL swap(tempgr, pold, tgroup(ig), pnew)
      CALL ci_unipdct(xcri,tgroup,tbond,ngr,brch,yuni,ig, &
                      idsyn,synnod,ntreesyn,treesyn,idanti,antinod,&
                      np,s,p,nref,ref)
      tgroup(:)=group(:) ; tbond(:,:)=bond(:,:)
    
    ENDIF
  ENDIF

END SUBROUTINE hot_criegee


! ======================================================================
! PURPOSE: return products (as s & p) of the C1 hot criegee (CI) 
! provided as input.
! ======================================================================
SUBROUTINE c1_hotcriegee(xcri,brch,ycri,ysci,np,s,p,nref,ref)
  USE dictstacktool, ONLY: bratio
  USE keyflag, ONLY: screenfg
  IMPLICIT NONE
  
  CHARACTER(LEN=*),INTENT(IN) :: xcri    ! formula of the hot criegee
  REAL,INTENT(IN)    :: brch             ! branching ratio of the hot criegee
  REAL,INTENT(IN)    :: ycri             ! yield of the hot criegee provided as input
  REAL,INTENT(IN)    :: ysci             ! yield of SCI from hot criegee
  INTEGER,INTENT(INOUT) :: np            ! # of product in the s and p list
  REAL,INTENT(INOUT)   :: s(:)           ! stoi. coef. of product in p(:) list
  CHARACTER(LEN=*),INTENT(INOUT) :: p(:) ! product list (short names)
  INTEGER,INTENT(INOUT) :: nref          ! # of references added in the reference list
  CHARACTER(LEN=*),INTENT(INOUT):: ref(:)! list of references 

  CHARACTER(LEN=LEN(xcri)) :: tempfo
  REAL           :: brtio

  CHARACTER(LEN=12),PARAMETER :: progname='hot_criegee'
  CHARACTER(LEN=80) :: mesg

  IF (screenfg) WRITE(6,*)progname

! -----------
! CH2.(OO.)*. 
! -----------
  IF (xcri(1:10)=='CH2.(OO.)*') THEN

! SCI channel
    tempfo='CH2.(OO.)' 
    CALL add1tonp(progname,xcri,np) ;  s(np)= ycri*ysci 
    brtio=brch*s(np)  ;  CALL bratio(tempfo,brtio,p(np),nref,ref)
  
! decomposition channel (H2+CO2 and H2O+CO, 50 % each)
    CALL add1tonp(progname,xcri,np) ; s(np)=0.5*ycri*(1.-ysci) ; p(np)='CO2 '
    CALL add1tonp(progname,xcri,np) ; s(np)=0.5*ycri*(1.-ysci) ; p(np)='H2  '
    CALL add1tonp(progname,xcri,np) ; s(np)=0.5*ycri*(1.-ysci) ; p(np)='CO  '
    CALL add1tonp(progname,xcri,np) ; s(np)=0.5*ycri*(1.-ysci) ; p(np)='H2O '

    RETURN

! -----------
! CH(OH).(OO.) => assume dioxirane C1H(OO1)(OH) decompose to CO2+OH+H
! ------------
  ELSE IF (xcri(1:13)=='CH(OH).(OO.)*') THEN
    CALL add1tonp(progname,xcri,np) ; s(np)=ycri ; p(np)='CO2  '
    CALL add1tonp(progname,xcri,np) ; s(np)=ycri ; p(np)='HO2  ' 
    CALL add1tonp(progname,xcri,np) ; s(np)=ycri ; p(np)='HO   ' 
    RETURN

! --------------
! CH(ONO2).(OO.) => assume dioxirane decompose to CO2+OH+NO2  
! --------------
  ELSE IF (xcri(1:15)=='CH(ONO2).(OO.)*') THEN
    CALL add1tonp(progname,xcri,np) ; s(np)=ycri ; p(np)='CO2  '
    CALL add1tonp(progname,xcri,np) ; s(np)=ycri ; p(np)='HO   ' 
    CALL add1tonp(progname,xcri,np) ; s(np)=ycri ; p(np)='NO2  ' 
    RETURN

! --------------
! CH(OOH).(OO.) =>  assume decompose to CO+OH+HO2
! --------------
  ELSE IF (xcri(1:14)=='CH(OOH).(OO.)*') THEN
    CALL add1tonp(progname,xcri,np) ; s(np)=ycri ; p(np)='CO   '
    CALL add1tonp(progname,xcri,np) ; s(np)=ycri ; p(np)='HO   ' 
    CALL add1tonp(progname,xcri,np) ; s(np)=ycri ; p(np)='HO2  ' 
    RETURN
      
! --------------
! CH(OOOH).(OO.) => assume decompose to CO2+OH+HO2
! --------------
  ELSE IF (xcri(1:15)=='CH(OOOH).(OO.)*') THEN
    CALL add1tonp(progname,xcri,np) ; s(np)=1.0 ; p(np)='CO2  '
    CALL add1tonp(progname,xcri,np) ; s(np)=1.0 ; p(np)='HO   ' 
    CALL add1tonp(progname,xcri,np) ; s(np)=1.0 ; p(np)='HO2  ' 
    RETURN
      
! --------------
! C(OOOH)(OOOH).(OO.) => assume decompose to CO2+2OH+2O2
! --------------
  ELSE IF (xcri(1:20)=='C(OOOH)(OOOH).(OO.)*') THEN
    CALL add1tonp(progname,xcri,np) ; s(np)=ycri   ; p(np)='CO2  '
    CALL add1tonp(progname,xcri,np) ; s(np)=2*ycri ; p(np)='HO   ' 
    CALL add1tonp(progname,xcri,np) ; s(np)=2*ycri ; p(np)='O2   ' 
    RETURN
      
! --------------
! CO.(OO.)"      => assume decompose to CO+O2
! --------------
  ELSE IF (xcri(1:9)=='CO.(OO.)*') THEN
    CALL add1tonp(progname,xcri,np) ; s(np)=ycri ; p(np)='CO   '
    CALL add1tonp(progname,xcri,np) ; s(np)=ycri ; p(np)='O2   '
    RETURN

  ELSE
    WRITE(mesg,*) "C1 hot_criegee not found. Change the program accordingly"
    CALL stoperr(progname,mesg,xcri)
  ENDIF
  
END SUBROUTINE c1_hotcriegee

END MODULE hotcriegeechem
