MODULE stbcriegeechem
USE keyparameter, ONLY: mxpd,mxnr,mxcopd,mecu,mxlco
USE reactool, ONLY: swap,rebond
USE dictstacktool, ONLY: bratio
USE normchem, ONLY: stdchm
USE rxwrttool, ONLY:rxwrit,rxinit
USE references, ONLY: mxlcod
USE toolbox, ONLY: stoperr,add1tonp,addref
IMPLICIT NONE

CONTAINS

! ======================================================================
! PURPOSE: performs the reactions of the stabilized criegee 
! intermediate (SCI). The following reactions are considered for the
! stabilized criegee (see Newland et al., 2022):
! - Unimolecular decomposition 
! - Reaction with water (monomer and dimer)
! - Reaction with SO2, NO, NO2, CO, HCl, HNO3, O3
! - Self/cross SCI reaction - need to be developped if needed (counters
!   might be introduced to count total SCI concentration)
! The routine avoids reactions with organics for  now. No reaction with
! NH3 because it makes amines (not handled currently).
! ======================================================================
SUBROUTINE stab_criegee(idnam,chem,bond,group,brch)
  USE keyparameter, ONLY: mxcp,mxnode
  USE mapping, ONLY: abcde_map,ciptree
  USE toolbox, ONLY: kval
  USE criegeetool, ONLY: init_done,idrchoo,idrrcoo, &
                init_stab_criegee,get_synanti,add2p,syn_index,anti_index
  USE keyflag, ONLY: bimolecrx4criegee,screenfg
  USE ringtool, ONLY: findring
  IMPLICIT NONE

  CHARACTER(LEN=*),INTENT(IN) :: idnam    ! name of current species
  CHARACTER(LEN=*),INTENT(IN) :: chem     ! formula of current species
  CHARACTER(LEN=*),INTENT(IN) :: group(:) ! group in chem
  INTEGER,INTENT(IN) :: bond(:,:)         ! bond matrix
  REAL,INTENT(IN)    :: brch              ! max yield of the input species

  CHARACTER(LEN=LEN(group(1))) :: tgroup(SIZE(group))      ! cp of group
  CHARACTER(LEN=LEN(group(1))) :: nozegroup(SIZE(group))   ! cp of group without Z/E 
  CHARACTER(LEN=LEN(group(1))) :: pold,pnew,tempgr
  INTEGER :: tbond(SIZE(bond,1),SIZE(bond,2))
  INTEGER :: ngr,ip,i,idci,zecase
  LOGICAL :: loze

  INTEGER :: noda,nodb                    ! node "a" and "b" next to criegee group
  INTEGER :: ntreea(mxnode)               ! # of distinct nodes at depth j from noda
  INTEGER :: treea(mxnode,mxcp,mxnode)    ! node tree starting from noda
  INTEGER :: ntreeb(mxnode)               ! # of distinct nodes at depth j from nodb
  INTEGER :: treeb(mxnode,mxcp,mxnode)    ! node tree starting from nodb
  
  INTEGER :: synnod,antinod               ! node # in syn and anti position
  INTEGER :: ntreesyn(mxnode)             ! # of distinct nodes at depth j from synnod
  INTEGER :: treesyn(mxnode,mxcp,mxnode)  ! node tree starting from synnod
  INTEGER :: ntreeanti(mxnode)            ! # of distinct nodes at depth j from antinod
  INTEGER :: treeanti(mxnode,mxcp,mxnode) ! node tree starting from antinod
  INTEGER :: idsyn,idanti                 ! SAR index of the syn and anti branch
  INTEGER :: rngflg                       ! flag for ring
  INTEGER :: ring(mxnode)                 ! nodes included in a ring

  INTEGER,PARAMETER :: mxcip=50           ! CIP: Cahn-Ingold-Prelog (CIP) priority rules
  INTEGER :: cipa(mxcip,mxcip)            ! CIP tree, starting from 1st atom in branch "a" 
  INTEGER :: cipb(mxcip,mxcip)            ! CIP tree, starting from 1st atom in branch "b"

  CHARACTER(LEN=12)     :: progname='stab_criegee'
  CHARACTER(LEN=70)     :: mesg
  INTEGER,PARAMETER     :: mxcom=10

  IF (screenfg) WRITE(6,*)progname

! load data to handle the chemistry of criegee (1st call only)
  IF (.NOT. init_done) CALL init_stab_criegee()

  nozegroup(:)=group(:)  ;  tbond(:,:)=bond(:,:)
  ntreea(:)=0   ; treea(:,:,:)=0   ; ntreeb(:)=0    ; treeb(:,:,:)=0
  ntreesyn(:)=0 ; treesyn(:,:,:)=0 ; ntreeanti(:)=0 ; treeanti(:,:,:)=0
  cipa(:,:)=0   ; cipb(:,:)=0

  ngr=COUNT(group/=' ')
  IF (ngr<2) THEN
    mesg="At least a 2 carbon species is expected as input "
    CALL stoperr(progname,mesg,chem)
  ENDIF

! check if species is allowed in this routine
  loze=.FALSE. ; ip=0 ; zecase=0
  DO i=1,ngr
    IF (INDEX(group(i),'.(OO.)')/=0)  ip=i
    IF (INDEX(group(i),'.(ZOO.)')/=0) ip=i
    IF (INDEX(group(i),'.(EOO.)')/=0) ip=i
  ENDDO
  IF (ip==0) THEN
    mesg="No criegee found "
    CALL stoperr(progname,mesg,chem)
  ENDIF

  IF (INDEX(group(ip),'*')/=0) THEN   ! hot criegee not allowed.
    mesg="Hot_criegee not allowed "
    CALL stoperr(progname,mesg,chem)
  ENDIF

! Save Z/E information and remove Z/E in group, if any
  IF (INDEX(group(ip),'.(ZOO.)')/=0) THEN
    loze=.TRUE. ; zecase=1
    pold='(ZOO.)' ;  pnew='(OO.)' ; tempgr=group(ip)
    CALL swap(tempgr,pold,nozegroup(ip),pnew)
  ELSE IF (INDEX(group(ip),'.(EOO.)')/=0) THEN
    loze=.TRUE. ; zecase=2
    pold='(EOO.)' ;  pnew='(OO.)' ; tempgr=group(ip)
    CALL swap(tempgr,pold,nozegroup(ip),pnew)
  ENDIF
  tgroup(:)=nozegroup(:)

! -----------------------
! CRIEGEE TYPE 
! -----------------------

! get nodes bounded to CI group and create a node tree for each branch
  noda=0 ; nodb=0
  DO i=1,ngr
    IF (bond(ip,i)/=0) THEN
      IF      (noda==0) THEN ; noda=i ; tbond(ip,noda)=0 ; tbond(noda,ip)=0
      ELSE IF (nodb==0) THEN ; nodb=i ; tbond(ip,nodb)=0 ; tbond(nodb,ip)=0
      ELSE
        mesg="more than 2 nodes bonded to a criegee group"
        CALL stoperr(progname,mesg,chem)
      ENDIF
    ENDIF
  ENDDO
  IF (noda/=0) CALL abcde_map(tbond,noda,ngr,ntreea,treea)
  IF (nodb/=0) CALL abcde_map(tbond,nodb,ngr,ntreeb,treeb)

  IF (noda==0) THEN
    mesg="expected branch on criegee not found"
    CALL stoperr(progname,mesg,chem)
  ENDIF
  IF (nodb==0) THEN ; idci=idrchoo
  ELSE              ; idci=idrrcoo
  ENDIF

  CALL findring(ip,noda,ngr,bond,rngflg,ring)

! External criegee
  IF (idci==idrchoo) THEN 
    IF (zecase==1) THEN                                   ! cis (ZOO.) configuration 
      idanti=1                                            ! H (index 1) is the anti branch  
      synnod=noda ; treesyn=treea ; ntreesyn=ntreea 
      CALL syn_index(chem,tgroup,bond,ring,ngr,ip,synnod,ntreesyn,treesyn,idsyn)  ! get syn index
    ELSE IF (zecase==2) THEN                              ! trans (EOO.) configuration
      idsyn=1                                             ! H (index 1) is the syn branch
      antinod=noda ; treeanti=treea ; ntreeanti=ntreea
      CALL anti_index(chem,tgroup,antinod,idanti)         ! get anti index
    ELSE 
      mesg="expected Z/E information not found on criegee"
      CALL stoperr(progname,mesg,chem)
    ENDIF
  
! Internal criegee
  ELSE IF (idci==idrrcoo) THEN    
    IF (loze) THEN                           ! identify which branch is syn or anti
      CALL ciptree(chem,tbond,tgroup,ngr,ntreea,treea,cipa)                 ! get 1st CIP tree 
      CALL ciptree(chem,tbond,tgroup,ngr,ntreeb,treeb,cipb)                 ! get 2nd CIP tree
      CALL get_synanti(chem,zecase,ntreea,treea,ntreeb,treeb,cipa,cipb, &   ! identify syn & ...
                       synnod,antinod,ntreesyn,treesyn,ntreeanti,treeanti)  ! ... anti branches
      CALL syn_index(chem,tgroup,bond,ring,ngr,ip,synnod,ntreesyn,treesyn,idsyn)    ! get syn index
      CALL anti_index(chem,tgroup,antinod,idanti)                           ! get anti index

    ELSE                                     ! 2 identical branches - set syn/anti arbitrarily
      synnod=noda ; treesyn=treea ; ntreesyn=ntreea
      CALL syn_index(chem,tgroup,bond,ring,ngr,ip,synnod,ntreesyn,treesyn,idsyn)    ! get syn index
      antinod=nodb ; treeanti=treeb ; ntreeanti=ntreeb
      CALL anti_index(chem,tgroup,antinod,idanti)                           ! get anti index
    ENDIF

  ENDIF
  tbond(:,:)=bond(:,:)

! ------------------------
! START REACTIONS        
! ------------------------

! Unimolecular decomposition
  CALL rx_sci_uni(idnam,chem,nozegroup,bond,ngr,brch,ip, &
                  idsyn,synnod,ntreesyn,treesyn,idanti,antinod)

! Reaction with water (monomer & dimer)
  CALL rx_sci_water(idnam,chem,nozegroup,bond,brch,ip,idsyn,idanti)

! Other bimolecular reactions
  IF (bimolecrx4criegee) THEN
    CALL rx_sci_so2(idnam,chem,nozegroup,bond,ip,idci,zecase,brch)
    CALL rx_sci_no2(idnam,chem,nozegroup,bond,ip,brch)
    CALL rx_sci_hno3(idnam,chem,nozegroup,bond,ip,brch)
    CALL rx_sci_rcooh(idnam,chem,nozegroup,bond,ngr,ip,idci,zecase,brch)
  ENDIF

END SUBROUTINE stab_criegee


! ======================================================================
! Purpose: perform the reaction of SCI with RCOOH. Rate are decribed in   
! Newland et al., 2022. Rate depend on the criegee ID # (idci):
! 1: CH2OO,   2: RCH.(ZOO.),   3: RCH.(ZOO.),   4: RC.(OO.)R
! The reaction is written to the output file within the subroutine.
! Two carboxylic acid are considered: HCOOH and CH3COOH. 
! ======================================================================
SUBROUTINE rx_sci_rcooh(idnam,chem,group,bond,ngr,cnod,idci,zecase,brch)
  USE keyparameter, ONLY: mxlco
  USE keyflag, ONLY: screenfg
  USE dictstacktool, ONLY: add_topvocstack
  IMPLICIT NONE

  CHARACTER(LEN=*),INTENT(IN) :: idnam      ! name of current criegee
  CHARACTER(LEN=*),INTENT(IN) :: chem       ! formula of current criegee
  CHARACTER(LEN=*),INTENT(IN) :: group(:)   ! group in chem
  INTEGER,INTENT(IN) :: bond(:,:)           ! bond matrix
  INTEGER,INTENT(IN) :: ngr                 ! # of nodes in chem
  INTEGER,INTENT(IN) :: cnod                ! criegee node
  INTEGER,INTENT(IN) :: idci                ! idci of current criegee
  INTEGER,INTENT(IN) :: zecase              ! conformer (1=Z ; 2=E)
  REAL,INTENT(IN)    :: brch                ! max yield of the input species

  CHARACTER(LEN=LEN(group(1))) :: tgroup(SIZE(group))      
  CHARACTER(LEN=LEN(group(1))) :: pold,pnew,tempgr
  CHARACTER(LEN=LEN(chem))     :: tchem,tchem2
  INTEGER :: tbond(SIZE(bond,1),SIZE(bond,2))

  CHARACTER(LEN=LEN(idnam)) :: p(mxpd), r(3)
  REAL    :: s(mxpd), arrh(3)
  INTEGER :: idreac, nlabel
  REAL    :: xlabel,folow(3),fotroe(4)
  REAL    :: brtio
  INTEGER :: tempring,io,tngr,icarb
  
  CHARACTER(LEN=mxlco) :: namacetic=' '   ! intialize at 1st call only
  
  CHARACTER(LEN=12)     :: progname='rx_sci_rcooh'
  CHARACTER(LEN=70)     :: mesg

  INTEGER,PARAMETER     :: mxcom=10
  INTEGER               :: nref
  CHARACTER(LEN=mxlcod) :: ref(mxcom)

  IF (screenfg) WRITE(6,*)progname

!----------------------
! REACTION WITH HCOOH
!----------------------
  tngr=ngr ; tgroup(:)=group(:) ; tbond(:,:)=bond(:,:)
  nref=0 ; ref(:)=' ' 

! make product: R1-C.(OO.)-R2 +HCOOH => HCO-O-C(OOH)R1R2
! -------------
  pold='.(OO.)' ; pnew='(OOH)' 
  tempgr=tgroup(cnod) ; CALL swap(tempgr,pold,tgroup(cnod),pnew)

  tngr=tngr+1 ; tgroup(tngr)='-O-' ; io=tngr
  tbond(tngr,cnod)=3 ; tbond(cnod,tngr)=3   

  tngr=tngr+1 ; tgroup(tngr)='CHO'
  tbond(tngr,io)=3 ; tbond(io,tngr)=3   
  CALL rebond(tbond,tgroup,tchem,tempring)
  CALL stdchm(tchem)

! make reaction
! -------------
  CALL rxinit(r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe)
  CALL addref(progname,'MN22SAR2',nref,ref,chem)
  
! set reactant
  r(1)=idnam  ;  r(2)='HCOOH '

! set product
  s(1)=1. ; brtio=brch*s(1)  
  CALL bratio(tchem,brtio,p(1),nref,ref)
 
! set rate constant
  SELECT CASE (idci)
  CASE(1)               ; arrh(1)=1.2E-10       ! CH2OO
  CASE(3)               ; arrh(1)=3.1E-10       ! RC.(OO.)R 
  CASE(2)  
    IF (zecase==1) THEN ; arrh(1)=2.1E-10       ! RCH.(ZOO.)
    ELSE                ; arrh(1)=3.8E-10       ! RCH.(EOO.)
    ENDIF 
  CASE DEFAULT 
    mesg="criegee group not identified "
    CALL stoperr(progname,mesg,chem)
  END SELECT

! write reaction
  CALL rxwrit(mecu,r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe,ref)

!----------------------
! REACTION WITH CH3COOH
!----------------------
  tngr=ngr ; tgroup(:)=group(:) ; tbond(:,:)=bond(:,:)
  nref=0 ; ref(:)=' ' 

! make product: R1-C.(OO.)-R2 + CH3COOH => CH3CO-O-C(OOH)R1R2
! -------------
  pold='.(OO.)' ; pnew='(OOH)' 
  tempgr=tgroup(cnod) ; CALL swap(tempgr,pold,tgroup(cnod),pnew)

  tngr=tngr+1 ; tgroup(tngr)='-O-' ; io=tngr
  tbond(tngr,cnod)=3 ; tbond(cnod,tngr)=3   

  tngr=tngr+1 ; tgroup(tngr)='CO' ; icarb=tngr
  tbond(tngr,io)=3 ; tbond(io,tngr)=3  
   
  tngr=tngr+1 ; tgroup(tngr)='CH3' 
  tbond(tngr,icarb)=1 ; tbond(icarb,tngr)=1  

  CALL rebond(tbond,tgroup,tchem,tempring)
  CALL stdchm(tchem)

! make reaction
! -------------
  CALL rxinit(r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe)
  CALL addref(progname,'MN22SAR2',nref,ref,chem)
  
! set reactant
  r(1)=idnam  

! add CH3COOH - if still unknown in the mechanism, add the species on top of 
! the VOC stack (with same generation # as the CI and an overall yield of 1)
  IF (namacetic==' ') THEN
    tchem2='CH3CO(OH) ' ; brtio=1.  ! CH3COOH assumed as a primary species
    CALL add_topvocstack(tchem2,brtio,namacetic,nref,ref)
  ENDIF
  r(2)=namacetic

! set product
  s(1)=1. ; brtio=brch*s(1)  
  CALL bratio(tchem,brtio,p(1),nref,ref)
 
! set rate constant
  SELECT CASE (idci)
  CASE(1)               ; arrh(1)=1.2E-10       ! CH2OO
  CASE(3)               ; arrh(1)=3.1E-10       ! RC.(OO.)R 
  CASE(2)  
    IF (zecase==1) THEN ; arrh(1)=2.1E-10       ! RCH.(ZOO.)
    ELSE                ; arrh(1)=3.8E-10       ! RCH.(EOO.)
    ENDIF 
  CASE DEFAULT 
    mesg="criegee group not identified "
    CALL stoperr(progname,mesg,chem)
  END SELECT

! write reaction
  CALL rxwrit(mecu,r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe,ref)

END SUBROUTINE rx_sci_rcooh

! ======================================================================
! Purpose: perform the reaction of SCI with HNO3. Rate are decribed in   
! Newland et al., 2022. Rate depend on the criegee ID # (idci):
! 1: CH2OO,   2: RCH.(ZOO.),   3: RCH.(ZOO.),   4: RC.(OO.)R
! The reaction is written to the output file within the subroutine.
! ======================================================================
SUBROUTINE rx_sci_hno3(idnam,chem,group,bond,cnod,brch)
  USE keyflag, ONLY: screenfg
  IMPLICIT NONE

  CHARACTER(LEN=*),INTENT(IN) :: idnam      ! name of current criegee
  CHARACTER(LEN=*),INTENT(IN) :: chem       ! formula of current criegee
  CHARACTER(LEN=*),INTENT(IN) :: group(:)   ! group in chem
  INTEGER,INTENT(IN) :: bond(:,:)           ! bond matrix
  INTEGER,INTENT(IN) :: cnod                ! criegee node
  REAL,INTENT(IN)    :: brch                ! max yield of the input species

  CHARACTER(LEN=LEN(group(1))) :: tgroup(SIZE(group))      
  CHARACTER(LEN=LEN(group(1))) :: pold,pnew,tempgr
  CHARACTER(LEN=LEN(chem))     :: tchem

  CHARACTER(LEN=LEN(idnam)) :: p(mxpd), r(3)
  REAL    :: s(mxpd), arrh(3)
  INTEGER :: idreac, nlabel
  REAL    :: xlabel,folow(3),fotroe(4)
  REAL    :: brtio
  INTEGER :: tempring
  
  CHARACTER(LEN=12)     :: progname='rx_sci_hno3'
!  CHARACTER(LEN=70)     :: mesg

  INTEGER,PARAMETER     :: mxcom=10
  INTEGER               :: nref
  CHARACTER(LEN=mxlcod) :: ref(mxcom)

  IF (screenfg) WRITE(6,*)progname

  tgroup(:)=group(:)
  nref=0 ; ref(:)=' ' 

! make product: switch -C.(OO.)- into -C(OOH)(ONO2)-
! -------------
  pold='.(OO.)' ; pnew='(OOH)(ONO2)' 
  tempgr=tgroup(cnod) ; CALL swap(tempgr,pold,tgroup(cnod),pnew)
  CALL rebond(bond,tgroup,tchem,tempring)
  CALL stdchm(tchem)

! make reaction
! -------------
  CALL rxinit(r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe)
  CALL addref(progname,'MN22SAR2',nref,ref,chem)
  
! set reactant
  r(1)=idnam  ;  r(2)='HNO3 '

! set product
  s(1)=1. ; brtio=brch*s(1)  
  CALL bratio(tchem,brtio,p(1),nref,ref)
 
! set rate constant
  arrh(1)=5.4E-10

! write reaction
  CALL rxwrit(mecu,r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe,ref)

END SUBROUTINE rx_sci_hno3 


! ======================================================================
! Purpose: perform the reaction of SCI with NO2. Rate are decribed in   
! Newland et al., 2022. Rate depend on the criegee ID # (idci):
! 1: CH2OO,   2: RCH.(ZOO.),   3: RCH.(ZOO.),   4: RC.(OO.)R
! The reaction is written to the output file within the subroutine.
! ======================================================================
SUBROUTINE rx_sci_no2(idnam,chem,group,bond,cnod,brch)
  USE keyflag, ONLY: screenfg
  IMPLICIT NONE

  CHARACTER(LEN=*),INTENT(IN) :: idnam      ! name of current criegee
  CHARACTER(LEN=*),INTENT(IN) :: chem       ! formula of current criegee
  CHARACTER(LEN=*),INTENT(IN) :: group(:)   ! group in chem
  INTEGER,INTENT(IN) :: bond(:,:)           ! bond matrix
  INTEGER,INTENT(IN) :: cnod                ! criegee node
  REAL,INTENT(IN)    :: brch                ! max yield of the input species

  CHARACTER(LEN=LEN(group(1))) :: tgroup(SIZE(group))      
  CHARACTER(LEN=LEN(group(1))) :: pold,pnew,tempgr
  CHARACTER(LEN=LEN(chem))     :: tchem

  CHARACTER(LEN=LEN(idnam)) :: p(mxpd), r(3)
  REAL    :: s(mxpd), arrh(3)
  INTEGER :: idreac, nlabel
  REAL    :: xlabel,folow(3),fotroe(4)
  REAL    :: brtio
  INTEGER :: tempring
  
  CHARACTER(LEN=12)     :: progname='rx_sci_no2'
!  CHARACTER(LEN=70)     :: mesg

  INTEGER,PARAMETER     :: mxcom=10
  INTEGER               :: nref
  CHARACTER(LEN=mxlcod) :: ref(mxcom)

  IF (screenfg) WRITE(6,*)progname

  tgroup(:)=group(:)
  nref=0 ; ref(:)=' ' 

! make product: switch -C.(OO.)- into -C(NO2)(OO.)-
! -------------
  pold='.(OO.)' ; pnew='(NO2)(OO.)' 
  tempgr=tgroup(cnod) ; CALL swap(tempgr,pold,tgroup(cnod),pnew)
  CALL rebond(bond,tgroup,tchem,tempring)
  CALL stdchm(tchem)

! make reaction
! -------------
  CALL rxinit(r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe)
  CALL addref(progname,'MN22SAR2',nref,ref,chem)
  
! set reactant
  r(1)=idnam  ;  r(2)='NO2 '

! set product
  s(1)=1. ; brtio=brch*s(1)  
  CALL bratio(tchem,brtio,p(1),nref,ref)
 
! set rate constant
  arrh(1)=2.0E-12

! write reaction
  CALL rxwrit(mecu,r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe,ref)

END SUBROUTINE rx_sci_no2 

! ======================================================================
! Purpose: perform the reaction of SCI with SO2. Rate are decribed in   
! Newland et al., 2022. Rate depend on the criegee ID number (idci):
! 1: CH2OO,   2: RCH.(ZOO.),   3: RCH.(ZOO.),   4: RC.(OO.)R
! The reaction is written to the output file within the subroutine.
! ======================================================================
SUBROUTINE rx_sci_so2(idnam,chem,group,bond,cnod,idci,zecase,brch)
  USE keyflag, ONLY: screenfg
  IMPLICIT NONE

  CHARACTER(LEN=*),INTENT(IN) :: idnam      ! name of current criegee
  CHARACTER(LEN=*),INTENT(IN) :: chem       ! formula of current criegee
  CHARACTER(LEN=*),INTENT(IN) :: group(:) ! group in chem
  INTEGER,INTENT(IN) :: bond(:,:)         ! bond matrix
  INTEGER,INTENT(IN) :: cnod              ! criegee node
  INTEGER,INTENT(IN) :: idci                ! idci of current criegee
  INTEGER,INTENT(IN) :: zecase              ! conformer (1=Z ; 2=E)
  REAL,INTENT(IN)    :: brch                ! max yield of the input species

  CHARACTER(LEN=LEN(group(1))) :: tgroup(SIZE(group))      
  CHARACTER(LEN=LEN(group(1))) :: pold,pnew,tempgr
  CHARACTER(LEN=LEN(chem))     :: p_carbonyl

  CHARACTER(LEN=LEN(idnam)) :: p(mxpd), r(3)
  REAL    :: s(mxpd), arrh(3)
  INTEGER :: idreac, nlabel
  REAL    :: xlabel,folow(3),fotroe(4)
  INTEGER :: tempring
  REAL    :: brtio
  
  CHARACTER(LEN=12)     :: progname='rx_sci_so2'
  CHARACTER(LEN=70)     :: mesg

  INTEGER,PARAMETER     :: mxcom=10
  INTEGER               :: nref
  CHARACTER(LEN=mxlcod) :: ref(mxcom)

  IF (screenfg) WRITE(6,*)progname

  tgroup(:)=group(:)

! make the carbonyl
  p_carbonyl=' '
  IF     (group(cnod)=='CH.(OO.)') THEN ; pold='CH.(OO.)' ; pnew='CHO'
  ELSEIF (group(cnod)=='C.(OO.)' ) THEN ; pold='C.(OO.)'  ; pnew='CO'
  ELSE
    mesg="criegee group not identified "
    CALL stoperr(progname,mesg,chem)
  ENDIF
  tempgr=tgroup(cnod)
  CALL swap(tempgr,pold,tgroup(cnod),pnew)
  CALL rebond(bond,tgroup,p_carbonyl,tempring)
  CALL stdchm(p_carbonyl)
  tgroup(:)=group(:)

  CALL rxinit(r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe)
  nref=0 ; ref(:)=' ' 
  CALL addref(progname,'MN22SAR2',nref,ref,chem)
 
! set reactant
  r(1)=idnam  ;  r(2)='SO2 '

! set products (carbonyl) and co-products
  s(1)=1. ; brtio=brch*s(1)  
  CALL bratio(p_carbonyl,brtio,p(1),nref,ref)
  s(2)=1.  ;  p(2)='SULF '

! set rate constant
  SELECT CASE (idci)
  CASE(1)               ; arrh(1)=3.7E-11       ! CH2OO
  CASE(3)               ; arrh(1)=1.6E-10       ! RC.(OO.)R 
  CASE(2)  
    IF (zecase==1) THEN ; arrh(1)=2.6E-11       ! RCH.(ZOO.)
    ELSE                ; arrh(1)=1.4E-10       ! RCH.(EOO.)
    ENDIF 
  CASE DEFAULT 
    mesg="criegee group not identified "
    CALL stoperr(progname,mesg,chem)
  END SELECT

  CALL rxwrit(mecu,r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe,ref)
  
END SUBROUTINE rx_sci_so2 

! ======================================================================
! Purpose: perform the reaction of SCI with water (monomer and dimer).
! Rate are from the Vereecken et al., 2017 SAR, as decribed in 
! Newland et al., 2022. The reaction is written to the output
! file within the subroutine.
! ======================================================================
SUBROUTINE rx_sci_water(idnam,chem,group,bond,brch,cnod,idsyn,idanti)
  USE criegeetool, ONLY:  k_mh2o,k_dh2o
  USE keyflag, ONLY: screenfg
  IMPLICIT NONE

  CHARACTER(LEN=*),INTENT(IN) :: idnam    ! name of current species
  CHARACTER(LEN=*),INTENT(IN) :: chem     ! formula of current species
  CHARACTER(LEN=*),INTENT(IN) :: group(:) ! group in chem
  INTEGER,INTENT(IN) :: bond(:,:)         ! bond matrix
  REAL,INTENT(IN)    :: brch              ! max yield of the input species
  INTEGER,INTENT(IN) :: cnod              ! criegee node
  INTEGER,INTENT(IN) :: idsyn             ! SAR index of the syn branch
  INTEGER,INTENT(IN) :: idanti            ! SAR index of the anti branch

  CHARACTER(LEN=LEN(group(1))) :: tgroup(SIZE(group))      
  CHARACTER(LEN=LEN(group(1))) :: pold,pnew,tempgr
  CHARACTER(LEN=LEN(chem))     :: p_hahp,p_carbonyl,p_acid
  REAL    :: brtio
  INTEGER :: np,tempring

  CHARACTER(LEN=LEN(idnam)) :: p(mxpd), r(3)
  REAL    :: s(mxpd), arrh(3)
  INTEGER :: idreac, nlabel
  REAL    :: xlabel,folow(3),fotroe(4)

  CHARACTER(LEN=15)     :: progname='rx_sci_water'
  CHARACTER(LEN=70)     :: mesg
  INTEGER,PARAMETER     :: mxcom=10
  INTEGER               :: nref
  CHARACTER(LEN=mxlcod) :: ref(mxcom)

  IF (screenfg) WRITE(6,*)progname

  tgroup(:)=group(:)
  nref=0 ; ref(:)=' ' ; CALL addref(progname,'MN22SAR2',nref,ref,chem)

! make the carbonyl
  p_carbonyl=' '
  IF     (group(cnod)=='CH.(OO.)') THEN ; pold='CH.(OO.)' ; pnew='CHO'
  ELSEIF (group(cnod)=='C.(OO.)' ) THEN ; pold='C.(OO.)'  ; pnew='CO'
  ELSE
    mesg="criegee group not identified "
    CALL stoperr(progname,mesg,chem)
  ENDIF
  tempgr=tgroup(cnod)
  CALL swap(tempgr,pold,tgroup(cnod),pnew)
  CALL rebond(bond,tgroup,p_carbonyl,tempring)
  CALL stdchm(p_carbonyl)
  tgroup(:)=group(:)

! make the carboxylic acid
  p_acid=' '
  IF (group(cnod)=='CH.(OO.)') THEN 
    pold='CH.(OO.)' ; pnew='CO(OH)'
    tempgr=tgroup(cnod)
    CALL swap(tempgr,pold,tgroup(cnod),pnew)
    CALL rebond(bond,tgroup,p_acid,tempring)
    CALL stdchm(p_acid)
    tgroup(:)=group(:)
  ENDIF

! make the hydroxy alkyl hydroperoxide (hahp)
  IF      (group(cnod)=='CH.(OO.)') THEN ; pold='CH.(OO.)' ; pnew='CH(OH)(OOH)'
  ELSE IF (group(cnod)=='C.(OO.)' ) THEN ; pold='C.(OO.)'  ; pnew='C(OH)(OOH)'
  ELSE
    mesg="criegee group not identified "
    CALL stoperr(progname,mesg,chem)
  ENDIF
  CALL swap(group(cnod),pold,tgroup(cnod),pnew)
  CALL rebond(bond,tgroup,p_hahp,tempring)
  CALL stdchm(p_hahp)
  tgroup(:)=group(:)

! set stoi. coef. and products (same for H2O monomer and dimer)

! reaction with water monomer
! ---------------------------
  CALL rxinit(r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe)
  r(1)=idnam  ;  r(2)='EXTRA '  ;  np=0  

  CALL add1tonp(progname,chem,np)  ; s(np)=0.55 ; brtio=brch*s(np)
  CALL bratio(p_hahp,brtio,p(np),nref,ref)
  IF (p_acid==' ') THEN
    CALL add1tonp(progname,chem,np)  ; s(np)=0.45 ; brtio=brch*s(np)
    CALL bratio(p_carbonyl,brtio,p(np),nref,ref)
    CALL add1tonp(progname,chem,np)  ; s(np)=0.45 ; p(np)='H2O2'
  ELSE
    CALL add1tonp(progname,chem,np)  ; s(np)=0.40 ; brtio=brch*s(np)
    CALL bratio(p_carbonyl,brtio,p(np),nref,ref)
    CALL add1tonp(progname,chem,np)  ; s(np)=0.40 ; p(np)='H2O2'
    CALL add1tonp(progname,chem,np)  ; s(np)=0.05 ; brtio=brch*s(np)
    CALL bratio(p_acid,brtio,p(np),nref,ref)
  ENDIF

  idreac=2  ;   nlabel=500
  arrh(:) = k_mh2o(idsyn,idanti,:)
  CALL rxwrit(mecu,r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe,ref)
  
! reaction with water dimer
! ---------------------------
  CALL rxinit(r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe)
  r(1)=idnam  ;  r(2)='EXTRA '  ; np=0

  CALL add1tonp(progname,chem,np)  ; s(np)=0.55 ; brtio=brch*s(np)
  CALL bratio(p_hahp,brtio,p(np),nref,ref)
  IF (p_acid==' ') THEN
    CALL add1tonp(progname,chem,np)  ; s(np)=0.45 ; brtio=brch*s(np)
    CALL bratio(p_carbonyl,brtio,p(np),nref,ref)
    CALL add1tonp(progname,chem,np)  ; s(np)=0.45 ; p(np)='H2O2'
  ELSE
    CALL add1tonp(progname,chem,np)  ; s(np)=0.40 ; brtio=brch*s(np)
    CALL bratio(p_carbonyl,brtio,p(np),nref,ref)
    CALL add1tonp(progname,chem,np)  ; s(np)=0.40 ; p(np)='H2O2'
    CALL add1tonp(progname,chem,np)  ; s(np)=0.05 ; brtio=brch*s(np)
    CALL bratio(p_acid,brtio,p(np),nref,ref)
  ENDIF

  idreac=2  ;  nlabel=502
  arrh(:) = k_dh2o(idsyn,idanti,:)
  CALL rxwrit(mecu,r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe,ref)

END SUBROUTINE rx_sci_water

! ======================================================================
! Purpose: perform the unimolecular decomposition of the stabilized
! criegee intermediate (SCI). Reaction to be performed is controled
! by the index of the syn node, following Vereecken et al. (2017) as
! described in Newland et al. 2022. The reaction is written to 
! the output file within the subroutine.
! ======================================================================
SUBROUTINE rx_sci_uni(idnam,chem,group,bond,ngr,brch,cnod, &
                      idsyn,synnod,ntreesyn,treesyn,idanti,antinod)
  USE criegeetool, ONLY: ci_unipdct,k_uni
  USE keyflag, ONLY: screenfg
  IMPLICIT NONE
  
  CHARACTER(LEN=*),INTENT(IN) :: idnam    ! name of current species
  CHARACTER(LEN=*),INTENT(IN) :: chem     ! formula of current species
  CHARACTER(LEN=*),INTENT(IN) :: group(:) ! group in chem
  INTEGER,INTENT(IN) :: bond(:,:)         ! bond matrix
  INTEGER,INTENT(IN) :: ngr               ! # of group/nod in chem
  REAL,INTENT(IN)    :: brch              ! max yield of the input species
  INTEGER,INTENT(IN) :: cnod              ! criegee node
  INTEGER,INTENT(IN) :: idsyn             ! SAR index of the syn branch
  INTEGER,INTENT(IN) :: synnod            ! syn node 
  INTEGER,INTENT(IN) :: ntreesyn(:)       ! # of distinct nodes at depth j (syn side)
  INTEGER,INTENT(IN) :: treesyn(:,:,:)    ! node tree starting from the syn node
  INTEGER,INTENT(IN) :: idanti            ! SAR index of the anti branch
  INTEGER,INTENT(IN) :: antinod           ! anti node

  INTEGER :: np
  REAL    :: yield 

  CHARACTER(LEN=LEN(idnam)) :: p(mxpd), r(3)
  REAL    :: s(mxpd), arrh(3)
  INTEGER :: idreac, nlabel
  REAL    :: xlabel,folow(3),fotroe(4)

  CHARACTER(LEN=15)     :: progname='rx_sci_uni'
  CHARACTER(LEN=70)     :: mesg
  INTEGER,PARAMETER     :: mxcom=10
  INTEGER               :: nref
  CHARACTER(LEN=mxlcod) :: ref(mxcom)

  IF (screenfg) WRITE(6,*)progname

  np=0
  CALL rxinit(r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe)
  nref=0 ; ref(:)=' ' ; CALL addref(progname,'MN22SAR2',nref,ref,chem)

! make unimolecular decomposition product (return as s & p)
  yield=1.
  CALL ci_unipdct(chem,group,bond,ngr,brch,yield,cnod, &
                  idsyn,synnod,ntreesyn,treesyn,idanti,antinod,&
                  np,s,p,nref,ref)

! Write reaction (if found) 
  IF (np==0) THEN
    mesg="No product found for the unimolecular decomposition of the SCI"
    CALL stoperr(progname,mesg,chem)
  ENDIF

  r(1)=idnam
  arrh(:) = k_uni(idsyn,idanti,:)
  CALL rxwrit(mecu,r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe,ref)
  
END SUBROUTINE rx_sci_uni

END MODULE stbcriegeechem
