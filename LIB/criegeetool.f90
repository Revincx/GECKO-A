MODULE criegeetool
USE toolbox, ONLY: stoperr,add1tonp,addref
USE reactool, ONLY: swap,rebond
IMPLICIT NONE

! initialization flag (raised after 1st call to "init_stab_criegee")
LOGICAL              :: init_done=.FALSE.  

! Criegee structures index
INTEGER,PARAMETER    :: idch2oo=1          ! CH2OO
INTEGER,PARAMETER    :: idrchoo=2          ! -CHOO : external CI
INTEGER,PARAMETER    :: idrrcoo=3          !  >COO : internal CI
INTEGER,PARAMETER    :: ncistruct=idrrcoo  ! # of CI structures

! -----------------------------------------
! Data arrays for Vereecken et al. 2017 SAR (unimolecular & water
! -----------------------------------------

! index for syn branches 
INTEGER,PARAMETER :: r1CH3 = 2          ! R1 => -CH3
INTEGER,PARAMETER :: r1CH2 = 3          ! R1 => -CH2Ra
INTEGER,PARAMETER :: r1CH  = 4          ! R1 => -CHRaRb
INTEGER,PARAMETER :: r1C = 5            ! R1 => -CRaRbRc
INTEGER,PARAMETER :: r1H = 1            ! R1 => -H
INTEGER,PARAMETER :: r1CH2Cd = 6        ! R1 => -CH2-CR3=CR4R5
INTEGER,PARAMETER :: r1CHCd = 7         ! R1 => -CHRa-CR3=CR4R5
INTEGER,PARAMETER :: r1CCd  = 8         ! R1 => -CRaRb-CR3=CR4R5
INTEGER,PARAMETER :: r1Cd2CH3 = 9       ! R1 => -(CR3=CR4CH3)
INTEGER,PARAMETER :: r1Cd2CH2 = 10      ! R1 => -CR3=CR4CH2Ra
INTEGER,PARAMETER :: r1Cd2CH = 11       ! R1 => -CR3=CR4CHRaRb
INTEGER,PARAMETER :: r1Cd2C = 12        ! R1 => -CR3=CR4R'
INTEGER,PARAMETER :: r1CHO = 13         ! R1 => -CHO
INTEGER,PARAMETER :: r1CO = 14          ! R1 => -C(O)Ra
INTEGER,PARAMETER :: r1asH = 15         ! R1 => -c and as H grp
INTEGER,PARAMETER :: nr1 = r1asH        ! # of id for syn branch

! index for anti branches
INTEGER,PARAMETER :: r2CH3 = 2          ! R2 = -CH3
INTEGER,PARAMETER :: r2CH2 = 3          ! R2 = -CH2Ra
INTEGER,PARAMETER :: r2CH  = 4          ! R2 = -CHRaRb
INTEGER,PARAMETER :: r2C = 5            ! R2 = -CRaRbRc
INTEGER,PARAMETER :: r2H = 1            ! R2 = -H
INTEGER,PARAMETER :: r2Cd2 = 6          ! R2 = vinyl (-CdR=CdRR)
INTEGER,PARAMETER :: r2CO = 7           ! R2 = oxo
INTEGER,PARAMETER :: r2asH = 8          ! R2 = -c and as H grp
INTEGER,PARAMETER :: nr2 = r2asH        ! # of id for anti branch

! index for various unimolecular decomposition pathways
INTEGER,PARAMETER :: id_13cycl=1, id_14Hshift=2, id_allyl16Hshift=3 
INTEGER,PARAMETER :: id_15cycl=4, id_allyl14Hshift=5               
INTEGER,PARAMETER :: n_uni_pathways=id_allyl14Hshift ! # of pathways

! rates for unimolecular decomposition, H2O monomer and dimer reactions
INTEGER,DIMENSION(nr1)    :: path_uni   ! pathway index for SCI decomposition
REAL,DIMENSION(nr1,nr2,3) :: k_uni      ! unimolecular rates
REAL,DIMENSION(nr1,nr2,3) :: k_mh2o     ! criegee + H2O rates
REAL,DIMENSION(nr1,nr2,3) :: k_dh2o     ! criegee + (H2O)2 rates

CONTAINS

! ======================================================================
! Purpose: perform the unimolecular decomposition of the criegee 
! intermediate (CI) following Vereecken et al. (2017) to identify the
! major decomposition pathway . Reaction to be performed is controled
! by the index of the syn node.
! ======================================================================
SUBROUTINE ci_unipdct(chem,group,bond,ngr,brch,yield,cnod, &
                      idsyn,synnod,ntreesyn,treesyn,idanti,antinod,&
                      np,s,p,nref,ref)
  USE toolbox, ONLY: kval
  USE fragmenttool, ONLY: fragm  
  USE ringtool, ONLY: ring_data
  USE dictstacktool, ONLY: bratio
  USE normchem, ONLY: stdchm
  USE cdtool, ONLY: switchenol

  IMPLICIT NONE
  
  CHARACTER(LEN=*),INTENT(IN) :: chem     ! formula of current species
  CHARACTER(LEN=*),INTENT(IN) :: group(:) ! group in chem
  INTEGER,INTENT(IN) :: bond(:,:)         ! bond matrix
  INTEGER,INTENT(IN) :: ngr               ! # of group/nod in chem
  REAL,INTENT(IN)    :: yield             ! yield of the current pathway
  REAL,INTENT(IN)    :: brch              ! max yield of the input species
  INTEGER,INTENT(IN) :: cnod              ! criegee node
  INTEGER,INTENT(IN) :: idsyn             ! SAR index of the syn branch
  INTEGER,INTENT(IN) :: synnod            ! syn node 
  INTEGER,INTENT(IN) :: ntreesyn(:)       ! # of distinct nodes at depth j (syn side)
  INTEGER,INTENT(IN) :: treesyn(:,:,:)    ! node tree starting from the syn node
  INTEGER,INTENT(IN) :: idanti            ! SAR index of the anti branch
  INTEGER,INTENT(IN) :: antinod           ! anti node
  INTEGER,INTENT(INOUT):: np              ! # of product (incremented here)
  REAL,INTENT(INOUT)   :: s(:)            ! stoi. coef. of the products in p
  CHARACTER(LEN=*),INTENT(INOUT):: p(:)   ! list of the products (short name)
  INTEGER,INTENT(INOUT):: nref            ! # of references
  CHARACTER(LEN=*),INTENT(INOUT):: ref(:) ! reac_info code 

  CHARACTER(LEN=LEN(group(1))) :: tgroup(SIZE(group))      
  CHARACTER(LEN=LEN(group(1))) :: pold,pnew
  CHARACTER(LEN=LEN(chem))     :: pchem,pchem2
  CHARACTER(LEN=LEN(chem))     :: fragpd(2)
  INTEGER :: tbond(SIZE(bond,1),SIZE(bond,2))

  INTEGER :: idpath                       ! id for the decomposition channel
  INTEGER :: nc,tempring,jj,nring
  INTEGER :: aa,bb,cc,ifrag,ia,wknod
  LOGICAL :: lodicarb,loester,lomolar
  REAL    :: brtio
  REAL    :: yroute, yall

  INTEGER,PARAMETER :: mxirg=6             ! max # of distinct rings 
  INTEGER  :: ndrg                         ! # of disctinct rings
  INTEGER  :: trackrg(mxirg,SIZE(tgroup))  ! (a,:)== track (node #) belonging ring a
  LOGICAL  :: lorgnod(mxirg,SIZE(tgroup))  ! (a,b)==true if node b belong to ring a

  CHARACTER(LEN=15)     :: progname='ci_unipdct'
  CHARACTER(LEN=70)     :: mesg

  tgroup(:)=group(:)  ;  tbond=bond(:,:)
  idpath=path_uni(idsyn)
  
! ------------
! 1-4 H SHIFT:  >CH-C.(OO.)- => >C(.)CO- + OH
! ------------
  IF ( (idpath==id_14Hshift) .OR. (idpath==id_allyl14Hshift) ) THEN
    IF     (group(synnod)(1:3)=='CH3') THEN ; pold='CH3' ; pnew='CH2'
    ELSEIF (group(synnod)(1:3)=='CH2') THEN ; pold='CH2' ; pnew='CH'
    ELSEIF (group(synnod)(1:3)=='CHO') THEN ; pold='CHO' ; pnew='CO'
    ELSEIF (group(synnod)(1:2)=='CH')  THEN ; pold='CH'  ; pnew='C'
    ELSE
      mesg= "--error-- Group not allowed in 1-4 H shift"
      CALL stoperr(progname,mesg,chem)
    ENDIF
    CALL swap(group(synnod),pold,tgroup(synnod),pnew)
    nc=INDEX(tgroup(synnod),' ')  ;  tgroup(synnod)(nc:nc)='.'
   
    IF      (group(cnod)=='CH.(OO.)') THEN ; pold='CH.(OO.)' ;  pnew='CHO'
    ELSE IF (group(cnod)=='C.(OO.)' ) THEN ; pold='C.(OO.)'  ;  pnew='CO'
    ENDIF
    CALL swap(group(cnod),pold,tgroup(cnod),pnew)

    CALL rebond(tbond,tgroup,pchem,tempring)
    CALL add2p(chem,pchem,yield,brch,np,s,p,nref,ref)
    CALL add1tonp(progname,chem,np) ; s(np)=yield ; p(np)='HO '           
    tgroup(:)=group(:)
    
! ------------
! 1-6 H SHIFT:  >CH-C=C-C.(OO.)- => >C=C-C(.)-CO- + OH
! ------------
  ELSE IF (idpath==id_allyl16Hshift) THEN  
    DO jj=1,ntreesyn(3)
      aa=treesyn(3,jj,1) ; bb=treesyn(3,jj,2) ; cc=treesyn(3,jj,3)
      IF (group(aa)(1:2)/='Cd') CYCLE
      IF (group(bb)(1:2)/='Cd') CYCLE
      IF (group(cc)(1:3)=='CHO') CYCLE        ! no info regarding aldehyde H
      IF (group(cc)(1:2)=='CO')  CYCLE        ! no H on this path
      IF (group(cc)(1:2)=='C(')  CYCLE        ! no H on this path
      IF (group(cc)(1:2)=='C ')  CYCLE        ! no H on this path
      IF (group(cc)(1:2)=='-O')  CYCLE        ! no H on this path
      IF (INDEX(group(cc),'(ONO2)')/=0) CYCLE ! do not form Cd(ONO2), structure not treated in GECKO-A
        
      IF      (group(cnod)=='CH.(OO.)') THEN ; pold='CH.(OO.)' ;  pnew='CHO'
      ELSE IF (group(cnod)=='C.(OO.)' ) THEN ; pold='C.(OO.)'  ;  pnew='CO'
      ENDIF
      CALL swap(group(cnod),pold,tgroup(cnod),pnew)

      pold='Cd'  ; pnew='C'    ; CALL swap(group(aa),pold,tgroup(aa),pnew)
      nc=INDEX(tgroup(aa),' ')  ;  tgroup(aa)(nc:nc)='.'

      IF      (group(cc)(1:3)=='CH3') THEN ; pold='CH3' ; pnew='CdH2'
      ELSE IF (group(cc)(1:3)=='CH2') THEN ; pold='CH2' ; pnew='CdH' 
      ELSE IF (group(cc)(1:2)=='CH')  THEN ; pold='CH'  ; pnew='Cd'
      ENDIF
      CALL swap(group(cc),pold,tgroup(cc),pnew)

      tbond(aa,bb)=1 ; tbond(bb,aa)=1 ; tbond(bb,cc)=2 ; tbond(cc,bb)=2
      CALL rebond(tbond,tgroup,pchem,tempring)
      CALL add2p(chem,pchem,yield,brch,np,s,p,nref,ref)
      CALL add1tonp(progname,chem,np) ; s(np)=yield ; p(np)='HO '           
      tgroup(:)=group(:) ; tbond(:,:)=bond(:,:)
    ENDDO

! ------------
! 1-5 RING CLOSURE:  >C=C-C.(OO.)- => >C1-C=C(-O-O1)R => products
! ------------
  ELSE IF (idpath==id_15cycl) THEN
    DO jj=1,ntreesyn(2)
      aa=treesyn(2,jj,1) ; bb=treesyn(2,jj,2) 
      IF (group(aa)(1:2)/='Cd') CYCLE
      IF (group(bb)(1:2)/='Cd') CYCLE
      lodicarb=.FALSE.

! make the carbonyl-epoxyde (Newland et al., 2022):  >C1-O-C1-CO-R
      pold='Cd' ; pnew='C' 
      CALL swap(group(bb),pold,tgroup(bb),pnew)
      CALL swap(group(aa),pold,tgroup(aa),pnew)
      tbond(aa,bb)=1 ; tbond(bb,aa)=1

      IF      (group(cnod)=='CH.(OO.)') THEN ; pold='CH.(OO.)' ;  pnew='CHO'
      ELSE IF (group(cnod)=='C.(OO.)' ) THEN ; pold='C.(OO.)'  ;  pnew='CO'
      ENDIF
      CALL swap(group(cnod),pold,tgroup(cnod),pnew)

      tgroup(ngr+1)='-O-'
      tbond(aa,ngr+1)=3 ; tbond(ngr+1,aa)=3
      tbond(bb,ngr+1)=3 ; tbond(ngr+1,bb)=3

      CALL rebond(tbond,tgroup,pchem,tempring)
      CALL stdchm(pchem)
      tgroup(:)=group(:) ; tbond(:,:)=bond(:,:)

! make the dicarbonyl (Newland et al., 2022):  -CO-CH-CO-R
      IF (group(bb)(1:3)=='CdH') THEN     ! check that H if available on beta Cd
        lodicarb=.TRUE.

        IF      (group(cnod)=='CH.(OO.)') THEN ; pold='CH.(OO.)' ;  pnew='CHO'
        ELSE IF (group(cnod)=='C.(OO.)' ) THEN ; pold='C.(OO.)'  ;  pnew='CO'
        ENDIF
        CALL swap(group(cnod),pold,tgroup(cnod),pnew)

        IF (group(bb)(1:4)=='CdH2') THEN ; pold='CdH2' ; pnew='CHO' 
        ELSE                             ; pold='CdH'  ; pnew='CO'
        ENDIF
        CALL swap(group(bb),pold,tgroup(bb),pnew)
        IF (group(aa)(1:3)=='CdH') THEN ; pold='CdH' ; pnew='CH2' 
        ELSE                            ; pold='Cd'  ; pnew='CH'
        ENDIF
        CALL swap(group(aa),pold,tgroup(aa),pnew)
        tbond(aa,bb)=1 ; tbond(bb,aa)=1

        CALL rebond(tbond,tgroup,pchem2,tempring)
        CALL stdchm(pchem2)
        tgroup(:)=group(:) ; tbond(:,:)=bond(:,:)
      ENDIF
      
! add the species in the stack and as reaction products
      CALL add1tonp(progname,chem,np)
      IF (lodicarb) THEN ; s(np)=yield*0.5 ; ELSE ;  s(np)=1. ; ENDIF
      brtio=brch*s(np)
      CALL bratio(pchem,brtio,p(np),nref,ref)
      IF (lodicarb) THEN
        CALL add1tonp(progname,chem,np) ; s(np)=yield*0.5
        brtio=brch*s(np)
        CALL bratio(pchem2,brtio,p(np),nref,ref)      
      ENDIF
    ENDDO

! ------------
! 1-3 RING CLOSURE:  RC.(OO.)R' -> R(COO)R' -> -> CO2+R+R'
! ------------
  ELSE IF  (idpath==id_13cycl) THEN       

! check whether ester (stabilization) or radicals (decomposition)  
    loester=.FALSE.                 ! default is decomposition 
    IF (ngr>8) THEN
      loester=.TRUE.                ! big molecule expected to stabilize
      DO ia=1,ngr
        IF (bond(cnod,ia)/=0) THEN
          IF (group(ia)(1:2)=='CO')  loester=.FALSE.
          IF (group(ia)(1:3)=='CHO') loester=.FALSE.
          IF (group(ia)(1:3)=='-O-') loester=.FALSE.
        ENDIF
      ENDDO
    ENDIF
    CALL ring_data(cnod,ngr,tbond,tgroup,ndrg,lorgnod,trackrg)
    nring=COUNT(lorgnod(:,cnod))    ! count # of rings involving node ia
    IF (nring/=0) loester=.TRUE.    ! ester channel for cyclic species (no diradical!)

! ester channel
! -------------
    IF (loester) THEN

! 1st ester break C-C bond between criegee and syn side
      IF  (group(cnod)=='CH.(OO.)') THEN ; pold='CH.(OO.)' 
        IF (idsyn==1) THEN ; pnew='CO(OH)' ; ELSE ; pnew='CHO' ; ENDIF
      ELSE IF (group(cnod)=='C.(OO.)' ) THEN ; pold='C.(OO.)'  ;  pnew='CO'
      ELSE
        mesg="Invalid criegee for ester channel"
        CALL stoperr(progname,mesg,chem)
      ENDIF
      CALL swap(group(cnod),pold,tgroup(cnod),pnew)
      IF (idsyn/=1) THEN                           ! ester does not exit
        tbond(cnod,synnod)=0 ; tbond(synnod,cnod)=0
        tgroup(ngr+1)='-O-'
        tbond(cnod,ngr+1)=3   ; tbond(ngr+1,cnod)=3
        tbond(synnod,ngr+1)=3 ; tbond(ngr+1,synnod)=3
      ENDIF
      CALL rebond(tbond,tgroup,pchem,tempring)
      CALL stdchm(pchem)

      CALL add1tonp(progname,chem,np)
      s(np)=yield*0.5 ; brtio=brch*s(np)
      CALL bratio(pchem,brtio,p(np),nref,ref)    ! add 1st ester to the stack
      tgroup(:)=group(:) ; tbond(:,:)=bond(:,:)
      
! 2nd ester break C-C bond between criegee and anti side
      IF  (group(cnod)=='CH.(OO.)') THEN ; pold='CH.(OO.)' 
        IF (idanti==1) THEN ; pnew='CO(OH)' ; ELSE ; pnew='CHO' ; ENDIF
      ELSE IF (group(cnod)=='C.(OO.)' ) THEN ; pold='C.(OO.)'  ;  pnew='CO'
      ELSE
        mesg="Invalid criegee for ester channel"
        CALL stoperr(progname,mesg,chem)
      ENDIF
      CALL swap(group(cnod),pold,tgroup(cnod),pnew)
      IF (idanti/=1) THEN                           ! ester does not exit
        tbond(cnod,antinod)=0 ; tbond(antinod,cnod)=0
        tgroup(ngr+1)='-O-'
        tbond(cnod,ngr+1)=3   ; tbond(ngr+1,cnod)=3
        tbond(antinod,ngr+1)=3 ; tbond(ngr+1,antinod)=3
      ENDIF
      CALL rebond(tbond,tgroup,pchem2,tempring)
      CALL stdchm(pchem2)

      CALL add1tonp(progname,chem,np)
      s(np)=yield*0.5 ; brtio=brch*s(np)
      CALL bratio(pchem2,brtio,p(np),nref,ref)   ! add 2nd ester to the stack
      tgroup(:)=group(:) ; tbond(:,:)=bond(:,:)

! decomposition of the "hot" ester
! --------------------------------
    ELSE 

! identify whether internal or external (wknod/=0) criegee 
      wknod=0
      IF      (idanti==1) THEN ; wknod=synnod   ! wknod is the 1st node ... 
      ELSE IF (idsyn==1)  THEN ; wknod=antinod  ! ... on the remaining C chain
      ENDIF

!---- EXTERNAL CRIEGEE:
      IF (wknod/=0) THEN
        tgroup(cnod)=' '
        tbond(cnod,wknod)=0 ; tbond(wknod,cnod)=0
        lomolar=.TRUE.    ! cancel molecular channel if create mess
        IF (INDEX(group(wknod),'O')/=0) THEN
          IF (group(wknod)(1:3)/='CO ') lomolar=.FALSE.
        ENDIF 
        
!       Route 1: => CO2 + R1. + HO2 (40 %) [or 100 % if molecular channel is cancel]
        IF (lomolar) THEN ; yroute=0.4 ; ELSE ; yroute=1. ; ENDIF
          yall=yroute*yield 
          nc=INDEX(tgroup(wknod),' ')  ;  tgroup(wknod)(nc:nc)='.'
          CALL rebond(tbond,tgroup,pchem,tempring)
          CALL add2p(chem,pchem,yall,brch,np,s,p,nref,ref)
          CALL add1tonp(progname,chem,np) ;  s(np)=yall ; p(np)='CO2 '
          CALL add1tonp(progname,chem,np) ;  s(np)=yall ; p(np)='HO2 ' 
     
!         Route 2: => CO2 + RH (20 %)
          IF (lomolar) THEN
          yroute=0.2   ;  yall=yroute*yield 
          tgroup(wknod)=group(wknod)
          IF      (group(wknod)(1:3)=='CH3') THEN ; pold='CH3' ; pnew='CH4'
          ELSE IF (group(wknod)(1:3)=='CH2') THEN ; pold='CH2' ; pnew='CH3'
          ELSE IF (group(wknod)(1:2)=='CH')  THEN ; pold='CH'  ; pnew='CH2'
          ELSE IF (group(wknod)(1:3)=='CO ') THEN ; pold='CO'  ; pnew='CHO'
          ELSE IF (group(wknod)(1:3)=='CdH') THEN ; pold='CdH' ; pnew='CdH2'
          ELSE IF (group(wknod)(1:2)=='Cd')  THEN ; pold='Cd'  ; pnew='CdH'
          ELSE IF (group(wknod)(1:1)=='C')   THEN ; pold='C'   ; pnew='CH' 
          ELSE IF (group(wknod)(1:1)=='c')   THEN ; pold='c'   ; pnew='cH' 
          ELSE
            mesg="Can not find appropriate group for 'CO2+RH' route hot ester channel"
            CALL stoperr(progname,mesg,chem)
          ENDIF
          CALL swap(group(wknod),pold,tgroup(wknod),pnew)
          CALL rebond(tbond,tgroup,pchem,tempring)
          CALL stdchm(pchem)
          CALL add1tonp(progname,chem,np) ; s(np)=yall ; brtio=brch*s(np)
          CALL bratio(pchem,brtio,p(np),nref,ref)   
          CALL add1tonp(progname,chem,np) ;  s(np)=yall ; p(np)='CO2 '

!         Route 3: => CO + R(OH) (40 %)
          yroute=0.4   ;  yall=yroute*yield 
          tgroup(wknod)=group(wknod)
          nc=INDEX(tgroup(wknod),' ')  ;  tgroup(wknod)(nc:nc+3)='(OH)'
          CALL rebond(tbond,tgroup,pchem,tempring)
          CALL stdchm(pchem)
          CALL add1tonp(progname,chem,np) ; s(np)=yall ; brtio=brch*s(np)
          CALL bratio(pchem,brtio,p(np),nref,ref)   
          CALL add1tonp(progname,chem,np) ;  s(np)=yall ; p(np)='CO '

        ENDIF
        tgroup(:)=group(:) ; tbond(:,:)=bond(:,:)
  
!---- INTERNAL CRIEGEE: 
      ELSE 

! route: =>   R1. + R2. + CO2   
        tgroup(cnod)=' '
        tbond(cnod,synnod)=0  ; tbond(synnod,cnod)=0
        tbond(antinod,cnod)=0 ; tbond(antinod,cnod)=0
        nc=INDEX(tgroup(synnod),' ')  ; tgroup(synnod)(nc:nc)='.'
        nc=INDEX(tgroup(antinod),' ') ; tgroup(antinod)(nc:nc)='.'
        CALL fragm(tbond,tgroup,fragpd(1),fragpd(2))

! loop over the 2 fragmentation products
        DO ifrag=1,2
          pchem=fragpd(ifrag)
          CALL add2p(chem,pchem,yield,brch,np,s,p,nref,ref)
        ENDDO
        CALL add1tonp(progname,chem,np) ;  s(np)=yield ; p(np)='CO2 '

        tgroup(:)=group(:) ; tbond(:,:)=bond(:,:)
      ENDIF
    ENDIF

! No path found - stop error
! ----------------------
  ELSE
    mesg="Invalid unimolecular decomposition pathway"
    CALL stoperr(progname,mesg,chem)
  ENDIF

END SUBROUTINE ci_unipdct


! ======================================================================
! PURPOSE: the routine alphao_criegee break the R-O-C.(OO.)-R  
! hot-criegee provided as input into R(O.)+CO2+R. and return the product
! and related stoi. coef.
! Note: crude approximation is performed here. The structure is not 
! considered in the Vereecken SAR and is considered as a special case.
! ======================================================================
SUBROUTINE alphao_criegee(xcri,group,bond,ngr,brch,cnod,noda,nodb, &
                          np,s,p,nref,ref)
  USE fragmenttool, ONLY: fragm  
  USE ringtool, ONLY: ring_data
  USE outtool, ONLY: wrt_dict
  IMPLICIT NONE

  CHARACTER(LEN=*),INTENT(IN) :: xcri     ! formula of the hot criegee
  CHARACTER(LEN=*),INTENT(IN) :: group(:) ! group in xcri
  INTEGER,INTENT(IN) :: bond(:,:)         ! bond matrix
  INTEGER,INTENT(IN) :: ngr               ! # of group/nod in xcri
  REAL,INTENT(IN)    :: brch              ! branching ratio of the hot criegee
  INTEGER,INTENT(IN) :: cnod              ! criegee node
  INTEGER,INTENT(IN) :: noda              ! node of the 1st branch (next to criegee) 
  INTEGER,INTENT(IN) :: nodb              ! node of the 2nd branch (next to criegee) 
  INTEGER,INTENT(INOUT) :: np             ! # of product in the s and p list
  REAL,INTENT(INOUT)    :: s(:)           ! stoi. coef. of product in p(:) list
  CHARACTER(LEN=*),INTENT(INOUT) :: p(:)  ! product list (short names)
  INTEGER,INTENT(INOUT) :: nref           ! # of references added in the reference list
  CHARACTER(LEN=*),INTENT(INOUT):: ref(:) ! list of references 

  CHARACTER(LEN=LEN(group(1))) :: tgroup(SIZE(group))      
  INTEGER :: tbond(SIZE(bond,1),SIZE(bond,2))
  CHARACTER(LEN=LEN(xcri))     :: fragpd(2),pchem
  INTEGER :: i,j,iao  
  INTEGER :: ifrag,nc,nring
  REAL    :: yield

  INTEGER,PARAMETER :: mxirg=6             ! max # of distinct rings 
  INTEGER  :: ndrg                         ! # of disctinct rings
  INTEGER  :: trackrg(mxirg,SIZE(tgroup))  ! (a,:)== track (node #) belonging ring a
  LOGICAL  :: lorgnod(mxirg,SIZE(tgroup))  ! (a,b)==true if node b belong to ring a

  CHARACTER(LEN=15),PARAMETER :: progname='alphao_criegee'
  CHARACTER(LEN=80) :: mesg

  tgroup(:)=group(:)  ;  tbond=bond(:,:)

! criegee should not be part of ring (no biradicals !)
  CALL ring_data(cnod,ngr,tbond,tgroup,ndrg,lorgnod,trackrg)
  nring=COUNT(lorgnod(:,cnod))    ! count # of rings involving node cnod
  IF (nring/=0) THEN ! ester channel (??) needed for cyclic species (no diradical!)
    RETURN
    !mesg= "--error-- for -O-C.(OO.)- structure in a ring - not accounted"
    !CALL wrt_dict()
    !CALL stoperr(progname,mesg,xcri)
  ENDIF

! external criegee R-O-CH.(OO.)
! ----------------
  IF (nodb==0) THEN
  
    yield=1. ! assume 100 % decomposition (overwrite sciyield, no stabilization !)
    IF (group(cnod)(1:8)=='CH.(OO.)') THEN
      CALL add1tonp(progname,xcri,np) ; s(np)=yield ; p(np)='CO2 '  ! add CO2
      CALL add1tonp(progname,xcri,np) ; s(np)=yield ; p(np)='HO2 '  ! add HO2
    ELSE IF (group(cnod)(1:12)=='C(OOH).(OO.)') THEN
      CALL add1tonp(progname,xcri,np) ; s(np)=yield ; p(np)='CO2 '  ! add CO2
      CALL add1tonp(progname,xcri,np) ; s(np)=yield ; p(np)='HO2 '  ! add HO2
    ELSE
      mesg= "--error-- chemistry must be added for external criegee -O-C.(OO.)"
      CALL stoperr(progname,mesg,xcri)
    ENDIF
    tgroup(cnod)=' ' ; tbond(cnod,noda)=0 ; tbond(noda,cnod)=0
  
    IF (tgroup(noda)/='-O-') THEN  
      mesg= "--error-- the expected external criegee -O-C.(OO.) not identified"
      CALL stoperr(progname,mesg,xcri)
    ENDIF
    iloop: DO i=1,ngr
      IF ((tbond(i,noda)/=0).AND.(i/=cnod)) THEN
         tbond(i,noda)=0 ; tbond(noda,i)=0 ; tgroup(noda)=' '
         nc=INDEX(tgroup(i),' ') ; tgroup(i)(nc:nc+3)='(O.)'

         DO j=1,ngr                           ! check for peroxide
           IF ((tbond(i,j)/=0).AND.(j/=noda)) THEN
             IF (group(i)=='-O-') THEN
               tbond(i,j)=0 ; tbond(j,i)=0 ; tgroup(i)=' '
               nc=INDEX(tgroup(j),' ') ; tgroup(j)(nc:nc+4)='(OO.)'
               EXIT iloop
             ENDIF
           ENDIF
         ENDDO
         
      ENDIF
    ENDDO iloop
    CALL rebond(tbond,tgroup,pchem,nring)
    CALL add2p(xcri,pchem,yield,brch,np,s,p,nref,ref)
    tgroup(:)=group(:) ; tbond(:,:)=bond(:,:)

! internal criegee R-O-C.(OO.)-R
! ----------------
  ELSE
  
    yield=1.  ! assume 100 % decomposition
    CALL add1tonp(progname,xcri,np) ; s(np)=yield ; p(np)='CO2 '  ! add CO2
    tgroup(cnod)=' '
  
    tbond(cnod,noda)=0 ; tbond(noda,cnod)=0
    tbond(cnod,nodb)=0 ; tbond(nodb,cnod)=0
    iao=0
    IF (tgroup(noda)=='-O-') THEN  
      DO i=1,ngr
        IF ((tbond(i,noda)/=0).AND.(i/=cnod)) THEN
          tbond(i,noda)=0 ; tbond(noda,i)=0 ; tgroup(noda)=' '
          nc=INDEX(tgroup(i),' ') ; tgroup(i)(nc:nc+3)='(O.)'
          iao=iao+1

          DO j=1,ngr      
            IF ((tbond(i,j)/=0).AND.(j/=noda)) THEN
              IF (group(i)=='-O-') THEN
                tbond(i,j)=0 ; tbond(j,i)=0 ; tgroup(i)=' '
                nc=INDEX(tgroup(j),' ') ; tgroup(j)(nc:nc+4)='(OO.)'
              ENDIF
            ENDIF
          ENDDO

        ENDIF
      ENDDO
    ELSE
      nc=INDEX(tgroup(noda),' ')  ; tgroup(noda)(nc:nc)='.'
    ENDIF
    
    IF (tgroup(nodb)=='-O-') THEN
      DO i=1,ngr
        IF ((tbond(i,nodb)/=0).AND.(i/=cnod)) THEN
          tbond(i,nodb)=0 ; tbond(nodb,i)=0 ; tgroup(nodb)=' '
          nc=INDEX(tgroup(i),' ') ; tgroup(i)(nc:nc+3)='(O.)'
          iao=iao+1

          DO j=1,ngr                           ! check for peroxide
            IF ((tbond(i,j)/=0).AND.(j/=nodb)) THEN
              IF (group(i)=='-O-') THEN
                tbond(i,j)=0 ; tbond(j,i)=0 ; tgroup(i)=' '
                nc=INDEX(tgroup(j),' ') ; tgroup(j)(nc:nc+4)='(OO.)'
              ENDIF
            ENDIF
          ENDDO

        ENDIF
      ENDDO
    ELSE
      nc=INDEX(tgroup(nodb),' ')  ; tgroup(nodb)(nc:nc)='.'
    ENDIF
    
    IF (iao==0) THEN
      mesg= "--error-- the expected internal criegee -O-C.(OO.)- not identified"
      CALL stoperr(progname,mesg,xcri)
    ENDIF
  
    CALL fragm(tbond,tgroup,fragpd(1),fragpd(2))

! loop over the 2 fragmentation products
    DO ifrag=1,2
      pchem=fragpd(ifrag)
      CALL add2p(xcri,pchem,yield,brch,np,s,p,nref,ref)
    ENDDO
    tgroup(:)=group(:) ; tbond(:,:)=bond(:,:)
    
  ENDIF

! add reac_info for hot criegee decomposition
  CALL addref(progname,'OCRIEGEE',nref,ref,xcri)

END SUBROUTINE alphao_criegee


! ======================================================================
! PURPOSE: the routine ADD2P checks the radical (chem) provided as 
! input, send the species to bratio (addition to the stack) and load the 
! short name species to the p(:) list, with the corresponding  
! stoichiometric coefficient in s(:). The coproducts of the species are 
! also added to the p(:) as well as the species that may be produced   
! after electron delocalisation.
! ======================================================================
SUBROUTINE add2p(xcri,chem,yield,brch,np,s,p,nref,ref)
  USE keyparameter, ONLY: mxcopd
  USE dictstacktool, ONLY: bratio
  USE radchktool, ONLY: radchk
  USE normchem, ONLY: stdchm      
  USE toolbox, ONLY: add1tonp
  IMPLICIT NONE

  CHARACTER(LEN=*),INTENT(IN) :: xcri ! parent species that produced chem
  CHARACTER(LEN=*),INTENT(IN) :: chem ! radical to be managed
  REAL,INTENT(IN)      :: brch        ! branching ratio of the parent species
  REAL,INTENT(IN)      :: yield       ! yield of the radical
  INTEGER,INTENT(INOUT):: np          ! # of product (incremented here)
  REAL,INTENT(INOUT)   :: s(:)        ! stoi. coef. of the products in p
  CHARACTER(LEN=*),INTENT(INOUT):: p(:)   ! list of the products (short name)
  INTEGER,INTENT(INOUT):: nref            ! # of references
  CHARACTER(LEN=*),INTENT(INOUT):: ref(:) ! reac_info code 
  
  CHARACTER(LEN=LEN(chem)) :: tempfo
  INTEGER :: j
  REAL :: brtio

  INTEGER,PARAMETER :: mxrpd=2  ! max # of products returned by radchk sub
  CHARACTER(LEN=LEN(chem))  :: rdckpd(mxrpd)
  CHARACTER(LEN=LEN(p(1))) :: rdckcopd(mxrpd,mxcopd)
  REAL :: sc(mxrpd)
  INTEGER :: nip
  CHARACTER(LEN=20),PARAMETER :: progname='add2p_from_xcrieg'

! send to radchk to get product and coproduct
  CALL radchk(chem,rdckpd,rdckcopd,nip,sc,nref,ref)

! 1st products from radchk 
  tempfo=rdckpd(1)
  CALL stdchm(tempfo)
  CALL add1tonp(progname,xcri,np)
  s(np)=yield*sc(1)
  brtio=brch*s(np)
  CALL bratio(tempfo,brtio,p(np),nref,ref)
  DO j=1,SIZE(rdckcopd,2)
    IF (rdckcopd(1,j)(1:1)==' ') CYCLE
    CALL add1tonp(progname,xcri,np)
    s(np)=yield*sc(1) ; p(np)=rdckcopd(1,j)
  ENDDO
  
! 2nd products from radchk 
  IF (nip==2) THEN
    tempfo=rdckpd(2)
    CALL stdchm(tempfo)
    CALL add1tonp(progname,xcri,np)
    s(np)=yield*sc(2)
    brtio=brch*s(np)
    CALL bratio(tempfo,brtio,p(np),nref,ref)
    DO j=1,SIZE(rdckcopd,2)
      IF (rdckcopd(2,j)(1:1)==' ') CYCLE
      CALL add1tonp(progname,xcri,np)
      s(np)=yield*sc(2) ; p(np)=rdckcopd(2,j)
    ENDDO
  ENDIF
END SUBROUTINE add2p

!=======================================================================
! PURPOSE: Return the index (as defined for the Vereecken 2017 SAR) for 
! the SYN branch bounded to the Criegee group.
!
! Branches considered for the syn configuration:
!  (1) r1H      => -H                  (2) r1CH3    => -CH3            
!  (3) r1CH2    => -CH2Ra              (4) r1CH     => -CHRaRb         
!  (5) r1C      => -CRaRbRc            (6) r1CH2Cd  => -CH2-CR3=CR4R5  
!  (7) r1CHCd   => -CHRa-CR3=CR4R5     (8) r1CCd    => -CRaRb-CR3=CR4R5
!  (9) r1Cd2CH3 => -(CR3=CR4CH3)      (10) r1Cd2CH2 => -CR3=CR4CH2Ra   
! (11) r1Cd2CH  => -CR3=CR4CHRaRb     (12) r1Cd2C   => -CR3=CR4R'      
! (13) r1CHO    => -CHO               (14) r1CO     => -C(O)Ra
! (15) r1asH    => -c and H like grp
! ======================================================================
SUBROUTINE syn_index(chem,group,bond,ring,ngr,icri,ia,ntreenod,treenod,brindex)
  USE keyparameter, ONLY: mxcp, mxring
  USE mapping, ONLY: abcde_map
  USE ringtool, ONLY: ring_data
  IMPLICIT NONE

  CHARACTER(LEN=*),INTENT(IN) :: chem     ! formula of current species
  CHARACTER(LEN=*),INTENT(IN) :: group(:) ! group in chem
  INTEGER,INTENT(IN) :: bond(:,:)         ! bond in chem (without break, for ring count)
  INTEGER,INTENT(IN) :: ring(:)           ! nodes belonging to ring
  INTEGER,INTENT(IN) :: ngr               ! # of node in chem
  INTEGER,INTENT(IN) :: icri              ! node bearing the criegee group
  INTEGER,INTENT(IN) :: ia                ! node next to criegee group (examined branch)  
  INTEGER,INTENT(IN) :: ntreenod(:)       ! # of distinct nodes at depth j
  INTEGER,INTENT(IN) :: treenod(:,:,:)    ! node tree starting from ia
  INTEGER,INTENT(OUT):: brindex           ! index for the branch bounded to the Criegee grp

  INTEGER,PARAMETER :: mxirg=6            ! max # of distinct rings 
  INTEGER  :: ndrg                        ! # of distinct rings 
  INTEGER  :: trackrg(mxirg,SIZE(group))  ! (a,:)== track (node #) belonging ring a
  LOGICAL  :: lorgnod(mxirg,SIZE(group))  ! (a,b)==true if node b belong to ring a

  INTEGER :: jt,nring
  LOGICAL :: notfound

  CHARACTER(LEN=24) :: progname='syn_index'
  CHARACTER(LEN=70) :: mesg

  brindex=0       ! a value must be set before returning
  
! check for simple "unique" structure
  IF      (group(ia)(1:3) == 'CH3') THEN ; brindex = r1CH3     ! -CH3
  ELSE IF (group(ia)(1:3) == 'CHO') THEN ; brindex = r1CHO     ! -CHO
  ELSE IF (group(ia)(1:2) == 'CO')  THEN ; brindex = r1CO      ! -CO-
  ELSE IF (group(ia)(1:1) == 'c')   THEN ; brindex = r1asH     ! aromatic ring 
 
  ELSE IF (group(ia)(1:3) == 'CH2') THEN ; brindex = r1CH2     ! -CH2R
    DO jt=1,ntreenod(2)   
      IF(group(treenod(2,jt,2))(1:2)=='Cd') brindex = r1CH2Cd  ! -CH2-Cd=Cd
    ENDDO
 
  ELSE IF (group(ia)(1:2) == 'CH')  THEN ; brindex = r1CH      ! -CH<
    DO jt=1,ntreenod(2)   
      IF(group(treenod(2,jt,2))(1:2)=='Cd') brindex = r1CHCd   ! -CH-Cd=Cd
    ENDDO    
    CALL ring_data(ia,ngr,bond,group,ndrg,lorgnod,trackrg)
    nring=COUNT(lorgnod(:,ia))    ! count # of rings involving node ia
    IF (nring>2)  brindex = r1C   ! No H shift if multiple rings on ia
    
! check -Cd=Cd-R
  ELSE IF (group(ia)(1:2) == 'Cd') THEN                        ! vinyl group
    notfound=.TRUE.

    IF (ring(icri)/=0) THEN
      IF (ring(ia)/=0) THEN
        DO jt=1,ntreenod(2)   
          IF (group(treenod(2,jt,2))(1:2)=='Cd') THEN
            IF (ring(treenod(2,jt,2))/=0) THEN
              brindex = r1asH                                  ! -C=C-C in ring
              notfound=.FALSE. ; EXIT
            ENDIF
          ENDIF
        ENDDO    
      ENDIF
    ENDIF

    IF (notfound) THEN
      DO jt=1,ntreenod(3)   ! scroll all groups for CH3
        IF (group(treenod(3,jt,3))(1:3)=='CH3') THEN
          IF (group(treenod(3,jt,2))(1:2)=='Cd') THEN
            brindex = r1Cd2CH3                                   ! -Cd=Cd-CH3
            notfound=.FALSE. ; EXIT
          ENDIF
        ENDIF
      ENDDO
    ENDIF
    
    IF (notfound) THEN
      DO jt=1,ntreenod(3)   ! scroll all groups for CH2
        IF ((group(treenod(3,jt,3))(1:3)=='CH2').AND.    &
            (INDEX(group(treenod(3,jt,3)),'(ONO2)')==0)) THEN  ! to avoid forming Cd(ONO2), structure not treated in GECKO-A
          IF (group(treenod(3,jt,2))(1:2)=='Cd') THEN
            brindex = r1Cd2CH2                                 ! -Cd=Cd-CH2
            notfound=.FALSE. ; EXIT
          ENDIF
        ENDIF
      ENDDO
    ENDIF
 
    IF (notfound) THEN
      DO jt=1,ntreenod(3)   ! scroll all groups for CH
        IF ((group(treenod(3,jt,3))(1:2)=='CH') .AND.    &
            (group(treenod(3,jt,3))(1:3)/='CHO').AND.    &     ! no info regarding aldehyde H
            (INDEX(group(treenod(3,jt,3)),'(ONO2)')==0)) THEN  ! to avoid forming Cd(ONO2), structure not treated in GECKO-A
          IF (group(treenod(3,jt,2))(1:2)=='Cd') THEN
            brindex = r1Cd2CH                                  ! -Cd=Cd-CH
            notfound=.FALSE. ; EXIT
          ENDIF
        ENDIF
      ENDDO
    ENDIF
           
    IF (notfound)   brindex = r1Cd2C                           ! -Cd=Cd

! check for -C
  ELSE IF (group(ia)(1:1) == 'C') THEN ;  brindex = r1C        ! -C
    DO jt=1,ntreenod(2)   
      IF(group(treenod(2,jt,2))(1:2)=='Cd') brindex = r1CCd    ! -C-Cd=Cd
    ENDDO
  ENDIF

! check that a value was set to the branch index
  IF (brindex==0) THEN
    mesg="no index found for a criegee branch "
    CALL stoperr(progname,mesg,chem)
  ENDIF 

END SUBROUTINE syn_index

!=======================================================================
! PURPOSE: Return the index (as defined for the Vereecken 2017 SAR) for 
! the ANTI branch bounded to the Criegee group.
!
! Branches considered for the anti configuration:
!  (1) r2H      => -H                  (2) r2CH3    => -CH3            
!  (3) r2CH2    => -CH2Ra              (4) r2CH     => -CHRaRb         
!  (5) r2C      => -CRaRbRc            (6) r2Cd2    => -CdR=CdRR 
!  (7) r2CO     => -CO-                (8) r2asH    => -c & as H grp
! ======================================================================
SUBROUTINE anti_index(chem,group,ia,brindex)
  USE keyparameter, ONLY: mxcp, mxring
  USE mapping, ONLY: abcde_map
  IMPLICIT NONE

  CHARACTER(LEN=*),INTENT(IN) :: chem     ! formula of current species
  CHARACTER(LEN=*),INTENT(IN) :: group(:) ! group in chem
  INTEGER,INTENT(IN) :: ia                ! node next to criegee group (examined branch)  
  INTEGER,INTENT(OUT):: brindex           ! index for the branch bounded to the Criegee grp

  CHARACTER(LEN=24) :: progname='anti_index'
  CHARACTER(LEN=70) :: mesg

  brindex=0  ! a value must be set before returning
  
! check for simple "unique" structure
  IF      (group(ia)(1:3) == 'CH3') THEN ; brindex = r2CH3      ! -CH3
  ELSE IF (group(ia)(1:3) == 'CH2') THEN ; brindex = r2CH2      ! -CH2R
  ELSE IF (group(ia)(1:3) == 'CHO') THEN ; brindex = r2CO       ! -CHO
  ELSE IF (group(ia)(1:2) == 'CO')  THEN ; brindex = r2CO       ! -CO-
  ELSE IF (group(ia)(1:2) == 'CH')  THEN ; brindex = r2CH       ! -CH<
  ELSE IF (group(ia)(1:2) == 'Cd')  THEN ; brindex = r2Cd2      ! -Cd=Cd
  ELSE IF (group(ia)(1:1) == 'C')   THEN ; brindex = r2C        ! -C
  ELSE IF (group(ia)(1:1) == 'c')   THEN ; brindex = r2asH      ! -c & as H grp
  ENDIF
  
! check that a value was set to the branch index
  IF (brindex==0) THEN
    mesg="no index found for a criegee anti branch "
    CALL stoperr(progname,mesg,chem)
  ENDIF 

END SUBROUTINE anti_index

! ======================================================================
! PURPOSE: given 2 branches "a" and "b" provided as input, identify
! which is branch is syn, which is anti (and return the corresponding
! trees).
! ======================================================================
SUBROUTINE get_synanti(chem,zecase,ntreea,treea,ntreeb,treeb,cipa,cipb, &
                       synnod,antinod,ntreesyn,treesyn,ntreeanti,treeanti)
  IMPLICIT NONE

  CHARACTER(LEN=*),INTENT(IN) :: chem   ! current species         
  INTEGER,INTENT(IN)  :: zecase         ! 1 = Z(cis) config ; 2 = E(trans) config         
  INTEGER,INTENT(IN)  :: ntreea(:)      ! # of distinct nodes at depth j (a side)
  INTEGER,INTENT(IN)  :: treea(:,:,:)   ! node tree starting from node a
  INTEGER,INTENT(IN)  :: ntreeb(:)      ! # of distinct nodes at depth j (b side)
  INTEGER,INTENT(IN)  :: treeb(:,:,:)   ! node tree starting from node b
  INTEGER,INTENT(IN)  :: cipa(:,:)      ! CIP tree starting from 1st atom (a side)
  INTEGER,INTENT(IN)  :: cipb(:,:)      ! CIP tree starting from 1st atom (b side)
  INTEGER,INTENT(OUT) :: synnod         ! syn node 
  INTEGER,INTENT(OUT) :: antinod        ! anti node
  INTEGER,INTENT(OUT) :: ntreesyn(:)    ! # of distinct nodes at depth j (syn side)
  INTEGER,INTENT(OUT) :: treesyn(:,:,:) ! node tree starting from the syn node
  INTEGER,INTENT(OUT) :: ntreeanti(:)   ! # of distinct nodes at depth j (anti side)
  INTEGER,INTENT(OUT) :: treeanti(:,:,:)! node tree starting from the anti node

  INTEGER :: i,j
  LOGICAL :: lofound 
  CHARACTER(LEN=15) :: progname='get_synanti'
  CHARACTER(LEN=70) :: mesg

  lofound=.FALSE.
 
  linelp: DO i=1,SIZE(cipa,1)      ! scroll the CIP table line by line
    IF (cipa(i,1)==0) EXIT linelp  ! nothing left to explore
    poslp: DO j=1,SIZE(cipa,2)     ! scroll element in each line
      IF (cipa(i,j)/=cipb(i,j)) THEN     ! find a difference in CIP numbers

        IF (cipa(i,j)>cipb(i,j)) THEN    ! rank of "a" branch is greater than "b" 
          IF (zecase==1) THEN                  ! syn case: "a" is syn & "b" is anti 
            synnod=treea(1,1,1)  ; antinod=treeb(1,1,1)
            ntreesyn=ntreea  ; treesyn=treea
            ntreeanti=ntreeb ; treeanti=treeb 
          ELSE                           ! anti case: "a" is anti & "b" is syn
            synnod=treeb(1,1,1)  ; antinod=treea(1,1,1)
            ntreesyn=ntreeb  ; treesyn=treeb
            ntreeanti=ntreea ; treeanti=treea 
          ENDIF
        ELSE                             ! rank of "b" branch is greater than "a" 
          IF (zecase==1) THEN                  ! syn case: "b" is syn & "a" is anti
            synnod=treeb(1,1,1)  ; antinod=treea(1,1,1)
            ntreesyn=ntreeb  ; treesyn=treeb
            ntreeanti=ntreea ; treeanti=treea 
          ELSE                           ! anti case: "b" is anti & "a" is syn
            synnod=treea(1,1,1)  ; antinod=treeb(1,1,1)
            ntreesyn=ntreea  ; treesyn=treea
            ntreeanti=ntreeb ; treeanti=treeb 
          ENDIF
        ENDIF
        lofound=.TRUE.
        EXIT linelp
      ENDIF
      IF (cipa(i,j)==0) CYCLE linelp ! scroll next line
    ENDDO poslp
  ENDDO linelp

  IF (.NOT. lofound) THEN
    mesg="Unexpected identical CIP rank found for 2 criegee branches"
    CALL stoperr(progname,mesg,chem)
  ENDIF
  
END SUBROUTINE get_synanti

! ======================================================================
! PURPOSE: Check whether 2 CIP trees provided are identical (return 0)
! or not (return 1)
! ======================================================================
INTEGER FUNCTION check_ze(cipa,cipb)
  IMPLICIT NONE
  INTEGER,INTENT(IN)  :: cipa(:,:)      ! CIP tree starting from 1st atom (a side)
  INTEGER,INTENT(IN)  :: cipb(:,:)      ! CIP tree starting from 1st atom (b side)
  
  INTEGER :: i,j
  
  check_ze=0

  linelp: DO i=1,SIZE(cipa,1)      ! scroll the CIP table line by line
    IF (cipa(i,1)==0) EXIT linelp  ! nothing left to explore
    poslp: DO j=1,SIZE(cipa,2)     ! scroll element in each line
      IF (cipa(i,j)/=cipb(i,j)) THEN     ! find a difference in CIP numbers
        check_ze=1 ; RETURN
      ENDIF
      IF (cipa(i,j)==0) CYCLE linelp ! scroll next line
    ENDDO poslp
  ENDDO linelp
  
END FUNCTION check_ze

! ======================================================================
! PURPOSE: return the stabilized criegee (SCI) yield for the hot 
! criegee (CI*) provided as input.
! ======================================================================
REAL FUNCTION sci_yield(chem,pozatom,idsyn)
  USE atomtool, ONLY:getatoms
  IMPLICIT NONE
  CHARACTER(LEN=*),INTENT(IN) :: chem   ! current species (i.e. the CI*)         
  INTEGER,INTENT(IN) :: pozatom         ! number O, C, N atoms in parent POZ
  INTEGER,INTENT(IN) :: idsyn           ! SAR index of the syn branch

  INTEGER :: idpath                     ! id for the decomposition channel
  INTEGER :: ciatom                     ! number O, C, N atoms in CI
  REAL    :: f_zpath                    ! F*Z_pathway
  INTEGER :: xxc,xxh,xxn,xxo,xxr,xxs,xxfl,xxbr,xxcl

  CHARACTER(LEN=15) :: progname='sci_yield'
  CHARACTER(LEN=70) :: mesg

! compute the number of atoms (C,N,O) in CI 
  CALL getatoms(chem,xxc,xxh,xxn,xxo,xxr,xxs,xxfl,xxbr,xxcl)
  ciatom=xxc+xxn+xxo   

! compute F*z_path
  idpath=path_uni(idsyn)
  IF  (idpath==id_13cycl) THEN  ; f_zpath=0.950*(5./REAL(ciatom+2)) ! 5-3 is 2 (see paper) 
  ELSE                          ; f_zpath=1.242*(5./REAL(ciatom))   ! 5-5 is 0 (see paper)
  ENDIF

! compute SCI yield
  sci_yield=1.-(f_zpath*ciatom/pozatom)
  IF (sci_yield>1.) THEN
    mesg="sci yield is greater than 1."
    CALL stoperr(progname,mesg,chem)
  ENDIF
  IF (sci_yield<0.) sci_yield=0.   ! this happen in CdH2=CdHCH3
  
END FUNCTION sci_yield

! ======================================================================
! PURPOSE: "init_stab_criegee" initializes the data used for the 
! stabilized criegee reactions.
! ======================================================================
SUBROUTINE init_stab_criegee()
  IMPLICIT NONE

! ------------------------------------------------------
! RATE FOR THE VEREECKEN SAR: H2O monomer, H2O dimer & unimolecular
! ------------------------------------------------------

! rate constant for unimolecular decomposition

! r1CH3 = 2       ! R1 => -CH3
path_uni(r1CH3) = id_14Hshift
k_uni(r1CH3, r2CH3, :)  = (/ 7.64E-60 , 23.59 , -2367. /)  ! R2 = -CH3
k_uni(r1CH3, r2CH2, :)  = (/ 8.63E-61 , 23.94 , -2390. /)  ! R2 = -CH2Ra
k_uni(r1CH3, r2CH , :)  = (/ 4.18E-59 , 23.38 , -2276. /)  ! R2 = -CHRaRb
k_uni(r1CH3, r2C  , :)  = (/ 9.26E-61 , 23.92 , -2362. /)  ! R2 = -CRaRbRc
k_uni(r1CH3, r2H  , :)  = (/ 3.11E-64 , 24.95 , -2685. /)  ! R2 = -H
k_uni(r1CH3, r2Cd2, :)  = (/ 4.36E-67 , 25.90 , -2737. /)  ! R2 = vinyl (-CdR=CdRR)
k_uni(r1CH3, r2CO , :)  = (/ 7.31E-80 , 30.08 , -3644. /)  ! R2 = oxo
k_uni(r1CH3, r2asH, :)  = (/ 3.11E-64 , 24.95 , -2685. /)  ! R2 = -c & H like grp (use H)
                                  
! r1CH2 = 3       ! R1 => -CH2Ra
path_uni(r1CH2) = id_14Hshift                   
k_uni(r1CH2, r2CH3, :)  = (/ 5.91E-58 , 22.91 , -2331. /)  ! R2 = -CH3
k_uni(r1CH2, r2CH2, :)  = (/ 4.72E-59 , 23.29 , -2358. /)  ! R2 = -CH2Ra
k_uni(r1CH2, r2CH , :)  = (/ 4.83E-61 , 23.93 , -2420. /)  ! R2 = -CHRaRb
k_uni(r1CH2, r2C  , :)  = (/ 2.13E-74 , 28.47 , -2964. /)  ! R2 = -CRaRbRc
k_uni(r1CH2, r2H  , :)  = (/ 2.41E-62 , 24.33 , -2571. /)  ! R2 = -H
k_uni(r1CH2, r2Cd2, :)  = (/ 1.42E-66 , 25.62 , -2780. /)  ! R2 = vinyl (-CdR=CdRR)
k_uni(r1CH2, r2CO , :)  = (/ 2.38E-79 , 29.80 , -3686. /)  ! R2 = oxo
k_uni(r1CH2, r2asH, :)  = (/ 2.41E-62 , 24.33 , -2571. /)  ! R2 = -c & H like grp (use H)
                                  
! r1CH  = 4       ! R1 => -CHRaRb
path_uni(r1CH) = id_14Hshift                   
k_uni(r1CH, r2CH3, :)  = (/ 4.45E-76 , 27.97 , -2967. /)  ! R2 = -CH3
k_uni(r1CH, r2CH2, :)  = (/ 8.54E-73 , 27.89 , -2858. /)  ! R2 = -CH2Ra
k_uni(r1CH, r2CH , :)  = (/ 1.62E-70 , 27.01 , -2912. /)  ! R2 = -CHRaRb
k_uni(r1CH, r2C  , :)  = (/ 8.80E-81 , 30.47 , -3249. /)  ! R2 = -CRaRbRc
!k_uni(r1CH, r2H  , :)  = (/ 8.54E-63 , 24.47 , -2490. /)  ! R2 = -H   ! published values are errors
k_uni(r1CH, r2H  , :)  = (/ 1.52E-75 , 28.72 , -3150. /)  ! R2 = -H
k_uni(r1CH, r2Cd2, :)  = (/ 2.55E-80 , 30.21 , -3387. /)  ! R2 = vinyl (-CdR=CdRR)
k_uni(r1CH, r2CO , :)  = (/ 4.28E-93 , 34.39 , -4294. /)  ! R2 = oxo
k_uni(r1CH, r2asH, :)  = (/ 1.52E-75 , 28.72 , -3150. /)  ! R2 = -c & H like grp (use H)
                                  
! r1C = 5         ! R1 => -CRaRbRc
path_uni(r1C) = id_13cycl                   
k_uni(r1C, r2CH3, :)  = (/ 2.55E10 , 1.02  ,  9399. /)  ! R2 = -CH3
k_uni(r1C, r2CH2, :)  = (/ 5.45E11 , 0.55  ,  9555. /)  ! R2 = -CH2Ra
k_uni(r1C, r2CH , :)  = (/ 6.83E10 , 0.80  ,  8759. /)  ! R2 = -CHRaRb
k_uni(r1C, r2C  , :)  = (/ 1.12E11 , 0.85  ,  9795. /)  ! R2 = -CRaRbRc
k_uni(r1C, r2H  , :)  = (/ 2.58E6  , 2.32  ,  9710. /)  ! R2 = -H
k_uni(r1C, r2Cd2, :)  = (/ 3.76E7  , 1.87  ,  8977. /)  ! R2 = vinyl (-CdR=CdRR)
k_uni(r1C, r2CO , :)  = (/ 1.92E2  , 3.63  , 10331. /)  ! R2 = oxo
k_uni(r1C, r2asH, :)  = (/ 2.58E6  , 2.32  ,  9710. /)  ! R2 = -c & H like grp (use H)
                                  
! r1H = 1         ! R1 => -H
path_uni(r1H) = id_13cycl                   
k_uni(r1H, r2CH3, :)  = (/ 1.69E9  , 1.35  ,  7445. /)  ! R2 = -CH3
k_uni(r1H, r2CH2, :)  = (/ 1.57E10 , 1.03  ,  7464. /)  ! R2 = -CH2Ra
k_uni(r1H, r2CH , :)  = (/ 9.22E9  , 1.13  ,  7387. /)  ! R2 = -CHRaRb
k_uni(r1H, r2C  , :)  = (/ 8.51E9  , 1.15  ,  7357. /)  ! R2 = -CRaRbRc
k_uni(r1H, r2H  , :)  = (/ 1.66E1  , 4.02  ,  8024. /)  ! R2 = -H
k_uni(r1H, r2Cd2, :)  = (/ 1.68E10 , 1.02  ,  7732. /)  ! R2 = vinyl (-CdR=CdRR)
k_uni(r1H, r2CO , :)  = (/ 8.58E4  , 2.78  ,  9085. /)  ! R2 = oxo
k_uni(r1H, r2asH, :)  = (/ 1.66E1  , 4.02  ,  8024. /)  ! R2 = -c & H like grp (use H)
                                  
! r1CH2Cd = 6     ! R1 => -CH2-CR3=CR4R5
path_uni(r1CH2Cd) = id_allyl14Hshift                   
k_uni(r1CH2Cd, r2CH3, :)  = (/ 9.00E-42 , 17.45 , -1390. /)  ! R2 = -CH3
k_uni(r1CH2Cd, r2CH2, :)  = (/ 7.18E-43 , 17.83 , -1417. /)  ! R2 = -CH2Ra
k_uni(r1CH2Cd, r2CH , :)  = (/ 7.35E-45 , 18.47 , -1479. /)  ! R2 = -CHRaRb
k_uni(r1CH2Cd, r2C  , :)  = (/ 3.24E-58 , 23.01 , -2023. /)  ! R2 = -CRaRbRc
k_uni(r1CH2Cd, r2H  , :)  = (/ 3.66E-46 , 18.87 , -1630. /)  ! R2 = -H
k_uni(r1CH2Cd, r2Cd2, :)  = (/ 2.16E-50 , 20.16 , -1839. /)  ! R2 = vinyl (-CdR=CdRR)
k_uni(r1CH2Cd, r2CO , :)  = (/ 3.63E-63 , 24.35 , -2746. /)  ! R2 = oxo
k_uni(r1CH2Cd, r2asH, :)  = (/ 3.66E-46 , 18.87 , -1630. /)  ! R2 = -c & H like grp (use H)
                                  
! r1CHCd = 7      ! R1 => -CHRa-CR3=CR4R5
path_uni(r1CHCd) = id_allyl14Hshift                   
k_uni(r1CHCd, r2CH3, :) = (/ 2.60E-55 , 21.92 , -2087. /)  ! R2 = -CH3
k_uni(r1CHCd, r2CH2, :) = (/ 4.99E-55 , 21.85 , -1978. /)  ! R2 = -CH2Ra
k_uni(r1CHCd, r2CH , :) = (/ 9.47E-53 , 20.97 , -2032. /)  ! R2 = -CHRaRb
k_uni(r1CHCd, r2C  , :) = (/ 5.14E-63 , 24.42 , -2369. /)  ! R2 = -CRaRbRc
k_uni(r1CHCd, r2H  , :) = (/ 8.89E-58 , 22.68 , -2269. /)  ! R2 = -H
k_uni(r1CHCd, r2Cd2, :) = (/ 1.49E-62 , 24.16 , -2507. /)  ! R2 = vinyl (-CdR=CdRR)
k_uni(r1CHCd, r2CO , :) = (/ 2.50E-75 , 28.34 , -3413. /)  ! R2 = oxo
k_uni(r1CHCd, r2asH, :) = (/ 8.89E-58 , 22.68 , -2269. /)  ! R2 = -c & H like grp (use H)
                                  
! r1CCd  = 8      ! R1 => -CRaRb-CR3=CR4R5
path_uni(r1CCd) = id_13cycl                   
k_uni(r1CCd, r2CH3, :)  = (/ 9.73E10 , 0.79 ,  9314. /)  ! R2 = -CH3
k_uni(r1CCd, r2CH2, :)  = (/ 2.08E12 , 0.32 ,  9471. /)  ! R2 = -CH2Ra
k_uni(r1CCd, r2CH , :)  = (/ 2.61E11 , 0.57 ,  8674. /)  ! R2 = -CHRaRb
k_uni(r1CCd, r2C  , :)  = (/ 4.27E11 , 0.62 ,  9711. /)  ! R2 = -CRaRbRc
k_uni(r1CCd, r2H  , :)  = (/ 9.87E6  , 2.09 ,  9625. /)  ! R2 = -H
k_uni(r1CCd, r2Cd2, :)  = (/ 1.44E8  , 1.63 ,  8893. /)  ! R2 = vinyl (-CdR=CdRR)
k_uni(r1CCd, r2CO , :)  = (/ 7.32E2  , 3.39 , 10246. /)  ! R2 = oxo
k_uni(r1CCd, r2asH, :)  = (/ 9.87E6  , 2.09 ,  9625. /)  ! R2 = -c & H like grp (use H)
                                  
!r1Cd2CH3 = 9    ! R1 => -(CR3=CR4CH3)
path_uni(r1Cd2CH3) = id_allyl16Hshift                   
k_uni(r1Cd2CH3, r2CH3, :)  = (/ 1.30E-5  , 5.52  , 1106. /)  ! R2 = -CH3
k_uni(r1Cd2CH3, r2CH2, :)  = (/ 2.63E-6  , 5.77  , 1117. /)  ! R2 = -CH2Ra
k_uni(r1Cd2CH3, r2CH , :)  = (/ 2.31E-5  , 5.42  , 1131. /)  ! R2 = -CHRaRb
k_uni(r1Cd2CH3, r2C  , :)  = (/ 5.83E-8  , 6.24  , 1063. /)  ! R2 = -CRaRbRc
k_uni(r1Cd2CH3, r2H  , :)  = (/ 1.40E-9  , 6.77  ,  838. /)  ! R2 = -H
k_uni(r1Cd2CH3, r2Cd2, :)  = (/ 2.79E-12 , 7.67  ,  807. /)  ! R2 = vinyl (-CdR=CdRR)
k_uni(r1Cd2CH3, r2CO , :)  = (/ 3.29E-25 , 11.90 , -121. /)  ! R2 = oxo
k_uni(r1Cd2CH3, r2asH, :)  = (/ 1.40E-9  , 6.77  ,  838. /)  ! R2 = -c & H like grp (use H)
                                  
! r1Cd2CH2 = 10   ! R1 => -CR3=CR4CH2Ra
path_uni(r1Cd2CH2) = id_allyl16Hshift                   
k_uni(r1Cd2CH2, r2CH3, :)  = (/ 1.03E-4  , 5.22  , 1163. /)  ! R2 = -CH3
k_uni(r1Cd2CH2, r2CH2, :)  = (/ 2.07E-5  , 5.48  , 1174. /)  ! R2 = -CH2Ra
k_uni(r1Cd2CH2, r2CH , :)  = (/ 1.82E-4  , 5.12  , 1188. /)  ! R2 = -CHRaRb
k_uni(r1Cd2CH2, r2C  , :)  = (/ 4.60E-7  , 5.95  , 1120. /)  ! R2 = -CRaRbRc
k_uni(r1Cd2CH2, r2H  , :)  = (/ 1.57E-8  , 6.42  ,  917. /)  ! R2 = -H
k_uni(r1Cd2CH2, r2Cd2, :)  = (/ 2.20E-11 , 7.37  ,  864. /)  ! R2 = vinyl (-CdR=CdRR)
k_uni(r1Cd2CH2, r2CO , :)  = (/ 3.69E-24 , 11.55 ,  -42. /)  ! R2 = oxo
k_uni(r1Cd2CH2, r2asH, :)  = (/ 1.57E-8  , 6.42  ,  917. /)  ! R2 = -c & H like grp (use H)
                                  
! r1Cd2CH = 11    ! R1 => -CR3=CR4CHRaRb
path_uni(r1Cd2CH) = id_allyl16Hshift                   
k_uni(r1Cd2CH, r2CH3, :)  = (/ 9.92E-6  , 5.62  , 1456. /)  ! R2 = -CH3
k_uni(r1Cd2CH, r2CH2, :)  = (/ 2.00E-6  , 5.88  , 1467. /)  ! R2 = -CH2Ra
k_uni(r1Cd2CH, r2CH , :)  = (/ 1.76E-5  , 5.53  , 1481. /)  ! R2 = -CHRaRb
k_uni(r1Cd2CH, r2C  , :)  = (/ 4.44E-8  , 6.35  , 1413. /)  ! R2 = -CRaRbRc
k_uni(r1Cd2CH, r2H  , :)  = (/ 1.51E-9  , 6.83  , 1210. /)  ! R2 = -H
k_uni(r1Cd2CH, r2Cd2, :)  = (/ 2.13E-12 , 7.77  , 1157. /)  ! R2 = vinyl (-CdR=CdRR)
k_uni(r1Cd2CH, r2CO , :)  = (/ 3.56E-25 , 11.95 ,  251. /)  ! R2 = oxo
k_uni(r1Cd2CH, r2asH, :)  = (/ 1.51E-9  , 6.83  , 1210. /)  ! R2 = -c & H like grp (use H)

! r1Cd2C = 12     ! R1 => -CR3=CR4R'
path_uni(r1Cd2C) = id_15cycl                   
k_uni(r1Cd2C, r2CH3, :)  = (/ 9.75E8 , 1.03 , 5081. /)  ! R2 = -CH3       ! original values corrrected 
k_uni(r1Cd2C, r2CH2, :)  = (/ 9.75E8 , 1.03 , 5081. /)  ! R2 = -CH2Ra     ! original values corrrected
k_uni(r1Cd2C, r2CH , :)  = (/ 9.75E8 , 1.03 , 5081. /)  ! R2 = -CHRaRb    ! original values corrrected
k_uni(r1Cd2C, r2C  , :)  = (/ 9.75E8 , 1.03 , 5081. /)  ! R2 = -CRaRbRc   ! original values corrrected
k_uni(r1Cd2C, r2H  , :)  = (/ 2.58E9 , 0.87 , 5090. /)  ! R2 = -H
k_uni(r1Cd2C, r2Cd2, :)  = (/ 9.75E8 , 1.03 , 5081. /)  ! R2 = vinyl (-CdR=CdRR)   ! values corrrected
k_uni(r1Cd2C, r2CO , :)  = (/ 9.75E8 , 1.03 , 5081. /)  ! R2 = oxo        ! original values corrrected
k_uni(r1Cd2C, r2asH, :)  = (/ 2.58E9 , 0.87 , 5090. /)  ! R2 = -c & H like grp (use H)

! r1CHO = 13       ! R1 => -CHO
path_uni(r1CHO) = id_13cycl                   
k_uni(r1CHO, r2CH3, :)  = (/ 2.76E10 , 0.78 , 5162. /)  ! R2 = -CH3
k_uni(r1CHO, r2CH2, :)  = (/ 2.01E10 , 0.82 , 5212. /)  ! R2 = -CH2Ra
k_uni(r1CHO, r2CH , :)  = (/ 7.28E6  , 1.84 , 4329. /)  ! R2 = -CHRaRb
k_uni(r1CHO, r2C  , :)  = (/ 2.95E6  , 1.96 , 5099. /)  ! R2 = -CRaRbRc
k_uni(r1CHO, r2H  , :)  = (/ 2.49E6  , 2.21 , 6656. /)  ! R2 = -H
k_uni(r1CHO, r2Cd2, :)  = (/ 2.33E10 , 0.71 , 5303. /)  ! R2 = vinyl (-CdR=CdRR)
k_uni(r1CHO, r2CO , :)  = (/ 4.06E4  , 2.70 , 6730. /)  ! R2 = oxo
k_uni(r1CHO, r2asH, :)  = (/ 2.49E6  , 2.21 , 6656. /)  ! R2 = -c & H like grp (use H)
       
! r1CO = 14        ! R1 => -C(O)Ra
path_uni(r1CO) = id_13cycl                   
k_uni(r1CO, r2CH3, :)  = (/ 3.95E9 , 1.02 , 5448. /)  ! R2 = -CH3
k_uni(r1CO, r2CH2, :)  = (/ 2.86E9 , 1.05 , 5489. /)  ! R2 = -CH2Ra
k_uni(r1CO, r2CH , :)  = (/ 1.50E6 , 2.02 , 4620. /)  ! R2 = -CHRaRb
k_uni(r1CO, r2C  , :)  = (/ 1.29E5 , 2.36 , 5315. /)  ! R2 = -CRaRbRc
k_uni(r1CO, r2H  , :)  = (/ 5.55E5 , 2.41 , 7141. /)  ! R2 = -H
k_uni(r1CO, r2Cd2, :)  = (/ 2.25E9 , 1.00 , 5560. /)  ! R2 = vinyl (-CdR=CdRR)
k_uni(r1CO, r2CO , :)  = (/ 3.92E3 , 2.99 , 6987. /)  ! R2 = oxo
k_uni(r1CO, r2asH, :) = (/ 5.55E5 , 2.41 , 7141. /)  ! R2 = -c & H like grp (use H)

! r1asH = 15        ! R1 => -c & other H like (Use -H as surrogate)
path_uni(r1asH) = id_13cycl                   
k_uni(r1asH, r2CH3, :)  = (/ 1.69E9  , 1.35  ,  7445. /)  ! R2 = -CH3
k_uni(r1asH, r2CH2, :)  = (/ 1.57E10 , 1.03  ,  7464. /)  ! R2 = -CH2Ra
k_uni(r1asH, r2CH , :)  = (/ 9.22E9  , 1.13  ,  7387. /)  ! R2 = -CHRaRb
k_uni(r1asH, r2C  , :)  = (/ 8.51E9  , 1.15  ,  7357. /)  ! R2 = -CRaRbRc
k_uni(r1asH, r2H  , :)  = (/ 1.66E1  , 4.02  ,  8024. /)  ! R2 = -H
k_uni(r1asH, r2Cd2, :)  = (/ 1.68E10 , 1.02  ,  7732. /)  ! R2 = vinyl (-CdR=CdRR)
k_uni(r1asH, r2CO , :)  = (/ 8.58E4  , 2.78  ,  9085. /)  ! R2 = oxo
k_uni(r1asH, r2asH, :)  = (/ 1.66E1  , 4.02  ,  8024. /)  ! R2 = -c & H like grp (use H)


! -----------------------------
! Rate constant with water monomer
! -----------------------------

! r1CH3 = 2       ! R1 => -CH3
k_mh2o(r1CH3, r2CH3, :)  = (/ 3.87E-20 , 1.91 , 1677. /)    ! R2 = -CH3
k_mh2o(r1CH3, r2CH2, :)  = (/ 1.41E-20 , 2.03 , 1524. /)    ! R2 = -CH2Ra
k_mh2o(r1CH3, r2CH , :)  = (/ 2.06E-19 , 1.60 , 1813. /)    ! R2 = -CHRaRb
k_mh2o(r1CH3, r2C  , :)  = (/ 4.94E-19 , 1.44 , 1174. /)    ! R2 = -CRaRbRc
k_mh2o(r1CH3, r2H  , :)  = (/ 2.19E-19 , 1.68 , 2513. /)    ! R2 = -H
k_mh2o(r1CH3, r2Cd2, :)  = (/ 7.07E-19 , 1.46 , 3132. /)    ! R2 = vinyl (-CdR=CdRR)
k_mh2o(r1CH3, r2CO , :)  = (/ 1.35E-22 , 2.75 , 3393. /)    ! R2 = oxo
k_mh2o(r1CH3, r2asH, :)  = (/ 2.19E-19 , 1.68 , 2513. /)    ! R2 = -c & H like grp (use H)
                  
! r1CH2 = 3       ! R1 => -CH2Ra                            
k_mh2o(r1CH2, r2CH3, :)  = (/ 4.19E-20 , 1.78 , 1265. /)    ! R2 = -CH3
k_mh2o(r1CH2, r2CH2, :)  = (/ 2.05E-19 , 1.45 ,  964. /)    ! R2 = -CH2Ra
k_mh2o(r1CH2, r2CH , :)  = (/ 4.21E-20 , 1.75 ,  970. /)    ! R2 = -CHRaRb
k_mh2o(r1CH2, r2C  , :)  = (/ 2.76E-18 , 1.15 ,  961. /)    ! R2 = -CRaRbRc
k_mh2o(r1CH2, r2H  , :)  = (/ 1.23E-18 , 1.39 , 2300. /)    ! R2 = -H
k_mh2o(r1CH2, r2Cd2, :)  = (/ 1.76E-19 , 1.62 , 2697. /)    ! R2 = vinyl (-CdR=CdRR)
k_mh2o(r1CH2, r2CO , :)  = (/ 1.51E-22 , 2.71 , 3227. /)    ! R2 = oxo
k_mh2o(r1CH2, r2asH, :)  = (/ 1.23E-18 , 1.39 , 2300. /)    ! R2 = -c & H like grp (use H)
                  
! r1CH  = 4       ! R1 => -CHRaRb                           
k_mh2o(r1CH, r2CH3, :)  = (/ 1.66E-20 , 1.71 , 1005. /)     ! R2 = -CH3
k_mh2o(r1CH, r2CH2, :)  = (/ 8.76E-20 , 1.36 ,  806. /)     ! R2 = -CH2Ra
k_mh2o(r1CH, r2CH , :)  = (/ 1.80E-20 , 1.66 ,  812. /)     ! R2 = -CHRaRb
k_mh2o(r1CH, r2C  , :)  = (/ 1.18E-18 , 1.06 ,  802. /)     ! R2 = -CRaRbRc
k_mh2o(r1CH, r2H  , :)  = (/ 5.76E-18 , 1.11 , 1740. /)     ! R2 = -H
k_mh2o(r1CH, r2Cd2, :)  = (/ 8.28E-19 , 1.34 , 2136. /)     ! R2 = vinyl (-CdR=CdRR)
k_mh2o(r1CH, r2CO , :)  = (/ 7.10E-22 , 2.43 , 2666. /)     ! R2 = oxo
k_mh2o(r1CH, r2asH, :)  = (/ 5.76E-18 , 1.11 , 1740. /)     ! R2 = -c & H like grp (use H)
                  
! r1C = 5         ! R1 => -CRaRbRc                          
k_mh2o(r1C, r2CH3, :)   = (/ 7.64E-21 , 1.98 , 1104. /)     ! R2 = -CH3
k_mh2o(r1C, r2CH2, :)   = (/ 4.02E-20 , 1.63 ,  905. /)     ! R2 = -CH2Ra
k_mh2o(r1C, r2CH , :)   = (/ 8.24E-21 , 1.93 ,  911. /)     ! R2 = -CHRaRb
k_mh2o(r1C, r2C  , :)   = (/ 5.40E-19 , 1.33 ,  902. /)     ! R2 = -CRaRbRc
k_mh2o(r1C, r2H  , :)   = (/ 2.40E-19 , 1.57 , 2241. /)     ! R2 = -H
k_mh2o(r1C, r2Cd2, :)   = (/ 3.45E-20 , 1.80 , 2638. /)     ! R2 = vinyl (-CdR=CdRR)
k_mh2o(r1C, r2CO , :)   = (/ 2.96E-23 , 2.89 , 3168. /)     ! R2 = oxo
k_mh2o(r1C, r2asH, :)   = (/ 2.40E-19 , 1.57 , 2241. /)     ! R2 = -c & H like grp (use H)
                  
! r1H = 1         ! R1 => -H                                
k_mh2o(r1H, r2CH3, :)   = (/ 4.32E-18 , 1.27 , -405. /)     ! R2 = -CH3
k_mh2o(r1H, r2CH2, :)   = (/ 9.64E-18 , 0.92 , -652. /)     ! R2 = -CH2Ra
k_mh2o(r1H, r2CH , :)   = (/ 6.69E-18 , 1.16 , -581. /)     ! R2 = -CHRaRb
k_mh2o(r1H, r2C  , :)   = (/ 1.83E-17 , 0.99 , -641. /)     ! R2 = -CRaRbRc
k_mh2o(r1H, r2H  , :)   = (/ 8.13E-18 , 1.23 ,  698. /)     ! R2 = -H
k_mh2o(r1H, r2Cd2, :)   = (/ 2.93E-19 , 1.66 ,  973. /)     ! R2 = vinyl (-CdR=CdRR)
k_mh2o(r1H, r2CO , :)   = (/ 1.59E-18 , 1.53 , 2258. /)     ! R2 = oxo
k_mh2o(r1H, r2asH, :)   = (/ 8.13E-18 , 1.23 ,  698. /)     ! R2 = -c & H like grp (use H)
                         
! r1CH2Cd = 6     ! R1 => -CH2-CR3=CR4R5
k_mh2o(r1CH2Cd, r2CH3, :)  = (/ 8.45E-21 , 2.05 , 1633. /)  ! R2 = -CH3
k_mh2o(r1CH2Cd, r2CH2, :)  = (/ 8.43E-21 , 1.89 ,  784. /)  ! R2 = -CH2Ra
k_mh2o(r1CH2Cd, r2CH , :)  = (/ 1.73E-21 , 2.20 ,  790. /)  ! R2 = -CHRaRb
k_mh2o(r1CH2Cd, r2C  , :)  = (/ 1.13E-19 , 1.60 ,  781. /)  ! R2 = -CRaRbRc
k_mh2o(r1CH2Cd, r2H  , :)  = (/ 5.03E-20 , 1.84 , 2120. /)  ! R2 = -H
k_mh2o(r1CH2Cd, r2Cd2, :)  = (/ 7.23E-21 , 2.07 , 2516. /)  ! R2 = vinyl (-CdR=CdRR)
k_mh2o(r1CH2Cd, r2CO , :)  = (/ 6.20E-24 , 3.16 , 3047. /)  ! R2 = oxo
k_mh2o(r1CH2Cd, r2asH, :)  = (/ 5.03E-20 , 1.84 , 2120. /)  ! R2 = -c & H like grp (use H)
                         
! r1CHCd = 7      ! R1 => -CHRa-CR3=CR4R5
k_mh2o(r1CHCd, r2CH3, :)  = (/ 4.48E-23 , 2.48 ,  749. /)   ! R2 = -CH3
k_mh2o(r1CHCd, r2CH2, :)  = (/ 2.36E-22 , 2.13 ,  550. /)   ! R2 = -CH2Ra
k_mh2o(r1CHCd, r2CH , :)  = (/ 4.84E-23 , 2.44 ,  556. /)   ! R2 = -CHRaRb
k_mh2o(r1CHCd, r2C  , :)  = (/ 3.17E-21 , 1.84 ,  547. /)   ! R2 = -CRaRbRc
k_mh2o(r1CHCd, r2H  , :)  = (/ 1.55E-20 , 1.88 , 1484. /)   ! R2 = -H
k_mh2o(r1CHCd, r2Cd2, :)  = (/ 2.23E-21 , 2.11 , 1881. /)   ! R2 = vinyl (-CdR=CdRR)
k_mh2o(r1CHCd, r2CO , :)  = (/ 1.91E-24 , 3.20 , 2411. /)   ! R2 = oxo
k_mh2o(r1CHCd, r2asH, :)  = (/ 1.55E-20 , 1.88 , 1484. /)   ! R2 = -c & H like grp (use H)
                         
! r1CCd  = 8      ! R1 => -CRaRb-CR3=CR4R5
k_mh2o(r1CCd, r2CH3, :)  = (/ 2.25E-20 , 1.79 , 1013. /)    ! R2 = -CH3
k_mh2o(r1CCd, r2CH2, :)  = (/ 1.18E-19 , 1.44 ,  813. /)    ! R2 = -CH2Ra
k_mh2o(r1CCd, r2CH , :)  = (/ 2.43E-20 , 1.75 ,  819. /)    ! R2 = -CHRaRb
k_mh2o(r1CCd, r2C  , :)  = (/ 1.59E-18 , 1.15 ,  810. /)    ! R2 = -CRaRbRc
k_mh2o(r1CCd, r2H  , :)  = (/ 7.07E-19 , 1.39 , 2150. /)    ! R2 = -H
k_mh2o(r1CCd, r2Cd2, :)  = (/ 1.02E-19 , 1.62 , 2546. /)    ! R2 = vinyl (-CdR=CdRR)
k_mh2o(r1CCd, r2CO , :)  = (/ 8.71E-23 , 2.71 , 3076. /)    ! R2 = oxo
k_mh2o(r1CCd, r2asH, :)  = (/ 7.07E-19 , 1.39 , 2150. /)    ! R2 = -c & H like grp (use H)
                         
! r1Cd2CH3 = 9    ! R1 => -(CR3=CR4CH3)
k_mh2o(r1Cd2CH3, r2CH3, :)  = (/ 2.21E-21 , 2.27 , 1858. /) ! R2 = -CH3
k_mh2o(r1Cd2CH3, r2CH2, :)  = (/ 2.13E-20 , 1.83 , 1589. /) ! R2 = -CH2Ra
k_mh2o(r1Cd2CH3, r2CH , :)  = (/ 4.37E-21 , 2.14 , 1595. /) ! R2 = -CHRaRb
k_mh2o(r1Cd2CH3, r2C  , :)  = (/ 2.87E-19 , 1.4  , 1586. /) ! R2 = -CRaRbRc
k_mh2o(r1Cd2CH3, r2H  , :)  = (/ 2.24E-19 , 1.65 , 2989. /) ! R2 = -H
k_mh2o(r1Cd2CH3, r2Cd2, :)  = (/ 3.22E-20 , 1.88 , 3385. /) ! R2 = vinyl (-CdR=CdRR)
k_mh2o(r1Cd2CH3, r2CO , :)  = (/ 2.76E-23 , 2.97 , 3915. /) ! R2 = oxo
k_mh2o(r1Cd2CH3, r2asH, :)  = (/ 2.24E-19 , 1.65 , 2989. /) ! R2 = -c & H like grp (use H)
                         
! r1Cd2CH2 = 10   ! R1 => -CR3=CR4CH2Ra
k_mh2o(r1Cd2CH2, r2CH3, :)  = (/ 2.21E-21 , 2.27 , 1858. /) ! R2 = -CH3
k_mh2o(r1Cd2CH2, r2CH2, :)  = (/ 2.13E-20 , 1.83 , 1589. /) ! R2 = -CH2Ra
k_mh2o(r1Cd2CH2, r2CH , :)  = (/ 4.37E-21 , 2.14 , 1595. /) ! R2 = -CHRaRb
k_mh2o(r1Cd2CH2, r2C  , :)  = (/ 2.87E-19 , 1.4  , 1586. /) ! R2 = -CRaRbRc
k_mh2o(r1Cd2CH2, r2H  , :)  = (/ 2.24E-19 , 1.65 , 2989. /) ! R2 = -H
k_mh2o(r1Cd2CH2, r2Cd2, :)  = (/ 3.22E-20 , 1.88 , 3385. /) ! R2 = vinyl (-CdR=CdRR)
k_mh2o(r1Cd2CH2, r2CO , :)  = (/ 2.76E-23 , 2.97 , 3915. /) ! R2 = oxo
k_mh2o(r1Cd2CH2, r2asH, :)  = (/ 2.24E-19 , 1.65 , 2989. /) ! R2 = -c & H like grp (use H)
                         
! r1Cd2CH = 11    ! R1 => -CR3=CR4CHRaRb
k_mh2o(r1Cd2CH, r2CH3, :)  = (/ 2.21E-21 , 2.27 , 1858. /)  ! R2 = -CH3
k_mh2o(r1Cd2CH, r2CH2, :)  = (/ 2.13E-20 , 1.83 , 1589. /)  ! R2 = -CH2Ra
k_mh2o(r1Cd2CH, r2CH , :)  = (/ 4.37E-21 , 2.14 , 1595. /)  ! R2 = -CHRaRb
k_mh2o(r1Cd2CH, r2C  , :)  = (/ 2.87E-19 , 1.4  , 1586. /)  ! R2 = -CRaRbRc
k_mh2o(r1Cd2CH, r2H  , :)  = (/ 2.24E-19 , 1.65 , 2989. /)  ! R2 = -H
k_mh2o(r1Cd2CH, r2Cd2, :)  = (/ 3.22E-20 , 1.88 , 3385. /)  ! R2 = vinyl (-CdR=CdRR)
k_mh2o(r1Cd2CH, r2CO , :)  = (/ 2.76E-23 , 2.97 , 3915. /)  ! R2 = oxo
k_mh2o(r1Cd2CH, r2asH, :)  = (/ 2.24E-19 , 1.65 , 2989. /)  ! R2 = -c & H like grp (use H)

! r1Cd2C = 12     ! R1 => -CR3=CR4R'
k_mh2o(r1Cd2C, r2CH3, :)  = (/ 2.21E-21 , 2.27 , 1858. /)   ! R2 = -CH3
k_mh2o(r1Cd2C, r2CH2, :)  = (/ 2.13E-20 , 1.83 , 1589. /)   ! R2 = -CH2Ra
k_mh2o(r1Cd2C, r2CH , :)  = (/ 4.37E-21 , 2.14 , 1595. /)   ! R2 = -CHRaRb
k_mh2o(r1Cd2C, r2C  , :)  = (/ 2.87E-19 , 1.4  , 1586. /)   ! R2 = -CRaRbRc
k_mh2o(r1Cd2C, r2H  , :)  = (/ 2.24E-19 , 1.65 , 2989. /)   ! R2 = -H
k_mh2o(r1Cd2C, r2Cd2, :)  = (/ 3.22E-20 , 1.88 , 3385. /)   ! R2 = vinyl (-CdR=CdRR)
k_mh2o(r1Cd2C, r2CO , :)  = (/ 2.76E-23 , 2.97 , 3915. /)   ! R2 = oxo
k_mh2o(r1Cd2C, r2asH, :)  = (/ 2.24E-19 , 1.65 , 2989. /)   ! R2 = -c & H like grp (use H)
                         
! r1CHO = 13       ! R1 => -CHO
k_mh2o(r1CHO, r2CH3, :)  = (/ 7.81E-20 , 1.68 , -757. /)    ! R2 = -CH3
k_mh2o(r1CHO, r2CH2, :)  = (/ 5.10E-20 , 1.65 , -966. /)    ! R2 = -CH2Ra
k_mh2o(r1CHO, r2CH , :)  = (/ 1.05E-20 , 1.95 , -960. /)    ! R2 = -CHRaRb
k_mh2o(r1CHO, r2C  , :)  = (/ 6.86E-19 , 1.35 , -969. /)    ! R2 = -CRaRbRc
k_mh2o(r1CHO, r2H  , :)  = (/ 3.04E-19 , 1.59 ,  370. /)    ! R2 = -H
k_mh2o(r1CHO, r2Cd2, :)  = (/ 4.37E-20 , 1.82 ,  767. /)    ! R2 = vinyl (-CdR=CdRR)
k_mh2o(r1CHO, r2CO , :)  = (/ 3.75E-23 , 2.91 , 1297. /)    ! R2 = oxo
k_mh2o(r1CHO, r2asH, :)  = (/ 3.04E-19 , 1.59 ,  370. /)    ! R2 = -c & H like grp (use H)
                  
! r1CO = 14        ! R1 => -C(O)Ra                          
k_mh2o(r1CO, r2CH3, :)  = (/ 2.18E-19 , 1.43 , -1268. /)    ! R2 = -CH3
k_mh2o(r1CO, r2CH2, :)  = (/ 1.15E-18 , 1.08 , -1468. /)    ! R2 = -CH2Ra
k_mh2o(r1CO, r2CH , :)  = (/ 2.35E-19 , 1.39 , -1462. /)    ! R2 = -CHRaRb
k_mh2o(r1CO, r2C  , :)  = (/ 1.54E-17 , 0.79 , -1471. /)    ! R2 = -CRaRbRc
k_mh2o(r1CO, r2H  , :)  = (/ 6.84E-18 , 1.03 ,  -131. /)    ! R2 = -H
k_mh2o(r1CO, r2Cd2, :)  = (/ 9.83E-19 , 1.26 ,   265. /)    ! R2 = vinyl (-CdR=CdRR)
k_mh2o(r1CO, r2CO , :)  = (/ 8.43E-22 , 2.35 ,   795. /)    ! R2 = oxo
k_mh2o(r1CO, r2asH, :)  = (/ 6.84E-18 , 1.03 ,  -131. /)    ! R2 = -c & H like grp (use H)

! r1asH= 15         ! R1 => -c & other H like (Use -H as surrogate)                                
k_mh2o(r1asH, r2CH3, :)   = (/ 4.32E-18 , 1.27 , -405. /)     ! R2 = -CH3
k_mh2o(r1asH, r2CH2, :)   = (/ 9.64E-18 , 0.92 , -652. /)     ! R2 = -CH2Ra
k_mh2o(r1asH, r2CH , :)   = (/ 6.69E-18 , 1.16 , -581. /)     ! R2 = -CHRaRb
k_mh2o(r1asH, r2C  , :)   = (/ 1.83E-17 , 0.99 , -641. /)     ! R2 = -CRaRbRc
k_mh2o(r1asH, r2H  , :)   = (/ 8.13E-18 , 1.23 ,  698. /)     ! R2 = -H
k_mh2o(r1asH, r2Cd2, :)   = (/ 2.93E-19 , 1.66 ,  973. /)     ! R2 = vinyl (-CdR=CdRR)
k_mh2o(r1asH, r2CO , :)   = (/ 1.59E-18 , 1.53 , 2258. /)     ! R2 = oxo
k_mh2o(r1asH, r2asH, :)   = (/ 8.13E-18 , 1.23 ,  698. /)     ! R2 = -c & H like grp (use H)

! -----------------------------
! Rate constant with water dimer
! -----------------------------

! r1CH3 = 2       ! R1 => -CH3
k_dh2o(r1CH3, r2CH3, :)  = (/ 3.90E-20 , 1.91 ,  -643. /)  ! R2 = -CH3   
k_dh2o(r1CH3, r2CH2, :)  = (/ 1.41E-20 , 2.03 ,  -786. /)  ! R2 = -CH2Ra   
k_dh2o(r1CH3, r2CH , :)  = (/ 2.05E-19 , 1.60 ,  -512. /)  ! R2 = -CHRaRb   
k_dh2o(r1CH3, r2C  , :)  = (/ 5.17E-19 , 1.43 , -1085. /)  ! R2 = -CRaRbRc  
k_dh2o(r1CH3, r2H  , :)  = (/ 2.27E-19 , 1.67 ,   121. /)  ! R2 = -H    
k_dh2o(r1CH3, r2Cd2, :)  = (/ 7.63E-19 , 1.45 ,   675. /)  ! R2 = vinyl (-CdR=CdRR)    
k_dh2o(r1CH3, r2CO , :)  = (/ 2.28E-22 , 2.68 ,   929. /)  ! R2 = oxo    
k_dh2o(r1CH3, r2asH, :)  = (/ 2.27E-19 , 1.67 ,   121. /)  ! R2 = -c & H like grp (use H)    
                        
! r1CH2 = 3       ! R1 => -CH2Ra
k_dh2o(r1CH2, r2CH3, :)  = (/ 4.27E-20 , 1.77 , -1022. /)  ! R2 = -CH3  
k_dh2o(r1CH2, r2CH2, :)  = (/ 2.35E-19 , 1.43 , -1286. /)  ! R2 = -CH2Ra  
k_dh2o(r1CH2, r2CH , :)  = (/ 5.70E-20 , 1.72 , -1255. /)  ! R2 = -CHRaRb  
k_dh2o(r1CH2, r2C  , :)  = (/ 2.88E-18 , 1.15 , -1274. /)  ! R2 = -CRaRbRc  
k_dh2o(r1CH2, r2H  , :)  = (/ 1.26E-18 , 1.39 ,   -68. /)  ! R2 = -H    
k_dh2o(r1CH2, r2Cd2, :)  = (/ 2.50E-19 , 1.57 ,   292. /)  ! R2 = vinyl (-CdR=CdRR)    
k_dh2o(r1CH2, r2CO , :)  = (/ 4.79E-22 , 2.55 ,   808. /)  ! R2 = oxo    
k_dh2o(r1CH2, r2asH, :)  = (/ 1.26E-18 , 1.39 ,   -68. /)  ! R2 = -c & H like grp (use H)    
                        
! r1CH  = 4       ! R1 => -CHRaRb
k_dh2o(r1CH, r2CH3, :)  = (/ 1.66E-20 , 1.71 , -1254. /)  ! R2 = -CH3  
k_dh2o(r1CH, r2CH2, :)  = (/ 8.47E-20 , 1.36 , -1441. /)  ! R2 = -CH2Ra  
k_dh2o(r1CH, r2CH , :)  = (/ 2.06E-20 , 1.65 , -1410. /)  ! R2 = -CHRaRb  
k_dh2o(r1CH, r2C  , :)  = (/ 1.04E-18 , 1.08 , -1429. /)  ! R2 = -CRaRbRc  
k_dh2o(r1CH, r2H  , :)  = (/ 5.93E-18 , 1.10 ,  -568. /)  ! R2 = -H   
k_dh2o(r1CH, r2Cd2, :)  = (/ 1.17E-18 , 1.29 ,  -207. /)  ! R2 = vinyl (-CdR=CdRR)   
k_dh2o(r1CH, r2CO , :)  = (/ 2.25E-21 , 2.27 ,   309. /)  ! R2 = oxo    
k_dh2o(r1CH, r2asH, :)  = (/ 5.93E-18 , 1.10 ,  -568. /)  ! R2 = -c & H like grp (use H)  
                        
! r1C = 5         ! R1 => -CRaRbRc
k_dh2o(r1C, r2CH3, :)  = (/ 8.81E-21 , 1.96 , -1155. /)  ! R2 = -CH3  
k_dh2o(r1C, r2CH2, :)  = (/ 4.49E-20 , 1.61 , -1342. /)  ! R2 = -CH2Ra  
k_dh2o(r1C, r2CH , :)  = (/ 1.09E-20 , 1.90 , -1311. /)  ! R2 = -CHRaRb  
k_dh2o(r1C, r2C  , :)  = (/ 5.51E-19 , 1.33 , -1329. /)  ! R2 = -CRaRbRc  
k_dh2o(r1C, r2H  , :)  = (/ 2.42E-19 , 1.57 ,  -124. /)  ! R2 = -H   
k_dh2o(r1C, r2Cd2, :)  = (/ 4.78E-20 , 1.76 ,   237. /)  ! R2 = vinyl (-CdR=CdRR)    
k_dh2o(r1C, r2CO , :)  = (/ 9.17E-23 , 2.74 ,   753. /)  ! R2 = oxo    
k_dh2o(r1C, r2asH, :)  = (/ 2.42E-19 , 1.57 ,  -124. /)  ! R2 = -c & H like grp (use H)   
                        
! r1H = 1         ! R1 => -H
k_dh2o(r1H, r2CH3, :)  = (/ 4.17E-18 , 1.27 , -2510. /)  ! R2 = -CH3  
k_dh2o(r1H, r2CH2, :)  = (/ 9.97E-18 , 0.91 , -2739. /)  ! R2 = -CH2Ra  
k_dh2o(r1H, r2CH , :)  = (/ 6.79E-18 , 1.16 , -2659. /)  ! R2 = -CHRaRb  
k_dh2o(r1H, r2C  , :)  = (/ 1.81E-17 , 0.99 , -2715. /)  ! R2 = -CRaRbRc  
k_dh2o(r1H, r2H  , :)  = (/ 7.95E-18 , 1.24 , -1510. /)  ! R2 = -H  
k_dh2o(r1H, r2Cd2, :)  = (/ 3.24E-19 , 1.65 , -1271. /)  ! R2 = vinyl (-CdR=CdRR)  
k_dh2o(r1H, r2CO , :)  = (/ 1.56E-18 , 1.53 ,  -100. /)  ! R2 = oxo   
k_dh2o(r1H, r2asH, :)  = (/ 7.95E-18 , 1.24 , -1510. /)  ! R2 = -c & H like grp (use H)  
                        
! r1CH2Cd = 6     ! R1 => -CH2-CR3=CR4R5
k_dh2o(r1CH2Cd, r2CH3, :)  = (/ 1.75E-21 , 2.22 , -1187. /)  ! R2 = -CH3      ! original values corrrected
k_dh2o(r1CH2Cd, r2CH2, :)  = (/ 9.59E-21 , 1.87 , -1452. /)  ! R2 = -CH2Ra    ! original values corrrected
k_dh2o(r1CH2Cd, r2CH , :)  = (/ 2.33E-21 , 2.16 , -1421. /)  ! R2 = -CHRaRb   ! original values corrrected
k_dh2o(r1CH2Cd, r2C  , :)  = (/ 1.18E-19 , 1.59 , -1440. /)  ! R2 = -CRaRbRc  ! original values corrrected
k_dh2o(r1CH2Cd, r2H  , :)  = (/ 5.17E-20 , 1.83 ,  -234. /)  ! R2 = -H   
k_dh2o(r1CH2Cd, r2Cd2, :)  = (/ 1.02E-20 , 2.02 ,   126. /)  ! R2 = vinyl (-CdR=CdRR)  ! values corrrected
k_dh2o(r1CH2Cd, r2CO , :)  = (/ 1.96E-23 , 3.00 ,   642. /)  ! R2 = oxo       ! original values corrrected
k_dh2o(r1CH2Cd, r2asH, :)  = (/ 5.17E-20 , 1.83 ,  -234. /)  ! R2 = -c & H like grp (use H)   
                        
! r1CHCd = 7      ! R1 => -CHRa-CR3=CR4R5
k_dh2o(r1CHCd, r2CH3, :)  = (/ 3.64E-23 , 2.51 , -1511. /)  ! R2 = -CH3      ! original values corrrected
k_dh2o(r1CHCd, r2CH2, :)  = (/ 1.86E-22 , 2.16 , -1698. /)  ! R2 = -CH2Ra    ! original values corrrected
k_dh2o(r1CHCd, r2CH , :)  = (/ 4.52E-23 , 2.45 , -1667. /)  ! R2 = -CHRaRb   ! original values corrrected 
k_dh2o(r1CHCd, r2C  , :)  = (/ 2.29E-21 , 1.88 , -1686. /)  ! R2 = -CRaRbRc  ! original values corrrected  
k_dh2o(r1CHCd, r2H  , :)  = (/ 1.55E-20 , 1.88 ,  -819. /)  ! R2 = -H      
k_dh2o(r1CHCd, r2Cd2, :)  = (/ 3.07E-21 , 2.07 ,  -458. /)  ! R2 = vinyl (-CdR=CdRR)  ! values corrrectedCdR=CdRR)   
k_dh2o(r1CHCd, r2CO , :)  = (/ 5.89E-24 , 3.05 ,    58. /)  ! R2 = oxo       ! original values corrrected
k_dh2o(r1CHCd, r2asH, :)  = (/ 1.55E-20 , 1.88 ,  -819. /)  ! R2 = -c & H like grp (use H)      
                        
! r1CCd  = 8      ! R1 => -CRaRb-CR3=CR4R5
k_dh2o(r1CCd, r2CH3, :)  = (/ 2.64E-20 , 1.77 , -1231. /)  ! R2 = -CH3    
k_dh2o(r1CCd, r2CH2, :)  = (/ 1.34E-19 , 1.42 , -1418. /)  ! R2 = -CH2Ra  
k_dh2o(r1CCd, r2CH , :)  = (/ 3.26E-20 , 1.71 , -1387. /)  ! R2 = -CHRaRb  
k_dh2o(r1CCd, r2C  , :)  = (/ 1.65E-18 , 1.14 , -1405. /)  ! R2 = -CRaRbRc  
k_dh2o(r1CCd, r2H  , :)  = (/ 7.24E-19 , 1.38 ,  -200. /)  ! R2 = -H   
k_dh2o(r1CCd, r2Cd2, :)  = (/ 1.43E-19 , 1.7  ,   161. /)  ! R2 = vinyl (-CdR=CdRR)    
k_dh2o(r1CCd, r2CO , :)  = (/ 2.74E-22 , 2.55 ,   677. /)  ! R2 = oxo    
k_dh2o(r1CCd, r2asH, :)  = (/ 7.24E-19 , 1.38 ,  -200. /)  ! R2 = -c & H like grp (use H)   
                        
! r1Cd2CH3 = 9    ! R1 => -(CR3=CR4CH3)
k_dh2o(r1Cd2CH3, r2CH3, :)  = (/ 2.25E-21 , 2.27 , -493. /)  ! R2 = -CH3   
k_dh2o(r1Cd2CH3, r2CH2, :)  = (/ 2.16E-20 , 1.83 , -738. /)  ! R2 = -CH2Ra   
k_dh2o(r1Cd2CH3, r2CH , :)  = (/ 5.24E-21 , 2.12 , -707. /)  ! R2 = -CHRaRb   
k_dh2o(r1Cd2CH3, r2C  , :)  = (/ 2.65E-19 , 1.55 , -726. /)  ! R2 = -CRaRbRc   
k_dh2o(r1Cd2CH3, r2H  , :)  = (/ 2.42E-19 , 1.64 ,  548. /)  ! R2 = -H    
k_dh2o(r1Cd2CH3, r2Cd2, :)  = (/ 4.78E-20 , 1.82 ,  909. /)  ! R2 = vinyl (-CdR=CdRR)    
k_dh2o(r1Cd2CH3, r2CO , :)  = (/ 9.16E-23 , 2.80 , 1424. /)  ! R2 = oxo   
k_dh2o(r1Cd2CH3, r2asH, :)  = (/ 2.42E-19 , 1.64 ,  548. /)  ! R2 = -c & H like grp (use H)    
                        
! r1Cd2CH2 = 10   ! R1 => -CR3=CR4CH2Ra
k_dh2o(r1Cd2CH2, r2CH3, :)  = (/ 2.25E-21 , 2.27 , -493. /)  ! R2 = -CH3   
k_dh2o(r1Cd2CH2, r2CH2, :)  = (/ 2.16E-20 , 1.83 , -738. /)  ! R2 = -CH2Ra   
k_dh2o(r1Cd2CH2, r2CH , :)  = (/ 5.24E-21 , 2.12 , -707. /)  ! R2 = -CHRaRb   
k_dh2o(r1Cd2CH2, r2C  , :)  = (/ 2.65E-19 , 1.55 , -726. /)  ! R2 = -CRaRbRc   
k_dh2o(r1Cd2CH2, r2H  , :)  = (/ 2.42E-19 , 1.64 ,  548. /)  ! R2 = -H    
k_dh2o(r1Cd2CH2, r2Cd2, :)  = (/ 4.78E-20 , 1.82 ,  909. /)  ! R2 = vinyl (-CdR=CdRR)    
k_dh2o(r1Cd2CH2, r2CO , :)  = (/ 9.16E-23 , 2.80 , 1424. /)  ! R2 = oxo   
k_dh2o(r1Cd2CH2, r2asH, :)  = (/ 2.42E-19 , 1.64 ,  548. /)  ! R2 = -c & H like grp (use H)    
                        
! r1Cd2CH = 11    ! R1 => -CR3=CR4CHRaRb
k_dh2o(r1Cd2CH, r2CH3, :)  = (/ 2.25E-21 , 2.27 , -493. /)  ! R2 = -CH3   
k_dh2o(r1Cd2CH, r2CH2, :)  = (/ 2.16E-20 , 1.83 , -738. /)  ! R2 = -CH2Ra   
k_dh2o(r1Cd2CH, r2CH , :)  = (/ 5.24E-21 , 2.12 , -707. /)  ! R2 = -CHRaRb   
k_dh2o(r1Cd2CH, r2C  , :)  = (/ 2.65E-19 , 1.55 , -726. /)  ! R2 = -CRaRbRc   
k_dh2o(r1Cd2CH, r2H  , :)  = (/ 2.42E-19 , 1.64 ,  548. /)  ! R2 = -H    
k_dh2o(r1Cd2CH, r2Cd2, :)  = (/ 4.78E-20 , 1.82 ,  909. /)  ! R2 = vinyl (-CdR=CdRR)    
k_dh2o(r1Cd2CH, r2CO , :)  = (/ 9.16E-23 , 2.80 , 1424. /)  ! R2 = oxo   
k_dh2o(r1Cd2CH, r2asH, :)  = (/ 2.42E-19 , 1.64 ,  548. /)  ! R2 = -c & H like grp (use H)    

! r1Cd2C = 12     ! R1 => -CR3=CR4R'
k_dh2o(r1Cd2C, r2CH3, :)  = (/ 2.25E-21 , 2.27 , -493. /)  ! R2 = -CH3   
k_dh2o(r1Cd2C, r2CH2, :)  = (/ 2.16E-20 , 1.83 , -738. /)  ! R2 = -CH2Ra   
k_dh2o(r1Cd2C, r2CH , :)  = (/ 5.24E-21 , 2.12 , -707. /)  ! R2 = -CHRaRb   
k_dh2o(r1Cd2C, r2C  , :)  = (/ 2.65E-19 , 1.55 , -726. /)  ! R2 = -CRaRbRc   
k_dh2o(r1Cd2C, r2H  , :)  = (/ 2.42E-19 , 1.64 ,  548. /)  ! R2 = -H    
k_dh2o(r1Cd2C, r2Cd2, :)  = (/ 4.78E-20 , 1.82 ,  909. /)  ! R2 = vinyl (-CdR=CdRR)    
k_dh2o(r1Cd2C, r2CO , :)  = (/ 9.16E-23 , 2.80 , 1424. /)  ! R2 = oxo   
k_dh2o(r1Cd2C, r2asH, :)  = (/ 2.42E-19 , 1.64 ,  548. /)  ! R2 = -c & H like grp (use H)    
                        
! r1CHO = 13      ! R1 => -CHO
k_dh2o(r1CHO, r2CH3, :)  = (/ 8.07E-20 , 1.67 , -2828. /)  ! R2 = -CH3  
k_dh2o(r1CHO, r2CH2, :)  = (/ 5.65E-20 , 1.63 , -3022. /)  ! R2 = -CH2Ra  
k_dh2o(r1CHO, r2CH , :)  = (/ 1.37E-20 , 1.92 , -2991. /)  ! R2 = -CHRaRb  
k_dh2o(r1CHO, r2C  , :)  = (/ 6.94E-19 , 1.35 , -3009. /)  ! R2 = -CRaRbRc  
k_dh2o(r1CHO, r2H  , :)  = (/ 3.04E-19 , 1.59 , -1804. /)  ! R2 = -H         ! error in paper, corrected
k_dh2o(r1CHO, r2Cd2, :)  = (/ 6.02E-20 , 1.78 , -1443. /)  ! R2 = vinyl (-CdR=CdRR)  
k_dh2o(r1CHO, r2CO , :)  = (/ 1.15E-22 , 2.76 ,  -927. /)  ! R2 = oxo   
k_dh2o(r1CHO, r2asH, :)  = (/ 3.04E-19 , 1.59 , -1804. /)  ! R2 = -c & H like grp (use H)
       
! r1CO = 14       ! R1 => -C(O)Ra
k_dh2o(r1CO, r2CH3, :)  = (/ 2.26E-19 , 1.43 , -3279. /)  ! R2 = -CH3  
k_dh2o(r1CO, r2CH2, :)  = (/ 1.15E-18 , 1.08 , -3466. /)  ! R2 = -CH2Ra  
k_dh2o(r1CO, r2CH , :)  = (/ 2.80E-19 , 1.37 , -3435. /)  ! R2 = -CHRaRb  
k_dh2o(r1CO, r2C  , :)  = (/ 1.41E-17 , 0.80 , -3454. /)  ! R2 = -CRaRbRc  
k_dh2o(r1CO, r2H  , :)  = (/ 6.20E-18 , 1.04 , -2248. /)  ! R2 = -H  
k_dh2o(r1CO, r2Cd2, :)  = (/ 1.23E-18 , 1.23 , -1887. /)  ! R2 = vinyl (-CdR=CdRR)  
k_dh2o(r1CO, r2CO , :)  = (/ 2.35E-21 , 2.21 , -1372. /)  ! R2 = oxo  
k_dh2o(r1CO, r2asH, :)  = (/ 6.20E-18 , 1.04 , -2248. /)  ! R2 = -c & H like grp (use H)  

! r1asH = 15         ! R1 => -c & other H like (Use -H as surrogate)
k_dh2o(r1asH, r2CH3, :)  = (/ 4.17E-18 , 1.27 , -2510. /)  ! R2 = -CH3  
k_dh2o(r1asH, r2CH2, :)  = (/ 9.97E-18 , 0.91 , -2739. /)  ! R2 = -CH2Ra  
k_dh2o(r1asH, r2CH , :)  = (/ 6.79E-18 , 1.16 , -2659. /)  ! R2 = -CHRaRb  
k_dh2o(r1asH, r2C  , :)  = (/ 1.81E-17 , 0.99 , -2715. /)  ! R2 = -CRaRbRc  
k_dh2o(r1asH, r2H  , :)  = (/ 7.95E-18 , 1.24 , -1510. /)  ! R2 = -H  
k_dh2o(r1asH, r2Cd2, :)  = (/ 3.24E-19 , 1.65 , -1271. /)  ! R2 = vinyl (-CdR=CdRR)  
k_dh2o(r1asH, r2CO , :)  = (/ 1.56E-18 , 1.53 ,  -100. /)  ! R2 = oxo   
k_dh2o(r1asH, r2asH, :)  = (/ 7.95E-18 , 1.24 , -1510. /)  ! R2 = -c & H like grp (use H)  

init_done = .TRUE.

END SUBROUTINE init_stab_criegee

END MODUlE criegeetool
