MODULE gromhetool
IMPLICIT NONE
CONTAINS

!SUBROUTINE gromhe(chem,Keff)
!SUBROUTINE nitrofol(ng,group,bond,onitrofol)
!SUBROUTINE haloic(ng,group,bond,haloica)
!SUBROUTINE caox(ng,group,bond,caoxa,caoxb)
!SUBROUTINE hydol(ng,group,bond,noh15,noh16)

! ======================================================================
! PURPOSE: This subroutine returns the effective Henry's law coefficient 
!          (Keff) for the species given as input (chem).
! 
! For details about the Gromhe SAR, see T. Raventos-Duran et al., 
! Atmospheric Chemistry and Physics, 7643-7654, 2010.
! ======================================================================
SUBROUTINE gromhe(chem,Keff)
  USE keyparameter, ONLY: mxnode, mxlgr, mxlfo, mxring, mxcp, mxhyd, mxhiso
  USE rjtool, ONLY:rjgrm
  USE searching, ONLY: srh5
  USE stdgrbond, ONLY: grbond
  USE atomtool, ONLY: getatoms
  USE mapping, ONLY: abcde_map, chemmap
  USE khydtool, ONLY: get_hydrate
  USE database, ONLY: nhlcdb,hlcdb_chem,hlcdb_dat 
  USE toolbox, ONLY: stoperr
  IMPLICIT NONE

  CHARACTER(LEN=*),INTENT(IN) :: chem  ! formula of the input species
  REAL,INTENT(OUT) :: keff             ! Effective Henry's law coefficient

! chem tables
  CHARACTER(LEN=LEN(chem)) :: tchem  ! cp of chem
  CHARACTER(LEN=mxlgr) :: group(mxnode)
  INTEGER              :: bond(mxnode,mxnode),node,dbflg,nring
  INTEGER              :: rjg(mxring,2)       ! ring-join group pairs

! group contributors :
  REAL    :: kest           ! estimated intrinsic Henry's law coef.
  REAL    :: sigma,caoxa,caoxb
  INTEGER :: cato,hato,onitrofol,haloica,nogrp,noh16,noh15

! info returned by chemmap
  INTEGER          :: ngrp
  CHARACTER(LEN=1) :: nodetype(mxnode)
  REAL             :: alifun(21),cdfun(21),arofun(21)
  REAL             :: mapfun(mxnode,3,21)
  INTEGER          :: funflg(mxnode)
  INTEGER          :: tabester(4,2)  ! 1= -O- side, 2= CO side
  INTEGER          :: nfcd,nfaro,ierr

! tracks (deep=9)
  INTEGER  :: nabcde(9), tabcde(9,mxcp,mxnode)

! returned by hydration (to compute Keff)
  INTEGER  ::  nwa
  CHARACTER(LEN=mxlfo) chemhyd(mxhyd,mxhiso) ! new formula with hydrates
  REAL     :: yhyd(mxhyd,mxhiso)
  INTEGER  :: nhyd(mxhyd)
  REAL     :: khydstar

  REAL    mult,tsig,dist
  INTEGER i,j,k,l,ipos,it,i1,i2,ncd,inte
  INTEGER nato,oato,ir,is,flato,brato,clato

  CHARACTER(LEN=7)      :: progname='gromhe '
  CHARACTER(LEN=70)     :: mesg

! Tafta sigma values, alifatic compounds:
!  1= -OH    ;  2= -NO2     ;  3= -ONO2   ; 4= -OOH    ;  5= -F     
!  6= -Cl    ;  7=  -Br     ;  8= -I      ; 9= -CHO    ; 10= -CO-   
! 11= -COOH  ; 12= -CO(OOH) ; 13= -PAN    ; 14= -O-    ; 15= R-COO-R
! 16 = HCO-O-R; 17= -CO(F) ; 18= -CO(Cl) ;  19= -CO(Br) ; 20= -CO(I)
! 21 = -CO(O-)
! Tafta sigma values for aliphatic compounds. sigma for ester depends
! whether ester is connected from the -CO- side or the -O- side. 
! Sigma=2.56 for -O- side and 2.00 for CO side.
! CMV 16/06/14 : add sigma taft for CO(O-), n°21
  REAL,DIMENSION(21),PARAMETER :: alisig(21)=  &
    (/ 0.62, 1.47, 1.38, 0.62, 1.10, 0.94, 1.00, 1.00, 2.15, 1.81, &
       2.08, 2.08, 2.00, 1.81, 2.56, 2.90, 2.44, 2.37, 2.37, 2.37, &
      -1.06/)
  REAL, PARAMETER :: sigester=2.00 ! sigma for ester, CO side

  tchem=chem
  IF (chem(1:3)=='#mm' ) tchem(1:)=chem(4:)   

! Check if henry's law constant is already known in the database
!-------------------------------------------------------------
  ipos=srh5(tchem,hlcdb_chem,nhlcdb)
  IF (ipos>0) THEN
    Keff = log10(hlcdb_dat(ipos,1))
    RETURN
  ENDIF
!  DO i=1,nhlcdb
!    IF (tchem == hlcdb_chem(i)) THEN
!      Keff = log10(hlcdb_dat(i,1))
!      RETURN
!    ENDIF
!  ENDDO

! build the group and bond matrix for chem
  CALL grbond(tchem,group,bond,dbflg,nring)
  DO i=1,mxnode  ;  bond(i,i)=0  ;  ENDDO
  CALL rjgrm(nring,group,rjg)  ! rm ring index and get ring closing nodes

  node=0
  DO i=1,mxnode  ;  IF (group(i)(1:1)/=' ') node = node + 1  ;  ENDDO

! Get hydration constant
! ----------------------
  CALL get_hydrate(tchem,nwa,nhyd,chemhyd,yhyd,khydstar)

! ---------------------------------
! Count functionalities
! ---------------------------------

! get the number of atom - C and H counters
  CALL getatoms(tchem,cato,hato,nato,oato,ir,is,flato,brato,clato)
      
  ierr = 0
  CALL chemmap(tchem,node,group,bond,ngrp,nodetype, &
              alifun,cdfun,arofun,mapfun,funflg,   &
              tabester,nfcd,nfaro,ierr)
  IF (ierr/=0) THEN
    mesg="chemmap returned problems"
    CALL stoperr(progname,mesg,chem)
  ENDIF

! ---------------------------------
! compute taft sigmas
! ---------------------------------
  sigma=0.

! 1-1 interaction
! ----------------
  DO i=1,node
    IF (funflg(i)>1) THEN
      DO k=1,8

        DO it=1,2 ! check for both saturated and unsaturated aliphatic
! self interaction
          IF (mapfun(i,it,k)/=0) THEN
            IF (mapfun(i,it,k)>1) THEN
              sigma=sigma+(alisig(k)/0.4)*mapfun(i,it,k)*(mapfun(i,it,k)-1)
            ENDIF
          ENDIF
! cross interaction
          DO l=1,8
            IF ((mapfun(i,it,k)>0.).and.(l/=k)) THEN
               sigma=sigma+(alisig(k)/0.4)*mapfun(i,it,k)*mapfun(i,it,l)
            ENDIF
          ENDDO

        ENDDO
      ENDDO
    ENDIF
  ENDDO   ! end do loop for the 1-1 interaction


! multiple node interaction
! -------------------------
  DO i=1,node    ! start loop among node
    IF (funflg(i)/=0) THEN   ! manage node i
      CALL abcde_map(bond,i,node,nabcde,tabcde)

! ester functionality is on 2 nodes - remove the path that would lead
! to count some interaction twice. This is performed by "killing" the
! "second" node of the functionality.
      DO k=1,4
        IF (tabester(k,1)/=0) THEN

          DO ipos=2,9
            DO it=1,nabcde(ipos)
              DO j=1,ipos-1

                IF (tabcde(ipos,it,j)==tabester(k,1)) THEN
                  IF (tabcde(ipos,it,j+1)==tabester(k,2)) tabcde(ipos,it,j+1)=0  ! set the node to 0
                ENDIF

                IF (tabcde(ipos,it,j)==tabester(k,2)) THEN
                  IF (tabcde(ipos,it,j+1)==tabester(k,1)) tabcde(ipos,it,j+1)=0  ! set the node to 0
                ENDIF

              ENDDO
            ENDDO
          ENDDO

        ENDIF
      ENDDO     ! end of loop to kill irrelevant ester nodes

! start loop for sigma effect
      DO j=1,20
        IF ( (mapfun(i,1,j)/=0).OR. (mapfun(i,2,j)/=0)    )  THEN  !start sigma - taft
          DO inte=2,9
            DO k=1,nabcde(inte)
              ncd=0
              IF (tabcde(inte,k,inte)== 0) CYCLE
              IF (funflg(tabcde(inte,k,inte))== 0)  CYCLE
 
!CMV! rule out case when tabcde(int,k,int) == 0 because
!CMV! associated funflg(0) (invisible seg fault!)
!CMV! can be /= 0, depending on memory issues

              DO l=1,inte-1 ! start the loop to see if Cd are in between
                i1=tabcde(inte,k,l)   ;  IF (i1==0) CYCLE
                i2=tabcde(inte,k,l+1) ;  IF (i2==0) CYCLE
                IF (bond(i1,i2)==2) ncd=ncd+1
              ENDDO
              mult=REAL(funflg(tabcde(inte,k,inte)))
              IF (tabcde(inte,k,2)==0) mult=0. ! kill second node on ester
              tsig=alisig(j)
              IF (j==15) THEN ! check ester
                IF (group(i)(1:2)=='CO') tsig=sigester   ! CO instead of -O-
              ENDIF
              dist=REAL(inte)-REAL(ncd)-2.
              IF (group(tabcde(inte,k,inte))(1:3)=='-O-') dist=dist-1.  ! decrease the distance
              IF (nodetype(i)=='o') dist=dist-1.
              sigma=sigma+(tsig*0.4**(dist)*mult)
            ENDDO
          ENDDO
        ENDIF    ! end sigma - taft

      ENDDO   ! end loop for sigma effect

    ENDIF  ! manage node i
  ENDDO   ! end loop among node

! lump group on Cd and saturated aliphatic
  DO i=1,20
    IF (cdfun(i)>0)  alifun(i)=alifun(i)+cdfun(i)
  ENDDO

! lump group on aromatic and saturated aliphatic
  DO i=1,20
    IF (arofun(i)>0)  alifun(i)=alifun(i)+arofun(i)
  ENDDO

! -------------------------------------
! compute other group-group descriptor
! -------------------------------------

! compute correcting factors
  nogrp=1
  IF (ngrp/=0) nogrp=0

! check for ortho nitro phenols
  CALL nitrofol(node,group,bond,onitrofol)

! check for halogen next to a carboxylic group
  CALL haloic(node,group,bond,haloica)

! check for H bounding from an alcohol
  noh15=0  ;  noh16=0
  IF (alifun(1)>0) CALL hydol(node,group,bond,noh15,noh16)

! check for -CO-C(X) and -CO-C-C(X) structure
  CALL caox(node,group,bond,caoxa,caoxb)


! ALL DATA
  kest = - 1.51807             &
         - 0.13764*sigma       &
         - 2.66308*onitrofol   &
         + 0.973418*haloica    &
         - 1.76731*caoxa       &
         - 1.09538*caoxb       &
         - 1.02636*noh16       &
         - 0.596601*noh15      &
         + 0.499757*cato       &
         - 0.31011*hato        &
         - 0.27637*nogrp       &
         - 0.524482*nfcd       &
         - 1.12707*nfaro       &
         + 4.56593*alifun(1)   &
         + 3.01781*alifun(2)   &
         + 2.03719*alifun(3)   &
         + 4.8687*alifun(4)    &
         + 0.59747*alifun(5)   &
         + 0.870225*alifun(6)  &
         + 1.05894*alifun(7)   &
         + 1.22831*alifun(8)   &
         + 2.59268*alifun(9)   &
         + 3.16535*alifun(10)  &
         + 5.09321*alifun(11)  &
         + 4.683*alifun(12)    &
         + 1.93179*alifun(13)  &
         + 2.4008*alifun(14)   &
         + 2.78457*alifun(15)  &
         + 2.35694*alifun(16)  

! keff=10**kest
! keff=keff*(1+khydstar)
! keff=log10(keff)
! Keff=kest! + log10(1+khydstar)
  Keff=kest + log10(1+khydstar)
  
END SUBROUTINE gromhe


!=======================================================================
! Purpose: Count the number of ortho nitro-phenols group in a molecule
! ======================================================================
SUBROUTINE nitrofol(ng,group,bond,onitrofol)
  IMPLICIT NONE
      
  INTEGER,INTENT(IN) :: ng
  CHARACTER(LEN=*),INTENT(IN) :: group(:)
  INTEGER,INTENT(IN) :: bond(:,:)
  INTEGER,INTENT(OUT):: onitrofol
      
  INTEGER :: tempoh16
  INTEGER :: i,k

  onitrofol=0

! only hydroxy group (whether alkohol or carboxylic) is seek of H bonding
  grloop: DO i=1,ng
    IF (INDEX(group(i),'(OH)')/=0) THEN  ! OPEN 'OH'
      IF (INDEX(group(i),'CO(OH)')/=0) CYCLE grloop ! EXCLUDE ACID H 
      IF (group(i)(1:1)/='c') CYCLE grloop ! EXCLUDE non aromatic OH

! search for nitro group on the alpha carbon
      tempoh16=0     ! a phenol is found
      inloop: DO k=1,ng
        IF (group(k)(1:1)/='c') CYCLE inloop ! EXCLUDE non aromatic groups
        IF (bond(i,k)/=0) THEN
          IF (INDEX(group(k),'(NO2)')/=0) THEN
           tempoh16=tempoh16+1
          ENDIF
        ENDIF
      ENDDO inloop

! increment is maximum 1 per OH (a given OH can only be involded in a single bond)
       IF (tempoh16>0) THEN
          onitrofol=onitrofol+1
       ENDIF
    ENDIF   ! CLOSE 'OH'
  ENDDO grloop

END SUBROUTINE nitrofol

!=======================================================================
! Purpose: Count the # of -C(X)-CO(OH) structure, where X is an halogen 
!=======================================================================
SUBROUTINE haloic(ng,group,bond,haloica)
  IMPLICIT NONE
      
  INTEGER,INTENT(IN) :: ng
  CHARACTER(LEN=*),INTENT(IN) :: group(:)
  INTEGER,INTENT(IN) :: bond(:,:)
  INTEGER,INTENT(OUT):: haloica
 
  INTEGER tempoh15
  INTEGER i,k

  haloica=0

! start loop - search for carboxylic group
  DO i=1,ng
    IF (INDEX(group(i),'CO(OH)')/=0) THEN 
      tempoh15=0
      
! H-bond with a functional group on the alpha carbon
      DO k=1,ng
        IF (bond(i,k)/=0) THEN
          IF (INDEX(group(k),'(F)')/=0) THEN
            tempoh15=tempoh15+1
          ENDIF
          IF (INDEX(group(k),'(Cl)')/=0) THEN
            tempoh15=tempoh15+1
          ENDIF
          IF (INDEX(group(k),'(Br)')/=0) THEN
            tempoh15=tempoh15+1
          ENDIF
          IF (INDEX(group(k),'(I)')/=0) THEN
            tempoh15=tempoh15+1
          ENDIF
        ENDIF
      ENDDO

! increment is maximum 1 per CO(OH) (a given OH can only be involded in a
! single bond)
      IF (tempoh15 > 0) THEN
        haloica=haloica+1
      ENDIF
    ENDIF   
  ENDDO

END SUBROUTINE haloic

!=======================================================================
! Purpose: count the # of -CO-C(X) -and CO-C-C(X) structures, where X is
! an oxygenated group (carbonyl,alkohol, nitro, hydroperoxide, ether) 
!=======================================================================
SUBROUTINE caox(ng,group,bond,caoxa,caoxb)
  IMPLICIT NONE
      
  INTEGER,INTENT(IN) :: ng
  CHARACTER(LEN=*),INTENT(IN) :: group(:)
  INTEGER,INTENT(IN) :: bond(:,:)
  REAL,INTENT(OUT)   :: caoxa  ! caoxa : # of -CO-C(X)
  REAL,INTENT(OUT)   :: caoxb  ! caoxb : # of -CO-C-C(X)
      
  INTEGER tempoh15, tempoh16
  INTEGER i,k,l,m
  REAL    ncoco15,ncoco16

  caoxa=0.  ;  caoxb=0.  ;  ncoco15=0.  ;  ncoco16=0.
      
! start loop - search for carbonyl
  DO i=1,ng
    IF ((group(i)(1:2)=='CO').OR.(group(i)(1:3)=='CHO')) THEN
      tempoh15=0  ;  tempoh16=0

! look for group on alpha position 
      DO k=1,ng
        IF (bond(i,k)/=0) THEN
          IF (INDEX(group(k),'(NO2)')/=0) tempoh15=tempoh15+1
          IF (group(k)(1:2)=='CO') THEN
            tempoh15=tempoh15+1  ;  ncoco15=ncoco15+1
          ENDIF
          IF (INDEX(group(k),'(OH)')/=0) tempoh15=tempoh15+1
          IF (INDEX(group(k),'CHO') /=0) THEN
            tempoh15=tempoh15+1  ;  ncoco15=ncoco15+1
          ENDIF
          IF (INDEX(group(k),'(OOH)')/=0) tempoh15=tempoh15+1

! look for group on beta position 
          betaloop: DO l=1,ng
            IF ((group(i)(1:1)=='c') .AND. (group(l)(1:1)=='c') ) CYCLE betaloop 
            IF ((bond(k,l)/=0).and.(l/=i)) THEN
              IF (INDEX(group(l),'(NO2)')/=0) tempoh16=tempoh16+1
              IF (group(l)(1:2)=='CO') THEN 
                tempoh16=tempoh16+1  ;  ncoco16=ncoco16+1
              ENDIF
              IF (group(l)(1:3)=='-O-') tempoh15=tempoh15+1
              IF (INDEX(group(l),'(OH)') /=0) tempoh16=tempoh16+1
              IF (INDEX(group(l),'(OOH)')/=0) tempoh16=tempoh16+1
              IF (INDEX(group(l),'CHO')  /=0) THEN
                tempoh16=tempoh16+1  ;  ncoco16=ncoco16+1
              ENDIF

              DO m=1,ng     ! check for ether on position gamma
                IF ((bond(l,m)/=0).and.(m/=k)) THEN
                  IF (group(m)(1:3)=='-O-') tempoh16=tempoh16+1
                ENDIF
              ENDDO

            ENDIF
          ENDDO betaloop  

        ENDIF
      ENDDO

      IF (tempoh15 > 0) caoxa=caoxa+1.
      IF (tempoh16 > 0) caoxb=caoxb+1.

    ENDIF
  ENDDO

! correct the value for dicarbonyl (otherwise counted twice)
  IF (ncoco15>1) THEN
    ncoco15=ncoco15/2.  ;  caoxa=caoxa-ncoco15
  ENDIF
  IF (ncoco16>1) THEN
    ncoco16=ncoco16/2.  ;  caoxb=caoxb-ncoco16
  ENDIF

END SUBROUTINE caox

!=======================================================================
! This subroutine count the number of alcohol moiety leading to 
! intramolecular H bounding thrue a 6 or 5 member ring. Priority
! is given to 6 member ring. A alcohol moiety can only be counted once.
!=======================================================================
SUBROUTINE hydol(ng,group,bond,noh15,noh16)
  USE keyparameter, ONLY:mxcp
  USE mapping,  ONLY: abcde_map
  USE ringtool, ONLY: findring
  IMPLICIT NONE
      
! input 
  INTEGER,INTENT(IN) :: ng
  CHARACTER(LEN=*),INTENT(IN) :: group(:)
  INTEGER,INTENT(IN) :: bond(:,:)
  INTEGER,INTENT(OUT):: noh15
  INTEGER,INTENT(OUT):: noh16
      
! local
  INTEGER tempoh15, tempoh16
  INTEGER nohoh15,nohoh16
  INTEGER i,j,k,l,i1,i2
  INTEGER tohoh15,tohoh16

  INTEGER,PARAMETER :: ndeep=4
  INTEGER  nabcde(ndeep), tabcde(ndeep,mxcp,SIZE(group))

  INTEGER    rngflg       ! 0 = 'no', 1 = 'yes'
  INTEGER    ring(SIZE(group))    ! =1 if node participates in current ring
  INTEGER    tring

! initialize
  noh15=0  ;  noh16=0  ;  nohoh15=0  ;  nohoh16=0

! start loop - only hydroxy group (whether alkohol or carboxylic) is seek of H bonding
  grloop: DO i=1,ng
    IF (INDEX(group(i),'(OH)')/=0) THEN              ! OPEN 'OH'
      IF (INDEX(group(i),'CO(OH)')/=0) CYCLE grloop  ! EXCLUDE ACID H 

      tohoh16=0  ;  tohoh15=0  ;  tempoh16=0  ;  tempoh15=0
      CALL abcde_map(bond,i,ng,nabcde,tabcde)

! gamma position - H bonding occurs only if the nodes does not belong to a cycle
      DO k=1,nabcde(4)
        l=tabcde(4,k,4)
        IF (group(l)(1:3)=='-O-') THEN 
           i1=i  ;  i2=tabcde(4,k,2)
           CALL findring(i1,i2,ng,bond,rngflg,ring)
           tring=0
           IF (rngflg/=0) THEN
             DO j=1,4
               IF (ring(tabcde(4,k,j))/=0) tring=tring+1
             ENDDO
           ENDIF
           IF (tring<=2) tempoh16=tempoh16+1
        ENDIF
      ENDDO

! beta position
      DO k=1,nabcde(3)
        l=tabcde(3,k,3)
        i1=i  ;  i2=tabcde(3,k,2)
        CALL findring(i1,i2,ng,bond,rngflg,ring)
        tring=0
        IF (rngflg/=0) THEN
          DO j=1,3
            IF (ring(tabcde(3,k,j))/=0) tring=tring+1
          ENDDO
        ENDIF

        IF (tring<=2) THEN  ! exit if nodes belong to the same ring
          IF (INDEX(group(l),'(ONO2)')/=0) tempoh16=tempoh16+1
          IF (INDEX(group(l),'(F)')/=0)    tempoh16=tempoh16+1
          IF (INDEX(group(l),'(Cl)')/=0)   tempoh16=tempoh16+1
          IF (INDEX(group(l),'(Br)')/=0)   tempoh16=tempoh16+1
          IF (INDEX(group(l),'(I)')/=0)    tempoh16=tempoh16+1
          IF (INDEX(group(l),'(OOH)')/=0)  tempoh16=tempoh16+1
          IF (INDEX(group(l),'CHO')/=0)    tempoh16=tempoh16+1
          IF (INDEX(group(l)(1:3),'CO ')/=0) tempoh16=tempoh16+1
          IF (INDEX(group(l),'(OH)')/=0) THEN
            tempoh16=tempoh16+1
            tohoh16=tohoh16+1       ! 2 hydroxy can only make 1 bond 
          ENDIF
          IF (group(l)(1:3)=='-O-') tempoh15=tempoh15+1
        ENDIF
      ENDDO

! alpha position
      DO k=1,nabcde(2)
        l=tabcde(2,k,2)
        IF (INDEX(group(l),'(NO2)')/=0)  tempoh16=tempoh16+1
        IF (INDEX(group(l),'(ONO2)')/=0) tempoh15=tempoh15+1
        IF (INDEX(group(l),'(F)')/=0)    tempoh15=tempoh15+1
        IF (INDEX(group(l),'(Cl)')/=0)   tempoh15=tempoh15+1
        IF (INDEX(group(l),'(Br)')/=0)   tempoh15=tempoh15+1
        IF (INDEX(group(l),'(I)')/=0)    tempoh15=tempoh15+1
        IF (INDEX(group(l),'(OOH)')/=0)  tempoh15=tempoh15+1
        IF (INDEX(group(l),'CHO')/=0)    tempoh15=tempoh15+1
        IF (INDEX(group(l)(1:3),'CO ')/=0) tempoh15=tempoh15+1
        IF (INDEX(group(l),'(OH)')/=0) THEN
          tempoh15=tempoh15+1  ;  tohoh15=tohoh15+1
        ENDIF
      ENDDO

! curent position
      IF (INDEX(group(i),'(OOH)')/=0) tempoh15 = tempoh15+1

! Analyse H bonds. A given -OH make only one bond. Priority is given to 
! 6 member ring, then 5 member ring. Care must be taken for dihydroxy 
! species, which can only make one H-bond
      IF (tempoh16>0) THEN
        noh16=noh16+1
        IF ((tohoh16>0).AND.(tempoh16==1)) nohoh16=nohoh16+1   ! Hbond through dihydroxy group only
      ELSE IF (tempoh15>0) THEN
        noh15=noh15+1
        IF ((tohoh15>0).AND.(tempoh15==1)) nohoh15=nohoh15+1   ! Hbond through dihydroxy group only
      ENDIF
    ENDIF 
  ENDDO grloop

  nohoh16=nohoh16/2       ! integer division on purpose
  IF (nohoh16/=0) noh16=noh16-nohoh16
  nohoh15=nohoh15/2       ! integer division on purpose
  IF (nohoh15/=0) noh15=noh15-nohoh15

END SUBROUTINE hydol


END MODULE gromhetool
