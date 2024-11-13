MODULE poztool
USE toolbox, ONLY: stoperr,add1tonp,addref
IMPLICIT NONE

CONTAINS
!=======================================================================
! PURPOSE: make the products from POZ decomposition according to the 
! nodes given as input for the C atom bearing the Criegee group (icri)
! and the C atom bearing the carbonyl (or acyl) group (icarb). The 
! subroutine also return the ratio for the Z conformer of the criegee 
! species (but may not be relevent). Indeed, the Z/E conformers may not  
! may exist (e.g. for identical branches). This will be tested in 
! the subroutines devoted to CI* chemistry.
!=======================================================================
SUBROUTINE pozprod(chem,bond,group,zebond,icarb,icri, &
                   chemcarb,chemacyl,chemcri,zconform)
  USE reactool, ONLY: swap,rebond
  USE ringtool, ONLY: findring,rejoin
  USE normchem, ONLY: stdchm
  USE toolbox, ONLY: addrx,stoperr,addref
  USE fragmenttool, ONLY: fragm
  USE stdgrbond, ONLY: grbond
  IMPLICIT NONE

  CHARACTER(LEN=*),INTENT(IN)  :: chem          ! chemical formula
  CHARACTER(LEN=*),INTENT(IN)  :: group(:)      ! group matrix
  INTEGER,INTENT(IN)           :: bond(:,:)     ! node matrix
  INTEGER,INTENT(IN)           :: zebond(:,:)   ! cis/trans info on C=C bond
  INTEGER,INTENT(IN)           :: icarb         ! node # bearing >C=O 
  INTEGER,INTENT(IN)           :: icri          ! node # bearing >C.OO.* 
  CHARACTER(LEN=*),INTENT(OUT) :: chemcarb      ! formula of the carbonyl product
  CHARACTER(LEN=*),INTENT(OUT) :: chemacyl      ! formula of the acyl radical product 
  CHARACTER(LEN=*),INTENT(OUT) :: chemcri       ! formula of the hot criegee CI* product 
  REAL,INTENT(OUT)             :: zconform      ! ratio of the CI* Z conformer

!  INTEGER :: alkene_isom=2                   ! Parent alkene is a cis isomer
  CHARACTER(LEN=LEN(group(1))) :: tgroup(SIZE(group))
  CHARACTER(LEN=LEN(group(1))) :: pold, pnew,tempgr
  INTEGER :: tbond(SIZE(bond,1),SIZE(bond,2))
  INTEGER :: ngr,rngflg,nc,i,ig
  INTEGER :: mm,nn,tempnring,tempdbflg
  INTEGER :: ring(SIZE(group))
  CHARACTER(LEN=LEN(chem)) :: chema,chemb       
  LOGICAL :: loze

  CHARACTER(LEN=10),PARAMETER :: progname='pozprod'
  CHARACTER(LEN=70)           :: mesg

  tgroup(:)=group(:) ; tbond(:,:)=bond(:,:)
  ngr=COUNT(tgroup/=' ')
  chemcarb=' ' ; chemcri=' ' ; chemacyl=' '

! check if bond that will be brocken belong to a cycle.
  CALL findring(icarb,icri,ngr,tbond,rngflg,ring)

! break the C=C bond and change "Cd" to "C"
  tbond(icarb,icri)=0 ; tbond(icri,icarb)=0
  pold='Cd' ; pnew='C' 
  tempgr=group(icarb) ; CALL swap(tempgr,pold,tgroup(icarb),pnew)
  tempgr=group(icri)  ; CALL swap(tempgr,pold,tgroup(icri),pnew)

! make carbonyl on icarb node
! ---------------------------
  IF     (tgroup(icarb)(1:3)=='CH2')THEN; pold='CH2'; pnew='CH2O'
  ELSE IF(tgroup(icarb)(1:2)=='CH') THEN; pold='CH' ; pnew='CHO'
  ELSE IF(tgroup(icarb)(1:2)=='CO') THEN; pold='CO' ; pnew='CO2' ! ketene structure
  ELSE IF(tgroup(icarb)(1:1)=='C')  THEN; pold='C'  ; pnew='CO'
  ELSE
    mesg="group not identified while trying to make carbonyl formula"
    CALL stoperr(progname,mesg,chem)
  ENDIF
  tempgr=tgroup(icarb) ; CALL swap(tempgr,pold,tgroup(icarb),pnew)

! make hot criegee on icri node
! -----------------------------
  nc=INDEX(tgroup(icri),' ') ; tgroup(icri)(nc:nc+6)='.(OO.)*'

! check the # of distinct molecules. If no ring => fragments into 2 products.
  IF (rngflg==1) THEN
    CALL rejoin(ngr,icarb,icri,mm,nn,tbond,tgroup)
    CALL rebond(tbond,tgroup,chemcri,tempnring)
    CALL stdchm(chemcri)

  ELSE IF (rngflg==0) THEN
    CALL fragm(tbond,tgroup,chema,chemb)
    CALL stdchm(chema) ; CALL stdchm(chemb)
    IF      (INDEX(chema,'(OO.)*')/=0) THEN ; chemcri=chema ; chemcarb=chemb
    ELSE IF (INDEX(chemb,'(OO.)*')/=0) THEN ; chemcri=chemb ; chemcarb=chema
    ELSE
      mesg="expected criegee not found in decomposition products"
      CALL stoperr(progname,mesg,chem)
    ENDIF

  ELSE
    mesg="expected number for rngflg (0 or 1) "
    CALL stoperr(progname,mesg,chem)
  ENDIF

  tgroup(:)=group(:) ; tbond(:,:)=bond(:,:)

! Compute Z/E conformer yield for the hot criegee (which might not exist !)
  loze=.TRUE.  ; zconform=0.5    ! default value is 50/50 for the Z and E ratio
  IF      (group(icri)(1:4)=='CdH2') THEN ; loze=.FALSE. ! no Z/E for CH2.(OO.) 
  ELSE IF (group(icri)(1:3)/='CdH')  THEN ; loze=.FALSE. ! Internal criegee (Z/E = 50/50)
  ENDIF 

! Check for C=C-CH.(OO.)
  IF (loze) THEN
    DO i=1,ngr
      IF (bond(icri,i)/=0) THEN
        IF (i==icarb) CYCLE
        IF (group(i)(1:2)=='Cd') THEN
          zconform=0.2 ; loze=.FALSE. 
        ENDIF
      ENDIF
    ENDDO
  ENDIF  

! Check for cis/trans RCH=CHR - change default value if cis alkene
  IF (loze) THEN
    IF (group(icarb)(1:3)=='CdH') THEN
      IF (zebond(icri,icarb)==1) THEN    ! cis RCH=CHR case here
        zconform=0.3 ; loze=.FALSE. 
      ENDIF
    ENDIF
  ENDIF
  
! make acyl on the criegee side (Carbonyl-hydroperoxide route) 
! -----------------------------
  CALL grbond(chemcri,tgroup,tbond,tempdbflg,tempnring) ! use chemcri here
  ngr=COUNT(tgroup/=' ')  ;  ig=0 
  DO i=1,ngr 
    IF (INDEX(tgroup(i),'.(OO.)*')/=0) ig=i  
  ENDDO
  IF (ig==0) THEN
    mesg="expected hot_criegee functional group not found"
    CALL stoperr(progname,mesg,chem)
  ENDIF

! make CO. on node bearing the criegee group
  IF (tgroup(ig)(1:11)=='CH2.(OO.)* ') THEN   ! ignore other C1 Criegee
    chemacyl='CHO.'
  ELSE IF (tgroup(ig)(1:10)=='CH.(OO.)* ') THEN
    tgroup(ig)='CO. '
    CALL rebond(tbond,tgroup,chemacyl,tempnring)
    CALL stdchm(chemacyl)
  ENDIF
  
END SUBROUTINE pozprod

!=======================================================================
! PURPOSE: Get CI* yield (yci) for each Cd side of the POZ structure. 
! Note that the carbonyl-hydroperoxide route is not considered here and
! yci(1)+yci(2)=1. CI* yields are based on Newland et al., 2022.  
!=======================================================================
SUBROUTINE poz_decomp(chem,group,bond,poznod,cdsub,ring,yci)
  USE toolbox, ONLY: stoperr,addref
  IMPLICIT NONE

  CHARACTER(LEN=*),INTENT(IN)  :: chem         ! chemical formula
  CHARACTER(LEN=*),INTENT(IN)  :: group(:)     ! group matrix
  INTEGER,INTENT(IN)           :: bond(:,:)    ! node matrix
  INTEGER,INTENT(IN)           :: poznod(:)    ! 2 Cd nodes # bearing the poz
  INTEGER,INTENT(IN)           :: cdsub(:)     ! # of C (including CO) bonded to Cd in cdtable(i)
  INTEGER,INTENT(IN)           :: ring(:)      ! nodes part of rings
  REAL,INTENT(OUT)             :: yci(:)       ! CI* yield on the 1st and 2nd Cd nodes
  
  INTEGER  :: ngr,i,j,k
  INTEGER  :: poz1,poz2,naa,nbb
  REAL     :: fval(SIZE(yci))

  CHARACTER(LEN=15),PARAMETER  :: progname='poz_decomp'
  CHARACTER(LEN=70)            :: mesg

  fval(:)=0.
  ngr=COUNT(group/=' ')

! compute the factors both sides of the C=C (POZ) bond
  DO i=1,2
    poz1=poznod(i)
    IF (i==1) THEN ; poz2=poznod(2) ; ELSE ;  poz2=poznod(1) ; ENDIF
    
    ! check for exocyclic C=C
    IF (ring(poz1)/=0) THEN
      IF (ring(poz2)==0) fval(i)=fval(i)+0.62 ! exocyclic C=C 
    ENDIF
    
    IF (cdsub(poz1)==0) CYCLE ! no nodes on Cd (BA, Warning OH from enols)

    subloop: DO j=1,ngr       ! scroll for neighbors
      IF (bond(poz1,j)==0) CYCLE subloop
      DO k=1,2
        IF (j==poznod(k)) CYCLE subloop ! exclude C=C nodes
      ENDDO

      IF (group(j)(1:3)=='CHO') THEN ; fval(i)=fval(i)+0.127 ; CYCLE subloop ; ENDIF  ! aldehyde groups
      IF (group(j)(1:2)=='-O')  THEN ; fval(i)=fval(i)-0.655 ; CYCLE subloop ; ENDIF  ! -O-R groups
      IF (group(j)(1:2)=='Cd')  THEN ; fval(i)=fval(i)-0.28  ; CYCLE subloop ; ENDIF  ! conjugated -C=C 
      IF (ring(poz1)==0) THEN
        IF (group(j)(1:1)=='c') THEN ; fval(i)=fval(i)-0.25  ; CYCLE subloop ; ENDIF  ! -c (aromatic) groups
        IF (ring(j)==1)         THEN ; fval(i)=fval(i)-0.25  ; CYCLE subloop ; ENDIF  ! cyclic groups
      ENDIF
      IF (group(j)(1:2)=='CO')  THEN                          ! carbonyl groups
        IF (group(j)(1:3)=='CO(') THEN                        ! CO(OH), CO(OOH), CO(OONO2)
          fval(i)=fval(i) + 0.00 ; CYCLE subloop                                     
        ELSE                                                  ! other 
          DO k=1,ngr
            IF (bond(k,j)/=0) THEN
              IF (group(k)(1:2)=='-O') THEN
                fval(i)=fval(i)+0.00 ; CYCLE subloop          ! CO-OR
              ENDIF
            ENDIF
          ENDDO
          fval(i)=fval(i)+0.127 ; CYCLE subloop                ! -CO-R ("simple" ketone)
        ENDIF
      ENDIF

      naa=COUNT(bond(j,:)/=0)  ! count the number of neighbors
      IF (naa==1)      THEN ; fval(i)=fval(i)+0.218 ! CH3 or CH2(X)
      ELSE IF (naa==3) THEN ; fval(i)=fval(i)-0.069 ! CHR1R2
      ELSE IF (naa==4) THEN ; fval(i)=fval(i)-0.386 ! CR1R2R3
      ELSE IF (naa==2) THEN                         ! CH2R
        DO k=1,ngr
          IF (bond(k,j)==0) CYCLE
          IF (k==poz1) CYCLE
          nbb = COUNT(bond(k,:)/=0)
          IF (nbb==1) THEN ; fval(i)=fval(i)+0.107
          ELSE             ; fval(i)=fval(i)+0.00
          ENDIF
          EXIT
        ENDDO
      ELSE
        mesg="unexpected value for cdsub"
        CALL stoperr(progname,mesg,chem)        
      ENDIF
    ENDDO subloop
  ENDDO

! compute the fragmentation yield for the production of the hot criegee
  yci(1)=(fval(1)-fval(2)+1.) / 2. 
  yci(2)= 1. - yci(1)
 
! overwrite for out of range values (negative yields)
  IF (yci(1) < 0.) THEN   
    IF (yci(2) < 0.) THEN
      mesg="unexpected values for yci (criegee yield) < 0 both side"
      CALL stoperr(progname,mesg,chem)        
    ENDIF
    yci(2)=1. ; yci(1)=0. 
  ELSE IF (yci(2) < 0.) THEN
    yci(1)=1. ; yci(2)=0.
  ENDIF   

! check consistencies
  DO i=1,2
    IF ((yci(i) > 1.) .OR. (yci(i) < 0.)) THEN
      mesg="unexpected value for yci (criegee yield)"
      CALL stoperr(progname,mesg,chem)        
    ENDIF
  ENDDO 
  IF (ABS(1.-(yci(1)+yci(2))) > 0.001) THEN
    mesg="Decomposition POZ pathways is not 1."
    CALL stoperr(progname,mesg,chem)        
  ENDIF
END SUBROUTINE poz_decomp

!=======================================================================
! PURPOSE: For C=C-C=C structure, the subroutine return the ratio of 
! O3 addition on each C=C bond. The ratios are based on the SAR by 
! Jenkin et al., 2020
!=======================================================================
SUBROUTINE get_pozratio(cdtable,cdsub,pozratio)
  IMPLICIT NONE

  INTEGER,INTENT(IN) :: cdtable(:)   ! node # bearing a "Cd"
  INTEGER,INTENT(IN) :: cdsub(:)     ! # of C (including CO) bonded to Cd in cdtable(i)
  REAL,INTENT(OUT)   :: pozratio(:)  ! branching ratio for POZ on 1st and 2nd C=C bonds 

  INTEGER :: nsub12,nsub34

  pozratio(:)=0.
  nsub12 = cdsub(cdtable(1)) + cdsub(cdtable(2))
  nsub34 = cdsub(cdtable(3)) + cdsub(cdtable(4))
  IF      (nsub12==nsub34)     THEN ; pozratio(1)=0.5 ; pozratio(2)=0.5
  ELSE IF (nsub12==(nsub34-1)) THEN ; pozratio(1)=0.6 ; pozratio(2)=0.4
  ELSE IF (nsub12==(nsub34-2)) THEN ; pozratio(1)=0.7 ; pozratio(2)=0.3
  ELSE IF (nsub12==(nsub34-3)) THEN ; pozratio(1)=0.8 ; pozratio(2)=0.2
  ELSE IF (nsub34==(nsub12-1)) THEN ; pozratio(1)=0.4 ; pozratio(2)=0.6
  ELSE IF (nsub34==(nsub12-2)) THEN ; pozratio(1)=0.3 ; pozratio(2)=0.7
  ELSE IF (nsub34==(nsub12-3)) THEN ; pozratio(1)=0.2 ; pozratio(2)=0.8
  ENDIF  

END SUBROUTINE get_pozratio  

END MODULE poztool
