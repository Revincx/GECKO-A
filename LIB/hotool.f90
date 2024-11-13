MODULE hotool
IMPLICIT NONE
CONTAINS

!=======================================================================
! PURPOSE :                                                           
!   Find rate constants for H-atom abstraction from VOC by HO based on
!   Jenkin et al., ACP, 2018
!=======================================================================
SUBROUTINE rabsoh(tbond,tgroup,ig,arrhc,nring)
  USE keyparameter, ONLY: mxcp,mxlfo,saru
  USE keyflag, ONLY: wrtsarinfo
  USE ho_aromtool, ONLY: arom_data
  USE ringtool, ONLY: ring_data
  IMPLICIT NONE

  INTEGER,INTENT(IN)         :: tbond(:,:) ! node matrix for the current species
  CHARACTER(LEN=*),INTENT(IN):: tgroup(:)  ! groups for the current species
  INTEGER,INTENT(IN)         :: ig         ! group bearing the leaving H 
  INTEGER,INTENT(IN)         :: nring      ! # of rings
  REAL,INTENT(OUT)           :: arrhc(:)   ! arrh. coef. for H abstraction

  REAL     :: aalk                ! alkyl contribution to A factor 
  REAL     :: balk                ! alkyl contribution to B factor 
  LOGICAL  :: loalk               ! alkyl contribution required (flag) 
  LOGICAL  :: loester             ! ester contribution required (flag) 
  INTEGER  :: nfun                ! # of contribution to substituent factor
  INTEGER,PARAMETER :: mxfact=10
  REAL     :: af(mxfact)          ! contribution to A substituent factor
  REAL     :: bf(mxfact)          ! contribution to B substituent factor
  REAL     :: ax                  ! overall substituent A factor
  REAL     :: bx                  ! overall substituent B factor

  INTEGER  :: i,j,k
  INTEGER  :: nca
  INTEGER  :: nbF                          ! # of "correction factors" applied
  INTEGER  :: o_sub(2),m_sub(2),p_sub(2)   ! branches in o,m,p position - if any
  INTEGER  :: nCd                          ! # of Cd bounded to ig
  INTEGER  :: alf(3)                       ! index of Cd neighbors
  INTEGER  :: Cdtable(4)                   ! Cd nodes (bounded to ig)

! initialize
! -----------
  arrhc(:)=0. 
  nca=COUNT(tgroup/=' ')
  IF (ig>nca) STOP "ig > nca in rabsoh" 
      
  IF (wrtsarinfo) THEN
    WRITE(saru,*) '   '                                  !! debug
    WRITE(saru,*) '============= RABSOH ======'          !! debug
    DO i=1,nca ;  WRITE(saru,*) i,tgroup(i) ; ENDDO      !! debug
    WRITE(saru,*) 'ig,tgroup(ig)',ig,tgroup(ig)          !! debug
  ENDIF
            
! --------------------
! --- FIND K(0) VALUE (Tables 1, 4 and 7 in Jenkin et al., 2018)
! ------------------ -
  IF (tgroup(ig)(1:3)=='CH3') THEN
    arrhc(1)=2.90E-12 ; arrhc(3)=925.

  ELSE IF(tgroup(ig)(1:3)=='CH2') THEN
    arrhc(1)=4.95E-12 ; arrhc(3)=555.

  ELSE IF(tgroup(ig)(1:3)=='CHO') THEN
    DO i=1,nca
      IF (tbond(ig,i)==3) THEN                     ! -O-CHO group
        arrhc(1)=1.70E-12 ; arrhc(3)=910.
        IF (wrtsarinfo) CALL wrtkabs(saru,tgroup(ig),arrhc) ; RETURN !--- Out
      ELSE IF (tbond(ig,i)==1) THEN
        DO j=1,nca
          IF ((tbond(j,i)==1).AND.(j/=ig)) THEN    ! C(OH)-C(polar)-CHO
            IF ((INDEX(tgroup(i),'O')==0).AND.(INDEX(tgroup(j),'(OH)')/=0)) THEN   
              arrhc(1)=11.7E-12 ; arrhc(3)=78.
              IF (wrtsarinfo) CALL wrtkabs(saru,tgroup(ig),arrhc) ; RETURN !--- Out
            ENDIF
          ENDIF
        ENDDO
        IF (tgroup(i)(1:3)=='CH3')       THEN ; arrhc(1)=4.60E-12 ; arrhc(3)=-350. ! CH3CHO
        ELSE IF (tgroup(i)(1:2)=='Cd')   THEN ; arrhc(1)=1.30E-11 ; arrhc(3)=0.    ! C=CCHO
        ELSE IF (tgroup(i)(1:4)=='CH2 ') THEN ; arrhc(1)=5.08E-12 ; arrhc(3)=-420. ! -CH2CHO
        ELSE IF (tgroup(i)(1:4)=='CHO ') THEN ; arrhc(1)=1.55E-12 ; arrhc(3)=-340. ! CHOCHO
        ELSE IF (tgroup(i)(1:5)=='CH2(O')THEN ; arrhc(1)=11.7E-12 ; arrhc(3)=140.  ! CH2(polar)CHO
        ELSE IF (tgroup(i)(1:2)=='CO')   THEN ; arrhc(1)=1.78E-12 ; arrhc(3)=-590. ! -COCHO           
        ELSE IF ((tgroup(i)(1:4)=='CH(O').OR.(tgroup(i)(1:3)=='C(O')) THEN         ! >C(polar)CHO
          arrhc(1)=11.7E-12 ; arrhc(3)=-25.            
        ELSE IF (tgroup(i)(1:1)=='c') THEN                                         ! phi-CHO
          arrhc(1)=1.21E-11 ; arrhc(3)=0.  
          CALL arom_data(i,tgroup,tbond,nca,o_sub,m_sub,p_sub)
          DO k=1,2
            IF (o_sub(k)/=0) THEN 
              IF (INDEX(tgroup(o_sub(k)),'O')==0) arrhc(3)=arrhc(3)-115.
            ENDIF
            IF (m_sub(k)/=0) THEN 
              IF (INDEX(tgroup(m_sub(k)),'O')==0) arrhc(3)=arrhc(3)-78.
            ENDIF
            IF (p_sub(k)/=0) THEN 
              IF (INDEX(tgroup(p_sub(k)),'O')==0) arrhc(3)=arrhc(3)-78.
            ENDIF
          ENDDO
        ELSE IF (tgroup(i)(1:1)=='C')   THEN ; arrhc(1)=5.22E-12 ;  arrhc(3)=-490. ! CCHO           
        ELSE
          WRITE(6,*) 'No constant found in rabsoh for group: ', tgroup(i)
          STOP "in rabsoh, expected reaction not found"
        ENDIF
        IF (wrtsarinfo) CALL wrtkabs(saru,tgroup(ig),arrhc) ; RETURN !--- Out
      ENDIF
    ENDDO
    
  ELSE IF(tgroup(ig)(1:2)=='CH') THEN
    arrhc(1)=3.17E-12 ; arrhc(3)=225.

  ELSE
    WRITE(6,*) '--error--, in rabsoh. No reaction found for group: ', tgroup(ig)
    STOP "in rabsoh"
  ENDIF
      
! K0 for alpha ether
! -----------------
  eloop: DO i=1,nca
    IF (tbond(ig,i)==3) THEN
      DO j=1,nca                                  ! skip -CO-O-CH-
        IF ( (tbond(i,j)==3).AND.(j/=ig).AND.&   
             ((tgroup(j)(1:3)=='CO ').OR.(tgroup(j)(1:3)=='CHO')) ) THEN        
          CYCLE eloop
        ENDIF
      ENDDO

      IF (tgroup(ig)(1:3)=='CH3') THEN
        arrhc(1)=2.22E-12 ; arrhc(3)=160.         ! -O-CH3
        IF (wrtsarinfo) CALL wrtkabs(saru,tgroup(ig),arrhc) ; RETURN !--- Out
      ELSE IF (tgroup(ig)(1:2)=='CH') THEN        ! -O-CH<
        DO j=1,nca
          IF (tbond(ig,j)==1) THEN
            DO k=1,nca
              IF ((tbond(j,k)==3).AND.(nring==0)) THEN  ! -O-CH-C-O- 
                arrhc(1)=1.17E-12 ;  arrhc(3)=-760.    
                IF (wrtsarinfo) CALL wrtkabs(saru,tgroup(ig),arrhc) ; RETURN !--- Out
              ENDIF
            ENDDO
          ENDIF
        ENDDO
        arrhc(1)=1.20E-12 ; arrhc(3)=-460.        ! -O-CH< "regular" 
        IF (nring/=0) THEN
          CALL ringfac(ig,nca,tbond,tgroup,ax,bx)  
          arrhc(1)=arrhc(1)*ax                     
          arrhc(3)=arrhc(3)+bx                     
        ENDIF
        IF (wrtsarinfo) CALL wrtkabs(saru,tgroup(ig),arrhc) ; RETURN !--- Out
      ENDIF
    ENDIF
  ENDDO eloop
  
! K0 for Aromatic branches (Fph1(1)=8.6, Fph1(3)=345, Fph2(1)=7.0, Fph2(3)=580)
  aloop: DO i=1,nca
    IF (tbond(ig,i)/=1) CYCLE
    IF (tgroup(i)(1:1)=='c') THEN
      IF ((tgroup(ig)(1:4)=='CH3 ').OR.(tgroup(ig)(1:4)=='CH2(')) THEN
        arrhc(1)=arrhc(1)*8.6 ; arrhc(3)=arrhc(3)+345.
      ELSE IF (tgroup(ig)(1:2)=='CH') THEN
        arrhc(1)=arrhc(1)*7.0 ; arrhc(3)=arrhc(3)+580.
      ENDIF
! scaling factor due to alkyl group in orto and para positions
      CALL arom_data(i,tgroup,tbond,nca,o_sub,m_sub,p_sub)
      DO k=1,2
        IF (o_sub(k)/=0) THEN 
          IF (INDEX(tgroup(o_sub(k)),'O')==0) arrhc(3)=arrhc(3)-140.
        ENDIF
        IF (p_sub(k)/=0) THEN 
          IF (INDEX(tgroup(p_sub(k)),'O')==0) arrhc(3)=arrhc(3)-140.
        ENDIF
      ENDDO
      EXIT aloop
    ENDIF
  ENDDO aloop
      
  IF (wrtsarinfo) WRITE(saru,*) ' 1) rate constant =',arrhc(1:3)

      
! -----------------------
! --- F(X) CORRECTION 
! -----------------------

! FIND SUBSTITUENTS ON SAME CARBON (Tables 5 and 9 in Jenkin et al., 2018)
  nfun=0 ; af(:)=0. ; bf(:)=0.                                                                 
  IF (INDEX(tgroup(ig),'(NO2)')/=0)  THEN ; arrhc(1)=0. ; arrhc(3)=0.    ; RETURN         ; ENDIF 
  IF (INDEX(tgroup(ig),'(OH)')/=0)   THEN ; nfun=nfun+1 ; af(nfun)=0.497 ; bf(nfun)=-590. ; ENDIF
  IF (INDEX(tgroup(ig),'(ONO2)')/=0) THEN ; nfun=nfun+1 ; af(nfun)=0.127 ; bf(nfun)=-70.  ; ENDIF
  IF (INDEX(tgroup(ig),'(OOH)') /=0) THEN ; nfun=nfun+1 ; af(nfun)=0.497 ; bf(nfun)=-590. ; ENDIF
  IF (INDEX(tgroup(ig),'(OOOH)')/=0) THEN ; nfun=nfun+1 ; af(nfun)=0.497 ; bf(nfun)=-590. ; ENDIF
  DO j=1,nfun
    arrhc(1)=arrhc(1)*af(j) ; arrhc(3)=arrhc(3)+bf(j)
    IF (wrtsarinfo) THEN
      WRITE(saru,*) ' 2) fact on carbon bearing leaving H:',af(j),'*exp(-',bf(j),'/T' !! debug
    ENDIF
  ENDDO

  IF (wrtsarinfo) THEN
    WRITE(saru,*) ' 3) Substituents factors:'                  !! debug
  ENDIF
            
! C=C neighbors - Seek the 2 neighbors (if any) in case of superallyl 
! resonant structures
  alf(:)=0 ; Cdtable(:)=0 ; nCd=0
  DO i=1,nca
    IF ((tgroup(i)(1:2)=='Cd').AND.(tbond(ig,i)==1)) THEN
      nCd=nCd+1 ; alf(nCd)=i
    ENDIF
  ENDDO

  IF (nCd/=0) THEN           ! CH-Cd=Cd
    Cdtable(1)=alf(1)        ! look for the Cd on the chain
    DO i=1,nca
      IF (tbond(Cdtable(1),i)==2) THEN ; Cdtable(2)=i ; EXIT ; ENDIF
    ENDDO
    DO i=1,nca
      IF ((tbond(Cdtable(2),i)==1).AND.(tgroup(i)(1:2)=='Cd')) THEN
        Cdtable(3)=i ; EXIT
      ENDIF
    ENDDO
    IF (Cdtable(3)/=0) THEN  ! CH-Cd=Cd-Cd=Cd-
      DO i=1,nca
        IF (tbond(Cdtable(3),i)==2) THEN ; Cdtable(4)=i ; EXIT ; ENDIF
      ENDDO
    ENDIF
  ENDIF

  IF (nCd==2) THEN           ! C=C at both side of the C-H node (overwrite if used above)   
    Cdtable(3)=alf(2)        ! look for the Cd on the chain
    DO i=1,nca
      IF (tbond(Cdtable(3),i)==2) THEN ; Cdtable(4)=i ; EXIT; ENDIF
    ENDDO
  ENDIF

  IF (nCd>2) STOP "unexpected 3 Cd attached to the C-H node"      

! FIND SUBSTITUENTS ON ALPHA CARBONS such as F(X)=A(FX) exp(-B(FX)/T)
  alphaloop: DO i=1,nca
    aalk=0  ; balk=0   ; loalk=.FALSE. ; loester=.FALSE.
    nfun=0  ; af(:)=0. ; bf(:)=0.
    IF (tbond(ig,i)/=0) THEN

      IF (wrtsarinfo) THEN
        WRITE(saru,*) ' 3.1) alpha group factors'           !! debug
        WRITE(saru,*) '  >>> group_neighbor=',i,tgroup(i)   !! debug
      ENDIF

! simple alkyl (Table 2 in Jenkin et al., 2018)
      IF (tgroup(i)(1:3)=='CH3') THEN ; aalk=1. ; balk= 0.  ; loalk=.TRUE. ; ENDIF
      IF (tgroup(i)(1:4)=='CH2 ')THEN ; aalk=1. ; balk=-89. ; loalk=.TRUE. ; ENDIF
      IF (tgroup(i)(1:3)=='CH ') THEN ; aalk=1. ; balk=-89. ; loalk=.TRUE. ; ENDIF
      IF (tgroup(i)(1:2)=='C ')  THEN ; aalk=1. ; balk=-89. ; loalk=.TRUE. ; ENDIF

! contribution of functional groups

! alcohol (Table 5 in Jenkin et al., 2018)
      IF ((tgroup(i)(1:2)/='CO').AND.(INDEX(tgroup(i),'(OH)')/=0)) THEN
        nfun=nfun+1 ; af(nfun)=0.119 ; bf(nfun)=-930.
      ENDIF

! carbonyl (Tables 5, 8 and 9 in Jenkin et al., 2018)
      IF (tgroup(i)(1:3)=='CO ') THEN
        DO j=1,nca
          IF (tbond(i,j)==3) THEN                  ! ester -CH-CO-O-R
            nfun=nfun+1 ; af(nfun)=0.783 ; bf(nfun)=200.
            loester=.TRUE.
            EXIT
          ENDIF
        ENDDO
        DO j=1,nca
          IF ((tbond(i,j)/=0).AND.(j/=ig)) THEN     ! CO(OH)-CO-CH-
            IF (tgroup(j)=='CO(OH)') THEN          
              nfun=1 ; af(nfun)=0 ; bf(nfun)=0.     ! reaction only occurs on -COOH side
              GOTO 10
            ENDIF
          ENDIF
        ENDDO
        IF (.NOT. loester) THEN
          nfun=nfun+1 ; af(nfun)=0.309 ; bf(nfun)=-350
        ENDIF
      ENDIF

      IF (tgroup(i)(1:6)=='CO(OH)')    THEN; nfun=nfun+1 ; af(nfun)= 0.783 ; bf(nfun)= 200. ; ENDIF
      IF (INDEX(tgroup(i),'(ONO2)')/=0)THEN; nfun=nfun+1 ; af(nfun)=0.0397 ; bf(nfun)=-640. ; ENDIF
      IF (INDEX(tgroup(i),'(NO2)')/=0) THEN; nfun=nfun+1 ; af(nfun)=0.31   ; bf(nfun)= 0.   ; ENDIF
      IF (tgroup(i)(1:10)=='CO(OONO2)')THEN; nfun=nfun+1 ; af(nfun)=0.1    ; bf(nfun)= 0.   ; ENDIF


! Esters (Table 8 in Jenkin et al., 2018)
      IF (INDEX(tgroup(i),'-O-')/=0) THEN          
        DO j=1,nca
          IF ((tbond(i,j)==3).AND.(j/=ig)) THEN
            IF (tgroup(j)=='CHO') THEN              ! -CH-O-CHO
              nfun=nfun+1 ; af(nfun)=0.0251 ; bf(nfun)=-1050.
              GOTO 10
            ELSE IF (tgroup(j)=='CO') THEN          ! -CH-O-CO-
              nfun=nfun+1 ; af(nfun)=0.0310 ; bf(nfun)=-1270.
              GOTO 10
            ENDIF
          ENDIF                 
        ENDDO
      ENDIF

! FIND SUBSTITUENTS ON BETA CARBONS (Tables 5, 8 in Jenkin et al., 2018):
! for -CH2CO- (but not -CO-O-), use MULT=3.4
! for -CH2CO-O- , use MULT=2.7 
! for -CH2-O-, >CH-O-,>C(-O-R)- , use MULT=3.5 (but not for -CH2-O-CO-)
      betaloop: DO j=1,nca
        IF ((tbond(i,j)/=0 ).AND.(j/=ig)) THEN
          IF (tgroup(j)(1:6)=='CO(OH)') THEN             ! CH-C-COOH
            IF (loalk) loalk=.FALSE. 
            nfun=nfun+1 ; af(nfun)=0.0215 ; bf(nfun)=-1440.
            CYCLE betaloop
          
          ELSE IF (tgroup(j)(1:2)=='CO') THEN            ! CH-C-CO
            DO k=1,nca
              IF (k==i) CYCLE 
              IF ((tbond(j,k)==1).AND.(k/=i)) THEN       ! CH-C-CO-C
                IF (loalk) loalk=.FALSE.             
                IF (tgroup(i)(1:3)/='CO ') THEN          ! CH-C-CO-C (but not CH-CO-CO-C) 
                  nfun=nfun+1 ; af(nfun)=0.0253 ; bf(nfun)=-1460. 
                  CYCLE betaloop
                ENDIF
              ELSE IF (tbond(j,k)==3) THEN               ! CH-C-CO-O
                IF (loalk) loalk=.FALSE. 
                nfun=nfun+1 ; af(nfun)=0.0215 ; bf(nfun)=-1440. 
                CYCLE betaloop
              ENDIF
            ENDDO    

          ELSE IF (tgroup(j)(1:3)=='CHO') THEN           ! CH-C-CHO
            IF (loalk) loalk=.FALSE.
            nfun=nfun+1 ; af(nfun)=0.0253 ; bf(nfun)=-1460.
            CYCLE betaloop               
          
          ELSE IF (tgroup(j)=='-O-') THEN                ! CH-C-O-
            DO k=1,nca
              IF ((tbond(j,k)==3).AND.(k/=i)) THEN
                IF ((tgroup(i)(1:3)/='CO ').AND.(tgroup(k)(1:3)/='CO ').AND. & 
                    (tgroup(k)(1:3)/='CHO')) THEN ! CH-C-O-C
                  IF (loalk) loalk=.FALSE.         
                  nfun=nfun+1 ; af(nfun)=0.122 ; bf(nfun)=-1000.
                  CYCLE betaloop
                ENDIF
              ENDIF
            ENDDO
          ENDIF
        ENDIF
      ENDDO betaloop
        
10    CONTINUE

! select the right factor related to the H-atom abstraction reactions of OH adjacent to C=C bonds (Table 12 in Jenkin et al., 2018)
! first IF : only count contribution if the neighbor is the 1st Cd
! to avoid to count twice the contribution for Cd=Cd-CH-Cd=Cd structures
      IF (i==Cdtable(1)) THEN 
        IF (Cdtable(3)==0) THEN     ! only one C=C double bond
          IF (tgroup(Cdtable(2))(1:4)=='CdH2') THEN
            nfun=nfun+1 ; af(nfun)=1. ; bf(nfun)=-275.
          ELSE 
            nfun=nfun+1 ; af(nfun)=1. ; bf(nfun)=-545.
          ENDIF
        ELSE                          ! two C=C double bonds
          IF ((tgroup(Cdtable(2))(1:4)=='CdH2').AND.(tgroup(Cdtable(4))(1:4)=='CdH2')) THEN
            nfun=nfun+1 ; af(nfun)=1. ; bf(nfun)=-480.
          ELSE IF ((tgroup(Cdtable(2))(1:4)=='CdH2').OR.(tgroup(Cdtable(4))(1:4)=='CdH2')) THEN
            nfun=nfun+1 ; af(nfun)=1. ; bf(nfun)=-645.
          ELSE 
            nfun=nfun+1 ; af(nfun)=1. ; bf(nfun)=-750.
          ENDIF
        ENDIF
      ENDIF

! add the alkyl correction if needed 
      IF (loalk) THEN
        nfun=nfun+1 ; af(nfun)=aalk ; bf(nfun)=balk
      ENDIF

! check size (just in case)
      IF (nfun>mxfact) THEN
        WRITE(6,*) '--error--, in rabsoh. Too many contributions'
        STOP "in rabsoh"
      ENDIF

! compute AF(x) and BF(x) 
      IF (nfun==0) THEN 
        ax=1. ; bx = 0.
      ELSE IF (nfun==1) THEN
        ax=af(1) ; bx=bf(1)
      ELSE IF (nfun > 1) THEN
        ax=1. ; bx=0.
        DO j=1,nfun
          ax=ax*af(j) ; bx=bx+bf(j)
        ENDDO
        ax=ax**(1./REAL(nfun)) ; bx=bx/REAL(nfun)
      ENDIF

! add correction due to the current alpha group
      arrhc(1)=arrhc(1)*ax
      arrhc(3)=arrhc(3)+bx

      IF (wrtsarinfo) THEN !! debug
        IF (nCd>0) WRITE(saru,*) 'nCd=',nCd              
        IF (nCd>0) WRITE(saru,*) 'Cdtable=',Cdtable(:)
        IF (nbF>1) WRITE(saru,*) ' nbF=',nbF          
        WRITE(saru,*) ' ax=',ax ; WRITE(saru,*) ' bx=',bx
        WRITE(saru,*) ' arrhc(1,3)=',arrhc(1),arrhc(3)
        DO j=1,nfun
          WRITE(saru,*) 'A(X',j,')=',af(j) ; WRITE(saru,*) 'B(X',j,')=',bf(j)                    
          WRITE(saru,*) 'F(X',j,')298=',af(j)*EXP(-bf(j)/298.)
        ENDDO
      ENDIF
    ENDIF
  ENDDO alphaloop

! ADD RING FACTOR 
  IF (nring/=0) THEN
    CALL ringfac(ig,nca,tbond,tgroup,ax,bx)
    arrhc(1)=arrhc(1)*ax
    arrhc(3)=arrhc(3)+bx
  ENDIF

  IF (wrtsarinfo) CALL wrtkabs(saru,tgroup(ig),arrhc)

END SUBROUTINE rabsoh

!=======================================================================
! Purpose: compute the ring factor (ring strength) for H-abstraction 
! on a node belonging to rings. Data are from tables 3 and 6 in Jenkin et al., ACPs, 2018 
!=======================================================================
SUBROUTINE ringfac(ig,nca,tbond,tgroup,ax,bx)
  USE keyparameter, ONLY: saru
  USE ringtool, ONLY: ring_data
  USE keyflag, ONLY: wrtsarinfo
  IMPLICIT NONE

  INTEGER,INTENT(IN)         :: ig         ! group bearing the leaving H 
  INTEGER,INTENT(IN)         :: nca        ! # of nodes 
  INTEGER,INTENT(IN)         :: tbond(:,:) ! node matrix for the current species
  CHARACTER(LEN=*),INTENT(IN):: tgroup(:)  ! groups for the current species
  REAL,INTENT(OUT)           :: ax         ! overall ring A factor
  REAL,INTENT(OUT)           :: bx         ! overall ring B factor

  INTEGER,PARAMETER :: mxirg=6             ! max # of distinct rings 
  INTEGER  :: nring_ind                    ! # of distinct rings
  INTEGER  :: trackrg(mxirg,SIZE(tgroup))  ! (a,:)== track (node #) belonging ring a
  LOGICAL  :: ring_ind(mxirg,SIZE(tgroup)) ! (a,b)==true if node b belong to ring a
  INTEGER  :: i,j,k
  INTEGER  :: rgord,nCO,nether

  INTEGER  :: nrf                 ! # of contribution to ring factor
  REAL     :: af(mxirg)           ! contribution to A ring factor
  REAL     :: bf(mxirg)           ! contribution to B ring factor

  nrf = 0 ; af(:)=0. ; bf(:)=0.

  CALL ring_data(ig,nca,tbond,tgroup,nring_ind,ring_ind,trackrg) 
  IF (wrtsarinfo) WRITE(saru,*) '# of distinct cycles:',nring_ind
  ringloop: DO i=1,nring_ind
    rgord=COUNT(ring_ind(i,:).EQV..TRUE.)   ! get ring size                    

    SELECT CASE (rgord)
      CASE (3)
        nrf = nrf+1 ; af(nrf)=0.395 ; bf(nrf)=920.
        DO k=1,rgord
          IF (tgroup(trackrg(i,k))=='-O-') THEN
            af(nrf)=1. ; bf(nrf)=-298*log(0.0079)
          ENDIF
        ENDDO
        
      CASE (4)
        nether=0 ; nCO=0
        DO k=1,rgord
          IF (tgroup(trackrg(i,k))(1:3)=='-O-') nether=nether+1
          IF (tgroup(trackrg(i,k))(1:3)=='CO ') nCO=nCO+1
        ENDDO

        IF ((nCO==0).AND.(nether==0))      THEN
          nrf = nrf+1 ; af(nrf)=0.634 ; bf(nrf)=130.                  ! factor for cyclobutane (table 3)
        ELSE IF ((nCO/=0).AND.(nether==1)) THEN
          nrf = nrf+1 ; af(nrf)=1. ; bf(nrf)=-298*log(sqrt(0.5*0.08)) ! square root of the ring factors
        ELSE IF (nCO/=0)                   THEN
          nrf = nrf+1 ; af(nrf)=1. ; bf(nrf)=-298*log(0.08)
        ELSE IF (nether==1)                THEN
          nrf = nrf+1 ; af(nrf)=1. ; bf(nrf)=-298*log(0.5)
        ENDIF                 
        
      CASE (5)
        nether=0 ; nCO=0
        DO k=1,rgord
          IF (tgroup(trackrg(i,k))(1:3)=='-O-') nether=nether+1
          IF (tgroup(trackrg(i,k))(1:3)=='CO ') nCO=nCO+1
        ENDDO
  
        IF ((nCO==0).AND.(nether==0))      THEN
          nrf = nrf+1 ; af(nrf)=0.873 ; bf(nrf)=70.                    ! factor for cyclopentane (table 3)
        ELSE IF ((nCO/=0).AND.(nether==2)) THEN
          nrf = nrf+1 ; af(nrf)=1. ; bf(nrf)=-298*log(sqrt(0.32*0.59)) ! square root of the ring factors
        ELSE IF ((nCO/=0).AND.(nether==1)) THEN
          nrf = nrf+1 ; af(nrf)=sqrt(1.20) ; bf(nrf)=((-298*log(0.32))+55.)/2.
        ELSE IF (nCO/=0)                   THEN
          nrf = nrf+1 ; af(nrf)=1. ; bf(nrf)=-298*log(0.32)    
        ELSE IF (nether==1)                THEN
          nrf = nrf+1 ; af(nrf)=1.20 ; bf(nrf)=55.  
        ELSE IF (nether==2)                THEN 
          nrf = nrf+1 ; af(nrf)=1. ; bf(nrf)=-298*log(0.59)    
        ENDIF          


      CASE (6)
        nether=0 ; nCO=0
        DO k=1,rgord
          IF (tgroup(trackrg(i,k))(1:3)=='-O-') nether=nether+1
          IF (tgroup(trackrg(i,k))(1:3)=='CO ') nCO=nCO+1
        ENDDO   
 
        IF ((nCO==0).AND.(nether==0))      THEN 
          nrf = nrf+1 ; af(nrf)=0.95 ; bf(nrf)=0.                      ! factor for cycloalkane (table 3)
        ELSE IF ((nCO/=0).AND.(nether==2)) THEN 
          nrf = nrf+1 ; af(nrf)=sqrt(2.53) ; bf(nrf)=((-298*log(0.61))+535.)/2. 
        ELSE IF ((nCO/=0).AND.(nether==1)) THEN
          nrf = nrf+1 ; af(nrf)=sqrt(1.59) ; bf(nrf)=((-298*log(0.61))+290.)/2. 
        ELSE IF (nCO/=0)                   THEN 
          nrf = nrf+1 ; af(nrf)=1. ; bf(nrf)=-298*log(0.61)
        ELSE IF (nether==1)                THEN
          nrf = nrf+1 ; af(nrf)=1.59 ; bf(nrf)=290. 
        ELSE IF (nether>1)                 THEN   
          nrf = nrf+1 ; af(nrf)=2.53 ; bf(nrf)=535.
        ENDIF        

                       
      CASE (7)
        nether=0
        DO k=1,rgord
          IF (tgroup(trackrg(i,k))(1:3)=='-O-') nether=nether+1
        ENDDO

        IF (nether==1) THEN
          nrf = nrf+1 ; af(nrf)=1.61 ; bf(nrf)=190.
        ELSE IF (nether==2) THEN
          nrf = nrf+1 ; af(nrf)=1.30 ; bf(nrf)=240.
        ENDIF

    END SELECT
    
    IF (wrtsarinfo) THEN ; WRITE(saru,*) 'CYCLE ',i,' of size:',rgord ; ENDIF

  ENDDO ringloop

! check size (just in case)                                        
  IF (nrf>mxirg) THEN                                              
    WRITE(6,*) '--error--, in ringfac. Too many ring contributions'
    STOP "in ringfac"                                              
  ENDIF

! compute AF(x) and BF(x)                                          
  IF (nrf == 0) THEN
    ax=1. ; bx = 0.
  ELSE IF (nrf == 1) THEN
    ax=af(1) ; bx=bf(1)  
      IF (wrtsarinfo) THEN
        WRITE(saru,*) 'A(1)=',af(1) ; WRITE(saru,*) 'B(1)=',bf(1)    
        WRITE(saru,*) 'Fring298=',af(1)*EXP(-bf(1)/298.) 
      ENDIF                                              
  ELSE IF (nrf > 1) THEN                                 
    ax=1. ; bx=0.                                        
    DO j=1,nrf                                           
      ax=ax*af(j) ; bx=bx+bf(j)                          
      IF (wrtsarinfo) THEN                              
        WRITE(saru,*) 'A(',j,')=',af(j) ; WRITE(saru,*) 'B(',j,')=',bf(j)
        WRITE(saru,*) 'Fring298=',af(j)*EXP(-bf(j)/298.) 
      ENDIF                                              
    ENDDO                                                
  ENDIF                                                  

END SUBROUTINE ringfac 

!=======================================================================
! PURPOSE: small routine to write basic info
!=======================================================================
SUBROUTINE wrtkabs(lout,grp,arrhc)
  IMPLICIT NONE

  INTEGER :: lout
  CHARACTER(LEN=*) :: grp
  REAL :: arrhc(:)

  WRITE(lout,*)' '
  WRITE(lout,*)' 4) TOTAL RATE for abstraction at:',TRIM(grp)
  WRITE(lout,*)'  H abs by OH rate 298:',arrhc(1)*EXP(-arrhc(3)/298.)
END SUBROUTINE wrtkabs

END MODULE hotool
