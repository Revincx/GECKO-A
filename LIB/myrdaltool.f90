MODULE myrdaltool
IMPLICIT NONE
CONTAINS

! SUBROUTINE jobakgr(group,bond,nca,nring,rjg,JobGroup)
! SUBROUTINE myrdalprop(chem,bond,group,nring,rjg,weight,Tb,logPvap,deltaHvap)

!=======================================================================
! PURPOSE : compute the parameters involved in the Myrdal & Yalkowsky 
! group contribution method for vapor pressure estimates.  
!
! Joback group contributions for boiling point are given in          
! Reid et al., 1986, except for -ONO2 and -COOONO2  (see Camredon    
! and Aumont, 2006). Joback groups are picked with the following     
! index in JdeltaTb (table of DeltaTb) and Jobgroup (number of Joback
! group in chem):                                                    
!                                                                    
!  1: CH3                  2: CH2-chain           3: CH-chain        
!  4: C-chain              5: CdH2                6: CdH-chain       
!  7: Cd-chain (=Cd<)      8: F                   9: Cl              
! 10: Br                  11: (OH)               12: -O- chain       
! 13: -CO- chain          14: CHO                15: CO(OH)          
! 16: CO-O-               17: S                  18: ONO2            
! 19: COOONO2             20: CH2- ring          21: CH-ring         
! 22: C- ring             23: CdH- ring          24: Cd-ring         
! 25: -O- ring            26: CO- ring           27: phenolic OH     
! 28: NO2                 29: OOH (-O-+OH)                           
!=======================================================================
SUBROUTINE myrdalprop(chem,bond,group,nring,rjg,weight,&
                      Tb,logPvap,deltaHvap)
  USE ringtool, ONLY: findring
  USE toolbox, ONLY: stoperr
  IMPLICIT NONE

  CHARACTER(LEN=*),INTENT(IN) :: chem
  CHARACTER(LEN=*),INTENT(IN) :: group(:) ! group, without ring joining char
  INTEGER,INTENT(IN) :: bond(:,:) 
  INTEGER,INTENT(IN) :: nring  
  INTEGER,INTENT(IN) :: rjg(:,:)          ! ring-join group pairs
  REAL,INTENT(IN)    :: weight
  REAL,INTENT(OUT)   :: Tb
  REAL,INTENT(OUT)   :: logPvap
  REAL,INTENT(OUT)   :: deltaHvap

  REAL,PARAMETER :: temp=298.   ! ref temperature for Psat and latent heat
  REAL    :: tau                ! effective # of torsional bonds
  REAL    :: HBN                ! Hydrogen bond number
  INTEGER :: i, j, k, nca
  INTEGER :: sp2, sp3, indring
  INTEGER :: Jobgroup(29) 
  REAL    :: OH, COOH
  INTEGER :: rgallpath(SIZE(rjg,1),SIZE(group))
  INTEGER :: rgpath(SIZE(group)), nshare
  INTEGER :: begrg, endrg
  INTEGER :: rngflg             ! 0 = 'no', 1 = 'yes'

  CHARACTER(LEN=12),PARAMETER :: progname='myrdalprop '
  CHARACTER(LEN=70)          :: mesg

! Joback data for boiling points (see index in the subroutine header) 
  REAL, DIMENSION(29), PARAMETER :: JdeltaTb=  &
    (/23.58,  22.88,  21.74,  18.25,  18.18,  24.96,  24.14,  -0.03,  &
      38.13,  66.86,  92.88,  22.42,  76.75,  72.24, 169.09,  81.10,  &
      68.78, 112.10, 157.42,  27.15,  21.78,  21.32,  26.73,  31.01,  &
      31.22,  94.97,  76.34, 152.54, 115.30 /)

! count the number of nodes
  nca=0
  DO i=1,SIZE(group)
    IF (group(i)(1:1)/=' ') nca=nca+1
  ENDDO

! --------------------------------------------
! get Joback group      
! --------------------------------------------
  CALL jobakgr(group,bond,nca,nring,rjg,JobGroup)

! ----------------------------------   
!  Compute boiling point Tb (K), by Joback's method
! ----------------------------------
  Tb=0.
  DO i = 1,29
    Tb = Tb + JdeltaTb(i)*Jobgroup(i)
  ENDDO
  Tb = Tb + 198.

! ---------------------------------------	   
! Compute the effective number of torsional bonds (tau)
! (requires non ring sp3 and sp2 and the number of independent rings) 
! ---------------------------------------

! non ring non terminal SP3 :
! ===========================
! -O- (12), -CH2- (2), >CH- (3), >C< (4), -ONO2 (18),
!  -PAN (19x2), -CO-O- (16), -OOH (29)
  sp3 = Jobgroup(12) + Jobgroup(2)  + Jobgroup(3)  + Jobgroup(4) + &
        Jobgroup(18) + Jobgroup(19) + Jobgroup(19) + Jobgroup(16)+ Jobgroup(29) 

! non ring non terminal SP2 :
! ===========================
! -CHO (14) , -COOH (15), -CO- (13),  -CdH=(6), -CO-O- (16),
! -ONO2 (18), -PAN (19x2), >Cd= (7), -NO2(28)
  sp2 = Jobgroup(14) + Jobgroup(15) + Jobgroup(13) + Jobgroup(6) + Jobgroup(16) + &
        Jobgroup(18) + Jobgroup(19) + Jobgroup(19) + Jobgroup(7) + Jobgroup(28)

! number of independent rings in the system 
! =========================================
  indring=nring

! If more than 2 rings, stop (need more checks ....)
  IF (nring > SIZE(rjg,1)) THEN
    mesg="Ring number > mri"
    CALL stoperr(progname,mesg,chem)
  ENDIF 

! If more than 1 ring, check whether 2 rings are independent or not
  IF (nring==2) THEN
    rgallpath(:,:)=0

! find the nodes that belong to a ring
    DO i=1,nring
       begrg=rjg(i,1)  ;  endrg=rjg(i,2)
       CALL findring(begrg,endrg,nca,bond,rngflg,rgpath)
       DO j=1,nca
         IF (rgpath(j)==1) rgallpath(i,j)=1
       ENDDO
    ENDDO

! if 2 nodes are shared by 2 rings => rings are not independent
    indring=1
    DO i=1,nring-1
      DO j=i+1,nring
        nshare=0
        DO k=1,nca
          IF (rgallpath(i,k)+rgallpath(i,k)==2) nshare=nshare+1
        ENDDO
        IF (nshare<=1) indring=indring+1
      ENDDO
    ENDDO
  ENDIF

! BA: need to check if indring if fine when more than 2 rings
  IF (nring>2) THEN
    PRINT*, TRIM(chem)
    PRINT*, 'indring=', indring
  ENDIF

! Compute tau 
! ============
  tau=sp3+0.5*sp2+0.5*indring-1.
  IF (tau<0) tau=0

! ---------------------------------------	   
! Compute the Hydrogen Bond Number (HBN) 
! OH = -OH(11) + -OOH(29); COOH = -COOH(15) ; NH2 = 0
! ---------------------------------------
  OH   = Jobgroup(11) + Jobgroup(29)
  COOH = Jobgroup(15)
  HBN=(sqrt(OH+COOH))/weight

! ---------------------------------------
!  compute the vapor pressure (log10) and latent heat 
! ---------------------------------------
  logPvap = -(86. + 0.4*tau + 1421*HBN)*(Tb-temp)/(19.1*temp) + &
               (-90.0-2.1*tau)/19.1*((Tb-temp)/temp-LOG(Tb/temp))
  deltaHvap = Tb*(86. + 0.4*tau + 1421*HBN) +  (-90.0-2.1*tau)*(temp-Tb)
END SUBROUTINE myrdalprop

! ======================================================================
! PURPOSE: Find the number of group increments for the Joback method 
! for Tb estimation.          
!                                                                    
! Joback group contribution for boiling point (Tb) are given in      
! Reid et al., 1986, except for -ONO2 and PAN (see Camredon & Aumont, 
! atmos. env. 2006). Joback groups are picked with the   
! following index in JdeltaTb (table of DeltaTb) and Jobgroup (number
! of Joback group in chem):                                          
!                                                                    
!  1: CH3                  2: CH2-chain           3: CH-chain        
!  4: C-chain              5: CdH2                6: CdH-chain       
!  7: Cd-chain (=Cd<)      8: F                   9: Cl              
! 10: Br                  11: (OH)               12: -O- chain       
! 13: -CO- chain          14: CHO                15: CO(OH)          
! 16: CO-O-               17: S                  18: ONO2            
! 19: COOONO2             20: CH2- ring          21: CH-ring         
! 22: C- ring             23: CdH- ring          24: Cd-ring         
! 25: -O- ring            26: CO- ring           27: phenolic OH     
! 28: NO2                 29: OOH (-O-+OH)                           
! ======================================================================
SUBROUTINE jobakgr(group,bond,nca,nring,rjg,JobGroup)
  USE ringtool, ONLY: findring
  IMPLICIT NONE
   
  CHARACTER(LEN=*),INTENT(IN) :: group(:) ! group, without ring joining char
  INTEGER,INTENT(IN) :: bond(:,:) 
  INTEGER,INTENT(IN) :: nca               ! # of nodes  
  INTEGER,INTENT(IN) :: nring  
  INTEGER,INTENT(IN) :: rjg(:,:)          ! ring-join group pairs
  INTEGER,INTENT(OUT):: Jobgroup(29) 
 
  INTEGER    i, j, k
  INTEGER    begrg,endrg
  INTEGER    ring(SIZE(group))     ! rg index for nodes, current ring (0=no, 1=yes)
  INTEGER    indexrg(SIZE(group))  ! rg index for nodes, any ring (0=no, 1=yes)
  INTEGER    rngflg        ! 0 = 'no ring', 1 = 'yes ring'
  INTEGER    nc, ibeg, eflg
  INTEGER    tbond(SIZE(group),SIZE(group))

  Jobgroup(:)=0  ;  indexrg(:)=0  ;  tbond(:,:)=bond(:,:)

! Stop if the number of rings is greater than mxring! (need additional checks)
  IF (nring>SIZE(rjg,1)) THEN
    WRITE(6,*) '--error-- in jobakgr. nring > mxring'
    STOP "in jobakgr"
  ENDIF

! IF RINGS EXIST, THEN FIND THE NODES THAT BELONG TO THE RINGS
! --------------------------------------------------------------
  IF (nring>0) THEN

! find the nodes that belong to a ring
    DO i=1,nring
       begrg=rjg(i,1)  ;  endrg=rjg(i,2)
       CALL findring(begrg,endrg,nca,bond,rngflg,ring)
       DO j=1,nca
         IF (ring(j)==1) indexrg(j)=1
       ENDDO
    ENDDO
  ENDIF

! LOOP OVER THE GROUPS AND FIND JOBACK'S GROUPS
! ----------------------------------------------
! Corresponding index for the joback group
! 1: CH3                  2: CH2-chain           3: CH-chain         *
! 4: C-chain              5: CdH2                6: CdH-chain
! 7: Cd-chain (=Cd<)      8: F                   9: Cl
!10: Br                  11: (OH)               12: -O- chain
!13: -CO- chain          14: CHO                15: CO(OH)
!16: CO-O-               17: S                  18: ONO2
!19: COOONO2             20: CH2- ring          21: CH-ring
!22: C- ring             23: CdH- ring          24: Cd-ring
!25: -O- ring            26: CO- ring           27: phenolic OH
!28: NO2                 29: OOH (=OH+-O-)

  grloop: DO i=1,nca 

! oxygen increments -CH=O,-COOH,CO(OONO2),-O-,>C=O
    IF (group(i)(1:4)=='CHO ') THEN
      Jobgroup(14)=Jobgroup(14)+1 
      CYCLE grloop
    ELSE IF (group(i)(1:7)=='CO(OH) ') THEN
      Jobgroup(15)=Jobgroup(15)+1              
      CYCLE grloop
    ELSE IF (group(i)(1:10)=='CO(OONO2) ') THEN
      Jobgroup(19)=Jobgroup(19)+1              
      CYCLE grloop

    ELSE IF (group(i)(1:4)=='-O- ') THEN
      IF (indexrg(i)==0) THEN
        eflg=0 
        DO j=1,nca
          IF (tbond(i,j)==3) THEN
            IF (group(j)(1:3)=='CO ') THEN
              eflg=eflg+1
              DO k=1,nca  ! loop to avoid double count for -O-CO-O-
                IF (tbond(j,k)==3) THEN
                  tbond(j,k)=0  ;  tbond(k,j)=0
                ENDIF 
              ENDDO
            ENDIF
          ENDIF
        ENDDO
        IF (eflg>=1) THEN 
          Jobgroup(16)=Jobgroup(16)+1
          Jobgroup(13)=Jobgroup(13)-1
        ELSE 
          Jobgroup(12)=Jobgroup(12)+1 
        ENDIF 
      ELSE            
         Jobgroup(25)=Jobgroup(25)+1
      ENDIF  
      CYCLE grloop

    ELSE IF (group(i)(1:2)=='CO') THEN
      IF (indexrg(i)==0) THEN 
        Jobgroup(13)=Jobgroup(13)+1           
      ELSE            
        Jobgroup(26)=Jobgroup(26)+1
      ENDIF   
      IF (group(i)(1:3)=='CO ') CYCLE grloop

! alkane-alkene-aromatic -CH3,-CH2-,=CH2-,-CH<,=CH-,>C<,=C< 
    ELSE IF (group(i)(1:4)=='CH3 ') THEN
      Jobgroup(1)=Jobgroup(1)+1   
      CYCLE grloop

    ELSE IF (group(i)(1:3)=='CH2') THEN
      IF (indexrg(i)==0) THEN 
        Jobgroup(2)=Jobgroup(2)+1           
      ELSE            
        Jobgroup(20)=Jobgroup(20)+1
      ENDIF
      IF (group(i)(1:4)=='CH2 ') CYCLE grloop

    ELSE IF (group(i)(1:5)=='CdH2 ') THEN
      Jobgroup(5)=Jobgroup(5)+1            
      CYCLE grloop

    ELSE IF (group(i)(1:2)=='CH') THEN         
      IF (indexrg(i)==0) THEN 
        Jobgroup(3)=Jobgroup(3)+1           
      ELSE            
        Jobgroup(21)=Jobgroup(21)+1
      ENDIF
      IF (group(i)(1:3)=='CH ') CYCLE grloop

    ELSE IF (group(i)(1:3)=='CdH') THEN
      IF (indexrg(i)==0) THEN 
        Jobgroup(6)=Jobgroup(6)+1           
      ELSE            
        Jobgroup(23)=Jobgroup(23)+1
      ENDIF           
      IF (group(i)(1:4)=='CdH ') CYCLE grloop

    ELSE IF (group(i)(1:3)=='Cd ') THEN
      IF (indexrg(i)==0) THEN 
        Jobgroup(7)=Jobgroup(7)+1           
      ELSE            
        Jobgroup(24)=Jobgroup(24)+1
      ENDIF        

    ELSE IF (group(i)(1:1)=='C') THEN
      IF (indexrg(i)==0) THEN 
        Jobgroup(4)=Jobgroup(4)+1           
      ELSE            
        Jobgroup(22)=Jobgroup(22)+1
      ENDIF          
      IF (group(i)(1:2)=='C ') CYCLE grloop

! aromatics
    ELSE IF (group(i)(1:2)=='cH') THEN
      Jobgroup(23)=Jobgroup(23)+1
      IF (group(i)(1:3)=='cH ') CYCLE grloop

    ELSE IF (group(i)(1:1)=='c') THEN
      Jobgroup(24)=Jobgroup(24)+1
      IF (group(i)(1:5)=='c(OH)') THEN
        Jobgroup(27)=Jobgroup(27)+1
        CYCLE grloop
      ENDIF               
    ENDIF
   
! SEARCH FOR FUNCTIONALTIES IN ()
    nc=INDEX(group(i),' ')

! alcohol increments (OH) but not carboxylic acid. Distinguish
! between primary and secondary alcohol for Stein and Brown correction
    ibeg=INDEX(group(i),'(OH)')
    IF (ibeg/=0) THEN 
      DO j=ibeg,nc-3
        IF ((group(i)(j:j+3)=='(OH)').AND.(group(i)(1:6)/='CO(OH)')) THEN  
           Jobgroup(11)=Jobgroup(11)+1
        ENDIF
      ENDDO
    ENDIF   

! hydroperoxides (OOH) increments   
    ibeg=INDEX(group(i),'(OOH)')
    IF (ibeg/=0) THEN
      DO j=ibeg,nc-4
        IF (group(i)(j:j+4)=='(OOH)') THEN
          Jobgroup(29)=Jobgroup(29)+1
        ENDIF
      ENDDO
    ENDIF
  
! (OOOH) increments   
    ibeg=INDEX(group(i),'(OOOH)')
    IF (ibeg/=0) THEN
      DO j=ibeg,nc-5
        IF (group(i)(j:j+5)=='(OOOH)') THEN
          Jobgroup(29)=Jobgroup(29)+1
        ENDIF
      ENDDO
    ENDIF
  
! (ONO2) increments
    ibeg=INDEX(group(i),'(ONO2)')
    IF (ibeg/=0) THEN
      DO j=ibeg,nc-5
        IF (group(i)(j:j+5)=='(ONO2)') THEN
          Jobgroup(18)=Jobgroup(18)+1
        ENDIF
      ENDDO
    ENDIF

    ibeg=INDEX(group(i),'(NO2)')
    IF (ibeg/=0) THEN
      DO j=ibeg,nc-4
        IF (group(i)(j:j+5)=='(NO2)') THEN
          Jobgroup(28)=Jobgroup(28)+1
        ENDIF
      ENDDO
    ENDIF
  ENDDO grloop
     
END SUBROUTINE jobakgr

END MODULE myrdaltool
