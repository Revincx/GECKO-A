MODULE o3addtool
  IMPLICIT NONE
  CONTAINS

!=======================================================================
! PURPOSE: Estimate rate constant according to Jenkin 2020 report/paper 
! (MAGNIFY project) for a non conjugated C=C double-bond. 
!=======================================================================
SUBROUTINE o3rate_mono(bond,zebond,group,ngr,cdtable,cdsub,arrhc)
  USE keyflag, ONLY: wrtsarinfo
  USE keyparameter, ONLY: saru
  USE ringtool, ONLY: findring
  IMPLICIT NONE

  CHARACTER(LEN=*),INTENT(IN) :: group(:)   ! group matrix
  INTEGER,INTENT(IN)          :: bond(:,:)  ! node matrix
  INTEGER,INTENT(IN)          :: zebond(:,:) ! cis/trans info on C=C bond
  INTEGER,INTENT(IN)          :: ngr        ! # of groups
  INTEGER,INTENT(IN)          :: cdtable(:) ! node # bearing a "Cd"
  INTEGER,INTENT(IN)          :: cdsub(:)   ! # of C (including CO) bonded to Cd in cdtable(i)
  REAL,INTENT(OUT)            :: arrhc(:)

  INTEGER  :: i,j,ifct,Ci,Cf,cd1,cd2
  INTEGER  :: ring(SIZE(group)),rngflg
  REAL     :: mult,Fring
  REAL     :: tmparrhc(3),k298

  arrhc(:)=0 ; mult=1.0 ; ifct=0 ; Fring=1 ; k298=0.
  cd1=cdtable(1) ; cd2=cdtable(2)
      
  IF (wrtsarinfo) THEN
    WRITE(saru,*) '  ' ;  WRITE(saru,*) ' ======= O3RATE MONO ====== '
  ENDIF

  CALL findring(cd1,cd2,ngr,bond,rngflg,ring)

! Use reference rate coefficients for acyclic vinylic compounds  
! >C=C(R)-O- and >C=C(R)-CO- (see Table 6-8 in Jenkin et al., 2020)
! Loop over all vinylic structure in alpha position of the C=C bond
! and if more than one vinylic function, keep the less reactive rate.
! For vinylic oxygenated compounds, do not apply substituent factors (RETURN).
  IF (rngflg==0) THEN
    ifct=0 
    arrhc(1)=2E-10 ; arrhc(2:3)=0. ! maximum rate constant
    DO i=1,ngr ; DO j=1,2
      IF ( ((bond(cdtable(j),i)==1).AND.((group(i)(1:2)=='CO').OR.(group(i)(1:3)=='CHO')))    &
           .OR.(bond(cdtable(j),i)==3) ) THEN                ! a vinylic structure is found here
        IF (j==1) THEN ; Ci=cd1 ; Cf=cd2            ! set Ci as Cd node attached to CO or -O ...   
        ELSE           ; Ci=cd2 ; Cf=cd1            ! and Cf as the opposite node of the C=C
        ENDIF
        ifct=i                                               ! flag raised - node with vinylic O
        CALL ref_rate_mono(bond,zebond,group,ngr,cdsub,tmparrhc,Cf,Ci,ifct)
        arrhc(1)=MIN(arrhc(1),tmparrhc(1))                   ! keep only the smallest rate
      ENDIF
    ENDDO ; ENDDO
    IF (ifct/=0) RETURN                                      ! vinyl O found, return
  ENDIF   

! other monoalkenes - assign the rate constant (table 1 in Jenkin et al., 2020)
  Ci=cd1 ; Cf=cd2
  CALL ref_rate_mono(bond,zebond,group,ngr,cdsub,arrhc,Cf,Ci,ifct) 

! compute substituent and ring factors (table 2 and 3 in Jenkin et al., 2020)
  CALL subs_fact(group,bond,ngr,ring,cdtable,mult)    
  IF (rngflg==1) THEN
    CALL ring_fact_mono(group,bond,cd1,ngr,Fring)
    
    DO i=1,ngr ; DO j=1,2        ! for cyclic vinyl O, reduce by 50 for each -CO- subst.
      IF ((bond(i,cdtable(j))==1).AND.((group(i)=='CO').OR.(group(i)=='CHO'))) THEN
        Fring=Fring/50.
      ENDIF
    ENDDO ; ENDDO
  ENDIF

! Apply factors to the reference rate coefficient
  arrhc(1)=arrhc(1)*mult*Fring
     
! write to debug file
  IF (wrtsarinfo) THEN
    WRITE(saru,*) '=> Non-conjugated' ;
    WRITE(saru,*) 'mult=',mult
    IF (rngflg==1) WRITE(saru,*) 'Fring=',Fring
    WRITE(saru,*) 'arrhc=',arrhc(1:3)
    WRITE(saru,*) 'k298=',k298
    WRITE(saru,*) ' ' ;  WRITE(saru,*) ' '
  ENDIF

END SUBROUTINE o3rate_mono

!=======================================================================
! PURPOSE: Estimate rate constant according to Jenkin 2020 report/paper 
! (MAGNIFY project) for conjugated C=C-C=C bonds. 
!=======================================================================
SUBROUTINE o3rate_conj(bond,group,ngr,cdtable,cdsub,arrhc)
  USE keyflag, ONLY: wrtsarinfo
  USE keyparameter, ONLY: saru
  USE ringtool, ONLY: findring
  IMPLICIT NONE

  CHARACTER(LEN=*),INTENT(IN) :: group(:)  ! group matrix
  INTEGER,INTENT(IN)   :: bond(:,:)        ! node matrix
  INTEGER,INTENT(IN)   :: ngr              ! # of groups
  INTEGER,INTENT(IN)   :: cdtable(:)       ! node # bearing a "Cd"
  INTEGER,INTENT(IN)   :: cdsub(:)         ! # of C (including CO) bonded to Cd in cdtable(i)
  REAL,INTENT(OUT)     :: arrhc(:)

  INTEGER  :: i,j,nbF
  INTEGER  :: ring(SIZE(group)),rngflg
  REAL     :: mult,Fring

  arrhc(:)=0 ; mult=1.0 ; Fring=1
      
  IF (wrtsarinfo) THEN
    WRITE(saru,*) '  '
    WRITE(saru,*) ' ====== O3RATE CONJ ====== '
  ENDIF

! get reference rate constant 
  CALL ref_rate_conj(cdtable,cdsub,arrhc)
    
! look for rings
  CALL findring(cdtable(1),cdtable(2),ngr,bond,rngflg,ring)

! look for conjugated vinyl oxygenates
  nbF=0
  DO i=1,ngr ; DO j=1,4
    IF (((bond(i,cdtable(j))==1).AND.((group(i)(1:2)=='CO').OR.  &
         (group(i)=='CHO'))).OR.(bond(i,cdtable(j))==3)) THEN
      nbF=nbF+1
    ENDIF
  ENDDO ; ENDDO
    
! Substituent factors for beta groups (if no vinylic oxygenates)
  !IF (nbF==0) THEN
    CALL subs_fact(group,bond,ngr,ring,cdtable,mult)
    arrhc(1)=arrhc(1)*mult**0.5
  !ENDIF

! For conjugated vinylic O, reduce by a factor 50 for each ketone substituent
  IF (nbF/=0) THEN 
    DO i=1,ngr ; DO j=1,4 
      IF ((bond(i,cdtable(j))==1).AND.((group(i)=='CO').OR.(group(i)=='CHO'))) THEN
        arrhc(1)=arrhc(1)/50.
      ENDIF
    ENDDO ; ENDDO
  ENDIF            

! Use conjugated ring factors if all C=C-C=C structure is in a cycle ...
  IF ((ring(cdtable(1))==1).AND.(ring(cdtable(2))==1).AND. &
      (ring(cdtable(3))==1).AND.(ring(cdtable(4))==1)) THEN 
    CALL ring_fact_conj(group,bond,cdtable,ngr,Fring)

! ... otherwise use ring factors for C=C bond in a cycle
  ELSE IF ((ring(cdtable(1))==1).AND.(ring(cdtable(2))==1)) THEN
    CALL ring_fact_mono(group,bond,cdtable(1),ngr,Fring)       
  ELSE IF ((ring(cdtable(3))==1).AND.(ring(cdtable(4))==1)) THEN
    CALL ring_fact_mono(group,bond,cdtable(3),ngr,Fring)       
  ENDIF
  arrhc(1)=arrhc(1)*Fring
  
! write to debug file
  IF (wrtsarinfo) THEN
    WRITE(saru,*) '=> conjugated'
    DO i=1,ngr ;  WRITE(saru,*) 'cdsub(',i,')=',cdsub(i) ; ENDDO
    WRITE(saru,*) ' mult =',mult
    IF (rngflg==1) WRITE(saru,*) 'Fring=',Fring      
    WRITE(saru,*) 'arrhc=',arrhc(1:3)
    WRITE(saru,*) ' ' ;  WRITE(saru,*) ' '
  ENDIF

END SUBROUTINE o3rate_conj

!=======================================================================
! Purpose: Assign reference rate coefficients, see tables 1,6-8 in 
! Jenkin et al., 2020. If the C=C bond have vinylic O, then input is  
! such that the numbering  is >C[f]=C[i]-CO-,  >C[f]=C[i]-O-R (vinylic 
! side is Ci).
!=======================================================================
SUBROUTINE ref_rate_mono(bond,zebond,group,ngr,cdsub,arrhc,Cf,Ci,ifct)
  USE keyparameter, ONLY: saru,mxcp
  USE keyflag, ONLY: wrtsarinfo
  USE mapping, ONLY: gettrack
  IMPLICIT NONE
  
  CHARACTER(LEN=*),INTENT(IN) :: group(:)  ! group matrix
  INTEGER,INTENT(IN) :: bond(:,:)  ! node matrix
  INTEGER,INTENT(IN) :: zebond(:,:) ! cis/trans info on C=C bond
  INTEGER,INTENT(IN) :: ngr        ! # of groups
  INTEGER,INTENT(IN) :: cdsub(:)   ! # of C (including CO) bonded to Cd in cdtable(i)
  REAL,INTENT(OUT)   :: arrhc(:)
  INTEGER,INTENT(IN) :: Ci         ! 1st Cd node (bearing vinyl O if any)
  INTEGER,INTENT(IN) :: Cf         ! 2nd Cd node 
  INTEGER,INTENT(IN) :: ifct       ! index of node bearing vinyl O if any

  INTEGER  :: i,j,k,nb,ncarb
  REAL     :: alpha,fact
  INTEGER  :: track(mxcp,SIZE(bond,1))
  INTEGER  :: trlen(mxcp),ntr
  LOGICAL  :: ch_eff(SIZE(group))
  LOGICAL  :: ester
  
  arrhc(:)=0.

! check 1st for vinylic O
! -----------------------
  IF (ifct /=0) THEN
! Reference rate coefficients for acyclic vinyl aldehyde, ketone, esters, 
! acids and ethers. Table 6-8 in Jenkin et al., 2020
  
! Cd-CO- and Cd-CO-O
    IF (group(ifct)(1:3)=='CO ') THEN
      ester=.FALSE.
      DO i=1,ngr
        IF ((bond(ifct,i)==3).AND.(ifct/=Ci)) THEN 
          ester=.TRUE. ; EXIT
        ENDIF
      ENDDO

      IF (ester) THEN
        IF ((group(Ci)(1:3)=='CdH').AND.(group(Cf)(1:4)=='CdH2')) THEN
               arrhc(1)=0.15E-17 
        ELSE ; arrhc(1)=0.65E-17
        ENDIF
      ELSE
        IF (group(Cf)(1:4)=='CdH2') THEN
          IF (group(Ci)(1:3)=='CdH') THEN    ; arrhc(1)=0.52E-17 
          ELSE                               ; arrhc(1)=1.20E-17  
          ENDIF
        ELSE IF (group(Cf)(1:3)=='CdH') THEN ; arrhc(1)=3.90E-17 
        ELSE                                 ; arrhc(1)=0.83E-17
        ENDIF
      ENDIF
  
! Cd-CHO
    ELSE IF (group(ifct)(1:3)=='CHO') THEN
      IF (group(Ci)(1:3)=='CdH') THEN
        IF (group(Cf)(1:3)=='CdH') THEN  ; arrhc(1)=0.14E-17 
        ELSE                             ; arrhc(1)=0.18E-17
        ENDIF
      ELSE
        IF (group(Cf)(1:4)=='CdH2') THEN ; arrhc(1)=0.12E-17 
        ELSE                             ; arrhc(1)=0.57E-17
        ENDIF
      ENDIF
       
! Cd-O-R and Cd-O-CO
    ELSE IF (group(ifct)(1:3)=='-O-') THEN
      ester=.FALSE.
      DO i=1,ngr
        IF ((bond(ifct,i)==3).AND.(group(i)(1:3)=='CO ')) THEN 
          ester=.TRUE. ; EXIT
        ENDIF      
      ENDDO

      IF (ester) THEN
        IF (group(Ci)(1:3)=='CdH') THEN  ; arrhc(1)=0.32E-17 
        ELSE                             ; arrhc(1)=0.054E-17
        ENDIF
      ELSE
        IF (group(Ci)(1:4)=='CdH ') THEN 
          IF (group(Cf)(1:4)=='CdH2') THEN ; arrhc(1)=17E-17 
          ELSE                             ; arrhc(1)=42E-17
          ENDIF
        ELSE
          arrhc(1)=1.3E-17 
        ENDIF
      ENDIF
  
! CO(OH)      
    ELSE IF (group(ifct)(1:6)=='CO(OH)') THEN
      arrhc(1)=0.23E-17

! CO(OOH) and CO(OONO2) see footnotes in Table 7     
    ELSE IF ((group(ifct)(1:7)=='CO(OOH)').OR.(group(ifct)(1:9)=='CO(OONO2)')) THEN
      IF ((group(Ci)(1:3)=='CdH').AND.(group(Cf)(1:4)=='CdH2')) THEN
             arrhc(1)=0.14E-17 
      ELSE ; arrhc(1)=0.65E-17
      ENDIF
    ENDIF

! multiple vinylic oxygenated substituents
    IF ((group(ifct)=='CO ').OR.(group(ifct)=='CHO')) THEN
      DO i=1,ngr
        IF ((bond(Cf,i)==1).AND.((group(i)=='CO ').OR.(group(i)=='CHO')) ) THEN
          arrhc(1)=0.50E-17
        ENDIF
      ENDDO
    ENDIF
    IF (group(ifct)=='-O-') THEN
      DO i=1,ngr
        IF ((bond(Cf,i)==3).OR.((bond(Ci,i)==3).AND.(i/=ifct))) THEN  ! -O-C=C-O- and C=C(-O-)-O-
          arrhc(1)=48E-17
        ENDIF
      ENDDO
    ENDIF

! chain length effect: loop over the R substituents and count number of C
    alpha=0.19

! substituent at Cf side
    fact=1.
    DO j=1,ngr
      ncarb=0
      IF ((bond(Cf,j)==1).OR.(bond(Cf,j)==3)) THEN
        CALL gettrack(bond,j,ngr,ntr,track,trlen)
        ch_eff(:)=.FALSE.
        DO i=1,ntr
          IF (track(i,2)==Cf) CYCLE
          DO k=1,trlen(i)
            IF (group(track(i,k))(1:3)/='-O-') ch_eff(track(i,k))=.TRUE. 
          ENDDO
          IF (group(track(i,1))(1:3)=='CO ') ch_eff(track(i,1))=.FALSE.
        ENDDO  
        ncarb=COUNT(ch_eff.EQV. .TRUE.)
    
        IF (ncarb>=1) fact=fact+alpha*(ncarb-1)
        
        IF ((wrtsarinfo).AND.(ncarb/=0)) THEN
          WRITE(saru,*) 'group(ifct)= ',group(ifct)
          WRITE(saru,*) 'chain length effect of substituent of Cf:',group(Cf)
          WRITE(saru,*) 'arrhc(1)=',arrhc(1)
          WRITE(saru,*) 'fact=',fact
          WRITE(saru,*) 'alpha =',alpha
          WRITE(saru,*) 'ncarb=',ncarb
        ENDIF
      ENDIF
    ENDDO
    arrhc(1)=arrhc(1)*fact

! substituent at Ci side
    fact=1.
    DO j=1,ngr
      IF (j==ifct) CYCLE
      ncarb=0
      IF (bond(Ci,j)==1) THEN
        CALL gettrack(bond,j,ngr,ntr,track,trlen)
        ch_eff(:)=.FALSE.
        DO i=1,ntr
          IF (track(i,2)==Ci) CYCLE
          ch_eff(track(i,1:trlen(i)))=.TRUE. 
        ENDDO  
        ncarb=COUNT(ch_eff.EQV. .TRUE.)
    
        IF (ncarb>=1) fact=fact+alpha*(ncarb-1)
        
        IF ((wrtsarinfo).AND.(ncarb/=0)) THEN
          WRITE(saru,*) 'chain length effect of substituent of Ci:',group(Ci)
          WRITE(saru,*) 'fact=',fact
          WRITE(saru,*) 'alpha =',alpha
          WRITE(saru,*) 'ncarb=',ncarb
        ENDIF
      ENDIF
    ENDDO
    arrhc(1)=arrhc(1)*fact

! substituent of ifct
    fact=1.
    DO j=1,ngr
      IF (j==Ci) CYCLE
      ncarb=0
      IF ((bond(ifct,j)==1).OR.(bond(ifct,j)==3)) THEN
        CALL gettrack(bond,j,ngr,ntr,track,trlen)
        ch_eff(:)=.FALSE.
        DO i=1,ntr
          IF (track(i,2)==ifct) CYCLE
          DO k=1,trlen(i)
            IF (group(track(i,k))(1:3)/='-O-') ch_eff(track(i,k))=.TRUE. 
          ENDDO
        ENDDO  
        ncarb=COUNT(ch_eff.EQV. .TRUE.)
        
        IF ((bond(ifct,j)==3).AND.(group(j)(1:3)=='CO ')) ncarb=ncarb-1 ! don't count ester C in chain length
        IF (ncarb>=1) fact=fact+alpha*(ncarb-1)
        
        IF ((wrtsarinfo).AND.(ncarb/=0)) THEN
          WRITE(saru,*) 'chain length effect of substituent of ',group(ifct)
          WRITE(saru,*) 'arrhc(1)=',arrhc(1)
          WRITE(saru,*) 'fact=',fact
          WRITE(saru,*) 'alpha =',alpha
          WRITE(saru,*) 'ncarb=',ncarb
        ENDIF
      ENDIF
    ENDDO
    arrhc(1)=arrhc(1)*fact

! Monoalkene without vinylic O
! ----------------------------
  ELSE

! Reference rate coefficients for generic monoalkenes. Table 1 in Jenkin et al., 2020
    nb=cdsub(Ci)+cdsub(Cf)
    IF      (nb==0) THEN ; arrhc(1)=9.14E-15 ; arrhc(3)=2580.  ! ethene
    ELSE IF (nb==1) THEN ; arrhc(1)=2.91E-15 ; arrhc(3)=1690.  ! -CH=CH2
    ELSE IF (nb==2) THEN                                       ! -CH=CH- or >C=CH2
      IF (cdsub(Ci)==1) THEN 
        IF (zebond(Ci,Cf)==1) THEN
          arrhc(1)=3.39E-15 ; arrhc(3)=995.  ! trans value
        ELSEIF (zebond(Ci,Cf)==2) THEN
          arrhc(1)=7.29E-15 ; arrhc(3)=1120. ! cis value
        ELSE
          arrhc(1)=1.50E-16 ; arrhc(3)=0.    ! k298 mean value between cis and trans if not specified
        ENDIF
      ELSE                   ; arrhc(1)=4.00E-15 ; arrhc(3)=1685.
      ENDIF
    ELSE IF (nb==3) THEN ; arrhc(1)=7.61E-15 ; arrhc(3)=830.   ! >C=CH- 
    ELSE IF (nb==4) THEN ; arrhc(1)=3.00E-15 ; arrhc(3)=300.   ! >C=C< 
    ELSE
      WRITE(6,*) '-error1-- in ref_rate, no rate constant found '
      WRITE(6,*) 'ref_rate, nb=',nb 
      STOP "in ref_rate"
    ENDIF
  ENDIF

END SUBROUTINE ref_rate_mono

!=======================================================================
! Compute substituent factors to apply on C=C and C=C-C=C structure for
! O3+alkene reaction rates as described by Jenkin et al., 2020 (see 
! factors in table 2). 
!=======================================================================
SUBROUTINE subs_fact(group,bond,ngr,ring,cdtable,mult)
  IMPLICIT NONE
  CHARACTER(LEN=*),INTENT(IN) :: group(:)   ! group matrix
  INTEGER,INTENT(IN)          :: bond(:,:)  ! node matrix
  INTEGER,INTENT(IN)          :: ngr        ! # of groups
  INTEGER,INTENT(IN)          :: ring(:)    ! =1 if the carbon belong to a ring
  INTEGER,INTENT(IN)          :: cdtable(:) ! node # bearing a "Cd"
  REAL,INTENT(OUT)            :: mult       ! substituent factor to apply

  INTEGER   :: i,j,k,l,nalkyl,nalkyl_cyc

  mult=1.0
  DO i=1,SIZE(cdtable)  
    IF (cdtable(i)==0) CYCLE
    DO j=1,ngr
      IF (bond(cdtable(i),j)==1) THEN
        IF ((INDEX(group(j),'(OH)')/=0).AND.(group(j)/='CO(OH)')) mult=mult*1.4
        IF (INDEX(group(j),'(ONO2)') /=0)  mult=mult*0.044
        nalkyl=0 ; nalkyl_cyc=0
        betaloop : DO k=1,ngr
          IF ((bond(k,j)==1).AND.(k/=(cdtable(i)))) THEN
            IF ((group(k)(1:4)=='CH3 ').OR.(group(k)(1:3)=='CH2').OR.  &
                (group(k)(1:3)=='CH ').OR.(group(k)(1:3)=='CH(').OR.(group(k)(1:2)=='C ' )) THEN
              IF (ring(k)==0) THEN 
                nalkyl=nalkyl+1
              ELSE
                nalkyl_cyc = nalkyl_cyc + 1
              ENDIF
            ELSE IF (group(k)(1:3) =='CO ') THEN ; mult=mult*0.32
            ELSE IF (group(k)(1:3) =='CHO') THEN ; mult=mult*0.32
            ENDIF
          ELSE IF (bond(k,j)==3) THEN
            DO l=1,ngr
              IF ((bond(l,k)==3).AND.(group(l)(1:3)=='CO ')) THEN
                mult=mult*0.25 ; EXIT betaloop
              ENDIF
            ENDDO
            mult=mult*0.6
          ENDIF
        ENDDO betaloop
        IF (nalkyl>=1) THEN
          IF (nalkyl_cyc==0) THEN 
            mult=mult*(0.54**(nalkyl-1))
          ELSE
            mult=mult*(0.54**(nalkyl))
          ENDIF
        ENDIF
      ENDIF    
    ENDDO
  ENDDO
END SUBROUTINE subs_fact

!=======================================================================
! Compute ring factors to apply for monoalkenes structure for the 
! O3+alkene reaction rates as described by Jenkin et al., 2020 (see 
! factors in table 3).
!=======================================================================
SUBROUTINE ring_fact_mono(group,bond,cd1,ngr,Fring)
  USE ringtool, ONLY: ring_data
  IMPLICIT NONE
  CHARACTER(LEN=*),INTENT(IN) :: group(:)             ! group matrix
  INTEGER,INTENT(IN)  :: bond(:,:)                    ! node matrix
  INTEGER,INTENT(IN)  :: cd1                          ! node # bearing a Cd
  INTEGER,INTENT(IN)  :: ngr                          ! # of groups
  REAL,INTENT(OUT)    :: Fring                        ! substituent factor to apply
											         
  INTEGER,PARAMETER   :: mxirg=6                      ! max # of distinct rings 
  INTEGER             :: i,rgord
  REAL                :: mult
  INTEGER             :: nring_ind                    ! # of disctinct rings
  INTEGER             :: trackrg(mxirg,SIZE(group))   ! (a,:)== track (node #) belonging ring a
  LOGICAL             :: ring_ind(mxirg,SIZE(group))  ! (a,b)==true if node b belong to ring a

  CALL ring_data(cd1,ngr,bond,group,nring_ind,ring_ind,trackrg) 
  
  Fring=1.0
  ringloop: DO i=1,nring_ind
    mult=1.
    rgord=COUNT(ring_ind(i,:).EQV..TRUE.)   ! get ring size                    
    IF (rgord==5)  mult=3.9   !  5 members
    IF (rgord==6)  mult=0.52  !  6 members
    IF (rgord==7)  mult=2.0   !  7 members
    IF (rgord==8)  mult=2.8   !  8 members
    IF (rgord==9)  mult=2.1   !  9 members
    IF (rgord==10) mult=0.24  ! 10 members
    IF (rgord==11) mult=12.   ! 11 members
    Fring=Fring*mult    
  ENDDO ringloop

END SUBROUTINE ring_fact_mono

!=======================================================================
! Purpose: Assign reference rate coefficients, see tables 4 in 
! Jenkin et al., 2020. 
!=======================================================================
SUBROUTINE ref_rate_conj(cdtable,cdsub,arrhc)
  USE keyparameter, ONLY: saru,mxcp
  USE mapping, ONLY: gettrack
  IMPLICIT NONE
  
  INTEGER,INTENT(IN) :: cdtable(:)         ! node # bearing a "Cd"
  INTEGER,INTENT(IN) :: cdsub(:)           ! # of C (including CO) bonded to Cd in cdtable(i)
  REAL,INTENT(OUT)   :: arrhc(:)

  INTEGER :: nb
  INTEGER :: ic1,ic2,ic3,ic4
  
  arrhc(1)=1E-14 ; arrhc(2:3)=0.

  ic1=cdtable(1) ; ic2=cdtable(2) ; ic3=cdtable(3) ; ic4=cdtable(4)
  nb=cdsub(ic1)+cdsub(ic2)+cdsub(ic3)+cdsub(ic4)

  SELECT CASE (nb)
    CASE (3) 
      IF ((cdsub(ic1)==1).OR.(cdsub(ic4)==1) ) THEN
         arrhc(3)=1677.
      ELSE IF ((cdsub(ic1)==0).AND.(cdsub(ic4)==0) ) THEN
         arrhc(3)=1980.
      ENDIF
    
    CASE (4) 
      IF ((cdsub(ic1)==0).AND.(cdsub(ic4)==0) .OR. &
          (cdsub(ic1)==0).AND.(cdsub(ic3)==0) .OR. &
          (cdsub(ic2)==0).AND.(cdsub(ic3)==0) .OR. &
          (cdsub(ic2)==0).AND.(cdsub(ic4)==0)) THEN
         arrhc(3)=1774.
      ENDIF
      IF (((cdsub(ic2)==2).AND.(cdsub(ic3)==1).AND.(cdsub(ic4)==1)) .OR. &
          ((cdsub(ic1)==1).AND.(cdsub(ic2)==1).AND.(cdsub(ic3)==2)) .OR. &
          ((cdsub(ic2)==1).AND.(cdsub(ic3)==2).AND.(cdsub(ic4)==1)) .OR. &
          ((cdsub(ic1)==1).AND.(cdsub(ic2)==2).AND.(cdsub(ic1)==1)) .OR. &
          ((cdsub(ic1)==0).AND.(cdsub(ic2)==1).AND.(cdsub(ic3)==1).AND.(cdsub(ic4)==2))) THEN
         arrhc(3)=1439.
      ENDIF
      IF ((cdsub(ic1)==1).AND.(cdsub(ic2)==1).AND.(cdsub(ic3)==1).AND.(cdsub(ic4)==1)) arrhc(3)=1008.

    CASE (5) 
      IF ((cdsub(ic1)==0).OR.(cdsub(ic4)==0))  arrhc(3)=1214.
      IF ((cdsub(ic1)/=0).AND.(cdsub(ic4)/=0)) arrhc(3)=686.

    CASE (6) 
      IF ((cdsub(ic1)==0).OR.(cdsub(ic4)==0))  arrhc(3)=887.
      IF ((cdsub(ic1)/=0).AND.(cdsub(ic4)/=0)) arrhc(3)=359.

    CASE (7) ; arrhc(3)=142.
    CASE (8) ; arrhc(3)=0.
  END SELECT 

END SUBROUTINE ref_rate_conj

!=======================================================================
! Compute ring factors to apply for conjugated structure for the 
! O3+alkene reaction rates as described by Jenkin et al., 2020 (see 
! factors in table 5).
!=======================================================================
SUBROUTINE ring_fact_conj(group,bond,cdtable,ngr,Fring)
  USE ringtool, ONLY: ring_data
  IMPLICIT NONE
  CHARACTER(LEN=*),INTENT(IN) :: group(:)             ! group matrix
  INTEGER,INTENT(IN)  :: bond(:,:)                    ! node matrix
  INTEGER,INTENT(IN)  :: cdtable(:)                   ! node # bearing a "Cd"
  INTEGER,INTENT(IN)  :: ngr                          ! # of groups
  REAL,INTENT(OUT)    :: Fring                        ! substituent factor to apply

  INTEGER,PARAMETER   :: mxirg=6                      ! max # of distinct rings 
  INTEGER             :: i,rgord
  REAL                :: mult
  INTEGER             :: nring_ind                    ! # of disctinct rings
  INTEGER             :: trackrg(mxirg,SIZE(group))   ! (a,:)== track (node #) belonging ring a
  LOGICAL             :: ring_ind(mxirg,SIZE(group))  ! (a,b)==true if node b belong to ring a

  CALL ring_data(cdtable(1),ngr,bond,group,nring_ind,ring_ind,trackrg) 

  Fring=1.0
  ringloop: DO i=1,nring_ind
    mult=1.
    rgord=COUNT(ring_ind(i,:).EQV..TRUE.)   ! get ring size                    
    IF (rgord==6) mult=4.50  ! 6 members
    IF (rgord==7) mult=0.44  ! 7 members
    IF (rgord==8) mult=0.06  ! 8 members
    Fring=Fring*mult
  ENDDO ringloop

END SUBROUTINE ring_fact_conj

END MODULE o3addtool
