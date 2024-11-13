MODULE no3addtool
  IMPLICIT NONE
  CONTAINS
!SUBROUTINE no3add_c1(chem,bond,group,ncd,conjug,cdtable,cdsub,&
!                     nr,flag,tarrhc,pchem,coprod,flag_del,pchem_del,&
!                     coprod_del,sc_del)
!SUBROUTINE no3add_c2(chem,bond,group, &
!                     ncd,conjug,cdtable,cdsub,cdcarbo, &
!                     nr,flag,tarrhc,pchem,coprod,flag_del, &
!                     pchem_del,coprod_del,sc_del)
!SUBROUTINE no3add_c3(chem,bond,group, &
!                     ncd,conjug,cdtable,cdsub,cdcarbo, &
!                     nr,flag,tarrhc,pchem,coprod,flag_del, &
!                     pchem_del,coprod_del,sc_del)
!SUBROUTINE no3add_c4(chem,bond,group,ncd,cdtable,conjug,cdsub,cdcarbo, &
!                     nr,flag,tarrhc,pchem,coprod,flag_del,&
!                     pchem_del,coprod_del,sc_del)
!SUBROUTINE no3add_c5(chem,bond,group,ncd,conjug,cdtable,cdsub,cdcarbo,&
!                     nr,flag,tarrhc,pchem,coprod,flag_del,pchem_del,&
!                     coprod_del,sc_del)
!SUBROUTINE no3add_c6(chem,bond,group,nr,flag,tarrhc,pchem,coprod)
!SUBROUTINE no3add_c7(chem,bond,group,ncd,conjug,cdtable,&
!                     cdeth,nr,flag,tarrhc,pchem,coprod)
!SUBROUTINE kaddno3(tbond,tgroup,i1,i2,arrhc)
!=======================================================================
!=======================================================================
! PURPOSE : compute the reaction rate for NO3 addition on >C=C< bond for
! "case 1", i.e. species without C=O group conjugated with the C=C.   
! Part I performs the reaction for simple >C=C< bond ; Part II for  
!   the conjugated >C=C-C=C< bonds.                           
!=======================================================================
!=======================================================================
SUBROUTINE no3add_c1(chem,bond,group,cdtable,cdsub,&
                     nr,flag,tarrhc,pchem,coprod,flag_del,pchem_del,&
                     coprod_del,sc_del,nrxref,rxref)
  USE references, ONLY:mxlcod
  USE reactool, ONLY: swap, rebond
  USE normchem, ONLY: stdchm
  USE radchktool, ONLY: radchk
  USE toolbox, ONLY: stoperr,addrx,addref
  IMPLICIT NONE

  CHARACTER(LEN=*),INTENT(IN) :: chem      ! chemical formula
  CHARACTER(LEN=*),INTENT(IN) :: group(:)  ! group matrix
  INTEGER,INTENT(IN) :: bond(:,:)          ! node matrix
  INTEGER,INTENT(IN) :: cdtable(:)         ! node # bearing a "Cd"
  INTEGER,INTENT(IN) :: cdsub(:)           ! # of C (including CO) bonded to Cd in cdtable(i)

  INTEGER,INTENT(INOUT) :: nr              ! # of reactions attached to chem (incremented here)   
  INTEGER,INTENT(INOUT) :: flag(:)         ! flag to activate reaction i
  REAL,INTENT(INOUT)    :: tarrhc(:,:)     ! arrhenius coefficient for reaction i
  CHARACTER(LEN=*),INTENT(INOUT) :: pchem(:)        ! "main" product of reaction channel i
  CHARACTER(LEN=*),INTENT(INOUT) :: coprod(:,:)     ! coproducts (short names) for reaction channel i
  INTEGER,INTENT(INOUT)          :: flag_del(:)     ! flag for a 2nd channel of reaction i (delocalisation)
  CHARACTER(LEN=*),INTENT(INOUT) :: pchem_del(:)    ! "main" product for the 2nd channel
  CHARACTER(LEN=*),INTENT(INOUT) :: coprod_del(:,:) ! coproducts for 2nd channel
  REAL,INTENT(INOUT)             :: sc_del(:,:)     ! ratio od the 2 channels (for delocalisation)
  INTEGER,INTENT(INOUT)          :: nrxref(:)       ! # of references added for each reaction
  CHARACTER(LEN=*),INTENT(INOUT) :: rxref(:,:)      ! references (short codes) for the reactions

  CHARACTER(LEN=LEN(group(1))) :: tgroup(SIZE(group)),pold,pnew,tempgr
  CHARACTER(LEN=LEN(chem)) :: tempkc
  INTEGER :: tbond(SIZE(bond,1),SIZE(bond,1))
  INTEGER :: tempnring
  INTEGER :: i,j1,nc,ic1,ic2
  REAL    :: fract(4)
  INTEGER :: ialpha, ibeta, igamma, idelta
  INTEGER :: nip
  REAL    :: arrhc(3) 
  CHARACTER(LEN=9),PARAMETER :: progname= 'no3add_c1'
  CHARACTER(LEN=70)          :: mesg

  INTEGER,PARAMETER :: mxrpd=2  ! max # of products returned by radchk sub
  CHARACTER(LEN=LEN(chem))  :: rdckpd(mxrpd)
  CHARACTER(LEN=LEN(coprod(1,1))) :: rdckcopd(mxrpd,SIZE(coprod,2))
  REAL :: sc(mxrpd)

  tgroup(:)=group(:)  ;  tbond(:,:)=bond(:,:)

! ------------------------------------------------------
! PART I : non conjugated Cds (simple >C=C<) structures
! ------------------------------------------------------
!  IF (conjug==0) THEN
  IF (cdtable(3)==0) THEN
  
! Loop over the 2 possible double bonds. If only 1 C=C bond, exit after 1st
    ic1=1 ; ic2=2
    CALL kaddno3(tbond,tgroup,cdtable(ic1),cdtable(ic2),arrhc)

! assign the branching ratio
    IF      (cdsub(ic1)==cdsub(ic2)) THEN ; fract(ic1)=0.5 ; fract(ic2)=0.5
    ELSE IF (cdsub(ic1)>cdsub(ic2))  THEN ; fract(ic2)=1. ; fract(ic1)=0.
    ELSE                                  ; fract(ic2)=0. ; fract(ic1)=1.
    ENDIF

! do the reaction
    DO i=ic1,ic2
      IF (fract(i)>0.) THEN
        CALL addrx(progname,chem,nr,flag)
        tarrhc(nr,1)=arrhc(1)*fract(i)
        tarrhc(nr,2:3)=arrhc(2:3)
        CALL addref(progname,'JK10KEV000',nrxref(nr),rxref(nr,:),chem)

! convert to single carbon bonds
        tbond(cdtable(ic1),cdtable(ic2))=1
        tbond(cdtable(ic2),cdtable(ic1))=1
	  
        pold='Cd' ;  pnew='C' 
        tempgr=group(cdtable(ic1))
        CALL swap(tempgr,pold,tgroup(cdtable(ic1)),pnew) ;
        tempgr=group(cdtable(ic2))
        CALL swap(tempgr,pold,tgroup(cdtable(ic2)),pnew)

! add (ONO2) to cdtable(i) carbon,
        nc=INDEX(tgroup(cdtable(i)),' ') ;
        tgroup(cdtable(i))(nc:nc+5)='(ONO2)'
             
! add radical dot to the other carbon:
        IF (i==ic1) THEN ; j1=ic2
        ELSE             ; j1=ic1
        ENDIF
        nc=INDEX(tgroup(cdtable(j1)),' ')
        tgroup(cdtable(j1))(nc:nc)='.'

! rebuild, check, and rename:
        CALL rebond(tbond,tgroup,tempkc,tempnring)
        CALL radchk(tempkc,rdckpd,rdckcopd,nip,sc,nrxref(nr),rxref(nr,:))
        pchem(nr)=rdckpd(1)
        sc_del(nr,1)=sc(1)
        IF (nip==2) THEN
          mesg="nip=2 unexpected "
          CALL stoperr(progname,mesg,chem)
        ENDIF
        CALL stdchm(pchem(nr))
        coprod(nr,:)=rdckcopd(1,:)

! reset groups and bonds:
        tgroup(:)=group(:)
        tbond(:,:)=bond(:,:)
      ENDIF
    ENDDO
  
! -------------------------------------------------
! PART II : conjugated Cds (>C=C-C=C<) structures
! -------------------------------------------------
!  ELSE IF (conjug==1) THEN
  ELSE IF (cdtable(3)/=0) THEN

! make reaction: NO3 adds at terminal position (ialpha positions), 
! radical dot is at beta position (ibeta) or delta position (idelta).
! Electron delocalisation considered in radchk 

    DO i=1,4,3    
      ialpha=i
      IF (ialpha==1) THEN ; ibeta=2 ; igamma=3 ; idelta=4
      ELSE                ; ibeta=3 ; igamma=2 ; idelta=1
      ENDIF

      CALL kaddno3(tbond,tgroup,cdtable(ialpha),cdtable(ibeta),arrhc)

      CALL addrx(progname,chem,nr,flag)
      tarrhc(nr,1)=arrhc(1)
      tarrhc(nr,2:3)=arrhc(2:3)
      CALL addref(progname,'JK10KEV000',nrxref(nr),rxref(nr,:),chem)

! convert to single carbon bonds:
      tbond(cdtable(ialpha),cdtable(ibeta))=1
      tbond(cdtable(ibeta),cdtable(ialpha))=1

      pold='Cd' ; pnew='C'      
      tempgr=group(cdtable(ialpha))
      CALL swap(tempgr,pold,tgroup(cdtable(ialpha)),pnew)
      tempgr=group(cdtable(ibeta))
      CALL swap(tempgr,pold,tgroup(cdtable(ibeta)),pnew)

! add (ONO2) to cdtable(i) carbon,
      nc=INDEX(tgroup(cdtable(ialpha)),' ')
      tgroup(cdtable(ialpha))(nc:nc+5)='(ONO2)'
             
! add radical dot to the other carbon:
      nc=INDEX(tgroup(cdtable(ibeta)),' ')
      tgroup(cdtable(ibeta))(nc:nc)='.'

! rebuild, check, and rename:
      CALL rebond(tbond,tgroup,tempkc,tempnring)
      CALL radchk(tempkc,rdckpd,rdckcopd,nip,sc,nrxref(nr),rxref(nr,:))
      pchem(nr)=rdckpd(1)
      sc_del(nr,1)=sc(1)
      IF (nip==2) THEN
        flag_del(nr)=1
        pchem_del(nr)=rdckpd(2)
        sc_del(nr,2)=sc(2)
        CALL stdchm(pchem_del(nr))
        coprod_del(nr,:)=rdckcopd(2,:)
      ENDIF
      CALL stdchm(pchem(nr))
      coprod(nr,:)=rdckcopd(1,:)

! reset groups and bonds:
      tgroup(:)=group(:)
      tbond(:,:)=bond(:,:)
    ENDDO

! unexpected value for conjug - stop  
  ELSE
    mesg="conjug not equal 0 or 1 "
    CALL stoperr(progname,mesg,chem)
  ENDIF

END SUBROUTINE no3add_c1

!=======================================================================
!=======================================================================
! PURPOSE : computes the reaction rate for NO3 addition on >C=C-C=O bond 
! (case 2). The -CO-C=C-C=C-CO- structure is not taken into account by 
! this routine (see case 3).                   
!=======================================================================
!=======================================================================
SUBROUTINE no3add_c2(chem,bond,group,cdtable,cdcarbo, &
                     nr,flag,tarrhc,pchem,coprod,flag_del, &
                     pchem_del,coprod_del,sc_del,nrxref,rxref)
  USE references, ONLY:mxlcod
  USE reactool, ONLY: swap, rebond
  USE normchem, ONLY: stdchm
  USE radchktool, ONLY: radchk
  USE toolbox, ONLY: stoperr,addrx,addref
  IMPLICIT NONE

  CHARACTER(LEN=*),INTENT(IN) :: chem      ! chemical formula
  INTEGER,INTENT(IN) :: bond(:,:)          ! bond matrix
  CHARACTER(LEN=*),INTENT(IN) :: group(:)  ! group matrix
  INTEGER,INTENT(IN) :: cdtable(:)         ! node bearing a "Cd", ordered by bonds
  INTEGER,INTENT(IN) :: cdcarbo(:,:)       ! node # of CO bonded to the Cd in cdtable(i)

  INTEGER,INTENT(INOUT) :: nr              ! # of reactions attached to chem (incremented here)   
  INTEGER,INTENT(INOUT) :: flag(:)         ! flag to activate reaction i
  REAL,INTENT(INOUT)    :: tarrhc(:,:)     ! arrhenius coefficient for reaction i
  CHARACTER(LEN=*),INTENT(INOUT) :: pchem(:)        ! "main" product of reaction channel i
  CHARACTER(LEN=*),INTENT(INOUT) :: coprod(:,:)     ! coproducts (short names) for reaction channel i
  INTEGER,INTENT(INOUT)          :: flag_del(:)     ! flag for a 2nd channel of reaction i (delocalisation)
  CHARACTER(LEN=*),INTENT(INOUT) :: pchem_del(:)    ! "main" product for the 2nd channel
  CHARACTER(LEN=*),INTENT(INOUT) :: coprod_del(:,:) ! coproducts for 2nd channel
  REAL,INTENT(INOUT)             :: sc_del(:,:)     ! ratio od the 2 channels (for delocalisation)
  INTEGER,INTENT(INOUT)          :: nrxref(:)       ! # of references added for each reaction
  CHARACTER(LEN=*),INTENT(INOUT) :: rxref(:,:)      ! references (short codes) for the reactions

  CHARACTER(LEN=LEN(group(1))) :: tgroup(SIZE(group)),pold,pnew
  CHARACTER(LEN=LEN(chem)) :: tempkc
  INTEGER :: tbond(SIZE(bond,1),SIZE(bond,1))
  INTEGER :: tempnring

  INTEGER  :: i,j,k,nc
  INTEGER  :: Ci,Cf
  REAL     :: kstruct
  REAL     :: rat(2)
  REAL     :: fac
  INTEGER  :: nbcarbi,nbcarbf
  INTEGER  :: posf, posi
  INTEGER  :: ratei,ratej
  INTEGER  :: nip
  REAL     :: arrhc(3)
  
  CHARACTER(LEN=9),PARAMETER :: progname='no3add_c2'

  INTEGER,PARAMETER :: mxrpd=2  ! max # of products returned by radchk sub
  CHARACTER(LEN=LEN(chem))  :: rdckpd(mxrpd)
  CHARACTER(LEN=LEN(coprod(1,1))) :: rdckcopd(mxrpd,SIZE(coprod,2))
  REAL :: sc(mxrpd)

  tgroup(:)=group(:) ; tbond(:,:)=bond(:,:)

! In cdtable, 1-2 and 3-4 indexes of cdtable are double bonded. Note
! that only one single bond must be "active" in this subroutine and 
! either cdtable(1-2) or cdtable(3-4) must be 0.
! find Ci and Cf such as Cdf=Cdi-C=O and count total number of 
! carbonyls bonded to the double bond
  iloop: DO i=1,3,2
    nbcarbi=0 ; nbcarbf=0
    IF (cdtable(i)/=0) THEN
      IF ( (cdcarbo(i,1)/=0).AND.(cdcarbo(i+1,1)==0) ) THEN
        Ci=cdtable(i) ; Cf=cdtable(i+1)
        posi=i  ; posf=i+1
        DO j=1,2
          IF (cdcarbo(posi,j)/=0) nbcarbi=nbcarbi+1
          IF (cdcarbo(posf,j)/=0) nbcarbf=nbcarbf+1
        ENDDO
        EXIT iloop
      ELSE IF ((cdcarbo(i+1,1)/=0).AND.(cdcarbo(i,1)==0))THEN
        Ci=cdtable(i+1) ; Cf=cdtable(i)
        posi=i+1  ;  posf=i
        DO j=1,2
          IF (cdcarbo(posi,j)/=0) nbcarbi=nbcarbi+1
          IF (cdcarbo(posf,j)/=0) nbcarbf=nbcarbf+1
        ENDDO
        EXIT iloop
! patch added to remove the problem of reaction depending on whether
! chem is written up side down or the opposite. More work needed to
! clean that
      ELSE IF ((cdcarbo(i,1)/=0).AND.(cdcarbo(i+1,1)/=0)) THEN
        IF      (group(cdtable(i))(1:4)=='CdH2') THEN ;  ratei=0
        ELSE IF (group(cdtable(i))(1:3)=='CdH')  THEN ;  ratei=1
        ELSE                                          ;  ratei=2
        ENDIF

        IF      (group(cdtable(i+1))(1:4)=='CdH2') THEN  ; ratej=0
        ELSE IF (group(cdtable(i+1))(1:3)=='CdH')  THEN  ; ratej=1
        ELSE                                             ; ratej=2
        ENDIF

        IF (ratei>=ratej) THEN
          Ci=cdtable(i)  ;  Cf=cdtable(i+1)
          posi=i  ;  posf=i+1
          DO j=1,2
            IF (cdcarbo(posi,j)/=0) nbcarbi=nbcarbi+1
            IF (cdcarbo(posf,j)/=0) nbcarbf=nbcarbf+1
          ENDDO
          EXIT iloop
        ELSE
          Ci=cdtable(i+1) ; Cf=cdtable(i)
          posi=i+1  ;  posf=i
          DO j=1,2
            IF (cdcarbo(posi,j)/=0) nbcarbi=nbcarbi+1
            IF (cdcarbo(posf,j)/=0) nbcarbf=nbcarbf+1
          ENDDO
          EXIT iloop
        ENDIF            
        
      ENDIF
    ENDIF
  ENDDO iloop

! double bond with carbonyl : k=kstruct*fac               
! RCH=CHCHO kstruct =2.11E-15  !  CH2=CRCHO kstruct =1.60E-15 
! RCH=CRCHO kstruct =7.60E-14  !  RRC=CHCHO kstruct =1.07E-13 
! RRC=CRCHO kstruct =2.30E-13  !  
! R=CHO or COR   fac=0.01      !  R=-C=C-C=O     fac=0.073  

! set rate and branching ratio (rat(1) is NO3 addition at Ci ; rat(2) is 
! NO3 addition at Cf)

! CH2=CR-CO- or CH2=CH-CO- structures. For acrolein, reaction is 
! assumed to proceed only by H-atom abstraction from -CHO 
  IF (group(Cf)(1:4)=='CdH2') THEN
    IF (group(Ci)(1:3)=='CdH') THEN
      kstruct=1.2E-16 ; rat(1)=0. ; rat(2)=1.
    ELSE
      kstruct=1.6E-15 ; rat(1)=0. ; rat(2)=1.
    ENDIF

! -CH=CR-CO- or -CH=CH-CO- structures        
  ELSE IF (group(Cf)(1:3)=='CdH') THEN
    IF (group(Ci)(1:3)=='CdH') THEN 
      kstruct=2.11E-15 ; rat(1)=0.5 ; rat(2)=0.5
    ELSE
      kstruct=7.6E-14  ; rat(1)=0.  ; rat(2)=1.
    ENDIF

! >C=CR-CO- or >C=CH-CO- structures
  ELSE IF (group(Cf)(1:2)=='Cd') THEN
    IF (group(Ci)(1:3)=='CdH') THEN 
      kstruct=1.07E-13 ; rat(1)=0.  ; rat(2)=1.
    ELSE 
      kstruct=2.3E-13  ; rat(1)=0.5 ; rat(2)=0.5
    ENDIF
  ENDIF

! correct kstruct as a function of the functional group on Ci and Cf
  fac=1.
  DO k=1,nbcarbf ; fac=fac*0.01 ; ENDDO
  IF (nbcarbi==2) fac=fac*0.01

! estimate rate constant
  CALL kaddno3(tbond,tgroup,Ci,Cf,arrhc)

! make the two NO3 additions to the >C=C<. Ci&Cf are swapped between j=1 and j=2.
  DO j=1,2
    CALL addrx(progname,chem,nr,flag)

! assign rate constant for ONO2 addition in Ci position
    tarrhc(nr,1)=arrhc(1)*rat(j)
    CALL addref(progname,'JK10KEV000',nrxref(nr),rxref(nr,:),chem)

! convert Cf and Ci to single bond carbon
    tbond(Cf,Ci)=1 ; tbond(Ci,Cf)=1
    pold='Cd' ; pnew='C' ;  CALL swap(group(Cf),pold,tgroup(Cf),pnew)
    pold='Cd' ; pnew='C' ;  CALL swap(group(Ci),pold,tgroup(Ci),pnew)

! add NO3 to Ci carbon and add radical dot to Cf
    nc=INDEX(tgroup(Ci),' ') ; tgroup(Ci)(nc:nc+5)='(ONO2)'
    nc=INDEX(tgroup(Cf),' ') ; tgroup(Cf)(nc:nc)='.'

! rebuild, check, and find coproducts:
    CALL rebond(tbond,tgroup,tempkc,tempnring)
    CALL radchk(tempkc,rdckpd,rdckcopd,nip,sc,nrxref(nr),rxref(nr,:))
    pchem(nr)=rdckpd(1)
    coprod(nr,:)=rdckcopd(1,:)
    CALL stdchm(pchem(nr))
    IF (nip/=1) THEN
      flag_del(nr)=1
      sc_del(nr,1)=sc(1) ; sc_del(nr,2)=sc(2)
      pchem_del(nr)=rdckpd(2)
      CALL stdchm(pchem_del(nr))
      coprod_del(nr,:)=rdckcopd(2,:)
    ENDIF

! reset:
    tgroup(:)=group(:)  ;  tbond(:,:)=bond(:,:)
                 
! convert Ci to Cf and Cf to Ci to form the other radical          
    k=Ci ; Ci=Cf ; Cf=k
  ENDDO
END SUBROUTINE no3add_c2


!=======================================================================
! PURPOSE: compute the reaction rate for NO3 addition on -CO-C=C-C=C-C=O 
! bond (case 3) only. 
!=======================================================================
SUBROUTINE no3add_c3(chem,bond,group,cdtable,cdcarbo, &
                     nr,flag,tarrhc,pchem,coprod,flag_del, &
                     pchem_del,coprod_del,sc_del,nrxref,rxref)
  USE references, ONLY:mxlcod
  USE reactool, ONLY: swap, rebond
  USE normchem, ONLY: stdchm
  USE radchktool, ONLY: radchk
  USE toolbox, ONLY: stoperr,addrx,addref
  IMPLICIT NONE

  CHARACTER(LEN=*),INTENT(IN) :: chem      ! chemical formula
  CHARACTER(LEN=*),INTENT(IN) :: group(:)  ! group matrix
  INTEGER,INTENT(IN) :: bond(:,:)          ! node matrix
  INTEGER,INTENT(IN) :: cdtable(:)         ! node # bearing a "Cd"
  INTEGER,INTENT(IN) :: cdcarbo(:,:)       ! node # of CO bonded to the Cd in cdtable(i)

  INTEGER,INTENT(INOUT) :: nr         ! # of reactions attached to chem (incremented here)
  INTEGER,INTENT(INOUT) :: flag(:)    ! flag to activate reaction i
  REAL,INTENT(INOUT) :: tarrhc(:,:)   ! arrhenius coefficient for reaction i
  CHARACTER(LEN=*),INTENT(INOUT) :: pchem(:)    ! "main" product of reaction channel i
  CHARACTER(LEN=*),INTENT(INOUT) :: coprod(:,:) ! coproducts (short names) for reaction channel i
  INTEGER,INTENT(INOUT)          :: flag_del(:)     ! flag for a 2nd channel of reaction i (delocalisation)
  CHARACTER(LEN=*),INTENT(INOUT) :: pchem_del(:)    ! "main" product for the 2nd channel
  CHARACTER(LEN=*),INTENT(INOUT) :: coprod_del(:,:) ! coproducts for 2nd channel
  REAL,INTENT(INOUT)             :: sc_del(:,:)     ! ratio od the 2 channels (for delocalisation)
  INTEGER,INTENT(INOUT)          :: nrxref(:)       ! # of references added for each reaction
  CHARACTER(LEN=*),INTENT(INOUT) :: rxref(:,:)      ! references (short codes) for the reactions

  CHARACTER(LEN=LEN(group(1))) :: tgroup(SIZE(group)), pold, pnew
  CHARACTER(LEN=LEN(chem)) :: tempkc
  INTEGER :: tbond(SIZE(bond,1),SIZE(bond,1))

  REAL    :: arrhc(3)
  INTEGER :: i,j,k,nc
  INTEGER :: Ci,Cf
  INTEGER :: tempnring
  REAL    :: kstruct
  REAL    :: rat(2)
  REAL    :: fac
  INTEGER :: nbcarbi,nbcarbf
  INTEGER :: posf, posi
  INTEGER :: nip

  INTEGER,PARAMETER :: mxrpd=2  ! max # of products returned by radchk sub
  CHARACTER(LEN=LEN(chem))  :: rdckpd(mxrpd)
  CHARACTER(LEN=LEN(coprod(1,1))) :: rdckcopd(mxrpd,SIZE(coprod,2))
  REAL :: sc(mxrpd)

  CHARACTER(LEN=9),PARAMETER :: progname='no3add_c3'
  CHARACTER(LEN=70)          :: mesg
 
! -----------
! Initialize
! -----------
  tgroup(:)=group(:) ; tbond(:,:)=bond(:,:)

! In cdtable, 1-2 and 3-4 indexes of cdtable are double bonded.
! find Ci and Cf such that Cdf=Cdi-C=O and count total number of 
! carbonyls bonded to the double bond. Each double bond is treated as
! if it is not conjugated. A weighting factor is used to take into
! account the effect of the reactivity of the "second" -C=C-CO- 
! on the "first" -C=C-CO- bond
  DO i=1,4,3
    nbcarbi=0 ; nbcarbf=0
    IF (i==1) THEN
      IF (cdcarbo(1,1)/=0) THEN        
        Ci=cdtable(1) ; Cf=cdtable(2) ;  posi=1 ; posf=2
      ELSE IF (cdcarbo(2,1)/=0) THEN
        Ci=cdtable(2) ; Cf=cdtable(1) ;  posi=2 ; posf=1
      ENDIF
    ELSE IF (i==4) THEN
      IF (cdcarbo(4,1)/=0) THEN        
        Ci=cdtable(4) ; Cf=cdtable(3) ; posi=4 ; posf=3
      ELSE IF (cdcarbo(3,1)/=0) THEN
        Ci=cdtable(3) ; Cf=cdtable(4) ; posi=3 ; posf=4
      ENDIF
    ENDIF
    DO j=1,2
      IF (cdcarbo(posi,j)/=0) nbcarbi=nbcarbi+1
      IF (cdcarbo(posf,j)/=0) nbcarbf=nbcarbf+1
    ENDDO

! set the rate constant

! -CH=CR-CO- or -CH=CH-CO- structures
    IF (group(Cf)(1:3)=='CdH') THEN
      IF (group(Ci)(1:3)=='CdH') THEN 
        kstruct=2.11E-15 ; rat(1)=0.5 ; rat(2)=0.5
      ELSE
        kstruct= 7.6E-14   ; rat(1)=0.  ; rat(2)=1.
      ENDIF

! >C=CR-CO- or >C=CH-CO- structures
    ELSE IF (group(Cf)(1:2)=='Cd') THEN
      IF (group(Ci)(1:3)=='CdH') THEN 
        kstruct=1.07E-13 ; rat(1)=0.  ; rat(2)=1.
      ELSE 
        kstruct=2.3E-13  ; rat(1)=0.5 ; rat(2)=0.5
      ENDIF

    ELSE
      mesg=" -CO-C=C-C=C-CO- expected "
      CALL stoperr(progname,mesg,chem)
    ENDIF

! correct kstruct as a function of the functional group on Ci and Cf
! and multiply by 0.073 due to -C=C-CO-
    fac=0.073
    DO k=1,nbcarbf ;  fac=fac*0.01 ; ENDDO
    IF (nbcarbi==2) fac=fac*0.01

! estimate rate constant
    CALL kaddno3(tbond,tgroup,Ci,Cf,arrhc)

! make the two NO3 additions to the >C=C<. Ci&Cf are swapped between j=1 and j=2.
    DO j=1,2
      CALL addrx(progname,chem,nr,flag)

! assign rate constant for ONO2 addition in Ci position
      tarrhc(nr,1)=arrhc(1)*rat(j)
      CALL addref(progname,'JK10KEV000',nrxref(nr),rxref(nr,:),chem)

! convert Cf to single bond carbon
      tbond(Cf,Ci)=1 ; tbond(Ci,Cf)=1
      pold='Cd' ; pnew='C' ;
      CALL swap(group(Cf),pold,tgroup(Cf),pnew)
      pold='Cd' ;  pnew='C'         
      CALL swap(group(Ci),pold,tgroup(Ci),pnew)

! add NO3 to Ci carbon and add radical dot to Cf
      nc=INDEX(tgroup(Ci),' ') ; tgroup(Ci)(nc:nc+5)='(ONO2)'
      nc=INDEX(tgroup(Cf),' ') ; tgroup(Cf)(nc:nc)='.'

! rebuild, check, and find coproducts:
      CALL rebond(tbond,tgroup,tempkc,tempnring)
      CALL radchk(tempkc,rdckpd,rdckcopd,nip,sc,nrxref(nr),rxref(nr,:))
      pchem(nr)=rdckpd(1)
      coprod(nr,:)=rdckcopd(1,:)
      CALL stdchm(pchem(nr))
      IF (nip/=1) THEN
        flag_del(nr)=1
        sc_del(nr,1)=sc(1)
        sc_del(nr,2)=sc(2)
        pchem_del(nr)=rdckpd(2)
        CALL stdchm(pchem_del(nr))
        coprod_del(nr,:)=rdckcopd(2,:)
      ENDIF

! reset:
      tgroup(:)=group(:) ;  tbond(:,:)=bond(:,:)
                 
! Convert Ci to Cf and Cf to Ci to form the other radical          
      k=Ci ; Ci=Cf ; Cf=k
    ENDDO
  ENDDO
END SUBROUTINE no3add_c3

!=======================================================================
!=======================================================================
SUBROUTINE no3add_c4(chem,bond,group,cdtable,cdsub,cdcarbo, &
                     nr,flag,tarrhc,pchem,coprod,flag_del,&
                     pchem_del,coprod_del,sc_del,nrxref,rxref)
  USE toolbox, ONLY: stoperr
  IMPLICIT NONE

  CHARACTER(LEN=*),INTENT(IN) :: chem      ! chemical formula
  INTEGER,INTENT(IN) :: bond(:,:)          ! node matrix
  CHARACTER(LEN=*),INTENT(IN) :: group(:)  ! group matrix
  INTEGER,INTENT(IN) :: cdtable(:)         ! node # bearing a "Cd"
  INTEGER,INTENT(IN) :: cdsub(:)           ! # of C (including CO) bonded to Cd in cdtable(i)
  INTEGER,INTENT(IN) :: cdcarbo(:,:)       ! node # of CO bonded to the Cd in cdtable(i)

  INTEGER,INTENT(INOUT) :: nr              ! # of reactions attached to chem (incremented here)   
  INTEGER,INTENT(INOUT) :: flag(:)         ! flag to activate reaction i
  REAL,INTENT(INOUT)    :: tarrhc(:,:)     ! arrhenius coefficient for reaction i
  CHARACTER(LEN=*),INTENT(INOUT) :: pchem(:)        ! "main" product of reaction channel i
  CHARACTER(LEN=*),INTENT(INOUT) :: coprod(:,:)     ! coproducts (short names) for reaction channel i
  INTEGER,INTENT(INOUT)          :: flag_del(:)     ! flag for a 2nd channel of reaction i (delocalisation)
  CHARACTER(LEN=*),INTENT(INOUT) :: pchem_del(:)    ! "main" product for the 2nd channel
  CHARACTER(LEN=*),INTENT(INOUT) :: coprod_del(:,:) ! coproducts for 2nd channel
  REAL,INTENT(INOUT)             :: sc_del(:,:)     ! ratio od the 2 channels (for delocalisation)
  INTEGER,INTENT(INOUT)          :: nrxref(:)       ! # of references added for each reaction
  CHARACTER(LEN=*),INTENT(INOUT) :: rxref(:,:)      ! references (short codes) for the reactions

  INTEGER :: i
  INTEGER :: cpcdtable(SIZE(cdtable))
  CHARACTER(LEN=9),PARAMETER :: progname='no3add_c5'
  CHARACTER(LEN=70)          :: mesg
  
! perform reaction for each C=C, ignoring the other C=C bond 
  DO i=1,3,2

! check everything is ok
    IF (cdtable(i)==0) THEN
      mesg="Cd not found. "
      CALL stoperr(progname,mesg,chem)
    ENDIF

! make a copy of cdtable and "kill" the Cd bond not considered
    cpcdtable(:)=0
    cpcdtable(i)=cdtable(i) ; cpcdtable(i+1)=cdtable(i+1)

! send to reaction subroutine
    IF (cdcarbo(i,1)/=0) THEN
      CALL no3add_c2(chem,bond,group,cpcdtable,cdcarbo, &
                     nr,flag,tarrhc,pchem,coprod,flag_del, &
                     pchem_del,coprod_del,sc_del,nrxref,rxref) 
    ELSE IF (cdcarbo(i+1,1)/=0) THEN
      CALL no3add_c2(chem,bond,group,cpcdtable,cdcarbo, &
                     nr,flag,tarrhc,pchem,coprod,flag_del, &
                     pchem_del,coprod_del,sc_del,nrxref,rxref)
    ELSE
      CALL no3add_c1(chem,bond,group,cpcdtable,cdsub, &
                     nr,flag,tarrhc,pchem,coprod,flag_del,   &
                     pchem_del,coprod_del,sc_del,nrxref,rxref)
    ENDIF
  ENDDO
END SUBROUTINE no3add_c4

!=======================================================================
!=======================================================================
! PURPOSE: computes the reaction rate for OH addition to >C=C-C=C-C=O 
! bond (case 5) only.
!=======================================================================
!=======================================================================
SUBROUTINE no3add_c5(chem,bond,group,cdtable,cdsub,cdcarbo,&
                     nr,flag,tarrhc,pchem,coprod,flag_del,pchem_del,&
                     coprod_del,sc_del,nrxref,rxref)
  USE references, ONLY:mxlcod
  USE reactool, ONLY: swap,rebond
  USE normchem, ONLY:stdchm
  USE radchktool, ONLY: radchk
  USE toolbox, ONLY: stoperr,addrx,addref
  IMPLICIT NONE

  CHARACTER(LEN=*),INTENT(IN) :: chem      ! chemical formula
  CHARACTER(LEN=*),INTENT(IN) :: group(:)  ! group matrix
  INTEGER,INTENT(IN) :: bond(:,:)          ! node matrix
  INTEGER,INTENT(IN) :: cdtable(:)         ! node # bearing a "Cd"
  INTEGER,INTENT(IN) :: cdsub(:)           ! # of C (including CO) bonded to Cd in cdtable(i)
  INTEGER,INTENT(IN) :: cdcarbo(:,:)       ! node # of CO bonded to the Cd in cdtable(i)

  INTEGER,INTENT(INOUT) :: nr         ! # of reactions attached to chem (incremented here)
  INTEGER,INTENT(INOUT) :: flag(:)    ! flag to activate reaction i
  REAL,INTENT(INOUT) :: tarrhc(:,:)   ! arrhenius coefficient for reaction i
  CHARACTER(LEN=*),INTENT(INOUT) :: pchem(:)    ! "main" product of reaction channel i
  CHARACTER(LEN=*),INTENT(INOUT) :: coprod(:,:) ! coproducts (short names) for reaction channel i
  INTEGER,INTENT(INOUT)          :: flag_del(:)     ! flag for a 2nd channel of reaction i (delocalisation)
  CHARACTER(LEN=*),INTENT(INOUT) :: pchem_del(:)    ! "main" product for the 2nd channel
  CHARACTER(LEN=*),INTENT(INOUT) :: coprod_del(:,:) ! coproducts for 2nd channel
  REAL,INTENT(INOUT)             :: sc_del(:,:)     ! ratio od the 2 channels (for delocalisation)
  INTEGER,INTENT(INOUT)          :: nrxref(:)       ! # of references added for each reaction
  CHARACTER(LEN=*),INTENT(INOUT) :: rxref(:,:)      ! references (short codes) for the reactions

  CHARACTER(LEN=LEN(group(1))) :: tgroup(SIZE(group)),pold,pnew,tempgr
  CHARACTER(LEN=LEN(chem)) :: tempkc
  INTEGER :: tbond(SIZE(bond,1),SIZE(bond,1))
  INTEGER :: nc,nb
  INTEGER :: tempnring
  REAL    :: fract(4)
  REAL    :: rtot
  INTEGER :: ialpha, ibeta
  INTEGER :: nip
  REAL    :: arrhc(3)

  INTEGER,PARAMETER :: mxrpd=2  ! max # of products returned by radchk sub
  CHARACTER(LEN=LEN(chem))  :: rdckpd(mxrpd)
  CHARACTER(LEN=LEN(coprod(1,1))) :: rdckcopd(mxrpd,SIZE(coprod,2))
  REAL :: sc(mxrpd)

  CHARACTER(LEN=9),PARAMETER :: progname='no3add_c5'
  CHARACTER(LEN=70)          :: mesg
 
  tgroup(:)=group(:) ; tbond(:,:)=bond(:,:)
  ialpha=0 ; 

! =================================================================
!  conjugated Cds (>C=C-C=C<) structures (note from no3add_c1)
! =================================================================
!
! alkyl group at terminal position (i.e. alpha and delta) seems to have
! a greater impact on rate constant than alkyl group at internal
! position (i.e beta and gamma). Therefore, rate constant is asigned
! by considering first the alkyl group at terminal position. A correction
! is then made to account for the contribution of the internal alkyl group.
! Rate constants are assigned as follows :
! CH2=CH-CH=CH2   : 1.00E-13
! C-CH=CH-CH=CH2  : 1.50E-12 (from CH3-CH=CH-CH=CH2)
! C-CH=CH-CH=CH-C : 1.60E-11 (from CH3-CH=CH-CH=CH-CH3)
! and then multiply by 6.8 if there is 1 internal alkyl group (ratio 
! isoprene / 1,3 butadiene) and 21 if there are 2 internal alkyl groups
! (ratio  2,3 dimethyl / 1,3 butadiene and butadiene)
! ALL THIS IS OF COURSE A VERY CRUDE ESTIMATE. (more work should be done in
! the future to better estimate these rate constants)
!
! For branching ratio : NO3 is assumed to add at the less substitued terminal
! position. If the 2 terminal positions are equivalent, then the dominant
! position is assumed to be the one yielding the highest substitued radical

! --------------------------------------
! set rate constant and branching ratio
! --------------------------------------

! set the rate constant - contribution of external alkyl group
  nb=cdsub(1)+cdsub(4)
  IF      (nb==0) THEN ; rtot=1.00E-13
  ELSE IF (nb==1) THEN ; rtot=1.50E-12
  ELSE IF (nb>1)  THEN ; rtot=1.60E-11
  ENDIF
  rtot=rtot*0.5

! contribution of internal alkyl group (nb=3 means 1 internal alkyl group)
  nb=cdsub(2)+cdsub(3)
  IF (nb==3) rtot=rtot*6.8
  IF (nb >3) rtot=rtot*21.

! ----------------
! define ialpha position (furthest from carbonyl group)
! ----------------
  IF      ((cdcarbo(1,1)/=0) .OR. (cdcarbo(2,1)/=0)) THEN ; ialpha=4
  ELSE IF ((cdcarbo(4,1)/=0) .OR. (cdcarbo(3,1)/=0)) THEN ; ialpha=1
  ELSE
    mesg="no carbonyl found "
    CALL stoperr(progname,mesg,chem)
  ENDIF

! define branching ratio for NO3 addition as 100% at ialpha
  fract(:)=0.  ;  fract(ialpha)=1.

! ----------------
! do the reaction 
! ----------------
! NO3 adds at terminal position (ialpha position), radical dot is at
! beta position (ibeta) or delta position (idelta). 
! !! reduce cdsub(ibeta) by 1 to avoid considering conjugation bond !!
  IF (ialpha==1) THEN ; ibeta=2
  ELSE                ; ibeta=3
  ENDIF

! do the beta addition only, delocalisation is done in radchk
! estimate rate constant
  CALL kaddno3(tbond,tgroup,cdtable(ialpha),cdtable(ibeta),arrhc)
  CALL addrx(progname,chem,nr,flag)
  tarrhc(nr,1)=arrhc(1)
  CALL addref(progname,'JK10KEV000',nrxref(nr),rxref(nr,:),chem)

! convert to single carbon bonds:
  tbond(cdtable(ialpha),cdtable(ibeta))=1
  tbond(cdtable(ibeta),cdtable(ialpha))=1
      
  pold='Cd' ; pnew='C'      
  tempgr=group(cdtable(ialpha))
  CALL swap(tempgr,pold,tgroup(cdtable(ialpha)),pnew)
  tempgr=group(cdtable(ibeta))
  CALL swap(tempgr,pold,tgroup(cdtable(ibeta)),pnew)

! add (ONO2) to cdtable(i) carbon,
  nc=INDEX(tgroup(cdtable(ialpha)),' ')
  tgroup(cdtable(ialpha))(nc:nc+5)='(ONO2)'
             
! add radical dot to the other carbon:
  nc=INDEX(tgroup(cdtable(ibeta)),' ')
  tgroup(cdtable(ibeta))(nc:nc)='.'

! rebuild, check, and rename:
  CALL rebond(tbond,tgroup,tempkc,tempnring)
  CALL radchk(tempkc,rdckpd,rdckcopd,nip,sc,nrxref(nr),rxref(nr,:))
  pchem(nr)=rdckpd(1)
  sc_del(nr,1)=sc(1)
  IF (nip==2) THEN
    flag_del(nr)=1
    pchem_del(nr)=rdckpd(2)
    sc_del(nr,2)=sc(2)
    CALL stdchm(pchem_del(nr))
    coprod_del(nr,:)=rdckcopd(2,:)
  ENDIF
  CALL stdchm(pchem(nr))
  coprod(nr,:)=rdckcopd(1,:)

! reset groups and bonds:
  tgroup(:)=group(:) ; tbond(:,:)=bond(:,:)
END SUBROUTINE no3add_c5


!=======================================================================
!=======================================================================
! PURPOSE: computes the reaction rate for NO3 addition to >C=C=O bond 
! (ketene) (case 6) only. 
!=======================================================================
!=======================================================================
SUBROUTINE no3add_c6(chem,bond,group,nr,flag,tarrhc,pchem,coprod,nrxref,rxref)
  USE references, ONLY:mxlcod
  USE reactool
  USE normchem
  USE radchktool, ONLY: radchk
  USE toolbox, ONLY: stoperr,addrx,addref
  IMPLICIT NONE

  CHARACTER(LEN=*),INTENT(IN) :: chem      ! chemical formula
  CHARACTER(LEN=*),INTENT(IN) :: group(:)  ! group matrix
  INTEGER,INTENT(IN) :: bond(:,:)          ! node matrix

  INTEGER,INTENT(INOUT) :: nr         ! # of reactions attached to chem (incremented here)
  INTEGER,INTENT(INOUT) :: flag(:)    ! flag to activate reaction i
  REAL,INTENT(INOUT) :: tarrhc(:,:)   ! arrhenius coefficient for reaction i
  CHARACTER(LEN=*),INTENT(INOUT) :: pchem(:)    ! "main" product of reaction channel i
  CHARACTER(LEN=*),INTENT(INOUT) :: coprod(:,:) ! coproducts (short names) for reaction channel i
  INTEGER,INTENT(INOUT)          :: nrxref(:)       ! # of references added for each reaction
  CHARACTER(LEN=*),INTENT(INOUT) :: rxref(:,:)      ! references (short codes) for the reactions

  CHARACTER(LEN=LEN(group(1))) :: tgroup(SIZE(group)),pold,pnew
  CHARACTER(LEN=LEN(chem)) :: tempkc
  INTEGER :: tbond(SIZE(bond,1),SIZE(bond,1))
  INTEGER :: i,j,k,nca,tempnring
  INTEGER :: nip

  INTEGER,PARAMETER :: mxrpd=2  ! max # of products returned by radchk sub
  CHARACTER(LEN=LEN(chem))  :: rdckpd(mxrpd)
  CHARACTER(LEN=LEN(coprod(1,1))) :: rdckcopd(mxrpd,SIZE(coprod,2))
  REAL :: sc(mxrpd)
  
  CHARACTER(LEN=9),PARAMETER :: progname='no3add_c6'
  CHARACTER(LEN=70)          :: mesg

  tgroup(:)=group(:) ; tbond(:,:)=bond(:,:)
  nca=COUNT(tgroup/=' ')

! add NO3 only at carbonyl group
  DO i=1,nca
    IF (INDEX(group(i),'CdO')==0) CYCLE
    CALL addrx(progname,chem,nr,flag)
   
    pold='CdO'  ;  pnew='CO.'
    CALL swap(group(i),pold,tgroup(i),pnew)
    DO j=1,nca
      IF (bond(i,j)==2) THEN
        tbond(i,j)=1 ; tbond(j,i)=1
        IF (INDEX(group(j),'Cd')/=0) THEN
          k=INDEX(group(j),' ')
          tgroup(j)(2:k-1)=group(j)(3:k)
        ELSE
          mesg="problem encountered for: "
          CALL stoperr(progname,mesg,chem)
        ENDIF
        k=INDEX(tgroup(j),' ') ; tgroup(j)(k:k+5)='(ONO2)'
      ENDIF
    ENDDO

! rebuild, check, and find co-products:
    CALL rebond(tbond,tgroup,tempkc,tempnring)
    CALL radchk(tempkc,rdckpd,rdckcopd,nip,sc,nrxref(nr),rxref(nr,:))
    pchem(nr)=rdckpd(1)
    coprod(nr,:)=rdckcopd(1,:)
    CALL stdchm(pchem(nr))
    IF (nip/=1) THEN
      mesg="unexpected 2 products "
      CALL stoperr(progname,mesg,chem)
    ENDIF

! reset groups:
    tgroup(:)=group(:)

! set rate constants based on CdH2=CdO
    tarrhc(nr,1)=1.1E-13   ! Seland et al., 1996
    CALL addref(progname,'JS96KMV000',nrxref(nr),rxref(nr,:),chem)
  ENDDO
END SUBROUTINE no3add_c6

!=======================================================================
!=======================================================================
! PURPOSE: compute the reaction rate for NO3 addition on O-C=C bond 
! (case 7). Species having a C=O group conjugated with the C=C are not 
! taken into account here.
!=======================================================================
!=======================================================================
SUBROUTINE no3add_c7(chem,bond,group,cdtable,cdeth,nr,flag,tarrhc,&
                     pchem,coprod,flag_del,pchem_del,&
                     coprod_del,sc_del,nrxref,rxref)
  USE references, ONLY:mxlcod
  USE reactool, ONLY:swap, rebond
  USE ringtool, ONLY: findring
  USE normchem, ONLY: stdchm
  USE radchktool, ONLY: radchk
  USE toolbox, ONLY: stoperr,addrx,addref
  IMPLICIT NONE

  CHARACTER(LEN=*),INTENT(IN) :: chem      ! chemical formula
  INTEGER,INTENT(IN) :: bond(:,:)          ! bond matrix
  CHARACTER(LEN=*),INTENT(IN) :: group(:)  ! group matrix
  INTEGER,INTENT(IN) :: cdtable(:)         ! node bearing a "Cd", ordered by bonds
  INTEGER,INTENT(IN) :: cdeth(:,:)         ! node # of -O- bonded to the Cd in cdtable(i)

  INTEGER,INTENT(INOUT) :: nr         ! # of reactions attached to chem (incremented here)
  INTEGER,INTENT(INOUT) :: flag(:)    ! flag to activate reaction i
  REAL,INTENT(INOUT) :: tarrhc(:,:)   ! arrhenius coefficient for reaction i
  CHARACTER(LEN=*),INTENT(INOUT) :: pchem(:)    ! "main" product of reaction channel i
  CHARACTER(LEN=*),INTENT(INOUT) :: coprod(:,:) ! coproducts (short names) for reaction channel i
  INTEGER,INTENT(INOUT)          :: flag_del(:)     ! flag for a 2nd channel of reaction i (delocalisation)
  CHARACTER(LEN=*),INTENT(INOUT) :: pchem_del(:)    ! "main" product for the 2nd channel
  CHARACTER(LEN=*),INTENT(INOUT) :: coprod_del(:,:) ! coproducts for 2nd channel
  REAL,INTENT(INOUT)             :: sc_del(:,:)     ! ratio od the 2 channels (for delocalisation)
  INTEGER,INTENT(INOUT)          :: nrxref(:)       ! # of references added for each reaction
  CHARACTER(LEN=*),INTENT(INOUT) :: rxref(:,:)      ! references (short codes) for the reactions

  CHARACTER(LEN=LEN(group(1))) :: tgroup(SIZE(group))
  CHARACTER(LEN=LEN(group(1))) :: pold, pnew,tempkg
  CHARACTER(LEN=LEN(chem))     :: tempkc
  INTEGER :: tbond(SIZE(bond,1),SIZE(bond,2))

  INTEGER :: i,j,nc,ic1,ic2,Cf,Ci,p
  INTEGER :: tempnring,ic
  REAL    :: fract(4), rtot
  INTEGER :: nip
  REAL    :: arrhc(3)

  INTEGER,PARAMETER :: mxrpd=2  ! max # of products returned by radchk sub
  CHARACTER(LEN=LEN(chem))  :: rdckpd(mxrpd)
  CHARACTER(LEN=LEN(coprod(1,1))) :: rdckcopd(mxrpd,SIZE(coprod,2))
  REAL :: sc(mxrpd)

  CHARACTER(LEN=9),PARAMETER :: progname='no3add_c7'

  Ci=0 ; Cf=0 ; ic1=0 ; ic2=0
  fract(:)=0.
  tgroup(:)=group(:)  ;  tbond(:,:)=bond(:,:)

! Find Ci and Cf (i.e. the node number in bond) such that Cdf=Cdi-O- 
  iloop: DO i=1,3,2
    IF (cdtable(i)/=0) THEN
      IF (cdeth(i,1)/=0) THEN
        Ci=cdtable(i)   ; Cf=cdtable(i+1) ; ic1=i   ; ic2=i+1

      ELSE IF (cdeth(i+1,1)/=0) THEN
        Ci=cdtable(i+1) ; Cf=cdtable(i)   ; ic1=i+1 ; ic2=i

      ELSE IF ((cdeth(i,1)==0).AND.(cdeth(i+1,1)==0)) THEN
        CYCLE iloop
      ENDIF

      fract(ic1)=0.5  ; fract(ic2)=0.5
      CALL kaddno3(tbond,tgroup,cdtable(ic1),cdtable(ic2),arrhc)
      rtot=arrhc(1)

! Reaction
!-----------

! set the step for the loop
      p = Cf-Ci

      DO j=Ci,Cf,p
! Change the double bond into a single bond
        tbond(Cf,Ci)=1  ;  tbond(Ci,Cf)=1
        pold='Cd' ; pnew='C' 
        tempkg=group(Ci)
        CALL swap(tempkg,pold,tgroup(Ci),pnew)
        tempkg=group(Cf)
        CALL swap(tempkg,pold,tgroup(Cf),pnew)
        ic=0
        IF      (j==Ci)              THEN ; ic=ic1
        ELSE IF (j==Cf)              THEN ; ic=ic2
        ELSE IF ((j/=Ci).OR.(j/=Cf)) THEN ; fract(ic)=0.
        ENDIF
        IF (fract(ic)>0.) THEN
          CALL addrx(progname,chem,nr,flag)
          tarrhc(nr,1)=rtot*fract(ic)
          CALL addref(progname,'JK10KEV000',nrxref(nr),rxref(nr,:),chem)

! Add NO3 to one carbon
          nc=INDEX(tgroup(j),' ') ; tgroup(j)(nc:nc+5)='(ONO2)' 

! Add '.' to the other one
          IF (j==Ci) THEN
             nc=INDEX(tgroup(Cf),' ') ; tgroup(Cf)(nc:nc)='.'
          ELSE IF (j==Cf) THEN
             nc=INDEX(tgroup(Ci),' ') ; tgroup(Ci)(nc:nc)='.'
          ENDIF
          CALL rebond(tbond,tgroup,tempkc,tempnring)
          CALL radchk(tempkc,rdckpd,rdckcopd,nip,sc,nrxref(nr),rxref(nr,:))
          pchem(nr)=rdckpd(1)
          CALL stdchm(pchem(nr))
          coprod(nr,:)=rdckcopd(1,:)
          sc_del(nr,1)=sc(1)
          IF (nip==2) THEN
            flag_del(nr)=1
            pchem_del(nr)=rdckpd(2)
            sc_del(nr,2)=sc(2)
            CALL stdchm(pchem_del(nr))
            coprod_del(nr,:)=rdckcopd(2,:)
          ENDIF

! Reset groups and bonds:
          tgroup(:)=group(:)  ;  tbond(:,:)=bond(:,:)
        ENDIF
      ENDDO
    ENDIF
  ENDDO iloop

END SUBROUTINE no3add_c7

!=======================================================================
! PURPOSE : Find rate constants for NO3 addition on VOC double bonds 
! based on Kerdouci et al., 2010, 2014                                        
!=======================================================================
SUBROUTINE kaddno3(tbond,tgroup,i1,i2,arrhc)
  USE keyparameter, ONLY: mxlfo,mxcp
  USE mapping, ONLY: gettrack
  USE reactool, ONLY: rebond      
  USE ringtool, ONLY: findring, ring_data
  IMPLICIT NONE

  INTEGER,INTENT(IN) :: tbond(:,:)           ! node matrix
  CHARACTER(LEN=*),INTENT(IN) :: tgroup(:)   ! group matrix
  INTEGER,INTENT(IN) :: i1, i2               ! group number bearing the Cd
  REAL,INTENT(OUT)   :: arrhc(:)             ! arrhenius coef. for NO3 addition

  INTEGER :: i,j,k
  REAL    :: mult
  INTEGER :: nether,nca,ngr,begrg
  INTEGER :: ring(SIZE(tgroup))     ! rg index for nodes, current ring (0=no, 1=yes)
  INTEGER :: rngflg        ! 0='no ring', 1='yes ring'
  INTEGER :: lring         ! length of the ring bearing the double bond
  INTEGER :: nring_ind,trackrg(12,SIZE(tgroup))
  LOGICAL :: ring_ind(12,SIZE(tgroup))

  arrhc(:)=0.
  
  ngr=COUNT(tgroup/=' ')
  nca=0 ; lring=0
  DO i=1,ngr
    IF (tgroup(i)(1:1)=='C') nca=nca+1
  ENDDO

! check ig value is ok
  IF (tbond(i1,i2)/=2) THEN
    WRITE(6,*) '--error--, in raddno3, called with a simple bond'
    STOP "in raddno3"
  ENDIF

! ------------------
! FIND k0 VALUE
! ------------------

! Values are from Kerdouci et al., 2014
  IF ( ((tgroup(i1)(1:4)=='CdH2') .AND. (tgroup(i2)(1:4)=='CdH ')) .OR. &
       ((tgroup(i1)(1:4)=='CdH ') .AND. (tgroup(i2)(1:4)=='CdH2')) .OR. &
       ((tgroup(i1)(1:4)=='CdH2') .AND. (tgroup(i2)(1:4)=='CdH(')) .OR. &
       ((tgroup(i1)(1:4)=='CdH(') .AND. (tgroup(i2)(1:4)=='CdH2')) ) THEN
    arrhc(1)=1.14E-14 ;  arrhc(2)=0.  ;  arrhc(3)=0.

! effect of the lenght of the carbon chain
    IF (nca>3) THEN
      mult= -0.64 + 3.14*(1.-EXP(-0.17*(nca))) ; arrhc(1)=arrhc(1)*mult
    ENDIF
  
  ELSE IF ( ((tgroup(i1)(1:4)=='CdH2').AND.(tgroup(i2)(1:3)=='Cd ')) .OR. &
            ((tgroup(i1)(1:3)=='Cd ') .AND.(tgroup(i2)(1:4)=='CdH2')).OR. &
            ((tgroup(i1)(1:4)=='CdH ').AND.(tgroup(i2)(1:4)=='CdH ')).OR. &
            ((tgroup(i1)(1:4)=='CdH ').AND.(tgroup(i2)(1:4)=='CdH(')).OR. &
            ((tgroup(i1)(1:4)=='CdH(').AND.(tgroup(i2)(1:4)=='CdH ')).OR. &
            ((tgroup(i1)(1:4)=='CdH2').AND.(tgroup(i2)(1:3)=='Cd(')) .OR. &
            ((tgroup(i1)(1:3)=='Cd(') .AND.(tgroup(i2)(1:4)=='CdH2')).OR. &
            ((tgroup(i1)(1:4)=='CdH(').AND.(tgroup(i2)(1:4)=='CdH(')) ) THEN
    arrhc(1)=3.40E-13 ; arrhc(2)=0. ; arrhc(3)=0.

  ELSE IF ( ((tgroup(i1)(1:4)=='CdH ').AND.(tgroup(i2)(1:3)=='Cd ')) .OR. &
            ((tgroup(i1)(1:4)=='CdH ').AND.(tgroup(i2)(1:3)=='Cd(')) .OR. &
            ((tgroup(i1)(1:4)=='CdH(').AND.(tgroup(i2)(1:3)=='Cd ')) .OR. &
            ((tgroup(i1)(1:4)=='CdH(').AND.(tgroup(i2)(1:3)=='Cd(')) .OR. &
            ((tgroup(i1)(1:3)=='Cd ') .AND.(tgroup(i2)(1:4)=='CdH ')).OR. &
            ((tgroup(i1)(1:3)=='Cd ') .AND.(tgroup(i2)(1:4)=='CdH(')).OR. &
            ((tgroup(i1)(1:3)=='Cd(') .AND.(tgroup(i2)(1:4)=='CdH ')).OR. &
            ((tgroup(i1)(1:3)=='Cd(') .AND.(tgroup(i2)(1:4)=='CdH(')) ) THEN
    arrhc(1)=8.35E-12 ; arrhc(2)=0. ; arrhc(3)=0.

  ELSE IF ( ((tgroup(i1)(1:3)=='Cd ').AND.(tgroup(i2)(1:3)=='Cd ')).OR. &
            ((tgroup(i1)(1:3)=='Cd(').AND.(tgroup(i2)(1:3)=='Cd(')).OR. &
            ((tgroup(i1)(1:3)=='Cd ').AND.(tgroup(i2)(1:3)=='Cd(')).OR. &
            ((tgroup(i1)(1:3)=='Cd(').AND.(tgroup(i2)(1:3)=='Cd ')) ) THEN
    arrhc(1)=4.50E-11 ; arrhc(2)=0. ; arrhc(3)=0.

  ELSE IF ((tgroup(i1)(1:4)=='CdH2').AND.(tgroup(i2)(1:4)=='CdH2')) THEN
    arrhc(1)=4.88E-18 ;  arrhc(2)=2. ; arrhc(3)=2282.

  ELSE
    WRITE(6,'(a)') '--error-- in raddno3, NO3 addition with groups:'
    WRITE(6,'(a)') '1= ',TRIM(tgroup(i1)),'  2=',TRIM(tgroup(i2))
    STOP "in raddno3"
  ENDIF

! effect of the lenght of the carbon chain
  IF ( (((tgroup(i1)(1:4)=='CdH2').AND.(tgroup(i2)(1:3)=='Cd ')).OR.&
       ((tgroup(i1)(1:3)=='Cd ') .AND.(tgroup(i2)(1:4)=='CdH2'))).AND.&
       (nca>4) ) THEN
    mult= -1.16 + 3.06*(1. -EXP(-0.32*nca)) ; arrhc(1)=arrhc(1)*mult
  ENDIF

! effect of the lenght of the cycle chain
  CALL findring(i1,i2,ngr,tbond,rngflg,ring)
  mult=1.
  DO j=1,2
    IF (rngflg==1) EXIT
    IF (j==1) begrg=i1
    IF (j==2) begrg=i2
    IF (ring(begrg)==0) CYCLE
    CALL ring_data(begrg,nca,tbond,tgroup,nring_ind,ring_ind,trackrg)
    mult=1.
    DO i=1,nring_ind
      lring=COUNT(ring_ind(i,:).EQV..TRUE.)  
      IF (lring==3) mult=mult*0.033
      IF (lring==4) mult=mult*0.73
      IF (lring==5) mult=mult*2.32
      IF (lring==6) mult=mult*0.93
      IF (lring==7) mult=mult*1.34
    ENDDO
  ENDDO    
  arrhc(1)=arrhc(1)*mult

! ---------------------------------------
! FIND MULTIPLIERS BASED ON SUBSTITUENTS
! ---------------------------------------

! conjugated double bonds factors

! on alpha carbons:
  nether=0
  DO i=1,ngr
    mult=1.
    IF ((tbond(i1,i)==1).OR.(tbond(i2,i)==1).OR. &
        (tbond(i1,i)==3).OR.(tbond(i2,i)==3) ) THEN

! simple alkyl:
      IF (tgroup(i)(1:3)=='CH3') THEN
        mult= 1. ; arrhc(1)=arrhc(1)*mult
  
      ELSE IF (tgroup(i)(1:4) == 'CH2 ') THEN
        j1loop: DO j=1,ngr 
          IF ((tbond(i,j)/=0).AND.(j/=i1).AND.(j/=i2)) THEN
            IF (tbond(i,j)==3) THEN
              DO k=1,ngr
                IF ((tbond(j,k)==3).AND.(tgroup(k)(1:2)=='CO')) THEN
                  mult=0.1 ;  EXIT j1loop
                ENDIF
              ENDDO
            ELSE IF ((tbond(i,j)==1).AND.(tgroup(j)=='CH2(OH)')) THEN
              mult=0.68
            ELSE IF ((tbond(i,j)==1).AND.(tgroup(j)(1:2)=='CO')) THEN
              mult=0.19
            ELSE 
              mult=1.15
            ENDIF
          ENDIF
        ENDDO j1loop
        arrhc(1)=arrhc(1)*mult
      
      ELSE IF(tgroup(i)(1:7) == 'CH2(OH)') THEN 
        mult= 0.63 ; arrhc(1)=arrhc(1)*mult
      
      ELSE IF(tgroup(i)(1:3) == 'CH ' ) THEN
        mult= 0.97 ; arrhc(1)=arrhc(1)*mult
      
      ELSE IF(tgroup(i)(1:6) == 'CH(OH)' ) THEN
        mult= 0.89 ; arrhc(1)=arrhc(1)*mult
      
      ELSE IF (tgroup(i)(1:2) == 'C('  ) THEN
        IF (tgroup(i)(1:5) == 'C(OH)'  ) THEN
          mult= 0.91 ; arrhc(1)=arrhc(1)*mult
        ELSE
          DO j=1,ngr
            IF ((tbond(i,j)==3).AND.(j/=i1).AND.(j/=i2)) THEN ; mult= 0.3
            ELSE                                              ; mult= 0.64
            ENDIF
          ENDDO
          arrhc(1)=arrhc(1)*mult
        ENDIF
      
      ELSE IF (tgroup(i)(1:2) == 'C '  ) THEN
        DO j=1,ngr
          IF ((tbond(i,j)==3).AND.(j/=i1).AND.(j/=i2)) THEN ; mult= 0.3
          ELSE                                              ; mult= 0.64
          ENDIF
        ENDDO
        arrhc(1)=arrhc(1)*mult
      ENDIF

! overwrite for carbonyls, esters and aldehydes values 
      IF (tgroup(i)(1:2)=='CO') THEN
        DO j=1,ngr
          IF      ((tbond(i,j)==3).AND.(j/=i1).AND.(j/=i2)) THEN ; mult= 0.008
          ELSE IF ((tbond(i,j)==1).AND.(j/=i1).AND.(j/=i2)) THEN ; mult= 0.02
          ENDIF
        ENDDO
        IF (tgroup(i)(1:9)=='CO(OONO2)') mult=0.0004
        arrhc(1)=arrhc(1)*mult

      ELSE IF (tgroup(i)(1:3)=='CHO') THEN
        mult= 0.002 ; arrhc(1)=arrhc(1)*mult
      ENDIF
! factor for beta-aldehyde, value not in the article, given by BPV          
      DO j=1,ngr
        IF ((tbond(i,j)==1).AND.(tgroup(j)(1:3)=='CHO').AND. &
            (INDEX(tgroup(i),'O')==0) ) THEN
          mult= 1.67  ;  arrhc(1)=arrhc(1)*mult
        ENDIF
      ENDDO

! Ethers and esters
! For acetal, consider only one -O- influence
      IF (INDEX(tgroup(i),'-O-')/=0) THEN
        DO j=1,ngr
          IF ((tbond(i,j)==3).AND.(j/=i1).AND.(j/=i2)) THEN
!! ether in alpha of ig
            IF      (tgroup(j)(1:3)=='CH3')  THEN ; mult=61.
            ELSE IF (tgroup(j)(1:4)=='CH2 ') THEN ; mult=127.
            ELSE IF (tgroup(j)(1:4)=='CH2(') THEN ; mult=127.
            ELSE IF (tgroup(j)(1:3)=='CH ' ) THEN ; mult=150.
            ELSE IF (tgroup(j)(1:3)=='CH(' ) THEN ; mult=150.
            ELSE IF (tgroup(j)(1:2)=='C('  ) THEN ; mult=294.
            ELSE IF (tgroup(j)(1:2)=='C '  ) THEN ; mult=294.
            ELSE IF (tgroup(j)(1:2)=='CO')   THEN ; mult=0.46  
            ENDIF
            arrhc(1)=arrhc(1)*mult
          ENDIF                 
        ENDDO
      ENDIF

! conjugated double bonds factors
      IF (tgroup(i)(1:4)=='CdH ') THEN
        mult=1
        DO j=1,ngr
          IF ((tbond(i,j)==2).AND.(j/=i1).AND.(j/=i2)) THEN
            IF      (tgroup(j)(1:4)=='CdH2') THEN ; mult=3.6
            ELSE IF (tgroup(j)(1:4)=='CdH ') THEN ; mult=18.   ! mean value between 12 for cis and 24 for trans
            ELSE IF (tgroup(j)(1:3)=='Cd ')  THEN ; mult=5.4
            ENDIF
          ENDIF
        ENDDO
        arrhc(1)=arrhc(1)*mult
      ENDIF

      IF (tgroup(i)(1:3)=='Cd ') THEN
        mult=1
        DO j=1,ngr
          IF ((tbond(i,j)==2).AND.(j/=i1).AND.(j/=i2)) THEN
            IF      (tgroup(j)(1:4)=='CdH2') THEN ; mult=2.1
            ELSE IF (tgroup(j)(1:4)=='CdH ') THEN ; mult=5.4      
            ENDIF
          ENDIF
        ENDDO
        arrhc(1)=arrhc(1)*mult
      ENDIF

    ENDIF
  ENDDO

END SUBROUTINE kaddno3

  
END  MODULE no3addtool
