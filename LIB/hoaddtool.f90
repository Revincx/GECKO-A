MODULE hoaddtool
IMPLICIT NONE 
CONTAINS
! SUBROUTINE hoadd_c1(chem,bond,group,ncd,cdtable,conjug, &
!                   nr,flag,tarrhc,pchem,coprod,flag_del,&
!                   pchem_del,coprod_del,sc_del)
! SUBROUTINE hoadd_c2(chem,bond,group, &
!                     ncd,conjug,cdtable,cdsub,cdcarbo, &
!                     nr,flag,tarrhc,pchem,coprod)
! SUBROUTINE hoadd_c3(chem,bond,group, &
!                     cdtable,cdsub,&
!                     nr,flag,tarrhc,pchem,coprod,flag_del,&
!                     pchem_del,coprod_del,sc_del)
! SUBROUTINE hoadd_c4(chem,bond,group,ncd,cdtable,conjug,cdsub,cdcarbo, &
!                   nr,flag,tarrhc,pchem,coprod,flag_del,&
!                   pchem_del,coprod_del,sc_del)
! SUBROUTINE hoadd_c5(chem,bond,group, &
!                     ncd,conjug,cdtable,cdcarbo, &
!                     nr,flag,tarrhc,pchem,coprod,flag_del, &
!                     pchem_del,coprod_del,sc_del)
! SUBROUTINE hoadd_c6(chem,bond,group,nr,flag,tarrhc,pchem,coprod)
! SUBROUTINE hoadd_c7(chem,bond,group,ncd,cdtable,cdeth,&
!                     nr,flag,tarrhc,pchem,coprod)
! SUBROUTINE hoadd_rate(bond,group,j0,j1,j2,j3,arrhc,w1,w2)
! SUBROUTINE hokadd_zie(ncdcase,group,j0,j1,j3,arrhc)
! SUBROUTINE hokaddc2(bond,group,cdtable,cdsub,ci,cf,arrhc)
! SUBROUTINE hokaddc3(bond,group,cdtable,cdsub,arrhc)
!
!=======================================================================
! PURPOSE : compute the reaction rate for OH addition on >C=C< bond for
! "case 1", i.e. species without C=O group conjugated with the C=C.   
!=======================================================================
SUBROUTINE hoadd_c1(chem,bond,group,cdtable, &
                  nr,flag,tarrhc,pchem,coprod,flag_del,&
                  pchem_del,coprod_del,sc_del,nrxref,rxref)
  USE keyflag, ONLY: kohadd_sar
  USE atomtool, ONLY: cnum
  USE reactool, ONLY: swap, rebond
  USE normchem, ONLY: stdchm
  USE radchktool, ONLY: radchk
  USE toolbox, ONLY: stoperr,addrx,addref
  IMPLICIT NONE

  CHARACTER(LEN=*),INTENT(IN) :: chem      ! chemical formula
  INTEGER,INTENT(IN) :: bond(:,:)          ! node matrix
  CHARACTER(LEN=*),INTENT(IN) :: group(:)  ! group matrix
  INTEGER,INTENT(IN) :: cdtable(:)         ! node # bearing a "Cd"

  INTEGER,INTENT(INOUT) :: nr              ! # of reactions attached to chem (incremented here)   
  INTEGER,INTENT(INOUT) :: flag(:)         ! flag to activate reaction i
  REAL,INTENT(INOUT)    :: tarrhc(:,:)     ! arrhenius coefficient for reaction i
  CHARACTER(LEN=*),INTENT(OUT)   :: pchem(:)        ! "main" product of reaction channel i
  CHARACTER(LEN=*),INTENT(OUT)   :: coprod(:,:)     ! coproducts (short names) for reaction channel i
  INTEGER,INTENT(INOUT)          :: flag_del(:)     ! flag for a 2nd channel of reaction i (delocalisation)
  CHARACTER(LEN=*),INTENT(INOUT) :: pchem_del(:)    ! "main" product for the 2nd channel
  CHARACTER(LEN=*),INTENT(INOUT) :: coprod_del(:,:) ! coproducts for 2nd channel
  REAL,INTENT(INOUT)             :: sc_del(:,:)     ! ratio od the 2 channels (for delocalisation)
  INTEGER,INTENT(INOUT)          :: nrxref(:)       ! # of references added for each reaction
  CHARACTER(LEN=*),INTENT(INOUT) :: rxref(:,:)      ! references (short codes) for the reactions

  CHARACTER(LEN=LEN(group(1))) :: tgroup(SIZE(group)), pold, pnew
  CHARACTER(LEN=LEN(chem)) :: tempkc
  INTEGER :: tbond(SIZE(bond,1),SIZE(bond,1))
  INTEGER :: nring
  INTEGER :: i,j0,j1,j2,j3,nc
  INTEGER :: nip

  INTEGER,PARAMETER :: mxrpd=2  ! max # of products returned by radchk sub
  CHARACTER(LEN=LEN(chem))  :: rdckpd(mxrpd)
  CHARACTER(LEN=LEN(coprod(1,1))) :: rdckcopd(mxrpd,SIZE(coprod,2))
  REAL :: sc(mxrpd)

  CHARACTER(LEN=9),PARAMETER :: progname='hoadd_c1'
  CHARACTER(LEN=70)          :: mesg

  tgroup(:)=group(:)  ;  tbond(:,:)=bond(:,:) 

! loop over each Cd
  cdloop: DO i=1,4  ! 4 is size of cdtable ... 
    IF (cdtable(i)==0) CYCLE cdloop 
    j0=cdtable(i)
    CALL addrx(progname,chem,nr,flag)
        
! find partner double-bond carbons (j1 is beta position, j2 is gamma 
! position and j3 is delta position in C=C-C=C structure)
    j1=0 ; j2=0 ; j3=0
    IF (i==1) THEN
      j1=cdtable(2)
      IF (cdtable(3)/=0) THEN
        IF (bond(cdtable(2),cdtable(3))==1) THEN
          j2=cdtable(3) ; j3=cdtable(4)
        ENDIF
      ENDIF
    ELSE IF (i==2) THEN
      j1=cdtable(1)
      IF (cdtable(3)/=0) THEN
        IF (bond(cdtable(1),cdtable(3))==1) THEN
          j2=cdtable(3) ; j3=cdtable(4)
        ENDIF
      ENDIF
    ELSE IF (i==3) THEN
      j1=cdtable(4)
    ELSE IF (i==4) THEN
      j1=cdtable(3)
      IF (bond(cdtable(2),cdtable(3))==1) THEN
        j2=cdtable(2) ; j3=cdtable(1)
      ELSE IF (bond(cdtable(1),cdtable(3))==1) THEN
        j2=cdtable(1) ; j3=cdtable(2)
      ENDIF
    ENDIF

! -----------------------------
! treat non conjugated C=C bond
! -----------------------------
    IF (j3==0) THEN

      IF (kohadd_sar==2) THEN         ! SAR based on Atkinson, Ziemann ...
        CALL hokadd_zie(1,group,j0,j1,j3,tarrhc(nr,:),&
                        nrxref(nr),rxref(nr,:))
     
      ELSE IF (kohadd_sar==3) THEN    ! magnify SAR - Jenkin et al. 2018
        CALL hoadd_rate(bond,group,j0,j1,j2,j3,tarrhc(nr,:))
        CALL addref(progname,'OHADDSAR1',nrxref(nr),rxref(nr,:),chem)
        CALL addref(progname,'MJ18KMV000',nrxref(nr),rxref(nr,:),chem)
          
! assign rate constant 
      ELSE        
        CALL addref(progname,'OHADDSAR1',nrxref(nr),rxref(nr,:),chem)
        CALL addref(progname,'DUMMYREF',nrxref(nr),rxref(nr,:),chem)
        IF (group(j1)(1:4)=='CdH2')     THEN ; tarrhc(nr,1)=0.45E-11
        ELSE IF (group(j1)(1:3)=='CdH') THEN ; tarrhc(nr,1)=3.0E-11
        ELSE IF (group(j1)(1:2)=='Cd')  THEN ; tarrhc(nr,1)=5.5E-11
        ELSE 
          mesg="A carbon in a C=C bond was found with no reactivity"
          CALL stoperr(progname,mesg,chem)
        ENDIF
      ENDIF

! convert I to single bond carbon and change Cd to C
      tbond(j0,j1)=1 ; tbond(j1,j0)=1
      pold='Cd' ; pnew='C'      
      CALL swap(group(j0),pold,tgroup(j0),pnew)
      pold='Cd' ; pnew='C'      
      CALL swap(group(j1),pold,tgroup(j1),pnew)

! add (OH) to j0 carbon, add radical dot to j1:
      nc=INDEX(tgroup(j0),' ') ; tgroup(j0)(nc:nc+3)='(OH)'
      nc=INDEX(tgroup(j1),' ') ; tgroup(j1)(nc:nc)='.'

! rebuild, check, rename and find co-products:
      CALL rebond(tbond,tgroup,tempkc,nring)
      CALL radchk(tempkc,rdckpd,rdckcopd,nip,sc,nrxref(nr),rxref(nr,:))
      pchem(nr)=rdckpd(1)
      CALL stdchm(pchem(nr))
      coprod(nr,:)=rdckcopd(1,:)
      IF (nip>1) THEN
        sc_del(nr,1)=sc(1) ; sc_del(nr,2)=sc(2) ; flag_del(nr)=1
        pchem_del(nr)=rdckpd(2)
        CALL stdchm(pchem_del(nr))
        coprod_del(nr,:)=rdckcopd(2,:)
      ENDIF

! reset groups,bonds:
      tbond(:,:)=bond(:,:)  ;  tgroup(:)=group(:)
    ENDIF
        
! -----------------------------
! treat conjugated C=C-C=C bond
! -----------------------------

! assign rate constant for conjugated C=C bond (i.e radical formed is a 
! C.-C=C structure) and "weight" the 2 positions for the radical
    IF (j3/=0) THEN
      IF (kohadd_sar==2) &    ! SAR based on Atkinson, Ziemann ...
        CALL hokadd_zie(1,group,j0,j1,j3,tarrhc(nr,:), &
                        nrxref(nr),rxref(nr,:))

      IF (kohadd_sar==3) THEN    ! magnify SAR - Jenkin et al. 2018
        CALL hoadd_rate(bond,group,j0,j1,j2,j3,tarrhc(nr,:))
        CALL addref(progname,'MJ18KMV000',nrxref(nr),rxref(nr,:),chem)
        CALL addref(progname,'OHADDSAR1',nrxref(nr),rxref(nr,:),chem)
      ENDIF
          
! check that a rate constant was set
      IF (tarrhc(nr,1)==0.) THEN
        mesg="A carbon in a C=C bond was found with no reactivity"
        CALL stoperr(progname,mesg,chem)
      ENDIF

! convert I to single bond carbon:
      tbond(j0,j1)=1 ; tbond(j1,j0)=1
      pold='Cd' ; pnew='C'      
      CALL swap(group(j0),pold,tgroup(j0),pnew)

! convert J1 to single bond C 
      pold='Cd' ; pnew='C'      
      CALL swap(group(j1),pold,tgroup(j1),pnew)

! add (OH) to I carbon, add radical dot to J1:
      nc=INDEX(tgroup(j0),' ') ; tgroup(j0)(nc:nc+3)='(OH)'
      nc=INDEX(tgroup(j1),' ') ; tgroup(j1)(nc:nc)='.'

! rebuild, check, and find co-products:
      CALL rebond(tbond,tgroup,tempkc,nring)
      CALL radchk(tempkc,rdckpd,rdckcopd,nip,sc,nrxref(nr),rxref(nr,:))
      pchem(nr)=rdckpd(1)
      CALL stdchm(pchem(nr))
      coprod(nr,:)=rdckcopd(1,:)
      IF (nip>1) THEN
        sc_del(nr,1)=sc(1) ; sc_del(nr,2)=sc(2) ; flag_del(nr)=1
        pchem_del(nr)=rdckpd(2)
        CALL stdchm(pchem_del(nr))
        coprod_del(nr,:)=rdckcopd(2,:)
      ENDIF

! reset groups,bonds:
      tgroup(:)=group(:) ; tbond(:,:)=bond(:,:)
    ENDIF
  ENDDO cdloop
END SUBROUTINE hoadd_c1

!=======================================================================
! PURPOSE : computes the reaction rate for OH addition on >C=C-C=O bond 
! (case 2). The -CO-C=C-C=C-CO- structure is not taken into account by 
! this routine (see case 3).                   
!=======================================================================
SUBROUTINE hoadd_c2(chem,bond,group,cdtable,cdsub,cdcarbo, &
                     nr,flag,tarrhc,pchem,coprod,nrxref,rxref)
  USE keyflag, ONLY: kohadd_sar
  USE reactool
  USE normchem
  USE radchktool, ONLY: radchk
  USE toolbox, ONLY: stoperr,addrx,addref
  IMPLICIT NONE

  CHARACTER(LEN=*),INTENT(IN) :: chem      ! chemical formula
  INTEGER,INTENT(IN) :: bond(:,:)          ! bond matrix
  CHARACTER(LEN=*),INTENT(IN) :: group(:)  ! group matrix
  INTEGER,INTENT(IN) :: cdtable(:)         ! node bearing a "Cd", ordered by bonds
  INTEGER,INTENT(IN) :: cdsub(:)           ! # of C (including CO) bonded to Cd in cdtable(i)
  INTEGER,INTENT(IN) :: cdcarbo(:,:)       ! node # of CO bonded to the Cd in cdtable(i)

  INTEGER,INTENT(INOUT) :: nr         ! # of reactions attached to chem (incremented here)
  INTEGER,INTENT(INOUT) :: flag(:)    ! flag to activate reaction i
  REAL,INTENT(INOUT) :: tarrhc(:,:)   ! arrhenius coefficient for reaction i
  CHARACTER(LEN=*),INTENT(INOUT) :: pchem(:)    ! "main" product of reaction channel i
  CHARACTER(LEN=*),INTENT(INOUT) :: coprod(:,:) ! coproducts (short names) for reaction channel i
  INTEGER,INTENT(INOUT)          :: nrxref(:)       ! # of references added for each reaction
  CHARACTER(LEN=*),INTENT(INOUT) :: rxref(:,:)      ! references (short codes) for the reactions

  CHARACTER(LEN=LEN(group(1))) :: tgroup(SIZE(group)), pold, pnew
  CHARACTER(LEN=LEN(chem))     :: tempkc
  INTEGER :: tbond(SIZE(bond,1),SIZE(bond,2))
  INTEGER :: i,j,k,ci,cf,nbcarbi,nbcarbf,nc,nring
  INTEGER :: posi,posf
  INTEGER :: nip
  INTEGER :: j0,j1,j2,j3
  REAL    :: arrhc(3)

  INTEGER,PARAMETER :: mxrpd=2  ! max # of products returned by radchk sub
  CHARACTER(LEN=LEN(chem))  :: rdckpd(mxrpd)
  CHARACTER(LEN=LEN(coprod(1,1))) :: rdckcopd(mxrpd,SIZE(coprod,2))
  REAL :: sc(mxrpd)

  CHARACTER(LEN=9),PARAMETER :: progname='hoadd_c2 '
  CHARACTER(LEN=70)          :: mesg

! initialize
! -----------
  tgroup(:)=group(:)  ;  tbond(:,:)=bond(:,:)
  ci=0  ;  cf=0 
      
! One double bond structures with carbonyl 
! -----------------------------------------
! Note: only one single bond must be "active" in this subroutine and 
! either cdtable(1-2) or  cdtable(3-4) must be 0.
! Find Ci and Cf such as Cdf=Cdi-C=O and count total number of 
! carbonyls bonded to the double bond
  cdloop: DO i=1,3,2
    nbcarbi=0  ;  nbcarbf=0
    IF (cdtable(i)/=0) THEN
      IF (cdcarbo(i,1)/=0) THEN
        Ci=cdtable(i)  ;  Cf=cdtable(i+1)
        posi=i  ;  posf=i+1
        DO j=1,2
          IF (cdcarbo(posi,j)/=0) nbcarbi=nbcarbi+1
          IF (cdcarbo(posf,j)/=0) nbcarbf=nbcarbf+1
        ENDDO
        EXIT cdloop 
      ELSE IF (cdcarbo(i+1,1)/=0) THEN
        Ci=cdtable(i+1)  ;  Cf=cdtable(i)
        posi=i+1  ;  posf=i
        DO j=1,2
          IF (cdcarbo(posi,j)/=0) nbcarbi=nbcarbi+1
          IF (cdcarbo(posf,j)/=0) nbcarbf=nbcarbf+1
        ENDDO
        EXIT cdloop 
      ENDIF
    ENDIF
  ENDDO cdloop

  IF (kohadd_sar==2) &  ! compute the arrh. coef. for the structure
    CALL hokaddc2(bond,group,cdtable,cdsub,ci,cf,arrhc)
  
! perform reaction
! ----------------
       
! loop over the two possibilities for the OH radical addition to 
! the >C=C<. Swap Ci and Cf at the end of the do loop.
  DO i=1,2
    CALL addrx(progname,chem,nr,flag)

! assign rate constant with original gecko for OH addition in Ci position
    IF (kohadd_sar==2) THEN
      IF (nbcarbf>nbcarbi)  tarrhc(nr,1)=arrhc(1)*0.8 ; tarrhc(nr,3)=0.
      IF (nbcarbi==nbcarbf) tarrhc(nr,1)=arrhc(1)*0.5 ; tarrhc(nr,3)=0.
      IF (nbcarbf<nbcarbi)  tarrhc(nr,1)=arrhc(1)*0.2 ; tarrhc(nr,3)=0.
      CALL addref(progname,'OHADDSAR1',nrxref(nr),rxref(nr,:),chem)
      CALL addref(progname,'EK95KEV000',nrxref(nr),rxref(nr,:),chem) ! xxx add the reference 
    ENDIF    

! assign rate with Jenkin 2018 SAR
    IF (kohadd_sar==3) THEN
      j0=Ci  ;  j1=Cf  ;  j2=0  ;  j3=0 
      CALL hoadd_rate(bond,group,j0,j1,j2,j3,tarrhc(nr,:))
      CALL addref(progname,'OHADDSAR1',nrxref(nr),rxref(nr,:),chem)
      CALL addref(progname,'MJ18KMV000',nrxref(nr),rxref(nr,:),chem)
    ENDIF

! convert Cf to single bond carbon
    tbond(Cf,Ci)=1  ;  tbond(Ci,Cf)=1
    pold='Cd'  ;  pnew='C'
    CALL swap(group(Cf),pold,tgroup(Cf),pnew)

! convert Ci to single bond carbon
    pold='Cd'  ;  pnew='C'         
    CALL swap(group(Ci),pold,tgroup(Ci),pnew)

! add OH to Ci carbon and add radical dot to Cf
    nc=INDEX(tgroup(Ci),' ')  ;  tgroup(Ci)(nc:nc+3)='(OH)'
    nc=INDEX(tgroup(Cf),' ')  ;  tgroup(Cf)(nc:nc)='.'
   
! rebuild, check rename and find coproducts:
    CALL rebond(tbond,tgroup,tempkc,nring)
    CALL radchk(tempkc,rdckpd,rdckcopd,nip,sc,nrxref(nr),rxref(nr,:))
    pchem(nr)=rdckpd(1)
    coprod(nr,:)=rdckcopd(1,:)
    CALL stdchm(pchem(nr))
    IF (nip/=1) THEN
      mesg="unexpected 2 pdcts from radchk in hoadd_c2.f"
      CALL stoperr(progname,mesg,chem)
    ENDIF

! reset:
    tgroup(:)=group(:)  ;  tbond(:,:)=bond(:,:)
    k=Ci ; Ci=Cf  ; Cf=k
    k=nbcarbi ;  nbcarbi=nbcarbf  ;  nbcarbf=k
         
  ENDDO
END SUBROUTINE hoadd_c2

!=======================================================================
! PURPOSE: compute the reaction rate for OH addition on -CO-C=C-C=C-C=O 
! bond (case 3) only. 
!=======================================================================
SUBROUTINE hoadd_c3(chem,bond,group,cdtable,cdsub,&
                    nr,flag,tarrhc,pchem,coprod,flag_del,&
                    pchem_del,coprod_del,sc_del,nrxref,rxref)
  USE keyflag, ONLY: kohadd_sar
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
  INTEGER :: i,j,ncoh,nco2,nc,nring
  INTEGER :: nip
  INTEGER :: j0,j1,j2,j3
  REAL    :: arrhc(3)

  INTEGER,PARAMETER :: mxrpd=2  ! max # of products returned by radchk sub
  CHARACTER(LEN=LEN(chem))  :: rdckpd(mxrpd)
  CHARACTER(LEN=LEN(coprod(1,1))) :: rdckcopd(mxrpd,SIZE(coprod,2))
  REAL :: sc(mxrpd)

  CHARACTER(LEN=9),PARAMETER :: progname='hoadd_c3 '

! initialize
! -----------
  tgroup(:)=group(:)  ;  tbond(:,:)=bond(:,:)

  IF (kohadd_sar==2)  THEN
    CALL hokaddc3(bond,group,cdtable,cdsub,arrhc)
  ENDIF
      
! perform OH addition (only 1-2 addition is considered, 1-4 addition is 
! neglected). The same branching ratio is assumed for each channel (kstruc/4).
  DO i=1,3,2
    DO j=1,2
      IF (j==1) THEN ;  ncoh=cdtable(i)   ; nco2=cdtable(i+1)
      ELSE           ;  ncoh=cdtable(i+1) ; nco2=cdtable(i)
      ENDIF

! add 1 to the channel counter
      CALL addrx(progname,chem,nr,flag)

! assign rate constant for OH addition in Ci position
      IF (kohadd_sar==2) THEN
        tarrhc(nr,1)=arrhc(1)*0.25 ; tarrhc(nr,2)=0.  ;  tarrhc(nr,3)=0.
        CALL addref(progname,'OHADDSAR1',nrxref(nr),rxref(nr,:),chem)
        CALL addref(progname,'EK95KEV000',nrxref(nr),rxref(nr,:),chem)
      ENDIF

! assign rate with Jenkin 2018 SAR
      IF (kohadd_sar==3) THEN
        IF ((i==1).AND.(j==1)) THEN
          j0=cdtable(1) ; j1=cdtable(2) ; j2=cdtable(3) ; j3=cdtable(4)
        ELSE IF ((i==1).AND.(j==2)) THEN
          j0=cdtable(2) ; j1=cdtable(1) ; j2=0 ; j3=0
        ELSE IF ((i==3).AND.(j==1)) THEN
          j0=cdtable(3) ; j1=cdtable(4) ; j2=0 ; j3=0
        ELSE IF ((i==3).AND.(j==2)) THEN
          j0=cdtable(4) ; j1=cdtable(3) ; j2=cdtable(2) ; j3=cdtable(1)
        ENDIF
        CALL hoadd_rate(bond,group,j0,j1,j2,j3,tarrhc(nr,:))
        CALL addref(progname,'OHADDSAR1',nrxref(nr),rxref(nr,:),chem)
        CALL addref(progname,'MJ18KMV000',nrxref(nr),rxref(nr,:),chem)  
      ENDIF
               
! convert Cf to single bond carbon
      tbond(ncoh,nco2)=1 ; tbond(nco2,ncoh)=1
      pold='Cd' ; pnew='C'
      CALL swap(group(ncoh),pold,tgroup(ncoh),pnew)

! convert Ci to single bond carbon
      pold='Cd' ; pnew='C'         
      CALL swap(group(nco2),pold,tgroup(nco2),pnew)

! add OH to ncoh carbon and add radical dot to nco2
      nc=INDEX(tgroup(ncoh),' ') ; tgroup(ncoh)(nc:nc+3)='(OH)'
      nc=INDEX(tgroup(nco2),' ') ; tgroup(nco2)(nc:nc)='.'         
      
! rebuild, check, rename and find coproducts:
      CALL rebond(tbond,tgroup,tempkc,nring)
      CALL radchk(tempkc,rdckpd,rdckcopd,nip,sc,nrxref(nr),rxref(nr,:))
      pchem(nr)=rdckpd(1)
      CALL stdchm(pchem(nr))
      coprod(nr,:)=rdckcopd(1,:)
      IF (nip/=1) THEN
        flag_del(nr)=1 ; sc_del(nr,1)=sc(1) ; sc_del(nr,2)=sc(2)
        pchem_del(nr)=rdckpd(2)
        CALL stdchm(pchem_del(nr))
        coprod_del(nr,:)=rdckcopd(2,:)
      ENDIF

! reset
      tgroup(:)=group(:) ; tbond(:,:)=bond(:,:)     
    ENDDO
  ENDDO  

END SUBROUTINE hoadd_c3


!=======================================================================
! PURPOSE: compute the reaction rate for OH addition on -C=C-R-C=C-C=O-
! non-conjugated C=C bonds
!=======================================================================
SUBROUTINE hoadd_c4(chem,bond,group,cdtable,cdsub,cdcarbo, &
                  nr,flag,tarrhc,pchem,coprod,flag_del,&
                  pchem_del,coprod_del,sc_del,nrxref,rxref)
  USE toolbox, ONLY: stoperr,addref
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
  CHARACTER(LEN=*),INTENT(OUT) :: pchem(:)        ! "main" product of reaction channel i
  CHARACTER(LEN=*),INTENT(OUT) :: coprod(:,:)     ! coproducts (short names) for reaction channel i
  INTEGER,INTENT(INOUT)          :: flag_del(:)     ! flag for a 2nd channel of reaction i (delocalisation)
  CHARACTER(LEN=*),INTENT(INOUT) :: pchem_del(:)    ! "main" product for the 2nd channel
  CHARACTER(LEN=*),INTENT(INOUT) :: coprod_del(:,:) ! coproducts for 2nd channel
  REAL,INTENT(INOUT)             :: sc_del(:,:)     ! ratio od the 2 channels (for delocalisation)
  INTEGER,INTENT(INOUT)          :: nrxref(:)       ! # of references added for each reaction
  CHARACTER(LEN=*),INTENT(INOUT) :: rxref(:,:)      ! references (short codes) for the reactions

  INTEGER :: i
  INTEGER :: cpcdtable(SIZE(cdtable))

  CHARACTER(LEN=9),PARAMETER :: progname='hoadd_c4 '
  CHARACTER(LEN=70)           :: mesg

! perform reaction for each C=C, ignoring the other C=C bond 
  DO i=1,3,2

! check everything is ok
    IF (cdtable(i)==0) THEN
      mesg="Cd not found"
      CALL stoperr(progname,mesg,chem)
    ENDIF

! make a copy of cdtable and "kill" the Cd bond not considered
    cpcdtable(:)=0
    cpcdtable(i)=cdtable(i) ; cpcdtable(i+1)=cdtable(i+1)

! send to reaction subroutine
    IF (cdcarbo(i,1)/=0) THEN
      CALL hoadd_c2(chem,bond,group,cpcdtable,cdsub,cdcarbo,&
                   nr,flag,tarrhc,pchem,coprod,nrxref,rxref)
    ELSE IF (cdcarbo(i+1,1)/=0) THEN
      CALL hoadd_c2(chem,bond,group,cpcdtable,cdsub,cdcarbo,&
                   nr,flag,tarrhc,pchem,coprod,nrxref,rxref)
    ELSE
      CALL hoadd_c1(chem,bond,group,cpcdtable,&
            nr,flag,tarrhc,pchem,coprod,flag_del,&
            pchem_del,coprod_del,sc_del,nrxref,rxref)
    ENDIF
  ENDDO

END SUBROUTINE hoadd_c4

!=======================================================================
! PURPOSE: computes the reaction rate for OH addition to >C=C-C=C-C=O 
! bond (case 5) only.
!=======================================================================
SUBROUTINE hoadd_c5(chem,bond,group,cdtable,cdcarbo, &
                    nr,flag,tarrhc,pchem,coprod,flag_del, &
                    pchem_del,coprod_del,sc_del,nrxref,rxref)
  USE keyflag, ONLY: kohadd_sar
  USE reactool
  USE normchem
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
  INTEGER :: i,j0,j1,j2,j3,nc,nring,ialpha
  INTEGER :: nip

  INTEGER,PARAMETER :: mxrpd=2  ! max # of products returned by radchk sub
  CHARACTER(LEN=LEN(chem))  :: rdckpd(mxrpd)
  CHARACTER(LEN=LEN(coprod(1,1))) :: rdckcopd(mxrpd,SIZE(coprod,2))
  REAL :: sc(mxrpd)
  
  CHARACTER(LEN=9),PARAMETER :: progname='hoadd_c5 '
  CHARACTER(LEN=70)          :: mesg

  tgroup(:)=group(:) ; tbond(:,:)=bond(:,:)
  ialpha=0

! add OH only at position furthest from carbonyl group (ialpha)
!   Note: Only 1-2 and 1-4 addition are considered,               
!         those are treated as per hoadd_c1                       
  IF      ((cdcarbo(1,1)/=0).OR.(cdcarbo(2,1)/=0)) THEN ; ialpha=4
  ELSE IF ((cdcarbo(4,1)/=0).OR.(cdcarbo(3,1)/=0)) THEN ; ialpha=1
  ELSE
    mesg="Expected carbonyl not found "
    CALL stoperr(progname,mesg,chem)
  ENDIF
  i=ialpha ; j0=cdtable(i)

  CALL addrx(progname,chem,nr,flag)

! from C1.find partner double-bond carbons (j1 is beta position, j2 is
! gamma position and j3 is delta position in C=C-C=C structure)
  j1=0 ; j2=0 ; j3=0
  IF (i==1) THEN      ; j1=cdtable(2) ; j2=cdtable(3) ; j3=cdtable(4)
  ELSE IF (i==4) THEN ; j1=cdtable(3) ; j2=cdtable(2) ; j3=cdtable(1)
  ENDIF

! treat conjugated C=C-C=C bond
! -----------------------------
! assign rate constant for conjug. C=C (i.e. radical is a C.-C=C structure) 
  IF (cdcarbo(i,1)/=0) RETURN
  
  IF (kohadd_sar/=3) THEN  
!  The method used is the SAR of Peeters et al., 1997, "kinetic studies of reactions 
!  of alkylperoxy and haloalkylperoxy radicals with NO. A SAR for reactions of OH with  
!  alkenes and polyalkenes", in Chemical processes in atmospheric oxidation, Eurotrac.
    IF (group(j1)(1:3)=='CdH') THEN
      IF (group(j3)(1:4)=='CdH2')     THEN ; tarrhc(nr,1)=3.00E-11 
      ELSE IF (group(j3)(1:3)=='CdH') THEN ; tarrhc(nr,1)=3.75E-11
      ELSE IF (group(j3)(1:2)=='Cd')  THEN ; tarrhc(nr,1)=5.05E-11
      ENDIF
    ELSE IF (group(j1)(1:2)=='Cd') THEN
      IF (group(j3)(1:4)=='CdH2')     THEN ; tarrhc(nr,1)=5.65E-11
      ELSE IF (group(j3)(1:3)=='CdH') THEN ; tarrhc(nr,1)=8.35E-11
      ELSE IF (group(j3)(1:2)=='Cd')  THEN ; tarrhc(nr,1)=9.85E-11
      ENDIF
    ENDIF
  ENDIF
  
! Use Jenkin et al., 2018 
  IF (kohadd_sar==3) THEN
    CALL hoadd_rate(bond,group,j0,j1,j2,j3,tarrhc(nr,:))
    CALL addref(progname,'OHADDSAR1',nrxref(nr),rxref(nr,:),chem)
    CALL addref(progname,'MJ18KMV000',nrxref(nr),rxref(nr,:),chem)  
  ENDIF

! check that a rate constant was set
  IF (tarrhc(nr,1)==0.) THEN
    PRINT*, 'Warning, no rate constant for OH addition for: ',TRIM(chem)
    !mesg="No rate constant set "
    !CALL stoperr(progname,mesg,chem)
  ENDIF

! convert i to single bond carbon:
  tbond(j0,j1)=1 ; tbond(j1,j0)=1
  pold='Cd' ; pnew='C'
  CALL swap(group(j0),pold,tgroup(j0),pnew)

! convert j1 to single bond C 
  pold='Cd' ; pnew='C'
  CALL swap(group(j1),pold,tgroup(j1),pnew)

! add (OH) to i carbon, add radical dot to j1:
  nc=INDEX(tgroup(j0),' ') ; tgroup(j0)(nc:nc+3)='(OH)'
  nc=INDEX(tgroup(j1),' ') ; tgroup(j1)(nc:nc)='.'

! rebuild, check, rename and find co-products:
  CALL rebond(tbond,tgroup,tempkc,nring)
  CALL radchk(tempkc,rdckpd,rdckcopd,nip,sc,nrxref(nr),rxref(nr,:))
  pchem(nr)=rdckpd(1)
  CALL stdchm(pchem(nr))
  coprod(nr,:)=rdckcopd(1,:)
  IF (nip>1) THEN
    sc_del(nr,1)=sc(1) ; sc_del(nr,2)=sc(2) ; flag_del(nr)=1
    pchem_del(nr)=rdckpd(2)
    CALL stdchm(pchem_del(nr))
    coprod_del(nr,:)=rdckcopd(2,:)
  ENDIF

! reset groups,bonds:
  tbond(:,:)=bond(:,:)  ;  tgroup(:)=group(:)

END SUBROUTINE hoadd_c5

!=======================================================================
! PURPOSE: computes the reaction rate for OH addition to >C=C=O bond 
! (ketene) (case 6) only. 
!=======================================================================
SUBROUTINE hoadd_c6(chem,bond,group,nr,flag,tarrhc,pchem,coprod,nrxref,rxref)
  USE reactool, ONLY: swap,rebond
  USE normchem, ONLY: stdchm
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

  CHARACTER(LEN=LEN(group(1))) :: tgroup(SIZE(group)), pold, pnew
  CHARACTER(LEN=LEN(chem)) :: tempkc
  INTEGER :: tbond(SIZE(bond,1),SIZE(bond,1))
  INTEGER :: i,j,k,nring
  INTEGER :: nip

  INTEGER,PARAMETER :: mxrpd=2  ! max # of products returned by radchk sub
  CHARACTER(LEN=LEN(chem))  :: rdckpd(mxrpd)
  CHARACTER(LEN=LEN(coprod(1,1))) :: rdckcopd(mxrpd,SIZE(coprod,2))
  REAL :: sc(mxrpd)

  CHARACTER(LEN=9),PARAMETER :: progname='hoadd_c6 '
  CHARACTER(LEN=70)           :: mesg

  tgroup(:)=group(:) ; tbond(:,:)=bond(:,:)
      
! add OH only at carbonyl group
  DO i=1,SIZE(group)
    IF (INDEX(group(i),'CdO')/=0) THEN
      CALL addrx(progname,chem,nr,flag)

      pold='CdO'  ;  pnew='CO.'
      CALL swap(group(i),pold,tgroup(i),pnew)
      DO j=1,SIZE(group)
        IF(bond(i,j)==2)THEN
          tbond(i,j)=1  ;  tbond(j,i)=1
          IF(INDEX(group(j),'Cd')/=0) THEN
            k=INDEX(group(j),' ') ; tgroup(j)(2:k-1)=group(j)(3:k)
          ELSE
            mesg="Unidentified species"
            CALL stoperr(progname,mesg,chem)
          ENDIF
          k=INDEX(tgroup(j),' ')  ;  tgroup(j)(k:k+3)='(OH)'
        ENDIF
      ENDDO

! rebuild, check, rename and find co-products:
      CALL rebond(tbond,tgroup,tempkc,nring)
      CALL radchk(tempkc,rdckpd,rdckcopd,nip,sc,nrxref(nr),rxref(nr,:))
      pchem(nr)=rdckpd(1)
      CALL stdchm(pchem(nr))
      coprod(nr,:)=rdckcopd(1,:)
      IF (nip/=1) THEN
        mesg="unexpected 2 pdcts from radchk"
        CALL stoperr(progname,mesg,chem)
      ENDIF

! set rate constants based on CdH2=CdO. Uses OH rate of Hatakeyama (1985),
! "Rate constants and mechanism for reactions of ketenes with OH radicals 
! in air at 299±2K", Hatakeyama S., Honda S., Washida N. & Akimoto H., 
! Bull. Chem. Soc. Japan, 58(8), 2157-2162, 10.1246/bcsj.58.2157, 1985.
      tarrhc(nr,1)=1.79E-11  ; tarrhc(nr,2)=0.  ;  tarrhc(nr,3)=0.
      CALL addref(progname,'SARKETENE',nrxref(nr),rxref(nr,:),chem)  

! reset groups:
      tgroup(:)=group(:)

    ENDIF
  ENDDO

END SUBROUTINE hoadd_c6

!=======================================================================
! PURPOSE: compute the reaction rate for OH addition on O-C=C bond 
! (case 7). Species having a C=O group conjugated with the C=C are not 
! taken into account here.
!=======================================================================
SUBROUTINE hoadd_c7(chem,bond,group,cdtable,cdeth, &
                     nr,flag,tarrhc,pchem,coprod,flag_del, &
                    pchem_del,coprod_del,sc_del,nrxref,rxref)
  USE reactool, ONLY: swap, rebond
  USE ringtool, ONLY: findring
  USE normchem, ONLY: stdchm
  USE radchktool, ONLY: radchk
  USE toolbox, ONLY: stoperr,addrx,addref
  IMPLICIT NONE

  CHARACTER(LEN=*),INTENT(IN) :: chem      ! chemical formula
  INTEGER,INTENT(IN) :: bond(:,:)          ! bond matrix
  CHARACTER(LEN=*),INTENT(IN) :: group(:)  ! group matrix
  INTEGER,INTENT(IN) :: cdtable(:)         ! node bearing a "Cd", ordered by bonds
  INTEGER,INTENT(IN) :: cdeth(:,:)         ! node # of CO bonded to the Cd in cdtable(i)
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
  INTEGER :: ring(SIZE(group))
  INTEGER :: rngflg,p

  INTEGER :: i,j,nc,nca,Ci,Cf,ic1,ic2
  INTEGER :: nrg,ic,nip
  REAL    :: fract(SIZE(cdtable))

  INTEGER,PARAMETER :: mxrpd=2  ! max # of products returned by radchk sub
  CHARACTER(LEN=LEN(chem))  :: rdckpd(mxrpd)
  CHARACTER(LEN=LEN(coprod(1,1))) :: rdckcopd(mxrpd,SIZE(coprod,2))
  REAL :: sc(mxrpd)

  CHARACTER(LEN=9),PARAMETER :: progname='hoadd_c7 '
  CHARACTER(LEN=70)           :: mesg

  nca=0 ; Ci=0 ; Cf=0 ; nc=0 ; tempkg=' ' ; ic1=0 ; ic2=0 ; nip=0

  fract(:)=0. 
  tgroup(:)=group(:) ; tbond(:,:)=bond(:,:)
  nca=COUNT(tgroup/=' ')

! Find Ci and Cf (i.e. the node number in bond) such that Cdf=Cdi-O- 
  cdloop: DO i=1,3,2
    IF (cdtable(i)/=0) THEN
      IF (cdeth(i,1)/=0) THEN
        Ci=cdtable(i)   ;  Cf=cdtable(i+1) ; ic1=i   ;  ic2=i+1
      ELSE IF (cdeth(i+1,1)/=0) THEN
        Ci=cdtable(i+1) ;  Cf=cdtable(i)   ; ic1=i+1 ;  ic2=i
      ELSE IF ((cdeth(i,1)==0).AND.(cdeth(i+1,1)==0)) THEN
        CYCLE cdloop
      ENDIF
         
! Check if there's a ring
      CALL findring(Ci,Cf,nca,bond,rngflg,ring)

! For the moment, 50% for each channel
      fract(ic1)=0.5  ;  fract(ic2)=0.5
     
! REACTION
!-----------

! set the step for the loop
      p = Cf-Ci
      DO j=Ci,Cf,p

! Change the double bond into a single bond
        tbond(Cf,Ci)=1  ;  tbond(Ci,Cf)=1
        pold='Cd'  ;  pnew='C'  ; tempkg=group(Ci)
        CALL swap(tempkg,pold,tgroup(Ci),pnew)
        tempkg=group(Cf)
        CALL swap(tempkg,pold,tgroup(Cf),pnew)
        ic=0
    
        IF (j==Ci) THEN      ; ic=ic1
        ELSE IF (j==Cf) THEN ; ic=ic2
        ENDIF
        IF (ic==0) THEN
          mesg="ic not set!"
          CALL stoperr(progname,mesg,chem)
        ENDIF

        IF (fract(ic)>0.) THEN
          CALL addrx(progname,chem,nr,flag)
    
          IF (j==Ci) THEN
            CALL hoadd_rate(bond,group,Ci,Cf,0,0,tarrhc(nr,:))
          ELSE IF (j==Cf) THEN
            CALL hoadd_rate(bond,group,Cf,Ci,0,0,tarrhc(nr,:))
          ENDIF
          CALL addref(progname,'OHADDSAR1',nrxref(nr),rxref(nr,:),chem)
          CALL addref(progname,'MJ18KMV000',nrxref(nr),rxref(nr,:),chem)  

! Add OH to one carbon
          nc=INDEX(tgroup(j),' ')   ;  tgroup(j)(nc:nc+3)='(OH)' 

! Add '.' to the other one
          IF (j==Ci) THEN
            nc=INDEX(tgroup(Cf),' ') ; tgroup(Cf)(nc:nc)='.'
          ELSE IF (j==Cf) THEN
            nc=INDEX(tgroup(Ci),' ') ; tgroup(Ci)(nc:nc)='.'
          ENDIF
          CALL rebond(tbond,tgroup,tempkc,nrg)
          CALL radchk(tempkc,rdckpd,rdckcopd,nip,sc,nrxref(nr),rxref(nr,:))
          pchem(nr)=rdckpd(1)
          coprod(nr,:)=rdckcopd(1,:)
          CALL stdchm(pchem(nr))
          IF (nip>1) THEN
            sc_del(nr,1)=sc(1) ; sc_del(nr,2)=sc(2) ; flag_del(nr)=1
            pchem_del(nr)=rdckpd(2)
            CALL stdchm(pchem_del(nr))
            coprod_del(nr,:)=rdckcopd(2,:)
          ENDIF

! reset groups,bonds:
          tbond(:,:)=bond(:,:)  ;  tgroup(:)=group(:)
        ENDIF
      ENDDO
    ENDIF
  ENDDO cdloop

END SUBROUTINE hoadd_c7


!=======================================================================
! PURPOSE: Compute the reaction rate for OH addition on >C=C< bonds
! based on the SAR described by Jenkin et al., ACP, 2018.
!=======================================================================
SUBROUTINE hoadd_rate(bond,group,j0,j1,j2,j3,arrhc)
  USE keyparameter, ONLY: saru
  USE keyflag, ONLY: wrtsarinfo
  USE toolbox, ONLY: stoperr
  IMPLICIT NONE

  CHARACTER(LEN=*),INTENT(IN) :: group(:)
  INTEGER,INTENT(IN) :: bond(:,:)
  INTEGER,INTENT(IN) :: j0,j1,j2,j3
  REAL,INTENT(OUT) :: arrhc(:)

  INTEGER ::  i,j,k,l,nca,ngr,xx
  INTEGER  :: nfun                ! # of contribution to substituent factor
  INTEGER,PARAMETER :: mxfact=10
  REAL     :: af(mxfact)          ! contribution to A substituent factor
  REAL     :: bf(mxfact)          ! contribution to B substituent factor
  REAL     :: ax                  ! overall substituent A factor
  REAL     :: bx                  ! overall substituent B factor
  LOGICAL  :: loester             ! ester contribution required (flag) 
  CHARACTER(LEN=11),PARAMETER :: progname='hoadd_rate '

  arrhc(:)=0.
  
  nca=0 ; ngr=0
  DO i=1,SIZE(group)
    IF (group(i)(1:1)=='C') nca=nca+1
    IF (group(i)(1:1)/=' ') ngr=ngr+1
  ENDDO

  IF (wrtsarinfo) THEN                                     ! debug
    WRITE(saru,*) '   '                                     ! debug
    WRITE(saru,*) '=========================='              ! debug
    WRITE(saru,*) '=====   HOADD_RATE   ====='              ! debug
    DO i=1,ngr                                              ! debug
      WRITE(saru,'(a6,i2,a2,a15)') 'group(',i,')=',group(i) ! debug
    ENDDO                                                   ! debug
    WRITE(saru,*) 'j0=',j0                                  ! debug
    WRITE(saru,*) 'j1=',j1                                  ! debug
    WRITE(saru,*) 'j2=',j2                                  ! debug
    WRITE(saru,*) 'j3=',j3                                  ! debug
    WRITE(saru,*) 'group(j0)=',group(j0)                    ! debug
    WRITE(saru,*) 'group(j1)=',group(j1)                    ! debug
  ENDIF                                                     ! debug

! simple C=C (Table 10 in Jenkin et al., 2018)
  IF (j3==0) THEN
    IF (group(j1)(1:4)=='CdH2') THEN 
      arrhc(1)=2.04E-12  ;  arrhc(3)=-215.
    ELSE IF ((group(j1)(1:4)=='CdH ').OR.(group(j1)(1:4)=='CdH(')) THEN
      arrhc(1)=4.30E-12  ;  arrhc(3)=-540.
    ELSE
      arrhc(1)=8.13E-12  ;  arrhc(3)=-550.
    ENDIF
  ELSE

! conjugated C=C-C=C (Table 10 in Jenkin et al., 2018)
    IF (wrtsarinfo) THEN                                 ! debug
      WRITE(saru,*) 'group(j2)=',group(j2)                ! debug
      WRITE(saru,*) 'group(j3)=',group(j3)                ! debug
    ENDIF                                                 ! debug

    IF (group(j1)(1:3)=='CdH') THEN
      IF      (group(j3)(1:4)=='CdH2') THEN ; arrhc(1)= 6.74E-12
      ELSE IF (group(j3)(1:3)=='CdH' ) THEN ; arrhc(1)= 8.99E-12
      ELSE IF (group(j3)(1:2)=='Cd'  ) THEN ; arrhc(1)=10.56E-12
      ENDIF
    ELSE IF (group(j1)(1:2)=='Cd') THEN
      IF      (group(j3)(1:4)=='CdH2') THEN ; arrhc(1)=13.70E-12
      ELSE IF (group(j3)(1:3)=='CdH' ) THEN ; arrhc(1)=16.62E-12
      ELSE IF (group(j3)(1:2)=='Cd'  ) THEN ; arrhc(1)=22.24E-12
      ENDIF
    ENDIF 
    arrhc(3)=-445.
  ENDIF       
  IF (wrtsarinfo) WRITE(saru,*) 'arrhc(1:3)',arrhc(1:3)  ! debug        

! Substituent factors (Table 11 in Jenkin et al., 2018)
  nfun=0 ; af(:)=0. ; bf(:)=0. ; loester=.FALSE.
  factloop : DO i=1,2
    IF (i==1) THEN ; xx=j1
    ELSE           ; xx=j3  
      IF (j3==0) EXIT factloop
    ENDIF
    
    IF (INDEX(group(xx),'(OH)')/=0)  THEN ; nfun=nfun+1 ; af(nfun)=0.249 ; bf(nfun)=-515. ; ENDIF
    IF (INDEX(group(xx),'(NO2)')/=0) THEN ; arrhc(1)=0. ; arrhc(3)=0.    ; RETURN         ; ENDIF
    
    DO l=1,ngr
      IF (bond(xx,l)==1) THEN
        IF (group(l)(1:4)=='CO(O')        THEN ; nfun=nfun+1 ; af(nfun)=0.094 ; bf(nfun)=-480. ; CYCLE ; ENDIF
        IF (INDEX(group(l),'(OH)')/=0)    THEN ; nfun=nfun+1 ; af(nfun)=0.951 ; bf(nfun)=-190. ; ENDIF
        IF (INDEX(group(l),'(OOH)')/=0)   THEN ; nfun=nfun+1 ; af(nfun)=1.    ; bf(nfun)=-298.*log(1.2) ; ENDIF
        IF (INDEX(group(l),'(NO2)')/=0)   THEN ; nfun=nfun+1 ; af(nfun)=1.    ; bf(nfun)=-298.*log(0.3) ; ENDIF
        IF (INDEX(group(l),'(ONO2)')/=0)  THEN ; nfun=nfun+1 ; af(nfun)=1.    ; bf(nfun)=-298.*log(0.26) ; ENDIF
        IF (INDEX(group(l),'CO(OONO2)')/=0) THEN ; nfun=nfun+1 ; af(nfun)=1.  ; bf(nfun)=-298.*log(0.47) ; CYCLE ; ENDIF
        
        IF (group(l)(1:3)=='CO ') THEN
          DO j=1,ngr
            IF (bond(l,j)==3) THEN
              nfun=nfun+1 ; af(nfun)=0.094 ; bf(nfun)=-480. ! -C(=O)OR  
              loester=.TRUE.
              EXIT
            ENDIF
          ENDDO
          IF (.NOT. loester) THEN ; nfun=nfun+1 ; af(nfun)=0.328 ; bf(nfun)=-180. ; ENDIF     ! -C(=O)-  
        ENDIF
        IF (group(l)(1:3)=='CHO') THEN ; nfun=nfun+1 ; af(nfun)=0.423 ; bf(nfun)=145. ; ENDIF ! -C(=O)H  
      
        DO j=1,ngr  ! effect of beta functions
          IF ((bond(l,j)==1).AND.(xx/=j)) THEN
            IF (INDEX(group(j),'(OH)')/=0)   THEN ; nfun=nfun+1 ; af(nfun)=0.951 ; bf(nfun)=-190. ; ENDIF
            IF (INDEX(group(j),'(ONO2)')/=0) THEN ; nfun=nfun+1 ; af(nfun)=1.    ; bf(nfun)=-298.*log(0.6) ; ENDIF
          ELSE IF (bond(l,j)==3) THEN
            DO k=1,ngr
              IF ((bond(k,j)==3).AND.(k/=l).AND.(group(k)(1:3)=='CO ')) THEN ! -C-OC(=O)R 
                nfun=nfun+1 ; af(nfun)=0.319 ; bf(nfun)=-230.
                EXIT
              ENDIF
            ENDDO
          ENDIF
        ENDDO
      ELSE IF (bond(xx,l)==3) THEN
        DO j=1,ngr
          IF ((bond(l,j)==3).AND.(j/=xx) ) THEN
            IF (group(j)(1:2)=='CO') THEN 
              nfun=nfun+1 ; af(nfun)=0.508 ; bf(nfun)=-180. ; EXIT ! -OC(=O)R
            ELSE
              nfun=nfun+1 ; af(nfun)=0.249 ; bf(nfun)=-515. ; EXIT ! -O- assume identical to F(-OH)            
            ENDIF
          ENDIF
        ENDDO
      ENDIF
    ENDDO
  ENDDO factloop

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
  ENDIF

  IF (j3/=0) THEN
    ax=ax**0.5 ; bx=bx/2.
  ENDIF

! add correction due to arrhenius coefficients
  arrhc(1)=arrhc(1)*ax
  arrhc(3)=arrhc(3)+bx

  IF (wrtsarinfo) THEN                                        ! debug 
! chain lenght factors ! need rewriting, only for substituent
! but can be optional (see Jenkin et al)
!      arrhc(1)=arrhc(1)*(1.+0.14*(1.-exp(-0.35*(nca-1))))
!      WRITE(saru,*)'fac longueur chaine=',(1.+0.14*(1-exp(-0.35*(nca-1))))
    DO j=1,nfun
      IF (wrtsarinfo) THEN                                    ! debug
        WRITE(saru,*) 'A(X',j,')=',af(j) ; WRITE(saru,*) 'B(X',j,')=',bf(j)                    
        WRITE(saru,*) 'F(X',j,')298=',af(j)*EXP(-bf(j)/298.)
      ENDIF
    ENDDO
    WRITE(saru,*) 'arrhc(1)',arrhc(1)                          ! debug
    WRITE(saru,*) 'arrhc(3)',arrhc(3)                          ! debug
    WRITE(saru,*) 'arrhc(298)',arrhc(1)*EXP(-arrhc(3)/298.)    ! debug
  ENDIF

END SUBROUTINE hoadd_rate

!=======================================================================
! PURPOSE: compute the reaction rate for OH addition for the case 1, 
! based on a "old" version of gecko (before the updating using the
! Jenkin, 2018 SAR).
!
!                                                                  
! The method assumes that addition to a given Cd depends only on   
! the structure (stability) of the addition products
!               
! FOR NON CONJUGATED ALKENE:
! The overall rate constant for OH addition is based on the SAR by 
! Nishino et al.: Rate constants for the gas-phase reactions of OH 
! radicals with a series of C6–C14alkenesat 299±2 K, 
! J. Phys. Chem. A, 113, 852–857, 2009.
! Rate depend on carbon chain length. Branching ratio are computed based
! Ziemann data when noted (see paper by La et al., ACP, 2016).
!
! FOR CONJUGATED ALKENE
! The method used is the SAR of Peeters et al., 1997, "kinetic studies 
! of reactions of alkylperoxy and haloalkylperoxy radicals with NO.
! A structure/reactivity relationship for reactions of OH with     
! alkenes and polyalkenes", in Chemical processes in atmospheric   
! oxidation, Eurotrac, G. Le Bras (edt).
! The method was also published in: Peeters, J., Boullart, W., 
! Pultau, V., Vandenberk, S., and Vereecken,L.: Structure-activity 
! relationship for the addition of OH to (poly)alkenes: site-specific 
! and total rate constants, J. Phys.Chem. A, 111, 1618–1631, 2007
!=======================================================================
SUBROUTINE hokadd_zie(ncdcase,group,j0,j1,j3,arrhc,nref,com)
  USE toolbox, ONLY: stoperr,addref
  IMPLICIT NONE

  CHARACTER(LEN=*),INTENT(IN) :: group(:)
  INTEGER,INTENT(IN) :: ncdcase
  INTEGER,INTENT(IN) :: j0,j1,j3
  REAL,INTENT(OUT)   :: arrhc(:)
  INTEGER,INTENT(INOUT) :: nref
  CHARACTER(LEN=*), INTENT(INOUT) :: com(:)
  
  REAL :: kadd
  INTEGER :: nca,ngr,i,xc
  CHARACTER(LEN=11),PARAMETER :: progname='hokadd_zie '
  CHARACTER(LEN=70)           :: mesg
  
  nca=0 ; ngr=0
  DO i=1,SIZE(group)
    IF (group(i)(1:1)=='C') nca=nca+1
  ENDDO

! CASE 1
! --------
  IF (ncdcase==1) THEN
    IF (j3==0) THEN   ! non conjugated Cd
      IF (group(j0)(1:4)=='CdH2') THEN 
        IF (group(j1)(1:4)== 'CdH2') THEN
          kadd=8.52E-12                            !Atkinson 2003, ktot(ethene)
          arrhc(1)=kadd*0.5
          CALL addref(progname,'OHETHENE',nref,com)
        ELSE IF (group(j1)(1:3)=='CdH') THEN 
          xc=nca-3
          kadd=2.8E-11+0.9E-11*(1.-exp(-0.35*xc))  !Nishino
          arrhc(1)=kadd*0.70 !Ziemann
          CALL addref(progname,'SARNISHINO',nref,com)
          CALL addref(progname,'BRCHZIE1',nref,com)
        ELSE IF (group(j1)(1:2)=='Cd') THEN
          xc=nca-4
          kadd=5.1E-11+1.6E-11*(1.-exp(-0.35*xc))  !Nishino
          arrhc(1)=kadd*0.70 !Ziemann
          CALL addref(progname,'SARNISHINO',nref,com)
          CALL addref(progname,'BRCHZIE1',nref,com)
        ENDIF
      ELSE IF(group(j0)(1:3)=='CdH') THEN
        IF (group(j1)(1:4)== 'CdH2') THEN
          xc=nca-3
          kadd=2.8E-11+0.9E-11*(1.-exp(-0.35*xc))  !Nishino
          arrhc(1)=kadd*0.30 !Ziemann
          CALL addref(progname,'SARNISHINO',nref,com)
          CALL addref(progname,'BRCHZIE1',nref,com)
        ELSE IF (group(j1)(1:3)=='CdH') THEN
          xc=nca-4
          kadd=6.3E-11+0.6E-11*(1.-exp(-0.35*xc)) !Nishino
          arrhc(1)=kadd *0.5 !Ziemann
          CALL addref(progname,'SARNISHINO',nref,com)
          CALL addref(progname,'BRCHZIE1',nref,com)
        ELSE IF (group(j1)(1:2)=='Cd') THEN
          xc=nca-4
          kadd=8.69E-11+0.1E-11*(1.-exp(-0.35*xc)) !3% de diff
          arrhc(1)=kadd *0.647              !meme rapport que Peeters
          CALL addref(progname,'SARNISHINO',nref,com)
          CALL addref(progname,'BRCHPEET1',nref,com)
        ENDIF
      ELSE IF (group(j0)(1:2)=='Cd') THEN   
        IF (group(j1)(1:4)=='CdH2') THEN
          xc=nca-4
          kadd=5.1E-11+1.6E-11*(1.-exp(-0.35*xc))
          arrhc(1)=kadd*0.30
          CALL addref(progname,'SARNISHINO',nref,com)
          CALL addref(progname,'BRCHZIE1',nref,com)
        ELSE IF (group(j1)(1:3)=='CdH') THEN
          xc=nca-4
          kadd=8.69E-11+0.1E-11*(1.-exp(-0.35*xc)) !3% diff
          arrhc(1)=kadd *0.353                !meme rapport que Peeters
          CALL addref(progname,'SARNISHINO',nref,com)
          CALL addref(progname,'BRCHPEET1',nref,com)
        ELSE IF (group(j1)(1:2)=='Cd') THEN
          xc=nca-4
          kadd=10.3E-11+1.0E-11*(1.-exp(-0.35*xc))!2% diff
          arrhc(1)=kadd *0.5
          CALL addref(progname,'SARNISHINO',nref,com)
          CALL addref(progname,'BRCHZIE1',nref,com)
        ENDIF
      ENDIF
    ENDIF

    IF (j3/=0) THEN ! conjugated Cd
      CALL addref(progname,'BRCHPEET2',nref,com)
      IF (group(j1)(1:3)=='CdH') THEN
        IF      (group(j3)(1:4)=='CdH2') THEN ; arrhc(1)=3.00E-11
        ELSE IF (group(j3)(1:3)=='CdH' ) THEN ; arrhc(1)=3.75E-11
        ELSE IF (group(j3)(1:2)=='Cd'  ) THEN ; arrhc(1)=5.05E-11
        ENDIF
      ELSE IF (group(j1)(1:2)=='Cd') THEN
        IF      (group(j3)(1:4)=='CdH2') THEN ; arrhc(1)=5.65E-11
        ELSE IF (group(j3)(1:3)=='CdH' ) THEN ; arrhc(1)=8.35E-11
        ELSE IF (group(j3)(1:2)=='Cd'  ) THEN ; arrhc(1)=9.85E-11
        ENDIF
      ENDIF
    ENDIF
  ELSE
    mesg="cdcase unknown in hoadd_rate_first"
    PRINT*, TRIM(mesg)
    PRINT*, "--error-- in ", TRIM(progname)
    STOP "in hokadd_zie"
  ENDIF  
  
END SUBROUTINE hokadd_zie

!=======================================================================
! PURPOSE: compute the reaction rate for OH addition for the case 2.
!
! The method used is the SAR of Kwok and Atkinson, Atmos. Environ., 
! 1685-1695, 1995. 
!
! Note : addition of OH is expected to occur mainly at the beta      
! position (with respect to the -CO group), i.e leading to the       
! >C(OH)-CR(.)-CO- structure. A yield of 80 % is set for this pathway.
!=======================================================================
SUBROUTINE hokaddc2(bond,group,cdtable,cdsub,ci,cf,arrhc)
  USE toolbox, ONLY: stoperr
  IMPLICIT NONE

  CHARACTER(LEN=*),INTENT(IN) :: group(:)
  INTEGER,INTENT(IN) :: bond(:,:)
  INTEGER,INTENT(IN) :: cdtable(:)
  INTEGER,INTENT(IN) :: cdsub(:)  
  INTEGER,INTENT(IN) :: ci,cf
  REAL,INTENT(OUT) :: arrhc(:)
  

  INTEGER :: i,l,posi,posf
  REAL    :: kstruct,fac
  
  kstruct=0. ; posi=0 ; posf=0 ; arrhc=0.

! Given Ci and Cf
  DO i=1,4
    IF (Ci==cdtable(i)) posi=i
    IF (Cf==cdtable(i)) posf=i
  ENDDO
  
! set the rate constant for the >C=C< structure
! ----------------------------------------------

! k struct for >C=CR-CO- or >C=CH-CO- structures
  IF (cdsub(posf)==2) THEN
    IF (cdsub(posi)==1) THEN  ;  kstruct=86.9E-12
    ELSE                      ;  kstruct=110E-12
    ENDIF

! kstruct for -CH=CR-CO- or -CH=CH-CO- structures        
  ELSE IF (cdsub(posf)==1) THEN
    IF (cdsub(posi)==1) THEN  ;   kstruct=60.2E-12
    ELSE                      ;   kstruct=86.9E-12
    ENDIF

! kstruct for CH2=CR-CO- or CH2=CH-CO- structures                
  ELSE IF (cdsub(posf)==0) THEN
    IF (cdsub(posi)==1) THEN  ;   kstruct=26.3E-12
    ELSE                      ;   kstruct=51.4E-12
    ENDIF
  ENDIF

! modify rate constant as a function of the substituent
! -----------------------------------------------------
  fac=1.
  DO l=1, SIZE(group)

! find group on Cf ('CO' and 'CO(O' do not have the same weight)
    IF ((bond(Cf,l)==1)) THEN
      IF (group(l)(1:3)=='CHO') fac=fac*0.34
      IF (group(l)(1:3)=='CO ') fac=fac*0.90                        
      IF (INDEX(group(l),'CH2(ONO2)')/=0) fac=fac*0.47
      IF (INDEX(group(l),'CH(ONO2)')/=0)  fac=fac*0.47
      IF (INDEX(group(l),'C(ONO2)')/=0)   fac=fac*0.47
      IF (INDEX(group(l),'CO(O')/=0)      fac=fac*0.25
      IF (INDEX(group(l),'CH2(OH)')/=0)   fac=fac*1.6
      IF (INDEX(group(l),'CH(OH)')/=0)    fac=fac*1.6
      IF (INDEX(group(l),'C(OH)')/=0)     fac=fac*1.6
    ENDIF 
          
! find group on Ci                
    IF ((bond(Ci,l)==1)) THEN
      IF (group(l)(1:3)=='CHO') fac=fac*0.34
      IF (group(l)(1:3)=='CO ') fac=fac*0.90
      IF (INDEX(group(l),'CH2(ONO2)')/=0) fac=fac*0.47
      IF (INDEX(group(l),'CH(ONO2)')/=0)  fac=fac*0.47
      IF (INDEX(group(l),'C(ONO2)')/=0)   fac=fac*0.47
      IF (INDEX(group(l),'CO(O')/=0)      fac=fac*0.25
      IF (INDEX(group(l),'CH2(OH)')/=0)   fac=fac*1.6
      IF (INDEX(group(l),'CH(OH)')/=0)    fac=fac*1.6
      IF (INDEX(group(l),'C(OH)')/=0)     fac=fac*1.6 
    ENDIF
  ENDDO
  arrhc(1)=kstruct*fac
END SUBROUTINE hokaddc2

!=======================================================================
! PURPOSE: compute the reaction rate for OH addition for the case 3.
!
! The method used is the SAR of Kwok and Atkinson, Atmos. Environ., 
! 1685-1695, 1995.
!                                                 
! Assign rate constant for OH addition to the double bond. The -C=C-C=C-
! structure is considered as a single structural unit. The rate 
! constant is a function of the number of substitutents, corrected by
! the group factor (assumed to be 1 for alkyl). Total number of 
! substituents is 6 (2 at both end+2 in the middle). 
! Note : The terminal positions of this structure must be position 1 
! and 4 in cdtable.
!=======================================================================
SUBROUTINE hokaddc3(bond,group,cdtable,cdsub,arrhc)
  USE toolbox, ONLY: stoperr
  IMPLICIT NONE

  CHARACTER(LEN=*),INTENT(IN) :: group(:)
  INTEGER,INTENT(IN) :: bond(:,:)
  INTEGER,INTENT(IN) :: cdtable(:)
  INTEGER,INTENT(IN) :: cdsub(:)  
  REAL,INTENT(OUT) :: arrhc(:)

  INTEGER :: i,j,nsub
  REAL    :: kstruct,fac
  CHARACTER(LEN=9),PARAMETER :: progname='hokaddc3'
  CHARACTER(LEN=70)          :: mesg
  
  kstruct=0. ; arrhc=0.

! count the number of substituents (remove 2 due to C-C bond in the 
! -C=C-C=C- structure)
  nsub=0
  DO i=1,4
    nsub=nsub+cdsub(i)
  ENDDO
  nsub=nsub-2

  IF (nsub==2) THEN      ; kstruct=142E-12
  ELSE IF (nsub==3) THEN ; kstruct=190E-12
  ELSE IF (nsub>3) THEN  ; kstruct=260E-12
  ELSE
    mesg="At least 2 substitents expected"
    PRINT*, TRIM(mesg)
    PRINT*, "--error-- in hokaddc3"
    STOP "in hokaddc3"
  ENDIF

! take group contributions into account
  fac=1.
  DO i=1,4
    DO j=1,SIZE(group)
      IF (bond(cdtable(i),j)==1) THEN
        IF (group(j)(1:3)=='CHO') fac=fac*0.34
        IF (group(j)(1:3)=='CO ') fac=fac*0.90 
        IF (INDEX(group(j),'CO(O')/=0)      fac=fac*0.25 
        IF (INDEX(group(j),'CH2(ONO2)')/=0) fac=fac*0.47
        IF (INDEX(group(j),'CH(ONO2)')/=0)  fac=fac*0.47
        IF (INDEX(group(j),'C(ONO2)')/=0)   fac=fac*0.47
        IF (INDEX(group(j),'CH2(OH)')/=0)   fac=fac*1.6 
        IF (INDEX(group(j),'CH(OH)')/=0)    fac=fac*1.6 
        IF (INDEX(group(j),'C(OH)')/=0)     fac=fac*1.6
      ENDIF
    ENDDO
  ENDDO
  arrhc(1)=kstruct*fac
END SUBROUTINE hokaddc3


END MODULE hoaddtool
