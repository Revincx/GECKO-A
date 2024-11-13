MODULE hochem
IMPLICIT NONE
CONTAINS
!=======================================================================                                                                      
! PURPOSE: Set the VOC+OH mechanism. Read first the database for known 
! reaction rates and mechanisms. If the mechanism is not provided in the
! database, then generate the VOC + OH reactions. All selected reactions
! are written in the output mechanism file.
!=======================================================================                                                                      
SUBROUTINE ho_voc(idnam,chem,bond,group,nring,brch,cut_off)
  USE keyparameter, ONLY: mxpd,mxnr,mecu,kohu,waru,ohu,mcou,mxcopd
  USE references, ONLY:mxlcod
  USE keyflag, ONLY: screenfg,wrtopeinfo,sar_only_fg,yldcut,TK
  USE toolbox, ONLY: addrx  
  USE searching, ONLY: srh5
  USE normchem, ONLY: stdchm
  USE database,ONLY:nkohdb,kohdb_chem,kohdb_arr,kohdb_com, &           ! k_OH database
      nkwoh,nkwoh_pd,kwoh_rct,kwoh_pd,kwoh_copd,kwoh_yld,kwoh_com      ! VOC+OH mechanism     
  USE dictstacktool, ONLY: bratio
  USE rxwrttool, ONLY:rxwrit,rxinit
  USE radchktool, ONLY: radchk
  USE toolbox, ONLY: stoperr,add1tonp,addref

  IMPLICIT NONE

  CHARACTER(LEN=*),INTENT(IN) :: idnam    ! name of current species
  CHARACTER(LEN=*),INTENT(IN) :: chem     ! formula of current species
  CHARACTER(LEN=*),INTENT(IN) :: group(:) ! group (without ring joining char)
  INTEGER,INTENT(IN) :: bond(:,:)     ! bond matrix
  INTEGER,INTENT(IN) :: nring         ! # of rings
  REAL,INTENT(IN)    :: brch          ! max yield of the input species
  REAL,INTENT(IN)    :: cut_off       ! ratio threshold below which rx is ignored
      
  CHARACTER(LEN=LEN(group(1))) :: tgroup(SIZE(group))
  INTEGER         :: tbond(SIZE(bond,1),SIZE(bond,2))
  CHARACTER(LEN=LEN(chem)) :: pchema(mxnr)

  INTEGER              :: flag(mxnr)
  CHARACTER(LEN=LEN(chem)) :: pchem(mxnr)
  CHARACTER(LEN=LEN(idnam)) :: coprod(mxnr,mxcopd)
  INTEGER              :: flag_del(mxnr)
  CHARACTER(LEN=LEN(chem)) :: pchem_del(mxnr)
  CHARACTER(LEN=LEN(idnam)) :: coprod_del(mxnr,mxcopd)
  REAL                 :: tarrhc(mxnr,3)
  
  REAL    :: smk298, k298(mxnr)
  REAL    :: smkT, kT(mxnr)
  REAL    :: ratio(mxnr)
  REAL    :: allakT, allak298, kmax
  INTEGER :: sumflg,npath
  REAL    :: brtio

  REAL    :: rk298, exparrhc(3)  ! known rate constant
  REAL    :: rkT  ! kT calculated with known rate constant coeff coeffs
  INTEGER :: nohprod

  CHARACTER(LEN=LEN(idnam)) :: p(mxpd), r(3)
  REAL    :: s(mxpd), arrh(3)
  INTEGER :: idreac, nlabel
  REAL    :: xlabel,folow(3),fotroe(4)
  
  CHARACTER(LEN=LEN(idnam)) :: pname(mxnr)
  CHARACTER(LEN=LEN(idnam)) :: pope
  REAL    :: a1,nt2,ae3, sope

  INTEGER :: nr,i,j,k,l,np,iloc,nca
  REAL    :: kscale

  INTEGER,PARAMETER  :: mxrpd=2  ! max # of products returned by radchk sub
  CHARACTER(LEN=LEN(chem)) :: rdckpd(mxrpd)
  CHARACTER(LEN=LEN(idnam)) :: rdckcopd(mxrpd,mxcopd)
  REAL    :: sc(mxrpd)
  REAL    :: sc_del(mxnr,mxrpd)
  INTEGER :: nip

  LOGICAL :: kdbflg,rxdbflg
  
  CHARACTER(LEN=6),PARAMETER :: progname='ho_voc'
  CHARACTER(LEN=70)          :: mesg
  
  INTEGER,PARAMETER     :: mxcom=10
  INTEGER               :: nrxref(mxnr)
  CHARACTER(LEN=mxlcod) :: rxref(mxnr,mxcom)
  CHARACTER(LEN=mxlcod) :: comkdb(SIZE(kohdb_com,2))  ! reference for koh in db

! for available exp. rate constant, give priority to SAR Ea (default) or
! to the experimental Ea (old version - not recommended) 
  INTEGER,PARAMETER :: expeaflg=0 ! set to 0 for default case

  IF (screenfg) WRITE(6,*)progname

! initialize
! -----------
  nr=0 ; rk298=0.
  pchem(:)=' ' ; pchema(:)=' '
  flag(:)=0 ; k298(:)=0. ; ratio(:)=0.
  coprod(:,:)=' ' ; tarrhc(:,:)=0. ; exparrhc(:)=0. 
  coprod_del(:,:)=' ' ; pchem_del(:)=' ' ; flag_del(:) =0 
  sc_del(:,:)=0. ;  pname(:)=' '
  nrxref(:)=0 ; rxref(:,:)=' ' ; comkdb(:)=' '

  kdbflg=.FALSE. ; rxdbflg=.FALSE.
  tgroup(:)=group(:) ; tbond(:,:)=bond(:,:)
      
  IF (LEN_TRIM(chem)<1) RETURN
  nca=COUNT(tgroup/=' ')

! ------------------------------------------------------------
! check if the species is known in the rate constant data base
! ------------------------------------------------------------
  iloc=srh5(chem,kohdb_chem,nkohdb)
  IF (iloc>0) THEN
    rk298=kohdb_arr(iloc,1)*298**kohdb_arr(iloc,2)*EXP(-kohdb_arr(iloc,3)/298.)
    rkT  =kohdb_arr(iloc,1)*TK**kohdb_arr(iloc,2) *EXP(-kohdb_arr(iloc,3)/TK)
    exparrhc(:)=kohdb_arr(iloc,:)
    comkdb(:)=kohdb_com(iloc,:)
    kdbflg=.TRUE.
  ENDIF

! -----------------------------------------------------------
! Check if rate or products are provided in the database. 
! Note: yields below threshold are not automatically eliminated.
! -----------------------------------------------------------
! SAR assessment - ignore the rate in the database
  IF (sar_only_fg) kdbflg=.FALSE.
    
  IF (kdbflg) THEN
    iloc=srh5(chem,kwoh_rct,nkwoh)
    IF (iloc>0) THEN
      rxdbflg=.TRUE.  ;  nohprod=nkwoh_pd(iloc)

! write reaction rate to kOH file 
      WRITE(kohu,'(a,2x,1pe10.2,1pe10.2,a)') idnam, rk298, rkT, " DATA, DATA"
      
! loop through known products
      DO j=1,nohprod
        CALL addrx(progname,chem,nr,flag)
        pchema(j)=kwoh_pd(iloc,j)
        ratio(j)=kwoh_yld(iloc,j)

! check the products
        IF (INDEX(pchema(j),'.')/=0) THEN
          CALL radchk(pchema(j),rdckpd,rdckcopd,nip,sc,nrxref(nr),rxref(nr,:))
          pchem(j)=rdckpd(1)
          coprod(j,:)=rdckcopd(1,:)   
          IF (nip==2) THEN
            flag_del(nr)=1
            pchem_del(nr)=rdckpd(2)
            sc_del(nr,1)=sc(1) ; sc_del(nr,2)=sc(2)
            coprod_del(nr,:)=rdckcopd(2,:)
            CALL stdchm(pchem_del(nr))
          ENDIF
        ELSE
          pchem(j)= pchema(j)
        ENDIF
        CALL stdchm(pchem(j))

! add known coproduct(s) if necessary
        kloop: DO k=1,SIZE(kwoh_copd,3)  ! loop possible coproducts
          IF (kwoh_copd(iloc,j,k)/=' ') THEN           
            DO l=1,mxcopd       
              IF (coprod(j,l)/=' ') CYCLE
              coprod(j,l)=kwoh_copd(iloc,j,k)
              CYCLE kloop
            ENDDO 
          ELSE
            EXIT kloop
          ENDIF           ! coproduct is non-zero
        ENDDO kloop       ! coproducts loop
        
! add reaction rate
        tarrhc(j,1)=exparrhc(1)*ratio(j) ! keep exp. value
        tarrhc(j,2)=exparrhc(2)          ! keep exp. value 
        tarrhc(j,3)=exparrhc(3)          ! keep exp. Ea 
        CALL addref(progname,'RXEXP',nrxref(nr),rxref(nr,:),chem)
        DO i=1,SIZE(comkdb)
          IF (comkdb(i)/=' ') CALL addref(progname,comkdb(i),nrxref(nr),rxref(nr,:),chem)
        ENDDO
! add references to the reaction
        DO i=1, SIZE(kwoh_com,2)
          IF (kwoh_com(iloc,i)/=' ') &
              CALL addref(progname,kwoh_com(iloc,i),nrxref(nr),rxref(nr,:),chem)
        ENDDO 
      ENDDO

    ENDIF  ! reactant has known products

  ENDIF

! =========================
! MAKE THE VOC+OH MECHANISM
! =========================

  IF (.NOT.rxdbflg) THEN

! perform all OH+VOC reactions
! ----------------------------
    CALL mkhorx(chem,bond,group,nca,nring, &
         nr,flag,pchem,coprod,flag_del,pchem_del,coprod_del,sc_del,tarrhc,&
         nrxref,rxref)
  
! select pathways
! ---------------

! compute the rate at 298 K and at requested T (might not be 298 K)
    smkT=0. ; smk298=0. ; k298(:)=0.; kT=0.
    DO i=1,nr
      IF (flag(i)==0) CYCLE
      k298(i)=tarrhc(i,1) * 298.**tarrhc(i,2)*EXP(-tarrhc(i,3)/298.) 
      k298(i)=AMIN1(k298(i),2.0E-10)        ! do not exceed collision rate
      smk298 =smk298+k298(i)
      kT(i)  =tarrhc(i,1) * TK**tarrhc(i,2)  *EXP(-tarrhc(i,3)/TK)
      kT(i)=AMIN1(kT(i),2.0E-10)        ! do not exceed collision rate
      smkT   =smkT+kT(i)
    ENDDO

! write reaction rate for SAR assessment
    IF (sar_only_fg) WRITE(ohu,'(1PE12.2,2x,a)') smk298, TRIM(chem)

! if no reactions, warning message (expected for a few species only)
    IF (smkT==0) THEN
      WRITE(waru,'(a)') '--Warning-- No OH reaction found for: ', TRIM(chem)
      RETURN
    ENDIF
    
! remove pathways below cutoff, at requested mechanism temperature TK
    WHERE (kT(:) < cut_off*smkT) flag(:)=0
    WHERE (kT(:) < (yldcut/brch)*smkT) flag(:)=0
    IF (SUM(flag)==0) THEN
      kmax=MAXVAL(kT)
      WHERE (kT(:) > kmax/2.) flag(:)=1
      sumflg=SUM(flag)
      IF (sumflg==0) THEN
        mesg="no path remains after cut_off"
        CALL stoperr(progname,mesg,chem)
      ENDIF
    ENDIF

! sum of all active rate constants (allak) and set branching ratio
    allakT=0. ; allak298=0.
    DO i=1,nr
      IF (flag(i)==1) allakT  =allakT  +kT(i)
      IF (flag(i)==1) allak298=allak298+k298(i)
    ENDDO

    ratio(:)=0.
    npath=0
    DO i=1,nr 
      IF (flag(i)==1) ratio(i)=kT(i)/allakT
      IF ((flag(i)==1).AND.(ratio(i)/=0.)) npath=npath+1 
    ENDDO

! SET THE RATE OF EACH CHANNEL

! If rate is not in the database, rescale tarrhc(i,1) to keep smk298
    IF (.NOT.kdbflg) THEN            
      kscale=smk298/allak298          ! scaling factor for sum active k = smk298
      DO i=1,nr
        IF (flag(i)==1) tarrhc(i,1)=tarrhc(i,1)*kscale
      ENDDO

! Write generated k298 & kT to the kivoci dataset
      WRITE(kohu,'(a,2x,2(1pe10.2),a)') idnam, smk298,smkT," SAR, SAR"

! If rate is in the database then rescale to experimental value needed.
    ELSE


      IF (expeaflg==0) THEN   ! use exp. value only at 298K (DEFAULT)

        IF ((npath==1).AND.(exparrhc(3)/=0.)) THEN 
        ! if only one pathway and exp Ea exist, keep experimental values
        ! Write reference k298 & kT to the kivoci output dataset 
          WRITE(kohu,'(a,2x,2(1pe10.2),a)') idnam, rk298,rkT," DATA, DATA"
          DO i=1,nr
            IF (ratio(i)/=0.) THEN
              tarrhc(i,1)=exparrhc(1)*ratio(i) ! keep exp. value
              tarrhc(i,2)=exparrhc(2)          ! keep exp. value 
              tarrhc(i,3)=exparrhc(3)          ! keep exp. Ea 
              CALL addref(progname,'KEXPBRSAR',nrxref(i),rxref(i,:),chem)
              DO j=1,SIZE(comkdb)    ! add ref for the exp. data
                IF (comkdb(j)/=' ') CALL addref(progname,comkdb(j),nrxref(i),rxref(i,:),chem)
              ENDDO
              EXIT
            ENDIF
          ENDDO
        ELSE
          
        ! Write reference k298 & calculated kT to the kivoci output dataset 
        WRITE(kohu,'(a,2x,2(1pe10.2),a)') idnam, rk298,smkT," DATA, SAR"

        kscale=rk298/allak298     ! scaling factor for SAR rate = exp. value @298K
        DO i=1,nr
          IF (flag(i)==1) THEN
            tarrhc(i,1)=tarrhc(i,1)*kscale ! keep SAR values
            CALL addref(progname,'KEXP298',nrxref(i),rxref(i,:),chem)
            DO j=1,SIZE(comkdb)    ! add ref for the exp. data
              IF (comkdb(j)/=' ') CALL addref(progname,comkdb(j),nrxref(i),rxref(i,:),chem)
            ENDDO
          ENDIF
        ENDDO
        ENDIF

      ELSE  ! full priority to exp. values (NOT RECOMMENDED)

        ! Write reference kT & k298 to the kivoci dataset
        WRITE(kohu,'(a,2x,1pe10.2,1pe10.2,a)') idnam, rk298, rkT, " DATA, DATA"

        DO i=1,nr
          IF (flag(i)==1) THEN
            tarrhc(i,1)=exparrhc(1)*ratio(i) ! keep exp. value, scaled for each active channel
            tarrhc(i,2)=exparrhc(2)          ! keep exp. value 
            tarrhc(i,3)=exparrhc(3)          ! keep exp. Ea 
            CALL addref(progname,'KEXPBRSAR',nrxref(i),rxref(i,:),chem)
            DO j=1,SIZE(comkdb)    ! add ref for the exp. data
              IF (comkdb(j)/=' ') CALL addref(progname,comkdb(j),nrxref(i),rxref(i,:),chem)
            ENDDO
          ENDIF
        ENDDO
      ENDIF
    ENDIF
  ENDIF ! end of generated reactions
    
! =================================
! WRITE REACTION
! =================================

  DO i=1,nr
    IF (flag(i)==0) CYCLE

! initialize reaction
    CALL rxinit(r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe)
    arrh(:)=tarrhc(i,:)
    r(1)=idnam  ; r(2)='HO   '
    s(1)=1.
    
    brtio=brch*ratio(i)
    CALL bratio(pchem(i),brtio,p(1),nrxref(nr),rxref(nr,:))
    pname(i)=p(1)          ! store for later use in operator information

! second product (if delocalisation allowed)
    np=1
    IF (flag_del(i)/=0) THEN
      np=np+1
      s(1)=sc_del(i,1) ; s(np)=sc_del(i,2)  ! overwrite the 1st coef. for s
      CALL bratio(pchem_del(i),brtio,p(np),nrxref(nr),rxref(nr,:))
    ENDIF

! write out coproducts
    DO j=1,mxcopd
      IF (coprod(i,j)(1:1)/=' ') THEN
        CALL add1tonp(progname,chem,np)  
        s(np)=1. ; p(np)=coprod(i,j)
      ENDIF
    ENDDO
    IF (flag_del(i)/=0) THEN
      DO j=1,mxcopd
        IF (coprod_del(i,j)(1:1)/=' ') THEN
          CALL add1tonp(progname,chem,np)  
          s(np)=sc_del(i,2) ; p(np)=coprod_del(i,j)
        ENDIF
      ENDDO
    ENDIF

! write out: reactions with OH are thermal reactions (idreac=0) and do not require labels
    CALL rxwrit(mecu,r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe,rxref(i,:))
  ENDDO

! =================================
! WRITE INFORMATION REQUIRED FOR OPERATOR
! =================================
  IF (wrtopeinfo) THEN
    nt2=0. ; ae3=0. ; a1=0. 
    opeloop2: DO i=1,nr
      IF (flag(i)/=0) THEN
        a1=a1+tarrhc(i,1)
        IF ((tarrhc(i,2)/=nt2).OR.(tarrhc(i,3)/=ae3)) THEN
          a1=allakT ; nt2=0. ; ae3=0.
          CYCLE opeloop2
        ELSE
          nt2=tarrhc(i,2) ; ae3=tarrhc(i,3)
        ENDIF 
      ENDIF
    ENDDO opeloop2
    WRITE(10,'(A15,A1,A6,A4,1X,ES10.3,1X,f4.1,1X,f7.0)') &
       '****INIT GHO +','G',idnam,'****',a1,nt2,ae3
      
    DO i=1,nr
      IF (flag(i)==0) CYCLE
      WRITE(10,'(f5.3,2X,A1,A6)') ratio(i), 'G',pname(i)

! write out coproducts
      DO j=1,mxcopd
        IF (coprod(i,j)(1:1)/=' ') THEN
          sope=1. ; pope=coprod(i,j)
          WRITE(10,'(f5.3,2X,A1,A6)') sope*ratio(i),'G',pope
        ENDIF
      ENDDO
      IF (flag_del(i)/=0) THEN
        DO j=1,mxcopd
          IF (coprod_del(i,j)(1:1)/=' ') THEN
            sope=sc_del(i,j) ; pope=coprod_del(i,j)
            WRITE(10,'(20X,A1,A6)') 'G',pope
          ENDIF
        ENDDO
      ENDIF
    ENDDO
    WRITE(10,*)'end',progname
  ENDIF

END SUBROUTINE ho_voc



!=======================================================================
!=======================================================================
! PURPOSE: Generate the VOC+OH mechanism. The routine identify: 
! (1) all possible H-abstraction (-CH, -OH, -OOH, -C(O)OH
! (2) all possible OH addition to C=C bond
! (3) all possible OH addition to aromatic rings
!=======================================================================
!=======================================================================
SUBROUTINE mkhorx(chem,bond,group,nca,nring, &
           nr,flag,pchem,coprod,flag_del,pchem_del,coprod_del,sc_del,tarrhc, &
           nrxref,rxref)
  USE references, ONLY:mxlcod
  USE toolbox, ONLY: addrx,stoperr,addref  
  USE reactool, ONLY: swap,rebond
  USE radchktool, ONLY: radchk
  USE normchem, ONLY: stdchm
  USE hotool, ONLY: rabsoh
  USE cdtool, ONLY: cdcase2
  USE hoaddtool, ONLY: hoadd_c1,hoadd_c2,hoadd_c3,hoadd_c4, hoadd_c5,hoadd_c6,hoadd_c7
  USE ho_aromtool, ONLY: hoadd_arom
  USE keyparameter, ONLY: mxtrk,mxlcd
  IMPLICIT NONE

  CHARACTER(LEN=*),INTENT(IN) :: chem ! formula of the VOC considered
  CHARACTER(LEN=*),INTENT(IN) :: group(:) ! group (without ring joining char)
  INTEGER,INTENT(IN) :: bond(:,:)     ! bond matrix
  INTEGER,INTENT(IN) :: nca           ! # number of nodes
  INTEGER,INTENT(IN) :: nring         ! # of rings
  INTEGER,INTENT(OUT):: nr            ! # of reaction channel
  INTEGER,INTENT(OUT):: flag(:)       ! flag for active channel
  CHARACTER(LEN=*),INTENT(OUT) :: pchem(:)    ! "main" product of reaction channel i
  CHARACTER(LEN=*),INTENT(OUT) :: coprod(:,:) ! coproducts (short names) for reaction channel i
  INTEGER,INTENT(OUT)          :: flag_del(:) ! flag for a 2nd channel of reaction i (delocalisation)
  CHARACTER(LEN=*),INTENT(OUT) :: pchem_del(:)    ! "main" product for the 2nd channel
  CHARACTER(LEN=*),INTENT(OUT) :: coprod_del(:,:) ! coproducts for 2nd channel
  REAL,INTENT(OUT)             :: sc_del(:,:)     ! ratio od the 2 channels (for delocalisation)
  REAL,INTENT(OUT)             :: tarrhc(:,:)     ! arrhenius coefficient for reaction i
  INTEGER,INTENT(OUT)          :: nrxref(:)       ! # of references added for each reaction
  CHARACTER(LEN=*),INTENT(OUT) :: rxref(:,:)      ! references for the reaction (1st index is rx #)

  CHARACTER(LEN=LEN(group(1))) :: tgroup(SIZE(group))
  CHARACTER(LEN=LEN(group(1))) :: pold,pnew
  INTEGER         :: tbond(SIZE(bond,1),SIZE(bond,2))
  CHARACTER(LEN=LEN(chem)) :: tempkc

  INTEGER :: carb
  INTEGER :: cdtable(4),tcdtable(4),cdsub(4),ncdcase
  INTEGER :: cdcarbo(4,2),cdeth(4,2)

  INTEGER :: i,j,k,nc,tempring,rxnflg,ngr

  INTEGER,PARAMETER  :: mxrpd=2  ! max # of products returned by radchk sub
  CHARACTER(LEN=LEN(pchem(1))) :: rdckpd(mxrpd)
  CHARACTER(LEN=LEN(coprod(1,1))) :: rdckcopd(mxrpd,SIZE(coprod,2))
  REAL    :: sc(mxrpd)
  INTEGER :: nip

  INTEGER :: ncdtrack                ! # of Cd tracks 
  INTEGER :: cdtrack(mxtrk,mxlcd)    ! cdtrack(i,j) Cd nodes for the ith track         
  INTEGER :: xcdconjug(mxtrk)        ! flag for conjugated tracks (C=C-C=C)
  INTEGER :: xcdsub(SIZE(group))     ! # of nodes bonded to a Cd group
  INTEGER :: xcdeth(SIZE(group),2)   ! -OR group # bounded to a Cd (max 2 -OR groups)
  INTEGER :: xcdcarbo(SIZE(group),2) ! -CO- group # bounded to a Cd (max 2 CO groups)  
  INTEGER :: xcdcase(mxtrk)          ! "case" of the track

  CHARACTER(LEN=6),PARAMETER :: progname='mkhorx'
  CHARACTER(LEN=70)          :: mesg

  nr=0  ;  flag(:)=0  ;  tarrhc(:,:)=0.
  pchem(:)=' '     ; coprod(:,:)=' ' 
  pchem_del(:)=' ' ; coprod_del(:,:)=' '
  flag_del(:)=0  ; sc_del(:,:)=0. 

  tgroup(:)=group(:) ; tbond(:,:)=bond(:,:)
  carb=0 ; cdtable(:)=0 ; tcdtable(:)=0 ; cdsub(:)=0
  ncdcase=0 ; cdcarbo(:,:)=0 ; cdeth(:,:)=0
  rxnflg=1 ; nrxref(:)=0 ; rxref(:,:)=' '
  ngr=COUNT(group/=' ')
  
! =================================
!  H ABSTRACTION                              
! =================================

! H-ABSTRACTION FROM CH
! ---------------------

! find all C-H bonds - loop over each C
  DO i=1,nca
    IF (INDEX(group(i),'CH')==0) CYCLE
    IF (INDEX(group(i),'(NO2)')/=0) CYCLE         ! no H abst (but OH addition below)
    CALL addrx(progname,chem,nr,flag)

! compute associated rate constant and add references
    CALL rabsoh(tbond,tgroup,i,tarrhc(nr,:),nring)
    CALL addref(progname,'HABSAR',nrxref(nr),rxref(nr,:),chem)
    CALL addref(progname,'MJ18KMV000',nrxref(nr),rxref(nr,:),chem)

! remove one H atom
    IF (group(i)(1:3)=='CH3')      THEN ; pold='CH3' ; pnew='CH2'
    ELSE IF (group(i)(1:3)=='CH2') THEN ; pold='CH2' ; pnew='CH'
    ELSE IF (group(i)(1:3)=='CHO') THEN ; pold='CHO' ; pnew='CO'
    ELSE IF (group(i)(1:2)=='CH')  THEN ; pold='CH'  ; pnew='C'            
    ENDIF
    CALL swap(group(i),pold,tgroup(i),pnew)

! add radical dot at end of group
    nc=INDEX(tgroup(i),' ') ; tgroup(i)(nc:nc)='.'

! re-build molecule, including ring if present
    CALL rebond(tbond,tgroup,tempkc,tempring)

! check radical and find co-products:
    CALL radchk(tempkc,rdckpd,rdckcopd,nip,sc,nrxref(nr),rxref(nr,:))
    pchem(nr)=rdckpd(1) ; sc_del(nr,1)=sc(1)
    coprod(nr,:)=rdckcopd(1,:)
    CALL stdchm(pchem(nr))
    IF (nip==2) THEN
      flag_del(nr)=1
      pchem_del(nr)=rdckpd(2) ; sc_del(nr,2)=sc(2)
      coprod_del(nr,:)=rdckcopd(2,:)
      CALL stdchm(pchem_del(nr))
    ENDIF

! reset group: 
    tgroup(:)=group(:)

  ENDDO

! H-ABSTRACTION FROM ROH
! ----------------------

! find all hydroxy group (but not carboxylic acid)
  DO i=1,nca
    IF (INDEX(group(i),'(OH)')==0) CYCLE
    IF (INDEX(group(i),'CO(OH)')/=0) CYCLE
    CALL addrx(progname,chem,nr,flag)

! assign rate constant by Jenkin 2015, MAGNIFY
    IF (group(i)(1:5)=='c(OH)') THEN                                 ! phenol
      tarrhc(nr,1)=2.6E-12 ; tarrhc(nr,2)=0. ; tarrhc(nr,3)=0.  
      CALL addref(progname,'HABSAR1',nrxref(nr),rxref(nr,:),chem)
      CALL addref(progname,'MJ18KMV001',nrxref(nr),rxref(nr,:),chem)
      DO j=1,nca
        IF ((tbond(i,j)==1).AND.(group(j)(1:6)=='c(NO2)')) THEN      ! nitrophenol
          tarrhc(nr,1)=1.28E-12 ; tarrhc(nr,2)=0. ; tarrhc(nr,3)=660.
          CALL addref(progname,'HABSAR1',nrxref(nr),rxref(nr,:),chem)
          CALL addref(progname,'MJ18KMV001',nrxref(nr),rxref(nr,:),chem)
          EXIT
        ENDIF
      ENDDO
    ELSE
      tarrhc(nr,1)=1.28E-12 ; tarrhc(nr,2)=0. ; tarrhc(nr,3)=660.    ! aliphatic -OH
      CALL addref(progname,'HABSAR1',nrxref(nr),rxref(nr,:),chem)
      CALL addref(progname,'MJ18KMV000',nrxref(nr),rxref(nr,:),chem)
    ENDIF

! replace hydroxy group by alkoxy group
    pold='(OH)' ; pnew='(O.)'
    CALL swap(group(i),pold,tgroup(i),pnew)
      
! rebuild, check, and find coproducts & rename
    CALL rebond(tbond,tgroup,tempkc,tempring)
    CALL radchk(tempkc,rdckpd,rdckcopd,nip,sc,nrxref(nr),rxref(nr,:))
    pchem(nr)=rdckpd(1)
    IF (nip/=1) THEN
      mesg="unexpected 2 products for ROH+OH "
      CALL stoperr(progname,mesg,chem)
    ENDIF
    coprod(nr,:)=rdckcopd(1,:)
    CALL stdchm(pchem(nr))

! reset group: 
    tgroup(i)=group(i)
  ENDDO

! H ABSTRACTION FROM ROOH
! ------------------------

! find hydroperoxide groups (peracid is not treated separately!)
  DO i=1,nca
    IF (INDEX(group(i),'(OOH)')==0) CYCLE
    CALL addrx(progname,chem,nr,flag)

! assign rate constants for alkyl hydroperoxide - Mike Jenkin report for MAGNIFY
    tarrhc(nr,1)=3.68E-13 ; tarrhc(nr,2)=0. ; tarrhc(nr,3)=-635.
    CALL addref(progname,'HABSAR2',nrxref(nr),rxref(nr,:),chem)
    CALL addref(progname,'MJ18KMV000',nrxref(nr),rxref(nr,:),chem)

! change (OOH) to (OO.)
    pold='(OOH)' ; pnew='(OO.)'
    CALL swap(group(i),pold,tgroup(i),pnew)

! rebuild, check, assign coproducts & rename:
    CALL rebond(tbond,tgroup,tempkc,tempring)
    CALL radchk(tempkc,rdckpd,rdckcopd,nip,sc,nrxref(nr),rxref(nr,:))
    pchem(nr)=rdckpd(1)
    IF (nip/=1) THEN
      mesg="unexpected 2 products for ROOH+OH "
      CALL stoperr(progname,mesg,chem)
    ENDIF
    coprod(nr,:)=rdckcopd(1,:)
    CALL stdchm(pchem(nr))

! reset group: 
    tgroup(i)=group(i)
  ENDDO

! H ABSTRACTION from -COOH 
!-------------------------

! find carboxylic groups: 
  DO i=1,nca
    IF (INDEX(group(i),'CO(OH)')==0) CYCLE
    CALL addrx(progname,chem,nr,flag)

! RCOOH + OH -> RCOO + H2O 
    ratecooh: DO j=1,nca
      IF (tbond(i,j)/=0) THEN
        IF (tgroup(j)(1:3)=='CO ') THEN
          tarrhc(nr,1)=4.77E-14 ; tarrhc(nr,2)=0. ; tarrhc(nr,3)=-275.   
          CALL addref(progname,'HABSAR3',nrxref(nr),rxref(nr,:),chem)
          CALL addref(progname,'MJ18KMV000',nrxref(nr),rxref(nr,:),chem)
        ELSE
          tarrhc(nr,1)=2.87E-14 ; tarrhc(nr,2)=0. ; tarrhc(nr,3)=-880.             
          CALL addref(progname,'HABSAR3',nrxref(nr),rxref(nr,:),chem)
          CALL addref(progname,'MJ18KMV000',nrxref(nr),rxref(nr,:),chem)
        ENDIF
        EXIT ratecooh
      ENDIF
    ENDDO ratecooh

! change CO(OH) to CO(O.)
    pold='CO(OH)' ; pnew='CO(O.)'
    CALL swap(group(i),pold,tgroup(i),pnew)

! rebuild, check, and assign coproducts & rename
    CALL rebond(tbond,tgroup,tempkc,tempring)
    CALL radchk(tempkc,rdckpd,rdckcopd,nip,sc,nrxref(nr),rxref(nr,:))
    pchem(nr)=rdckpd(1)
    sc_del(nr,1)=sc(1)
    coprod(nr,:)=rdckcopd(1,:)
    CALL stdchm(pchem(nr))
    IF (nip==2) THEN
      flag_del(nr)=1
      pchem_del(nr)=rdckpd(2)
      sc_del(nr,2)=sc(2)
      coprod_del(nr,:)=rdckcopd(2,:)
    ENDIF
 
! reset group: 
    tgroup(i)=group(i)
  ENDDO

! H-ATOM ABSTRACTION FROM (OOOH)
!------------------------------------------------

! find hydrotrioxides groups: (Peroxy acid is not treated separately!)
  DO i=1,nca
    IF (INDEX(group(i),'(OOOH)') == 0) CYCLE
    CALL addrx(progname,chem,nr,flag)

! Rate from Jenkin et al., 2019 (rate from Anglada and Sole, 2018)
    tarrhc(nr,1)=1.46E-12 ; tarrhc(nr,2)=0. ; tarrhc(nr,3)=-1037.
    CALL addref(progname,'HABSAR5',nrxref(nr),rxref(nr,:),chem)
    CALL addref(progname,'JA18KMV000',nrxref(nr),rxref(nr,:),chem)
 
    pold = '(OOOH)' ; pnew = '(O.)'
    CALL swap(group(i),pold,tgroup(i),pnew)
    
    CALL rebond(tbond,tgroup,tempkc,tempring)
    CALL radchk(tempkc,rdckpd,rdckcopd,nip,sc,nrxref(nr),rxref(nr,:))
    pchem(nr) = rdckpd(1)
    sc_del(nr,1)=sc(1)
    coprod(nr,:)=rdckcopd(1,:)
    IF (nip/=1) WRITE(6,*) '2 products in hochem - H abs from OOOH'
    IF (nip/=1) STOP
    
    CALL stdchm(pchem(nr))
 
    tgroup(i) = group(i)
  ENDDO


            
! =================================
! OH ADDITION TO -NO2
! =================================

! RC(NO2) + OH -> RC(OH) + NO2 (reactions differs from Jenkin et al., 2018)
! based on personnal communication from MJ.
! find nitro groups: 
  DO i=1,nca
    IF (INDEX(group(i),'(NO2)')==0) CYCLE
    IF (INDEX(group(i),'c')/=0) CYCLE
    IF (INDEX(group(i),'Cd')/=0) CYCLE
    CALL addrx(progname,chem,nr,flag)

! rates from Jenkin, 2018 - but RC(NO2) + OH -> RC(OH) + NO2
    tarrhc(nr,1)=1.1E-13  ; tarrhc(nr,2)=0. ; tarrhc(nr,3)=0.             
    CALL addref(progname,'HABSAR4',nrxref(nr),rxref(nr,:),chem)
    CALL addref(progname,'ON89KMV000',nrxref(nr),rxref(nr,:),chem)
    CALL addref(progname,'MJ18KMV000',nrxref(nr),rxref(nr,:),chem)

! change >C(NO2)- to >C(OH)-
    pold='(NO2)' ; pnew='(OH)'
    CALL swap(group(i),pold,tgroup(i),pnew)

! rebuild, add NO2 as coproduct & rename
    CALL rebond(tbond,tgroup,tempkc,tempring)
    pchem(nr)=tempkc
    coprod(nr,1)='NO2'
    CALL stdchm(pchem(nr))

! reset group: 
    tgroup(:)=group(:)
  ENDDO

! =================================
!  OH ADDITION TO >C=C< BOND 
! =================================

! The parameterization used to assign a rate constant for OH addition to
! C=C bonds depends on whether the C=C bond is conjugated or not with a 
! C=O bond (i.e. structure of type -C=C-C=O). The subroutine return the 
! "case" the species belong. Five cases are considered:     
! CASE 1: regular Cd molecule: only >C=C< and >C=C-C=C< bonds in the 
!         molecule (i.e, without conjugated C=C-C=O)
! CASE 2: are for structure containing the >C=C-C=O structure but 
!         no C=C-C=C 
! CASE 3: are for the -CO-C=C-C=C-C=O structure only (i.e. containing
!         carbonyl at both side of the conjugated C=C-C=C)
! CASE 4: Two double bonds non conjugated (i.e. C=C-C-C=C) but one 
!         containing at least one C=C-C=O  
! CASE 5: are for the -CO-C=C-C=C< structure (i.e. containing carbonyl 
!         at only one side of the conjugated C=C-C=C 
! CASE 6: -C=C=O
! CASE 7: addition of OH to C=C bond of vinyl ether (-O-C=C)
!         and dihydrofurans

! find double-bonded carbons in the chain (label 600 => end C=C)
  IF (INDEX(chem,'Cd')/=0) THEN

! get the case for each C=C bonds in chem (xcdcase) 
    CALL cdcase2(chem,bond,group,rxnflg,ncdtrack,cdtrack,&
               xcdconjug,xcdsub,xcdeth,xcdcarbo,xcdcase)

! loop over the various C=C group of bonds (either mono or conjug. bonds)
! Note that cdcarbo, cdeth, cdsub (etc) are old version of the xcdcarbo
! xcdeth, xcdsub (etc) provided by cdcase2. Link between both set of table
! is performed here (this part must be rewritten!). 
    DO i=1,ncdtrack
      ncdcase=xcdcase(i)
      cdtable(:)=0 ; cdtable(1:SIZE(cdtrack,2))=cdtrack(i,:)
      cdsub(:)=0 ; cdcarbo(:,:)=0 ; cdeth(:,:)=0
      DO j=1,4
        IF (cdtrack(i,j)==0) CYCLE
        cdsub(j)=xcdsub(cdtrack(i,j))
        cdcarbo(j,1)=xcdcarbo(cdtrack(i,j),1)
        cdcarbo(j,2)=xcdcarbo(cdtrack(i,j),2)
        cdeth(j,1)=xcdeth(cdtrack(i,j),1)
        cdeth(j,2)=xcdeth(cdtrack(i,j),2)
      ENDDO

! case 1 : C=C bonds (conjugated or not)
! -------
      IF (ncdcase==1) THEN
        CALL hoadd_c1(chem,bond,group,cdtable,&
                      nr,flag,tarrhc,pchem,coprod,flag_del,&
                      pchem_del,coprod_del,sc_del,nrxref,rxref)

! case 2 : C=C-C=O bonds, only 1 C=C (i.e no conjugated Cds)
! -------
      ELSE IF (ncdcase==2) THEN
        CALL hoadd_c2(chem,bond,group,cdtable,cdsub,cdcarbo,&
                      nr,flag,tarrhc,pchem,coprod,nrxref,rxref)

! case 3 : -C=O-C=C-C=C-C=O- structure 
! -------
      ELSE IF (ncdcase==3) THEN 
        CALL hoadd_c3(chem,bond,group,cdtable,cdsub,&
                      nr,flag,tarrhc,pchem,coprod,flag_del,&
                      pchem_del,coprod_del,sc_del,nrxref,rxref)
        
! case 4 : -C=C-C-C=C-C=O-
! -------
      ELSE IF (ncdcase==4) THEN
        CALL hoadd_c4(chem,bond,group,cdtable,cdsub,cdcarbo, &
                    nr,flag,tarrhc,pchem,coprod,flag_del,&
                    pchem_del,coprod_del,sc_del,nrxref,rxref)
      
! case 5 : -CO-C=C-C=C<  
! -------
      ELSE IF (ncdcase==5) THEN
        CALL hoadd_c5(chem,bond,group,cdtable,cdcarbo,&
                      nr,flag,tarrhc,pchem,coprod,flag_del,&
                      pchem_del,coprod_del,sc_del,nrxref,rxref)
    
! case 6 : -C=C=O
! -------
      ELSE IF (ncdcase==6) THEN
        CALL hoadd_c6(chem,bond,group,nr,flag,tarrhc,pchem,coprod,nrxref,rxref)

! Case 7 : -C=C-O-
! ----------------
      ELSE IF (ncdcase==7) THEN
        CALL hoadd_c7(chem,bond,group,cdtable,cdeth, &
                      nr,flag,tarrhc,pchem,coprod,flag_del, &
                      pchem_del,coprod_del,sc_del,nrxref,rxref)
      ENDIF
    ENDDO
  ENDIF ! end of OH addition to C=C

! =================================
! OH ADDITION ON AROMATIC RING
! =================================
  IF (INDEX(chem,'c')/=0) &
    CALL hoadd_arom(chem,bond,group,nca,nr,flag,tarrhc,pchem,coprod,nrxref,rxref)

! ================================
! LUMP IDENTICAL REACTION
! ================================
! lump what is identical in the product, coproduct, and rate constant
  IF (nr>1) THEN
    DO i=1,nr-1
      jloop: DO j=i+1,nr
        IF (pchem(i)==pchem(j)) THEN
          DO k=1,SIZE(coprod,2)
            IF (coprod(i,k)/=coprod(j,k)) CYCLE jloop
          ENDDO
          IF (flag(i)==0) CYCLE jloop
          IF (ABS(tarrhc(i,3)-tarrhc(j,3)) > 0.001*ABS(tarrhc(i,3))) CYCLE jloop
          IF (ABS(tarrhc(i,2)-tarrhc(j,2)) > 0.001*ABS(tarrhc(i,2))) CYCLE jloop
          tarrhc(i,1)=tarrhc(i,1)+tarrhc(j,1)   ! same reaction if reach that point
          tarrhc(j,:)=0 ; flag(j)=0
        ENDIF
      ENDDO jloop
    ENDDO
  ENDIF

END SUBROUTINE mkhorx


END MODULE hochem
