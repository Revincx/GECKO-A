MODULE no3chem
IMPLICIT NONE
CONTAINS
!=======================================================================                                                                      
! PURPOSE : Set the VOC+OH mechanism. Read first the database for known 
! reaction rates and mechanisms. If the mechanism is not provided in the
! database, then generate the VOC + OH reactions. All selected reactions
! are written in the output mechanism file.
!=======================================================================                                                                      
SUBROUTINE no3_voc(idnam,chem,bond,group,brch,cut_off)
  USE keyparameter, ONLY: mxpd,mxnr,mecu,no3u,kno3u,mxcopd
  USE database, ONLY:nkno3db,kno3db_chem,kno3db_arr,kno3db_com, &         ! K_NO3 database
      nkwno3,nkwno3_pd,kwno3_rct,kwno3_pd,kwno3_copd,kwno3_yld,kwno3_com  ! VOC+NO3 mechanism
  USE references, ONLY:mxlcod
  USE keyflag, ONLY: wrtopeinfo,screenfg,sar_only_fg,yldcut,TK,Tref
  USE searching, ONLY: srh5
  USE reactool, ONLY: swap,rebond
  USE normchem, ONLY: stdchm
  USE dictstacktool, ONLY: bratio
  USE rxwrttool, ONLY:rxwrit,rxinit
  USE radchktool, ONLY: radchk
  USE toolbox, ONLY: stoperr,add1tonp,addrx,addref
  IMPLICIT NONE

  CHARACTER(LEN=*),INTENT(IN) :: idnam    ! name of current species
  CHARACTER(LEN=*),INTENT(IN) :: chem     ! formula of current species
  CHARACTER(LEN=*),INTENT(IN) :: group(:) ! group (without ring joining char)
  INTEGER,INTENT(IN) :: bond(:,:)    ! bond matrix
  REAL,INTENT(IN)    :: brch         ! max yield of the input species
  REAL,INTENT(IN)    :: cut_off      ! ratio threshold below which rx is ignored

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

  REAL    :: smkT, kT(mxnr)
  REAL    :: smk298, k298(mxnr)
  REAL    :: ratio(mxnr)
  REAL    :: allakT,allak298, kmax
  INTEGER :: sumflg
  REAL    :: brtio

  REAL    :: rk298, exparrhc(3)  ! know rate constant
  REAL    :: rkT  ! kT calculated with known rate constant coeffs
  INTEGER :: nno3prod

  CHARACTER(LEN=LEN(idnam)) :: p(mxpd), r(3)
  REAL    :: s(mxpd), arrh(3)
  INTEGER :: idreac, nlabel
  REAL    :: xlabel,folow(3),fotroe(4)

  CHARACTER(LEN=LEN(idnam)) :: pname(mxnr)
  CHARACTER(LEN=LEN(idnam)) :: pope
  REAL    :: a1,nt2,ae3, sope

  INTEGER :: nr,i,j,k,l,np,iloc,ngr
  REAL    :: kscale

  INTEGER,PARAMETER  :: mxrpd=2  ! max # of products returned by radchk sub
  CHARACTER(LEN=LEN(chem)) :: rdckpd(mxrpd)
  CHARACTER(LEN=LEN(idnam)) :: rdckcopd(mxrpd,mxcopd)
  REAL    :: sc(mxrpd)
  REAL    :: sc_del(mxnr,mxrpd)
  INTEGER :: nip

  LOGICAL :: kdbflg,rxdbflg
  
  CHARACTER(LEN=8),PARAMETER :: progname='no3_voc'
  CHARACTER(LEN=70)          :: mesg
  
  INTEGER,PARAMETER     :: mxcom=10
  INTEGER               :: nrxref(mxnr)
  CHARACTER(LEN=mxlcod) :: rxref(mxnr,mxcom)
  CHARACTER(LEN=mxlcod) :: comkdb(SIZE(kno3db_com,2))  ! reference for koh in db

! for available exp. rate constant, give priority to SAR Ea (not available for
! current SAR) or to the experimental Ea (recommended here for NO3 SAR) 
  INTEGER,PARAMETER :: expeaflg=1 ! set to 0 for default case

  IF (screenfg) WRITE(6,*)progname

!-------------
! initialize:
!-------------
  nr=0  ;  rk298=0.
  pchem(:)=' ' ; pchem_del(:)=' ' ; pchema(:)=' '
  flag(:)=0 ; flag_del(:)=0
  coprod(:,:)=' ' ; coprod_del(:,:)=' '
  tarrhc(:,:)=0.
  nrxref(:)=0 ; rxref(:,:)=' ' ; comkdb(:)=' '

  kdbflg=.FALSE. ; rxdbflg=.FALSE.
  tgroup(:)=group(:) ; tbond(:,:)=bond(:,:)

  IF (LEN_TRIM(chem)<1) RETURN
  ngr=COUNT(tgroup/=' ')

! -----------------------------------------------------------
! check if the species is known in the rate constant data base
! -----------------------------------------------------------
  iloc=srh5(chem,kno3db_chem,nkno3db)
  IF (iloc>0) THEN
    rk298=kno3db_arr(iloc,1)*298**kno3db_arr(iloc,2)*EXP(-kno3db_arr(iloc,3)/298.)
    rkT  =kno3db_arr(iloc,1)*TK**kno3db_arr(iloc,2) *EXP(-kno3db_arr(iloc,3)/TK)
    exparrhc(:)=kno3db_arr(iloc,:)
    comkdb(:)=kno3db_com(iloc,:)
    kdbflg=.TRUE.
  ENDIF

! -----------------------------------------------------------
! Check file no3_prod.dat for any known products & branching ratios 
! WARNING: for those reactions, yields < 5% are not automatically eliminated
! -----------------------------------------------------------
! SAR assessment - ignore the rate in the database
  IF (sar_only_fg) kdbflg=.FALSE.    

  IF (kdbflg) THEN
    iloc=srh5(chem,kwno3_rct,nkwno3)
    IF (iloc>0) THEN
      rxdbflg=.TRUE.  ;  nno3prod=nkwno3_pd(iloc)

! loop through known products
      DO j=1,nno3prod
        CALL addrx(progname,chem,nr,flag)
        pchema(j)=kwno3_pd(iloc,j)
        ratio(j)=kwno3_yld(iloc,j)

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
        kloop: DO k=1,SIZE(kwno3_copd,3)  ! loop possible coproducts
          IF (kwno3_copd(iloc,j,k)/=' ') THEN           
            DO l=1,mxcopd       
              IF (coprod(j,l)/=' ') CYCLE
              coprod(j,l)=kwno3_copd(iloc,j,k)
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
        DO i=1, SIZE(kwno3_com,2)
          IF (kwno3_com(iloc,i)/=' ') &
              CALL addref(progname,kwno3_com(iloc,i),nrxref(nr),rxref(nr,:),chem)
        ENDDO 

      ENDDO

    ENDIF  ! reactant has known products
  ENDIF

! ==========================
! MAKE THE VOC+NO3 MECHANISM
! ==========================

  IF (.NOT.rxdbflg) THEN

! perform all all NO3+VOC reaction
! --------------------------------
    CALL mkno3rx(chem,bond,group,ngr,nr,flag,pchem,coprod,&
                 flag_del,pchem_del,coprod_del,sc_del,tarrhc, &
                 nrxref,rxref)
    IF (nr==0) RETURN

! select pathways
! ---------------

! compute the rate at 298 K and at TK (requested T might not be 298 K)
    smkT=0. ; smk298=0. ; k298(:)=0.; kT=0.
    DO i=1,nr
      IF (flag(i)==0) CYCLE
      k298(i)=tarrhc(i,1) * 298.**tarrhc(i,2)*EXP(-tarrhc(i,3)/298.) 
!!      k298(i)=AMIN1(k298(i),2.0E-10)        ! do not exceed collision rate
      smk298=smk298+k298(i)
      kT(i) =tarrhc(i,1) * TK**tarrhc(i,2)*EXP(-tarrhc(i,3)/TK)
      smkT  =smkT+kT(i)
    ENDDO
    IF (smkT==0.) RETURN 

! write reaction rate for SAR assessment
    IF (sar_only_fg) WRITE(no3u,'(1PE12.2,2x,a)') smk298, TRIM(chem)

! remove pathways below cutoff at mechanism temperature
    WHERE (kT(:) < cut_off*smkT) flag(:)=0
    WHERE (kT(:) < (yldcut/brch)*smkT) flag(:)=0
    IF (SUM(flag)==0) THEN
      kmax=MAXVAL(kT)
      WHERE (kT(:) > kmax/2.) flag(:)=1
      sumflg=SUM(flag)
      IF (sumflg==0) THEN
        mesg="no path remains after cut_off "
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
    IF (allakT==0.) THEN
      mesg="allakT is 0 before division "
      CALL stoperr(progname,mesg,chem)
    ENDIF
    DO i=1,nr
      IF (flag(i)==1) ratio(i)=kT(i)/allakT
    ENDDO

! SET THE RATE OF EACH CHANNEL

! If rate not in the database, rescale tarrhc(i,1) to keep smk298
    IF (.NOT.kdbflg) THEN            
! Write kT to the kjvocj dataset
      WRITE(kno3u,'(a,2x,2(1pe10.2),a)') idnam, smk298, smkT," SAR, SAR"

      kscale=smk298/allak298          ! scaling factor for sum active k=smk298
      DO i=1,nr
        IF (flag(i)==1) tarrhc(i,1)=tarrhc(i,1)*kscale
      ENDDO

! If rate in the database then rescaling to experimental value needed.
! NOTE: The NO3 SAR does not provide Ea. Therefore, recommendation is expeaflg=1
    ELSE

      IF (expeaflg==0) THEN  ! use exp. value only at 298 (not recommended)
! Write k298 & kT to the kjvocj output dataset
        WRITE(kno3u,'(a,2x,2(1pe10.2),a)') idnam, rk298,smkT," DATA, SAR"

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

      ELSE                         ! full priority to exp. values (RECOMMENDED)
! Write DATA k298 & kT to the kjvocj output dataset
        WRITE(kno3u,'(a,2x,2(1pe10.2),a)') idnam, rk298,rkT," DATA, DATA"

        DO i=1,nr
          IF (flag(i)==1) THEN
            tarrhc(i,1)=exparrhc(1)*ratio(i) ! keep exp. value, scaled for each active channel
            tarrhc(i,2)=exparrhc(2)          ! keep exp. value 
            tarrhc(i,3)=exparrhc(3)          ! keep exp. Ea 
            CALL addref(progname,'KEXPBRSAR',nrxref(i),rxref(i,:),chem)
            DO j=1,SIZE(comkdb)              ! add ref for the exp. data
              IF (comkdb(j)/=' ') CALL addref(progname,comkdb(j),nrxref(i),rxref(i,:),chem)
            ENDDO
         ENDIF
        ENDDO
      ENDIF
    ENDIF
  ENDIF ! end of generated reactions

! =================================
! WRITE REACTIONS
! =================================
     
  DO i=1,nr
    IF (flag(i)==0) CYCLE

! initialize reaction
    CALL rxinit(r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe)
    arrh(:)=tarrhc(i,:)
    r(1)=idnam  ; r(2)='NO3 '
    s(1)=1.
    
    brtio=brch*ratio(i)
    CALL bratio(pchem(i),brtio,p(1),nrxref(i),rxref(i,:))
    pname(i)=p(1)          ! store for later use in operator information

! second product (if delocalisation allowed)
    np=1
    IF (flag_del(i)/=0) THEN
      np=np+1
      s(1)=sc_del(i,1) ; s(np)=sc_del(i,2)  ! overwrite the 1st coef. for s
      CALL bratio(pchem_del(i),brtio,p(np),nrxref(i),rxref(i,:))
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

! write out: reactions with NO3 are thermal reactions (idreac=0) and do not require labels
    CALL rxwrit(mecu,r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe,rxref(i,:))
  ENDDO

! =================================
! WRITE INFORMATION REQUIRED FOR OPERATOR
! =================================
  IF (wrtopeinfo) THEN !#Write information required for operator
    nt2=0. ; ae3=0. ; a1=0. 
    opeloop2: DO i=1,nr
      IF (flag(i)/=0) THEN
        a1=a1+tarrhc(i,1)
        IF ((tarrhc(i,2)/=nt2).OR.(tarrhc(i,3)/=ae3)) THEN
          a1=allakT ; nt2=0. ; ae3=0.  ! need revision, allak undefined if known chem.
          CYCLE opeloop2
        ELSE
          nt2=tarrhc(i,2) ; ae3=tarrhc(i,3)
        ENDIF 
      ENDIF
    ENDDO opeloop2
    WRITE(10,'(A16,A1,A6,A4,1X,ES10.3,1X,f4.1,1X,f7.0)') &
       '****INIT GNO3 + ','G',idnam,'****',a1,nt2,ae3

    DO i=1,nr
      IF (flag(i)==0) CYCLE
      WRITE(10,'(f5.3,2X,A1,A6)') ratio(i), 'G',pname(i)
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
    WRITE(10,*) 'end'
  ENDIF

END SUBROUTINE no3_voc

!=======================================================================
!=======================================================================
! PURPOSE: Generate the VOC+NO3 mechanism. The routine identify: 
! (1) all possible H-abstraction (-CH, -OH, -OOH, -C(O)OH
! (2) all possible NO3 addition to C=C bond
!=======================================================================
!=======================================================================
SUBROUTINE mkno3rx(chem,bond,group,ngr,nr,flag,pchem,coprod,&
                   flag_del,pchem_del,coprod_del,sc_del,tarrhc,&
                   nrxref,rxref)
  USE references, ONLY:mxlcod
  USE toolbox, ONLY: stoperr,addrx,addref
  USE reactool, ONLY: swap,rebond
  USE radchktool, ONLY: radchk
  USE normchem, ONLY: stdchm
  USE cdtool, ONLY: cdcase2
  USE no3addtool, ONLY: no3add_c1,no3add_c2,no3add_c3,no3add_c4,&
                        no3add_c5,no3add_c6,no3add_c7
  USE no3abstool, ONLY: rabsno3
  USE keyparameter, ONLY: mxtrk,mxlcd
  IMPLICIT NONE

  CHARACTER(LEN=*),INTENT(IN) :: chem ! formula of the VOC considered
  CHARACTER(LEN=*),INTENT(IN) :: group(:) ! group (without ring joining char)
  INTEGER,INTENT(IN) :: bond(:,:)     ! bond matrix
  INTEGER,INTENT(IN) :: ngr           ! # number of nodes
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

  INTEGER :: i,j,k,nc,tempring,rxnflg

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

  CHARACTER(LEN=8),PARAMETER :: progname='mkno3rx'
  CHARACTER(LEN=70)          :: mesg

  nr=0  ;  flag(:)=0  ;  tarrhc(:,:)=0.
  pchem(:)=' '     ; coprod(:,:)=' ' 
  pchem_del(:)=' ' ; coprod_del(:,:)=' '
  flag_del(:)=0  ; sc_del(:,:)=0. 

  tgroup(:)=group(:) ; tbond(:,:)=bond(:,:)
  carb=0 ; cdtable(:)=0 ; tcdtable(:)=0 ; cdsub(:)=0
  ncdcase=0 ; cdcarbo(:,:)=0 ; cdeth(:,:)=0
  rxnflg=1 ; nrxref(:)=0 ; rxref(:,:)=' '
  
! =================================
!  H ABSTRACTION                              
! =================================

! H-ABSTRACTION FROM CH
! ---------------------
  DO i=1,ngr
    IF (INDEX(group(i),'CH')==0) CYCLE
    CALL addrx(progname,chem,nr,flag)

! compute associated rate constant: k=A*T**2*exp(-E/T)
    CALL rabsno3(tbond,tgroup,i,tarrhc(nr,:))
    CALL addref(progname,'HABSAR',nrxref(nr),rxref(nr,:),chem)
    CALL addref(progname,'JK10KEV000',nrxref(nr),rxref(nr,:),chem)
!    tarrhc(nr,:)=arrh(:)

! remove one H atom
    IF (group(i)(1:3)=='CH3')      THEN ; pold='CH3' ; pnew='CH2'
    ELSE IF (group(i)(1:3)=='CH2') THEN ; pold='CH2' ; pnew='CH'
    ELSE IF (group(i)(1:3)=='CHO') THEN ; pold='CHO' ; pnew='CO'
    ELSE IF (group(i)(1:2)=='CH')  THEN ; pold='CH'  ; pnew='C'  
    ENDIF
    CALL swap(group(i),pold,tgroup(i),pnew)

! add radical dot at end of group:
    nc=INDEX(tgroup(i),' ') ; tgroup(i)(nc:nc)='.'

! re-build molecule, including ring if present:
    CALL rebond(tbond,tgroup,tempkc,tempring)

! check radical and find co-products:
    CALL radchk(tempkc,rdckpd,rdckcopd,nip,sc,nrxref(nr),rxref(nr,:))
    pchem(nr)=rdckpd(1) ; sc_del(nr,1)=sc(1)
    CALL stdchm(pchem(nr))
    coprod(nr,:)=rdckcopd(1,:)
    DO j=1,SIZE(coprod,2)
      IF (coprod(nr,j)==' ') THEN 
        coprod(nr,j)='HNO3' ; EXIT 
      ENDIF
      IF (j==SIZE(coprod,2)) THEN
        mesg="No slot available to add HNO3 "
        CALL stoperr(progname,mesg,chem)
      ENDIF  
    ENDDO
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
  DO i=1,ngr
    IF (INDEX(group(i),'(OH)')==0) CYCLE
    IF (INDEX(group(i),'CO(OH)')/=0) CYCLE
    CALL addrx(progname,chem,nr,flag)

! assign rate constant by Kerdouci et al., 2010
     tarrhc(nr,1)=2.00E-17 ; tarrhc(nr,2)=0. ; tarrhc(nr,3)=0.  
    CALL addref(progname,'HABSAR1',nrxref(nr),rxref(nr,:),chem)
    CALL addref(progname,'JK10KEV000',nrxref(nr),rxref(nr,:),chem)

! replace hydroxy group by alkoxy group:
    pold='(OH)' ; pnew='(O.)'
    CALL swap(group(i),pold,tgroup(i),pnew)
      
! rebuild, check, and find coproducts:
    CALL rebond(tbond,tgroup,tempkc,tempring)
    CALL radchk(tempkc,rdckpd,rdckcopd,nip,sc,nrxref(nr),rxref(nr,:))
    pchem(nr)=rdckpd(1)
    IF (nip/=1) THEN
      mesg="unexpected 2 products for ROH+OH "
      CALL stoperr(progname,mesg,chem)
    ENDIF
    coprod(nr,:)=rdckcopd(1,:)
    DO j=1,SIZE(coprod,2)
      IF (coprod(nr,j)==' ') THEN 
        coprod(nr,j)='HNO3' ; EXIT 
      ENDIF
      IF (j==SIZE(coprod,2)) THEN
        mesg="no slot available to add HNO3 "
        CALL stoperr(progname,mesg,chem)
      ENDIF
    ENDDO
    CALL stdchm(pchem(nr))

! reset group: 
    tgroup(i)=group(i)
  ENDDO

! =================================
!  NO3 ADDITION TO >C=C< BOND 
! =================================

! The parameterization used to assign a rate constant for NO3 addition to
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
      
! check if the molecule contains C=C bond 
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

! case 1 : C=C bonds conjugated or not, without any C=C-C=O structure
! -------
    IF (ncdcase==1) THEN
      CALL no3add_c1(chem,bond,group,cdtable,cdsub, &
                    nr,flag,tarrhc,pchem,coprod,flag_del,        &
                    pchem_del,coprod_del,sc_del,nrxref,rxref)
      
! case 2 : C=C-C=O bonds, only 1 C=C (i.e no conjugated Cds)
! -------
    ELSE IF (ncdcase==2) THEN
      CALL no3add_c2(chem,bond,group,cdtable,cdcarbo,   &
                     nr,flag,tarrhc,pchem,coprod,flag_del, &
                     pchem_del,coprod_del,sc_del,nrxref,rxref)

! case 3 : -CO-C=C-C=C-C=O- structure 
! -------
    ELSE IF (ncdcase==3) THEN
      CALL no3add_c3(chem,bond,group,cdtable,cdcarbo,   &
                     nr,flag,tarrhc,pchem,coprod,flag_del, &
                     pchem_del,coprod_del,sc_del,nrxref,rxref)

! case 4 : -C=C-C-C=C-C=O-
! -------
    ELSE IF (ncdcase==4) THEN
      CALL no3add_c4(chem,bond,group,cdtable,cdsub,cdcarbo, &
                     nr,flag,tarrhc,pchem,coprod,flag_del, &
                     pchem_del,coprod_del,sc_del,nrxref,rxref)

! case 5 : -CO-C=C-C=C< 
! -------
    ELSE IF (ncdcase==5) THEN
      CALL no3add_c5(chem,bond,group,cdtable,cdsub,cdcarbo, &
                     nr,flag,tarrhc,pchem,coprod,flag_del,pchem_del, &
                     coprod_del,sc_del,nrxref,rxref)

! case 6 : -C=C=O
! -------
    ELSE IF (ncdcase==6) THEN
      CALL no3add_c6(chem,bond,group,nr,flag,tarrhc,pchem,coprod,nrxref,rxref)

! case 7 : -C=C-O-
! ----------------
    ELSE IF (ncdcase==7) THEN
      CALL no3add_c7(chem,bond,group,cdtable,cdeth,nr,flag,tarrhc, &
                     pchem,coprod,flag_del,pchem_del,&
                     coprod_del,sc_del,nrxref,rxref)
    ENDIF
    
    ENDDO  
! end of NO3 addition to C=C
  ENDIF

! lump what is identical in the product, coproduct, and rate constant
! --------------------------------------------------------------------
  IF (nr>1) THEN
    DO i=1,nr-1
      jloop: DO j=i+1,nr
        IF (pchem(i)==pchem(j))THEN 
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

END SUBROUTINE mkno3rx

END MODULE no3chem
