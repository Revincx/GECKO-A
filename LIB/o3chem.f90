MODULE o3chem
IMPLICIT NONE
CONTAINS

! SUBROUTINE o3_voc(chem,rdct,bond,group,brch,cut_off)
! SUBROUTINE mko3rx(chem,bond,group,nr,flag,pchem,tarrhc,parent_bond)

! ======================================================================
! PURPOSE: Set the VOC+O3 mechanism. Read first the database for known 
! reaction rates and mechanisms. If the mechanism is not provided in the
! database, then generate the VOC+O3 reactions. All selected reactions
! are written in the output mechanism file.
! 
! NOTE: Decomposition/stabilisation of hot criegee is performed before
! writing the reaction in the output file (see mko3rx routine).
! Reaction products include stabilized criegee only.
! ======================================================================
SUBROUTINE o3_voc(idnam,chem,bond,group,zebond,brch,cut_off)
  USE keyparameter, ONLY: mxpd,mxnr,mxcopd,mecu,waru,o3u
  USE keyflag, ONLY: sar_only_fg,TK,screenfg
  USE references, ONLY:mxlcod
  USE searching, ONLY: srh5
  USE normchem, ONLY: stdchm
  USE database, ONLY:nko3db,ko3db_chem,ko3db_arr,ko3db_com, &           ! K_O3 database
           nkwo3,nkwo3_pd,kwo3_rct,kwo3_pd,kwo3_copd,kwo3_yld,kwo3_com  ! VOC+O3 mechanism
  USE dictstacktool, ONLY: bratio
  USE rxwrttool, ONLY:rxwrit,rxinit
  USE radchktool, ONLY: radchk
  USE toolbox, ONLY: kval,stoperr,add1tonp,addrx,addref
  USE reac_infotool, ONLY: fullref
  IMPLICIT NONE

  CHARACTER(LEN=*),INTENT(IN) :: idnam        ! name of current species
  CHARACTER(LEN=*),INTENT(IN) :: chem         ! formula of current species
  CHARACTER(LEN=*),INTENT(IN) :: group(:)     ! group (without ring joining char)
  INTEGER,INTENT(IN) :: bond(:,:)             ! bond matrix
  INTEGER,INTENT(IN) :: zebond(:,:)           ! cis/trans info on C=C bond
  REAL,INTENT(IN)    :: brch                  ! max yield of the input species
  REAL,INTENT(IN)    :: cut_off               ! ratio threshold below which rx is ignored

  CHARACTER(LEN=LEN(group(1))) :: tgroup(SIZE(group))
  INTEGER              :: tbond(SIZE(bond,1),SIZE(bond,2))
  CHARACTER(LEN=LEN(chem)) :: tempkc
  CHARACTER(LEN=LEN(chem)) :: pchem(mxnr,2)
  INTEGER :: flag(mxnr)
  INTEGER :: nr,np,i,j,k,l,ngr,iloc,no3prod,icode
  INTEGER :: ncdtrack
  REAL    :: ratio(mxnr)
  REAL    :: tarrhc(mxnr,3),exparrhc(3)
  REAL    :: rk298,smk298,k298(mxnr)
  REAL    :: rkT,smkT,kT(mxnr)
  REAL    :: kscale
  LOGICAL :: kdbflg,rxdbflg
  CHARACTER(LEN=LEN(idnam)) :: coprod(mxnr,mxcopd)

  INTEGER                   :: rxnp(mxnr)
  CHARACTER(LEN=LEN(idnam)) :: rxp(mxnr,mxpd)
  REAL                      :: rxs(mxnr,mxpd)

  CHARACTER(LEN=LEN(idnam)) :: p(mxpd), r(3)
  REAL    :: s(mxpd), arrh(3)
  INTEGER :: idreac, nlabel
  REAL    :: xlabel,folow(3),fotroe(4)

  INTEGER,PARAMETER  :: mxrpd=2         ! max # of products returned by radchk sub
  CHARACTER(LEN=LEN(chem))  :: rdckpd(mxrpd)
  CHARACTER(LEN=LEN(idnam)) :: rdckcopd(mxrpd,SIZE(coprod,2))
  REAL    :: sc(mxrpd)
  INTEGER :: nip
  
  INTEGER,PARAMETER     :: mxcom=10
  INTEGER               :: nrxref(mxnr)
  CHARACTER(LEN=mxlcod) :: rxref(mxnr,mxcom)
  CHARACTER(LEN=mxlcod) :: comkdb(SIZE(ko3db_com,2))  ! reference for ko3 in db

  CHARACTER(LEN=6),PARAMETER :: progname='o3_voc'
  CHARACTER(LEN=80) :: mesg

! for available exp. rate constant, give priority to SAR Ea (default) or
! to the  experimental Ea (not recommended) 
  INTEGER,PARAMETER :: expeaflg=0 ! set to 0 for default case

  IF (screenfg) WRITE(6,*)progname

  nrxref(:)=0 ; rxref(:,:)=' ' ; comkdb(:)=' ' ; kdbflg=.FALSE.
  nr=0 ;  rk298=0. ; rxdbflg=.FALSE. ;  coprod(:,:)=' ' ; pchem(:,:)=' ' 
  rkT = 0. ; flag(:)=0 ; tarrhc(:,:)=0. ;  exparrhc(:)=0.  ; ratio(:)=0.

  tgroup(:)=group(:) ; tbond(:,:)=bond(:,:)

  ngr=COUNT(group/=' ')

! check that calling "ozone reaction" make sense
  IF (INDEX(chem,'Cd')==0) THEN
    mesg='the O3_VOC routine was called without C=C bond in chem'
    CALL stoperr(progname,mesg,chem)
  ENDIF

! ------------------------------------------------------------
! check if the species is known in the rate constant data base
! ------------------------------------------------------------
  iloc=srh5(chem,ko3db_chem,nko3db)
  IF (iloc>0) THEN
    rk298=ko3db_arr(iloc,1)*298**ko3db_arr(iloc,2)*exp(-ko3db_arr(iloc,3)/298.)
    rkT=ko3db_arr(iloc,1)*TK**ko3db_arr(iloc,2)*exp(-ko3db_arr(iloc,3)/TK)
    !IF (rk298<1.0E-20) RETURN     ! if database rate constant < 1.0E-20, no reaction
    IF (rkT<1.0E-20) RETURN     ! if database rate expression < 1.0E-20, no reaction
    exparrhc(:)=ko3db_arr(iloc,:)
    comkdb(:)=ko3db_com(iloc,:)
    kdbflg=.TRUE.
  ENDIF
      
! -----------------------------------------------------------
! Check if rate or products are provided in the database. 
! Note: yields below threshold are not automatically eliminated.                                                     
! -----------------------------------------------------------

! SAR assessment - ignore the rate in the database
  IF (sar_only_fg) kdbflg=.FALSE.

  IF (kdbflg) THEN
    iloc=srh5(chem,kwo3_rct,nkwo3)
    IF (iloc>0) THEN
      rxdbflg=.TRUE.  ;  no3prod=nkwo3_pd(iloc)
  
! loop through known products
      DO j=1,no3prod
        CALL addrx(progname,chem,nr,flag)
        tempkc=kwo3_pd(iloc,j)
        ratio(j)=kwo3_yld(iloc,j)
              
! check the products
        IF (INDEX(tempkc,'.')/=0) THEN
          CALL radchk(tempkc,rdckpd,rdckcopd,nip,sc,nrxref(nr),rxref(nr,:))
          pchem(j,1)=rdckpd(1)
          coprod(j,:)=rdckcopd(1,:)   
          IF (nip==2) THEN
            mesg="nip/=1 unexpected "
            CALL stoperr(progname,mesg,chem)
          ENDIF
        ELSE
          pchem(j,1)=tempkc
        ENDIF
        CALL stdchm(pchem(j,1))

! add known coproduct(s) if necessary
        kloop: DO k=1,mxcopd  ! loop possible coproducts
          IF (kwo3_copd(iloc,j,k)/=' ') THEN                   
            DO l=1,SIZE(coprod,2)           
              IF (coprod(j,l)/=' ') CYCLE  ! find available slot
              coprod(j,l)=kwo3_copd(iloc,j,k)
              CYCLE kloop
            ENDDO
            IF (l==SIZE(coprod,2)) THEN
              mesg="l==mxcopd, no slot available "
              CALL stoperr(progname,mesg,chem)
            ENDIF   
          ELSE
            EXIT kloop
          ENDIF          ! coproduct is non-zero
        ENDDO kloop      ! coproducts loop

! add reaction rate
        tarrhc(j,1)=exparrhc(1)*ratio(j) ! keep exp. value
        tarrhc(j,2)=exparrhc(2)          ! keep exp. value 
        tarrhc(j,3)=exparrhc(3)          ! keep exp. Ea 
        CALL addref(progname,'RXEXP',nrxref(nr),rxref(nr,:),chem)
        DO i=1,SIZE(comkdb)
          IF (comkdb(i)/=' ') CALL addref(progname,comkdb(i),nrxref(nr),rxref(nr,:),chem)
        ENDDO

! add references to the reaction
        DO i=1, SIZE(kwo3_com,2)
          IF (kwo3_com(iloc,i)/=' ') &
              CALL addref(progname,kwo3_com(iloc,i),nrxref(nr),rxref(nr,:),chem)
        ENDDO 

      ENDDO 

    ENDIF  ! reactant has known products
  ENDIF

! =========================
! MAKE THE VOC+O3 MECHANISM
! =========================
 
  IF (.NOT.rxdbflg) THEN

! perform all O3+VOC reaction
! --------------------------------
    CALL mko3rx(chem,bond,group,zebond,brch,cut_off,ncdtrack,nr, &
                flag,tarrhc,rxnp,rxp,rxs,icode,nrxref,rxref)

    IF (nr==0) THEN   ! icode must be 1 if no reaction returned 
      IF (icode==1) THEN
        RETURN        ! reactions are ignored for some species
      ELSE
        mesg="expected O3 reaction not found"
        CALL stoperr(progname,mesg,chem)
      ENDIF
    ENDIF

! all reactions must be considered (filtered in mko3rx)
! -----------------------------------------------------

    ! compute the rate at 298 K and TK
    smk298=0. ; k298(:)=0.
    smkT=0. ; kT(:)=0.
    DO i=1,nr
      IF (flag(i)==0) CYCLE
      k298(i)=kval(tarrhc(i,:),298.)
      kT(i)  =kval(tarrhc(i,:),TK)
      smk298=smk298+k298(i)
      smkT  =smkT  +kT(i)
    ENDDO

    ! write reaction rate for SAR assessment
    IF (sar_only_fg) THEN
      IF (smk298/=0.) THEN ; WRITE(o3u,'(1ES12.3,2x,a)') smk298, TRIM(chem)
      ELSE               ; WRITE(o3u,'(1ES12.3,2x,a)')  1E-35, TRIM(chem)
      ENDIF
    ENDIF
    
! If rate in the database then rescale to experimental value
! -----------------------------------------------------------
    IF (kdbflg) THEN            
    
      ! set branching ratio
      ratio(:)=0.
      DO i=1,nr
        !IF (flag(i)==1) ratio(i)=k298(i)/smk298
        IF (flag(i)==1) ratio(i)=kT(i)/smkT
      ENDDO

      ! only one Cd track simple (C=C or C=C-C=C) => use experimental values 
      IF (ncdtrack==1) THEN  
        DO i=1,nr            ! nr could either be 1 (simple C=C) or 2 (C=C-C=C)
          IF (flag(i)==1) THEN
            tarrhc(i,1)=exparrhc(1)*ratio(i)
            tarrhc(i,2)=exparrhc(2)  ;  tarrhc(i,3)=exparrhc(3)
            CALL addref(progname,'KEXPBRSAR',nrxref(i),rxref(i,:),chem)
            DO j=1,SIZE(comkdb)    ! add ref for the exp. data
              IF (comkdb(j)/=' ') CALL addref(progname,comkdb(j),nrxref(i),rxref(i,:),chem)
            ENDDO
          ENDIF
        ENDDO

      ! multiple tracks - scaling required
      ELSE                   
        IF (expeaflg==0) THEN        ! use exp. value at 298 only (DEFAULT)
          kscale=rk298/smk298        ! scaling factor for SAR rate = exp. value @298K 
          DO i=1,nr
            IF (flag(i)==1) THEN
              tarrhc(i,1)=tarrhc(i,1)*kscale ! keep SAR values
              CALL addref(progname,'KEXP298',nrxref(i),rxref(i,:),chem)
              DO j=1,SIZE(comkdb)    ! add ref for the exp. data
                IF (comkdb(j)/=' ') CALL addref(progname,comkdb(j),nrxref(i),rxref(i,:),chem)
              ENDDO
            ENDIF
          ENDDO

        ELSE                         ! full priority to exp. values (NOT RECOMMANDED)
          DO i=1,nr
            IF (flag(i)==1) THEN
              tarrhc(i,1)=exparrhc(1)*ratio(i) ! keep exp. value, scaled for each active channel
              tarrhc(i,2)=exparrhc(2) ; tarrhc(i,3)=exparrhc(3)         ! keep exp. value 
              CALL addref(progname,'KEXPBRSAR',nrxref(i),rxref(i,:),chem)
              DO j=1,SIZE(comkdb)    ! add ref for the exp. data
                IF (comkdb(j)/=' ') CALL addref(progname,comkdb(j),nrxref(i),rxref(i,:),chem)
              ENDDO
            ENDIF
          ENDDO
        ENDIF
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
    r(1)=idnam ;  r(2)='O3 '
       
    IF (.NOT.rxdbflg) THEN
      s(:)=rxs(i,:)  ;  p(:)=rxp(i,:)

    ELSE 
      DO j=1,SIZE(coprod,2)
        IF (coprod(i,j)(1:1)/=' ') THEN
          CALL add1tonp(progname,chem,np)   
          s(np)=1.0 ; p(np)=coprod(i,j)
        ENDIF
      ENDDO
    ENDIF

! write out: reaction with O3 are thermal reaction (idreac=0 and does 
! not require labels)
    CALL rxwrit(mecu,r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe,rxref(i,:))
  ENDDO  
  
END SUBROUTINE o3_voc

!=======================================================================
!=======================================================================
! PURPOSE: Generate the VOC+O3 mechanism. The routine identifies all 
! possible O3 addition to C=C bonds. A specific reaction is considered 
! for each C=C bond on the carbon skeleton (with related arrhenius 
! parameter). A primary ozononide (POZ) is associated to each C=C. 
! Each POZ decompose according to (1) a biradical peroxy-alkoxy route 
! (minor) and (2) carbonyl + hot criegee (2 routes), as described in 
! Newland et al. 2022. The decomposition of hot Criegee is then
! considered. Products returned are lists of products (short names) and
! related stoi. coefficients (rxp and rxs).   
!=======================================================================
!=======================================================================
SUBROUTINE mko3rx(chem,bond,group,zebond,brch,cut_off,ncdtrack,nr, &
                  flag,tarrhc,rxnp,rxp,rxs,icode,nrxref,rxref)
  USE references, ONLY:mxlcod
  USE keyparameter, ONLY: mxtrk,mxlcd,mxpd,mxlco,mxnr,waru
  USE keyflag, ONLY: TK,screenfg
  USE cdtool, ONLY: cdcase2
  USE o3addtool, ONLY: o3rate_mono,o3rate_conj
  USE toolbox, ONLY: kval,add1tonp,addrx,stoperr,addref
  USE hotcriegeechem, ONLY: hot_criegee
  USE atomtool, ONLY:getatoms
  USE dictstacktool, ONLY: bratio
  USE criegeetool, ONLY: add2p
  USE poztool, ONLY:pozprod,poz_decomp,get_pozratio
  USE ringtool, ONLY: findring
  IMPLICIT NONE

  CHARACTER(LEN=*),INTENT(IN) :: chem        ! formula of the VOC considered
  CHARACTER(LEN=*),INTENT(IN) :: group(:)    ! group (without ring joining char)
  INTEGER,INTENT(IN)  :: bond(:,:)           ! bond matrix
  INTEGER,INTENT(IN)  :: zebond(:,:)         ! cis/trans info on C=C bond
  REAL,INTENT(IN)     :: brch                ! max yield of the input species
  REAL,INTENT(IN)     :: cut_off             ! ratio threshold below which rx is ignored
  INTEGER,INTENT(OUT) :: ncdtrack            ! # of Cd track (i.e. C=C or  C=C-C=C structure) 
  INTEGER,INTENT(OUT) :: nr                  ! # of reaction channel
  INTEGER,INTENT(OUT) :: flag(:)             ! flag for active rxn channel
  REAL,INTENT(OUT)    :: tarrhc(:,:)         ! arrhenius coefficient for reaction i (i.e. the ith POZ)
  INTEGER,INTENT(OUT) :: rxnp(:)             ! # of pdcts in each the ith rxn channel
  CHARACTER(LEN=*),INTENT(OUT) :: rxp(:,:)   ! list of pdcts (short names) in the ith rxn channel
  REAL,INTENT(OUT)    :: rxs(:,:)            ! stoi. coef. of pdcts in the ith rxn channel
  INTEGER,INTENT(OUT) :: icode               ! info code (check for errors) 
  INTEGER,INTENT(OUT) :: nrxref(:)           ! # of references added for each reaction
  CHARACTER(LEN=*),INTENT(OUT) :: rxref(:,:) ! references for the reaction (1st index is rx #)

  INTEGER :: cdtrack(mxtrk,mxlcd)            ! cdtrack(i,j), Cd nodes for the ith track         
  INTEGER :: xcdconjug(mxtrk)                ! flag for conjugated tracks (C=C-C=C)
  INTEGER :: xcdsub(SIZE(group))             ! # of nodes bonded to a Cd group
  INTEGER :: xcdeth(SIZE(group),2)           ! -OR group # bounded to a Cd (max 2 -OR groups)
  INTEGER :: xcdcarbo(SIZE(group),2)         ! -CO- group # bounded to a Cd (max 2 CO groups)  
  INTEGER :: xcdcase(mxtrk)                  ! "case" of the track

  INTEGER :: npoz                            ! number of POZ in O3+chem (i.e. # of C=C bonds)
  INTEGER :: pozcc(mxnr,2)                   ! C nodes making the ith POZ 
  REAL    :: pozratio(2)                     ! for conjugated C=C, POZ ratio of each [C=C]-[C=C]
  REAL    :: pozyci(2)                       ! hot criegee yield at each C node (i.e. side) of the poz
  INTEGER :: pozatom                         ! number O, C, N atoms in POZ
  REAL    :: pozarrhc(mxnr,3)                ! arrhenius coefficient for POZ i formation
  REAL    :: k298(mxnr),kT(mxnr)             ! rate for reaction [i] 
  REAL    :: smk298,smkT                     ! overall rate
  LOGICAL :: loendo                          ! true if C=C bond is endo cyclic
  INTEGER :: rxpozflg(mxnr)                  ! flag rxn to be considered for a given poz
  LOGICAL :: lorescale                       ! true if rate constant must be rescaled
  REAL    :: kscale, allak298                ! scaling factor and all k
  REAL    :: ratio, allakT                   ! scaling factor and all k
  
  INTEGER              :: np                 ! # of product in a given channel (in a given rxn)
  CHARACTER(LEN=mxlco) :: p(SIZE(rxs,2))     ! list of the pdcts in a given channel
  REAL                 :: s(SIZE(rxs,2))     ! stoi. coef. of the pdcts
  INTEGER              :: npchannel(2)             ! # of products in channel i (i=[1,2])
  CHARACTER(LEN=mxlco) :: pchannel(2,SIZE(rxs,2))  ! list of the pdcts in channel i (i=[1,2])
  REAL                 :: schannel(2,SIZE(rxs,2))  ! stoi. coef. in channel i (i=[1,2])
  
  CHARACTER(LEN=LEN(chem)) :: chemcri        ! formula of the criegee product
  CHARACTER(LEN=LEN(chem)) :: chemcarb       ! formula of the carbonyl product
  CHARACTER(LEN=LEN(chem)) :: chemacyl       ! formula of the acyl radical product
  REAL    :: chbratio                        ! overall yield of carbonylhydroperoxide route
  REAL    :: iayield                         ! yield of the channel from POZ decomp.
  REAL    :: chyield                         ! yield of carbonylhydroperoxide route
  REAL    :: yci                             ! criegee yield (from POZ decomp & Carb-OOH channel)
  REAL    :: zconform                        ! yield of the Z conformer 
  INTEGER :: ia,icri,icarb
  REAL    :: brtio
  INTEGER :: rngflg
  INTEGER :: ring(SIZE(group))

  INTEGER :: xxc,xxh,xxn,xxo,xxr,xxs,xxfl,xxbr,xxcl
  INTEGER :: i,j,rxnflg,ngr,ni,nj
  REAL    :: arrhc(3)

  CHARACTER(LEN=7),PARAMETER :: progname='mko3rx'
  CHARACTER(LEN=70) :: mesg

  IF (screenfg) WRITE(6,*)progname

  nr=0  ;  flag(:)=0 ; pozarrhc(:,:)=0. ; tarrhc(:,:)=0. ; arrhc(:)=0.
  nrxref(:)=0  ;  rxref(:,:)=' '
  npoz=0  ; pozcc(:,:)=0 ; pozratio(:)=0. 
  s(:)=0. ; p(:)=' '
  rxnp(:)=0 ; rxp(:,:)=' ' ; rxs(:,:)=0. 
  k298(:)=0. ; kT(:)=0.
  npchannel(:)=0 ; schannel(:,:)=0. ; pchannel(:,:)=' '   

! set carbonylhydroperoxide route to ethene only
  IF (chem=='CdH2=CdH2 ') THEN ; chbratio=0.12  ; ELSE ; chbratio=0. ; ENDIF 
  
  rxnflg=1 ; icode=0
  ngr=COUNT(group/=' ')

! compute the number of C, N, O atoms in poz (needed to compute SCI yield next)
  CALL getatoms(chem,xxc,xxh,xxn,xxo,xxr,xxs,xxfl,xxbr,xxcl)
  pozatom=xxc+xxn+xxo+3                     ! 3 O added (not in chem but in POZ)

! NITRO-ALKENE: are assumed to have a O3 reactivity slow compared to OH
! reactivity and are ignored here (waiting for additional information).
  IF (INDEX(chem,'Cd(NO2)')/=0)  icode=1 
  IF (INDEX(chem,'CdH(NO2)')/=0) icode=1
  IF (icode==1) RETURN 

! get Cd tracks and nodes
  CALL cdcase2(chem,bond,group,rxnflg,ncdtrack,cdtrack,&
               xcdconjug,xcdsub,xcdeth,xcdcarbo,xcdcase)

! COMPUTE ARRHENIUS PARAMETER FOR EACH POZ  
! ----------------------------------------
  DO i=1,ncdtrack
    IF (xcdconjug(i)==0) THEN
      CALL o3rate_mono(bond,zebond,group,ngr,cdtrack(i,:),xcdsub,arrhc)  ! non conjugated C=C
      npoz=npoz+1
      pozarrhc(npoz,:)=arrhc(:)
      pozcc(npoz,1:2)=cdtrack(i,1:2)
    ELSE IF (xcdconjug(i)==1) THEN
      CALL o3rate_conj(bond,group,ngr,cdtrack(i,:),xcdsub,arrhc)  ! conjugated C=C
      CALL get_pozratio(cdtrack(i,:),xcdsub,pozratio)             ! fraction of O3 added on each C=C
      DO j=1,2                                                    ! loop over each C=C
        IF (pozratio(j)==0.) THEN
          mesg="O3 addition not found on a C=C bond in a C=C-C=C structure"
          CALL stoperr(progname,mesg,chem)
        ENDIF
        npoz=npoz+1
        pozarrhc(npoz,1)=arrhc(1)*pozratio(j)                       ! set rate for each rxn
        pozarrhc(npoz,2:3)=arrhc(2:3)
        IF (j==1) pozcc(npoz,1:2)=cdtrack(i,1:2)                  ! poz located at 1st C=C
        IF (j==2) pozcc(npoz,1:2)=cdtrack(i,3:4)                  ! poz located at 2nd C=C
      ENDDO
    ENDIF
  ENDDO

! select pathway to keep
! ----------------------
  smk298=0. ; smkT=0. ; rxpozflg(:)=0
  DO i=1,npoz
    rxpozflg(i)=1 ; k298(i)=kval(pozarrhc(i,:),298.)  ;  smk298=smk298+k298(i)
    rxpozflg(i)=1 ; kT(i)=kval(pozarrhc(i,:),TK)  ;  smkT=smkT+kT(i)
  ENDDO

  ! if total rate constant < 1.0E-20, then no reaction
  !IF (smk298 < 1.0E-20) THEN
  IF (smkT < 1.0E-20) THEN
    WRITE(waru,'(a,a)') '--warning-- k_O3<1E-20. Ignore reaction for: ',TRIM(chem)
    icode=1 ; RETURN
  ENDIF

  ! remove pathways below cutoff
  lorescale=.FALSE.
  DO i=1,npoz
    !IF (k298(i) < cut_off*smk298) THEN 
    IF (kT(i) < cut_off*smkT) THEN 
      rxpozflg(i)=0 ; lorescale=.TRUE.
    ENDIF
  ENDDO
  IF (SUM(rxpozflg)==0) THEN
    mesg="no path remains after cut_off"
    CALL stoperr(progname,mesg,chem)
  ENDIF

  ! rescale rate constant if needed
  IF (lorescale) THEN
    allakT=0.
    allak298=0.
    DO i=1,npoz
      IF (rxpozflg(i)==1) THEN
        allak298=allak298+k298(i)
        allakT=allakT+kT(i)
      ENDIF
    ENDDO
 
    !kscale=smk298/allak          ! scaling factor for sum active k
    ratio=smkT/allakT         ! scaling factor for sum active k
    DO i=1,npoz
      !IF (rxpozflg(i)==1) pozarrhc(i,1)=pozarrhc(i,1)*kscale
      IF (rxpozflg(i)==1) pozarrhc(i,1)=pozarrhc(i,1)*ratio
    ENDDO
  ENDIF

! MAKE PRODUCTS FOR EACH POZ
! --------------------------
  pozloop: DO i=1,npoz
  
    ! chek if C=C is endocyclic (needed for sci yield) 
    CALL findring(pozcc(i,1),pozcc(i,2),ngr,bond,rngflg,ring)
    IF (rngflg==1) THEN ; loendo=.TRUE. ; ELSE ; loendo=.FALSE. ; ENDIF
    
    ! ignore poz chemistry for reaction rate below threshold
    IF (rxpozflg(i)==0) CYCLE pozloop

    ! add one reaction per POZ
    CALL addrx(progname,chem,nr,flag)     
    npchannel(:)=0 ; schannel(:,:)=0. ; pchannel(:,:)=' '   
    CALL addref(progname,'MJ20KMV000',nrxref(nr),rxref(nr,:),chem)  ! ref to k
    CALL addref(progname,'MN22SAR1  ',nrxref(nr),rxref(nr,:),chem)  ! ref to POZ and CI* pdcts

    ! POZ decomposition - get hot CI yield for each Cd of the POZ structure 
    CALL poz_decomp(chem,group,bond,pozcc(i,:),xcdsub,ring,pozyci(:))

    ! make 2 channels per poz (i.e. >Ca.OO. + >Cb=O  and  >Ca=O + >Cb.OO.)
    ialoop: DO ia=1,2
      np=0 ; p(:)=' ' ; s(:)=0.           ! reset pdct list for the channel
      iayield=pozyci(ia)                  ! set yield of the current channel
      chyield=0.                          ! reset carbonyl hydroperoxide route yield 
      IF (iayield==0.) CYCLE ialoop       ! no channel if hot criegee yield is 0

      ! identify node # bearing the criegee and carbonyl grp for the current channel
      IF      (ia==1) THEN ; icri=pozcc(i,1) ; icarb=pozcc(i,2)  
      ELSE IF (ia==2) THEN ; icri=pozcc(i,2) ; icarb=pozcc(i,1)
      ENDIF
      
      ! break current POZ and return formula of the CI*, carbonyl and acyl (if any)
      CALL pozprod(chem,bond,group,zebond,icarb,icri, &
                   chemcarb,chemacyl,chemcri,zconform)

      ! rescale decomposition yield for the carbonyl-hydroperoxide route
      chyield=chbratio/2. ;  iayield=iayield*(1.-chbratio)
      IF (iayield<0.) THEN
        mesg="yield below 0 after considering the carbonyl hydroperoxyde channel"
        CALL stoperr(progname,mesg,chem)
      ENDIF

      ! add the carbonyl in the product list (might not exist for cyclic species)
      IF (chemcarb/=' ') THEN  
        brtio=brch*(iayield+chyield)
        CALL add1tonp(progname,chem,np) ; s(np)=iayield+chyield   ! 2 routes: acyl & criegee   
        CALL bratio(chemcarb,brtio,p(np),nrxref(nr),rxref(nr,:))
      ENDIF

      ! add products generated by the hot criegee products in the list
      yci=iayield ; IF (chemacyl==' ') yci=yci+chyield            ! 2 routes if no acyl: diradical & criegee
      CALL hot_criegee(chemcri,yci,zconform,brch,pozatom,loendo, &
                       np,s,p,nrxref(nr),rxref(nr,:))

      ! add product of the carbonyl-hydroperoxide route (ie. acyl species)
      IF ( (chemacyl/=' ') .AND. (chyield/=0.) ) THEN
        CALL add2p(chem,chemacyl,chyield,brch,np,s,p,nrxref(nr),rxref(nr,:))
        CALL add1tonp(progname,chem,np) ; s(np)=chyield ; p(np)='HO '
      ENDIF
      npchannel(ia)=np ; schannel(ia,:)=s(:) ; pchannel(ia,:)=p  

    ENDDO ialoop
    
    ! gather the products of the current reaction (i.e. POZ) into a single list
    CALL gather_pdct(chem,npchannel,schannel,pchannel,np,s,p)

    ! save the list of products for the reaction (i.e. current poz) 
    rxnp(nr)=np ; rxp(nr,:)=p(:) ; rxs(nr,:)=s(:) 
    
    ! save rate constant of current poz for nr reaction
    tarrhc(nr,:)=pozarrhc(i,:)

  ENDDO pozloop
  

! LUMP IDENTICAL REACTION
! -----------------------
  IF (npoz>1) THEN
    poziloop: DO i=1,npoz-1
      pozjloop :DO j=i+1,npoz
        IF (ABS(tarrhc(i,3)-tarrhc(j,3))>0.001*ABS(tarrhc(j,3))) CYCLE pozjloop
        IF (ABS(tarrhc(i,2)-tarrhc(j,2))>0.001*ABS(tarrhc(j,2))) CYCLE pozjloop

        ! check if species in rx i is also in rx j (with same stoi.coef.)
        niloop: DO ni=1,rxnp(i)
          DO nj=1,rxnp(j)
            IF (rxp(i,ni)==rxp(j,nj)) THEN
              IF (ABS(rxs(i,ni)-rxs(j,nj)) < 0.0001)  CYCLE niloop
            ENDIF 
          ENDDO
          CYCLE pozjloop                          ! species in rx i is not in rx j
        ENDDO niloop

        ! check if species in rx j is also in rx i (with same stoi.coef.)
        njloop: DO nj=1,rxnp(j)
          DO ni=1,rxnp(i)
            IF (rxp(j,nj)==rxp(i,ni)) THEN
              IF (ABS(rxs(j,nj)-rxs(i,ni)) < 0.0001)  CYCLE njloop
            ENDIF
          ENDDO
          CYCLE pozjloop                          ! species in rx i is not in rx j
        ENDDO njloop

        ! if that point is reached then rxn i is identical to rxn j
        tarrhc(i,1)=tarrhc(i,1)+tarrhc(j,1) 
        flag(j)=0 ; rxp(j,:)=' ' ; rxs(j,:)=0. ; rxnp(j)=0 ; tarrhc(j,:)=0. 
      ENDDO pozjloop
    ENDDO poziloop
  ENDIF

END SUBROUTINE mko3rx

!=======================================================================
! PURPOSE: gather the products and respective stoi. coef. found in 2 
! lists. Duplicate products can happen within each list and among 
! the 2 lists.   
!=======================================================================
SUBROUTINE gather_pdct(chem,nplist,slist,plist,np,s,p)
  USE toolbox, ONLY: add1tonp
  USE keyflag, ONLY: screenfg
  IMPLICIT NONE
  
  CHARACTER(LEN=*),INTENT(IN)  :: chem       ! chemical formula
  INTEGER,INTENT(IN)           :: nplist(:)  ! # of product in list i
  REAL,INTENT(IN)              :: slist(:,:) ! stoi coef of pdct j in i
  CHARACTER(LEN=*),INTENT(IN)  :: plist(:,:) ! pdct j in list i
  INTEGER,INTENT(OUT)          :: np         ! overall # of distinct pdct
  REAL,INTENT(OUT)             :: s(:)       ! stoi coef of each pdct
  CHARACTER(LEN=*),INTENT(OUT) :: p(:)       ! list of distinct pdct
  
  CHARACTER(LEN=LEN(p(1))) :: pwklist(SIZE(plist,1),SIZE(plist,2)) ! working list 
  REAL    :: swklist(SIZE(slist,1),SIZE(slist,2))                  ! working list 
  INTEGER :: npwklist(SIZE(nplist))                                ! working list 
  INTEGER :: ia,kk,j                                             

  CHARACTER(LEN=15),PARAMETER :: progname='gather_pdct'

  IF (screenfg) WRITE(6,*)progname

  np=0 ; p(:)=' ' ; s(:)=0.

! copy lists into working tables
  npwklist(:)=nplist(:) ; swklist(:,:)=slist(:,:) ; pwklist(:,:)=plist(:,:)
  
! gather duplicate products within each wklist of products
  DO ia=1,2
    DO kk=1,npwklist(ia)-1
      DO j=kk+1,npwklist(ia)
        IF (pwklist(ia,kk)==pwklist(ia,j)) THEN
          swklist(ia,kk)=swklist(ia,kk)+swklist(ia,j)
          pwklist(ia,j)=' '  ;  swklist(ia,j)=0.
        ENDIF
      ENDDO
    ENDDO
  ENDDO
        
! merge the products of the 2 lists into a single list of products
  DO kk=1,npwklist(1)
    IF (pwklist(1,kk)==' ') CYCLE
    DO j=1,npwklist(2)
      IF (pwklist(2,j)==' ') CYCLE
      IF (pwklist(1,kk)==pwklist(2,j)) THEN
        swklist(1,kk)=swklist(1,kk)+swklist(2,j)
        pwklist(2,j)=' ' ; pwklist(2,j)=' '
      ENDIF
    ENDDO
  ENDDO

  DO ia=1,2
    DO kk=1,npwklist(ia)
      IF (pwklist(ia,kk)/=' ') THEN
        CALL add1tonp(progname,chem,np)
        s(np)=swklist(ia,kk)  ;  p(np)=pwklist(ia,kk)
      ENDIF
    ENDDO    
  ENDDO
END SUBROUTINE gather_pdct


END MODULE o3chem
