MODULE rochem
IMPLICIT NONE
CONTAINS
!SUBROUTINE ro(idnam,chem,bond,group,nring,brch,cut_off)
!SUBROUTINE roselector(cut_off,chem,&
!                      arrhest,arrho2,arrhno2,arrhiso,FSD,arrhdec,&
!                      flest,flo2,flno2,niso,fliso,ndec,fldec,&
!                      ksum,kest,ko2,kiso,kdec)

! ======================================================================
! PURPOSE: manage the chemistry of alkoxy -R(O.)- radicals.
! Four reaction types are considered: (1) reaction with O2, 
! (2) decomposition, (3) H-shift isomerisation and (4) ester 
! rearrangement. Some particular cases are considered first, for
! structure having an expected "unique" (dominant) pathways.
! ======================================================================
SUBROUTINE ro(idnam,chem,bond,group,nring,rjg,brch,cut_off)
  USE keyparameter, ONLY: mxpd,mxring,mecu,refu,mxcopd
  USE references, ONLY:mxlcod
  USE keyflag, ONLY: wrtopeinfo,screenfg,kisom_sar,rx_ro_no2,wrtref    ! select SAR for reaction rate  
  USE database, ONLY: nkwro, kwro_arrh,kwro_stoi,kwro_rct,kwro_prd,kwro_com
  USE dictstacktool, ONLY: bratio
  USE rxwrttool, ONLY:rxwrit,rxinit
  USE roisotool, ONLY: roiso
  USE rorxtool, ONLY: rooxy,roester,rodecnitro,roaro,rono2
  USE rodectool, ONLY: rospedec, rodec, rodecO
  USE toolbox, ONLY: stoperr,add1tonp,addref
  IMPLICIT NONE
      
  CHARACTER(LEN=*),INTENT(IN) :: idnam    ! name of current species
  CHARACTER(LEN=*),INTENT(IN) :: chem     ! formula of current species
  CHARACTER(LEN=*),INTENT(IN) :: group(:) ! group (without ring joining char)
  INTEGER,INTENT(IN) :: bond(:,:)         ! bond matrix
  INTEGER,INTENT(IN) :: nring             ! # of rings
  INTEGER,INTENT(IN) :: rjg(:,:)          ! ring joining groups (if any)
  REAL,INTENT(IN)    :: brch              ! max yield of the input species
  REAL,INTENT(IN)    :: cut_off           ! ratio threshold below which rx is ignored

  CHARACTER(LEN=LEN(group(1))) :: tgroup(SIZE(group))
  INTEGER                      :: tbond(SIZE(bond,1),SIZE(bond,2))
  LOGICAL ::  loknown                     ! .TRUE. if RO chem is available in the database
  LOGICAL ::  loregro                     ! .TRUE. if regular "RO" (all reactions considered)

  INTEGER,PARAMETER    :: mxref=10        !--- max # of refs/comments per reaction
  INTEGER              :: flest           !--- ESTER CHANNEL, flag
  CHARACTER(LEN=LEN(chem)) :: pdest(2)            ! main products for ester channel
  CHARACTER(LEN=LEN(idnam)) :: copdest(mxcopd)    ! coproduct, ester channel
  REAL                 :: arrhest(3), kest   ! arrhenius coef. & rate constant, ester channel
  INTEGER              :: nrfest             ! # of references in the reference list
  CHARACTER(LEN=mxlcod):: rfest(mxref)       ! ref/comment, ester channel
  INTEGER              :: flo2         !--- O2 CHANNEL, flag 
  CHARACTER(LEN=LEN(chem)) :: po2(2)           ! main products for O2 channel
  REAL                 :: ypo2(2)              ! yield (ratio) of each po2 (if > 2)
  CHARACTER(LEN=LEN(idnam)) :: copdo2(mxcopd)  ! coproduct, O2 channel
  REAL                 :: ycopdo2(mxcopd)      ! yield of each coproduct in O2 channel
  REAL                 :: arrho2(3),ko2        ! arrhenius coef. & rate constant, O2 channel
  INTEGER              :: nrfo2                ! # of references in the reference list
  CHARACTER(LEN=mxlcod):: rfo2(mxref)          ! ref/comment, O2 channel
  INTEGER              :: flno2        !--- NO2 CHANNEL, flag
  CHARACTER(LEN=LEN(chem)) :: pno2(2)          ! product(s) for NO2 channel
  REAL                 :: ypno2(2)             ! yield (ratio) of each pno2 (if > 2)
  CHARACTER(LEN=LEN(idnam)) :: copdno2(mxcopd) ! coproduct, NO2 channel
  REAL                 :: ycopdno2(mxcopd)     ! yield of each coproduct in NO2 channel
  REAL                 :: arrhno2(3),kno2      ! arrhenius coef. & rate constant, NO2 channel
  INTEGER              :: nrfno2               ! # of references in the reference list
  CHARACTER(LEN=mxlcod):: rfno2(mxref)         ! ref/comment, NO2 channel
  INTEGER,PARAMETER    :: mxiso=6      !--- ISOMERISATION, max # of reactions 
  INTEGER              :: fliso(mxiso)         ! flag for isomerisation
  INTEGER              :: fliso2(mxiso)        ! 2nd flag for isomerisation (2 pdcts)
  INTEGER              :: niso                 ! # of isomerisation reactions
  CHARACTER(LEN=LEN(chem)) :: piso(mxiso)          ! pdct of the isomerisation reactions
  CHARACTER(LEN=LEN(chem)) :: piso2(mxiso)         ! 2nd pdct of the isomerisation reactions
  CHARACTER(LEN=LEN(idnam)) :: copiso(mxiso,mxcopd) ! co-pdct of the product
  CHARACTER(LEN=LEN(idnam)) :: copiso2(mxiso,mxcopd)! co-pdct for the 2nd product
  REAL                 :: siso(mxiso,2)        ! ratio of 1st vs 2nd pdct
  REAL                 :: arrhiso(mxiso,3)     ! arrhenius coef., isomerisation
  REAL                 :: FSD(mxiso,5)         ! polynomial coef. for Vereecken SAR
  REAL                 :: kiso(mxiso)          ! rate constant, isomerisation
  INTEGER              :: nrfiso(mxiso)        ! # of ref/comment in the list (isomerisation)
  CHARACTER(LEN=mxlcod):: rfiso(mxiso,mxref)   ! ref/comment, isomerisation channel
  INTEGER,PARAMETER    :: mxdec=5      !--- DECOMPOSITION, max # of  reactions 
  INTEGER              :: fldec(mxdec)         ! flag for decomposition
  INTEGER              :: fldec2(mxdec)        ! 2nd flag for decomposition
  INTEGER              :: ndec                 ! # of decomposition reactions
  CHARACTER(LEN=LEN(chem)) :: pdec(mxdec,2)        ! pdct, need size 2 (2 dictinct pdcts for acyclic molecules)
  CHARACTER(LEN=LEN(chem)) :: pdec2(mxdec)         ! 2nd pdct
  CHARACTER(LEN=LEN(idnam)) :: copdec(mxdec,mxcopd) ! co-pdct of the product
  CHARACTER(LEN=LEN(idnam)) :: copdec2(mxdec,mxcopd)! co-pdct of the 2nd product
  REAL                 :: sdec(mxdec,2)        ! ratio of 1st vs 2nd pdct
  REAL                 :: arrhdec(mxdec,3)     ! arrhenius coef., decomposition
  REAL                 :: kdec(mxdec)          ! rate constant, decomposition
  INTEGER              :: nrfdec(mxdec)        ! # of ref/comment in the list (decomposition)
  CHARACTER(LEN=mxlcod):: rfdec(mxdec,mxref)  ! ref/comment, isomerisation channel

  INTEGER              :: nref
  CHARACTER(LEN=mxlcod):: ref(mxref)
  
  INTEGER :: nca,np,ia,i,j,ib,sumflg
  REAL    :: brtio
  REAL    :: ksum

  CHARACTER(LEN=LEN(idnam)) :: p(mxpd), r(3)
  REAL    :: s(mxpd), arrh(3)
  INTEGER :: idreac, nlabel
  REAL    :: xlabel,folow(3),fotroe(4)

  CHARACTER(LEN=5),PARAMETER :: progname='ro '
  CHARACTER(LEN=70)          :: mesg

  IF (screenfg) WRITE(6,*)progname

! check if species is allowed in this routine
  IF (LEN_TRIM(chem)<1) RETURN
  IF (INDEX(chem,'(O.)')==0) THEN
    mesg="The routine was called with no (O.) grp "
    CALL stoperr(progname,mesg,chem)
  ENDIF

! initialize reaction line and transcribe groups and bonds
  flest=0 ; pdest(:)=' ' ; copdest(:)=' ' ; arrhest(:)=0. ; kest=0.    ! ester
  nrfest=0 ; rfest(:)=' '                                              ! ester
  flo2=0  ; po2(:)=' ' ; ypo2(:)=0. ; copdo2(:)=' ' ; ycopdo2(:)=0.    ! O2
  arrho2(:)=0. ; ko2=0. ; nrfo2=0 ; rfo2(:)=' '                        ! O2
  flno2=0  ; pno2(:)=' '                                               ! NO2
  arrhno2(:)=0. ; kno2=0. ; nrfno2=0 ; rfno2(:)=' '                    ! NO2
  fliso(:)=0 ; fliso2(:)=0 ; niso=0 ; piso(:)=' ' ; piso2(:)=' '       ! isomerisation
  copiso(:,:)=' ' ; copiso2(:,:)=' ' ; siso(:,:)=0. ; arrhiso(:,:)=0.  ! isomerisation
  kiso(:)=0. ; FSD(:,:)=0. ; nrfiso(:)=0 ; rfiso(:,:)=' '              ! isomerisation
  fldec(:)=0 ; fldec2(:)=0 ; ndec=0 ; pdec(:,:)=' ' ; pdec2(:)=' '     ! decomposition
  copdec(:,:)=' ' ; copdec2(:,:)=' ' ; sdec(:,:)=0. ; arrhdec(:,:)=0.  ! decomposition   
  kdec(:)=0 ; nrfdec(:)=0  ;  rfdec(:,:)=' '                           ! decomposition
  loregro=.TRUE.  ! default is "regular" RO
  
! ----------------------------------------------
! check if the reaction is known in the database
! ----------------------------------------------
  loknown=.FALSE.              ! default is unknown reactions !
  DO i=1,nkwro
    IF (chem==kwro_rct(i,1)) THEN
      loknown = .TRUE.
      CALL rxinit(r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe)
      nref=0 ; ref(:)=' '
      DO j=1,SIZE(kwro_com,2)
        IF (kwro_com(i,j)/=' ') THEN
          CALL addref(progname,kwro_com(i,j),nref,ref,chem)
        ENDIF
      ENDDO
 
! read reactants & reaction rate
      r(1) = idnam 
      IF (kwro_rct(i,2)/=' ')  r(2) = kwro_rct(i,2)(1:6)
      arrh(1:3)=kwro_arrh(i,1:3)

! read products
      s(1) = kwro_stoi(i,1)          ! 1st product
      brtio = brch*s(1)
      CALL bratio(kwro_prd(i,1),brtio,p(1),nref,ref)
      np=1

      IF (kwro_prd(i,2)/=' ') THEN   ! 2nd product
        np=np+1  ; s(np) = kwro_stoi(i,2)
        brtio = brch*s(np)
        CALL bratio(kwro_prd(i,2),brtio,p(np),nref,ref)
      ENDIF

      IF (kwro_prd(i,3)/=' ') THEN   ! 3rd product
        np=np+1 ; s(np) = kwro_stoi(i,3)
        brtio = brch * s(np)
        CALL bratio(kwro_prd(i,3),brtio,p(np),nref,ref)
      ENDIF

! write reaction
      CALL rxwrit(mecu,r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe,ref)

    ENDIF
  ENDDO
  IF (loknown) RETURN ! exit if reactions are known

! -----------------
! set and copy working tables
! -----------------
  tgroup(:) = group(:)  ;  tbond(:,:) = bond(:,:)                   ! cp bond & groups
  nca=COUNT(tgroup/=' ')
  DO i=1,nca                                                        ! identify the R(O.) node
    IF (INDEX(group(i),'(O.)')/=0) THEN
      ia = i ; EXIT                             
    ENDIF    
  ENDDO

! -----------------
! Particular cases:
! -----------------

! phenoxy radicals: write reaction and return
  IF (group(ia)(1:1)=='c'  ) THEN
    CALL roaro(idnam,bond,group,nca,ia)
    RETURN
  ENDIF

! If R(NO2)(O.): decompose to R=O + NO2 only
  IF (INDEX(group(ia),'(NO2)')/=0) THEN
    CALL rodecnitro(bond,group,ia,ndec,fldec,pdec,copdec,arrhdec,nrfdec,rfdec)
    sdec(ndec,:)=1.
    loregro=.FALSE.   ! jump to reaction rates
  ENDIF
      
! If  >C=CR-CR(polar)-C(O.)< or >C=CR-CO-C(O.) decompose only
  IF (loregro) THEN
    DO ib=1,nca 
      IF (bond(ib,ia)/=0) THEN
        DO  j=1,nca
          IF ((bond(ib,j)==1).AND.(group(j)(1:2)=='Cd').AND.(j/=ia)) THEN
            IF ((INDEX(tgroup(ib),'(OH)')/=0).OR. &
               (INDEX(tgroup(ib),'(ONO2)')/=0).OR. &
               (INDEX(tgroup(ib)(1:3),'CO ')/=0).OR. &
               (INDEX(tgroup(ib),'(OOH)')/=0)) THEN
               CALL rospedec(bond,group,nca,ia,ib,nring,ndec,fldec,pdec,copdec,&
                             fldec2,pdec2,copdec2,sdec,arrhdec,nrfdec,rfdec)
               loregro=.FALSE.  ! jump to reaction rates
            ENDIF
          ENDIF
        ENDDO
      ENDIF
    ENDDO
  ENDIF
  
! --------------------------
! alpha-ester rearrangement 
! --------------------------
  IF (loregro) THEN
    CALL roester(bond,group,nca,ia,nring,flest,pdest,copdest,arrhest,&
                 loregro,nrfest,rfest)
  ENDIF

  IF (loregro) THEN
! -------------------
! reaction with O2 
! -------------------
    CALL rooxy(chem,bond,group,nca,ia,nring,flo2,po2,ypo2,copdo2, &
               ycopdo2,arrho2,nrfo2,rfo2)

! -------------------
! reaction with NO2
! -------------------
    IF(rx_ro_no2) THEN
    CALL rono2(chem,bond,group,nca,ia,nring,flno2,pno2,ypno2,copdno2, &
               ycopdno2,arrhno2,nrfno2,rfno2)
    ENDIF

! -------------------
! isomerization (for non-ring nodes only):
! -------------------
    CALL roiso(bond,group,nca,ia,nring,rjg, &
           niso,fliso,piso,copiso,fliso2,piso2,copiso2,siso, &
           arrhiso,FSD,nrfiso,rfiso)

! -----------------------------
! decomposition / ring opening
! -----------------------------
    CALL rodec(chem,bond,group,nca,ia,nring,ndec,fldec,pdec,copdec, &
               fldec2,pdec2,copdec2,sdec,arrhdec,nrfdec,rfdec)

! -----------------------------
! if no path, search for other decomposition (e.g. C(O.)(OH) -> C=O + OH)
! -----------------------------
    sumflg=flest+flo2+SUM(fliso(1:niso))+SUM(fldec(1:ndec))
    IF (sumflg==0) THEN  
      CALL rodecO(bond,group,nca,ia,ndec,fldec,pdec,copdec,sdec,arrhdec,nrfdec,rfdec)
    ENDIF
  ENDIF

! ------------------------------------------------
! compare rate constants, keep only if > cutoff
! ------------------------------------------------
  CALL roselector(cut_off,chem,brch,&
                      arrhest,arrho2,arrhno2,arrhiso,FSD,arrhdec,&
                      flest,flo2,flno2,niso,fliso,ndec,fldec,&
                      ksum,kest,ko2,kno2,kiso,kdec)
  IF (wrtopeinfo) WRITE(10,*)'***********! alkoxy ','G',idnam,'********'

! -------------------
! WRITE ALL REACTIONS
! -------------------

!  alpha ester rearrangement
! ---------------------------
  IF (flest==1) THEN
    CALL rxinit(r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe)
    brtio = brch*kest/ksum
!    brtio = brch

! 1st & 2nd products 
    s(1)=1.  ;  CALL bratio(pdest(1),brtio,p(1),nrfest,rfest)
    s(2)=1.  ;  CALL bratio(pdest(2),brtio,p(2),nrfest,rfest)
    IF (wrtopeinfo)  WRITE(10,21)  'ester rearrgt /',brch,'G',p(1),' + ','G',p(2)

! other products
    np=2
    DO j=1,mxcopd      
      IF (copdest(j)(1:1)/=' ') THEN
        CALL add1tonp(progname,chem,np)  
        s(np)=1.  ;  p(np)=copdest(j)
        IF (wrtopeinfo) WRITE(10,'(20X,A1,A6)') 'G',p(np)
      ENDIF
    ENDDO

! reactant
    r(1) = idnam
    arrh(1:3) = arrhest(1:3)
    CALL rxwrit(mecu,r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe,rfest)
  ENDIF

! O2 reaction:
! -------------
  IF (flo2==1) THEN
    np=0
    CALL rxinit(r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe)
    
    ! add the expected carbonyl
    CALL add1tonp(progname,chem,np) ; s(np) = ypo2(1)
    brtio = (brch*ko2/ksum)*ypo2(1)
!    brtio = brch*ypo2(1)
    CALL bratio(po2(1),brtio,p(np),nrfo2,rfo2)
    IF (wrtopeinfo) WRITE(10,22)'O2 -> HO2 /',ko2/ksum,'G',p(1)

    ! add the second product (if any)
    IF (po2(2)/=' ') THEN
      CALL add1tonp(progname,chem,np) ; s(np) = ypo2(2)
      brtio = (brch*ko2/ksum)*ypo2(2)
!      brtio = brch*ypo2(2)
      CALL bratio(po2(2),brtio,p(np),nrfo2,rfo2)
    ENDIF
    
    ! add HO2 
    CALL add1tonp(progname,chem,np) ; s(np) = 1.  ;  p(np) = 'HO2  '

    ! add other products if any)
    copdloop: DO j=1,SIZE(copdo2)      
      IF (copdo2(j)==' ') EXIT copdloop
      CALL add1tonp(progname,chem,np)  
      s(np)=ycopdo2(j)  ;  p(np)=copdo2(j)
      IF (wrtopeinfo) WRITE(10,'(20X,A1,A6)') 'G',p(np)
    ENDDO copdloop

    r(1)=idnam  ;  r(2)='OXYGEN'
    arrh(:)=arrho2(:)

! write out - extra reaction => idreac=2 (nlabel was set above)
    idreac=0
    CALL rxwrit(mecu,r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe,rfo2)
  ENDIF

! NO2 reaction (usually: addition, single product, no coproducts):
!              ... But R-C(=O)ONO2 decomposes =>  R. + CO2 + NO2
! -------------
  IF (flno2==1) THEN
    np=0
    CALL rxinit(r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe)

    ! add the expected nitrate
    CALL add1tonp(progname,chem,np) ; s(np) = 1
    brtio = (brch*kno2/ksum)*ypno2(1)
    CALL bratio(pno2(1),brtio,p(np),nrfno2,rfno2)
    IF (wrtopeinfo) WRITE(10,22)'NO2 addition /',kno2/ksum,'G',p(1)

    ! add other products if any)
    copdno2loop: DO j=1,SIZE(copdno2)
      IF (copdno2(j)==' ') EXIT copdno2loop
      CALL add1tonp(progname,chem,np)
      s(np)=ycopdno2(j)  ;  p(np)=copdno2(j)
      IF (wrtopeinfo) WRITE(10,'(20X,A1,A6)') 'G',p(np)
    ENDDO copdno2loop

    r(1)=idnam  ;  r(2)='NO2'
    arrh(:)=arrhno2(:)

! write out - extra reaction => idreac=2 (nlabel was set above)
    idreac=0
    CALL rxwrit(mecu,r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe,rfno2)
  ENDIF

! isomerization reactions
! ----------------------- 
  irxloop: DO i=1,niso
    IF (fliso(i)==0) CYCLE irxloop
    CALL rxinit(r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe)
    s(1) = 1.
    brtio = brch*kiso(i)/ksum
!    brtio=brch
    CALL bratio(piso(i),brtio,p(1),nrfiso(i),rfiso(i,:))
    IF (wrtopeinfo)  WRITE(10,22) 'isomerisa° /',kiso(i)/ksum,'G',p(1)

! second product (if delocalisation allowed)
    np=1
    IF (fliso2(i)/=0) THEN
      np=np+1
      s(1)=siso(i,1)  ;  s(np)=siso(i,2)
      CALL bratio(piso2(i),brtio,p(np),nrfiso(i),rfiso(i,:))
    ENDIF

! other products
    DO j=1,mxcopd
      IF (copiso(i,j)(1:1)/=' ') THEN
        CALL add1tonp(progname,chem,np)  
        s(np)=1.  ;  p(np)=copiso(i,j)
        IF (wrtopeinfo) WRITE(10,'(20X,A1,A6)') 'G',p(np)
      ENDIF
    ENDDO

    r(1)=idnam
    IF (kisom_sar==2) r(2)='ISOM '
    arrh(1)=arrhiso(i,1) ; arrh(2)=arrhiso(i,2) ; arrh(3)=arrhiso(i,3)

    CALL rxwrit(mecu,r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe,rfiso(i,:))

    IF (kisom_sar==2) THEN
      WRITE(mecu,'(A9,5(ES10.3,1x),A1)')            "  ISOM / ",(FSD(i,j),j=1,5),"/"
      IF (wrtref) WRITE(refu,'(A9,5(ES10.3,1x),A1)')"  ISOM / ",(FSD(i,j),j=1,5),"/"
    ENDIF
  ENDDO irxloop

! decomposition reaction
! -----------------------
  drxloop: DO i=1,ndec
    IF (fldec(i)==0) CYCLE drxloop
    CALL rxinit(r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe)
      
! first & other products
    s(1)=1.
    brtio=brch*kdec(i)/ksum
!    brtio=brch
    CALL bratio(pdec(i,1),brtio,p(1),nrfdec(i),rfdec(i,:))
    np=1

    IF (INDEX(pdec(i,2),' ')>1) THEN  ! 2nd pdct
      np=np+1  ;  s(np)=1.
      CALL bratio(pdec(i,2),brtio,p(np),nrfdec(i),rfdec(i,:))
      IF (wrtopeinfo) &
        WRITE(10,21)'decomp /      ',kdec(i)/ksum,'G',p(1),' + ','G',p(2)
    ELSE
      IF (wrtopeinfo) &
        WRITE(10,21)'decomp /      ',kdec(i)/ksum,'G',p(1)
    ENDIF

    IF (fldec2(i)/=0) THEN        ! additional pdcts from delocalisation
      IF (INDEX(pdec(i,1),'.')/=0) THEN ; s(1) =sdec(i,1)
      ELSE                              ; s(np)=sdec(i,1)
      ENDIF
      np=np+1  ;  s(np)=sdec(i,2)
      CALL bratio(pdec2(i),brtio,p(np),nrfdec(i),rfdec(i,:))
    ENDIF

! other (coproducts)
    DO j=1,mxcopd
      IF (copdec(i,j)(1:1)/=' ') THEN
        CALL add1tonp(progname,chem,np)  
        s(np)=sdec(i,1)  ;  p(np)=copdec(i,j)
        IF (wrtopeinfo) WRITE(10,'(20X,A1,A6)') 'G',p(np)
      ENDIF
    ENDDO
    DO j=1,mxcopd
      IF ((copdec2(i,j)(1:1)/=' ').AND.(fldec2(i)/=0)) THEN
        CALL add1tonp(progname,chem,np)  
        s(np)=sdec(i,2)  ;  p(np)=copdec2(i,j)
        IF (wrtopeinfo) WRITE(10,'(20X,A1,A6)') 'G',p(np)
      ENDIF
    ENDDO

    r(1)=idnam
    arrh(:)=arrhdec(i,:)
    CALL rxwrit(mecu,r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe,rfdec(i,:))

  ENDDO drxloop
  
!  IF (wrtopeinfo) WRITE(6,*)'end of ro'

21  FORMAT(A15,5X,f5.3,2X,A1,A6,A3,A1,A6)         
22  FORMAT(A14,4X,f5.3,2X,A1,A6)         
END SUBROUTINE ro

!=======================================================================
! PURPOSE: select reaction pathways for RO chemistry. 1st: all pathways
! below threshold are killed. If no reaction remains, keep only major
! reaction (> kmax/2, where kmax is the highest rate constant). 
!=======================================================================
SUBROUTINE roselector(cut_off,chem,brch,  &
                      arrhest,arrho2,arrhno2,arrhiso,FSD,arrhdec,&
                      flest,flo2,flno2,niso,fliso,ndec,fldec,&
                      ksum,kest,ko2,kno2,kiso,kdec)
  USE keyflag, ONLY: TK,Mref,C_NO2,kisom_sar,yldcut        ! select SAR for reaction rate
  USE toolbox, ONLY: stoperr
  IMPLICIT NONE
  
  CHARACTER(LEN=*),INTENT(IN):: chem
  REAL, INTENT(IN) :: cut_off
  REAL,INTENT(IN)  :: brch              ! max yield of the input species
  REAL, INTENT(IN) :: arrhest(:)
  REAL, INTENT(IN) :: arrho2(:)
  REAL, INTENT(IN) :: arrhno2(:)
  REAL, INTENT(IN) :: arrhiso(:,:)
  REAL, INTENT(IN) :: FSD(:,:)
  REAL, INTENT(IN) :: arrhdec(:,:)
  INTEGER,INTENT(INOUT) :: flest
  INTEGER,INTENT(INOUT) :: flo2
  INTEGER,INTENT(INOUT) :: flno2
  INTEGER,INTENT(IN)    :: niso
  INTEGER,INTENT(INOUT) :: fliso(:)
  INTEGER,INTENT(IN)    :: ndec
  INTEGER,INTENT(INOUT) :: fldec(:) 
  REAL,INTENT(OUT)      :: ksum
  REAL,INTENT(OUT)      :: kest
  REAL,INTENT(OUT)      :: ko2
  REAL,INTENT(OUT)      :: kno2
  REAL,INTENT(OUT)      :: kiso(:)
  REAL,INTENT(OUT)      :: kdec(:)

  INTEGER :: sumflg,i
  REAL    :: kmax,kcut1,kcut2

  CHARACTER(LEN=12),PARAMETER :: progname='roselector '
  CHARACTER(LEN=70)           :: mesg

  ksum=0. 
  kest=0 ; ko2=0. ; kiso(:)=0. ; kdec(:)=0. 

! compute rate constant for ester rearrangement
  IF (flest==1) THEN
    kest = arrhest(1)*(TK**arrhest(2))*exp(-arrhest(3)/TK)
    ksum = ksum + kest
  ENDIF
  
! compute rate constant for O2 reaction
  IF (flo2==1) THEN
    ko2 = arrho2(1)*(TK**arrho2(2))*exp(-arrho2(3)/TK)*Mref*0.209
    ksum = ksum + ko2
  ENDIF 
  
! compute rate constant for NO2 reaction
  IF (flno2==1) THEN
    kno2 = arrhno2(1)*(TK**arrhno2(2))*exp(-arrhno2(3)/TK)*C_NO2
    ksum = ksum + kno2
  ENDIF 
  
! compute rate constant for isomerization reaction
  DO i=1,niso
    IF (fliso(i)==1) THEN
      kiso(i) = arrhiso(i,1)*(TK**arrhiso(i,2))*exp(-arrhiso(i,3)/TK)
      IF (kisom_sar==2) THEN
        kiso(i)=kiso(i)*(FSD(i,1)*TK**4 + FSD(i,2)*TK**3 + &
                FSD(i,3)*TK**2 + FSD(i,4)*TK + FSD(i,5))
      ENDIF
      ksum = ksum + kiso(i)
    ENDIF
  ENDDO

! compute rate constant for decomposition reaction
  DO i=1,ndec
    IF (fldec(i)==1) THEN
      kdec(i) = arrhdec(i,1)*(TK**arrhdec(i,2))*exp(-arrhdec(i,3)/TK)
      ksum = ksum + kdec(i)
   ENDIF
  ENDDO
  
! remove reaction below threshold
  kcut1 = cut_off*ksum
  IF (kest<kcut1) flest=0
  IF (ko2<kcut1)  flo2=0
  IF (kno2<kcut1) flno2=0
  DO i=1,niso  ;  IF (kiso(i)<kcut1) fliso(i)=0  ;  ENDDO
  DO i=1,ndec  ;  IF (kdec(i)<kcut1) fldec(i)=0  ;  ENDDO

! remove reactions with overall production yield below threshold
  kcut2 = (yldcut/brch)*ksum
  IF (kest<kcut2) flest=0
  IF (ko2<kcut2)  flo2=0
  IF (kno2<kcut2) flno2=0
  DO i=1,niso  ;  IF (kiso(i)<kcut2) fliso(i)=0  ;  ENDDO
  DO i=1,ndec  ;  IF (kdec(i)<kcut2) fldec(i)=0  ;  ENDDO

! Check that a reaction remains
  sumflg=flest+flo2+flno2+SUM(fliso(1:niso))+SUM(fldec(1:ndec))

! If no reaction remains, find the fatest and keep major pathways
  IF (sumflg==0) THEN
    kmax=MAX(kest,ko2)
    kmax=MAX(kmax,kno2)
    kmax=MAX(kmax,MAXVAL(kiso(1:niso)))
    kmax=MAX(kmax,MAXVAL(kdec(1:ndec)))
    
    IF (kest>kmax/2.) flest=1
    IF (ko2>kmax/2.)  flo2=1
    IF (kno2>kmax/2.) flno2=1
    WHERE (kiso(:)>kmax/2.) fliso(:)=1
    WHERE (kdec(:)>kmax/2.) fldec(:)=1
  ENDIF

! If no reaction still remains, stop error
  sumflg=flest+flo2+flno2+SUM(fliso(1:niso))+SUM(fldec(1:ndec))
  IF (sumflg==0) THEN
    mesg="no reaction path found"
    CALL stoperr(progname,mesg,chem)
  ENDIF
  
! compute sum of the rate constant at TK
  ksum=0.
  IF (flest==1) ksum=ksum+kest
  IF (flo2==1)  ksum=ksum+ko2
  IF (flno2==1) ksum=ksum+kno2
  DO i=1,niso
    IF (fliso(i)==1) ksum=ksum+kiso(i)
  ENDDO
  DO i=1,ndec
    IF (fldec(i)==1) ksum=ksum+kdec(i)
  ENDDO
END SUBROUTINE roselector


END MODULE rochem
