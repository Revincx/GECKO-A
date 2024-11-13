MODULE ro2chem
IMPLICIT NONE
CONTAINS
! ======================================================================
! PURPOSE: manage the chemistry of peroxy -R(OO.)- radicals.
! First, the database is explored to check if the peroxy provided as
! input is "known". If not, ro2 chemistry is then generated. Reactions
! considered are : (1) ring closure (aromatic only), (2) RO2+NO 
! reaction,  (3) RO2+HO2 reaction, (4) RO2+NO3 reaction, (5) RO2+OH
! reaction and (6) RO2+RO2 reaction.
! Protocol and SARs from Jenkin et al., ACP, 2019
! ======================================================================
SUBROUTINE ro2(idnam,chem,bond,group,nring,brch,cut_off)
  USE keyparameter, ONLY: mxpd,mxnr,mecu,mxcopd
  USE keyflag, ONLY: rx_ro2_oh, highnoxfg, screenfg
  USE references, ONLY:mxlcod
  USE atomtool, ONLY:getatoms,cnum, onum
  USE normchem, ONLY:stdchm
  USE database, ONLY: nkwro2,kwro2_arrh,kwro2_stoi,kwro2_rct,kwro2_prd,kwro2_com
  USE dictstacktool, ONLY: bratio
  USE rxwrttool, ONLY:rxwrit,rxinit
  USE radchktool, ONLY: radchk
  USE toolbox, ONLY: stoperr,add1tonp,addrx
  USE ro2tool, ONLY: ro2no,ro2ho2,ro2no3,ro2ro2,ro2ho,ro2aro,ro2no2
  IMPLICIT NONE

  CHARACTER(LEN=*),INTENT(IN) :: idnam    ! name of current species
  CHARACTER(LEN=*),INTENT(IN) :: chem     ! formula of current species
  CHARACTER(LEN=*),INTENT(IN) :: group(:) ! group (without ring joining char)
  INTEGER,INTENT(IN) :: bond(:,:)         ! bond matrix
  INTEGER,INTENT(IN) :: nring             ! # of rings
  REAL,INTENT(IN)    :: brch              ! max yield of the input species
  REAL,INTENT(IN)    :: cut_off           ! ratio threshold below which rx is ignored

  CHARACTER(LEN=LEN(idnam)) :: pname      ! short name of the parent compound 
  INTEGER :: np,ip,i,j,k,ngr,ncd
  REAL    :: brtio, sctemp
  LOGICAL :: loaro                   ! true if chem is aromatic RO2
  LOGICAL :: loknown                 ! true if RO2 chem is available in the database

  INTEGER :: ic,ih,in,io,ir,is,ifl,icl,ibr,Ncon
  INTEGER :: ncopd, npd2
  REAL    :: scpd2(SIZE(kwro2_prd,2))      ! stoi. coef. delocalisation products 
  CHARACTER(LEN=LEN(chem)) :: delpd2(SIZE(kwro2_prd,2)) ! formula of deloc pdcts

  CHARACTER(LEN=LEN(idnam)) :: coprod(mxcopd)
  REAL                 :: sccopd(mxcopd)

  INTEGER,PARAMETER    :: mxrpd=2   ! max # of products returned by radchk sub
  CHARACTER(LEN=LEN(chem)) :: rdckpd(mxrpd)
  CHARACTER(LEN=LEN(idnam)) :: rdckcopd(mxrpd,mxcopd)
  REAL    :: sc(mxrpd)
  INTEGER :: nip

  CHARACTER(LEN=LEN(idnam)) :: p(mxpd), r(3)
  REAL    :: s(mxpd), arrh(3)
  INTEGER :: idreac, nlabel
  REAL    :: xlabel,folow(3),fotroe(4)

  INTEGER,PARAMETER     :: mxref=10      !--- max # of ref/comment per reaction
  INTEGER               :: nref
  CHARACTER(LEN=mxlcod) :: ref(mxref)
     
  CHARACTER(LEN=8) :: progname='ro2'
  CHARACTER(LEN=70) :: mesg

  IF (screenfg) WRITE(6,*)progname

  loknown=.FALSE.
  pname=idnam

! check if species is allowed in this routine
  IF (INDEX(chem,'(OO.)')==0) THEN
    mesg="Routine for ro2 called but no peroxy found in chem provided as input"
    CALL stoperr(progname,mesg,chem)
  ENDIF
  IF (INDEX(chem,'CO(OO.)')/=0) RETURN


!====== CHECK IF THE REACTION IS KNOWN IN THE DATABASE
  DO i=1,nkwro2
    IF (chem==kwro2_rct(i,1)) THEN
      loknown=.TRUE.

      CALL rxinit(r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe)
      nref=0 ; ref(:)=' '
      DO j=1,SIZE(kwro2_com,2)
        IF (kwro2_com(i,j)/=' ') THEN
          nref=nref+1 ; ref(nref)=kwro2_com(i,j)
        ENDIF
      ENDDO
      np=0 ; arrh(:)=kwro2_arrh(i,:)

      r(1)=pname                                          ! 1st reactant
      IF (kwro2_rct(i,2)/=' ')  r(2)=kwro2_rct(i,2)(1:6)  ! 2nd reactant

! send the products in radchk (only carbon/radicals species)
      ncopd=0 ; npd2=0 ; sccopd(:)=0. ; coprod(:)=' '
      scpd2(:)=0. ; delpd2(:)=' '
      DO j=1,SIZE(kwro2_prd,2)           
        IF (kwro2_prd(i,j)(1:1)==' ') EXIT
        IF (INDEX(kwro2_prd(i,j),'.')==0) CYCLE 
        IF (INDEX(kwro2_prd(i,j),'C')==0) CYCLE 

        sctemp=kwro2_stoi(i,j)                    ! store stoi. coef.
        CALL radchk(kwro2_prd(i,j),rdckpd,rdckcopd,nip,sc,nref,ref)
        kwro2_prd(i,j)=rdckpd(1)                  ! keep the 1st product
        CALL stdchm(kwro2_prd(i,j))
        kwro2_stoi(i,j)= sctemp*sc(1)    ! set coef. by delocalisation fraction
        DO k=1,SIZE(rdckcopd,2)
          IF (rdckcopd(1,k)(1:1)/=' ') THEN
            ncopd=ncopd+1 ; sccopd(ncopd)=sctemp*sc(1)
            coprod(ncopd)=rdckcopd(1,k)
          ENDIF
        ENDDO

        IF (nip==2) THEN            ! if there is a delocalisation product then ...
          npd2=npd2+1               ! add one product
          scpd2(npd2)=sctemp*sc(2)  ! store the initial stochiometric coeff of the product
          CALL stdchm(rdckpd(2))
          delpd2(npd2)=rdckpd(2)    ! store the delocalisation product
          DO k=1,SIZE(rdckcopd,2)
            IF (rdckcopd(2,k)(1:1)/=' ') THEN
              ncopd=ncopd+1 ; sccopd(ncopd)=sctemp*sc(2)
              coprod(ncopd)=rdckcopd(2,k)
            ENDIF
          ENDDO
        ENDIF
      ENDDO

! write products
      DO j=1,SIZE(kwro2_prd,2)
        IF (kwro2_prd(i,j)/=' ') THEN
          CALL add1tonp(progname,chem,np)  
          s(np)=kwro2_stoi(i,j) ; brtio=brch*s(np)
          CALL bratio(kwro2_prd(i,j),brtio,p(np),nref,ref)
        ENDIF 
      ENDDO

! write delocalisation products if any
      DO j=1,npd2
        CALL add1tonp(progname,chem,np) 
        s(np)=scpd2(j) ; brtio=brch*s(np) 
        CALL bratio(delpd2(j),brtio,p(np),nref,ref)
      ENDDO

! write co-products if any
      DO j=1,ncopd
        CALL add1tonp(progname,chem,np) 
        s(np)=sccopd(j) ; p(np)=coprod(j)
      ENDDO

! write reaction
      CALL rxwrit(mecu,r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe,ref)

    ENDIF
  ENDDO
  IF (loknown) RETURN

!====== MAKE CHEMISTRY OF THE RO2

! locate peroxy group and count number of "Cds"
  ncd=0 ; ngr=COUNT(group/=' ')
  DO i=1,ngr
    IF (INDEX(group(i),'(OO.)')/=0)  ip=i
    IF (group(i)(1:2)=='Cd')  ncd=ncd+1
  ENDDO
  CALL getatoms(chem,ic,ih,in,io,ir,is,ifl,ibr,icl)
  Ncon=ic+io-2+in
  
! aromatic RO2 - ring closure
! ----------------------------
  IF ((ncd>3) .AND. (nring/=0)) THEN
    CALL ro2aro(chem,pname,bond,group,ngr,brch,ip,loaro)
    IF (loaro) RETURN
  ENDIF

! reaction with NO
! -------------------
  CALL ro2no(chem,pname,bond,group,ngr,nring,brch,ip,Ncon)
  IF (highnoxfg) RETURN   ! mechanism for high NOx only

! reaction with NO2
! -------------------
  CALL ro2no2(chem,pname,bond,group,brch,ip)

! reaction with HO2
! -------------------
  CALL ro2ho2(chem,pname,bond,group,ngr,brch,ip,Ncon)

! reaction with NO3
! -------------------
  CALL ro2no3(chem,pname,bond,group,brch,ip)

! reaction with OH
! -------------------
  IF (rx_ro2_oh)  CALL ro2ho(chem,pname,bond,group,brch,ip)

! RO2+RO2 reactions
! -------------------------------
  CALL ro2ro2(chem,pname,bond,group,ngr,brch,ip,Ncon)

END SUBROUTINE ro2
END MODULE ro2chem
