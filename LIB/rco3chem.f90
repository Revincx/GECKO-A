MODULE rco3chem
IMPLICIT NONE
CONTAINS

! ======================================================================
! PURPOSE: manage the chemitsry of peroxy -RCO(OO.)- radicals.
! First, the database is explored to check if the peroxy acyl provided 
! as input is "known". If not, rco3 chemistry is then generated. 
! Reactions considered are : (1) RCO3+NO reaction, 
! (2) RCO3+NO2 reaction, (3) RO2+NO3 reaction, (4) RO2+HO2 reaction, 
! (5) RCO3+RO2 reaction.
! Protocol and SARs from Jenkin et al., ACP, 2019
! ======================================================================
SUBROUTINE rco3(idnam,chem,bond,group,brch,cut_off)
  USE keyparameter, ONLY: mxpd,mxnr,mecu,mxcopd
  USE references, ONLY:mxlcod
  USE keyflag, ONLY: highnoxfg,screenfg
  USE database, ONLY: nkwrco3,kwrco3_arrh,kwrco3_stoi,kwrco3_rct,kwrco3_prd,kwrco3_com
  USE dictstacktool, ONLY: bratio
  USE normchem, ONLY: stdchm
  USE rxwrttool, ONLY:rxwrit,rxinit
  USE radchktool, ONLY: radchk
  USE rco3tool, ONLY: rco3no3,rco3no,rco3no2,rco3ho2,rco3ro2
  USE toolbox, ONLY: stoperr,add1tonp,addrx,addref
  IMPLICIT NONE

  CHARACTER(LEN=*),INTENT(IN) :: idnam    ! name of current species
  CHARACTER(LEN=*),INTENT(IN) :: chem     ! formula of current species
  CHARACTER(LEN=*),INTENT(IN) :: group(:) ! group (without ring joining char)
  INTEGER,INTENT(IN) :: bond(:,:)         ! bond matrix
  REAL,INTENT(IN)    :: brch              ! max yield of the input species
  REAL,INTENT(IN)    :: cut_off           ! ratio threshold below which rx is ignored

  CHARACTER(LEN=LEN(idnam)) :: pname      ! short name of the parent compound 

  INTEGER :: ip,i,j,k,ngr,np
  REAL    :: brtio

  CHARACTER(LEN=LEN(idnam)) :: p(mxpd), r(3)
  REAL    :: s(mxpd), arrh(3)
  INTEGER :: idreac, nlabel
  REAL    :: xlabel,folow(3),fotroe(4)

  INTEGER :: ncopd, npd2
  LOGICAL :: loknow                       ! true if RCO3 chem is available in database
  REAL    :: scpd2(SIZE(kwrco3_prd,2))    ! stoi. coef. delocalisation products 
  REAL    :: sctemp
  CHARACTER(LEN=LEN(chem)) :: delpd2(SIZE(kwrco3_prd,2)) ! formula of deloc pdcts

  CHARACTER(LEN=LEN(idnam)) :: coprod(mxcopd)
  REAL                 :: sccopd(mxcopd)

  INTEGER,PARAMETER    :: mxrpd=2   ! max # of products returned by radchk sub
  CHARACTER(LEN=LEN(chem)) :: rdckpd(mxrpd)
  CHARACTER(LEN=LEN(idnam)):: rdckcopd(mxrpd,mxcopd)
  REAL    :: sc(mxrpd)
  INTEGER :: nip

  CHARACTER(LEN=7)  :: progname='rco3'
  CHARACTER(LEN=70) :: mesg
  INTEGER,PARAMETER     :: mxref=10      !--- max # of ref/comment per reaction
  INTEGER               :: nref
  CHARACTER(LEN=mxlcod) :: ref(mxref)

  IF (screenfg) WRITE(6,*)progname

!  nref=0 ; com(:)=' '
  p(:)=' ' ; s(:)=0.
  loknow=.FALSE.
  pname=idnam

! check if species is allowed in this routine
  IF (INDEX(chem,'CO(OO.)')==0) THEN
    mesg="Routine for rco3 called but not found in chem provided as input"
    CALL stoperr(progname,mesg,chem)
  ENDIF

!====== CHECK IF THE REACTION IS KNOWN IN THE DATABASE

  DO i=1,nkwrco3
    IF (chem==kwrco3_rct(i,1)) THEN
      loknow=.TRUE.

      CALL rxinit(r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe)
      nref=0 ; ref(:)=' '
      DO j=1,SIZE(kwrco3_com,2)
        IF (kwrco3_com(i,j)/=' ') THEN
          CALL addref(progname,kwrco3_com(i,j),nref,ref,chem)
        ENDIF
      ENDDO

      np=0  ; arrh(:)=kwrco3_arrh(i,:)

      r(1)=pname                                             ! 1st reactant
      IF (kwrco3_rct(i,2)/=' ') r(2)=kwrco3_rct(i,2)(1:6)    ! 2nd reactant

! send the products in radchk (only carbon/radicals species)
      ncopd=0 ; npd2=0 ; sccopd(:)=0. ; coprod(:)=' '
      scpd2(:)=0. ; delpd2(:)=' '
      DO j=1,SIZE(kwrco3_prd,2)      ! only 3 products
        IF (kwrco3_prd(i,j)(1:1)==' ') EXIT
        IF (INDEX(kwrco3_prd(i,j),'.')==0) CYCLE 
        IF (INDEX(kwrco3_prd(i,j),'C')==0) CYCLE 

        sctemp=kwrco3_stoi(i,j)                    ! store stoi. coef.
        CALL radchk(kwrco3_prd(i,j),rdckpd,rdckcopd,nip,sc,nref,ref)
        kwrco3_prd(i,j)=rdckpd(1)                  ! keep the 1st product
        CALL stdchm(kwrco3_prd(i,j))
        kwrco3_stoi(i,j)= sctemp*sc(1)    ! set coef. by delocalisation fraction
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
      DO j=1,SIZE(kwrco3_prd,2)
        IF (kwrco3_prd(i,j)/=' ') THEN
          CALL add1tonp(progname,chem,np)  
          s(np)=kwrco3_stoi(i,j) ; brtio=brch*s(np)
          CALL bratio(kwrco3_prd(i,j),brtio,p(np),nref,ref)
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
  IF (loknow) RETURN

!====== MAKE CHEMISTRY OF THE RCO3

! locate peroxy_acyl group and count number of carbons:
  ngr=COUNT(group/=' ')
  DO i=1,ngr
    IF (INDEX(group(i),'CO(OO.)')/=0) ip=i
  ENDDO

! reaction with NO
! -------------------
  CALL rco3no(chem,pname,bond,group,brch,ip)

! reaction with NO2
! -------------------
  CALL rco3no2(chem,pname,bond,group,ngr,brch,ip)
  IF (highnoxfg) RETURN

! reaction with NO3
! -------------------
  CALL rco3no3(chem,pname,bond,group,brch,ip)

! reaction with HO2
! -------------------
  CALL rco3ho2(chem,pname,bond,group,ngr,brch,ip)

! RCO3+RO2 reactions
!-------------------------------
  CALL rco3ro2(chem,pname,bond,group,brch,ip)
      
END SUBROUTINE rco3

END MODULE rco3chem
