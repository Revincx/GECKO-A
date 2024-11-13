MODULE panchem
IMPLICIT NONE
CONTAINS
! SUBROUTINE pandec(idnam,chem,bond,group,brch)
!=======================================================================
! PURPOSE: In pandec, each PAN group is treated independently and all 
! possible break-up of the PAN moiety are considered, i.e. 
! RCO(OONO2)-> RCO(OO.)+NO2     
!=======================================================================
SUBROUTINE pandec(idnam,chem,bond,group,brch)
  USE keyparameter, ONLY: mxpd,mxnr,mecu,mxcopd
  USE keyflag, ONLY: wrtopeinfo,screenfg,rx_ro2_no2
  USE references, ONLY:mxlcod
  USE atomtool, ONLY: cnum
  USE reactool, ONLY: swap, rebond
  USE ro2tool, ONLY: ro2no2_fo_reverse
  USE normchem, ONLY: stdchm
  USE toolbox, ONLY: stoperr,add1tonp,addrx,addref
  USE dictstacktool, ONLY: bratio
  USE rxwrttool, ONLY:rxwrit,rxinit
  USE radchktool, ONLY: radchk
  IMPLICIT NONE

  CHARACTER(LEN=*),INTENT(IN) :: idnam    ! name of current species
  CHARACTER(LEN=*),INTENT(IN) :: chem     ! formula of current species
  CHARACTER(LEN=*),INTENT(IN) :: group(:) ! group (without ring joining char)
  INTEGER,INTENT(IN) :: bond(:,:)         ! bond matrix
  REAL,INTENT(IN)    :: brch              ! max yield of the input species

  CHARACTER(LEN=LEN(group(1))) :: tgroup(SIZE(group))
  CHARACTER(LEN=LEN(group(1))) :: pold,pnew
  INTEGER         :: tbond(SIZE(bond,1),SIZE(bond,2))
  CHARACTER(LEN=LEN(chem)) :: tempkc

  CHARACTER(LEN=LEN(idnam)) :: p(mxpd), r(3)
  REAL    :: s(mxpd), arrh(3)
  INTEGER :: idreac, nlabel
  REAL    :: xlabel,folow(3),fotroe(4)

  INTEGER :: nr,i,j,np,nca,ngr,tempnring
  REAL    :: brtio,aa
  CHARACTER(LEN=LEN(idnam)) :: pname(mxnr)

  INTEGER              :: flag(mxnr)
  CHARACTER(LEN=LEN(chem)) :: pchem(mxnr)
  CHARACTER(LEN=LEN(idnam)) :: coprod(mxnr,mxcopd)
  REAL                 :: tarrhc(mxnr,3), lotarrhc(mxnr,3), fc(mxnr)

  INTEGER,PARAMETER  :: mxrpd=2  ! max # of products returned by radchk sub
  CHARACTER(LEN=LEN(chem)) :: rdckpd(mxrpd)
  CHARACTER(LEN=LEN(idnam)) :: rdckcopd(mxrpd,mxcopd)
  REAL    :: sc(mxrpd)
  INTEGER :: nip

  CHARACTER(LEN=7),PARAMETER :: progname='PANdec'
  CHARACTER(LEN=70)          :: mesg
  INTEGER,PARAMETER     :: mxcom=10
  INTEGER               :: nrxref(mxnr)
  CHARACTER(LEN=mxlcod) :: rxref(mxnr,mxcom)

  IF (screenfg) WRITE(6,*)progname

  nrxref(:)=0 ; rxref(:,:)=' '
  tgroup(:)=group(:) ; tbond(:,:)=bond(:,:)
  tarrhc(:,:)=0. ; lotarrhc(:,:)=0. ; fc(:)=0.
  flag(:)=0  ; pchem(:)=' ' ; coprod(:,:)=' '
  pname(:)=' '

  ngr=COUNT(tgroup/=' ')
  nca=cnum(chem)   ! count number of carbon atoms

! FIND ALL POSSIBLE CHANNELS
! --------------------------
  nr=0
  panloop: DO i=1,ngr
    IF (INDEX(group(i),'(OONO2)')==0) CYCLE panloop
    CALL addrx(progname,chem,nr,flag)

! change (OONO2) to (OO.)
    pold='(OONO2)' ; pnew='(OO.)'
    CALL swap(group(i),pold,tgroup(i),pnew)

! rebuild, check and rename:
    CALL rebond(tbond,tgroup,tempkc,tempnring)
    CALL radchk(tempkc,rdckpd,rdckcopd,nip,sc,nrxref(nr),rxref(nr,:))
    pchem(nr)=rdckpd(1)
    CALL stdchm(pchem(nr))
    coprod(nr,:)=rdckcopd(1,:)
    IF (nip/=1) THEN
      mesg="unexpected 2nd product "
      CALL stoperr(progname,mesg,chem)
    ENDIF

! add NO2 as coproduct
    jloop: DO j=1,SIZE(coprod,2)
      IF (coprod(nr,j)==' ') THEN
        coprod(nr,j)='NO2 ' ; EXIT jloop
      ENDIF
      IF (j==SIZE(coprod,2)) THEN
        mesg="no slot find to add NO2"
        CALL stoperr(progname,mesg,chem)
      ENDIF
    ENDDO jloop
     
! PAN like species
! ----------------
    IF (INDEX(group(i),'CO(OONO2)')/=0) THEN

! Rate constant for PAN is used for every species having 2 carbons. This
! reaction is in the fall-off regime. Rate constant are from IUPAC (last access 2019)
! and is consistent with the values taken for the
! CH3CO(OO.)+NO and +NO2 reactions. For species with C>2, the value is 
! taken in Jenkin et al., ACP, 2019. 
! For the fall off reaction, k is given with the expression: 
! k=arrh1*(T/300)**arrh2*exp(-arrh3/T)
      IF (nca==2) THEN
        lotarrhc(nr,1)=1.1E-5 ; lotarrhc(nr,3)=10100.  ! low pressure rate constant
        tarrhc(nr,1)=1.9E17   ; tarrhc(nr,3)=14100.    ! high pressure rate constant
        fc(nr)= 0.3                                    ! broadening factor            
        CALL addref(progname,'IUPAKMV000',nrxref(nr),rxref(nr,:),chem)
      ELSE 
        tarrhc(nr,1)=5.2E16   ; tarrhc(nr,3)=13850.
        DO j=1,nca
          IF (bond(i,j)==3) THEN
            tarrhc(nr,1)=2*tarrhc(nr,1)
            EXIT
          ENDIF 
        ENDDO
        CALL addref(progname,'MJ19KMV000',nrxref(nr),rxref(nr,:),chem)
      ENDIF

! peroxy nitrate (OONO2) not PAN like forbidden
! UNLESS we are allowing RO2+NO2 reversible reactions (JPL 2006)
! Added by JMLT & Orlando, NCAR, Sept 2022
! ---------------------------------------------
    ELSE
      IF(rx_ro2_no2)THEN
! decompose RO2NO2 into RO2 + NO2 (activated if rx_ro2_no2 = .TRUE.)
! and write the reaction. No other chemistry is allowed.
       CALL ro2no2_fo_reverse(chem,idnam,bond,group,brch,i)
       RETURN
      ELSE
        mesg="unexpected (OONO2) "
        CALL stoperr(progname,mesg,chem)
      ENDIF
    ENDIF

! reset:
    tgroup(:)=group(:)
  ENDDO panloop

! ALL REACTIONS FOUND - TREAT REACTIONS
! ---------------------------------------

! collapse identical products:
! NOTE: if the reaction rates are in the fall off regime (nca<3), reaction
! rates can not easily be added - reaction will be written 2 times.
  IF (nr>1) THEN
    IF (nca>2) THEN
      DO i=1,nr-1
        DO j=i+1,nr
          IF (pchem(i)==pchem(j)) THEN
            flag(j)=0 ; tarrhc(i,1)=tarrhc(i,1)+tarrhc(j,1)
          ENDIF
        ENDDO
      ENDDO
    ENDIF
  ENDIF      

! WRITE REACTIONS
! ---------------
              
! loop over each reaction pathway
  DO i=1,nr
    IF (flag(i)==0) CYCLE
    CALL rxinit(r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe)
    r(1)=idnam

    s(1)=1. ; brtio=brch*s(1)
    CALL bratio(pchem(i),brtio,p(1),nrxref(i),rxref(i,:))
    pname(i)=p(1)          ! store for later use in operator information
       
    np=1
    DO j=1,mxcopd
      IF (coprod(i,j)(1:1)/=' ') THEN
        CALL add1tonp(progname,chem,np)
        s(np)=1.0 ; p(np)=coprod(i,j)
      ENDIF
    ENDDO

! Case 1: add keyword (+M) and write fall off reaction (idreac=3)
    IF (nca==2) THEN
      r(2)='(+M)'
      idreac=3
      arrh(:)=tarrhc(i,:) ; folow(:)=lotarrhc(i,:) ; fotroe(1)=fc(i)
      CALL rxwrit(mecu,r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe,com=rxref(i,:))
     
! Case 2: regular thermal reaction (C>2)
    ELSE
      arrh(:)=tarrhc(i,:)
      CALL rxwrit(mecu,r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe,com=rxref(i,:))
    ENDIF
  ENDDO

! -----------------------------------------
! WRITE INFORMATION REQUIRED FOR OPERATOR
! -----------------------------------------
! following block needed for "postprocessing" only (operator reduction)
  IF (wrtopeinfo) THEN  
    IF (nca/=2) THEN
      aa=0.
      DO i=1,nr
       IF (flag(i)==1) aa=aa+tarrhc(i,1)
      ENDDO
      WRITE(10,'(A10,A1,A6,5X,A4,1X,ES10.3,1X,f4.1,1X,f7.0)') &
        '**** INIT ','G',idnam,'****',aa,tarrhc(1,2),tarrhc(1,3)
    ELSE
      WRITE(10,*)'**** INIT ','G',idnam,' (+M) ','****'
    ENDIF

    DO i=1,nr
      IF (flag(i)==0) CYCLE
      IF (nca==2) THEN   
        WRITE(10,'(A6,E8.2,F5.1,F7.0,A1)') &
              'LOW  /',lotarrhc(i,1),lotarrhc(i,2),lotarrhc(i,3),'/'
        WRITE(10,'(A6,F3.1,A10)')'TROE /',fc(i), ' 0. 0. 0./'
        WRITE(10,'(E10.4,F5.1,F7.0)') tarrhc(i,1),tarrhc(i,2),tarrhc(i,3)
      ENDIF
      s(1)=1.
      WRITE(10,'(f5.3,2X,A1,A6)') s(1), 'G',pname(i)

      np=1
      DO j=1,mxcopd
        IF (coprod(i,j)(1:1)/=' ') THEN
          np=np+1 ; s(np)=1.0 ; p(np)=coprod(i,j)
        ENDIF
      ENDDO
      np=np+1 ; s(np)=1. ; p(np)='NO2  '
      WRITE(10,'(f5.3,2X,A1,A6)') s(np)*s(1),'G',p(np)
    ENDDO
    WRITE(10,*)'end'
  ENDIF 

END SUBROUTINE pandec

END MODULE panchem
