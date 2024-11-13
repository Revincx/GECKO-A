MODULE rooohchem
IMPLICIT NONE
CONTAINS
! SUBROUTINE rooohdec(idnam,chem,bond,group,brch)
!=======================================================================
! PURPOSE: 
! R(OOOH)-> R(O.)+HO2     
! R(OOOH)-> R(OH)+ O2     
!=======================================================================
SUBROUTINE rooohdec(idnam,chem,bond,group,brch)
  USE keyparameter, ONLY: mxpd,mxnr,mecu
  USE keyflag, ONLY: screenfg
  USE references, ONLY:mxlcod
  USE atomtool, ONLY: cnum
  USE reactool, ONLY: swap, rebond
  USE normchem, ONLY: stdchm
  USE toolbox, ONLY: stoperr,add1tonp,addrx,addref
  USE dictstacktool, ONLY: bratio
  USE rxwrttool, ONLY:rxwrit,rxinit
  USE radchktool, ONLY: radchk
  USE tempo, ONLY: chknadd2p
  IMPLICIT NONE

  CHARACTER(LEN=*),INTENT(IN) :: idnam    ! name of current species
  CHARACTER(LEN=*),INTENT(IN) :: chem     ! formula of current species
  CHARACTER(LEN=*),INTENT(IN) :: group(:) ! group (without ring joining char)
  INTEGER,INTENT(IN) :: bond(:,:)         ! bond matrix
  REAL,INTENT(IN)    :: brch              ! max yield of the input species

  CHARACTER(LEN=LEN(group(1))) :: tgroup(SIZE(group))
  CHARACTER(LEN=LEN(group(1))) :: pold,pnew
  INTEGER         :: tbond(SIZE(bond,1),SIZE(bond,2))
  CHARACTER(LEN=LEN(chem)) :: tempfo

  CHARACTER(LEN=LEN(idnam)) :: p(mxpd), r(3)
  REAL    :: s(mxpd), arrh(3)
  INTEGER :: idreac, nlabel
  REAL    :: xlabel,folow(3),fotroe(4)

  INTEGER :: i,np,nca,ngr,tempnring
  REAL    :: brtio,yield

  INTEGER,PARAMETER    :: mxcopd=4     !--- max # of copdct per reaction
  REAL                 :: arrhc(3)

  INTEGER,PARAMETER  :: mxrpd=2  ! max # of products returned by radchk sub
  CHARACTER(LEN=LEN(chem)) :: rdckpd(mxrpd)
  CHARACTER(LEN=LEN(idnam)) :: rdckcopd(mxrpd,mxcopd)
  REAL    :: sc(mxrpd)
  INTEGER :: nip

  CHARACTER(LEN=8),PARAMETER :: progname='ROOOHdec'
  INTEGER,PARAMETER     :: mxref=10
  INTEGER               :: nref
  CHARACTER(LEN=mxlcod) :: ref(mxref)

  IF (screenfg) WRITE(6,*)progname

  tgroup(:)=group(:) ; tbond(:,:)=bond(:,:) ; arrhc(:)=0. 
  
  ngr=COUNT(tgroup/=' ')
  nca=cnum(chem)   ! count number of carbon atoms

! Rates constant from Assaf et al., 2018
    arrhc(1) = 1.9E10 ; arrhc(2) = 1.35 ; arrhc(3) = 12000
    yield=1.

! FIND ALL POSSIBLE CHANNEL
! -------------------------
  rooohloop: DO i=1,ngr
    IF (INDEX(group(i),'(OOOH)')==0) CYCLE rooohloop

! == First pathway : ROOOH -> RO + HO2
    CALL rxinit(r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe)
    nref=0  ;  ref(:)=' ' 
    CALL addref(progname,'MJ19KMV000',nref,ref,chem)
    CALL addref(progname,'EA18KMR000',nref,ref,chem)
    CALL addref(progname,'ROOOH_DEC ',nref,ref,chem)
    np=0

    pold='(OOOH)' ; pnew='(O.)'
    CALL swap(group(i),pold,tgroup(i),pnew)
    CALL rebond(tbond,tgroup,tempfo,tempnring)
    CALL radchk(tempfo,rdckpd,rdckcopd,nip,sc,nref,ref)
    tempfo=rdckpd(1)
    brtio=brch
    IF (INDEX(chem,'.')/=0) THEN
      CALL chknadd2p(chem,tempfo,yield,brch,np,s,p,nref,ref)
    ELSE
      CALL stdchm(tempfo)
      CALL add1tonp(progname,chem,np)   
      s(np)=1.0 
      brtio=s(np)*brch
      CALL bratio(tempfo,brtio,p(np),nref,ref)
    ENDIF

    CALL add1tonp(progname,chem,np)   
    s(np)=1. ;  p(np)='HO2 '

    r(1)=idnam
    arrh(1)=arrhc(1)*0.8 ; arrh(2:3)=arrhc(2:3)
    CALL rxwrit(mecu,r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe,ref)
! reset:
    tgroup(:)=group(:)

! == Second pathway : ROOOH -> ROH + O2
    CALL rxinit(r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe)
    nref=0  ;  ref(:)=' '
    CALL addref(progname,'MJ19KMV000',nref,ref,chem)
    CALL addref(progname,'EA18KMR000',nref,ref,chem)
    CALL addref(progname,'ROOOH_DEC ',nref,ref,chem)
    np=0

    pold='(OOOH)' ; pnew='(OH)'
    CALL swap(group(i),pold,tgroup(i),pnew)
     
    CALL rebond(tbond,tgroup,tempfo,tempnring)
    CALL stdchm(tempfo)
    CALL add1tonp(progname,chem,np)   
    s(np)=1.0 
    brtio=s(np)*brch
    CALL bratio(tempfo,brtio,p(np),nref,ref)

    r(1)=idnam
    arrh(1)=arrhc(1)*0.2 ; arrh(2:3)=arrhc(2:3)
    CALL rxwrit(mecu,r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe,ref)
! reset:
    tgroup(:)=group(:)

  ENDDO rooohloop

END SUBROUTINE rooohdec

END MODULE rooohchem
