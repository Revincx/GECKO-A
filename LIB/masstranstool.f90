MODULE masstranstool 
IMPLICIT NONE
CONTAINS

!SUBROUTINE changephase(ncha,chatab)
!SUBROUTINE gas2aero(lrea,ldic,name,mweight)
!SUBROUTINE gas2wall(lrea,ldic,chem,name,mweight)
! ======================================================================
! Purpose: scroll the dictionary and write the mass transfer pseudo
! reactions for any non radical species identified. 
!
! Note: CHA (Cyclic Hemi Acetal) might be considered has non-volatile
! species. 
! ======================================================================
SUBROUTINE changephase()
  USE keyparameter, ONLY: mxlfo,mxlco,mecu,prtu,walu,tfu1,dirout
  USE dictstackdb, ONLY: ncha,chatab,nwpspe
  USE keyflag, ONLY: g2pfg,g2wfg,chafg
  USE dictstackdb, ONLY: nrec,dict
  USE atomtool, ONLY: cnum, molweight  
  IMPLICIT NONE
    
  INTEGER :: i,j
  CHARACTER(LEN=mxlfo) :: chem
  CHARACTER(LEN=mxlco) :: idnam
  REAL :: molmass
  
  IF (chafg) OPEN(tfu1, FILE=TRIM(dirout)//'chaname.dat')
  
! scroll the dictionary
  recloop: DO i=2,nrec
    chem=dict(i)(10:129)
    IF (INDEX(chem,'.')/=0) CYCLE recloop           ! rm radicals
    IF (cnum(chem) < 2)     CYCLE recloop           ! rm C1 and inorg. 

! check if '#' can be managed
    IF (chem(1:1)/='C'.AND.chem(1:1)/='c'.AND.chem(1:2)/='-O' ) THEN
      IF (chem(1:3)=='#mm' ) THEN
        chem(1:)=chem(4:)
      ELSE IF (chem(1:1)=='#' ) THEN
        IF (INDEX(chem,'CH2OO')/=0)  CYCLE recloop  ! rm crieggee
        IF (INDEX(chem,'.(OO.)')/=0) CYCLE recloop  ! rm crieggee
        chem(1:)=chem(2:)
      ELSE
        CYCLE recloop                               ! rm unexpected species 
      ENDIF
    ENDIF
    CALL molweight(chem,molmass)
    idnam=dict(i)(1:mxlco)

! manage CHA species (assumed non volatile if chafg true: write dict only) 
    IF (chafg) THEN  
      DO j=1,ncha
        IF (idnam==chatab(j)) THEN 
          WRITE(tfu1,'(A1,A6)') "G",chatab(j)  ! saved here for info only
          IF (g2pfg) THEN
            nwpspe=nwpspe+1 
            WRITE(prtu,'(A1,A6,10X,A1,F6.1,A1)') "A",idnam,"/",molmass ,"/"
          ENDIF
          IF (g2wfg) THEN
            nwpspe=nwpspe+1 
            WRITE(walu,'(A1,A6,10X,A1,F6.1,A1)') "W",idnam,"/",molmass ,"/"
          ENDIF
          CYCLE recloop   ! no mass tranfer !
        ENDIF 
      ENDDO
    ENDIF

! write mass transfer equations 
    IF (g2pfg) CALL gas2aero(mecu,prtu,idnam,molmass)
    IF (g2wfg) CALL gas2wall(mecu,walu,chem,idnam,molmass)
    
  ENDDO recloop
  IF (chafg) CLOSE(tfu1)

END SUBROUTINE changephase

! ======================================================================
! Purpose: Add the particle phase species in the list of species and 
! write the 2 pseudo reactions for the gas <-> particle mass transfer 
! for the species provided as input (idnam). 
! ======================================================================
SUBROUTINE gas2aero(lrea,ldic,idnam,mweight)
  USE rxwrttool, ONLY: rxinit, rxwrit_dyn
  USE dictstackdb, ONLY: nwpspe

  CHARACTER(LEN=*),INTENT(IN) :: idnam ! name of species undergoing phase transfer
  INTEGER,INTENT(IN)   :: lrea         ! unit file to write the pseudo reaction
  INTEGER,INTENT(IN)   :: ldic         ! unit file to write new condensed phase species
  REAL,INTENT(IN)      :: mweight      ! molecular weigth of the species

  INTEGER,PARAMETER :: mxprod=4  ! (only 1 product used here) 
  INTEGER,PARAMETER :: mxreac=3  ! (only 2 reactants here (chem + keywd) 
  CHARACTER(LEN=LEN(idnam)) :: r(mxreac), p(mxprod)
  REAL            :: s(mxprod),arrh(3)
  INTEGER         :: idreac
  REAL            :: auxinfo(9)  
  CHARACTER*1     :: charfrom,charto

! initialize 
  r(:)=' '   ; s(:)=0.  ;  p(:)=' ' ;  arrh(:)=0. ; idreac=0
  auxinfo(:)=0.  

! write the aerosol species in the dictionary
  nwpspe=nwpspe+1
  WRITE(ldic,'(A1,A6,10X,A1,F6.1,A1)')  'A',idnam,'/',mweight ,'/'

! ---- from gas to aerosol
  r(1)= idnam ;   r(2)= 'AIN ' ;  p(1)= idnam  ;  s(1) = 1.
  arrh(1) = 1.0
  idreac=1
  charfrom='G'  ;  charto='A'
  CALL rxwrit_dyn(lrea,r,s,p,arrh,idreac,auxinfo,charfrom,charto)

! ---- from aerosol to gas
  r(1)= idnam  ;   r(2)= 'AOU '  ;  p(1)= idnam  ;  s(1) = 1.
  arrh(1) = 1.0
  idreac=2
  charfrom='A'  ;   charto='G'
  CALL rxwrit_dyn(lrea,r,s,p,arrh,idreac,auxinfo,charfrom,charto)

END SUBROUTINE gas2aero

! ======================================================================
! Purpose: Add the wall phase species in the list of species and
! write the 2 pseudo reactions for the gas <-> wall mass transfer 
! for the species provided as input (idnam). 
! ======================================================================
SUBROUTINE gas2wall(lrea,ldic,chem,idnam,mweight)
  USE dictstackdb, ONLY: nwpspe
  USE rxwrttool, ONLY: rxinit, rxwrit_dyn

  CHARACTER(LEN=*),INTENT(IN) :: chem   ! formula of species undergoing phase transfer
  CHARACTER(LEN=*),INTENT(IN) :: idnam   ! name of species undergoing phase transfer
  INTEGER,INTENT(IN) :: lrea     ! unit file to write the pseudo reaction
  INTEGER,INTENT(IN) :: ldic     ! unit file to write new condensed phase species
  REAL,INTENT(IN)    :: mweight  ! molecular weigth of the species

  INTEGER,PARAMETER :: mxprod=4  ! (only 1 product used here) 
  INTEGER,PARAMETER :: mxreac=3  ! (only 2 reactants here (chem + keywd) 
  CHARACTER(LEN=LEN(idnam)) :: r(mxreac), p(mxprod)
  REAL            :: s(mxprod),arrh(3)
  INTEGER         :: idreac
  REAL            :: auxinfo(9)  
  CHARACTER*1     :: charfrom,charto

! initialize
  r(:)=' '   ; s(:)=0.  ;  p(:)=' ' ;  arrh(:)=0. ; idreac=0
  auxinfo(:)=0.  

! write the aerosol species in the dictionary
  nwpspe=nwpspe+1
  WRITE(ldic,'(A1,A6,10X,A1,F6.1,A1)') 'W',idnam,'/',mweight ,'/'

! ---- from gas to wall
  r(1)=idnam  ;  r(2)='WIN '  ;  p(1)=idnam  ;  s(1)=1.
  arrh(1) = 1.111E-3 ! measures 2-ketones Matsunaga et al.2010 AST
  idreac=3
  charfrom='G'  ;  charto='W'
  IF (INDEX(chem,'O')==0) THEN
    auxinfo(1)=2.E-5  ! measures 1-alkenes Matsunaga et al.2010 AST
  ELSE
    auxinfo(1)=1.2E-4 ! measures 2-ketones Matsunaga et al.2010 AST
  ENDIF 
  CALL rxwrit_dyn(lrea,r,s,p,arrh,idreac,auxinfo,charfrom,charto)

! ---- from wall to gas
  r(1)=idnam  ;   r(2)='WOU '  ;  p(1)=idnam  ;  s(1) = 1.
  arrh(1)=1.111E-3 
  idreac=4
  charfrom='W'  ;  charto='G'
  IF (INDEX(chem,'O')==0) THEN
    auxinfo(1)=2.E-5  ! measures 1-alkenes Matsunaga et al.2010 AST
  ELSE
    auxinfo(1)=1.2E-4 ! measures 2-ketones Matsunaga et al.2010 AST
  ENDIF 
  CALL rxwrit_dyn(lrea,r,s,p,arrh,idreac,auxinfo,charfrom,charto)

END SUBROUTINE gas2wall

END MODULE masstranstool 
