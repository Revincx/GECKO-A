MODULE atomtool
  IMPLICIT NONE
  CONTAINS
  
!INTEGER FUNCTION cnum(chem)
!INTEGER FUNCTION onum(chem)
!SUBROUTINE molweight(chem,weight)
!SUBROUTINE getatoms(chem,ica,ih,init,io,ir,isul,ifl,ibr,icl)

! ======================================================================
! Return the # of carbon atoms in a molecule                    
! ======================================================================
INTEGER FUNCTION cnum(chem)
  IMPLICIT NONE

  CHARACTER(LEN=*), INTENT(IN) :: chem   !  chemical formula
  INTEGER :: ic, i, nc

  cnum = 0 ; ic = 0 
  nc=INDEX(chem,' ')
  DO i=1,nc
    IF  (chem(i:i)=='C' .AND. chem(i:i+1)/='Cl') ic=ic+1
    IF  (chem(i:i)=='c') ic=ic+1
  ENDDO
  cnum = ic
  
END FUNCTION cnum

! ======================================================================
! Return the # of C-O-C (ether) in a molecule.
! (e.g. to estimate the # of of nodes in chem  
! ======================================================================
INTEGER FUNCTION onum(chem)
  IMPLICIT NONE

  CHARACTER(LEN=*),INTENT(IN) :: chem  ! chemical formula
  INTEGER :: io, i,nc

  onum = 0  ;  io   = 0
  nc=INDEX(chem,' ')-1
  DO i=1,nc
    IF (chem(i:i+1)=='-O') io=io+1
  ENDDO
  onum = io

END FUNCTION onum

! ======================================================================
! PURPOSE: compute the molecular weight of the species provided as input
! ======================================================================
SUBROUTINE molweight(chem,weight)
  IMPLICIT NONE

  CHARACTER(LEN=*),INTENT(IN) :: chem
  REAL,INTENT(OUT) :: weight

  INTEGER  ica, iha, ina, ioa, ira, isa, ifl, ibr, icl

  CALL getatoms(chem,ica,iha,ina,ioa,ira,isa,ifl,ibr,icl)
 
  weight = ica*12.01 + iha*1.01 + ina*14.01 + ioa*16 + isa*32.07 + &
           ifl*19 + ibr*79.9 + icl*35.45

END SUBROUTINE molweight

! ======================================================================
! Find the number of each possible atom in the species CHEM,
! e.g. C's, H's, N's, O's, S's, F's, Cl's, Br's, and .'s  
! ======================================================================
SUBROUTINE getatoms(chem,ica,ih,init,io,ir,isul,ifl,ibr,icl)
  USE keyparameter, ONLY:mxring
  USE rjtool
  IMPLICIT NONE
	 
  CHARACTER(LEN=*),INTENT(IN) ::  chem ! formula of the species 
  INTEGER,INTENT(OUT) :: ica  ! # of carbon atoms
  INTEGER,INTENT(OUT) :: ih   ! # of hydrogen atoms
  INTEGER,INTENT(OUT) :: init ! # of nitrogen atoms
  INTEGER,INTENT(OUT) :: io   ! # of oxygen atoms
  INTEGER,INTENT(OUT) :: ir   ! # of radical dot
  INTEGER,INTENT(OUT) :: isul ! # of sulfur atoms
  INTEGER,INTENT(OUT) :: ifl  ! # of fluorine atoms
  INTEGER,INTENT(OUT) :: ibr  ! # of bromine atoms
  INTEGER,INTENT(OUT) :: icl  ! # of chlorine atoms

  INTEGER         nring, i, n, nc
  CHARACTER(LEN=LEN(chem)) :: tchem  ! a copy of chem
  INTEGER         rjs(mxring,2)      ! ring-join group pairs
  LOGICAL         loring

! initialize
! ----------
  tchem=chem  ; nc=INDEX(tchem,' ')
  ifl = 0  ;  ibr = 0  ;  icl = 0 
  ica = 0  ;  ih  = 0  ;  init= 0  ;  io  = 0
  isul= 0  ;  ir  = 0
  n  = 1

  loring=.FALSE.
  IF (INDEX(tchem,'C1')/=0)  loring=.TRUE.
  IF (INDEX(tchem,'Cd1')/=0)  loring=.TRUE.
  IF (INDEX(tchem,'c1')/=0)  loring=.TRUE.
  IF (INDEX(tchem,'-O1')/=0) loring=.TRUE.
  IF (loring) THEN
     nring=2  ! set max # of rings (since # of rings is unknown here)
    CALL rjsrm(nring,tchem,rjs)
  ENDIF
  
! start counting
! --------------
  i=nc
  loopchar: DO
    i=i-1
    IF (i<1) EXIT loopchar
    IF (tchem(i:i)=='2') n = 2
    IF (tchem(i:i)=='3') n = 3
    IF (tchem(i:i)=='4') n = 4
    IF (i>1) THEN
      IF (tchem(i-1:i)=='Cl') THEN
        icl = icl + n
        n = 1  
        i = i-1
      ENDIF
    ENDIF
    IF (tchem(i:i)=='F') THEN      ; ifl  = ifl+n  ; n=1
    ELSE IF (tchem(i:i+1)=='Br') THEN ; ibr  = ibr+n  ; n=1
    ELSE IF (tchem(i:i)=='H') THEN ; ih   = ih+n   ; n=1
    ELSE IF (tchem(i:i)=='N') THEN ; init = init+n ; n=1
    ELSE IF (tchem(i:i)=='O') THEN ; io   = io+n   ; n=1
    ELSE IF (tchem(i:i)=='S') THEN ; isul = isul+1 ; n=1
    ELSE IF (tchem(i:i)=='.') THEN ; ir = ir+n     ; n=1
    ELSE IF ((tchem(i:i+1)/='Cl') .AND. &
             ((tchem(i:i)=='C').OR.(tchem(i:i)=='c')) ) THEN
       ica = ica+1 ; n  = 1
    ENDIF
  ENDDO loopchar

END SUBROUTINE getatoms

SUBROUTINE oxnum(chem,numoxC)
  USE keyparameter, ONLY:mxring,mxnode,mxlgr
  USE rjtool
  USE stdgrbond
  USE toolbox     , ONLY:countstring
  IMPLICIT NONE

  CHARACTER(LEN=*),INTENT(IN) :: chem   ! formula of the species
  REAL,INTENT(OUT)            :: numoxC ! Mean number of oxidation of carbon

  INTEGER                     :: noC(mxnode),nC,ngr
  INTEGER                     :: bond(mxnode,mxnode),dbflg,nring
  CHARACTER(LEN=mxlgr)        :: group(mxnode)
  INTEGER                     :: rjg(mxring,2) ! ring-join group pairs
  INTEGER                     :: i,j,nfunc

  numoxC=0. ; noC(:)=0 ; nC=0

  CALL grbond(chem,group,bond,dbflg,nring)
  ngr=COUNT(group/=' ')
  CALL rjgrm(nring,group,rjg)

  DO i=1,ngr
    IF ((group(i)(1:1)/='C').AND.(group(i)(1:1)/='c')) CYCLE
    nC=nC+1
    IF      (group(i)(1:3)=='CH3')  THEN ; noC(i)= -3
    ELSE IF (group(i)(1:3)=='CH2')  THEN ; noC(i)= -2
    ELSE IF (group(i)(1:2)=='CH')   THEN ; noC(i)= -1
    ELSE IF (group(i)(1:2)=='cH')   THEN ; noC(i)= -1
    ELSE IF (group(i)(1:3)=='CHO')  THEN ; noC(i)=  1
    ELSE IF (group(i)(1:2)=='CO')   THEN ; noC(i)=  2
    ENDIF

    nfunc=countstring(group(i),'(O') ; IF (nfunc/=0) noC(i)=noC(i)+1
    nfunc=countstring(group(i),'(N') ; IF (nfunc/=0) noC(i)=noC(i)+1

    DO j=1,ngr ; IF (bond(i,j)==3) noC(i)=noC(i)+1 ; ENDDO
  ENDDO

  numoxC=(SUM(noC(:)))/REAL(nC)

END SUBROUTINE oxnum

END MODULE atomtool
