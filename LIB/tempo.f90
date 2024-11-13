MODULE tempo
IMPLICIT NONE
CONTAINS
! ======================================================================
! PURPOSE: the routine chknadd2p checks the radical (radchem) provided  
! as input, send the species to bratio (addition to the stack) and load  
! the short name species to the p(:) list, with the corresponding  
! stoichiometric coefficient in s(:). The coproducts of the species are 
! also added to the p(:) as well as the species that may be produced   
! after electron delocalisation.
! ======================================================================
SUBROUTINE chknadd2p(chem,radchem,yield,brch,np,s,p,nref,com)
  USE keyparameter, ONLY: mxcopd
  USE dictstacktool, ONLY: bratio
  USE radchktool, ONLY: radchk
  USE normchem, ONLY: stdchm      
  USE toolbox, ONLY: stoperr,add1tonp
  IMPLICIT NONE

  CHARACTER(LEN=*),INTENT(IN) :: chem ! parent species that produced chem
  CHARACTER(LEN=*),INTENT(IN) :: radchem  ! radical to be managed
  REAL,INTENT(IN)      :: brch        ! branching ratio of the parent species
  REAL,INTENT(IN)      :: yield       ! yield of the radical
  INTEGER,INTENT(INOUT):: np          ! # of product (incremented here)
  REAL,INTENT(INOUT)   :: s(:)        ! stoi. coef. of the products in p
  CHARACTER(LEN=*),INTENT(INOUT):: p(:)   ! list of the products (short name)
  INTEGER,INTENT(INOUT):: nref            ! # of references
  CHARACTER(LEN=*),INTENT(INOUT):: com(:) ! reac_info code 
  
  CHARACTER(LEN=LEN(chem)) :: tempfo
  INTEGER :: j
  REAL :: brtio

  INTEGER,PARAMETER :: mxrpd=2  ! max # of products returned by radchk sub
  CHARACTER(LEN=LEN(chem))  :: rdckpd(mxrpd)
  CHARACTER(LEN=LEN(p(1))) :: rdckcopd(mxrpd,mxcopd)
  REAL :: sc(mxrpd)
  INTEGER :: nip
  CHARACTER(LEN=20),PARAMETER :: progname='chknadd2p'

! send to radchk to get product and coproduct
  CALL radchk(radchem,rdckpd,rdckcopd,nip,sc,nref,com)

! 1st products from radchk 
  tempfo=rdckpd(1)
  CALL stdchm(tempfo)
  CALL add1tonp(progname,chem,np)
  s(np)=yield*sc(1)
  brtio=brch*s(np)
  CALL bratio(tempfo,brtio,p(np),nref,com)
  DO j=1,SIZE(rdckcopd,2)
    IF (rdckcopd(1,j)(1:1)==' ') CYCLE
    CALL add1tonp(progname,chem,np)
    s(np)=yield*sc(1) ; p(np)=rdckcopd(1,j)
   ENDDO

! 2nd products from radchk 
   IF (nip==2) THEN
     tempfo=rdckpd(2)
     CALL stdchm(tempfo)
     CALL add1tonp(progname,chem,np)
     s(np)=yield*sc(2)
     brtio=brch*s(np)
     CALL bratio(tempfo,brtio,p(np),nref,com)
     DO j=1,SIZE(rdckcopd,2)
       IF (rdckcopd(2,j)(1:1)==' ') CYCLE
       CALL add1tonp(progname,chem,np)
       s(np)=yield*sc(2) ; p(np)=rdckcopd(2,j)
     ENDDO
   ENDIF
END SUBROUTINE chknadd2p

END MODULE tempo
