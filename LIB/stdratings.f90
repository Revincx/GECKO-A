MODULE stdratings
  IMPLICIT NONE
  CONTAINS

!=======================================================================
! PURPOSE: Rates groups in molecule based on node bond and group types.  
! Small rating gives higher priority. Default rating = '999....' so is 
! LOW priority.    
!=======================================================================
SUBROUTINE ratings(nca,ogroup,bond,nring,rank)
  USE keyparameter, ONLY : digit,mxring
  USE rjtool
  USE primetool
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nca   ! # of nodes in molecule 
  INTEGER, INTENT(IN) :: nring ! # of rings in molecule
  CHARACTER(LEN=*), INTENT(IN) :: ogroup(:) ! string (group) for each node
  INTEGER, INTENT(IN) :: bond(:,:)         ! node-node bond matrix  
  INTEGER, INTENT(OUT):: rank(:)           ! integer ranks for groups

  CHARACTER(LEN=LEN(ogroup(1))) :: group(SIZE(ogroup)) ! cp of group
  INTEGER         :: i,j,cntr
  INTEGER         :: ncx(SIZE(ogroup))  ! n_connections on node
  INTEGER         :: rjg(mxring,2)      ! ring-join group pairs
  INTEGER,PARAMETER  :: nrat=24
  CHARACTER(LEN=nrat):: nines, maxm
  CHARACTER(LEN=nrat):: rating(SIZE(ogroup))

! initialize
  rank(:)=0
! ba - cp needed to work on group to handle ring (conflict with "intent(in)")      
  group(:)=ogroup(:)  

! find # of adjacent nodes, convert to character
  DO i=1,nca
    ncx(i)=0
    rating(i)(1:1) = '0'
    DO j=1,nca
      IF (bond(i,j)>0) THEN
        ncx(i)=ncx(i)+1
        rating(i)(1:1) =  digit(ncx(i))
      ENDIF
    ENDDO
  ENDDO

! remove ring-joining characters from groups if rings exist
  IF (nring > 0) CALL rjgrm(nring,group,rjg)

  nines = '999999999999999999999999'
  DO i=1,nca
    rating(i)(2:) = nines
!rating(i)(1:1) reserved for # non-H bonds
    IF(INDEX(group(i),'.' )/=0)         rating(i)(2:2) = '1'
    IF(INDEX(group(i),'CdO' )/=0)       rating(i)(4:4) = '1'
    IF(INDEX(group(i),'Cd' )/=0)        rating(i)(5:5) = '3'
    IF(INDEX(group(i),'CdH' )/=0)       rating(i)(5:5) = '2'
    IF(INDEX(group(i),'CdH2' )/=0)      rating(i)(5:5) = '1'
    IF(INDEX(group(i),'-O-' )/=0)       rating(i)(6:6) = '1'
    IF(INDEX(group(i),'CHO' )  /=0)     rating(i)(7:7) = '1'
    IF(INDEX(group(i),'CO(OONO2)' )/=0) rating(i)(8:8) = '1'
    IF(INDEX(group(i),'CO(OOH)' )/=0)   rating(i)(9:9) = '1'
    IF(INDEX(group(i),'(OOH)' )/=0)     rating(i)(10:10) = '3'
    IF(INDEX(group(i),'(OOH)(OOH)' )/=0) rating(i)(10:10) = '2'
    IF(INDEX(group(i),'(OOH)(OOH)(OOH)')/=0)rating(i)(10:10) = '1'
    IF(INDEX(group(i),'(ONO2)')/=0)     rating(i)(11:11) = '3'
    IF(INDEX(group(i),'(ONO2)(ONO2)')/=0)rating(i)(11:11) = '2'
    IF(INDEX(group(i),'(ONO2)(ONO2)(ONO2)')/=0) rating(i)(11:11) = '1'
    IF(INDEX(group(i),'CO(OH)' )/=0)    rating(i)(12:12) = '1'
    IF(INDEX(group(i),'CO' )/=0)        rating(i)(13:13) = '1'
    IF(INDEX(group(i),'F' )/=0)         rating(i)(14:14) = '1'
    IF(INDEX(group(i),'Br' )/=0)        rating(i)(15:15) = '1'
    IF(INDEX(group(i),'Cl' )/=0)        rating(i)(16:16) = '1'
    IF(INDEX(group(i),'S' )/=0)         rating(i)(17:17) = '1'
    IF(INDEX(group(i),'NH' )/=0)        rating(i)(18:18) = '1'
    IF(INDEX(group(i),'(NO2)' )/=0)     rating(i)(19:19) = '3'
    IF(INDEX(group(i),'(NO2)(NO2)' )/=0) rating(i)(19:19) = '2' 
    IF(INDEX(group(i),'(NO2)(NO2)(NO2)')/=0)rating(i)(19:19) = '1'
    IF(INDEX(group(i),'NO' )/=0)        rating(i)(20:20) = '1'
    IF(INDEX(group(i),'(OH)' )/=0)      rating(i)(21:21) = '3'
    IF(INDEX(group(i),'(OH)(OH)' )/=0)  rating(i)(21:21) = '2'
    IF(INDEX(group(i),'(OH)(OH)(OH)' )/=0) rating(i)(21:21) = '1'
    IF(INDEX(group(i),'C' )/=0)         rating(i)(22:22) = '4'
    IF(INDEX(group(i),'CH' )/=0)        rating(i)(22:22) = '3'
    IF(INDEX(group(i),'CH2' )/=0)       rating(i)(22:22) = '2'
    IF(INDEX(group(i),'CH3' )/=0)       rating(i)(22:22) = '1'
    IF(INDEX(group(i),'c' )/=0)         rating(i)(23:23) = '2'
    IF(INDEX(group(i),'cH' )/=0)        rating(i)(23:23) = '1'
  ENDDO

  cntr = 0
  nines = '999999999999999999999999'
  grloop:DO i = 1, nca
    maxm = nines
    DO j = 1, nca
      IF (rating(j)<maxm) maxm = rating(j)
    ENDDO
    IF (maxm==nines) EXIT grloop
    cntr = cntr + 1
    DO j = 1, nca
      IF (rating(j) <= maxm) THEN
         rating(j) = nines
         rank(j) = cntr
      ENDIF
    ENDDO
  ENDDO grloop

  DO i = 1, nca
    IF(rank(i)==0)THEN
      WRITE(6,*) '--error--, in ratings, rank not assigned for group:',i
      STOP "in ratings"
    ENDIF
  ENDDO

! add ring-joining characters if rings exist
  IF (nring>0) CALL rjgadd(nring,group,rjg)

! use product of adjacent primes to find unique ranks
  CALL primes(nca,bond,rank)

END SUBROUTINE ratings
  
END MODULE stdratings
