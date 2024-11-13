MODULE rjtool
  IMPLICIT NONE
  CONTAINS

!---------------------------------------------------------------------
! Add ring-joining numerical characters to group strings     
!          at nodes given by rjg.                                     
!---------------------------------------------------------------------
SUBROUTINE rjgadd(nring,group,rjg)
  USE keyparameter, ONLY : digit
  IMPLICIT NONE

  INTEGER,INTENT(IN) :: nring     ! # of rings in chem 
  INTEGER,INTENT(IN) :: rjg(:,:)  ! group # of ring-joining pairs for each ring
  CHARACTER(LEN=*),INTENT(INOUT) :: group(:)  ! groups at position (node) i

  INTEGER :: n,i,ii,j

! loop in reverse order => chars in numerical order if >1 exist at any node
  DO n=nring,1,-1
    DO ii=1,2
      i = rjg(n,ii)
      IF(i/=0)THEN
        IF (group(i)(1:2)=='-O' .OR. group(i)(2:2)=='d') THEN
          j=3
        ELSE
          j=2
        ENDIF
        group(i)(j+1:)=group(i)(j:)
        group(i)(j:j) = digit(n)
      ENDIF
    ENDDO
  ENDDO
  
  END SUBROUTINE rjgadd

!-------------------------------------------------------------
! Strip ring-joining numerical characters from group strings. 
! Subroutine rjgadd puts them back.                           
!-------------------------------------------------------------
SUBROUTINE rjgrm(nring,group,rjg)
  USE keyparameter, ONLY : digit
  IMPLICIT NONE

  INTEGER,INTENT(IN)  :: nring                ! # of rings in chem
  CHARACTER(LEN=*),INTENT(INOUT) :: group(:)  ! groups at position (node) i
  INTEGER,INTENT(OUT) ::   rjg(:,:)  ! group # of ring-joining pairs for each ring

  INTEGER :: n,i,j

  rjg(:,:)=0
  DO n=1,nring
    DO i=1,SIZE(group)
      IF(group(i)(1:2)=='-O' .OR. group(i)(2:2)=='d') THEN
        j=3
      ELSE
        j=2
      ENDIF
      IF(group(i)(j:j)==digit(n)) THEN
        group(i)(j:)=group(i)(j+1:)
        IF (rjg(n,1)==0) THEN
          rjg(n,1) = i
        ELSE
          rjg(n,2) = i
        ENDIF
      ENDIF
    ENDDO
  ENDDO

END SUBROUTINE rjgrm

!-----------------------------------------------------------------
! PURPOSE: Add ring-joining numerical characters to group strings 
!          at nodes given by rj.                                  
!-----------------------------------------------------------------
SUBROUTINE rjsadd(nring,chem,rjs)
  USE keyparameter, ONLY:digit
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nring
  INTEGER, INTENT(IN) :: rjs(:,:)  ! character # of ring-joining pairs for each ring
  CHARACTER(LEN=*),INTENT(INOUT) :: chem  ! chemical formula

  INTEGER ::   n,i,j

! loop in reverse order => chars in numerical order if >1 exist at any node
  loopr: DO n=nring,1,-1
    loopp: DO j=2,1,-1
      i = rjs(n,j)
      IF (i==0) CYCLE loopp
        chem(i+1:) = chem(i:)
        chem(i:i)=digit(n)
    ENDDO loopp
  ENDDO loopr

END SUBROUTINE rjsadd

!----------------------------------------------------------------------
! PURPOSE: Strip ring-joining numerical characters from chem strings.  
!          Subroutine rjsadd puts them back.                           
!----------------------------------------------------------------------
SUBROUTINE rjsrm(nring,chem,rjs)
  USE keyparameter, ONLY: digit
  IMPLICIT NONE

  INTEGER, INTENT(IN) ::  nring          ! 
  CHARACTER(LEN=*),INTENT(INOUT) :: chem ! chemical name 
  INTEGER, INTENT(OUT) :: rjs(:,:)  ! character # of ring-joining pairs for each ring

  INTEGER  ::  n,i,j,ptr,lenchem

  rjs(:,:)=0  ; lenchem=INDEX(chem,' ')-1
  ringloop : DO n=1,nring
    pairloop: DO j=1,2
      ptr = 1
      DO i=1,lenchem
        IF(chem(i:i)=='C' .OR. chem(i:i)=='c' .OR. chem(i:i)=='-')THEN
          ptr=i+1 
          IF (chem(ptr:ptr)=='d') ptr=ptr+1
          IF (chem(ptr-1:ptr)=='-O') ptr=ptr+1
          IF (chem(ptr:ptr)==digit(n)) THEN ! found ring-joining character
            rjs(n,j)=ptr
            chem(ptr:) = chem(ptr+1:)
            CYCLE pairloop ! exit
          ENDIF
        ENDIF
      ENDDO
    ENDDO  pairloop
  ENDDO ringloop

END SUBROUTINE rjsrm

END MODULE rjtool

