MODULE searching
  IMPLICIT NONE
  CONTAINS

!=======================================================================
! PURPOSE: Binary tree search for string      
! The value returned (srh5) is: 
!    <0 if not found (|srh5|==pointer where aseek should be inserted) 
!    >0 if found (|srh5|==pointer where aseek is inserted in alist)   
!=======================================================================
INTEGER FUNCTION srh5(aseek,alist,nlist)
  IMPLICIT NONE

  INTEGER,INTENT(IN)  ::  nlist           ! # of records in alist
  CHARACTER(LEN=*),INTENT(IN) :: aseek    ! string to look for in alist
  CHARACTER(LEN=*),INTENT(IN) :: alist(:) ! list of strings (e.g. dictionary)

! internal:
  INTEGER :: jhi    ! highest index value in alist                        
  INTEGER :: jlo    ! lowest index value in alist                         
  INTEGER :: jold   ! temporary storage of actual pointer in alist        
  INTEGER :: j      ! actual pointer in alist                             


! initialize:
  srh5 = 0  ;  jold = 0  ;  jlo  = 1  ;  jhi  = nlist + 1

! search loop
  searchloop: DO
    j = (jhi+jlo)/2
    IF (j == jold) EXIT searchloop 

    jold = j
    IF(aseek > alist(j)) THEN
      jlo  = j
      CYCLE searchloop
    ENDIF
    IF(aseek == alist(j)) THEN 
      srh5 = j
      RETURN
    ENDIF  
    jhi  = j

  ENDDO searchloop

! string not found 
  srh5 = -j
  RETURN
  
END FUNCTION srh5

!=======================================================================
! PURPOSE: Binary tree search for formula in dictionary      
! The value returned (srch) is: 
!    <0 if not found (|srch|==pointer where chem should be inserted) 
!    >0 if found (|srch|==pointer where chem is inserted in dict)
! WARNING (BA) : "same" routine as srh5 => merge later    
!=======================================================================
INTEGER FUNCTION srch(nrec,chem,dict)
  IMPLICIT NONE

  INTEGER,INTENT(IN)  ::  nrec           ! # of records in dict
  CHARACTER(LEN=*),INTENT(IN) :: chem    ! formula to be searched in dict.
  CHARACTER(LEN=*),INTENT(IN) :: dict(:) ! dictionary lines

! internal:
  INTEGER :: jhi     ! highest index value in alist                        
  INTEGER :: jlo     ! lowest index value in alist                         
  INTEGER :: jold    ! temporary storage of actual pointer in alist        
  INTEGER :: j       ! actual pointer in alist                             


! initialize:
  srch = 0  ;  jold = 0  ;  jlo  = 1  ;  jhi  = nrec + 1

! search loop
  searchloop: DO
    j = (jhi+jlo)/2
    IF (j == jold) EXIT searchloop 

    jold = j
    IF(chem > dict(j)(10:129)) THEN
      jlo  = j
      CYCLE searchloop
    ENDIF
    IF(chem == dict(j)(10:129)) THEN 
      srch = j
      RETURN
    ENDIF  
    jhi  = j

  ENDDO searchloop

! string not found 
  srch = -j
  RETURN

END FUNCTION srch

!=======================================================================
! PURPOSE: return the rank (position) of an integer provided as input 
! (iseek) in a sorted list of integer (from large to small number). The
! list provided as input (ilist) is expected to be short: a simple 
! search is performed. Turn the function using binary tree search for 
! long list.  
!=======================================================================
INTEGER FUNCTION search_ipos(iseek,ilist)
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: iseek    ! inteer to look for in list
  INTEGER,INTENT(IN) :: ilist(:) ! list of integer

  INTEGER :: i

  search_ipos=0
  IF (iseek<1) THEN
    WRITE(*,*) "--error-- in search_ipos, requested search using <1 int."
    STOP "in search_ipos"
  ENDIF
  DO i=1,SIZE(ilist)
   IF (iseek>ilist(i)) THEN
     search_ipos=i
     RETURN
   ENDIF
 ENDDO
 
  IF (search_ipos<1) THEN
    WRITE(*,*) "--error-- in search_ipos, no slot found"
    STOP "in search_ipos"
  ENDIF

END FUNCTION search_ipos

END MODULE searching
