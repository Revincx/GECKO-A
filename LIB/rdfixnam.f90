!*****************************************************************
! PURPOSE: Read the file prescibing (forcing) the name of some
! particular species. Needed to systematically set the same name
! for a given species.
!                                                                 
! Data are store in the fixed name (fn) database:                                                         
!   - nfn       : total nb. of species having a fixed name   
!   - namfn(i)  : table of the fixed name (6 character)       
!   - chemfn(i) : formula corresponding the species  withfixed name                               
!*****************************************************************
SUBROUTINE rdfixnam()
  USE keyparameter, ONLY: mxldi,mxlco,mxlfo,tfu1,dirgecko
  USE normchem, ONLY: stdchm
  USE database, ONLY: mxfn,nfn,namfn,chemfn
  USE searching, ONLY: srh5
  USE sortstring, ONLY:sort_string   ! to sort chemical
  
  IMPLICIT NONE

  INTEGER :: i, ilin
  INTEGER :: ierr
  CHARACTER(LEN=512) :: filename
  CHARACTER(LEN=mxldi) :: line
  CHARACTER(LEN=mxlco) :: tempnam(mxfn)
  CHARACTER(LEN=mxlfo) :: tempchem(mxfn)

  nfn=0 ; namfn(:)=' ' ; chemfn(:)=' ' ; tempnam(:)=' ' ; tempchem(:)=' '

  filename=TRIM(dirgecko)//'DATA/fixedname.dat'
  OPEN(tfu1,FILE=filename,STATUS='OLD')

  rdloop: DO 
    READ(tfu1,'(a)',IOSTAT=ierr) line
    IF (ierr /= 0) THEN
      WRITE (6,*) '--error-- in rdfixnam. Missing keyword END ?'
      STOP "in rdfixnam, while reading inputs"
    ENDIF
    IF (line(1:1)=='!') CYCLE rdloop
    IF (line(1:3)=='END') EXIT rdloop

    nfn = nfn + 1
    IF (nfn>=mxfn) THEN
      WRITE(6,'(a)') '--error--, in rdfixnam: nfn > mxfn'
      STOP "in rdfixnam"
    ENDIF  
!    READ (line(1:n-1),'(a6,3x,a)') tempnam(nfn), tempchem(nfn)
    READ (line,*,IOSTAT=ierr) tempnam(nfn), tempchem(nfn)
    IF (ierr /= 0) THEN
      WRITE (6,*) '--error-- in rdfixnam, while reading data:',TRIM(line)
      STOP "in rdfixnam, while reading inputs"
    ENDIF
    CALL stdchm(tempchem(nfn))
  ENDDO rdloop
  CLOSE(tfu1)

! sort the formula
  chemfn(:)=tempchem(:)
  CALL sort_string(chemfn(1:nfn)) 

! check for duplicate 
  ierr=0
  IF (nfn>1) THEN
    DO i=1,nfn-1
      IF (chemfn(i)==chemfn(i+1)) THEN
        WRITE(6,*) '--error--, while reading file: ',TRIM(filename)
        WRITE(6,*)  'Following species identified 2 times: ',TRIM(chemfn(i))
        ierr=1
      ENDIF
    ENDDO
  ENDIF
  IF (ierr/=0) STOP "in rdfixnam" 
  
! sort the tables according to formula
  DO i=1,nfn
    ilin=srh5(tempchem(i),chemfn,nfn)
    IF (ilin <= 0) THEN
      WRITE(6,*) '--error--, while sorting species in: ',TRIM(filename)
      WRITE(6,*) 'species "lost" after sorting the list: ', TRIM(tempchem(i))
      STOP "in rdfixnam" 
    ENDIF
    namfn(ilin)=tempnam(i)
  ENDDO

END SUBROUTINE rdfixnam
