MODULE reac_infotool
IMPLICIT NONE
CONTAINS

! ======================================================================
SUBROUTINE rdreac_info()
  USE keyparameter, ONLY: tfu1,dirgecko
  USE references, ONLY: mxreac_info,mxlcod,mxlreac_info,nreac_info,code,reac_info
  USE sortstring              ! to sort the reac_info codes
  USE searching, ONLY: srh5   ! to search in a sorted list
  IMPLICIT NONE

  CHARACTER(LEN=256) :: filename
  INTEGER       :: i,n1,ierr,ilin,lcode,cblank
  CHARACTER(LEN=420) :: line                   ! string to be read
  CHARACTER(LEN=mxlcod+10) :: leftline         ! code side string
  CHARACTER(LEN=mxlreac_info+5):: rightline        ! reac_info side string
  CHARACTER(LEN=mxlcod)    :: tpcode(mxreac_info)  ! temporary copy (unsorted) of code
  CHARACTER(LEN=mxlreac_info)  :: tpreac_info(mxreac_info) ! temporary copy (unsorted) of code
  
  nreac_info=0 ; tpcode(:)=' ' ; tpreac_info(:)=' '
  
! open the file
  filename=TRIM(dirgecko)//'DATA/references.dat'
  OPEN(tfu1,FILE=filename, FORM='FORMATTED',STATUS='OLD',IOSTAT=ierr)
  IF (ierr/=0) THEN
    WRITE(6,*) '--error--, while trying to open: ',TRIM(filename)
    STOP "in rdreac_info" 
  ENDIF

! Read the file
! ---------------
  ilin=0
  rdloop: DO
    leftline=' ' ; rightline=' '
    ilin = ilin+1
    READ(tfu1,'(a)',IOSTAT=ierr) line
    IF (ierr/=0) THEN
      WRITE(6,*) '--error--, while reading file: ',TRIM(filename)
      WRITE(6,*) 'at line number: ',ilin
      WRITE(6,*) 'keyword "END" missing ?'
      STOP "in rdreac_info" 
    ENDIF
    IF (line(1:1)=='!') CYCLE rdloop
    IF (line(1:3)=='END') EXIT rdloop
    
! read line
    line=ADJUSTL(line) ; n1=INDEX(line,':') 
    IF (n1<=1) THEN 
      WRITE(6,*) '--error--, while reading file: ',TRIM(filename)
      WRITE(6,*) 'at line number: ',ilin
      WRITE(6,*) 'Expect ":" separator in line: ',TRIM(line)
      STOP "in rdreac_info" 
    ENDIF

! read the code
    leftline=line(1:n1-1) ; rightline=ADJUSTL(line(n1+1:))
    lcode=LEN_TRIM(leftline)              ! check the length of the code
    IF (lcode > mxlcod) THEN
      WRITE(6,*) '--error--, while reading file: ',TRIM(filename)
      WRITE(6,*) 'at line number: ',ilin
      WRITE(6,*) 'Length of the code exceed mxlcod in: ',TRIM(line)
      STOP "in rdreac_info" 
    ENDIF
    cblank=INDEX(leftline(1:lcode),' ')   ! check for " " char
    IF (cblank /=0 ) THEN
      WRITE(6,*) '--error--, while reading file: ',TRIM(filename)
      WRITE(6,*) 'at line number: ',ilin
      WRITE(6,*) 'Unexpected " " char in code: ',TRIM(line)
      STOP "in rdreac_info" 
    ENDIF
    nreac_info=nreac_info+1
    IF (nreac_info > mxreac_info) THEN             ! check # of reac_infos
      WRITE(6,*) '--error--, while reading file: ',TRIM(filename)
      WRITE(6,*) 'too many reac_infos (check table size). mxreac_info= ',mxreac_info
      STOP "in rdreac_info" 
    ENDIF
    tpcode(nreac_info)=leftline(1:mxlcod)

! read the reac_info
    IF (LEN_TRIM(leftline)>mxlreac_info) THEN
      WRITE(6,*) '--error--, while reading file: ',TRIM(filename)
      WRITE(6,*) 'at line number: ',ilin
      WRITE(6,*) 'Length of the code exceed mxlreac_info in: ',TRIM(line)
      STOP "in rdreac_info" 
    ENDIF
    tpreac_info(nreac_info)=rightline(1:mxlreac_info)
  ENDDO rdloop
  CLOSE(tfu1)

! load sorted tables tables
! -------------------------

! sort the codes
  code(:)=tpcode(:)
  CALL sort_string(code(1:nreac_info)) 

! check for duplicate 
  ierr=0
  DO i=1,nreac_info-1
    IF (code(i)==code(i+1)) THEN
      WRITE(6,*) '--error--, while reading file: ',TRIM(filename)
      WRITE(6,*)  'Following code identified 2 times: ',TRIM(code(i))
      ierr=1
    ENDIF
  ENDDO
  IF (ierr/=0) STOP "in rdreac_info" 
  
! sort the tables according to formula
  DO i=1,nreac_info
    ilin=srh5(tpcode(i),code,nreac_info)
    IF (ilin <= 0) THEN
      WRITE(6,*) '--error--, while sorting reac_info in: ',TRIM(filename)
      WRITE(6,*) 'reac_info "lost" after sorting the list: ', TRIM(tpcode(i))
      STOP "in rdreac_info"
    ENDIF
    reac_info(ilin)=tpreac_info(i)
  ENDDO

END SUBROUTINE rdreac_info ! -----------------------------------------------

CHARACTER*(mxlreac_info) FUNCTION fullref(shortcode)
 USE references
 IMPLICIT NONE
 CHARACTER*(mxlcod) :: shortcode
 INTEGER            :: i
 
 fullref=' '
 DO i=1,mxreac_info
   IF (code(i)==shortcode) THEN
     fullref=reac_info(i)
     RETURN
   ENDIF 
 ENDDO
END FUNCTION

END MODULE reac_infotool
