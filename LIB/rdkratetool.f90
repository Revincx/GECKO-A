MODULE rdkratetool
IMPLICIT NONE
CONTAINS

! SUBROUTINE rdratedb(filename,ndat,oxchem,k298,arrh,oxcom)
! SUBROUTINE rdkoxfiles()

! ======================================================================
! Purpose: 
! ======================================================================
SUBROUTINE rdkoxfiles()
  USE keyparameter, ONLY: mxlfo,dirgecko
  USE database, ONLY: &
    mxkdb, & ! max # of element in the database
    nkohdb,kohdb_chem,kohdb_298,kohdb_arr,kohdb_com, &   ! K_OH database
    nko3db,ko3db_chem,ko3db_298,ko3db_arr,ko3db_com, &   ! K_O3 database
    nkno3db,kno3db_chem,kno3db_298,kno3db_arr,kno3db_com ! K_NO3 database
  IMPLICIT NONE

  CHARACTER(LEN=512) :: filename
  
  WRITE (6,*) '  ...reading OH rate constants'
  filename=TRIM(dirgecko)//'DATA/koh_rate.dat'
  CALL rdratedb(filename,nkohdb,kohdb_chem,kohdb_298,kohdb_arr,kohdb_com)

  WRITE (6,*) '  ...reading O3 rate constants'
  filename=TRIM(dirgecko)//'DATA/ko3_rate.dat'
  CALL rdratedb(filename,nko3db,ko3db_chem,ko3db_298,ko3db_arr,ko3db_com)

  WRITE (6,*) '  ...reading NO3 rate constants'
  filename=TRIM(dirgecko)//'DATA/kno3_rate.dat'
  CALL rdratedb(filename,nkno3db,kno3db_chem,kno3db_298,kno3db_arr,kno3db_com)
END SUBROUTINE rdkoxfiles

! ======================================================================
! Purpose: read a database file for a given set of reactions (e.g. 
! VOC+OH, VOC+O3, ...)
! ======================================================================
SUBROUTINE rdratedb(filename,ndat,oxchem,k298,arrh,oxcom)
  USE keyparameter, ONLY: tfu1
  USE sortstring              ! to sort chemical
  USE searching, ONLY: srh5   ! to seach in a sorted list
  USE normchem, ONLY:stdchm   ! to standardize formulas 
  USE references, ONLY: nreac_info, code  ! available code
  IMPLICIT NONE

  CHARACTER(LEN=*),INTENT(IN)  :: filename  ! file to be read
  INTEGER,INTENT(OUT) :: ndat               ! # of data in the database (bd)
  CHARACTER(LEN=*),INTENT(OUT) :: oxchem(:) ! name od the species in the db
  REAL,INTENT(OUT) :: k298(:)               ! rate constant k @ 298K
  REAL,INTENT(OUT) :: arrh(:,:)             ! arrhenius parameter 
  CHARACTER(LEN=*),INTENT(OUT) :: oxcom(:,:)! comment (source) of the rate const.

! local :
  CHARACTER(LEN=200) :: line                ! string in the db
  INTEGER       :: i,j,n1,n2,ierr,ilin
  INTEGER       :: ipc(7)                   ! expected # of ";" char in string
  INTEGER       :: maxdat                   ! max # of data allowed
  CHARACTER(LEN=LEN(oxcom(1,1))) :: tpcom   ! a temporary comment
! temporary copy of the output tables (need to sort tables)  
  CHARACTER(LEN=LEN(oxchem(1))) :: tpoxchem(SIZE(oxchem))
  REAL :: tpk298(SIZE(arrh,1))
  REAL :: tparrh(SIZE(arrh,1),SIZE(arrh,2))
  CHARACTER(LEN=LEN(oxcom(1,1))) :: tpoxcom(SIZE(oxcom,1),SIZE(oxcom,2))
  
! initialize
  maxdat=SIZE(arrh)
  ndat=0 ; k298(:)=0. ; arrh(:,:)=0.
  oxchem(:)=' ' ; oxcom(:,:)=' '
  tpoxchem(:)=oxchem ; tpk298(:)=k298(:) ; tparrh(:,:)=arrh(:,:)
  tpoxcom(:,:)=oxcom(:,:)
  
! open the file
  OPEN(tfu1,FILE=filename, FORM='FORMATTED',STATUS='OLD')
  
! read the file
  ilin=0
  rdloop: DO
    ilin = ilin+1
    READ(tfu1,'(a)',IOSTAT=ierr) line
    IF (ierr/=0) THEN
      WRITE(6,*) '--error--, while reading file: ',TRIM(filename)
      WRITE(6,*) 'at line number: ',ilin
      WRITE(6,*) 'keyword "END" missing ?'
      STOP "in rdratedb" 
    ENDIF
    IF (line(1:1)=='!') CYCLE rdloop
    IF (line(1:3)=='END') EXIT rdloop

! remove anything after "!" (not necessary as 1st character)
    n1=INDEX(line,"!")
    IF (n1>0) line(n1:)=' '
    
! get the position of the ";" separating the field
    n1=0
    DO i=1,LEN_TRIM(line)
      IF (line(i:i)==';') THEN ; n1=n1+1 ; ipc(n1)=i ; ENDIF
      IF (n1==7) EXIT
    ENDDO
    IF (n1 /= 7) THEN 
      WRITE(6,*) '--error--, while reading file: ',TRIM(filename)
      WRITE(6,*) 'at line number: ',ilin
      WRITE(6,*) 'missing ";" separator. n1= ', n1
      STOP "in rdratedb" 
    ENDIF
    IF (INDEX(line(ipc(7)+1:),";")/=0) THEN
      WRITE(6,*) '--error--, while reading file: ',TRIM(filename)
      WRITE(6,*) 'at line number: ',ilin
      WRITE(6,*) 'too many ";" separator. n1 > 7'
      STOP "in rdratedb" 
    ENDIF

! read the species
    ndat=ndat+1
    IF (ndat > maxdat) THEN
      WRITE(6,*) '--error--, while reading file: ',TRIM(filename)
      WRITE(6,*) 'too many species (check table size). maxdat= ',maxdat
      STOP "in rdratedb" 
    ENDIF
    READ(line(1:ipc(1)-1),'(a)',IOSTAT=ierr) tpoxchem(ndat)
    IF (ierr/=0) THEN
      WRITE(6,*) '--error--, while reading file: ',TRIM(filename)
      WRITE(6,*) 'at line number: ',ilin
      WRITE(6,*) 'while reading species: ',line(1:ipc(1)-1)
      STOP "in rdratedb" 
    ENDIF

! read the rate constant @ 298K
    READ(line(ipc(1)+1:ipc(2)-1),*,IOSTAT=ierr) tpk298(ndat)
    IF (ierr/=0) THEN
      WRITE(6,*) '--error--, while reading file: ',TRIM(filename)
      WRITE(6,*) 'while reading k298 at: ',line(ipc(1)+1:ipc(2)-1)
      WRITE(6,*) 'at line number: ',ilin
      STOP "in rdratedb" 
    ENDIF

! read arrhenius coefficients
    DO i=1,3
      READ(line(ipc(i+1)+1:ipc(i+2)-1),*,IOSTAT=ierr) tparrh(ndat,i)
      IF (ierr/=0) THEN
        WRITE(6,*) '--error--, while reading file: ',TRIM(filename)
        WRITE(6,*) 'while reading rate coef. at: ', line(ipc(i+1)+1:ipc(i+2)-1)
        WRITE(6,*) 'at line number: ',ilin
        STOP "in rdratedb" 
      ENDIF
    ENDDO

! read the comments
    DO i=1,3
      tpcom=ADJUSTL(line(ipc(i+4)+1:)) 
      n2=INDEX(tpcom,';')
      IF (n2/=0) tpcom(n2:)='           '  ! remove string not part of comment
      IF (tpcom(1:1)==' ') CYCLE           ! nothing to read
      READ(tpcom,*,IOSTAT=ierr) tpoxcom(ndat,i)
      IF (ierr/=0) THEN
        WRITE(6,*) '--error--, while reading file: ',TRIM(filename)
        WRITE(6,*) 'while reading comment (code): ', tpcom
        WRITE(6,*) 'at line number: ',ilin
        WRITE(6,*) TRIM(line)
        STOP "in rdratedb" 
      ENDIF
    ENDDO
  ENDDO rdloop  
  CLOSE(tfu1)

! check and normalize the formula
  DO i=1,ndat  ;  CALL stdchm(tpoxchem(i))  ;  ENDDO

! sort the formula
  oxchem(:)=tpoxchem(:)
  CALL sort_string(oxchem(1:ndat)) 

! check for duplicate 
  ierr=0
  DO i=1,ndat-1
    IF (oxchem(i)==oxchem(i+1)) THEN
      WRITE(6,*) '--error--, while reading file: ',TRIM(filename)
      WRITE(6,*)  'Following species identified 2 times: ',TRIM(oxchem(i))
      ierr=1
    ENDIF
  ENDDO
  IF (ierr/=0) STOP "in rdratedb" 
  
! sort the tables according to formula
  DO i=1,ndat
    ilin=srh5(tpoxchem(i),oxchem,ndat)
    IF (ilin <= 0) THEN
      WRITE(6,*) '--error--, while sorting species in: ',TRIM(filename)
      WRITE(6,*) 'species "lost" after sorting the list: ', TRIM(tpoxchem(i))
      STOP "in rdratedb" 
    ENDIF
    k298(ilin)=tpk298(i) ; arrh(ilin,:)=tparrh(i,:) ; oxcom(ilin,:)=tpoxcom(i,:)
  ENDDO

! check the availability of the codes
  ierr=0
  DO i=1,ndat
    DO j=1,SIZE(oxcom,2)
      IF (oxcom(i,j)(1:1) /= ' ') THEN 
        n1=srh5(oxcom(i,j),code,nreac_info)
        IF (n1 <= 0) THEN
          WRITE(6,*) '--error--, while reading file ',TRIM(filename)
          WRITE(6,*) 'unknow comment: ', oxcom(i,j)
          ierr=1
        ENDIF
      ENDIF
    ENDDO
  ENDDO
  IF (ierr/=0) STOP "in rdratedb" 
  
END SUBROUTINE rdratedb

END MODULE rdkratetool
