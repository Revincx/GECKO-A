MODULE rdchemprop
IMPLICIT NONE
CONTAINS

!SUBROUTINE rdhenry()
!SUBROUTINE rdhydrat()
!SUBROUTINE rddatabase(filename,ncol,ndat,chemframe,dataframe,comframe)

! ======================================================================
! Read the effective Henry's law coefficient and store the data in the
! "hlcdb" database (see module "database").
! ======================================================================
SUBROUTINE rdhenry()
  USE database, ONLY: nhlcdb,hlcdb_chem,hlcdb_dat,hlcdb_com
  USE keyparameter, ONLY: dirgecko
  IMPLICIT NONE
  CHARACTER(LEN=200) filename

  nhlcdb=0 ; hlcdb_chem(:)=' ' ; hlcdb_dat(:,:)=0.0 ; hlcdb_com(:,:)=' '
  filename=TRIM(dirgecko)//'DATA/henry_const.dat'
  CALL rddatabase(filename,3, nhlcdb, hlcdb_chem, hlcdb_dat, hlcdb_com)
END SUBROUTINE rdhenry

! ======================================================================
! Read the hydration constant and store the data in the "khydb" 
! database (see module "database").
! ======================================================================

SUBROUTINE rdhydrat()
  USE database, ONLY: nkhydb, khydb_chem, khydb_dat, khydb_com
  USE keyparameter, ONLY: dirgecko
  IMPLICIT NONE
  CHARACTER(LEN=200) filename

  nkhydb=0 ; khydb_chem(:)=' ' ;  khydb_dat(:,:)=0. ; khydb_com(:,:)=' '
  filename=TRIM(dirgecko)//'DATA/hydrat_const.dat'
  CALL rddatabase(filename,1, nkhydb, khydb_chem, khydb_dat, khydb_com)
END SUBROUTINE rdhydrat

! ======================================================================
! Purpose: read a database file for a given set of species (e.g. 
! Henry's law coefficient, K_hydration ...)
! ======================================================================
SUBROUTINE rddatabase(filename,ncol,ndat,chemframe,dataframe,comframe)
  USE keyparameter, ONLY: tfu1
  USE sortstring,ONLY:sort_string   ! to sort chemical
  USE searching, ONLY: srh5   ! to seach in a sorted list
  USE normchem, ONLY:stdchm   ! to standardize formulas 
  USE references, ONLY: nreac_info, code  ! available code
  IMPLICIT NONE

  CHARACTER(LEN=*),INTENT(IN)  :: filename  ! file to be read
  INTEGER,INTENT(IN) :: ncol                ! # of column to be read (also size(dataframe,2))
  INTEGER,INTENT(OUT) :: ndat               ! # of data in the database (bd)
  CHARACTER(LEN=*),INTENT(OUT) :: chemframe(:) ! name od the species in the db
  REAL,INTENT(OUT) :: dataframe(:,:)           ! data to be read 
  CHARACTER(LEN=*),INTENT(OUT) :: comframe(:,:)! comment (source) of the data

! local :
  CHARACTER(LEN=200) :: line                ! string in the db
  INTEGER       :: i,j,n1,n2,ierr,ilin      !
  INTEGER       :: ipc(ncol+3)              ! position of the expected ";"  in string
  INTEGER       :: maxdat                   ! max # of data allowed
  CHARACTER(LEN=LEN(comframe(1,1))) :: tpcom   ! a temporary comment

! temporary copy of the output tables (need to sort tables)  
  CHARACTER(LEN=LEN(chemframe(1))) :: tpchemframe(SIZE(chemframe))
  REAL :: tpdataframe(SIZE(dataframe,1),SIZE(dataframe,2))
  CHARACTER(LEN=LEN(comframe(1,1))) :: tpcomframe(SIZE(comframe,1),SIZE(comframe,2))
  INTEGER :: nsep
  
! initialize
  maxdat=SIZE(dataframe,1)  ;  ndat=0 ; dataframe(:,:)=0.
  chemframe(:)=' ' ; comframe(:,:)=' '
  tpchemframe(:)=chemframe(:)  ;  tpdataframe(:,:)=dataframe(:,:)
  tpcomframe(:,:)=comframe(:,:)
  nsep=ncol+3
  
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
      STOP "in rddatabase" 
    ENDIF
    IF (line(1:1)=='!') CYCLE rdloop
    IF (line(1:3)=='END') EXIT rdloop
    IF (INDEX(line,'!')/=0) line=line(1:INDEX(line,'!')-1) ! rm txt after "!"
    
! get the position of the ";" separating the field
    n1=0
    DO i=1,LEN_TRIM(line)
      IF (line(i:i)==';') THEN ; n1=n1+1 ; ipc(n1)=i ; ENDIF
      IF (n1==nsep) EXIT
    ENDDO
    IF (n1 /= nsep) THEN 
      WRITE(6,*) '--error--, while reading file: ',TRIM(filename)
      WRITE(6,*) 'at line number: ',ilin
      WRITE(6,*) 'missing ";" separator. n1= ', n1
      STOP "in rddatabase" 
    ENDIF
    IF (INDEX(line(ipc(nsep)+1:),";")/=0) THEN
      WRITE(6,*) '--error--, while reading file: ',TRIM(filename)
      WRITE(6,*) 'at line number: ',ilin
      WRITE(6,*) 'too many ";" separator. n1 > ',nsep
      STOP "in rddatabase" 
    ENDIF

! read the species
    ndat=ndat+1
    IF (ndat > maxdat) THEN
      WRITE(6,*) '--error--, while reading file: ',TRIM(filename)
      WRITE(6,*) 'too many species (check table size). maxdat= ',maxdat
      STOP "in rddatabase" 
    ENDIF
    READ(line(1:ipc(1)-1),'(a)',IOSTAT=ierr) tpchemframe(ndat)
    IF (ierr/=0) THEN
      WRITE(6,*) '--error--, while reading file: ',TRIM(filename)
      WRITE(6,*) 'at line number: ',ilin
      WRITE(6,*) 'while reading species: ',line(1:ipc(1)-1)
      STOP "in rddatabase" 
    ENDIF

! read dataframe coefficients
    DO i=1,ncol
      READ(line(ipc(i)+1:ipc(i+1)-1),*,IOSTAT=ierr) tpdataframe(ndat,i)
      IF (ierr/=0) THEN
        WRITE(6,*) '--error--, while reading file: ',TRIM(filename)
        WRITE(6,*) 'while reading rate coef. at: ', line(ipc(i)+1:ipc(i+1)-1)
        WRITE(6,*) 'at line number: ',ilin
        STOP "in rddatabase" 
      ENDIF
    ENDDO

! read the comments
    DO i=1,3
      tpcom=ADJUSTL(line(ipc(i+ncol)+1:)) 
      n2=INDEX(tpcom,';')
      IF (n2/=0) tpcom(n2:)='           '  ! remove string not part of comment
      IF (tpcom(1:1)==' ') CYCLE           ! nothing to read
      READ(tpcom,*,IOSTAT=ierr) tpcomframe(ndat,i)
      IF (ierr/=0) THEN
        WRITE(6,*) '--error--, while reading file: ',TRIM(filename)
        WRITE(6,*) 'while reading comment (code): ', tpcom
        WRITE(6,*) 'at line number: ',ilin
        WRITE(6,*) TRIM(line)
        STOP "in rddatabase" 
      ENDIF
    ENDDO
  ENDDO rdloop  
  CLOSE(tfu1)

! check and normalize the formula
  DO i=1,ndat  ;  CALL stdchm(tpchemframe(i))  ;  ENDDO

! sort the formula
  chemframe(:)=tpchemframe(:)
  CALL sort_string(chemframe(1:ndat)) 

! check for duplicate 
  ierr=0
  DO i=1,ndat-1
    IF (chemframe(i)==chemframe(i+1)) THEN
      WRITE(6,*) '--error--, while reading file: ',TRIM(filename)
      WRITE(6,*)  'Following species identified 2 times: ',TRIM(chemframe(i))
      ierr=1
    ENDIF
  ENDDO
  IF (ierr/=0) STOP "in rddatabase" 
  
! sort the tables according to formula
  DO i=1,ndat
    ilin=srh5(tpchemframe(i),chemframe,ndat)
    IF (ilin <= 0) THEN
      WRITE(6,*) '--error--, while sorting species in: ',TRIM(filename)
      WRITE(6,*) 'species "lost" after sorting the list: ', TRIM(tpchemframe(i))
      STOP "in rddatabase" 
    ENDIF
    dataframe(ilin,:)=tpdataframe(i,:) ; comframe(ilin,:)=tpcomframe(i,:)
  ENDDO

! check the availability of the codes
  ierr=0
  DO i=1,ndat
    DO j=1,SIZE(comframe,2)
      IF (comframe(i,j)(1:1) /= ' ') THEN 
        n1=srh5(comframe(i,j),code,nreac_info)
        IF (n1 <= 0) THEN
          WRITE(6,*) '--error--, while reading file ',TRIM(filename)
          WRITE(6,*) 'unknown comment: ', comframe(i,j)," ",n1
          ierr=1
        ENDIF
      ENDIF
    ENDDO
  ENDDO
  IF (ierr/=0) STOP "in rddatabase" 
  
END SUBROUTINE rddatabase

END MODULE rdchemprop
