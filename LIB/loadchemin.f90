MODULE loadchemin
IMPLICIT NONE
CONTAINS

!=======================================================================
! PURPOSE: read the list of primary species for which the chemical
! scheme must be generated.                                
!=======================================================================
SUBROUTINE rdchemin(filename,ninp,input)
  USE keyparameter, ONLY: mxlfo,mxnode,tfu1
  USE atomtool, ONLY: cnum, onum
  USE spsptool, ONLY: chcksp
  IMPLICIT NONE

  CHARACTER(LEN=*),INTENT(IN) :: filename   ! name of the file to be read 
  CHARACTER(LEN=*),INTENT(OUT):: input(:)   ! table of the primary species 
  INTEGER,INTENT(OUT) :: ninp               ! # of "primary" species 

  INTEGER,PARAMETER :: lenline=200
  CHARACTER(LEN=lenline) :: line
  CHARACTER(LEN=mxlfo) :: chem,tpchem
  INTEGER :: nc,ierr,ilin
  INTEGER :: ic, io

  ninp=0  ;  input(:)=' '

  OPEN(tfu1,FILE=filename, FORM='FORMATTED',STATUS='OLD')

! read the data
  ilin=0
  rdloop: DO
    ilin=ilin+1
    READ (tfu1,'(a)',IOSTAT=ierr) line
    IF (ierr/=0) THEN
      WRITE(6,*) '--error--, while reading file: ',TRIM(filename)
      WRITE(6,*) 'at line number: ',ilin
      WRITE(6,*) 'keyword "END" missing ?'
      STOP "in rddict"
    ENDIF
    IF (line(1:3)=='END') EXIT rdloop
    IF (line(1:1)=='!') CYCLE rdloop

! For SAR evaluation: store species with "*", not to be treated but 
! to keep original sorting
    IF (line(1:1)=='*') THEN
      ninp = ninp+1 ; input(ninp) = TRIM(line)
      CYCLE rdloop
    ENDIF


! check that the species can be handled in the generator
    nc=INDEX(line,' ')
    IF (nc > mxlfo) THEN
      WRITE(6,'(a)') '--error-- in dchemin, species > lfo for :'
      WRITE(6,'(a)') TRIM(line)
      STOP "in rdchemin"
    ENDIF
    chem=line(1:nc)

    CALL chcksp(chem)
       
! If the species is a special species (must start with #) ...
    IF (chem(1:1)=='#') THEN
      ninp = ninp+1
      input(ninp) = chem
      CYCLE rdloop
 
! If not a "special species", then must start with a C 
    ELSE IF ( (chem(1:1)=='C')  .OR. (chem(1:4)=='=Cd1') .OR. &
             (chem(1:2)=='-O') .OR. (chem(1:1)=='c') ) THEN

! get the number of group (C) in the molecule
      ic = cnum(chem)  ;  io = onum(chem)
      ic = ic + io
      IF (ic > mxnode) THEN
        WRITE(6,*) '--error-- in dchemin, node in species > mxnode: '
        WRITE(6,*) TRIM(chem)
        STOP "in rdchemin"
      ENDIF

      IF (ic==1) THEN
        tpchem=chem ; nc=INDEX(tpchem,' ')
        chem='* '; chem(2:nc+1)=tpchem(1:nc)        
!        CYCLE rdloop    ! cycle if only 1 carbon
      ENDIF 
      ninp = ninp+1
      input(ninp) = chem

! Otherwise, species can not be managed by gecko
    ELSE 
      WRITE(6,*) '--error-- in rdchemin, formula must start with C or #: '
      WRITE(6,*) TRIM(chem)
      STOP "in rdchemin"
    ENDIF

  ENDDO rdloop
  CLOSE(tfu1)

! check that # of input does not exceed max #
  IF (ninp > SIZE(input)) THEN
    WRITE(6,'(a)') '--error-- in rdchemin. Too many input species'
    WRITE(6,'(a)') 'max # of input species = ', SIZE(input)
    STOP "in rdchemin"
  ENDIF

END SUBROUTINE rdchemin

!=======================================================================
! PURPOSE: Handle primary species. The subroutine checks the species, 
! gives it a name (if not already in the dictionary), updates the 
! dictionary and puts the species at the beginning of the stack (i.e. 
! in holdvoc(1)). Stack variables are updated.
!=======================================================================
SUBROUTINE in1chm(chem,xname)
  USE keyparameter, ONLY: mxlfl
  USE keyflag, ONLY: sar_only_fg,isomerfg
  USE atomtool, ONLY: cnum
  USE searching, ONLY: srh5, srch
  USE normchem, ONLY: stdchm
  USE database, ONLY:  nspsp, dictsp  ! special species input
  USE dictstackdb, ONLY: nrec,dict,namlst,stabl,level,dbrch,lotopstack,diccri
  USE dictstacktool, ONLY:loader
  USE namingtool, ONLY: naming
  USE toolbox, ONLY: stoperr
  IMPLICIT NONE

  CHARACTER(LEN=*),INTENT(INOUT) :: chem  ! formula of the primary species to be added in the stack
  CHARACTER(LEN=*),INTENT(OUT)   :: xname ! short name of chem

  CHARACTER(LEN=mxlfl) :: fgrp
  INTEGER  :: i,nca, dicptr, namptr
  INTEGER  :: ipos

  CHARACTER(LEN=7),PARAMETER :: progname='in1chm '
  CHARACTER(LEN=70)          :: mesg

  xname = ' ' ; lotopstack=.FALSE.  ! default is fifo stack

! THE SPECIES IS A SPECIAL SPECIES 
! ------------------------------------
  IF (chem(1:1)=='#') THEN
    dicptr = srch(nrec,chem,dict)

! already recorded
    IF (dicptr > 0) THEN
      level = -1
      WRITE(6,*) '::cheminput found in dictionary: ',TRIM(chem)
      xname=dict(dicptr)(1:6)
      IF (.NOT.sar_only_fg) RETURN  ! for SAR assessment, need to run all entries
    ENDIF

! search the name
    ipos=0
    nsploop: DO i=1,nspsp
      IF (dictsp(i)(10:130)==chem) THEN
        ipos=i
        xname=dictsp(i)(1:6)
        namptr = srh5(xname,namlst,nrec)
        fgrp=dictsp(i)(132:)
        namptr=-namptr
        EXIT nsploop
      ENDIF
    ENDDO nsploop

    IF (ipos==0) THEN
      PRINT*, "-error- species not found in in1chm"
      STOP "in in1chm"
    ENDIF

  ELSE 
! THE SPECIES IS A REGULAR SPECIES
! -----------------------------------

! get the number of carbons in species - must be greater than 1
    nca = cnum(chem)
    IF (nca <= 1) THEN
      mesg="Species with < 2 carbons "
      CALL stoperr(progname,mesg,chem)
    ENDIF
      
! standardize chem
    CALL stdchm(chem)
        
! search if chem already exists in the dictionary 
    dicptr = srch(nrec,chem,dict)

! If already known, return
    IF (dicptr > 0) THEN
      level = -1
      WRITE(6,*) 'cheminput found in dictionary: ',TRIM(chem)
      xname=dict(dicptr)(1:6)
      IF (.NOT.sar_only_fg) RETURN  ! for SAR assessment, need to run all entries
    ENDIF

! if new species, get the short name for the species (xname) and 
! position after which it must be added in the namlst table (namptr)
    CALL naming(chem,namptr,xname,fgrp)

  ENDIF

! ADD SPECIES TO THE STACK AND RAISE TABLES
! -----------------------------------------
! raise the record counters
  nrec = nrec + 1  
  IF (nrec > SIZE(dict)) THEN
    WRITE (6,*) '--error-- in in1chem. Too many species (nrec > mxspe)'
    STOP "in in1chm"
  ENDIF 

! add new name and raise the name array, insert new name
  namptr = namptr + 1
  namlst(namptr+1:nrec+1)=namlst(namptr:nrec)
  namlst(namptr) = xname

! raise upper part of dictionary arrays and branching array 
! insert new line for the new species 
  dicptr = ABS(dicptr) + 1
  dict(dicptr+1:nrec+1)=dict(dicptr:nrec)
  dbrch(dicptr+1:nrec+1)=dbrch(dicptr:nrec)
  WRITE(dict(dicptr),'(a6,3x,a120,2x,a15)')  xname, chem, fgrp
  dbrch(dicptr)  = 1.   ! initialize dbrch for the parent

  IF (isomerfg) THEN
    diccri(dicptr+1:nrec+1,:) = diccri(dicptr:nrec,:)
    diccri(dicptr,:)=0
  ENDIF

! load species in the stack
  level = -1                              ! loader adds 1 to level
  IF (INDEX(chem,'.')== 0)THEN ; stabl=-1 ! loader adds 1 for non radical ...
  ELSE                         ; stabl=0  ! ... but 0 to radicals
  ENDIF
  CALL loader(chem,xname)

  WRITE(6,*) '::cheminput added to dictionary::',TRIM(chem)
END SUBROUTINE in1chm

END MODULE loadchemin
