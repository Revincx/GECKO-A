MODULE maintool
  IMPLICIT NONE
  CONTAINS
!SUBROUTINE initdictstack()
!SUBROUTINE cleanstack()
!SUBROUTINE sortstack(nrec,dict,dbrch,nhldvoc,holdvoc)
!INTEGER FUNCTION nphotogrp(chem)
!SUBROUTINE wrtscreen(nrec,nhldvoc,stabl,chem)
!REAL FUNCTION vaporpressure(chem,bond,group,nring,rjg)
!SUBROUTINE InOrOut(chem,parent,brch,stabl)
!SUBROUTINE cut_path(brch,cut_default,cut_OH,cut_O3,cut_NO3,cut_PAN,&
!                    cut_HV, cut_RO,cut_RO2,cut_RCOO2) 

!=======================================================================
! PURPOSE: Initialise dictionaries and chemical stacks
!=======================================================================
SUBROUTINE initdictstack()
  USE dictstackdb
  IMPLICIT NONE

! initialize dictionaries
  nrec=0      ;  ninorg = 0   ; nwpspe = 0  
  dict(:)=' ' ; namlst(:)=' ' ; inorglst(:)=' '  
  dbrch(:)=0. 
  
! Initialize isomer stack 
  diccri(:,:)=0

! initialize data in the stack
  nhldvoc=0      ;  nhldrad=0
  holdvoc(:)=' ' ;  holdrad(:)=' '
  level=0        ;  stabl=0   
     
! initialize tetrahydrofuran stack
  ncha=0  ; chatab(:)=' ' 
END SUBROUTINE initdictstack  

!=======================================================================
! PURPOSE: Clean the chemical stacks
!=======================================================================
SUBROUTINE cleanstack()
  USE dictstackdb, ONLY: nhldvoc,holdvoc,nhldrad,holdrad
  
  holdvoc(1:nhldvoc+100)=' '  ! big table, useless to erase over all index
  holdrad(:)=' '  ;  nhldvoc=0  ;  nhldrad=0  
END SUBROUTINE cleanstack

!=======================================================================
! PURPOSE: sort the stack
!=======================================================================
SUBROUTINE sortstack(nrec,dict,dbrch,nhldvoc,holdvoc)
  USE sortstring, ONLY: sort_string
  USE searching, ONLY: srch
  IMPLICIT NONE

  INTEGER,INTENT(IN)          :: nrec     ! # of species recorded in dict      
  CHARACTER(LEN=*),INTENT(IN) :: dict(:)  ! dictionary line (code+formula+fg)
  REAL,INTENT(IN)             :: dbrch(:) ! yield attached to a formula in dict
  INTEGER,INTENT(IN)          :: nhldvoc  ! # of (stable) VOC in the stack
  CHARACTER(LEN=*),INTENT(INOUT) :: holdvoc(:)  ! VOCs in the stack (name[a6]+formula[a120]+stabl[i3]+level[i3])

  CHARACTER(LEN=20+LEN(holdvoc(1))), ALLOCATABLE :: tmpstack(:)
  INTEGER :: ios,i,j, ind
  REAL    :: yld1, yld2

! create a temporary copy of the stack (including species yield)
  ALLOCATE(tmpstack(nhldvoc),STAT=ios)
  IF (ios/=0) THEN
    WRITE(6,*) '-error- in allocating tmpstack in sortstack'
    STOP "in sortstack"
  ENDIF

! create the string to be sorted (yield+holdvoc) and store in tmpstack 
  DO i=1,nhldvoc
    ind=srch(nrec,holdvoc(i)(7:126),dict)
    IF (ind<1) THEN
      WRITE(6,*) '-error- in sortstack, species in stack not found in dict'
      STOP "in sortstack"
    ENDIF
    WRITE(tmpstack(i),'(f20.10,a132)') dbrch(ind), holdvoc(i)
  ENDDO
  
! sort tmpstack according to yield(1st) and formula (next)
  CALL sort_string(tmpstack(1:nhldvoc)) 
  
! restore the stack holdvoc (from greater to smaller yield) 
  j=0 ; yld1=100.
  DO i=nhldvoc,1,-1
    j=j+1
    READ(tmpstack(i),'(f20.10,a132)') yld2, holdvoc(j)  
    IF (yld2>yld1) THEN
      PRINT*, 'yld1=',yld1 ; PRINT*, 'yld2=',yld2
      WRITE(6,*) '-error-, unexpected yld2>yld1 while sorting stack'
      STOP "in sortstack"
    ENDIF
    yld1=yld2
  ENDDO
  DEALLOCATE(tmpstack)

END SUBROUTINE sortstack

!=======================================================================
! PURPOSE: return the number of chromophore found in chem
!=======================================================================
INTEGER FUNCTION nphotogrp(chem)
  IMPLICIT NONE
  CHARACTER(LEN=*),INTENT(IN) :: chem  ! formula
  
  nphotogrp=0
  IF (INDEX(chem,'COC')/=0)      nphotogrp = nphotogrp+1 
  IF (INDEX(chem,'CO(OOH)')/=0)  nphotogrp = nphotogrp+1 
  IF (INDEX(chem,'CO(ONO2)')/=0) nphotogrp = nphotogrp+1
  IF ( (INDEX(chem,'CO(')/=0) .AND. (INDEX(chem,'CO(O')==0) ) nphotogrp=nphotogrp+1
  IF (INDEX(chem,'CdO')/=0)      nphotogrp = nphotogrp+1 
  IF (INDEX(chem,'CHO')/=0)      nphotogrp = nphotogrp+1 
  IF (INDEX(chem,'(OOH)')/=0)    nphotogrp = nphotogrp+1 
  IF (INDEX(chem,'(ONO2)')/=0)   nphotogrp = nphotogrp+1 
  IF (INDEX(chem,'CO(OONO2)')/=0)nphotogrp = nphotogrp+1 
  IF (INDEX(chem,'CO(OH)')/=0)   nphotogrp = nphotogrp+1 
END FUNCTION nphotogrp

!=======================================================================
! PURPOSE: Write (or not) the species under progress to the screen. Idea
! is to avoid unecessary screen writing to save running time. Species
! are written/erased in a temporary file during mechanism generation. 
!=======================================================================
SUBROUTINE wrtscreen(nrec,nhldvoc,stabl,code,chem)
  USE keyparameter, ONLY: scru
  USE keyflag, ONLY: screenfg
  IMPLICIT NONE
  
  INTEGER,INTENT(IN) :: nrec     ! # of species recorded in dict      
  INTEGER,INTENT(IN) :: nhldvoc  ! # of (stable) VOC in the stack
  INTEGER,INTENT(IN) :: stabl    ! # of generations needed to produce current species
  CHARACTER(LEN=*),INTENT(IN) :: code  ! formula
  CHARACTER(LEN=*),INTENT(IN) :: chem  ! formula

  INTEGER,SAVE :: irewind=0 ! static variable, initialised at 1st call only
  INTEGER,PARAMETER :: nrewind=100 ! number of species before rewind

! write the current "species line" in the file
  WRITE(scru,'(i8,1x,i5,1x,i3,2(1x,a))') nrec,nhldvoc,stabl,TRIM(code),TRIM(chem)

  IF (screenfg) THEN
    WRITE(6,*) '..........................................'
    WRITE(6,'(i8,1x,i6,1x,i3,1x,2(1x,a))') nrec,nhldvoc,stabl,TRIM(code),TRIM(chem)
  ELSE
    IF (nhldvoc==1 .AND. stabl==0) THEN
      WRITE(6,*) '..........................................'
      WRITE(6,'(i8,1x,i6,1x,i3,2(1x,a))') nrec,nhldvoc,stabl,TRIM(code),TRIM(chem)
      irewind=0
    ENDIF
  ENDIF

! every 'nrewind'th species
  irewind=irewind+1
  IF (irewind==nrewind) THEN
    irewind=0  ;  REWIND(scru)
    WRITE(scru,'(i8,1x,i5,1x,i3,2(1x,a))') nrec,nhldvoc,stabl,TRIM(code),TRIM(chem) 
    IF (.NOT.screenfg) WRITE(6,'(i8,1x,i6,1x,i3,2(1x,a))') nrec,nhldvoc,stabl,TRIM(code),TRIM(chem)
  ENDIF
END SUBROUTINE wrtscreen

!=======================================================================
! PURPOSE: Return the vapor pressure of the species provided as input.
!=======================================================================
REAL FUNCTION vaporpressure(chem,bond,group,nring,rjg,sarfg)
  USE simpoltool, ONLY: simpol          ! for simpol SAR
  USE nannoolaltool, ONLY:nannoolalprop ! for nannoonlal SAR
  USE myrdaltool, ONLY: myrdalprop      ! for Myrdal & Yalkowsky SAR
  USE atomtool, ONLY: molweight         ! needed in MY SAR
  USE keyflag, ONLY: pvap_sar             ! SAR to be used
  IMPLICIT NONE
  
  CHARACTER(LEN=*),INTENT(IN) :: chem
  INTEGER,INTENT(IN) :: bond(:,:)
  CHARACTER(LEN=*),INTENT(IN) :: group(:)
  INTEGER,INTENT(IN) :: nring
  INTEGER,INTENT(IN) :: rjg(:,:)
  INTEGER,INTENT(IN),OPTIONAL :: sarfg

  REAL :: latentheat, Tb, mw
  INTEGER :: ifg

! check if SAR is provided
  IF (PRESENT(sarfg)) THEN; ifg=sarfg ; ELSE ; ifg=pvap_sar ; ENDIF

  IF (ifg==1) THEN
    CALL molweight(chem,mw)
    CALL myrdalprop(chem,bond,group,nring,rjg,mw,Tb,vaporpressure,latentheat)
  ELSE IF (ifg==2) THEN
    CALL nannoolalprop(chem,bond,group,nring,rjg,Tb,vaporpressure,latentheat)
  ELSE IF (ifg==3) THEN
    CALL simpol(chem,bond,group,nring,vaporpressure,latentheat)
  ELSE
    WRITE(6,*) 'no method selected in vaporpressure function. Check pvap_sar'
    STOP "in vaporpressure"
  ENDIF
END FUNCTION vaporpressure

!=======================================================================
! PURPOSE: check if the chemistry of the species provided as input must
! be treated of "flushed". If flushed, then reaction still occurs (to
! remove the species) but no products are made. When flushing is wanted,
! the "iflost" flag (in module tempflag) is raised.
!=======================================================================
SUBROUTINE InOrOut(chem,parent,brch,stabl)
  USE keyflag, ONLY: maxgen,rxloss
  USE atomtool, ONLY: cnum, getatoms
  USE tempflag, ONLY: iflost,xxc,xxh,xxn,xxo,xxr,xxs,xxfl,xxbr,xxcl

  CHARACTER(LEN=*),INTENT(IN) :: chem
  CHARACTER(LEN=*),INTENT(IN) :: parent
  REAL, INTENT(IN) :: brch
  INTEGER,INTENT(IN) :: stabl

  CHARACTER(LEN=LEN(chem)) :: tchem ! copy of chem, without '#'
  INTEGER :: ih, in, io, ir, is, ifl, ibr, icl
  INTEGER :: tc
  INTEGER :: loout
  LOGICAL,PARAMETER :: losize=.FALSE.  ! account for parent size or not

  iflost=0 ! default value (make full chemistry)
  loout=0

! flush out species based on max # of generation or max yield
  IF (stabl>=maxgen) loout=1
  IF (brch<rxloss) loout=1

! flush out species based on max # of generation of a given parent size
  IF (losize) THEN
    tchem=parent
    IF (parent(1:3)=='#mm')   THEN ; tchem(1:)=parent(4:)
    ELSEIF (parent(1:1)=='#') THEN ; tchem(1:)=parent(2:)
    ENDIF
    CALL getatoms(tchem,tc,ih,in,io,ir,is,ifl,ibr,icl)
    
    IF      ( (tc>15).AND.(stabl>0) ) THEN ; loout=1
    ELSE IF ( (tc>9) .AND.(stabl>1) ) THEN ; loout=1
    ELSE IF ( (tc>7) .AND.(stabl>2) ) THEN ; loout=1
    ELSE IF ( (tc>5) .AND.(stabl>4) ) THEN ; loout=1
    ENDIF
  ENDIF
  
! count # of carbon atom to be removed (store in tempflag module)
  IF (loout==1) THEN
    iflost=1  ;  tchem=chem
    IF (chem(1:3)=='#mm' )    THEN ; tchem(1:)=chem(4:)
    ELSE IF (chem(1:1)=='#' ) THEN ; tchem(1:)=chem(2:)
    ENDIF
    CALL getatoms(tchem,xxc,xxh,xxn,xxo,xxr,xxs,xxfl,xxbr,xxcl)
  ENDIF

END SUBROUTINE InOrOut

!=======================================================================
! PURPOSE: set cutoff for the various reaction types.
!=======================================================================
SUBROUTINE cut_path(brch,cut_default,cut_OH,cut_O3,cut_NO3,cut_PAN,&
                    cut_HV, cut_RO,cut_RO2,cut_RCOO2) 
  USE keyflag, ONLY: brcut
  IMPLICIT NONE

  REAL,INTENT(IN)  :: brch
  REAL,INTENT(OUT) :: cut_default
  REAL,INTENT(OUT) :: cut_OH 
  REAL,INTENT(OUT) :: cut_O3 
  REAL,INTENT(OUT) :: cut_NO3 
  REAL,INTENT(OUT) :: cut_PAN 
  REAL,INTENT(OUT) :: cut_HV 
  REAL,INTENT(OUT) :: cut_RO 
  REAL,INTENT(OUT) :: cut_RO2 
  REAL,INTENT(OUT) :: cut_RCOO2 
  
  REAL :: cut

!  IF      (brch > 0.5)   THEN ; cut=0.01
!  ELSE IF (brch > 0.1)   THEN ; cut=0.05
!  ELSE IF (brch > 0.01)  THEN ; cut=0.10
!  ELSE IF (brch > 0.001) THEN ; cut=0.20
!  ELSE                        ; cut=0.30
!  ENDIF

!  IF      (brch > 0.5)   THEN ; cut=0.05
!  ELSE IF (brch > 0.1)   THEN ; cut=0.05
!  ELSE IF (brch > 0.01)  THEN ; cut=0.05
!  ELSE IF (brch > 0.001) THEN ; cut=0.05
!  ELSE                        ; cut=0.05
!  ENDIF
  
!  cut=0.05
  cut=brcut
  
  cut_default=cut
  cut_OH=cut
  cut_O3=cut
  cut_NO3=cut
  cut_PAN=cut
  cut_HV=cut
  cut_RO=cut
  cut_RO2=cut
  cut_RCOO2=cut
 
END SUBROUTINE cut_path

!=======================================================================
! PURPOSE: check that functions between parenthesis is allowed.
!=======================================================================
SUBROUTINE check_parenthesis(chem)
  USE keyparameter, ONLY: mxnode,mxlgr
  USE stdgrbond, ONLY: grbond
  USE toolbox, ONLY: countstring,stoperr
  IMPLICIT NONE

  CHARACTER(LEN=*),INTENT(IN) :: chem

  CHARACTER(LEN=mxlgr)        :: group(mxnode),funct
  INTEGER                     :: bond(mxnode,mxnode)
  INTEGER                     :: i,j,i1,i2,ngr,dbflg,nring,nf
  LOGICAL                     :: locheck

  CHARACTER(LEN=19),PARAMETER :: prog='*check_parenthesis*'
  CHARACTER(LEN=70)           :: mesg
  CHARACTER(LEN=LEN(chem))    :: tchem

  IF (chem(1:1)=='#') THEN 
    tchem=chem(2:)
  ELSE
    tchem=chem
  ENDIF
  
  CALL grbond(tchem,group,bond,dbflg,nring)
  ngr=COUNT(group/=' ')

  DO i=1,ngr
    nf=countstring(group(i),'(') ; IF (nf==0) CYCLE

    funloop : DO j=1,nf
      locheck=.FALSE. ; funct=' ' 
      i1=INDEX(group(i),'(') ; i2=INDEX(group(i),')')
      funct=group(i)(i1+1:i2-1)
      group(i)(i1:i1)=' ' ; group(i)(i2:i2)=' '

      IF (TRIM(funct)=='OH')    locheck=.TRUE.
      IF (TRIM(funct)=='OOH')   locheck=.TRUE.
      IF (TRIM(funct)=='OOOH')  locheck=.TRUE.
      IF (TRIM(funct)=='NO2')   locheck=.TRUE.
      IF (TRIM(funct)=='ONO2')  locheck=.TRUE.
      IF (TRIM(funct)=='OONO2') locheck=.TRUE.

      IF (TRIM(funct)=='OO.')   locheck=.TRUE.
      IF (TRIM(funct)=='ZOO.')  locheck=.TRUE.
      IF (TRIM(funct)=='EOO.')  locheck=.TRUE.
      IF (TRIM(funct)=='O.')    locheck=.TRUE.
      
      IF (locheck) CYCLE funloop
      
      mesg= '/!\ Function not allowed or mispelled'
      CALL stoperr(prog,mesg,chem)
    ENDDO funloop
  ENDDO
  RETURN

END SUBROUTINE check_parenthesis

END MODULE maintool
  
