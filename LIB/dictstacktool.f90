MODULE dictstacktool
USE toolbox, ONLY: stoperr,addref
IMPLICIT NONE
CONTAINS
! SUBROUTINE bratio(pchem,brtio,pname)
! SUBROUTINE add_topvocstack(inchem,brtio,pname,ncom,com)
! SUBROUTINE loader(chem,name)
! SUBROUTINE scrollc1stack(c1stack)
! SUBROUTINE addc1chem(pchem,pname,nstck,c1stack,ierr)
! SUBROUTINE addc1dict(pchem,pname,fgrp)

!=======================================================================
! PURPOSE : gives the short name of the species (pchem) given as input 
! and updates the stack and dictionary arrays if necessary.
!  If the species is already known, the routine just returns its short 
!  name. If the species is new, then:                                   
!    1- a short name is given to the species (see naming)          
!    2- the "dictionary" arrays are updated (dict, namlst, dbrch).            
!       Species is added in such a way that the tables remain sorted   
!    3- the species is added to stack for future reactions (see loader)
!=======================================================================
SUBROUTINE bratio(inchem,brtio,pname,ncom,com)
  USE keyparameter, ONLY: mxlfl,mxlco,mxlfo  ! length of funct. list, name, formula
  USE references, ONLY:mxlcod
  USE keyflag, ONLY: enolflg,isomerfg
  USE cdtool, ONLY: switchenol
  USE atomtool, ONLY: cnum, onum
  USE searching, ONLY: srch, srh5
  USE database, ONLY: nspsp, dictsp     ! special species input
  USE dictstackdb, ONLY: nrec,ninorg,inorglst,dict,namlst,dbrch,diccri
  USE namingtool, ONLY: naming
  USE switchisomer, ONLY: isomer  
  USE tempflag, ONLY: iflost
  IMPLICIT NONE

  REAL,INTENT(IN)             :: brtio   ! yield of the species
  CHARACTER(LEN=*),INTENT(IN) :: inchem  ! formula for which a name must be given 
  CHARACTER(LEN=*),INTENT(OUT):: pname   ! name corresponding to formula
  INTEGER,INTENT(INOUT) :: ncom          ! # of ref/com in the current list (tpref)
  CHARACTER(LEN=*), INTENT(INOUT) :: com(:)  ! list of references/comments

  CHARACTER(LEN=mxlfl) :: fgrp
  CHARACTER(LEN=LEN(inchem)) :: pchem
  INTEGER :: dicptr,namptr,nca,chg,i,ipos,ierr
  INTEGER :: tabinfo(SIZE(diccri,2))
  LOGICAL :: loswitch                    ! switch for keto/enol change
  
  INTEGER, PARAMETER :: mx1c=15                  ! size of the C1 stack
  CHARACTER(LEN=LEN(inchem)) :: c1stack(mx1c)    ! C1 stack 
  INTEGER :: nstck                               ! # of element in C1 stack
  

  CHARACTER(LEN=7)  :: progname='bratio'
  CHARACTER(LEN=70) :: mesg

  IF (inchem==' ') RETURN   ! return if no species 

! check keyword (can be given by special mechanisms)
  IF      (inchem(1:5)=='EXTRA')  THEN ; pname='EXTRA ' ; RETURN
  ELSE IF (inchem(1:2)=='HV')     THEN ; pname='HV    ' ; RETURN
  ELSE IF (inchem(1:6)=='OXYGEN') THEN ; pname='OXYGEN' ; RETURN
  ELSE IF (inchem(1:6)=='ISOM  ') THEN ; pname='ISOM  ' ; RETURN
  ELSE IF (inchem(1:4)=='(+M)')   THEN ; pname='(+M)  ' ; RETURN !! "FALLOFF"=7
  ENDIF
  
  
  pname=' '
  pchem=inchem                          ! make a working copy of inchem
  nca=cnum(pchem)+onum(pchem)
  tabinfo(:)=0

! special name (inorganics, formulae that cannot be held, and C1)
! =============================================================

! no carbon - check if known in the list of inorganics or keywords. 
  IF (nca==0) THEN
    DO i=1,ninorg
      IF (INDEX(inorglst(i)(10:mxlfo+11),pchem(1:LEN_TRIM(pchem)+1)).EQ.1) THEN
        pname = inorglst(i)(1:mxlco)
        RETURN
      ENDIF
    ENDDO
  ENDIF

!catcher for "COO"
  IF (inchem(1:4)=='COO ') pchem(1:4)='CO2 '

! Check if "special species" (formula starts with '#')
! =============================================================
  IF (pchem(1:1)=='#') THEN  
    ipos=0
    dicptr = srch(nrec,pchem,dict)

! if already recorded ...
    IF (dicptr>0) THEN
      pname = dict(dicptr)(1:mxlco)
      dbrch(dicptr) = MAX(dbrch(dicptr),brtio)
      RETURN
    ENDIF

! If not recorded, then search in the special list 
    sploop: DO i=1,nspsp
      IF (INDEX(dictsp(i)(10:mxlfo+11),pchem(1:INDEX(pchem," "))).EQ.1) THEN
        ipos=i
        pname=dictsp(i)(1:mxlco)
        namptr=srh5(pname,namlst,nrec)
        namptr=-namptr
! if species is already in the namlist, we don't need to add it.
        IF (namptr<0) RETURN
        fgrp=dictsp(i)(mxlfo+12:)
        EXIT sploop
      ENDIF
    ENDDO sploop

    IF (ipos==0) THEN
      mesg="Special species cannot be managed. Check formula."
      CALL stoperr(progname,mesg,pchem)
    ENDIF     

! C1 species, must already be known !
! =============================================================
  ELSE IF (nca==1) THEN
    dicptr = srch(nrec,pchem,dict)
    IF (dicptr>0) THEN                          ! species already included
      pname = dict(dicptr)(1:mxlco)
      RETURN

! not yet known ...
    ELSE IF (iflost==0) THEN                                       
      
! check if species & chemistry is "available" in C1 routines (add all then & return pname) 
      ierr=1 ; nstck=0 ; c1stack(:)=' '
      CALL addc1chem(pchem,pname,nstck,c1stack,ierr)

! species unknown - stop
      IF (ierr/=0) THEN                         
        mesg="Following C1 species not in the dictionary"
        CALL stoperr(progname,mesg,pchem)
      ENDIF

! If additional species was produced by pchem, then add also those species      
      IF (nstck/=0) CALL scrollc1stack(c1stack) 

      RETURN                                    ! all clear - return with pname
    ENDIF

! Regular species 
! =============================================================
  ELSE 

! Regular formula and must start with a "C,c, or O"  
    IF ((pchem(1:1)/='C').AND.(pchem(1:2)/='-O').AND.(pchem(1:1)/='c')) THEN
      mesg="Species cannot be managed. Check formula."
      CALL stoperr(progname,mesg,pchem)
    ENDIF

! check enols and switch to keto form
   IF (enolflg) THEN
     IF (INDEX(pchem,'Cd')/=0) THEN
       CALL switchenol(pchem,loswitch)  ! pchem is returned as std keto if enol
       IF (loswitch) CALL addref(progname,'KETOENOL',ncom,com,inchem)
     ENDIF
   ENDIF

! Search if pchem is already recorded. If yes, return the short name
    dicptr = srch(nrec,pchem,dict)
    IF (dicptr>0) THEN
      pname = dict(dicptr)(1:mxlco)
      dbrch(dicptr) = MAX(dbrch(dicptr),brtio)
      RETURN
    ENDIF

! the species is unknown 
! --------------------------------------

! LUMP SECONDARY SPECIES - not available. See old version of geckoa 
! including subroutine "lump_sec" and the corresponding bratio version. 

! ISOMER SWITCH - search if an isomer is already known and replace formula 
    IF ((isomerfg).AND.(nca>3).AND.(INDEX(pchem,'.')==0)) THEN
      CALL isomer(pchem,brtio,chg,tabinfo)
      IF (chg==1) THEN     ! pchem was switched - get corresponding name
        dicptr = srch(nrec,pchem,dict)
        IF (dicptr<=0)  THEN
          mesg="expected species after isomer switch not found"
          CALL stoperr(progname,mesg,pchem)
        ENDIF
        dbrch(dicptr) = MAX(dbrch(dicptr),brtio)
        pname = dict(dicptr)(1:mxlco)
        CALL addref(progname,'SWAPISOMER',ncom,com,inchem)
        RETURN 
      ENDIF
    ENDIF 

! Get name (pname) for pchem and position (namptr) in namlst for addition
    CALL naming(pchem,namptr,pname,fgrp)

  ENDIF

! update stack and dictionary array (species may come be a '#' species)
! ====================================================================

! if the flag to stop the chemistry is raised, then return
  IF (iflost==1) RETURN

! raise the counters
  nrec=nrec+1  
  IF (nrec > SIZE(dict)) THEN
    mesg="Too many species added to dictionary"
    CALL stoperr(progname,mesg,inchem)
  ENDIF 

! add new name and raise the name array
  namptr=namptr+1
  namlst(namptr+1:nrec+1)=namlst(namptr:nrec)
  namlst(namptr)=pname

! raise dictionary and branching table. Add new line for the new species 
  dicptr = ABS(dicptr) + 1
  dict(dicptr+1:nrec+1)=dict(dicptr:nrec)
  dbrch(dicptr+1:nrec+1)=dbrch(dicptr:nrec)
  WRITE(dict(dicptr),'(a6,3x,a120,2x,a15)') pname, pchem, fgrp
  dbrch(dicptr)  = brtio

! store information required to search for an isomer (info saved in tabinfo)
  IF (isomerfg) THEN
    diccri(dicptr+1:nrec+1,:) = diccri(dicptr:nrec,:)
    diccri(dicptr,:)=tabinfo(:)
  ENDIF

  CALL loader(pchem,pname)
    
END SUBROUTINE bratio

!=======================================================================
! PURPOSE: raise flag for an addition on top of the VOC stack before
! calling bratio. If the species is already known in the mechanism, then
! nothing should happen: the routine simply return the short name 
! (pname) corresponding to the input species (inchem). If the species 
! is new, then it will be added on top of the VOC stack, with a # of 
! generation identical to the current value.
!=======================================================================
SUBROUTINE add_topvocstack(inchem,brtio,pname,ncom,com)
  USE dictstackdb, ONLY: stabl,lotopstack
  IMPLICIT NONE

  REAL,INTENT(IN)             :: brtio   ! yield of the species
  CHARACTER(LEN=*),INTENT(IN) :: inchem  ! formula for which a name must be given 
  CHARACTER(LEN=*),INTENT(OUT):: pname   ! name corresponding to formula
  INTEGER,INTENT(INOUT) :: ncom          ! # of ref/com in the current list (tpref)
  CHARACTER(LEN=*), INTENT(INOUT) :: com(:)  ! list of references/comments
  
  INTEGER :: savstabl
  
  savstabl=stabl 

! change data in dictstackdb
  stabl=stabl-1          ! because loader will add 1 to stabl once called
  lotopstack=.TRUE.      ! raise flag for top addition
  
! return short name and add species to the stack (if new)
  CALL bratio(inchem,brtio,pname,ncom,com)

! restore data in dictstackdb 
  stabl=savstabl      ! curretnt value if stabl
  lotopstack=.FALSE.  ! default value for stack management

END SUBROUTINE add_topvocstack

!=======================================================================
! PURPOSE: Load molecules in the stack for further reactions. 
! Two stacks are considered, one for the non radical VOC and one for 
! the radicals. Reactions for radical are written first (see main.f).
! New species are added to the bottom of the stack (first in first out)
! as the default situation.
!                                                                     
! Information loaded in the stack comprises (holdvoc & holdrad):
!   holdvoc(1:6) = short name of the chemical                       
!   holdvoc(7:126) = chemical formula                               
!   holdvoc(127:129) = # of stable generations
!   holdvoc(130:132) = # of levels (including radicals)
!                                                                     
! INPUT from module dictstackdb:
!  - level      : # of level (stable+radical) needed to produce "chem"
!  - stabl      : # of generation (stable) needed to produce "chem"
!  - lotopstack : if true, then VOC is added on top on the stack (lifo)  
! INPUT/OUTPUT from module dictstackdb:
!  - nhldvoc    : number of (stable) VOC in the stack 
!  - holdvoc(i) : list of the VOC in the stack
!  - nhldrad    : number of radical in the stack
!  - holdrad(i) : list of the radicals in the stack
!=======================================================================
SUBROUTINE loader(chem,idnam)
  USE dictstackdb, ONLY: nhldvoc,holdvoc,nhldrad,holdrad,stabl, &
                         level,lotopstack
  IMPLICIT NONE

  CHARACTER(LEN=*),INTENT(IN) :: chem  ! formula to be added to the stack
  CHARACTER(LEN=*),INTENT(IN) :: idnam ! name to be added to the stack

  INTEGER           :: tplevel, tpstabl
  CHARACTER(LEN=7)  :: progname='loader'
  CHARACTER(LEN=70) :: mesg

! raise the counters:
  tplevel=level+1 
  IF (INDEX(chem,'.') == 0) THEN ; tpstabl=stabl+1
  ELSE                           ; tpstabl=stabl
  ENDIF

! VOC (non radical) species
! --------------------------
  IF (INDEX(chem,'.')==0) THEN     ! usual VOC
    nhldvoc = nhldvoc + 1
    IF (nhldvoc > SIZE(holdvoc)) THEN
      mesg="Stack full. # of VOC species > VOC stack size"
      CALL stoperr(progname,mesg,chem)
    ENDIF
    IF (lotopstack) THEN   ! special case - add species on top of the stack (lifo)
      holdvoc(2:nhldvoc)=holdvoc(1:nhldvoc-1)     
      WRITE(holdvoc(1),'(a6,a120,i3,i3)') idnam,chem,tpstabl,tplevel
    ELSE                   ! default (fifo)
      WRITE(holdvoc(nhldvoc),'(a6,a120,i3,i3)') idnam,chem,tpstabl,tplevel 
    ENDIF
      
! radical species
! ---------------
  ELSE                                    ! usual radical
    nhldrad = nhldrad + 1
    IF (nhldrad > SIZE(holdrad)) THEN
      mesg="Stack full. # of radical species > radical stack size"
      CALL stoperr(progname,mesg,chem)
    ENDIF
    WRITE(holdrad(nhldrad),'(a6,a120,i3,i3)') idnam,chem,tpstabl,tplevel
  ENDIF

END SUBROUTINE loader

! ======================================================================
! PURPOSE: Write the species in the "C1stack" to the dictionaries and 
! the corresponding chemistry to the mechanism. That chemistry may lead
! to the addition of more C1 species to the stack. The routine scrolls
! all species in the stack and stops when the stack is empty.
!
! For future development: all this C1 chemistry/species is "hard coded".
! Some future development might be useful to make it "soft" (reading
! input files)
! ======================================================================
SUBROUTINE scrollc1stack(c1stack)
  USE keyparameter, ONLY: mxlco
  USE dictstackdb, ONLY: nrec,dict
  USE searching, ONLY: srch
  IMPLICIT NONE
  CHARACTER(LEN=*),INTENT(OUT) :: c1stack(:)  ! new C1 species not considered in mechanism 
  
  CHARACTER(LEN=LEN(c1stack(1))) :: chem
  INTEGER :: icount,nstck,ierr,dicptr
  CHARACTER(LEN=14),PARAMETER :: progname='scrollc1stack'
  CHARACTER(LEN=70) :: mesg
  
  CHARACTER(LEN=mxlco) :: dumnam

  nstck=COUNT(c1stack/=' ')
  icount=0

! loop over C1 species in stack
  DO                                        ! just in case, count the ...   
    icount=icount+1                         ! ... number of iteration to ...
    IF (icount==100) THEN                   ! ... infinite loop.
      mesg="infinite loop when adding C1 species"
      CALL stoperr(progname,mesg,c1stack(1))
    ENDIF

! get next species from stack (exit if empty)
    chem=c1stack(1)
    !PRINT*, 'next from stack:',chem 
    IF (chem==' ') EXIT                     ! no more species in stack, exit

! check if the species already added to dictionary 
    dicptr = srch(nrec,chem,dict)
    !PRINT*, 'dictpr=',dicptr
    IF  (dicptr > 0) THEN 
      IF (nstck==1) EXIT                    ! exit - nothing left in stack 
      c1stack(1:nstck-1)=c1stack(2:nstck)   ! move -1 stack list
      nstck=nstck-1

! add new species
    ELSE                                      
      ierr=1
      CALL addc1chem(chem,dumnam,nstck,c1stack,ierr)
      IF (ierr/=0) THEN                                    ! species not available
        mesg="Following C1 species not in the dictionary"
        CALL stoperr(progname,mesg,c1stack(1))
      ENDIF

      IF (nstck==1) EXIT                    ! exit - nothing left in stack 
      c1stack(1:nstck-1)=c1stack(2:nstck)   ! move -1 stack list
      nstck=nstck-1
    ENDIF
  ENDDO

END SUBROUTINE scrollc1stack

! ======================================================================
! PURPOSE: check if the formula provided as input (pchem) is "known" as
! a "regular" species. If so, then the species is added in the
! dictionaries and the chemistry (if any) is added in the mechanism.
! The subroutine returns the "regular" short name (pname) for the input
! species. 
! New C1 species (provided from the chemistry of pchem) that may not
! already be considered are added to a C1 stack for further addition in
! dictionaries and mechanisms.
! NB: if new C1 species are added here, they MUST ALSO be added to
!     array "optional_c1" in spsptool.f90
! For future development: all this C1 chemistry/species is "hard coded".
! Some further development might be useful to make it "soft" (reading
! input files).
! The list of allowed formulae and corresponding names is:
! CH2OO  : CH2.(OO.)
! CH3O3H : CH3(OOOH)
! COHOOH : CH2(OH)(OOH)
! COHO3H : CH2(OH)(OOOH)
! CHONO2 : CHO(NO2)
! C10001 : CH2(ONO2)(OO.)
! C10002 : CO(NO2)(ONO2)
! C10003 : CO(OH)(NO2)
! C10004 : CO(OH)(ONO2)
! C10005 : CO(OH)(OOH)
! C10006 : CO(ONO2)(ONO2)
! C10007 : CO(ONO2)(OOH)
! C10008 : CO(OOH)(OO.)
! C10009 : CO(NO2)(OOOH)
! C10010 : CO(OOOH)(OOOH)
! H2CO3  : CO(OH)(OH)
! HCOO3H : CHO(OOOH)                   
! NO1001 : CH2(OH)(ONO2)
! NH1001 : CH2(ONO2)(OOH)
! 2V1001 : CH2(NO2)(OO.)
! C1ND00 : CHO(ONO2)
! N01003 : CH3(ONO)        ! still in dic_c1.dat and mch_singlec.dat
!
! to be continued ...
! ======================================================================
SUBROUTINE addc1chem(pchem,pname,nstck,c1stack,ierr)
  USE keyparameter, ONLY: mxpd,mecu,mxlfl,mxlco
  USE normchem, ONLY: stdchm
  USE references, ONLY: mxlcod
  USE dictstackdb, ONLY: nrec,namlst
  USE toolbox, ONLY: add1tonp
  USE searching, ONLY: srh5
  USE rxwrttool, ONLY:rxwrit,rxinit
  IMPLICIT NONE

  CHARACTER(LEN=*),INTENT(IN) :: pchem    ! species and chemistry to be added 
  CHARACTER(LEN=*),INTENT(OUT) :: pname   ! name of chem to be returned (if found)
  CHARACTER(LEN=*),INTENT(INOUT) :: c1stack(:)  ! stack of C1 species to be added to mechanism
  INTEGER,INTENT(INOUT) :: nstck          ! # of species in C1 stack  
  INTEGER,INTENT(OUT)   :: ierr           ! error flag 

  CHARACTER(LEN=mxlfl) :: fgrp
  CHARACTER(LEN=LEN(pchem)) :: tchem
  INTEGER :: np,namptr
  CHARACTER(LEN=LEN(pname)) :: p(mxpd), r(3)
  REAL    :: s(mxpd), arrh(3)
  INTEGER :: idreac, nlabel
  REAL    :: xlabel,folow(3),fotroe(4)
  REAL    :: brtio
  CHARACTER(LEN=mxlco)  :: namacetic=' '   ! intialize at 1st call only

  INTEGER,PARAMETER     :: mxcom=10
  INTEGER               :: nref
  CHARACTER(LEN=mxlcod) :: ref(mxcom)
  CHARACTER(LEN=10),PARAMETER :: progname='addc1chem'
  CHARACTER(LEN=70) :: mesg

! rate constants for CH2.(OO). criegee radical
  REAL,DIMENSION(3),PARAMETER :: kuni_ch2oo  =(/ 1.66E+01, 4.02,  8024. /)
  REAL,DIMENSION(3),PARAMETER :: kmh2o_ch2oo =(/ 8.13E-18, 1.23,   698. /)
  REAL,DIMENSION(3),PARAMETER :: kdh2o_ch2oo =(/ 7.95E-18, 1.24, -1510. /)
  REAL,DIMENSION(3),PARAMETER :: kso2_ch2oo  =(/ 3.70E-11, 0.  ,     0. /)
  REAL,DIMENSION(3),PARAMETER :: kno2_ch2oo  =(/ 3.00E-12, 0.  ,     0. /)
  REAL,DIMENSION(3),PARAMETER :: khno3_ch2oo =(/ 5.40E-10, 0.  ,     0. /)
  REAL,DIMENSION(3),PARAMETER :: krcooh_ch2oo=(/ 1.20E-10, 0.  ,     0. /)
  
  pname=' '  ;  ierr=1   ! raise flag as default 

! ----------
! CH2.(OO).
! ----------
  IF (pchem=='CH2.(OO.)') THEN

!-- add the species to the dictionary
    pname='CH2OO'  ;   fgrp='4 '
    CALL addc1dict(pchem,pname,fgrp)
    ierr=0               ! flag down, species found

!-- add chemistry

    ! unimolecular decomposition (CO2+H2, CO+H2O, 50% each)
    np=0
    CALL rxinit(r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe)
    nref=0 ; ref(:)=' ' ; CALL addref(progname,'MN22SAR2',nref,ref,pchem)

    CALL add1tonp(progname,pchem,np) ; s(np)=0.5 ; p(np)='CO2  '
    CALL add1tonp(progname,pchem,np) ; s(np)=0.5 ; p(np)='H2   '
    CALL add1tonp(progname,pchem,np) ; s(np)=0.5 ; p(np)='CO   ' 
    CALL add1tonp(progname,pchem,np) ; s(np)=0.5 ; p(np)='H2O  ' 

    r(1)=pname
    arrh(:) = kuni_ch2oo(:)
    CALL rxwrit(mecu,r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe,ref)

    ! reaction with water
    np=0
    CALL rxinit(r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe)
    nref=0 ; ref(:)=' ' ; CALL addref(progname,'MN22SAR2',nref,ref,pchem)

    CALL add1tonp(progname,pchem,np) ; s(np)=0.55 ; p(np)='COHOOH'
    CALL add1tonp(progname,pchem,np) ; s(np)=0.40 ; p(np)='CH2O  ' 
    CALL add1tonp(progname,pchem,np) ; s(np)=0.40 ; p(np)='H2O2  ' 
    CALL add1tonp(progname,pchem,np) ; s(np)=0.05 ; p(np)='HCOOH ' 

    r(1)=pname
    r(2)='EXTRA'
    idreac=2  ;   nlabel=500
    arrh(:) = kmh2o_ch2oo(:)
    CALL rxwrit(mecu,r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe,ref)

   ! reaction with water dimer
    np=0
    CALL rxinit(r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe)
    nref=0 ; ref(:)=' ' ; CALL addref(progname,'MN22SAR2',nref,ref,pchem)

    CALL add1tonp(progname,pchem,np) ; s(np)=0.55 ; p(np)='COHOOH'
    CALL add1tonp(progname,pchem,np) ; s(np)=0.40 ; p(np)='CH2O  ' 
    CALL add1tonp(progname,pchem,np) ; s(np)=0.40 ; p(np)='H2O2  ' 
    CALL add1tonp(progname,pchem,np) ; s(np)=0.05 ; p(np)='HCOOH ' 

    r(1)=pname
    r(2)='EXTRA'
    idreac=2  ;  nlabel=502
    arrh(:) = kdh2o_ch2oo(:)
    CALL rxwrit(mecu,r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe,ref)

    ! new C1 species that may need to be added    
    namptr=srh5('COHOOH',namlst,nrec)
    IF (namptr<0) THEN
      nstck=nstck+1 
      IF (nstck > SIZE(c1stack)) THEN
        mesg="too many species added to C1 stack"
        CALL stoperr(progname,mesg,pchem)
      ENDIF
      c1stack(nstck)='CH2(OH)(OOH) '
    ENDIF

    ! reaction with SO2
    np=0
    CALL rxinit(r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe)
    nref=0 ; ref(:)=' ' ; CALL addref(progname,'MN22SAR2',nref,ref,pchem)

    CALL add1tonp(progname,pchem,np) ; s(np)=1. ; p(np)='CH2O  ' 
    CALL add1tonp(progname,pchem,np) ; s(np)=1. ; p(np)='SULF  ' 

    r(1)=pname  ;  r(2)='SO2 '  ;  arrh(:) = kso2_ch2oo(:)
    CALL rxwrit(mecu,r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe,ref)

    ! reaction with NO2 (temporary version, product expected is CH2(NO2)(OO.)
    np=0
    CALL rxinit(r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe)
    nref=0 ; ref(:)=' ' ; CALL addref(progname,'MN22SAR2',nref,ref,pchem)

    CALL add1tonp(progname,pchem,np) ; s(np)=1. ; p(np)='2V1001' 

    r(1)=pname ;  r(2)='NO2 ' ;   arrh(:) = kno2_ch2oo(:)
    CALL rxwrit(mecu,r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe,ref)

    ! new C1 species that may need to be added    
    namptr=srh5('2V1001',namlst,nrec)
    IF (namptr<0) THEN
      nstck=nstck+1 
      IF (nstck > SIZE(c1stack)) THEN
        mesg="too many species added to C1 stack"
        CALL stoperr(progname,mesg,pchem)
      ENDIF
      c1stack(nstck)='CH2(NO2)(OO.) '
    ENDIF

! reaction with HNO3 (temporary version, product expected is CH2(ONO2)(OOH)
    np=0
    CALL rxinit(r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe)
    nref=0 ; ref(:)=' ' ; CALL addref(progname,'MN22SAR2',nref,ref,pchem)

    CALL add1tonp(progname,pchem,np) ; s(np)=1. ; p(np)='NH1001' 

    r(1)=pname ;  r(2)='HNO3 ' ;   arrh(:) = khno3_ch2oo(:)
    CALL rxwrit(mecu,r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe,ref)

    ! new C1 species that may need to be added    
    namptr=srh5('NH1001',namlst,nrec)
    IF (namptr<0) THEN
      nstck=nstck+1 
      IF (nstck > SIZE(c1stack)) THEN
        mesg="too many species added to C1 stack"
        CALL stoperr(progname,mesg,pchem)
      ENDIF
      c1stack(nstck)='CH2(ONO2)(OOH) '
    ENDIF

! reaction with HCOOH
    np=0
    CALL rxinit(r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe)
    nref=0 ; ref(:)=' ' ; CALL addref(progname,'MN22SAR2',nref,ref,pchem)

    CALL add1tonp(progname,pchem,np) ; s(np)=1.
    tchem='CHO-O-CH2(OOH)'  ;   CALL stdchm(tchem)
    brtio=1.  ; CALL bratio(tchem,brtio,p(1),nref,ref)

    r(1)=pname ;  r(2)='HCOOH ' ;   arrh(:) = krcooh_ch2oo(:)
    CALL rxwrit(mecu,r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe,ref)

! reaction with CH3COOH
    np=0
    CALL rxinit(r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe)
    nref=0 ; ref(:)=' ' ; CALL addref(progname,'MN22SAR2',nref,ref,pchem)

    CALL add1tonp(progname,pchem,np) ; s(np)=1.
    tchem='CH3CO-O-CH2(OOH)'  ;   CALL stdchm(tchem)
    brtio=1.  ; CALL bratio(tchem,brtio,p(1),nref,ref)

    r(1)=pname 
    IF (namacetic==' ') THEN
      tchem='CH3CO(OH) ' ; brtio=1.  ! CH3COOH assumed as a primary species
      CALL add_topvocstack(tchem,brtio,namacetic,nref,ref)
    ENDIF
    r(2)=namacetic
    arrh(:) = krcooh_ch2oo(:)
    CALL rxwrit(mecu,r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe,ref)

! ------------
! CH3(OOOH)
! ------------
  ELSE IF (pchem=='CH3(OOOH)') THEN
! these species are produced from the RO2 + OH reaction (see Jenkin et al., 2019)
! use rate for CH3OOOH + HO => CH3O. + H2O + O2 from Anglada and Solé (2018)
! use decompostion rate from Assaf et al., 2018 (80% path 1, 20% path 2)

!-- add the species to the dictionary
    pname='CH3O3H'  ;   fgrp='   '
    CALL addc1dict(pchem,pname,fgrp)
    ierr=0               ! flag down, species found

!-- add chemistry
! reaction with OH
    np=0
    CALL rxinit(r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe)
    nref=0 ; ref(:)=' ' ; CALL addref(progname,'JA18KMV000',nref,ref,pchem)

    CALL add1tonp(progname,pchem,np) ; s(np)=1.0 ; p(np)='CH3O  '

    r(1)=pname
    r(2)='HO '
    arrh(1) = 1.46E-12 ; arrh(2)=0. ; arrh(3)=-1037
    CALL rxwrit(mecu,r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe,ref)

! decomposition
    np=0
    CALL rxinit(r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe)
    nref=0 ; ref(:)=' ' ; CALL addref(progname,'EA18KMR000',nref,ref,pchem)

    CALL add1tonp(progname,pchem,np) ; s(np)=0.8 ; p(np)='CH3O  '
    CALL add1tonp(progname,pchem,np) ; s(np)=0.8 ; p(np)='HO2   '
    CALL add1tonp(progname,pchem,np) ; s(np)=0.2 ; p(np)='CH3OH '

    r(1)=pname
    arrh(1)=1.52E10 ; arrh(2)=1.35 ; arrh(3)=12000
    CALL rxwrit(mecu,r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe,ref)

! ------------
! CHO(OOOH)
! ------------
  ELSE IF (pchem=='CHO(OOOH)') THEN
! these species are produced from the RO2 + OH reaction (see Jenkin et al., 2019)
! use rate for CH3OOOH + HO => CH3O. + H2O + O2 from Anglada and Solé (2018)
! use decompostion rate from Assaf et al., 2018 (80% path 1, 20% path 2)

!-- add the species to the dictionary
    pname='HCOO3H'  ;   fgrp='   '
    CALL addc1dict(pchem,pname,fgrp)
    ierr=0               ! flag down, species found

!-- add chemistry
! reaction with OH
    np=0
    CALL rxinit(r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe)
    nref=0 ; ref(:)=' ' ; CALL addref(progname,'JA18KMV000',nref,ref,pchem)

    CALL add1tonp(progname,pchem,np) ; s(np)=1.0 ; p(np)='CO    '
    CALL add1tonp(progname,pchem,np) ; s(np)=1.0 ; p(np)='HO2   '

    r(1)=pname
    r(2)='HO '
    arrh(1) = 1.46E-12 ; arrh(2)=0. ; arrh(3)=-1037
    CALL rxwrit(mecu,r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe,ref)

! decomposition
    np=0
    CALL rxinit(r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe)
    nref=0 ; ref(:)=' ' ; CALL addref(progname,'EA18KMR000',nref,ref,pchem)

    CALL add1tonp(progname,pchem,np) ; s(np)=0.8 ; p(np)='CO    '
    CALL add1tonp(progname,pchem,np) ; s(np)=1.6 ; p(np)='HO2   '
    CALL add1tonp(progname,pchem,np) ; s(np)=0.2 ; p(np)='HCOOH '

    r(1)=pname
    arrh(1)=1.52E10 ; arrh(2)=1.35 ; arrh(3)=12000
    CALL rxwrit(mecu,r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe,ref)

! ------------
! CH2(OH)(OOH)
! ------------
  ELSE IF (pchem=='CH2(OH)(OOH)') THEN

!-- add the species to the dictionary
    pname='COHOOH'  ;   fgrp='HO '
    CALL addc1dict(pchem,pname,fgrp)
    ierr=0               ! flag down, species found

!-- add chemistry
! reaction with OH
    np=0
    CALL rxinit(r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe)
    nref=0 ; ref(:)=' ' ; CALL addref(progname,'HA18KMV000',nref,ref,pchem)

    CALL add1tonp(progname,pchem,np) ; s(np)=0.55 ; p(np)='CH2O  '
    CALL add1tonp(progname,pchem,np) ; s(np)=0.55 ; p(np)='HO2   '
    CALL add1tonp(progname,pchem,np) ; s(np)=0.45 ; p(np)='HCOOH '
    CALL add1tonp(progname,pchem,np) ; s(np)=0.45 ; p(np)='HO    '

    r(1)=pname
    r(2)='HO '
    arrh(1) = 7.1E-12 ; arrh(2)=0. ; arrh(3)=0.
    CALL rxwrit(mecu,r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe,ref)

! ------------
! CH2(OH)(OOOH)
! ------------
  ELSE IF (pchem=='CH2(OH)(OOOH)') THEN

!-- add the species to the dictionary
    pname='COHO3H'  ;   fgrp='O '
    CALL addc1dict(pchem,pname,fgrp)
    ierr=0               ! flag down, species found

!-- add chemistry
!       ..... WRITE THE REACTIONS (INCLUDING PHASE TRANSFER) ...
    np=0
    CALL rxinit(r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe)
    nref=0 ; ref(:)=' ' ; CALL addref(progname,'JO23PRSCOM',nref,ref,pchem)

    CALL add1tonp(progname,pchem,np) ; s(np)=1.0 ; p(np)='CH2O  '
    CALL add1tonp(progname,pchem,np) ; s(np)=1.0 ; p(np)='HO2   '

    r(1) = pname
    arrh(1)=0.1 ; arrh(2)=0. ; arrh(3)=0.
    r(2) = 'HV   ' ; idreac=1 ; nlabel=4 ; xlabel=arrh(1)
    CALL rxwrit(mecu,r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe,ref)

! ------------
! CHO(NO2)
! ------------
  ELSE IF (pchem=='CHO(NO2)') THEN

!-- add the species to the dictionary
    pname='CHONO2'  ;   fgrp='DV '
    CALL addc1dict(pchem,pname,fgrp)
    ierr=0               ! flag down, species found

!-- add chemistry
!       ..... WRITE THE REACTIONS (INCLUDING PHASE TRANSFER) ...
    np=0
    CALL rxinit(r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe)
    nref=0 ; ref(:)=' ' ; CALL addref(progname,'JO23PRSCOM',nref,ref,pchem)

    CALL add1tonp(progname,pchem,np) ; s(np)=1.0 ; p(np)='CO    '
    CALL add1tonp(progname,pchem,np) ; s(np)=1.0 ; p(np)='HO2   '
    CALL add1tonp(progname,pchem,np) ; s(np)=1.0 ; p(np)='NO2   '

    r(1) = pname
    arrh(1)=0.1 ; arrh(2)=0. ; arrh(3)=0.
    r(2) = 'HV   ' ; idreac=1 ; nlabel=4 ; xlabel=arrh(1)
    CALL rxwrit(mecu,r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe,ref)

! ------------
! CH2(ONO2)(OO.)
! ------------
  ELSE IF (pchem=='CH2(ONO2)(OO.)') THEN

!-- add the species to the dictionary
    pname='C10001'  ;   fgrp='2N  '
    CALL addc1dict(pchem,pname,fgrp)
    ierr=0               ! flag down, species found

!-- add chemistry
! decomposition
    np=0
    CALL rxinit(r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe)
    nref=0 ; ref(:)=' ' !; CALL addref(progname,'XXXXXXXXXX',nref,ref,pchem)

    CALL add1tonp(progname,pchem,np) ; s(np)=1.0 ; p(np)='CH2O  '
    CALL add1tonp(progname,pchem,np) ; s(np)=1.0 ; p(np)='NO2   '

    r(1)=pname
    arrh(1)=1.00E6 ; arrh(2)=0. ; arrh(3)=0.
    CALL rxwrit(mecu,r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe,ref)

! ------------
! CO(NO2)(ONO2)
! ------------
  ELSE IF (pchem=='CO(NO2)(ONO2)') THEN

!-- add the species to the dictionary
    pname='C10002'  ;   fgrp='KNV '
    CALL addc1dict(pchem,pname,fgrp)
    ierr=0               ! flag down, species found

!-- add chemistry
!       ..... WRITE THE REACTIONS (INCLUDING PHASE TRANSFER) ...

    np=0
    CALL rxinit(r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe)
    nref=0 ; ref(:)=' ' ; CALL addref(progname,'JO23PRSCOM',nref,ref,pchem)

    CALL add1tonp(progname,pchem,np) ; s(np)=1.0 ; p(np)='CO2  '
    CALL add1tonp(progname,pchem,np) ; s(np)=2.0 ; p(np)='NO2   '

    r(1) = pname
    arrh(1)=0.1 ; arrh(2)=0. ; arrh(3)=0.
    r(2) = 'HV   ' ; idreac=1 ; nlabel=4 ; xlabel=arrh(1)
    CALL rxwrit(mecu,r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe,ref)

! ------------
! CO(OH)(NO2)
! ------------
  ELSE IF (pchem=='CO(OH)(NO2)') THEN

!-- add the species to the dictionary
    pname='C10003'  ;   fgrp='AV '
    CALL addc1dict(pchem,pname,fgrp)
    ierr=0               ! flag down, species found

!-- add chemistry
!       ..... WRITE THE REACTIONS (INCLUDING PHASE TRANSFER) ...

    np=0
    CALL rxinit(r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe)
    nref=0 ; ref(:)=' ' ; CALL addref(progname,'JO23PRSCOM',nref,ref,pchem)

    CALL add1tonp(progname,pchem,np) ; s(np)=1.0 ; p(np)='HO2   '
    CALL add1tonp(progname,pchem,np) ; s(np)=1.0 ; p(np)='CO2   '
    CALL add1tonp(progname,pchem,np) ; s(np)=1.0 ; p(np)='NO2   '

    r(1) = pname
    arrh(1)=0.1 ; arrh(2)=0. ; arrh(3)=0.
    r(2) = 'HV   ' ; idreac=1 ; nlabel=4 ; xlabel=arrh(1)
    CALL rxwrit(mecu,r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe,ref)
! ------------
! CO(OH)(ONO2)
! ------------
  ELSE IF (pchem=='CO(OH)(ONO2)') THEN

!-- add the species to the dictionary
    pname='C10004'  ;   fgrp='AN '
    CALL addc1dict(pchem,pname,fgrp)
    ierr=0               ! flag down, species found

!-- add chemistry
!       ..... WRITE THE REACTIONS (INCLUDING PHASE TRANSFER) ...

    np=0
    CALL rxinit(r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe)
    nref=0 ; ref(:)=' ' ; CALL addref(progname,'JO23PRSCOM',nref,ref,pchem)

    CALL add1tonp(progname,pchem,np) ; s(np)=1.0 ; p(np)='HO2   '
    CALL add1tonp(progname,pchem,np) ; s(np)=1.0 ; p(np)='CO2   '
    CALL add1tonp(progname,pchem,np) ; s(np)=1.0 ; p(np)='NO3   '

    r(1) = pname
    arrh(1)=0.1 ; arrh(2)=0. ; arrh(3)=0.
    r(2) = 'HV   ' ; idreac=1 ; nlabel=4 ; xlabel=arrh(1)
    CALL rxwrit(mecu,r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe,ref)

! ------------
! CO(OH)(OOH)
! ------------
  ELSE IF (pchem=='CO(OH)(OOH)') THEN

!-- add the species to the dictionary
    pname='C10005'  ;   fgrp='AH '
    CALL addc1dict(pchem,pname,fgrp)
    ierr=0               ! flag down, species found

!-- add chemistry
!       ..... WRITE THE REACTIONS (INCLUDING PHASE TRANSFER) ...

    np=0
    CALL rxinit(r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe)
    nref=0 ; ref(:)=' ' ; CALL addref(progname,'JO23PRSCOM',nref,ref,pchem)

    CALL add1tonp(progname,pchem,np) ; s(np)=2.0 ; p(np)='HO    '
    CALL add1tonp(progname,pchem,np) ; s(np)=1.0 ; p(np)='CO2   '

    r(1) = pname
    arrh(1)=0.1 ; arrh(2)=0. ; arrh(3)=0.
    r(2) = 'HV   ' ; idreac=1 ; nlabel=4 ; xlabel=arrh(1)
    CALL rxwrit(mecu,r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe,ref)

! ------------
! CO(ONO2)(ONO2)
! ------------
  ELSE IF (pchem=='CO(ONO2)(ONO2)') THEN

!-- add the species to the dictionary
    pname='C10006'  ;   fgrp='NN '
    CALL addc1dict(pchem,pname,fgrp)
    ierr=0               ! flag down, species found

!-- add chemistry
!       ..... WRITE THE REACTIONS (INCLUDING PHASE TRANSFER) ...

    np=0
    CALL rxinit(r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe)
    nref=0 ; ref(:)=' ' ; CALL addref(progname,'JO23PRSCOM',nref,ref,pchem)

    CALL add1tonp(progname,pchem,np) ; s(np)=1.0 ; p(np)='NO3   '
    CALL add1tonp(progname,pchem,np) ; s(np)=1.0 ; p(np)='CO2   '
    CALL add1tonp(progname,pchem,np) ; s(np)=1.0 ; p(np)='NO2   '

    r(1) = pname
    arrh(1)=0.1 ; arrh(2)=0. ; arrh(3)=0.
    r(2) = 'HV   ' ; idreac=1 ; nlabel=4 ; xlabel=arrh(1)
    CALL rxwrit(mecu,r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe,ref)

! ------------
! CO(ONO2)(OOH)
! ------------
  ELSE IF (pchem=='CO(ONO2)(OOH)') THEN

!-- add the species to the dictionary
    pname='C10007'  ;   fgrp='NH '
    CALL addc1dict(pchem,pname,fgrp)
    ierr=0               ! flag down, species found

!-- add chemistry
!       ..... WRITE THE REACTIONS (INCLUDING PHASE TRANSFER) ...

    np=0
    CALL rxinit(r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe)
    nref=0 ; ref(:)=' ' ; CALL addref(progname,'JO23PRSCOM',nref,ref,pchem)

    CALL add1tonp(progname,pchem,np) ; s(np)=1.0 ; p(np)='HO    '
    CALL add1tonp(progname,pchem,np) ; s(np)=1.0 ; p(np)='CO2   '
    CALL add1tonp(progname,pchem,np) ; s(np)=1.0 ; p(np)='NO3   '

    r(1) = pname
    arrh(1)=0.1 ; arrh(2)=0. ; arrh(3)=0.
    r(2) = 'HV   ' ; idreac=1 ; nlabel=4 ; xlabel=arrh(1)
    CALL rxwrit(mecu,r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe,ref)

! ------------
! CO(OOH)(OO.)
! ------------
  ELSE IF (pchem=='CO(OOH)(OO.)') THEN

!-- add the species to the dictionary
    pname='C10008'  ;   fgrp='3H  '
    CALL addc1dict(pchem,pname,fgrp)
    ierr=0               ! flag down, species found

!-- add chemistry
! decomposition
    np=0
    CALL rxinit(r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe)
    nref=0 ; ref(:)=' ' !; CALL addref(progname,'XXXXXXXXXX',nref,ref,pchem)

    CALL add1tonp(progname,pchem,np) ; s(np)=1.0 ; p(np)='CO2   '
    CALL add1tonp(progname,pchem,np) ; s(np)=1.0 ; p(np)='HO    '

    r(1)=pname
    arrh(1)=1.00E6 ; arrh(2)=0. ; arrh(3)=0.
    CALL rxwrit(mecu,r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe,ref)

! ------------
! CO(NO2)(OOOH)
! ------------
  ELSE IF (pchem=='CO(NO2)(OOOH)') THEN

!-- add the species to the dictionary
    pname='C10009'  ;   fgrp='KNV '
    CALL addc1dict(pchem,pname,fgrp)
    ierr=0               ! flag down, species found

!-- add chemistry
!       ..... WRITE THE REACTIONS (INCLUDING PHASE TRANSFER) ...
    np=0
    CALL rxinit(r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe)
    nref=0 ; ref(:)=' ' ; CALL addref(progname,'JO23PRSCOM',nref,ref,pchem)

    CALL add1tonp(progname,pchem,np) ; s(np)=1.0 ; p(np)='HO2   '
    CALL add1tonp(progname,pchem,np) ; s(np)=1.0 ; p(np)='CO2   '
    CALL add1tonp(progname,pchem,np) ; s(np)=1.0 ; p(np)='NO2   '

!-- assign photolysis @ 0.1 x jNO2
    r(1)=pname
    arrh(1)=0.1 ; arrh(2)=0. ; arrh(3)=0.
    r(2) = 'HV   ' ; idreac=1 ; nlabel=4 ; xlabel=arrh(1)
    CALL rxwrit(mecu,r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe,ref)

! --------------
! CO(OOOH)(OOOH)
! --------------
  ELSE IF (pchem=='CO(OOOH)(OOOH)') THEN

!-- add the species to the dictionary
    pname='C10010'  ;   fgrp='KNV '
    CALL addc1dict(pchem,pname,fgrp)
    ierr=0               ! flag down, species found

!-- add chemistry
!       ..... WRITE THE REACTIONS (INCLUDING PHASE TRANSFER) ...
    np=0
    CALL rxinit(r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe)
    nref=0 ; ref(:)=' ' ; CALL addref(progname,'JO23PRSCOM',nref,ref,pchem)

    CALL add1tonp(progname,pchem,np) ; s(np)=1.0 ; p(np)='HO    '
    CALL add1tonp(progname,pchem,np) ; s(np)=1.0 ; p(np)='CO2   '
    CALL add1tonp(progname,pchem,np) ; s(np)=1.0 ; p(np)='HO2   '

!-- assign photolysis @ 0.1 x jNO2
    r(1)=pname
    arrh(1)=0.1 ; arrh(2)=0. ; arrh(3)=0.
    r(2) = 'HV   ' ; idreac=1 ; nlabel=4 ; xlabel=arrh(1)
    CALL rxwrit(mecu,r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe,ref)
              
! ------------
! CO(OH)(OH)
! ------------
  ELSE IF (pchem=='CO(OH)(OH)') THEN

!-- add the species to the dictionary
    pname='H2CO3'  ;   fgrp='D  '
    CALL addc1dict(pchem,pname,fgrp)
    ierr=0               ! flag down, species found

!-- add chemistry
!       ..... WRITE THE REACTIONS (INCLUDING PHASE TRANSFER) ...

    np=0
    CALL rxinit(r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe)
    nref=0 ; ref(:)=' ' ; CALL addref(progname,'JO23PRSCOM',nref,ref,pchem)

    CALL add1tonp(progname,pchem,np) ; s(np)=1.0 ; p(np)='HO    '
    CALL add1tonp(progname,pchem,np) ; s(np)=1.0 ; p(np)='CO2   '
    CALL add1tonp(progname,pchem,np) ; s(np)=1.0 ; p(np)='HO2   '

    r(1) = pname
    arrh(1)=0.1 ; arrh(2)=0. ; arrh(3)=0.
    r(2) = 'HV   ' ; idreac=1 ; nlabel=4 ; xlabel=arrh(1)
    CALL rxwrit(mecu,r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe,ref)

! ------------
! CH2(OH)(ONO2)
! ------------
  ELSE IF (pchem=='CH2(OH)(ONO2)') THEN
! The chemistry of this compound has to be checked. It's produced during 
! the oxidation of alkanes. For the moment, the constants are the same as CH3OH.
! It's assumed that only a reaction with OH occurs.

!-- add the species to the dictionary
    pname='NO1001'  ;   fgrp='NO  '
    CALL addc1dict(pchem,pname,fgrp)
    ierr=0               ! flag down, species found

!-- add chemistry
! decomposition
    np=0
    CALL rxinit(r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe)
    nref=0 ; ref(:)=' ' !; CALL addref(progname,'XXXXXXXXXX',nref,ref,pchem)

    CALL add1tonp(progname,pchem,np) ; s(np)=1.0 ; p(np)='HCOOH '
    CALL add1tonp(progname,pchem,np) ; s(np)=1.0 ; p(np)='NO2   '

    r(1)=pname ; r(2)='HO    '
    arrh(1)=3.10E-12 ; arrh(2)=0. ; arrh(3)=0.
    CALL rxwrit(mecu,r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe,ref)

! ------------
! CH2(ONO2)(OOH)
! ------------
  ELSE IF (pchem=='CH2(ONO2)(OOH)') THEN
! The chemistry of this compound has to be checked. It's produced by 
! the CH2(.OO.)+HNO3 reaction. The rate constant is based on Jenkin et
! al., 2018 SAR and assume H abstraction from the CH2 group only. 
! Reaction products are preliminary. 

!-- add the species to the dictionary
    pname='NH1001'  ;   fgrp='NH  '
    CALL addc1dict(pchem,pname,fgrp)
    ierr=0               ! flag down, species found

!-- add chemistry
    np=0
    CALL rxinit(r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe)
    nref=0 ; ref(:)=' ' ; 
    CALL addref(progname,'MJ18KMV000',nref,ref,pchem)
    CALL addref(progname,'NH1001COM1',nref,ref,pchem)

    CALL add1tonp(progname,pchem,np) ; s(np)=1.0 ; p(np)='HCOO2H '
    CALL add1tonp(progname,pchem,np) ; s(np)=1.0 ; p(np)='NO3   '
    CALL add1tonp(progname,pchem,np) ; s(np)=1.0 ; p(np)='H2O   '

    r(1)=pname ; r(2)='HO    '
    arrh(1)=4.95E-12 ; arrh(2)=0. ; arrh(3)=719.
    CALL rxwrit(mecu,r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe,ref)

! ------------
! CH2(NO2)(OO.)
! ------------
! The chemistry of this compound still need to be developped. The 
! species is produced by the CH2(.OO.)+NO2 reaction. To avoid the 
! accumulation of the radical, decomposition is performed with an 
! arbitrary high rate constant.  
  ELSE IF (pchem=='CH2(NO2)(OO.)') THEN

!-- add the species to the dictionary
    pname='2V1001'  ;   fgrp='2V  '
    CALL addc1dict(pchem,pname,fgrp)
    ierr=0               ! flag down, species found

!-- add chemistry
! decomposition
    np=0
    CALL rxinit(r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe)
    nref=0 ; ref(:)=' ' ; CALL addref(progname,'RO2TEMPO',nref,ref,pchem)

    CALL add1tonp(progname,pchem,np) ; s(np)=1.0 ; p(np)='CH2O  '
    CALL add1tonp(progname,pchem,np) ; s(np)=1.0 ; p(np)='NO3   '

    r(1)=pname
    arrh(1)=1.00E6 ; arrh(2)=0. ; arrh(3)=0.
    CALL rxwrit(mecu,r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe,ref)

! ------------
! CHO(ONO2)
! ------------
  ELSE IF (pchem=='CHO(ONO2)') THEN
! this species is produced as a minor product during the oxidation of 
! some cyclic HC. We assume that this species behaves like R-CO(ONO2)

!-- add the species to the dictionary
    pname='C1ND00'  ;   fgrp='N  '
    CALL addc1dict(pchem,pname,fgrp)
    ierr=0               ! flag down, species found

!-- add chemistry
! decomposition
    np=0
    CALL rxinit(r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe)
    nref=0 ; ref(:)=' ' !; CALL addref(progname,'XXXXXXXXXX',nref,ref,pchem)

    CALL add1tonp(progname,pchem,np) ; s(np)=1.0 ; p(np)='CO2   '
    CALL add1tonp(progname,pchem,np) ; s(np)=1.0 ; p(np)='HO2   '
    CALL add1tonp(progname,pchem,np) ; s(np)=1.0 ; p(np)='NO2   '

    r(1)=pname
    arrh(1)=1.00E1 ; arrh(2)=0. ; arrh(3)=0.
    CALL rxwrit(mecu,r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe,ref)

! ------------
! CH3(ONO)
! ------------
!  ELSE IF (pchem=='CH3(ONO)') THEN
!
!!-- add the species to the dictionary
!    pname='N01003'  ;   fgrp='   '
!    CALL addc1dict(pchem,pname,fgrp)
!    ierr=0               ! flag down, species found
!
!!-- add chemistry
!! reaction with OH
!    np=0
!    CALL rxinit(r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe)
!    nref=0 ; ref(:)=' ' ; !CALL addref(progname,'XXXXXXXXX',nref,ref,pchem)
!
!    CALL add1tonp(progname,pchem,np) ; s(np)=0.55 ; p(np)='CH2O  '
!    CALL add1tonp(progname,pchem,np) ; s(np)=0.45 ; p(np)='NO    '
!
!    r(1)=pname ; r(2)='HO '
!    arrh(1) = 0.301E-12  ; arrh(2)=0. ; arrh(3)=0.
!    CALL rxwrit(mecu,r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe,ref)
!    
!! photolysis
!    np=0
!    CALL rxinit(r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe)
!    nref=0 ; ref(:)=' ' ; !CALL addref(progname,'XXXXXXXXX',nref,ref,pchem)
!
!    CALL add1tonp(progname,pchem,np) ; s(np)=0.55 ; p(np)='CH3O  '
!    CALL add1tonp(progname,pchem,np) ; s(np)=0.45 ; p(np)='NO    '
!
!    r(1)=pname
!    idreac=1 ; nlabel=40000   
!    arrh(1) = 1.00E00  ; arrh(2)=0. ; arrh(3)=0.
!    CALL rxwrit(mecu,r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe,ref)
!    

  ENDIF

END SUBROUTINE addc1chem

! ======================================================================
! PURPOSE: Add the formula (pchem) and name (pname) provided as input 
! in the dictionaries.
! ======================================================================
SUBROUTINE addc1dict(pchem,pname,fgrp)
  USE dictstackdb, ONLY: nrec,dict,namlst,dbrch,diccri
  USE keyflag, ONLY: isomerfg
  USE searching, ONLY: srch, srh5
  IMPLICIT NONE

  CHARACTER(LEN=*),INTENT(IN) :: pchem  ! formula to be added to dict 
  CHARACTER(LEN=*),INTENT(IN) :: pname  ! name to be added in namlst 
  CHARACTER(LEN=*),INTENT(IN) :: fgrp   ! functionalities for pchem to be added 

  INTEGER :: dicptr,namptr
  CHARACTER(LEN=10) :: progname='addc1dict'
  CHARACTER(LEN=70) :: mesg

! get pointer to add the species in namlst and dict
  dicptr=srch(nrec,pchem,dict)
  namptr=srh5(pname,namlst,nrec)
  IF (dicptr>0) THEN
    mesg="formula for a C1 species to be added to dict already exists"
    CALL stoperr(progname,mesg,pchem)
  ENDIF
  IF (namptr>0) THEN
    mesg="short name assigned to a C1 species already exists"
    CALL stoperr(progname,mesg,pchem)
  ENDIF

! raise the counters
  nrec=nrec+1  
  IF (nrec > SIZE(dict)) THEN
    mesg="Too many species added to dictionary"
    CALL stoperr(progname,mesg,pchem)
  ENDIF 

! add new name and raise the name array
  namptr=ABS(namptr)+1
  namlst(namptr+1:nrec+1)=namlst(namptr:nrec)
  namlst(namptr)=pname

! raise dictionary and branching table. Add new line for the new species 
  dicptr = ABS(dicptr) + 1
  dict(dicptr+1:nrec+1)=dict(dicptr:nrec)
  dbrch(dicptr+1:nrec+1)=dbrch(dicptr:nrec)
  WRITE(dict(dicptr),'(a6,3x,a120,2x,a15)') pname, pchem, fgrp
  dbrch(dicptr)  = 1.0  ! set max  yield to 1 (C1 species, number not needed)

! raise index of isomer table
  IF (isomerfg) THEN
    diccri(dicptr+1:nrec+1,:) = diccri(dicptr:nrec,:)
    diccri(dicptr,:)=0
  ENDIF

END SUBROUTINE addc1dict

END MODULE dictstacktool
