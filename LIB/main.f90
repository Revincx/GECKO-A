! =======================================================================
! GECKO-A program
! 
! PURPOSE: Generate the oxidation mechanism for organic compounds under
! tropospheric conditions.
!
! DEVELOPERS:
!
! COPYRIGHTS: All rights reserved.
!
! DISCLAIMER:
!=======================================================================
PROGRAM main
  USE keyparameter, ONLY:mxring,mxlgr,mxlfo,mxlco,mxps,mxnode,dirout,&
                         logu,mecu,scru,kohu,kno3u,waru,gasu,prtu,walu,&
                         tfu1,ohu,no3u,o3u,dhfu,saru,refu,prmu,rszu
  USE keyflag,    ONLY: critvp,g2pfg,g2wfg, &
                        dhffg,rx_ro2_no2,wrtsarinfo,sar_only_fg,pvap_sar, &
                        wrtref,wrtdhf,wrtpvap,wrthenry,wrttg,wrtdepo, &
                        wrtkivoci,wrtmaxyield
  USE rjtool,     ONLY: rjgrm
  USE stdgrbond,  ONLY: grbond
  USE searching,  ONLY: srch
  USE sortstring, ONLY: sort_string
  USE reac_infotool,  ONLY: rdreac_info            ! to read comments/references in databases
  USE spsptool,   ONLY: rdoutgene, spreact ! for special species mechanism (i.e. #species)
  USE loaddbtool, ONLY: loaddb             ! to load the database
  USE loadc1tool, ONLY: loadc1mch          ! to load inorganic & C1 chemistry
  USE dictstackdb,ONLY: nrec,dict,namlst,dbrch, &                      !  dictionary data
                        nhldvoc,holdvoc,nhldrad,holdrad,stabl,level, & ! stack data
                        lotopstack                                     ! fifo (default) or lifo
  USE loadchemin, ONLY: rdchemin,in1chm  ! to load the input (primary) species
  USE tempflag,   ONLY: iflost
  USE masstranstool, ONLY:changephase
  USE outtool,    ONLY: dictelement,wrt_henry,wrt_depo,wrt_heatf,wrt_Tg,&
                        wrt_dict,wrt_psat,wrt_mxyield,wrt_ro2,&
                        wrt_size,diffusion_vol
  USE rochem,     ONLY: ro
  USE hochem,     ONLY: ho_voc
  USE no3chem,    ONLY: no3_voc
  USE o3chem,     ONLY: o3_voc
  USE panchem,    ONLY: pandec
  USE rooohchem,  ONLY: rooohdec
  USE stbcriegeechem, ONLY: stab_criegee
  USE ro2chem,    ONLY: ro2
  USE rco3chem,   ONLY: rco3
  USE toolbox,    ONLY: stoperr
  USE dhftool,    ONLY: dhf, dhf_thf
  USE logtool,    ONLY: wrtlog
  USE maintool,   ONLY: initdictstack,sortstack,nphotogrp,wrtscreen, &
                        vaporpressure,InOrOut,cut_path,cleanstack,check_parenthesis
  USE rdnml, ONLY:rd_nml  
  IMPLICIT NONE

! current species info
  CHARACTER(LEN=mxlfo) :: chem          ! chemical formula of the species being processed
  CHARACTER(LEN=mxlfo) :: parent        ! parent compound formula 
  CHARACTER(LEN=mxlco) :: idnam         ! ID name of the species in the mechanism
  CHARACTER(LEN=mxlgr) :: group(mxnode) ! groups in chem
  INTEGER :: bond(mxnode,mxnode)        ! node bond matrix of chem
  INTEGER :: zebond(mxnode,mxnode)      ! cis/trans info on C=C bond in chem 
  REAL    :: brch
  INTEGER :: rjg(mxring,2)              ! ring-join group pairs
  INTEGER :: dbflg, nring
  LOGICAL :: lorad

  INTEGER  :: ninp                           ! # of "primary" species 
  CHARACTER(LEN=mxlfo) :: input(mxps)        ! table of the primary species (formula)  
  REAL    :: log10Psat                       ! data for Pvap estimates
  INTEGER :: ired
  INTEGER :: ind,i,k
  CHARACTER(LEN=6)   :: progname='main '
  CHARACTER(LEN=256) :: mesg
  CHARACTER(LEN=256) :: filename

! branching ratio below which reaction pathways are not considered
  REAL :: cut_off,cut_default,cut_OH,cut_O3,cut_NO3,cut_PAN
  REAL :: cut_HV,cut_RO,cut_RO2,cut_RCOO2

! Get default values fpr user-settable parameters & flags
  CALL define_defaults

! Read namelist
  CALL rd_nml

! initialize data in the dictionaries and stacks
! ----------------------------------------------
  CALL initdictstack()
  brch=1. 

! ------------------------------------------
! OPEN SOME OUTPUT FILES
! ------------------------------------------
  OPEN(logu,FILE=TRIM(dirout)//'scheme.log')  ;  CALL wrtlog() ! operating conditions and flags
  OPEN(prmu,FILE=TRIM(dirout)//'listprimary.dat')
  OPEN(scru,FILE=TRIM(dirout)//'screeninfo.out')               !
  IF (wrtkivoci) THEN
    OPEN(kohu,FILE=TRIM(dirout)//'kivoci.dat')                   !
    OPEN(kno3u,FILE=TRIM(dirout)//'kjvocj.dat')                  !
  ELSE
    OPEN(kohu,FILE='/dev/null')                   !
    OPEN(kno3u,FILE='/dev/null')                  !
  ENDIF
  OPEN(waru,FILE=TRIM(dirout)//'warning.out')                  !
  OPEN(mecu,FILE=TRIM(dirout)//'reactions.dum')                !
  IF (wrtref) OPEN(refu,FILE=TRIM(dirout)//'reactionswithcom.dum')         !
  OPEN(rszu,FILE=TRIM(dirout)//'stack_size_record.dat')        ! info file on #species & stack
  WRITE(mecu,'(a)') "END SPECIES"; WRITE(mecu,'(a)') "!----------------"
  WRITE(mecu,'(a)') "REACTIONS " ; WRITE(mecu,'(a)') "!----------------"

  IF (sar_only_fg) THEN
    OPEN(ohu, FILE=TRIM(dirout)//'rateoh.dat')                ! OH rate constant 
    OPEN(o3u, FILE=TRIM(dirout)//'rateo3.dat')                ! O3 rate constant
    OPEN(no3u,FILE=TRIM(dirout)//'rateno3.dat')               ! NO3 rate constant
  ENDIF
  IF (dhffg) OPEN(dhfu,FILE=TRIM(dirout)//'dhfrx.dat')
  IF (wrtsarinfo) OPEN(saru,FILE=TRIM(dirout)//'wrtsarinfo.dat')

! ------------------------------------------
! READ ALL DATA
! ------------------------------------------
! read comment/references dictionnary and relate code
  CALL rdreac_info()

! Read the dictionaries and mechanisms for C1, inorg ... and sort
  CALL loadc1mch()      ! read inorganic and C1 mechanism
  CALL wrt_kc1()        ! ... write koh for C1 species

! read all data to be stored in the database module
  CALL loaddb()

! load species from cheminput.dat 
  WRITE (6,*) "  ...reading the list of the primary species"
  filename='./cheminput.dat'
  CALL rdchemin(filename,ninp,input)

! -----------------------------------------------
! LOAD INPUT SPECIES - START LOOP OVER PARENTS
! -----------------------------------------------
  psloop: DO k=1,ninp  ! primary species loop

    chem=input(k)  ;  parent=chem  ;  iflost=0  ;  lotopstack=.FALSE.
    WRITE(logu,*) 'Precursor: ',TRIM(parent)

! cycle if "*" as 1st character (note: chem might be simple string)
    IF (chem(1:1) == "*") THEN
      IF (sar_only_fg) WRITE(ohu, '(1PE12.3,2x,a)') 0., TRIM(chem)
      IF (sar_only_fg) WRITE(o3u, '(1PE12.3,2x,a)') 0., TRIM(chem)
      IF (sar_only_fg) WRITE(no3u,'(1PE12.3,2x,a)') 0., TRIM(chem)
      CYCLE psloop
    ENDIF

! check and name the species, update the dictionary and put the species
! at the beginning of the stack (i.e. in holdvoc(1))
    CALL check_parenthesis(chem)  
    CALL in1chm(chem,idnam)
! write primary species (short name) in listprimary.dat
    WRITE(prmu,*) idnam,' ',TRIM(chem)

! sort namlst
    CALL sort_string(namlst(1:nrec))

! highest yield species treated 1st (required for isomer substitution)
    ired=1   ! target to sort the stack (when stabl==ired) 

! ------------------------------------
! LOOP OVER SPECIES IN THE STACKS
! ------------------------------------
! Give the priority to the stack holding radical species (holdrad). If this 
! stack is empty, then take the next stable VOC in the stack (holdvoc). 
    stackloop: DO 

! load 1st radical in stack
      IF (holdrad(1)/=' ') THEN   
        READ (holdrad(1),'(a6,a120,i3,i3)') idnam,chem,stabl,level
        lorad=.TRUE. ; lotopstack=.FALSE.

! load 1st VOC in stack
      ELSE                        
        IF (holdvoc(1)==' ') EXIT stackloop ! empty stack, load next primary species
        READ (holdvoc(1),'(a6,a120,i3,i3)') idnam,chem,stabl,level
        lorad=.FALSE. ; lotopstack=.FALSE.
     
        IF (stabl==ired) THEN       ! new generation - sort the stack (highest yield at top)
          ! record the stack and dictionary size of the current generation (information only)
          WRITE(rszu,'(i3,3(i8))') ired,nrec,nhldvoc,nhldrad
          ired=ired+1
          CALL sortstack(nrec,dict,dbrch,nhldvoc,holdvoc)
          READ (holdvoc(1),'(a6,a120,i3,i3)') idnam,chem,stabl,level
        ENDIF
      ENDIF

! remove the current species in the stack 
      IF (lorad) THEN ; holdrad(1:nhldrad)=holdrad(2:nhldrad+1) ; nhldrad=nhldrad-1 
      ELSE            ; holdvoc(1:nhldvoc)=holdvoc(2:nhldvoc+1) ; nhldvoc=nhldvoc-1
      ENDIF

! screen writing of current species (see possible options in the routine)
      CALL wrtscreen(nrec,nhldvoc,stabl,idnam,chem)

! special chemistry (don't call reaction subroutine)
      IF (chem(1:1)=='#') THEN
        CALL spreact(idnam,chem,brch)   ! perform special reactions
        CYCLE stackloop                 ! load next species
      ENDIF

! rebuild species tables (bond, group, ...) and check that the species 
! is effectively known in the dictionary (just in case).
      ind=srch(nrec,chem,dict)
      IF (ind<0) THEN
        DO i=1,nrec ;  WRITE(18,'(a)') dict(i)  ;  ENDDO
        mesg ="species loaded from stack is unknown in dictionary"
        CALL stoperr(progname,mesg,chem)
      ENDIF
      brch=dbrch(ind)
      IF (brch<1E-20) brch=1E-20
      CALL grbond(chem,group,bond,dbflg,nring,zebond)
      IF (nring>0)  CALL rjgrm(nring,group,rjg)       ! rm ring characters in group

! SEND SPECIES TO REACTION SUBROUTINES
! -------------------------------------

! set threshold for reaction paths and restore reaction flag (iflost)
      CALL cut_path(brch,cut_default,cut_OH,cut_O3,cut_NO3,cut_PAN,&
                    cut_HV, cut_RO,cut_RO2,cut_RCOO2) 
      iflost=0              ! restore to default, make full reaction
     
! STABLE VOC
! =============
      IF (.not.lorad) THEN

! compute the vapor pressure - exit if below threshold
        log10Psat=vaporpressure(chem,bond,group,nring,rjg)
        IF (log10Psat<critvp) CYCLE stackloop        ! load next species

! check if chemistry wanted (see iflost flag: VOC+Ox -> LostCarbon)
        CALL InOrOut(chem,parent,brch,stabl)

! Conversion to dihydrofuran
        IF ((dhffg).AND.(stabl<2)) THEN
          IF (((INDEX(chem,'CHO')/=0).OR.(INDEX(chem,'CO')/=0)).AND.(INDEX(chem,'OH')/=0)) THEN
            cut_off=cut_default
            IF (g2pfg .OR. g2wfg) THEN
              CALL dhf_thf(idnam,chem,bond,group,nring,brch,cut_off)
            ELSE 
              CALL dhf(idnam,chem,bond,group,nring,brch,cut_off)
            ENDIF
          ENDIF
        ENDIF

! PAN decomposition 
        IF (INDEX(chem,'OONO2')/=0) THEN
          CALL pandec(idnam,chem,bond,group,brch)
          IF (INDEX(chem,'CO(OONO2)')==0) CYCLE stackloop  ! alkyl, not acyl moiety
          IF (rx_ro2_no2) CYCLE stackloop  ! no further chemistry after RO2NO2 decomp
        ENDIF

! ROOOH decomposition 
        IF (INDEX(chem,'(OOOH)')/=0) THEN
          CALL rooohdec(idnam,chem,bond,group,brch)
        ENDIF

! reaction with OH 
        IF ((INDEX(chem,'H')/=0).OR.(INDEX(chem,'Cd')/=0)) THEN
          CALL ho_voc(idnam,chem,bond,group,nring,brch,cut_OH)
        ENDIF
            
! reaction with O3
        IF (INDEX(chem,'Cd')/=0) THEN
          CALL o3_voc(idnam,chem,bond,group,zebond,brch,cut_O3)
        ENDIF

! reaction with NO3
        IF((INDEX(chem,'H')/=0).OR.(INDEX(chem,'Cd')/=0))THEN
          CALL no3_voc(idnam,chem,bond,group,brch,cut_NO3)
        ENDIF
        
! test for chromophores and perform reaction with HV
        IF (nphotogrp(chem)/=0) THEN
          CALL hvdiss2(idnam,chem,bond,group,nring,brch,cut_HV)
        ENDIF
 
      ENDIF
         
! RADICAL VOC
! =============
      IF (lorad) THEN

! alkoxy radical
        IF (idnam(1:1)=='1') THEN
          CALL ro(idnam,chem,bond,group,nring,rjg,brch,cut_RO)

! peroxy radical
        ELSE IF (idnam(1:1)=='2')  THEN
          CALL ro2(idnam,chem,bond,group,nring,brch,cut_RO2)

! acyl peroxy radical
        ELSE IF (idnam(1:1)=='3') THEN
          CALL rco3(idnam,chem,bond,group,brch,cut_RCOO2)

! criegee radical
        ELSE IF (idnam(1:1)=='4') THEN
          CALL stab_criegee(idnam,chem,bond,group,brch)

! unidentified radical
        ELSE
          mesg ="Radical type not found."
          CALL stoperr(progname,mesg,chem)
        ENDIF 

      ENDIF

! for SAR assessment, clean stack to jump to the next parent (psloop)
      IF (sar_only_fg) CALL cleanstack() 
      
    ENDDO stackloop

! load the next "primary" species in the input table
  ENDDO psloop

! ----------------------------------
! END OF REACTIONS - WRITE DATA OUT
! ----------------------------------
  DO i=1,nrec-1
    IF (dict(i)(10:120)==dict(i+1)(10:120)) THEN
      WRITE(*,*) 'duplicate=', dict(i)(1:70)
      WRITE(*,*) 'duplicate=', dict(i+1)(1:70)
      IF (.NOT.sar_only_fg) STOP "duplicate in chem"
    ENDIF
    IF (namlst(i)==namlst(i+1)) THEN
      STOP "duplicate in namlst"
    ENDIF
  ENDDO

! close files - kivoci, kjvocj
  WRITE(kohu,'(a)') "END   "  ; CLOSE (kohu)
  WRITE(kno3u,'(a)') "END   "  ; CLOSE (kno3u)
  IF (sar_only_fg) THEN ; CLOSE(ohu) ; CLOSE(o3u) ; CLOSE(no3u) ; ENDIF
      
! compute molar mass and atoms for species (stored in dctmw & dctatom)
  PRINT*, '... compute molar masses'
  CALL dictelement()

! write the dictionaries (gas phase mechanism species and dictionary)
  PRINT*, '... write dictionary'
  OPEN(gasu,FILE=TRIM(dirout)//'gasspe.dum')  
  WRITE(gasu,'(a)') "SPECIES"  ;  WRITE(gasu,'(a)') "PHASE: START GAS"
  CALL wrt_dict()
  WRITE(gasu,'(a)') "PHASE: END GAS"  ;  CLOSE(gasu) 
      
! write max yield for each species in the mechanism
  IF (wrtmaxyield) THEN
    PRINT*, '... write max yields'
    CALL wrt_mxyield()
  ENDIF

! write peroxy species in counting files
  CALL wrt_ro2()

! write mass transport equations and species in the various phase
  PRINT*, '... write mass transfer equation (if any)'
  OPEN(prtu,FILE=TRIM(dirout)//'partspe.dum') ; WRITE(prtu,'(a)') "PHASE: START PART."
  OPEN(walu,FILE=TRIM(dirout)//'wallspe.dum') ; WRITE(walu,'(a)') "PHASE: START WALL"
  IF (g2pfg .OR. g2wfg) CALL changephase()
  WRITE(prtu,'(a)') "PHASE: END PART."  ;  CLOSE(prtu)
  WRITE(walu,'(a)') "PHASE: END WALL"   ;  CLOSE(walu)

! write diffusion volume
  IF (g2pfg) CALL diffusion_vol()

! write the vapor pressure of the species in dict.
  IF (wrtpvap) THEN
    PRINT*, '... compute data for Psat evaluation'
    IF      (pvap_sar==1) THEN ; CALL wrt_psat(1,0,0)
    ELSE IF (pvap_sar==2) THEN ; CALL wrt_psat(0,1,0)
    ELSE IF (pvap_sar==3) THEN ; CALL wrt_psat(0,0,1)
    ENDIF
  ENDIF

! write Henry's law coef. for species in dict.
  IF (wrthenry) THEN
    PRINT*, '... compute data for Henry parameters'
    CALL wrt_henry()
  ENDIF

! write deposition parameters for species in dict.
  IF (wrtdepo) THEN
    PRINT*, '... compute data for Deposition'
    CALL wrt_depo()
  ENDIF

! write the heat of formation for species in dict.
  IF (wrtdhf) THEN
    PRINT*, '... compute heats of formation'
    CALL wrt_heatf()
  ENDIF

! Write Tg for viscosity (calculated in the boxmod)
  IF (wrttg) THEN
    PRINT*, '... compute Tg'
    CALL wrt_Tg()
  ENDIF

! write numbers providing mechanism size (species, numbers ...)
  CALL wrt_size()
  
! close files still open
  WRITE(mecu,'(a3)') "END"  ; CLOSE(mecu)
  CLOSE(logu) ; CLOSE(prmu) ; CLOSE(rszu)
  WRITE(waru,*) "" ; CLOSE (waru) 
  mesg="cat "//TRIM(dirout)//"warning.out"  ; CALL SYSTEM(mesg)

END PROGRAM main
