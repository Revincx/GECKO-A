MODULE spsptool
IMPLICIT NONE
CONTAINS

! ======================================================================
! PURPOSE: read the list of known 'species' for which the stucture  
! can not be handled by the generator (e.g. cyclic molecule like furan 
! made by some photolytic reaction). Species are in the database as 
! dictionary line (short name, formula, funct. grp).
! ======================================================================
SUBROUTINE rdspsp(filename)
  USE keyparameter, ONLY: tfu1
  USE database, ONLY:  mxspsp, nspsp, lospsp, dictsp  
  IMPLICIT NONE

  CHARACTER(LEN=*),INTENT(IN):: filename      ! file to be read

  CHARACTER(LEN=LEN(dictsp(1))) :: line
  INTEGER :: ierr

! Initialize the database
  nspsp=0  ;  dictsp(:)=' ' ; lospsp(:)=.FALSE.  

! read the list of the special species.
  OPEN(tfu1,FILE=filename, FORM='FORMATTED', STATUS='OLD')

  rdloop: DO 
    READ(tfu1,'(a)',IOSTAT=ierr) line
    IF (ierr /= 0) THEN
      WRITE (6,*) '--error-- in rdspsp. Missing keyword END ?'
      STOP "in rdspsp, while reading inputs"
    ENDIF
    IF (line(1:1)=='!') CYCLE rdloop
    IF (line(1:3)=='END') EXIT rdloop

    nspsp = nspsp + 1
    IF (nspsp > mxspsp) THEN
      WRITE(6,*) '--error-- in rdspsp subroutine. The # of species in the "special" dictionary'
      WRITE(6,*) 'exceed maximum size allowed. Change the "mxspsp" parameter in database'
      STOP "in rdspsp"
    ENDIF
    READ (line,'(a)') dictsp(nspsp)
  ENDDO rdloop
  CLOSE(tfu1)

END SUBROUTINE rdspsp

! ======================================================================
! PURPOSE: read the chemical scheme that cannot be managed by the  
! generator (e.g. alkyne chemistry). The program checks if the species 
! in a reaction are "known" in the generator (list of special species spsp).
! If not, the program stops. Data of each reaction are store in the "oge"
! database (out of gecko):                                 
!  - noge           ! # of reactions given (oge = "out generator")
!  - ogelab(:)      ! label for the reaction (if EXTRA or HV)
!  - ogertve(:,3)   ! reactants for reaction i
!  - ogeprod(:,mxpd)! products for reaction i
!  - ogearh(:,3)    ! arrehnius coefficients (A, n, Ea)
!  - ogestoe(:,mxpd)! stochiometric coefficients for product j 
!  - ogeaux(:,7)    ! aux. info. for reaction i (e.g. falloff react.)
! ======================================================================
SUBROUTINE rdoutgene(filename)
  USE keyparameter, ONLY: tfu1
  USE database, ONLY: noge,ogelab,ogertve,ogeprod,ogearh,ogestoe,ogeaux ! out generator database 

  IMPLICIT NONE

  CHARACTER(LEN=*),INTENT(IN):: filename      ! file to be read

  INTEGER,PARAMETER :: lenlin=200
  CHARACTER(LEN=lenlin)            :: line,line1,line2,line3
  CHARACTER(LEN=LEN(ogertve(1,1))) :: chem1,chem2 
  INTEGER  :: n1,n2,i,j,nlin,ierr,ifo
  LOGICAL  :: locheck, lohv

! -------------
! INITIALIZE
! -------------
  noge=0 ; ogelab(:)=0 ; ogearh(:,:)=0. ; ogeaux(:,:)=0. ; ogestoe(:,:)=0.
  ogertve(:,:)=' ' ;  ogeprod(:,:)=' '

! -------------
! OPEN THE FILE
! -------------
  OPEN(tfu1,FILE=filename, FORM='FORMATTED',STATUS='OLD',IOSTAT=ierr)
  IF (ierr /= 0) THEN
    WRITE (6,*) '--error-- in rdoutgene, when opening input file'
    STOP "in rdgeneout (error to open file)"
  ENDIF

! -------------
! READ THE FILE
! -------------
  readloop: DO
    READ(tfu1,'(a)',IOSTAT=ierr) line
    IF (ierr /= 0) THEN
      WRITE (6,*) '--error-- in rdoutgene, while reading line:',TRIM(line)
      STOP "in rdgeneout, while reading inputs"
    ENDIF
  
    IF (line(1:3)=='END') EXIT readloop
    IF (line(1:1)=='!') CYCLE readloop

! NEW REACTION - GET THE NUMBER OF LINE TO BE READ
! ------------------------------------------------
    IF (line(1:4)=='REAC') THEN
      noge=noge+1

! get the total number of lines necessary to read the reaction
      READ(line(5:),*,IOSTAT=ierr) nlin
      IF (ierr/=0) THEN
        WRITE(6,*) ' --error--, in rdoutgene. Cannot read # of line :'
        WRITE(6,*) TRIM(line)
        STOP "in rdoutgene"
      ENDIF
      IF (nlin <= 0) THEN
        WRITE(6,*) '--error--, in rdoutgene. # of line < 1 in:'
        WRITE(6,*) TRIM(line)
        STOP "in rdoutgene"
      ENDIF

! READ THE REACTANTS, ARRHENIUS PARAMETERS, LABEL
! --------------------------------------------------
      READ (tfu1,'(a)',IOSTAT=ierr) line
      IF (ierr /= 0) THEN
        WRITE (6,*) '--error-- in rdoutgene, while reading line:',TRIM(line)
        STOP 'in rdgeneout, while reading inputs'
      ENDIF
      n1=INDEX(line,';') ; n2=INDEX(line,'&') 
      IF ((n1==0).OR.(n2==0)) THEN
        WRITE(6,*) ' --error--, in rdoutgene. ";" and "&" expected in:'
        WRITE(6,*) TRIM(line)
        STOP "in rdoutgene"
      ENDIF

! split the line into 3 parts
      line1=' ' ; line2=' ' ; line3=' '
      line1=line(1:n1-1)     ! 2 reactants
      line2=line(n1+1:n2-1)  ! 3 Arrhenius coefficients
      line3=line(n2+1:)      ! auxiliary info (for fall-off reaction only)

! get the 2 reactants - first the special species then the coreactant
      chem1=' ' ;  chem2=' '
      IF (line1(1:1)==' ') THEN
        WRITE(6,*) ' --error--, in rdoutgene. 1st pos not a char in line :'
        WRITE(6,*) TRIM(line1)
        STOP  "in rdoutgene"
      ENDIF
   
      n2=INDEX(line1,' + ') ; ifo=INDEX(line,'(+M)')
      IF (n2/=0) THEN
        chem1=line1(1:n2)
      ELSE 
        chem1=line1(1:n1-1)
      ENDIF
      
      IF (n2/=0) THEN
        chem2=line1(n2+3:)
        IF (chem2(1:1)==' ') THEN
          WRITE(6,*) ' --error--, in rdoutgene, while reading 2nd reactant'
          WRITE(6,*) ' separation between the 2 species must be " + "'
          WRITE(6,*) TRIM(line)
          STOP "in rdoutgene"
        ENDIF
      ENDIF

! check for fall-off reaction - remove the flag '(+M)' stored in chem
      IF (ifo/=0) THEN
        i=INDEX(chem1,'(+M)')
        IF (i/=0) chem1(i:i+4)='    '
        i=INDEX(chem2,'(+M)')
        IF (i/=0) chem2(i:i+4)='    '
      ENDIF

! check that the special species are known
      CALL chcksp(chem1)
      IF (n2/=0) CALL chcksp(chem2)

! put the species in the table of reactants
      ogertve(noge,1)=chem1  ;  ogertve(noge,2)=chem2

! add the fall off flag if required
      IF (ifo/=0) THEN
        IF (ogertve(noge,2)(1:1)==' ') THEN 
           ogertve(noge,2)='(+M)'
        ELSE
           ogertve(noge,3)='(+M)'
        ENDIF
      ENDIF

! get the arrhenius coefficients (must read 3 values)
      READ(line2,*,IOSTAT=ierr) (ogearh(noge,i),i=1,3) 
      IF (ierr/=0) THEN
        WRITE(6,*) '--error-- in the routine rdoutgene.'
        WRITE(6,*) 'Cannot read the 3 arrhenius coefficients at line:'
        WRITE(6,*) TRIM(line)
        STOP "in rdoutgene"
      ENDIF

! get the label (only if EXTRA or HV is the coreactant)
      locheck=.false.  ;  lohv=.false.
      IF (chem2(1:6)=='EXTRA ') locheck=.true.
      IF (chem2(1:3)=='HV ') THEN  
        locheck=.true.
        lohv=.true.
      ENDIF
      IF (locheck) THEN
        IF (lohv) THEN
          READ(line3,*,IOSTAT=ierr) ogelab(noge), ogeaux(noge,1)
          IF (ierr/=0) THEN
           WRITE(6,*) '--error-- in the routine rdoutgene.'
           WRITE(6,*) 'Cannot read label or scal. factor at line :'
           WRITE(6,*) TRIM(line)
           STOP "in rdoutgene"
          ENDIF
        ELSE 
          READ(line3,*,IOSTAT=ierr) ogelab(noge)
          IF (ierr/=0) THEN
           WRITE(6,*) '--error-- in the routine rdoutgene.'
           WRITE(6,*) 'Cannot read label at line :'
           WRITE(6,*) TRIM(line)
           STOP "in rdoutgene"
          ENDIF
        ENDIF
      ENDIF

! get auxiliary information - for fall off reaction only
      IF (ifo/=0) THEN
        READ(tfu1,'(a)',IOSTAT=ierr) line
        IF (ierr/=0) THEN
          WRITE(6,*) '--error-- in rdoutgene. Cannot read line :',TRIM(line)
          STOP "in rdoutgene"
        ENDIF
        READ(line,*,IOSTAT=ierr) (ogeaux(noge,i),i=1,7)
        IF (ierr/=0) THEN
          WRITE(6,*) '--error-- in the routine rdoutgene.'
          WRITE(6,*) 'Cannot read the 7 coef. at line: '
          WRITE(6,*) TRIM(line)
          STOP "in rdoutgene"
        ENDIF
      ENDIF

! GET THE PRODUCTS AND STOICHIOMETRIC COEFFICIENTS
! ------------------------------------------------
      nlin=nlin-1
      IF (ifo/=0) nlin=nlin-1

      DO j=1,nlin
        READ(tfu1,'(a)',IOSTAT=ierr) line
        IF (ierr/=0) THEN
          WRITE(6,*) '--error-- in rdoutgene. Cannot read line :',TRIM(line)
          STOP "in rdoutgene"
        ENDIF
        n1=INDEX(line,' ')
        chem1=line(1:n1-1)

! check the product (#, C1 or inorganic species)
        CALL chcksp(chem1)
        
! set product and get the stoechiometric coefficient
        ogeprod(noge,j)=chem1
        READ(line(n1:),*,IOSTAT=ierr) ogestoe(noge,j)
        IF (ierr/=0) THEN
          WRITE(6,*) '--error-- in the routine rdoutgene.'
          WRITE(6,*) 'Cannot read stoechiometric coefficient at line:'
          WRITE(6,*) TRIM(line)
          STOP "in rdoutgene"
        ENDIF

      ENDDO

! error - line does not start with *, END or REAC
    ELSE
      WRITE (6,*) '--error--, in rdoutgene. keyword unknown at line :'
      WRITE (6,*) TRIM(line)
      STOP "in rdoutgene"
    ENDIF

! read next reaction
  ENDDO readloop
  
  CLOSE(tfu1)

END SUBROUTINE rdoutgene

! ======================================================================
! PURPOSE: perform the reactions that cannot be managed by the
! generator for "#xxx" species.     
! ======================================================================
SUBROUTINE spreact(idnam,chem,brch)
  USE keyparameter, ONLY: mxlco,mxnr,mxpd,mecu,waru,kohu,kno3u
  USE keyflag, ONLY: TK
  USE references, ONLY:mxlcod
  USE database, ONLY: noge,ogelab,ogertve,ogeprod,ogearh,ogestoe,ogeaux ! out generator database 
  USE dictstacktool, ONLY: bratio
  USE rxwrttool, ONLY:rxwrit,rxinit
  USE toolbox, ONLY: stoperr,add1tonp,addrx,addref
  IMPLICIT NONE

  CHARACTER(LEN=*),INTENT(IN) :: idnam  ! name of the species     
  CHARACTER(LEN=*),INTENT(IN) :: chem   ! formula     
  REAL,INTENT(IN) :: brch               ! yield of the species
    
  CHARACTER(LEN=mxlco) :: r(3), p(mxpd)
  INTEGER  :: i,j, nr, nhv, nex, ipos, nfo
  REAL     :: s(mxpd),arrh(3)
  REAL     :: brtio
  INTEGER  :: idreac, nlabel
  REAL     :: xlabel,folow(3),fotroe(4)

  INTEGER,PARAMETER     :: mxcom=10
  INTEGER               :: nref
  CHARACTER(LEN=mxlcod) :: ref(mxcom)
  CHARACTER(LEN=7),PARAMETER :: progname='spreact'
  CHARACTER(LEN=80) :: mesg
  INTEGER :: idummy(mxnr)
! net rxn rates @ 298K, TK with OH, NO3
  REAL :: kOH_298, kOH_T, kNO3_298, kNO3_T 
  LOGICAL :: kOHfg,kNO3fg

  nr=0 ; idummy(:)=0 
  kOH_298 = 0. ; kOH_T = 0. ;  kNO3_298 = 0. ; kNO3_T = 0.
  kOHfg = .FALSE. ; kNO3fg = .FALSE.

! find reaction 
  DO i=1,noge
    DO ipos=1,2
      IF (chem==ogertve(i,ipos)) THEN
        CALL addrx(progname,chem,nr,idummy)
        CALL rxinit(r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe)
        nref=0 ; ref(:)=' ' ; CALL addref(progname,'MCMRX',nref,ref)

! set arrhenius coefficient
        arrh(1)=ogearh(i,1) ; arrh(2)=ogearh(i,2) ; arrh(3)=ogearh(i,3)
        brtio=brch

! set the reactants of the reaction
        DO j=1,3
          IF (j==ipos) THEN
            r(j)=idnam
          ELSE IF (ogertve(i,j)(1:1)/=' ') THEN
            CALL bratio(ogertve(i,j),brtio,r(j),nref,ref)
          ENDIF             
          IF (TRIM(r(j)).EQ."HO") kOHfg = .TRUE.
          IF (TRIM(r(j)).EQ."NO3") kNO3fg = .TRUE.
        ENDDO
    
! if reactant = OH or NO3, accumulate k(298) and k(TK)
! (do not exceed collision rate for individual ks)
        IF(kOHfg)THEN
          kOH_298 = kOH_298 + amin1(2.0E-10, &
                 arrh(1) * 298.**arrh(2)*EXP(-arrh(3)/298.))
          kOH_T = kOH_T + amin1(2.0E-10, &
                 arrh(1) * TK**arrh(2)*EXP(-arrh(3)/TK))
        ENDIF
        IF(kNO3fg)THEN
          kNO3_298 = kNO3_298 + amin1(2.0E-10, &
                 arrh(1) * 298.**arrh(2)*EXP(-arrh(3)/298.))
          kNO3_T = kNO3_T + amin1(2.0E-10, &
                 arrh(1) * TK**arrh(2)*EXP(-arrh(3)/TK))
        ENDIF

! set the products of the reaction
        DO j=1,mxpd
          IF (ogeprod(i,j)(1:1)/=' ') THEN
            CALL bratio(ogeprod(i,j),brtio,p(j),nref,ref)
            s(j)=ogestoe(i,j)
          ENDIF             
        ENDDO

! set id and label for the reaction (if necessary). For HV reaction,
! weighting factor is set as the first value of ogeaux (ogeaux(i,1)).
        idreac=0 ; nlabel=0 ; nhv=0 ; nex=0 ; nfo=0
    
        DO j=1,3
          IF (ogertve(i,j)(1:3)=='HV ')    nhv=nhv+1
          IF (ogertve(i,j)(1:6)=='EXTRA ') nex=nex+1
          IF (ogertve(i,j)(1:6)=='(+M) ')  nfo=nfo+1
        ENDDO
    
        IF (nhv+nex+nfo>1) THEN
          mesg='Too many keywords.'
          CALL stoperr(progname,mesg,chem)
        ENDIF
    
! photolysis
        IF (nhv==1) THEN
          idreac=1  ;  xlabel=ogeaux(i,1)  ;  nlabel=ogelab(i)
        ENDIF

! extra
        IF (nex==1) THEN
          idreac=2  ;  nlabel=ogelab(i)
        ENDIF
   
! Fall off reaction
        IF (nfo==1) THEN
          idreac=3 ; folow(1:3)=ogeaux(i,1:3) ; fotroe(1:4)=ogeaux(i,4:7)
        ENDIF

! write the reaction
        CALL rxwrit(mecu,r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe,ref)
      ENDIF
    ENDDO
  ENDDO

! if kOH or kNO3 found, write the net rates
  IF(kOHfg)THEN
      WRITE(kohu,'(a,2x,2(1pe10.2),a)') idnam, kOH_298,kOH_T," HAND, HAND"
  ENDIF
  IF(kNO3fg)THEN
      WRITE(kno3u,'(a,2x,2(1pe10.2),a)') idnam, kNO3_298,kNO3_T," HAND, HAND"
  ENDIF

  IF (nr==0) THEN
    WRITE (waru,*) '--WARNING--, the following species has no sink:'
    WRITE (waru,*) TRIM(chem)
  ENDIF

END SUBROUTINE spreact

!=======================================================================
! PURPOSE: check if the species makes sense in the generator. The 
! program stops if the species is not OK. 
!=======================================================================
SUBROUTINE chcksp(chem)
  USE atomtool, ONLY: cnum
  USE searching, ONLY: srch
  USE normchem, ONLY: stdchm
  USE database, ONLY:  nspsp, dictsp  ! special species input
  USE dictstackdb, ONLY: nrec,ninorg,inorglst,dict
  USE toolbox, ONLY: stoperr
  IMPLICIT NONE

  CHARACTER(LEN=*),INTENT(INOUT):: chem       ! the formula to be checked

! internal 
  INTEGER         :: i, nca, dicptr
! These c1 are not in dic_c1.dat but are added if necessary 
! during generation, in dictstacktool.f90  
  CHARACTER(LEN=15),DIMENSION(22) :: optional_c1=(/'CH2.(OO.)      ', &
               'CH3(OOOH)      ','CH2(OH)(OOH)   ','CH2(OH)(OOOH)  ', &
               'CHO(NO2)       ','CH2(ONO2)(OO.) ','CO(NO2)(ONO2)  ', &
               'CO(OH)(NO2)    ','CO(OH)(ONO2)   ','CO(OH)(OOH)    ', &
               'CO(ONO2)(ONO2) ','CO(ONO2)(OOH)  ','CO(OOH)(OO.)   ', &
               'CO(NO2)(OOOH)  ','CO(OOOH)(OOOH) ','CO(OH)(OH)     ', &
               'CHO(OOOH)      ','CH2(OH)(ONO2)  ','CH2(ONO2)(OOH) ', &
               'CH2(NO2)(OO.)  ','CHO(ONO2)      ','CH3(ONO)       '/)

  CHARACTER(LEN=7),PARAMETER :: progname='chcksp '
  CHARACTER(LEN=70)          :: mesg

! check if the species is a keyword
  IF (chem(1:6)=='EXTRA ') RETURN
  IF (chem(1:3)=='HV ') RETURN
  IF (chem(1:4)=='PERO') RETURN
  IF (chem(1:7)=='NOTHING') RETURN
  IF (chem(1:7)=='OXYGEN ') RETURN
  IF (chem(1:7)=='ISOM   ') RETURN

! check if the species is a special species
  IF (chem(1:1)=='#') THEN
     DO i=1,nspsp
       IF (dictsp(i)(10:131)==chem) RETURN
     ENDDO
     mesg="Species starting with # not in 'dictsp' table (see 'special_dic.dat')"
     CALL stoperr(progname,mesg,chem)

! check if the species start with a C. If the species is a C1, then must  
! already be known, otherwise check the name using stdchm
  ELSE IF ((chem(1:1)=='C').OR.(chem(1:1)=='c').OR. &
           (chem(1:4)=='=Cd1').OR.(chem(1:4)=='-O1-')) THEN
    nca=cnum(chem)
    IF (nca==1) THEN
      dicptr = srch(nrec,chem,dict)
      IF (dicptr>0) RETURN  ! known in dic_c1.dat
      DO i=1,SIZE(optional_c1) 
        IF (TRIM(chem)==TRIM(optional_c1(i))) RETURN ! known in dictstacktool.f90
      ENDDO
      mesg="C1 species not found in dictionary (see 'dic_c1.dat') or in dictstacktool.f90"
      CALL stoperr(progname,mesg,chem)
    ELSE
      CALL stdchm(chem)
    ENDIF

! If the test above failed, then the species must be inorganic
  ELSE
    DO i=1,ninorg
      IF (chem==inorglst(i)(10:129)) RETURN
    ENDDO
    WRITE(6,*)'ninorg=',ninorg
    mesg="Species not found in the inorganic dictionary"
    CALL stoperr(progname,mesg,chem)
  ENDIF

END SUBROUTINE chcksp


END MODULE spsptool
