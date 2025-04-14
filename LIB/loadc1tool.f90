MODULE loadc1tool
IMPLICIT NONE
CONTAINS

!SUBROUTINE loadc1mch()
!SUBROUTINE rddict(filename,ninorg,iname,ichem,ifgrp,nrec,oname,ochem,ofgrp)
!SUBROUTINE rdfixmch(filename)
!SUBROUTINE rxmeo2ro2()
!SUBROUTINE getcomcod(filename,comline,comtab)

!=======================================================================
! PURPOSE: load the various dictionaries and mechanisms for inorganic
! species, C1 and any particular "external" mechanisms. 
!=======================================================================
SUBROUTINE loadc1mch()
  USE keyflag, ONLY: pamfg
  USE keyparameter, ONLY: mxlco,mxlfo,mxlfl,dirgecko
  USE dictstackdb, ONLY: nrec,ninorg,inorglst,dict,namlst
  USE sortstring, ONLY: sort_string ! to sort chemical
  IMPLICIT NONE

  INTEGER,PARAMETER :: mxtprec=200     ! max size of the inorg/C1 dictionary
  CHARACTER(LEN=mxlco) :: oname(mxtprec)
  CHARACTER(LEN=mxlfo) :: ochem(mxtprec)
  CHARACTER(LEN=mxlfl) :: ofgrp(mxtprec) 
  CHARACTER(LEN=mxlco) :: iname(mxtprec)
  CHARACTER(LEN=mxlfo) :: ichem(mxtprec)
  CHARACTER(LEN=mxlfl) :: ifgrp(mxtprec) 
  CHARACTER(LEN=256) :: filename
  INTEGER :: i  

  oname(:)=' ' ; ochem(:)=' '; ofgrp(:)=' '
  iname(:)=' ' ; ichem(:)=' '; ifgrp(:)=' '
  nrec = 1 ; oname(1)='######'  ; ochem(1)='###################'

! read dictionaries
! -----------------
  filename=TRIM(dirgecko)//"DATA/dic_inorg.dat"
  CALL rddict(filename,ninorg,iname,ichem,ifgrp,nrec,oname,ochem,ofgrp)

  filename=TRIM(dirgecko)//"DATA/dic_c1.dat"
  CALL rddict(filename,ninorg,iname,ichem,ifgrp,nrec,oname,ochem,ofgrp)

  IF (pamfg) THEN
    filename=TRIM(dirgecko)//"DATA/dic_OFR.dat"
    CALL rddict(filename,ninorg,iname,ichem,ifgrp,nrec,oname,ochem,ofgrp)
  ENDIF

  IF (nrec>SIZE(dict)) THEN
    WRITE (6,*) 'error in loadc1mch, number of records exceeds mni'
    STOP "in rddict"
  ENDIF  
  IF (ninorg>SIZE(inorglst)) THEN
    WRITE (6,*) 'error in loadc1mch, number of records exceeds mni'
    STOP "in loadc1mch"
  ENDIF  

  DO i=1,ninorg
    WRITE(inorglst(i),'(a6,3x,a120,2x,a15)') iname(i), ichem(i), ifgrp(i)
  ENDDO

  DO i=1,nrec
    WRITE(dict(i),'(a6,3x,a120,2x,a15)') oname(i), ochem(i), ofgrp(i)
  ENDDO  

! set the sorted for shortnames 
  namlst(1:nrec)=oname(1:nrec)
  CALL sort_string(namlst(1:nrec))
  
! read the mechanisms
! ------------------- 
  WRITE (6,*) '  ...reading inorganic reactions'
  filename=TRIM(dirgecko)//"DATA/mch_inorg.dat"
  CALL rdfixmch(filename)

  IF (pamfg) THEN
    filename=TRIM(dirgecko)//'DATA/mch_OFR.dat'
    CALL rdfixmch(filename)
  ENDIF   
  
  WRITE (6,*) '  ...reading CH4 chemistry'
  filename=TRIM(dirgecko)//"DATA/mch_singlec.dat"
  CALL rdfixmch(filename)

! write CH3O2+counters reactions
  WRITE (6,*) '  ...writing CH3O2+counters reactions'
  CALL rxmeo2ro2

END SUBROUTINE loadc1mch

!=======================================================================
! PURPOSE: read the species that are not "self generated" (e.g. 
! inorganic, C1 species) and add the species to the tables related
! to the dictionary (either inorg or org). Species are sorted according
! to formula.                                          
!=======================================================================
SUBROUTINE rddict(filename,ninorg,iname,ichem,ifgrp,nrec,oname,ochem,ofgrp)
  USE keyparameter, ONLY: tfu1
  USE atomtool, ONLY: cnum
  USE sortstring, ONLY: sort_string ! to sort chemical
  USE searching, ONLY: srh5   ! to seach in a sorted list
  IMPLICIT NONE

  CHARACTER(LEN=*),INTENT(IN)    :: filename
  INTEGER,INTENT(INOUT)          :: ninorg   ! # of inorganic species 
  INTEGER,INTENT(INOUT)          :: nrec     ! # of species (carbonaceous)
  CHARACTER(LEN=*),INTENT(INOUT) :: oname(:) ! name of the org. species
  CHARACTER(LEN=*),INTENT(INOUT) :: ochem(:) ! formula of the org. species
  CHARACTER(LEN=*),INTENT(INOUT) :: ofgrp(:) ! fun. grp. info of the org. species
  CHARACTER(LEN=*),INTENT(INOUT) :: iname(:) ! name of the inorg. species
  CHARACTER(LEN=*),INTENT(INOUT) :: ichem(:) ! formula of the inorg. species
  CHARACTER(LEN=*),INTENT(INOUT) :: ifgrp(:) ! fun. grp. the inorg. species

  INTEGER,PARAMETER :: lenlin=200
  CHARACTER(LEN=lenlin) :: line, tpline
  INTEGER  :: i,ilin,ierr,ipos,loerr,nca
  CHARACTER(LEN=LEN(oname(1))) :: tpname
  CHARACTER(LEN=LEN(ochem(1))) :: tpchem
  CHARACTER(LEN=LEN(ofgrp(1))) :: tpfgrp
  CHARACTER(LEN=LEN(oname(1))) :: tponame(SIZE(oname))
  CHARACTER(LEN=LEN(ochem(1))) :: tpochem(SIZE(ochem))
  CHARACTER(LEN=LEN(ofgrp(1))) :: tpofgrp(SIZE(ofgrp))
  CHARACTER(LEN=LEN(iname(1))) :: tpiname(SIZE(iname))
  CHARACTER(LEN=LEN(ichem(1))) :: tpichem(SIZE(ichem))
  CHARACTER(LEN=LEN(ifgrp(1))) :: tpifgrp(SIZE(ifgrp))

  tponame(:)=oname(:) ; tpochem(:)=ochem(:) ; tpofgrp(:)=ofgrp(:)
  tpiname(:)=iname(:) ; tpichem(:)=ichem(:) ; tpifgrp(:)=ifgrp(:)

! --------------------------------------------------------
! read dictionary - store the data in temporary tables (tp)
! -------------------------------------------------------
  OPEN(tfu1,FILE=filename,FORM='FORMATTED', STATUS='OLD')

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

! read the line
    loerr=0
    tpline=ADJUSTL(line)
    READ(tpline,*,IOSTAT=ierr) tpname
    IF (ierr/=0) loerr=1
    ipos=INDEX(tpline,' ') ; tpline=tpline(ipos+1:) ; tpline=ADJUSTL(tpline)
    READ(tpline,*,IOSTAT=ierr) tpchem
    IF (ierr/=0) loerr=1
    ipos=INDEX(tpline,' ') ; tpline=tpline(ipos+1:) 
    IF (tpline/=' ') THEN 
      READ(tpline,*,IOSTAT=ierr) tpfgrp
      IF (ierr/=0) loerr=1
    ELSE
      tpfgrp=' '
    ENDIF
    IF (loerr/=0) THEN
      WRITE(6,*) '--error--, while reading file: ',TRIM(filename)
      WRITE(6,*) 'at line number: ',ilin
      WRITE(6,*) 'can not read parameters in :', TRIM(line)   
      STOP "in rddict"
    ENDIF     

! get the number of C and check that it is not greater than 1
    nca = cnum(tpchem)
    IF (nca>1) THEN
      WRITE (6,*) '--error--, while reading file: ',TRIM(filename)
      WRITE (6,*) 'following species have more than 1 C:'
      WRITE (6,*) TRIM(line)
      STOP "in rddict"
    ENDIF

! if C1 species, store data in the "organic" tables
    IF (nca==1) THEN
      IF (tpchem(1:1)/='C') THEN
        WRITE (6,*) '-error-, while reading file: ',TRIM(filename)
        WRITE (6,*) 'organic species does not start with C:', tpchem
        WRITE (6,*) TRIM(line)
        STOP "in rddict"
      ENDIF
      nrec = nrec + 1
      tponame(nrec)=tpname ; tpochem(nrec)=tpchem ; tpofgrp(nrec)=tpfgrp
      
! else inorganic species - store data the "inorganic" tables
    ELSE
      ninorg = ninorg + 1
      tpiname(ninorg)=tpname ; tpichem(ninorg)=tpchem ; tpifgrp(ninorg)=tpfgrp
    ENDIF
  ENDDO rdloop  
  CLOSE(tfu1)

! check the size
  IF ((nrec>SIZE(oname)) .OR. (ninorg>SIZE(iname))) THEN
    WRITE (6,*) '--error--, while reading file: ',TRIM(filename)
    WRITE (6,*) 'too many species in the input file (nrec,ninorg > mxtprec)'
    STOP "in rddict"
  ENDIF

! -------------------------------------
! sort the formula of carbonaceous species
! -------------------------------------
  ochem(:)=tpochem(:)
  CALL sort_string(ochem(1:nrec)) 

! check for duplicate 
  ierr=0
  DO i=1,nrec-1
    IF (ochem(i)==ochem(i+1)) THEN
      WRITE(6,*) '--error--, while reading file: ',TRIM(filename)
      WRITE(6,*)  'Following species identified 2 times: ',TRIM(ochem(i))
      ierr=1
    ENDIF
  ENDDO
  IF (ierr/=0) STOP "in rddict" 

! sort the tables according to formula
  DO i=1,nrec
    ilin=srh5(tpochem(i),ochem,nrec)
    IF (ilin <= 0) THEN
      WRITE(6,*) '--error--, while sorting species in: ',TRIM(filename)
      WRITE(6,*) 'species "lost" after sorting the list: ', TRIM(tpochem(i))
      STOP "in rddict" 
    ENDIF
    oname(ilin)=tponame(i)  ;  ofgrp(ilin)=tpofgrp(i)
  ENDDO

! -------------------------------------
! sort the formula of inorganic species
! -------------------------------------
  ichem(:)=tpichem(:)
  CALL sort_string(ichem(1:ninorg)) 

! check for duplicate 
  ierr=0
  DO i=1,ninorg-1
    IF (ichem(i)==ichem(i+1)) THEN
      WRITE(6,*) '--error--, while reading file: ',TRIM(filename)
      WRITE(6,*)  'Following species identified 2 times: ',TRIM(ichem(i))
      ierr=1
    ENDIF
  ENDDO
  IF (ierr/=0) STOP "in rddict" 

! sort the tables according to formula
  DO i=1,ninorg
    ilin=srh5(tpichem(i),ichem,ninorg)
    IF (ilin <= 0) THEN
      WRITE(6,*) '--error--, while sorting species in: ',TRIM(filename)
      WRITE(6,*) 'species "lost" after sorting the list: ', TRIM(tpichem(i))
      STOP "in rddict" 
    ENDIF
    iname(ilin)=tpiname(i)  ;  ifgrp(ilin)=tpifgrp(i)
  ENDDO

END SUBROUTINE rddict

!=======================================================================
! PURPOSE: read fixed chemical scheme (e.g. inorganic & C1 chemistry).      
! The reactions in "filename" are just read and copied in mechanism 
! given in the output (unit lout)  
! 
! NOTE: Species names in filename may not exist in the dictionary.
! This point need to be adressed in future versions.
!=======================================================================
SUBROUTINE rdfixmch(filename)
  USE keyparameter, ONLY: tfu1,mecu,refu
  USE keyflag , ONLY: wrtref
  USE references, ONLY: mxlcod   ! max length of comment's code in database
  USE rxwrttool, ONLY: count4rxn
  IMPLICIT NONE
      
  CHARACTER(LEN=*),INTENT(IN) :: filename

  INTEGER,PARAMETER :: lenlin=100
  CHARACTER(LEN=lenlin) :: line, line1, line2
  CHARACTER(LEN=lenlin) :: reaction,info,comline
  CHARACTER(LEN=mxlcod) :: comtab(3)
  INTEGER :: label,ilin
  REAL    :: A, n, E, m, fact
  REAL    :: one, zero
  REAL    :: F0_300, Finf_300,E0,Einf
  REAL    :: Fc1, Fc2, Fc3, Fc4
  INTEGER :: n1, n2, loerr, ierr,idrx
              
  one=1.  ;  zero=0.  ;  loerr=0

! open the file
  OPEN(tfu1,FILE=filename,STATUS='OLD',FORM='FORMATTED')

! read next line
  ilin=0
  rdloop: DO
    ilin=ilin+1
    READ (tfu1,'(a)',IOSTAT=ierr) line
    IF (ierr/=0) THEN
      WRITE(6,*) '--error--, while reading file: ',TRIM(filename)
      WRITE(6,*) 'at line number: ',ilin
      WRITE(6,*) 'keyword "END" missing ?'
      STOP "in rdfixmch"
    ENDIF     
    IF (line(1:3)=='END') EXIT rdloop
    IF (line(1:1)=='!') CYCLE rdloop

! remove anything after "!" (not necessary as 1st character)
    n1=INDEX(line,"!")
    IF (n1>0) line(n1:)=' '

! Any new reaction must have ";" and ":" character
    n1 = INDEX(line,':')  ;  n2 = INDEX(line,';')
    IF (n1==0) loerr=1    ;  IF (n2==0) loerr=1
    IF (loerr/=0) THEN
      WRITE(6,*) '--error--, missing ":" or ";" while reading: ', TRIM(filename)
      WRITE(6,*) 'at line number: ',ilin
      WRITE(6,*) TRIM(line)
      STOP "in rdfixmch"
    ENDIF
 
! check fall off reaction
! -----------------------
    IF (INDEX(line,'FALLOFF')/=0) THEN
! check that the reaction is correctly formatted
      READ(line(n1+1:n2-1),*,IOSTAT=ierr)  Fc1, Fc2, Fc3, Fc4
      IF (ierr/=0) THEN
        WRITE(6,*) '--error--, while reading in ', TRIM(filename)
        WRITE(6,*) 'at line number: ',ilin
        WRITE(6,*) 'while reading # in: ',TRIM(line)
        STOP "in rdfixmch"
      ENDIF
      reaction = line(1:n1-1)
      comline=line(n2:)   ! define comment line with ";" as 1st char
      line(n2:n2)='/'
      info=line(n1+1:n2)
      CALL getcomcod(filename,comline,comtab)

      ilin=ilin+1
      READ(tfu1,'(a)') line1
      READ(line1,*,IOSTAT=ierr) F0_300, n, E0
      IF (ierr/=0) THEN
        WRITE(6,*) '--error--, while reading in ', TRIM(filename)
        WRITE(6,*) 'at line number: ',ilin
        WRITE(6,*) 'while reading # in: ',TRIM(line1)
        STOP "in rdfixmch"
      ENDIF

      ilin=ilin+1
      READ(tfu1,'(a)') line2
      READ(line2,*,IOSTAT=ierr) Finf_300, m, Einf
      IF (ierr/=0) THEN
        WRITE(6,*) '--error--, while reading:', TRIM(filename)
        WRITE(6,*) 'at line number: ', ilin
        WRITE(6,*) 'while reading # in: ',TRIM(line2)
        STOP "in rdfixmch"
      ENDIF

! write reaction to the output file
      line1=' '
      WRITE(line1,'(1X,E10.3,1X,f4.1,1X,f7.0)') F0_300,-n,E0
      line='  FALLOFF /'//TRIM(line1)//TRIM(info)

      CALL count4rxn(3)   ! index 3 is falloff type
      WRITE(mecu,'(A70,19X,ES10.3,1X,f4.1,1X,f7.0)') reaction,Finf_300,-m,Einf
      WRITE(mecu,'(a)') TRIM(line)
      IF (wrtref) THEN
        WRITE(refu,'(A70,19X,ES10.3,1X,f4.1,1X,f7.0)') reaction,Finf_300,-m,Einf
        WRITE(refu,'(a)') TRIM(line)
      ENDIF
 
! check extra reaction
! --------------------
    ELSE IF (INDEX(line,'EXTRA')/=0) THEN
        
! check that the reaction is correctly formatted
      reaction = line(1:n1-1)
      READ(line(n1+1:n2-1),*,IOSTAT=ierr)  A, n, E
      IF (ierr/=0) THEN
        WRITE(6,*) '--error--, while reading:', TRIM(filename)
        WRITE(6,*) 'at line number: ',ilin
        WRITE(6,*) 'while reading # in: ',TRIM(line)
        STOP "in rdfixmch"
      ENDIF
      comline=line(n2:)
      CALL getcomcod(filename,comline,comtab)

      ilin=ilin+1
      READ(tfu1,'(a)') line1
      n1 = INDEX(line1,':')  ;  n2 = INDEX(line1,';')
      IF (n1==0) loerr=1     ;  IF (n2==0) loerr=1
      IF (loerr/=0) THEN
        WRITE(6,*) '--error--, missing ":" or ";" while reading: ', TRIM(filename)
        WRITE(6,*) TRIM(line1)
        STOP "in rdfixmch"
      ENDIF
      READ (line1(1:n1-1),*,IOSTAT=ierr) label
      IF (ierr/=0) THEN
        WRITE(6,*) '--error--, while reading:', TRIM(filename)
        WRITE(6,*) 'at line number: ',ilin
        WRITE(6,*) 'while reading label in: ',TRIM(line1)
        STOP "in rdfixmch"
      ENDIF
      line1(n2:n2)='/'
      info=' '
      info=line1(n1+1:n2)

! write reaction to the output file
      CALL count4rxn(2) ! index 2 is extra type
      WRITE(mecu,'(A70,19X,ES10.3,1X,f4.1,1X,f7.0)') reaction,A,n,E
      WRITE(mecu,'(2X,A7,i4,1X,A50)') 'EXTRA /',label,info
      IF (wrtref) THEN
        WRITE(refu,'(A70,19X,ES10.3,1X,f4.1,1X,f7.0)') reaction,A,n,E
        WRITE(refu,'(2X,A7,i4,1X,A50)') 'EXTRA /',label,info
      ENDIF

! check photolytic reaction
! -------------------------
    ELSE IF (INDEX(line,'HV') /= 0) THEN

! check that the reaction is correctly formatted
      reaction = line(1:n1-1)
      READ(line(n1+1:n2-1),*,IOSTAT=ierr)  label, fact
      IF (ierr/=0) THEN
        WRITE(6,*) '--error--, while reading:', TRIM(filename)
        WRITE(6,*) 'while reading label & fact in: ',TRIM(line)
        STOP "in rdfixmch"
      ENDIF
      comline=line(n2:)
      CALL getcomcod(filename,comline,comtab)

! write reaction to the output file
      CALL count4rxn(1) ! index 1 is hv type
      WRITE(mecu,'(A70,19X,ES10.3,1X,f4.1,1X,f7.0)') reaction,one,zero,zero
      WRITE(mecu,'(A7,i5,1x,f5.2,A2)') '  HV / ',label,fact,' /'
      IF (wrtref) THEN
        WRITE(refu,'(A70,19X,ES10.3,1X,f4.1,1X,f7.0)') reaction,one,zero,zero
        WRITE(refu,'(A7,i5,1x,f5.2,A2)') '  HV / ',label,fact,' /'
      ENDIF
            
! otherwise thermal reaction
! ---------------------------
    ELSE 
      reaction = line(1:n1-1)
      READ(line(n1+1:n2-1),*,IOSTAT=ierr) A,n,E
      IF (ierr/=0) THEN
        WRITE(6,*) '--error--, while reading:', TRIM(filename)
        WRITE(6,*) 'while reading # in: ',TRIM(line)
        STOP "in rdfixmch"
      ENDIF
      comline=line(n2:)
      CALL getcomcod(filename,comline,comtab)

! write reaction to the output file
      idrx=0                                  ! default label (0 is simple thermal rxn)
      IF (INDEX(reaction,'TBODY')/=0) idrx= 4 ! label for third body reaction
      CALL count4rxn(idrx)   
      WRITE(mecu,'(A70,19X,ES10.3,1X,f4.1,1X,f7.0)') reaction,A,n,E
      IF (wrtref) WRITE(refu,'(A70,19X,ES10.3,1X,f4.1,1X,f7.0)') reaction,A,n,E
    ENDIF

  ENDDO rdloop
  CLOSE(tfu1)

END SUBROUTINE rdfixmch

!=======================================================================
! PURPOSE: write the CH3O2+RO2 reactions in the the output file.
! The reactions are just written in the mechanism output file
!                 
! NOTE: Data in this file must be consistent with those in ro2tool and    
! rco3tool
!=======================================================================
SUBROUTINE rxmeo2ro2()
  USE keyparameter, ONLY: mxlco, mxpd,mecu  ! max name length & prod. per reaction
  USE rxwrttool, ONLY:rxwrit,rxinit
  USE keyflag, ONLY: rx_ro2_multiclass      ! If true, consider all classes of RO2
  IMPLICIT NONE

  CHARACTER(LEN=mxlco) :: p(mxpd), r(3)
  REAL     :: xlabel,folow(3),fotroe(4)
  REAL     :: s(mxpd), arrh(3)
  REAL     :: rrad,rmol1,rmol2
  INTEGER  :: itype, idreac, nlabel

  CHARACTER(LEN=mxlco),DIMENSION(9), PARAMETER :: kwdclass= &
    (/'PERO1 ' , 'PERO2 ' , 'PERO3 ' , 'PERO4 ' , 'PERO5 ', &
      'MEPERO' , 'PERO7 ' , 'PERO8 ' , 'PERO9 '/)

  REAL,DIMENSION(9,3),PARAMETER :: ro2dat=RESHAPE( &  ! reshape here for 2 dim. array
   (/ 1.00E-13,  0.,  627. ,  &  ! PERO1- unsubstituted ter RO2
      1.00E-13,  0.,  341. ,  &  ! PERO2- iC3H7O2
      1.00E-13,  0., -257. ,  &  ! PERO3- ter RO2 with alpha or beta O or N
      1.00E-13,  0., -289. ,  &  ! PERO4- C2H5O2 & unsubstituted sec RO2
      1.00E-13,  0., -635. ,  &  ! PERO5- subst. ter RO2 and allyl or benzyl
      2.06E-13,  0., -365. ,  &  ! MEPERO- CH3O2
      1.00E-13,  0., -702. ,  &  ! PERO7- unsubst. prim RO2 & subst. sec RO2
      1.00E-13,  0., -936. ,  &  ! PERO8- subst. prim RO2 & subst. sec RO2 with allyl
      2.00E-12,  0., -508./), &  ! PERO9- acyl RO2
      SHAPE(ro2dat), ORDER=(/2,1/) )  !- order (/2,1/) is to switch row & column from reshape

  REAL,DIMENSION(9,3),PARAMETER :: stoi_m=reshape( &  ! stoi. coef. for CH3O2+RO2
      (/ 0.8 , 0.2 , 0.0 ,  & ! PERO1- unsubstituted ter RO2
         0.6 , 0.2 , 0.2 ,  & ! PERO2- iC3H7O2
         0.8 , 0.2 , 0.0 ,  & ! PERO3- ter RO2 with alpha or beta O or N
         0.6 , 0.2 , 0.2 ,  & ! PERO4- C2H5O2 & unsubstituted sec RO2
         0.8 , 0.2 , 0.0 ,  & ! PERO5- subst. ter RO2 and allyl or benzyl
         0.370, 0.315, 0.315,  & ! MEPERO- CH3O2
         0.6 , 0.2 , 0.2 ,  & ! PERO7- unsubst. prim RO2 & subst. sec RO2
         0.6 , 0.2 , 0.2 ,  & ! PERO8- subst. prim RO2 & subst. sec RO2 with allyl
         0.8 , 0.2 , 0.0 /),& ! PERO9- acyl RO2
       SHAPE(stoi_m), ORDER=(/2,1/) )  !- order (/2,1/) is to switch row & column from reshape

  r(:)=' ' ; folow(:)=0. ; fotroe(:)=0.
  s(:)=0.  ; p(:)=' '    ; idreac=0  ;  nlabel=0

! -------------------------------------------
! ALL RO2 class: loop over the various RO2 class
! -------------------------------------------
  IF (rx_ro2_multiclass) THEN 

    ro2loop: DO itype=1,9
      rrad =stoi_m(itype,1) 
      rmol1=stoi_m(itype,2) 
      rmol2=stoi_m(itype,3)

! initialize reaction
      CALL rxinit(r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe)

! set reactants
      r(1)='CH3O2'  ;  r(2)=kwdclass(itype)

! set rate constant
      arrh(:)=ro2dat(itype,:)

      s(1)=rrad   ;  p(1) = 'CH3O '  ! radical channel
      s(2)=rmol1  ;  p(2) = 'CH2O '  ! H given by CH3O2
      s(3)=rmol2  ;  p(3) = 'CH3OH'  ! H taken by CH3O2

      CALL rxwrit(mecu,r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe)
    ENDDO ro2loop 

! -------------------------------------------
! SINGLE RO2 class only
! -------------------------------------------
  ELSE  
    itype=6
    rrad =stoi_m(itype,1) 
    rmol1=stoi_m(itype,2) 
    rmol2=stoi_m(itype,3)
      
! initialize reaction
    CALL rxinit(r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe)

! set reactants
    r(1)='CH3O2'  ;  r(2)='PERO1'

! set rate constant :
    arrh(:)=ro2dat(6,:)

    s(1)=rrad   ;  p(1) = 'CH3O '  ! radical channel
    s(2)=rmol1  ;  p(2) = 'CH2O '  ! H given by CH3O2
    s(3)=rmol2  ;  p(3) = 'CH3OH'  ! H taken by CH3O2

    CALL rxwrit(mecu,r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe)

  ENDIF

END SUBROUTINE rxmeo2ro2

!=======================================================================
! Purpose: read the line (comline) provided as input and return the
! table (comtab) with the various codes.
!=======================================================================
SUBROUTINE getcomcod(filename,comline,comtab)
  USE searching, ONLY: srh5   ! to seach in a sorted list
  USE references, ONLY: nreac_info, code  ! available code
  IMPLICIT NONE
  CHARACTER(LEN=*),INTENT(IN) :: filename  ! file currently read 
  CHARACTER(LEN=*),INTENT(IN) :: comline   ! line (piece) to be read for comment code
  CHARACTER(LEN=*),INTENT(OUT):: comtab(:) ! table of the codes
  
  INTEGER :: maxcom,ncpbeg,ncpend,ierr
  INTEGER :: i,n1
  
  maxcom=SIZE(comtab) ; ncpbeg=1 ; comtab(:)=' '

  IF (comline(1:1)/=";") THEN
    WRITE(6,*) '--error--, while reading file: ',TRIM(filename)
    WRITE(6,*) 'expected 1st char (;) not found in:',TRIM(comline)
    STOP "in getcomcod" 
  ENDIF
  
! read the comment code
  DO i=1,maxcom
    ncpend=INDEX(comline(ncpbeg+1:),';')
    IF (ncpend==0) THEN
      IF (comline(ncpbeg+1:)/=" ") THEN
        READ(comline(ncpbeg+1:),*,IOSTAT=ierr) comtab(i)
        IF (ierr/=0) THEN
          WRITE(6,*) '--error--, while reading file: ',TRIM(filename)
          WRITE(6,*) 'while reading reference code at comline: ',TRIM(comline)
          STOP "in getcomcod" 
        ENDIF
      ENDIF
      EXIT
    ELSE
      ncpend=ncpbeg+ncpend  ! count ncpend from the 1st char in comline 
      IF (comline(ncpbeg+1:ncpend-1)/=" ") THEN
        READ(comline(ncpbeg+1:ncpend-1),*,IOSTAT=ierr) comtab(i)
        IF (ierr/=0) THEN
          WRITE(6,*) '--error--, while reading file: ',TRIM(filename)
          WRITE(6,*) 'while reading reference code at comline: ',TRIM(comline)
          STOP "in getcomcod" 
        ENDIF
      ENDIF
      ncpbeg=ncpend
      CYCLE       
    ENDIF 
  ENDDO
  
! check the availability of the codes
  ierr=0
  DO i=1,maxcom
    IF (comtab(i)(1:1) /= ' ') THEN 
      n1=srh5(comtab(i),code,nreac_info)
      IF (n1 <= 0) THEN
        WRITE(6,*) '--error--, while reading file ',TRIM(filename)
        WRITE(6,*) 'unknow comment: ', comtab(i)
        ierr=1
      ENDIF
    ENDIF
  ENDDO
  IF (ierr/=0) STOP "in getcomcod" 
  
END SUBROUTINE getcomcod

END MODULE loadc1tool
