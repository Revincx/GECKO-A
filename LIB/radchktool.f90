MODULE radchktool
  IMPLICIT NONE
  CONTAINS
!SUBROUTINE radchk(chem,prod,coprod,nip,sc,nref,com)
!SUBROUTINE multip(inchem,prod,coprod,np,nref,com)
!SUBROUTINE single(tchem,prod,coprod)
!SUBROUTINE deloc(chem,tchem,nchem,sc)
!SUBROUTINE ooalkdec(chem,tchem,nchem,sc,dcom)
!SUBROUTINE epoxalkdec(inchem,outchem,ecom)
!=======================================================================
! PURPOSE: Drives radical rearrangement. First, check if electron can
! happen. For each mesomere form, then check the formula: if the number 
! of carbons is 1, then calls "single" otherwise calls "multip". 
!=======================================================================
SUBROUTINE radchk(chem,prod,coprod,nip,sc,nref,com)
  USE references, ONLY:mxlcod
  USE cdtool, ONLY: alkcheck
  USE atomtool, ONLY: cnum
  USE normchem, ONLY: stdchm
  USE sortstring, ONLY:sort_string
  USE toolbox, ONLY: stoperr
  IMPLICIT NONE

  CHARACTER(LEN=*),INTENT(IN) :: chem         ! radical species to be checked
  CHARACTER(LEN=*),INTENT(OUT):: prod(:)      ! new forms of chem (after e- deloc and rearrangt)
  CHARACTER(LEN=*),INTENT(OUT):: coprod(:,:)  ! coproduct list associated with each form of chem
  INTEGER,INTENT(OUT)   :: nip                ! # of possible formula (if e- delocalisation)
  REAL,INTENT(OUT)      :: sc(:)              ! "stoi. coef" of the formula
  INTEGER,INTENT(INOUT) :: nref               ! # of used ref in com table
  CHARACTER(LEN=*),INTENT(INOUT):: com(:)     ! ref/comment for the transformation operated
      
  CHARACTER(LEN=LEN(chem)) :: tchem(SIZE(prod)),dchem
  CHARACTER(LEN=LEN(coprod(1,1))) :: mcopd(SIZE(coprod,2))  ! co-products of multip
  CHARACTER(LEN=LEN(coprod(1,1))) :: scopd(4)               ! co-products of single 
  CHARACTER(LEN=LEN(coprod(1,1))) :: tcopd
  INTEGER :: nc,nca,i,ip
  INTEGER :: npd,nmco
  CHARACTER(LEN=mxlcod) :: acom,dcom,ecom                   ! comment code 
  CHARACTER(LEN=9),PARAMETER :: progname='radchk '
  CHARACTER(LEN=70)          :: mesg

! return if no species:
  nc = INDEX(chem,' ') - 1 ; IF (nc<1) RETURN

! stop if no radical:
  IF (INDEX(chem,'.') == 0) THEN
    mesg="non radical species as input "
    CALL stoperr(progname,mesg,chem)
  ENDIF

! initialize:
  tchem(:) = ' '    ;  prod(:) = ' '  ;  mcopd(:) = ' ' ;  sc(:) = 0
  coprod(:,:) = ' ' ;  scopd(:) = ' '
  nip = 1   ;   sc(1) = 1.  ;   tchem(1) = chem
  
! check for electron delocalisation
  IF (INDEX(chem,'=')/=0)  THEN
    CALL deloc(chem,tchem,nip,sc)
  ENDIF

! check for -O--O-C. structures
  IF ((nip==1).AND.(INDEX(tchem(1),'-O-')/=0)) THEN
    dchem=tchem(1)
    CALL ooalkdec(dchem,tchem,nip,sc,dcom)
    IF (dcom/=' ') THEN ; nref=nref+1 ; com(nref)=dcom ; ENDIF
  ENDIF  
  
  DO i=1,nip
    dchem=tchem(i)
    CALL epoxalkdec(dchem,tchem(i),ecom)
    IF (ecom/=' ') THEN ; nref=nref+1 ; com(nref)=ecom ; ENDIF
  ENDDO
  
! loop on the different products formed from delocalisation      
  DO ip=1,nip
    nca=cnum(tchem(ip))  ;  npd=0
    IF (INDEX(tchem(ip),'.')==0) CYCLE

! if multiple carbons, call multip:
! ---------------------------------
    IF (nca > 1) THEN
      CALL multip(tchem(ip),prod(ip),mcopd,nmco,nref,com)
      npd=npd+nmco
! check if product is a substituted alkene and rearrange if necessary
      IF ((INDEX(prod(ip),'Cd(O')/=0).OR.(INDEX(prod(ip),'CdH(O')/=0)) THEN
        CALL alkcheck(prod(ip),tcopd,acom)
        IF (tcopd/=' ') THEN ;  npd=npd+1 ; mcopd(npd)=tcopd ; ENDIF
        IF (acom/=' ') THEN ; nref=nref+1 ; com(nref)=acom ; ENDIF
      ENDIF

! update tchem and nca
      tchem(ip) = prod(ip)  ;  nca = cnum(tchem(ip))
    ENDIF

! if single carbon, call subroutine single & collect co-products
! -------------------------------------------------
    IF (nca==1) THEN
      CALL single(tchem(ip),prod(ip),scopd)

! check size & add scopd to mcopd
      DO i=1,SIZE(scopd)
        IF (scopd(i)/=' ') THEN 
          npd=npd+1
          IF (npd > SIZE(mcopd)) THEN
            mesg="Too many copdcts created for "
            CALL stoperr(progname,mesg,chem)
          ENDIF
          mcopd(npd)=scopd(i) 
        ENDIF
      ENDDO
    ENDIF                

! make list of coproduct & standardize prod
    IF (npd>0) coprod(ip,1:npd)=mcopd(1:npd)   
    CALL stdchm(prod(ip))
  ENDDO
      
! check if the two products are not the same (i.e. symmetrical products)
  IF (nip > 1) THEN
    IF (prod(1)==prod(2)) THEN
      nip = 1  ;  sc(1)=1. ;  sc(2)=0.  ;  prod(2)=' '
    ENDIF
  ENDIF
  
END SUBROUTINE radchk


! ======================================================================
! ======================================================================
! PURPOSE: check organic radical species (C>1). Perform rearrangement 
! when a unique transformation is identified (e.g. R. (+O2) -> RO2) or 
! RCO(O.) -> R. + CO2 etc ..). No e- delocalisation here (see deloc 
! subroutine).
! ======================================================================
! ======================================================================
SUBROUTINE multip(inchem,prod,coprod,np,nref,com)
  USE keyparameter, ONLY: mxnode,mxlgr,mxring  
  USE rjtool, ONLY: rjgrm
  USE atomtool, ONLY: cnum,onum
  USE stdgrbond, ONLY: grbond
  USE reactool, ONLY: swap,rebond
  USE toolbox, ONLY: stoperr
  IMPLICIT NONE

  CHARACTER(LEN=*),INTENT(IN)  :: inchem     ! species in 
  CHARACTER(LEN=*),INTENT(OUT) :: prod       ! possibly rearranged species out 
  CHARACTER(LEN=*),INTENT(OUT) :: coprod(:)  ! coproducts (if any)
  INTEGER,INTENT(OUT) :: np                  ! # of coproducts
  INTEGER,INTENT(INOUT) :: nref              ! # of used ref in com table
  CHARACTER(LEN=*),INTENT(INOUT):: com(:)    ! ref/comment for the transformation operated

  CHARACTER(LEN=LEN(inchem)) :: tchem
  CHARACTER(LEN=mxlgr)       :: tgroup(mxnode), tempkg, pold, pnew
  INTEGER :: tbond(mxnode,mxnode)
  INTEGER :: rjg(mxring,2)
  INTEGER :: nc,nca,j1,j2,k,keep
  INTEGER :: i,j,nh,dbflg,nring
  LOGICAL :: locriegee
  CHARACTER(LEN=9),PARAMETER :: progname='multip '
  CHARACTER(LEN=70)          :: mesg

! initialize:
  tchem=inchem ; prod = ' ' ; coprod(:) = ' ' ;  np = 0 

! WARNING: large "rdloop" loop (up to the end of the subroutine).  Exit loop
! when no "." remains on tchem or additional transformation is not needed.
rdloop: DO
  j1=0  ;  j2=0
  nca = cnum(tchem)+onum(tchem)
  IF (nca==1) EXIT rdloop

  CALL grbond(tchem,tgroup,tbond,dbflg,nring)
  IF (nring > 0) CALL rjgrm(nring,tgroup,rjg)

! find group which contains '.' (if none, escape)
  DO i=1,nca
    IF (INDEX(tgroup(i),'.')/=0) THEN ; j1 = i ; EXIT ; ENDIF
  ENDDO
  IF (j1==0) EXIT rdloop

! check for various criegee configuration
  locriegee=.FALSE.
  IF (INDEX(tgroup(j1),'.(OO.)')/=0)  locriegee=.TRUE.
  IF (INDEX(tgroup(j1),'.(ZOO.)')/=0) locriegee=.TRUE.
  IF (INDEX(tgroup(j1),'.(EOO.)')/=0) locriegee=.TRUE.

!---------------------------------------------------------
!  '(O.)' and (OO.) on a non 'Cd' carbon
!---------------------------------------------------------
  IF (tgroup(j1)(1:2)/='Cd') THEN

! plain '(O.)' exits but 'CO(O.)' group continues
    IF (INDEX(tgroup(j1),'(O.)')/=0) THEN
      IF (tgroup(j1)(1:2)/='CO') EXIT rdloop 
    ENDIF

! peroxy (but not criegee: that's treated later) 
    IF ( (INDEX(tgroup(j1),'(OO.)')/=0).AND. (.NOT. locriegee) ) THEN

! If needed, transform: X-O-CO(OO.) -> X. + CO2 (+O2)
      IF (INDEX(tgroup(j1),'CO(OO.)')/=0)  THEN
        DO i=1,nca
          IF (tbond(i,j1)==3)  THEN
            tgroup(j1) = ' ' 
            np = np+1 ;  coprod(np) = 'CO2 '
            tbond(j1,i) = 0  ; tbond(i,j1) = 0
            DO j=1,nca
              IF (tbond(i,j)==3) THEN
                nc = INDEX(tgroup(j),' ')  ;  tgroup(j)(nc:nc) = '.'
                tgroup(i) = ' '
                tbond(j,i) = 0  ;  tbond(i,j) = 0
! if X is another -O- center then : R(-O-.) => R(O.)
                IF (tgroup(j)(1:4)=='-O-.') THEN
                  DO k=1,nca
                    IF (tbond(k,j)/=0) THEN
                      nc = INDEX(tgroup(k),' ') 
                      tgroup(k)(nc:nc+3) = '(O.)'
                      tgroup(j) = ' '
                      tbond(j,k) = 0  ;  tbond(k,j) = 0
                    ENDIF
                  ENDDO
                ENDIF 
              ENDIF
            ENDDO
            CALL rebond(tbond,tgroup,tchem,nring)
            nref=nref+1 ; com(nref)='XOCO(OO)'
            CYCLE rdloop
          ENDIF
        ENDDO
      ENDIF

! C(OOH)(OO.) => O2 + C.(OOH) - overwrite if C(OOH)(OO.)
      IF (INDEX(tgroup(j1),'(OOH)')/=0)  THEN
        nc = INDEX(tgroup(j1),'(OO.)')
        nh = INDEX(tgroup(j1),' ')
        tgroup(j1) = tgroup(j1)(1:nc-1) // '.' // tgroup(j1)(nc+5:nh) // ' ' ! rm (OO.)
        CALL rebond(tbond,tgroup,tchem,nring)
        nref=nref+1 ; com(nref)='C(OOH)(OO)'
        CYCLE rdloop
      ENDIF

      EXIT rdloop   ! regular RO2: exit
    ENDIF
  ENDIF

!---------------------------------------------------------
!  '(.)' on a Cd carbon: C=C.
!---------------------------------------------------------
  IF (tgroup(j1)(1:2) == 'Cd') THEN

! find partner of the C(j2)=C(j1)(.) bond:
    DO i=1,nca
      IF (tbond(i,j1)==2) THEN ; j2 = i ; EXIT; ENDIF
    ENDDO

! check if j2 has another bouble bond (Cd=Cd=Cd). if yes, keep Cd notation for j2
    keep = 0
    DO i=1,nca
      IF (tbond(i,j2)==2 .AND. i/=j1) keep = 1
    ENDDO

! two radicals on carbon: >Cd=Cd..  -->  >Cd=Cd.(OO.)
    IF (tgroup(j1)(1:5)== 'Cd.. ') THEN
      tgroup(j1) = 'Cd.(OO.)' ; locriegee=.TRUE.
    ENDIF

! '.(OO.)' radicals on C=C: R2Cd=Cd.(OO.) --> RCOR + CO
!    IF (INDEX(tgroup(j1),'.(OO.)')/=0) THEN
    IF (locriegee) THEN
      np = np+1  ;  coprod(np) = 'CO '
      tgroup(j1) = ' '
      tbond(j1,j2) = 0  ;  tbond(j2,j1) = 0
      IF (keep==1) THEN
        tgroup(j2) = 'CdO'
      ELSE
        IF (tgroup(j2)(1:4)=='CdH2') THEN
          np = np+1 ; tchem = ' ' ; coprod(np) = 'CH2O ' ; EXIT rdloop
        ELSE IF (tgroup(j2)(1:3)=='CdH') THEN
          pold = 'CdH' ; pnew = 'CHO'
          CALL swap(tgroup(j2),pold,tempkg,pnew)
          tgroup(j2) = tempkg
        ELSE IF (tgroup(j2)(1:3)=='CdO') THEN
          np = np+1 ; tchem = ' ' ; coprod(np) = 'CO2 ' ; EXIT rdloop 
        ELSE IF (tgroup(j2)(1:2)=='Cd') THEN
          pold = 'Cd'  ;  pnew = 'CO'
          CALL swap(tgroup(j2),pold,tempkg,pnew)
          tgroup(j2) = tempkg
        ELSE IF (tgroup(j2)(1:2)/='Cd') THEN
          mesg="Molecule not identified "
          CALL stoperr(progname,mesg,inchem)
        ENDIF
      ENDIF
      CALL rebond(tbond,tgroup,tchem,nring)

! vinyl-peroxy  Cd=Cd(OO.)  -->  C(O.)-CO
    ELSE IF (INDEX(tgroup(j1),'(OO.)')/=0) THEN
      tbond(j1,j2) = 1  ;  tbond(j2,j1) = 1
      IF (tgroup(j1)(1:3)=='CdH') THEN
        tgroup(j1) = 'CHO'
      ELSE
        pold = '(OO.)'  ;  pnew = ' '
        CALL swap(tgroup(j1),pold,tempkg,pnew)
        pold = 'Cd'  ;  pnew = 'CO'
        CALL swap(tempkg,pold,tgroup(j1),pnew)
      ENDIF
      IF (keep== 0) THEN
        pold = 'Cd'  ;  pnew = 'C'
        CALL swap(tgroup(j2),pold,tempkg,pnew)
        tgroup(j2) = tempkg
      ENDIF
      nc = INDEX(tgroup(j2),' ')  ;  tgroup(j2)(nc:nc+3) = '(O.)'
      CALL rebond(tbond,tgroup,tchem,nring)
      nref=nref+1 ; com(nref)='CdCd(OO)'

! vinyl-oxy radicals:  Cd=Cd(O.)  -->  C.-CO
    ELSE IF(INDEX(tgroup(j1),'(O.)')/=0) THEN
      tbond(j1,j2) = 1  ;  tbond(j2,j1) = 1
      IF (tgroup(j1)(1:3)=='CdH') THEN
        tgroup(j1) = 'CHO'
      ELSE
        pold = '(O.)'  ;  pnew = ' '
        CALL swap(tgroup(j1),pold,tempkg,pnew)
        pold = 'Cd'  ;  pnew = 'CO'
        CALL swap(tempkg,pold,tgroup(j1),pnew)
      ENDIF
      IF (keep==0) THEN
        pold = 'Cd'  ;  pnew = 'C'
        CALL swap(tgroup(j2),pold,tempkg,pnew)
        tgroup(j2) = tempkg
      ENDIF
      nc = INDEX(tgroup(j2),' ')  ;  tgroup(j2)(nc:nc) = '.'
      CALL rebond(tbond,tgroup,tchem,nring)
      nref=nref+1 ; com(nref)='CdCd(O)'

! vinyl radicals:  Cd=Cd(.)  -->  Cd=Cd(OO.) 
!                  Cd=Cd(X)(.)  -->  Cd=CdO + rad
    ELSE IF (tgroup(j1)(1:4)=='CdH.') THEN
      tgroup(j1) = 'CdH(OO.)'
      CALL rebond(tbond,tgroup,tchem,nring)
    ELSE IF (tgroup(j1)(1:4)=='Cd. ') THEN
      tgroup(j1) = 'Cd(OO.)'
      CALL rebond(tbond,tgroup,tchem,nring)
    ELSE IF (tgroup(j1)(1:2) == 'Cd') THEN
      IF (INDEX(tgroup(j1),'(OH)')/=0) THEN
        tgroup(j1) = 'CdO'  ;  np = np+1  ;  coprod(np) = 'HO2 '
      ELSE IF (INDEX(tgroup(j1),'(OOH)')/=0) THEN
        tgroup(j1) = 'CdO'  ;  np = np+1  ;  coprod(np) = 'HO '
      ELSE IF (INDEX(tgroup(j1),'(OOOH)')/=0) THEN
        tgroup(j1) = 'CdO'  ;  np = np+1  ;  coprod(np) = 'HO2 '
      ELSE IF (INDEX(tgroup(j1),'(ONO2)')/=0) THEN
        tgroup(j1) = 'CdO'  ;  np = np+1  ;  coprod(np) = 'NO2  '
      ELSE IF (INDEX(tgroup(j1),'Cd(NO2).')/=0) THEN
        tgroup(j1)='Cd(NO2)(OO.)'
      ELSE
        mesg="Molecule not identified (1): "
        CALL stoperr(progname,mesg,inchem)
      ENDIF
      CALL rebond(tbond,tgroup,tchem,nring)
      nref=nref+1 ; com(nref)='CdCd(X)(.)'
    ENDIF
         
!---------------------------------------------------------
! radical on 'CO' group :
!---------------------------------------------------------
  ELSE IF (tgroup(j1)(1:2)=='CO') THEN
                  
! cannot be a criegee:
    IF (locriegee) THEN
      mesg="Molecule not identified (2) "
      CALL stoperr(progname,mesg,inchem)
    ENDIF

! acyl oxy radicals: R-CO(O.) --> R.  +  CO2
    IF (tgroup(j1)(1:6)=='CO(O.)') THEN
      tgroup(j1) = ' '
      np = np+1 ; coprod(np) = 'CO2 '
      loopa: DO i=1,nca
        IF (tbond(j1,i)==0) CYCLE loopa
        IF (tbond(j1,i)==1) THEN 
          j2 = i
          nc = INDEX(tgroup(j2),' ')  ;  tgroup(j2)(nc:nc) = '.'
          tbond(j1,j2) = 0  ;  tbond(j2,j1) = 0
          CALL rebond(tbond,tgroup,tchem,nring)
          EXIT loopa
          
        ELSE IF (tbond(j1,i)==3) THEN                ! R-C-O-CO(O.) --> R-C(O.) + CO2
          tgroup(j1) = ' '  ;  tbond(i,j1) = 0  ;  tbond(j1,i) = 0
          DO j= 1,nca
            IF (tbond(i,j)==3) THEN
              tgroup(i) = ' '  ;  tbond(i,j) = 0  ;  tbond(j,i) = 0
              nc = INDEX(tgroup(j),' ')
              tgroup(j)(nc:nc+3) = '(O.)'

              DO k=1,nca
                IF ((tbond(j,k)==3).AND.(tgroup(j)(1:3)=='-O-')) THEN  ! R-C-O--O-CO(O.) --> R-C(OO.) + CO2
                  nc = INDEX(tgroup(k),' ')  ;  tgroup(k)(nc:nc+4) = '(OO.)'
                  tgroup(j) = ' '  ;  tbond(j,k) = 0  ;  tbond(k,j) = 0
                  IF (tgroup(k)(1:3)=='CHO') THEN          !  CHO-O--O-CO(O.) --> CO2 + CO + HO2
                    tgroup(k)=' ' 
                    np = np+1  ;  coprod(np) = 'CO '
                    np = np+1  ;  coprod(np) = 'HO2 '
                    tchem=' '  ; EXIT loopa ! need to exit here, nothing left in chem
                  ELSE IF (tgroup(k)(1:6)=='CO(OH)') THEN  !  CO(OH)-O--O-CO(O.) --> CO2 + CO2 + HO2
                    tgroup(k)=' '  
                    np = np+1  ;  coprod(np) = 'CO2 '
                    np = np+1  ;  coprod(np) = 'HO2 '
                    tchem=' '  ; EXIT loopa ! need to exit here, nothing left in chem
                  ELSE IF (tgroup(k)(1:9)=='CO(OONO2)') THEN  !  
                    tgroup(k)=' '  
                    np = np+1  ;  coprod(np) = 'CO2 '
                    np = np+1  ;  coprod(np) = 'NO3 '
                    tchem=' '  ; EXIT loopa ! need to exit here, nothing left in chem
                  ELSE IF (tgroup(k)(1:8)=='CO(ONO2)') THEN  !  
                    tgroup(k)=' '  
                    np = np+1  ;  coprod(np) = 'CO2 '
                    np = np+1  ;  coprod(np) = 'NO2 '
                    tchem=' '  ; EXIT loopa ! need to exit here, nothing left in chem
                  ENDIF
                ENDIF
              ENDDO
              CALL rebond(tbond,tgroup,tchem,nring)
              EXIT loopa
            ENDIF
          ENDDO
        ENDIF
      ENDDO loopa      

! acyl radicals: C-C-CO --> C-C-CO(OO.) but  C-CO-CO. -->  C-CO. + CO
    ELSE IF (tgroup(j1)(1:3)=='CO.') THEN
      DO i=1,nca
        IF (tbond(j1,i) /=0) THEN ; j2 = i  ;  EXIT  ;  ENDIF
      ENDDO

      IF (tgroup(j2)(1:2)=='CO') THEN
        tgroup(j1) = ' '  ;  np = np+1  ;  coprod(np) = 'CO   '
        tbond(j1,j2) = 0  ;  tbond(j2,j1) = 0
        nc = INDEX(tgroup(j2),' ')  ;  tgroup(j2)(nc:nc) = '.'
        CALL rebond(tbond,tgroup,tchem,nring)
      ELSE
        tgroup(j1) = 'CO(OO.)'
        CALL rebond(tbond,tgroup,tchem,nring)
      ENDIF
    ENDIF

!---------------------------------------------------------
! radical on alkyl carbon : criegee
!---------------------------------------------------------
!  ELSE IF (INDEX(tgroup(j1),'.(OO.)')/=0) THEN
  ELSE IF (locriegee) THEN

! functional groups on '.(OO.)' carbon:
    IF (INDEX(tgroup(j1),'(OOH)')/=0) THEN
      tgroup(j1) = 'CO(OO.)'  ;  np = np+1  ;  coprod(np) = 'HO '
      CALL rebond(tbond,tgroup,tchem,nring)
    ELSE IF (INDEX(tgroup(j1),'(OH)')/=0) THEN
      tgroup(j1) = 'CO(O.)'  ;  np = np+1  ;  coprod(np) = 'HO '
      CALL rebond(tbond,tgroup,tchem,nring)
    ELSE IF (INDEX(tgroup(j1),'(ONO2)')/=0) THEN
      tgroup(j1) = 'CO(OO.)'  ;  np = np+1  ;  coprod(np) = 'NO2 '
      CALL rebond(tbond,tgroup,tchem,nring)
    ELSE 
      EXIT rdloop 
    ENDIF

!---------------------------------------------------------
! radical on alkyl carbon : di-radical
!---------------------------------------------------------
  ELSE IF (INDEX(tgroup(j1),'..')/=0) THEN
    pold = '..'  ;  pnew = '.(OO.)' ; locriegee=.TRUE.
    CALL swap(tgroup(j1),pold,tempkg,pnew)
    tgroup(j1) = tempkg
    CALL rebond(tbond,tgroup,tchem,nring)

!---------------------------------------------------------
! radical on a CH2. : RCH2. --> RCH2(OO.)
!---------------------------------------------------------
  ELSE IF (INDEX(tgroup(j1),'CH2')/=0) THEN
    tgroup(j1) = 'CH2(OO.)'
    CALL rebond(tbond,tgroup,tchem,nring)

!---------------------------------------------------------
! radical on a >CH.
!---------------------------------------------------------
  ELSE IF (INDEX(tgroup(j1),'CH')/=0) THEN

    IF (INDEX(tgroup(j1),'(OOH)')/=0) THEN       ! -CH(OOH). --> -CHO + OH
      tgroup(j1) = 'CHO'  ;  np = np+1  ;  coprod(np) = 'HO   '
    ELSE IF (INDEX(tgroup(j1),'(OH)')/=0) THEN   ! -CH(OH). --> -CHO + HO2
      tgroup(j1) = 'CHO'  ;  np = np+1  ;  coprod(np) = 'HO2  '
    ELSE IF (INDEX(tgroup(j1),'(ONO2)')/=0) THEN ! -CH(ONO2). --> -CHO + NO2
      tgroup(j1) = 'CHO'  ;  np = np+1  ;  coprod(np) = 'NO2  '
    ELSE IF (INDEX(tgroup(j1),'(NO2)')/=0) THEN  ! -CH(NO2). --> -CH(NO2)(OO.)
        tgroup(j1)='CH(NO2)(OO.)'
    ELSE IF (INDEX(tgroup(j1),'(OOOH)')/=0) THEN   ! -CH(OOOH). --> -CHO + HO2
      tgroup(j1) = 'CHO'  ;  np = np+1  ;  coprod(np) = 'HO2  '
    ELSE                                         ! >CH. --> >CH(OO.)
      tgroup(j1) = 'CH(OO.)'
    ENDIF         
    CALL rebond(tbond,tgroup,tchem,nring)

!---------------------------------------------------------
! radical on alkyl carbon : tertiary alkyl
!---------------------------------------------------------
  ELSE IF ( (INDEX(tgroup(j1),'Cd'    )==0).AND.(INDEX(tgroup(j1),'CO' )==0).AND. &
            (.NOT. locriegee)              .AND.(INDEX(tgroup(j1),'..' )==0).AND. &
            (INDEX(tgroup(j1),'CH'    )==0) ) THEN

    IF ((INDEX(tgroup(j1),'(OOH)')/=0).AND.(INDEX(tgroup(j1),'(ONO2)')/=0)) THEN
      pold = '.'  ;  pnew = ' '
      CALL swap(tgroup(j1),pold,tempkg,pnew) ; tgroup(j1) = tempkg  
      pold = 'C'  ;  pnew = 'CO'
      CALL swap(tgroup(j1),pold,tempkg,pnew) ; tgroup(j1) = tempkg
      pold = '(ONO2)'  ;   pnew = ' '
      CALL swap(tgroup(j1),pold,tempkg,pnew) ; tgroup(j1) = tempkg
      np = np+1  ;  coprod(np) = 'NO2  '

    ELSE IF (INDEX(tgroup(j1),'(OOH)')/=0) THEN
      pold = '.'  ;  pnew = ' '
      CALL swap(tgroup(j1),pold,tempkg,pnew) ; tgroup(j1) = tempkg  
      pold = 'C'  ;  pnew = 'CO'
      CALL swap(tgroup(j1),pold,tempkg,pnew) ; tgroup(j1) = tempkg
      pold = '(OOH)'  ;   pnew = ' '
      CALL swap(tgroup(j1),pold,tempkg,pnew) ; tgroup(j1) = tempkg
      np = np+1  ;  coprod(np) = 'HO   '

    ELSE IF ((INDEX(tgroup(j1),'(OH)')/=0).AND.(INDEX(tgroup(j1),'(ONO2)')/=0))THEN
      pold = '.'  ;  pnew = ' '
      CALL swap(tgroup(j1),pold,tempkg,pnew) ; tgroup(j1) = tempkg
      pold = '(ONO2)'  ;  pnew = ' '
      CALL swap(tgroup(j1),pold,tempkg,pnew) ; tgroup(j1) = tempkg
      pold = 'C'  ;  pnew = 'CO'
      CALL swap(tgroup(j1),pold,tempkg,pnew) ; tgroup(j1) = tempkg
      np = np+1 ;  coprod(np) = 'NO2  '

    ELSE IF (INDEX(tgroup(j1),'(OH)')/=0) THEN
      pold = '.'  ;  pnew = ' '
      CALL swap(tgroup(j1),pold,tempkg,pnew) ;  tgroup(j1) = tempkg
      pold = 'C'  ;  pnew = 'CO'
      CALL swap(tgroup(j1),pold,tempkg,pnew) ;  tgroup(j1) = tempkg
      pold = '(OH)'  ;  pnew = ' '
      CALL swap(tgroup(j1),pold,tempkg,pnew) ;  tgroup(j1) = tempkg
      np = np+1  ;  coprod(np) = 'HO2  '

    ELSE IF (INDEX(tgroup(j1),'(ONO2)')/=0) THEN
      pold = '.'  ;   pnew = ' '
      CALL swap(tgroup(j1),pold,tempkg,pnew) ;  tgroup(j1) = tempkg
      pold = 'C'  ;  pnew = 'CO'
      CALL swap(tgroup(j1),pold,tempkg,pnew) ;  tgroup(j1) = tempkg
      pold = '(ONO2)'  ;  pnew = ' '
      CALL swap(tgroup(j1),pold,tempkg,pnew) ;  tgroup(j1) = tempkg
      np = np+1  ;  coprod(np) = 'NO2  '

    ELSE IF (INDEX(tgroup(j1),'(OOOH)')/=0) THEN
      pold = '.'  ;  pnew = ' '
      CALL swap(tgroup(j1),pold,tempkg,pnew) ;  tgroup(j1) = tempkg
      pold = 'C'  ;  pnew = 'CO'
      CALL swap(tgroup(j1),pold,tempkg,pnew) ;  tgroup(j1) = tempkg
      pold = '(OOOH)'  ;  pnew = ' '
      CALL swap(tgroup(j1),pold,tempkg,pnew) ;  tgroup(j1) = tempkg
      np = np+1  ;  coprod(np) = 'HO2  '

    ELSE
      pold = '.'  ;  pnew='(OO.)'
      CALL swap(tgroup(j1),pold,tempkg,pnew) ;  tgroup(j1) = tempkg
    ENDIF

    CALL rebond(tbond,tgroup,tchem,nring)
  ENDIF 

! exit loop if tchem is not a radical species anymore
  IF (INDEX(tchem,'.')==0) EXIT rdloop  
ENDDO rdloop

  prod = tchem
! check that np is small enough before returnig value
  IF (np > SIZE(coprod)) THEN
    mesg="Too many coprod."
    CALL stoperr(progname,mesg,inchem)
  ENDIF

END SUBROUTINE multip

! ======================================================================
! Purpose: for a C1 species provided as input, return a C1 species that 
! can be handle by gecko (i.e. in the set of allowed C1 species) and the 
! various "coproducts" linked to the "transformation".
! ======================================================================
SUBROUTINE single(tchem,prod,coprod)
  USE keyparameter, ONLY: mxlgr
  USE reactool, ONLY: swap
  USE toolbox, ONLY: stoperr
  IMPLICIT NONE

  CHARACTER(LEN=*),INTENT(IN) :: tchem      ! C1 species as input
  CHARACTER(LEN=*),INTENT(OUT):: prod       ! C1 species as output
  CHARACTER(LEN=*),INTENT(OUT):: coprod(:)  ! coproducts 

  CHARACTER(LEN=mxlgr) :: tgroup, pnew, pold, tempkg
  INTEGER :: np
  CHARACTER(LEN=9),PARAMETER :: progname='single '
  CHARACTER(LEN=70)          :: mesg

! initialise
  prod = ' ' ; tempkg = tchem(1:mxlgr) ; coprod(:)=' '

! double bond not allowed:
  IF (INDEX(tchem,'Cd')/=0) THEN
    mesg="illegal single C species (Cd): "
    CALL stoperr(progname,mesg,tchem)
  ENDIF 

! --------------------
! CH2.. ---> CH2.(OO.)
! --------------------
  IF (tchem(1:5)=='CH2..') THEN
    prod = 'CH2.(OO.)'  ;  RETURN
  ENDIF

! --------------------
! Break single carbon peroxy
! --------------------
  IF ((INDEX(tchem,'(OO.)')/=0).AND.(INDEX(tchem,'.(OO.)')==0))  THEN
    IF (tchem(1:2)=='CO') THEN
      IF      (INDEX(tchem,'(OH)')   /=0) THEN ; prod='CO2' ; coprod(1) = 'HO2 '
      ELSE IF (INDEX(tchem,'(OOH)')  /=0) THEN ; prod='CO2' ; coprod(1) = 'HO '
      ELSE IF (INDEX(tchem,'(ONO2)') /=0) THEN ; prod='CO2' ; coprod(1) = 'NO2 '
      ELSE IF (INDEX(tchem,'(OONO2)')/=0) THEN ; prod='CO'  ; coprod(1) = 'NO2 '
      ELSE IF (INDEX(tchem,'(NO2)')  /=0) THEN ; prod='CO'  ; coprod(1) = 'NO2 '
      ELSE
        mesg="illegal single C species (1A) "
        CALL stoperr(progname,mesg,tchem)
      ENDIF

    ELSE IF (tchem(1:3)=='CH2') THEN  ;  prod = 'CH2O'
      IF      (INDEX(tchem,'(OOH)') /=0)   THEN ; coprod(1) = 'HO  '
      ELSE IF (INDEX(tchem,'(OOOH)')  /=0) THEN ; coprod(1) = 'HO2 '
      ELSE IF (INDEX(tchem,'(OH)')  /=0)   THEN ; coprod(1) = 'HO2 '
      ELSE IF (INDEX(tchem,'(ONO2)')/=0)   THEN ; coprod(1) = 'NO2 '
      ELSE IF (INDEX(tchem,'(NO2)') /=0)   THEN
        prod = 'CH2.(OO.)'  ;  coprod(1) = 'NO2 '
      ELSE
        mesg="illegal single C species (1b) "
        CALL stoperr(progname,mesg,tchem)
      ENDIF

    ELSE IF (tchem(1:3)=='CHO') THEN
      prod=tchem
      RETURN

    ELSE IF (tchem(1:2)=='CH') THEN ! need additionnal work
      IF (INDEX(tchem,'(OOH)')/=0) THEN
         prod = 'CO ' ; coprod(1) = 'HO2 ' ; coprod(2) = 'HO '
      ELSE IF (INDEX(tchem,'(ONO2)')/=0) THEN
         prod = 'CO ' ; coprod(1) = 'HO2 ' ; coprod(2) = 'NO2 '
      ELSE IF (INDEX(tchem,'(NO2)') /=0) THEN
         prod = 'CO ' ; coprod(1) = 'HO ' ;  coprod(2) = 'NO2 '
      ELSE IF (INDEX(tchem,'(OH)') /=0) THEN
         prod = 'CO2 ' ; coprod(1) = 'HO2 ' ; coprod(2) = 'HO '
      ELSE
        mesg="illegal single C species (1c) "
        CALL stoperr(progname,mesg,tchem)
      ENDIF

    ELSE
      mesg="illegal single C species (1d) "
      CALL stoperr(progname,mesg,tchem)
    ENDIF

! --------------------
! Break single carbon (O.)(OOOH)
! --------------------
  ELSE IF ((INDEX(tchem,'(O.)')/=0).AND.(INDEX(tchem,'(OOOH)')/=0))  THEN
    IF (tchem(1:3)=='CH2') THEN  ;  prod = 'CH2O'
      IF (INDEX(tchem,'(OOOH)')  /=0) THEN ; coprod(1) = 'HO  '
      ELSE
        mesg="illegal single C species (1d) "
        CALL stoperr(progname,mesg,tchem)
      ENDIF
    ENDIF
    
! --------------------
! carbonyl: formyl
! --------------------
  ELSE IF ((tchem(1:3)=='CHO').OR.(tchem(1:4)=='CHO.')) THEN
     prod = 'CO ' ; coprod(1) = 'HO2 '

     IF (tchem=='CHO(O.)') THEN
       prod = 'CO2 ' ; coprod(1) = 'HO2  '
     ENDIF
              
! --------------------
! carbonyl
! --------------------
  ELSE IF (tchem(1:2)=='CO') THEN

! carbonyl-carbene:
    IF ( (INDEX(tchem,'..')/=0 ) .OR. (tchem=='CO.')) THEN ; prod = 'CO '

! carbonyl-criegee:
    ELSE IF (INDEX(tchem,'.(OO.)')/=0) THEN ; prod = 'CO '

! acyl-oxy  and  radical on carbonyl:
    ELSE 
      IF (INDEX(tchem,'(O.)')/=0) THEN ; prod = 'CO2 '
      ELSE                             ; prod = 'CO  '
      ENDIF
     
      IF      (INDEX(tchem,'(OH)')   /=0) THEN ; coprod(1) = 'HO '
      ELSE IF (INDEX(tchem,'(OOH)')  /=0) THEN ; coprod(1) = 'HO2 '
      ELSE IF (INDEX(tchem,'(ONO2)') /=0) THEN ; coprod(1) = 'NO3 '
      ELSE IF (INDEX(tchem,'(OONO2)')/=0) THEN ; coprod(1) = 'NO2 '
      ELSE IF (INDEX(tchem,'(NO2)')  /=0) THEN ; coprod(1) = 'NO2 '
      ELSE IF (INDEX(tchem,'H')      /=0) THEN ; coprod(1) = 'HO2 '
      ELSE
        mesg="illegal single C species (2) "
        CALL stoperr(progname,mesg,tchem)
      ENDIF
    ENDIF  

  ELSE IF (INDEX(tchem,'(O.)')/=0) THEN
    IF (tchem=='CH2(OOH)(O.)') THEN
      prod= 'CHO(OOH)' ; coprod(1) = 'HO2'
    ELSE IF (tchem=='CH2(ONO2)(O.)') THEN
! BA : chge decomposition to conserve H.
      prod= 'CH2O' ; coprod(1) = 'NO3'
    ELSE IF (tchem=='CH2(OH)(O.)') THEN
      prod= 'CHO(OH)' ; coprod(1) = 'HO2'
    ELSE IF (tchem(1:3)=='CH(') THEN
      prod= 'CO' ; coprod(1) = 'HO2' ; np=1
      IF (INDEX(tchem,'(OH)')/=0) THEN
        np=np+1 ; coprod(np) = 'HO '
      ELSE IF (INDEX(tchem,'(OOH)')/=0) THEN
        np=np+1 ; coprod(np) = 'HO2 '
      ELSE IF (INDEX(tchem,'(ONO2)')/=0) THEN
        np=np+1 ; coprod(np) = 'NO3 '
      ELSE IF (INDEX(tchem,'(OONO2)')/=0) THEN 
        np=np+1 ; coprod(np) = 'NO2 '
      ELSE IF (INDEX(tchem, '(NO2)') /= 0) THEN
        np=np+1 ; coprod(np) = 'NO2  '
      ENDIF         
    ELSE  
      prod = tchem
!      RETURN
    ENDIF

! --------------------
! alkyl carbon: carbenes = C..
! --------------------
  ELSE IF (INDEX(tchem,'..')/=0) THEN
    IF (tchem(1:2)=='CH') THEN
      IF (INDEX(tchem,'(ONO2)')/=0) THEN
         prod = 'CO ' ; coprod(1) = 'HO2 ' ; coprod(2) = 'NO2 '
      ELSE IF (INDEX(tchem,'(OOH)')/=0) THEN
         prod = 'CO ' ; coprod(1) = 'HO2 ' ; coprod(2) = 'HO  '
      ELSE IF (INDEX(tchem,'(OH)')/=0) THEN
         prod = 'CO2 ' ; coprod(1) = 'HO2 ' ; coprod(2) = 'HO  '
      ELSE
        mesg="illegal single C species (3) "
        CALL stoperr(progname,mesg,tchem)
      ENDIF

    ELSE IF (tchem(1:1)=='C') THEN
      IF (INDEX(tchem,'(ONO2)')/=0) THEN
        pold = '(ONO2)' ; pnew = ' '
        CALL swap(tempkg,pold,tgroup,pnew)
        IF (INDEX(tgroup,'(ONO2)')/=0) THEN
          prod = 'CO2 ' ; coprod(1) = 'NO2 ' ;  coprod(2) = 'NO2 '
        ELSE IF (INDEX(tgroup,'(OOH)')/=0) THEN
          prod = 'CO2 ' ; coprod(1) = 'NO2 ' ; coprod(2) = 'HO '
        ELSE IF (INDEX(tgroup,'(OH)')/=0) THEN
          prod = 'CO2 ' ; coprod(1) = 'NO2 ' ; coprod(2) = 'HO2 '
        ELSE
          mesg="illegal single C species (4) "
          CALL stoperr(progname,mesg,tchem)
        ENDIF

      ELSE IF (INDEX(tchem,'(OOH)')/=0) THEN
         pold = '(OOH)' ; pnew = ' '
         CALL swap(tempkg,pold,tgroup,pnew)
         IF (INDEX(tgroup,'(OOH)')/=0) THEN
           prod = 'CO2 ' ;  coprod(1) = 'HO ' ;  coprod(2) = 'HO '
         ELSE IF (INDEX(tgroup,'(OH)')/=0) THEN
           prod = 'CO2 ' ;  coprod(1) = 'HO ' ; coprod(2) = 'HO2  '
         ELSE
           mesg="illegal single C species (5) "
           CALL stoperr(progname,mesg,tchem)
         ENDIF

      ELSE IF (INDEX(tchem,'(OH)')/=0) THEN
        pold = '(OH)' ; pnew = ' '
        CALL swap(tempkg,pold,tgroup,pnew)
        IF (INDEX(tgroup,'(OH)')/=0) THEN
          prod = 'CO2 ' ; coprod(1) = 'HO ' ; coprod(2) = 'HO   '
        ELSE
         mesg="illegal single C species (6) "
         CALL stoperr(progname,mesg,tchem)
        ENDIF

      ELSE
        mesg="illegal single C species (7) "
        CALL stoperr(progname,mesg,tchem)
      ENDIF

    ELSE
      mesg="illegal single C species (7b) "
      CALL stoperr(progname,mesg,tchem)
    ENDIF

! --------------------
! alkyl carbon: criegee diradicals
! --------------------
  ELSE IF (INDEX(tchem,'.(OO.)')/=0) THEN

    IF (tchem(1:2)=='CH') THEN
      IF (INDEX(tchem,'(OOH)')/=0) THEN
        prod = 'CO ' ; coprod(1) = 'HO2 ' ; coprod(2) = 'HO '
      ELSE IF (INDEX(tchem,'(ONO2)')/=0) THEN
        prod = 'CO ' ; coprod(1) = 'HO2 ' ; coprod(2) = 'NO2 '
      ELSE IF (INDEX(tchem,'(OH)')/=0) THEN
        prod = 'CO2 ' ; coprod(1) = 'HO2 ' ; coprod(2) = 'HO '
      ELSE
        mesg="illegal single C species (7c) "
        CALL stoperr(progname,mesg,tchem)
      ENDIF
    
    ELSE IF (tchem(1:1)=='C') THEN
      IF (INDEX(tchem,'(ONO2)')/=0) THEN
        pold = '(ONO2)' ; pnew = ' '
        CALL swap(tempkg,pold,tgroup,pnew)
        IF (INDEX(tgroup,'(ONO2)')/=0) THEN
           prod = 'CO2 ' ; coprod(1) = 'NO2 ' ; coprod(2) = 'NO2 '
        ELSE IF (INDEX(tgroup,'(OOH)')/=0) THEN
           prod = 'CO2 ' ; coprod(1) = 'NO2 ' ; coprod(2) = 'HO '
        ELSE IF (INDEX(tgroup,'(OH)')/=0) THEN
           prod = 'CO2 ' ; coprod(1) = 'NO2 ' ; coprod(2) = 'HO2 '
        ELSE
          mesg="illegal single C species (8) "
          CALL stoperr(progname,mesg,tchem)
        ENDIF
    
      ELSE IF (INDEX(tchem,'(OOH)')/=0) THEN
        pold = '(OOH)'  ; pnew = ' '
        CALL swap(tempkg,pold,tgroup,pnew)
        IF (INDEX(tgroup,'(OOH)')/=0) THEN
           prod = 'CO2 ' ; coprod(1) = 'HO ' ; coprod(2) = 'HO '
        ELSE IF (INDEX(tgroup,'(OH)')/=0) THEN
           prod = 'CO2 ' ; coprod(1) = 'HO ' ; coprod(2) = 'HO2 '
        ELSE
          mesg="illegal single C species (9) "
          CALL stoperr(progname,mesg,tchem)
        ENDIF
    
      ELSE IF (INDEX(tchem,'(OH)')/=0) THEN
        pold = '(OH)' ; pnew = ' '
        CALL swap(tempkg,pold,tgroup,pnew)
        IF (INDEX(tgroup,'(OH)')/=0) THEN
           prod = 'CO2 ' ; coprod(1) = 'HO ' ; coprod(2) = 'HO '
        ELSE
          mesg="illegal single C species (10) "
          CALL stoperr(progname,mesg,tchem)
        ENDIF
    
      ELSE
        mesg="illegal single C species (11) "
        CALL stoperr(progname,mesg,tchem)
      ENDIF
    
    ELSE
      mesg="illegal single C species (12) "
      CALL stoperr(progname,mesg,tchem)
    ENDIF

! -----------------------------
! alkyl carbon: methyl
! ----------------------------
  ELSE IF (tchem(1:3)=='CH3') THEN
    prod = "CH3(OO.)"

! -----------------------------
! alkyl carbon: primary
! ----------------------------
  ELSE IF (tchem(1:3)=='CH2') THEN
    prod = 'CH2O'
    IF (INDEX(tchem,'(OOH)')/=0) THEN
       coprod(1) = 'HO '
    ELSE IF (INDEX(tchem,'(OH)') /=0) THEN
       coprod(1) = 'HO2 '
    ELSE IF (INDEX(tchem,'(OOOH)')/=0) THEN
       coprod(1) = 'HO2 '
    ELSE IF (INDEX(tchem,'(ONO2)')/=0) THEN
       coprod(1) = 'NO2 '
    ELSE IF (INDEX(tchem,'(NO2)')/=0) THEN
       prod = 'CH2.(OO.)'
       coprod(1) = 'NO2 '
!       RETURN
    ELSE
      mesg="illegal single C species (13) "
      CALL stoperr(progname,mesg,tchem)
    ENDIF

! -----------------------------
! alkyl carbon: secondary
! ----------------------------
  ELSE IF (tchem(1:2)=='CH') THEN 
    pold = '.' ; pnew = ' ' ;
    CALL swap(tempkg,pold,tgroup,pnew)
    tempkg = tgroup
    pold = 'CH' ; pnew = 'CHO'
    CALL swap(tempkg,pold,tgroup,pnew)
    IF (INDEX(tchem,'(OOH)')/=0) THEN
      coprod(2) = 'HO   '
      pold = '(OOH)' ; pnew = ' '
      CALL swap(tgroup,pold,prod(1:mxlgr),pnew)
    ELSE IF (INDEX(tchem,'(OH)')/=0) THEN
      coprod(2) = 'HO2  '
      pold = '(OH)' ; pnew = ' '
      CALL swap(tgroup,pold,prod(1:mxlgr),pnew)
    ELSE IF (INDEX(tchem,'(OOOH)')/=0) THEN
      coprod(2) = 'HO2  '
      pold = '(OOOH)' ; pnew = ' '
      CALL swap(tgroup,pold,prod(1:mxlgr),pnew)
    ELSE IF (INDEX(tchem,'(ONO2)')/=0) THEN
      coprod(2) = 'NO2  '
      pold = '(ONO2)' ; pnew = ' '
      CALL swap(tgroup,pold,prod(1:mxlgr),pnew)
    ELSE IF (INDEX(tchem,'(NO2)')/=0) THEN
      coprod(2) = 'NO2 '
      pold = '(NO2)' ; pnew = ' '
      CALL swap(tgroup,pold,prod(1:mxlgr),pnew)
    ELSE
      mesg="illegal single C species (14) "
      CALL stoperr(progname,mesg,tchem)
    ENDIF

! -----------------------------
! alkyl carbon: C
! ----------------------------
  ELSE IF (tchem(1:1)=='C') THEN 
    pold = '.' ;  pnew = ' '   
    CALL swap(tempkg,pold,tgroup,pnew)
    tempkg = tgroup
    pold = 'C' ; pnew = 'CO'
    CALL swap(tempkg,pold,tgroup,pnew)
    IF (INDEX(tchem,'(OOH)')/=0) THEN
      coprod(2) = 'HO   '
      pold = '(OOH)' ;  pnew = ' '
      CALL swap(tgroup,pold,prod(1:mxlgr),pnew)
    ELSE IF (INDEX(tchem,'(OH)')/=0) THEN
      coprod(2) = 'HO2  '
      pold = '(OH)' ; pnew = ' '
      CALL swap(tgroup,pold,prod(1:mxlgr),pnew)
    ELSE IF (INDEX(tchem,'(ONO2)')/=0) THEN
       coprod(2) = 'NO2 '
       pold = '(ONO2)' ;  pnew = ' '
       CALL swap(tgroup,pold,prod(1:mxlgr),pnew)
    ELSE IF (INDEX(tchem,'(NO2)')/=0) THEN
       coprod(2) = 'NO2 '
       pold = '(NO2)' ;   pnew = ' '
       CALL swap(tgroup,pold,prod(1:mxlgr),pnew)
    ELSE
      mesg="illegal single C species (15) "
      CALL stoperr(progname,mesg,tchem)
    ENDIF

! -----------------------------
! error
! ----------------------------
  ELSE
    mesg="illegal single C species (16) "
    CALL stoperr(progname,mesg,tchem)
  ENDIF

 END SUBROUTINE single

  
!=======================================================================
! Purpose: this subroutine checks if there is a possible electron 
! delocalisation in the case of a C=C bond in alpha position of the 
! alkyl radical: -C.-C=C<  <->  -C=C-C.<
! The program returns the most stable radical form. If the various 
! radical forms have similar probabilities, the output is a table 
! with the differents formula. 
!=======================================================================
SUBROUTINE deloc(chem,tchem,nchem,sc)
  USE keyparameter, ONLY: mxnode,mxlgr,mxcp,mxring
  USE rjtool, ONLY:rjgrm
  USE stdgrbond, ONLY: grbond
  USE mapping, ONLY: gettrack
  USE reactool, ONLY: swap,rebond
  USE toolbox, ONLY: stoperr 
  IMPLICIT NONE

  CHARACTER(LEN=*),INTENT(IN) :: chem     ! input species
  CHARACTER(LEN=*),INTENT(OUT):: tchem(:) ! output species
  INTEGER,INTENT(OUT) :: nchem            ! # of species as output
  REAL,INTENT(OUT)    :: sc(:)            ! "stoi.coef." of the various form

  CHARACTER(LEN=LEN(chem)) :: temp
  CHARACTER(LEN=mxlgr) :: group(mxnode)
  CHARACTER(LEN=mxlgr) :: tgr, pold, pnew
  INTEGER :: bond(mxnode,mxnode)
  INTEGER :: i,j,k,n1,n3,nring,dbflg
  INTEGER :: j1,j2,j3,ngr
  INTEGER :: rjg(mxring,2)
  REAL    :: sctmp
  INTEGER :: track(mxcp,mxnode)
  INTEGER :: trlen(mxcp)
  INTEGER :: ntr
  REAL    :: kbcadd1,kbcadd2
  REAL    :: k1,k3  
  CHARACTER(LEN=mxlgr) :: tgroup(mxnode)
  INTEGER :: tbond(mxnode,mxnode)
  CHARACTER(LEN=9),PARAMETER :: progname='deloc '
  CHARACTER(LEN=70)          :: mesg  
  
 
!-------------
! initialize:
!-------------
  j1=0  ;  j2=0  ;  j3=0  ;  ngr=0
  tchem(:)=' '   ;  tchem(1) = chem
  nchem=1  ;  sc(:)=0.  ;  sc(1) = 1.

! find groups and bond matrix from chem      
  CALL grbond(chem,group,bond,dbflg,nring)
  IF (nring > 0) CALL rjgrm(nring,group,rjg)
  tbond(:,:)=bond ; tgroup(:)=group

  ngr=COUNT(group/=' ')

! j1 = radical group (if none, escape)
  DO i=1,ngr
    IF (INDEX(group(i),'.')/=0) THEN ; j1=i ; EXIT ; ENDIF
  ENDDO
  IF (j1==0) RETURN

! only treat uncomplicated radical groups
  IF ( (INDEX(group(j1),'O.')==0).AND.(INDEX(group(j1),'Cd')==0) )  THEN

! "Paulot" case: disallow ').' except CH(OH).Cd with no conjugation.
    IF ((INDEX(group(j1),').')==0).OR.(INDEX(group(j1),'CH(OH).')/=0)) THEN

! j2 = double-bonded group attached to radical group (if none, escape)
      DO i=1,ngr
        IF ((bond(i,j1)/=0).AND.group(i)(1:2)=='Cd') THEN
          j2 = i ; EXIT
        ENDIF
      ENDDO
      IF (j2==0) RETURN

      DO i=1,ngr ! No delocalisation if j1 and j2 are part of an epoxide
        IF ((bond(i,j1)==3).AND.(bond(i,j2)==3).AND.(group(i)(1:3)=='-O-')) THEN
          RETURN
        ENDIF
      ENDDO

! j3 = group at other end of double bond
      DO i=1,ngr
        IF (bond(i,j2)==2) THEN ; j3=i ; EXIT ; ENDIF
      ENDDO

! only treat if not conjugated to C=O
      DO i=1,ngr
        IF (i/=j2) THEN
          IF (bond(i,j3)==1) THEN
            IF ((INDEX(group(i),'CO')/=0).OR.(INDEX(group(i),'CHO')/=0)) RETURN
          ENDIF
        ENDIF
      ENDDO

! n1: number of bonds to radical group j1
! n3: number of bonds to double-bonded group j3
      n1=0  ;  n3=0  
      DO i=1,ngr
        IF (bond(i,j1)/=0) n1=n1+1 
        IF (bond(i,j3)/=0) n3=n3+1 
      ENDDO

! move double bond one carbon over, to radical center 
      tbond(j1,j2)=2 ; tbond(j2,j1)=2
      tbond(j2,j3)=1 ; tbond(j3,j2)=1

! convert radical group j1 to double-bonded group (allow CH(OH). structures)
      IF (tgroup(j1)(1:2)=='C.')   tgroup(j1)(1:2)='Cd'
      IF (tgroup(j1)(1:3)=='CH.')  tgroup(j1)(1:3)='CdH'
      IF (tgroup(j1)(1:4)=='CH2.') tgroup(j1)(1:4)='CdH2'
      IF (tgroup(j1)(1:7)=='CH(OH).') tgroup(j1)(1:7)='CdH(OH)'

! convert 'double' group j3 to radical group
      IF      (tgroup(j3)(1:4)=='CdH2') THEN ; tgroup(j3)(1:4)='CH2.'
      ELSE IF (tgroup(j3)(1:3)=='CdH')  THEN ; tgroup(j3)(1:3)='CH.'
      ELSE IF (tgroup(j3)(1:2)=='Cd')   THEN
        pold='Cd'  ;  pnew='C' ;
        CALL swap(tgroup(j3),pold,tgr,pnew) ; tgroup(j3)=tgr
        k=INDEX(tgroup(j3),' ')  ;  tgroup(j3)(k:k)='.'
      ENDIF

! write (non-standard) formula for new product
      CALL rebond(tbond,tgroup,tchem(2),nring)

! prioritise possible products
! BA: original version= Ratio of mesomeric forms to 0/1 based on brch #  
!      IF (n1 > n3)      THEN ; sc(1) = 1.  ; sc(2) = 0.  ; nchem=1
!      ELSE IF (n3 > n1) THEN ; sc(1) = 0.  ; sc(2) = 1.  ; nchem=1
!      ELSE IF (n3==n1)  THEN ; sc(1) = 0.5 ; sc(2) = 0.5 ; nchem=2
!      ENDIF
! BA (May 2020): Ratio of mesomeric forms based on ratio of OH addition to C=C bond  
! according to the approach described by Peeter et al. (1007) for OH+diene (see also
! Szopa thesis). Rate used here are from Jenkin et al., 2018 (in 1E-12, page 9310)
      IF      (n1==1) THEN ; k1= 4.2   ! -C=CH2       
      ELSE IF (n1==2) THEN ; k1=26.3   ! -C=CHCH3     
      ELSE IF (n1==3) THEN ; k1=51.5   ! -C=C(CH3)CH3 
      ELSE                                            
        mesg="Too many bond at the n1 radical center."
        CALL stoperr(progname,mesg,chem)              
      ENDIF                                           
      IF      (n3==1) THEN ; k3= 4.2   ! -C=CH2       
      ELSE IF (n3==2) THEN ; k3=26.3   ! -C=CHCH3     
      ELSE IF (n3==3) THEN ; k3=51.5   ! -C=C(CH3)CH3 
      ELSE                                            
        mesg="Too many bond at the n3 center."        
        CALL stoperr(progname,mesg,chem)              
      ENDIF                                           
      sc(1)=k1/(k1+k3) ; sc(2)=1.-sc(1)  ; nchem=2    

! allow OH to attract electron density, as per Paulot ('09) by re-setting value of sc.
      IF (INDEX(tgroup(j1),'OH')/=0) THEN
        sc(1) = 0.65  ;  sc(2) = 0.35  ;  nchem = 2
      ENDIF

! post aromatic chemistry: for radical formed after -O-O- ring closure within an aromatic ring
      tbond(:,:)=bond ; tgroup(:)=group
      CALL gettrack(bond,j1,ngr,ntr,track,trlen)
      kbcadd1=4E-16  ;  kbcadd2=4E-16
      DO i=1,ntr
        IF (trlen(i) < 6) CYCLE
        IF ((track(i,2)==j2).AND.(track(i,3)==j3).AND. &
            (INDEX(group(track(i,5)),'(OH)')/=0) .AND. &
            (bond(track(i,6),j1)/=0)) THEN
          DO j=1,ngr
            IF (bond(track(i,4),j)==3) THEN
              DO k=1,ngr
                IF ((bond(track(i,6),k)==3).AND.(k/=j)) THEN
                  IF (group(j1)(1:3)=='C. ') kbcadd1=kbcadd1*1000.
                  IF (group(j3)(1:3)=='Cd ') kbcadd2=kbcadd2*1000.
                  IF (group(track(i,4))(1:2)=='C ') kbcadd1=kbcadd1*3.
                  IF (group(track(i,6))(1:2)=='C ') kbcadd2=kbcadd2*3.
                  sc(1)=kbcadd1 / (kbcadd1 + kbcadd2)
                  sc(2)=kbcadd2 / (kbcadd1 + kbcadd2)
                  nchem=2
                  EXIT
                ENDIF
              ENDDO
            ENDIF
          ENDDO
        ENDIF
      ENDDO

! if C=C-C.-O--O- and C.-C=C-O--O-, only keep the 1st one
      IF ((nchem==2).AND.(INDEX(tchem(1),'-O-')/=0)) THEN
        coo_loop : DO k=1,2
          CALL grbond(tchem(k),group,bond,dbflg,nring)
          IF (nring > 0) CALL rjgrm(nring,group,rjg)
          DO i=1,ngr
            IF (INDEX(group(i),'.')/=0) THEN ; j1=i ; EXIT ; ENDIF
          ENDDO
         
          DO i=1,ngr 
            IF (bond(j1,i)==3) THEN
              DO j=1,ngr
                IF ((bond(i,j)==3).AND.(j/=j1).AND.(group(j)=='-O-')) THEN
                   IF (k==1) THEN
                     sc(1)=1. ; sc(2)=0. ; tchem(2)=' ' ; nchem=1
                   ELSE
                     sc(2)=1. ; sc(1)=0. ; tchem(1)=' ' ; nchem=1
                   ENDIF
                   EXIT coo_loop
                ENDIF
              ENDDO
            ENDIF
          ENDDO
        ENDDO coo_loop
      ENDIF

! sort tchem according to yield (tchem(1) is with the greater yield)
      IF (sc(1) < sc(2)) THEN
        temp=tchem(1)      ;  sctmp=sc(1)
        tchem(1)=tchem(2)  ;  sc(1)=sc(2)
        tchem(2)=temp      ;  sc(2)=sctmp
      ENDIF
            
    ENDIF
  ENDIF

END SUBROUTINE deloc

! ======================================================================  
! PURPOSE : merge the list of coproducts (i.e. short names) returned 
! by radchck into a single list of products. Yield of each species
! in the merged list are also provided.
! ======================================================================  
SUBROUTINE lump_copd(chem,copdin,scin,nall,copdout,ycopd) 
  USE toolbox, ONLY: stoperr
  IMPLICIT NONE
  
  CHARACTER(LEN=*),INTENT(IN) :: chem         ! current species
  CHARACTER(LEN=*),INTENT(IN) :: copdin(:,:)  ! [2,j] list of copd in each list
  REAL,INTENT(IN)             :: scin(:)      ! yield of each list of copd
  INTEGER,INTENT(OUT)         :: nall         ! number of coproducts in the merged list
  CHARACTER(LEN=*),INTENT(OUT):: copdout(:)   ! merged list of copd
  REAL,INTENT(OUT)            :: ycopd(:)     ! yield of each copd in the merged list

  INTEGER :: i,j,n1,n2,mxcopd
  CHARACTER(LEN=12),PARAMETER :: progname='lump_copd '
  CHARACTER(LEN=70)          :: mesg
 
  n1=0 ; n2=0 ; nall=0
! check sizes 
  IF (SIZE(scin)/=2) THEN
    mesg="size of the imput table (sc) is out of range"
    CALL stoperr(progname,mesg,chem)
  ENDIF
  IF (SIZE(copdin,1)/=2) THEN
    mesg="size of the imput table (copd) is out of range"
    CALL stoperr(progname,mesg,chem)
  ENDIF
  IF (SIZE(ycopd)/=SIZE(copdout)) THEN
    mesg="size of the copdout and ycopd are not identical"
    CALL stoperr(progname,mesg,chem)
  ENDIF
  mxcopd=SIZE(copdout)
  
! add species in list 1
  DO i=1,SIZE(copdin,2)
    IF (copdin(1,i)==' ') EXIT
    nall=nall+1 ; copdout(nall)=copdin(1,i) ; ycopd(nall)=scin(1)
  ENDDO
  n1=nall
  
! add species in list 2
  lst2: DO i=1,SIZE(copdin,2)
    IF (copdin(2,i)==' ') EXIT lst2
    DO j=1,n1                              
      IF (copdin(2,i)==copdout(j)) THEN
        ycopd(j)=ycopd(j)+scin(2) ; CYCLE lst2
      ENDIF
    ENDDO
    nall=nall+1 ; copdout(nall)=copdin(2,i) ; ycopd(nall)=scin(2)
    IF (nall>mxcopd) THEN
      mesg="number of coproducts exceed the size allocated for copdout"
      CALL stoperr(progname,mesg,chem)
    ENDIF
  ENDDO lst2
  
END SUBROUTINE lump_copd 

!=======================================================================
! Purpose: this subroutine breaks the following structure: 
! -C.-O--O-R  ->  -CO +R(O.)
!=======================================================================
SUBROUTINE ooalkdec(chem,tchem,nchem,sc,dcom)
  USE keyparameter, ONLY: mxnode,mxlgr,mxcp,mxring
  USE rjtool, ONLY:rjgrm
  USE stdgrbond, ONLY: grbond
  USE mapping, ONLY: gettrack
  USE reactool, ONLY: swap,rebond
  USE toolbox, ONLY: stoperr,setbond
  USE ringtool, ONLY: findring
  USE fragmenttool, ONLY: fragm
  IMPLICIT NONE

  CHARACTER(LEN=*),INTENT(IN)  :: chem     ! input species
  CHARACTER(LEN=*),INTENT(OUT) :: tchem(:) ! output species
  INTEGER,INTENT(OUT)          :: nchem    ! # of species as output
  REAL,INTENT(OUT)             :: sc(:)    ! "stoi.coef." of the various form
  CHARACTER(LEN=*),INTENT(OUT) :: dcom     ! code for reference

  CHARACTER(LEN=mxlgr) :: group(mxnode)
  CHARACTER(LEN=mxlgr) :: tgr, pold, pnew
  INTEGER :: bond(mxnode,mxnode)
  INTEGER :: i,j,k,nring,dbflg
  INTEGER :: j1,j2,j3,j4,ngr,nc
  INTEGER :: rjg(mxring,2)
  INTEGER :: track(mxcp,mxnode)
  INTEGER :: trlen(mxcp)
  INTEGER :: ntr
  CHARACTER(LEN=mxlgr) :: tgroup(mxnode)
  INTEGER :: tbond(mxnode,mxnode)
  INTEGER :: ring(SIZE(bond,1)),rngflg
  LOGICAL :: ooalkfg          ! flag set to true if decomp is performed
  CHARACTER(LEN=9),PARAMETER :: progname='ooalkdec '
  CHARACTER(LEN=70)          :: mesg  
  
!-------------
! initialize:
!-------------
  j1=0  ;  j2=0  ;  j3=0  ;  ngr=0 ; ring(:)=0
  tchem(:)=' '   ;  tchem(1) = chem
  nchem=1  ;  sc(:)=0.  ;  sc(1) = 1.
  ooalkfg=.FALSE. ; dcom=' '

! find groups and bond matrix from chem      
  CALL grbond(chem,group,bond,dbflg,nring)
  IF (nring > 0) CALL rjgrm(nring,group,rjg)
  tbond(:,:)=bond ; tgroup(:)=group

  ngr=COUNT(group/=' ')

! j1 = radical group (if none, escape)
  DO i=1,ngr
    IF (INDEX(group(i),'.')/=0) THEN ; j1=i ; EXIT ; ENDIF
  ENDDO
  IF (j1==0) RETURN 
  IF (INDEX(group(j1),'O.')/=0) RETURN
  IF (INDEX(group(j1),'Cd.')/=0) RETURN
  
! check if -O--O-C. structure exist in chem (if none, escape)
  DO i=1,ngr ; IF (bond(j1,i)==3) THEN
    DO j=1,ngr ; IF ((bond(i,j)==3).AND.(tgroup(j)=='-O-').AND.(j/=j1)) THEN
      ooalkfg=.TRUE.
      EXIT
    ENDIF ; ENDDO
  ENDIF ; ENDDO
  IF (.NOT.ooalkfg) RETURN

  CALL gettrack(bond,j1,ngr,ntr,track,trlen)
  IF (nring > 0) CALL findring(1,2,ngr,tbond,rngflg,ring)
  trloop : DO i=1,ntr
    j2 = track(i,2) ; j3 = track(i,3) ; j4 = track(i,4)
    IF ((j2/=0).AND.(j3/=0)) THEN
      IF ((group(j2)=='-O-').AND.(group(j3)=='-O-')) THEN
	  
        IF      (tgroup(j1)=='CH2.')      THEN ; pold='CH2'      ; pnew='CH2O'
        ELSE IF (tgroup(j1)=='CH.')       THEN ; pold='CH'       ; pnew='CHO'
        ELSE IF (tgroup(j1)=='C.')        THEN ; pold='C'        ; pnew='CO'
        ELSE IF (tgroup(j1)(1:2)=='C(')   THEN ; pold='C('       ; pnew='CO('
        ELSE IF (tgroup(j1)(1:3)=='CH(')  THEN ; pold='CH('      ; pnew='CHO('
        ELSE  
          mesg='structure not handled in oolakdec'
          CALL stoperr(progname,mesg,chem)
        ENDIF
	  
        nc=INDEX(tgroup(j1),'.') ; tgroup(j1)(nc:nc)=' '
	  
        CALL swap(tgroup(j1),pold,tgr,pnew)
        tgroup(j1)=tgr
        tgroup(j2)=' '
        tgroup(j3)=' '
        k=INDEX(tgroup(j4),' ')  ;  tgroup(j4)(k:k+3)='(O.)'
        
        CALL setbond(tbond,j1,j2,0)
        CALL setbond(tbond,j2,j3,0)
        CALL setbond(tbond,j3,j4,0)
        
! write (non-standard) formula for new product
! Ring opening
        IF ((ring(j2)==1).AND.(ring(j3)==1)) THEN
          CALL rebond(tbond,tgroup,tchem(1),nring)
! linear decomposition
        ELSE
          CALL fragm(tbond,tgroup,tchem(1),tchem(2))
          sc(2)=1. ; nchem=2
        ENDIF
        EXIT trloop
      ENDIF
    ENDIF
  ENDDO trloop

  dcom='OOC._DEC'

END SUBROUTINE ooalkdec

!=======================================================================
! Purpose: this subroutine breaks the following structure: 
! -C.C1-O-C1-  ->  -C=C-C(O.)-
!=======================================================================
SUBROUTINE epoxalkdec(inchem,outchem,ecom)
  USE keyparameter, ONLY: mxnode,mxlgr,mxcp,mxring
  USE rjtool, ONLY:rjgrm
  USE stdgrbond, ONLY: grbond
  USE mapping, ONLY: gettrack
  USE reactool, ONLY: swap,rebond
  USE toolbox, ONLY: stoperr,setbond
  USE ringtool, ONLY: findring
  USE fragmenttool, ONLY: fragm
  USE normchem, ONLY: stdchm ! FOR DEBUG ONLY
  IMPLICIT NONE

  CHARACTER(LEN=*),INTENT(IN)  :: inchem   ! input species
  CHARACTER(LEN=*),INTENT(OUT) :: outchem  ! output species
  CHARACTER(LEN=*),INTENT(OUT) :: ecom     ! code for reference

  CHARACTER(LEN=mxlgr) :: group(mxnode)
  CHARACTER(LEN=mxlgr) :: tgr, pold, pnew
  INTEGER :: bond(mxnode,mxnode)
  INTEGER :: i,j,k,nring,dbflg
  INTEGER :: j1,ngr,nc
  INTEGER :: ec1,ec2,eo
  INTEGER :: rjg(mxring,2)
  CHARACTER(LEN=mxlgr) :: tgroup(mxnode)
  INTEGER :: tbond(mxnode,mxnode)
  LOGICAL :: epoxalkfg          ! flag set to true if decomp is performed
  CHARACTER(LEN=9),PARAMETER :: progname='ooalkdec '
  !CHARACTER(LEN=70)          :: mesg  
  
!-------------
! initialize:
!-------------
  j1=0 ; ngr=0 ; outchem=inchem ; epoxalkfg=.FALSE. ; ecom=' '

! find groups and bond matrix from chem      
  CALL grbond(inchem,group,bond,dbflg,nring)
  IF (nring > 0) CALL rjgrm(nring,group,rjg)
  tbond(:,:)=bond ; tgroup(:)=group

  ngr=COUNT(group/=' ')

! j1 = radical group (if none, escape)
  DO i=1,ngr
    IF (INDEX(group(i),'.')/=0) THEN ; j1=i ; EXIT ; ENDIF
  ENDDO
  IF (j1==0) RETURN 
  IF (INDEX(group(j1),'O.')/=0) RETURN
  IF (INDEX(group(j1),'Cd.')/=0) RETURN
  
! check if C1-O-C1-C. structure exist in chem (if none, escape)
  griloop : DO i=1,ngr 
    IF (bond(j1,i)==1) THEN
      DO j=1,ngr 
        IF ((bond(i,j)==3).AND.(tgroup(j)=='-O-')) THEN
          DO k=1,ngr
            IF ((bond(j,k)==3).AND.(tgroup(k)/='-O-').AND.(bond(i,k)==1).AND.(k/=j1)) THEN
              ec1=i ; eo=j ; ec2=k
              epoxalkfg=.TRUE.
              EXIT griloop
            ENDIF
          ENDDO
        ENDIF
      ENDDO
    ENDIF
  ENDDO griloop
  IF (.NOT.epoxalkfg) RETURN
  IF (INDEX(tgroup(j1) ,'(ONO2)')/=0) RETURN
  IF (INDEX(tgroup(ec1),'(ONO2)')/=0) RETURN

  pold='C'      ; pnew='Cd'
  CALL swap(tgroup(j1),pold,tgr,pnew)  ; tgroup(j1)=tgr
  nc=INDEX(tgroup(j1),'.') ; tgroup(j1)(nc:nc)=' '
  CALL swap(tgroup(ec1),pold,tgr,pnew) ; tgroup(ec1)=tgr
  tgroup(eo)=' ' 
  nc=INDEX(tgroup(ec2),' ') ; tgroup(ec2)(nc:nc+3)='(O.)'

  tbond(j1,ec1)=2 ; tbond(ec1,j1)=2
  tbond(ec1,eo)=0 ; tbond(eo,ec1)=0
  tbond(ec2,eo)=0 ; tbond(eo,ec2)=0
  
  CALL rebond(tbond,tgroup,outchem,nring)
  ecom='EPOX_DEC'
  !PRINT*, TRIM(inchem),' > ',TRIM(outchem)

END SUBROUTINE epoxalkdec


END MODULE radchktool
