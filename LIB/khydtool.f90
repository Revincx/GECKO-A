MODULE khydtool
IMPLICIT NONE
CONTAINS

! SUBROUTINE get_hydrate(chem,nwa,nhyd,chemhyd,yhyd,khydstar,nhydratdat,hydratdat,chemhydrat)
! SUBROUTINE khydration(chem,yield,ndat,chemdat,ydat)
! SUBROUTINE get_tsig(bond,group,ndeep,nabcde,tabcde,mapfun,funflg,nodetype,sigma)
! SUBROUTINE get_hsig(nabcde,tabcde,mapfun,funflg,nodetype,hsigma)
! SUBROUTINE kill_ester(tabester,nabcde,tabcde)

! ======================================================================
! PURPOSE: get all hydrates that can be made from a given formula and 
! related hydration constants. 
!
! For details about the Gromhe SAR and Khydration, see 
! T. Raventos-Duran et al., Atmospheric Chemistry and Physics, 7643-7654, 2010.
! C. Mouchel-Vallon et al., Geoscientific Model Development, 1339-1362, 2017 
! https://doi.org/10.5194/gmd-10-1339-2017, 
! 
! OUTPUT SUMMARY:
! - nwa: number of water molecule that can be added to chem. (e.g. if 
!     the molecule bear 3 ketones, the nwa is return as 3).
! - nhyd(i): number of distinct hydrate that can be made for i molecule
!     of water added.
! - chemhyd(j,i): formula of the jth hydrate bearing ith added water 
!     molecule
! - yhyd(j,i): hydration constant (with respect to the fully non-hydrated 
!     molec) of the jth hydrate bearing ith added water molecules (log scale)
! - khydstar: hydration constant taking all possible hydrate into
!     account (apparent constant K_star). NOT A LOG SCALE !
! ======================================================================
SUBROUTINE get_hydrate(chem,nwa,nhyd,chemhyd,yhyd,khydstar)
  USE database, ONLY: nkhydb,khydb_chem,khydb_dat ! database for hydrate
  USE searching, ONLY: srh5
  IMPLICIT NONE

  CHARACTER(LEN=*),INTENT(IN) :: chem
  INTEGER,INTENT(OUT) :: nwa                  ! # of water molecule that can be added to chem (= # of carbonyl).
  CHARACTER(LEN=*),INTENT(OUT) :: chemhyd(:,:)  ! (j:i) formula of the jth hydrate bearing ith added water molecule
  REAL,INTENT(OUT)    :: yhyd(:,:)   ! hydration constant (log scale)
  INTEGER,INTENT(OUT) :: nhyd(:)     ! # of distinct hydrate for i molecule of water added.
  REAL,INTENT(OUT)    :: khydstar    ! hydration constant taking all possible hydrate into account

  INTEGER :: i,j
  CHARACTER(LEN=LEN(chem)) :: t1chemhyd(SIZE(yhyd,2))
  REAL    :: t1yhyd(SIZE(yhyd,2))
  INTEGER :: t1nhyd

  CHARACTER(LEN=LEN(chem)) :: chemdat(SIZE(nhyd))  ! formula of the i'th hydrate
  CHARACTER(LEN=LEN(chem)) :: tempkc
  REAL    :: yield,numy,sumy
  REAL    :: ydat(SIZE(nhyd))  ! hydration constant of the i'th hydrate
  INTEGER :: ndat         ! # of distinct hydrate made from chem (one H2O added only !)
  INTEGER :: jp
  INTEGER :: hlev

! table of functions - Index of functionalities
! --------------------------
!  1= -OH    ;  2= -NO2     ;  3= -ONO2   ; 4= -OOH    ;  5= -F     
!  6= -Cl    ;  7=  -Br     ;  8= -I      ; 9= -CHO    ; 10= -CO-   
! 11= -COOH  ; 12= -CO(OOH) ; 13= -PAN    ; 14= -O-    ; 15= R-COO-R
! 16 = HCO-O-R; 17= -CO(F) ; 18= -CO(Cl) ;  19= -CO(Br) ; 20= -CO(I)

  nwa=0 ; khydstar=0. ; nhyd(:)=0 ; chemhyd(:,:)=' ' ; yhyd(:,:)=0.
      
! -----------
! monohydrate
! -----------
  t1nhyd=0  ;  t1chemhyd(:)=' '  ;  t1yhyd(:)=0.

  yield=0. ! parent compound
  CALL khydration(chem,yield,ndat,chemdat,ydat)
  IF (ndat==0) RETURN   ! no hydrate can be made from chem
  t1nhyd=ndat
  t1chemhyd(1:ndat)=chemdat(1:ndat)  ;  t1yhyd(1:ndat)=ydat(1:ndat)
      
! collapse identical formula
  DO i=1,t1nhyd-1
    numy=1  ;  sumy=t1yhyd(i)
    DO j=i+1,t1nhyd
      IF (t1chemhyd(i)==t1chemhyd(j)) THEN
        numy=numy+1  ;  sumy=sumy+t1yhyd(j)
        t1chemhyd(j)=' '  ;  t1yhyd(j)=0.
      ENDIF
    ENDDO
    IF (numy>1) t1yhyd(i)=sumy/REAL(numy)
  ENDDO      

! write the table for monohydrate
  DO i=1,t1nhyd
    IF (t1chemhyd(i)(1:1)/=' ') THEN
      nhyd(1)=nhyd(1)+1
      chemhyd(1,nhyd(1))=t1chemhyd(i)
      yhyd(1,nhyd(1))=t1yhyd(i)
    ENDIF
  ENDDO

! -----------
! multi-hydrate
! -----------
  multiloop: DO hlev=2,SIZE(nhyd)
    IF (nhyd(hlev-1)==0) THEN  ! all hydrate were found
      nwa=hlev-1
      EXIT multiloop
    ENDIF

    t1nhyd=0  ;  t1chemhyd(:)=' '  ;  t1yhyd(:)=0.
    DO jp=1,nhyd(hlev-1)
      yield=yhyd(hlev-1,jp) ! parent compound
      tempkc=chemhyd(hlev-1,jp)
      CALL khydration(tempkc,yield,ndat,chemdat,ydat)
      DO i=1,ndat
        t1chemhyd(t1nhyd+i)=chemdat(i)
        t1yhyd(t1nhyd+i)=ydat(i)
      ENDDO
      t1nhyd=t1nhyd+ndat
    ENDDO

! collapse identical formula
    DO i=1,t1nhyd-1
      numy=1  ;  sumy=t1yhyd(i)
      DO j=i+1,t1nhyd
        IF (t1chemhyd(i)==t1chemhyd(j)) THEN
          numy=numy+1  ;  sumy=sumy+t1yhyd(j)
          t1chemhyd(j)=' '  ;  t1yhyd(j)=0.
        ENDIF
      ENDDO
      IF (numy>1) t1yhyd(i)=sumy/REAL(numy)
    ENDDO      

! write the table for monohydrate
    DO i=1,t1nhyd
      IF (t1chemhyd(i)(1:1)/=' ') THEN
        nhyd(hlev)=nhyd(hlev)+1
        chemhyd(hlev,nhyd(hlev))=t1chemhyd(i)
        yhyd(hlev,nhyd(hlev))=t1yhyd(i)
      ENDIF
    ENDDO

  ENDDO multiloop

! Adjust nwa if exit did not occur in the formula loop
  IF (hlev>=SIZE(nhyd)) nwa=SIZE(nhyd)   

! compute Kstar
  DO i=1,nwa
    DO j=1,nhyd(i)
      khydstar=khydstar+10**(yhyd(i,j))
    ENDDO
  ENDDO

!-------------------------------------------------------------
!Check if hydration constant is already known in the database
!-------------------------------------------------------------
  i=srh5(chem,khydb_chem,nkhydb)
  IF (i>0) THEN
    khydstar = 10**khydb_dat(i,1)
  ENDIF

!  DO i=1,nkhydb
!    IF (chem == khydb_chem(i)) khydstar = 10**khydb_dat(i,1)
!  ENDDO

END  SUBROUTINE get_hydrate

! ======================================================================
! PURPOSE: return the list of hydrate (1 water molecule added) that
! can be made from the formula given as input
!
! OUTPUT SUMMARY:
! - ndat : number of distinct molecule (hydrate) that can be made from
!     chem (1 water molecule added only !)
! - chemdat(i) : formula of the i'th hydrate
! - ydat(i) : hydration constant of the i'th hydrate (with respect
!             to the fully non hydrated molecule !)
!
! INTERNAL
! - nabcde(k): number of distinct pathways that end up at a position k
!     relative to top (e.g. nabcde(3) gives the number of distinct 
!     pathways finishing in a beta position relative to top.
! - tabcde(k,i,j): give the pathways (nodes j), for the track number i
!     to reach the position k (k=2 is beta position ...). For example, 
!     tabcde(4,1,j) give the first track to reach a gamma position (node
!     given by tabcde(4,1,4), using the track given by tabcde(4,1,*).
! - mapfun(a,b,c): provide the number of function of type 'c' at position
!     (node) 'a'. index 'b' if for node type with 1=aliphatic, 2=cd and
!     3=aromatic. For example, the molecule CH2(OH)CdH=CdHCHO should 
!     have non zero values at the positions : mapfun(1,1,1)=1 and 
!     mapfun(4,2,9)=1
! - funflg(a): get the number of functional group at node a. For the 
!     example above, non-zero values are found at position 1 and 4, 
!     where it is set to 1.
! - nodetype : table of character for type node:
!      'y' = carbonyl      'r' = aromatic        'o'= -O- node
!      'd' = Cd            'n' = others (i.e. normal)
!
! table of functions - Index of functionalities (in mapfun)
! -----------------------------------------------------------------
!  1= -OH    ;  2= -NO2     ;  3= -ONO2   ; 4= -OOH    ;  5= -F
!  6= -Cl    ;  7=  -Br     ;  8= -I      ; 9= -CHO    ; 10= -CO-
! 11= -COOH  ; 12= -CO(OOH) ; 13= -PAN    ; 14= -O-    ; 15= R-COO-R
! 16 = HCO-O-R; 17= -CO(F) ; 18= -CO(Cl) ;  19= -CO(Br) ; 20= -CO(I)
! This list is used in alisig which provide sigma for a given group
! ======================================================================
SUBROUTINE khydration(chem,yield,ndat,chemdat,ydat)
  USE keyparameter, ONLY: mxnode,mxlgr,mxring,mxcp  ! need to create new bond and group matrix
  USE rjtool, ONLY: rjgrm
  USE stdgrbond, ONLY: grbond
  USE mapping, ONLY: abcde_map, chemmap
  USE reactool, ONLY: swap, rebond
  USE normchem, ONLY: stdchm
  USE toolbox, ONLY: stoperr
  IMPLICIT NONE

  CHARACTER(LEN=*),INTENT(IN) :: chem
  REAL,INTENT(IN) :: yield             ! Ratio (based on the non hydrated molecule)
  INTEGER,INTENT(OUT) :: ndat
  CHARACTER(LEN=LEN(chem)),INTENT(OUT) :: chemdat(:)
  REAL,INTENT(OUT) :: ydat(:)

! grbond and rjgrm variables
  INTEGER :: bond(mxnode,mxnode)
  CHARACTER(LEN=mxlgr) :: group(mxnode)
  INTEGER  :: nring, dbflg
  INTEGER  :: rjg(mxring,2)       ! ring-join group pairs

! local
  CHARACTER(LEN=LEN(chem)) :: tchemdat(SIZE(chemdat))
  REAL    :: tydat(SIZE(ydat))
  INTEGER :: tndat
  INTEGER :: node
  INTEGER :: i,j,k
  CHARACTER(LEN=LEN(chem)) :: tempkc
  CHARACTER(LEN=mxlgr) :: tgroup(mxnode), pold, pnew
  INTEGER :: tbond(mxnode,mxnode)
  REAL    :: sigma, hsigma
  REAL    :: khyd, sumy
  INTEGER :: numy

! chemmap variables 
  CHARACTER*1     nodetype(mxnode)
  REAL            alifun(21),cdfun(21),arofun(21)
  REAL            mapfun(mxnode,3,21)
  INTEGER         funflg(mxnode)
  INTEGER         tabester(4,2)  ! 1= -O- side, 2= CO side
  INTEGER         ngrp
  INTEGER         nfcd,nfcr
  INTEGER         ierr

! abcde_map variables
  INTEGER,PARAMETER :: mxdeep=7 ! sigma effect scroll up to position 7 (para group in aromatic)
  INTEGER :: nabcde(mxdeep), tabcde(mxdeep,mxcp,mxnode)

  CHARACTER(LEN=12),PARAMETER :: progname='khydration '
  CHARACTER(LEN=70)           :: mesg

  nfcr=0 ;  nfcd=0 ; ndat=0 ; tndat=0 ;  chemdat(:)=' ' 
  ydat(:)=0. ; tchemdat(:)=' ' ; tydat(:)=0.  

! build the group and bond matrix for chem      
  CALL grbond(chem,group,bond,dbflg,nring)
  DO i=1,mxnode  ;  bond(i,i)=0  ;  ENDDO
  CALL rjgrm(nring,group,rjg)  ! rm ring index and get ring closing nodes

! make copy 
  tgroup(:)=group(:)  ;  tbond(:,:)=bond(:,:)

! count the number of group in chem and get the map of the functional 
! group (see index number related the organic function)
  node=0
  DO  i=1,mxnode  ; IF (group(i)(1:1)/=' ') node = node + 1  ; ENDDO 
  CALL chemmap(chem,node,group,bond,ngrp,nodetype, &
               alifun,cdfun,arofun,mapfun,funflg,  &
               tabester,nfcd,nfcr,ierr)
  IF (ierr/=0) THEN
    mesg="--error--, after a chemmap call "
    CALL stoperr(progname,mesg,chem)
  ENDIF

! ---------------------------------
! hydration of aliphatic aldehyde 
! ---------------------------------
   IF ((alifun(9)/=0).OR.(cdfun(9)/=0) )THEN

! find aldehyde position
    alidloop: DO i=1,node
      IF ( (mapfun(i,1,9)==1.).OR.(mapfun(i,2,9)==1.) ) THEN
        sigma=0.
        CALL abcde_map(bond,i,node,nabcde,tabcde)  ! get the neighboors.

! check that the aldehyde is not an hidden anhydride (-CO-O-CO-) 
        DO k=1,nabcde(2)
          IF (nodetype(tabcde(2,k,2))=='o') CYCLE alidloop 
        ENDDO   

! Remove the path that would lead to count the ester twice. 
! Performed by "killing" the "second" node of the functionality.
        CALL kill_ester(tabester,nabcde,tabcde)

! get taft sigma (look for groups up to delta position)
        CALL get_tsig(bond,group,5,nabcde,tabcde,mapfun, &
                      funflg,nodetype,sigma)

! compute the hydration constant for aliphatic aldehyde 
        khyd=0.0818 + 1.2663*sigma

! If the species is a 'Cd' aldehyde, then add the correction applied
! for aromatic (no data available)
        IF  (mapfun(i,2,9)==1.) khyd=khyd-1.5817

! make the hydrate, rebuild, rename and reset tgroup & tbond
        pold = 'CHO '
        pnew = 'CH(OH)(OH) '
        CALL swap(group(i),pold,tgroup(i),pnew)
        CALL rebond(tbond,tgroup,tempkc,nring)
        CALL stdchm(tempkc)
        tgroup(:)=group(:)  ;  tbond(:,:)=bond(:,:)

! store the data 
        tndat=tndat+1
        IF (tndat>SIZE(ydat)) THEN
          mesg="# exceed max hydrate "
          CALL stoperr(progname,mesg,chem)
        ENDIF
        tydat(tndat)=khyd+yield
        tchemdat(tndat)=tempkc

      ENDIF
    ENDDO alidloop
  ENDIF  ! end of aliphatic aldehyde

! ---------------------------------
! hydration of aliphatic ketone 
! ---------------------------------
  IF (alifun(10)/=0) THEN

! find ketone position
    alikloop: DO i=1,node
      IF (mapfun(i,1,10)==1.) THEN
        sigma=0.

! get the neighboors.
        CALL abcde_map(bond,i,node,nabcde,tabcde)  

! check that the ketone is not an hidden anhydride (-CO-O-CO-) 
        DO k=1,nabcde(2)
          IF (nodetype(tabcde(2,k,2))=='o') CYCLE alikloop
        ENDDO   

! Remove the path that would lead to count the ester twice. 
! Performed by "killing" the "second" node of the functionality.
        CALL kill_ester(tabester,nabcde,tabcde)

! get taft sigma (look for groups up to delta position)
        CALL get_tsig(bond,group,5,nabcde,tabcde,mapfun,&
                      funflg,nodetype,sigma)

! compute the hydration constant for aliphatic ketone 
        khyd=0.0818 + 1.2663*sigma - 2.5039

! make the hydrate, rebuild, rename and reset tgroup & tbond
        pold = 'CO '  ;  pnew = 'C(OH)(OH) '
        CALL swap(group(i),pold,tgroup(i),pnew)
        CALL rebond(tbond,tgroup,tempkc,nring)
        CALL stdchm(tempkc)
        tgroup(:)=group(:)  ;  tbond(:,:)=bond(:,:)
       
! store the data 
        tndat=tndat+1
        IF (tndat>SIZE(ydat)) THEN
          mesg="# exceed max hydrate "
          CALL stoperr(progname,mesg,chem)
        ENDIF
        tydat(tndat)=khyd+yield
        tchemdat(tndat)=tempkc

      ENDIF
    ENDDO alikloop  
  ENDIF  ! end of aliphatic ketone

      
! ---------------------------------
! hydration of aromatic aldehyde 
! ---------------------------------
  IF (arofun(9)/=0) THEN

! find aldehyde position
    DO i=1,node
      IF (mapfun(i,3,9)==1.) THEN
        sigma=0.

        CALL abcde_map(bond,i,node,nabcde,tabcde)  ! get the neighboors

! Remove the path that would lead to count the ester twice. 
! Performed by "killing" the "second" node of the functionality.
        CALL kill_ester(tabester,nabcde,tabcde)

! get hammet sigma (look for groups up to delta position)
        CALL get_hsig(nabcde,tabcde,mapfun,funflg,nodetype,hsigma)

! compute the hydration constant for alihpatic aldehyde 
        khyd=0.0818 + 0.4969*hsigma - 1.5817

! make the hydrate, rebuild, rename and reset tgroup & tbond
        pold = 'CHO '  ;  pnew = 'CH(OH)(OH) '
        CALL swap(group(i),pold,tgroup(i),pnew)
        CALL rebond(tbond,tgroup,tempkc,nring)
        CALL stdchm(tempkc)
        tgroup(:)=group(:)  ;  tbond(:,:)=bond(:,:)
        
! store the data 
        tndat=tndat+1
        IF (tndat>SIZE(ydat)) THEN
          mesg="# exceed max hydrate "
          CALL stoperr(progname,mesg,chem)
        ENDIF
        tydat(tndat)=khyd+yield
        tchemdat(tndat)=tempkc
      ENDIF
    ENDDO       
  ENDIF  ! end of aromatic aldehyde

! ----------------------------------------------
! collapse identical product 
! ----------------------------------------------
! note : must add khyd, not the log of it
  DO i=1,tndat
    numy=1
    sumy=10**(tydat(i))
    DO j=i+1,tndat
      IF (tchemdat(i)==tchemdat(j)) THEN
        numy=numy+1
        sumy=sumy+10**(tydat(j))
        tchemdat(j)=' '
      ENDIF
    ENDDO
    IF (numy>1) THEN
      tydat(i)=log10(sumy)
    ENDIF
  ENDDO

  DO i=1,tndat
    IF (tchemdat(i)(1:1)/=' ') THEN
      ndat=ndat+1
      chemdat(ndat)=tchemdat(i)
      ydat(ndat)=tydat(i)
    ENDIF
  ENDDO
      
END SUBROUTINE khydration

! ======================================================================
! Purpose: Return the taft sigma value at for the positions given by the 
! of the functional group given in "tabcde". Note : the position for 
! which sigma is computed is given by tabcde(1,1,1).
!
! - ndeep: deepest position for which sigma must be considered (ndeep=5 
!   is delta position)
! - nabcde(k): number of distinct pathways that end up at a position k
!     relative to top (e.g. nabcde(3) gives the number of distinct 
!     pathways finishing in a beta position relative to top.
! - tabcde(k,i,j): give the pathways (nodes j), for the track number i
!     to reach the position k (k=2 is beta position ...). For example, 
!     tabcde(4,1,j) give the first track to reach a gamma position (node
!     given by tabcde(4,1,4), using the track given by tabcde(4,1,*).
! - mapfun(a,b,c): provide the number of function of type 'c' at position
!     (node) 'a'. index 'b' if for node type with 1=aliphatic, 2=cd and
!     3=aromatic. For example, the molecule CH2(OH)CdH=CdHCHO should 
!     have non zero values at the positions : mapfun(1,1,1)=1 and 
!     mapfun(4,2,9)=1
! - funflg(a): get the number of functional group at node a. For the 
!     example above, non-zero values are found at position 1 and 4, 
!     where it is set to 1.
! - nodetype: table of character for type node:
!      'y' = carbonyl      'r' = aromatic        'o'= -O- node
!      'd' = Cd            'n' = others (i.e. normal)
! - sigma: is the taft sigma due to the neighboors relative to the
!     position being considered (and given by tabcde(1,1,1).
!
! table of functions - Index of functionalities (in mapfun)
! -----------------------------------------------------------------
!  1= -OH    ;  2= -NO2     ;  3= -ONO2   ; 4= -OOH    ;  5= -F     
!  6= -Cl    ;  7=  -Br     ;  8= -I      ; 9= -CHO    ; 10= -CO-   
! 11= -COOH  ; 12= -CO(OOH) ; 13= -PAN    ; 14= -O-    ; 15= R-COO-R
! 16 = HCO-O-R; 17= -CO(F) ; 18= -CO(Cl) ;  19= -CO(Br) ; 20= -CO(I)
! 21 = -CO(O-)
! This list is used in alisig which provide sigma for a given group
! ======================================================================
SUBROUTINE get_tsig(bond,group,ndeep,nabcde,tabcde,mapfun, &
                    funflg,nodetype,sigma)
  IMPLICIT NONE

  INTEGER,INTENT(IN) :: bond(:,:)
  CHARACTER(LEN=*),INTENT(IN) :: group(:)
  INTEGER,INTENT(IN) :: ndeep
  INTEGER,INTENT(IN) :: nabcde(:)
  INTEGER,INTENT(IN) :: tabcde(:,:,:)
  REAL,INTENT(IN)    :: mapfun(:,:,:)
  INTEGER,INTENT(IN) :: funflg(:)
  CHARACTER*1,INTENT(IN) :: nodetype(:)
  REAL,INTENT(OUT) :: sigma

  INTEGER         idp, k, ncd, locn, l, i1, i2, ifun
  REAL            dist, tsig
  INTEGER         aroflg

! Tafta sigma values for aliphatic compounds. sigma for ester depends
! whether ester is connected from the -CO- side or the -O- side. 
! Sigma=2.56 for -O- side and 2.00 for CO side.
! CMV 16/06/14 : add sigma taft for CO(O-), n°21
  REAL,DIMENSION(21),PARAMETER :: alisig(21)=  &
    (/ 0.62, 1.47, 1.38, 0.62, 1.10, 0.94, 1.00, 1.00, 2.15, 1.81, &
       2.08, 2.08, 2.00, 1.81, 2.56, 2.90, 2.44, 2.37, 2.37, 2.37, &
      -1.06/)
  REAL, PARAMETER :: sigester=2.00 ! sigma for ester, CO side

  IF (ndeep>SIZE(nabcde)) THEN
    WRITE(6,*) '--error-- in get_tsig, ndeep>SIZE(nabcde)=',SIZE(nabcde)
    STOP "in get_tsig"
  ENDIF
  
  sigma=0.
! scroll the neighboors up to ndeep position 
  DO idp=2,ndeep
    DO k=1,nabcde(idp)
      ncd=0
      locn=tabcde(idp,k,idp)
      IF (locn == 0) CYCLE
      IF (funflg(locn)/=0) THEN

        DO l=1,idp-1 ! start the loop to see if Cd are in between
          i1=tabcde(idp,k,l)   ; IF (i1==0) CYCLE;
          i2=tabcde(idp,k,l+1) ; IF (i2==0) CYCLE
          IF (bond(i1,i2)==2) ncd=ncd+1
          IF ( (nodetype(i1)=='d').AND.(nodetype(i2)=='d') ) ncd=ncd+1
        ENDDO
        dist=REAL(idp)-REAL(ncd)-2.
        IF (group(locn)(1:3)=='-O-') dist=dist-1.  ! decrease the distance
   
        aroflg=0               ! check for aromatic nodes between groups
        DO l=2,idp
          IF (tabcde(idp,k,l)==0) CYCLE
          IF (nodetype(tabcde(idp,k,l))=='r') aroflg=1
        ENDDO

        DO ifun=1,21       ! find the function on node locn (local node)
          IF (mapfun(locn,1,ifun)/=0) THEN
            tsig=mapfun(locn,1,ifun)*alisig(ifun)
            IF (ifun==15) THEN                           ! check ester
              IF (group(locn)(1:2)=='CO') tsig=sigester  ! CO instead of -O-
            ENDIF
            sigma=sigma+(tsig*0.4**(dist)) 
          ENDIF

          IF (mapfun(locn,2,ifun)/=0) THEN     ! =CdC(X)CHO be must counted
            tsig=mapfun(locn,2,ifun)*alisig(ifun)
            IF (ifun==15) THEN                           ! check ester
              IF (group(locn)(1:2)=='CO') tsig=sigester  ! CO instead of -O-
            ENDIF
            sigma=sigma+(tsig*0.4**(dist)) 
          ENDIF

          IF (mapfun(locn,3,ifun)/=0) THEN    ! Phi-COCHO must be counted
            IF (aroflg==0) THEN  
              tsig=mapfun(locn,3,ifun)*alisig(ifun)
              IF (ifun==15) THEN                           ! check ester
                IF (group(locn)(1:2)=='CO') tsig=sigester  ! CO instead of -O-
              ENDIF
              sigma=sigma+(tsig*0.4**(dist)) 
            ENDIF
          ENDIF
        ENDDO

      ENDIF
    ENDDO
  ENDDO

END SUBROUTINE get_tsig

! ======================================================================
! Purpose: Return the hammet sigma (see header for get_tsig above for)
! group numbers. arorto, arometa, aropara give the hammet sigma for 
! each group in o,m,p.
! ======================================================================
SUBROUTINE get_hsig(nabcde,tabcde,mapfun,funflg,nodetype,hsigma)
  IMPLICIT NONE

  INTEGER,INTENT(IN) :: nabcde(:)
  INTEGER,INTENT(IN) :: tabcde(:,:,:)
  REAL,INTENT(IN)    :: mapfun(:,:,:)
  INTEGER,INTENT(IN) :: funflg(:)
  CHARACTER(LEN=*),INTENT(IN) :: nodetype(:)
  REAL,INTENT(OUT)   :: hsigma

  INTEGER :: idp,k,locn,l,ifun,funindex
  REAL    :: dist

! CMV - 16/06/14. Add hammet values for -CO(O-)
! benzoic as a reference for substituents in ortho
  REAL, DIMENSION(21),PARAMETER :: arometa = &
    (/ 0.13, 0.74, 0.55, 0.00, 0.34, 0.37, 0.39, 0.35, 0.36, 0.36, &
       0.35, 0.00, 0.00, 0.11, 0.32, 0.00, 0.55, 0.53, 0.53, 0.53, &
       0.09 /)
  REAL, DIMENSION(21),PARAMETER :: aropara = &
    (/ -0.38, 0.78, 0.70, 0.00, 0.06, 0.24, 0.22, 0.21, 0.44, 0.47, &
        0.44, 0.00, 0.00,-0.28, 0.39, 0.00, 0.70, 0.69, 0.69, 0.69, &
       -0.05 /)
  REAL, DIMENSION(21),PARAMETER :: arorto = &
    (/  1.22, 1.99, 0.00, 0.00, 0.93, 1.28, 1.35, 1.34, 0.72, 0.07, &
        0.95, 0.00, 0.00, 0.12, 0.63, 0.00, 0.00, 0.00, 0.00, 0.00, &
       -0.91 /)
  
  hsigma=0.  ;  funindex = 0
! scroll the neighboors up to 7 position (i.e. a carbonyl on para)
  DO idp=2,7
    kloop: DO k=1,nabcde(idp)

! exit if already taken into account by another pathway (it happens at 
! every para position)
      IF (k>=2) THEN
        DO l=1,k-1
          IF (tabcde(idp,k,idp)==tabcde(idp,l,idp) ) CYCLE kloop 
        ENDDO
      ENDIF

      locn=tabcde(idp,k,idp)
      IF (funflg(locn)/=0) THEN

! find the function on node locn (local node)
        DO ifun=1,21  
!CMV - 20/06/16. Look for function on aromatic AND aliphatic carbons
          IF ((mapfun(locn,3,ifun)/=0).OR.(mapfun(locn,1,ifun)/=0)) funindex=ifun
        ENDDO
            
!CMV - 20/06/16. funindex should not be 0 at this stage
        IF (funindex == 0) THEN
            PRINT*,"error, funindex = 0 in get_hsig"
            STOP "in get_hsig"
        ENDIF

! find the distance (ortho, meta, para) on node locn (local node)
        dist=0
        DO l=1,idp
          IF (nodetype(tabcde(idp,k,l))=='r') dist=dist+1
        ENDDO

        IF (dist==2) THEN 
          hsigma = hsigma+arorto(funindex)
        ELSE IF (dist==3) THEN
          hsigma = hsigma+arometa(funindex)
        ELSE IF (dist==4) THEN
          hsigma = hsigma+aropara(funindex)
        ENDIF
        IF (tabcde(idp,k,2)==0) hsigma=0. ! kill second node on ester
      ENDIF
    ENDDO kloop
  ENDDO

END SUBROUTINE get_hsig

! ======================================================================
! Purpose: Remove the path that would lead to count the ester twice. 
! Performed by "killing" the "second" node of the functionality.
!
! tabester: provide the position of ester "couple" (i.e. the O and CO 
! nodes. For example, the molecule CH3CO-O-CH2-O-COCH3 has the following
! values: tabester(1,1)=3,tabester(1,2)=2 
!         tabester(2,1)=5,tabester(2,2)=6
! nabcde(k): number of distinct pathways that end up at a position k 
! relative to top (e.g. nabcde(3) gives the number of distinct pathways
! finishing in a beta position relative to top.                   
! tabcde(k,i,j) : give the pathways (node j), for the track number i  
! to reach the position k (k=2 is beta position ...). For example, 
! tabcde(4,1,j) give the first track to reach a gamma position (node 
! given by tabcde(4,1,4), using the track given by tabcde(4,1,*).
! ======================================================================
SUBROUTINE kill_ester(tabester,nabcde,tabcde)
  IMPLICIT NONE

  INTEGER,INTENT(IN) :: nabcde(:)
  INTEGER,INTENT(IN) :: tabester(:,:)  ! (:,2) 1= -O- side, 2= CO side
  INTEGER,INTENT(INOUT) :: tabcde(:,:,:)

  INTEGER k, ipos, it, j

  DO k=1,SIZE(tabester,1)
    IF (tabester(k,1)/=0) THEN

      DO ipos=2,SIZE(nabcde)
        DO it=1,nabcde(ipos)
          DO j=1,ipos-1

            IF (tabcde(ipos,it,j)==tabester(k,1)) THEN
              IF (tabcde(ipos,it,j+1)==tabester(k,2)) THEN
                tabcde(ipos,it,j+1)=0  ! set the node to 0
              ENDIF 
            ENDIF 

            IF (tabcde(ipos,it,j)==tabester(k,2)) THEN
              IF (tabcde(ipos,it,j+1)==tabester(k,1)) THEN
                tabcde(ipos,it,j+1)=0  ! set the node to 0
              ENDIF 
            ENDIF 

          ENDDO
        ENDDO
      ENDDO

    ENDIF
  ENDDO

END SUBROUTINE kill_ester

END MODULE khydtool
