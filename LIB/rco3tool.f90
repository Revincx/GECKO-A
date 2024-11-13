MODULE rco3tool
IMPLICIT NONE
CONTAINS

!=======================================================================
! PURPOSE: perform the RCO3+NO reaction and write the reaction. Only
! one reaction is considered : RCO3+NO => RCO(O.) + NO2
!=======================================================================
SUBROUTINE rco3no(chem,pname,bond,group,brch,ip)
  USE keyparameter, ONLY: mxlco,mxpd,mecu,kno3u
  USE references, ONLY:mxlcod
  USE reactool, ONLY: swap,rebond
  USE rxwrttool, ONLY:rxwrit,rxinit
  USE toolbox, ONLY: stoperr,add1tonp,addrx,addref
  USE tempo, ONLY: chknadd2p

  CHARACTER(LEN=*),INTENT(IN) :: chem  ! formula of the VOC considered
  CHARACTER(LEN=*),INTENT(IN) :: pname ! short name of chem (parent compound) 
  CHARACTER(LEN=*),INTENT(IN) :: group(:) ! group (without ring joining char)
  INTEGER,INTENT(IN) :: bond(:,:)     ! bond matrix
  REAL,INTENT(IN)    :: brch          ! max yield of the input species
  INTEGER,INTENT(IN) :: ip            ! node bearing the peroxy (OO.) group

  CHARACTER(LEN=LEN(group(1))) :: tgroup(SIZE(group))
  CHARACTER(LEN=LEN(group(1))) :: pold,pnew
  INTEGER         :: tbond(SIZE(bond,1),SIZE(bond,2))
  CHARACTER(LEN=LEN(chem)) :: tempfo

  CHARACTER(LEN=mxlco) :: p(mxpd), r(3)
  REAL    :: s(mxpd), arrh(3)
  INTEGER :: idreac, nlabel
  REAL    :: xlabel,folow(3),fotroe(4)

  INTEGER :: tempring, np
  REAL    :: brtio, yield

  INTEGER,PARAMETER     :: mxref=10
  INTEGER               :: nref
  CHARACTER(LEN=mxlcod) :: ref(mxref)
  CHARACTER(LEN=7)      :: progname='rco3no'

! set the reaction rate (Jenkin et al., ACP, 2019, MJ19KMV000 eqn (1)) 
  REAL,DIMENSION(3),PARAMETER :: kno=(/7.5E-12, 0., -290. /)

  tgroup(:)=group(:) ;  tbond(:,:)=bond(:,:)

  CALL rxinit(r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe)
  nref=0  ;  ref(:)=' ' ;  CALL addref(progname,'MJ19KMV000',nref,ref,chem)
  np=0
  
! change (OO.) to (O.), rebuild, check and rename
  pold='(OO.)' ;  pnew='(O.)'
  CALL swap(group(ip),pold,tgroup(ip),pnew)
  CALL rebond(tbond,tgroup,tempfo,tempring)
  yield=1.;   brtio=brch

! check new formula, standardize and add products and coproducts to p(:)   
  CALL chknadd2p(chem,tempfo,yield,brtio,np,s,p,nref,ref)

! add NO2 to the product list
  CALL add1tonp(progname,chem,np)   
  s(np)=1.  ;  p(np)='NO2 '

! write the reaction
  r(1)=pname  ;  r(2)='NO ' ; arrh(:)=kno(:)
  CALL rxwrit(mecu,r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe,ref)

END SUBROUTINE rco3no


!=======================================================================
! PURPOSE: perform the RCO3+NO2 reaction and write the reaction. The
! reaction leads to PAN (RCO(OONO2)), except if a carbonyl is next to
! peroxy acyl group: RCOCO3+NO2 => RCOCO(O.)+NO3 (as described in 
! Jenkin et al., ACP, 2019).
!=======================================================================
SUBROUTINE rco3no2(chem,pname,bond,group,ngr,brch,ip)
  USE keyparameter, ONLY: mxlco,mxpd,mecu
  USE references, ONLY:mxlcod
  USE atomtool, ONLY: cnum
  USE reactool, ONLY: swap,rebond
  USE dictstacktool, ONLY: bratio
  USE normchem, ONLY: stdchm  
  USE rxwrttool, ONLY:rxwrit,rxinit
  USE toolbox, ONLY: stoperr,add1tonp,addrx,addref
  USE tempo, ONLY: chknadd2p

  CHARACTER(LEN=*),INTENT(IN) :: chem  ! formula of the VOC considered
  CHARACTER(LEN=*),INTENT(IN) :: pname ! short name of chem (parent compound) 
  CHARACTER(LEN=*),INTENT(IN) :: group(:) ! group (without ring joining char)
  INTEGER,INTENT(IN) :: bond(:,:)     ! bond matrix
  INTEGER,INTENT(IN) :: ngr           ! # of nodes
  REAL,INTENT(IN)    :: brch          ! max yield of the input species
  INTEGER,INTENT(IN) :: ip            ! node bearing the peroxy (OO.) group

  CHARACTER(LEN=LEN(group(1))) :: tgroup(SIZE(group))
  CHARACTER(LEN=LEN(group(1))) :: pold,pnew
  INTEGER         :: tbond(SIZE(bond,1),SIZE(bond,2))
  CHARACTER(LEN=LEN(chem)) :: tempfo

  CHARACTER(LEN=mxlco) :: p(mxpd), r(3)
  REAL    :: s(mxpd), arrh(3)
  INTEGER :: idreac, nlabel
  REAL    :: xlabel,folow(3),fotroe(4)

  INTEGER :: tempring, np, i, nca
  REAL    :: brtio, yield
  LOGICAL :: alkpath

  INTEGER,PARAMETER     :: mxref=10
  INTEGER               :: nref
  CHARACTER(LEN=mxlcod) :: ref(mxref)
  CHARACTER(LEN=7)      :: progname='rco3no2'

! Set the reaction rate (Jenkin et al., ACP, 2019). 
! For the fall off reaction, rate constant are given with the following
! expression : k=arrh1*(T/300)**arrh2*exp(-arrh3/T)
  REAL,DIMENSION(3),PARAMETER :: kno2high=(/1.125E-11 , -1.105 , 0. /)                  ! fall off exp.
  REAL,DIMENSION(3),PARAMETER :: kno2=(/1.125E-11*((1./300.)**(-1.105)), -1.105 , 0. /) ! usual arrh exp.

  tgroup(:)=group(:) ;  tbond(:,:)=bond(:,:)
  nca=cnum(chem)

! check for R-CO-CO(OO.) structure
  alkpath=.FALSE.
  DO i=1,ngr
    IF (bond(i,ip)==1) THEN
      IF ((group(i)(1:3)=='CO ').OR.(group(i)(1:3)=='CHO')) THEN
        alkpath=.TRUE.
        EXIT
      ENDIF
    ENDIF
  ENDDO
  
! R-CO-CO(OO.) + NO2 -> R-CO-CO(O.) + NO3  (see Jenkin et al., 2019)
  IF (alkpath) THEN

    CALL rxinit(r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe)
    nref=0  ;  ref(:)=' ' ;  CALL addref(progname,'MJ19KMV000',nref,ref,chem)
    CALL addref(progname,'RCOCO3_NO2',nref,ref,chem)
    np=0

! change (OO.) to (O.), rebuild, check and rename
    pold='(OO.)' ;  pnew='(O.)'
    CALL swap(group(ip),pold,tgroup(ip),pnew)
    CALL rebond(tbond,tgroup,tempfo,tempring)

! check new formula, standardize and add products and coproducts to p(:)   
    yield=1.;   brtio=brch
    CALL chknadd2p(chem,tempfo,yield,brtio,np,s,p,nref,ref)

! add NO2 to the product list
    CALL add1tonp(progname,chem,np)   
    s(np)=1.  ;  p(np)='NO3 '

! write the reaction
    r(1)=pname  ;  r(2)='NO2 '
    arrh(:)=kno2(:)
    CALL rxwrit(mecu,r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe,ref)
     
  ELSE

    CALL rxinit(r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe)
    nref=0  ;  ref(:)=' '
    np=0

! change (OO.) to (OONO2), rebuild and rename:
    pold='(OO.)'  ;   pnew='(OONO2)'
    CALL swap(group(ip),pold,tgroup(ip),pnew)
    CALL rebond(tbond,tgroup,tempfo,tempring)
    CALL stdchm(tempfo)

    s(1)=1.0  ;  CALL bratio(tempfo,brch,p(1),nref,ref)

    r(1)=pname  ;  r(2)='NO2 '      ! reactant
   
    IF (nca==2) THEN
! rate constant (fall off regime) from IUPAC, last access 2019
      idreac=3                                           ! fall off reaction (idreac=3)
      folow(1)=3.28E-28  ; folow(2)=-6.87  ; folow(3)=0. ! low pressure rate constant
      arrh(1)=1.125E-11 ; arrh(2)=-1.105 ; arrh(3)=0.    ! high pressure rate constant
      fotroe(1)=0.3                                      ! broadening factor

      r(3)='(+M)'
      CALL addref(progname,'IUPAKMV000',nref,ref,chem)
      CALL rxwrit(mecu,r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe,ref)
      
    ELSE
      arrh(:)=kno2(:)
      CALL addref(progname,'MJ19KMV000',nref,ref,chem)
      CALL rxwrit(mecu,r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe,ref)

    ENDIF
  ENDIF
END SUBROUTINE rco3no2

!=======================================================================
! PURPOSE: perform the RCO3+NO3 reaction and write the reaction. Only
! one reaction is considered: RCO3+NO3 => RCO(O.) + NO2 (+O2)
!=======================================================================
SUBROUTINE rco3no3(chem,pname,bond,group,brch,ip)
  USE keyparameter, ONLY: mxlco,mxpd,mecu
  USE references, ONLY:mxlcod
  USE reactool, ONLY: swap,rebond
  USE rxwrttool, ONLY:rxwrit,rxinit
  USE toolbox, ONLY: stoperr,add1tonp,addrx,addref
  USE tempo, ONLY: chknadd2p

  CHARACTER(LEN=*),INTENT(IN) :: chem  ! formula of the VOC considered
  CHARACTER(LEN=*),INTENT(IN) :: pname ! short name of chem (parent compound) 
  CHARACTER(LEN=*),INTENT(IN) :: group(:) ! group (without ring joining char)
  INTEGER,INTENT(IN) :: bond(:,:)      ! bond matrix
  REAL,INTENT(IN)    :: brch           ! max yield of the input species
  INTEGER,INTENT(IN) :: ip             ! node bearig the peroxy (OO.) group

  CHARACTER(LEN=LEN(group(1))) :: tgroup(SIZE(group))
  CHARACTER(LEN=LEN(group(1))) :: pold,pnew
  INTEGER         :: tbond(SIZE(bond,1),SIZE(bond,2))
  CHARACTER(LEN=LEN(chem)) :: tempfo

  CHARACTER(LEN=mxlco) :: p(mxpd), r(3)
  REAL    :: s(mxpd), arrh(3)
  INTEGER :: idreac, nlabel
  REAL    :: xlabel,folow(3),fotroe(4)

  INTEGER :: tempring, np
  REAL    :: brtio, yield

  INTEGER,PARAMETER     :: mxref=10
  INTEGER               :: nref
  CHARACTER(LEN=mxlcod) :: ref(mxref)
  CHARACTER(LEN=7)      :: progname='rco3no3'

! set the reaction rate (Jenkin et al., ACPD, 2019) 
  REAL,DIMENSION(3),PARAMETER :: kno3=(/8.9E-12, 0., 305. /)

  tgroup(:)=group(:) ;  tbond(:,:)=bond(:,:)

  CALL rxinit(r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe)
  nref=0  ;  ref(:)=' ' ;  CALL addref(progname,'MJ19KMV000',nref,ref,chem)
  np=0
  
! change (OO.) to (O.), rebuild, check and rename
  pold='(OO.)' ;  pnew='(O.)'
  CALL swap(group(ip),pold,tgroup(ip),pnew)
  CALL rebond(tbond,tgroup,tempfo,tempring)

! check new formula, standardize and add products and coproducts to p(:)   
  yield=1.;   brtio=brch*0.01
  CALL chknadd2p(chem,tempfo,yield,brtio,np,s,p,nref,ref)

! add NO2 to the product list
  CALL add1tonp(progname,chem,np)   
  s(np)=1.  ;  p(np)='NO2 '

! write the reaction
  r(1)=pname  ;  r(2)='NO3 ' ; arrh(:)=kno3(:)
  CALL rxwrit(mecu,r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe,ref)

END SUBROUTINE rco3no3


!=======================================================================
! PURPOSE: perform and write the RCO3+HO2 reactions. 
! 
! Following channels are considered for the reactions of RCO3 with HO2, 
! as discussed in Jenkin et al., ACP, 2019:
! RCO3+HO2 -> RCO(OOH)+O2     (a) INDEX 1 
! RCO3+HO2 -> RCO(OH)+O3      (b) INDEX 2 
! RO2+HO2  -> R(-H)O+H2O+O2   (c) INDEX 3  (only for RO2, not used here)
! RCO3+HO2 -> RCO(O.)+OH+O2   (d) INDEX 4
! RO2+HO2  -> R(-H)O+OH+HO2   (e) INDEX 5  (only for RO2, not used here)
!
! Four distinct structures are considered to set branching ratios:
!   "regular" RCO3 and "aromatic RCO3  
!=======================================================================
SUBROUTINE rco3ho2(chem,pname,bond,group,ngr,brch,ip)
  USE keyparameter, ONLY: mxlco,mxpd,mecu
  USE references, ONLY:mxlcod
  USE atomtool, ONLY: getatoms
  USE reactool, ONLY: swap,rebond
  USE dictstacktool, ONLY: bratio
  USE normchem, ONLY: stdchm  
  USE rxwrttool, ONLY:rxwrit,rxinit
  USE toolbox, ONLY: stoperr,add1tonp,addrx,addref
  USE tempo, ONLY: chknadd2p

  CHARACTER(LEN=*),INTENT(IN) :: chem  ! formula of the VOC considered
  CHARACTER(LEN=*),INTENT(IN) :: pname ! short name of chem (parent compound) 
  CHARACTER(LEN=*),INTENT(IN) :: group(:) ! group (without ring joining char)
  INTEGER,INTENT(IN) :: bond(:,:)      ! bond matrix
  INTEGER,INTENT(IN) :: ngr            ! # of nodes
  REAL,INTENT(IN)    :: brch           ! max yield of the input species
  INTEGER,INTENT(IN) :: ip             ! node bearig the peroxy (OO.) group

  CHARACTER(LEN=LEN(group(1))) :: tgroup(SIZE(group))
  CHARACTER(LEN=LEN(group(1))) :: pold,pnew
  INTEGER         :: tbond(SIZE(bond,1),SIZE(bond,2))
  CHARACTER(LEN=LEN(chem)) :: tempfo

  CHARACTER(LEN=mxlco) :: p(mxpd), r(3)
  REAL    :: s(mxpd), arrh(3)
  INTEGER :: idreac, nlabel
  REAL    :: xlabel,folow(3),fotroe(4)

  INTEGER :: np,i,tempring
  INTEGER :: ic,ih,in,io,ir,is,ifl,icl,ibr,Ncon
  REAL    :: brtio, yield
  REAL    :: kho2(3)
  REAL    :: patratio(5)               ! path ratio for channel (a)-(e) 

  INTEGER,PARAMETER     :: mxref=10
  INTEGER               :: nref
  CHARACTER(LEN=mxlcod) :: ref(mxref)
  CHARACTER(LEN=7)      :: progname='rco3ho2'

  tgroup(:)=group(:) ;  tbond(:,:)=bond(:,:)

! Rates from Jenkin et al. 2019
  CALL getatoms(chem,ic,ih,in,io,ir,is,ifl,ibr,icl)
  Ncon=ic+io-2+in
! set the rate constant (see Jenkin, ACP, 2019)
  kho2(1)=3.5E-12*(1-exp(-0.23*REAL(Ncon))) ; kho2(2)=0. ; kho2(3)=-730.

  CALL rxinit(r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe)
  nref=0  ;  ref(:)=' ' ;  CALL addref(progname,'MJ19KMV000',nref,ref,chem)
  np=0

! Following channels are considered for the reactions of RO2 with HO2, 
! as discussed in Jenkin et al., ACP, 2019:
! RCO(OO.)+HO2 -> RCO(OOH)+O2     (a) INDEX 1 
! RCO(OO.)+HO2 -> RCO(OH)+O3      (b) INDEX 2 
! RCO(OO.)+HO2 -> R(-H)O+H2O+O2   (c) INDEX 3  ! unavailable for RCO3
! RCO(OO.)+HO2 -> RCO(O.)+OH+O2   (d) INDEX 4
! RCO(OO.)+HO2 -> R(-H)O+OH+HO2   (e) INDEX 5  ! unavailable for RCO3

  patratio(:)=0.
  patratio(1)=0.37 ;  patratio(2)=0.13 ;  patratio(4)=0.50
  
! overwrite for aromatic R-CO(OO.)
  DO i=1,ngr
    IF ((tbond(i,ip) == 1).AND.(tgroup(i)(1:1) == 'c')) THEN
      patratio(1)=0.65 ;  patratio(2)=0.15 ;  patratio(4)=0.20
    ENDIF
  ENDDO
 
! --- channel (a) - index 1: => RCO(OOH)+O2
  pold='(OO.)'  ;  pnew='(OOH)'
  CALL swap(group(ip),pold,tgroup(ip),pnew)
  CALL rebond(tbond,tgroup,tempfo,tempring)
  CALL stdchm(tempfo)
  CALL add1tonp(progname,chem,np)   
  s(np)=patratio(1)  ;  brtio=brch*s(np)
  CALL bratio(tempfo,brtio,p(np),nref,ref)
  tgroup(:)=group(:)

! --- channel (b) - index 2: => RCO(OH) + O3
  pold='(OO.)'  ;  pnew='(OH)'
  CALL swap(group(ip),pold,tgroup(ip),pnew)
  CALL rebond(tbond,tgroup,tempfo,tempring)
  CALL stdchm(tempfo)
  CALL add1tonp(progname,chem,np)   
  s(np)=patratio(2)  ;   brtio=brch*s(np)
  CALL bratio(tempfo,brtio,p(np),nref,ref)
  CALL add1tonp(progname,chem,np)   
  s(np)=patratio(2)  ;  p(np)='O3 '
  tgroup(:)=group(:)

! --- channel (b) - index 2: => RCO(O.) + OH (+O2)
  pold='(OO.)'  ;  pnew='(O.)'
  CALL swap(group(ip),pold,tgroup(ip),pnew)
  CALL rebond(tbond,tgroup,tempfo,tempring)

! check new formula, standardize and add products and coproducts to p(:)   
  yield=patratio(4);   brtio=brch
  CALL chknadd2p(chem,tempfo,yield,brtio,np,s,p,nref,ref)

  CALL add1tonp(progname,chem,np)   
  s(np)=patratio(4)  ;  p(np)='HO '
  tgroup(:)=group(:)

! reactant
  r(1)=pname    
  r(2)='HO2 '
 
  arrh(:)=kho2(:)
  CALL rxwrit(mecu,r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe,ref)

END SUBROUTINE rco3ho2


!=======================================================================
! PURPOSE: perform and write the RCO3+RO2 permutation reactions.
! Two cases are considered: 
!   - take the 9 classes of RO2 into account (if rx_ro2_multiclass==.TRUE.)
!   - take only CH3O2 account (if rx_ro2_multiclass==.FALSE.)
!
! PROTOCOL FROM M. JENKIN ET AL., ACP, 2019
!=======================================================================
SUBROUTINE rco3ro2(chem,pname,bond,group,brch,ip)
  USE keyparameter, ONLY: mxlco,mxpd,mecu,mxcopd
  USE keyflag, ONLY: rx_ro2_multiclass      ! If true, consider all classes of RO2
  USE references, ONLY:mxlcod
  USE atomtool, ONLY: getatoms
  USE reactool, ONLY: swap,rebond
  USE dictstacktool, ONLY: bratio
  USE normchem, ONLY: stdchm  
  USE radchktool, ONLY: radchk
  USE rxwrttool, ONLY:rxwrit,rxinit
  USE toolbox, ONLY: stoperr,add1tonp,addrx,addref
  USE tempo, ONLY: chknadd2p

  CHARACTER(LEN=*),INTENT(IN) :: chem  ! formula of the VOC considered
  CHARACTER(LEN=*),INTENT(IN) :: pname ! short name of chem (parent compound) 
  CHARACTER(LEN=*),INTENT(IN) :: group(:) ! group (without ring joining char)
  INTEGER,INTENT(IN) :: bond(:,:)      ! bond matrix
  REAL,INTENT(IN)    :: brch           ! max yield of the input species
  INTEGER,INTENT(IN) :: ip             ! node bearig the peroxy (OO.) group

  CHARACTER(LEN=LEN(group(1))) :: tgroup(SIZE(group))
  CHARACTER(LEN=LEN(group(1))) :: pold,pnew
  INTEGER         :: tbond(SIZE(bond,1),SIZE(bond,2))
  CHARACTER(LEN=LEN(chem)) :: tempfo,p_ol,p_rad(2)

  CHARACTER(LEN=mxlco) :: p(mxpd), r(3)
  REAL    :: s(mxpd), arrh(3)
  INTEGER :: idreac, nlabel
  REAL    :: xlabel,folow(3),fotroe(4)

  INTEGER :: np,i,j,tempring
  REAL    :: brtio, wf

  INTEGER,PARAMETER    :: mxrpd=2   ! max # of products returned by radchk sub
  CHARACTER(LEN=mxlco) :: coprod(mxcopd),coprod2(mxcopd)
  CHARACTER(LEN=LEN(chem)) :: rdckpd(mxrpd) 
  CHARACTER(LEN=mxlco) :: rdckcopd(mxrpd,mxcopd)
  REAL    :: sc(mxrpd)
  INTEGER :: nip

  INTEGER,PARAMETER     :: mxref=10
  INTEGER               :: nref
  CHARACTER(LEN=mxlcod) :: ref(mxref)
  INTEGER               :: ntpcom
  CHARACTER(LEN=mxlcod) :: tpcom(mxref)
  CHARACTER(LEN=7)      :: progname='rco3ro2'

! set the reaction rate (Jenkin et al., ACPD, 2019) 
  REAL,DIMENSION(3),PARAMETER :: krco3=(/2.0E-12, 0., -508. /)

  CHARACTER(LEN=mxlco),DIMENSION(9), PARAMETER :: kwdclass= &
    (/'PERO1 ' , 'PERO2 ' , 'PERO3 ' , 'PERO4 ' , 'PERO5 ', &
      'MEPERO' , 'PERO7 ' , 'PERO8 ' , 'PERO9 '/)

  REAL,DIMENSION(9,3),PARAMETER :: stoi_a=RESHAPE( &  ! acyl RO2
      (/  1.0 , 0.0 , 0.0 ,  & ! PERO1- unsubstituted ter RO2
          0.8 , 0.0 , 0.2 ,  & ! PERO2- iC3H7O2
          1.0 , 0.0 , 0.0 ,  & ! PERO3- ter RO2 with alpha or beta O or N
          0.8 , 0.0 , 0.2 ,  & ! PERO4- C2H5O2 & unsubstituted sec RO2
          1.0 , 0.0 , 0.0 ,  & ! PERO5- subst. ter RO2 and allyl or benzyl
          0.8 , 0.0 , 0.2 ,  & ! MEPERO- CH3O2 
          0.8 , 0.0 , 0.2 ,  & ! PERO7- unsubst. prim RO2 & subst. sec RO2 
          0.8 , 0.0 , 0.2 ,  & ! PERO8- subst. prim RO2 & subst. sec RO2 with allyl 
          1.0 , 0.0 , 0.0 /),& ! PERO9- acyl RO2
       SHAPE(stoi_a), ORDER=(/2,1/) )  !- order (/2,1/) is to switch row & column from reshape

  ntpcom=0 ; tpcom(:)=' '
  tgroup(:)=group(:) ;  tbond(:,:)=bond(:,:)
  coprod(:)=' ' ;  coprod2(:)=' '

  wf=0.02
       
! find the various products linked to the species (alkoxy, acid)
  p_rad=' ' ;   p_ol=' '
  
  pold='(OO.)' ;  pnew='(O.)'                        !--- alkoxy product
  CALL swap(group(ip),pold,tgroup(ip),pnew)
  CALL rebond(tbond,tgroup,tempfo,tempring)
  CALL radchk(tempfo,rdckpd,rdckcopd,nip,sc,ntpcom,tpcom)
  p_rad(1)=rdckpd(1)   ;   CALL stdchm(p_rad(1))
  coprod(:)=rdckcopd(1,:) 
  IF (nip==2) THEN
    p_rad(2)=rdckpd(2) ;  CALL stdchm(p_rad(2))
    coprod2(:)=rdckcopd(2,:)
  ENDIF
  
  pold='(OO.)' ;  pnew='(OH)'                        !--- acid product
  CALL swap(group(ip),pold,tgroup(ip),pnew)
  CALL rebond(tbond,tgroup,p_ol,tempring)
  CALL stdchm(p_ol)
  
! multi-RO2-class parameterization
! --------------------------------
  IF (rx_ro2_multiclass) THEN 
    DO j=1,9
      CALL rxinit(r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe)
      nref=0  ;  ref(:)=' ' ;  CALL addref(progname,'MJ19KMV000',nref,ref,chem)
      IF (ntpcom/=0) THEN    ! add reference linked to R(O.) 
        !DO i=1,ntpcom  ;  nref=nref+1  ;  ref(nref)=tpcom(i)  ;  ENDDO
        DO i=1,ntpcom ; CALL addref(progname,tpcom(i),nref,ref,chem) ; ENDDO
      ENDIF
      np=0

! radical channel, 1st product
      CALL add1tonp(progname,chem,np)
      s(np)=stoi_a(j,1)*sc(1)  ;  brtio=brch*wf 
      CALL bratio(p_rad(1),brtio,p(np),nref,ref)
      DO i=1,SIZE(coprod)
        IF (coprod(i)(1:1)/=' ') THEN
          CALL add1tonp(progname,chem,np)
          s(np)=sc(1)*stoi_a(j,1) ; p(np)=coprod(i)
        ENDIF
      ENDDO

! radical channel, 2nd product (if any)
      IF (nip==2) THEN
        CALL add1tonp(progname,chem,np)
        s(np)=sc(2)*stoi_a(j,1)
        CALL bratio(p_rad(2),brch,p(np),nref,ref)
        DO i=1,SIZE(coprod2)
          IF (coprod2(i)(1:1)/=' ') THEN
            CALL add1tonp(progname,chem,np)
            s(np)=sc(2)*stoi_a(j,1) ;  p(np)=coprod2(i)
          ENDIF
        ENDDO
      ENDIF

! "H received from counter" channel (make carboxylic acid)
      IF (stoi_a(j,3)>0.) THEN
        CALL add1tonp(progname,chem,np)
        s(np)=stoi_a(j,3)  ;  brtio=brch*s(np)*wf
        CALL bratio(p_ol,brtio,p(np),nref,ref)
      ENDIF

! write out
      r(1)=pname  ;  r(2)=kwdclass(j)    ! reactant
      arrh(:)=krco3(:)
      CALL rxwrit(mecu,r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe,ref)
    ENDDO

! single RO2 class parameterization
! ---------------------------------
  ELSE 

    CALL rxinit(r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe)
    nref=0  ;  ref(:)=' ' ;  CALL addref(progname,'MJ19KMV000',nref,ref,chem)
    np=0
            
! radical channel, 1st product
    CALL add1tonp(progname,chem,np)
    s(np)=stoi_a(6,1)*sc(1) ;  brtio=brch*wf 
    CALL bratio(p_rad(1),brtio,p(np),nref,ref)
    DO i=1,SIZE(coprod)
      IF (coprod(i)(1:1)/=' ') THEN
        CALL add1tonp(progname,chem,np)
        s(np)=sc(1)*stoi_a(6,1) ; p(np)=coprod(i)
      ENDIF
    ENDDO

! radical channel, 2nd product (if any)
    IF (nip==2) THEN
      CALL add1tonp(progname,chem,np)
      s(np)=sc(2)*stoi_a(6,1)
      CALL bratio(p_rad(2),brch,p(np),nref,ref)
      DO i=1,SIZE(coprod2)
        IF (coprod2(i)(1:1)/=' ') THEN
          CALL add1tonp(progname,chem,np)
          s(np)=sc(2)*stoi_a(6,1) ;  p(np)=coprod2(i)
        ENDIF
      ENDDO
    ENDIF

! "H received" from counter channel
    IF (stoi_a(6,3)>0.) THEN
      CALL add1tonp(progname,chem,np)
      s(np)=stoi_a(6,3) ;  brtio=brch*s(np)*wf
      CALL bratio(p_ol,brtio,p(np),nref,ref)
    ENDIF

! write the reaction
    r(1)=pname ; r(2)='PERO1'       ! reactant
    arrh(:)=krco3(:)
    CALL rxwrit(mecu,r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe,ref)
  ENDIF

END SUBROUTINE rco3ro2

END MODULE rco3tool
