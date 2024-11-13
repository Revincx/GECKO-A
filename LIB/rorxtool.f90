MODULE rorxtool
IMPLICIT NONE
CONTAINS
!SUBROUTINE rooxy(bond,group,nca,ia,nring,flo2,po2,copo2,arrho2)
!SUBROUTINE rono2(bond,group,nca,ia,nring,flno2,pno2,arrno2)
!SUBROUTINE roester(bond,group,nca,ia,nring,flest,pdest,copdest,arrhest)
!SUBROUTINE roaro(idnam,bond,group,nca,ia)
!SUBROUTINE rodecnitro(bond,group,ia,ndec,fldec,pdec,copdec, arrhdec)

! ======================================================================
! PURPOSE: perform the R(O.)+O2 reaction. Reaction rates are based on 
! the SAPRC99 chemical scheme, updated from Atkinson 2007. 
! Rate constants are computed using kO2=2.5E-14 exp(-300/T) for both
! primary and secondary alkoxy radicals.
! ======================================================================
SUBROUTINE rooxy(chem,bond,group,nca,ia,nring,flo2,po2,scpo2,copdo2, &
                 ycopdo2,arrho2,nrfo2,rfo2)
  USE keyparameter, ONLY: mxcp
  USE ringtool, ONLY: findring 
  USE reactool, ONLY: swap, rebond
  USE normchem, ONLY: stdchm
  USE radchktool, ONLY: radchk,lump_copd  
  USE mapping, ONLY: abcde_map    
  USE toolbox, ONLY: stoperr,addref
  IMPLICIT NONE

  CHARACTER(LEN=*),INTENT(IN):: chem
  INTEGER,INTENT(IN)  :: bond(:,:)
  CHARACTER(LEN=*),INTENT(IN) :: group(:)
  INTEGER,INTENT(IN)  :: nca     ! # of nodes
  INTEGER,INTENT(IN)  :: ia      ! index for node bearing (O.) group
  INTEGER,INTENT(IN)  :: nring   
  INTEGER,INTENT(OUT) :: flo2    ! O2 reaction flag  
  CHARACTER(LEN=*),INTENT(OUT):: po2(:)     ! main product of O2 reaction (size 2)
  REAL,INTENT(OUT)      :: scpo2(:)         ! stoi. coef of O2 products (size 2)
  CHARACTER(LEN=*),INTENT(OUT):: copdo2(:)  ! list of coproducts (if decomposition)
  REAl,INTENT(OUT)      :: ycopdo2(:)       ! yield of the stoi; coef. in the copdo2 list
  REAL,INTENT(OUT)      :: arrho2(:)
  INTEGER,INTENT(INOUT) :: nrfo2            ! # of references added in the reference list
  CHARACTER(LEN=*),INTENT(INOUT):: rfo2(:)  ! reference for O2 reaction

  INTEGER :: tbond(SIZE(bond,1),SIZE(bond,2))
  CHARACTER(LEN=LEN(group(1))):: tgroup(SIZE(group))
  CHARACTER(LEN=LEN(group(1))):: pold, pnew,tempkg
  CHARACTER(LEN=LEN(po2))     :: tpchem
  INTEGER :: ring(SIZE(bond,1))
  INTEGER :: i,j,k,rngflg,tpnring,nip,nc,tnode,ncopd

  INTEGER :: nabcde(5)                   ! deepest path is 5 here
  INTEGER :: tabcde(5,mxcp,SIZE(bond,1)) ! deepest path is 5 here

  INTEGER,PARAMETER :: mxrpd=2  ! max # of product returned by radchk sub
  CHARACTER(LEN=LEN(copdo2(1))) :: rdckcopd(mxrpd,SIZE(copdo2))
  REAL :: sc(mxrpd)

  CHARACTER(LEN=6),PARAMETER :: progname='rooxy '
  CHARACTER(LEN=70)          :: mesg

  flo2=0  ; po2(:)=' ' ; scpo2(:)=0. ; 
  copdo2(:)=' ' ; ycopdo2(:)=0. ; arrho2(:)=0.
  tbond(:,:)=bond(:,:) ; tgroup(:)=group(:)  
  
! primary RO: R-CH2(O.)
  IF (group(ia)(1:3)=='CH2') THEN
    flo2=1  ;  tgroup(ia)='CHO'

! secondary RO: R-CH(O.)-R
  ELSE IF (group(ia)(1:2)=='CH') THEN
    flo2=1
    IF (nring>0) THEN  ! disallow reaction for rings < 5 members (only decomposition)
      DO i=1,nca
        IF (tbond(i,ia)/=0) THEN
          CALL findring(i,ia,nca,tbond,rngflg,ring)
          IF (rngflg/=0) EXIT
        ENDIF
      ENDDO
      IF (rngflg/=0) THEN
        CALL abcde_map(tbond,ia,nca,nabcde,tabcde)
        rgloop: DO j=3,5
          DO k=1,nabcde(j)
            tnode=tabcde(j,k,j)
            IF (tbond(ia,tnode)/=0) THEN
              flo2=0
              EXIT rgloop
            ENDIF
          ENDDO
        ENDDO rgloop
      ENDIF
    ENDIF
    pold='CH'   ; pnew='CO' ; CALL swap(group(ia),pold,tempkg,pnew)
    pold='(O.)' ; pnew=' '  ; CALL swap(tempkg,pold,tgroup(ia),pnew)
  ENDIF

! rebuild and rename:
  IF (flo2==1) THEN
    CALL rebond(tbond,tgroup,po2(1),tpnring)
    CALL stdchm(po2(1)) ; scpo2(1)=1.
    arrho2(1) =2.5E-14  ;  arrho2(3)=300.
    CALL addref(progname,'RA07KER000',nrfo2,rfo2,chem)
  ENDIF
  
! If po2 is R-C(=O)ONO2 , assume decomposition =>  R. + CO2 + NO2 
! nb: doesn't involve breaking a ring bond
  IF (INDEX(po2(1),'CO(ONO2)')/=0) THEN
    grloop : DO i=1,nca
      IF (INDEX(tgroup(i),'CO(ONO2)')/=0) THEN
        tgroup(i)=' '
        DO j=1,nca
          IF (bond(i,j)==1) THEN
            tbond(i,j)=0 ;  tbond(j,i)=0
            nc=INDEX(tgroup(j),' ') ; tgroup(j)(nc:nc)='.'
            CALL rebond(tbond,tgroup,tpchem,tpnring)
            CALL radchk(tpchem,po2,rdckcopd,nip,sc,nrfo2,rfo2)
            scpo2(:)=sc(:) ; ncopd=0
            CALL lump_copd(chem,rdckcopd,sc,ncopd,copdo2,ycopdo2)
            ncopd=ncopd+1 ; copdo2(ncopd)='CO2' ; ycopdo2(ncopd)=1.
            ncopd=ncopd+1 ; copdo2(ncopd)='NO2' ; ycopdo2(ncopd)=1.
            IF (ncopd > SIZE(copdo2)) THEN
              mesg="Too many copdct for a 'CO(ONO2)' species "
              CALL stoperr(progname,mesg,chem)
            ENDIF
          ENDIF
        ENDDO  
        nrfo2=nrfo2+1 ; rfo2(nrfo2)='RO_CO(NO3)'
        CALL addref(progname,'RO_CO(NO3)',nrfo2,rfo2,chem)
        EXIT grloop ! Exit after one decomposition to avoid diradical
      ENDIF
    ENDDO grloop
  ENDIF    
  
END SUBROUTINE rooxy

! ======================================================================
! PURPOSE: perform the R(O.)+NO2 reaction. Reaction rates are based on
! an average of IUPAC rates (Atkinson et al 2006)
! Rate constants average 3.3e-11 molec-1 s-1.
! ======================================================================
SUBROUTINE rono2(chem,bond,group,nca,ia,nring,flno2,pno2,scpno2,copdno2, &
                 ycopdno2,arrhno2,nrfno2,rfno2)
  USE keyparameter, ONLY: mxcp
  USE ringtool, ONLY: findring
  USE reactool, ONLY: swap, rebond
  USE normchem, ONLY: stdchm
  USE radchktool, ONLY: radchk,lump_copd
  USE mapping, ONLY: abcde_map
  USE toolbox, ONLY: stoperr,addref
  IMPLICIT NONE

  CHARACTER(LEN=*),INTENT(IN):: chem
  INTEGER,INTENT(IN)  :: bond(:,:)
  CHARACTER(LEN=*),INTENT(IN) :: group(:)
  INTEGER,INTENT(IN)  :: nca     ! # of nodes
  INTEGER,INTENT(IN)  :: ia      ! index for node bearing (O.) group
  INTEGER,INTENT(IN)  :: nring
  INTEGER,INTENT(OUT) :: flno2   ! NO2 reaction flag
  CHARACTER(LEN=*),INTENT(OUT):: pno2(:)   ! main product of NO2 reaction (size 2)
  REAL,INTENT(OUT)      :: scpno2(:)         ! stoi. coef of NO2 products (size 2)
  CHARACTER(LEN=*),INTENT(OUT):: copdno2(:)  ! list of coproducts (if decomposition)
  REAl,INTENT(OUT)      :: ycopdno2(:)       ! yields of the coproducts 
  REAL,INTENT(OUT)      :: arrhno2(:)
  INTEGER,INTENT(INOUT) :: nrfno2           ! # of references added to the reference list
  CHARACTER(LEN=*),INTENT(INOUT):: rfno2(:) ! reference for O2 reaction

  INTEGER :: tbond(SIZE(bond,1),SIZE(bond,2))
  CHARACTER(LEN=LEN(group(1))):: tgroup(SIZE(group))
  CHARACTER(LEN=LEN(group(1))):: pold, pnew,tempkg
  CHARACTER(LEN=LEN(pno2))     :: tpchem
  INTEGER :: ring(SIZE(bond,1))
  INTEGER :: i,j,k,rngflg,tpnring,nip,nc,tnode,ncopd

  INTEGER :: nabcde(5)                   ! deepest path is 5 here
  INTEGER :: tabcde(5,mxcp,SIZE(bond,1)) ! deepest path is 5 here

  INTEGER,PARAMETER :: mxrpd=2  ! max # of product returned by radchk sub
  CHARACTER(LEN=LEN(copdno2(1))) :: rdckcopd(mxrpd,SIZE(copdno2))
  REAL :: sc(mxrpd)

  CHARACTER(LEN=6),PARAMETER :: progname='rono2 '
  CHARACTER(LEN=70)          :: mesg

  flno2=0  ; pno2(:)=' ' ; arrhno2(:)=0.
  copdno2(:)=' ' ; ycopdno2(:)=0.  ; scpno2(:)=0. ;
  tbond(:,:)=bond(:,:) ; tgroup(:)=group(:)

! primary RO: R-CH2(O.)
  IF (group(ia)(1:3)=='CH2') THEN
    flno2=1  ;  tgroup(ia)='CH2(ONO2)'

! secondary RO: R-CH(O.)-R
  ELSE IF (group(ia)(1:2)=='CH') THEN
    flno2=1
    IF (nring>0) THEN  ! disallow reaction for rings < 5 members (only decomposition)
      DO i=1,nca
        IF (tbond(i,ia)/=0) THEN
          CALL findring(i,ia,nca,tbond,rngflg,ring)
          IF (rngflg/=0) EXIT
        ENDIF
      ENDDO
      IF (rngflg/=0) THEN
        CALL abcde_map(tbond,ia,nca,nabcde,tabcde)
        rgloop: DO j=3,5
          DO k=1,nabcde(j)
            tnode=tabcde(j,k,j)
            IF (tbond(ia,tnode)/=0) THEN
              flno2=0
              EXIT rgloop
            ENDIF
          ENDDO
        ENDDO rgloop
      ENDIF
    ENDIF
    pold='CH'   ; pnew='CH' ; CALL swap(group(ia),pold,tempkg,pnew)
    pold='(O.)' ; pnew='(ONO2)'  ; CALL swap(tempkg,pold,tgroup(ia),pnew)
  ENDIF

! rebuild and rename:
  IF (flno2==1) THEN
    CALL rebond(tbond,tgroup,pno2(1),tpnring)
    CALL stdchm(pno2(1))
    arrhno2(1) = 3.3E-11
! average of IUPAC 2006 values in Atkinson et al, 2006
    CALL addref(progname,'RA06KMR000',nrfno2,rfno2,chem)
  ENDIF

! If pno2 is R-C(=O)ONO2 , assume decomposition =>  R. + CO2 + NO2
! nb: doesn't involve breaking a ring bond

  IF (INDEX(pno2(1),'CO(ONO2)')/=0) THEN
    grloop : DO i=1,nca
      IF (INDEX(tgroup(i),'CO(ONO2)')/=0) THEN
        tgroup(i)=' '
        DO j=1,nca
          IF (bond(i,j)==1) THEN
            tbond(i,j)=0 ;  tbond(j,i)=0
            nc=INDEX(tgroup(j),' ') ; tgroup(j)(nc:nc)='.'
            CALL rebond(tbond,tgroup,tpchem,tpnring)
            CALL radchk(tpchem,pno2,rdckcopd,nip,sc,nrfno2,rfno2)
            scpno2(:)=sc(:) ; ncopd=nip
            CALL lump_copd(chem,rdckcopd,sc,ncopd,copdno2,ycopdno2)
            ncopd=ncopd+1 ; copdno2(ncopd)='CO2' ; ycopdno2(ncopd)=1.
            ncopd=ncopd+1 ; copdno2(ncopd)='NO2' ; ycopdno2(ncopd)=1.
            IF (ncopd > SIZE(copdno2)) THEN
              mesg="Too many copdct for a 'CO(ONO2)' species "
              CALL stoperr(progname,mesg,chem)
            ENDIF
          ENDIF
        ENDDO
        nrfno2=nrfno2+1 ; rfno2(nrfno2)='RO_CO(NO3)'
        CALL addref(progname,'RO_CO(NO3)',nrfno2,rfno2,chem)
        EXIT grloop ! Exit after one decomposition to avoid diradical
      ENDIF
    ENDDO grloop
  ENDIF

END SUBROUTINE rono2

! ======================================================================
! PURPOSE: manage the alpha-ester rearrangement, ie. the transformation:
!    R-CO-O-CH(O.)R' => R-CO(OH) + .CO-R'
! Rate constants are based on the SAPRC99 chemical scheme, using an
! arrhenius expression, where A=8.0E+10 and Ea/R=3723, based on the
! qualitative data of Tuazon (98). Only 1-4 shift is considered ; 1-5 
! assumed to be not significant. No reaction considered for formate (no 
! data available)
! ======================================================================
SUBROUTINE roester(bond,group,nca,ia,nring,flest,pdest,copdest,arrhest,&
                   loregro,nrfest,rfest)
  USE ringtool, ONLY: findring 
  USE reactool, ONLY: swap, rebond
  USE normchem, ONLY: stdchm
  USE radchktool, ONLY: radchk    
  USE fragmenttool, ONLY: fragm
  USE toolbox, ONLY: addref
  IMPLICIT NONE
  
  INTEGER,INTENT(IN)  :: bond(:,:)
  CHARACTER(LEN=*),INTENT(IN) :: group(:)
  INTEGER,INTENT(IN)  :: nca     ! # of nodes
  INTEGER,INTENT(IN)  :: ia      ! index for node bearing (O.) group
  INTEGER,INTENT(IN)  :: nring   
  INTEGER,INTENT(OUT) :: flest   ! ester rearrangement flag  
  CHARACTER(LEN=*),INTENT(OUT)  :: pdest(:)   ! major product (size 2)
  CHARACTER(LEN=*),INTENT(OUT)  :: copdest(:) ! coproducts (in short names)
  REAL,INTENT(OUT)      :: arrhest(:)
  LOGICAL,INTENT(INOUT) :: loregro            ! .TRUE. if regular "RO" (all reactions considered)
  INTEGER,INTENT(INOUT) :: nrfest             ! # of references added in the reference list
  CHARACTER(LEN=*),INTENT(INOUT):: rfest(:)   ! ref/comment, ester channel

  CHARACTER(LEN=LEN(group(1))) :: tgroup(SIZE(group))
  CHARACTER(LEN=LEN(group(1))) :: pold, pnew,tempgr
  INTEGER :: tbond(SIZE(bond,1),SIZE(bond,1))
  INTEGER :: ib,j,k,tpnring,nip,rngflg
  INTEGER :: ring(SIZE(group))

  CHARACTER(LEN=LEN(pdest(1)))   :: tpchem
  CHARACTER(LEN=LEN(pdest(1)))   :: rdckpd(SIZE(pdest))
  CHARACTER(LEN=LEN(copdest(1))) :: rdckcopd(SIZE(pdest),SIZE(copdest))
  REAL                           :: sc(SIZE(pdest))

  CHARACTER(LEN=8),PARAMETER :: progname='roester '
  
  pdest(:)=' ' ; copdest(:)=' ' ; flest=0 ; ring(:)=0 ; arrhest(:)=0.
  tbond(:,:)=bond(:,:)  ;  tgroup(:)=group(:)  
    
  DO ib=1,nca
    IF ((tbond(ib,ia)==3).AND.(tgroup(ia)(1:2)=='CH')) THEN
      IF (nring>0) CALL findring(ib,ia,nca,tbond,rngflg,ring)
      DO j=1,nca
        IF ((tbond(j,ib)/=0).AND.(tgroup(j)=='CO')) THEN
          flest = 1  
          pold = 'CO' ;  pnew = 'CO(OH)'
          CALL swap (group(j),pold,tgroup(j),pnew)
          pold = '-O-'  ;  pnew = ' '
          CALL swap (group(ib),pold,tgroup(ib),pnew)
          IF (group(ia)(1:3)=='CH2') THEN 
            pold='CH2(O.)' ; pnew='CHO.' 
            CALL swap(group(ia),pold,tgroup(ia),pnew)
          ELSE IF (group(ia)(1:2)=='CH') THEN 
            pold='CH'  ; pnew='CO'
            CALL swap(group(ia),pold,tgroup(ia),pnew)
            pold='(O.)'  ; pnew='.' ; tempgr=tgroup(ia)
            CALL swap(tempgr,pold,tgroup(ia),pnew)
          ENDIF
          tbond(ib,j)=0 ; tbond(j,ib)=0 ; tbond(ia,ib)=0 ; tbond(ib,ia)=0
      
          IF (ring(ia)/=0) THEN  ! ring case (1 product)
            CALL rebond(tbond,tgroup,tpchem,tpnring)
            CALL radchk(tpchem,rdckpd,rdckcopd,nip,sc,nrfest,rfest)
            pdest(1) = rdckpd(1)
            IF (nip/=1) STOP "A: in roester, unexpected 2 products"
            copdest(:) = rdckcopd(1,:)
            CALL stdchm(pdest(1))
            pdest(2)=' '
          ELSE                   ! no ring (2 products)
            CALL fragm(tbond,tgroup,pdest(1),pdest(2))            
            DO k=1,2
              IF (INDEX(pdest(k),'.')/=0) THEN
                tpchem = pdest(k)
                CALL radchk(tpchem,rdckpd,rdckcopd,nip,sc,nrfest,rfest)
                pdest(k) = rdckpd(1)
                IF (nip/=1) STOP "B: in roester, unexpected 2 products"
                copdest(:) = rdckcopd(1,:)
              ENDIF
              CALL stdchm (pdest(k))
            ENDDO
          ENDIF
          arrhest(1)=8.E+10  ;  arrhest(2)=0.  ;  arrhest(3)=3.723E+03
          CALL addref(progname,'RA07KER000',nrfest,rfest)
          tgroup(:)=group(:)  ;  tbond(:,:)=bond(:,:) 
        ENDIF
      ENDDO
    ENDIF 
  ENDDO

! If >C=CR-CO-O-C(O.)- : assumed to be rearranged by the alpha-ester only
! raise flag
  IF (flest/=0) THEN
    DO ib=1,nca
      IF ((bond(ib,ia)==3).AND.(group(ia)/='C(O.)')) THEN
        DO  j=1,nca
          IF ((bond(ib,j)==3).AND.(group(j)(1:2)=='CO')) THEN
            DO k=1,nca
              IF ((bond(k,j)==1).AND.(group(k)(1:2)=='Cd')) THEN
                loregro = .FALSE.
                nrfest=nrfest+1 ; rfest(nrfest)='RO_ESTER_U'
                CALL addref(progname,'RO_ESTER_U',nrfest,rfest)
              ENDIF
            ENDDO        
          ENDIF
        ENDDO
      ENDIF
    ENDDO
  ENDIF          
END SUBROUTINE roester

!=======================================================================
! PURPOSE: manage phenoxy [phi-(O.)] radical chemistry. 
! Three pathways considered: (1) reaction with O3 to form a peroxy 
! radical, (2) reaction with NO2 to form a nitro phenol (if H available 
! close next to O. and (3) reaction with HO2 to form a phenolic 
! product + O2
! Note: the reactions are directly written to the output file!
!=======================================================================
SUBROUTINE roaro(idnam,bond,group,nca,ia)
  USE references, ONLY:mxlcod
  USE keyparameter, ONLY: mxlco,mxlfo,mxpd,mecu
  USE reactool, ONLY: rebond
  USE normchem, ONLY: stdchm
  USE radchktool, ONLY: radchk    
  USE rxwrttool, ONLY:rxwrit,rxinit
  USE dictstacktool, ONLY: bratio
  USE toolbox, ONLY: addref
  IMPLICIT NONE

  CHARACTER(LEN=*) :: idnam       !  
  INTEGER,INTENT(IN)  :: bond(:,:)
  CHARACTER(LEN=*),INTENT(IN) :: group(:)
  INTEGER,INTENT(IN)  :: nca     ! # of nodes
  INTEGER,INTENT(IN)  :: ia      ! index for node bearing (O.) group

  INTEGER :: tbond(SIZE(bond,1),SIZE(bond,1))
  CHARACTER(LEN=LEN(group(1))) :: tgroup(SIZE(group))
  INTEGER :: i,j,tpnring,nip
  CHARACTER(LEN=mxlfo) :: tpchem, tchem(2)

! variable required for radchk
  INTEGER,PARAMETER :: mxrpd=2  ! max # of products expected to be returned by radchk sub
  CHARACTER(LEN=mxlfo) :: rdckpd(mxrpd)
  CHARACTER(LEN=mxlco) :: rdckcopd(mxrpd,SIZE(group)) ! no copdct involved
  REAL    :: sc(mxrpd)

! variable required to write the reactions
  CHARACTER(LEN=mxlco) :: p(mxpd), r(3)
  REAL    :: s(mxpd), arrh(3)
  INTEGER :: idreac, nlabel
  REAL    :: xlabel,folow(3),fotroe(4)
  REAL    :: brtio

  INTEGER,PARAMETER     :: mxref=10     !--- max # of ref/comment per reaction
  INTEGER,PARAMETER     :: mxrx=5       !--- max # of reaction per phenoxy
  INTEGER               :: nref(mxrx),nrx
  CHARACTER(LEN=mxlcod) :: com(mxrx,mxref)
  CHARACTER(LEN=6),PARAMETER :: progname='roaro '

  tgroup(:)=group(:)  ;  tbond(:,:)=bond(:,:)
  nref(:)=0 ; nrx=0 ; com(:,:)=' '

  IF (tgroup(ia)/='c(O.)') THEN 
    PRINT*, 'error in roaro - phenoxy + O3 '
    STOP "in roaro"
  ENDIF
        
! 1st pathway : ozone reaction 
  nrx=nrx+1 
  tgroup(ia)='c(OO.)' ;  CALL rebond(tbond,tgroup,tpchem,tpnring)
  CALL radchk(tpchem,rdckpd,rdckcopd,nip,sc,nref(nrx),com(nrx,:))
  tpchem = rdckpd(1)
  IF (nip/=1) THEN
    PRINT*, "1 products only expected for phenoxy in roaro"
    STOP "in roaro"
  ENDIF
  CALL stdchm(tpchem)
  tgroup(:)=group(:) ; tbond(:,:)=bond(:,:)

  CALL rxinit(r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe)
  r(1) = idnam  ;  r(2) = 'O3'  ;  s(1) = 1. 
  arrh(1) = 2.86E-13
  CALL addref(progname,'ZT99KMR000',nref(nrx),com(nrx,:))
  brtio=1.  ;  CALL bratio(tpchem,brtio,p(1),nref(nrx),com(nrx,:))
  CALL rxwrit(mecu,r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe,com(nrx,:))

! 2nd pathaway : nitro formation
  tchem(:)=' '  ;  j=0
  DO i=1,nca
    IF ((tbond(i,ia)==1).AND.(tgroup(i)(1:3)=='cH ')) THEN
      j=j+1 ; tgroup(ia)='c(OH)' ;  tgroup(i)='c(NO2)'
      CALL rebond(tbond,tgroup,tpchem,tpnring)
      CALL stdchm(tpchem)
      tchem(j)=tpchem
      tgroup(:)=group(:)  ; tbond(:,:)=bond(:,:)
    ENDIF
  ENDDO
  IF (j>2) THEN
    PRINT*, "in roaro. Too many alpha cH for a phenoxy"
    STOP "in roaro"
  ENDIF
  IF (tchem(1)==tchem(2)) tchem(2)=' '

  DO i=1,2
    IF (tchem(i)/=' ') THEN
      CALL rxinit(r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe)
      r(1) = idnam ; r(2) = 'NO2' ; s(1) = 1.
      arrh(1) = 2.08E-12
      nrx=nrx+1;  CALL addref(progname,'JP98KMR000',nref(nrx),com(nrx,:))
      brtio=1. ;  CALL bratio(tchem(i),brtio,p(1),nref(nrx),com(nrx,:))
      CALL rxwrit(mecu,r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe,com(nrx,:))
    ENDIF
  ENDDO

! 3rd reaction with HO2 : c(O.) + HO2 -> c(OH)
  tgroup(ia)='c(OH)'  ;  CALL rebond(tbond,tgroup,tpchem,tpnring)
  CALL stdchm(tpchem)
  tgroup(:)=group(:)  ;  tbond(:,:)=bond(:,:)

  CALL rxinit(r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe)
  r(1) = idnam  ;  s(1) = 1.
  r(2)='HO2 '
  arrh(1) = 2.3E-13 
  nrx=nrx+1 ;  CALL addref(progname,'SM11KER000',nref(nrx),com(nrx,:))
  brtio=1.  ;  CALL bratio(tpchem,brtio,p(1),nref(nrx),com(nrx,:))
  CALL rxwrit(mecu,r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe,com(nrx,:))

END SUBROUTINE roaro

!=======================================================================
! PURPOSE: perform the decomposition of nitro alkoxy radicals, i.e. 
!    R(NO2)(O.) -> R=O + NO2 only
! No competition expected (other pathways should be set to 0).
!=======================================================================
SUBROUTINE rodecnitro(bond,group,ia,ndec,fldec,pdec,copdec,arrhdec,nrfdec,rfdec)
  USE reactool, ONLY: swap, rebond
  USE normchem, ONLY: stdchm
  USE toolbox, ONLY: addref
  USE dectool, ONLY: coono2_dec
  IMPLICIT NONE

  INTEGER,INTENT(IN)  :: bond(:,:)
  CHARACTER(LEN=*),INTENT(IN) :: group(:)
  INTEGER,INTENT(IN)  :: ia      ! index for node bearing (O.) group
  INTEGER,INTENT(OUT) :: ndec          ! # of decomposition reactions
  INTEGER,INTENT(OUT) :: fldec(:)    ! decomposition reaction flag  
  CHARACTER(LEN=*),INTENT(OUT) :: pdec(:,:)
  CHARACTER(LEN=*),INTENT(OUT) :: copdec(:,:)
  REAL,INTENT(OUT)    :: arrhdec(:,:)
  INTEGER,INTENT(INOUT) :: nrfdec(:)          ! # of reference added for each reactions
  CHARACTER(LEN=*),INTENT(INOUT):: rfdec(:,:) ! ref/comments for each reactions

  INTEGER :: tbond(SIZE(bond,1),SIZE(bond,2))
  CHARACTER(LEN=LEN(group(1))) :: tgroup(SIZE(group))
  CHARACTER(LEN=LEN(group(1))) :: pold, pnew,tempkg
  INTEGER :: tpnring,nc,ngr

  CHARACTER(LEN=11),PARAMETER :: progname='rodecnitro '

  ndec=0  ; arrhdec(:,:)=0.  
  fldec(:)=0   ;  pdec(:,:)=' ' ; copdec(:,:)=' '
  tbond(:,:)=bond(:,:) ; tgroup(:)=group(:)  
  ngr=COUNT(group/=' ')
  
  ndec=1  ;  fldec(ndec)=1

  
  pold='(O.)' ;  pnew=' '
  CALL swap(tgroup(ia),pold,tempkg,pnew)  ;  tgroup(ia)=tempkg
  IF (INDEX(tgroup(ia),'CH')/=0) THEN ; pold='CH' ;  pnew='CHO'
  ELSE                                ; pold='C'  ;  pnew='CO'
  ENDIF

  CALL swap(tgroup(ia),pold,tempkg,pnew)
  nc=INDEX(tempkg,'(NO2)') ; tempkg(nc:)=tempkg(nc+5:)   ! rm NO2
  tgroup(ia)=tempkg
  CALL rebond(tbond,tgroup,pdec(ndec,1),tpnring)
  CALL stdchm(pdec(ndec,1))

  IF (INDEX(pdec(ndec,1),'CO(ONO2)')/=0) THEN !For RCO(ONO2), assume decomposition to R. + CO2 + NO2
    CALL coono2_dec(pdec(ndec,1),copdec(ndec,:),ngr,nrfdec(ndec),rfdec(ndec,:))
  ENDIF

  pdec(ndec,2)='NO2'
   
  arrhdec(ndec,1)=1E12 ; arrhdec(ndec,2)=0. ;  arrhdec(ndec,3)=0.
  CALL addref(progname,'RO_NITRO',nrfdec(ndec),rfdec(ndec,:))

END SUBROUTINE rodecnitro

END MODULE rorxtool
