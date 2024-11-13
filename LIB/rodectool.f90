MODULE rodectool
IMPLICIT NONE
CONTAINS
!SUBROUTINE rodec(chem,bond,group,nca,ia,nring, &
!           ndec,fldec,pdec,copdec,fldec2,pdec2,copdec2,sdec,arrhdec)
!SUBROUTINE rospedec(bond,group,nca,ia,ib,nring, &
!           ndec,fldec,pdec,copdec,fldec2,pdec2,copdec2,sdec, arrhdec)
!SUBROUTINE krodec_ver(tbond,tgroup,nca,ia,ib,nring,arrh)
!SUBROUTINE krodec_atk(chem,pdct1,pdct2,lgrp,copd1,oflg,arrh)

!=======================================================================
! PURPOSE: operate alkoxy decomposition. All C-C(O.) bonds are examined. 
! Note for particular cases: 
!  1. If Cd next to C(O.) (i.e. C=C-C(O.)-), then ignore dissociation
!  2. If phenyl (i.e. >c-C(O.)-), then ignore dissociation
! NOTE: Reation rate constant can be estimated using:
!   - Atkinson, 2007 (Atmospheric Environment)
!   - Vereecken and Peeters (2009) 
!=======================================================================
SUBROUTINE rodec(chem,bond,group,nca,ia,nring,ndec,fldec,pdec,copdec, &
                 fldec2,pdec2,copdec2,sdec,arrhdec,nrfdec,rfdec)
  USE keyparameter, ONLY: mxcp           
  USE cdtool, ONLY: alkcheck
  USE keyflag, ONLY: kdiss_sar          ! flag for the estimation method   
  USE reactool, ONLY: swap, rebond
  USE normchem, ONLY: stdchm
  USE mapping, ONLY: gettrack  
  USE radchktool, ONLY: radchk, single
  USE ringtool, ONLY: findring 
  USE stdgrbond, ONLY: grbond  
  USE fragmenttool, ONLY: fragm
  USE toolbox, ONLY: setbond,stoperr,addref
  USE dectool, ONLY: coono2_dec
  IMPLICIT NONE

  CHARACTER(LEN=*),INTENT(IN) :: chem
  INTEGER,INTENT(IN)  :: bond(:,:)
  CHARACTER(LEN=*),INTENT(IN) :: group(:)
  INTEGER,INTENT(IN)  :: nca     ! # of nodes
  INTEGER,INTENT(IN)  :: ia      ! index for node bearing (O.) group
  INTEGER,INTENT(IN)  :: nring   
  INTEGER,INTENT(OUT) :: ndec        ! # of decomp reactions
  INTEGER,INTENT(OUT) :: fldec(:)    ! decomp reaction flag  
  CHARACTER(LEN=*),INTENT(OUT) :: pdec(:,:)
  CHARACTER(LEN=*),INTENT(OUT) :: copdec(:,:)
  INTEGER,INTENT(OUT) :: fldec2(:)    ! decomp reaction flag  
  CHARACTER(LEN=*),INTENT(OUT) :: pdec2(:)
  CHARACTER(LEN=*),INTENT(OUT) :: copdec2(:,:)
  REAL,INTENT(OUT) :: sdec(:,:)
  REAL,INTENT(OUT) :: arrhdec(:,:)
  INTEGER,INTENT(INOUT) :: nrfdec(:)          ! # of reference added for each reactions
  CHARACTER(LEN=*),INTENT(INOUT):: rfdec(:,:) ! ref/comments for each reactions

  INTEGER :: tbond(SIZE(bond,1),SIZE(bond,2))
  CHARACTER(LEN=LEN(group(1))) :: tgroup(SIZE(group))
  CHARACTER(LEN=LEN(group(1))) :: pold, pnew,tempkg
  CHARACTER(LEN=LEN(pdec(1,1))):: tpchem
  INTEGER :: ring(SIZE(group))  
  INTEGER :: fldecether
  INTEGER :: i,ib,j,k,l,rngflg,tpnring,nip,nc
  
  INTEGER,PARAMETER :: mxrpd=2  ! max # of products returned by radchk sub
  CHARACTER(LEN=LEN(pdec(1,1)))  :: rdckpd(mxrpd)
  CHARACTER(LEN=LEN(copdec(1,1))):: rdckcopd(mxrpd,SIZE(copdec,2))
  CHARACTER(LEN=LEN(copdec(1,1))):: acopd
  REAL :: sc(mxrpd)

  CHARACTER(LEN=LEN(rfdec(1,1))) :: acom   ! a single code for comment 
  CHARACTER(LEN=12) :: resu

  CHARACTER(LEN=6),PARAMETER :: progname='rodec '
  CHARACTER(LEN=70)          :: mesg

  ndec=0  ;  sdec(:,:)=0.  ; arrhdec(:,:)=0. 
  fldec(:)=0   ;  pdec(:,:)=' ' ; copdec(:,:)=' '
  fldec2(:)=0  ;  pdec2(:)=' '  ; copdec2(:,:)=' '

  tbond(:,:)=bond(:,:) ; tgroup(:)=group(:)  
  fldecether=0  ;   ring(:)=0  

! scroll for beta group (ib) 
  grloop: DO ib=1,nca                       ! long loop
    IF (bond(ib,ia)==0) CYCLE grloop
    IF (group(ib)(1:2)=='Cd') CYCLE grloop
    IF (group(ib)(1:1)=='c') CYCLE grloop   ! ignore phenyl formation (magnify 2016)

    ndec=ndec+1
    IF (ndec>SIZE(fldec)) THEN
      mesg="ndec > SIZE(fldec) "
      CALL stoperr(progname,mesg,chem)
    ENDIF
    fldec(ndec)=1

! check if ring opening required
    rngflg=0                ! overwrite next if the bond belong to ring
    IF (nring>0) CALL findring(ib,ia,nca,tbond,rngflg,ring)

! ------------------------------   
! LINEAR MOLECULE DECOMPOSITION   
! ------------------------------   
    IF (rngflg==0) THEN

! change (O.) to carbonyl:
      IF (group(ia)(1:3)=='CH2') THEN ; tgroup(ia)='CH2O'
      ELSE IF (group(ia)(1:2)=='CH') THEN
        pold='(O.)' ;  pnew=' ' ;  CALL swap(group(ia),pold,tempkg,pnew)
        pold='CH' ; pnew='CHO'  ;  CALL swap(tempkg,pold,tgroup(ia),pnew)
      ELSE IF (group(ia)(1:1)=='C' ) THEN
        pold='(O.)' ;  pnew=' ' ;  CALL swap(group(ia),pold,tempkg,pnew)
        IF (group(ia)(1:2)/='Cd') THEN ; pold='C'  ; pnew='CO'
        ELSE                           ; pold='Cd' ; pnew='CdO'
        ENDIF
        CALL swap(tempkg,pold,tgroup(ia),pnew)
      ENDIF      
      tbond(ib,ia)=0  ;  tbond(ia,ib)=0   ! break bond (2 species !)

! linear decomposition: case of R-O-C(O.)<
      IF (tgroup(ib)(1:3)=='-O-') THEN
        oloop: DO j=1,nca
          IF ((tbond(ib,j)==3).AND.(j/=ia)) THEN
            nc=INDEX(tgroup(j),' ') ; tgroup(j)(nc:nc+3)='(O.)'
            tgroup(ib)=' '
            tbond(ib,j)=0  ;  tbond(j,ib)=0
! If 2 ether groups => must be converted as RO2: R-O--O-C(O.)< --> >C=O + R(OO.)
            IF (tgroup(j)=='-O-(O.)') THEN
              DO k=1,nca
                IF ((tbond(j,k)/=0).AND.(k/=ib)) THEN
                  tbond(k,j)=0  ; tbond(j,k)=0  ;  tgroup(j)=' '
                  nc=INDEX(tgroup(k),' ')  ;  tgroup(k)(nc:nc+4)='(OO.)'
                ENDIF
              ENDDO
            ENDIF
            EXIT oloop
          ENDIF
        ENDDO oloop  
! linear decomposition: other cases
      ELSE
        nc=INDEX(tgroup(ib),' ') ;  tgroup(ib)(nc:nc)='.'   ! add dot to R fragment
      ENDIF

      CALL fragm(tbond,tgroup,pdec(ndec,1),pdec(ndec,2))

! check if product is substituted alkene (and decompose it if necessary)
      DO k=1,2
        IF (INDEX(pdec(ndec,k),'.')==0) THEN
          IF (INDEX(pdec(ndec,k),'Cd(O')/=0  .OR. INDEX(pdec(ndec,k),'CdH(O')/=0) THEN
            CALL alkcheck(pdec(ndec,k),acopd,acom)
            IF (acom/=' ') CALL addref(progname,acom,nrfdec(ndec),rfdec(ndec,:))
            IF (acopd/=' ') THEN
              mesg="acopd not empty "
              CALL stoperr(progname,mesg,chem)
            ENDIF  
          ENDIF
        ENDIF
      ENDDO

! ------------------------------   
! CYCLIC MOLECULE DECOMPOSITION   
! ------------------------------   
    ELSE
! regular cyclic opening
!------------------------
      tbond(ib,ia)=0  ;  tbond(ia,ib)=0  ! break bond matrix, opening the ring

! if group(ib) is -O-, it must be converted to (O.)
      IF (tgroup(ib)=='-O-') THEN
        DO j=1,nca
          IF ((tbond(ib,j)/=0).AND.(j/=ia)) THEN
            tbond(ib,j)=0  ;  tbond(j,ib)=0  ;  tgroup(ib)=' '
            nc=INDEX(tgroup(j),' ') ;  tgroup(j)(nc:nc+3)='(O.)'
! if group(j) is another -O- => convert to (OO.) 
            IF (tgroup(j)=='-O-(O.)') THEN
              DO k=1,nca
                IF ((tbond(j,k)/=0).AND.(k/=ib)) THEN
                  tbond(k,j)=0  ;  tbond(j,k)=0  ;  tgroup(j)=' '
                  nc=INDEX(tgroup(k),' ') ; tgroup(k)(nc:nc+4)='(OO.)'
                ENDIF
              ENDDO
            ENDIF
            fldecether=1
          ENDIF
        ENDDO
      ELSE  
        nc=INDEX(tgroup(ib),' ')  ;  tgroup(ib)(nc:nc)='.'
      ENDIF

! change (O.) to carbonyl:
      IF (group(ia)(1:2)=='CH' ) THEN
        pold='(O.)' ; pnew=' '   ; CALL swap(tgroup(ia),pold,tempkg,pnew)
        pold='CH'   ; pnew='CHO' ; CALL swap(tempkg,pold,tgroup(ia),pnew)
      ELSE IF (group(ia)(1:1)=='C'  ) THEN
        pold='(O.)' ; pnew=' '   ; CALL swap(tgroup(ia),pold,tempkg,pnew)
        pold='C'  ;  pnew='CO'   ; CALL swap(tempkg,pold,tgroup(ia),pnew)
      ENDIF      

! standardise formula of new linear molecule
      CALL rebond(tbond,tgroup,tpchem,tpnring)
      CALL stdchm(tpchem)
      pdec(ndec,1)=tpchem
! radchk gets called later (a couple of screens-worth down)
    ENDIF

! -------------------------
! compute the rate constant
! -------------------------
    IF (kdiss_sar==2) THEN
      CALL krodec_ver(bond,group,nca,ia,ib,nring,arrhdec(ndec,:))
      CALL addref(progname,'LV09KER000',nrfdec(ndec),rfdec(ndec,:))
    ELSE
      CALL krodec_atk(chem,pdec(ndec,1),pdec(ndec,2), &
                   tgroup(ib),copdec(ndec,1),fldecether,arrhdec(ndec,:))      
      CALL addref(progname,'RA07KER000',nrfdec(ndec),rfdec(ndec,:))
    ENDIF
 
! -------------------------
! standardize the products
! -------------------------
    DO k=1,2
      IF (pdec(ndec,k)==' ') CYCLE
      IF (INDEX(pdec(ndec,k),'.')/=0) THEN
        tpchem=pdec(ndec,k)
        CALL radchk(tpchem,rdckpd,rdckcopd,nip,sc,nrfdec(ndec),rfdec(ndec,:))
        pdec(ndec,k)=rdckpd(1)  ;  sdec(ndec,1)=sc(1)
        IF (nip==2) THEN
          fldec2(ndec)=1
          pdec2(ndec)=rdckpd(2)  ;  sdec(ndec,2)=sc(2)
          copdec2(ndec,:)=rdckcopd(2,:)
        ENDIF

        DO j=1,SIZE(copdec,2)
          IF (copdec(ndec,j)==' ') EXIT
        ENDDO
        DO l=1,SIZE(rdckcopd,2)
          IF (rdckcopd(1,l)/=' ') THEN
            copdec(ndec,j)=rdckcopd(1,l)  ; j=j+1
          ENDIF
        ENDDO
      ENDIF
      CALL stdchm(pdec(ndec,k))  ! required to std non radical products
    ENDDO


! -------------------------
! check non-radical fragment
! -------------------------
! If pdec is R-C(=O)ONO2, assume decomposition to R. + CO2 + NO2 
    DO i=1,2
      IF (INDEX(pdec(ndec,i),'.')/=0) CYCLE
      CALL coono2_dec(pdec(ndec,i),copdec(ndec,:),nca,nrfdec(ndec),rfdec(ndec,:))
    ENDDO 

! restore tables
    tgroup(:)=group(:)  ;  tbond(:,:)=bond(:,:)

  ENDDO grloop

! ------------------------------
! check for duplicate reactions:
! ------------------------------

! decomposition (keep only 5 numbers after comma before comparing ar3dec)
  DO i=1,ndec
    WRITE(resu,'(e12.5)') arrhdec(i,3)  ;  READ(resu,'(e12.5)') arrhdec(i,3)
  ENDDO

  DO i=1,ndec-1
    dflloop: DO j=i+1,ndec
      IF ( (pdec(j,1)==pdec(i,1) .AND. pdec(j,2)==pdec(i,2)).OR. &
           (pdec(j,1)==pdec(i,2) .AND. pdec(j,2)==pdec(i,1))) THEN
        DO k=1,SIZE(copdec,2)
           IF (copdec(i,k)/=copdec(j,k)) CYCLE dflloop
        ENDDO
        IF (arrhdec(i,1)/=arrhdec(j,1)) CYCLE dflloop
        IF (arrhdec(i,3)/=arrhdec(j,3)) CYCLE dflloop
        IF (fldec(i)==0) CYCLE dflloop
        fldec(j)=0  ;  fldec(i)=fldec(i)+1     ! identical reactions
      ENDIF
    ENDDO dflloop
  ENDDO
  DO i=1,ndec
    IF (fldec(i)>1) THEN
      arrhdec(i,1)=REAL(fldec(i))*arrhdec(i,1)  ;  fldec(i)=1
    ENDIF
  ENDDO
END SUBROUTINE rodec

!=======================================================================
! PURPOSE: manage the decomposition of alkoxy when no pathway was found: 
! look for C-O bond breaking such as C(O.)OX -> C=O + OX
!=======================================================================
SUBROUTINE rodecO(bond,group,nca,ia,ndec,fldec,pdec,copdec, &
               sdec,arrhdec,nrfdec,rfdec)
  USE reactool, ONLY: swap, rebond
  USE normchem, ONLY: stdchm
  USE toolbox, ONLY: addref
  USE dectool, ONLY: coono2_dec
  IMPLICIT NONE

  INTEGER,INTENT(IN)            :: bond(:,:)
  CHARACTER(LEN=*),INTENT(IN)   :: group(:)
  INTEGER,INTENT(IN)            :: nca         ! # of nodes
  INTEGER,INTENT(IN)            :: ia          ! index for node bearing (O.) group
  INTEGER,INTENT(OUT)           :: ndec        ! # of decomp reactions
  INTEGER,INTENT(OUT)           :: fldec(:)    ! decomp reaction flag  
  CHARACTER(LEN=*),INTENT(OUT)  :: pdec(:,:)
  CHARACTER(LEN=*),INTENT(OUT)  :: copdec(:,:)
  REAL,INTENT(OUT)              :: sdec(:,:)
  REAL,INTENT(OUT)              :: arrhdec(:,:)
  INTEGER,INTENT(INOUT)         :: nrfdec(:)     ! # of reference added for each reactions
  CHARACTER(LEN=*),INTENT(INOUT):: rfdec(:,:)    ! ref/comments for each reactions

  INTEGER :: tbond(SIZE(bond,1),SIZE(bond,2))
  CHARACTER(LEN=LEN(group(1)))  :: tgroup(SIZE(group))
  CHARACTER(LEN=LEN(group(1)))  :: pold, pnew,tempkg
  CHARACTER(LEN=LEN(pdec(1,1))) :: tpchem
  INTEGER :: tpnring

  CHARACTER(LEN=9),PARAMETER :: progname='rodecO '

  ndec=0  ; arrhdec(:,:)=0.  ;  sdec(:,:)=0.
  fldec(:)=0   ;  pdec(:,:)=' ' ; copdec(:,:)=' '
  tbond(:,:)=bond(:,:) ; tgroup(:)=group(:) 
 
  IF (INDEX(group(ia),'(O')==0) RETURN

  ndec=1 ; fldec(ndec) = 1 ; nrfdec(:)=0

! change (O.) to carbonyl:
    IF (group(ia)(1:2)=='CH') THEN
      pold= '(O.)' ; pnew=' '   ; CALL swap(group(ia),pold,tempkg,pnew)
      pold='CH'    ; pnew='CHO' ; CALL swap(tempkg,pold,tgroup(ia),pnew)
    ELSE IF (group(ia)(1:1)=='C'  ) THEN
      pold='(O.)'  ; pnew=' '   ; CALL swap(group(ia),pold,tempkg,pnew)
      pold='C'     ; pnew='CO'  ; CALL swap(tempkg,pold,tgroup(ia),pnew)
    ELSE
      PRINT*, "error in rodecO, unknow ro type"
      STOP "in rodecO"
    ENDIF      

! add coproduct
    IF (INDEX(tgroup(ia),'(OH)')/=0) THEN 
      pold='(OH)' ; copdec(ndec,1)='HO'
    ELSE IF (INDEX(tgroup(ia),'(OOH)')/=0) THEN 
      pold='(OOH)' ; copdec(ndec,1)='HO2'
    ELSE IF (INDEX(tgroup(ia),'(NO2)')/=0) THEN 
      pold='(NO2)' ; copdec(ndec,1)='NO2'
    ELSE IF (INDEX(tgroup(ia),'(ONO2)')/=0) THEN 
      pold='(ONO2)' ; copdec(ndec,1)='NO3'
    ENDIF
    pnew=' '
    
    tempkg=tgroup(ia)
    CALL swap(tempkg,pold,tgroup(ia),pnew)
    CALL rebond(tbond,tgroup,tpchem,tpnring)
    
    pdec(ndec,1)=tpchem  ;  pdec(ndec,2)=' '
    sdec(ndec,1)=1.

    arrhdec(ndec,1)=1E10 ; arrhdec(ndec,2)=0. ; arrhdec(ndec,3)=0.

    CALL stdchm(pdec(ndec,1))

    CALL addref(progname,'RO_DEC_C-O',nrfdec(ndec),rfdec(ndec,:))

! If pdec is R-C(=O)ONO2, assume decomposition to R. + CO2 + NO2 
    CALL coono2_dec(pdec(ndec,1),copdec(ndec,:),nca,nrfdec(ndec),rfdec(ndec,:))
    
END SUBROUTINE rodecO

!=======================================================================
! PURPOSE: manage the decomposition of some beta unsaturated alkoxy 
! radicals, in particular: 
!     >C=CR-CR(polar)-C(O.)< or >C=CR-CO-C(O.) 
! This R(O.) are assumed to decompose only.
!=======================================================================
SUBROUTINE rospedec(bond,group,nca,ia,ib,nring,ndec,fldec,pdec,copdec, &
                    fldec2,pdec2,copdec2,sdec,arrhdec,nrfdec,rfdec)
  USE keyflag, ONLY: kdiss_sar          ! flag for the estimation method   
  USE ringtool, ONLY: findring 
  USE reactool, ONLY: swap, rebond
  USE normchem, ONLY: stdchm
  USE radchktool, ONLY: radchk    
  USE fragmenttool, ONLY: fragm
  USE toolbox, ONLY: addref
  USE dectool, ONLY: coono2_dec
  IMPLICIT NONE

  INTEGER,INTENT(IN)  :: bond(:,:)
  CHARACTER(LEN=*),INTENT(IN) :: group(:)
  INTEGER,INTENT(IN)  :: nca     ! # of nodes
  INTEGER,INTENT(IN)  :: ia      ! index for node bearing (O.) group
  INTEGER,INTENT(IN)  :: ib      ! index for beta node bearing polar group 
  INTEGER,INTENT(IN)  :: nring   !     
  INTEGER,INTENT(OUT) :: ndec        ! # of decomp reactions
  INTEGER,INTENT(OUT) :: fldec(:)    ! decomp reaction flag  
  INTEGER,INTENT(OUT) :: fldec2(:)   ! decomp reaction flag (2nd spe) 
  CHARACTER(LEN=*),INTENT(OUT) :: pdec(:,:)    ! decomposition products 
  CHARACTER(LEN=*),INTENT(OUT) :: copdec(:,:)  ! decomposition co-products
  CHARACTER(LEN=*),INTENT(OUT) :: pdec2(:)     ! decomposition products (2nd spe)
  CHARACTER(LEN=*),INTENT(OUT) :: copdec2(:,:) ! decomposition co-products (2nd spe)
  REAL,INTENT(OUT)    :: sdec(:,:)             ! "stoi. coef." (spe1 vs spe2)
  REAL,INTENT(OUT)    :: arrhdec(:,:)
  INTEGER,INTENT(INOUT) :: nrfdec(:)          ! # of reference added for each reactions
  CHARACTER(LEN=*),INTENT(OUT) :: rfdec(:,:)  ! ref/comments

  INTEGER :: tbond(SIZE(bond,1),SIZE(bond,2))
  CHARACTER(LEN=LEN(group(1)))  :: tgroup(SIZE(group))
  CHARACTER(LEN=LEN(group(1)))  :: pold, pnew,tempkg
  CHARACTER(LEN=LEN(pdec(1,1))) :: tpchem, tpchemb
  INTEGER :: ring(SIZE(group))
  INTEGER :: rngflg,tpnring,nip,nc,irad,imol,i

  INTEGER,PARAMETER :: mxrpd=2  ! max # of product returned by radchk sub
  CHARACTER(LEN=LEN(pdec(1,1)))   :: rdckpd(mxrpd)
  CHARACTER(LEN=LEN(copdec(1,1))) :: rdckcopd(mxrpd,SIZE(copdec,2))
  REAL :: sc(mxrpd)

  CHARACTER(LEN=9),PARAMETER :: progname='rospedec '

  ndec=0  ; arrhdec(:,:)=0.  ;  sdec(:,:)=0.  ; ring(:)=0
  fldec(:)=0   ;  pdec(:,:)=' ' ; copdec(:,:)=' '
  fldec2(:)=0  ;  pdec2(:)=' '  ; copdec2(:,:)=' '
  tbond(:,:)=bond(:,:) ; tgroup(:)=group(:) 
 
  ndec=1  ;  fldec(ndec) = 1  ; nrfdec(:)=0

! change (O.) to carbonyl:
    IF (group(ia)(1:3)=='CH2') THEN
      tgroup(ia)='CH2O'
    ELSE IF (group(ia)(1:2)=='CH') THEN
      pold= '(O.)' ; pnew=' '   ; CALL swap(group(ia),pold,tempkg,pnew)
      pold='CH'    ; pnew='CHO' ; CALL swap(tempkg,pold,tgroup(ia),pnew)
    ELSE IF (group(ia)(1:1)=='C'  ) THEN
      pold='(O.)'  ; pnew=' '   ; CALL swap(group(ia),pold,tempkg,pnew)
      pold='C'     ; pnew='CO'  ; CALL swap(tempkg,pold,tgroup(ia),pnew)
    ELSE
      PRINT*, "error in rospedec, unknow ro type"
      STOP "in rospedecs"
    ENDIF      

! check if bond that will be broken belongs to a ring & break bonds
    CALL findring(ib,ia,nca,tbond,rngflg,ring)
    tbond(ib,ia)=0  ;  tbond(ia,ib)=0

! add dot to alkyl radical fragment:
    nc=INDEX(group(ib),' ')  ;  tgroup(ib)(nc:nc)='.'

    IF (rngflg==1) CALL rebond(tbond,tgroup,tpchem,tpnring)
    IF (rngflg==0) THEN
      CALL fragm(tbond,tgroup,pdec(ndec,1),pdec(ndec,2))
    ELSE 
      pdec(ndec,1)=tpchem  ;  pdec(ndec,2)=' '
    ENDIF            

    IF (kdiss_sar==2) THEN
      CALL krodec_ver(bond,group,nca,ia,ib,nring,arrhdec(ndec,:))
    ELSE
      arrhdec(ndec,1)=2.0E14 ; arrhdec(ndec,2)=0.  ! rate constant setto k=2E14exp(-Ea/RT), Ea = 0.75
      arrhdec(ndec,3)=0.75   ; arrhdec(ndec,3)=(arrhdec(ndec,3)*1000.)/1.9872
    ENDIF

! identify radical fragment, run through radical checker routine
    IF (INDEX(pdec(ndec,1),'.')/=0) THEN ; irad = 1 ; imol = 2
    ELSE                                 ; irad = 2 ; imol = 1
    ENDIF

    tpchem=pdec(ndec,irad)
    CALL radchk(tpchem,rdckpd,rdckcopd,nip,sc,nrfdec(ndec),rfdec(ndec,:))
    pdec(ndec,irad)=rdckpd(1)  ;  sdec(ndec,1)=sc(1)
    IF (nip==2) THEN
      fldec2(ndec)=1
      pdec2(ndec)=rdckpd(2)  ;  sdec(ndec,2)=sc(2)
      copdec2(ndec,:)=rdckcopd(2,:)
    ENDIF

!    CALL stdchm(pdec(ndec,irad))  ! useless: radchk return std formula
    CALL stdchm(pdec(ndec,imol))

! always put radical species in 2nd, if multiple delocalisation pdcts
    tpchem=pdec(ndec,irad)  ;  tpchemb=pdec(ndec,imol) 
    pdec(ndec,1)=tpchemb    ;  pdec(ndec,2)=tpchem
    copdec(ndec,:)=rdckcopd(1,:)
    CALL addref(progname,'RO_BPOL_U',nrfdec(ndec),rfdec(ndec,:))

! If pdec is R-C(=O)ONO2, assume decomposition to R. + CO2 + NO2 
    DO i=1,2
      IF (INDEX(pdec(ndec,i),'.')/=0) CYCLE
      CALL coono2_dec(pdec(ndec,i),copdec(ndec,:),nca,nrfdec(ndec),rfdec(ndec,:))
    ENDDO
    
END SUBROUTINE rospedec

!=======================================================================
! PURPOSE: Find rate constants for R(O.) decomposition based on 
! Vereecken and Peeters 2009 SAR. 
!=======================================================================
SUBROUTINE krodec_ver(tbond,tgroup,nca,ia,ib,nring,arrh)
  USE keyparameter, ONLY:mxcp
  USE mapping, ONLY: gettrack
  USE ringtool, ONLY: findring
  IMPLICIT NONE

  INTEGER,INTENT(IN) :: tbond(:,:)
  CHARACTER(LEN=*),INTENT(IN) :: tgroup(:)
  INTEGER,INTENT(IN) :: nca   ! # of nodes
  INTEGER,INTENT(IN) :: ia    ! group bearing the (O.) center
  INTEGER,INTENT(IN) :: ib    ! bond to be broken is (ia,ib) 
  INTEGER,INTENT(IN) :: nring ! # of rings in the parent compound
  REAL,INTENT(OUT)   :: arrh(3)

  INTEGER :: i,j,k
  INTEGER :: ring(SIZE(tgroup))   ! rg index for nodes, current ring (0=no, 1=yes)
  INTEGER :: rngflg               ! 0 = 'no ring', 1 = 'yes ring'
  INTEGER :: lring(5)
  INTEGER :: track(mxcp,SIZE(tgroup))
  INTEGER :: trlen(mxcp)
  INTEGER :: ntr
  REAL    :: Eb                   ! barrier height (kcal.mol-1) for RO decomposition

  lring(:)=0  ;  ring(:)=0
  Eb= 17.9

! IF RINGS EXIST, THEN FIND THE NODES THAT BELONG TO THE RINGS
! --------------------------------------------------------------
! check if ring opening required
  IF (nring>0) THEN
    CALL findring(ia,ib,nca,tbond,rngflg,ring)

! if ia-ib bond is not part of a cycle, check if ia or ib is a member of a cycle
    IF (rngflg==0) THEN
      DO i=1,nca
        IF ((tbond(i,ia)==0).OR.(i==ib)) CYCLE
        CALL findring(i,ia,nca,tbond,rngflg,ring)
        IF (rngflg==1) EXIT
      ENDDO
      IF (rngflg/=1) THEN 
        DO i=1,nca
          IF ((tbond(i,ib)==0).OR.(i==ia)) CYCLE
          CALL findring(i,ib,nca,tbond,rngflg,ring)
          IF (rngflg==1) EXIT
        ENDDO
      ENDIF 
      rngflg=0   ! reset rngflg to 0 as ia-ib not part of a cycle
    ENDIF
  ENDIF

  IF (nring>0) THEN
    grloop1: DO i=1,nca
      k=1
      CALL gettrack(tbond,i,nca,ntr,track,trlen)
      DO j=1,ntr
        IF ((trlen(j)>2).AND.(tbond(i,track(j,trlen(j)))>0)) THEN
          IF (trlen(j)/=lring(k)) THEN 
            lring(k)=trlen(j)
            k=k+1
            IF (k==nring) EXIT grloop1  ! to be updated for species having 3 cycles
          ENDIF
        ENDIF
      ENDDO
    ENDDO grloop1
  ENDIF

! SUBSTITUENTS ON THE CARBON BEARING THE ALKOXY RADICAL 
! -----------------------------------------------------
  IF (tgroup(ia)(1:6)=='CO(O.)')     Eb=Eb-12.7
  IF (INDEX(tgroup(ia),'(OH)')/=0)   Eb=Eb-8.9
  IF (INDEX(tgroup(ia),'(OOH)')/=0)  Eb=Eb-8.9
  IF (INDEX(tgroup(ia),'(NO2)')/=0)  Eb=Eb-2.2
  IF (INDEX(tgroup(ia),'(ONO)')/=0)  Eb=Eb-4.2
  IF (INDEX(tgroup(ia),'(ONO2)')/=0) Eb=Eb-3.8

  DO j=1,nca
    IF ((tbond(ia,j)/=0).AND.(j/=ib)) THEN                      
      IF (tgroup(j)(1:3)=='-O-')  THEN
        DO k=1,nca
          IF ((tbond(k,j)/=0).AND.(k/=ia)) THEN
            IF (tgroup(k)(1:3)=='-O-') THEN ; Eb=Eb-6.4 ! alpha-O-O-R
            ELSE                            ; Eb=Eb-9.2 ! alpha-O-R
            ENDIF
          ENDIF
        ENDDO
      ENDIF       
    ENDIF   
  ENDDO

  DO j=1,nca
    IF (ring(j)==1) CYCLE   ! cycle nodes: contribution included in cyclic factors
    IF ((tbond(ia,j)==2).AND.(j/=ib)) THEN                      
      IF (tgroup(ia)(1:2)=='Cd')  Eb=Eb + 21.5    ! alpha=C

    ELSE IF ((tbond(ia,j)==1).AND.(j/=ib)) THEN                 
      IF (tgroup(j)(1:2)=='Cd')   Eb=Eb-4.9       ! alpha-C=C
      IF (tgroup(j)(1:2)=='CH')   Eb=Eb-2.3       ! alpha-alkyl
      IF (tgroup(j)(1:2)=='C ')   Eb=Eb-2.3       ! alpha-alkyl
      IF (tgroup(j)(1:3)=='C(O')  Eb=Eb-2.3       ! alpha-alkyl

! correction factor if there is only one substituent to the alpha carbon
      IF (tgroup(j)(1:4)=='CHO ') THEN
        IF (tgroup(ia)(1:6)=='CH(O.)')    Eb=Eb-0.7+2.3  ! R1-CH(O.)CHO
      ENDIF
      IF (tgroup(j)(1:8)=='CH2(OOH)') THEN
        IF (tgroup(ia)(1:6)=='CH(O.)')    Eb=Eb-0.7+2.3  ! R1-CH(O.)CH2(OOH)
      ENDIF
      IF (tgroup(j)(1:3)=='CO ') THEN
        IF (tgroup(ia)(1:6)=='CH(O.)') THEN ; Eb=Eb-0.7  ! R1-CH(O.)CO-R2
        ELSE                                ; Eb=Eb-2.3  ! R1-C(O.)R3CO-R2
        ENDIF
      ENDIF
      DO k=1,nca
        IF ((tbond(k,j)/=0).AND.(k/=ia)) THEN
          IF ((tgroup(j)(1:4)=='CH2 ').AND.(tgroup(k)(1:3)=='-O-')) THEN
            Eb=Eb-0.7+2.3                                ! alpha-CH2-O-R
          ENDIF
        ENDIF
      ENDDO
    ENDIF
  ENDDO

! SUBSTITUENTS ON THE LEAVING CARBON 
! ----------------------------------
  IF (tgroup(ib)(1:3)=='CO ')        Eb=Eb-8.5
  IF (tgroup(ib)(1:4)=='CHO ')       Eb=Eb-8.5
  IF (INDEX(tgroup(ib),'(OH)')/=0)   Eb=Eb-7.5
  IF (INDEX(tgroup(ib),'(OOH)')/=0)  Eb=Eb-9.3
  IF (INDEX(tgroup(ib),'(NO)')/=0)   Eb=Eb-16.0
  IF (INDEX(tgroup(ib),'(NO2)')/=0)  Eb=Eb+0.4
  IF (INDEX(tgroup(ib),'(ONO)')/=0)  Eb=Eb-6.0
  IF (INDEX(tgroup(ib),'(ONO2)')/=0) Eb=Eb-2.8

  DO j=1,nca
    IF ((tbond(ib,j)/=0).AND.(j/=ia)) THEN                      
      IF (tgroup(j)(1:3)=='-O-')  THEN
        DO k=1,nca
          IF ((tbond(k,j)/=0).AND.(k/=ib)) THEN
            IF (tgroup(k)(1:3)=='-O-') THEN ; Eb=Eb-7.2  ! alpha-O-O-R
            ELSE                            ; Eb=Eb-9.1  ! alpha-O-R
            ENDIF
          ENDIF
        ENDDO
      ENDIF  
    ENDIF
  ENDDO        

  DO j=1,nca
    IF (ring(j)==1) CYCLE  ! cycle nodes: contribution included in cyclic factors
    IF ((tbond(ib,j)==2).AND.(j/=ia)) THEN                      
      IF (tgroup(ib)(1:2)=='Cd') Eb=Eb+5    ! alpha=C
    ELSE IF ((tbond(ib,j)==1).AND.(j/=ia)) THEN                 
      IF (tgroup(j)(1:2)=='Cd')  Eb=Eb-9.6  ! beta-C=C
      IF (tgroup(j)(1:2)=='CH')  Eb=Eb-3.4  ! beta-alkyl
      IF (tgroup(j)(1:2)=='C ')  Eb=Eb-3.4  ! beta-alkyl
      IF (tgroup(j)(1:3)=='C(O') Eb=Eb-3.4  ! beta-alkyl
    ENDIF
  ENDDO

! CORRECTIONS DUE TO RINGS 
! -------------------------
  IF (rngflg==0) THEN
    IF (ring(ia)==1) Eb=Eb-2.0
    DO j=1,nring
      IF ((ring(ib)==1).AND.(lring(j)==3)) Eb=Eb+2.4
      IF ((ring(ib)==1).AND.(lring(j)==4)) Eb=Eb-4.2
      IF ((ring(ib)==1).AND.(lring(j)==5)) Eb=Eb-7.0
      IF ((ring(ib)==1).AND.(lring(j)>5))  Eb=Eb-7.0
    ENDDO
  ELSE IF (rngflg==1) THEN
    DO j=1,nring
      IF (lring(j)==3) Eb=Eb-24.6
      IF (lring(j)==4) Eb=Eb-17.1
      IF (lring(j)==5) Eb=Eb-8.7
      IF (lring(j)>5)  Eb=Eb-6.3
    ENDDO
  ENDIF

! correction for Eb<7kcal.mol-1
  IF (Eb<7.0) Eb=19*exp(-((Eb-22.)*(Eb-22.))/225.)

! Eb/R  : Eb in cal.mol-1 and R in Kcal.mol-1
  Eb=(Eb*1000.)/1.9872    

! set arrhenius coefficient
  arrh(1)=1.12E09  ;  arrh(2)=1.7  ;  arrh(3)=Eb
END SUBROUTINE krodec_ver


!=======================================================================
! PURPOSE: compute rate constants for R(O.) decomposition based on the 
! SAPRC99 SAR, as updated by Atkinson 2007. Rate is k=2E14exp(-Ea/RT), 
! where Ea is given as: Ea=Eaa+0.44DHr=-8.73+2.35IP+0.44DHr with 
! (1) DHr: heat of the reaction (in kcal) et (2) IP: ionization 
! potential (in eV). Values for Eaa (in kcal) are taken from the SAPRC99 
! documentation with updated Eea taken from Atkinson 2007.
!=======================================================================
SUBROUTINE krodec_atk(chem,pdct1,pdct2,lgrp,copd1,oflg,arrh)
  USE bensontool, ONLY: heat
  IMPLICIT NONE
  
  CHARACTER(LEN=*),INTENT(IN) :: chem   ! formula of the alkoxy
  CHARACTER(LEN=*),INTENT(IN) :: pdct1  ! formula of 1st product
  CHARACTER(LEN=*),INTENT(IN) :: pdct2  ! formula of 2nd product
  CHARACTER(LEN=*),INTENT(IN) :: lgrp   ! leaving grp in bond breaking
  CHARACTER(LEN=*),INTENT(IN) :: copd1  ! 1st coproduct (check for CO2)
  INTEGER,INTENT(IN) :: oflg            ! ether flag (leaving grp is R-O-)
  REAL,INTENT(OUT) :: arrh(3)

  REAL :: dhfp1,dhfp2,dhfro,dhnet 
  REAL :: Eaa, ea
          
! compute heat of formation of each product
  dhfp1=heat(pdct1)  ;   dhfp2=heat(pdct2)  ;  dhfro=heat(chem)
  dhnet=dhfp1+dhfp2-dhfro

! set the Eaa factor
  IF (lgrp(1:4)=='CH3.') THEN        ; Eaa=12.9
  ELSE IF(lgrp(1:3)=='CH2') THEN
    IF (INDEX(lgrp,'(OH)')/=0) THEN  ; Eaa=6.8
    ELSE                             ; Eaa=9.5
    ENDIF
  ELSE IF(lgrp(1:2)=='CH') THEN
    IF (lgrp(1:3)=='CHO') THEN             ; Eaa=9.99
    ELSE IF (INDEX(lgrp,'(OH)')/=0) THEN   ; Eaa=5.2
    ELSE IF (INDEX(lgrp,'(ONO2)')/=0) THEN ; Eaa=7.46
    ELSE IF (INDEX(lgrp,'(OOH)')/=0) THEN  ; Eaa=7.46
    ELSE                                   ; Eaa=8.2
    ENDIF
  ELSE IF(lgrp(1:2)=='CO') THEN            ; Eaa=5.
  ELSE IF(lgrp(1:1)=='C') THEN
    IF (INDEX(lgrp,'(OH)')/=0) THEN        ; Eaa=4.8
    ELSE IF (INDEX(lgrp,'(ONO2)')/=0) THEN ; Eaa=5.58
    ELSE IF (INDEX(lgrp,'(OOH)')/=0) THEN  ; Eaa=5.58
    ELSE                                   ; Eaa=7.1
    ENDIF
  ELSE IF (oflg==1) THEN                   ; Eaa=7.5
  ELSE IF ((lgrp(1:1)==' ').AND.(copd1=='CO2')) THEN   ; Eaa=5.
  ELSE
    WRITE(6,'(a)') 'lgrp=',TRIM(lgrp)
    WRITE(6,'(a)') 'copdec=',copd1
    WRITE(6,'(a)') '--error--, in krodec_atk. Eaa factor not attributed'
    WRITE(6,'(a)') 'for the decomposition  of the species :'
    WRITE(6,'(a)') TRIM(chem)
    STOP "in krodec_atk"
  ENDIF
  
! set the arrhenius parameter
  arrh(1)=5.0E13   ;  arrh(2)=0.  
  ea=(Eaa+0.40*dhnet)  ; ea=MAX(0.75,ea)  ;  arrh(3)=(ea*1000.)/1.9872

END SUBROUTINE krodec_atk

END MODULE rodectool
