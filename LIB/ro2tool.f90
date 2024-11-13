MODULE ro2tool
IMPLICIT NONE
CONTAINS
! SUBROUTINE ro2aro(chem,pname,bond,group,ngr,nring,brch,ip,Ncon,loaro)
! SUBROUTINE ro2no(chem,pname,bond,group,ngr,nring,brch,ip,Ncon)
! REAL FUNCTION nityield(bond,group,ngr,nring,Ncon,ip,tflg)
! SUBROUTINE ro2ho2(chem,pname,bond,group,ngr,brch,ip,Ncon)
! SUBROUTINE ro2no3(chem,pname,bond,group,brch,ip)
! SUBROUTINE ro2ro2(chem,pname,bond,group,ngr,brch,ip,Ncon)
! SUBROUTINE ro2ho(chem,pname,bond,group,brch,ip)
!
!=======================================================================
! PURPOSE: perform ring closure for aromatic RO2. It is assumed that 
! ring closure is fast enough to not compete with other RO2 reactions.
!
! Reaction rates needs to be set (although not a key issue, competing
! reactions being ignored. 
!
! See Jenkin et al., ACP, 2018 and Jenkin et al., ACP 2019 for 
! discussion about the ring closure reactions
!=======================================================================
SUBROUTINE ro2aro(chem,pname,bond,group,ngr,brch,ip,loaro)
  USE keyparameter, ONLY: mxlco,mxpd,mxcp,mecu,mxcopd
  USE references, ONLY:mxlcod
  USE reactool, ONLY: swap,rebond
  USE dictstacktool, ONLY: bratio
  USE mapping, ONLY: gettrack
  USE ringtool, ONLY: findring
  USE radchktool, ONLY: radchk
  USE normchem, ONLY: stdchm
  USE rxwrttool, ONLY:rxwrit,rxinit
  USE toolbox, ONLY: stoperr,add1tonp,addrx,setbond,addref
  USE tempo, ONLY: chknadd2p

  CHARACTER(LEN=*),INTENT(IN) :: chem  ! formula of the VOC considered
  CHARACTER(LEN=*),INTENT(IN) :: pname ! short name of chem (parent compound) 
  CHARACTER(LEN=*),INTENT(IN) :: group(:) ! group (without ring joining char)
  INTEGER,INTENT(IN) :: bond(:,:)     ! bond matrix
  INTEGER,INTENT(IN) :: ngr           ! # of nodes
  REAL,   INTENT(IN) :: brch          ! max yield of the input species
  INTEGER,INTENT(IN) :: ip            ! node bearing the peroxy (OO.) group
  LOGICAL,INTENT(OUT):: loaro         ! true if chem is aromatic RO2

  CHARACTER(LEN=LEN(group(1))) :: tgroup(SIZE(group))
  CHARACTER(LEN=LEN(group(1))) :: pold,pnew
  INTEGER         :: tbond(SIZE(bond,1),SIZE(bond,2))
  CHARACTER(LEN=LEN(chem)) :: tempfo,p1chem,p2chem

  INTEGER :: i,j,k,np,tempring,nc,ngrp
  REAL    :: brtio

  INTEGER :: ring(SIZE(bond,1)),rngflg
  INTEGER :: track(mxcp,SIZE(bond,1))
  INTEGER :: trlen(mxcp)
  INTEGER :: ntr

  INTEGER,PARAMETER    :: mxrpd=2   ! max # of products returned by radchk sub
  CHARACTER(LEN=LEN(chem)) :: rdckpd(mxrpd)
  CHARACTER(LEN=mxlco) :: rdckcopd(mxrpd,mxcopd)
  REAL    :: sc(mxrpd)
  INTEGER :: nip

  CHARACTER(LEN=mxlco) :: p(mxpd), r(3)
  REAL    :: s(mxpd), arrh(3)
  INTEGER :: idreac, nlabel
  REAL    :: xlabel,folow(3),fotroe(4)

  INTEGER,PARAMETER     :: mxref=10
  INTEGER               :: nref
  CHARACTER(LEN=mxlcod) :: ref(mxref)

  CHARACTER(LEN=7)      :: progname='ro2aro'
  CHARACTER(LEN=70)     :: mesg

  tgroup(:)=group(:) ;  tbond(:,:)=bond(:,:)
  loaro=.FALSE.
  
! check if the peroxy belong to "aromatic RO2" class (C1(OO.)Cd2Cd3Cd4Cd5C6(OH))
  CALL gettrack(tbond,ip,ngr,ntr,track,trlen)
  ngrp=ngr
  DO i=1,ngr 
    IF ((bond(ip,i)/=0).AND.(INDEX(group(i),'(OH)')/=0)) THEN
      CALL findring(ip,i,ngr,tbond,rngflg,ring)
      IF (rngflg==1) THEN
        DO j=1,ntr
          IF (trlen(j)<6) CYCLE
          IF (tgroup(track(j,2))(1:2)/='Cd') CYCLE
          IF (tgroup(track(j,3))(1:2)/='Cd') CYCLE
          IF (tgroup(track(j,4))(1:2)/='Cd') CYCLE
          IF (tgroup(track(j,5))(1:2)/='Cd') CYCLE
          IF (track(j,6)/=i) CYCLE

! === 1st reaction: make -O--O-bridge ring closure and stabilize
          CALL rxinit(r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe)
          nref=0  ;  ref(:)=' ' ;  CALL addref(progname,'MJ19KMV000',nref,ref,chem)
          np=0
       
! add and bridge new groups
          tgroup(ngrp+1)='-O-'  ;  tgroup(ngrp+2)='-O-'              
          CALL setbond(tbond,ngrp+1,ip,3)                
          CALL setbond(tbond,ngrp+1,ngrp+2,3)
          CALL setbond(tbond,ngrp+2,track(j,5),3)

! rm C=C bond at position 4 and 5 and peroxy at position 1
          pold='(OO.)'  ;  pnew=' '
          CALL swap(group(ip),pold,tgroup(ip),pnew)
          pold='Cd'  ;  pnew='C'
          CALL swap(group(track(j,4)),pold,tgroup(track(j,4)),pnew)
          nc=INDEX(tgroup(track(j,4)),' ') ;
          tgroup(track(j,4))(nc:nc)='.'
          pold='Cd'  ;  pnew='C'
          CALL swap(group(track(j,5)),pold,tgroup(track(j,5)),pnew)
          CALL setbond(tbond,track(j,4),track(j,5),1)
          
! rebond and normalize new products (2 products expected from delocalisation)
          CALL rebond(tbond,tgroup,tempfo,tempring)
          CALL radchk(tempfo,rdckpd,rdckcopd,nip,sc,nref,ref)
          CALL stdchm(rdckpd(1))
          tempfo=rdckpd(1)
          CALL add1tonp(progname,chem,np)   
          s(np)=sc(1) ; brtio=brch         ! force to highest yield
          CALL bratio(tempfo,brtio,p(np),nref,ref)
          DO k=1,SIZE(rdckcopd,2)
            IF (rdckcopd(1,k)==' ') EXIT
            CALL add1tonp(progname,chem,np)
            s(np)=sc(1) ; p(np)=rdckcopd(1,k)
          ENDDO

          IF (nip==2) THEN
            CALL stdchm(rdckpd(2))
            tempfo=rdckpd(2)
            CALL add1tonp(progname,chem,np)   
            s(np)=sc(2) ; brtio=brch         ! force to highest yield
            CALL bratio(tempfo,brtio,p(np),nref,ref)
            DO k=1,SIZE(rdckcopd,2)
              IF (rdckcopd(2,k)==' ') EXIT
              CALL add1tonp(progname,chem,np)
              s(np)=sc(2) ; p(np)=rdckcopd(2,k)
            ENDDO
          ENDIF
     
          r(1)=pname ; r(2)=' '           ! reactant
          arrh(1)=0.75 ; arrh(2)=0. ; arrh(3)=0.
          CALL rxwrit(mecu,r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe,ref)
          tbond(:,:)=bond(:,:)  ;  tgroup(:)=group(:)

!==== make an epoxy-alkoxy (break -O-O- bridge). 2 species can be formed 
          nref=0  ;  ref(:)=' '

! species one 
          tgroup(ngrp+1)='-O-'
          pold='(OO.)'  ;  pnew='(O.)'
          CALL swap(group(ip),pold,tgroup(ip),pnew)
          pold='Cd' ; pnew='C'
          CALL swap(group(track(j,4)),pold,tgroup(track(j,4)),pnew)
          pold='Cd'  ;  pnew='C'
          CALL swap(group(track(j,5)),pold,tgroup(track(j,5)),pnew)
          CALL setbond(tbond,track(j,4),track(j,5),1)
          CALL setbond(tbond,track(j,4),ngrp+1,3)
          CALL setbond(tbond,track(j,5),ngrp+1,3) 
          CALL rebond(tbond,tgroup,tempfo,tempring)
          CALL radchk(tempfo,rdckpd,rdckcopd,nip,sc,nref,ref)
          CALL stdchm(rdckpd(1))
          p1chem=rdckpd(1)
          IF ( (nip/=1) .OR. (rdckcopd(1,1)/=' ')) THEN
            mesg="unexpected 2nd product or coproducts in epoxy route"
            CALL stoperr(progname,mesg,chem)
          ENDIF
          tbond(:,:)=bond(:,:)  ; tgroup(:)=group(:)

! species two
          tgroup(ngrp+1)='-O-'
          pold='(OO.)'  ;  pnew=' '
          CALL swap(group(ip),pold,tgroup(ip),pnew)
          pold='Cd'  ;  pnew='C'
          CALL swap(group(track(j,2)),pold,tgroup(track(j,2)),pnew)
          pold='Cd'  ;  pnew='C'
          CALL swap(group(track(j,5)),pold,tgroup(track(j,5)),pnew)
          nc=INDEX(tgroup(track(j,5)),' ')
          tgroup(track(j,5))(nc:nc+3)='(O.)'
          CALL setbond(tbond,track(j,2),track(j,3),1)
          CALL setbond(tbond,track(j,3),track(j,4),2)
          CALL setbond(tbond,track(j,4),track(j,5),1)
          CALL setbond(tbond,ip,ngrp+1,3)
          CALL setbond(tbond,track(j,2),ngrp+1,3)
          CALL rebond(tbond,tgroup,tempfo,tempring)
          CALL radchk(tempfo,rdckpd,rdckcopd,nip,sc,nref,ref)
          CALL stdchm(rdckpd(1))
          p2chem=rdckpd(1)
          IF ( (nip/=1) .OR. (rdckcopd(1,1)/=' ')) THEN
            mesg="unexpected 2nd product or coproducts in epoxy route"
            CALL stoperr(progname,mesg,chem)
          ENDIF

! write the reactions
          CALL rxinit(r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe)
          CALL addref(progname,'MJ19KMV000',nref,ref,chem)
          np=0

          CALL add1tonp(progname,chem,np)   
          s(np)=0.5 ; brtio=brch         ! force to highest yield
          CALL bratio(p1chem,brtio,p(np),nref,ref)

          CALL add1tonp(progname,chem,np)   
          s(np)=0.5 ; brtio=brch         ! force to highest yield
          CALL bratio(p2chem,brtio,p(np),nref,ref)
           
          arrh(1)=0.25 ; arrh(2)=0. ; arrh(3)=0.
            
          r(1)=pname ;  r(2)='  '   ! reactant
          CALL rxwrit(mecu,r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe,ref)

          tbond(:,:)=bond(:,:) ; tgroup(:)=group(:)
          loaro=.TRUE.
          RETURN
        ENDDO

      ENDIF
    ENDIF
  ENDDO
  tbond(:,:)=bond(:,:)  ;  tgroup(:)=group(:)
END SUBROUTINE ro2aro

!=======================================================================
! PURPOSE: perform RO2+NO reactions. Two channels are considered:
! - RO2+NO => RO +NO2
! - RO2+NO => RONO2 (nitrate)
!=======================================================================
SUBROUTINE ro2no(chem,pname,bond,group,ngr,nring,brch,ip,Ncon)
  USE keyparameter, ONLY: mxlco,mxpd,mecu
  USE references, ONLY:mxlcod
  USE reactool, ONLY: swap,rebond
  USE dictstacktool, ONLY: bratio
  USE normchem, ONLY: stdchm
  USE rxwrttool, ONLY:rxwrit,rxinit
  USE toolbox, ONLY: stoperr,add1tonp,addrx,addref
  USE fragmenttool, ONLY: fragm
  USE tempo, ONLY: chknadd2p

  CHARACTER(LEN=*),INTENT(IN) :: chem  ! formula of the VOC considered
  CHARACTER(LEN=*),INTENT(IN) :: pname ! short name of chem (parent compound) 
  CHARACTER(LEN=*),INTENT(IN) :: group(:) ! group (without ring joining char)
  INTEGER,INTENT(IN) :: bond(:,:)     ! bond matrix
  INTEGER,INTENT(IN) :: ngr           ! # of nodes
  INTEGER,INTENT(IN) :: nring         ! # of ring in parent
  REAL,   INTENT(IN) :: brch          ! max yield of the input species
  INTEGER,INTENT(IN) :: ip            ! node bearing the peroxy (OO.) group
  INTEGER,INTENT(IN) :: Ncon          ! # of C, N, O atoms in peroxy "R"

  CHARACTER(LEN=LEN(group(1))) :: tgroup(SIZE(group))
  CHARACTER(LEN=LEN(group(1))) :: pold,pnew,tempkg
  INTEGER         :: tbond(SIZE(bond,1),SIZE(bond,2))
  CHARACTER(LEN=LEN(chem)) :: tempfo,p1chem,p2chem

  INTEGER :: i,np,tempring,npath,flag,nc
  REAL    :: brtio, yield,rad,rab,rodec

  CHARACTER(LEN=mxlco) :: p(mxpd), r(3)
  REAL    :: s(mxpd), arrh(3)
  INTEGER :: idreac, nlabel
  REAL    :: xlabel,folow(3),fotroe(4)

  INTEGER,PARAMETER     :: mxref=10
  INTEGER               :: nref
  CHARACTER(LEN=mxlcod) :: ref(mxref)
  CHARACTER(LEN=7)      :: progname='ro2no'

! set the reaction rate (Jenkin et al., ACP, 2019, MJ19KMV000 eqn (2)) 
  REAL,DIMENSION(3),PARAMETER :: kno=(/2.7E-12, 0., -360. /)

  tgroup(:)=group(:)  ;  tbond(:,:)=bond(:,:)

! nitrate yield
  flag=0          ! use std T and P provided as input
  IF (ngr>2) THEN
    rad=nityield(bond,group,ngr,nring,Ncon,ip,flag)
    rab=(1.-rad)
  ELSE
    rad=0. ; rab=1.
  ENDIF

! special case - remove nitrate
  IF (group(ip)(1:6)=='c(OO.)') THEN    ! No nitrate on aromatic ring
    rad=0. ; rab=1.
  ENDIF
  IF (INDEX(group(ip),'Cd')/=0) THEN    ! no nitrate on >Cd(OO.) 
    rad=0. ; rab=1.
  ENDIF   
  IF (INDEX(group(ip),'(ONO2)')/=0) THEN  ! no double nitrate 
    rad=0. ; rab=1.
  ENDIF   

! fraction of "hot" alkoxy that decompose once produced 
  rodec=0. ; npath=0
  IF ((Ncon>3).AND.(Ncon<10)) THEN
    DO i=1,ngr
      IF ( (bond(i,ip)/=0).AND.(group(ip)(1:1)/='c').AND. &
           ((INDEX(group(i),'(OH)')/=0).OR.(group(i)(1:3)=='CO ')) ) THEN
          rodec=1.34-0.145*Ncon
          npath=npath+1
      ENDIF
    ENDDO
  ENDIF

! RO2+NO => RO+NO2 reaction
! -------------------------
  CALL rxinit(r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe)
  nref=0  ;  ref(:)=' ' ;  CALL addref(progname,'MJ19KMV000',nref,ref,chem)
  np=0

! add NO2 as product
  CALL add1tonp(progname,chem,np)   
  s(np)=1. ;  p(np)='NO2 '

! change (OO.) to (O.), rebuild, check and rename
  pold='(OO.)'  ;  pnew='(O.)'
  CALL swap(group(ip),pold,tgroup(ip),pnew)
  CALL rebond(tbond,tgroup,tempfo,tempring)
  yield=1.-rodec
  brtio=brch*rab  
  brtio=brch  ! comparison old scheme
  CALL chknadd2p(chem,tempfo,yield,brch,np,s,p,nref,ref)
  tgroup(:)=group(:)  ;  tbond(:,:)=bond(:,:)
  
! Check if the formed RO* decompose
  IF (npath/=0) THEN
    DO i=1,ngr
      IF ( (bond(i,ip)/=0).AND.(group(ip)(1:1)/='c').AND. &
         ((INDEX(group(i),'(OH)')/=0).OR.(group(i)(1:3)=='CO ')) ) THEN
        tbond(i,ip)=0 ; tbond(ip,i)=0
        pold='(OO.)' ; pnew=' '
        CALL swap(tgroup(ip),pold,tempkg,pnew)
        tgroup(ip)=tempkg
        IF      (tgroup(ip)(1:3)=='CH2') THEN ; pold='CH2' ; pnew='CH2O'
        ELSE IF (tgroup(ip)(1:2)=='CH')  THEN ; pold='CH'  ; pnew='CHO'
        ELSE                                  ; pold='C'   ; pnew='CO'
        ENDIF
        CALL swap(tgroup(ip),pold,tempkg,pnew)
        tgroup(ip)=tempkg

        nc=INDEX(tgroup(i),' ') ; tgroup(i)(nc:nc)='.'
 
        CALL fragm(tbond,tgroup,p1chem,p2chem)
 
        yield=rodec/REAL(npath)
        brtio=brch*rab  
        brtio=brch  ! comparison old scheme

        IF (INDEX(p1chem,'.')/=0) THEN
          CALL chknadd2p(chem,p1chem,yield,brch,np,s,p,nref,ref)
          CALL stdchm(p2chem)
          CALL add1tonp(progname,chem,np)
          s(np)=yield ; brtio=brch*s(np)
          CALL bratio(p2chem,brtio,p(np),nref,ref)
        ELSE
          CALL stdchm(p1chem)
          CALL add1tonp(progname,chem,np)
          s(np)=yield ; brtio=brch*s(np)
          CALL bratio(p1chem,brtio,p(np),nref,ref)
          CALL chknadd2p(chem,p2chem,yield,brch,np,s,p,nref,ref)
        ENDIF             
        tgroup(:)=group(:)  ;  tbond(:,:)=bond(:,:)
        
      ENDIF
    ENDDO  
  ENDIF

  r(1)=pname  ;  r(2)='NO '  ! reactant
  arrh(1)=kno(1)*rab ; arrh(2:3)=kno(2:3)
  CALL rxwrit(mecu,r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe,ref)

! RO2+NO => RONO2 reaction
! ------------------------
  IF (rad/=0.) THEN

    CALL rxinit(r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe)
    nref=0  ;  ref(:)=' ' ;  CALL addref(progname,'MJ19KMV000',nref,ref,chem)
    np=0
    
! change (OO.) to (ONO2), rebuild and rename:
    pold='(OO.)' ;  pnew='(ONO2)'
    CALL swap(group(ip),pold,tgroup(ip),pnew)
    CALL rebond(tbond,tgroup,tempfo,tempring)
    CALL stdchm(tempfo)
    CALL add1tonp(progname,chem,np)   
    s(np)=1. ;  brtio=rad*brch
    CALL bratio(tempfo,brtio,p(1),nref,ref)

    DO i=1,ngr
      IF (INDEX(tgroup(i),'Cd')/=0) THEN
        IF (INDEX(tgroup(i),'(ONO2)')/=0) THEN
          PRINT*, "Cd(ONO2) found in: ",TRIM(tempfo)
          PRINT*, "parent is: ",TRIM(chem)
          STOP "in ro2no"
        ENDIF
      ENDIF
    ENDDO

    r(1)=pname ; r(2)='NO '    ! reactant
    arrh(1)=kno(1)*rad ; arrh(2:3)=kno(2:3)
    CALL rxwrit(mecu,r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe,ref)

  ENDIF
END SUBROUTINE ro2no

!=======================================================================
! PURPOSE: perform RO2+NO2 -> RO + NO3 reactions.
! Only for benzyl peroxy, see Jenkin et al., 2019
! EXCEPTION: if flag rx_ro2_no2 is activated,
!          perform reactions RO2+NO2 -> ROONO2
! Reverse reactions RO2+NO2 <- ROONO2 are called from panchem.f90
!=======================================================================
SUBROUTINE ro2no2(chem,pname,bond,group,brch,ip)
  USE keyparameter, ONLY: mxlco,mxpd,mecu
  USE keyflag, ONLY: rx_ro2_no2
  USE references, ONLY:mxlcod
  USE reactool, ONLY: swap,rebond
  USE rxwrttool, ONLY:rxwrit, rxinit
  USE toolbox, ONLY: stoperr,add1tonp,addref
  USE tempo, ONLY: chknadd2p

  CHARACTER(LEN=*),INTENT(IN) :: chem  ! formula of the VOC considered
  CHARACTER(LEN=*),INTENT(IN) :: pname ! short name of chem (parent compound) 
  CHARACTER(LEN=*),INTENT(IN) :: group(:) ! group (without ring joining char)
  INTEGER,INTENT(IN) :: bond(:,:)     ! bond matrix
  REAL,   INTENT(IN) :: brch          ! max yield of the input species
  INTEGER,INTENT(IN) :: ip            ! node bearing the peroxy (OO.) group

  CHARACTER(LEN=LEN(group(1))) :: tgroup(SIZE(group))
  CHARACTER(LEN=LEN(group(1))) :: pold,pnew
  INTEGER         :: tbond(SIZE(bond,1),SIZE(bond,2))
  CHARACTER(LEN=LEN(chem)) :: tempfo

  INTEGER :: np,tempring
  REAL    :: brtio, yield

  CHARACTER(LEN=mxlco) :: p(mxpd), r(3)
  REAL    :: s(mxpd), arrh(3)
  INTEGER :: idreac, nlabel
  REAL    :: xlabel,folow(3),fotroe(4)

  INTEGER,PARAMETER     :: mxref=10
  INTEGER               :: nref
  CHARACTER(LEN=mxlcod) :: ref(mxref)
  CHARACTER(LEN=7)      :: progname='ro2no2'

! Set the reaction rate (Jenkin et al., ACP, 2019).  
  REAL,DIMENSION(3),PARAMETER :: kno2=(/1.125E-11*((1./300.)**(-1.105)), -1.105 , 0. /) ! usual arrh exp.

  tgroup(:)=group(:)  ;  tbond(:,:)=bond(:,:)

! DO RO2+NO2->ROONO2 FALLOFF reaction, if requested by flag
  IF(rx_ro2_no2) CALL ro2no2_fo(chem,pname,bond,group,brch,ip)

! Following reaction is only for benzyl peroxy
  IF (tgroup(ip)/='c(OO.)') RETURN

! RO2+NO2 => RO+NO3 reaction
! -------------------------
  CALL rxinit(r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe)
  nref=0  ;  ref(:)=' ' ;  CALL addref(progname,'MJ19KMR000',nref,ref,chem)
  np=0

! add NO3 as product
  CALL add1tonp(progname,chem,np)   
  s(np)=1. ;  p(np)='NO3 '

! change (OO.) to (O.), rebuild, check and rename
  pold='(OO.)'  ;  pnew='(O.)'
  CALL swap(group(ip),pold,tgroup(ip),pnew)
  CALL rebond(tbond,tgroup,tempfo,tempring)
  
  yield=1.  
  brtio=brch  ! comparison old scheme
  CALL chknadd2p(chem,tempfo,yield,brch,np,s,p,nref,ref)
  
  r(1)=pname  ;  r(2)='NO2 '  ! reactant
  arrh(:)=kno2(:)
  CALL rxwrit(mecu,r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe,ref)

END SUBROUTINE ro2no2

!=======================================================================
! PURPOSE: write RO2+NO2 -> R(OONO2) FALLOFF reaction per NASA-JPL-2006
!          (optional, activated with flag rx_ro2_no2)
!=======================================================================
SUBROUTINE ro2no2_fo(chem,pname,bond,group,brch,ip)
  USE keyparameter, ONLY: mxlco,mxpd,mecu
  USE references, ONLY:mxlcod
  USE reactool, ONLY: swap,rebond
  USE normchem, ONLY: stdchm
  USE rxwrttool, ONLY:rxwrit, rxinit
  USE dictstacktool, ONLY: bratio
  USE toolbox, ONLY: stoperr,add1tonp,addref
  USE tempo, ONLY: chknadd2p

  CHARACTER(LEN=*),INTENT(IN) :: chem  ! formula of the VOC considered
  CHARACTER(LEN=*),INTENT(IN) :: pname ! short name of chem (parent compound)
  CHARACTER(LEN=*),INTENT(IN) :: group(:) ! group (without ring joining char)
  INTEGER,INTENT(IN) :: bond(:,:)     ! bond matrix
  REAL,   INTENT(IN) :: brch          ! max yield of the input species
  INTEGER,INTENT(IN) :: ip            ! node bearing the peroxy (OO.) group

  CHARACTER(LEN=LEN(group(1))) :: tgroup(SIZE(group))
  CHARACTER(LEN=LEN(group(1))) :: pold,pnew
  INTEGER         :: tbond(SIZE(bond,1),SIZE(bond,2))
  CHARACTER(LEN=LEN(chem)) :: tempfo

  INTEGER :: np,tempring
  REAL    :: brtio

  CHARACTER(LEN=mxlco) :: p(mxpd), r(3)
  REAL    :: s(mxpd), arrh(3)
  INTEGER :: idreac, nlabel
  REAL    :: xlabel,folow(3),fotroe(4)
  INTEGER,PARAMETER     :: mxref=10
  INTEGER,PARAMETER     :: mxrx=5
  INTEGER               :: nref
  CHARACTER(LEN=mxlcod) :: ref(mxref)
  CHARACTER(LEN=9)      :: progname='ro2no2_fo'
  CHARACTER(LEN=mxlcod) :: com(mxrx,mxref)

  tgroup(:)=group(:)  ;  tbond(:,:)=bond(:,:) ; p(:) = " "

! Exclude PAN-like RO2 : RO(OO.)
! -------------------------
  IF (tgroup(ip).EQ.'CO(OO.)') RETURN

! RO2+NO2 => ROONO2 reaction
! rates based on Table 2-1, JPL 2006
! -------------------------
  nref=0  ;  ref(:)=' ' ;  np=0
  CALL addref(progname,'SS06KMR000',nref,ref,chem)

! change (OO.) to (OONO2), rebuild, check and rename
  pold='(OO.)'  ;  pnew='(OONO2)'
  CALL swap(group(ip),pold,tgroup(ip),pnew)
  CALL rebond(tbond,tgroup,tempfo,tempring)

  CALL stdchm(tempfo)
  brtio = 1.
  CALL bratio(tempfo,brtio,p(1),nref,com(1,:))

  CALL rxinit(r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe)
  r(1) = pname  ;  r(2)='NO2 ' ; r(3)='(+M)' ! reactants
  s(1) = 1.

  idreac = 3
  arrh(1)= 7.2e-12
  arrh(2)= -2.1
  arrh(3)= 0.
  folow(1)= 1.0E-30
  folow(2)= -4.8
  folow(3)= 0.
  fotroe(1) = 0.6

  CALL rxwrit(mecu,r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe,ref)

END SUBROUTINE ro2no2_fo

!=======================================================================
! PURPOSE: write RO2+NO2 <- R(OONO2) FALLOFF REVERSE reaction (optional)
!=======================================================================
SUBROUTINE ro2no2_fo_reverse(chem,pname,bond,group,brch,ip)
  USE keyparameter, ONLY: mxlco,mxpd,mecu
  USE references, ONLY:mxlcod
  USE reactool, ONLY: swap,rebond
  USE normchem, ONLY: stdchm
  USE rxwrttool, ONLY:rxwrit, rxinit
  USE dictstacktool, ONLY: bratio
  USE toolbox, ONLY: stoperr,add1tonp,addref
  USE tempo, ONLY: chknadd2p

  CHARACTER(LEN=*),INTENT(IN) :: chem  ! formula of the VOC considered
  CHARACTER(LEN=*),INTENT(IN) :: pname ! short name of chem (parent compound)
  CHARACTER(LEN=*),INTENT(IN) :: group(:) ! group (without ring joining char)
  INTEGER,INTENT(IN) :: bond(:,:)     ! bond matrix
  REAL,   INTENT(IN) :: brch          ! max yield of the input species
  INTEGER,INTENT(IN) :: ip            ! node to receive the peroxy (OO.) group

  CHARACTER(LEN=LEN(group(1))) :: tgroup(SIZE(group))
  CHARACTER(LEN=LEN(group(1))) :: pold,pnew
  INTEGER         :: tbond(SIZE(bond,1),SIZE(bond,2))
  CHARACTER(LEN=LEN(chem)) :: tempfo

  INTEGER :: np,tempring
  REAL    :: brtio

  CHARACTER(LEN=mxlco) :: p(mxpd), r(3)
  REAL    :: s(mxpd), arrh(3)
  INTEGER :: idreac, nlabel
  REAL    :: xlabel,folow(3),fotroe(4)
  INTEGER,PARAMETER     :: mxref=10
  INTEGER,PARAMETER     :: mxrx=5
  INTEGER               :: nref
  CHARACTER(LEN=mxlcod) :: ref(mxref)
  CHARACTER(LEN=17)      :: progname='ro2no2_fo_reverse'
  CHARACTER(LEN=mxlcod) :: com(mxrx,mxref)

  tgroup(:)=group(:)  ;  tbond(:,:)=bond(:,:) ; p(:) = " "
  !PRINT*,progname," ",chem

! ROONO2 + M => RO2+NO2 reaction
! rates based on Tables 2-1 & 3-1, JPL 2006
! -------------------------
  CALL rxinit(r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe)
  nref=0  ;  ref(:)=' ' ;  np=0
  CALL addref(progname,'SS06KMR000',nref,ref,chem)

! change (OO.) to (OONO2), rebuild, check and rename
  pold='(OONO2)' ; pnew='(OO.)  '
  CALL swap(group(ip),pold,tgroup(ip),pnew)
  CALL rebond(tbond,tgroup,tempfo,tempring)
  CALL stdchm(tempfo)

  CALL rxinit(r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe)
  brtio = 1.
  r(1) = pname  ; r(2)='(+M)' ! reactants
  p(2)='NO2 '
  s(1:2) = 1.
  CALL bratio(tempfo,brtio,p(1),nref,com(1,:))

  idreac = 3
  arrh(1)= 7.58e16
  arrh(2)= -2.1
  arrh(3)= 11234
  folow(1)= 1.05E-2
  folow(2)= -4.8
  folow(3)= 11234
  fotroe(1) = 0.6

  CALL rxwrit(mecu,r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe,ref)
!  PRINT*,pname," + ",r(2)," -> ",p(1)," + ",p(2)

END SUBROUTINE ro2no2_fo_reverse

!=======================================================================
REAL FUNCTION nityield(bond,group,ngr,nring,Ncon,ip,tflg)
  USE keyparameter, ONLY: mxcp           
  USE keyflag, ONLY: TK, Mref           ! requested temperature, default pressure
  USE mapping, ONLY: gettrack
  IMPLICIT NONE

  CHARACTER(LEN=*),INTENT(IN) :: group(:) ! group (without ring joining char)
  INTEGER,INTENT(IN) :: ngr           ! # of nodes
  INTEGER,INTENT(IN) :: bond(:,:)     ! bond matrix
  INTEGER,INTENT(IN) :: ip            ! node bearing the peroxy (OO.) group
  INTEGER,INTENT(IN) :: Ncon          ! # of C, N, O atoms in peroxy "R"
  INTEGER,INTENT(IN) :: tflg          ! flag (0=prescribed T,P - 1=std T,P) 
  INTEGER,INTENT(IN) :: nring         ! # of rings in parent RO2
      
  REAL    :: r0,ri,Z,rapk,fa,fb,t,xm
  INTEGER :: i,j,k
  INTEGER :: track(mxcp,SIZE(bond,1))
  INTEGER :: trlen(mxcp)
  INTEGER :: ntr

  REAL,PARAMETER :: F=0.410
  REAL,PARAMETER :: alpha=2E-22
  REAL,PARAMETER :: beta=1.
  REAL,PARAMETER :: r_inf=0.430
  REAL,PARAMETER :: xm_inf=8.

  fa=1. ;  fb=1. ;  nityield=0.

  IF (tflg==1) THEN ; xm=2.46E+19 ; t=298.
  ELSE              ; xm=Mref       ; t=TK
  ENDIF

! COMPUTE NITRATE YIELD (before scaling with fa & fb)
! ---------------------------------------------------
  r0=alpha*EXP(beta*Ncon)*(t/300.)
  ri=r_inf*(t/300.)**(-xm_inf)
  Z=1.+(LOG10(r0*xm/ri))**2.
  Z=1./Z
  rapk=( (r0*xm)/(1.+(r0*xm/ri)) )*(F**Z)
  nityield=rapk/(1.+rapk)

! COMPUTE fa AND fb SCALING FACTORS
! --------------------------------- 
  IF  (group(ip)=='CH2(OO.)') fa=0.65
  IF  (group(ip)=='CH(OO.)')  fa=1.0
  IF  (group(ip)=='C(OO.)')   fa=1.0
  
! search for structures from aromatic oxidation
  IF (nring==2) THEN
    CALL gettrack(bond,ip,ngr,ntr,track,trlen)
 
    tracksloop: DO i=1,ntr
      IF (trlen(i) > 6) THEN
      IF ((group(track(i,3))=='-O-').AND.(group(track(i,4))=='-O-').AND. &
          (group(track(i,6))(1:2)=='Cd').AND.(group(track(i,7))(1:2)=='Cd')) THEN
        DO j=1,ngr
          IF ((INDEX(group(j),'(OH)')/=0).AND.(bond(track(i,2),j)/=0) ) THEN
            IF (group(ip)=='CH(OO.)') THEN 
              IF (group(track(i,2))(1:3)=='CH ') THEN ; fa=1.0
              ELSE                                    ; fa=0.43
              ENDIF
            ELSE
              IF (group(track(i,2))(1:3)=='CH ') THEN ; fa=0.13
              ELSE                                    ; fa=0.06
              ENDIF
            ENDIF
            fb=0.33
            nityield=nityield*fa*fb
            RETURN
          ENDIF
        ENDDO
      ENDIF
      ENDIF
    ENDDO tracksloop
  ENDIF

! check additional structures for fb factor

! alpha substituents
  IF (group(ip)=='C(OH)(OO.)') fb=fb*0.65
  IF (group(ip)=='CO(OO.)')    fb=0.
      
! beta substituents
  DO i=1,ngr
    IF (bond(ip,i)/=0) THEN
      IF (group(i)=='-O-') THEN
        fb=fb*0.65
        DO j=1,ngr
          IF ((bond(j,i)==3).AND.(group(j)(1:3)=='CO ')) fb=0.
        ENDDO
      ENDIF
 
      IF (INDEX(group(i),'(ONO2)')/=0) fb=fb*0.65
      IF (INDEX(group(i),'(OH)')  /=0) fb=fb*0.65
      IF (INDEX(group(i),'(OOH)') /=0) fb=fb*0.65
      
      IF (group(i)(1:3)=='CO ') fb=fb*0.3

! delta OH
      deltaloop : DO j=1,ngr
        IF ((bond(i,j)/=0).AND.(j/=ip)) THEN
          IF (group(j)=='-O-') fb=fb*0.65
          DO k=1,ngr
            IF ((bond(j,k)/=0).AND.(k/=i)) THEN
              IF (INDEX(group(k),'(OH)')/=0) THEN 
                fb=fb*0.8
                EXIT deltaloop
              ENDIF
            ENDIF
          ENDDO
        ENDIF
      ENDDO deltaloop   
      
    ENDIF
  ENDDO

! return the nitrate yield
  nityield=nityield*fa*fb

END FUNCTION nityield


!=======================================================================
! PURPOSE: perform and write the RO2+HO2 reactions. 
! 
! Following channels are considered for the reactions of RO2 with HO2, 
! as discussed in Jenkin et al., ACPD, 2019:
! RO2+HO2 -> ROOH+O2       (a) INDEX 1 
! RO2+HO2 -> ROH+O3        (b) INDEX 2 (only for RCOO2, not used in this routine)
! RO2+HO2 -> R(-H)O+H2O+O2 (c) INDEX 3
! RO2+HO2 -> RO+OH+O2      (d) INDEX 4
! RO2+HO2 -> R(-H)O+OH+HO2 (e) INDEX 5
!
! Four disctinct structures are considered to set branching ratios:
!   1: -CO-C(OO.)<            ,  2: R-O-C(OO.)<  
!   3: >C=C-C(OO.)-C(ONO2)<   ,  4: >C=C-C(OO.)-C(OH)< 
! If multiple "structure" for the same species are identified, an 
! average of the branching ratio of the different cases is used.      
! Default case assume channel (a) is 100% (i.e. ROOH formation).
!=======================================================================
SUBROUTINE ro2ho2(chem,pname,bond,group,ngr,brch,ip,Ncon)
  USE keyparameter, ONLY: mxlco,mxpd,mecu
  USE references, ONLY:mxlcod
  USE reactool, ONLY: swap,rebond
  USE dictstacktool, ONLY: bratio
  USE normchem, ONLY: stdchm
  USE rxwrttool, ONLY:rxwrit,rxinit
  USE toolbox, ONLY: stoperr,add1tonp,addrx,addref
  USE tempo, ONLY: chknadd2p

  CHARACTER(LEN=*),INTENT(IN) :: chem  ! formula of the VOC considered
  CHARACTER(LEN=*),INTENT(IN) :: pname ! short name of chem (parent compound) 
  CHARACTER(LEN=*),INTENT(IN) :: group(:) ! group (without ring joining char)
  INTEGER,INTENT(IN) :: ngr           ! # of nodes
  INTEGER,INTENT(IN) :: bond(:,:)     ! bond matrix
  REAL,INTENT(IN)    :: brch          ! max yield of the input species
  INTEGER,INTENT(IN) :: ip            ! node bearing the peroxy (OO.) group
  INTEGER,INTENT(IN) :: Ncon          ! # of C, N, O atoms in peroxy "R"

  CHARACTER(LEN=LEN(group(1))) :: tgroup(SIZE(group))
  CHARACTER(LEN=LEN(group(1))) :: pold,pnew
  INTEGER         :: tbond(SIZE(bond,1),SIZE(bond,2))
  CHARACTER(LEN=LEN(chem)) :: tempfo, p_onyl

  INTEGER :: i,j,np,tempring
  REAL    :: brtio, yield
  REAL    :: kho2(3)
  REAL    :: patratio(5)               ! path ratio for channel (a)-(e) 
  REAL    :: s1(5),s2(5),s3(5),s4(5) 
  INTEGER :: av(4)

  CHARACTER(LEN=mxlco) :: p(mxpd), r(3)
  REAL    :: s(mxpd), arrh(3)
  INTEGER :: idreac, nlabel
  REAL    :: xlabel,folow(3),fotroe(4)

  INTEGER,PARAMETER     :: mxref=10
  INTEGER               :: nref
  CHARACTER(LEN=mxlcod) :: ref(mxref)

  CHARACTER(LEN=7)      :: progname='ro2ho2'
  CHARACTER(LEN=70)     :: mesg

  tgroup(:)=group(:) ;  tbond(:,:)=bond(:,:)

! set the rate constant (see Jenkin, ACPD, 2019)
  kho2(1)=2.8E-13*(1-exp(-0.23*REAL(Ncon))) ; kho2(2)=0. ; kho2(3)=-1300.

  CALL rxinit(r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe)
  nref=0  ;  ref(:)=' ' ;  CALL addref(progname,'MJ19KMV000',nref,ref,chem)
  np=0
  
! set branching ratio according to the structure of parent RO2     
  patratio(1)=1.0                 ! default => ROOH (100% yield)
  patratio(2:5)=0.0  ; av(:)=0

! beta-oxoalkyl: -CO-C(OO.)<   (s1(:) group) 
  s1(:)=0.0 
  DO i=1,ngr
    IF ((tbond(i,ip) /= 0).AND.(tgroup(i)=='CO ')) THEN
      av(1)=av(1)+1
      IF (tgroup(ip)=='CH2(OO.)') THEN ; s1(1)=0.82 ; s1(4)=0.18   ! prim
      ELSE                             ; s1(1)=0.52 ; s1(4)=0.48   ! sec, ter
      ENDIF
      EXIT
    ENDIF      
  ENDDO

! alpha-alkoxyalkyl: R-O-C(OO.)<   (s2(:) group)  
  s2(:)=0.0
  DO i=1,ngr
    IF (tbond(i,ip)==3) THEN
      av(2)=av(2)+1
      IF (tgroup(ip)(1:2)=='CH') THEN ; s2(1)=0.54 ; s2(3)=0.26 ; s2(5)=0.20 ! prim, sec
      ELSE                            ; s2(1)=1.00                           ! ter
      ENDIF
      EXIT
    ENDIF      
  ENDDO

! beta-nitrooxy-allyl: >C=C-C(OO.)-C(ONO2)<  (s3(:) group)
  s3(:)=0.0
  betanitrooxyallylloop : DO i=1,ngr
    IF ((tbond(i,ip)==1).AND.(tgroup(i)(1:2)=='Cd')) THEN
      DO j=1,ngr
        IF ((tbond(j,ip)==1).AND.(INDEX(tgroup(j),'(ONO2)')/=0).AND.(i/=j)) THEN
          av(3)=av(3)+1 ; s3(1)=0.50 ; s3(4)=0.50
          EXIT betanitrooxyallylloop
        ENDIF
      ENDDO
    ENDIF
  ENDDO betanitrooxyallylloop

! beta-hydroxy-allyl: >C=C-C(OO.)-C(OH)<  (s4(:) group)
  s4(:)=0.0
  betahydroxyallylloop : DO i=1,ngr
    IF ((tbond(i,ip)==1).AND.(tgroup(i)(1:2)=='Cd')) THEN
      DO j=1,ngr
        IF ((tbond(j,ip)==1).AND.(INDEX(tgroup(j),'(OH)')/=0).AND.(i/=j)) THEN
          av(4)=av(4)+1 ; s4(1)=0.92 ; s4(4)=0.08
          EXIT betahydroxyallylloop
        ENDIF
      ENDDO
    ENDIF
  ENDDO betahydroxyallylloop
      
  IF (SUM(av(:))/=0) THEN
    patratio(1:5)=(s1(:)*av(1)+s2(:)*av(2)+s3(:)*av(3)+ s4(:)*av(4))/REAL(sum(av(:)))
  ENDIF

  IF (ABS(SUM(patratio(1:5))-1.)>0.001) THEN
    mesg="problem in ro2ho2 while seting the rate coefficient."
    CALL stoperr(progname,mesg,chem)
  ENDIF
       
! --- channel (a) - index 1: => ROOH+O2
  pold='(OO.)' ;  pnew='(OOH)'
  CALL swap(group(ip),pold,tgroup(ip),pnew)
  CALL rebond(tbond,tgroup,tempfo,tempring)
  CALL stdchm(tempfo)
  CALL add1tonp(progname,chem,np)   
  s(np)=patratio(1) ; brtio=brch*s(np)
  CALL bratio(tempfo,brtio,p(np),nref,ref)
  tgroup(:)=group(:)

! --- Channel (c) - index 3, => R(-H)O+H2O+O2
  IF (patratio(3) /= 0.) THEN
    pold='(OO.)' ; pnew=' '
    CALL swap(group(ip),pold,tgroup(ip),pnew)

    IF      (tgroup(ip)(1:3)=='CH2') THEN ; tgroup(ip)='CHO     '
    ELSE IF (tgroup(ip)(1:2)=='CH')  THEN ; tgroup(ip)(1:2)='CO'
    ENDIF
 
    CALL rebond(tbond,tgroup,p_onyl,tempring)
    CALL stdchm(p_onyl)
    CALL add1tonp(progname,chem,np)   
    s(np)=patratio(3) ; brtio=brch*s(np)
    CALL bratio(p_onyl,brtio,p(np),nref,ref)
    tgroup(:)=group(:)
  ENDIF

! --- Channel (d), index 4: => RO+OH+O2
  IF (patratio(4) /= 0.) THEN
    pold='(OO.)' ;  pnew='(O.)'
    CALL swap(group(ip),pold,tgroup(ip),pnew)
    CALL rebond(tbond,tgroup,tempfo,tempring)
    yield=patratio(4)
    CALL chknadd2p(chem,tempfo,yield,brch,np,s,p,nref,ref)
    CALL add1tonp(progname,chem,np)   
    s(np)=patratio(4) ; p(np)='HO '
    tgroup(:)=group(:)
  ENDIF
      
! --- Channel (e), index 5: => R(-H)O+OH+HO2
  IF (patratio(5) /= 0.) THEN
    CALL add1tonp(progname,chem,np)   
    s(np)=patratio(5) ; brtio=brch*s(np)
    IF (p_onyl==' ') THEN
      mesg="missing p_onyl "
      CALL stoperr(progname,mesg,chem)
    ENDIF
    CALL bratio(p_onyl,brtio,p(np),nref,ref)
    CALL add1tonp(progname,chem,np)   
    s(np)=patratio(5) ; p(np)='HO '
    CALL add1tonp(progname,chem,np)   
    s(np)=patratio(5) ; p(np)='HO2 '
    tgroup(:)=group(:)
  ENDIF

! Write the reaction
  r(1)=pname  
  r(2)='HO2 '
  arrh(:)=kho2(:)
  CALL rxwrit(mecu,r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe,ref)

END SUBROUTINE ro2ho2


!=======================================================================
! PURPOSE: perform the RO2+NO3 reaction and write the reaction. Only
! one reaction is considered : RO2+NO3 => RO + NO2 (+O2)
!=======================================================================
SUBROUTINE ro2no3(chem,pname,bond,group,brch,ip)
  USE keyparameter, ONLY: mxlco,mxpd,mecu,kno3u
  USE keyflag, ONLY: TK
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
  REAL    :: k298, kT
  INTEGER :: idreac, nlabel
  REAL    :: xlabel,folow(3),fotroe(4)

  INTEGER :: tempring, np
  REAL    :: brtio, yield

  INTEGER,PARAMETER     :: mxref=10
  INTEGER               :: nref
  CHARACTER(LEN=mxlcod) :: ref(mxref)

  CHARACTER(LEN=7)      :: progname='ro2no3'

! set the reaction rate (Jenkin et al., ACPD, 2019) 
  REAL,DIMENSION(3),PARAMETER :: kno3=(/8.9E-12, 0., 390. /)

  tgroup(:)=group(:) ;  tbond(:,:)=bond(:,:)

  CALL rxinit(r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe)
  nref=0  ;  ref(:)=' ' ;  CALL addref(progname,'MJ19KMV000',nref,ref,chem)
  np=0
  
! change (OO.) to (O.), rebuild, check and rename
  pold='(OO.)' ;  pnew='(O.)'
  CALL swap(group(ip),pold,tgroup(ip),pnew)
  CALL rebond(tbond,tgroup,tempfo,tempring)
  yield=1.;   brtio=brch*0.01

! check new formula, standardize and add products and coproducts to p(:)   
  CALL chknadd2p(chem,tempfo,yield,brch,np,s,p,nref,ref)

! add NO2 to the product list
  CALL add1tonp(progname,chem,np)   
  s(np)=1.  ;  p(np)='NO2 '

! write the reaction
  r(1)=pname  ;  r(2)='NO3 ' ; arrh(:)=kno3(:)
  CALL rxwrit(mecu,r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe,ref)

! write kNO3 to kjvocj
  k298 = arrh(1)*298.**arrh(2)*EXP(-arrh(3)/298.)
  kT = arrh(1)*TK**arrh(2)*EXP(-arrh(3)/TK)
  WRITE(kno3u,'(a,2x,2(1pe10.2),a)') r(1), k298, kT," SAR, SAR"

END SUBROUTINE ro2no3

!=======================================================================
! PURPOSE: perform and write the RO2+RO2 permutation reactions.
! Two cases are considered: 
!   - take the 9 classes of RO2 into account (if rx_ro2_multiclass==.TRUE.)
!   - take only CH3O2 into account (if rx_ro2_multiclass==.FALSE.)
!
! PROTOCOL FROM M. JENKIN ET AL., ACP, 2019
!=======================================================================
SUBROUTINE ro2ro2(chem,pname,bond,group,ngr,brch,ip,Ncon)
  USE keyparameter, ONLY: mxlco,mxpd,mecu,mxcopd
  USE keyflag, ONLY: rx_ro2_multiclass   ! If true, consider all classes of RO2
  USE references, ONLY:mxlcod
  USE reactool, ONLY: swap,rebond
  USE dictstacktool, ONLY: bratio
  USE normchem, ONLY: stdchm
  USE radchktool, ONLY: radchk
  USE rxwrttool, ONLY:rxwrit,rxinit
  USE toolbox, ONLY: stoperr,add1tonp,addrx,addref
  USE tempo, ONLY: chknadd2p
  USE dectool, ONLY: coono2_dec

  CHARACTER(LEN=*),INTENT(IN) :: chem  ! formula of the VOC considered
  CHARACTER(LEN=*),INTENT(IN) :: pname ! short name of chem (parent compound) 
  CHARACTER(LEN=*),INTENT(IN) :: group(:) ! group (without ring joining char)
  INTEGER,INTENT(IN) :: bond(:,:)     ! bond matrix
  INTEGER,INTENT(IN) :: ngr           ! # of nodes
  REAL,   INTENT(IN) :: brch          ! max yield of the input species
  INTEGER,INTENT(IN) :: ip            ! node bearing the peroxy (OO.) group
  INTEGER,INTENT(IN) :: Ncon          ! # of C, N, O atoms in peroxy "R"

  CHARACTER(LEN=LEN(group(1))) :: tgroup(SIZE(group))
  CHARACTER(LEN=LEN(group(1))) :: pold,pnew,tempkg
  INTEGER         :: tbond(SIZE(bond,1),SIZE(bond,2))
  CHARACTER(LEN=LEN(chem)) :: tempfo,p_onyl,p_ol,p_rad(2)

  INTEGER :: i,j,np,tempring,nsub,nsub1,nsub2
  REAL    :: brtio,wf
  REAL    :: k0ro2ro2,kro2ro2,fsub(9),kro2(9),s_ro2(9,3)
  REAL    :: f1(9),f2(9),fsub1,fsub2,ktest,fRO2

  INTEGER,PARAMETER    :: mxrpd=2   ! max # of products returned by radchk sub
  CHARACTER(LEN=mxlco) :: coprod(mxcopd),coprod2(mxcopd)
  CHARACTER(LEN=mxlco) :: coprod_onyl(mxcopd),coprod_ol(mxcopd)
  CHARACTER(LEN=LEN(chem)) :: rdckpd(mxrpd)
  CHARACTER(LEN=mxlco) :: rdckcopd(mxrpd,mxcopd)
  REAL    :: sc(mxrpd)
  INTEGER :: nip

  CHARACTER(LEN=mxlco) :: p(mxpd), r(3)
  REAL    :: s(mxpd), arrh(3)
  INTEGER :: idreac, nlabel
  REAL    :: xlabel,folow(3),fotroe(4)

  INTEGER,PARAMETER     :: mxref=10
  INTEGER               :: nref,tpnref
  CHARACTER(LEN=mxlcod) :: ref(mxref),tpref(mxref)
  INTEGER               :: ntpcom
  CHARACTER(LEN=mxlcod) :: tpcom(mxref)

  CHARACTER(LEN=7)      :: progname='ro2ro2'
  CHARACTER(LEN=70)     :: mesg

  CHARACTER(LEN=mxlco),DIMENSION(9), PARAMETER :: kwdclass= &
    (/'PERO1 ' , 'PERO2 ' , 'PERO3 ' , 'PERO4 ' , 'PERO5 ', &
      'MEPERO' , 'PERO7 ' , 'PERO8 ' , 'PERO9 '/)

  REAL, DIMENSION(9), PARAMETER :: kclass=  &  
    (/ 0.0021E-14, &  ! 1- unsubstituted ter RO2
       0.10E-14  , &  ! 2- iC3H7O2
       0.79E-14  , &  ! 3- ter RO2 with alpha or beta O or N
       6.90E-14  , &  ! 4- C2H5O2 & unsubstituted sec RO2
      10.00E-14  , &  ! 5- ter RO2 with alpha or beta O or N and allyl or benzyl
      35.00E-14  , &  ! 6- CH3O2
     110.00E-14  , &  ! 7- unsubstituted prim RO2 & sec RO2 with alpha or beta O or N
     530.00E-14  , &  ! 8- prim RO2 with alpha or beta O or N & sec RO2 with allyl and O, N
    1400.00E-14  /)   ! 9- acyl RO2

  REAL,DIMENSION(9,3),PARAMETER :: stoi_sp=RESHAPE( &  ! primary, secondary RO2
       (/ 0.8 , 0.2 , 0.0 ,  & ! PERO1- unsubstituted ter RO2
          0.6 , 0.2 , 0.2 ,  & ! PERO2- iC3H7O2
          0.8 , 0.2 , 0.0 ,  & ! PERO3- ter RO2 with alpha or beta O or N
          0.6 , 0.2 , 0.2 ,  & ! PERO4- C2H5O2 & unsubstituted sec RO2
          0.8 , 0.2 , 0.0 ,  & ! PERO5- subst. ter RO2 and allyl or benzyl
          0.6 , 0.2 , 0.2 ,  & ! MEPERO- CH3O2
          0.6 , 0.2 , 0.2 ,  & ! PERO7- unsubst. prim RO2 & subst. sec RO2
          0.6 , 0.2 , 0.2 ,  & ! PERO8- subst. prim RO2 & subst. sec RO2 with allyl
          0.8 , 0.2 , 0.0 /),& ! PERO9- acyl RO2
       SHAPE(stoi_sp), ORDER=(/2,1/) )  !- order (/2,1/) is to switch row & column from reshape

  REAL,DIMENSION(9,3),PARAMETER :: stoi_t=RESHAPE( &  ! tertiary RO2
      (/  1.0 , 0.0 , 0.0 ,  & ! PERO1- unsubstituted ter RO2
          0.8 , 0.0 , 0.2 ,  & ! PERO2- iC3H7O2
          1.0 , 0.0 , 0.0 ,  & ! PERO3- ter RO2 with alpha or beta O or N
          0.8 , 0.0 , 0.2 ,  & ! PERO4- C2H5O2 & unsubstituted sec RO2
          1.0 , 0.0 , 0.0 ,  & ! PERO5- subst. ter RO2 and allyl or benzyl
          0.8 , 0.0 , 0.2 ,  & ! MEPERO- CH3O2 
          0.8 , 0.0 , 0.2 ,  & ! PERO7- unsubst. prim RO2 & subst. sec RO2 
          0.8 , 0.0 , 0.2 ,  & ! PERO8- subst. prim RO2 & subst. sec RO2 with allyl 
          1.0 , 0.0 , 0.0 /),& ! PERO9- acyl RO2
       SHAPE(stoi_t), ORDER=(/2,1/) )  !- order (/2,1/) is to switch row & column from reshape

  ntpcom=0 ; tpcom(:)=' '
  tgroup(:)=group(:) ;  tbond(:,:)=bond(:,:)
  
  wf=0.02 ! weighting factor applied to bratio for the RO2+RO2 reactions

! /1/ find the various products linked to the species (alkoxy, alcohol, carbonyl)
  p_rad(1) =' '  ;  p_rad(2) =' ' ;  p_ol  =' ' ;  p_onyl=' '

  pold='(OO.)' ;  pnew='(O.)'                        !--- alkoxy product
  CALL swap(group(ip),pold,tgroup(ip),pnew)
  CALL rebond(tbond,tgroup,tempfo,tempring)
  CALL radchk(tempfo,rdckpd,rdckcopd,nip,sc,ntpcom,tpcom)
  p_rad(1)=rdckpd(1)
  CALL stdchm(p_rad(1))
  IF (nip==2) THEN
    p_rad(2)=rdckpd(2)
    CALL stdchm(p_rad(2))
  ENDIF
  coprod(:)=rdckcopd(1,:) ;  coprod2(:)=rdckcopd(2,:)

  pold='(OO.)' ;  pnew='(OH)'                        !--- alcohol product
  CALL swap(group(ip),pold,tgroup(ip),pnew)
  CALL rebond(tbond,tgroup,p_ol,tempring)
  CALL stdchm(p_ol)

  IF (group(ip)=='CH2(OO.) ') THEN                   !--- carbonyl product
    tgroup(ip)='CHO'
    CALL rebond(tbond,tgroup,p_onyl,tempring)
    CALL stdchm(p_onyl)
  ELSE IF (INDEX(group(ip),'CH')/=0) THEN
    pold='CH' ;  pnew='CO'
    CALL swap(group(ip),pold,tgroup(ip),pnew)
    pold='(OO.)' ;  pnew=' '
    tempkg=tgroup(ip)
    CALL swap(tempkg,pold,tgroup(ip),pnew)
    CALL rebond(tbond,tgroup,p_onyl,tempring)
    CALL stdchm(p_onyl)
  ENDIF

  coprod_onyl(:)=' ' ; coprod_ol(:)=' ' ; tpnref=0 ; tpref(:)=' '
  CALL coono2_dec(p_onyl,coprod_onyl,ngr,tpnref,tpref)
  CALL coono2_dec(p_ol,coprod_ol,ngr,tpnref,tpref)
        
! /2/ set RO2 type (primary, secondary or tertiary)
  nsub=0
  DO i=1,ngr
    IF (bond(ip,i)/=0) nsub=nsub+1
  ENDDO

  SELECT CASE (nsub)
    CASE(1) ; k0ro2ro2 = 10**(-11.7-(3.2*EXP(-0.55*(Ncon-0.52))))  ! primary RO2
    CASE(2) ; k0ro2ro2 = 10**(-12.9-(3.2*EXP(-0.64*(Ncon-2.30))))  ! secondary RO2
    CASE(3) ; k0ro2ro2 = 2.1E-17                                   ! tertiary RO2
    CASE DEFAULT
      mesg="problem in RO2+RO2 rate estimation - unindentified RO2 type"
      CALL stoperr(progname,mesg,chem)
  END SELECT

  kro2ro2=k0ro2ro2
            
! /3/ add substituent factor
  IF (group(ip)(1:1)/='c') THEN ! ignore scaling factor for aromatic RO2
    nsub2=0 ; nsub1=0 ; fsub(:)=0.  ; f1(:)=0. ; f2(:)=0.
    DO i=1,ngr
      IF (bond(ip,i)/=0) THEN
    
! 1st group: allyl or benzyl
        IF (group(i)(1:2)=='Cd')  THEN
          kro2ro2=kro2ro2*(4E-2/(k0ro2ro2**0.15))
          nsub1=nsub1+1 ; f1(nsub1)=(4E-2/(k0ro2ro2**0.15))
        ELSE IF (group(i)(1:1)=='c')   THEN
          nsub1=nsub1+1 ; f1(nsub1)=(5.8E-2/(k0ro2ro2**0.15))
          kro2ro2=kro2ro2*(5.8E-2/(k0ro2ro2**0.15))
        ENDIF
      
! 2nd group: oxygenated functional group
        IF ((INDEX(group(i),'(OH)')/=0).AND.(group(i)(1:2)/='CO')) THEN ! beta_hydroxy
          nsub2=nsub2+1  ;  fsub(nsub2)=(8E-5/(k0ro2ro2**0.4))
          f2(nsub2)=(8E-5/(k0ro2ro2**0.4))
        ENDIF
    
        IF (group(i)(1:3)=='-O-') THEN                                  ! alpha alkoxy
          nsub2=nsub2+1  ;  fsub(nsub2)=(7E-5/(k0ro2ro2**0.4))
          f2(nsub2)=(7E-5/(k0ro2ro2**0.4))
        ENDIF
    
        IF (group(i)(1:3)=='CO ') THEN                                  ! beta oxo
          nsub2=nsub2+1 ;  fsub(nsub2)=(1.6E-4/(k0ro2ro2**0.4))
          f2(nsub2)=(1.6E-4/(k0ro2ro2**0.4))
        ENDIF
    
        DO j=1,ngr
          IF ((bond(i,j)/=0).AND.(j/=ip)) THEN                          ! gamma oxo
            IF (group(j)(1:3)=='CO ') THEN
              nsub2=nsub2+1 ;  fsub(nsub2)= (5.3E-5/(k0ro2ro2**0.4))
              f2(nsub2)=(5.3E-5/(k0ro2ro2**0.4))
            ENDIF
          ENDIF
        ENDDO
      ENDIF
    ENDDO
  
! keep the most important activating factors only
    IF (nsub1/=0) THEN ; fsub1=MAXVAL(f1) ; ELSE ; fsub1=1.; ENDIF
    IF (nsub2/=0) THEN ; fsub2=MAXVAL(f2) ; ELSE ; fsub2=1.; ENDIF
    ktest=kro2ro2
    IF (MAXVAL(fsub)/=0) ktest=kro2ro2*MAXVAL(fsub)
    kro2ro2=k0ro2ro2*fsub1*fsub2
    IF (kro2ro2==0) THEN 
      mesg="problem in RO2+RO2 rate estimation - k_RO2=0"
      CALL stoperr(progname,mesg,chem)
    ENDIF
    
  ENDIF
  
! /4/ Compute the pseudo 1st order rate
  kro2(:)=0.  ;  s_ro2(:,:)=0.0

! multi-RO2-class parameterization
! --------------------------------
  IF (rx_ro2_multiclass) THEN   

    fRO2=1.
! MJ19KMV000 eqn (22)
    DO i=1,8
      kro2(i)=fRO2*2.*(kro2ro2*kclass(i))**0.5
    ENDDO

! scaling factor for cross reaction between (prim,sec) and (ter) RO2
! nsub refers to the species named in the reaction,
! pero class refers to the lumped class with which it is reacting
    fRO2=2.                   ! scaling factor set to 2.0
    IF (nsub/=3) THEN         ! non tertiary RO2 with tertiary class {1,3,5}
      kro2(1)=fRO2*kro2(1)    
      kro2(3)=fRO2*kro2(3)    
      kro2(5)=fRO2*kro2(5)    
    ELSE IF (nsub==3) THEN    ! tertiary RO2 with non tertiary class {2,4,6,7,8}
      kro2(2)=fRO2*kro2(2)
      kro2(4)=fRO2*kro2(4)
      kro2(6)=fRO2*kro2(6)
      kro2(7)=fRO2*kro2(7)
      kro2(8)=fRO2*kro2(8)
    ENDIF

! set branching ratio for radical & molecular (H shift) channel
    IF (group(ip)(1:1)=='c') nsub=3   
    IF (group(ip)(1:2)=='Cd') nsub=nsub+1   ! problem
    IF (p_onyl==' ') nsub=3                 ! no H available
    IF ((nsub==1).OR.(nsub==2)) THEN ; s_ro2(:,:)=stoi_sp(:,:)
    ELSE IF (nsub==3)           THEN ; s_ro2(:,:)=stoi_t(:,:)
    ENDIF

! /5/ loop over the 9 classes and write the reaction     
    DO j=1,9
      CALL rxinit(r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe)
      nref=0  ;  ref(:)=' ' ;  CALL addref(progname,'MJ19KMV000',nref,ref,chem)
      IF (tpnref/=0) THEN 
        DO i=1,tpnref ; CALL addref(progname,tpref(i),nref,ref,chem) ; ENDDO 
      ENDIF
      IF (ntpcom/=0) THEN    ! add reference linked to R(O.) 
        DO i=1,ntpcom ; CALL addref(progname,tpcom(i),nref,ref,chem) ; ENDDO
      ENDIF
      np=0

      IF (j<9) THEN      
! general case: MJ19KMV000, eqn (24): kro2 from eqn (22)
        arrh(1)=1E-13 ; arrh(2)=0. ; arrh(3)=-298*log(kro2(j)/1E-13)
      ELSE IF (j==9) THEN
! acyl peroxy: MJ19KMV000, eqn (23)
        arrh(1)=2E-12 ; arrh(2)=0. ; arrh(3)=-508
      ENDIF

! reactant
      r(1)=pname  ;  r(2)=kwdclass(j)

! radical channel
      CALL add1tonp(progname,chem,np)   
      s(np)=s_ro2(j,1)*sc(1)
      brtio=brch*wf 
      CALL bratio(p_rad(1),brtio,p(np),nref,ref)
      DO i=1,mxcopd                         ! add coproduct (if any)
        IF (coprod(i)(1:1)/=' ') THEN
          CALL add1tonp(progname,chem,np)   
           s(np)=s_ro2(j,1)*sc(1) ; p(np)=coprod(i)
        ENDIF
      ENDDO
      IF (nip==2) THEN
        CALL add1tonp(progname,chem,np)   
        s(np)=s_ro2(j,1)*sc(2)
        CALL bratio(p_rad(2),brtio,p(np),nref,ref)
        DO i=1,mxcopd                      ! add coproduct (if any)
          IF (coprod2(i)(1:1)/=' ') THEN
            CALL add1tonp(progname,chem,np)   
            s(np)=s_ro2(j,1)*sc(2)  ;  p(np)=coprod2(i)
          ENDIF
        ENDDO
      ENDIF

! "H given to counter" channel
      IF (s_ro2(j,2)>0.) THEN
        CALL add1tonp(progname,chem,np)   
        s(np)=s_ro2(j,2)
        brtio=brch*wf
        CALL bratio(p_onyl,brtio,p(np),nref,ref)
        DO i=1,mxcopd                         ! add coproduct (if any)
          IF (coprod_onyl(i)(1:1)/=' ') THEN
            CALL add1tonp(progname,chem,np)   
             s(np)=s_ro2(j,2) ; p(np)=coprod_onyl(i)
          ENDIF
        ENDDO
      ENDIF

! "H received from counter" channel
      IF (s_ro2(j,3)>0.) THEN
        CALL add1tonp(progname,chem,np)   
        s(np)=s_ro2(j,3)
        brtio=brch*0.5*wf
        CALL bratio(p_ol,brtio,p(np),nref,ref)
        DO i=1,mxcopd                         ! add coproduct (if any)
          IF (coprod_ol(i)(1:1)/=' ') THEN
            CALL add1tonp(progname,chem,np)   
             s(np)=s_ro2(j,3) ; p(np)=coprod_ol(i)
          ENDIF
        ENDDO
      ENDIF

      idreac=0  ;  nlabel=0
      CALL rxwrit(mecu,r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe,ref)
    ENDDO

! single RO2 class parameterization
! ---------------------------------
  ELSE  

    CALL rxinit(r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe)
    nref=0  ;  ref(:)=' ' ; CALL addref(progname,'MJ19KMV000',nref,ref,chem)
    IF (tpnref/=0) THEN 
      DO i=1,tpnref ; CALL addref(progname,tpref(i),nref,ref,chem) ; ENDDO 
    ENDIF
    IF (ntpcom/=0) THEN    ! add reference linked to R(O.) 
      DO i=1,ntpcom ; CALL addref(progname,tpcom(i),nref,ref,chem) ; ENDDO
    ENDIF
    np=0     

    fRO2=1.
    IF (nsub==3) fRO2=2.
! MJ19KMV000 eqn (22): 3.5E−13 = k298(CH3O2+CH3O2)
    kro2(1)=fRO2*2.*(kro2ro2*3.5E-13)**0.5
! MJ19KMV000 eqn (24): using kro2 from eqn (22)
    arrh(1)=1E-13 ; arrh(2)=0. ; arrh(3)=-298*log(kro2(1)/1E-13)
    
    IF (group(ip)(1:1)=='c') nsub=3   
    IF (group(ip)(1:2)=='Cd') nsub=nsub+1   
    IF (p_onyl==' ') nsub=3                 ! no H available
    s_ro2(1,1)=0.6   ; s_ro2(1,2)=0.2 ; s_ro2(1,3)=0.2 
    IF (nsub==3) THEN
      s_ro2(1,1)=0.8 ; s_ro2(1,2)=0.0 ; s_ro2(1,3)=0.2 
    ENDIF

! reactant
    r(1)=pname ; r(2)='PERO1'

! radical channel
    CALL add1tonp(progname,chem,np)   
    s(np)=s_ro2(1,1)*sc(1)
    brtio=brch*wf 
    CALL bratio(p_rad(1),brtio,p(np),nref,ref)
    DO i=1,mxcopd
      IF (coprod(i)(1:1)/=' ') THEN
        CALL add1tonp(progname,chem,np)   
        s(np)=s_ro2(1,1)*sc(1)  ;  p(np)=coprod(i)
      ENDIF
    ENDDO
    IF (nip==2) THEN
      CALL add1tonp(progname,chem,np)   
      CALL bratio(p_rad(2),brtio,p(np),nref,ref)
      s(np)=s_ro2(1,1)*sc(2)
      DO i=1,mxcopd
         IF (coprod2(i)(1:1)/=' ') THEN
           CALL add1tonp(progname,chem,np)   
           s(np)=s(2)  ;  p(np)=coprod2(i)
         ENDIF
      ENDDO
    ENDIF

! "H given to counter" channel
    IF (s_ro2(1,2)>0.) THEN
      CALL add1tonp(progname,chem,np)   
      s(np)=s_ro2(1,2)
      brtio=brch*wf
      CALL bratio(p_onyl,brtio,p(np),nref,ref)
      DO i=1,mxcopd                         ! add coproduct (if any)
        IF (coprod_onyl(i)(1:1)/=' ') THEN
          CALL add1tonp(progname,chem,np)   
           s(np)=s_ro2(j,2) ; p(np)=coprod_onyl(i)
        ENDIF
      ENDDO
    ENDIF

! "H received from counter" channel
    IF (s_ro2(1,3)>0.) THEN
      CALL add1tonp(progname,chem,np)   
      s(np)=s_ro2(1,3)
      brtio=brch*0.5*wf
      CALL bratio(p_ol,brtio,p(np),nref,ref)
      DO i=1,mxcopd                         ! add coproduct (if any)
        IF (coprod_ol(i)(1:1)/=' ') THEN
          CALL add1tonp(progname,chem,np)   
           s(np)=s_ro2(j,3) ; p(np)=coprod_ol(i)
        ENDIF
      ENDDO
    ENDIF

    idreac=0  ;  nlabel=0
    CALL rxwrit(mecu,r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe,ref)
 
  ENDIF
END SUBROUTINE ro2ro2

!=======================================================================
! PURPOSE: Perform and write the RO2+OH reaction. TO BE CHANGED,
! since based on the very first draft paper by Jenkin et al., 2019
! NEW ROOOH TO BE CONSIDERED HERE  - reaction should be changed !
!=======================================================================
SUBROUTINE ro2ho(chem,pname,bond,group,brch,ip)
  USE keyparameter, ONLY: mxlco,mxpd,mecu
  USE references, ONLY:mxlcod
  USE reactool, ONLY: swap,rebond
  USE rxwrttool, ONLY:rxwrit,rxinit
  USE toolbox, ONLY: stoperr,add1tonp,addrx,addref
  USE tempo, ONLY: chknadd2p
  USE atomtool, ONLY: getatoms
  USE dictstacktool, ONLY: bratio
  USE normchem, ONLY: stdchm
  USE namingtool, ONLY: codefg

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
  INTEGER :: ic,ih,in,io,ir,is,ifl,ibr,icl,Ncon

  INTEGER :: tempring,np
  REAL    :: brtio, yield

  INTEGER,PARAMETER     :: mxref=10
  INTEGER               :: nref
  CHARACTER(LEN=mxlcod) :: ref(mxref)
  CHARACTER(LEN=7)      :: progname='ro2ho'

! set the reaction rate (Jenkin et al., ACPD, 2019) 
  REAL,DIMENSION(3),PARAMETER :: kho=(/3.7E-11, 0., -350. /)

  tgroup(:)=group(:) ;  tbond(:,:)=bond(:,:)

  CALL rxinit(r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe)
  nref=0  ;  ref(:)=' ' ;  CALL addref(progname,'MJ19KMV000',nref,ref,chem)
  np=0

! Two pathways depending on Ncon
  CALL getatoms(chem,ic,ih,in,io,ir,is,ifl,ibr,icl)
  Ncon=ic+io-2+in

! === First pathway : RO2 + OH -> ROOOH
  pold ='(OO.)' ; pnew ='(OOOH)'
  CALL swap(group(ip),pold,tgroup(ip),pnew)
  
  CALL rebond(tbond,tgroup,tempfo,tempring)
  CALL stdchm(tempfo)
  CALL add1tonp(progname,chem,np)   
  s(np)=1.0 ; IF (Ncon==2) s(np) = 0.8
  brtio=s(np)*brch
  CALL bratio(tempfo,brtio,p(1),nref,ref)

! === Second pathway : RO2 + OH -> RO + HO2
  IF (Ncon==2) THEN
    tgroup(:)=group(:) ;  tbond(:,:)=bond(:,:)
     
! change (OO.) to (O.), rebuild, check and rename
    pold='(OO.)'  ;  pnew='(O.)'
    CALL swap(group(ip),pold,tgroup(ip),pnew)
    CALL rebond(tbond,tgroup,tempfo,tempring)
    yield=0.2;   brtio=s(1)*brch
    CALL chknadd2p(chem,tempfo,yield,brch,np,s,p,nref,ref)

! second product is HO2
    CALL add1tonp(progname,chem,np)   
    s(np)=0.2  ;  p(np)='HO2 '
  ENDIF

! reactant
  r(1)=pname  ;  r(2)='HO '  ;  arrh(:)=kho(:)
  CALL rxwrit(mecu,r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe,ref)

END SUBROUTINE ro2ho

END MODULE ro2tool
