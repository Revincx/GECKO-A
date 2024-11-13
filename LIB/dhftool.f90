MODULE dhftool
IMPLICIT NONE
CONTAINS
!SUBROUTINE dhf(rname,chem,bond,group,nring,brch,cut_off)
!SUBROUTINE dhf_thf(rname,chem,bond,group,nring,brch,cut_off)

!=======================================================================
! Purpose: Conversion of hydroxyketone (1,4HC) into dihydrofuran (DHF). 
! The reaction is represented as a first order process, without 
! considering the condensed phase processes (i.e. 1,4HC and DHF are 
! gas phase product only). For a more detailed representation, see
! the dhf_thf subroutine. 
!=======================================================================
SUBROUTINE dhf(rname,chem,bond,group,nring,brch,cut_off)
  USE keyparameter, ONLY: mxnode,mxcp,mxnr,mxpd,mecu,dhfu
  USE references, ONLY:mxlcod
  USE mapping, ONLY: abcde_map
  USE reactool, ONLY: swap,rebond
  USE ringtool, ONLY: findring
  USE normchem, ONLY: stdchm
  USE dictstacktool, ONLY: bratio,add_topvocstack
  USE rxwrttool, ONLY: rxwrit,rxinit
  USE tempflag, ONLY: iflost
  USE toolbox, ONLY: addrx,addref  
  IMPLICIT NONE

  CHARACTER(LEN=*),INTENT(IN) :: rname     ! short name of the input species
  CHARACTER(LEN=*),INTENT(IN) :: chem      ! formula of the input species
  CHARACTER(LEN=*),INTENT(IN) :: group(:)  ! group (without ring joining char)
  INTEGER,INTENT(IN)    :: bond(:,:)       ! bond matrix
  INTEGER,INTENT(IN)    :: nring           ! # of rings 
  REAL,INTENT(IN)       :: brch            ! max yield of the input species
  REAL,INTENT(IN)       :: cut_off         ! ratio below which a pathway is ignored 

! internal
  CHARACTER(LEN=LEN(chem)) :: pchem(mxnr)
  CHARACTER(LEN=LEN(group(1))) :: tgroup(SIZE(group)), pold, pnew
  INTEGER         :: ic,j,i,ia,igam,ngr,ibet
  INTEGER         :: itr, org
  INTEGER         :: rngflg,ring(SIZE(group))
  INTEGER         :: nr,ncd
  INTEGER         :: tbond(SIZE(bond,1),SIZE(bond,2))
  INTEGER         :: flag(mxnr)

! abcde_map variables
  INTEGER,PARAMETER :: mxdeep=5 ! max expected is gamma position (i.e. deep=4)
  INTEGER :: nabcde(mxdeep), tabcde(mxdeep,mxcp,mxnode)

  REAL            :: brtio
  REAL            :: ratio(mxnr)
  REAL            :: tarrhc(mxnr,3)

  CHARACTER(LEN=LEN(rname)) :: r(3), p(mxpd)
  REAL            :: s(mxpd), arrh(3)
  INTEGER         :: idreac, nlabel
  REAL            :: xlabel,folow(3),fotroe(4)

  INTEGER         :: tempiflost
  INTEGER         :: tempring

  INTEGER,PARAMETER     :: mxref=10
  INTEGER               :: nref
  CHARACTER(LEN=mxlcod) :: ref(mxref)
  CHARACTER(LEN=7),PARAMETER :: progname='dhf'

! Initialize
! ----------
  tempiflost=iflost
  ic=0  ; itr=0  ;  igam=0  ;  ia=0  ;  ibet=0  ;  ngr=0
  rngflg=0  ; org=0  ;  nr=0
  flag(:)=0 ; ratio(:)=0     ;  pchem(:)=' '
  tgroup(:)=group(:)       ;  tbond(:,:)=bond(:,:)  ;  tarrhc(:,:)=0.

  ngr=COUNT(group/=' ')

  IF (nring>0) RETURN

! Check the number of double bonds
  ncd=0
  DO i=1,ngr
    IF (tgroup(i)(1:2)=='Cd') ncd=ncd+1
  ENDDO
  IF (ncd>2) RETURN

! search for the right structure
! -------------------------------
  DO ic=1,ngr
    IF ((INDEX(group(ic)(1:3),'CO ')/=0).OR.(INDEX(group(ic)(1:4),'CHO ')/=0)) THEN

! search of an OH group in gamma
! check in the futur that both node does not belong to a cyle
      CALL abcde_map(bond,ic,ngr,nabcde,tabcde)

      igloop: DO itr=1,nabcde(4)
        igam=tabcde(4,itr,4)
        ia=tabcde(4,itr,2)
        ibet=tabcde(4,itr,3)
        IF ((INDEX(group(igam),'(OH)')/=0).AND.(group(igam)(1:2)/='CO')) THEN ! find one !

! Check structure 
          IF (group(igam)(1:1)/='C') CYCLE igloop
          IF (group(igam)(1:1)/='C')   CYCLE igloop ! no phenol OH
          IF (group(ia)(1:3)=='-O-')   CYCLE igloop ! no ether allowed on alpha
          IF (group(ibet)(1:3)=='-O-') CYCLE igloop ! no ether allowed on beta
          DO i=1, ngr        ! checks neighbors to -CO- and C-OH are not ethers
            IF (bond(i,ic)/=0) THEN
              IF (group(i)=='-O-') CYCLE igloop
            ENDIF
            IF (bond(i,igam)/=0) THEN
              IF (group(i)=='-O-') CYCLE igloop
            ENDIF
          ENDDO
          IF (group(ia)  /='CH2') CYCLE igloop ! alpha must be "CH2"
          IF (group(ibet)/='CH2') CYCLE igloop ! beta must be "CH2"

! If there's a ring with an ether function, the conversion can't occur
          DO i=1, ngr
            CALL findring(i,ia,ngr,bond,rngflg,ring)
            IF ((group(i)(1:3)=='-O-').AND.(ring(i)==1)) org=1
            IF ((rngflg==1).AND.(org==1)) CYCLE igloop
          ENDDO        

! find one - add a reaction
          CALL addrx(progname,chem,nr,flag)

! structure is change form RCO(1)-C2-C3-C4(OH)-R to 
! RCd(1)=Cd(2)-C3-C4[-O-(ngr+1)]-R

! group 1 (carbonyl to Cd)
          IF (group(ic)(1:2)=='CO')       THEN ; pold='CO'  ; pnew='Cd'
          ELSE IF (group(ic)(1:3)=='CHO') THEN ; pold='CHO' ; pnew='CdH'
          ENDIF
          CALL swap(group(ic),pold,tgroup(ic),pnew)

! group 2 (alpha to carbonyl)
          IF      (group(ia)(1:3)=='CH2') THEN ; pold='CH2' ; pnew='CdH'
          ELSE IF (group(ia)(1:2)=='CH')  THEN ; pold='CH'  ; pnew='Cd'
          ENDIF
          CALL swap(group(ia),pold,tgroup(ia),pnew)
                 
! group 4 (gamma to carbonyl bearing OH)
          pold='(OH)'  ;  pnew=''
          CALL swap(group(igam),pold,tgroup(igam),pnew)

! add group for -O- bond and create the bonds
          tgroup(ngr+1)='-O-'
          tbond(ngr+1,ic)=3    ;  tbond(ic,ngr+1)=3
          tbond(igam,ngr+1)=3  ;  tbond(ngr+1,igam)=3
          tbond(ic,ia)=2       ;  tbond(ia,ic)=2

! rebuild, check and rename:
          CALL rebond(tbond,tgroup,pchem(nr),tempring)
          CALL stdchm(pchem(nr))

! restore tgroup
          tgroup(:)=group(:)   ;   tbond(:,:)= bond(:,:)

        ENDIF
      ENDDO igloop
    ENDIF
  ENDDO              

! -------------------
! Write the reaction
!--------------------

! assign a rate constant
  DO i=1,nr
!    tarrhc(i,1)=1.7E-3
    tarrhc(i,1)=1E-3
    tarrhc(i,2)=0. ; tarrhc(i,3)=0.
    ratio(i)=1./nr  ! since all channel all equal
  ENDDO

! flag down for duplicate products
  IF (nr>1) THEN
    DO i=1,nr-1
      DO j= i+1,nr
        IF (pchem(i)==pchem(j)) THEN
           flag(j)=0
           tarrhc(i,1)=tarrhc(i,1) + tarrhc(j,1)
           ratio(i)=ratio(i)+ratio(j)
        ENDIF
      ENDDO
    ENDDO
  ENDIF

  iflost=0 ! cancel potential loss at the dhf level but reset below
  DO i=1,nr
    IF (flag(i)==0) CYCLE

! initialize reaction
    CALL rxinit(r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe)
    nref=0  ;  ref(:)=' ' ;  CALL addref(progname,'GASDHF',nref,ref,chem)
    arrh(:)=tarrhc(i,:)
    r(1)=rname ; r(2)=' ' ; r(3)=' '

! DHF and 1,4HC belong to the same generation. Add DHF on top
! of the VOC stack with appropriate generation number.
    brtio=brch * ratio(i)
    CALL add_topvocstack(pchem(i),brtio,p(1),nref,ref)
    CALL bratio(pchem(i),brtio,p(1),nref,ref)
    s(1)=1.

! write out: reaction set as thermal reaction (idreac=0, no labels)
    CALL rxwrit(mecu,r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe,ref)
    CALL rxwrit(dhfu,r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe)
  ENDDO
  iflost=tempiflost ! reset iflost to the correct value
                
END SUBROUTINE dhf
                 
! ======================================================================
! Purpose: Write the reaction in the condensed phase for the cyclization 
! of 1,4 hydroxycarbonyl (1,4HC) into cyclic hemiacetals (CHAs) that  
! next dehydrate to form dihydrofurans (DHF):  1,4HC -> CHA <-> DHF
! See papers by La et al., ACP, 1417–1431, 2016
! Note: Mass transfer between the gas phase and condensed phase are not
! treated here (see masstranstool)
! ======================================================================
SUBROUTINE dhf_thf(rname,chem,bond,group,nring,brch,cut_off)
  USE keyparameter, ONLY: mxnode,mxcp,mxnr,mxpd,mecu,dhfu
  USE references, ONLY:mxlcod
  USE keyflag, ONLY: g2pfg, g2wfg
  USE mapping, ONLY: abcde_map
  USE reactool, ONLY: swap,rebond
  USE ringtool, ONLY: findring
  USE normchem, ONLY: stdchm
  USE dictstackdb, ONLY: mxcha,ncha,chatab
  USE dictstacktool, ONLY: bratio,add_topvocstack
  USE rxwrttool, ONLY: rxinit,rxwrit
  USE tempflag, ONLY: iflost
  USE toolbox, ONLY: addrx,addref  
  IMPLICIT NONE

  CHARACTER(LEN=*),INTENT(IN) :: rname     ! short name of the input species
  CHARACTER(LEN=*),INTENT(IN) :: chem      ! formula of the input species
  CHARACTER(LEN=*),INTENT(IN) :: group(:)  ! group (without ring joining char)
  INTEGER,INTENT(IN)    :: bond(:,:)       ! bond matrix
  INTEGER,INTENT(IN)    :: nring           ! # of rings 
  REAL,INTENT(IN)       :: brch            ! max yield of the input species
  REAL,INTENT(IN)       :: cut_off         ! ratio below which a pathway is ignored 

! Internal
  CHARACTER(LEN=LEN(chem)) :: pchem(mxnr,2)
  CHARACTER(LEN=LEN(group(1))) :: tgroup(SIZE(group)), pold, pnew
  INTEGER         :: ic,j,i,ia,igam,ngr,ibet
  INTEGER         :: itr, org
  INTEGER         :: rngflg,ring(SIZE(group))
  INTEGER         :: nr,ncd
  INTEGER         :: tbond(SIZE(bond,1),SIZE(bond,2))
  INTEGER         :: flag(mxnr)
! abcde_map variables
  INTEGER,PARAMETER :: mxdeep=5 ! max expected is gamma position (i.e. deep=4)
  INTEGER :: nabcde(mxdeep), tabcde(mxdeep,mxcp,mxnode)

  REAL            :: brtio
  REAL            :: ratio(mxnr)
  INTEGER         :: tempiflost

  CHARACTER(LEN=LEN(rname)) :: r(3), p(mxpd)
  REAL            :: s(mxpd), arrh(3)
  INTEGER         :: idreac, nlabel
  REAL            :: xlabel,folow(3),fotroe(4)

  CHARACTER*1     :: fphase
  CHARACTER(LEN=LEN(chatab(1))) :: chaname,dhfname
  INTEGER         :: tempring

  INTEGER,PARAMETER     :: mxref=10
  INTEGER               :: nref
  CHARACTER(LEN=mxlcod) :: ref(mxref)
  CHARACTER(LEN=7),PARAMETER :: progname='dhf_thf'

! Initialize
! ----------
  tempiflost=iflost  ! saved to restore iflost
  ic=0  ; itr=0  ;  igam=0  ;  ia=0  ;  ibet=0  ;  ngr=0
  rngflg=0  ; org=0  ;  nr=0
  flag(:)=0 ; ratio(:)=0     ;  pchem(:,:)=' '
  tgroup(:)=group(:)  ;  tbond(:,:)=bond(:,:)  

  ngr=COUNT(group/=' ')

  IF (nring>0) RETURN

! Check the number of double bonds
  ncd=0
  DO i=1,ngr
    IF (tgroup(i)(1:2)=='Cd') ncd=ncd+1
  ENDDO
  IF (ncd>2) RETURN

! search for the right structure
! -------------------------------
  DO ic=1,ngr
    IF ((INDEX(group(ic)(1:3),'CO ')/=0).OR.(INDEX(group(ic)(1:4),'CHO ')/=0)) THEN

! search of an OH group in gamma
! check in the futur that both node does not belong to a cyle
      CALL abcde_map(bond,ic,ngr,nabcde,tabcde)

      igloop: DO itr=1,nabcde(4)
        igam=tabcde(4,itr,4)
        ia=tabcde(4,itr,2)
        ibet=tabcde(4,itr,3)
        IF ((INDEX(group(igam),'(OH)')/=0).AND.(group(igam)(1:2)/='CO')) THEN ! find one !

! Check structure 
          IF (group(igam)(1:1)/='C')   CYCLE igloop ! no phenol OH
          IF (group(ia)(1:3)=='-O-')   CYCLE igloop ! no ether allowed on alpha
          IF (group(ibet)(1:3)=='-O-') CYCLE igloop ! no ether allowed on beta
          DO i=1, ngr        ! checks neighbors to -CO- and C-OH are not ethers
            IF (bond(i,ic)/=0) THEN
              IF (group(i)=='-O-') CYCLE igloop
            ENDIF
            IF (bond(i,igam)/=0) THEN
              IF (group(i)=='-O-') CYCLE igloop
            ENDIF
          ENDDO
          IF (group(ia)  /='CH2') CYCLE igloop ! alpha must be "CH2"
          IF (group(ibet)/='CH2') CYCLE igloop ! beta must be "CH2"

! If there's a ring with an ether function, the conversion can't occur
          DO i=1, ngr
            CALL findring(i,ia,ngr,bond,rngflg,ring)
            IF ((group(i)(1:3)=='-O-').AND.(ring(i)==1)) org=1
            IF ((rngflg==1).AND.(org==1)) CYCLE igloop
          ENDDO        

! find one - add a reaction
          CALL addrx(progname,chem,nr,flag)

! --------------------
! make the DHF product
! --------------------
! change RCO(1)-C2-C3-C4(OH)-R  ==>  RCd(1)=Cd(2)-C3-C4[-O-(ngr+1)]-R

! group 1 (carbonyl to Cd)
          IF      (group(ic)(1:2)=='CO')  THEN ; pold='CO'  ; pnew='Cd'
          ELSE IF (group(ic)(1:3)=='CHO') THEN ; pold='CHO' ; pnew='CdH'
          ENDIF
          CALL swap(group(ic),pold,tgroup(ic),pnew)

! group 2 (alpha to carbonyl)
          IF      (group(ia)(1:3)=='CH2') THEN ; pold='CH2' ; pnew='CdH'
          ELSE IF (group(ia)(1:2)=='CH')  THEN ; pold='CH'  ; pnew='Cd'
          ENDIF
          CALL swap(group(ia),pold,tgroup(ia),pnew)
                 
! group 4 (gamma to carbonyl bearing OH)
          pold='(OH)'  ;  pnew=''
          CALL swap(group(igam),pold,tgroup(igam),pnew)

! add group for -O- bond and create the bonds
          tgroup(ngr+1)='-O-'
          tbond(ngr+1,ic)=3    ;  tbond(ic,ngr+1)=3
          tbond(igam,ngr+1)=3  ;  tbond(ngr+1,igam)=3
          tbond(ic,ia)=2       ;  tbond(ia,ic)=2

! rebuild, check and rename:
          CALL rebond(tbond,tgroup,pchem(nr,1),tempring)
          CALL stdchm(pchem(nr,1))

! restore tgroup
          tgroup(:)=group(:)  ;  tbond(:,:)= bond(:,:)

! --------------------
! make the THF product
! --------------------
! change form RCO(1)-C2-C3-C4(OH)-R  ==> RCHOH(1)-C2-C3-C4[-O-(ngr+1)]-R

! group 1 (carbonyl to Cd)
          IF     (group(ic)(1:2)=='CO') THEN ; pold='CO'  ; pnew='C(OH)'
          ELSE IF(group(ic)(1:3)=='CHO')THEN ; pold='CHO' ; pnew='CH(OH)'
          ENDIF
          CALL swap(group(ic),pold,tgroup(ic),pnew)
                 
! group 4 (gamma to carbonyl bearing OH)
          pold='(OH)'  ;  pnew=''
          CALL swap(group(igam),pold,tgroup(igam),pnew)

! add group for -O- bond and create the bonds
          tgroup(ngr+1)='-O-'
          tbond(ngr+1,ic)=3    ;  tbond(ic,ngr+1)=3
          tbond(igam,ngr+1)=3  ;  tbond(ngr+1,igam)=3

! rebuild, check and rename:
          CALL rebond(tbond,tgroup,pchem(nr,2),tempring)
          CALL stdchm(pchem(nr,2))

! restore tgroup
          tgroup(:)=group(:)   ;   tbond(:,:)= bond(:,:)

        ENDIF
      ENDDO igloop
    ENDIF
  ENDDO              

! -------------------
! Write the reaction
!--------------------

! assign a rate constant to make THF
  DO i=1,nr
    ratio(i)=1./nr  ! since all channel all equal
  ENDDO

! flag down for duplicate products
  IF (nr>1) THEN
    DO i=1,nr-1
      DO j= i+1,nr
        IF (pchem(i,1)==pchem(j,1)) THEN
           flag(j)=0
           ratio(i)=ratio(i)+ratio(j)
        ENDIF
      ENDDO
    ENDDO
  ENDIF

  iflost=0 ! cancel potential C loss here but reset below
  DO i=1,nr
    IF (flag(i)==0) CYCLE

! === 1-4 HC ---> CHA (1,4 hydroxy ketone to hydroxy tetrahydrofurans
    CALL rxinit(r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe)
    nref=0  ;  ref(:)=' ' ;  CALL addref(progname,'CP_HC2CHA',nref,ref,chem)
    arrh(1)=1E+7 ! arbitrary, set to a high value
    arrh(2)=0. ; arrh(3)=0.
    r(1)=rname

! THF, DHF and 1,4HC belong to the same generation.Add DHF on top
! of the VOC stack with appropriate generation number.
    brtio=brch*ratio(i)
    CALL add_topvocstack(pchem(i,2),brtio,p(1),nref,ref)
    s(1)=1.

! add the tetrahydrofuran (CHA) to the cha table (might be considered non-volatile)
    ncha=ncha+1
    IF (ncha > mxcha) THEN
      PRINT*, "max cha reached in dhf_thf" 
      STOP "in dhf_thf" 
    ENDIF
    chatab(ncha)=p(1)  ;  chaname=p(1)

! write out: reaction set as thermal reaction (idreac=0 and does 
! not require labels)
    IF (g2pfg) THEN
      fphase='A'
      CALL rxwrit(mecu,r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe,ref,phase=fphase)
      CALL rxwrit(dhfu,r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe,phase=fphase)
    ENDIF

    IF (g2wfg) THEN
      fphase='W'
      CALL rxwrit(mecu,r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe,ref,phase=fphase)
      CALL rxwrit(dhfu,r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe,phase=fphase)
    ENDIF

! === CHA ---> DHF (dihydration to dihydrofurans)
    CALL rxinit(r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe)
    nref=0  ;  ref(:)=' ' ;  CALL addref(progname,'CP_CHA2DHF',nref,ref,chem)
    arrh(1)=1E-3  ! La et al., 2016 (figure 3) 
!    arrh(1)=1E+7 ! arbitrary fast rate (full conversion 14HC -> DHF)
    arrh(2)=0. ; arrh(3)=0.
    r(1)=chaname

! THF, DHF and 1,4HC belong to the same generation. Add DHF on top
! of the VOC stack with appropriate generation number.
    brtio=brch*ratio(i)
    CALL add_topvocstack(pchem(i,1),brtio,p(1),nref,ref)
    CALL bratio(pchem(i,1),brtio,p(1),nref,ref)
    dhfname=p(1)
    s(1)=1.

! write out: reaction set as thermal reaction (idreac=0 and does not require labels)
    IF (g2pfg) THEN
      fphase='A'
      CALL rxwrit(mecu,r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe,ref,phase=fphase)
      CALL rxwrit(dhfu,r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe,phase=fphase)
    ENDIF

    IF (g2wfg) THEN
      fphase='W'
      CALL rxwrit(mecu,r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe,ref,phase=fphase)
      CALL rxwrit(dhfu,r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe,phase=fphase)
    ENDIF

! === DHF ---> CHA (hydration of dihydrofurans)
    CALL rxinit(r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe)
!    arrh(1)=1.E-7  ! ignore hydratation of DHF (full conversion 14HC -> DHF)
    nref=0  ;  ref(:)=' ' ;  CALL addref(progname,'CP_DHF2CHA',nref,ref,chem)
    arrh(1)=0.15  ! La et al., 2016 (figure 3) 
    arrh(2)=0. ; arrh(3)=0.
    r(1)=dhfname
    s(1)=1. ;  p(1)=chaname

! Note: the DHF -> CHA reaction not considered in La et al. (2016) on the wall 
    IF (g2pfg) THEN
      fphase='A'
      CALL rxwrit(mecu,r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe,ref,phase=fphase)
      CALL rxwrit(dhfu,r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe,phase=fphase)
    ENDIF
    
  ENDDO
  iflost=tempiflost ! reset iflost to the correct value                
END SUBROUTINE dhf_thf
!=======================================================================
!
!=======================================================================
! returns true if the bond i1-i2 is in a dhf cycle
LOGICAL FUNCTION in_dhf(group,bond,i1,i2)
  USE keyparameter, ONLY: mxcp
  USE mapping, ONLY: gettrack  
  IMPLICIT NONE
  
  CHARACTER(LEN=*),INTENT(IN) :: group(:)  ! group (without ring joining char)
  INTEGER,INTENT(IN)          :: bond(:,:) ! bond matrix
  INTEGER,INTENT(IN)          :: i1,i2     ! indexes of the nodes of the bond to test

  INTEGER :: track(mxcp,SIZE(bond,1))
  INTEGER :: trlen(mxcp)
  INTEGER :: ntr,ngr,nCd,nO,i,j
  
  in_dhf=.FALSE.
  
  ngr=COUNT(group/=' ')
  CALL gettrack(bond,i1,ngr,ntr,track,trlen)
  DO i=1,ntr
    IF ((track(i,5)==i2).AND.(bond(i1,i2)/=0)) THEN ! if bond(i1,i2) belongs to a 5 member ring
      nCd=0 ; nO=0
      DO j=1,5     ! count # of Cd and -O- in the 5-ring
        IF (group(track(i,j))(1:2)=='Cd')  nCd=nCd+1
        IF (group(track(i,j))(1:3)=='-O-') nO=nO+1
      ENDDO
      IF ((nCd==2).AND.(nO==1)) THEN
        in_dhf=.TRUE.
        RETURN
      ENDIF
    ENDIF
  ENDDO
  
END FUNCTION in_dhf
                   
END MODULE dhftool
