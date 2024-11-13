MODULE stuff4hvdiss2
IMPLICIT NONE
CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! OLD SUBROUTINE USED IN HVDISS2 - should be removed after updating
!!! the protocol for photolysis (BA, April 2019)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!SUBROUTINE setchrom2(rdct,bond,group,chromtab)
!SUBROUTINE openr(inbond,ingroup,innring,chem1,chem2,coprod)
!SUBROUTINE oprad(ingold,gnew,tcop)
!SUBROUTINE xcrieg(xcri,brch,s,p,cut_off)
!****************************************************************
!   MASTER MECHANISM V.3.1 ROUTINE NAME - SETCHROM              *
!   sets the chromophore table for a compound with a known      *
!   bond matrix and group list                                  *
!****************************************************************
SUBROUTINE setchrom2(chem,bond,group,chromtab)
  IMPLICIT NONE

  CHARACTER(LEN=*),INTENT(IN) :: chem
  CHARACTER(LEN=*),INTENT(IN) :: group(:)
  INTEGER,INTENT(IN)         :: bond(:,:)
  CHARACTER*1, INTENT(OUT)   :: chromtab(:,:)

! internal:
  INTEGER         :: i,j,k
  CHARACTER(LEN=LEN(group(1))) :: tgroup(SIZE(group))

  DO i=1,SIZE(group)
    tgroup(i)=group(i)
    DO j=1,4
      chromtab(i,j)=' '
    ENDDO
  ENDDO

! SET THE CHROMOPHORE TABLE  
  DO 10 i=1,SIZE(group)

! carbonyl chromophore
    IF (tgroup(i)(1:6)=='CO(OH)') THEN
       chromtab(i,1)='a'
       GOTO 10
    ENDIF
    IF (tgroup(i)(1:7)=='CO(OOH)') THEN
       chromtab(i,1)='g'
       GOTO 10
    ENDIF
    IF (tgroup(i)(1:8)=='CO(OOOH)') THEN
       !chromtab(i,1)='z'
       GOTO 10
    ENDIF
    IF (tgroup(i)(1:9)=='CO(OONO2)') THEN
       chromtab(i,1)='p'
       GOTO 10
    ENDIF
    IF (tgroup(i)(1:2)=='CO' .AND.INDEX(tgroup(i),'(ONO2)')/=0)THEN
       chromtab(i,1)='q'
       GOTO 10
    ENDIF
    IF (tgroup(i)(1:2)=='CO' .AND.INDEX(tgroup(i),'(NO2)')/=0)THEN
       chromtab(i,1)='v'
       GOTO 10
    ENDIF
    IF (tgroup(i)(1:4)=='CHO ') THEN
      DO j=1,SIZE(group)
        IF (bond(i,j)==3) GOTO 10
      ENDDO
       chromtab(i,1)='d'
       GOTO 10
    ENDIF
    IF (tgroup(i)(1:3)=='CO ') THEN
 
! No photolysis for ester function
      DO j=1,SIZE(group)
        IF (bond(i,j)==3) GOTO 10
      ENDDO
      chromtab(i,1)='k'
      GOTO 10
    ENDIF
    IF ((tgroup(i)(1:2)=='CO').AND.(tgroup(i)(3:3)/=' ')) THEN
       WRITE(6,'(a)') '--error in setchrom'
       WRITE(6,'(a)') 'carbonyl tgroup found but not identified'
       WRITE(6,'(a)') 'in the tgroup : '
       WRITE(6,'(a)') tgroup(i)
       WRITE(6,'(a)') 'in the molecule :'
       WRITE(6,'(a)') chem
       STOP
    ENDIF

! other chromophore
    j=0
20  j=j+1
  
    k=INDEX(tgroup(i),'(ONO2)')
    IF (k/=0) THEN
        chromtab(i,j)='n'
        tgroup(i)(k:k+5)='xxxxxx'
        GOTO 20
    ENDIF
  
    k=INDEX(tgroup(i),'(NO2)')
    IF (k/=0) THEN
        chromtab(i,j)='t'
        tgroup(i)(k:k+2)='xxx'
        GOTO 20
    ENDIF
  
    k=INDEX(tgroup(i),'(OOH)')
    IF (k/=0) THEN
        chromtab(i,j)='h'
        tgroup(i)(k:k+4)='xxxxx'
        GOTO 20
    ENDIF
  
    k=INDEX(tgroup(i),'(OH)')
    IF (k/=0) THEN
        chromtab(i,j)='o'
        tgroup(i)(k:k+3)='xxxx'
        GOTO 20
    ENDIF
    IF (INDEX(tgroup(i),'CO(OONO2)')/=0) THEN
       WRITE(6,'(a)') '--error in hvdisss2'
       WRITE(6,'(a)') 'peroxy nitrate found in the molecule :'
       WRITE(6,'(a)') chem
       WRITE(6,'(a)') 'chromophore not treated in hvdiss2'
       STOP
    ENDIF
  
10 CONTINUE

END SUBROUTINE setchrom2

!***********************************************************************
!  This subroutine defines the products of a ring-opening reaction.
!  It is assumed that the ring-opening is already done in the calling
!  program, by setting the corresponding element of BOND to zero
!  This subroutine expects to receive a bi-radical (e.g. R.-CO.)
!  with a bond matrix that makes sense (rejoin has already been called).
!  Reactions are the dominant pathway (only) from Calvert et al (2007)
!***********************************************************************
SUBROUTINE openr(inbond,ingroup,innring,chem1,chem2,coprod)
  USE keyparameter, ONLY: mxnode,mxlgr,waru
  USE references, ONLY:mxlcod
  USE cdtool, ONLY: alkcheck
  USE reactool
  USE ringtool, ONLY: findring
  USE normchem
  USE radchktool, ONLY: radchk
  IMPLICIT NONE

  INTEGER,INTENT(IN)          :: inbond(:,:)
  CHARACTER(LEN=*),INTENT(IN) :: ingroup(:)
  INTEGER,INTENT(IN)          :: innring
  CHARACTER(LEN=*),INTENT(OUT) :: chem1, chem2
  CHARACTER(LEN=*),INTENT(OUT) :: coprod(:)

! internal
  INTEGER                   :: tbnd1(SIZE(inbond,1),SIZE(inbond,2))
  INTEGER                   :: tbnd2(SIZE(inbond,1),SIZE(inbond,2))
  CHARACTER(LEN=LEN(chem1)) :: tchem
  CHARACTER(LEN=LEN(ingroup(1))) :: tgrp1(SIZE(ingroup))
  CHARACTER(LEN=LEN(ingroup(1))) :: tgrp2(SIZE(ingroup))
  CHARACTER(LEN=LEN(ingroup(1))) :: pold,pnew
  CHARACTER(LEN=LEN(coprod(1)))  :: tcop(SIZE(coprod))
  INTEGER         i,j,ii,jj,k,ia,ib,ic,ix,iy,iz
  INTEGER         ncr,nrad,nca,ind
  INTEGER         rngflg,ring(SIZE(inbond,1))
  CHARACTER(LEN=LEN(chem1)) rdckprod(mxnode)
  CHARACTER(LEN=LEN(coprod(1))) rdcktprod(mxnode,mxnode)
  INTEGER         nip
  REAL            sc(mxnode)
  INTEGER           :: nref
  CHARACTER(LEN=mxlcod) :: com(5)
  CHARACTER(LEN=mxlcod) :: acom                             ! comment code 

! save inputs      
  INTEGER :: bond(SIZE(inbond,1),SIZE(inbond,2))
  CHARACTER(LEN=LEN(ingroup(1))) :: group(SIZE(ingroup))
  INTEGER :: nring
  
  bond(:,:)=inbond(:,:)  ;  group(:)=ingroup(:)  ;  nring=innring
  nref=0 ; com(:)=' '
! Initialise
  nca=0 ; ncr=0 ; nrad=0
  chem1=' '  ; chem2=' ' ; tchem=' '  
  ia=0 ; ib=0 ; ix=0 ; iy=0 ; iz=0 ; ind= 0

  coprod(:)=' ' ;  tcop(:)=' ' ;  tgrp1(:)=' ' ;  tgrp2(:)=' '
  tbnd1(:,:)=0  ;  tbnd2(:,:)=0


! Count groups, identify radical groups (ia,iz)
! Define iz as furthest right CO., or as 
!           furthest right radical if no CO. exists
! Define ia as 'other' radical, count groups
  DO i=1,SIZE(inbond,1)
    IF(group(i)/=' ') nca=i
    IF(INDEX(group(i),'.')/=0) nrad=nrad+1
    IF(INDEX(group(i),'CO.')/=0) iz = i 
  ENDDO
  IF(iz==0)THEN
    DO i=nca,1,-1
      IF(INDEX(group(i),'.')/=0) iz = i
    ENDDO
  ENDIF
! if other radical end is substituted with (OOH), (OH) or (ONO2), 
! perform decomposition, treat resulting single radical, 
! and return to calling program
  DO i=nca,1,-1
    IF(INDEX(group(i),'.')/=0.AND.i/=iz) THEN
      IF(INDEX(group(i),'(OOH)')/=0 .OR. &
         INDEX(group(i),'(OH)')/=0 .OR.        &
         INDEX(group(i),'(ONO2)')/=0) THEN
        CALL oprad(group(i),tgrp1(i),coprod(i))
        group(i)=tgrp1(i)
        CALL rebond(bond,group,tchem,nring)
        CALL radchk(tchem,rdckprod,rdcktprod,nip,sc,nref,com)
        IF (nip==1) chem1 = rdckprod(1)
        IF (nip/=1) STOP 'openr.f'
        DO j=1,mxnode
          tcop(j) = rdcktprod(1,j)
        ENDDO
        DO j=1,nca
          ind=0
          DO k=j,nca
            IF(tcop(j)(1:1)/=' '.AND.coprod(k)(1:1)==' '.AND. ind==0)THEN
              coprod(k)=tcop(j)
              ind=1
            ENDIF
          ENDDO
        ENDDO
        RETURN 
      ENDIF

! otherwise, continue with cyclic bi-radical
      ia = i
    ENDIF
  ENDDO

  IF(nrad<2)THEN
    print*,'ERROR in subroutine openr!'
    print*,'bi-radical not found'
    WRITE(waru,*) 'openr',(group(i),i=1,nca) !STOP
  ELSE IF(nrad>2)THEN
    print*,'ERROR in subroutine openr!'
    print*,'too many radical centers'
    WRITE(waru,*) 'openr',(group(i),i=1,nca) !STOP
  ENDIF
     
! Find which carbons were previously on ring
  bond(ia,iz)=1
  bond(iz,ia)=1
  CALL findring(ia,iz,nca,bond,rngflg,ring)
  bond(ia,iz)=0
  bond(iz,ia)=0

! Find # of carbons previously on ring
  DO i=1,nca
    IF(ring(i)==1) ncr=ncr+1
  ENDDO

! identify groups in chain adjacent to radicals
! ia = '.'; ib = ia+1; ix = ib+1; iz = 'CO.'; iy = iz-1
  IF(iz>ia)THEN
    DO i=iz-1,ia+1,-1
      IF(ring(i)==1) ib = i
    ENDDO
    DO i=ia+1,iz-1
      IF(ring(i)==1) iy = i
    ENDDO
    DO i=iy-1,ib+1,-1
      IF(ring(i)==1) ic = i
    ENDDO
    DO i=ib+1,iy-1
      IF(ring(i)==1) ix = i
    ENDDO
  ELSE
    DO i=ia-1,iz+1,-1
      IF(ring(i)==1) iy = i
    ENDDO
    DO i=iz+1,ia-1
      IF(ring(i)==1) ib = i
    ENDDO
    DO i=ib-1,iy+1,-1
      IF(ring(i)==1) ic = i
    ENDDO
    DO i=iy+1,ib-1
      IF(ring(i)==1) ix = i
    ENDDO
  ENDIF

! ----------------------------------------
! a,b) bi-radical rearrangements (following photolysis of cyclic ketones)
! Pathways given here are the dominant among several channels.
! reference: Calvert et al, Chemistry of Alkanes (2007)
! ----------------------------------------
  IF(ncr>5) THEN
! C6+ =(a,b,c,x,y,z) ===========================

! a) acyl-to-aldehyde rearrangement (CH_. at ia, CH/CH2 at ib, CO. at iz)
    IF(INDEX(group(iz),'CO.')/=0)THEN
      IF(INDEX(group(ia),'CH')/=0 .AND. INDEX(group(ib),'CH')/=0)THEN
        !print*,'cC6: acyl-to-ketone rearrangement'
        pold='CO.'
        pnew='CHO'
        CALL swap(group(iz),pold,tgrp1(iz),pnew)
        group(iz)=tgrp1(iz)

        k=INDEX(group(ia),'.')
        group(ia)(3:k)=group(ia)(2:k-1)
        group(ia)(1:2)='Cd'
        IF(INDEX(group(ib),'CH2')==1)THEN
!         hydrogen is abstracted from ib
          group(ib)='CdH'
        ELSE
          group(ib)(1:2)='Cd'
          IF(INDEX(group(ib),'Cd(ONO2)')/=0)THEN
            pold='(ONO2)'
            pnew='(O.)'
            CALL swap(group(ib),pold,tgrp1(ib),pnew)
            group(ib)=tgrp1(ib)
            coprod(ib)='NO2   '
            CALL rebond(bond,group,tchem,nring)
          ENDIF
        ENDIF
        bond(ia,ib)=2
        bond(ib,ia)=2
        CALL rebond(bond,group,tchem,nring)
! jettison leaving groups if appropriate
        CALL alkcheck(tchem,coprod(ia),acom)
        IF (INDEX(tchem,'.')/=0)THEN
          CALL radchk(tchem,rdckprod,rdcktprod,nip,sc,nref,com)
          IF (nip==1) chem1 = rdckprod(1)
          IF (nip/=1) STOP 'openr.f'
          DO j=1,mxnode
            tcop(j) = rdcktprod(1,j)
          ENDDO
          CALL stdchm(chem1)
        ELSE
          chem1=tchem
        ENDIF
        DO i=1,nca
          IF(tcop(i)==' ') coprod(ib)=tcop(i)
        ENDDO

        GO TO 900
! b) CO elimination + new ring formation (for R.-C(O). if no H at both ia, ib)
      ELSE
        !print*,'cC6 : CO elimination'
        bond(ia,iy)=1
        bond(iy,ia)=1
        bond(iz,iy)=0
        bond(iy,iz)=0
        k=INDEX(group(ia),'.')
        group(ia)(3:k)=group(ia)(2:k-1)
        group(ia)(1:2)='C1'
        k=INDEX(group(iy),' ')
        group(iy)(3:k)=group(iy)(2:k-1)
        group(iy)(1:2)='C1'
        group(iz) = '      '
        coprod(iz) = 'CO    '
        CALL rebond(bond,group,chem1,nring)
        IF(INDEX(chem1,'.')/=0) THEN
          CALL radchk(chem1,rdckprod,rdcktprod,nip,sc,nref,com)
          IF (nip==1) tchem = rdckprod(1)
          IF (nip/=1) STOP 'openr.f'
          DO j=1,mxnode
            tcop(j) = rdcktprod(1,j)
          ENDDO
          DO j=1,nca
            IF(tcop(j)(1:1)/=' ') coprod(ib)=tcop(j)
          ENDDO
          chem1=tchem
          CALL stdchm(chem1)
        ENDIF
        GO TO 900
      ENDIF
    ENDIF

    print*,'ERROR in openr: C6 compound not found'
    DO i=1,nca
      print*,group(i)
    ENDDO
    !STOP
    WRITE(waru,*) 'openr',(group(i),i=1,nca) !STOP

  ELSE IF(ncr==5) THEN
! C5 =(a,b,x,y,z)====================================

! CO elimination and 2 x ethene formation (for R.-C(O).)
    IF(INDEX(group(iz),'CO.')/=0)THEN
      !print*,'cC5 : CO elimination'
      bond(iz,iy)=0
      bond(iy,iz)=0

      bond(iz,ia)=0
      bond(ia,iz)=0

      bond(ib,ix)=0
      bond(ix,ib)=0

! add radical dots (ia and iz are already radicals)
      DO j=1,3
        IF(j==1)i=ib
        IF(j==2)i=ix
        IF(j==3)i=iy
        k=INDEX(group(i),' ')
        group(i)(k:k) = '.'
      ENDDO

! first fragment: ethene (ia=ib) 
      IF(iz>ia)THEN ! ia-ib is already at start of array
        DO i=1,ix-1
          tgrp1(i)=group(i)
          DO j=1,ix-1
            tbnd1(i,j)=bond(i,j)
          ENDDO 
        ENDDO 
        ii=ia
        jj=ib
      ELSE ! iz < ia : shift indices of 2nd half of molecule
        DO i=ib,nca
          ii=i-ib+1
          tgrp1(ii)=group(i)
          DO j=ib,nca
            jj=j-ib+1
            tbnd1(ii,jj)=bond(i,j)
          ENDDO 
        ENDDO 
        ii=ia-ib+1
        jj=1
      ENDIF
! make ethene 
      DO j=1,2
        IF(j==1)i=ii
        IF(j==2)i=jj
        k=INDEX(tgrp1(i),'.')
        tgrp1(i)(3:k)=tgrp1(i)(2:k-1)
        tgrp1(i)(1:2)='Cd'
      ENDDO
      tbnd1(ii,jj)=2
      tbnd1(jj,ii)=2
      CALL rebond(tbnd1,tgrp1,chem1,nring)
! jettison leaving groups if appropriate
      CALL alkcheck(chem1,coprod(ia),acom)
      IF(INDEX(chem1,'.')/=0) THEN
        CALL radchk(chem1,rdckprod,rdcktprod,nip,sc,nref,com)
        IF (nip==1) tchem = rdckprod(1)
        IF (nip/=1) STOP 'openr.f'
        DO j=1,mxnode
          tcop(j) = rdcktprod(1,j)
        ENDDO
        DO j=1,nca
          IF(tcop(j)(1:1)/=' ') coprod(ib)=tcop(j)
        ENDDO
        chem1=tchem
      ENDIF

! second fragment: ethene (ix=iy) 
      IF(iz>ia)THEN  ! shift ix-iy to start of array 
        DO i=ix,iz-1
          ii=i-ix+1
          tgrp2(ii)=group(i)
          DO j=ix,iz-1
            jj=j-ix+1
            tbnd2(ii,jj)=bond(i,j)
          ENDDO 
        ENDDO 
        ii=iy-ix+1
        jj=1 
      ELSE ! iz < ia ; iz-iy is at start of array, shift 1 place to left
        DO i=2,ib-1
          ii = i-1
          tgrp2(ii)=group(i)
          DO j=2,ib-1
            jj = j-1
            tbnd2(ii,jj)=bond(i,j)
          ENDDO 
        ENDDO 
        jj=ix-iy+1
        ii=1
      ENDIF
! make ethene
      DO j=1,2
        IF(j==1)i=ii
        IF(j==2)i=jj
        k=INDEX(tgrp2(i),'.')
        tgrp2(i)(3:k)=tgrp2(i)(2:k-1)
        tgrp2(i)(1:2)='Cd'
      ENDDO
      tbnd2(ii,jj)=2
      tbnd2(jj,ii)=2
      CALL rebond(tbnd2,tgrp2,chem2,nring)
! jettison leaving groups if appropriate
      CALL alkcheck(chem2,coprod(ix),acom)
      IF(INDEX(chem2,'.')/=0) THEN
        CALL radchk(chem2,rdckprod,rdcktprod,nip,sc,nref,com)
        IF (nip==1) tchem = rdckprod(1)
        IF (nip/=1) STOP 'openr.f'
        DO j=1,mxnode
          tcop(j) = rdcktprod(1,j)
        ENDDO
        DO j=1,nca
          IF(tcop(j)(1:1)/=' ') coprod(iy)=tcop(j)
        ENDDO
        chem2=tchem
      ENDIF

! third fragment: CO coproduct
      coprod(iz) = 'CO    '

      GO TO 900
    ENDIF

    print*,'ERROR in openr: C5 compound not found'
    DO i=1,nca
      print*,group(i)
    ENDDO
    !STOP
    WRITE(waru,*) 'openr',(group(i),i=1,nca) !STOP

  ELSE IF(ncr==4) THEN
! C4 =(a,b,y,z)====================================

! R.-C(O). fragments to ketene + ethene (Calvert et al.)
    IF(INDEX(group(iz),'CO.')/=0)THEN
      bond(iz,ia)=0
      bond(ia,iz)=0

      bond(iy,ib)=0
      bond(ib,iy)=0

! add radical dots (ia and iz are already radicals)
      DO j=1,2
        IF(j==1)i=ib
        IF(j==2)i=iy
        k=INDEX(group(i),' ')
        group(i)(k:k) = '.'
      ENDDO

! ethene (ia=ib) 
      IF(iz>ia)THEN ! ia-ib is already at start of array
        DO i=1,iy-1
          tgrp1(i)=group(i)
          DO j=1,iy-1
            tbnd1(i,j)=bond(i,j)
          ENDDO 
        ENDDO 
        ii=ia
        jj=ib
      ELSE ! iz < ia : shift indices of 2nd half of molecule
        DO i=ib,nca
          ii=i-ib+1
          tgrp1(ii)=group(i)
          DO j=1,iy-1
            jj=i-ib+1
            tbnd1(ii,jj)=bond(i,j)
          ENDDO 
        ENDDO 
        ii=ia-ib+1
        jj=1
      ENDIF
! ii = radical center, jj = other center from ring
      DO j=1,2
        IF(j==1)i=ii
        IF(j==2)i=jj
        k=INDEX(tgrp1(i),'.')
        tgrp1(i)(3:k)=tgrp1(i)(2:k-1)
        tgrp1(i)(1:2)='Cd'
      ENDDO
      tbnd1(ii,jj)=2
      tbnd1(jj,ii)=2

      CALL rebond(tbnd1,tgrp1,chem1,nring)
      CALL alkcheck(chem1,coprod(ia),acom)
      IF(INDEX(chem1,'.')/=0) THEN
        CALL radchk(chem1,rdckprod,rdcktprod,nip,sc,nref,com)
        IF (nip==1) tchem = rdckprod(1)
        IF (nip/=1) STOP 'openr.f'
        DO j=1,mxnode
          tcop(j) = rdcktprod(1,j)
        ENDDO
         DO j=1,nca
          IF(tcop(j)(1:1)/=' ') coprod(ib)=tcop(j)
        ENDDO
        chem1=tchem
      ENDIF

! ketene (iy=iz=O) 
      IF(iz>ia)THEN  ! shift iy-iz to start of array 
        DO i=iy,nca
          ii=i-iy+1
          tgrp2(ii)=group(i)
          DO j=iy,iz
            jj=i-iy+1
            tbnd2(ii,jj)=bond(i,j)
          ENDDO 
        ENDDO 
        ii=nca-iy+1
        jj=1 
      ELSE ! iz < ia ; iz-iy is already at start of array
        DO i=1,nca-ib+1
          tgrp2(i)=group(i)
          DO j=1,iy-1
            tbnd2(i,j)=bond(i,j)
          ENDDO 
        ENDDO 
        jj=iy
        ii=iz
      ENDIF
      DO j=1,2
        IF(j==1)i=ii
        IF(j==2)i=jj
        k=INDEX(tgrp2(i),'.')
        tgrp2(i)(3:k)=tgrp2(i)(2:k-1)
        tgrp2(i)(1:2)='Cd'
      ENDDO
      tbnd2(ii,jj)=2
      tbnd2(jj,ii)=2

      CALL rebond(tbnd2,tgrp2,chem2,nring)
      CALL alkcheck(chem2,coprod(iy),acom)
      IF(INDEX(chem2,'.')/=0) THEN
        CALL radchk(chem2,rdckprod,rdcktprod,nip,sc,nref,com)
        IF (nip==1) tchem = rdckprod(1)
        IF (nip/=1) STOP 'openr.f'
        DO j=1,mxnode
          tcop(j) = rdcktprod(1,j)
        ENDDO
        DO j=1,nca
          IF(tcop(j)(1:1)/=' ') coprod(iz)=tcop(j)
        ENDDO
        chem2=tchem
      ENDIF

      GO TO 900
    ENDIF

    print*,'ERROR in openr: C4 compound not found'
    DO i=1,nca
      print*,group(i)
    ENDDO
    WRITE(waru,*) 'openr',(group(i),i=1,nca) !STOP
    !STOP

  ELSE IF(ncr==3) THEN
! C3 =(a,b,z)===================================

! CO elimination and ethene formation (for R.-C(O). or C(O).-R.)
    IF(INDEX(group(iz),'CO.')/=0)THEN
      bond(ia,ib)=2
      bond(ib,ia)=2
      bond(iz,ib)=0
      bond(ib,iz)=0
! add radical dot (ia and iz are already radicals)
      i=ib
      k=INDEX(group(i),' ')
      group(i)(k:k) = '.'

      DO j=1,2
        IF(j==1)i=ia
        IF(j==2)i=ib
        k=INDEX(group(i),'.')
        group(i)(3:k)=group(i)(2:k-1)
        group(i)(1:2)='Cd'
      ENDDO

      CALL rebond(bond,group,chem1,nring)
      CALL alkcheck(chem1,coprod(ia),acom)
      IF(INDEX(chem1,'.')/=0) THEN
        CALL radchk(chem1,rdckprod,rdcktprod,nip,sc,nref,com)
        IF (nip==1) tchem = rdckprod(1)
        IF (nip/=1) STOP 'openr.f'
        DO j=1,mxnode
          tcop(j) = rdcktprod(1,j)
        ENDDO
        DO j=1,nca
          IF(tcop(j)(1:1)/=' ') coprod(ib)=tcop(j)
        ENDDO
        chem1=tchem
      ENDIF
      group(iz) = '      '
      coprod(iz) = 'CO    '
      GO TO 900
    ENDIF

    print*,'ERROR in openr: C3 compound not found'
    DO i=1,nca
      print*,group(i)
    ENDDO
    WRITE(waru,*) 'openr',(group(i),i=1,nca) !STOP

  ENDIF

900 CONTINUE      
  IF(chem2(1:3)=='CO ')THEN
    chem2(1:3)='   '
    coprod(nca+1)(1:3)='CO '
  ENDIF

END SUBROUTINE openr

!=======================================================================
! Based on Multip, but treats only ONE group: the substituted end of a di-
! radical formed in ring-opening 
!=======================================================================
SUBROUTINE oprad(ingold,gnew,tcop)
  USE reactool
  IMPLICIT NONE

  CHARACTER(LEN=*),INTENT(IN)  :: ingold
  CHARACTER(LEN=*),INTENT(OUT) :: gnew
  CHARACTER(LEN=*),INTENT(OUT) :: tcop

! internal
  CHARACTER(LEN=LEN(ingold))  :: pold, pnew, gold
  
! initialize:
  gold=ingold  ;  gnew=' '  ;  tcop=' '

! ----------------------------------------
! OPTIONS FROM SUBROUTINE MULTIP
! ----------------------------------------
! radical on double-TBOND carbon TGROUP     : not likely in openr
! ----------------------------------------
! radical on alkyl carbon : criegge         : not likely in openr
! ----------------------------------------
! radical on alkyl carbon : di-radical      : not likely in openr
! ----------------------------------------
! radical on carbonyl-carbon TGROUP         : deal with later in openr
! ----------------------------------------
! radical on alkyl carbon : primary alkyl   : deal with later in openr
! ----------------------------------------
! radical on alkyl carbon : secondary alkyl : allow decomposition
! -----------------------------------------
! CH(OOH).  / CH(OH).  / CH(ONO2). -> CHO

  IF (INDEX(gold,'CH')/=0) THEN

     IF (INDEX(gold,'(OOH)')/=0) THEN
        gnew = 'CHO'
        tcop = 'HO   '
     ELSE IF (INDEX(gold,'(OH)')/=0) THEN
        gnew = 'CHO'
        tcop = 'HO2  '
     ELSE IF (INDEX(gold,'(ONO2)')/=0) THEN
        gnew = 'CHO'
        tcop = 'NO2  '
     ELSE IF ( (INDEX(gold,'(OOH)')==0) .AND. &
              (INDEX(gold,'(OH)')==0) .AND.        &
              (INDEX(gold,'(ONO2)')==0) ) THEN
        print*,"ERROR: subroutine oprad invoked for group"
        print*,gold
        STOP
     ENDIF         

! -----------------------------------------
! radical on alkyl carbon : tertiary alkyl  : allow decomposition
! -----------------------------------------
! C(OOH)(R).  / C(OH)(R).  / C(ONO2)(R). -> CO(R)

  ELSE IF ( (INDEX(gold,'Cd'    )==0) .AND. &
            (INDEX(gold,'CO'    )==0) .AND. &
            (INDEX(gold,'.(OO.)')==0) .AND. &
            (INDEX(gold,'..'    )==0) .AND. &
            (INDEX(gold,'CH'    )==0) ) THEN

     IF (INDEX(gold,'(OOH)')/=0) THEN
        pold = '.'
        pnew = ' '
        CALL swap(gold,pold,gnew,pnew)
        gold = gnew
        pold = 'C'
        pnew = 'CO'
        CALL swap(gold,pold,gnew,pnew)
        gold = gnew
        pold = '(OOH)'
        pnew = ' '
        CALL swap(gold,pold,gnew,pnew)
        gold = gnew
        tcop = 'HO   '
     ELSE IF (INDEX(gold,'(OH)')/=0) THEN
        pold = '.'
        pnew = ' '
        CALL swap(gold,pold,gnew,pnew)
        gold = gnew
        pold = 'C'
        pnew = 'CO'
        CALL swap(gold,pold,gnew,pnew)
        gold = gnew
        pold = '(OH)'
        pnew = ' '
        CALL swap(gold,pold,gnew,pnew)
        gold = gnew
        tcop = 'HO2  '
     ELSE IF (INDEX(gold,'(ONO2)')/=0) THEN
        pold = '.'
        pnew = ' '
        CALL swap(gold,pold,gnew,pnew)
        gold = gnew
        pold = 'C'
        pnew = 'CO'
        CALL swap(gold,pold,gnew,pnew)
        gold = gnew
        pold = '(ONO2)'      
        pnew = ' '
        CALL swap(gold,pold,gnew,pnew)
        gold = gnew
        tcop = 'NO2  '
     ELSE
        print*,"ERROR: subroutine oprad invoked for group"
        print*,gold
        STOP
     ENDIF
  ENDIF 
END SUBROUTINE oprad


!=======================================================================
! MASTER MECHANISM - ROUTINE NAME : xcrieg                           
!                                                                    
!                                                                    
! PURPOSE :                                                          
!    perform the reactions (decomposition, stabilization) of hot     
!    criegge produced by O3+alkene reactions. The routine returns the
!    list of products and associated stoi. coef. in tables p (short  
!    name) and s, respectively. The dictionary, stack and related    
!    tables are updated accordingly.                                 
!                                                                    
! INPUT:                                                             
!  - xcri        : formula of the hot criegge                        
!  - cut_off     : ratio below which a pathway is not considered     
!                                                                    
! OUTPUT:                                                            
!  - s(i)        : stoichiometric coef. associated with product i    
!  - p(i)        : product i                                         
!=======================================================================
SUBROUTINE xcrieg(xcri,brch,s,p,cut_off,nref,ref)
  USE keyparameter, ONLY: mxnode, mxlgr,waru,mxring
  USE references, ONLY:mxlcod
  USE atomtool
  USE stdgrbond
  USE reactool
  USE ringtool, ONLY: findring
  USE normchem
  USE dictstacktool, ONLY: bratio
  USE radchktool, ONLY: radchk
  USE hotool, ONLY: rabsoh
  USE fragmenttool, ONLY: fragm
  USE criegeetool, ONLY: add2p
  USE rjtool,     ONLY: rjgrm

  IMPLICIT NONE
     
  CHARACTER(LEN=*),INTENT(IN) :: xcri
  REAL,INTENT(IN)             :: brch
  REAL,INTENT(IN)             :: cut_off
  CHARACTER(LEN=*),INTENT(OUT):: p(:)
  REAL,INTENT(OUT)            :: s(:)
  INTEGER,INTENT(INOUT) :: nref          ! # of references added in the reference list
  CHARACTER(LEN=*),INTENT(INOUT):: ref(:)! list of references 

! internal
  INTEGER        :: l, mnp, mnr
  CHARACTER(LEN=LEN(xcri)) :: pchem(SIZE(p)),pchem1(SIZE(p)),pchem2(SIZE(p))
  CHARACTER(LEN=LEN(xcri)) :: prod,prod2,tempfo
  CHARACTER(LEN=mxlgr)     :: tgroup(mxnode),group(mxnode),tempgr,pold, pnew
  CHARACTER(LEN=LEN(p(1))) :: coprod(mxnode),coprod_del(mxnode)
  INTEGER         tbond(mxnode,mxnode),bond(mxnode,mxnode)
  INTEGER         dbflg,nring
  INTEGER         np,nc,nca,ig,j1,j2,i,j,nold,ngr
  INTEGER         posj1,posj2
  REAL            ftherm,fdec1,fmol,yld1,yld2
  REAL            brtio
  REAL            garrhc(3),c1,c2
  CHARACTER(LEN=LEN(p(1)))  :: copchem
  REAL              rdtcopchem 
  INTEGER         rngflg       ! 0 = 'no', 1 = 'yes'
  INTEGER         ring(mxnode)    ! =1 if node participates in current ring
  CHARACTER(LEN=LEN(xcri)) :: rdckprod(mxnode)
  CHARACTER(LEN=LEN(p(1)))  :: rdcktprod(mxnode,mxnode)
  INTEGER         nip
  REAL            sc(mxnode),zeyield,sciyield
  CHARACTER(LEN=LEN(xcri)) :: pchem_del
  INTEGER                  :: rjg(mxring,2)              ! ring-join group pairs
    
! initialize:
  np=0  ;  ig=0  ;  rdtcopchem=0.

  pchem(:)=' '  ;  p(:)=' '  ;  s(:)=0.
  mnp=SIZE(p) ; mnr=SIZE(p)

! calling function to define number of carbons in hot Criegee
      nc = INDEX(xcri,' ') - 1
      nca = cnum(xcri)
      ngr=0
      DO i=1,mxnode
        IF (group(i)/=' ') ngr=ngr+1
      ENDDO

! define functional groups and bond-matrix of hot Criegee and
! store the data into bond and group
      CALL grbond(xcri,group,bond,dbflg,nring)
      IF (nring>0)  CALL rjgrm(nring,group,rjg)  
      tgroup(:)=group(:)
      tbond(:,:)=bond(:,:)

! find hot_criegge groups:
      DO i=1,mxnode
         IF (INDEX(tgroup(i),'.(OO.)*')/=0) ig = i
      ENDDO
      IF (ig==0) THEN
         WRITE(6,'(a)') '--warning--(stop) in xcrieg'
         WRITE(6,'(a)') 'hot_criegge functional group not found in :'
         WRITE(6,'(a)') xcri
         WRITE(waru,*) 'xcrieg',xcri !STOP
      ENDIF

! -----------------------------------------
!  1st section: single carbon hot criegees:
! ------------------------------------------

      IF (nca==1) THEN

! branching ratio for CH2.(OO.)* are taken in Atkinson, 1999, J. Chem. Ref. 
! Data. The stabilized criegee radical is expected to produce a carboxylic 
! acid and the "cold criegee" is no longer treated. If criegee radical is 
! not specifically adressed below then the program stops. 
         IF (xcri(1:10)=='CH2.(OO.)*') THEN
           s(1) = 0.37
           pchem(1) = 'CHO(OH)'
           brtio = brch * s(1)
           CALL bratio(pchem(1),brtio,p(1),nref,ref)
           s(2) = 0.13
           p(2) = 'CO2  '
           s(3) = 0.13
           p(3) = 'H2   '
           s(4) = 0.50
           p(4) = 'CO   '
           s(5) = 0.38
           p(5) = 'H2O  '
           s(6) = 0.12
           p(6) = 'HO2  '
           s(7) = 0.12
           p(7) = 'HO   '
           RETURN
! escape routes for minor products of substituted ring-opening
! if criegee is "CH(OH).(OO.)" then assume NO conversion to acid, 100% yield.
         ELSE IF (xcri(1:13)=='CH(OH).(OO.)*') THEN
           s(1) = 1.0
           p(1) = 'CHO(OH)  '
           s(2) = 1.0
           p(2) = 'NO2  '
           RETURN
! if criegee is "CH(ONO2).(OO.)" then assume NO2 elimination, 100% yield.
         ELSE IF (xcri(1:15)=='CH(ONO2).(OO.)*') THEN
           pchem(1) = 'CHO(OO.) '
           s(1) = 1.0
           brtio = brch * s(1)
           CALL bratio(pchem(1),brtio,p(1),nref,ref)
           s(2) = 1.0
           p(2) = 'NO2  '
           RETURN
! if criegee is "CH(OOH).(OO.)" then assume OH elimination, 100% yield.
         ELSE IF (xcri(1:14)=='CH(OOH).(OO.)*') THEN
           pchem(1) = 'CHO(OO.) '
           s(1) = 1.0
           brtio = brch * s(1)
           CALL bratio(pchem(1),brtio,p(1),nref,ref)
           s(2) = 1.0
           p(2) = 'HO   '
           RETURN
! if criegee is "CO.(OO.)" then self-destruct, 100% yield.
         ELSE IF (xcri(1:9)=='CO.(OO.)*') THEN
           s(1) = 1.0
           p(1) = 'CO   '
           s(2) = 1.0
           p(2) = 'O2   '
           RETURN

         ELSE
           WRITE(6,'(a)') '--warning--(stop) in xcrieg'
           WRITE(6,'(a)') 'following C1 hot_criegge not found'
           WRITE(6,'(a)') 'in the routine. Change the program'
           WRITE(6,'(a)') 'accordingly.'
           WRITE(6,'(a)') xcri
           WRITE(waru,*) 'xcrieg',xcri !STOP
         ENDIF
       ENDIF

! -------------------------------------------
! 2nd section:  multi-carbon hot criegees:
! -------------------------------------------

       j1 = 0
       j2 = 0

! see if external or internal Criegee (internal:J2>0, external:J2=0):
       DO 210 i=1,mxnode
          IF (tbond(i,ig)==0) GOTO 210
          IF (j1==0) THEN
             j1 = i
          ELSE
             j2 = i
          ENDIF
 210   CONTINUE

! --------------------------------
!  open chapter "external criegee"  
! --------------------------------

! First the program checks if the criegee radical is an alpha
! carbonyl or alpha unsaturated carbon. If not, then perform
! the reaction for the "regular" R-CH.(OO.)* criegee.

       IF (j2==0) THEN
         IF (tgroup(ig)(1:9)=='CH.(OO.)*') THEN

! if criegee is "CHO-CH.(OO.)" then assume decomposition with 100% yield.
! The 2 channels (CHO+CO2+H, HCO+CO+OH) are arbitrarly set with equal 
! probability. This scheme is borrowed from the SAPRC99 mechanism and 
! is based on the fact that HCO-C.(OO.) has weaker bond energy than the
! corresponding CH3-CH.(OO.) radical
           IF (tgroup(j1)(1:3)=='CHO') THEN
             s(1) = 1.5
             p(1) = 'CO   '
             s(2) = 1.5
             p(2) = 'HO2  '
             s(3) = 0.5
             p(3) = 'HO   '
             s(4) = 0.50
             p(4) = 'CO2  '
             RETURN
           ENDIF

! if criegee is "R-CO-CH.(OO.)" then assume O shift to -CO- group with 
! 100% yield. This scheme is borrowed from the SAPRC99 mechanism and 
! is based on the fact that the resulting criegee is more stable than
! the first one. This reaction changes the external criegee into an 
! internal which is treated in the next chapter
           IF (tgroup(j1)(1:3)=='CO ') THEN

! swap criegee into aldehyde
             pold='CH.(OO.)*'
             pnew='CHO'
             tempgr=tgroup(ig)
             CALL swap(tempgr,pold,tgroup(ig),pnew)

! swap carbonyl into criegee
             pold='CO'
             pnew='C.(OO.)*'
             tempgr=tgroup(j1)
             CALL swap(tempgr,pold,tgroup(j1),pnew)

! store new position and new group matrix
             ig=j1
             DO i=1,mxnode
               group(i)=tgroup(i)
             ENDDO
             j1 = 0
             j2 = 0
             DO 220 i=1,mxnode
               IF (tbond(i,ig)==0) GOTO 220
               IF (j1==0) THEN
                 j1 = i
               ELSE
                 j2 = i
               ENDIF
220          CONTINUE

! jump to next chapter
             GOTO 456
           ENDIF

! if criegee is "C=C-CH.(OO.)" then assume H shift to give the
! corresponding alkene with 25 % yield. The remaining fraction is 
! expected to be stabilized (and hence yield the corresponding
! carboxylic acid). This scheme is borrowed from the SAPRC99 mechanism  
! and is based on detailed isoprene chemistry.
           IF (tgroup(j1)(1:2)=='Cd') THEN

! H shift
! -------

! add CO2
             np=np+1
             IF (np>mnr) THEN
                WRITE(6,'(a)') '--error-- in xcrieg'
                WRITE(6,'(a)') 'np is greater than mnr'
                WRITE(waru,*) 'xcrieg',xcri !STOP
             ENDIF
             s(np) = 0.25
             p(np) = 'CO2  '

! add corresponding alkene
             np=np+1
             IF (np>mnr) THEN
                WRITE(6,'(a)') '--error-- in xcrieg'
                WRITE(6,'(a)') 'np is greater than mnr'
                WRITE(waru,*) 'xcrieg',xcri !STOP
             ENDIF
             s(np) = 0.25

             tbond(ig,j1)=0
             tbond(j1,ig)=0

             IF (tgroup(j1)(1:3)=='CdH') THEN
               pold='CdH'
               pnew='CdH2'
             ELSE
               pold='Cd'
               pnew='CdH'
             ENDIF
             tempgr=tgroup(j1)
             CALL swap(tempgr,pold,tgroup(j1),pnew)
             CALL rebond(tbond,tgroup,pchem(np),nring)

! reset
             tgroup(:)=group(:)
             tbond(:,:)=bond(:,:)

! stabilization
! -------------
             sciyield=0.75
! ====== Stabilization: assume 50/50 for the cis (Z) and trans (E) isomer
             zeyield=sciyield/2.

! CIS (Z) isomer - replace hot_criegee by cold criegee and reset 
             pold='.(OO.)*' ;  pnew='.(ZOO.)' ; tempgr=tgroup(ig)  
             CALL swap(tempgr, pold, tgroup(ig), pnew)
             CALL rebond(tbond,tgroup,tempfo,nring)
! check radical, add species and coproducts in reaction products (s & p)
             CALL add2p(xcri,tempfo,zeyield,brch,np,s,p,nref,ref)
             tgroup(:)=group(:) ; tbond(:,:)=bond(:,:)

! TRANS (E) isomer - replace hot_criegee by cold criegee and reset 
             pold='.(OO.)*' ;  pnew='.(EOO.)' ; tempgr=tgroup(ig)  
             CALL swap(tempgr, pold, tgroup(ig), pnew)
             CALL rebond(tbond,tgroup,tempfo,nring)
! check radical, add species and coproducts in reaction products (s & p)
             CALL add2p(xcri,tempfo,zeyield,brch,np,s,p,nref,ref)
             tgroup(:)=group(:) ; tbond(:,:)=bond(:,:)

! jump to end criegee reaction
             GOTO 920             
           ENDIF

!  R-CH.(OO.) 
! If the program has reached this point, then the criegee radical react
! has a regular R-CH.(OO.)
! Atkinson, 1999, J. Phys. Chem. Ref. Data recommanded for CH3-CH.(OO.):
! -> CH3CH.(OO.) (stab.) : 15%
! -> CH3+CO+OH           : 54%
! -> CH3+CO2+H           : 17%
! -> CH4+CO              : 14%
! The above scheme was found to produce a too large amount
! of radical to match observation in smog chambers (see Carter,1999,
! SAPRC99 scheme). Measurements of Paulson et al 99 confirm the decrease
!  of OH yield with increasing length of radical.
! Therefore, we choose to affect for CH3-CH.(OO.) (Carter SAPRC99 in
! agreement with Rickard 99 OH measurements) :
! -> CH3CH.(OO.) (stab.) : 34%   (A)
! -> CH3+CO+OH           : 52%   (B)  
! -> CH4+CO              : 14%   (C)
! for RCH.(OO.), b ratios fixed to fit OH measurement from Paulson 99
! if R = 2C => (A) 0.54 ;(B) 0.46 ;(C) 0 ;
! if R = 3C => (A) 0.64 ;(B) 0.36 ;(C) 0 ;
! if R = 4C => (A) 0.76 ;(B) 0.24 ;(C) 0 ; 
! if R = 5C => (A) 0.84 ;(B) 0.16 ;(C) 0 ; 
! if R = 6C => (A) 0.92 ;(B) 0.08 ;(C) 0 ; 
! if R > 6C => (A) 1.   ;(B) 0    ;(C) 0 ; 
 
! Stabilization
! ------------- 
           IF (nca==2) ftherm = 0.34
           IF (nca==3) ftherm = 0.54
           IF (nca==4) ftherm = 0.64
           IF (nca==5) ftherm = 0.76
           IF (nca==6) ftherm = 0.84
           IF (nca==7) ftherm = 0.92
           IF (nca>7)  ftherm = 1.
           sciyield=ftherm
           
! ====== Stabilization: assume 50/50 for the cis (Z) and trans (E) isomer
           zeyield=sciyield/2.

! CIS (Z) isomer - replace hot_criegee by cold criegee and reset 
           pold='.(OO.)*' ;  pnew='.(ZOO.)' ; tempgr=tgroup(ig)  
           CALL swap(tempgr, pold, tgroup(ig), pnew)
           CALL rebond(tbond,tgroup,tempfo,nring)
! check radical, add species and coproducts in reaction products (s & p)
           CALL add2p(xcri,tempfo,zeyield,brch,np,s,p,nref,ref)
           tgroup(:)=group(:) ; tbond(:,:)=bond(:,:)

! TRANS (E) isomer - replace hot_criegee by cold criegee and reset 
           pold='.(OO.)*' ;  pnew='.(EOO.)' ; tempgr=tgroup(ig)  
           CALL swap(tempgr, pold, tgroup(ig), pnew)
           CALL rebond(tbond,tgroup,tempfo,nring)
! check radical, add species and coproducts in reaction products (s & p)
           CALL add2p(xcri,tempfo,zeyield,brch,np,s,p,nref,ref)
           tgroup(:)=group(:) ; tbond(:,:)=bond(:,:)

! H shift (molecular channel):
! ----------------------------
! This channel occurs only for C2 species.
! Old version was setting reaction products in the stack even
! if not produced.
           IF (nca==2) THEN
!c           IF (nca==2) fmol = 0.14
!c           IF (nca.ge.3) fmol = 0.
             fmol=0.14
! add CO2
             np=np+1
             IF (np>mnr) THEN
                WRITE(6,'(a)') '--error-- in xcrieg'
                WRITE(6,'(a)') 'np is greater than mnr'
                WRITE(waru,*) 'xcrieg',xcri !STOP
             ENDIF
             s(np) = fmol
             p(np) = 'CO2  '

! add corresponding molecule
             np    = np + 1
             IF (np>mnr) THEN
                WRITE(6,'(a)') '--error-- in xcrieg'
                WRITE(6,'(a)') 'np is greater than mnr'
                WRITE(waru,*) 'xcrieg',xcri !STOP
             ENDIF
             s(np) = fmol

             tgroup(ig) = ' '

             tbond(ig,j1) = 0
             tbond(j1,ig) = 0
            
             IF (tgroup(j1)(1:3)=='CH3') THEN
                p(np) = 'CH4  '
             ELSE
               IF (tgroup(j1)(1:3)=='CH2') THEN
                  pold = 'CH2'
                  pnew = 'CH3'
               ELSE IF (tgroup(j1)(1:2)=='CH') THEN
! JMLT Nov'15: decomposition of multi-substituted C1 coproducts
!              similar to hvdiss2.f
! NB: xcrieg 0.28 chan transfers H from adjacent =CdH, eliminates CO2
!     We might prefer to eliminate CO + OH (but do not do so here)
                  IF (INDEX(tgroup(j1),'(OOH)(ONO2)')/=0)THEN
                    pold='CH(OOH)(ONO2)'
                    pnew='CH2(OH)(OOH)'
                    p(np) = 'NO2  '
                    np = np+1
                    s(np) = fmol
                  ELSE IF (INDEX(tgroup(j1),'(ONO2)(OOH)')/=0)THEN
                    pold='CH(ONO2)(OOH)'
                    pnew='CH2(OH)(OOH)'
                    p(np) = 'NO2  '
                    np = np+1
                    s(np) = fmol
                  ELSE IF (INDEX(tgroup(j1),'(ONO2)(ONO2)')/=0)THEN
                    pold='CH(ONO2)(ONO2)'
                    pnew='CH2(OH)(ONO2)'
                    p(np) = 'NO2  '
                    np = np+1
                    s(np) = fmol
                  ELSE IF (INDEX(tgroup(j1),'(OOH)(OOH)')/=0)THEN
                    pold='CH(OOH)(OOH)'
                    pnew='CH2(OH)(OOH)'
                    p(np) = 'OH   '
                    np = np+1
                    s(np) = fmol
                  ELSE
                    pold = 'CH'
                    pnew = 'CH2'
                  ENDIF
               ELSE IF (tgroup(j1)(1:1)=='C') THEN
                  pold = 'C'
                  pnew = 'CH'
               ELSE
                 WRITE(6,'(a)') '--error-- in xcrieg'
                 WRITE(6,'(a)') 'group not found to perform'
                 WRITE(6,'(a)') 'the molecular channel for :'
                 WRITE(6,'(a)') xcri
                 WRITE(waru,*) 'xcrieg',xcri !STOP
               ENDIF

               tempgr=tgroup(j1)
               CALL swap(tempgr,pold,tgroup(j1),pnew)
               CALL rebond(tbond,tgroup,pchem(np),nring)
             ENDIF
! reset
             DO i=1,mxnode
               tgroup(i)=group(i)
               DO j=1,mxnode
                 tbond(i,j)=bond(i,j)
               ENDDO
             ENDDO

           ENDIF

! decomposition 1 channel (R. + CO + HO)
! --------------------------------------
! This channel occurs only for C<8 species.
! Old version was setting reaction products in the stack even
! if not produced.
           IF (nca>7) GOTO 42
           IF (nca==2) fdec1 = 0.52
           IF (nca==3) fdec1 = 0.46
           IF (nca==4) fdec1 = 0.36
           IF (nca==5) fdec1 = 0.24
           IF (nca==6) fdec1 = 0.16
           IF (nca==7) fdec1 = 0.08
           
! add CO
           np=np+1
           IF (np>mnr) THEN
              WRITE(6,'(a)') '--error-- in xcrieg'
              WRITE(6,'(a)') 'np is greater than mnr'
              WRITE(waru,*) 'xcrieg',xcri !STOP
           ENDIF
           s(np) = fdec1
           p(np) = 'CO   '

! add HO
           np=np+1
           IF (np>mnr) THEN
              WRITE(6,'(a)') '--error-- in xcrieg'
              WRITE(6,'(a)') 'np is greater than mnr'
              WRITE(waru,*) 'xcrieg',xcri !STOP
           ENDIF
           s(np) = fdec1
           p(np) = 'HO   '

! add corresponding radical
           np    = np + 1
           IF (np>mnr) THEN
              WRITE(6,'(a)') '--error-- in xcrieg'
              WRITE(6,'(a)') 'np is greater than mnr'
              WRITE(waru,*) 'xcrieg',xcri !STOP
           ENDIF
           s(np) = fdec1

           tgroup(ig) = ' '
           tbond(ig,j1) = 0
           tbond(j1,ig) = 0
            
! if the group next to the criegee radical is an ether,
! it must be converted into an alkoxy radical (ric nov 2008)
           IF (tgroup(j1)=='-O-') THEN
             DO i=1,mxnode
               IF ((tbond(i,j1)/=0).AND.(i/=ig)) THEN
                 tbond(i,j1) = 0
                 tbond(j1,i) = 0
                 tgroup(j1) = ' '
                 nc = INDEX(tgroup(i),' ')
                 tgroup(i)(nc:nc+3)='(O.)' 
               ENDIF
             ENDDO
           ELSE  
             nc = INDEX(tgroup(j1),' ')
             tgroup(j1)(nc:nc) = '.'
           ENDIF

           CALL rebond(tbond,tgroup,pchem(np),nring)
! reset
           DO i=1,mxnode
             tgroup(i)=group(i)
             DO j=1,mxnode
               tbond(i,j)=bond(i,j)
             ENDDO
           ENDDO

42         CONTINUE

! decomposition 2 channel (R. + CO2 + H)
! --------------------------------------
!c           fdec2=0.17
!* add CO2
!c           np=np+1
!c           IF (np>mnr) THEN
!c              WRITE(6,'(a)') '--error-- in xcrieg'
!c              WRITE(6,'(a)') 'np is greater than mnr'
!c              WRITE(waru,*) 'xcrieg',xcri !STOP
!c           ENDIF
!c           s(np) = fdec2
!c           p(np) = 'CO2  '

!* add HO2
!c           np=np+1
!c           IF (np>mnr) THEN
!c              WRITE(6,'(a)') '--error-- in xcrieg'
!c              WRITE(6,'(a)') 'np is greater than mnr'
!c              WRITE(waru,*) 'xcrieg',xcri !STOP
!c           ENDIF
!c           s(np) = fdec2
!c           p(np) = 'HO2  '

!* add corresponding radical
!c           np    = np + 1
!c           IF (np>mnr) THEN
!c              WRITE(6,'(a)') '--error-- in xcrieg'
!c              WRITE(6,'(a)') 'np is greater than mnr'
!c              WRITE(waru,*) 'xcrieg',xcri !STOP
!c           ENDIF
!c           s(np) = fdec2

!c           tgroup(ig) = ' '
!c           tbond(ig,j1) = 0
!c           tbond(j1,ig) = 0
            
!c           nc = INDEX(tgroup(j1),' ')
!c           tgroup(j1)(nc:nc) = '.'
!c           CALL rebond(tbond,tgroup,pchem(np),nring)

!* reset
!c           DO i=1,mxnode
!c             tgroup(i)=group(i)
!c             DO j=1,mxnode
!c               tbond(i,j)=bond(i,j)
!c             ENDDO
!c           ENDDO

!* go to assignment 920, where the products are checked to either be 
!* or not be loaded in the stack.
           GO TO 920

!* escape routes for products of ring-opening (John Orlando, no reference)
!* if criegee is "R-C(OH).(OO.)*" then assume conversion to acid with NO, 100% yield.
         ELSE IF (tgroup(ig)(1:12)=='C(OH).(OO.)*') THEN
           s(2) = 1.0
           p(2) = 'NO2  '
           s(1) = 1.0
           pold='C(OH).(OO.)*'
           pnew='CO(OH)'
           tempgr=tgroup(ig)
           CALL swap(tempgr,pold,tgroup(ig),pnew)
           CALL rebond(tbond,tgroup,pchem(1),nring)
           brtio = brch * s(1)
           CALL bratio(pchem(1),brtio,p(1),nref,ref)
           GOTO 920             
!* if criegee is "R-C(ONO2).(OO.)*" then assume NO2 elimination, 100% yield.
         ELSE IF (tgroup(ig)(1:14)=='C(ONO2).(OO.)*') THEN
           s(2) = 1.0
           p(2) = 'NO2  '
           s(1) = 1.0
           pold='C(ONO2).(OO.)*'
           pnew='CO(OO.)'
           tempgr=tgroup(ig)
           CALL swap(tempgr,pold,tgroup(ig),pnew)
           CALL rebond(tbond,tgroup,pchem(1),nring)
           CALL stdchm(pchem(1))
           brtio = brch * s(1)
           CALL bratio(pchem(1),brtio,p(1),nref,ref)
           GOTO 920             
!* if criegee is "R-C(OOH).(OO.)*" then assume OH elimination, 100% yield.
         ELSE IF (tgroup(ig)(1:13)=='C(OOH).(OO.)*') THEN
           s(2) = 1.0
           p(2) = 'HO   '
           s(1) = 1.0
           pold='C(OOH).(OO.)*'
           pnew='CO(OO.)'
           tempgr=tgroup(ig)
           CALL swap(tempgr,pold,tgroup(ig),pnew)
           CALL rebond(tbond,tgroup,pchem(1),nring)
           CALL stdchm(pchem(1))
           brtio = brch * s(1)
           CALL bratio(pchem(1),brtio,p(1),nref,ref)
           GOTO 920             
!* if criegee is "R-C(OOOH).(OO.)*" then assume OOH elimination, 100% yield.
         ELSE IF (tgroup(ig)(1:14)=='C(OOOH).(OO.)*') THEN
           s(2) = 1.0
           p(2) = 'HO2  '
           s(1) = 1.0
           pold='C(OOOH).(OO.)*'
           pnew='CO(OO.)'
           tempgr=tgroup(ig)
           CALL swap(tempgr,pold,tgroup(ig),pnew)
           CALL rebond(tbond,tgroup,pchem(1),nring)
           CALL stdchm(pchem(1))
           brtio = brch * s(1)
           CALL bratio(pchem(1),brtio,p(1),nref,ref)
           GOTO 920             

         ELSE
!* If the program continues up to this point for external criegee then
!* the criegee is not a -CH.(OO.) criegee. The program will stop. Add
!* some additional program lines if required.
           WRITE(6,'(a)') '--warning--(stop) in xcrieg'
           WRITE(6,'(a)') 'following external hot_criegge cannot be'
           WRITE(6,'(a)') 'computed. Please add some'
           WRITE(6,'(a)') 'additional fortran lines'
           WRITE(6,'(a)') xcri
           WRITE(waru,*) 'xcrieg',xcri 
           STOP
         ENDIF

!* end external criegee
       ENDIF



!* -----------------------------------------------
!* else open chapter "internal criegees" 
!* -----------------------------------------------

456    CONTINUE

!* For -CH-C.(OO.)-C criegee radical, major evolution pathway is
!* expected to be the hydroperoxide channel (e.g. see Atkinson, 1999,
!* J. Chem. Ref. Data) :
!* -CH-C.(OO.)-C => -C=C(OOH)-C => -C.-CO-C + OH 
!* This pathway is set with 100% yield. This mechanism requires an
!* H for the group in alpha position. Note that for H in -CHO, the 
!* above mechanism is expected to not occur, because of the formation
!* of a strained transition state (see Carter, 1999, SAPRC99 mechanism).
!* Therefore, the program first checks the two group at alpha position 
       IF (j2/=0) THEN

!* check that that criegee make sense
         IF (tgroup(ig)(1:8)/='C.(OO.)*') THEN
           WRITE(6,'(a)') '--warning--(stop) in xcrieg'
           WRITE(6,'(a)') 'something wrong for the following'
           WRITE(6,'(a)') 'criegee radical (was first expected'
           WRITE(6,'(a)') 'to be C.(OO.)*'
           WRITE(6,'(a)') xcri
           WRITE(waru,*) 'xcrieg',xcri !STOP
         ENDIF

!* check for CH groups
         posj1=1
         posj2=1
         IF (INDEX(tgroup(j1),'CH')==0) posj1=0
         IF (INDEX(tgroup(j2),'CH')==0) posj2=0

!* check for CH groups (but not CHO)
! I'd put this in for cyclic compounds: removed later for test
!         posj1=1
!         posj2=1
!         IF (INDEX(tgroup(j1),'CH')==0 .OR.
!     &       INDEX(tgroup(j1),'CHO')/=0) posj1=0
!         IF (INDEX(tgroup(j2),'CH')==0 .OR.
!     &       INDEX(tgroup(j2),'CHO')/=0) posj2=0

!* no CH group.
!* ------------ 

! Carter (SAPRC99) arbitrarily assume that 90% of this type of Criegee
! is stabilized and 10% decomposes to CO2 and 2 R.
! We consider that the thermalised radical only reacts with H2O to 
! product H2O2 + ketone (Baker 2002 assumed it is the only way for  
! R1R2C(OOH)(OH) )

         IF ( (posj1==0).AND.(posj2==0) ) THEN

! 10% decomposition
! if the carbon bearing the criegee belongs to a ring, this channel
! lead to a di-radical species, which is not currently treated in the
! generator. need update. (ric 2008)
           CALL findring(ig,j1,ngr,tbond,rngflg,ring)
           IF (ring(ig)==1) GOTO 200
           np=np+1 ; s(np)=0.1

           tbond(ig,j1) = 0
           tbond(j1,ig) = 0
           tbond(ig,j2) = 0
           tbond(j2,ig) = 0
           tgroup(ig) = ' '
           nc = INDEX(group(j1),' ')
           tgroup(j1)(nc:nc) = '.'
           nc = INDEX(group(j2),' ')
           tgroup(j2)(nc:nc) = '.'
! fragment and write in correct format:
           CALL fragm(tbond,tgroup,pchem1(np),pchem2(np))
           CALL stdchm(pchem1(np))
           CALL stdchm(pchem2(np))
           prod2 = pchem2(np)
! check pchem1, write in standard format:
!c          CALL radchk(pchem1(np),prod,coprod)
           CALL radchk(pchem1(np),rdckprod,rdcktprod,nip,sc,nref,ref)
           prod = rdckprod(1)
           IF (nip/=1) WRITE(6,*) 'xcrieg.f coprod à revoir'
           IF (nip/=1) WRITE(6,*) prod
           IF (nip/=1) STOP
           DO l = 1,mxnode
             coprod(l) = rdcktprod(1,l)
           ENDDO
           CALL stdchm(prod)
           DO l=1,mxnode
             IF (coprod(l)/=' ') THEN
               np = np + 1
               s(np) = 0.1
               p(np) = coprod(l)
             ENDIF
           ENDDO
           IF (np+1>mnr) THEN
             WRITE(6,'(a)') '--error-- in xcrieg'
             WRITE(6,'(a)') 'np is greater than mnr'
             WRITE(waru,*) 'xcrieg',xcri !STOP
           ENDIF
           
           np = np + 1
           s(np) = 0.1
           brtio = brch * s(np)

           CALL bratio(prod,brtio,p(np),nref,ref)
           IF (rdtcopchem>0.) THEN
             np = np + 1
             IF (np>mnp) then
               WRITE(6,'(a)') '--error-- in xcrieg'
               WRITE(6,'(a)') 'np is greater than mnp'
               WRITE(6,'(a)') '(too much product in the reaction)'
               WRITE(waru,*) 'xcrieg',xcri !STOP
             ENDIF
             s(np) = rdtcopchem
             p(np) = copchem
           ENDIF
           IF (np+2>mnr) THEN
              WRITE(6,'(a)') '--error-- in xcrieg'
              WRITE(6,'(a)') 'np is greater than mnr'
              WRITE(waru,*) 'xcrieg',xcri !STOP
           ENDIF
           
           
! check prod2 = pchem2, write in standard format:

           CALL radchk(prod2,rdckprod,rdcktprod,nip,sc,nref,ref)
           prod = rdckprod(1)
           IF (nip/=1) WRITE(6,*) 'xcrieg.f coprod à revoir 2'
           IF (nip/=1) WRITE(6,*) prod
           IF (nip/=1) STOP

           coprod(:) = rdcktprod(1,:)

           CALL stdchm(prod)
           DO l=1,mxnode
             IF (coprod(l)/=' ') THEN
             np = np + 1
             s(np) = 0.1
             p(np) = coprod(l)
             ENDIF
           ENDDO
           
           np = np + 1
           s(np) = 0.1
           brtio = brch * s(np)

           CALL bratio(prod,brtio,p(np),nref,ref)
           np=np+1
           s(np) = 0.1
           p(np) = 'CO2'
           IF (rdtcopchem>0.) THEN
             np = np + 1
             IF (np>mnp) then
               WRITE(6,'(a)') '--error-- in xcrieg'
               WRITE(6,'(a)') 'np is greater than mnp'
               WRITE(6,'(a)') '(too much product in the reaction)'
               WRITE(waru,*) 'xcrieg',xcri !STOP
             ENDIF
             s(np) = rdtcopchem
             p(np) = copchem
           ENDIF
! reset
           DO i=1,mxnode
             tgroup(i)=group(i)
             DO j=1,mxnode
               tbond(i,j)=bond(i,j)
             ENDDO
           ENDDO
           
200        CONTINUE
           
! 90% product H2O2 + ketone
           sciyield = 0.9
           IF (ring(ig)==1) sciyield = 1.

! ====== Stabilization: assume 50/50 for the cis (Z) and trans (E) isomer
           zeyield=sciyield/2.

! CIS (Z) isomer - replace hot_criegee by cold criegee and reset 
           pold='.(OO.)*' ;  pnew='.(ZOO.)' ; tempgr=tgroup(ig)  
           CALL swap(tempgr, pold, tgroup(ig), pnew)
           CALL rebond(tbond,tgroup,tempfo,nring)
! check radical, add species and coproducts in reaction products (s & p)
           CALL add2p(xcri,tempfo,zeyield,brch,np,s,p,nref,ref)
           tgroup(:)=group(:) ; tbond(:,:)=bond(:,:)

! TRANS (E) isomer - replace hot_criegee by cold criegee and reset 
           pold='.(OO.)*' ;  pnew='.(EOO.)' ; tempgr=tgroup(ig)  
           CALL swap(tempgr, pold, tgroup(ig), pnew)
           CALL rebond(tbond,tgroup,tempfo,nring)
! check radical, add species and coproducts in reaction products (s & p)
           CALL add2p(xcri,tempfo,zeyield,brch,np,s,p,nref,ref)
           tgroup(:)=group(:) ; tbond(:,:)=bond(:,:)
        
! reset
           tgroup(:)=group(:)
           tbond(:,:)=bond(:,:)
           RETURN
         ENDIF

! CH in alpha-position
! --------------------
! put yield for the various channels. If the reaction can 
! occur at the two positions, then ratio is arbitrarily set
! using rate constant for OH reaction at the given site. A
! cut off of 5% is used

         IF ( (posj1==1).AND.(posj2==0) ) THEN
           yld1=1.
           yld2=0.
         ENDIF
         IF ( (posj1==0).AND.(posj2==1) ) THEN
           yld1=0.
           yld2=1.
         ENDIF
         IF ( (posj1==1).AND.(posj2==1) ) THEN
           CALL rabsoh(tbond,tgroup,j1,garrhc,nring)
           c1 = garrhc(1) * 298.**garrhc(2) * exp(-garrhc(3)/298.)
           c1 = AMIN1(c1,2.0E-10)
           CALL rabsoh(tbond,tgroup,j2,garrhc,nring)
           c2 = garrhc(1) * 298.**garrhc(2) * exp(-garrhc(3)/298.)
           c2 = AMIN1(c2,2.0E-10)

           yld1=c1/(c1+c2)
           yld2=c2/(c1+c2)
           IF (yld1.LT.0.05) THEN
             yld1=0.
             yld2=1.
           ENDIF
           IF (yld2.LT.0.05) THEN
             yld1=1.
             yld2=0.
           ENDIF
         ENDIF

! CH at position 1
! ----------------
         IF (yld1>0.) THEN

! add OH
           np=np+1
           IF (np>mnr) THEN
              WRITE(6,'(a)') '--error-- in xcrieg'
              WRITE(6,'(a)') 'np is greater than mnr'
              WRITE(waru,*) 'xcrieg',xcri !STOP
           ENDIF
           s(np) = yld1
           p(np) = 'HO   '

! add corresponding molecule
           np    = np + 1
           IF (np>mnr) THEN
              WRITE(6,'(a)') '--error-- in xcrieg'
              WRITE(6,'(a)') 'np is greater than mnr'
              WRITE(waru,*) 'xcrieg',xcri !STOP
           ENDIF
           s(np) = yld1

! swap C.(OO).* into CO and add radical dot at j1
           tgroup(ig) = 'CO'
           IF (tgroup(j1)(1:3)=='CH3') THEN
                pold = 'CH3'
                pnew = 'CH2'
           ELSE IF (tgroup(j1)(1:3)=='CH2') THEN
                pold = 'CH2'
                pnew = 'CH'
           ELSE IF (tgroup(j1)(1:2)=='CH') THEN
                pold = 'CH'
                pnew = 'C'
           ELSE
               WRITE(6,'(a)') '--error-- in xcrieg'
               WRITE(6,'(a)') 'group not found to perform'
               WRITE(6,'(a)') 'the molecular channel for :'
               WRITE(6,'(a)') xcri
               WRITE(waru,*) 'xcrieg',xcri !STOP
           ENDIF

           tempgr=tgroup(j1)
           CALL swap(tempgr,pold,tgroup(j1),pnew)
           nc = INDEX(tgroup(j1),' ')
           tgroup(j1)(nc:nc) = '.'
           CALL rebond(tbond,tgroup,pchem(np),nring)

! reset
           DO i=1,mxnode
             tgroup(i)=group(i)
             DO j=1,mxnode
               tbond(i,j)=bond(i,j)
             ENDDO
           ENDDO
         ENDIF

! CH at position 2
! ----------------
         IF (yld2>0.) THEN

! add OH
           np=np+1
           IF (np>mnr) THEN
              WRITE(6,'(a)') '--error-- in xcrieg'
              WRITE(6,'(a)') 'np is greater than mnr'
              WRITE(waru,*) 'xcrieg',xcri !STOP
           ENDIF
           s(np) = yld2
           p(np) = 'HO   '

! add corresponding molecule
           np    = np + 1
           IF (np>mnr) THEN
              WRITE(6,'(a)') '--error-- in xcrieg'
              WRITE(6,'(a)') 'np is greater than mnr'
              WRITE(waru,*) 'xcrieg',xcri !STOP
           ENDIF
           s(np) = yld2

! swap C.(OO).* into CO and add radical dot at j1
           tgroup(ig) = 'CO'
           IF (tgroup(j2)(1:3)=='CH3') THEN
                pold = 'CH3'
                pnew = 'CH2'
           ELSE IF (tgroup(j2)(1:3)=='CH2') THEN
                pold = 'CH2'
                pnew = 'CH'
           ELSE IF (tgroup(j2)(1:2)=='CH') THEN
                pold = 'CH'
                pnew = 'C'
           ELSE
               WRITE(6,'(a)') '--error-- in xcrieg'
               WRITE(6,'(a)') 'group not found to perform'
               WRITE(6,'(a)') 'the molecular channel for :'
               WRITE(6,'(a)') xcri
               WRITE(waru,*) 'xcrieg',xcri !STOP
           ENDIF

           tempgr=tgroup(j2)
           CALL swap(tempgr,pold,tgroup(j2),pnew)
           nc = INDEX(tgroup(j2),' ')
           tgroup(j2)(nc:nc) = '.'
           CALL rebond(tbond,tgroup,pchem(np),nring)

! reset
           DO i=1,mxnode
             tgroup(i)=group(i)
             DO j=1,mxnode
               tbond(i,j)=bond(i,j)
             ENDDO
           ENDDO
         ENDIF

! end internal criegee
       ENDIF

! ------------------------------------------------------
! check the various product and load species in the
! stack (if required).
! ------------------------------------------------------
920   nold = np
      
      DO 950 i=1,nold
         IF (pchem(i)(1:1)==' ') GO TO 950
         prod = pchem(i)

         pchem_del=' '
         coprod_del(:) = ' '

         IF (INDEX(pchem(i),'.')/=0) THEN
            CALL radchk(prod,rdckprod,rdcktprod,nip,sc,nref,ref)
            pchem(i) = rdckprod(1)
            IF (nip==2) THEN
              pchem_del = rdckprod(2)
              DO j=1,mxnode
                coprod_del(j) = rdcktprod(2,j)
              ENDDO
            ENDIF
            DO l = 1,mxnode
              coprod(l) = rdcktprod(1,l)
            ENDDO
            DO j=1,mxnode
               IF (coprod(j)(1:1)/=' ') THEN
                  np = np + 1
                  IF (np>mnr) THEN
                     WRITE(6,'(a)') '--error-- in xcrieg'
                     WRITE(6,'(a)') 'np is greater than mnr'
                     WRITE(waru,*) 'xcrieg',xcri !STOP
                  ENDIF
                  s(np)=s(i)
                  p(np) = coprod(j)
               ENDIF
            ENDDO
            DO j=1,mxnode
               IF (coprod_del(j)(1:1)/=' ') THEN
                   WRITE(6,*) pchem_del
                   WRITE(6,*) coprod_del(j)
                  np = np + 1
                  IF (np>mnr) THEN
                     WRITE(6,'(a)') '--error-- in xcrieg'
                     WRITE(6,'(a)') 'np is greater than mnr'
                     WRITE(waru,*) 'xcrieg',xcri !STOP
                  ENDIF
                  s(np)=sc(2)
                  p(np) = coprod_del(j)
               ENDIF
            ENDDO
         ENDIF

         CALL stdchm(pchem(i))
         IF (nip==2) s(i)=s(i)*sc(1)
         brtio = brch * s(i)
         CALL bratio(pchem(i),brtio,p(i),nref,ref)
         IF (nip==2) THEN
           np=np+1
           s(np)=(s(i)*sc(2))/sc(1)
           brtio = brch * s(np)
           CALL bratio(pchem_del,brtio,p(np),nref,ref)
         ENDIF
         IF (rdtcopchem>0.) THEN
           np = np + 1
           IF (np>mnp) then
             WRITE(6,'(a)') '--error-- in ho_rad'
             WRITE(6,'(a)') 'np is greater than mnp'
             WRITE(6,'(a)') '(too much product in the reaction)'
             WRITE(waru,*) 'xcrieg',xcri !STOP
           ENDIF
           s(np) = rdtcopchem
           p(np) = copchem
         ENDIF

!      IF(wtflag/=0) WRITE(*,*) '*xcrieg* : product : ',i,pchem(i)
950   CONTINUE

END SUBROUTINE xcrieg

END MODULE stuff4hvdiss2
