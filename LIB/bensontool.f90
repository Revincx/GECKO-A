MODULE bensontool
IMPLICIT NONE
CONTAINS

! - SUBROUTINE rdbenson()
! - REAL FUNCTION heat(chem)
! - SUBROUTINE nameben(base,point,ic,icdot,io,ico,icodot,icd,init,ialko,ipero,bengrp)
! - SUBROUTINE getben(bengrp,bvalue,check)

!=======================================================================
! PURPOSE: read heat of formation of the benson groups
! The data are stored in module database as:     
!   - nbson        : total number of benson group              
!   - bsongrp(:)   : table of benson group                     
!   - bsonval(:)   : heat of formation corresponding to group i
!=======================================================================
SUBROUTINE rdbenson()
  USE keyparameter, ONLY: tfu1,dirgecko
  USE database, ONLY: mxbg,lbg,nbson,bsonval,bsongrp
  IMPLICIT NONE

  INTEGER :: ierr
  CHARACTER(LEN=100) line

! initialize
  nbson=0  ;  bsongrp(:)=' '  ;  bsonval(:)=0.

! open the file
  OPEN(UNIT=tfu1,FILE=TRIM(dirgecko)//'DATA/benson.dat',STATUS='OLD',FORM='FORMATTED')

! read benson group
  rdloop: DO
    READ(tfu1,'(a)',IOSTAT=ierr) line
    IF (ierr /= 0) THEN
      WRITE (6,*) '--error-- in rdbenson. Missing keyword END ?'
      STOP "in rdbenson, while reading inputs"
    ENDIF
    IF (line(1:3)=='END') EXIT rdloop
    IF (line(1:1)=='*') CYCLE rdloop

    IF (line(25:25)/=' ') THEN
      WRITE(6,*)'--error--, reading bensongrp.dat. Group > 24 characters'
      WRITE(6,*) TRIM(line)
      STOP "in rdbenson"
    ENDIF
    nbson=nbson+1
    IF (nbson >= mxbg) THEN
      WRITE (6,'(a)') '--error--, while reading bensongrp.dat'
      WRITE (6,'(a)') 'number of benson grp is greater than mbg'
      STOP "in rdbenson"
    ENDIF
    READ(line,'(a24,1x,f10.3)',IOSTAT=ierr) bsongrp(nbson),bsonval(nbson)
    IF (ierr /= 0) THEN
      WRITE (6,*) '--error-- in rdbenson, while reading bensongrp.dat:'
      WRITE(6,*) TRIM(line)
      STOP "in rdbenson, , while reading inputs"
    ENDIF
  ENDDO rdloop
  CLOSE (tfu1)
END SUBROUTINE rdbenson

!=======================================================================
! PURPOSE: return the heat (enthalpie) formation for the chemical 
! provided as input based on the Benson group contribution method.
!=======================================================================
REAL FUNCTION heat(chem)
  USE keyparameter, ONLY: mxnode,mxlgr,mxring
  USE database, ONLY:lbg
  USE rjtool
  USE stdgrbond
  IMPLICIT NONE

  CHARACTER(LEN=*),INTENT(IN) :: chem

! internal:
  CHARACTER(LEN=3)   :: base      ! base of the benson group
  CHARACTER(LEN=lbg) :: bengrp    ! benson group
  CHARACTER(LEN=mxlgr) :: tgroup(mxnode) ! group (gecko) table
  CHARACTER(LEN=mxlgr) :: tempkg         ! temp group
  INTEGER  :: tbond(mxnode,mxnode)       ! bond matrix
  INTEGER  :: nca                   ! # of nodes
  INTEGER  :: point, numo, nring
  INTEGER  :: dbflg
  REAL     :: bvalue
  INTEGER  :: i,j,ic,init,ir,icd,io,ico,nc
  INTEGER  :: ib,ie,fg,op(4),cp(4),icdot,icodot,ialko,ipero
  INTEGER  :: check
  INTEGER      :: rjg(mxring,2)
  INTEGER,PARAMETER  :: wtflag=0

  IF (wtflag/=0) WRITE(6,*) '*heat*'

! initialize:
  heat = 0.  ;  bvalue  = 0.  ;  point=0  ;  numo=0
  base = ' ' ;  bengrp = ' '  ;  tempkg= ' '
  tbond(:,:)=0  ; tgroup(:)=' '

  nc = INDEX(chem,' ')-1
  IF (nc < 1) RETURN
  CALL grbond(chem,tgroup,tbond,dbflg,nring)
  IF(nring/=0) CALL rjgrm(nring,tgroup,rjg)

! count the number of node in chem
  nca=0
  DO i=1, mxnode ; IF (tgroup(i)(1:1)/=' ') nca=nca+1 ; ENDDO

! treat some one-carbon molecules as special cases:
  IF (nca < 2) THEN
    IF (tgroup(1)=='CH2O') THEN 
       heat = -26.0 ; RETURN
    ELSE IF ((tgroup(1)=='CO  ').AND.(tgroup(2)(1:1)/='C')) THEN
       heat = -26.4 ; RETURN
    ELSE IF (tgroup(1)=='CH3O') THEN
       heat = 4.62 ; RETURN
    ELSE IF (tgroup(1)=='CH3.') THEN
       heat = 34.9 ; RETURN
    ELSE IF (tgroup(1)=='CHO.') THEN
       heat = 9.99 ; RETURN
    ENDIF
  ENDIF

! treat multinodes:
  grloop: DO i=1,nca
    tempkg = tgroup(i)
    IF (tempkg(1:2)=='  ') CYCLE grloop

!  -O- node
! -------------
 
    IF (tempkg(1:3)=='-O-') THEN
      ic=0 ; icdot=0 ; ico=0 ; io=0 ; icodot=0 ; init=0 ; ir=0
      base = ' '

! find bonds for C based groups, 1st check that the C is not a radical.
      point = 2
      base  = 'O'

      DO j=1,mxnode   ! find O-C  bonds:
        ir = 0
        IF (tbond(i,j)/=0) THEN
          nc = INDEX(tgroup(j),' ') - 1
          IF (tgroup(j)(nc:nc)=='.') ir=1
      
          IF ((tgroup(j)(1:2)=='CO') .OR. (tgroup(j)(1:3)=='CHO')) THEN
            IF (ir==1) THEN ; icodot=icodot+1 ; ELSE ; ico=ico+1 ; ENDIF
          ELSE
            IF (ir==1) THEN ; icdot=icdot+1 ; ELSE ; ic=ic+1 ; ENDIF
          ENDIF
        ENDIF 
      ENDDO
      ialko = 0 ; ipero = 0 ; icd = 0
      CALL nameben(base,point,ic,icdot,io,ico,icodot,icd,init,ialko,ipero,bengrp)
      CALL getben(bengrp,bvalue,check)
      IF (check/=0) THEN 
        WRITE (6,*) '--error-- in heat. Benson group not found:',TRIM(bengrp)
        WRITE (6,*) 'for chemical: ',TRIM(chem)
        STOP "in heat"
      ENDIF
      heat=heat+bvalue
    ENDIF    !'-O-'

 ! carbon center
 ! -------------
    IF (tempkg(1:1)=='C') THEN
      ic=0 ; icdot=0 ; ico=0 ; io=0 ; icodot=0 ; init=0 ; ir=0
      icd=0 ; ialko=0 ; ipero=0 ; base = ' '

! treat some special groups. 
! ketene: DHf C2H2O = -47.6 => DHf ketene group = -47.6 - Cd_(Cd) [6.26] = -53.86
! reference for ketene: Cohen, N. (1996) J.Phys.Chem.Ref.Data, 25(6), 1411-1481
      IF (tempkg(1:3)=='CdO') THEN
        bvalue = -53.86 ; heat = heat + bvalue
        CYCLE grloop
      ENDIF

! set the base 
      nc = INDEX(tempkg,' ') - 1
      IF (tempkg(nc:nc)=='.') ir=1
      IF (tempkg(1:2)=='CO') THEN
        IF (ir==1) THEN
           point=4 ; base = 'CO*'
        ELSE
           point=3 ; base = 'CO'
        ENDIF

      ELSE IF (tempkg(1:3)=='CHO') THEN 
        IF (ir==1) THEN
          WRITE (6,*) '--error--, Benson group HCO* unexpected in: ', TRIM(chem)
          STOP "in heat"
        ELSE
           point=3 ; base = 'CO'
        ENDIF

      ELSE IF (tempkg(1:2)=='Cd') THEN
        IF (ir==1) THEN
          WRITE (6,*) '--error--, Benson group Cd* unexpected in: ', TRIM(chem)
          STOP "in heat"
        ELSE
          point=3 ; base = 'Cd'
        ENDIF

      ELSE
        IF (ir==1) THEN
          point=3 ; base = 'C*'
        ELSE
          point=2 ; base = 'C'
        ENDIF
      ENDIF

! find node-node bond :
      DO j=1,mxnode
        ir = 0
        IF (tbond(i,j)/=0) THEN
          nc = INDEX(tgroup(j),' ') - 1
          IF (tgroup(j)(nc:nc)=='.') ir=1
		  
          IF ((tgroup(j)(1:2)=='CO') .OR. (tgroup(j)(1:3)=='CHO')) THEN
            IF (ir==1) THEN ; icodot=icodot+1 ; ELSE ; ico=ico+1 ; ENDIF
		  
          ELSE IF (tgroup(j)(1:3)=='-O-') THEN ; io  = io+1
          ELSE IF (tgroup(j)(1:2)=='Cd') THEN  ; icd = icd+1
          ELSE
            IF (ir==1) THEN ; icdot=icdot+1 ; ELSE ; ic=ic+1 ; ENDIF
          ENDIF
        ENDIF 
      ENDDO

! find carbon - functional group bonds:
! first the program finds the open and closing parenthesis (since functional
! groups are included into parenthesis) and count number of groups.
      fg=0 ; op(:)=0 ; cp(:)=0
      DO j=point-1,mxlgr
        IF (tempkg(j:j)=='(') THEN ; fg=fg+1 ; op(fg)=j ; ENDIF
        IF (tempkg(j:j)==')') cp(fg)=j
      ENDDO

! functional groups (in parenthesis)
      IF (fg > 0) THEN
        DO j=1,fg
          ib=op(j) ; ie=cp(j)
          IF (tempkg(ib:ie)=='(O.)') THEN         !--- alkoxy
            ialko=ialko+1

          ELSE IF (tempkg(ib:ie)=='(OO.)') THEN   !--- peroxy 
            ipero = ipero + 1

          ELSE IF (tempkg(ib:ie)=='(OH)') THEN    !--- alcohols
            io = io + 1
            bengrp='O_('//base(1:point-1)
            nc=INDEX(bengrp,' ')
            bengrp(nc:nc)=')'
            CALL getben(bengrp,bvalue,check)
            IF (check/=0) THEN 
              WRITE (6,*) '--error-- in heat. Benson group not found:',TRIM(bengrp)
              WRITE (6,*) 'for chemical: ',TRIM(chem)
              STOP "in heat"
            ENDIF
            heat=heat+bvalue

          ELSE IF ((tempkg(ib:ie)=='(OOH)').OR.(tempkg(ib:ie)=='(OOOH)')) THEN   !--- hydroperoxides
            io = io + 1
            bengrp='O_('//base(1:point-1)
            nc=INDEX(bengrp,' ')
            bengrp(nc:nc+3)=')(O)'
            CALL getben(bengrp,bvalue,check)
            IF (check/=0) THEN 
              WRITE (6,*) '--error-- in heat. Benson group not found:',TRIM(bengrp)
              WRITE (6,*) 'for chemical: ',TRIM(chem)
              STOP "in heat"
            ENDIF
            heat=heat+bvalue
            bengrp='O_(O)'
            CALL getben(bengrp,bvalue,check)
            IF (check/=0) THEN 
              WRITE (6,*) '--error-- in heat. Benson group not found:',TRIM(bengrp)
              WRITE (6,*) 'for chemical: ',TRIM(chem)
              STOP "in heat"
            ENDIF
            heat=heat+bvalue

          ELSE IF (tempkg(ib:ie)=='(ONO2)') THEN  !--- nitrate
            init = init + 1
            bengrp='ONO2_('//base(1:point-1)
            nc=INDEX(bengrp,' ')
            bengrp(nc:nc)=')'
            CALL getben(bengrp,bvalue,check)
            IF (check/=0) THEN 
              WRITE (6,*) '--error-- in heat. Benson group not found:',TRIM(bengrp)
              WRITE (6,*) 'for chemical: ',TRIM(chem)
              STOP "in heat"
            ENDIF
            heat=heat+bvalue

          ELSE IF (tempkg(ib:ie)=='(OONO2)') THEN  !--- pernitrate
            io = io + 1
            bengrp='O_('//base(1:point-1)
            nc=INDEX(bengrp,' ')
            bengrp(nc:nc+7)=')(ONO2)'
            CALL getben(bengrp,bvalue,check)
            IF (check/=0) THEN 
              WRITE (6,*) '--error-- in heat. Benson group not found:',TRIM(bengrp)
              WRITE (6,*) 'for chemical: ',TRIM(chem)
              STOP "in heat"
            ENDIF
            heat=heat+bvalue
            bengrp='ONO2_(O)'
            CALL getben(bengrp,bvalue,check)
            IF (check/=0) THEN 
              WRITE (6,*) '--error-- in heat. Benson group not found:',TRIM(bengrp)
              WRITE (6,*) 'for chemical: ',TRIM(chem)
              STOP "in heat"
            ENDIF
            heat=heat+bvalue

          ELSE IF (tempkg(ib:ie)=='(NO2)') THEN    !--- nitro group
            io = io + 1
            bengrp='O_('//base(1:point-1)
            nc=INDEX(bengrp,' ')
            bengrp(nc:nc+6)=')(NO2)'
            CALL getben(bengrp,bvalue,check)
            IF (check/=0) THEN 
              WRITE (6,*) '--error-- in heat. Benson group not found:',TRIM(bengrp)
              WRITE (6,*) 'for chemical: ',TRIM(chem)
              STOP "in heat"
            ENDIF
            heat=heat+bvalue
            bengrp='NO2_(O)'
            CALL getben(bengrp,bvalue,check)
            IF (check/=0) THEN 
              WRITE (6,*) '--error-- in heat. Benson group not found:',TRIM(bengrp)
              WRITE (6,*) 'for chemical: ',TRIM(chem)
              STOP "in heat"
            ENDIF
            heat=heat+bvalue

          ELSE
            WRITE (6,*) '--error-- in heat. Benson group contribution can not'
            WRITE (6,*) 'be estimated for the group: ',TRIM(tempkg)
            WRITE (6,*) 'in the species: ',TRIM(chem)
            STOP "in heat"
          ENDIF
        ENDDO  ! functional group (fg) loop
      ENDIF 

      CALL nameben(base,point,ic,icdot,io,ico,icodot,icd,init,ialko,ipero,bengrp)
      CALL getben(bengrp,bvalue,check)
      IF (check/=0) THEN 
        WRITE (6,*) '--error-- in heat. Benson group not found:',TRIM(bengrp)
        WRITE (6,*) 'for chemical: ',TRIM(chem)
        STOP "in heat"
      ENDIF
      heat=heat+bvalue
    ENDIF ! node "C"
  ENDDO grloop  

END FUNCTION heat

!=======================================================================
! NAMEBEN is used in "heat" to name benson groups
!=======================================================================
SUBROUTINE nameben(base,point,ic,icdot,io,ico,icodot, &
                   icd,init,ialko,ipero,bengrp)
  IMPLICIT NONE

  INTEGER,INTENT(IN)    :: ic,icdot,io,ico,icodot,icd,init,ialko,ipero
  INTEGER,INTENT(INOUT) :: point
  CHARACTER(LEN=*),INTENT(IN) :: base
  CHARACTER(LEN=*),INTENT(OUT) :: bengrp

  INTEGER :: i

  bengrp = ' '

! form benson groups : start with the base:
  bengrp = base
  bengrp(point:point) = '_'
  point=point+1

  IF (icodot/=0) THEN    !--- first (CO*)
    DO i=1,icodot
      bengrp(point:) = '(CO*)' ; point = point + 5
    ENDDO 
  ENDIF

  IF (ic/=0) THEN        !--- then (C):
    DO i=1,ic
      bengrp(point:) = '(C)' ; point = point + 3
    ENDDO 
  ENDIF

  IF (icdot/=0) THEN     !--- then (C*):
    DO i=1,icdot
      bengrp(point:) = '(C*)' ; point = point + 4
    ENDDO 
  ENDIF

  IF (icd/=0) THEN       !--- then (Cd):
    DO i=1,icd
      bengrp(point:) = '(Cd)' ; point = point + 4
    ENDDO 
  ENDIF

  IF (ico/=0) THEN       !--- then (CO): 
    DO i=1,ico
      bengrp(point:) = '(CO)' ; point = point + 4
    ENDDO 
  ENDIF

  IF (io/=0) THEN        !--- then (O): 
    DO i=1,io
      bengrp(point:) = '(O)' ; point = point + 3
    ENDDO 
  ENDIF

  IF (ipero/=0) THEN     !--- then (OO*):
    DO i=1,ipero
      bengrp(point:) = '(OO*)' ; point = point + 5
    ENDDO 
  ENDIF

  IF (ialko/=0) THEN     !--- then (O*):
    DO i=1,ialko
      bengrp(point:) = '(O*)' ; point = point + 4
    ENDDO 
  ENDIF

  IF (init/=0) THEN      !--- then (ONO2):
    DO i=1,init
      bengrp(point:) = '(ONO2)' ; point = point + 6
    ENDDO 
  ENDIF

END SUBROUTINE nameben

! ======================================================================
! Return the benson value for the group provided as input
! ======================================================================
SUBROUTINE getben(bengrp,bvalue,check)
  USE database, ONLY:lbg,nbson,bsonval,bsongrp
  IMPLICIT NONE

  CHARACTER(LEN=lbg),INTENT(IN) :: bengrp
  REAL,INTENT(OUT)    :: bvalue
  INTEGER,INTENT(OUT) :: check

  INTEGER :: i

  bvalue=0.  ;  check=1
  DO i=1,nbson
    IF (bengrp==bsongrp(i)) THEN
      check=0  ;  bvalue=bsonval(i)
      RETURN
    ENDIF
  ENDDO
END SUBROUTINE getben


END MODULE bensontool
