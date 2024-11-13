SUBROUTINE alphao_criegee(xcri,group,bond,ngr,brch,cnod,noda,nodb, &
                          np,s,p,nref,ref)
  USE fragmenttool, ONLY: fragm  
  USE ringtool, ONLY: ring_data
  IMPLICIT NONE

  CHARACTER(LEN=*),INTENT(IN) :: xcri     ! formula of the hot criegee
  CHARACTER(LEN=*),INTENT(IN) :: group(:) ! group in chem
  INTEGER,INTENT(IN) :: bond(:,:)         ! bond matrix
  INTEGER,INTENT(IN) :: ngr               ! # of group/nod in chem
  REAL,INTENT(IN)    :: brch              ! branching ratio of the hot criegee
!  REAL,INTENT(IN)    :: yield             ! yield of the current pathway
  INTEGER,INTENT(IN) :: cnod              ! criegee node
  INTEGER,INTENT(IN) :: noda              ! node of the 1st branch (next to criegee) 
  INTEGER,INTENT(IN) :: nodb              ! node of the 2nd branch (next to criegee) 
  INTEGER,INTENT(INOUT) :: np             ! # of product in the s and p list
  REAL,INTENT(INOUT)    :: s(:)           ! stoi. coef. of product in p(:) list
  CHARACTER(LEN=*),INTENT(INOUT) :: p(:)  ! product list (short names)
  INTEGER,INTENT(INOUT) :: nref           ! # of references added in the reference list
  CHARACTER(LEN=*),INTENT(INOUT):: ref(:) ! list of references 

  CHARACTER(LEN=LEN(group(1))) :: tgroup(SIZE(group))      
  INTEGER :: tbond(SIZE(bond,1),SIZE(bond,2))
  CHARACTER(LEN=LEN(xcri))     :: fragpd(2), pchem
  INTEGER :: i  
  INTEGER :: ifrag,nc,nring

  INTEGER,PARAMETER :: mxirg=6             ! max # of distinct rings 
  INTEGER  :: ndrg                         ! # of disctinct rings
  INTEGER  :: trackrg(mxirg,SIZE(tgroup))  ! (a,:)== track (node #) belonging ring a
  LOGICAL  :: lorgnod(mxirg,SIZE(tgroup))  ! (a,b)==true if node b belong to ring a

  tgroup(:)=group(:)  ;  tbond=bond(:,:)

! criegee should not be part of ring (no biradicals !)
  CALL ring_data(cnod,ngr,tbond,tgroup,ndrg,lorgnod,trackrg)
  nring=COUNT(lorgnod(:,cnod))    ! count # of rings involving node cnod
  IF (nring/=0) THEN
    STOP "ring in -O-C.(OO.)- structure"    ! ester channel for cyclic species (no diradical!)
  ENDIF
    
! external criegee R-O-CH.(OO.)
! ----------------
  IF (nodb==0) THEN
  
    yield=1. ! assume 100 % decomposition
    IF (group(cnod)(1:8)=='CH.(OO.)') THEN
      CALL add1tonp(progname,xcri,np) ; s(np)=yield ; p(np)='CO2 '  ! add CO2
      CALL add1tonp(progname,xcri,np) ; s(np)=yield ; p(np)='HO2 '  ! add HO2
    ELSE
      STOP 'write group external for R-O-C(OO.)'  
    ENDIF
    tgroup(cnod)=' ' ; tbond(cnod,noda)=0 ; tbond(noda,cnod)=0
  
    IF (tgroup(noda)/='-O-') THEN  
      PRINT*, 'error not criegee'
      STOP
    ENDIF
    DO i=1,ngr
      IF ((tbond(i,noda)/=0).AND.(i/=cnod)) THEN
         tbond(i,noda)=0 ; tbond(noda,i)=0 ; tgroup(noda)=' '
         nc=INDEX(tgroup(i),' ') ; tgroup(i)(nc:nc+3)='(O.)'
      ENDIF
    ENDDO
    CALL rebond(tbond,tgroup,tempfo,nring)
    CALL add2p(xcri,tempfo,yield,brch,np,s,p,nref,ref)
    tgroup(:)=group(:) ; tbond(:,:)=bond(:,:)

! internal criegee R-O-CH.(OO.)
! ----------------
  ELSE
    PRINT*,'inside internal'
  
    yield=1.  ! assume 100 % decomposition
    CALL add1tonp(progname,xcri,np) ; s(np)=yield ; p(np)='CO2 '  ! add CO2
    tgroup(cnod)=' '
  
    tbond(cnod,noda)=0 ; tbond(noda,cnod)=0
    tbond(cnod,nodb)=0 ; tbond(nodb,cnod)=0
    IF (tgroup(noda)=='-O-') THEN  
      DO i=1,ngr
        IF ((tbond(i,noda)/=0).AND.(i/=cnod)) THEN
          tbond(i,noda)=0 ; tbond(noda,i)=0 ; tgroup(noda)=' '
          nc=INDEX(tgroup(i),' ') ; tgroup(i)(nc:nc+3)='(O.)'
        ENDIF
      ENDDO
      nc=INDEX(tgroup(nodb),' ')  ; tgroup(nodb)(nc:nc)='.'
  
    ELSE IF (tgroup(nodb)=='-O-') THEN
      DO i=1,ngr
        IF ((tbond(i,nodb)/=0).AND.(i/=cnod)) THEN
          tbond(i,nodb)=0 ; tbond(nodb,i)=0 ; tgroup(nodb)=' '
          nc=INDEX(tgroup(i),' ') ; tgroup(i)(nc:nc+3)='(O.)'
        ENDIF
      ENDDO
      nc=INDEX(tgroup(noda),' ')  ; tgroup(noda)(nc:nc)='.'
  
    ELSE
      STOP 'error no -O- in internal criegee'
    ENDIF
    PRINT*,'before frgm'
  
    CALL fragm(tbond,tgroup,fragpd(1),fragpd(2))
    PRINT*,'fragment 1:',fragpd(1)
    PRINT*,'fragment 2:',fragpd(2)

! loop over the 2 fragmentation products
    DO ifrag=1,2
      pchem=fragpd(ifrag)
      CALL add2p(xcri,pchem,yield,brch,np,s,p,nref,ref)
    ENDDO
    tgroup(:)=group(:) ; tbond(:,:)=bond(:,:)
    
    DO i=1,10
      WRITE(*,'(i3,f6.3,2x,a)') i,s(i),p(i)
    ENDDO
    DO i=1,nref
      WRITE(*,*) ref(i)
    ENDDO
  ENDIF

END SUBROUTINE alphaocriegee
