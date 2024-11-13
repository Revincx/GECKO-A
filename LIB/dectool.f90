MODULE dectool
IMPLICIT NONE
CONTAINS
!SUBROUTINE coono2_dec(chem,copdec,nca,nrf,rf)

!=======================================================================
! PURPOSE: For RCO(ONO2), assume decomposition to R. + CO2 + NO2
!=======================================================================
SUBROUTINE coono2_dec(chem,copdec,nca,nrf,rf)
  USE keyparameter, ONLY: mxnode,mxlgr
  USE reactool, ONLY: rebond
  USE radchktool, ONLY: radchk
  USE stdgrbond, ONLY: grbond  
  USE toolbox, ONLY: stoperr,addref
  IMPLICIT NONE
 
  CHARACTER(LEN=*),INTENT(INOUT) :: chem
  CHARACTER(LEN=*),INTENT(INOUT) :: copdec(:)
  INTEGER,INTENT(IN)             :: nca     ! # of nodes
  INTEGER,INTENT(INOUT)          :: nrf          ! # of reference added for each reactions
  CHARACTER(LEN=*),INTENT(INOUT) :: rf(:)        ! ref/comments for each reactions

  INTEGER                        :: tbond(mxnode,mxnode)
  CHARACTER(LEN=mxlgr)           :: tgroup(mxnode)
  INTEGER                        :: dbf,tpnring,nc,j,k,l,m,nip,idmx
  CHARACTER(LEN=LEN(chem))       :: tpchem
  
  INTEGER,PARAMETER              :: mxrpd=2  ! max # of products returned by radchk sub
  CHARACTER(LEN=LEN(chem))       :: rdckpd(mxrpd)
  CHARACTER(LEN=LEN(copdec(1)))  :: rdckcopd(mxrpd,SIZE(copdec,1))
  REAL                           :: sc(mxrpd)

  CHARACTER(LEN=10),PARAMETER    :: progname='coono2_dec'
  CHARACTER(LEN=70)              :: mesg

  IF (INDEX(chem,'.')/=0) RETURN
  
  IF (INDEX(chem,'CO(ONO2)')/=0) THEN
    CALL grbond(chem,tgroup,tbond,dbf,tpnring)
    kloop: DO k=1,nca
      IF (INDEX(tgroup(k),'CO(ONO2)')/=0) THEN
        DO l=1,nca
          IF (tbond(l,k)==1) THEN
            tbond(l,k)=0  ;  tbond(k,l)=0
            nc=INDEX(tgroup(l),' ')  ;  tgroup(l)(nc:nc)='.'

            tgroup(k)=' '     ! erase CO(ONO2) grp
            CALL rebond(tbond,tgroup,tpchem,tpnring)
            CALL radchk(tpchem,rdckpd,rdckcopd,nip,sc,nrf,rf)
            idmx=MAXLOC(sc,DIM=1) ! keep only the species with max yield
            chem=rdckpd(idmx)     ! species standardized in radchk

            DO m=1,SIZE(copdec,1) 
              IF (copdec(m)==' ') EXIT   ! 1st available slot
            ENDDO  
            DO j=1,SIZE(rdckcopd,2)           ! cp the copdct
              IF (rdckcopd(idmx,j)==' ') EXIT 
              IF (rdckcopd(idmx,j)==' ') EXIT 
              copdec(m)=rdckcopd(idmx,j) ; m=m+1
            ENDDO
            copdec(m)='CO2'  ;  copdec(m+1)='NO2'  ! add CO2+NO2

            IF (m+1 > SIZE(copdec,1)) THEN
              mesg="m > SIZE(copdec,1) in decomposition "
              CALL stoperr(progname,mesg,chem)
            ENDIF         
            EXIT kloop
            
          ENDIF
        ENDDO
      ENDIF
    ENDDO kloop
    CALL addref(progname,'RCO(NO3)',nrf,rf)

  ELSE IF (INDEX(chem,'CHO(ONO2)')/=0) THEN
    chem='CO2'
    mloop: DO m=1,SIZE(copdec,1)
      IF (copdec(m)(1:1)==' ')THEN
        copdec(m)='NO2'  ;  copdec(m+1)='HO2'
        EXIT mloop
      ENDIF
    ENDDO  mloop
    CALL addref(progname,'CHO(ONO2)',nrf,rf)
  ENDIF

END SUBROUTINE

END MODULE
