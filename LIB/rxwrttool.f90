MODULE rxwrttool
IMPLICIT NONE

! various reaction counters (set to 0 the 1st call)
INTEGER, SAVE :: nrx_all=0   ! counter for all reactions 
INTEGER, SAVE :: nrx_n=0     ! counter for thermal reaction
INTEGER, SAVE :: nrx_hv=0    ! counter for hv reaction
INTEGER, SAVE :: nrx_extra=0 ! counter for extra reaction
INTEGER, SAVE :: nrx_fo=0    ! counter for fall-off eaction
INTEGER, SAVE :: nrx_tb=0    ! counter for third body reaction
INTEGER, SAVE :: nrx_o2=0    ! counter for oxygen reaction
INTEGER, SAVE :: nrx_meo2=0  ! counter for CH3O2 reaction
INTEGER, SAVE :: nrx_ro2=0   ! counter for RO2 reaction
INTEGER, SAVE :: nrx_isom=0  ! counter for isomerization reaction
INTEGER, SAVE :: nrx_tabcf=0 ! counter for tabulated stoi. coef. (not used)
INTEGER, SAVE :: nrx_ain=0   ! counter for gas   -> part. reaction
INTEGER, SAVE :: nrx_aou=0   ! counter for part. -> gas   reaction
INTEGER, SAVE :: nrx_win=0   ! counter for gas  -> wall reaction
INTEGER, SAVE :: nrx_wou=0   ! counter for wall -> gas  reaction

CONTAINS
! SUBROUTINE rxinit(r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe)
! SUBROUTINE rxwrit(lout,r,s,p,arrh,idreac,nlab,xlabel,folow,fotroe,phase)
! SUBROUTINE rxwrit_dyn(lout,r,s,p,arrh,idreac,auxinfo,charfrom,charto)
! SUBROUTINE count4rxn(idrx)

!=======================================================================
! PURPOSE:  -  Initialize the reaction line (see also rxwrit).    
!=======================================================================
SUBROUTINE rxinit(r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe)
  IMPLICIT NONE

  CHARACTER(LEN=*),INTENT(OUT) :: r(:) ! array of reagents
  CHARACTER(LEN=*),INTENT(OUT) :: p(:) ! array of products
  INTEGER,INTENT(OUT) :: idreac        ! ID for the reaction type (thermal, HV ...)
  INTEGER,INTENT(OUT) :: nlabel        ! label if the reaction is not thermal
  REAL,INTENT(OUT)    :: s(:)          ! array of stoichiometry coefficients
  REAL,INTENT(OUT)    :: arrh(:)       ! arrhenius coefficients
  REAL,INTENT(OUT)    :: xlabel        ! weighting factor for HV reaction
  REAL,INTENT(OUT)    :: folow(:)      ! low pressure fall off arrhenius coefficient
  REAL,INTENT(OUT)    :: fotroe(:)     ! troe parameter for fall off reaction

  r(:)=' '   ; s(:)=0.  ; p(:)=' ' ; folow(:)=0. ; fotroe(:)=0.
  arrh(:)=0. ; idreac=0 ; nlabel=0 ; xlabel=0.
  
END SUBROUTINE rxinit

!=======================================================================
! PURPOSE: write the reaction to a file in a standard format.                                                                 
! This subroutine computes the number of required splits (npmax 
! products per line), adjusts stoichiometry coefficients and rate
! constants accordingly, collapses identical products and blanks out
! unity stoichiometry coefficients. Assume that the reaction to be
! written is a gas phase reaction (default) - if "phase" is provided
! as input, then the reaction is written for the phase provided (i.e. 
! using the 1st character for the considered phase. References are 
! optional and written in the output if the corresponding flag is 
! raised.
!
! NOTE: if "iflost" in "tempflag" is turned on, then no reaction 
! products are written (just add carbon "loss" for C balance).                                    
!                                
! ID number are :                
!    0 => simple thermal reaction
!    1 => HV reaction            
!    2 => EXTRA reaction         
!    3 => Fall Off reaction      
!=======================================================================
SUBROUTINE rxwrit(lout,r,s,p,arrh,idreac,nlab,xlabel,folow,fotroe,com,phase)
  USE keyparameter, ONLY: refu,mxnp,mxlco 
  USE keyflag, ONLY: wrtref  ! write reference flag
  USE references, ONLY: mxlcod
  USE tempflag, ONLY: iflost,xxc     ! rm rx products if "iflost" set to 1
  USE toolbox, ONLY: addref
  IMPLICIT NONE  

  CHARACTER(LEN=*),INTENT(IN) :: r(:)   ! array of reagents
  CHARACTER(LEN=*),INTENT(IN) :: p(:)   ! array of products
  REAL,INTENT(IN)    :: s(:)        ! array of stoi. coef.
  INTEGER,INTENT(IN) :: lout        ! unit number for the mechanism file
  INTEGER,INTENT(IN) :: idreac      ! ID for the reaction type (thermal, HV ...)
  INTEGER,INTENT(IN) :: nlab        ! label if the reaction is not thermal 
  REAL,INTENT(IN)    :: arrh(:)     ! arrhenius coefficient
  REAL,INTENT(IN)    :: xlabel      ! weighting factor for HV reaction
  REAL,INTENT(IN)    :: folow(:)    ! low pressure fall off arrhenius coef.
  REAL,INTENT(IN)    :: fotroe(:)   ! troe parameter for fall off reaction
  CHARACTER(LEN=*),INTENT(IN),OPTIONAL :: com(:) ! references for the reaction 
  CHARACTER(LEN=*),INTENT(IN),OPTIONAL :: phase  ! 1 optional char for phase (default is gas) 
  
  CHARACTER(LEN=6)   :: charstoi(SIZE(s))
  CHARACTER(LEN=1)   :: c1,c2,signc(SIZE(s))
  CHARACTER(LEN=1)   :: rphase
  CHARACTER(LEN=LEN(p(1))) :: p1(SIZE(s)), p2(SIZE(s))
  REAL          :: s1(SIZE(s)), s2(SIZE(s))
  INTEGER       :: nsplit
  REAL          :: rc1, xlab
  INTEGER       :: i,j,k,np,j1,j2,jj,mnp
  CHARACTER(LEN=mxlco+1) :: rg(SIZE(r)), pg(SIZE(s))
  LOGICAL       :: locheck
  CHARACTER(LEN=200) :: refline
  INTEGER       :: ipos,idrx
  INTEGER       :: nref,nref1
  CHARACTER(LEN=mxlcod), ALLOCATABLE :: ref(:),ref1(:)
  CHARACTER(LEN=7),PARAMETER :: progname='rxwrit'

  mnp=SIZE(s)                          ! max # of products in a reaction
  IF (PRESENT(phase)) THEN
    IF (LEN(phase)/=1) THEN
      PRINT*, "in rxwrit,unexpected phase string" 
      STOP "in rxwrit" 
    ENDIF
  ENDIF
! phase default = gas phase
  IF (PRESENT(phase)) THEN; rphase=phase ; ELSE ; rphase='G' ; ENDIF 

! handle references whether provided or not as input (size of com is unknow)
  IF (PRESENT(com)) THEN
    ALLOCATE(ref(SIZE(com))) ; ref(:)=com(:) ; nref=COUNT(ref/=' ')
  ELSE
    ALLOCATE(ref(1)) ; ref(:)=' ' ; nref=0
  ENDIF
  ALLOCATE(ref1(SIZE(ref))) ; ref1(:)=' '  ! cp of ref without duplicate

  
! initialize & duplicate variables which may change (i.e. preserve output) 
  rc1 = arrh(1)
  signc(:) = ' ' ;  charstoi(:) = ' '  ;   pg(:) = ' '
  s1(:) = s(:)   ;  p1(:) = p(:)       ;   xlab = xlabel
  s2(:) = 0.     ;  p2(:) = ' '

! prepare the reaction string (add keyword, sign, etc ...)
! ------------------------------------------------------- 

! remove the product if flag to stop chemistry is raised
  IF (iflost==1) THEN
    s1(:) = 0.         ;  p1(:)=' '
    s1(1) = REAL(xxc)  ;  p1(1) = 'XCLOST'
  ENDIF

! clean table of products and stoi coefs (just in case)
  WHERE (s1(:)==0.)  p1(:) = ' '
  WHERE (p1(:)==' ') s1(:) = 0.

! collapse identical products:
  DO i =1,mnp-1
    IF (p1(i)(1:1)/=' ') THEN
      DO j=i+1,mnp
        IF (p1(j)==p1(i)) THEN
          s1(i) = s1(i) + s1(j)  ;  s1(j) = 0. ;  p1(j) = ' '
        ENDIF
      ENDDO
    ENDIF
  ENDDO

! rm blank and count distinct products. Data are stored in s2 and p2.
  np = 0
  DO i=1,mnp
    IF (p1(i)/=' ') THEN
      np = np+1  ;  s2(np) = s1(i)  ;  p2(np) = p1(i)      
    ENDIF
  ENDDO

! transcribe stoi. coef. to character string, blanking out unity values
  DO i=1,np
    IF (s2(i)==0.) CYCLE
    IF (s2(i)==1.) CYCLE
    WRITE(charstoi(i),'(F6.3)') ABS(s2(i))
  ENDDO

! add "G" (for gas phase) or any other phase to each name, expect keywords
  rg(:)=' ' ; idrx=0
  DO i=1,SIZE(r)    ! reactants
    IF (r(i)(1:1)/=' ')  THEN
      IF      (r(i)(1:3)=='HV ')    THEN ; idrx=1 ; rg(i)='HV'
      ELSE IF (r(i)(1:6)=='EXTRA ') THEN ; idrx=2 ; rg(i)='EXTRA'
      ELSE IF (r(i)(1:4)=='(+M)')   THEN ; idrx=3 ; rg(i)='FALLOFF'
      ELSE IF (r(i)(1:2)=='M ')     THEN ; idrx=4 ; rg(i)='TBODY' 
      ELSE IF (r(i)(1:6)=='OXYGEN') THEN ; idrx=5 ; rg(i)='OXYGEN'
      ELSE IF (r(i)(1:4)=='PERO' )  THEN ; idrx=6 ; rg(i)=r(i)      ! warning pero1, ...
      ELSE IF (r(i)(1:6)=='MEPERO') THEN ; idrx=7 ; rg(i)='MEPERO'
      ELSE IF (r(i)(1:6)=='ISOM  ') THEN ; idrx=8 ; rg(i)='ISOM'
      ELSE IF (r(i)(1:6)=='TABCF ') THEN ; idrx=9 ; rg(i)='TABCF'   ! for tabulated stoi. coef.
      ELSE 
        rg(i)(1:1)=rphase  ;  rg(i)(2:7)=r(i)
      ENDIF  
    ENDIF  
  ENDDO

  DO i=1,np    ! products
    IF (p2(i)/=' ') THEN
      locheck=.FALSE.
      IF (p2(i)(1:4)=='(+M)') THEN
        locheck=.TRUE.
        p2(i)=' '
      ENDIF  
      IF (p2(i)=='EMPTY ') locheck=.TRUE.
      IF (locheck) THEN
        pg(i)(1:6)=p2(i)
        IF (p2(i)=='EMPTY ') pg(i)='NOTHING'
      ELSE
        pg(i)(1:1)=rphase
        pg(i)(2:7)=p2(i)
      ENDIF
    ENDIF
  ENDDO

! add the sign '+' between reactants
  c1=' '  ;  c2=' '
  IF (rg(2)/=' ')  c1='+'
  IF (rg(3)/=' ')  c2='+'

! add the sign '+' between products
  DO i=2,np
    j=i-1 
    IF (p2(i)/=' ') THEN
      IF (s2(i) < 0.) THEN ; signc(j)='-'
      ELSE                 ; signc(j)='+'
      ENDIF
    ENDIF 
  ENDDO

! write comment line (reac_info codes)
! --------------------------------
  IF (wrtref) THEN
    DO i=1,nref-1  ! kill duplicate reac_infos
      DO j=i+1,nref
        IF (ref(i)==ref(j)) ref(j)=' ' 
      ENDDO
    ENDDO
    nref1=0        ! make a copy of ref without duplicate
    DO i=1,nref
      IF (ref(i)==' ') CYCLE
      nref1=nref1 + 1 ; ref1(nref1)=ref(i)
    ENDDO
    
    refline=' '
    DO i=1,nref1
      ipos=13*(i-1)+1
      refline(ipos:)=" ; "//ref1(i)
    ENDDO
  ENDIF
  
! write reaction
! --------------
  CALL count4rxn(idrx)   ! update reactions counters 

  ! get # "line" required to write the reaction (mxnp pdct per line) 
  nsplit=1
  IF (np > mxnp)  nsplit = (np-1)/mxnp + 1

  ! reaction stand in 1 line
  IF (nsplit==1) THEN
    j1 = 1  ;  j2 = 4
    WRITE(lout,'(3(A7,1x,A1,1x),A2,4(A6,1X,A7,1X,A1),4X,ES10.3,1X,f4.1,1X,f7.0)') &
      rg(1),c1,rg(2),c2,rg(3),' ','=>',&
      (charstoi(jj),pg(jj),signc(jj),jj=j1,j2),rc1,arrh(2:3) 
    IF (wrtref) THEN
      WRITE(refu,'(3(A7,1x,A1,1x),A2,4(A6,1X,A7,1X,A1),4X,ES10.3,1X,0P,f4.1,1X,f7.0,2a)') &
      rg(1),c1,rg(2),c2,rg(3),' ','=>',&
      (charstoi(jj),pg(jj),signc(jj),jj=j1,j2),rc1,arrh(2:3),' ',TRIM(refline) 
    ENDIF

  ! reaction written using many lines
  ELSE
    DO k=1,nsplit
      j1 = 1 + mxnp*(k-1)
      j2 = j1 + mxnp - 1

      ! 1st line
      IF (k==1) THEN
        WRITE(lout,'(3(A7,1x,A1,1x),A2,4(A6,1X,A7,1X,A1))') &
          rg(1),c1,rg(2),c2,rg(3),' ','=>',(charstoi(jj),pg(jj),signc(jj),jj=j1,j2) 
        IF (wrtref) THEN
          WRITE(refu,'(3(A7,1x,A1,1x),A2,4(A6,1X,A7,1X,A1))') &
            rg(1),c1,rg(2),c2,rg(3),' ','=>',(charstoi(jj),pg(jj),signc(jj),jj=j1,j2) 
        ENDIF

      ! last line
      ELSE IF (k==nsplit) THEN
        WRITE(lout,'(32X,4(A6,1X,A7,1X,A1),4X,ES10.3,1X,f4.1,1X,f7.0)') &
          (charstoi(jj),pg(jj),signc(jj),jj=j1,j2),rc1,arrh(2:3) 
        IF (wrtref) THEN
          WRITE(refu,'(32X,4(A6,1X,A7,1X,A1),4X,ES10.3,1X,f4.1,1X,f7.0,2a)') &
            (charstoi(jj),pg(jj),signc(jj),jj=j1,j2),rc1,arrh(2:3),' ',TRIM(refline) 
        ENDIF

      ! middle line
      ELSE
        WRITE(lout,'(32X,4(A6,1X,A7,1X,A1))') &
          (charstoi(jj),pg(jj),signc(jj),jj=j1,j2)
        IF (wrtref) THEN
          WRITE(refu,'(32X,4(A6,1X,A7,1X,A1))') &
             (charstoi(jj),pg(jj),signc(jj),jj=j1,j2)
        ENDIF

      ENDIF
    ENDDO
  ENDIF

! write keyword and label if necessary
  IF (idreac>0) THEN
    IF (idreac==1) THEN
      WRITE(lout,'(A7,I5,2x,f6.3,A2)')             "  HV / ",nlab,xlab," /"
      IF (wrtref) WRITE(refu,'(A7,I5,2x,f6.3,A2)') "  HV / ",nlab,xlab," /"
  
    ELSE IF (idreac==3) THEN
      WRITE(lout,'(A11,E9.3,F5.1,F7.0,1X,4(F6.1,1x),a1)') &
           "  FALLOFF /",folow(1:3),fotroe(1:4),"/"
      IF (wrtref) WRITE(refu,'(A11,1PE9.2,0P,F5.1,F7.0,1X,4(F6.1,1x),a1)') &
           "  FALLOFF /",folow(1:3),fotroe(1:4),"/"
  
    ELSE IF (idreac==2) THEN
      WRITE(lout,'(A10,I5,A2)')              "  EXTRA / ",nlab," /"
      IF (wrtref) WRITE(refu,'(A10,I5,A2)')  "  EXTRA / ",nlab," /"
      
    ELSE
      WRITE(6,*) '--error--, in rxwrit. Idreac is > 3, in reaction :'
      WRITE(6,'(a)') rg(1),c1,rg(2),c2,rg(3),' =>'
      STOP "in rxwrit"
    ENDIF
  ENDIF
  DEALLOCATE(ref) ; DEALLOCATE(ref1)
            
END SUBROUTINE rxwrit

!=======================================================================
! PURPOSE: write the pseudo reaction describing mass transfer between
! phases (gas <-> particles or gas <-> wall).
!
! ID number are :
!   1 => AIN (aerosol in : gas -> particle)
!   2 => AOU (aerosol out: particle -> gas)
!   3 => WIN (wall in : gas -> wall)
!   4 => WOU (wall out: wall -> gas)
!=======================================================================
SUBROUTINE rxwrit_dyn(lout,r,s,p,arrh,idreac,auxinfo,charfrom,charto)
  USE keyparameter, ONLY: mxnp  ! max # of products per reaction
  IMPLICIT NONE  

  CHARACTER(LEN=*),INTENT(IN) :: r(:)   ! array of reagents
  CHARACTER(LEN=*),INTENT(IN) :: p(:)   ! array of products
  REAL,INTENT(IN)    :: s(:)         ! array of stoi. coef.
  REAL,INTENT(IN)    :: arrh(:)      ! arrhenius coefficient
  INTEGER,INTENT(IN) :: lout         ! unit number for the mechanism file
  INTEGER,INTENT(IN) :: idreac       ! ID for the reaction type (AIN, AOU, ...)
  REAL,INTENT(IN)    :: auxinfo(9)   ! only 2 needed here
  CHARACTER(LEN=*),INTENT(IN) :: charfrom  ! char "phase" of "reagant"
  CHARACTER(LEN=*),INTENT(IN) :: charto    ! char "phase" of "product"

  CHARACTER*1   ::  c1,c2,signc(SIZE(s))
  CHARACTER*(6) ::  charstoi(SIZE(s))
  INTEGER       ::  i,j1,j2,jj,mnp
  CHARACTER*(7) ::  rg(SIZE(r)),pg(SIZE(s))
  LOGICAL       ::  locheck

  mnp=SIZE(s)                          ! max # of products in a reaction 
      
! check that species is the same both side
  IF (r(1) /= p(1)) THEN
    WRITE(6,*) '--error--, in rxwrit_dyn. distinct reactant and product:'
    WRITE(6,'(a)') r(1), '=>', p(1)
    STOP "in rxwrit_dyn"
  ENDIF

! check that only p(1) is filled with s(1)=1
  IF (s(1)/=1.) THEN
    WRITE(6,*) '--error--, in rxwrit_dyn. Expect stoe. coef. equal 1'
    WRITE(6,'(a)') s(1)
    STOP "in rxwrit_dyn"
  ENDIF
  IF (p(1)(1:2)=='  ') THEN
    WRITE(6,*) '--error--, in rxwrit_dyn. No product found'
    WRITE(6,'(a)') p(1)
    STOP "in rxwrit_dyn"
  ENDIF
  DO i=2,mnp
    IF (p(i)(1:2)/=' ') THEN
      WRITE(6,*) '--error--, in rxwrit_dyn. More than one product found'
      WRITE(6,'(a)') p(i)
      STOP "in rxwrit_dyn"
    ENDIF
  ENDDO

! idreac must be between 1 and 4
  IF ((idreac < 1).OR.(idreac > 5)) THEN
    WRITE(6,*) '--error--, in rxwrit_dyn. Idreac out of bound'
    WRITE(6,'(a)') p(i)
    STOP "in rxwrit_dyn"
  ENDIF

! check that second is one of the expected keyword
  locheck=.FALSE.
  IF (r(2)(1:4)=='AIN ' ) locheck=.TRUE.
  IF (r(2)(1:4)=='AOU ' ) locheck=.TRUE.
  IF (r(2)(1:4)=='WIN ' ) locheck=.TRUE.
  IF (r(2)(1:4)=='WOU ' ) locheck=.TRUE.
  IF (.not.locheck) THEN
    WRITE(6,*) '--error--, in rxwrit_dyn. Unexpected keyword :'
    WRITE(6,'(a)') r(2)
    STOP "in rxwrit_dyn"
  ENDIF

! write reaction using the format in rxwrit.f90
! --------------------------------------------
  charstoi(:)=' ' ! no coef, since should be 1 or 0
  signc(:)=' '  ; pg(:) = ' '    ; rg(:)= ' '

  rg(1)(1:1)=charfrom  ;  rg(1)(2:7)=r(1)
  rg(2)(1:)=r(2)(1:)  
  pg(1)(1:1)=charto    ;  pg(1)(2:7)=p(1)

  j1 = 1  ; j2 = 4
  c1='+'  ; c2=' '

  CALL count4rxn(20+idreac)  ! count number of reaction type
  WRITE(lout,120)  rg(1),c1,rg(2),c2,rg(3),' ','=>', &
       (charstoi(jj),pg(jj),signc(jj),jj=j1,j2),arrh(1),arrh(2),arrh(3)

120   FORMAT(3(A7,A1),A2,4(A5,1X,A7,1X,A1),4X,ES10.3,1X,f4.1,1X,f7.0)

! write keyword and label if necessary
! -------------------------------------
  IF (idreac==4) THEN   
      WRITE(lout,'(A7,1P,E10.2,0P,A1)')  '  WOU/ ',auxinfo(1),'/'
  ELSE IF (idreac==3) THEN 
     WRITE(lout,'(A7,1P,E10.2,0P,A1)')   '  WIN/ ',auxinfo(1),'/'
  ENDIF

END SUBROUTINE rxwrit_dyn

!=======================================================================
! PURPOSE: count reaction type. Not used to generate the mechanisms but
! provided as output. Numbers might be needed to allocate size of 
! tables in the boxmodel.
!=======================================================================
SUBROUTINE count4rxn(idrx)
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: idrx         ! reaction type
    
  nrx_all=nrx_all+1
  SELECT CASE (idrx)
  CASE(0)  ; nrx_n=nrx_n+1         ! thermal reaction
  CASE(1)  ; nrx_hv=nrx_hv+1       ! HV reaction 
  CASE(2)  ; nrx_extra=nrx_extra+1 ! EXTRA reaction 
  CASE(3)  ; nrx_fo=nrx_fo+1       ! FALLOFF reaction
  CASE(4)  ; nrx_tb=nrx_tb+1       ! TBODY reaction
  CASE(5)  ; nrx_o2=nrx_o2+1       ! OXYGEN reaction
  CASE(6)  ; nrx_ro2 =nrx_ro2+1    ! RO2 reaction
  CASE(7)  ; nrx_meo2=nrx_meo2+1   ! CH3O2 reaction
  CASE(8)  ; nrx_isom=nrx_isom+1   ! ISOMerization reaction
  CASE(9)  ; nrx_tabcf=nrx_tabcf+1 ! tabulated stoi. coef. (not used)
  CASE(21) ; nrx_ain=nrx_ain+1     ! gas   -> part. reaction
  CASE(22) ; nrx_aou=nrx_aou+1     ! part. -> gas   reaction
  CASE(23) ; nrx_win=nrx_win+1     ! gas  -> wall reaction
  CASE(24) ; nrx_wou=nrx_wou+1     ! wall -> gas  reaction
  CASE DEFAULT 
     STOP "reaction ID not identified"
  END SELECT

END SUBROUTINE count4rxn

END MODULE rxwrttool
