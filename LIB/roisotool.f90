MODULE roisotool
IMPLICIT NONE
CONTAINS
!SUBROUTINE roiso(bond,group,nca,ia,nring,rjg,&
!           niso,fliso,piso,copiso,fliso2,piso2,copiso2,siso, &
!           arrhiso,FSD)
!SUBROUTINE roisom_ver(bond,group,nca,ia,ig,span,arrh,FSD)
!SUBROUTINE roisom_atk(tbond,tgroup,nca,ig,arrh)

!=======================================================================
! PURPOSE: perform the H-shift isomerisation reaction for R(O.) 
! radicals (products is output tables) and compute the reaction rates 
! (arrhiso, FSD tables).
! NOTE: reaction rate constant can be estimated with:
! - Atkinson, 2007 (Atmospheric Environment)
! - Vereecken and Peeters (2010) 
!=======================================================================
SUBROUTINE roiso(bond,group,nca,ia,nring,rjg,&
           niso,fliso,piso,copiso,fliso2,piso2,copiso2,siso, &
           arrhiso,FSD,nrfiso,rfiso)
  USE ringtool, ONLY: findring 
  USE reactool, ONLY: swap, rebond
  USE radchktool, ONLY: radchk    
  USE keyflag, ONLY: kisom_sar         ! select SAR for reaction rate  
  USE toolbox, ONLY: addref
  IMPLICIT NONE

  INTEGER,INTENT(IN)  :: bond(:,:)
  CHARACTER(LEN=*),INTENT(IN) :: group(:)
  INTEGER,INTENT(IN)  :: nca         ! # of nodes
  INTEGER,INTENT(IN)  :: ia          ! index for node bearing (O.) group
  INTEGER,INTENT(IN)  :: nring
  INTEGER,INTENT(IN)  :: rjg(:,:)
  INTEGER,INTENT(OUT) :: niso        ! # of isom reactions
  INTEGER,INTENT(OUT) :: fliso(:)    ! isom reaction flag  
  CHARACTER(LEN=*),INTENT(OUT) :: piso(:)     ! products of the isomerization
  CHARACTER(LEN=*),INTENT(OUT) :: copiso(:,:) ! co-products of the reaction
  INTEGER,INTENT(OUT) :: fliso2(:)            ! isom reaction flag (2nd set) 
  CHARACTER(LEN=*),INTENT(OUT) :: piso2(:)    ! products of the isomerization (2nd set)
  CHARACTER(LEN=*),INTENT(OUT) :: copiso2(:,:)! co-products of reaction (2nd set)
  REAL,INTENT(OUT)    :: siso(:,:)            ! stoe. coef. 
  REAL,INTENT(OUT)    :: arrhiso(:,:)         ! arrhenius coef.
  REAL,INTENT(OUT)    :: FSD(:,:)             ! coef. for polynomial fct. (Vereecken SAR)
  INTEGER,INTENT(INOUT) :: nrfiso(:)          ! # of reference added for each reactions
  CHARACTER(LEN=*),INTENT(INOUT):: rfiso(:,:) ! ref/comments  

  INTEGER :: tbond(SIZE(bond,1),SIZE(bond,1))
  CHARACTER(LEN=LEN(group(1))) :: tgroup(SIZE(group))
  CHARACTER(LEN=LEN(group(1))) :: pold, pnew
  CHARACTER(LEN=LEN(piso(1)))  :: tpchem
  CHARACTER(LEN=LEN(piso(1)))     :: rdckpd(SIZE(piso))
  CHARACTER(LEN=LEN(copiso(1,1))) :: rdckcopd(SIZE(copiso,1),SIZE(copiso,2))
  INTEGER :: ring(SIZE(bond,1)), temprg(SIZE(bond,1))
  INTEGER :: i,j,k,rngflg,tpnring,nip,nc,begrg,endrg
  REAL    :: sc(SIZE(rdckpd))

  CHARACTER(LEN=6),PARAMETER :: progname='roiso '

  niso=0       ;  siso(:,:)=0.
  fliso(:)=0   ;  piso(:)=' '    ; copiso(:,:)=' '
  fliso2(:)=0  ;  piso2(:)=' '   ; copiso2(:,:)=' '
  ring(:)=0    ; arrhiso(:,:)=0. ;  FSD(:,:)=0. ;
  tbond(:,:)=bond(:,:) ; tgroup(:)=group(:)  
  
! find nodes that are already part of any ring
  IF (nring>0) THEN
    DO i=1,nring
      begrg=rjg(i,1)  ;  endrg=rjg(i,2)  ;  temprg(:)=0
      CALL findring(begrg,endrg,nca,bond,rngflg,temprg)
      WHERE (temprg(:) /= 0) ring(:)=1
    ENDDO
  ENDIF

! look down-chain (if it is not already on a ring). Note: isomerisation 
! rate for 5 membered ring are expected to be to slow to have a 
! significant contribution and is not considered anymore below.
  aloop: DO i=1,nca
    IF (bond(ia,i)==0 .OR. ring(i)>0) CYCLE aloop
    bloop: DO j=1,nca
      IF ( (bond(i,j)==0) .OR. (j==ia) ) CYCLE bloop ! avoid step-back
      cloop: DO k=1,nca                              ! look for 6-membered rings
        IF ( (bond(k,j)==0) .OR. (k==i) ) CYCLE cloop

        IF (INDEX(group(k),'CH')/=0) THEN
          niso=niso+1
          IF (niso>SIZE(fliso)) THEN
            WRITE(6,'(a)') '--error--, in ro: niso > SIZE(fliso).'
            STOP "in ro - niso > SIZE(fliso)"
          ENDIF
          fliso(niso)=1
          IF (kisom_sar==2) THEN
            CALL roisom_ver(tbond,tgroup,nca,ia,k,5,arrhiso(niso,:),FSD(niso,:))
            CALL addref(progname,'LV10KER000',nrfiso(niso),rfiso(niso,:))
          ELSE
            CALL roisom_atk(tbond,tgroup,nca,k,arrhiso(niso,:))
            CALL addref(progname,'RA07KER000',nrfiso(niso),rfiso(niso,:))
          ENDIF

          IF      (group(k)(1:3)=='CH3') THEN ; pold='CH3' ; pnew='CH2'
          ELSE IF (group(k)(1:3)=='CH2') THEN ; pold='CH2' ; pnew='CH'
          ELSE IF (group(k)(1:3)=='CHO') THEN ; pold='CHO' ; pnew='CO'
          ELSE IF (group(k)(1:2)=='CH')  THEN ; pold='CH'  ; pnew='C'
          ELSE
            WRITE(6,'(a)') '--error-- in roisotool.f90. Group not allowed for isomerization:'
            WRITE(6,'(a)') tgroup(k)
            STOP "in roisotool"
          ENDIF
          CALL swap(group(k),pold,tgroup(k),pnew)
          nc=INDEX(tgroup(k),' ')  ;  tgroup(k)(nc:nc)='.'

          pold='(O.)'  ;   pnew='(OH)' 
          CALL swap(group(ia),pold,tgroup(ia),pnew)
          CALL rebond(tbond,tgroup,tpchem,tpnring)
          CALL radchk(tpchem,rdckpd,rdckcopd,nip,sc,nrfiso(niso),rfiso(niso,:))
          piso(niso)=rdckpd(1)  ;  siso(niso,1)=sc(1)
          copiso(niso,:)=rdckcopd(1,:)
          IF (nip==2) THEN             ! 2nd species 
            fliso2(niso)=1 
            piso2(niso)=rdckpd(2) ;  siso(niso,2)=sc(2)
            copiso2(niso,:)=rdckcopd(2,:)
          ENDIF

          tgroup(:) =group(:)   ! reset
        ENDIF
         
      ENDDO cloop
    ENDDO bloop
  ENDDO aloop
  
! check for duplicate reactions
  DO i=1,niso-1
    iflloop: DO j=i+1,niso
      IF (piso(j)==piso(i)) THEN
        IF (fliso(j)==0) CYCLE iflloop
        DO k=1,SIZE(copiso,2)
            IF (copiso(i,k)/=copiso(j,k)) CYCLE iflloop
        ENDDO
        IF (arrhiso(i,1)/=arrhiso(j,1)) CYCLE iflloop
        IF (arrhiso(i,3)/=arrhiso(j,3)) CYCLE iflloop
        fliso(j)=0  ;  fliso(i)=fliso(i)+1      ! identical reactions
      ENDIF
    ENDDO iflloop
  ENDDO
  DO i=1,niso
    IF (fliso(i)>1) THEN
      arrhiso(i,1)=REAL(fliso(i))*arrhiso(i,1)  ;  fliso(i)=1
    ENDIF
  ENDDO
  
END SUBROUTINE roiso

!=======================================================================
! PURPOSE: compute rate constants for alkoxy radical isomerization
! based on Vereecken and Peeters 2010.
! span: is the span of the migration, i.e. : 4 for a 1-4 migration,                   
! 5 for a 1-5 migration, etc ... 
! Fspan,tunnel=a_span,tunnel * exp(b_span_tunnel/T)
! Fsubs,tunnel=a_subs,tunnel * exp(b_subs_tunnel/T)
!=======================================================================
SUBROUTINE roisom_ver(bond,group,nca,ia,ig,span,arrh,FSD)
  USE keyparameter, ONLY: mxcp
  USE mapping, ONLY: abcde_map
  IMPLICIT NONE

  INTEGER,INTENT(IN)         :: bond(:,:)
  CHARACTER(LEN=*),INTENT(IN):: group(:)
  INTEGER,INTENT(IN) :: nca         ! # of nodes
  INTEGER,INTENT(IN) :: ia          ! index for node bearing (O.) group
  INTEGER,INTENT(IN) :: ig          ! group bearing the leaving H
  INTEGER,INTENT(IN) :: span        ! span of H migration (5 for 1-5 H shift)
  REAL,INTENT(OUT)   :: arrh(:)
  REAL,INTENT(OUT)   :: FSD(5)
  
  CHARACTER(LEN=LEN(group(1))) :: tgroup(SIZE(group))
  INTEGER :: tbond(SIZE(bond,1),SIZE(bond,2))
  INTEGER :: i,ipath
  INTEGER :: nabcde(span)                   ! deepest path is "span"
  INTEGER :: tabcde(span,mxcp,SIZE(bond,1)) ! deepest path is "span"
  REAL    :: a1,b1,c1,d1,a2,b2,c2,d2        ! polynome coef.
  REAL    :: ar1,ar2,ar3

  tgroup(:)=group(:) ; tbond(:,:)=bond(:,:)
  FSD(1:4)=0. ; FSD(5)=1. ; arrh(:)=0.  ;  ipath=1 
  a1=0 ; b1=1 ; c1=0 ; d1=1 ; a2=0 ; b2=1 ; c2=0 ; d2=1

! check that the span given as an input is authorized
  IF ((span < 4).OR.(span>8)) THEN
    WRITE(6,*) 'ERROR in isomrate.f: Span <4 or >8'
    STOP "in isomrate"
  ENDIF

! get the track between the O. and the leaving H
  CALL abcde_map(tbond,ia,nca,nabcde,tabcde)
! get the right pathway to reach the leaving H
  DO i=1,nabcde(span-1)
    IF (tabcde(span-1,i,span-1)==ig) THEN
      ipath=i  ;  EXIT
    ENDIF
  ENDDO

! reference rate coefficients and correction for F_span_tunnel:
  IF (tgroup(ig)(1:3)=='CH3') THEN
    ar1=3*4E10  ;  ar2=0.  ;  ar3=3825.
    IF      (span==6) THEN ;  ar1=ar1*1.13  ;  ar3=ar3+59.
    ELSE IF (span==7) THEN ;  ar1=ar1*1.09  ;  ar3=ar3+49.
    ELSE IF (span==8) THEN ;  ar1=ar1*0.85  ;  ar3=ar3-65.
    ELSE IF (span==4) THEN ;  ar1=ar1*0.19  ;  ar3=ar3-678.
    ENDIF
  ELSE IF (tgroup(ig)(1:3)=='CH2') THEN
    ar1=2*4E10  ;  ar2=0.  ;  ar3=3010.
    IF      (span==6) THEN ;  ar1=ar1*1.26  ;  ar3=ar3+117.
    ELSE IF (span==7) THEN ;  ar1=ar1*1.24  ;  ar3=ar3+116.
    ELSE IF (span==8) THEN ;  ar1=ar1*0.94  ;  ar3=ar3-16.
    ELSE IF (span==4) THEN ;  ar1=ar1*0.28  ;  ar3=ar3-556.
    ENDIF
  ELSE IF (tgroup(ig)(1:2)=='CH') THEN
    ar1=4E10  ;  ar2=0. ; ar3=2440.
    IF      (span==6) THEN ;  ar1=ar1*1.25  ;  ar3=ar3+123
    ELSE IF (span==7) THEN ;  ar1=ar1*0.97  ;  ar3=ar3-9
    ELSE IF (span==8) THEN ;  ar1=ar1*0.96  ;  ar3=ar3-15
    ELSE IF (span==4) THEN ;  ar1=ar1*0.36  ;  ar3=ar3-473
    ENDIF
  ENDIF

! Cstrain and Fspan(SD) Correction factors
  IF (span==4) THEN
    ar3=ar3+5440. ; a1=2.68E-3 ; b1=1.46 ; c1=2.4E-3 ; d1=0.498
  ELSE IF (span==6) THEN
    ar3=ar3+136.  ; a1=-4.88E-4 ; b1=0.6 ; c1=-6.12E-4 ; d1=0.854
  ELSE IF (span==7) THEN
    ar3=ar3+1243. ; a1=-3.09E-4 ; b1=0.37 ; c1=-2.51E-4 ; d1=0.562
  ELSE IF (span==8) THEN
    ar3=ar3+2179. ; a1=-2.48E-4 ; b1=0.19 ; c1=3.54E-4 ; d1=0.603
  ENDIF

! Fspan(SD) and Fsubst(SD) : these coefficent are expressed in the form: 
!   F(SD)=Frig*Fconf=(aT + b)*(cT + d)=acT^2 + (ad+bc)T + bd
! giving for Fspan(SD)*Fsubst(SD) something like:
!   aa'cc'T^4 + [(ac(a'd'+b'c')+a'c'(ad+bc)]*T^3 +
!   [(ad+bc)(a'd'+b'c')+(acb'd')+(bda'c')]*T^2 +
!   [b'd'(ad+bc)bd(a'd'+b'c')]*T + bdb'd'
! The coef. of the polynome have to be given in output to compute the 
! appropriate rate. Coef. of the polynome are stored in FSD(1:5).
    
! Substition correction factors

! leaving H in an aldehyde
  IF (group(ig)(1:4)=='CHO ') THEN
    IF (span==5) THEN
      ar1=ar1*1.5 ; ar3=ar3+327.+223. ; a2=7.74E-3 ; b2=2.34 ; c2=-1.03E-3 ; d2=2.72
    ELSE IF (span==6) THEN
      ar1=ar1*1.27 ; ar3=ar3+574.+134. ; a2=4.49E-4 ; b2=2.23 ; c2=-8.32E-4 ; d2=0.802
    ELSE IF (span==7) THEN
      ar1=ar1*1.64 ; ar3=ar3+287.+267. ; a2=1.55E-3 ; b2=4.92 ; c2=-1.03E-3 ; d2=2.72
    ENDIF
  ENDIF

! endo B oxo : i.e. R(O.)-...-CO-CH2-...
  IF (tgroup(tabcde(span-1,ipath,span-2))=='CO') THEN
    IF (tgroup(ig)(1:3)=='CH3') THEN
      IF (span==5) THEN
        ar1=ar1*0.57 ; ar3=ar3-237.+2159. ; b2=0.74
      ELSE IF (span==6) THEN
        ar1=ar1*0.83 ; ar3=ar3-76.+1283. ; a2=1.82E-3 ; b2=1.15
      ENDIF
    ELSEIF (tgroup(ig)(1:2)=='CH') THEN
      IF (span==5) THEN
        ar1=ar1*0.69 ; ar3=ar3-165.+1651. ; b2=0.74
      ELSE IF (span==6) THEN
        ar1=ar1*0.92 ; ar3=ar3-26.+1057. ; a2=1.82E-3 ; b2=1.15
      ENDIF
    ENDIF
  ENDIF
      
! search for generic endo -oxo : R(O.)-CX-..-CO-..-CX-CHX-...
  IF (span > 4) THEN
    DO i=2,span-3
      IF (tgroup(tabcde(span-1,ipath,i))=='CO') THEN
        IF (span==5) THEN
          ar1=ar1*0.92 ; ar3=ar3-34.+354. ; a2=4.02E-3 ; b2=1.1
        ELSE IF (span==6) THEN
          ar1=ar1*0.92 ; ar3=ar3-34.+851. ; a2=4.02E-3 ; b2=1.1
        ENDIF
      ENDIF
    ENDDO
  ENDIF

! and acyl-oxi R(O.)-CO-R-CHX- ... subst correction factor
  IF (tgroup(ia)=='CO(O.)') THEN
    IF (tgroup(ig)(1:3)=='CH3') THEN
      IF (span==5) THEN
        ar1=ar1*1.26 ; ar3=ar3+114.+667. ; a2=7.3E-4 ; b2=0.13
      ELSE IF (span==6) THEN
        ar1=ar1*1.26+549. ; ar3=ar3+114. ; a2=7.3E-4 ; b2=0.13
      ENDIF
    ELSEIF (tgroup(ig)(1:2)=='CH') THEN
      IF (span==5) THEN
        ar1=ar1*1.43+667. ; ar3=ar3+187. ; a2=4.27E-4 ; b2=6.03E-2
      ELSE IF (span==6) THEN
        ar1=ar1*1.43+549. ; ar3=ar3+187. ; a2=4.27E-4 ; b2=6.03E-2
      ENDIF
    ENDIF
  ENDIF

! alpha-OH subst correction factor
  IF (INDEX(tgroup(ig),'(OH)')/=0) THEN
    IF (span==5) THEN
      ar1=ar1*1.62 ; ar3=ar3+258.-310. ; a2=-9.21E-3 ; b2=6.00
    ELSE IF (span==6) THEN
      ar1=ar1*1.62 ; ar3=ar3+258.-539. ; a2=-3.98E-2 ; b2=18.2
    ENDIF
  ENDIF

! R(OH)(O.) subst correction factor
  IF (INDEX(tgroup(ia),'(OH)')/=0) THEN
    ar1=ar1*1.17;  ar3=ar3+68.-956. ; a2=-5.61E-4 ; b2=0.633 ; c2=1.16E-4 ; d2=0.313
  ENDIF

! endo beta OH
  IF (INDEX(tgroup(tabcde(span-1,ipath,span-2)),'(OH)')/=0) THEN
    IF (span==5) THEN
      ar1=ar1*0.89 ; ar3=ar3-47.+974. ;  a2=-3.64E-3 ; b2=6.70
    ELSE IF (span==6) THEN
      ar1=ar1*0.93 ; ar3=ar3-35.-30.
    ENDIF
  ENDIF

! exo beta OH
  DO i=1,nca
    IF ((tbond(i,ig)==1).AND.(i/=tabcde(span-1,ipath,span-2)).AND.(INDEX(tgroup(i),'(OH)')/=0)) THEN
      IF (span==5) THEN
        ar3=ar3-252. ; a2=-4.53E-3 ; b2=2.12
      ELSE IF (span==6) THEN
        ar1=ar1*1.18 ; ar3=ar3+83.-392. ; a2=-4.53E-3 ; b2=2.12
      ENDIF
      EXIT
    ENDIF
  ENDDO

  FSD(1)=a1*a2*c1*c2
  FSD(2)=(a1*c1*(a2*d2+b2*c2) + a2*c2*(a1*d1+b1*c1))
  FSD(3)=((a1*d1+b1*c1)*(a2*d2+b2*c2) + (a1*c1*b2*d2) + (a2*c2*b1*d1))
  FSD(4)=b2*d2*(a1*d1+b1*c1) + b1*d1*(a2*d2+b2*c2)
  FSD(5)=b1*d1*b2*d2

  arrh(1)=ar1  ;  arrh(2)=ar2  ;  arrh(3)=ar3

END SUBROUTINE roisom_ver

!=======================================================================
! PURPOSE: Find rate constants for H-atom abstraction in isomerization
! reactions from alkoxy. Based on SAPRC99 and Atkinson 2007 and uses
! Kwok & Atkinson SAR. 
!=======================================================================
SUBROUTINE roisom_atk(tbond,tgroup,nca,ig,arrh)
  IMPLICIT NONE

  INTEGER,INTENT(IN) :: tbond(:,:)
  CHARACTER(LEN=*),INTENT(IN) :: tgroup(:)
  INTEGER,INTENT(IN) :: nca         ! # of nodes
  INTEGER,INTENT(IN) :: ig
  REAL,INTENT(OUT)   :: arrh(:)

  INTEGER :: i,j
  REAL    :: mult,ar1,ar2,ar3

  arrh(:)=0.  ;  ar1=0. ;  ar2=0. ;  ar3=0.

! find k(0) value: Except when noted, value are from Carter 1999 and are
! described in the SAPRC99 documentation. UPDATE from Atkinson 2007 (ric)
  IF (tgroup(ig)(1:3)=='CH3')     THEN ; ar1=1.2E11 ; ar2=0. ; ar3=3825.
  ELSE IF(tgroup(ig)(1:3)=='CH2') THEN ; ar1=8.0E10 ; ar2=0. ; ar3=3010.
  ELSE IF(tgroup(ig)(1:3)=='CHO') THEN ; ar1=0.8E11 ; ar2=0. ; ar3=2890.
  ELSE IF(tgroup(ig)(1:2)=='CH')  THEN ; ar1=4.0E10 ; ar2=0. ; ar3=2440.
  ELSE
    WRITE(6,'(a)') '--error-- in rabsisom. Group not allowed for isomerisation:'
    WRITE(6,'(a)') TRIM(tgroup(ig))
    STOP "in roisom_atk"
  ENDIF

! find factors based on substituents: multipliers are identical to the 
! OH+VOC factors but to the power of 1.5 (update to 1.3 from 
! Atkinson 2007), i.e. F(X)isom=F(X)OH^1.3

! on same carbon:
  mult=1.
  IF (INDEX(tgroup(ig),'(OH)')/=0) THEN
    mult=3.50  ;  ar3=ar3-298.*1.3*log(mult)
  ENDIF
  IF (INDEX(tgroup(ig),'(ONO2)')/=0) THEN
     mult=0.04 ;  ar3=ar3-298.*1.3*log(mult)
  ENDIF
! group contribution for -OOH is based following the comparison between
! rate constant for CH3OOH+OH -> CH2OOH and the rate constant provided 
! by Kwok for the CH3 group. This lead to a factor 13 for the -OOH group
! at 298 K. On the other hand, comparison between the rate constant for
! the CH3OH+OH -> CH2OH provide a group contribution for -OH of 5.9 
! while the value given by Kwok is only 3.5. On that basis, the group 
! factor for the -OOH group was lowered by the same factor (i.e 1.7) 
! leading finally to F(-OOH)=7.6. 
  IF (INDEX(tgroup(ig),'(OOH)') /=0) THEN
    mult=7.6  ;  ar3=ar3-298.*1.3*log(mult)
  ENDIF

! on alpha carbons:
  aloop: DO i=1,nca
    mult=1.
    IF (tbond(ig,i)/=0) THEN

! simple alkyl:
      IF (tgroup(i)(1:3)=='CH3') THEN ; mult=1.   ; ar3=ar3-298.*1.3*log(mult) ; ENDIF
      IF(tgroup(i)(1:4)=='CH2 ') THEN ; mult=1.23 ; ar3=ar3-298.*1.3*log(mult) ; ENDIF
      IF(tgroup(i)(1:4)=='CH2(') THEN ; mult=1.23 ; ar3=ar3-298.*1.3*log(mult) ; ENDIF
      IF(tgroup(i)(1:3)=='CH ' ) THEN ; mult=1.23 ; ar3=ar3-298.*1.3*log(mult) ; ENDIF
      IF(tgroup(i)(1:3)=='CH(' ) THEN ; mult=1.23 ; ar3=ar3-298.*1.3*log(mult) ; ENDIF
      IF(tgroup(i)(1:2)=='C('  ) THEN ; mult=1.23 ; ar3=ar3-298.*1.3*log(mult) ; ENDIF
      IF(tgroup(i)(1:2)=='C '  ) THEN ; mult=1.23 ; ar3=ar3-298.*1.3*log(mult) ; ENDIF

! overwrite for 'CO's, carboxylic acid and PAN values 
      IF (tgroup(i)(1:2)=='CO')  THEN ; mult=0.75 ; ar3=ar3-298.*1.3*log(mult) ; ENDIF
      IF (tgroup(i)(1:3)=='CHO') THEN ; mult=0.75 ; ar3=ar3-298.*1.3*log(mult) ; ENDIF
      IF (tgroup(i)(1:7)=='CO(OH)') THEN ; mult=0.74 ; ar3=ar3-298.*1.3*log(mult) ; ENDIF
      IF (INDEX(tgroup(i),'(ONO2)')/=0) THEN ; mult=0.20 ; ar3=ar3-298.*1.3*log(mult) ; ENDIF
! For PAN, koh is close to the rate constant for CH3 (respectively
! 1.1E-13 and 1.36E-13). Therefore group factor for PAN was set to 1.
      IF (tgroup(i)(1:10)=='CO(OONO2)') THEN ; mult=1.00 ; ar3=ar3-298.*1.3*log(mult) ; ENDIF

! on beta carbons: for -CH2CO-, use MULT=3.9. Group contribution
! taken above for the CH2 or CH group must be removed. 
      IF (mult==1.23) THEN
        DO j=1,nca
          IF (tbond(i,j)/=0 .AND. j/=ig) THEN
             IF (tgroup(j)(1:2)=='CO') THEN
               mult=3.9 ;  ar3=ar3-298.*log(mult)+298.*1.3*log(1.23)
             ENDIF
             IF (tgroup(j)(1:3)=='CHO') THEN
               mult=3.9 ;  ar3=ar3-298.*log(mult)+298.*1.3*log(1.23)
             ENDIF
          ENDIF
        ENDDO
      ENDIF

    ENDIF
  ENDDO aloop

  arrh(1)=ar1  ;  arrh(2)=ar2  ;  arrh(3)=ar3

END SUBROUTINE roisom_atk

END MODULE roisotool
