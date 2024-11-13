MODULE no3abstool
IMPLICIT NONE
CONTAINS
!SUBROUTINE rabsno3(tbond,tgroup,ig,arrhc)
!=======================================================================
!=======================================================================
! PURPOSE : Find rate constants for H-atom abstraction by NO3 based on 
! Kerdouci et al., chemphyschem, 3909, 2010 and Atmos. Env., 363, 2014.
!=======================================================================
!=======================================================================
SUBROUTINE rabsno3(tbond,tgroup,ig,arrhc)
  IMPLICIT NONE

  CHARACTER(LEN=*),INTENT(IN) :: tgroup(:) ! group matrix
  INTEGER,INTENT(IN) ::         tbond(:,:) ! node matrix
  INTEGER,INTENT(IN) ::         ig         ! group bearing the leaving H
  REAL,INTENT(OUT)   ::         arrhc(:)   ! arrhenius coefficients

  INTEGER :: i,j,nether,nca,ngr
  REAL    :: mult
  LOGICAL :: saturated

  arrhc(:)=0.
  ngr=COUNT(tgroup/=' ')
  nca=0  ;  saturated=.TRUE.
  DO i=1,ngr 
    IF (tgroup(i)(1:1)=='C') nca=nca+1
    IF (tgroup(i)(1:2)=='Cd') saturated=.FALSE.
  ENDDO

! check ig value is ok
  IF (ig>ngr) THEN
    WRITE(6,*) '--error--, in rabsno3, => ig is greater than ngr'
    STOP "in rabsno3"
  ENDIF

! FIND K(0) VALUE
! ------------------
! Values are from Kerdouci et al., 2011
  IF (tgroup(ig)(1:3)=='CH3') THEN
    arrhc(1)=1.00E-18 ; arrhc(2)=0. ; arrhc(3)=0.
  ELSE IF(tgroup(ig)(1:3)=='CH2') THEN
    arrhc(1)=2.56E-17 ; arrhc(2)=0. ; arrhc(3)=0.
  ELSE IF(tgroup(ig)(1:3)=='CHO') THEN
    DO i=1,ngr
      IF (tbond(ig,i)==3) THEN
        arrhc(:)=0. ; RETURN
      ENDIF
    ENDDO
    arrhc(1)=2.415E-15 ; arrhc(2)=0. ; arrhc(3)=0.
  ELSE IF(tgroup(ig)(1:2)=='CH') THEN
    arrhc(1)=1.05E-16  ; arrhc(2)=0. ; arrhc(3)=0.
  ELSE
    WRITE(6,'(a)') '--error-- in rabsno3, no NO3 reaction for: ', TRIM(tgroup(ig))
    STOP "in rabsno3"
  ENDIF
       
! ---------------------------------------
! FIND MULTIPLIERS BASED ON SUBSTITUENTS
! ---------------------------------------

! on same carbon:
  mult=1.
  IF (INDEX(tgroup(ig),'(OH)')/=0) THEN
    mult=18. ;  arrhc(1)=arrhc(1)*mult
  ENDIF
  IF ((INDEX(tgroup(ig),'CHO')/=0).AND.(nca>3).AND.(saturated)) THEN
    mult= -14.5 + 21.4*(1-EXP(-0.43*nca)) ; arrhc(1)=arrhc(1)*mult
  ENDIF

! on alpha carbons:
  nether=0
  DO i=1,ngr
    mult=1.
    IF (tbond(ig,i)/=0) THEN
! simple alkyl:
      IF (tgroup(i)(1:3)=='CH3') THEN
        mult= 1. ; arrhc(1)=arrhc(1)*mult
      ENDIF

      IF (tgroup(i)(1:4) == 'CH2 ') THEN
        mult= 1.02
        DO j=1,ngr
          IF ((tbond(i,j)==3).AND.(j/=ig)) mult=1.
        ENDDO
        arrhc(1)=arrhc(1)*mult
      ENDIF
      IF (tgroup(i)(1:4) == 'CH2(') THEN
        mult= 1.02 ; arrhc(1)=arrhc(1)*mult
      ENDIF

      IF (tgroup(i)(1:3) == 'CH ') THEN
        DO j=1,ngr
          IF ((tbond(i,j)==3).AND.(j/=ig)) THEN ; mult= 10.  ; EXIT
          ELSE                                  ; mult= 1.61
          ENDIF
        ENDDO
        arrhc(1)=arrhc(1)*mult
      ENDIF
      IF(tgroup(i)(1:3) == 'CH(' ) THEN
        DO j=1,ngr
          IF ((tbond(i,j)==3).AND.(j/=ig)) THEN ; mult= 10.  ; EXIT
          ELSE                                  ; mult= 1.61
          ENDIF
        ENDDO
        arrhc(1)=arrhc(1)*mult
      ENDIF

      IF(tgroup(i)(1:2) == 'C('  ) THEN
        DO j=1,ngr
          IF ((tbond(i,j)==3).AND.(j/=ig)) THEN ; mult= 48.  ; EXIT
          ELSE                                  ; mult= 2.03
          ENDIF
        ENDDO
        arrhc(1)=arrhc(1)*mult
      ENDIF
      IF(tgroup(i)(1:2) == 'C '  ) THEN
        DO j=1,ngr
          IF ((tbond(i,j)==3).AND.(j/=ig)) THEN ; mult= 48.  ; EXIT
          ELSE                                  ; mult= 2.03
          ENDIF
        ENDDO
        arrhc(1)=arrhc(1)*mult
      ENDIF

! overwrite for carbonyls, alcohols and double bonds values 
      IF (tgroup(i)(1:2)=='CO') THEN 
        mult=0.64 ; arrhc(1)=arrhc(1)*mult
      ENDIF
      IF (tgroup(i)(1:3)=='CHO') THEN
        mult=143. ; arrhc(1)=arrhc(1)*mult
      ENDIF
      IF (tgroup(i)(1:2)=='Cd') THEN
        mult=1.   ; arrhc(1)=arrhc(1)*mult
      ENDIF

! Ethers and esters
! For acetal, consider only one -O- influence
      IF (INDEX(tgroup(i),'-O-')/=0) THEN
        DO j=1,ngr
          IF ((tbond(i,j)==3).AND.(j/=ig)) THEN
!! ether in alpha of ig
            IF (tgroup(j)(1:3)=='CH3') THEN
              mult= 130. ;  arrhc(1)=arrhc(1)*mult
            ELSE IF(tgroup(j)(1:4) == 'CH2 ') THEN
              mult= 58.  ;  arrhc(1)=arrhc(1)*mult
            ELSE IF(tgroup(j)(1:4) == 'CH2(') THEN
              mult= 58.  ;  arrhc(1)=arrhc(1)*mult
            ELSE IF(tgroup(j)(1:3) == 'CH ' ) THEN
              mult= 23.  ;  arrhc(1)=arrhc(1)*mult
            ELSE IF(tgroup(j)(1:3) == 'CH(' ) THEN
              mult= 23.  ;  arrhc(1)=arrhc(1)*mult
            ELSE IF(tgroup(j)(1:2) == 'C('  ) THEN
              mult= 495. ;  arrhc(1)=arrhc(1)*mult
            ELSE IF(tgroup(j)(1:2) == 'C '  ) THEN
              mult= 495. ;  arrhc(1)=arrhc(1)*mult
            ENDIF
          ENDIF                 
        ENDDO
      ENDIF
    ENDIF
  ENDDO

END SUBROUTINE rabsno3


END MODULE no3abstool
