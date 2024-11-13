MODULE minidict
  USE keyparameter, ONLY: mxlfo,mxlco  ! string length in the dictionary
  IMPLICIT NONE

  INTEGER, PARAMETER   :: mxminid=3000
  INTEGER              :: nmini
  CHARACTER(LEN=mxlco) :: nam_minid(mxminid)
  CHARACTER(LEN=mxlfo) :: fo_minid(mxminid)

CONTAINS

! ======================================================================
SUBROUTINE get_fo(nam,formula)
  IMPLICIT NONE
  CHARACTER(LEN=mxlco), INTENT(IN) :: nam
  CHARACTER(LEN=mxlfo), INTENT(OUT) :: formula
  INTEGER :: i

  formula=' '
  DO i=1,nmini
    IF (nam_minid(i)==nam) THEN 
      formula=fo_minid(i)  ; RETURN
    ENDIF
  ENDDO
  formula=nam  ! no formula found if that point is reached
END SUBROUTINE get_fo

! ======================================================================
SUBROUTINE add_fo(nam,formula)
  IMPLICIT NONE
  CHARACTER(LEN=mxlco), INTENT(IN) :: nam
  CHARACTER(LEN=mxlfo), INTENT(IN) :: formula
  INTEGER :: i
  LOGICAL :: toadd
  
  toadd=.TRUE.
  DO i=1,nmini
    IF (formula==fo_minid(i)) THEN
      toadd=.FALSE.  ; EXIT 
    ENDIF
  ENDDO

  IF (toadd) THEN
    nmini=nmini+1
    IF (nmini>mxminid) THEN
      WRITE(*,*) "maximum formula in minid reached - see module minidict"
      STOP "in minid"
    ENDIF
    fo_minid(nmini)=formula
    nam_minid(nmini)=nam
  ENDIF
END SUBROUTINE add_fo
  
! ======================================================================
SUBROUTINE clean_minid()
  IMPLICIT NONE
  nmini=0 ; nam_minid=' ' ; fo_minid=' '
END SUBROUTINE clean_minid


END MODULE minidict
