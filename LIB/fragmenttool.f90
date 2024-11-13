MODULE fragmenttool
IMPLICIT NONE
CONTAINS
!=======================================================================
! PURPOSE: fragments a species into 2 parts. The bond matrix provided
! as input contains the 2 fragments that must be converted as 2 new
! species provided as output (not standardized). It is thus assumed 
! that the bond breaking is already done in the calling program, by 
! setting the corresponding element of the bond matrix to zero
!=======================================================================
SUBROUTINE fragm(bond,group,chem1,chem2)
  USE keyparameter, ONLY: mxcp
  USE mapping, ONLY: gettrack
  USE reactool, ONLY: rebond     
  USE toolbox, ONLY: erase_blank 
  IMPLICIT NONE

  INTEGER,INTENT(IN)           :: bond(:,:) ! bond matrix (containing the 2 fragments)
  CHARACTER(LEN=*),INTENT(IN)  :: group(:)  ! group matrix
  CHARACTER(LEN=*),INTENT(OUT) :: chem1     ! 1st fragment
  CHARACTER(LEN=*),INTENT(OUT) :: chem2     ! 2nd fragment

  CHARACTER(LEN=LEN(group(1))) :: tgrp1(SIZE(group))
  INTEGER         :: tbnd1(SIZE(bond,1),SIZE(bond,2))
  INTEGER         :: trace1(SIZE(bond,1))
  
  CHARACTER(LEN=LEN(group(1))) :: tgrp2(SIZE(group))
  INTEGER         :: tbnd2(SIZE(bond,1),SIZE(bond,2))
  INTEGER         :: trace2(SIZE(bond,1))
  
  CHARACTER(LEN=LEN(group(1))) :: tgroup(SIZE(group))
  INTEGER         :: tbond(SIZE(bond,1),SIZE(bond,1))
  INTEGER         :: i,nring,ngr

  INTEGER         :: track(mxcp,SIZE(bond,1))
  INTEGER         :: trlen(mxcp)
  INTEGER         :: ntr
      
  ngr=COUNT(group/=' ')

! copy bond and groups
  tgroup(:)=group(:) ; tbond(:,:)=bond(:,:)
  CALL erase_blank(tbond,tgroup)

! initialize tracer arrays
  chem1=' '     ;  chem2=' '
  trace1(:)=0   ;  trace2(:)=0
  tgrp1(:)=' '  ;  tgrp2(:)=' '
  tbnd1(:,:)=0  ;  tbnd2(:,:)=0

! connect TRACE1 (if not connected, skip)
  trace1(1)=1
  CALL gettrack(tbond,1,ngr,ntr,track,trlen)
  DO i=1,ntr
    trace1(track(i,1:trlen(i)))=1
  ENDDO
 
! construct TRACE2 from non-blank remainder 
  DO i=1,ngr
    IF( (trace1(i)==0) .AND. (tgroup(i)/=' ') ) trace2(i) = 1
  ENDDO

! write new group1, group2, bond1 and bond2 elements:
  DO i=1,ngr
    IF (trace1(i)==1) THEN
      tgrp1(i)=tgroup(i)
      tbnd1(i,1:ngr)=tbond(i,1:ngr)
    ENDIF
  
    IF (trace2(i)==1) THEN
      tgrp2(i)=tgroup(i)
      tbnd2(i,1:ngr)=tbond(i,1:ngr)
    ENDIF
  ENDDO
  
  CALL rebond(tbnd1,tgrp1,chem1,nring)
  CALL rebond(tbnd2,tgrp2,chem2,nring)
      
END SUBROUTINE fragm

END MODULE fragmenttool
