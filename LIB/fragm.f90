! This subroutine fragments a chemical into two parts.  
! It is assumed that the bond breaking is already done in the calling
! program, by setting the corresponding element of BOND to zero
!*************************************************************************
SUBROUTINE fragm(bond,group,chem1,chem2)
  USE keyparameter, ONLY: mxcp,mxnode,mxlfo,mxlgr
  USE mapping, ONLY: gettrack
  USE reactool, ONLY: rebond      
  IMPLICIT NONE

  INTEGER,INTENT(IN)          :: bond(mxnode,mxnode)
  CHARACTER(LEN=mxlgr),INTENT(IN)  :: group(mxnode)
  CHARACTER(LEN=mxlfo),INTENT(OUT) :: chem1, chem2

! internal
  INTEGER         :: tbnd1(SIZE(bond,1),SIZE(bond,2))
  CHARACTER(LEN=mxlgr) :: tgrp1(SIZE(group))
  INTEGER         :: trace1(SIZE(bond,1))
  
  INTEGER         :: tbnd2(SIZE(bond,1),SIZE(bond,2))
  CHARACTER(LEN=mxlgr) :: tgrp2(SIZE(group))
  INTEGER         :: trace2(SIZE(bond,1))
  
  INTEGER         :: tbond(SIZE(bond,1),SIZE(bond,1))
  CHARACTER(LEN=mxlgr) :: tgroup(SIZE(group))
  INTEGER         :: i,j,nring,ngr

  INTEGER         :: track(mxcp,SIZE(bond,1))
  INTEGER         :: trlen(mxcp)
  INTEGER         :: ntr
      
!  IF(wtflag>0) WRITE(6,*) '*fragm*'
  ngr=COUNT(group/=' ')

! copy bond and groups
  tgroup(:)=group(:) ; tbond(:,:)=bond(:,:)
  CALL erase_blank(tbond,tgroup)

! initialize tracer arrays
  chem1=' '     ;  chem2=' '
  trace1(:)=0   ;  trace2(:)=0
  tgrp1(:)=' '  ;  tgrp2(:)=' '
  tbnd1(:,:)=0  ;  tbnd2(:,:)=0

! connect TRACE1
! if not connected, skip
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
!  IF(wtflag>0) print*,'output 1 from fragm: ',chem1
!  IF(wtflag>0) print*,'output 2 from fragm: ',chem2
      
END SUBROUTINE fragm
