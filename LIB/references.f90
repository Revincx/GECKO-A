MODULE references
  IMPLICIT NONE
  INTEGER,PARAMETER :: mxreac_info=1000 ! max # of comment/references 
  INTEGER,PARAMETER :: mxlcod=10    ! max length of comment's code in database
  INTEGER,PARAMETER :: mxlreac_info=400 ! max length of a comment/reference
  
  INTEGER,SAVE           :: nreac_info         ! # of available "reac_info" (references)  
  CHARACTER(LEN=mxlcod)  :: code(mxreac_info)  ! Code for each "reac_info" 
  CHARACTER(LEN=mxlreac_info):: reac_info(mxreac_info) ! Full text correxponding to the reac_info 

END MODULE references
