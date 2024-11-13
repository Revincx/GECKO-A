MODULE toolbox
IMPLICIT NONE
CONTAINS
!=======================================================================
! PURPOSE: increment and add a reference to the reference list provided
! as input.
!=======================================================================
SUBROUTINE addref(progname,reference,nref,reflist,chem)
  USE keyparameter, ONLY: mxlfo
  IMPLICIT NONE
  CHARACTER(LEN=*),INTENT(IN) :: progname
  CHARACTER(LEN=*),INTENT(IN) :: reference
  INTEGER,INTENT(INOUT) :: nref
  CHARACTER(LEN=*),INTENT(INOUT) :: reflist(:)
  CHARACTER(LEN=*),INTENT(IN),OPTIONAL :: chem

  CHARACTER(LEN=mxlfo) :: species
  INTEGER              :: i

! check if chem is provided
  IF (PRESENT(chem)) THEN; species=chem ; ELSE ; species=' ' ; ENDIF

  DO i=1,nref
    IF (reflist(i)==reference) RETURN ! ref already in list
  ENDDO

  nref=nref+1
  IF (nref> SIZE(reflist)) THEN
    PRINT*, "--error-- too many references added for the list"
    PRINT*, "          in subroutine:",TRIM(progname)
    PRINT*, "size of reflist is: ", SIZE(reflist)
    IF (species(1:1)/=' ') &
    PRINT*, "          for the species:",TRIM(species)
    PRINT*, "          reference was:",TRIM(reference)
    PRINT*,reflist(:)
    STOP "in addref"  
  ENDIF
  reflist(nref)=reference
END SUBROUTINE addref
  
!=======================================================================
! PURPOSE: increment number of reactions.
!=======================================================================
SUBROUTINE addrx(progname,chem,nr,flag)
  IMPLICIT NONE

  CHARACTER(LEN=*),INTENT(IN) :: progname
  CHARACTER(LEN=*),INTENT(IN) :: chem
  INTEGER,INTENT(INOUT) :: nr
  INTEGER,INTENT(INOUT) :: flag(:)
  CHARACTER(LEN=70)     :: mesg

  nr=nr+1
  IF (nr>SIZE(flag)) THEN
    mesg='-error- too many reactions created for species: '//TRIM(chem)//' in addrx'
    CALL stoperr(progname,mesg,chem)
!    PRINT*, "--error-- too many reactions created for species: ", TRIM(chem)
!    PRINT*, "          in subroutine: ",TRIM(progname)
!    STOP "in addrx"
  ENDIF
  flag(nr) = 1
END SUBROUTINE addrx

!=======================================================================
! PURPOSE: increment number of products in a reaction.
!=======================================================================
SUBROUTINE add1tonp(progname,chem,np)
  USE keyparameter, ONLY: mxpd  ! max # of possible products in a reaction
  IMPLICIT NONE

  CHARACTER(LEN=*),INTENT(IN) :: progname
  CHARACTER(LEN=*),INTENT(IN) :: chem
  INTEGER,INTENT(INOUT) :: np

  np=np+1
  IF (np>mxpd) THEN
    PRINT*, "--error-- too many products created, np > mxpd"
    PRINT*, "          for species: ", TRIM(chem)
    PRINT*, "          in subroutine: ",TRIM(progname)
    STOP "in add1tonp"
  ENDIF
END SUBROUTINE add1tonp

!=======================================================================
! PURPOSE: generates error message output
!=======================================================================
SUBROUTINE stoperr(prog,mesg,chem)
  IMPLICIT NONE
  CHARACTER(LEN=*),INTENT(IN) :: chem  ! formula of species triggering error
  CHARACTER(LEN=*),INTENT(IN) :: prog  ! name of calling prog routine
  CHARACTER(LEN=*),INTENT(IN) :: mesg  ! nature of problem

!  LOGICAL,PARAMETER :: kill=.FALSE.
  LOGICAL,PARAMETER :: kill=.TRUE.

  WRITE(6,'(a)') ' '
  WRITE(6,'(a)') '--error-- in: '//TRIM(prog)
  WRITE(6,'(a)') TRIM(mesg)
  WRITE(6,'(a)') 'for species: '//TRIM(chem)
  IF (kill) STOP "in stoperr"

END SUBROUTINE stoperr

!=======================================================================
! PURPOSE: set bond value
!=======================================================================
SUBROUTINE setbond(bondtb,x,y,bondvalue)
  IMPLICIT NONE

  INTEGER,INTENT(IN) :: x,y
  INTEGER,INTENT(IN) :: bondvalue
  INTEGER,INTENT(INOUT) :: bondtb(:,:)

  bondtb(x,y)=bondvalue
  bondtb(y,x)=bondvalue

END SUBROUTINE setbond

!===========================================================
! PURPOSE: remove empty nodes in group and bond
!===========================================================
SUBROUTINE erase_blank(bond,group)
  IMPLICIT NONE

  INTEGER,INTENT(INOUT)          :: bond(:,:)
  CHARACTER(LEN=*),INTENT(INOUT) :: group(:)

  INTEGER         :: i,j,k,last,n,nca,icheck
  INTEGER         :: tbond(SIZE(bond,1),SIZE(bond,2))
  CHARACTER(LEN=LEN(group(1))) :: tgroup(SIZE(group))

  tbond(:,:)=bond(:,:)  ;  tgroup(:)=group(:)

  nca=0
  DO n=1,SIZE(group)
    IF (tgroup(n)(1:1) /= ' ') THEN
      nca=nca+1
      last=n
    ENDIF
  ENDDO
 
  IF (nca<2) RETURN
  IF (last==nca) RETURN

! erase blank lines in bond and group
! icheck is just to know if we are in an infinite loop
  icheck=0  ;  i=1

  ckloop: DO
    IF (i>nca) EXIT ckloop
    IF (tgroup(i)==' ') THEN
      DO j=i,last
        tgroup(j)=tgroup(j+1)
      ENDDO
      tgroup(last)=' '
      DO j=1,last
        DO k=i,last-1
          tbond(k,j)=tbond(k+1,j)
        ENDDO
        tbond(last,j)=0
      ENDDO
      DO j=1,last
        DO k=i,last-1
          tbond(j,k)=tbond(j,k+1)
        ENDDO
        tbond(j,last)=0
      ENDDO
      
    ENDIF
     
    IF (tgroup(i)/=' ') i=i+1
    
    icheck=icheck+1
    IF (icheck>SIZE(bond)) THEN
      WRITE(*,*) "--error-- in erase_blank"
      WRITE(*,*) "infinite loop when erasing blanks lines"
      STOP "in erase_blank"
    ENDIF
  ENDDO ckloop

  bond(:,:)=tbond(:,:)
  group(:)=tgroup(:)
      
END SUBROUTINE erase_blank

! ======================================================================
! PURPOSE: count the number of occurence of a substring in a string.
! ======================================================================
INTEGER FUNCTION countstring(line, str) RESULT(n)
  CHARACTER(LEN=*), INTENT(IN) :: line      ! input line 
  CHARACTER(LEN=*), INTENT(IN) :: str       ! string to found in line
  INTEGER :: p, ipos
 
  n=0
  IF (LEN(str) == 0) RETURN
  p = 1
  DO 
    ipos = INDEX(line(p:), str)
    IF (ipos == 0) RETURN
    n = n+1
    p = p+ipos+LEN(str)-1
  ENDDO
END FUNCTION countstring

!=======================================================================
! PURPOSE: return the rate constant at T for the arrhenius parameter
! provide as input. 
!=======================================================================
REAL FUNCTION kval(arrh,T)
  IMPLICIT NONE
  REAl,INTENT(IN) :: arrh(:)                    ! arrhenius parameter
  REAl,INTENT(IN) :: T                          ! input T
  kval=arrh(1)*(T**arrh(2))*exp(-arrh(3)/T)
  RETURN
END FUNCTION  kval


END MODULE toolbox
