MODULE reactool
  IMPLICIT NONE
  CONTAINS

!SUBROUTINE rebond(bond,group,chem,nring) 
!SUBROUTINE swap(gold,pold,gnew,pnew)

! ======================================================================
! Purpose: return the (non-standardised) chemical formula "chem" 
! for the bond and group tables provided has input
! ======================================================================
SUBROUTINE rebond(bond,group,chem,nring) 
  USE keyparameter, ONLY: mxring
  USE stdtool, ONLY: lntree,ckgrppt,mkcopy,dwrite
  USE stdratings, ONLY: ratings
  USE ringtool, ONLY: uniqring
  IMPLICIT NONE

  INTEGER,INTENT(IN)          :: bond(:,:)
  CHARACTER(LEN=*),INTENT(IN) :: group(:)
  CHARACTER(LEN=*),INTENT(OUT):: chem
  INTEGER,INTENT(OUT)         :: nring

  INTEGER :: i,j,k,n,idbflg,beg
  INTEGER :: nca,ncx,last,icheck
  INTEGER :: rjg(mxring,2)
  INTEGER :: rank(SIZE(bond,1))
  INTEGER :: path(SIZE(bond,1),SIZE(bond,1),SIZE(bond,1)) 
  INTEGER :: clngth(SIZE(bond,1),SIZE(bond,1))
  INTEGER :: leaf, last2, ptr
  INTEGER :: ig, pg, ng
  LOGICAL :: lobond(SIZE(bond,1),SIZE(bond,2))
  INTEGER :: locat
  INTEGER :: tbond(SIZE(bond,1),SIZE(bond,2))
  CHARACTER(LEN=LEN(group(1))) :: tgroup(SIZE(group))
  INTEGER :: mxca

  chem=' ' ; clngth(:,:)=0 ; lobond(:,:)=.false. ; path(:,:,:)=0
  tgroup(:)=group(:)  ;  tbond(:,:)=bond(:,:) ; mxca=SIZE(bond,1)

  nca=0 ;  idbflg=0
  DO n=1,mxca
    IF (tgroup(n)(1:1) /= ' ') THEN
      nca=nca+1  ;  last=n
    ENDIF
    IF (INDEX(tgroup(n),'Cd')/=0) idbflg=1
  ENDDO

! check that the functional groups are correctly sorted for each group
  DO i=1,last
    locat=INDEX(tgroup(i),')(')
    IF (locat/=0) CALL ckgrppt(locat,tgroup(i))
  ENDDO
  IF (nca<1) RETURN

! If there is only one group in tgroup, chem=group and return
  IF (nca==1) THEN
    DO i=1,last
      IF (tgroup(i)/= ' ') THEN
        chem=tgroup(i) ; RETURN
      ENDIF
    ENDDO
  ENDIF

! Find if we have ring(s)
  ncx=0
  DO i=1,last
    DO j=i+1,last
      IF (tbond(i,j)>0) ncx=ncx+1
    ENDDO
  ENDDO
  nring=ncx-nca+1
      
! erase blank lines in bond and group
  icheck=0  ! check for infinite loop
  i=1
  rmloop: DO
   IF (i>nca) EXIT rmloop
   IF (tgroup(i)==' ') THEN
     DO j=i,last
       tgroup(j)=tgroup(j+1)
     ENDDO
     tgroup(last)=' '
     DO j=1,last
       DO k=i,last-1
         tbond(k,j)=tbond(k+1,j)
       ENDDO
     ENDDO
     DO j=1,last
       DO k=i,last-1
         tbond(j,k)=tbond(j,k+1)
       ENDDO
     ENDDO
   ENDIF
   
   IF (tgroup(i)/=' ') i=i+1

   icheck=icheck+1
   IF (icheck>last) THEN
     WRITE(*,*) '--error-- in rebond. Infinite loop when erasing blanks lines'
     STOP "in rebond"
   ENDIF
  ENDDO rmloop

! find 'ends' of rings, add ring-join characters.
! Since uniqring modify bond and group, a copy of the table are used
  CALL ratings(nca,tgroup,tbond,nring,rank)
  IF (nring>0) THEN
    CALL uniqring(nring,nca,tgroup,tbond,rank,rjg)
    DO k=1,nring
      i=rjg(k,1)   ; j=rjg(k,2)
      tbond(i,j)=0 ; tbond(j,i)=0
    ENDDO
  ENDIF

! make a logical copy of the bond matrix
  WHERE (tbond(:,:)/=0) lobond(:,:)=.TRUE.

! bond between node 1 and 2 might be broken.
  beg=0
  grloop: DO i=1,last
    IF (tgroup(i)/=' ') THEN
      beg=i
      DO j=i+1,mxca
        IF (lobond(i,j)) THEN
          CALL lntree(tbond,i,j,nca,clngth,path)
          EXIT grloop
        ENDIF
      ENDDO  
    ENDIF
  ENDDO grloop

! look down-top for the very longest tree ...
  wrtloop: DO i=1,nca
    IF (clngth(beg,i)/=0) THEN
      leaf=path(beg,i,clngth(beg,i))
      last2=path(beg,i,clngth(beg,i)-1)
      CALL lntree(tbond,leaf,last2,nca,clngth,path)
      DO j=1,nca
        IF (clngth(leaf,j)/=    0) THEN
          ptr=1
          DO k=1,nca
            ig=path(leaf,j,k)
            IF (ig/=0) THEN
             IF (k>1)   pg=path(leaf,j,k-1)
             IF (k<nca) ng=path(leaf,j,k+1)
             CALL mkcopy(lobond,tgroup,nca,rank,nring,ig,pg,ng,ptr,chem)
            ENDIF
          ENDDO
          EXIT wrtloop ! the first copy has been written in chem
        ENDIF
      ENDDO
    ENDIF
  ENDDO wrtloop

! add '=' id double bond
  IF (idbflg/=0) CALL dwrite(chem)

END SUBROUTINE rebond


! ======================================================================
! Purpose : swaps pieces pold, pnew in group (gold, gnew)
! ======================================================================
SUBROUTINE swap(gold,pold,gnew,pnew)
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN)  :: gold,pold,pnew
  CHARACTER(LEN=*), INTENT(OUT) :: gnew

  INTEGER :: lengr,lg,lold,lnew,ibeg,iend,n,i,j
  INTEGER :: lengr2

!- initialize
  lengr=LEN(gold)  ;  lengr2=LEN(gnew)
  IF (lengr /= lengr2) THEN
    PRINT*, "stop in swap - wrong size"
    STOP "in swap"
  ENDIF 
  
  gnew=' '
  lg  =INDEX(gold,' ')-1
  lold=INDEX(pold,' ')-1
  lnew=INDEX(pnew,' ')-1

  IF (lg+lnew-lold > lengr) THEN
    WRITE (6,'(a)') '--error-- in swap'
    WRITE (6,'(a)') 'new group more than lengr characters'
    STOP 'in swap'
  ENDIF

!- locate beginning and ending index of old piece:
  ibeg=INDEX(gold,pold(1:lold))
  IF (ibeg==0) THEN
    WRITE (6,'(a)') '--error-- in swap'
    WRITE (6,'(a)') pold,' not in ',gold
    STOP "in swap"
  ENDIF
  iend=ibeg+lold-1
      
!- write the new group:
  n=0
  charloop : DO i=1,lengr
    IF (n==lengr) EXIT charloop
      IF (i<ibeg .OR. i>iend) THEN
        n=n+1
        gnew(n:n)=gold(i:i)
      ENDIF
      IF (i==ibeg) THEN
        DO j=1,lnew
          n=n+1
          gnew(n:n)=pnew(j:j)
        ENDDO
      ENDIF
  ENDDO charloop  

END SUBROUTINE swap

END MODULE reactool

