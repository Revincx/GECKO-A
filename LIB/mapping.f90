MODULE mapping
  IMPLICIT NONE
  CONTAINS

!SUBROUTINE nodmap(bond,top,nca,npos,nnod,tnod)
!SUBROUTINE gettrack(bond,top,nca,ntr,track,trlen)
!SUBROUTINE abcde_map(bond,top,nca,nabcde,tabcde)
!SUBROUTINE ciptree(chem,bond,group,ngr,nabcde,tabcde,ciptr)
!SUBROUTINE alkenetrack(chem,bond,group,ngr,ncdtrack,cdtracklen,cdtrack)
!SUBROUTINE estertrack(chem,bond,group,ngr,netrack,etracklen,etrack)
!SUBROUTINE chemmap(chem,node,group,bond,ngrp,nodetype,alifun,cdfun,arofun,mapfun,funflg,tabester,nfcd,nfcr,ierr)

! ======================================================================
! PURPOSE: Browse the bond matrix to find all possible node at a given        
! position (npos), starting from a given node (top). Ring are allowed
! (track is stopped when a full loop along the circle is made).       
!                                                                     
! The subroutine call "gettrack" which give all the possible track,   
! starting from top.                                                  
!  track(*,2) give all nodes in alpha position regarding node top     
!  track(*,3) give all nodes in beta position regarding node top      
!  track(*,4) give all nodes in gamma position regarding node top     
!  etc...                                                             
! ======================================================================
SUBROUTINE nodmap(bond,top,nca,npos,nnod,tnod)
  USE keyparameter, ONLY:mxcp

  IMPLICIT NONE
  INTEGER,INTENT(IN) :: top     ! starting node number
  INTEGER,INTENT(IN) :: nca     ! number of nodes
  INTEGER,INTENT(IN) :: bond(:,:) ! node matrix
  INTEGER,INTENT(IN) :: npos    ! position seeked (n=2 is alpha relative to top)
  INTEGER,INTENT(OUT):: nnod    ! # of distinct nodes at "npos" relative to top
  INTEGER,INTENT(OUT):: tnod(:) ! table of all nodes at position "npos"

! internal:
  INTEGER  track(mxcp,SIZE(bond,1))
  INTEGER  trlen(mxcp)
  INTEGER  ntr
  INTEGER  i,j

  nnod=0  ;  tnod(:)=0
  CALL gettrack(bond,top,nca,ntr,track,trlen)

! avoid duplicate - set 0 to duplicate nodes
  IF (ntr>1) THEN
    trloop: DO i=1,ntr-1
      IF (track(i,npos)==0) CYCLE trloop
      DO j=i+1,ntr
        IF (track(j,npos)/=0) THEN
          IF (track(i,npos)==track(j,npos)) track(j,npos)=0
        ENDIF  
      ENDDO
    ENDDO trloop
  ENDIF
       
! get nodes
  DO i=1,ntr
    IF (track(i,npos)/=0) THEN
      nnod=nnod+1
      tnod(nnod)=track(i,npos)
    ENDIF
  ENDDO
END SUBROUTINE nodmap

! ======================================================================
! PURPOSE : Browse the bond matrix to find all possible tracks starting     
! from a given node (top). A track ends when a "ending" node (i.e. the
! last positions on the carbon skeleton) is reached. Ring are allowed 
! (track is stopped when a full loop along the circle is made).                                    *
!                                                                      
!  track(1,*) give all nodes in 1st track
!  track(2,*) give all nodes in 2nd track       
!  etc...                                                              
!
!  track(*,2) give all nodes in alpha position regarding node top      
!  track(*,3) give all nodes in beta position regarding node top       
!  track(*,4) give all nodes in gamma position regarding node top      
!  etc...                                                              
! ======================================================================
SUBROUTINE gettrack(bond,top,nca,ntr,track,trlen)
  IMPLICIT NONE
  INTEGER, INTENT(IN)  :: top        ! starting node number                                
  INTEGER, INTENT(IN)  :: nca        ! number of groups (nodes)
  INTEGER, INTENT(IN)  :: bond(:,:)  ! bond matrix
  INTEGER, INTENT(OUT) :: ntr        ! # of distinct tracks found starting at "top"
  INTEGER, INTENT(OUT) :: track(:,:) ! various track (see header)
  INTEGER, INTENT(OUT) :: trlen(:)   ! length of each track

! internal:
  INTEGER :: ptr,niv,nod
  INTEGER :: i,j,k
  INTEGER :: memo(SIZE(bond,1))
  INTEGER :: slope, mxtr

! -----------
! Initialize
! -----------
  ptr=0  ;  trlen(:)=0  ;  track(:,:)=0  ;  memo(:)=0
  mxtr = SIZE(track,1)

! initialize parameters to find the tracks 
  track(1,1)=top  ;  memo(1)=top
  niv=1     ! # of node since starting node (i.e. top)
  ptr=0     ! current "search" pointer (i.e. to find next node in the track)
  ntr=1     ! # of track found - current track
  nod=top   ! current node along the track
  slope=1   ! equal 1 when going forward along the track, otherwise 0

! -----------
! start loop
! -----------

! get next node - reentry point
  nextnode : DO    ! exit the loop when all possible nodes found 
    ptr=ptr+1

! end of possible nodes reached - must go backward or exit
    IF (ptr > nca) THEN      
      IF (niv == 1) THEN       ! all the possible track are found - exit
        EXIT nextnode
      ELSE                     ! go backward
        ptr=memo(niv)          ! set pointer to previous memo
        memo(niv)=0            ! reset memo   
        niv=niv-1              ! decrease niv (go backward)
        nod=memo(niv)          ! set new current node 

        IF (slope /= 0) THEN   ! make a new track (if required)
          ntr=ntr+1
          IF (ntr > mxtr) THEN
            PRINT*, '--error-- in gettrack. ntr is greater than mxtr'
            STOP "in gettrack"
          ENDIF
          track(ntr,1:niv)=track(ntr-1,1:niv)
        ENDIF
        track(ntr,niv+1)=0     ! remove previous track
        slope=0                ! remember ... I am now going backward
        CYCLE nextnode
      ENDIF
    ENDIF

! no bond between current node (nod) and next possible node (ptr)
    IF (bond(ptr,nod)==0) CYCLE nextnode  ! get next

! new bond found - must go forward
    IF (bond(ptr,nod) /= 0) THEN
      DO k=1,niv
        IF (ptr == memo(k)) CYCLE nextnode   ! end circle or previous node
      ENDDO
      niv=niv+1                       ! increase niv (go forward)
      track(ntr,niv)=ptr              ! keep track
      nod=ptr                         ! set current node to the new node
      memo(niv)=nod                   ! remember which track was used
      ptr=0                           ! set pointer to 0 (to find next node)
      slope=1                         ! remember ... I am going forward
      CYCLE nextnode                  ! go find next node
    ENDIF

  ENDDO nextnode

! remove the last "case" (fill with top only) and clean track
  track(ntr,:)=0  ;  ntr=ntr-1

! compute the length of each track 
  lentrack:&
  DO i=1,ntr
    DO j=nca,1,-1
      IF (track(i,j)/=0) THEN
        trlen(i)=j
        CYCLE lentrack
      ENDIF
    ENDDO
  ENDDO lentrack

! end of gettrack
  RETURN
END SUBROUTINE gettrack

! ======================================================================
! PURPOSE :                                                            
!  Browse the bond matrix to find all possible node at various position 
!  (i.e. alpha, beta, gamma ...), starting from node 'top'. Ring are   
!  allowed (track is stopped when a full loop along the circle is made).
!                                                                      
! The subroutine call "gettrack" which give all the possible track,    
! starting from top (see gettrack in the same module).                                                   
!                                                                      
! The subroutine returns 2 tables:                                                              
! -nabcde(k)    : number of distinct pathways that end up at a        
!                  position k relative to top (e.g. nabcde(3) gives    
!                  the number of distinct pathways finishing in a      
!                  beta position relative to top                        
! -tabcde(k,i,j): give the pathways (node j), for the track number i  
!                 to reach the position k (k=2 is beta position ...).
!                 For example, tabcde(4,1,j) give the first track to 
!                 reach a gamma position (node given by tabcde(4,1,4)),                 
!                 using the track given by tabcde(4,1,*)                                      
! ======================================================================
SUBROUTINE abcde_map(bond,top,nca,nabcde,tabcde)
  IMPLICIT NONE
  INTEGER, INTENT(IN)  :: top       ! starting node number
  INTEGER, INTENT(IN)  :: nca       ! number of groups (nodes)
  INTEGER, INTENT(IN)  :: bond(:,:) ! bond matrix
  INTEGER, INTENT(OUT) :: nabcde(:) ! # of tracks ending at position k (see header)
  INTEGER, INTENT(OUT) :: tabcde(:,:,:) ! tabcde(k,i,j) is pathways (node j), for the

! internal:
  INTEGER  :: track(SIZE(tabcde,2),SIZE(tabcde,3))
  INTEGER  :: trlen(SIZE(tabcde,2))
  INTEGER  :: ntr, mxtr, mxdeep
  INTEGER  :: ltabc(SIZE(tabcde,1),SIZE(tabcde,2),SIZE(tabcde,3))
  INTEGER  :: i,j,k,ii

! initialize
  mxdeep=SIZE(tabcde,1) ; mxtr = SIZE(tabcde,2)
  nabcde(:)=0  ;  tabcde(:,:,:)=0  ;  ltabc(:,:,:)=0

! get all tracks starting at top
  CALL gettrack(bond,top,nca,ntr,track,trlen)

! range the track for each length (up to mxdeep)
! Remember that 2=alpha position, 3=beta ... 
  IF (ntr > mxtr) THEN
    PRINT*, '--error-- in abcde_map, ntr > mxtr'
    STOP "in abcde_map"
  ENDIF  
  DO k=2,mxdeep
    DO i=1,ntr
      IF (trlen(i) >= k) THEN
        DO j=1,k
           ltabc(k,i,j)=track(i,j)
        ENDDO
      ENDIF        
    ENDDO        
  ENDDO      
 
! avoid duplicate - set 0 to duplicate nodes tracks
  DO k=2,mxdeep
    trackloop: DO i=1,ntr-1
      IF (ltabc(k,i,k)==0) CYCLE
      DO ii=i+1,ntr
        DO j=1,k
          IF (ltabc(k,i,j)/=ltabc(k,ii,j)) CYCLE trackloop
        ENDDO
        DO j=1,k             ! If that point is reached, track i= track ii
          ltabc(k,ii,j)=0
        ENDDO
      ENDDO
    ENDDO trackloop
  ENDDO
 
! get nodes
  DO k=2,mxdeep
    DO i=1,ntr
      IF (ltabc(k,i,k)/=0) THEN
        nabcde(k)=nabcde(k)+1     ! add one more pathway that end up at a k position
        DO j=1,k
         ii=nabcde(k)
         tabcde(k,ii,j)=ltabc(k,i,j)
        ENDDO
      ENDIF
    ENDDO
  ENDDO

! add top nod (needed top scroll group from top, included)
  nabcde(1)=1
  tabcde(1,1,1)=top
  
END SUBROUTINE abcde_map

!=======================================================================
! PURPOSE : Return a tree of the atom masses to rank branches according
! to the Cahn-Ingold-Prelog (CIP) priority rules. The node tree (tabcde) 
! starting from the top node in branch is provided as input.
! For example, the CIT tree returned for -C(ONO2)(CdH=CdH2)COCH2(OH) is:
! line 1:  12   
! line 2:  16  12  12   
! line 3:  16  16  14  12  12  12   1   
! line 4:  16  16  16   1   1   1   1   
! line 5:   1     
!=======================================================================
SUBROUTINE ciptree(chem,bond,group,ngr,nabcde,tabcde,ciptr)
  USE searching, ONLY: search_ipos
  USE toolbox, ONLY:countstring,stoperr
  IMPLICIT NONE

  CHARACTER(LEN=*),INTENT(IN) :: chem     ! formula of current species
  CHARACTER(LEN=*),INTENT(IN) :: group(:) ! group
  INTEGER,INTENT(IN)  :: bond(:,:)        ! bond (might be for a branch of chem only)
  INTEGER,INTENT(IN)  :: ngr              ! # of groups
  INTEGER,INTENT(IN)  :: nabcde(:)        ! number of distinct nodes distant from top
  INTEGER,INTENT(IN)  :: tabcde(:,:,:)    ! node map 
  INTEGER,INTENT(OUT) :: ciptr(:,:)       ! CIP tree (see header) 

  INTEGER :: nciptr(SIZE(ciptr,1))        ! # of atoms (masses) in a CIP tree line
  INTEGER, PARAMETER :: mxbd=5            ! max tree length in a group
  INTEGER :: nC(mxbd)                     ! # carbon in group at pos j
  INTEGER :: nO(mxbd)                     ! # oxygen in group at pos. j
  INTEGER :: nH(mxbd)                     ! # hydrogen in group at pos. j
  INTEGER :: nN(mxbd)                     ! # nitrogen in group at pos. j
  INTEGER :: i,idepth,ig,nnod,ipos
  INTEGER :: pvnod,nfun,ii,nel,irk,ilast

  CHARACTER(LEN=15) :: progname='ciptree'
  CHARACTER(LEN=70) :: mesg

  ciptr(:,:)=0 ; nciptr(:)=0 

! scroll from close to remote nodes (top, alpha, beta, gamma etc)
  DO idepth=1,ngr
    nnod=nabcde(idepth)    ! # of nodes having a distance "idepth" from top (rank 1)
    IF (nnod==0) EXIT      ! no more groups to consider

    nC(:)=0 ; nH(:)=0 ; nO(:)=0 ; nN(:)=0

! loop over all nodes distant idepth from top
    DO i=1,nnod
      ig=tabcde(idepth,i,idepth)  ! current group index

      IF      (group(ig)(1:3)=='CH3')  THEN ; nC(1)=nC(1)+1 ; nH(2)=nH(2)+3
      ELSE IF (group(ig)(1:3)=='CH2')  THEN ; nC(1)=nC(1)+1 ; nH(2)=nH(2)+2
      ELSE IF (group(ig)(1:3)=='CHO')  THEN ; nC(1)=nC(1)+1 ; nH(2)=nH(2)+1 ; nO(2)=nO(2)+2
      ELSE IF (group(ig)(1:2)=='CO')   THEN ; nC(1)=nC(1)+1 ; nO(2)=nO(2)+2        
      ELSE IF (group(ig)(1:2)=='CH')   THEN ; nC(1)=nC(1)+1 ; nH(2)=nH(2)+1 
      ELSE IF (group(ig)(1:4)=='CdH2') THEN ; nC(1)=nC(1)+1 ; nH(2)=nH(2)+2 
      ELSE IF (group(ig)(1:3)=='CdH')  THEN ; nC(1)=nC(1)+1 ; nH(2)=nH(2)+1 
      ELSE IF (group(ig)(1:3)=='CdO')  THEN ; nC(1)=nC(1)+1 ; nO(2)=nO(2)+2 
      ELSE IF (group(ig)(1:1)=='C')    THEN ; nC(1)=nC(1)+1
      ELSE IF (group(ig)(1:3)=='-O-')  THEN ; nO(1)=nO(1)+1 
      ELSE IF (group(ig)(1:2)=='cH')   THEN ; nC(1)=nC(1)+1 ; nH(2)=nH(2)+1 
      ELSE IF (group(ig)(1:1)=='c')    THEN ; nC(1)=nC(1)+1  
      ENDIF

! count twice the C if Cd is the "second" Cd in C=C
      IF (group(ig)(1:2)=='Cd')   THEN 
        IF (idepth>1) THEN
          pvnod=tabcde(idepth,i,idepth-1)   ! index of the previous node in path 
          IF (bond(pvnod,ig)==2) nC(1)=nC(1)+1  ! count 2 C for the 2nd Cd
        ENDIF
      ENDIF

! loop over functional group
      nfun=countstring(group(ig),'(OH)')
      IF (nfun/=0) THEN 
        nO(2)=nO(2)+nfun ; nH(3)=nH(3)+nfun 
      ENDIF

      nfun=countstring(group(ig),'(OOH)')
      IF (nfun/=0) THEN 
        nO(2)=nO(2)+nfun ; nO(3)=nO(3)+nfun ; nH(4)=nH(4)+nfun 
      ENDIF

      nfun=countstring(group(ig),'(OOOH)')
      IF (nfun/=0) THEN 
        nO(2)=nO(2)+nfun ; nO(3)=nO(3)+nfun ; nO(4)=nO(4)+nfun ;  nH(5)=nH(5)+nfun 
      ENDIF

      nfun=countstring(group(ig),'(ONO2)')
      IF (nfun/=0) THEN 
        nO(2)=nO(2)+nfun ; nN(3)=nN(3)+nfun ; nO(4)=nO(4)+2*nfun 
      ENDIF

      nfun=countstring(group(ig),'(OONO2)')
      IF (nfun/=0) THEN 
        nO(2)=nO(2)+nfun ; nO(3)=nO(3)+nfun ; nN(4)=nN(4)+nfun ; nO(5)=nO(5)+2*nfun 
      ENDIF

      nfun=countstring(group(ig),'(NO2)')
      IF (nfun/=0) THEN 
        nN(2)=nN(2)+nfun ; nO(3)=nO(3)+2*nfun
      ENDIF

      nfun=countstring(group(ig),'(OO.)')
      IF (nfun/=0) THEN 
        nO(2)=nO(2)+nfun ; nO(3)=nO(3)+nfun
      ENDIF

      nfun=countstring(group(ig),'(O.)')
      IF (nfun/=0) THEN 
        nO(2)=nO(2)+nfun
      ENDIF
    ENDDO     
    
! add C in CIP tree 
    DO ii=1,mxbd
      IF (nC(ii)==0) CYCLE
      nel=nC(ii) ; irk=idepth+ii-1               ! # of C (nel) to be added at rank irk 
      ipos=search_ipos(12,ciptr(irk,:))          ! ipos: index to add element
      ilast=nciptr(irk)                          ! last occupied pos. in rank
      IF (ipos<=ilast) &                         ! make room if needed
        ciptr(irk,ipos+nel:ipos+nel+ilast)=ciptr(irk,ipos:ipos+ilast) 
      ciptr(irk,ipos:ipos+nel-1)=12              ! add C
      nciptr(irk)=nciptr(irk)+nel                ! store length of the rank
    ENDDO

! add O in CIP tree 
    DO ii=1,mxbd
      IF (nO(ii)==0) CYCLE
      nel=nO(ii) ; irk=idepth+ii-1               ! # of O (nel) to be added at rank irk 
      ipos=search_ipos(16,ciptr(irk,:))          ! ipos: index to add element
      ilast=nciptr(irk)                          ! last occupied pos. in rank
      IF (ipos<=ilast) &                         ! make room if needed
        ciptr(irk,ipos+nel:ilast+nel)=ciptr(irk,ipos:ilast) 
      ciptr(irk,ipos:ipos+nel-1)=16              ! add O
      nciptr(irk)=nciptr(irk)+nel                ! store length of the rank
    ENDDO
 
! add N in CIP tree
    DO ii=1,mxbd
      IF (nN(ii)==0) CYCLE
      nel=nN(ii) ; irk=idepth+ii-1               ! # of N (nel) to be added at rank irk  
      ipos=search_ipos(14,ciptr(irk,:))          ! ipos: index to add element
      ilast=nciptr(irk)                          ! last occupied pos. in rank
      IF (ipos<=ilast) &                         ! make room if needed
        ciptr(irk,ipos+nel:ilast+nel)=ciptr(irk,ipos:ilast) 
      ciptr(irk,ipos:ipos+nel-1)=14              ! add N
      nciptr(irk)=nciptr(irk)+nel                ! store length of the rank
    ENDDO

! add H in CIP tree
    DO ii=2,mxbd                       
      IF (nH(ii)==0) CYCLE
      nel=nH(ii) ; irk=idepth+ii-1               ! # of H (nel) to be added at rank irk  
      ipos=search_ipos(1,ciptr(irk,:))           ! ipos: index to add element
      ilast=nciptr(irk)                          ! last occupied pos. in rank
      IF (ipos<=ilast) &                         ! make room if needed
        ciptr(irk,ipos+nel:ilast+nel)=ciptr(irk,ipos:ilast) ! make room if needed
      ciptr(irk,ipos:ipos+nel-1)=1               ! add H
      nciptr(irk)=nciptr(irk)+nel                ! store length of the rank
    ENDDO

  ENDDO 

! check size of CIP tree before returning
  DO i=1,SIZE(nciptr)
    IF (nciptr(i)==0) EXIT  ! no more lines
    IF (nciptr(i) > SIZE(ciptr,2)) THEN
      mesg="too many elements added in a line of CIP tree"
      CALL stoperr(progname,mesg,chem)
    ENDIF
  ENDDO
END SUBROUTINE ciptree

!=======================================================================
! PURPOSE: create tracks for Cd groups. For example, the molecule
! CH3CdH=CdH-CH2-CdH=CdH-CdH=CdH2 as 2 tracks: 1) for the C=C group and
! 2) for the C=C-C=C group.
! The routine output the number of track, the length of each track and
! the node belonging to each track.
!=======================================================================
SUBROUTINE alkenetrack(chem,bond,group,ngr,ncdtrack,cdtracklen,cdtrack)
  USE keyparameter, ONLY: mxcp           
  USE toolbox, ONLY: stoperr
  IMPLICIT NONE

  CHARACTER(LEN=*),INTENT(IN) :: chem      ! formula
  CHARACTER(LEN=*),INTENT(IN) :: group(:)  ! grp matrix
  INTEGER,INTENT(IN)  :: bond(:,:)         ! bond matrix
  INTEGER,INTENT(IN)  :: ngr               ! # of nodes in chem
  INTEGER,INTENT(OUT) :: ncdtrack          ! # of Cd tracks 
  INTEGER,INTENT(OUT) :: cdtrack(:,:)      ! (i,j) Cd "j" nodes for the ith track         
  INTEGER,INTENT(OUT) :: cdtracklen(:)     ! length of cdtrack i

  INTEGER :: cdbond(SIZE(group),SIZE(group))
  INTEGER :: cdpst(SIZE(group))            ! none, primary, secondary or tertiary Cd
  INTEGER :: i,j,ncd,nprim,mxcd

  INTEGER :: track(mxcp,SIZE(group))
  INTEGER :: trlen(mxcp)
  INTEGER :: ntr

  CHARACTER(LEN=12) :: progname='alkenetrack '
  CHARACTER(LEN=70) :: mesg

  ncdtrack=0    ; cdtrack(:,:)=0  ; cdtracklen(:)=0 

  mxcd=SIZE(cdtrack,2)
  cdbond(:,:)=0
  
! ------------------------------
! Create and check the Cd tracks
! ------------------------------

! create a bond matrix for Cd only
  ncd=0
  DO i=1,ngr
    IF (group(i)(1:2)=='Cd') THEN
      ncd=ncd+1
      DO j=1,ngr
        IF (j==i) CYCLE
        IF (bond(i,j)/=0) THEN
          IF (group(j)(1:2)=='Cd') THEN
            cdbond(i,j)=1 ; cdbond(j,i)=1
          ENDIF
        ENDIF
      ENDDO
    ENDIF
  ENDDO
  IF (ncd==0) RETURN  

! Count # of connection of each Cd to other Cds
  DO i=1,ngr ; cdpst(i)=SUM(cdbond(i,:)) ; ENDDO

! check primary and tertiary Cd nodes
  nprim=0
  DO i=1,ngr
    IF (cdpst(i)==1) nprim=nprim+1
    IF (cdpst(i)>2) THEN
      mesg="tertiary Cd structure identified (and not allowed)"
      CALL stoperr(progname,mesg,chem)
    ENDIF    
  ENDDO
  IF (nprim==0) THEN
    mesg="unexpected cyclic Cd structure identified"
    CALL stoperr(progname,mesg,chem)
  ENDIF    

! get Cd chains (must start from primary Cd)
  DO i=1,ngr
    IF (cdpst(i)==1) THEN
      CALL gettrack(cdbond,i,ngr,ntr,track,trlen)
      IF (ntr>1) THEN
        mesg="unexpected branching identified on Cd structure"
        CALL stoperr(progname,mesg,chem)
      ENDIF
      IF (trlen(1)>4) THEN
        mesg="More than 4 Cd identified (and not allowed)"
        CALL stoperr(progname,mesg,chem)
      ENDIF
      IF (trlen(1)==3) THEN       ! check for >C=C=C< structure
        mesg=">C=C=C< structure identified and not allowed"
        CALL stoperr(progname,mesg,chem)
      ENDIF
      ncdtrack=ncdtrack+1                      ! add a new Cd track
      IF (ncdtrack>SIZE(cdtrack,1)) THEN
        mesg="maximum number of Cd tracks reached."
        CALL stoperr(progname,mesg,chem)
      ENDIF
      cdtrack(ncdtrack,1:mxcd)=track(1,1:mxcd) ! save the Cd track
      cdtracklen(ncdtrack)=trlen(1)
      cdpst(track(1,trlen(1)))=0               ! kill reverse track 
    ENDIF
  ENDDO

END SUBROUTINE alkenetrack

!=======================================================================
! PURPOSE: create tracks for -CO-O- groups. For example, the molecule
! CH3CO-O-CH2-O-CO-O as 2 tracks: 1) for the ester group -CO-O- group 
! and 2) for the carbonate group -O-CO-O.
! The routine output the number of track, the length of each track and
! the node belonging to each track.
!=======================================================================
SUBROUTINE estertrack(chem,bond,group,ngr,netrack,etracklen,etrack)
  USE keyparameter, ONLY: mxcp           
  USE toolbox, ONLY: stoperr
  IMPLICIT NONE

  CHARACTER(LEN=*),INTENT(IN) :: chem     ! formula
  CHARACTER(LEN=*),INTENT(IN) :: group(:) ! group without ring joining char
  INTEGER,INTENT(IN)  :: bond(:,:)        ! bond matrix
  INTEGER,INTENT(IN)  :: ngr              ! # of nodes in chem
  INTEGER,INTENT(OUT) :: netrack          ! # of (CO-O-)x track 
  INTEGER,INTENT(OUT) :: etracklen(:)     ! length of "ester" track i
  INTEGER,INTENT(OUT) :: etrack(:,:)      ! (i,j) node j for track i         
  
  INTEGER :: ebond(SIZE(bond,1),SIZE(bond,1))
  INTEGER :: i,j,nest,ntr,nprim,mxenod
  INTEGER :: epst(SIZE(bond,1))

  INTEGER :: track(mxcp,SIZE(group))
  INTEGER :: trlen(mxcp)
  CHARACTER(LEN=12),PARAMETER  :: progname='estertrack'
  CHARACTER(LEN=70) :: mesg

  netrack=0 ; etrack(:,:)=0 ; etracklen(:)=0
  ebond(:,:)=0  ;  mxenod=SIZE(etrack,2)
  
! ------------------------------
! Create and check the CO-O- tracks
! ------------------------------

! create a bond matrix for ester only might be bounded
  nest=0
  DO i=1,ngr
    IF (group(i)(1:3)=='-O-') THEN
      DO j=1,ngr
        IF (j==i) CYCLE
        IF (bond(i,j)/=0) THEN
          IF (group(j)=='CO')  THEN 
            ebond(i,j)=1 ; ebond(j,i)=1 ; nest=nest+1
          ELSEIF (group(j)=='CHO') THEN 
            ebond(i,j)=1 ; ebond(j,i)=1 ; nest=nest+1
          ENDIF
        ENDIF
      ENDDO
    ENDIF
  ENDDO
  IF (nest==0) RETURN  

! Count # of connection of each ester nodes to other ester nodes
  DO i=1,ngr ; epst(i)=SUM(ebond(i,:)) ; ENDDO

! check primary and secondary ester nodes
  nprim=0
  DO i=1,ngr
    IF (epst(i)==1) nprim=nprim+1
    IF (epst(i)>2) THEN
      mesg="tertiary ester structure identified (and not allowed)"
      CALL stoperr(progname,mesg,chem)
    ENDIF    
  ENDDO
  IF (nprim==0) THEN
    mesg="unexpected cyclic ester structure identified"
    CALL stoperr(progname,mesg,chem)
  ENDIF    

! get ester chains (must start from primary Cd)
  DO i=1,ngr
    IF (epst(i)==1) THEN
      CALL gettrack(ebond,i,ngr,ntr,track,trlen)
      IF (ntr>1) THEN
        mesg="unexpected branching identified on ester structure"
        CALL stoperr(progname,mesg,chem)
      ENDIF
      IF (trlen(1) > mxenod) THEN
        mesg="More than mx possible ester identified (and not allowed)"
        CALL stoperr(progname,mesg,chem)
      ENDIF
      netrack=netrack+1                      ! add a new ester track
      IF (netrack>SIZE(etrack,1)) THEN
        mesg="maximum number of Cd tracks reached."
        CALL stoperr(progname,mesg,chem)
      ENDIF
      etrack(netrack,1:mxenod)=track(1,1:mxenod) ! save the Cd track
      etracklen(netrack)=trlen(1)
      epst(track(1,trlen(1)))=0               ! kill reverse track 
    ENDIF
  ENDDO
END SUBROUTINE estertrack

! ======================================================================
! Purpose: return a set of table about the carbon skeleton and 
! functional moieties for the species provided as input. The 
! returned arguments are :
! - ngrp : number of functional group in chem
! - nodetype: table of character for type node:
!      'y' = carbonyl      'r' = aromatic        'o'= -O- node
!      'd' = Cd            'n' = others (i.e. normal)
! - mapfun(a,b,c): provide the number of function of type 'c' at position
!   (node) 'a'. index 'b' if for node type with 1=aliphatic, 2=cd and
!   3=aromatic. For example, the molecule CH2(OH)CdH=CdHCHO should 
!   have non zero values at the positions : mapfun(1,1,1)=1 and 
!   mapfun(4,2,9)=1
! - funflg(a): get the number of functional group at node a. For the 
!     example above, non-zero values are found at position 1 and 4, 
!     where it is set to 1.
! - alifun(i) : number of group of type "i" on a aliphatic carbon
! - cdfun(i) : number of group of type "i" on a cd carbon
! - arofun(i) : number of group of type "i" on a aromatic carbon
!   the index in those tables are :
!     index  1 -  5 : 'ROH  ','RNO2 ','RONO2','ROOH ','RF   '  
!     index  6 - 10 : 'RCl  ','RBr  ','RI   ','RCHO ','RCOR '
!     index 11 - 15 : 'RCOOH','COOOH','PAN  ','ROR  ','RCOOR'
!     index 16 - 20 : 'HCOOR','RCOF ','RCOCl','RCOBr','RCOI '
!     index 21      : 'CO(ONO2)'
! - tabester : provide the position of ester "couple" (i.e
!   the O and CO nodes. For example, the molecule CH3CO-O-CH2-O-COCH3 
!   has the following values: 
!            tabester(1,1)=3, tabester(1,2)=2 
!            tabester(2,1)=5, tabester(2,2)=6
! - ierr     : if not set to 0, then an error occured in the
!             subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE chemmap(chem,node,group,bond,ngrp,nodetype, &
                 alifun,cdfun,arofun,mapfun,funflg,   &
                 tabester,nfcd,nfcr,ierr)
  IMPLICIT NONE

  CHARACTER(LEN=*),INTENT(IN) :: chem   ! chem : chemical formula
  CHARACTER(LEN=*),INTENT(IN) :: group(:) ! table of groups
  INTEGER,INTENT(IN)          :: bond(:,:) ! bond matrix
  INTEGER,INTENT(IN)          :: node ! # of nodes is chem
  INTEGER,INTENT(OUT)         :: ngrp  
  CHARACTER(LEN=*),INTENT(OUT):: nodetype(:)
  REAL,INTENT(OUT)            :: alifun(:),cdfun(:),arofun(:)
  REAL,INTENT(OUT)            :: mapfun(:,:,:)
  INTEGER,INTENT(OUT)         :: funflg(:)
  INTEGER,INTENT(OUT)         :: tabester(:,:)  ! 1= -O- side, 2= CO side
  INTEGER,INTENT(OUT)         :: nfcd,nfcr
  INTEGER,INTENT(OUT)         :: ierr

  INTEGER                    :: i,j,k,lgr
  INTEGER                    :: nnod,tnod(SIZE(bond,1))
  INTEGER                    :: nf,ialpha,ialpha2,iy,rflg,dflg,yflg
  INTEGER                    :: ichecko, ichecky
  INTEGER                    :: ytab(2)
  INTEGER                    :: nester

  ierr=0  ;  nfcd=0  ;  nfcr=0  ;  ngrp=0  ;  alifun(:)=0  ;  cdfun(:)=0
  arofun(:)=0  ;  mapfun(:,:,:)=0  ;  nodetype(:)=' '  ;  funflg(:)=0
  nester=0  ;  tabester(:,:)=0  ;  lgr=LEN(group(1))

! get the type of nodes in groups
  DO i=1,node
    IF (group(i)(1:2)=='CO') THEN       ; nodetype(i)='y'
    ELSE IF (group(i)(1:3)=='CHO') THEN ; nodetype(i)='y'
    ElSE IF (group(i)(1:1)=='c') THEN   ; nodetype(i)='r'
    ELSE IF (group(i)(1:3)=='-O-') THEN ; nodetype(i)='o'
    ELSE IF (group(i)(1:2)=='Cd') THEN  ; nodetype(i)='d'
    ELSE                                ; nodetype(i)='n'
    ENDIF 
  ENDDO

! feed table of functions 
! --------------------------
!  1= -OH    ;  2= -NO2     ;  3= -ONO2   ; 4= -OOH    ;  5= -F     
!  6= -Cl    ;  7=  -Br     ;  8= -I      ; 9= -CHO    ; 10= -CO-   
! 11= -COOH  ; 12= -CO(OOH) ; 13= -PAN    ; 14= -O-    ; 15= R-COO-H 
! 16 = HCO-O-R; 17= -CO(F) ; 18= -CO(Cl) ;  19= -CO(Br) ; 20= -CO(I)
! 21 = R-CO(ONO2)

! ------- Alkohols (index 1) and Acids (index 11) ----------
  IF (INDEX(chem,'(OH)')/=0) THEN
    DO i=1, node
      IF (INDEX(group(i),'(OH)')/=0) THEN
        nf=0            
        DO j=1,lgr-3
          IF (group(i)(j:j+3)=='(OH)') nf=nf+1 
        ENDDO
        IF (nodetype(i)=='n') THEN       ! alcohol aliphatic
          alifun(1)=alifun(1)+nf
          mapfun(i,1,1)=nf
          ngrp=ngrp+nf
        ELSE IF (nodetype(i)=='r') THEN  ! phenols 
          arofun(1)=arofun(1)+nf
          mapfun(i,3,1)=nf
          ngrp=ngrp+nf
          nfcr=nfcr+nf
        ELSE IF (nodetype(i)=='d') THEN  ! =Cd-OH (not expected, since enol) 
          cdfun(1)=cdfun(1)+nf
          mapfun(i,2,1)=nf
          ngrp=ngrp+nf
          nfcd=nfcd+nf
        ELSE IF (nodetype(i)=='y') THEN  ! Carboxylic acid
          CALL nodmap(bond,i,node,2,nnod,tnod)
          IF (nnod>1) THEN 
            DO j=1,node ; WRITE(*,*) 'group',j,'-',group(j) ; ENDDO
            WRITE(6,*) '-- error --, a unique C is expected in'
            WRITE(6,*) 'alpha position of a carboxylic group'
            STOP "in chemmap"
          ENDIF
          ialpha=tnod(1)
          IF (nodetype(ialpha)=='r')THEN         ! CO(OH) aromatic 
            arofun(11)=arofun(11)+1
            mapfun(i,3,11)=1
            ngrp=ngrp+1
            nfcr=nfcr+1
          ELSE IF (nodetype(ialpha)=='d')THEN    ! CO(OH) on Cd
            cdfun(11)=cdfun(11)+1
            mapfun(i,2,11)=1
            ngrp=ngrp+1
            nfcd=nfcd+1
          ELSE                                    ! CO(OH) aliphatic
            alifun(11)=alifun(11)+1
            mapfun(i,1,11)=1
            ngrp=ngrp+1
          ENDIF 
        ENDIF
      ENDIF
    ENDDO 
  ENDIF 

! ----------- Nitro (index 2) -----------------
  IF (INDEX(chem,'(NO2)')/=0) THEN
    DO i = 1, node 
      IF (INDEX(group(i),'(NO2)')/=0) THEN
        nf=0
        DO j=1,lgr-4
          IF (group(i)(j:j+4)=='(NO2)') nf=nf+1 
        ENDDO
        IF (nodetype(i)=='n') THEN     ! NO2 aliphatic
          alifun(2)=alifun(2)+nf
          mapfun(i,1,2)=nf
          ngrp=ngrp+nf
        ELSE IF (nodetype(i)=='y') THEN     ! NO2 aliphatic
          alifun(2)=alifun(2)+nf
          mapfun(i,1,2)=nf
          ngrp=ngrp+nf
        ELSE IF (nodetype(i)=='d') THEN     ! NO2 on Cd 
          cdfun(2)=cdfun(2)+nf
          mapfun(i,2,2)=nf
          ngrp=ngrp+nf
          nfcd=nfcd+nf
        ELSE IF (nodetype(i)=='r') THEN     ! NO2 on aromatic 
          arofun(2)=arofun(2)+nf
          mapfun(i,3,2)=nf
          ngrp=ngrp+nf
          nfcr=nfcr+nf
        ELSE
          WRITE(6,*) '-- error --, a (NO2) group is borne by an '
          WRITE(6,*) ' unexpected group in chem :'
          WRITE(6,*) TRIM(chem)
          STOP "in chemmap"
        ENDIF
      ENDIF
    ENDDO
  ENDIF

! ----------- Nitrate (index 3) -----------------
  IF (INDEX(chem,'(ONO2)')/=0) THEN
    DO i = 1, node 
      IF (INDEX(group(i),'(ONO2)')/=0) THEN

!        IF (INDEX(group(i),'CO(ONO2)')/=0) GOTO 120  
        IF (INDEX(group(i),'CO(ONO2)')/=0) CYCLE ! ba: unexpected but out of IF in original version
        nf=0
        DO j=1,lgr-5
          IF (group(i)(j:j+5)=='(ONO2)') nf=nf+1 
        ENDDO

        IF (nodetype(i)=='n') THEN          ! ONO2 aliphatic
          alifun(3)=alifun(3)+nf
          mapfun(i,1,3)=nf
          ngrp=ngrp+nf
        ELSE IF (nodetype(i)=='d') THEN     ! ONO2 on Cd 
          cdfun(3)=cdfun(3)+nf
          mapfun(i,2,3)=nf
          ngrp=ngrp+nf
          nfcd=nfcd+nf
        ELSE IF (nodetype(i)=='r') THEN     ! ONO2 on aromatic 
          arofun(3)=arofun(3)+nf
          mapfun(i,3,3)=nf
          ngrp=ngrp+nf
         nfcr=nfcr+nf
        ELSE
          WRITE(6,*) '-- error --, a (ONO2) group is borne by an '
          WRITE(6,*) 'unexpected group in chem :'
          WRITE(6,*) TRIM(chem)
          STOP "in chemmap"
        ENDIF
      ENDIF
    ENDDO
  ENDIF

! ------- hydroperoxydes (index 4) and peracids (index 12) ----------
  IF (INDEX(chem,'(OOH)')/=0) THEN
    DO i=1, node
      IF (INDEX(group(i),'(OOH)')/=0) THEN
        nf=0            
        DO j=1,lgr-4
          IF (group(i)(j:j+4)=='(OOH)') nf=nf+1 
        ENDDO
        IF (nodetype(i)=='n') THEN       ! hydroperoxyde aliphatic
          alifun(4)=alifun(4)+nf
          mapfun(i,1,4)=nf
          ngrp=ngrp+nf
        ELSE IF (nodetype(i)=='r') THEN  ! aromaric -OOH 
          arofun(4)=arofun(4)+nf
          mapfun(i,3,4)=nf
          ngrp=ngrp+nf
          nfcr=nfcr+nf
        ELSE IF (nodetype(i)=='d') THEN  ! =Cd-OOH (not expected) 
          cdfun(4)=cdfun(4)+nf
          mapfun(i,2,4)=nf
          ngrp=ngrp+nf
          nfcd=nfcd+nf
        ELSE IF (nodetype(i)=='y') THEN  ! peracid acid
          CALL nodmap(bond,i,node,2,nnod,tnod)
          IF (nnod>1) THEN 
            DO j=1,node  ;  WRITE(*,*) 'group',j,'-',group(j) ; ENDDO
            WRITE(6,*) '-- error --, a unique C is expected in'
            WRITE(6,*) 'alpha position of a peracid group'
            WRITE(6,*) TRIM(chem)
            STOP "in chemmap"
          ENDIF
          ialpha=tnod(1)
          IF (nodetype(ialpha)=='r')THEN         ! CO(OOH) aromatic 
            arofun(12)=arofun(12)+1
            mapfun(i,3,12)=1
            ngrp=ngrp+1
          ELSE IF (nodetype(ialpha)=='d')THEN    ! CO(OOH) on Cd
            cdfun(12)=cdfun(12)+1
            mapfun(i,2,12)=1
            ngrp=ngrp+1
          ELSE                                    ! CO(OOH) aliphatic
            alifun(12)=alifun(12)+1
            mapfun(i,1,12)=1
            ngrp=ngrp+1
          ENDIF 
        ENDIF
      ENDIF
    ENDDO 
  ENDIF

! ------- fluroro (index 5) and fluoro acyl (index 17) ----------
  IF (INDEX(chem,'(F)')/=0) THEN
    DO i=1, node
      IF (INDEX(group(i),'(F)')/=0) THEN
        nf=0            
        DO j=1,lgr-2
          IF (group(i)(j:j+2)=='(F)') nf=nf+1 
        ENDDO
        IF (nodetype(i)=='n') THEN       ! F aliphatic
          alifun(5)=alifun(5)+nf
          mapfun(i,1,5)=nf
          ngrp=ngrp+nf
        ELSE IF (nodetype(i)=='r') THEN  ! aromaric F 
          arofun(5)=arofun(5)+nf
          mapfun(i,3,5)=nf
          ngrp=ngrp+nf
          nfcr=nfcr+nf
        ELSE IF (nodetype(i)=='d') THEN  ! =Cd-F 
          cdfun(5)=cdfun(5)+nf
          mapfun(i,2,5)=nf
          ngrp=ngrp+nf
          nfcd=nfcd+nf
        ELSE IF (nodetype(i)=='y') THEN  ! fluoro acyl
          CALL nodmap(bond,i,node,2,nnod,tnod)
          IF (nnod>1) THEN 
            DO j=1,node ; WRITE(*,*) 'group',j,'-',group(j) ;  ENDDO
            WRITE(6,*) '-- error --, a unique C is expected in'
            WRITE(6,*) 'alpha position of a fluoro acyl group'
            WRITE(6,*) TRIM(chem)
            STOP "in chemmap"
          ENDIF
          ialpha=tnod(1)
          IF (nodetype(ialpha)=='r')THEN         ! CO(F) aromatic 
            arofun(17)=arofun(17)+1
            mapfun(i,3,17)=1
            ngrp=ngrp+1
          ELSE IF (nodetype(ialpha)=='d')THEN    ! CO(F) on Cd
            cdfun(17)=cdfun(17)+1
            mapfun(i,2,17)=1
            ngrp=ngrp+1
          ELSE                                    ! CO(F) aliphatic
            alifun(17)=alifun(17)+1
            mapfun(i,1,17)=1
            ngrp=ngrp+1
          ENDIF 
        ENDIF
      ENDIF
    ENDDO 
  ENDIF

! ------- chloro (index 6) and chloro acyl (index 18) ----------
  IF (INDEX(chem,'(Cl)')/=0) THEN
    DO i=1, node
      IF (INDEX(group(i),'(Cl)')/=0) THEN
        nf=0            
        DO j=1,lgr-3
          IF (group(i)(j:j+3)=='(Cl)') nf=nf+1 
        ENDDO
        IF (nodetype(i)=='n') THEN       ! Cl aliphatic
          alifun(6)=alifun(6)+nf
          mapfun(i,1,6)=nf
          ngrp=ngrp+nf
        ELSE IF (nodetype(i)=='r') THEN  ! aromaric Cl
          arofun(6)=arofun(6)+nf
          mapfun(i,3,6)=nf
          ngrp=ngrp+nf
          nfcr=nfcr+nf
        ELSE IF (nodetype(i)=='d') THEN  ! =Cd-Cl  
          cdfun(6)=cdfun(6)+nf
          mapfun(i,2,6)=nf
          ngrp=ngrp+nf
          nfcd=nfcd+nf
        ELSE IF (nodetype(i)=='y') THEN  ! Chloro acyl
          CALL nodmap(bond,i,node,2,nnod,tnod)
          IF (nnod>1) THEN 
            DO j=1,node ; WRITE(*,*) 'group',j,'-',group(j) ; ENDDO
            WRITE(6,*) '-- error --, a unique C is expected in'
            WRITE(6,*) 'alpha position of a chloro acyl group'
            WRITE(6,*) TRIM(chem)
            STOP "in chemmap"
          ENDIF
          ialpha=tnod(1)
          IF (nodetype(ialpha)=='r')THEN         ! CO(Cl) aromatic 
            arofun(18)=arofun(18)+1
            mapfun(i,3,18)=1
            ngrp=ngrp+1
          ELSE IF (nodetype(ialpha)=='d')THEN    ! CO(Cl) on Cd
            cdfun(18)=cdfun(18)+1
            mapfun(i,2,18)=1
            ngrp=ngrp+1
          ELSE                                    ! CO(Cl) aliphatic
            alifun(18)=alifun(18)+1
            mapfun(i,1,18)=1
            ngrp=ngrp+1
          ENDIF 
        ENDIF
      ENDIF
    ENDDO 
  ENDIF

! ------- Bromo (index 7) and bromo acyl (index 19) ----------
  IF (INDEX(chem,'(Br)')/=0) THEN
    DO i=1, node
      IF (INDEX(group(i),'(Br)')/=0) THEN
        nf=0            
        DO j=1,lgr-3
          IF (group(i)(j:j+3)=='(Br)') nf=nf+1 
        ENDDO
        IF (nodetype(i)=='n') THEN       ! Br aliphatic
          alifun(7)=alifun(7)+nf
          mapfun(i,1,7)=nf
          ngrp=ngrp+nf
        ELSE IF (nodetype(i)=='r') THEN  ! aromatic Br
          arofun(7)=arofun(7)+nf
          mapfun(i,3,7)=nf
          ngrp=ngrp+nf
          nfcr=nfcr+nf
        ELSE IF (nodetype(i)=='d') THEN  ! =Cd-Br  
          cdfun(7)=cdfun(7)+nf
          mapfun(i,2,7)=nf
          ngrp=ngrp+nf
          nfcd=nfcd+nf
        ELSE IF (nodetype(i)=='y') THEN  ! Bromo acyl
          CALL nodmap(bond,i,node,2,nnod,tnod)
          IF (nnod>1) THEN 
            DO j=1,node  ; WRITE(*,*) 'group',j,'-',group(j)  ;  ENDDO
            WRITE(6,*) '-- error --, a unique C is expected in'
            WRITE(6,*) 'alpha position of a bromo acyl group'
            WRITE(6,*) chem(1:70)
            STOP "in chemmap"
          ENDIF
          ialpha=tnod(1)
          IF (nodetype(ialpha)=='r')THEN         ! CO(Br) aromatic 
            arofun(19)=arofun(19)+1
            mapfun(i,3,19)=1
            ngrp=ngrp+1
          ELSE IF (nodetype(ialpha)=='d')THEN    ! CO(Br) on Cd
            cdfun(19)=cdfun(19)+1
            mapfun(i,2,19)=1
            ngrp=ngrp+1
          ELSE                                    ! CO(Br) aliphatic
            alifun(19)=alifun(19)+1
            mapfun(i,1,19)=1
            ngrp=ngrp+1
          ENDIF 
        ENDIF
      ENDIF
    ENDDO 
  ENDIF

! ------- iodo (index 8) and iodo acyl (index 20) ----------
  IF (INDEX(chem,'(I)')/=0) THEN
    DO i=1, node
      IF (INDEX(group(i),'(I)')/=0) THEN
        nf=0            
        DO j=1,lgr-2
          IF (group(i)(j:j+2)=='(I)') nf=nf+1 
        ENDDO
        IF (nodetype(i)=='n') THEN       ! I aliphatic
          alifun(8)=alifun(8)+nf
          mapfun(i,1,8)=nf
          ngrp=ngrp+nf
        ELSE IF (nodetype(i)=='r') THEN  ! aromaric I 
          arofun(8)=arofun(8)+nf
          mapfun(i,3,8)=nf
          ngrp=ngrp+nf
          nfcr=nfcr+nf
        ELSE IF (nodetype(i)=='d') THEN  ! =Cd-I 
          cdfun(8)=cdfun(8)+nf
          mapfun(i,2,8)=nf
          ngrp=ngrp+nf
          nfcd=nfcd+nf
        ELSE IF (nodetype(i)=='y') THEN  ! fluoro acyl
          CALL nodmap(bond,i,node,2,nnod,tnod)
          IF (nnod>1) THEN 
            DO j=1,node ; WRITE(*,*) 'group',j,'-',group(j)  ;  ENDDO
            WRITE(6,*) '-- error --, a unique C is expected in'
            WRITE(6,*) 'alpha position of a Iodo acyl group'
            WRITE(6,*) chem(1:70)
            STOP "in chemmap"
          ENDIF
          ialpha=tnod(1)
          IF (nodetype(ialpha)=='r')THEN         ! CO(I) aromatic 
            arofun(20)=arofun(20)+1
            mapfun(i,3,20)=1
            ngrp=ngrp+1
            nfcr=nfcr+1
          ELSE IF (nodetype(ialpha)=='d')THEN    ! CO(I) on Cd
            cdfun(20)=cdfun(20)+1
            mapfun(i,2,20)=1
            ngrp=ngrp+1
            nfcd=nfcd+1
          ELSE                                    ! CO(I) aliphatic
            alifun(20)=alifun(20)+1
            mapfun(i,1,20)=1
            ngrp=ngrp+1
          ENDIF 
        ENDIF
      ENDIF
    ENDDO 
  ENDIF

! ---------- ester (index 15 and 16) -------------
  IF (INDEX(chem,'-O')/=0) THEN
    esterloop: DO i = 1, node
      IF (group(i)(1:3)=='-O-') THEN
        CALL nodmap(bond,i,node,2,nnod,tnod)
        rflg=0  ;  dflg=0  ;  yflg=0
        DO j=1,nnod
          ialpha=tnod(j)
          IF (nodetype(ialpha)=='y') THEN   ! RCO-O-R function
            ichecky=0
            DO k=1,4
              IF (tabester(k,2)==ialpha) ichecky=1   ! carbonyl already used
            ENDDO
            IF (ichecky==0) THEN
              yflg=yflg+1  ;  ytab(yflg)=ialpha
            ENDIF
          ENDIF
        ENDDO
        IF (yflg==0) CYCLE  esterloop               ! simple ether

! simple ester 
        IF (yflg/=0) THEN
          nester=nester+1
          IF (nester>4) THEN
            PRINT*, "in chemmap, nester > 4"
            STOP "in chemmap"
          ENDIF  
          iy=ytab(1)
          tabester(nester,1)=i  ;  tabester(nester,2)=iy
          DO j=1,2
            IF (tnod(j)/=ytab(1)) THEN
              ialpha=tnod(j)   ! second side of the -O- (first is iy)
            ENDIF
          ENDDO
          IF (group(iy)(1:3)=='CHO') THEN    ! HCO-O-R
            IF (nodetype(ialpha)=='r') THEN       ! aromatic CHO-O-R
              arofun(16)=arofun(16)+1
              mapfun(i,3,16)=1
              mapfun(iy,3,16)=1
              ngrp=ngrp+1
              nfcr=nfcr+1
            ELSE IF (nodetype(ialpha)=='d') THEN  ! =C-O-CHO
              cdfun(16)=cdfun(16)+1
              mapfun(i,2,16)=1
              mapfun(iy,2,16)=1
              ngrp=ngrp+1
              nfcd=nfcd+1
            ELSE                                    ! R-O-CHO
              alifun(16)=alifun(16)+1
              mapfun(i,1,16)=1
              mapfun(iy,1,16)=1
              ngrp=ngrp+1
            ENDIF
          ELSE IF (group(iy)(1:3)=='CO ') THEN    ! RCO-O-R
            CALL nodmap(bond,iy,node,2,nnod,tnod)
            DO j=1,2
              IF (tnod(j)/=i) ialpha2=tnod(j)  ! ialpha2 is the node next to the CO of the ester
            ENDDO
            rflg=0  ;  dflg=0

! structure is ialpha2-CO-O-ialpha
            IF (nodetype(ialpha)=='r') rflg=rflg+1       
            IF (nodetype(ialpha)=='d') dflg=dflg+1       
            IF (nodetype(ialpha2)=='r') rflg=rflg+1       
            IF (nodetype(ialpha2)=='d') dflg=dflg+1       

            IF (rflg/=0) THEN                    ! aromatic ester
              arofun(15)=arofun(15)+1.
              mapfun(i,3,15)=1.
              mapfun(iy,3,15)=1.
              nfcr=nfcr+rflg
            ELSE IF (dflg/=0) THEN                    ! =C-CO-O-
              cdfun(15)=cdfun(15)+1.
              mapfun(i,2,15)=1.
              mapfun(iy,2,15)=1.
              nfcd=nfcd+dflg
            ELSE                       ! R-COO-R
              alifun(15)=alifun(15)+1.
              mapfun(i,1,15)=1.
              mapfun(iy,1,15)=1.
            ENDIF
            ngrp=ngrp+1
          ENDIF
        ENDIF

      ENDIF 
        
    ENDDO esterloop
  ENDIF

! ------------Aldehydes (index 9) --------- 
  IF (INDEX(chem,'CHO')/=0) THEN
    dloop: DO i = 1, node
      IF (group(i)(1:3)=='CHO') THEN
        CALL nodmap(bond,i,node,2,nnod,tnod)
        IF (nnod/=1) THEN
          WRITE(6,*) '-- error --, a unique C is expected in'
          WRITE(6,*) 'alpha position of a CHO  group'
          WRITE(6,*) TRIM(chem)
          STOP "in chemmap"
        ENDIF
        ialpha=tnod(1)
        IF (nodetype(ialpha)=='o') THEN   ! HCO-O-R function
           ichecko=0 ! check if the ether is already involved in an ester
           ichecky=0 ! check if the carbonyl is already involved in an ester
           DO k=1,4
             IF (tabester(k,1)==ialpha) ichecko=1   ! ether already used
             IF (tabester(k,2)==i) ichecky=1        ! carbonyl already used
           ENDDO
           IF (ichecky==1) CYCLE dloop ! ether already involved
           IF (ichecko==0) CYCLE dloop ! carbonyl that must be an ester
        ENDIF  ! if that point is reached then must be counted as aldehyde
	    
        IF (nodetype(ialpha)=='r') THEN       ! aromatic aldehyde
           arofun(9)=arofun(9)+1
           mapfun(i,3,9)=1
           ngrp=ngrp+1
           nfcr=nfcr+1
        ELSE IF (nodetype(ialpha)=='d') THEN  ! =C-CHO
           cdfun(9)=cdfun(9)+1
           mapfun(i,2,9)=1
           ngrp=ngrp+1
           nfcd=nfcd+1
        ELSE                                  ! R-CHO
           alifun(9)=alifun(9)+1
           mapfun(i,1,9)=1
           ngrp=ngrp+1
        ENDIF
      ENDIF
    ENDDO dloop
  ENDIF

! ---------- ketones (index 10) -------------
  IF (INDEX(chem,'CO')/=0) THEN
    kloop: DO i = 1, node
      IF (group(i)(1:3)=='CO ') THEN
        CALL nodmap(bond,i,node,2,nnod,tnod)
        IF (nnod/=2) THEN
          WRITE(6,*) '-- error --, only 2 C is expected in'
          WRITE(6,*) 'alpha position of a -CO-  group'
          WRITE(6,*) TRIM(chem),'  nnod=',nnod
          WRITE(6,*) bond(1:node,1:node)
          STOP "in chemmap"
        ENDIF
        rflg=0  ;  dflg=0
        DO j=1,nnod
          ialpha=tnod(j)
          IF (nodetype(ialpha)=='o') THEN   ! RCO-O-R function
            ichecko=0 ! check if the ether is already involved in an ester
            ichecky=0 ! check if the carbonyl is already involved in an ester
            DO k=1,4
              IF (tabester(k,1)==ialpha) ichecko=1   ! ether already used
              IF (tabester(k,2)==i) ichecky=1   ! carbonyl already used
            ENDDO
            IF (ichecky==1) CYCLE kloop  ! carbonyl already involved
            IF (ichecko==0) CYCLE kloop  ! ether that must be an ester
          ENDIF  ! if that point is reached then must be counted as ketone
          IF (nodetype(ialpha)=='r') rflg=rflg+1       
          IF (nodetype(ialpha)=='d') dflg=dflg+1       
        ENDDO
        IF (rflg/=0) THEN                    ! aromatic ketone
          arofun(10)=arofun(10)+1.
          mapfun(i,3,10)=1.
          nfcr=nfcr+rflg
        ELSE IF (dflg/=0) THEN                    ! =C-CO-R
          cdfun(10)=cdfun(10)+1.
          mapfun(i,2,10)=1.
          nfcd=nfcd+dflg
        ELSE                       ! R-CO-R
          alifun(10)=alifun(10)+1.
          mapfun(i,1,10)=1.
        ENDIF
        ngrp=ngrp+1
      ENDIF
    ENDDO kloop
  ENDIF

! ------------- PAN (index 13) -----------------
! BA- July 2020. The group CO(ONO2) is ignored in the current version.
! Some work was likely initiated (see index 21) but was apparently not
! finished - these lines being commented. All this need to be revisited.
! As a preliminary patch, I added CO(ONO2) to the PAN group. Chemmap is 
! only called by gromhe ... so that should not create issues. Note that
! work was done also on GECKO July 2020 to avoid the production of  
! these structures.

!baba  IF (INDEX(chem,'CO(OONO2')/=0) THEN
  IF ((INDEX(chem,'CO(OONO2')/=0).OR.(INDEX(chem,'CO(ONO2')/=0)) THEN
    DO i = 1, node
      IF ( (group(i)(1:9)=='CO(OONO2)').OR. &
           (group(i)(1:8)=='CO(ONO2)') ) THEN
       CALL nodmap(bond,i,node,2,nnod,tnod)
       IF (nnod/=1) THEN
         WRITE(6,*) '-- error --, a unique C is expected in'
         WRITE(6,*) 'alpha position of a CO(OONO2)  group'
         WRITE(6,*) TRIM(chem)
         STOP "in chemmap"
       ENDIF
       ialpha=tnod(1)
       IF (nodetype(ialpha)=='o') THEN   ! R-O-CO(OONO2) function
          WRITE(6,*) '-- error --,  -O-CO(OONO2) group is unexpected'
          WRITE(6,*) TRIM(chem)
!BABA july 2020           STOP "in chemmap"
          nodetype(ialpha)='n' ! patch BABA july 2020 (just to make it work)
       ENDIF
       IF (nodetype(ialpha)=='r') THEN       ! aromatic PAN
          arofun(13)=arofun(13)+1
          mapfun(i,3,13)=1
          ngrp=ngrp+1
          nfcr=nfcr+1
       ELSE IF (nodetype(ialpha)=='d') THEN  ! =C-CO(OONO2)
          cdfun(13)=cdfun(13)+1
          mapfun(i,2,13)=1
          ngrp=ngrp+1
          nfcd=nfcd+1
       ELSE                                  ! R-CO(OONO2)
          alifun(13)=alifun(13)+1
          mapfun(i,1,13)=1
          ngrp=ngrp+1
       ENDIF
      ENDIF
    ENDDO
  ENDIF

! ---------- ether (index 14) -------------
  IF (INDEX(chem,'-O')/=0) THEN
    etherloop: DO i = 1, node
      IF (group(i)(1:3)=='-O-') THEN
        CALL nodmap(bond,i,node,2,nnod,tnod)
        IF (nnod/=2) THEN
          WRITE(6,*) '-- error --, only 2 C is expected in'
          WRITE(6,*) 'alpha position of a -O-  group'
          WRITE(6,*) TRIM(chem)
          STOP "in chemmap"
        ENDIF
        rflg=0  ;  dflg=0
        DO j=1,nnod
          ialpha=tnod(j)
          IF (nodetype(ialpha)=='y') THEN   ! RCO-O-R function
            ichecko=0 ! check if the ether is already involved in an ester
            ichecky=0 ! check if the carbonyl is already involved in an ester
            DO k=1,4
              IF (tabester(k,1)==i) ichecko=1   ! ether already used
              IF (tabester(k,2)==ialpha) ichecky=1   ! carbonyl already used
            ENDDO
            IF (ichecko==1) CYCLE etherloop  ! ether already involved in an ester
            IF (ichecky==0) CYCLE etherloop  ! carbonyl that must be an ester
          ENDIF  ! if that point is reached then must be counted as ether

          IF (nodetype(ialpha)=='r') rflg=rflg+1       
          IF (nodetype(ialpha)=='d') dflg=dflg+1       
        ENDDO
        IF (rflg/=0) THEN                    ! aromatic ether
          arofun(14)=arofun(14)+1.
          mapfun(i,3,14)=1.
          nfcr=nfcr+rflg

          IF (rflg>1) THEN
            ierr=1
            RETURN                           ! return error
          ENDIF
        ELSE IF (dflg/=0) THEN               ! =C-O-R
          cdfun(14)=cdfun(14)+1.
          mapfun(i,2,14)=1.
          nfcd=nfcd+dflg
        ELSE                      ! R-O-R
          alifun(14)=alifun(14)+1.
          mapfun(i,1,14)=1.
        ENDIF
        ngrp=ngrp+1
      ENDIF
    ENDDO etherloop
  ENDIF

!baba !------------- CO(ONO2) (index 21) --------------
!baba       IF (INDEX(chem,'CO(ONO2')/=0) THEN
!baba          WRITE(6,*) '-- error --, a CO(ONO2) group is unexpected'
!baba          WRITE(6,*) TRIM(chem)
!baba !         STOP "in chemmap"
!baba !         ierr=1
!baba !         RETURN                           ! return error
!baba       
!baba !      DO i = 1, node
!baba !        IF (group(i)(1:8)=='CO(ONO2)') THEN
!baba !         CALL nodmap(bond,i,node,2,nnod,tnod)
!baba !         IF (nnod/=1) THEN
!baba !            WRITE(6,*) '-- error --, a unique C is expected in'
!baba !            WRITE(6,*) 'alpha position of a CO(ONO2)  group'
!baba !            WRITE(6,*) chem(1:70)
!baba !            STOP "in chemmap"
!baba !         ENDIF
!baba !         ialpha=tnod(1)
!baba !         IF (nodetype(ialpha)=='o') THEN   ! R-O-CO(ONO2) function
!baba !            WRITE(6,*) '-- error --,  in'
!baba !            WRITE(6,*) '-O-CO(ONO2) group is unexpected'
!baba !            WRITE(6,*) chem(1:70)
!baba !            STOP "in chemmap"
!baba !         ENDIF
!baba !         IF (nodetype(ialpha)=='r') THEN       ! aromatic
!baba !            arofun(21)=arofun(21)+1
!baba !            mapfun(i,3,21)=1
!baba !            ngrp=ngrp+1
!baba !         ELSE IF (nodetype(ialpha)=='d') THEN  ! =C-CO(ONO2)
!baba !            cdfun(21)=cdfun(21)+1
!baba !c            alifun(13)=alifun(13)+1
!baba !            mapfun(i,2,21)=1
!baba !c            mapfun(i,1,13)=1
!baba !            ngrp=ngrp+1
!baba !         ELSE                                    ! R-CO(ONO2)
!baba !            alifun(21)=alifun(21)+1
!baba !            mapfun(i,1,21)=1
!baba !            ngrp=ngrp+1
!baba !         ENDIF
!baba !        ENDIF
!baba !      ENDDO
!baba   ENDIF


! set the table telling if a "functional group" is available at a given node
  DO i=1,node
    DO j=1,3
      DO k=1,20
        IF (mapfun(i,j,k)/=0.) THEN
          IF (mapfun(i,j,k).lt.1.) THEN
            funflg(i)=funflg(i)+1
          ELSE
            funflg(i)=funflg(i)+INT(mapfun(i,j,k))
          ENDIF
        ENDIF
      ENDDO
    ENDDO
  ENDDO

END SUBROUTINE chemmap


END MODULE mapping
