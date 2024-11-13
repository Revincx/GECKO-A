MODULE simpoltool
IMPLICIT NONE
CONTAINS

! SUBROUTINE simpol(chem,bond,group,nring,logPvap,latentheat)
! SUBROUTINE load_simpoldat()
! SUBROUTINE simpolgrp(chem,group,bond,nring,simpgroup)

! ======================================================================
! Purpose : Compute the vapor pressures and enthalpies of vaporization
! of the species provided as input based on the simpol SAR.
! Ref. for the SAR: Pankow, J. F. and Asher, W. E.: SIMPOL.1: a simple 
! group contribution method for predicting vapor pressures and 
! enthalpies of vaporization of multifunctional organic compounds, 
! Atmos. Chem. Phys., 2773-2796, 2008. 
! https://doi.org/10.5194/acp-8-2773-2008
! https://www.atmos-chem-phys.net/8/2773/2008/
! ======================================================================
SUBROUTINE simpol(chem,bond,group,nring,logPvap,latentheat)
  USE database, ONLY: psimpgrp298,hsimpgrp298,p0simpgrp298,h0simpgrp298
  IMPLICIT NONE

  CHARACTER(LEN=*),INTENT(IN) :: chem
  CHARACTER(LEN=*),INTENT(IN) :: group(:) ! group, without ring joining char
  INTEGER,INTENT(IN) :: bond(:,:)
  INTEGER,INTENT(IN) :: nring
  REAL,INTENT(OUT)   :: logPvap
  REAL,INTENT(OUT)   :: latentheat

  REAL,PARAMETER :: temp=298. ! ref temperature for Psat and latent heat
  INTEGER :: j
  INTEGER :: simpgroup(30)
  
! get simpol groups
  CALL simpolgrp(chem,group,bond,nring,simpgroup)

! add constant term
  logPvap = p0simpgrp298  ;  latentheat = h0simpgrp298 

! add group contribution
  DO j=1,30
    IF (simpgroup(j)/=0) THEN
      logPvap = logPvap + simpgroup(j)*psimpgrp298(j) 
      latentheat = latentheat + simpgroup(j)*hsimpgrp298(j)
    ENDIF
  ENDDO

  latentheat = -2.303 * 8.314 * latentheat
END SUBROUTINE simpol

! ======================================================================
! Purpose: read and store the values of the groups for the Simpol SAR 
! ======================================================================
SUBROUTINE load_simpoldat()
  USE keyparameter, ONLY: tfu1,dirgecko
  USE database, ONLY: psimpgrp298,hsimpgrp298,p0simpgrp298,h0simpgrp298
  IMPLICIT NONE

  REAL,PARAMETER :: temp=298. ! ref temperature for Psat and latent heat
  INTEGER :: j,ilin,ierr
  CHARACTER(LEN=512) filename
  CHARACTER(LEN=150) line
  REAL :: bk0(4)        ! parameter of the simpol SAR ("constant" group)
  REAL :: bkgr(30,4)    ! parameter of the simpol SAR
  REAL :: simpdat(SIZE(bkgr,1)+1,4)  ! parameter of the simpol SAR
  
  filename=TRIM(dirgecko)//'DATA/simpol.dat'
  OPEN(tfu1,FILE=filename, FORM='FORMATTED', STATUS='OLD', IOSTAT=ierr)
  IF (ierr/=0) THEN
    WRITE(6,*) '--error--, in load_simpoldat while trying to open file:',TRIM(filename)
    STOP "in load_simpoldat"
  ENDIF

! read the data
  ilin=0
  rdloop: DO
    READ(tfu1,'(a)') line
    IF (line(1:1)=='!') CYCLE rdloop
    IF (line(1:3)=='END') EXIT rdloop
    
    ilin=ilin+1
    READ(line,*,IOSTAT=ierr) (simpdat(ilin,j),j=1,4)
    IF (ierr/=0) THEN
      WRITE(6,*) '--error--, in load_simpoldat while reading line:',TRIM(line)
      STOP "in load_simpoldat"
    ENDIF
  ENDDO rdloop
  IF (ilin/=SIZE(simpdat,1)) THEN
    WRITE(6,*) '--error--, unexpected # of data in :',TRIM(filename)
    STOP "in load_simpoldat"
  ENDIF
  CLOSE(tfu1)

  bk0(:)=simpdat(1,:)        ! constant term
  bkgr(1:,:)=simpdat(2:,:)   ! group term
  
! store the overall group contribution @ the reference T (here 298K)
  p0simpgrp298 = (bk0(1)/temp) + bk0(2) + &
                 (bk0(3)*temp) + bk0(4)*log(temp)
  h0simpgrp298 = bk0(1) - bk0(3)*temp*temp - bk0(4)*temp

  DO j=1,30
    psimpgrp298(j) = (bkgr(j,1)/temp) + bkgr(j,2) + &
                     (bkgr(j,3)*temp) + bkgr(j,4)*log(temp)
    hsimpgrp298(j) = bkgr(j,1) - bkgr(j,3)*temp*temp - bkgr(j,4)*temp
  ENDDO

END SUBROUTINE load_simpoldat

! ======================================================================
! Purpose: provide the simpol groups for the simpol SAR to estimate
! the vapor pressure of organic species.
! ======================================================================
SUBROUTINE simpolgrp(chem,group,bond,nring,simpgroup)
  USE keyparameter, ONLY: mxtrk,mxlcd,mxlest
  USE ringtool, ONLY: findring      
  USE mapping, ONLY: estertrack,alkenetrack
  USE toolbox, ONLY: stoperr,countstring
  IMPLICIT NONE

  CHARACTER(LEN=*),INTENT(IN) :: chem
  CHARACTER(LEN=*),INTENT(IN) :: group(:) ! group without ring joining char
  INTEGER,INTENT(IN)    :: bond(:,:)
  INTEGER,INTENT(IN)    :: nring
  INTEGER,INTENT(OUT)   :: simpgroup(:)   

! internal
  CHARACTER(LEN=LEN(group(1))) :: neigh(SIZE(bond,1),4)
  CHARACTER(LEN=LEN(group(1))) :: tgroup(SIZE(group))
  INTEGER :: i,j,nester,ifirst,ipero,inod,aronitro,nn,nitroester,icd,ico
  INTEGER :: ring(SIZE(bond,1))     ! rg index for nodes, current ring (0=no, 1=yes)
  INTEGER :: rngflg                 ! 0 = 'no ring', 1 = 'yes ring'
  INTEGER :: ngr                    ! number of nodes
  INTEGER :: ogr                    ! number of -O- nodes
  INTEGER :: ncd                    ! number of Cd node
  INTEGER :: nbnei(SIZE(bond,1)),nei_ind(SIZE(bond,1),4)

  INTEGER :: netrack                ! # of ester tracks 
  INTEGER :: etrack(mxtrk,mxlest)    ! etrack(i,j) ester nodes for the ith track         
  INTEGER :: etracklen(mxtrk)

  CHARACTER(LEN=12),PARAMETER  :: progname='simpolgrp'
  CHARACTER(LEN=70) :: mesg

  ring=0 ; nbnei=0 ; neigh=' ' ; simpgroup(:)=0  
  tgroup(:)=group(:)
  
! count the number of nodes
  ngr=COUNT(group/=' ')
  ogr=COUNT(group(1:ngr)=='-O-')
  ncd=COUNT(group(1:ngr)(1:2)=='Cd')

! if rings exist, then find the nodes that belong to the rings
  IF (nring > 0) CALL findring(1,2,ngr,bond,rngflg,ring)

! find neighbours
  DO i=1,ngr
    DO j=1,ngr
      IF (bond(i,j)/=0) THEN
        nbnei(i)=nbnei(i)+1                ! number of neighbours
        neigh(i,nbnei(i))=group(j)         ! groups of neighbours
        nei_ind(i,nbnei(i))=j
      ENDIF
    ENDDO
  ENDDO

! carbon number (GROUP 1) 
! -----------------------
  simpgroup(1)=ngr-ogr

! Rings (GROUP 3 & 4) - need revision but should be fine if aromatic has only 1 ring
! --------------------
  IF (nring/=0) THEN
    IF (INDEX(chem,'c1')/= 0) THEN ;  simpgroup(3)=nring   !--- GROUP 3
    ELSE                           ;  simpgroup(4)=nring   !--- GROUP 4    
    ENDIF    
  ENDIF

! C=C bonds (GROUP 5 >C=C< and 6 C=C-C=O in a ring)
! --------------------
  IF (ncd/=0) THEN
    simpgroup(5) = ncd/2                                           !--- GROUP 5
    DO i=1,ngr
      IF (group(i)(1:2)=='Cd') THEN
        ico=0 ; icd=0
        DO j=1,nbnei(i)
!          IF (neigh(i,j)(1:2)=='CO')   simpgroup(6)=simpgroup(6)+1 !--- GROUP 6  
!          IF (neigh(i,j)(1:3)=='CHO')  simpgroup(6)=simpgroup(6)+1 !--- GROUP 6  
          IF (neigh(i,j)(1:2)=='CO') ico=j
          IF (neigh(i,j)(1:2)=='Cd') icd=j
        ENDDO
        IF ((ico/=0).AND.(icd/=0)) THEN
          IF ((ring(i)==1).AND.(ring(nei_ind(i,ico))==1).AND.(ring(nei_ind(i,icd))==1)) THEN
            simpgroup(6)=simpgroup(6)+1 !--- GROUP 6  
          ENDIF
        ENDIF
      ENDIF     
    ENDDO
  ENDIF

! Ester (GROUP 11 & 30) - formate considered as ester
! --------------------
  IF (ogr/=0) THEN
  
! make ester tracks, i.e. (CO-O)x tracks
    CALL estertrack(chem,bond,group,ngr,netrack,etracklen,etrack)

    DO i=1,netrack                           ! scroll each identified tracks
      nester=etracklen(i)/2                  ! integer division here 
      ogr = ogr-nester                       ! keep record of remaining -O-

! "simple" CO-O, CO-O-CO-O, etc structures (odd # of nodes in track)
      IF (MOD(etracklen(i),2)==0) THEN   
      
!       check for nitroester >C(NO2)-CO-O- special case
        nitroester=0                         
        IF (tgroup(etrack(i,1))(1:2)=='CO') THEN ; inod=etrack(i,1)
        ELSE                                     ; inod=etrack(i,etracklen(i))
        ENDIF
        DO j=1,nbnei(inod)
          IF (INDEX(neigh(inod,j),'(NO2)')/=0) nitroester=1   
        ENDDO

        IF (nitroester==0) THEN 
          simpgroup(11)=simpgroup(11)+nester             !--- GROUP 11
        ELSE
          simpgroup(30)=simpgroup(30)+1                  !--- GROUP 30
          simpgroup(11)=simpgroup(11)+nester-1           !--- GROUP 11
        ENDIF
        
!       "erase" all groups in the track
        DO j=1,etracklen(i)
          tgroup(etrack(i,j))=' '        
        ENDDO

! "Incomplete" ester structure O-CO-O, CO-O-CO etc (even # of nodes in track)
! Here, one node is left "out" of the ester structure - need to select which
      ELSE                 
        simpgroup(11)=simpgroup(11)+nester             !--- GROUP 11 (ignore possible nitroester)
        ifirst=etrack(i,1)
        
!       -O-CO-O- structure  - carbonate
        IF (tgroup(ifirst)=='-O-') THEN  
          ipero=0
          DO j=1,nbnei(ifirst)
            IF (neigh(ifirst,j)(1:3)=='-O-') ipero=1  ! first is "ROOR"
          ENDDO
          IF (ipero==0) THEN                          ! kill from 1st, leave last 
            DO j=1,etracklen(i)-1
              tgroup(etrack(i,j))=' '     
            ENDDO
          ELSE                                        ! leave first (peroxide), kill others
            DO j=2,etracklen(i)
              tgroup(etrack(i,j))=' '     
            ENDDO
          ENDIF    

!       -CO-O-CO- structure - Anhydre
        ELSE                             
          IF (tgroup(ifirst)(1:3)=='CHO') THEN        ! keep aldehyde, kill others
            DO j=2,etracklen(i)
              tgroup(etrack(i,j))=' '     
            ENDDO
          ELSE
            DO j=1,etracklen(i)-1
              tgroup(etrack(i,j))=' '                 ! keep last, kill others
            ENDDO
          ENDIF
        ENDIF
      ENDIF
    ENDDO
  ENDIF

! Peroxide R-O-O-R (GROUP 26)
! ---------------------------
  IF (ogr/=0) THEN
    pero: DO i=1,ngr
      IF (tgroup(i)(1:3)=='-O-') THEN
        DO j=1,nbnei(i)
          inod=nei_ind(i,j)
          IF (tgroup(inod)(1:3)=='-O-') THEN
            tgroup(i)=' ' ; tgroup(inod)=' ' ; ogr=ogr-2
            simpgroup(26) = simpgroup(26)+1           !--- GROUP 26 
            CYCLE pero
          ENDIF
        ENDDO
      ENDIF
    ENDDO pero
  ENDIF

! Ether (GROUP 12, 13, 14)
! ---------------------------
  IF (ogr/=0) THEN
    ether: DO i=1,ngr
      IF (tgroup(i)(1:3)=='-O-') THEN
        IF (ring(i)==0) THEN
          simpgroup(12) = simpgroup(12)+1        !--- GROUP 12 - simple ether
          tgroup(i)=' '; ogr=ogr-1
        ELSE
          IF ((neigh(i,1)(1:1)=='c').AND.(neigh(i,2)(1:1)=='c')) THEN
            simpgroup(14) = simpgroup(14)+1      !--- GROUP 14 - aromatic ether
            tgroup(i)=' '; ogr=ogr-1
          ELSE
            simpgroup(13) = simpgroup(13)+1      !--- GROUP 13 - cyclic ether
            tgroup(i)=' '; ogr=ogr-1
          ENDIF
        ENDIF 
      ENDIF      
    ENDDO ether
  ENDIF
  
! check that all -O- group were considered
  IF (ogr/=0) THEN
    mesg="expect no more -O- groups remains but did not happen"
    CALL stoperr(progname,mesg,chem)
  ENDIF

! loop over ending - unique nodes
! -------------------------------
  DO i=1,ngr
    IF (tgroup(i)==' ') CYCLE
    IF     (tgroup(i)=='CHO')       THEN ; simpgroup(8)=simpgroup(8)+1   ; tgroup(i)=' '  !--- GROUP 8
    ELSEIF (tgroup(i)=='CO ')       THEN ; simpgroup(9)=simpgroup(9)+1   ; tgroup(i)=' '  !--- GROUP 9
    ELSEIF (tgroup(i)=='CO(OH)')    THEN ; simpgroup(10)=simpgroup(10)+1 ; tgroup(i)=' '  !--- GROUP 10
    ELSEIF (tgroup(i)=='CO(OOH)')   THEN ; simpgroup(28)=simpgroup(28)+1 ; tgroup(i)=' '  !--- GROUP 28
    ELSEIF (tgroup(i)=='CO(OONO2)') THEN ; simpgroup(25)=simpgroup(25)+1 ; tgroup(i)=' '  !--- GROUP 25
    ELSEIF (tgroup(i)=='CO(OOOH)')  THEN ; simpgroup(28)=simpgroup(28)+1 ; tgroup(i)=' '  !--- GROUP 28 - peracid like
    ELSEIF (tgroup(i)=='CO(ONO2)')  THEN ; simpgroup(25)=simpgroup(25)+1 ; tgroup(i)=' '  !--- GROUP 25 - PAN like
    ELSEIF (tgroup(i)=='CO(NO2)')   THEN 
      simpgroup(9)=simpgroup(9)+1 ; simpgroup(16)=simpgroup(16)+1 ; tgroup(i)=' '         !--- ketone (9) + nitro (16) ?
    ENDIF
  ENDDO

! loop over other moieties
! ------------------------
  aronitro=0
  DO i=1,ngr
    IF (tgroup(i)==' ') CYCLE
    
    nn=countstring(tgroup(i),'(OH)')                                       
    IF (nn/=0) THEN                                                        
      IF     (tgroup(i)(1:1)=='C') THEN ; simpgroup(7)=simpgroup(7)+nn      !--- GROUP 7
      ELSEIF (tgroup(i)(1:1)=='c') THEN ; simpgroup(17)=simpgroup(17)+nn    !--- GROUP 17
      ENDIF                                                                
    ENDIF                                                                  
   
    nn=countstring(tgroup(i),'(OOH)')                                      
    IF (nn/=0) simpgroup(27)=simpgroup(27)+nn                               !--- GROUP 27
   
    nn=countstring(tgroup(i),'(OOOH)')                                     
    IF (nn/=0) simpgroup(27)=simpgroup(27)+nn                               !--- GROUP 27, OOH like
   
    nn=countstring(tgroup(i),'(ONO2)')                                     
    IF (nn/=0) simpgroup(15)=simpgroup(15)+nn                               !--- GROUP 15
   
    nn=countstring(tgroup(i),'(NO2)')                                      
    IF (nn/=0) THEN 
      simpgroup(16)=simpgroup(16)+nn                                        !--- GROUP 16  
      IF (tgroup(i)(1:1)=='c') aronitro=1                                   ! raise flag for nitrophenols, GROUP 30
    ENDIF
  ENDDO
  
! check nitrophenol - switch phenol (17) to nitrophenol (29)
  IF ((simpgroup(16)/=0).AND.(simpgroup(17)/=0)) THEN
    simpgroup(29)=simpgroup(17)                                             !--- GROUP 29 & 17
    simpgroup(17)=0
  ENDIF

END SUBROUTINE simpolgrp

END MODULE simpoltool
