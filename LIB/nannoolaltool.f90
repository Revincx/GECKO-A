MODULE nannoolaltool
IMPLICIT NONE
CONTAINS

! SUBROUTINE nannoolalprop(chem,bond,group,nring,rjg,Tb,logPvap,deltaHvap)
! SUBROUTINE load_nandat()
! SUBROUTINE getnangrp(chem,group,bond,nring,rjg,nangroup,nbatom)

! ======================================================================
! Purpose: Compute the vapor pressure, enthalpie of vaporization and 
! the boiling point of the species provided as input based on the 
! nannoolal SARs.
!
! Data for Tb: 
! Nannoolal, Y., Rarey, J., Ramjugernath, D., and Cordes, W.: Estimation
! of pure component properties:  Part 1. Estimation of the normal 
! boiling point of non-electrolyte organic compounds via group contributions
! and group interactions, Fluid Phase Equilibria, 226, 45-63, 2004.
!
! Data for Psat: 
! Nannoolal, Y., Rarey, J., and Ramjugernath, D.: Estimation of pure
! component properties: Part 3. Estimation of the vapor pressure of
! non-electrolyte organic compounds via group contributions and group 
! interactions, Fluid Phase Equilibria, 269, 117-133, 2008.
! ======================================================================
SUBROUTINE nannoolalprop(chem,bond,group,nring,rjg,Tb,logPvap,deltaHvap)
  USE database, ONLY: nanweight  ! weight of the various nannoolal groups
  USE toolbox, ONLY: stoperr
  IMPLICIT NONE
   
  CHARACTER(LEN=*),INTENT(IN) :: chem
  CHARACTER(LEN=*),INTENT(IN) :: group(:) ! group, without ring joining char
  INTEGER,INTENT(IN) :: bond(:,:)
  INTEGER,INTENT(IN) :: nring  
  INTEGER,INTENT(IN) :: rjg(:,:)      ! ring-join group pairs
  REAL,INTENT(OUT)   :: Tb            ! boiling T
  REAL,INTENT(OUT)   :: logPvap       ! vapor pressure (log10)
  REAL,INTENT(OUT)   :: deltaHvap     ! vaporization enthalpie

  REAL,PARAMETER :: temp=298.       ! ref temperature for Psat and latent heat
  INTEGER :: Nangroup(219)          ! Nannoolal groups

  INTEGER :: i
  INTEGER :: sum1
  REAL    :: sum2, xsum
  REAL    :: GI
  REAL    :: Trb,dB
  INTEGER :: nbatom

  CHARACTER(LEN=15),PARAMETER :: progname='nannoolalprop '
  CHARACTER(LEN=256)           :: mesg

! get the Nannoolal group numbers for chem (without group interaction)
  CALL getnangrp(chem,group,bond,nring,rjg,Nangroup,nbatom)    

! ----------
! COMPUTE Tb
! ----------

! In Nannoolal 2004, group 134 (C=C-C=O) was assigned to # 118
  IF (Nangroup(118)/=0) THEN
    mesg="Group 118 in Nannoolal SAR (diketone) is expected to be unused"
    CALL stoperr(progname,mesg,chem)
  ENDIF
  Nangroup(118)=Nangroup(134) ; Nangroup(134)=0

! add group contribution 
  xsum=0.
  DO i=1,134  
    IF (Nangroup(i)/=0.) xsum=xsum+Nangroup(i)*nanweight(i,1)
  ENDDO

! add group contribution for -OOH, C(O)OOH, PAN (added here)
! BA Warning: in Nannoolal 2008, these ID # are used for something else
  DO i=213,219  
    IF (Nangroup(i)/=0.) xsum=xsum+Nangroup(i)*nanweight(i,1)
  ENDDO

! get group interaction for Tb SAR
  CALL gi4tb(Nangroup)

! add contribution from group interactions
  sum1=0  ;  sum2=0.
  DO i=135,212
    IF (Nangroup(i) > 0) THEN
      sum1=sum1+Nangroup(i)                  ! count # of interaction
      sum2=sum2+Nangroup(i)*nanweight(i,1)   ! add contribution
    ENDIF
  ENDDO
  sum1=sum1-1
  IF (sum1<1) THEN ; GI=0. ; ELSE ; GI=(1./nbatom)*(sum2/sum1) ; ENDIF

! add all and compute Tb
  xsum = xsum+GI
  Tb = (xsum/(nbatom**0.6583 + 1.6868))+84.3395

! retore group 134 # 118
  Nangroup(134)=Nangroup(118) ; Nangroup(118)=0

! -----------------------
! COMPUTE PVAP & VAPORIZATION ENTHALPIE
! -----------------------
 
! add group contribution 
  xsum=0.
  DO i=1,134
    IF (Nangroup(i)/=0) xsum = xsum+Nangroup(i)*nanweight(i,2)*0.001
  ENDDO

! add group contribution for -OOH, C(O)OOH, PAN (added here)
! BA Warning: in Nannoolal 2008, these ID # are used for something else
  DO i=213,219
    IF (Nangroup(i)/=0) xsum = xsum+Nangroup(i)*nanweight(i,2)*0.001
  ENDDO

! get group interaction for Tb SAR
  CALL gi4vp(Nangroup)

! add contribution from group interactions
  sum1=0 ;  sum2=0.
  DO i=135,212
    IF (Nangroup(i) > 0) THEN
      sum1=sum1+Nangroup(i)
      sum2=sum2+Nangroup(i)*nanweight(i,2)
    ENDIF
  ENDDO
  sum1 = sum1-1
  sum2=sum2*0.001 ! weight was 1E3
  IF (sum1<1) THEN ; GI=0. ; ELSE ; GI=(1./nbatom)*(sum2/sum1) ; ENDIF
  
! add all and compute Tb
  Trb=temp/Tb
  dB=(xsum + GI) - 0.176055
  logPvap = (4.1012+dB)*((Trb-1.)/(Trb-(1./8.)))  ! log10 value
  deltaHvap = (56*(4.1012+dB)*8.314*Tb*2.303) / ((Tb/temp)-8.)**2.

! check errors
! ------------
  IF (logPvap > 10) THEN
    PRINT*, 'chem=',TRIM(chem)
    DO i=1,219
      IF (Nangroup(i)/=0) PRINT*, 'Nangroup(',i,')=',Nangroup(i)
    ENDDO        
    PRINT*, 'Tb=',Tb
    PRINT*, 'sum2=',sum2
    PRINT*, 'm=',sum1+1
    PRINT*, 'nbatom=',nbatom
    PRINT*, 'GI=',GI
    PRINT*, 'dB=',dB
    mesg="logPvap > 10"
    CALL stoperr(progname,mesg,chem)
  ENDIF      
      
END SUBROUTINE nannoolalprop 

! ======================================================================
! Purpose: load the data for the nannoolal SAR, i.e. the group
! for Tb and Psat. Data are stored in the database module 
! ======================================================================
SUBROUTINE load_nandat()
  USE keyparameter, ONLY: tfu1,dirgecko
  USE database, ONLY: nanweight ! weight of the various nannoolal groups
  IMPLICIT NONE

  INTEGER :: j,ilin,ierr
  CHARACTER(LEN=512) filename
  CHARACTER(LEN=150) line
  
  filename=TRIM(dirgecko)//'DATA/nannoolal.dat'
  OPEN(tfu1,FILE=filename, FORM='FORMATTED', STATUS='OLD', IOSTAT=ierr)
  IF (ierr/=0) THEN
    WRITE(6,*) '--error--, in load_nandat while trying top open file:',TRIM(filename)
    STOP "in load_nandat"
  ENDIF

! read the data
  ilin=0
  rdloop: DO
    READ(tfu1,'(a)') line
    IF (line(1:1)=='!') CYCLE rdloop
    IF (line(1:3)=='END') EXIT rdloop
    
    ilin=ilin+1
    READ(line,*,IOSTAT=ierr) (nanweight(ilin,j),j=1,2)
    IF (ierr/=0) THEN
      WRITE(6,*) '--error--, in load_nandat while reading line:',TRIM(line)
      STOP "in load_nandat"
    ENDIF
  ENDDO rdloop
  CLOSE(tfu1)

  IF (ilin/=SIZE(nanweight,1)) THEN
    WRITE(6,*) '--error--, unexpected # of data in :',TRIM(filename)
    STOP "in load_nandat"
  ENDIF

END SUBROUTINE load_nandat

!=======================================================================
! Purpose: return the number of Nannoonal groups for the species 
! provided as input
!=======================================================================
SUBROUTINE getnangrp(chem,group,bond,nring,rjg,nangroup,nbatom)
  USE keyparameter, ONLY:mxring,mxcp,mxtrk,mxlcd,mxlest
  USE rjtool, ONLY: rjsrm
  USE mapping, ONLY: chemmap, gettrack
  USE atomtool
  USE ringtool, ONLY: findring,ring_data
  USE cdtool, ONLY: cdcase2
  USE toolbox, ONLY: stoperr,countstring
  USE mapping, ONLY: estertrack,alkenetrack
  IMPLICIT NONE

  CHARACTER(LEN=*),INTENT(IN) :: chem
  CHARACTER(LEN=*),INTENT(IN) :: group(:) ! group without ring joining char
  INTEGER,INTENT(IN)  :: bond(:,:)
  INTEGER,INTENT(IN)  :: nring  
  INTEGER,INTENT(IN)  :: rjg(mxring,2)    ! ring-join group pairs
  INTEGER,INTENT(OUT) :: nangroup(219)    ! Nannoonal groups  
  INTEGER,INTENT(OUT) :: nbatom           ! # of atoms in chem (exclude H atom)

! internal
  CHARACTER(LEN=LEN(group(1))) :: tgroup(SIZE(group))
  CHARACTER(LEN=LEN(group(1))) :: neigh(SIZE(group),4)
  INTEGER :: tbond(SIZE(bond,1),SIZE(bond,2))
  INTEGER :: i,j,k
  INTEGER :: begrg,endrg
  INTEGER :: ring(SIZE(group))     ! rg index for nodes, current ring (0=no, 1=yes)
  INTEGER :: rngflg                ! 0 = 'no ring', 1 = 'yes ring'
  INTEGER :: nbnei(SIZE(group)),nei_ind(SIZE(group),4)
  INTEGER :: xxc,xxh,xxn,xxo,xxr,xxs,xxfl,xxbr,xxcl
  INTEGER :: nodflg(SIZE(group))
  
  INTEGER :: ngr,ncd,ogr,nca
  INTEGER :: ifirst,ilast,nod1,nod2,r1,r2
  INTEGER :: iformate,iunsat,iaro,nn,b2ene,b2engr,b2ennod,irg,ene1
  INTEGER :: unsat1,unsat2,sat1,sat2
  INTEGER :: northo,nmeta,npara

  INTEGER :: netrack                   ! # of ester tracks 
  INTEGER :: etrack(mxtrk,mxlest)      ! etrack(i,j) ester nodes for the ith track         
  INTEGER :: etracklen(mxtrk)          ! length of the ith track

  INTEGER :: ncdtrack                  ! # of Cd tracks 
  INTEGER :: cdtrack(mxtrk,mxlcd)      ! (i,j) Cd nodes "j" for the ith track         
  INTEGER :: cdtracklen(mxtrk)         ! length of the ith track         

  INTEGER :: rgrec(mxcp,SIZE(group))     ! record of all rings (include duplicates)
  INTEGER :: nrec                        ! # of records in rgrec
  INTEGER :: lrec(mxcp)                  ! length of the record
  INTEGER :: fgrec(mxcp)                 ! record flag (duplicate or not)
  INTEGER,PARAMETER :: mxirg=6           ! max # of distinct rings 
  INTEGER :: ndrg                        ! # of distinct rings
  INTEGER :: trackrg(mxirg,SIZE(group))  ! (a,:)== track (node #) belonging ring a
  LOGICAL :: lorgnod(mxirg,SIZE(group))  ! (a,b)==true if node b belong to ring a
  LOGICAL :: allfound
  INTEGER :: lenring  

  INTEGER :: track(mxcp,SIZE(group))
  INTEGER :: trlen(mxcp)
  INTEGER :: ntr

  CHARACTER(LEN=11),PARAMETER :: progname='getnangrp '
  CHARACTER(LEN=70)           :: mesg

! initialize
  ring(:)=0 ;  nbnei(:)=0 ; neigh(:,:)=' ' ; nei_ind(:,:)=0 
  tgroup(:)=group(:) ; tbond(:,:)=bond(:,:)
  nangroup(:)=0  ; nodflg(:)=0

! count the number of nodes
  ngr=COUNT(group/=' ')
  ogr=COUNT(group(1:ngr)=='-O-')
  ncd=COUNT(group(1:ngr)(1:2)=='Cd')
  nca=ngr-ogr

! if rings exist, then find the nodes that belong to the rings
  IF (nring > 0) THEN
    DO i=1,nring                ! find the nodes that belong to a ring
      begrg=rjg(i,1)  ;  endrg=rjg(i,2)
      CALL findring(begrg,endrg,ngr,bond,rngflg,ring)
    ENDDO
  ENDIF

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
  
! ----------------------------------------------
! ESTER, ANHYDRE, CARBONATE
! ----------------------------------------------

  IF (ogr/=0) THEN
  
! make ester tracks, i.e. (CO-O)x tracks
    CALL estertrack(chem,bond,group,ngr,netrack,etracklen,etrack)

    DO i=1,netrack                           ! scroll each identified tracks
    
      SELECT CASE (etracklen(i))

 ! simple ester
      CASE(2)  
        iformate=0 
        ifirst=etrack(i,1)  ;  ilast=etrack(i,2)
        IF (tgroup(ifirst)(1:3)=='CHO') iformate=1        
        IF (tgroup(ilast)(1:3)=='CHO')  iformate=1        
        
        IF (iformate==1) THEN         ;  nangroup(46)=nangroup(46)+1  ! formic acid ester
        ELSEIF (ring(ifirst)/=0) THEN ;  nangroup(47)=nangroup(47)+1  ! lactone
        ELSE                          ;  nangroup(45)=nangroup(45)+1  ! ester  
        ENDIF 

        ! "erase" all groups in the track and update ogr
        ogr = ogr-1
        DO j=1,2 ; tgroup(etrack(i,j))=' ' ;  ENDDO

! anhydre or carbonate
      CASE(3)
        ifirst=etrack(i,1)
        
        ! -O-CO-O- structure  - carbonate
        IF (tgroup(ifirst)=='-O-') THEN  
          IF (ring(ifirst)/=0) THEN ;  nangroup(103) = nangroup(103)+1 ! cyclic carbonate
          ELSE                      ;  nangroup(79)  = nangroup(79)+1  ! non-cyclic carbonate
          ENDIF
          ogr=ogr-2
          
        ! -CO-O-CO- structure - Anhydre
        ELSE                             
          IF (ring(ifirst)/=0) THEN

            iunsat=0
            DO j=1,nbnei(ifirst)
              IF (INDEX(neigh(ifirst,j)(1:2),'Cd')/=0) THEN 
                iunsat=iunsat+1 ; nod1=nei_ind(ifirst,j)
              ELSEIF (INDEX(neigh(ifirst,j)(1:1),'c')/=0) THEN
                iunsat=iunsat+1 ; nod1=nei_ind(ifirst,j)
              ENDIF
            ENDDO   
            DO j=1,nbnei(etrack(i,3))
              IF (INDEX(neigh(etrack(i,3),j)(1:2),'Cd')/=0) THEN
                iunsat=iunsat+1 ; nod2=nei_ind(etrack(i,3),j)  
              ELSEIF (INDEX(neigh(etrack(i,3),j)(1:1),'c')/=0) THEN
                iunsat=iunsat+1 ; nod2=nei_ind(etrack(i,3),j)  
              ENDIF
            ENDDO

            IF (iunsat==2) THEN
              IF (tbond(nod1,nod2)/=0) THEN ; nangroup(96)=nangroup(96)+1   ! cyclic anhydride via Cd=Cd or aro
              ELSE                          ; nangroup(76)=nangroup(76)+1   ! anhydride
              ENDIF
            ELSE                            ; nangroup(76)=nangroup(76)+1   ! anhydride
            ENDIF

          ELSE
            nangroup(76)=nangroup(76)+1   ! anhydride
          ENDIF
          ogr = ogr-1
        ENDIF

        ! "erase" all groups in the track and update ogr
        DO j=1,3 ; tgroup(etrack(i,j))=' ' ; ENDDO

! CO-O-CO-O- : treated as anhydre+ether according to priority rules     
      CASE(4)
        ifirst=etrack(i,1)
        nangroup(76)=nangroup(76)+1   ! anhydride

        ! "erase" groups in the track but leave one -O-
        ogr=ogr-1
        IF (tgroup(ifirst)=='-O-') THEN  
          DO j=2,4  ;  tgroup(etrack(i,j))=' ' ;   ENDDO  ! erase grp 2-4
        ELSE
          DO j=1,3  ;  tgroup(etrack(i,j))=' ' ;   ENDDO  ! erase grp 1-3
        ENDIF
      
! O-CO-O-CO-O- or CO-O-CO-O-CO: priority to anhydre according to priority rules     
      CASE(5)
        ifirst=etrack(i,1)

        IF (tgroup(ifirst)=='-O-') THEN  
          nangroup(76)=nangroup(76)+1   ! anhydride
          ! "erase" groups in the track but leave two -O-
          ogr=ogr-1
          DO j=2,4  ;  tgroup(etrack(i,j))=' ' ;   ENDDO

        ELSE
          nangroup(76)=nangroup(76)+1   ! anhydride
          nangroup(45)=nangroup(45)+1   ! ester  
          ogr=ogr-2
          DO j=1,5  ;  tgroup(etrack(i,j))=' ' ;   ENDDO 
        ENDIF
         
! CO-O-CO-O-CO-O: anhydre+carbonate     
      CASE(6)
        nangroup(79)=nangroup(79)+1     ! non-cyclic carbonate
        nangroup(76)=nangroup(76)+1     ! anhydride
        ogr=ogr-3
        DO j=1,6  ;  tgroup(etrack(i,j))=' ' ;   ENDDO      
      
! ester track not found      
      CASE DEFAULT
        mesg="unidentified length track in ester track"
        CALL stoperr(progname,mesg,chem) 

      END SELECT

    ENDDO    ! loop ester track
  ENDIF      ! ester, anhydre and carbonate

! ----------------------------------------------
! PEROXYDE
! ----------------------------------------------
  IF (ogr/=0) THEN
    pero: DO i=1,ngr
      IF (tgroup(i)(1:3)=='-O-') THEN
        DO j=1,nbnei(i)
          nod1=nei_ind(i,j)
          IF (tgroup(nod1)(1:3)=='-O-') THEN
            tgroup(i)=' ' ; tgroup(nod1)=' ' ; ogr=ogr-2
            nangroup(94)=nangroup(94) + 1                   ! peroxide
            CYCLE pero
          ENDIF
        ENDDO
      ENDIF
    ENDDO pero
  ENDIF

! ----------------------------------------------
! ETHER - EPOXIDE 
! ----------------------------------------------
  IF (ogr/=0) THEN
    ether: DO i=1,ngr
      IF (tgroup(i)(1:3)=='-O-') THEN
        nod1=nei_ind(i,1)  ;  nod2=nei_ind(i,2)
        IF (bond(nod1,nod2)/=0) THEN
            nangroup(39)=nangroup(39) + 1        ! epoxide
            tgroup(i)=' ' ; ogr=ogr-1 
            nodflg(nod1)=1 ; nodflg(nod2)=1
        ELSEIF (ring(i)/=0) THEN 
          iaro=0
          ! check for aromatic cycle Cd(ring)-O-Cd(ring) assumed as aromatic behavior
          IF ((neigh(i,1)(1:1)=='c').AND.(neigh(i,2)(1:1)=='c'))   iaro=1
          IF ((neigh(i,1)(1:2)=='Cd').AND.(neigh(i,2)(1:2)=='Cd')) iaro=1 
          IF (iaro==1) THEN
            nangroup(65)=nangroup(65)+1          ! furan
            tgroup(i)=' ' ; ogr=ogr-1       
          ELSE
            nangroup(38)=nangroup(38)+1          ! ether
            tgroup(i)=' ' ; ogr=ogr-1       
          ENDIF
        ELSE
          nangroup(38)=nangroup(38)+1          ! ether
          tgroup(i)=' ' ; ogr=ogr-1       
        ENDIF 
      ENDIF      
    ENDDO ether
  ENDIF
  
! check that all -O- group were considered
  IF (ogr/=0) THEN
    mesg="expect no more -O- groups remains but did not happen"
    CALL stoperr(progname,mesg,chem)
  ENDIF

! -------------------------------
! CARBONYLS: -CO-, -CHO, -CO(OH), -CO(OONO2), etc ...
! -------------------------------
  DO i=1,ngr

! aldehyde
    IF (tgroup(i)=='CHO ') THEN 
      IF (neigh(i,1)(1:1)=='c' ) THEN ; nangroup(90)=nangroup(90)+1   ! aromatic
      ELSE                            ; nangroup(52)=nangroup(52)+1   ! aliphatic
      ENDIF
      tgroup(i)=' '

! ketone 
! BA july 2020: warning: diketone with priority 1 but no estimate provided  
! in Nannoolal 2008 - here the group is ignored and diketone are treated  
! as 2 ketones, waiting for additional information)
    ELSEIF (tgroup(i)=='CO ') THEN 
      nod1=nei_ind(i,1)  ;  nod2=nei_ind(i,2)
!      IF (tgroup(nod1)=='CO ' ) THEN       
!        nangroup(118)=nangroup(118)+1          ! diketone (kill 2 groups)
!        tgroup(i)=' ' ; tgroup(nod1)=' '
!      ELSEIF (tgroup(nod2)=='CO ' ) THEN
!        nangroup(118)=nangroup(118)+1          ! diketone (kill 2 groups)
!        tgroup(i)=' ' ; tgroup(nod2)=' '
!      ELSEIF (tgroup(nod1)(1:1)=='c' ) THEN
      IF (tgroup(nod1)(1:1)=='c' ) THEN
        nangroup(92)=nangroup(92)+1            ! ketone aromatic
        tgroup(i)=' '
      ELSEIF (tgroup(nod2)(1:1)=='c' ) THEN
        nangroup(92)=nangroup(92)+1            ! aromatic
        tgroup(i)=' '
      ELSE
        nangroup(51)=nangroup(51)+1            ! aliphatic
        tgroup(i)=' '
      ENDIF
      
! other carbonyls (ending groups)
    ELSEIF (tgroup(i)(1:6)=='CO(OH)') THEN 
      nangroup(44)=nangroup(44)+1   ; tgroup(i)=' '  

    ELSEIF (tgroup(i)=='CO(OOH)')   THEN 
      nangroup(219)=nangroup(219)+1 ; tgroup(i)=' '               ! compernolle
													   
    ELSEIF (tgroup(i)=='CO(OOOH)')  THEN               
      nangroup(219)=nangroup(219)+1 ; tgroup(i)=' '               ! assume as CO(OOH)
													   
    ELSEIF (tgroup(i)=='CO(OONO2)') THEN               
      nangroup(213)=nangroup(213)+1 ; tgroup(i)=' '               ! compernolle
       tgroup(i)=' '              

! BA, July 2020. Following the idea that PAN could be represented by
! >C< (7) + -O- (48) + ONO2 (72). Same hypothesis here for "CO(ONO2) and "CO(NO2)".
    ELSEIF (tgroup(i)=='CO(NO2)')   THEN           
      nangroup(7)=nangroup(7)+1                   
      nangroup(68)=nangroup(68)+1
      tgroup(i)=' '

    ELSEIF (tgroup(i)=='CO(ONO2)')  THEN 
      nangroup(7)=nangroup(7)+1 
      nangroup(72)=nangroup(72)+1
       tgroup(i)=' '               
    ENDIF
  ENDDO

! -------------------------------
! ALCOHOL (OH)
! -------------------------------
  DO i=1,ngr
    IF (tgroup(i)==' ') CYCLE
    nn=countstring(tgroup(i),'(OH)')                                       
    IF (nn/=0) THEN
      IF     (tgroup(i)(1:1)=='c')   THEN ; nangroup(37)=nangroup(37)+1  ! phenol
      ELSEIF (nca<5)                 THEN ; nangroup(36)=nangroup(36)+nn ! short chain
      ELSEIF (tgroup(i)(1:3)=='CH2') THEN ; nangroup(35)=nangroup(35)+nn ! primary OH
      ELSEIF (tgroup(i)=='CdH(OH)')  THEN ; nangroup(34)=nangroup(34)+nn ! secondary OH 
      ELSEIF (tgroup(i)(1:2)=='CH')  THEN ; nangroup(34)=nangroup(34)+nn ! secondary OH 
      ELSE                                ; nangroup(33)=nangroup(33)+nn ! tertiary
      ENDIF
    ENDIF
  ENDDO

! -------------------------------
! HYDROPEROXYDE (OOH) & (OOOH)
! -------------------------------
  DO i=1,ngr
    IF (tgroup(i)==' ') CYCLE
    nn=countstring(tgroup(i),'(OOH)')                                       
    IF (nn/=0) nangroup(214)=nangroup(214)+nn                ! Compernolle

    nn=countstring(tgroup(i),'(OOOH)')                                       
    IF (nn/=0) nangroup(214)=nangroup(214)+nn
  ENDDO

! -------------------------------
! NITROGEN MOIETIES (ONO2) & (NO2) & (ONO)
! -------------------------------
  DO i=1,ngr
    IF (tgroup(i)==' ') CYCLE

    nn=countstring(tgroup(i),'(ONO2)')         ! nitrate                              
    IF (nn/=0) nangroup(72)=nangroup(72)+nn

    nn=countstring(tgroup(i),'(NO2)')          ! nitro
    IF (nn/=0) THEN
      IF (tgroup(i)(1:1)=='c') THEN ; nangroup(69)=nangroup(69)+nn
      ELSE                          ; nangroup(68)=nangroup(68)+nn
      ENDIF
    ENDIF

    nn=countstring(tgroup(i),'(ONO)')          ! nitrite
    IF (nn/=0) nangroup(74)=nangroup(74)+nn

  ENDDO

! -----------------
! SATURATED CARBON SKELETON 
! -----------------
  skeloop: DO i=1,ngr
    IF (tgroup(i)=='-O-')     CYCLE skeloop
    IF (tgroup(i)(1:2)=='CO') CYCLE skeloop
    IF (tgroup(i)=='CHO')     CYCLE skeloop
    IF (nodflg(i)==1)         CYCLE skeloop

! raise flag if bonded to electronegative element (O,N,F,Cl) or aromatic chain   
    b2ene=0 ; b2engr=0 ;  b2ennod=0 
    IF (INDEX(tgroup(i),'O')/=0)  b2engr=1   ! bonded to EN in group (OH,NO2,ONO2 ...)
    IF (INDEX(tgroup(i),'F')/=0)  b2engr=1 
    IF (INDEX(tgroup(i),'Cl')/=0) b2engr=1 
    IF (neigh(i,1)=='-O-')        b2ennod=1  ! bonded to EN out group
    IF (neigh(i,2)=='-O-')        b2ennod=1
    b2ene=b2engr+b2ennod                     ! bonded to EN element

    iaro=0
    DO j=1,nbnei(i)
      IF (neigh(i,j)(1:1)=='c') iaro=1 
    ENDDO

! -CH3
    IF (tgroup(i)=='CH3') THEN
      IF     (b2ene/=0) THEN ; nangroup(2)=nangroup(2)+1     ! priority 103
      ELSEIF (iaro/=0)  THEN ; nangroup(3)=nangroup(3)+1     ! priority 104
      ELSE                   ; nangroup(1)=nangroup(1)+1     ! priority 105
      ENDIF

! -CH2-
    ELSEIF (tgroup(i)(1:3)=='CH2') THEN

      !--- inside ring
      IF (ring(i)/=0) THEN
        IF (iaro/=0) THEN     ; nangroup(14)=nangroup(14)+1  ! priority 107

        ELSEIF (b2ene/=0) THEN 
          IF (b2engr/=0) THEN ; nangroup(12)=nangroup(12)+1  ! Priority 110
          ELSE 
            irg=0 
            DO j=1,nbnei(i)
              IF ((neigh(i,j)=='-O-').AND.(ring(nei_ind(i,j))/=0))  irg=1 
            ENDDO
            IF (irg==1) THEN ; nangroup(13)=nangroup(13)+1   ! Priority 111
            ELSE             ; nangroup(9) =nangroup(9)+1    ! Priority 113
            ENDIF
          ENDIF

        ELSE                 ; nangroup(9)=nangroup(9)+1     ! Priority 113
        ENDIF

      !--- outside ring
      ELSE
        IF     (b2ene/=0) THEN ; nangroup(7)=nangroup(7)+1   ! Priority 108
        ELSEIF (iaro/=0)  THEN ; nangroup(8)=nangroup(8)+1   ! priority 109
        ELSE                   ; nangroup(4)=nangroup(4)+1   ! Priority 112
        ENDIF
      ENDIF 

! >CH-
    ELSEIF (tgroup(i)(1:2)=='CH') THEN

      !--- inside ring
      IF (ring(i)/=0) THEN
        IF (iaro/=0) THEN     ; nangroup(14)=nangroup(14)+1   ! priority 107

        ELSEIF (b2ene/=0) THEN 
          IF (b2engr/=0) THEN ; nangroup(12)=nangroup(12)+1   ! Priority 110
          ELSE 
            irg=0 
            DO j=1,nbnei(i)
              IF ((neigh(i,j)=='-O-').AND.(ring(nei_ind(i,j))/=0))  irg=1 
            ENDDO
            IF (irg==1) THEN ; nangroup(13)=nangroup(13)+1  ! Priority 111
            ELSE             ; nangroup(10)=nangroup(10)+1  ! Priority 118
            ENDIF 
          ENDIF

        ELSE                   ; nangroup(10)=nangroup(10)+1  ! Priority 118
        ENDIF

      !--- outside ring
      ELSE
        IF     (b2ene/=0)   THEN ; nangroup(7)=nangroup(7)+1  ! Priority 108
        ELSEIF (iaro/=0)    THEN ; nangroup(8)=nangroup(8)+1  ! priority 109
        ELSE                     ; nangroup(5)=nangroup(5)+1  ! Priority 119
        ENDIF
      ENDIF 

! >C<      
    ELSE IF ((tgroup(i)=='C').OR.(tgroup(i)(1:2)=='C(')) THEN

      !--- inside ring
      IF (ring(i)/=0) THEN
        IF (iaro/=0) THEN     ; nangroup(14)=nangroup(14)+1   ! priority 107

        ELSEIF (b2ene/=0) THEN 
          IF (b2engr/=0) THEN ; nangroup(12)=nangroup(12)+1   ! Priority 110
          ELSE 
            irg=0 
            DO j=1,nbnei(i)
              IF ((neigh(i,j)=='-O-').AND.(ring(nei_ind(i,j))/=0))  irg=1 
            ENDDO
            IF (irg==1) THEN ; nangroup(13)=nangroup(13)+1    ! Priority 111
            ELSE             ; nangroup(11)=nangroup(11)+1    ! Priority 120
            ENDIF
          ENDIF

        ELSE                 ; nangroup(11)=nangroup(11)+1    ! Priority 120
        ENDIF

      !--- outside ring
      ELSE
        IF     (b2ene/=0) THEN ; nangroup(7)=nangroup(7)+1    ! Priority 108
        ELSEIF (iaro/=0)  THEN ; nangroup(8)=nangroup(8)+1    ! priority 109
        ELSE                   ; nangroup(6)=nangroup(6)+1    ! Priority 121
        ENDIF
      ENDIF
    ENDIF 
  ENDDO skeloop

! -----------------
! AROMATIC CARBON SKELETON 
! -----------------
  IF (INDEX(chem,'c')/=0) THEN
    aroloop: DO i=1,ngr
      IF (tgroup(i)(1:2)=='cH') THEN 
        nangroup(15)=nangroup(15)+1                             ! priority 106
    
      ELSE IF (tgroup(i)(1:1)=='c') THEN
        IF (b2ene/=0) THEN
          nangroup(17)=nangroup(17)+1                           ! Priority 114
        ELSEIF ((nbnei(i)==3).AND. (neigh(i,1)(1:1)=='c') .AND. &
                (neigh(i,2)(1:1)=='c') .AND. (neigh(i,3)(1:1)=='c')) THEN 
          nangroup(18)=nangroup(18)+1                           ! Priority 115
        ELSE  
          nangroup(16)=nangroup(16)+1                           ! Priority 117
        ENDIF

        ! add correction factor for "-C=C-CO"
        IF (neigh(i,1)(1:2)=='CO') nangroup(134)=nangroup(134)+1
        IF (neigh(i,1)=='CHO')     nangroup(134)=nangroup(134)+1
      ENDIF
    
    ENDDO aroloop
  ENDIF

! -----------------
! >C=C< CARBON SKELETON  
! -----------------
  IF (ncd/=0) THEN
   CALL alkenetrack(chem,tbond,tgroup,ngr,ncdtrack,cdtracklen,cdtrack)
   DO i=1,ncdtrack

     ! conjugated alkene
     IF (cdtracklen(i)==4) THEN
       irg=0
       DO j=1,4 ; IF (ring(cdtrack(i,j))/=0) irg=irg+1 ;  ENDDO
       IF (irg==4) THEN ; nangroup(88)=nangroup(88)+1         ! Priority 6
       ELSE             ; nangroup(89)=nangroup(89)+1         ! Priority 7
       ENDIF

     ! "simple" alkene     
     ELSE
       ene1=0 ; b2ene=0 ; iaro=0 ; irg=0
       nod1=cdtrack(i,1) ; nod2=cdtrack(i,2)
       IF (tgroup(nod1)=='CdH2') ene1=ene1+1
       IF (tgroup(nod2)=='CdH2') ene1=ene1+1
       IF (INDEX(tgroup(nod1),'O')/=0) b2ene=1 ! enol like
       IF (INDEX(tgroup(nod2),'O')/=0) b2ene=1 ! enol like
       DO j=1,nbnei(nod1)
         IF (neigh(nod1,j)=='-O-') b2ene=1
         IF (neigh(nod1,j)(1:1)=='c') iaro=1
       ENDDO
       DO j=1,nbnei(nod2)
         IF (neigh(nod2,j)=='-O-') b2ene=1
         IF (neigh(nod2,j)(1:1)=='c') iaro=1
       ENDDO
       IF ((ring(nod1)/=0).AND.(ring(nod2)/=0)) irg=1
       
       IF     (ene1/=0)  THEN ; nangroup(61)=nangroup(61)+1    ! Priority 57 
       ELSEIF (irg/=0)   THEN ; nangroup(62)=nangroup(62)+1    ! Priority 60 
       ELSEIF (b2ene/=0) THEN ; nangroup(60)=nangroup(60)+1    ! Priority 58 
       ELSEIF (iaro/=0)  THEN ; nangroup(59)=nangroup(59)+1    ! Priority 59 
       ELSE                   ; nangroup(58)=nangroup(58)+1    ! Priority 62
       ENDIF         
     ENDIF
     
     ! add correction factor for "-C=C-CO"
     DO j=1,cdtracklen(i)
       nod1=cdtrack(i,j)
       DO k=1,nbnei(nod1)
         IF (neigh(nod1,k)(1:2)=='CO') nangroup(134)=nangroup(134)+1
         IF (neigh(nod1,k)=='CHO')     nangroup(134)=nangroup(134)+1
       ENDDO
     ENDDO
   ENDDO
 ENDIF

! ------------------------
! RING CORRECTION FACTORS 
! -----------------------
  tgroup(:)=group(:)    ! restore groups 

! find all rings shorter than 6 nodes
  rgrec(:,:)=0 ; nrec=0 ; lrec(:)=0 ; fgrec(:)=0
  DO i=1,ngr
    IF (ring(i)==0) CYCLE
    CALL ring_data(i,ngr,tbond,tgroup,ndrg,lorgnod,trackrg)
    DO j=1,ndrg
      lenring=COUNT(lorgnod(j,:))
      IF (lenring>5) CYCLE
      nrec=nrec+1
      lrec(nrec)=lenring
      rgrec(nrec,:)=trackrg(j,:)
    ENDDO
  ENDDO

! flag down duplicate rings
  fgrec(1:nrec)=1
  DO i=1,nrec-1
    IF (fgrec(i)==0) CYCLE
    DO j=i+1,nrec
      IF (fgrec(j)==0) CYCLE
      IF (lrec(i)/=lrec(j)) CYCLE
      allfound=.TRUE.
      r1loop: DO r1=1,lrec(i)
        DO r2=1,lrec(i)
          IF (rgrec(i,r1)==rgrec(j,r2)) CYCLE r1loop
        ENDDO
        allfound=.FALSE.
      ENDDO r1loop
      IF (allfound) fgrec(j)=0
    ENDDO
  ENDDO

! set group
  DO i=1,nrec
    IF (fgrec(i)==1) THEN
      IF (lrec(i)==3) nangroup(125) = nangroup(125)+1
      IF (lrec(i)==4) nangroup(125) = nangroup(125)+1
      IF (lrec(i)==5) nangroup(126) = nangroup(126)+1
    ENDIF
  ENDDO

! -----------------------
! STERIC EFFECTS
! -----------------------
  n1loop: DO i=1,ngr-1
    IF (nbnei(i)<3) CYCLE n1loop      ! not enough neighbors on nod1
    n2loop: DO j=i,ngr
      IF (nbnei(j)<3) CYCLE n2loop    ! not enough neighbors on nod2
      IF (bond(i,j)==1) THEN
        
        unsat1=0 ; unsat2=0 ; sat1=0 ; sat2=0
        DO k=1,nbnei(i)
          IF (nei_ind(i,k)==j)   CYCLE
          IF (neigh(i,k)=='-O-') CYCLE
          IF ( (neigh(i,k)(1:2)=='Cd').OR.(neigh(i,k)(1:1)=='c') ) THEN
            unsat1=unsat1+1
          ELSE
            sat1=sat1+1
          ENDIF
        ENDDO
        
        DO k=1,nbnei(j)
          IF (nei_ind(j,k)==i)   CYCLE
          IF (neigh(j,k)=='-O-') CYCLE
          IF ( (neigh(j,k)(1:2)=='Cd').OR.(neigh(j,k)(1:1)=='c') ) THEN
            unsat2=unsat2+1
          ELSE
            sat2=sat2+1
          ENDIF
        ENDDO

        IF     (sat1+sat2==6) THEN ; nangroup(133) = nangroup(133)+1
        ELSEIF (sat1+sat2==5) THEN ; nangroup(132) = nangroup(132)+1
        ELSEIF (sat1==3) THEN
          IF (sat2+unsat2==2) nangroup(130) = nangroup(130)+1
        ELSEIF (sat2==3) THEN
          IF (sat1+unsat1==2) nangroup(130) = nangroup(130)+1
        ELSEIF ( (sat1==2).AND.(sat2==2)   ) THEN 
          nangroup(131) = nangroup(131)+1
        ENDIF
      ENDIF
    ENDDO n2loop
  ENDDO n1loop

! -----------------------
! ORTHO, META, PARA
! -----------------------
  IF ((nring==1).AND.(INDEX(chem,'c')/=0)) THEN
    northo=0 ; nmeta=0 ; npara=0
    omploop: DO i=1,ngr
      IF (group(i)/='c') CYCLE omploop
      CALL gettrack(bond,i,ngr,ntr,track,trlen)
      DO j=1,ntr
        IF (trlen(j)<4) CYCLE
        IF (group(track(j,2))=='c') northo=1 ! raise flag
        IF (group(track(j,3))=='c') nmeta=1  ! raise flag
        IF (group(track(j,4))=='c') npara=1  ! raise flag
      ENDDO
      IF ((northo+nmeta+npara)==1) THEN
        IF     (northo==1) THEN ; nangroup(127) = nangroup(127)+1
        ELSEIF (nmeta==1)  THEN ; nangroup(128) = nangroup(128)+1
        ELSEIF (npara==1)  THEN ; nangroup(129) = nangroup(129)+1
        ENDIF
      ENDIF
      EXIT omploop
    ENDDO omploop
  ENDIF
 
! -------------------
! CORRECTIONS - miscellaneous
! ------------------
  CALL getatoms(chem,xxc,xxh,xxn,xxo,xxr,xxs,xxfl,xxbr,xxcl)
  nbatom=xxc+xxn+xxo+xxs+xxfl+xxbr+xxcl
 
  
  IF (xxh==0) nangroup(123) = 1    ! no H atom correction
  IF (xxh==1) nangroup(124) = 1    ! only 1 H atom correction
     
END SUBROUTINE getnangrp 

!=======================================================================
! PURPOSE: Compute group interaction for the Nannoolal SAR for Tb.
! Groupe definition for the Tb and Vapor Pressure are not strictly 
! identical - need to be computed specifically.
!=======================================================================
SUBROUTINE gi4tb(Nangroup)
  IMPLICIT NONE
  INTEGER,INTENT(INOUT) :: nangroup(:)
  
  INTEGER :: sumgrp(19)

! clean any previous GI data
  Nangroup(135:212)=0
  
! Group definition for Tb SAR (same as VP, but (37 only in B)
! ----------------------------------------
  sumgrp(1)=nangroup(34)+nangroup(35)+nangroup(36)+ & ! group A (OH) - (37) not in A for VP 
              nangroup(37)+nangroup(214)                ! OOH added to A 
  sumgrp(2)=nangroup(37)                              ! group B (OH aromatic)
  sumgrp(3)=nangroup(44)+nangroup(219)                ! group C (CO(OH) with CO(OOH) added)
  sumgrp(4)=nangroup(38)                              ! group D (Ether)
  sumgrp(5)=nangroup(39)                              ! group E (Epox)
  sumgrp(6)=nangroup(45)+nangroup(46)+nangroup(47)    ! group F (Ester)
  sumgrp(7)=nangroup(51)+nangroup(92)                 ! group G (ketone)
  sumgrp(8)=nangroup(52)+nangroup(90)                 ! group H (aldehyde)
  sumgrp(9)=nangroup(65)                              ! group I ether (aromatic)
  sumgrp(10)=nangroup(54)                             ! group J sulfide (not used)
  sumgrp(11)=nangroup(56)                             ! group K aromatic sufur (not used)
  sumgrp(12)=nangroup(53)                             ! group L (SH thiol) (not used)
  sumgrp(13)=nangroup(40)+nangroup(41)                ! group M primary amine (not used)
  sumgrp(14)=nangroup(42)+nangroup(97)                ! group N secondary amine (not used)
  sumgrp(15)=nangroup(80)                             ! group O isocyanate -OCN (not used)
  sumgrp(16)=nangroup(57)                             ! group P cyanide -CN (not used)
  sumgrp(17)=nangroup(69)                             ! group Q nitrate 
  sumgrp(18)=nangroup(66)                             ! group R aromatic N cyclic 5 (not used)
  sumgrp(19)=nangroup(67)                             ! group S aromatic N cyclic 6 (not used)

! ----------------------------
! GROUP INTERACTION FOR TB SAR
! ----------------------------
  IF (sumgrp(1)/=0) THEN    ! OH interactions
    IF (sumgrp(1) >1)  nangroup(135)=sumgrp(1)*(sumgrp(1)-1)  ! 135 = A-A 
    IF (sumgrp(3)/=0)  nangroup(139)=sumgrp(1)*sumgrp(3)*2    ! 139 = A-C 
    IF (sumgrp(4)/=0)  nangroup(140)=sumgrp(1)*sumgrp(4)*2    ! 140 = A-D
    IF (sumgrp(5)/=0)  nangroup(141)=sumgrp(1)*sumgrp(5)*2    ! 141 = A-E
    IF (sumgrp(6)/=0)  nangroup(142)=sumgrp(1)*sumgrp(6)*2    ! 142 = A-f
    IF (sumgrp(7)/=0)  nangroup(143)=sumgrp(1)*sumgrp(7)*2    ! 143 = A-G
    IF (sumgrp(9)/=0)  nangroup(146)=sumgrp(1)*sumgrp(9)*2    ! 146 = A-I
  ENDIF

  IF (sumgrp(2)/=0) THEN    ! OH AROMATIC interactions (NOTE 37 BELONG TO A AND B ?)
    IF (sumgrp(2) >1)  nangroup(148)=sumgrp(2)*(sumgrp(2)-1)  ! 148 = BB 
    IF (sumgrp(3)/=0)  nangroup(150)=sumgrp(2)*sumgrp(3)*2    ! 150 = B-C 
    IF (sumgrp(4)/=0)  nangroup(151)=sumgrp(2)*sumgrp(4)*2    ! 151 = A-D
    IF (sumgrp(6)/=0)  nangroup(152)=sumgrp(2)*sumgrp(6)*2    ! 152 = A-F
!    IF (sumgrp(7)/=0)  nangroup(153)=sumgrp(2)*sumgrp(7)*2    ! 153 = A-G
    IF (sumgrp(17)/=0) nangroup(155)=sumgrp(2)*sumgrp(17)*2   ! 155 = BQ
  ENDIF

  IF (sumgrp(3)/=0) THEN  ! CO(OH) interactions
    IF (sumgrp(3) > 1) nangroup(172)=sumgrp(3)*(sumgrp(3)-1)  ! 172 = C-C
    IF (sumgrp(4)/=0)  nangroup(173)=sumgrp(3)*sumgrp(4)*2    ! 173 = C-D
    IF (sumgrp(6)/=0)  nangroup(174)=sumgrp(3)*sumgrp(6)*2    ! 174 = C-F
    IF (sumgrp(7)/=0)  nangroup(175)=sumgrp(3)*sumgrp(7)*2    ! 175 = C-G
  ENDIF

  IF (sumgrp(4)/=0) THEN  ! -O- interactions
    IF (sumgrp(4) > 1) nangroup(178)=sumgrp(4)*(sumgrp(4)-1)  ! 178 = D-D
    IF (sumgrp(5)/=0)  nangroup(179)=sumgrp(4)*sumgrp(5)*2    ! 179 = D-E
    IF (sumgrp(6)/=0)  nangroup(180)=sumgrp(4)*sumgrp(6)*2    ! 180 = D-F
    IF (sumgrp(7)/=0)  nangroup(181)=sumgrp(4)*sumgrp(7)*2    ! 181 = D-G
    IF (sumgrp(8)/=0)  nangroup(182)=sumgrp(4)*sumgrp(8)*2    ! 182 = D-H
    IF (sumgrp(17)/=0) nangroup(184)=sumgrp(4)*sumgrp(17)*2   ! 184 = D-Q
    IF (sumgrp(9)/=0)  nangroup(186)=sumgrp(4)*sumgrp(9)*2    ! 186 = D-I
  ENDIF

  IF (sumgrp(5)/=0) THEN  ! -O- interactions
    IF (sumgrp(5) > 1) nangroup(187)=sumgrp(5)*(sumgrp(5)-1)  ! 187 = E-E
    IF (sumgrp(8)/=0)  nangroup(188)=sumgrp(5)*sumgrp(8)*2    ! 188 = E-H
  ENDIF  

  IF (sumgrp(6)/=0) THEN  ! ester interactions
    IF (sumgrp(6) > 1) nangroup(189)=sumgrp(6)*(sumgrp(6)-1)  ! 189 = F-F
    IF (sumgrp(7)/=0)  nangroup(190)=sumgrp(6)*sumgrp(7)*2    ! 190 = F-G
    IF (sumgrp(17)/=0) nangroup(197)=sumgrp(6)*sumgrp(17)*2   ! 191 = F-Q
    IF (sumgrp(9)/=0)  nangroup(193)=sumgrp(6)*sumgrp(9)*2    ! 193 = F-I
  ENDIF

  IF (sumgrp(7)/=0) THEN  ! ketones interactions
    IF (sumgrp(7) > 1) nangroup(194)=sumgrp(7)*(sumgrp(7)-1)  ! 194 = G-G
    IF (sumgrp(8)/=0)  nangroup(195)=sumgrp(8)*sumgrp(7)*2    ! 195 = G-H
    IF (sumgrp(17)/=0) nangroup(196)=sumgrp(17)*sumgrp(7)*2   ! 196 = G-Q
    IF (sumgrp(9)/=0)  nangroup(199)=sumgrp(9)*sumgrp(7)*2    ! 199 = G-I
  ENDIF

  IF (sumgrp(8)/=0) THEN  ! aldehydes interactions
    IF (sumgrp(8) > 1) nangroup(201)=sumgrp(8)*(sumgrp(8)-1)  ! 201 = H-H
    IF (sumgrp(17)/=0) nangroup(202)=sumgrp(8)*sumgrp(17)*2   ! 202 = H-Q
    IF (sumgrp(9)/=0)  nangroup(204)=sumgrp(8)*sumgrp(9)*2    ! 204 = H-I
  ENDIF
  
  IF (sumgrp(17)/=0) THEN  ! nitro interactions
    IF (sumgrp(17)>1) nangroup(206)=sumgrp(17)*(sumgrp(17)-1) ! 201 = H-H
  ENDIF

END SUBROUTINE gi4tb


!=======================================================================
! PURPOSE: Compute group interaction for the Nannoolal SAR for Vapor P.
! Groupe definition for the Tb and Vapor Pressure are not strictly 
! identical - need to be computed specifically.
!=======================================================================
SUBROUTINE gi4vp(Nangroup)
  IMPLICIT NONE
  INTEGER,INTENT(INOUT) :: nangroup(:)
  
  INTEGER :: sumgrp(19)

! clean any previous GI data
  Nangroup(135:212)=0

! Group definition for Vapor Pressure SAR (but 37 not only in B)
! ----------------------------------------
  sumgrp(1)=nangroup(34)+nangroup(35)+nangroup(36)+ &   ! group A (OH)
            nangroup(214)                               ! OOH added to A
  sumgrp(2)=nangroup(37)                                ! group B (OH aromatic)
  sumgrp(3)=nangroup(44)+nangroup(219)                  ! group C (CO(OH) with CO(OOH) added)
  sumgrp(4)=nangroup(38)                                ! group D (Ether)
  sumgrp(5)=nangroup(39)                                ! group E (Epox)
  sumgrp(6)=nangroup(45)+nangroup(46)+nangroup(47)      ! group F (Ester)
  sumgrp(7)=nangroup(51)+nangroup(92)                   ! group G (ketone)
  sumgrp(8)=nangroup(52)+nangroup(90)                   ! group H (aldehyde)
  sumgrp(9)=nangroup(65)                                ! group I ether (aromatic)
  sumgrp(10)=nangroup(54)                               ! group J sulfide (not used)
  sumgrp(11)=nangroup(56)                               ! group K aromatic sufur (not used)
  sumgrp(12)=nangroup(53)                               ! group L (SH thiol) (not used)
  sumgrp(13)=nangroup(40)+nangroup(41)                  ! group M primary amine (not used)
  sumgrp(14)=nangroup(42)+nangroup(97)                  ! group N secondary amine (not used)
  sumgrp(15)=nangroup(80)                               ! group O isocyanate -OCN (not used)
  sumgrp(16)=nangroup(57)                               ! group P cyanide -CN (not used)
  sumgrp(17)=nangroup(69)                               ! group Q nitrate 
  sumgrp(18)=nangroup(66)                               ! group R aromatic N cyclic 5 (not used)
  sumgrp(19)=nangroup(67)                               ! group S aromatic N cyclic 6 (not used)


! ----------------------------
! GROUP INTERACTION FOR VP SAR
! ----------------------------
  IF (sumgrp(1)/=0) THEN    ! OH interactions
    IF (sumgrp(1) >1)  nangroup(135)=sumgrp(1)*(sumgrp(1)-1)  ! 135 = A-A 
    IF (sumgrp(4)/=0)  nangroup(140)=sumgrp(1)*sumgrp(4)*2    ! 140 = A-D
    IF (sumgrp(5)/=0)  nangroup(141)=sumgrp(1)*sumgrp(5)*2    ! 141 = A-E
    IF (sumgrp(6)/=0)  nangroup(142)=sumgrp(1)*sumgrp(6)*2    ! 142 = A-F
    IF (sumgrp(7)/=0)  nangroup(143)=sumgrp(1)*sumgrp(7)*2    ! 143 = A-G
    IF (sumgrp(9)/=0)  nangroup(146)=sumgrp(1)*sumgrp(9)*2    ! 146 = A-I
  ENDIF

  IF (sumgrp(2)/=0) THEN    ! OH AROMATIC interactions (NOTE 37 BELONG TO A AND B ?)
    IF (sumgrp(2) >1)  nangroup(148)=sumgrp(2)*(sumgrp(2)-1)  ! 148 = B-B 
    IF (sumgrp(4)/=0)  nangroup(151)=sumgrp(2)*sumgrp(4)*2    ! 151 = B-D
    IF (sumgrp(6)/=0)  nangroup(152)=sumgrp(2)*sumgrp(6)*2    ! 152 = B-F
  ENDIF

  IF (sumgrp(3)/=0) THEN  ! CO(OH) interactions
    IF (sumgrp(3) > 1) nangroup(172)=sumgrp(3)*(sumgrp(3)-1)  ! 172 = C-C
    IF (sumgrp(7)/=0)  nangroup(175)=sumgrp(3)*sumgrp(7)*2    ! 175 = C-G
  ENDIF

  IF (sumgrp(4)/=0) THEN  ! -O- interactions
    IF (sumgrp(4) > 1) nangroup(178)=sumgrp(4)*(sumgrp(4)-1)  ! 178 = D-D
    IF (sumgrp(5)/=0)  nangroup(179)=sumgrp(4)*sumgrp(5)*2    ! 179 = D-E
    IF (sumgrp(6)/=0)  nangroup(180)=sumgrp(4)*sumgrp(6)*2    ! 180 = D-F
    IF (sumgrp(7)/=0)  nangroup(181)=sumgrp(4)*sumgrp(7)*2    ! 181 = D-G
    IF (sumgrp(8)/=0)  nangroup(182)=sumgrp(4)*sumgrp(8)*2    ! 182 = D-H
  ENDIF

  IF (sumgrp(5)/=0) THEN  ! -O- interactions
    IF (sumgrp(5) > 1) nangroup(187)=sumgrp(5)*(sumgrp(5)-1)  ! 187 = E-E
  ENDIF  

  IF (sumgrp(6)/=0) THEN  ! ester interactions
    IF (sumgrp(6) > 1) nangroup(189)=sumgrp(6)*(sumgrp(6)-1)  ! 189 = F-F
    IF (sumgrp(7)/=0)  nangroup(190)=sumgrp(6)*sumgrp(7)*2    ! 190 = F-G
    IF (sumgrp(9)/=0)  nangroup(193)=sumgrp(6)*sumgrp(9)*2    ! 193 = F-I
  ENDIF

  IF (sumgrp(7)/=0) THEN  ! ketones interactions
    IF (sumgrp(7) > 1) nangroup(194)=sumgrp(7)*(sumgrp(7)-1)  ! 194 = G-G
  ENDIF

  IF (sumgrp(8)/=0) THEN  ! aldehydes interactions
    IF (sumgrp(8) > 1) nangroup(201)=sumgrp(8)*(sumgrp(8)-1)  ! 201 = H-H
    IF (sumgrp(9)/=0)  nangroup(204)=sumgrp(8)*sumgrp(9)*2    ! 204 = H-I
  ENDIF
  
  IF (sumgrp(17)/=0) THEN  ! nitro interactions
    IF (sumgrp(17)>1) nangroup(206)=sumgrp(17)*(sumgrp(17)-1) ! 206 = Q-Q
  ENDIF


END SUBROUTINE gi4vp

END MODULE nannoolaltool
