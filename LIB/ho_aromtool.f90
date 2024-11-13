MODULE ho_aromtool
IMPLICIT NONE
CONTAINS
!=======================================================================
! PURPOSE: perform the OH+aromatic reaction based on Jenkin et al.,
! ACP, 2018
!=======================================================================
SUBROUTINE hoadd_arom(chem,bond,group,nca,nr,flag,tarrhc,pchem,coprod,nrxref,rxref)
  USE keyparameter, ONLY: mxcp,saru
  USE references, ONLY: mxlcod
  USE keyflag, ONLY: wrtsarinfo
  USE mapping, ONLY: gettrack
  USE reactool, ONLY: swap, rebond
  USE normchem, ONLY: stdchm
  USE radchktool, ONLY: radchk
  USE toolbox, ONLY: setbond
  USE toolbox, ONLY: stoperr,addrx,addref
  IMPLICIT NONE

  CHARACTER(LEN=*),INTENT(IN) :: chem      ! chemical formula
  INTEGER,INTENT(IN)          :: bond(:,:) ! bond matrix
  CHARACTER(LEN=*),INTENT(IN) :: group(:)  ! group matrix
  INTEGER,INTENT(IN)          :: nca       ! # of groups
  INTEGER,INTENT(INOUT) :: nr         ! # of reactions attached to chem (incremented here)
  INTEGER,INTENT(INOUT) :: flag(:)    ! flag to activate reaction i
  REAL,INTENT(INOUT) :: tarrhc(:,:)   ! arrhenius coefficient for reaction i
  CHARACTER(LEN=*),INTENT(INOUT) :: pchem(:)    ! "main" product of reaction channel i
  CHARACTER(LEN=*),INTENT(INOUT) :: coprod(:,:) ! coproducts (short names) for reaction channel i
  INTEGER,INTENT(INOUT)          :: nrxref(:)       ! # of references added for each reaction
  CHARACTER(LEN=*),INTENT(INOUT) :: rxref(:,:)      ! references (short codes) for the reactions

  INTEGER :: i,j,k,l,num_sub,jnext
  REAL    :: karom(3),kipso(3)
  REAL    :: kabs(3),kadd1(3),kadd2(3)
  REAL    :: kabs_298,kadd1_298,kadd2_298
  INTEGER :: o_sub(2),m_sub(2),p_sub(2)
  INTEGER :: path(6),nring
  REAL    :: frac(3),sum298,arrhc(3)
  INTEGER :: nb_alp,alp(4),nsub
  INTEGER :: track(mxcp,SIZE(group))
  INTEGER :: trlen(mxcp)
  INTEGER :: ntr,nc,nch3

  CHARACTER(LEN=LEN(group(1))) :: tgroup(SIZE(group)), pold, pnew
  CHARACTER(LEN=LEN(chem)) :: tempkc
  INTEGER :: tbond(SIZE(bond,1),SIZE(bond,1))

  INTEGER :: nip
  LOGICAL :: lo_ester

  INTEGER,PARAMETER :: mxrpd=2  ! max # of products returned by radchk sub
  CHARACTER(LEN=LEN(chem))  :: rdckpd(mxrpd)
  CHARACTER(LEN=LEN(coprod(1,1))) :: rdckcopd(mxrpd,SIZE(coprod,2))
  REAL :: sc(mxrpd)

  CHARACTER(LEN=12),PARAMETER :: progname='hoadd_arom '
  CHARACTER(LEN=70)          :: mesg

  karom(1)=0.378E-12 ; karom(2)=0. ; karom(3)=190.  ! Table 2 in Jenkin et al., 2018
  kipso(1)=0.378E-12 ; kipso(2)=0. ; kipso(3)=89.   ! Table 2 in Jenkin et al., 2018
  frac(:)=0. ; arrhc(:)=0.

  IF (wrtsarinfo) THEN                                           !! debug
    WRITE(saru,*) '   '                                           !! debug
    WRITE(saru,*) '   '                                           !! debug
    WRITE(saru,*) '============== OH ADD AROM  ============== '   !! debug
    WRITE(saru,*)                                                 !! debug
    DO i=1,nca ; WRITE(saru,*) i,group(i) ; ENDDO                 !! debug
  ENDIF

! compute addition rate constant on each carbon of the aromatic ring      
  DO i=1,nca
    IF (group(i)(1:1)=='c') THEN
      CALL arom_data(i,group,bond,nca,o_sub,m_sub,p_sub)
      CALL gettrack(bond,i,nca,ntr,track,trlen)
      CALL arom_path(i,ntr,track,group,bond,path)          
      
      nb_alp=0  ;  alp(:)=0
      DO j=1,nca
        IF ((group(j)(1:1)=='c').AND.(bond(i,j)==1)) THEN
          nb_alp=nb_alp+1  ;  alp(nb_alp)=j
        ENDIF
      ENDDO

!--------------------------------
!- COMPUTE RATE CONSTANTS 
!--------------------------------

      IF (group(i)(1:3)=='cH ') THEN ! addition on a non substitued C
        arrhc(1)=karom(1) ; arrhc(3)=karom(3)
      ELSE                           ! addition on a substitued C
        arrhc(1)=kipso(1) ; arrhc(3)=kipso(3)
      ENDIF

      IF (group(path(2))(1:2)=='c(') THEN
        IF (o_sub(1)==0) THEN ; o_sub(1)=path(2)
        ELSE                  ; o_sub(2)=path(2)
        ENDIF
      ENDIF
      IF (group(path(3))(1:2)=='c(') THEN
        IF (m_sub(1)==0) THEN ; m_sub(1)=path(3)
        ELSE                  ; m_sub(2)=path(3)
        ENDIF
      ENDIF
      IF (group(path(4))(1:2)=='c(') THEN
        IF (p_sub(1)==0) THEN ; p_sub(1)=path(4)
        ELSE                  ; p_sub(2)=path(4)
        ENDIF
      ENDIF
      IF (group(path(5))(1:2)=='c(') THEN
        IF (m_sub(1)==0) THEN ; m_sub(1)=path(5)
        ELSE                  ; m_sub(2)=path(5)
        ENDIF
      ENDIF
      IF (group(path(6))(1:2)=='c(') THEN
        IF (o_sub(1)==0) THEN ; o_sub(1)=path(6)
        ELSE                  ; o_sub(2)=path(6)
        ENDIF
      ENDIF

      num_sub=0
      DO j=1,2
        IF (o_sub(j)/=0) num_sub=num_sub+1
        IF (m_sub(j)/=0) num_sub=num_sub+1
        IF (p_sub(j)/=0) num_sub=num_sub+1
      ENDDO
      
      IF (wrtsarinfo) THEN                                              !! debug
        WRITE(saru,*) 'Nb subst:', num_sub                               !! debug
        DO j=1,2                                                         !! debug
          IF (o_sub(j)/=0) WRITE(saru,*) '-orto(',j,')=',group(o_sub(j)) !! debug
          IF (m_sub(j)/=0) WRITE(saru,*) '-meta(',j,')=',group(m_sub(j)) !! debug
          IF (p_sub(j)/=0) WRITE(saru,*) '-para(',j,')=',group(p_sub(j)) !! debug
        ENDDO                                                            !! debug
      ENDIF                                                              !! debug

!--- SUBSTITUENT FACTORS F() --- Table 3 in Jenkin et al., 2018                     
      IF (num_sub==1) THEN
        IF (m_sub(1)/=0) THEN ; arrhc(1)=arrhc(1)*0.7 ; arrhc(3)=arrhc(3)-207.
        ELSE                  ; arrhc(1)=arrhc(1)*0.8 ; arrhc(3)=arrhc(3)-659.
        ENDIF

      ELSE IF (num_sub==2) THEN
        IF ((m_sub(1)/=0).AND.(m_sub(2)/=0)) THEN
          arrhc(1)=arrhc(1)*1.9 ; arrhc(3)=arrhc(3)-409.
        ELSE IF ( ((o_sub(1)/=0).AND.(o_sub(2)/=0)) .OR. &
                  ((o_sub(1)/=0).AND.(p_sub(1)/=0)) ) THEN
          arrhc(1)=arrhc(1)*0.6 ; arrhc(3)=arrhc(3)-1203.
        ELSE
          arrhc(1)=arrhc(1)*2.6 ; arrhc(3)=arrhc(3)-416.
        ENDIF

      ELSE IF (num_sub==3) THEN
        IF ( (o_sub(1)/=0).AND.(o_sub(2)/=0).AND.(p_sub(1)/=0) )THEN
          arrhc(1)=arrhc(1)*6.8 ; arrhc(3)=arrhc(3)-760.
        ELSE IF ( ((o_sub(1)/=0).AND.(o_sub(2)/=0).AND.(m_sub(1)/=0)) .OR. & 
                  ((o_sub(1)/=0).AND.(p_sub(1)/=0).AND.(m_sub(1)/=0)) ) THEN
          arrhc(1)=arrhc(1)*0.5 ; arrhc(3)=arrhc(3)-1200.
        ELSE
          arrhc(1)=arrhc(1)*3.5 ; arrhc(3)=arrhc(3)-341.
        ENDIF

      ELSE IF (num_sub==4) THEN
        IF ( (o_sub(1)/=0).AND.(o_sub(2)/=0).AND. &
             (m_sub(1)/=0).AND.(p_sub(1)/=0) ) THEN
          arrhc(1)=arrhc(1)*2.0 ; arrhc(3)=arrhc(3)-998.
        ELSE
          arrhc(1)=arrhc(1)*0.3 ; arrhc(3)=arrhc(3)-1564.
        ENDIF

      ELSE IF (num_sub==5) THEN
        arrhc(1)=arrhc(1)*4.7 ; arrhc(3)=arrhc(3)-809.
      ENDIF
          
! reset substituent to only keep carbon substituents
      CALL arom_data(i,group,bond,nca,o_sub,m_sub,p_sub)
          
!--- SUBSTITUENT FACTORS R()---  Table 6 in Jenkin et al., 2018

! ortho position
      loop_o: DO j=1,2
        IF (o_sub(j)/=0) THEN
          IF (group(o_sub(j))(1:4)=='CH2 ') THEN
            DO k=1,nca
             IF ((bond(o_sub(j),k)==1).AND.(group(k)(1:3)=='CH3'))THEN
               arrhc(1)=arrhc(1)*0.029 ; arrhc(3)=arrhc(3)-1014 ;  CYCLE loop_o
             ELSE IF ((bond(o_sub(j),k)==1).AND.(group(k)(1:4)=='CH2 ')) THEN
               DO l=1,nca
                 IF ((bond(l,k)==1).AND.(group(l)(1:3)=='CH3')) THEN
                   arrhc(1)=arrhc(1)*0.029 ; arrhc(3)=arrhc(3)-1000  ;  CYCLE loop_o
                 ENDIF
               ENDDO
             ENDIF
            ENDDO
          ELSE IF (group(o_sub(j))(1:3)=='CH ') THEN
            nch3=0
            DO k=1,nca
              IF ((bond(o_sub(j),k)==1).AND.(group(k)(1:3)=='CH3')) THEN
                nch3=nch3+1
              ENDIF
            ENDDO
            IF (nch3==2) THEN
              arrhc(1)=arrhc(1)*0.029 ; arrhc(3)=arrhc(3)-1000
            ENDIF
          ELSE IF (group(o_sub(j))(1:2)=='C ') THEN
            nch3=0
            DO k=1,nca
              IF ((bond(o_sub(j),k)==1).AND.(group(k)(1:3)=='CH3'))THEN
                nch3=nch3+1
              ENDIF
            ENDDO
            IF (nch3==3) THEN
              arrhc(1)=arrhc(1)*0.029 ; arrhc(3)=arrhc(3)-957
            ENDIF
          ENDIF
        
          IF ( (INDEX(group(o_sub(j)),'(OH)')/=0).AND. &
               (group(o_sub(j))(1:6)/='CO(OH)') ) THEN
            arrhc(3)=arrhc(3)-390
          ENDIF
          
          IF (group(o_sub(j))(1:3)=='CHO') arrhc(3)=arrhc(3)+698
          IF (group(o_sub(j))(1:3)=='-O-') arrhc(3)=arrhc(3)-365
          IF (group(o_sub(j))(1:3)=='CO ') THEN
            lo_ester=.FALSE.
            DO k=1,nca
              IF (bond(k,o_sub(j))==3) THEN
                lo_ester=.TRUE. ;  EXIT
              ENDIF
            ENDDO
            IF (lo_ester) THEN ; arrhc(3)=arrhc(3)+400
            ELSE               ; arrhc(3)=arrhc(3)+698
            ENDIF
          ENDIF 
  
        ENDIF
      ENDDO loop_o
      IF ((group(path(2))(1:5)=='c(OH)').OR.(group(path(2))(1:6)=='c(OOH)')) THEN 
        arrhc(1)=arrhc(1)*0.69 ; arrhc(3)=arrhc(3)-395
      ENDIF
      IF ((group(path(6))(1:5)=='c(OH)').OR.(group(path(6))(1:6)=='c(OOH)')) THEN 
        arrhc(1)=arrhc(1)*0.69 ; arrhc(3)=arrhc(3)-395
        ENDIF
      IF ((group(path(2))(1:6)=='c(NO2)').OR.(group(path(2))(1:7)=='c(ONO2)')) arrhc(3)=arrhc(3)+1110
      IF ((group(path(6))(1:6)=='c(NO2)').OR.(group(path(6))(1:7)=='c(ONO2)')) arrhc(3)=arrhc(3)+1110
          
! para position
      IF (p_sub(1)/=0) THEN
        IF (group(p_sub(1))(1:4)=='CH2 ') THEN
          DO k=1,nca
            IF ((bond(p_sub(1),k)==1).AND.(group(k)(1:3)=='CH3'))THEN
             arrhc(1)=arrhc(1)*0.029 ; arrhc(3)=arrhc(3)-1014 ; EXIT
           ELSE IF ((bond(p_sub(1),k)==1).AND.(group(k)(1:4)=='CH2 ')) THEN
             DO l=1,nca
               IF ((bond(l,k)==1).AND.(group(l)(1:3)=='CH3')) THEN 
                 arrhc(1)=arrhc(1)*0.029 ; arrhc(3)=arrhc(3)-1000 ; EXIT
               ENDIF
             ENDDO                 
           ENDIF
          ENDDO
        ELSE IF (group(p_sub(1))(1:3)=='CH ') THEN
          nch3=0
          DO k=1,nca
            IF ((bond(p_sub(1),k)==1).AND.(group(k)(1:3)=='CH3')) THEN
              nch3=nch3+1
            ENDIF
          ENDDO
          IF (nch3==2) THEN
            arrhc(1)=arrhc(1)*0.029 ; arrhc(3)=arrhc(3)-1000
          ENDIF
        ELSE IF (group(p_sub(1))(1:2)=='C ') THEN
          nch3=0
          DO k=1,nca
            IF ((bond(p_sub(1),k)==1).AND.(group(k)(1:3)=='CH3')) THEN
              nch3=nch3+1
            ENDIF
          ENDDO
          IF (nch3==3) THEN
            arrhc(1)=arrhc(1)*0.029 ; arrhc(3)=arrhc(3)-957
          ENDIF
        ENDIF
     
        IF ( (INDEX(group(p_sub(1)),'(OH)')/=0).AND. &
             (group(p_sub(1))(1:6)/='CO(OH)') ) THEN
          arrhc(3)=arrhc(3)-390
        ENDIF
     
        IF (group(p_sub(1))(1:3)=='CHO') arrhc(3)=arrhc(3)+698
        IF (group(p_sub(1))(1:3)=='-O-') arrhc(3)=arrhc(3)-365
        IF (group(p_sub(1))(1:3)=='CO ') THEN
          lo_ester=.FALSE.
          DO k=1,nca
            IF (bond(k,p_sub(1))==3) THEN
              lo_ester=.TRUE. ;  EXIT
            ENDIF
          ENDDO
          IF (lo_ester) THEN ; arrhc(3)=arrhc(3)+400
          ELSE               ; arrhc(3)=arrhc(3)+698
          ENDIF
        ENDIF 
     
      ENDIF
      IF ((group(path(4))(1:5)=='c(OH)').OR.(group(path(4))(1:6)=='c(OOH)')) THEN 
        arrhc(1)=arrhc(1)*0.69 ; arrhc(3)=arrhc(3)-395
      ENDIF
      IF ((group(path(4))(1:6)=='c(NO2)').OR.(group(path(4))(1:7)=='c(ONO2)')) arrhc(3)=arrhc(3)+1110

! meta position
      loop_m : DO j=1,2
        IF (m_sub(j)/=0) THEN
 
          IF ( (INDEX(group(m_sub(j)),'(OH)')/=0).AND. &
               (group(m_sub(j))(1:6)/='CO(OH)')) THEN
            arrhc(3)=arrhc(3)-390
          ENDIF
        
          IF (group(m_sub(j))(1:3)=='CHO') arrhc(3)=arrhc(3)+698
          IF (group(m_sub(j))(1:3)=='-O-') arrhc(3)=arrhc(3)+70
          IF (group(m_sub(j))(1:3)=='CO ') THEN
            lo_ester=.FALSE.
            DO k=1,nca
              IF (bond(k,m_sub(j))==3) THEN
                lo_ester=.TRUE. ;  EXIT
              ENDIF
            ENDDO
            IF (lo_ester) THEN ; arrhc(3)=arrhc(3)+400
            ELSE               ; arrhc(3)=arrhc(3)+698
            ENDIF
          ENDIF 
 
        ENDIF
      ENDDO loop_m
      IF ((group(path(3))(1:5)=='c(OH)').OR.(group(path(3))(1:6)=='c(OOH)')) THEN 
        arrhc(1)=arrhc(1)*0.025 ; arrhc(3)=arrhc(3)-1360
      ENDIF
      IF ((group(path(5))(1:5)=='c(OH)').OR.(group(path(5))(1:6)=='c(OOH)')) THEN 
        arrhc(1)=arrhc(1)*0.025 ; arrhc(3)=arrhc(3)-1360
      ENDIF
      IF ((group(path(3))(1:6)=='c(NO2)').OR.(group(path(3))(1:7)=='c(ONO2)')) arrhc(3)=arrhc(3)+792
      IF ((group(path(5))(1:6)=='c(NO2)').OR.(group(path(5))(1:7)=='c(ONO2)')) arrhc(3)=arrhc(3)+792
          
! ipso position
      IF ((group(path(1))(1:5)=='c(OH)').OR.(group(path(1))(1:6)=='c(OOH)')) THEN 
        arrhc(1)=arrhc(1)*0.025 ; arrhc(3)=arrhc(3)-1360
      ENDIF
      IF ((group(path(1))(1:6)=='c(NO2)').OR.(group(path(1))(1:7)=='c(ONO2)')) arrhc(3)=arrhc(3)+792
      
      DO j=1,nca
        IF ((bond(j,path(1))/=0).AND.(group(j)(1:1)/='c')) THEN
          IF ( (INDEX(group(j),'(OH)')/=0).AND. &
               (group(j)(1:6)/='CO(OH)') ) THEN
            arrhc(3)=arrhc(3)-390
          ENDIF
          IF (group(j)(1:3)=='CHO') arrhc(3)=arrhc(3)+698
          IF (group(j)(1:3)=='-O-') arrhc(3)=arrhc(3)+70
          IF (group(j)(1:3)=='CO ') THEN
            lo_ester=.FALSE.
            DO k=1,nca
              IF (bond(k,j)==3) THEN
                lo_ester=.TRUE. ; EXIT
              ENDIF
            ENDDO
            IF (lo_ester) THEN ; arrhc(3)=arrhc(3)+400
            ELSE               ; arrhc(3)=arrhc(3)+698
            ENDIF
          ENDIF
        ENDIF
      ENDDO

! alk-1-enyl
      Cd_loop: DO j=1,6
        DO k=1,nca
          IF ((bond(path(j),k)==1).AND.(group(k)(1:4)=='CdH ')) THEN 
            DO l=1,nca
              IF ((bond(l,k)==2).AND.(group(l)(1:4)=='CdH2')) THEN
                arrhc(:)=0. ;  EXIT Cd_loop
              ENDIF
            ENDDO
          ENDIF              
        ENDDO
      ENDDO Cd_loop

      IF (wrtsarinfo) THEN                                     !! debug
        WRITE(saru,*) 'addition on group :',i                    !! debug
        WRITE(saru,*) 'k 298:', arrhc(1)*EXP(- arrhc(3)/298.)     !! debug
      ENDIF                                                     !! debug

!--------------------------------
!- REACTIONS 
!--------------------------------

! to this point, OH has been added on the ring, forming a delocalized
! radical. This radical can follow three reactions with O2: 
! 1- H abstraction if OH was added on a cH group
! 2 and 3 : addition of O2 in alpha position of the OH bearing carbon
! Rates in table 7 and substituent factors in table 8 in Jenkin et al., 2018
          
! CASE 1
      kabs(:)=0.
      IF (group(i)(1:3)=='cH ') THEN
        kabs(1)=1.75E-14 ; kabs(2)=0. ;  kabs(3)=1500.
      ENDIF

! CASE 2:
      CALL gettrack(bond,alp(1),nca,ntr,track,trlen)
      CALL arom_path(i,ntr,track,group,bond,path)

      kadd1(1)=1.50E-14
      kadd1(2)=0.
      kadd1(3)=1700.
      IF (group(path(2))(1:3)/='cH ') kadd1(3)=kadd1(3)-207.
      IF (group(path(3))(1:3)/='cH ') kadd1(3)=kadd1(3)-620.
      IF (group(path(4))(1:3)/='cH ') kadd1(3)=kadd1(3)-207.
      IF (group(path(5))(1:3)/='cH ') kadd1(3)=kadd1(3)-558.
      nsub=0
      IF (group(path(2))(1:3)/='cH ') nsub=nsub+1
      IF (group(path(3))(1:3)/='cH ') nsub=nsub+1
      IF (group(path(4))(1:3)/='cH ') nsub=nsub+1
      IF (group(path(5))(1:3)/='cH ') nsub=nsub+1
      IF (nsub>1) kadd1(1)=kadd1(1)/(nsub**0.5)

! CASE 3:
      CALL gettrack(bond,alp(2),nca,ntr,track,trlen)
      CALL arom_path(i,ntr,track,group,bond,path)

      kadd2(1)=1.50E-14 ; kadd2(2)=0. ; kadd2(3)=1700.
      IF (group(path(2))(1:3)/='cH ') kadd2(3)=kadd2(3)-207.
      IF (group(path(3))(1:3)/='cH ') kadd2(3)=kadd2(3)-620.
      IF (group(path(4))(1:3)/='cH ') kadd2(3)=kadd2(3)-207.
      IF (group(path(5))(1:3)/='cH ') kadd2(3)=kadd2(3)-558.
      nsub=0
      IF (group(path(2))(1:3)/='cH ') nsub=nsub+1
      IF (group(path(3))(1:3)/='cH ') nsub=nsub+1
      IF (group(path(4))(1:3)/='cH ') nsub=nsub+1
      IF (group(path(5))(1:3)/='cH ') nsub=nsub+1
      IF (nsub>1) kadd2(1)=kadd2(1)/(nsub**0.5)

! compute fraction for three cases
      kabs_298=kabs(1)*EXP(-kabs(3)/298.)
      kadd1_298=kadd1(1)*EXP(-kadd1(3)/298.)
      kadd2_298=kadd2(1)*EXP(-kadd2(3)/298.)

      sum298=kabs_298+kadd1_298+kadd2_298

      frac(1)=kabs_298/sum298
      frac(2)=kadd1_298/sum298
      frac(3)=kadd2_298/sum298

      IF (wrtsarinfo) WRITE(saru,*) 'frac(1:3) (%): ',INT(100*frac(:)) !! debug

!--- WRITE NEW SPECIES

! CASE 1  
      tgroup(:)=group(:) ;  tbond(:,:)=bond(:,:)
	  
      IF (tgroup(i)(1:3)=='cH ') THEN
        CALL addrx(progname,chem,nr,flag)
        tarrhc(nr,1)=arrhc(1)*frac(1) 
        tarrhc(nr,2)=0.
        tarrhc(nr,3)=arrhc(3)
        CALL addref(progname,'OHAROSAR1',nrxref(nr),rxref(nr,:),chem)
        CALL addref(progname,'MJ18KMV001',nrxref(nr),rxref(nr,:),chem)  
	  
        CALL gettrack(bond,alp(1),nca,ntr,track,trlen)
        CALL arom_path(i,ntr,track,group,bond,path)
        tgroup(i)='c(OH)'
	  
        CALL rebond(tbond,tgroup,pchem(nr),nring)
        CALL stdchm(pchem(nr))
        tgroup(:)=group(:) ;  tbond(:,:)=bond(:,:)

        coprod(nr,1)='HO2'
      ENDIF

! CASE 2
      DO k=1,2
        CALL addrx(progname,chem,nr,flag)
        tarrhc(nr,1)=arrhc(1)*frac(k+1)
        tarrhc(nr,2)=0
        tarrhc(nr,3)=arrhc(3)
        CALL addref(progname,'OHAROSAR1',nrxref(nr),rxref(nr,:),chem)
        CALL addref(progname,'MJ18KMV001',nrxref(nr),rxref(nr,:),chem)  

        IF (tgroup(i)(1:3)=='cH ') THEN 
          tgroup(i)='CH(OH)'
        ELSE
          nc=INDEX(tgroup(i),' ') ; tgroup(i)(nc:nc+3)='(OH)'
          tgroup(i)(1:1)='C'
        ENDIF
        CALL gettrack(bond,alp(k),nca,ntr,track,trlen)
        CALL arom_path(i,ntr,track,group,bond,path)
        nc=INDEX(tgroup(path(1)),' ') ; tgroup(path(1))(nc:nc+4)='(OO.)'
        tgroup(path(1))(1:1)='C'
        pold='c' ; pnew='Cd'
        DO j=2,5
          CALL swap(group(path(j)),pold,tgroup(path(j)),pnew)
        ENDDO
	  
        CALL setbond(tbond,i,path(1),1)
        CALL setbond(tbond,path(1),path(2),1)
        CALL setbond(tbond,path(2),path(3),2)
        CALL setbond(tbond,path(3),path(4),1)
        CALL setbond(tbond,path(4),path(5),2)
        CALL setbond(tbond,path(5),i,1)
        CALL rebond(tbond,tgroup,tempkc,nring)
        CALL radchk(tempkc,rdckpd,rdckcopd,nip,sc,nrxref(nr),rxref(nr,:))
        pchem(nr) = rdckpd(1)
        CALL stdchm(pchem(nr))
        IF (nip/=1) THEN
          mesg="unexpected 2 pdcts from radchk "
          CALL stoperr(progname,mesg,chem)
        ENDIF
        IF (rdckcopd(1,1)/=' ') THEN      ! add coproducts
          copdloop: DO j=1,SIZE(coprod,2)
            IF (coprod(1,j)==' ') THEN    ! j is the 1st available slot
              jnext=j
              DO l=1,SIZE(rdckcopd,2)
                IF (rdckcopd(1,l)/=' ') THEN
                  coprod(nr,jnext)=rdckcopd(1,l)
                  jnext=jnext+1
                ENDIF
              ENDDO
              EXIT copdloop
            ENDIF
          ENDDO copdloop
        ENDIF

        tgroup(:)=group(:) ; tbond(:,:)=bond(:,:)
      ENDDO

    ENDIF
  ENDDO

END SUBROUTINE hoadd_arom

!=======================================================================
! Purpose: analyze an aromatic cycle and returns the index of group 
! substituent in o, m and para position, beginning to top)
!=======================================================================
SUBROUTINE arom_data(top,group,bond,nca,o_sub,m_sub,p_sub)
  USE keyparameter, ONLY: mxcp
  USE mapping, ONLY: gettrack

  IMPLICIT NONE

  CHARACTER(LEN=*),INTENT(IN) :: group(:)
  INTEGER,INTENT(IN)  :: bond(:,:)
  INTEGER,INTENT(IN)  :: top,nca
  INTEGER,INTENT(OUT) :: o_sub(2),m_sub(2),p_sub(2)

  INTEGER :: no_sub,nm_sub,np_sub
  INTEGER :: track(mxcp,SIZE(group))
  INTEGER :: trlen(mxcp)
  INTEGER :: ntr,k
  
  o_sub(:)=0  ;  m_sub(:)=0  ;  p_sub(:)=0
  no_sub=0    ;  nm_sub=0    ;  np_sub=0

  CALL gettrack(bond,top,nca,ntr,track,trlen)

  DO k=1,ntr
    IF (group(track(k,2))(1:1)=='C') THEN
!      np_sub=np_sub+1
!      p_sub(np_sub)=track(k,2)
    ELSE IF (group(track(k,2))(1:1)=='c') THEN
      IF ((group(track(k,3))(1:1)=='C').OR.(group(track(k,3))(1:3)=='-O-')) THEN
        IF ((o_sub(1)/=track(k,3)).AND.(o_sub(2)/=track(k,3))) THEN
          no_sub=no_sub+1 ; o_sub(no_sub)=track(k,3)
        ENDIF
      ELSE IF (group(track(k,3))(1:1)=='c') THEN
        IF ((group(track(k,4))(1:1)=='C').OR.(group(track(k,4))(1:3)=='-O-')) THEN
          IF ((m_sub(1)/=track(k,4)).AND.(m_sub(2)/=track(k,4))) THEN
            nm_sub=nm_sub+1 ; m_sub(nm_sub)=track(k,4)
          ENDIF
        ELSE IF (group(track(k,4))(1:1)=='c') THEN
          IF ((group(track(k,5))(1:1)=='C').OR.(group(track(k,5))(1:3)=='-O-')) THEN
            IF ((p_sub(1)/=track(k,5)).AND.(p_sub(2)/=track(k,5))) THEN
              np_sub=np_sub+1 ; p_sub(np_sub)=track(k,5)
            ENDIF
          ENDIF
        ENDIF
      ENDIF
    ENDIF
  ENDDO 
END SUBROUTINE arom_data

!=======================================================================
! Purpose:
!=======================================================================
SUBROUTINE arom_path(i,ntr,track,group,bond,path)
  IMPLICIT NONE

  INTEGER,INTENT(IN)  :: track(:,:)
  INTEGER,INTENT(IN)  :: bond(:,:)
  INTEGER,INTENT(IN)  :: ntr,i
  CHARACTER(LEN=*),INTENT(IN) :: group(:)
  INTEGER,INTENT(OUT) :: path(6)

  INTEGER             :: k

  path(:)=0

  DO k=1,ntr
    IF (track(k,6)==0) CYCLE
    IF ( (group(track(k,1))(1:1)=='c').AND. &
         (track(k,2)/=i).AND.               &
         (group(track(k,2))(1:1)=='c').AND. &
         (group(track(k,3))(1:1)=='c').AND. &
         (group(track(k,4))(1:1)=='c').AND. &
         (group(track(k,5))(1:1)=='c').AND. &
         (group(track(k,6))(1:1)=='c').AND. &
         (bond(track(k,1),track(k,6))==1)  ) THEN
      path(1:6)=track(k,1:6)
    ENDIF
  ENDDO
END SUBROUTINE arom_path

END MODULE ho_aromtool
