MODULE rdmecatool
IMPLICIT NONE
CONTAINS

!SUBROUTINE rdmeca()
!SUBROUTINE rdvocmchdb(filename,nkwdat,nkwpd,kwrct,kwpd,kwcopd,kwyld,kwcom)
!SUBROUTINE rdradmchdb(filename,nkr,krct,kprd,arrh,kcost,kcom)
!SUBROUTINE rdhvdb()

! ======================================================================
! PURPOSE: read known mechanism and store the data in the 
! database module
! ======================================================================
SUBROUTINE rdmeca()
  USE database, ONLY: &
      nkwoh, nkwoh_pd, kwoh_rct, kwoh_pd, kwoh_copd, kwoh_yld, kwoh_com, &
      nkwo3, nkwo3_pd, kwo3_rct, kwo3_pd, kwo3_copd, kwo3_yld, kwo3_com, &
      nkwno3, nkwno3_pd, kwno3_rct, kwno3_pd, kwno3_copd, kwno3_yld, kwno3_com, &
      nkwro2, kwro2_arrh, kwro2_stoi, kwro2_rct, kwro2_prd, kwro2_com, &
      nkwrco3, kwrco3_arrh, kwrco3_stoi, kwrco3_rct, kwrco3_prd, kwrco3_com, &
      nkwro, kwro_arrh, kwro_stoi, kwro_rct, kwro_prd, kwro_com,  &
      nkwcri, kwcri_arrh, kwcri_stoi, kwcri_rct, kwcri_prd, kwcri_com
  USE keyparameter, ONLY: dirgecko

  IMPLICIT NONE

  CHARACTER(LEN=100) :: filename

! VOC+OH 
  WRITE (6,*) '  ...reading known VOC+OH reactions ...'
  filename=TRIM(dirgecko)//'DATA/oh_prod.dat'
  CALL rdvocmchdb(filename,nkwoh,nkwoh_pd,kwoh_rct,kwoh_pd,kwoh_copd,kwoh_yld,kwoh_com)

! VOC+NO3 
  WRITE (6,*) '  ...reading known VOC+NO3 reactions ...'
  filename=TRIM(dirgecko)//'DATA/no3_prod.dat'
  CALL rdvocmchdb(filename,nkwno3,nkwno3_pd,kwno3_rct,kwno3_pd,kwno3_copd,kwno3_yld,kwno3_com)

! VOC+O3 
  WRITE (6,*) '  ...reading known VOC+O3 reactions ...'
  filename=TRIM(dirgecko)//'DATA/o3_prod.dat'
  CALL rdvocmchdb(filename,nkwo3,nkwo3_pd,kwo3_rct,kwo3_pd,kwo3_copd,kwo3_yld,kwo3_com)

! RO2 chemistry
  WRITE (6,*) '  ...reading known RO2 reactions ...'
  filename=TRIM(dirgecko)//'DATA/ro2.dat'
  CALL rdradmchdb(filename,nkwro2,kwro2_rct,kwro2_prd,kwro2_arrh,kwro2_stoi,kwro2_com)

! RCOO2 chemistry
  WRITE (6,*) '  ...reading known RO2 reactions ...'
  filename=TRIM(dirgecko)//'DATA/rcoo2.dat'
  CALL rdradmchdb(filename,nkwrco3,kwrco3_rct,kwrco3_prd,kwrco3_arrh,kwrco3_stoi,kwrco3_com)

! RO chemistry
  WRITE (6,*) '  ...reading known RO reactions ...'
  filename=TRIM(dirgecko)//'DATA/ro.dat'
  CALL rdradmchdb(filename,nkwro,kwro_rct,kwro_prd,kwro_arrh,kwro_stoi,kwro_com)

! Criegee chemistry
  WRITE (6,*) '  ...reading known criegee reactions ...'
  filename=TRIM(dirgecko)//'DATA/criegee.dat'
  CALL rdradmchdb(filename,nkwcri,kwcri_rct,kwcri_prd,kwcri_arrh,kwcri_stoi,kwcri_com)

! read known species for photolysis reactions
  CALL rdhvdb()

END SUBROUTINE rdmeca

! ======================================================================
! PURPOSE: Read mechanism data in the file provided as input. The
! routine is in particular called to read known VOC+OH, VOC+NO3, VOC+O3
! reactions.
! ======================================================================
SUBROUTINE rdvocmchdb(filename,nkwdat,nkwpd,kwrct,kwpd,kwcopd,kwyld,kwcom)
  USE keyparameter, ONLY: tfu1
  USE normchem, ONLY: stdchm
  USE sortstring, ONLY: sort_string   ! to sort chemical
  USE searching, ONLY: srh5           ! to seach in a sorted list
  USE references, ONLY: nreac_info, code  ! available code
  IMPLICIT NONE

  CHARACTER(LEN=*),INTENT(IN):: filename
  INTEGER,INTENT(OUT)        :: nkwdat          ! # of species in the database
  INTEGER,INTENT(OUT)        :: nkwpd(:)        ! # of channel involved in mechanism
  REAL,INTENT(OUT)           :: kwyld(:,:)      ! yield of each channel
  CHARACTER(LEN=*),INTENT(OUT) :: kwrct(:)      ! formula of the reacting VOC 
  CHARACTER(LEN=*),INTENT(OUT) :: kwpd(:,:)     ! formula of the main product per channel
  CHARACTER(LEN=*),INTENT(OUT) :: kwcopd(:,:,:) ! coproducts per channel
  CHARACTER(LEN=*)             :: kwcom(:,:)    ! code for references

  INTEGER, PARAMETER :: lenlin=200  ! length (max) of the line to be read
  CHARACTER(LEN=lenlin) :: line
  INTEGER  :: i,j,k,n1,n2,ndat,ncop,istart,iend,ierr,ilin
  INTEGER  :: ncpbeg, ncpend
  INTEGER  :: maxprd,maxcoprd,maxcom,maxent

! temporary copies of the table (before sorting)
  INTEGER                           :: tpnkwpd(SIZE(nkwpd))
  REAL                              :: tpkwyld(SIZE(kwyld,1),SIZE(kwyld,2))
  CHARACTER(LEN=LEN(kwrct(1)))      :: tpkwrct(SIZE(kwrct))
  CHARACTER(LEN=LEN(kwpd(1,1)))     :: tpkwpd(SIZE(kwpd,1),SIZE(kwpd,2))
  CHARACTER(LEN=LEN(kwcopd(1,1,1))) :: tpkwcopd(SIZE(kwcopd,1),SIZE(kwcopd,2),SIZE(kwcopd,3))
  CHARACTER(LEN=LEN(kwcom(1,1)))    :: tpkwcom(SIZE(kwcom,1),SIZE(kwcom,2))

  nkwdat = 0 ; nkwpd(:)=0 ; kwrct(:)=' ' ; kwpd(:,:)=' ' ; kwyld(:,:)=0. 
  kwcopd(:,:,:)=' ' ; kwcom(:,:)=' '

  tpkwrct(:)=' ' ; tpkwpd(:,:)=' ' ; tpkwyld(:,:)=0. ; 
  tpkwcopd(:,:,:)=' ' ; tpkwcom(:,:)=' '

  maxent=SIZE(nkwpd)
  maxprd=SIZE(kwcopd,2) ; maxcoprd=SIZE(kwcopd,3) ; maxcom=SIZE(kwcom,2)
  
  OPEN(tfu1,FILE=filename, FORM='FORMATTED',STATUS='OLD')

! read the data
  ilin=0
  rdloop: DO
    ilin=ilin+1
    READ (tfu1,'(a)',IOSTAT=ierr) line
    IF (ierr/=0) THEN
      WRITE(6,*) '--error--, while reading file: ',TRIM(filename)
      WRITE(6,*) 'at line number: ',ilin
      WRITE(6,*) 'keyword "END" missing ?'
      STOP "in rdvocmchdb"
    ENDIF     
    IF (line(1:3)=='END') EXIT rdloop
    IF (line(1:1)=='!') CYCLE rdloop

! check that the line is correctly formatted
    n1=INDEX(line,':')       
    IF (n1==0) THEN
      WRITE(6,*) '--error--, while reading file: ',TRIM(filename)
      WRITE (6,*) 'in rdvocmchdb, ":" not found at line:', ilin
      STOP "in rdvocmchdb"
    ENDIF

    ncpbeg=INDEX(line,';')       
    IF (ncpbeg==0) THEN
      WRITE(6,*) '--error--, while reading file: ',TRIM(filename)
      WRITE (6,*) 'in rdvocmchdb, ";" not found at line:', ilin
      STOP "in rdvocmchdb"
    ENDIF

    nkwdat=nkwdat+1
    IF (nkwdat>maxent) THEN
      WRITE(6,*) '--error--, while reading file: ',TRIM(filename)
      WRITE(6,*) 'tnumber of reactions is greater than mxkwr:',maxent
      STOP "in rdvocmchdb"
    ENDIF

! read reactant(s) and the number of products to be read
    READ(line(1:n1-1),*,IOSTAT=ierr) tpkwrct(nkwdat)
    IF (ierr/=0) THEN
      WRITE(6,*) '--error--, while reading file: ',TRIM(filename)
      WRITE(6,*) 'while reading reactant at line number:',ilin
      WRITE(6,*) TRIM(line)
      STOP "in rdvocmchdb" 
    ENDIF
    READ(line(n1+1:ncpbeg-1),*,IOSTAT=ierr) ndat
    IF (ierr/=0) THEN
      WRITE(6,*) '--error--, while reading file: ',TRIM(filename)
      WRITE(6,*) 'while reading # of data at line number:',ilin
      WRITE(6,*) TRIM(line)
      STOP "in rdvocmchdb" 
    ENDIF
    IF (ndat > maxprd) THEN
      WRITE(6,*) '--error--, while reading file: ',TRIM(filename)
      WRITE(6,*) 'at line number: ',ilin
      WRITE(6,*) ' number of products is greater than maxprd:', maxprd
      STOP "in rdvocmchdb"
    ENDIF
    
! read the comment code
    DO i=1,maxcom
      ncpend=INDEX(line(ncpbeg+1:),';')
      IF (ncpend==0) THEN
        IF (line(ncpbeg+1:)/=" ") THEN
          READ(line(ncpbeg+1:),*,IOSTAT=ierr) tpkwcom(nkwdat,i)
          IF (ierr/=0) THEN
            WRITE(6,*) '--error--, while reading file: ',TRIM(filename)
            WRITE(6,*) 'while reading reference code at line:',ilin
            WRITE(6,*) TRIM(line)
            STOP "in rdvocmchdb" 
          ENDIF
        ENDIF
        EXIT
      ELSE
        ncpend=ncpbeg+ncpend  ! count ncpend from the 1st char in line 
        IF (line(ncpbeg+1:ncpend-1)/=" ") THEN
          READ(line(ncpbeg+1:ncpend-1),*,IOSTAT=ierr) tpkwcom(nkwdat,i)
          IF (ierr/=0) THEN
            WRITE(6,*) '--error--, while reading file: ',TRIM(filename)
            WRITE(6,*) 'while reading reference code at line:',ilin
            WRITE(6,*) TRIM(line)
            STOP "in rdvocmchdb" 
          ENDIF
        ENDIF
        ncpbeg=ncpend
        CYCLE       
      ENDIF 
    ENDDO

    CALL stdchm(tpkwrct(nkwdat))
    tpnkwpd(nkwdat) = ndat
        
! loop over the products
    pdloop: DO j=1,ndat
      ilin=ilin+1
      READ (tfu1,'(a)',IOSTAT=ierr) line
      IF (ierr/=0) THEN
        WRITE(6,*) '--error--, while reading file: ',TRIM(filename)
        WRITE(6,*) 'at line number: ',ilin
        WRITE(6,*) 'data missing ??? '
        STOP "in rdvocmchdb"
      ENDIF
      n1=INDEX(line,':')       
      IF (n1==0) THEN
        WRITE(6,*) '--error--, while reading file: ',TRIM(filename)
        WRITE (6,*) '--error-- in rdvocmchdb, ":" not found at line:', ilin
        STOP "in rdvocmchdb"
      ENDIF

! count the number of coproducts
      ncop=0
      DO k=n1+1,lenlin
        IF (line(k:k)=='+') ncop=ncop+1
      ENDDO
      IF (ncop>maxcoprd) THEN
        WRITE(6,*) '--error--, while reading file: ',TRIM(filename)
        WRITE(6,*) 'at line number: ',ilin
        WRITE(6,*) ' number of coproducts is greater than maxcoprd:', maxcoprd
        STOP "in rdvocmchdb"
      ENDIF

! read the yield and the 'main' product
      READ (line(1:n1-1),*,IOSTAT=ierr) tpkwyld(nkwdat,j)
      IF (ierr/=0) THEN
        WRITE(6,*) '--error--, while reading file: ',TRIM(filename)
        WRITE(6,*) 'while reading yield at line number: ',ilin
        WRITE(6,*) TRIM(line)
        STOP "in rdvocmchdb"
      ENDIF
      
      IF (ncop==0) THEN               ! only 1 product
        READ (line(n1+2:),*,IOSTAT=ierr) tpkwpd(nkwdat,j)
        IF (ierr/=0) THEN
          WRITE(6,*) '--error--, while reading file: ',TRIM(filename)
          WRITE(6,*) 'while reading product at line number: ',ilin
          WRITE(6,*) TRIM(line)
          STOP "in rdvocmchdb"
        ENDIF
      ELSE                             ! read the product and coproducts
        n2=INDEX(line,'+')       
        READ (line(n1+2:n2-1),'(a)',IOSTAT=ierr) tpkwpd(nkwdat,j)  ! product
        IF (ierr/=0) THEN
          WRITE(6,*) '--error--, while reading file: ',TRIM(filename)
          WRITE(6,*) 'while reading product at line number: ',ilin
          WRITE(6,*) TRIM(line)
          STOP "in rdvocmchdb"
        ENDIF
        istart=n2
        DO k=1,ncop
          iend=INDEX(line(istart+1:),'+')
          iend=istart+iend
          IF (iend==istart) iend=lenlin
          READ (line(istart+2:iend-1),*,IOSTAT=ierr) tpkwcopd(nkwdat,j,k)  ! coproducts
          IF (ierr/=0) THEN
            WRITE(6,*) '--error--, while reading file: ',TRIM(filename)
            WRITE(6,*) 'while reading co-products at line number: ',ilin
            WRITE(6,*) TRIM(line)
            STOP "in rdvocmchdb"
          ENDIF
          istart=iend
        ENDDO
      ENDIF

      CALL stdchm(tpkwpd(nkwdat,j))

      IF (ncop>0) THEN
        DO k=1, ncop
          IF (tpkwcopd(nkwdat,j,k)(1:1)==' ') THEN
            WRITE(6,*) '--error--, while reading file: ',TRIM(filename)
            WRITE(6,*) 'expected co-product not found in line #:',ilin
            WRITE(6,*) TRIM(line)
            STOP "in rdvocmchdb"
          ENDIF
        ENDDO
      ENDIF

    ENDDO pdloop

  ENDDO rdloop
  CLOSE(tfu1)

! sort the formula
  kwrct(:)=tpkwrct(:)
  CALL sort_string(kwrct(1:nkwdat)) 

! check for duplicate 
  ierr=0
  DO i=1,nkwdat-1
    IF (kwrct(i)==kwrct(i+1)) THEN
      WRITE(6,*) '--error--, while reading file: ',TRIM(filename)
      WRITE(6,*)  'Following species identified 2 times: ',TRIM(kwrct(i))
      ierr=1
    ENDIF
  ENDDO
  IF (ierr/=0) STOP "in rdvocmchdb" 
  
! sort the tables according to formula
  DO i=1,nkwdat
    ilin=srh5(tpkwrct(i),kwrct,nkwdat)
    IF (ilin <= 0) THEN
      WRITE(6,*) '--error--, while sorting species in: ',TRIM(filename)
      WRITE(6,*) 'species "lost" after sorting the list: ', TRIM(tpkwrct(i))
      STOP "in rdvocmchdb" 
    ENDIF
    nkwpd(ilin)=tpnkwpd(i)   ; kwyld(ilin,:)=tpkwyld(i,:) 
    kwpd(ilin,:)=tpkwpd(i,:) ; kwcopd(ilin,:,:)=tpkwcopd(i,:,:)
    kwcom(ilin,:)=tpkwcom(i,:)
  ENDDO

! check the availability of the codes
  ierr=0
  DO i=1,nkwdat
    DO j=1,maxcom
      IF (kwcom(i,j)(1:1) /= ' ') THEN 
        n1=srh5(kwcom(i,j),code,nreac_info)
        IF (n1 <= 0) THEN
          WRITE(6,*) '--error--, while reading file ',TRIM(filename)
          WRITE(6,*) 'unknow comment: ', kwcom(i,j)
          ierr=1
        ENDIF
      ENDIF
    ENDDO
  ENDDO
  IF (ierr/=0) STOP "in rdvocmchdb" 

END SUBROUTINE rdvocmchdb

! ======================================================================
! PURPOSE: Read mechanism data in the file provided as input. The
! routine is in particular called to read know RO2, RCO3 and RO 
! chemistry.
!
! Data are returned to be next stored the database module:  
! - nkr : # of reaction given as input in the file
! - krct(:,2) : formula of the 2 reactants for reaction i
! - kprd(:,3) : formula the main products (up to 3) 
! - arrh(:,3) : arrhenius para. for the reactions         
! - kcost(:,3): stoi. coef. for the main products
! - kcom(:3)  : comment's code for the reactions
! ======================================================================
SUBROUTINE rdradmchdb(filename,nkr,krct,kprd,arrh,kcost,kcom)
  USE keyparameter, ONLY: tfu1
  USE normchem, ONLY: stdchm
  USE searching, ONLY: srh5           ! to seach in a sorted list
  USE references, ONLY: nreac_info, code  ! available code
  IMPLICIT NONE

  CHARACTER(LEN=*),INTENT(IN) :: filename
  INTEGER,INTENT(OUT)         :: nkr
  CHARACTER(LEN=*),INTENT(OUT):: krct(:,:)
  CHARACTER(LEN=*),INTENT(OUT):: kprd(:,:)
  REAL,INTENT(OUT)            :: arrh(:,:)
  REAL,INTENT(OUT)            :: kcost(:,:)
  CHARACTER(LEN=*),INTENT(OUT):: kcom(:,:)

  INTEGER :: i,j,n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,nlast
  INTEGER :: ierr, loerr, ilin
  INTEGER :: ncpbeg, ncpend
  INTEGER, PARAMETER :: lenlin=300  ! length (max) of the line to be read
  CHARACTER(LEN=lenlin) :: line
  INTEGER :: maxent, maxcom

  nkr  = 0 ; arrh(:,:) = 0. ; kcost(:,:) = 0.
  krct(:,:) = ' ' ; kprd(:,:) = ' ' ; kcom(:,:)=' '
  maxent=SIZE(arrh,1) ; maxcom=SIZE(kcom,2) 
  
  OPEN (UNIT=tfu1,FILE=filename,STATUS='OLD')

! read data
  ilin=0
  rdloop: DO
    ilin=ilin+1
    READ (tfu1,'(a)',IOSTAT=ierr) line
    IF (ierr/=0) THEN
      WRITE(6,*) '--error--, while reading file: ',TRIM(filename)
      WRITE(6,*) 'at line number: ',ilin
      WRITE(6,*) 'keyword "END" missing ?'
      STOP "in rdradmchdb"
    ENDIF     
    IF (line(1:3)=='END') EXIT rdloop
    IF (line(1:1)=='!') CYCLE rdloop

    n1=INDEX(line,'|')
    n2=INDEX(line(n1+1:lenlin),'|')+ n1
    n3=INDEX(line(n2+1:lenlin),'|')+ n2
    n4=INDEX(line(n3+1:lenlin),'|')+ n3
    n5=INDEX(line(n4+1:lenlin),'|')+ n4
    n6=INDEX(line(n5+1:lenlin),'|')+ n5
    n7=INDEX(line(n6+1:lenlin),'|')+ n6
    n8=INDEX(line(n7+1:lenlin),'|')+ n7
    n9=INDEX(line(n8+1:lenlin),'|')+ n8
    nlast=INDEX(line(n9+1:lenlin),'|')
    n10=nlast+n9
    
! check that the line is correctly formatted
    IF (nlast==0) THEN
      WRITE(6,*) '--error--, while reading: ',TRIM(filename)
      WRITE(6,*) ' missing "|" at line: ', ilin
      WRITE(6,*) TRIM(line) 
      STOP "in rdradmchdb"
    ENDIF
    
    nkr= nkr+1
    IF (nkr >= maxent) THEN
      WRITE (6,'(a)') '--error--, while reading ro.dat'
      WRITE (6,'(a)') 'number of reaction greater than mkr'
      STOP "in rdradmchdb"
    ENDIF

! search for ';' (reference code) 
    ncpbeg=INDEX(line,';')       
    IF (ncpbeg==0) THEN
      WRITE(6,*) '--error--, while reading file: ',TRIM(filename)
      WRITE (6,*) 'in rdvocmchdb, ";" not found at line:', ilin
      STOP "in rdradmchdb"
    ENDIF
    
! read data
    loerr=0
    READ (line(1:n1-1),*,IOSTAT=ierr)    krct(nkr,1)
      IF (ierr/=0) loerr=loerr+1
    READ (line(n1+1:n2-1),*,IOSTAT=ierr) krct(nkr,2)
      IF (ierr/=0) loerr=loerr+1
    READ (line(n2+1:n3-1),*,IOSTAT=ierr) kcost(nkr,1)
      IF (ierr/=0) loerr=loerr+1
    READ (line(n3+1:n4-1),*,IOSTAT=ierr) kprd(nkr,1)
      IF (ierr/=0) loerr=loerr+1
    READ (line(n4+1:n5-1),*,IOSTAT=ierr) kcost(nkr,2)
      IF (ierr/=0) loerr=loerr+1
    READ (line(n5+1:n6-1),*,IOSTAT=ierr) kprd(nkr,2)
      IF (ierr/=0) loerr=loerr+1
    READ (line(n6+1:n7-1),*,IOSTAT=ierr) kcost(nkr,3)
      IF (ierr/=0) loerr=loerr+1
    READ (line(n7+1:n8-1),*,IOSTAT=ierr) kprd(nkr,3)
      IF (ierr/=0) loerr=loerr+1
    READ (line(n8+1:n9-1),*,IOSTAT=ierr) kcost(nkr,4)
      IF (ierr/=0) loerr=loerr+1
    READ (line(n9+1:n10-1),*,IOSTAT=ierr) kprd(nkr,4)
      IF (ierr/=0) loerr=loerr+1
    READ (line(n10+1:lenlin),*,IOSTAT=ierr)  (arrh(nkr,j),j=1,3)
      IF (ierr/=0) loerr=loerr+1
    
    IF (loerr/=0) THEN
      WRITE(6,*) '--error--, while reading: ',TRIM(filename)
      WRITE(6,*) 'An error identified while reading line #: ', ilin
      WRITE(6,*) TRIM(line) 
      STOP "in rdradmchdb"
    ENDIF

    IF (krct(nkr,2)(1:2)=='- ') krct(nkr,2)=' '
    IF (kprd(nkr,2)(1:2)=='- ') kprd(nkr,2)=' '
    IF (kprd(nkr,3)(1:2)=='- ') kprd(nkr,3)=' '
    IF (kprd(nkr,4)(1:2)=='- ') kprd(nkr,4)=' '

! check the formula
    CALL stdchm(krct(nkr,1))
    CALL stdchm(krct(nkr,2))
    CALL stdchm(kprd(nkr,1))
    CALL stdchm(kprd(nkr,2))
    CALL stdchm(kprd(nkr,3))
    CALL stdchm(kprd(nkr,4))

! read the comment code
    DO i=1,maxcom
      ncpend=INDEX(line(ncpbeg+1:),';')
      IF (ncpend==0) THEN
        IF (line(ncpbeg+1:)/=" ") THEN
          READ(line(ncpbeg+1:),*,IOSTAT=ierr) kcom(nkr,i)
          IF (ierr/=0) THEN
            WRITE(6,*) '--error--, while reading file: ',TRIM(filename)
            WRITE(6,*) 'while reading reference code at line:',ilin
            WRITE(6,*) TRIM(line)
            STOP "in rdradmchdb" 
          ENDIF
        ENDIF
        EXIT
      ELSE
        ncpend=ncpbeg+ncpend  ! count ncpend from the 1st char in line 
        IF (line(ncpbeg+1:ncpend-1)/=" ") THEN
          READ(line(ncpbeg+1:ncpend-1),*,IOSTAT=ierr) kcom(nkr,i)
          IF (ierr/=0) THEN
            WRITE(6,*) '--error--, while reading file: ',TRIM(filename)
            WRITE(6,*) 'while reading reference code at line:',ilin
            WRITE(6,*) TRIM(line)
            STOP "in rdradmchdb" 
          ENDIF
        ENDIF
        ncpbeg=ncpend
        CYCLE       
      ENDIF 
    ENDDO

  ENDDO rdloop
  CLOSE (tfu1)

! check the availability of the codes
  ierr=0
  DO i=1,nkr
    DO j=1,maxcom
      IF (kcom(i,j)(1:1) /= ' ') THEN 
        n1=srh5(kcom(i,j),code,nreac_info)
        IF (n1 <= 0) THEN
          WRITE(6,*) '--error--, while reading file ',TRIM(filename)
          WRITE(6,*) 'unknown comment: ', kcom(i,j)
          ierr=1
        ENDIF
      ENDIF
    ENDDO
  ENDDO
  IF (ierr/=0) STOP "in rdradmchdb" 
      
END SUBROUTINE rdradmchdb

! ======================================================================
! PURPOSE: (1) read the data in photo.dat file (i.e the photolytic 
! data for the species used as "reference") and (2) read typical J 
! values (at a zenithal angle of 40).
!
! Data are stored in the database module:                    
! - njdat : # of photolytic reaction.    
! - jlabel(:) : ID # of photolytic reaction 
! - jchem(:)  : formula of the species being photolyzed  
! - jprod(:,2): the two main photodissociation fragments
! - coprodj(:): additional inorganic coproduct
! - j_com(:,:): comment's code for the photolysis reactions
! - nj40      : # of J40 data in database
! - jlab40(i) : label for which J4O is provided    
! - j40(i)    : J40 values for various labels      
! ======================================================================
SUBROUTINE rdhvdb()
  USE keyparameter, ONLY: tfu1,dirgecko
  USE references, ONLY: nreac_info, code  ! available code
  USE searching, ONLY: srh5           ! to seach in a sorted list
  USE normchem, ONLY: stdchm
  USE database, ONLY: njdat,jlabel,jchem,jprod,coprodj,nj40,jlab40,j40,j_com
  IMPLICIT NONE

  INTEGER :: i,j,n1,n2,n3,n4,ilin,ierr,loerr
  INTEGER :: idat,jdat
  INTEGER :: maxj, maxj40, maxcom
  INTEGER :: ncpend, ncpbeg
  INTEGER,PARAMETER :: lenlin=300  ! length (max) of the line to be read
  CHARACTER(LEN=lenlin) :: line

  jchem(:)=' '   ; jlabel(:)=0 ; jprod(:,:)=' ' ; coprodj(:)=' ' 
  j_com(:,:)=' ' ; j40(:)=0.   ; jlab40(:)=0
  jdat = 0       ; idat=0 ; 
  maxj=SIZE(jlabel) ; maxj40=SIZE(jlab40) ; maxcom=SIZE(j_com,2)

! read the known reactions and their associated label
! ---------------------------------------------------
  OPEN (UNIT=tfu1,FILE=TRIM(dirgecko)//'DATA/photo.dat',STATUS='OLD')

! read photolysis data
  ilin=0
  rdloop: DO
    ilin=ilin+1
    READ (tfu1,'(a)',IOSTAT=ierr) line
    IF (ierr/=0) THEN
      WRITE(6,*) '--error--, while reading file: photo.dat'
      WRITE(6,*) 'at line number: ',ilin
      WRITE(6,*) 'keyword "END" missing ?'
      STOP "in rdhvdb"
    ENDIF     
    IF (line(1:3)=='END') EXIT rdloop
    IF (line(1:1)=='!') CYCLE rdloop

    n1=INDEX(line,'|')
    n2=INDEX(line(n1+1:lenlin),'|')+ n1
    n3=INDEX(line(n2+1:lenlin),'|')+ n2
    n4=INDEX(line(n3+1:lenlin),'|')

! check that the line is correctly formatted
    IF (n4==0) THEN
      WRITE(6,*) '--error--, while reading photo.dat'
      WRITE(6,*) ' missing "|" at line: ', ilin
      WRITE(6,*) TRIM(line) 
      STOP "in rdhvdb"
    ENDIF

    n4 = n4 + n3
 
! search for ';' (reference code) 
    ncpbeg=INDEX(line,';')       
    IF (ncpbeg==0) THEN
      WRITE(6,*) '--error--, while reading file: phot.dat'
      WRITE (6,*) 'in rdvocmchdb, ";" not found at line:', ilin
      STOP "in rdhvdb"
    ENDIF
   
    jdat = jdat+1
    IF (jdat >= maxj) THEN
      WRITE (6,'(a)') '--error--, while reading photo.dat'
      WRITE (6,'(a)') 'number of reaction greater than mkr'
      STOP
    ENDIF
   
    loerr=0
    READ (line(1:n1-1),*,IOSTAT=ierr)    jchem(jdat)
      IF (ierr/=0) loerr=loerr+1
    READ (line(n1+1:n2-1),*,IOSTAT=ierr) jlabel(jdat)
      IF (ierr/=0) loerr=loerr+1
    READ (line(n2+1:n3-1),*,IOSTAT=ierr) jprod(jdat,1)
      IF (ierr/=0) loerr=loerr+1
    READ (line(n3+1:n4-1),*,IOSTAT=ierr) jprod(jdat,2)
      IF (ierr/=0) loerr=loerr+1
    READ (line(n4+1:lenlin),*,IOSTAT=ierr)  coprodj(jdat)
      IF (ierr/=0) loerr=loerr+1
    IF (loerr/=0) THEN
      WRITE(6,*) '--error--, while reading photo.dat '
      WRITE(6,*) 'An error identified while reading line #: ', ilin
      WRITE(6,*) TRIM(line) 
      STOP "in rdhvdb"
    ENDIF
 
! check the formula (if C>1 only)
    CALL stdchm(jchem(jdat))
    CALL stdchm(jprod(jdat,1))
    CALL stdchm(jprod(jdat,2))

    IF (coprodj(jdat)(1:2)=='- ') coprodj(jdat)=' '

! read the comment code
    DO i=1,maxcom
      ncpend=INDEX(line(ncpbeg+1:),';')
      IF (ncpend==0) THEN
        IF (line(ncpbeg+1:)/=" ") THEN
          READ(line(ncpbeg+1:),*,IOSTAT=ierr) j_com(jdat,i)
          IF (ierr/=0) THEN
            WRITE(6,*) '--error--, while reading file: photo.dat'
            WRITE(6,*) 'while reading reference code at line:',ilin
            WRITE(6,*) TRIM(line)
            STOP "in rdhvdb" 
          ENDIF
        ENDIF
        EXIT
      ELSE
        ncpend=ncpbeg+ncpend  ! count ncpend from the 1st char in line 
        IF (line(ncpbeg+1:ncpend-1)/=" ") THEN
          READ(line(ncpbeg+1:ncpend-1),*,IOSTAT=ierr) j_com(jdat,i)
          IF (ierr/=0) THEN
            WRITE(6,*) '--error--, while reading file: photo.dat'
            WRITE(6,*) 'while reading reference code at line:',ilin
            WRITE(6,*) TRIM(line)
            STOP "in rdhvdb" 
          ENDIF
        ENDIF
        ncpbeg=ncpend
        CYCLE       
      ENDIF 
    ENDDO

  ENDDO rdloop
  CLOSE (tfu1)
  njdat=jdat
  
! check the availability of the codes
  ierr=0
  DO i=1,njdat
    DO j=1,maxcom
      IF (j_com(i,j)(1:1) /= ' ') THEN 
        n1=srh5(j_com(i,j),code,nreac_info)
        IF (n1 <= 0) THEN
          WRITE(6,*) '--error--, while reading file: photo.dat'
          WRITE(6,*) 'unknown comment: ', j_com(i,j)
          ierr=1
        ENDIF
      ENDIF
    ENDDO
  ENDDO
  IF (ierr/=0) STOP "in rdradmchdb" 

      
! Read the photolysis rates for a 40° zenith angle        
! ---------------------------------------------------
  OPEN (UNIT=tfu1,FILE=TRIM(dirgecko)//'DATA/j40.dat',STATUS='OLD')

  ilin=0
  rdloop2: DO
    ilin=ilin+1
    READ (tfu1,'(a)',IOSTAT=ierr) line
    IF (ierr/=0) THEN
      WRITE(6,*) '--error--, while reading file j40.dat' 
      WRITE(6,*) 'at line number: ',ilin
      WRITE(6,*) 'keyword "END" missing ?'
      STOP "in rdhvdb"
    ENDIF     
    IF (line(1:3)=='END') EXIT rdloop2
    IF (line(1:1)=='!') CYCLE rdloop2
 
    idat=idat+1
    IF (idat>=maxj40) THEN
      WRITE (6,'(a)') '--error--,  while reading j40.dat'
      WRITE (6,'(a)') ' number of data is greater than maxj40'
      STOP
    ENDIF
    READ (line,*,IOSTAT=ierr) jlab40(idat), j40(idat)
  ENDDO rdloop2
  CLOSE (tfu1)
  nj40=idat
  
END SUBROUTINE rdhvdb

END MODULE rdmecatool
