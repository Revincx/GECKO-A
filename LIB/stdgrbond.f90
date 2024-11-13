MODULE stdgrbond
  IMPLICIT NONE
  CONTAINS

!=======================================================================
! PURPOSE: Contruction of the group matrix and the carbon-carbon bond   
! matrix of a molecule (chem) given as input. Characters (/, \) used to 
! mark cis/trans C=C bonds in chem are removed in groups. The zebond 
! matrix (optional argument!) handle the cis/trans configuration.  
!                                                                   
! First groups are constructed. A group starts with a carbon "C", an
! aromatic carbon "c", or an ether bond "-O-", and ends with either 
! a next "C" or "c", or a double-bond "=", or an ether bond "-O-", 
! or  blank "NUL". In each group the trailing parentheses are deleted 
! if it is an opening with no closing and vice-versa. Next the carbon 
! skeleton is built with carbon, oxygen (if ether) and remaining 
! parentheses and double-bonds. The bond matrix is then built based on 
! the skeleton. The program checks the valence for each carbon center.            
!=======================================================================
SUBROUTINE grbond(chem,group,bond,dbflg,nring,zebond)
  USE keyparameter, ONLY: mxring, digit
  USE toolbox, ONLY: stoperr
  IMPLICIT NONE

  CHARACTER(LEN=*),INTENT(IN) :: chem       ! chemical formula
  CHARACTER(LEN=*),INTENT(OUT):: group(:)   ! groups at position (carbon) i
  INTEGER,INTENT(OUT)         :: bond(:,:)  ! carbon-carbon bond matrix of chem
  INTEGER,INTENT(OUT)         :: dbflg      ! double-bond flag
  INTEGER,INTENT(OUT)         :: nring      ! # of distinct rings
  INTEGER,OPTIONAL,INTENT(OUT):: zebond(:,:)! cis/trans info in C=C bond

! internal:
  CHARACTER(LEN=LEN(chem)) :: skelet
  INTEGER         :: ig,ng,ik,istop,nb,tp,p,nc
  INTEGER         :: ring(SIZE(bond,1),SIZE(bond,2))   ! size(ring)=size(bond)
  INTEGER         :: i,j,k,l,nok,ii,jj,kk,n,nn
  INTEGER         :: dbbrk                             ! flag where SMILES contains 'C=1'
  INTEGER         :: lengr, numgr
  LOGICAL         :: loeter
  INTEGER         :: dbtype(SIZE(bond,1),SIZE(bond,2)) ! cis/trans Ci=Cj type 
  LOGICAL         :: loze                              !

  CHARACTER(LEN=7),PARAMETER :: progname='grbond '
  CHARACTER(LEN=70)          :: mesg

! ----------
! INITIALIZE
! ----------
  ig    = 0   ;   ik = 0  ;  dbflg  = 0   ;  nring = 0  ;  dbbrk = 0
  istop = 0   ;    l = 0  ;  skelet=' '   ;  loeter=.FALSE.
  group(:) = ' '     ;  bond(:,:) = 0     ;  ring(:,:) = 0
  lengr=LEN(group(1))  ;  numgr=SIZE(group)
  nc=INDEX(chem,' ')-1

! initialize Z/E table and check for cis/trans C=C bonds
  dbtype(:,:)=0  ; loze=.FALSE.
  IF (PRESENT(zebond)) THEN 
    zebond(:,:)=0
    IF (SIZE(dbtype)/=SIZE(zebond)) THEN
      mesg="size of zebond table does not match size of bond"
      CALL stoperr(progname,mesg,chem)
    ENDIF
  ENDIF
  IF      (INDEX(chem,'/')/=0) THEN ; loze=.TRUE.
  ELSE IF (INDEX(chem,'\')/=0) THEN ; loze=.TRUE.
  ENDIF
 
  IF (chem(1:1)/='C' .AND. chem(1:1)/='c' .AND. chem(1:2)/='-O') THEN
    IF (chem(1:4) == '=Cd1') THEN
      dbbrk=1
    ELSE
      mesg="input to grbond must begin with C, c or -O "
      CALL stoperr(progname,mesg,chem)
    ENDIF
  ENDIF

! -------------------------------------------------------------
! CONSTRUCT CARBON-CENTERED GROUPS, SKELETON AND DIAGONAL TERM
! -------------------------------------------------------------
  charloop: DO i=1,nc

! identify group
! --------------

! find where the group starts
    nok=0
    IF (chem(i:i)=='C' .OR. chem(i:i)== 'c'  .OR. chem(i:i+1)=='-O') nok=1
    IF (chem(i:i+1)=='Cl') nok=0
    IF (nok==0) CYCLE charloop  
    ig = ig + 1 ! one more group

! find where the group ends
    loopgr : &
    DO j=i+1,i+lengr
      nok=0
      IF (chem(j:j)=='C' .OR. chem(j:j)=='c' .OR. chem(j:j+1)=='-O') nok=1
      IF (chem(j:j+1) == 'Cl')  nok=0
      IF (chem(j:j) == '=')  nok=1
      IF (chem(j:j) == ' ')  nok=1
      IF (nok==1) THEN
        l = j-i   ;  istop = j-1
        EXIT loopgr
      ENDIF
    ENDDO loopgr
    group(ig) = chem(i:istop)  ! write out group

! delete trailing parentheses from group
    delploop: DO
      IF (group(ig)(l:l) == '(') THEN
         group(ig)(l:l) = ' '  ; 
         l = l-1  ; CYCLE delploop
      ELSE IF (group(ig)(l:l) == ')') THEN
         p = 0
         DO j=1,l
           IF (group(ig)(j:j) == '(') p = p+1
           IF (group(ig)(j:j) == ')') p = p-1
         ENDDO
         IF (p < 0) THEN
            group(ig)(l:l) = ' '
            l = l-1 ;  CYCLE delploop
         ENDIF
      ENDIF
      EXIT delploop
    ENDDO delploop

! build skeleton
! --------------
    ik = ik+1

! oxygen skeleton
    IF (group(ig)(1:2)=='-O') THEN
      skelet(ik:ik) = 'O'
      j=3
! ring-joining oxygen (Oxygen is never a multiple center)
      DO n=1,mxring
        IF (group(ig)(3:3)==digit(n)) THEN
          nring = nring+1
          ik=ik+1
          skelet(ik:ik)=digit(n)
          j=4
        ENDIF
      ENDDO
      group(ig)(j:j) = '-'
      loeter=.TRUE.
    ELSE

! carbon (either flavor)
      IF (group(ig)(1:1)=='c' .OR. group(ig)(1:1)=='C') THEN
        skelet(ik:ik) = group(ig)(1:1)
        jj=2
        IF (group(ig)(2:2)=='d') jj=3

! ring-joining carbon
        rjloop: DO j=jj,l
          DO n=1,mxring
            IF (group(ig)(j:j) == digit(n)) THEN
              ik=ik+1
              skelet(ik:ik)=group(ig)(j:j)
              nring = nring+1
              CYCLE rjloop ! try next character (in case multiple center)
            ENDIF
          ENDDO
          EXIT rjloop  ! no more ring-joiners on group  
        ENDDO rjloop
      ENDIF

    ENDIF
        
! put parenthesis if any
    DO j=i+l,istop
      ik = ik+1  ;  skelet(ik:ik) = chem(j:j)
    ENDDO
    IF (chem(istop+1:istop+1) == '=') THEN
      ik = ik+1  ;  skelet(ik:ik) = '='
    ENDIF

! get the number of bonds at each node. Construct diagonal of bond matrix
! (aromaticity is regarded as valence 3 for bookkeeping)
! --------------------------------------------
    nb=0  ;  p=0
    valloop: DO j=2,l
! number of functional group inside parenthesis
      IF (group(ig)(j:j)=='(' .AND. p == 0) nb=nb+1
      IF (group(ig)(j:j)=='(') p = p+1
      IF (group(ig)(j:j)==')') p = p-1

! H-,O-atom and radical attached to a carbon:
      IF (p /= 0) CYCLE valloop 
      IF (group(ig)(j:j)=='H')  nb = nb + 1
      IF (group(ig)(j-1:j-1)=='H' .AND. group(ig)(j:j)=='3')  nb = nb+2
      IF (group(ig)(j-1:j-1)=='H' .AND. group(ig)(j:j)=='2')  nb = nb+1
      IF (group(ig)(j:j)=='O')  nb = nb+2
      IF (group(ig)(j:j)=='.')  nb = nb+1 
    ENDDO valloop

! assign to the diagonal term of the bond matrix
    bond(ig,ig) = nb

  ENDDO charloop

! store number of groups and rings
  ng = ig
  IF (MOD(nring,2)/=0) THEN
    mesg="ERROR in grbond: RING NOT CLOSED "
    CALL stoperr(progname,mesg,chem)
  ENDIF
  nring=nring/2

! ----------------------------------------------------------
! FIND TERMS OF BOND MATRIX (C-C, c-c, c-C, AND C-O-C BONDS)
! ----------------------------------------------------------
  ig = 0
  mainloop: DO i=1,ik-1
    nok=0
    IF (skelet(i:i)=='C'.OR.skelet(i:i)=='c'.OR.skelet(i:i)=='O') nok=1
    IF (nok==0) CYCLE mainloop 
 
    ig=ig+1  ;  k=ig  ;  p=0  ;  tp=0  ;  nb=1
    brchloop: DO j=i+1,ik
      IF (skelet(j:j)=='(') THEN
        tp = 1  ;  p = p+1

      ELSE IF (skelet(j:j)==')') THEN
        tp = 0  ;  p  = p-1
        IF (p<0) CYCLE mainloop

      ELSE IF (skelet(j:j)=='C' .OR. skelet(j:j)=='c') THEN
        k = k+1
        IF (p==0 .AND. tp==0) THEN
          bond(ig,k) = nb  ;  bond(k,ig) = nb
          CYCLE mainloop
        ELSE IF (p==1 .AND. tp==1) THEN
          bond(ig,k) = nb  ;  bond(k,ig) = nb
          nb = 1  ;  tp = 0 
          CYCLE brchloop 
        ELSE IF (p.NE.tp) THEN  ! double bonds up branches
          nb = 1
        ENDIF                       

      ELSE IF (skelet(j:j)=='O') THEN
        k = k + 1
        IF (p==0 .AND. tp==0) THEN
           bond(ig,k) = nb  ;  bond(k,ig) = nb
           CYCLE mainloop
        ELSE IF (p==1 .AND. tp==1) THEN
           bond(ig,k) = nb  ;  bond(k,ig) = nb
           nb = 1  ;  tp = 0 
           CYCLE brchloop 
        ENDIF                       

      ELSE IF (skelet(j-1:j)=='(='.AND. p==1) THEN
        nb = 2  ;  dbflg = 1

      ELSE IF (skelet(j:j)=='=' .AND. p==0) THEN
        nb = 2  ;  dbflg = 1
      ENDIF

    ENDDO brchloop
  ENDDO mainloop

! ----------------------------------------------------------
! FIND BOND MATRIX TERMS FOR RING-JOINING CARBONS
! ----------------------------------------------------------
! i  = group index of current node
! ii = group index of subsequent node 
! j  = skeleton index of current node
! jj = skeleton index of subsequent node
! k  = skeleton index of current ring-join character
! kk = skeleton index of subsequent ring-join character
! ----------------------------------------------------------
  IF (nring> 0) THEN
    i=0  ;  ii=0  ;  j=0  ;  jj=0  ;  k=0  ;  kk=0

    DO j=1,ik-1  ! loop skeleton
      IF (skelet(j:j)=='C'.OR.skelet(j:j)=='c'.OR.skelet(j:j)=='O') THEN  ! node
        k=j  ;  i=i+1  ;  nb = 0

60      CONTINUE
        k=k+1
        DO n=1,nring  ! loop digits
          IF (skelet(k:k)==digit(n)) THEN  ! ring start
            IF (ring(n,1)==1) THEN ! ring already dealt with
              GO TO 60  ! test for multi-ring center 
            ENDIF
            ii = i  
            nb = 1
            DO jj=k+1,ik  ! loop skeleton forwards
              IF (skelet(jj:jj)=='C'.OR.skelet(jj:jj)=='c'.OR.skelet(jj:jj)=='O')THEN  ! node
                ii=ii+1 
                kk = jj+1
                IF(skelet(kk:kk)==digit(n))THEN ! ring joined
                  ring(n,1)=1
                  bond(i,ii)=nb  ;  bond(ii,i)=nb
                  IF(digit(n)=='1'.AND.dbbrk==1)THEN
                    bond(i,ii)=bond(i,ii)+1  ;  bond(ii,i)=bond(ii,i)+1
                  ENDIF
                  GO TO 60  ! test for multi-ring center on 1st node
                ELSE ! test for multi-ring center on 2nd node
                  DO nn=1,nring
                    IF (skelet(kk:kk)==digit(nn)) THEN  ! ring char
                      kk=kk+1  ! RIC 2016 for 3 cycles species
                      IF(skelet(kk:kk)==digit(n))THEN ! ring joined
                        ring(n,1)=1
                        bond(i,ii)=nb  ;  bond(ii,i) = nb
                        GO TO 60  ! test for multi-ring center on 1st node
                      ENDIF  ! ring joined
                    ENDIF  ! ring char
                  ENDDO  ! loop digits for ring char to ignore
                ENDIF  ! ring joined
              ENDIF  ! node
            ENDDO  ! loop skeleton forwards
          ENDIF  ! ring start
        ENDDO  ! loop digits

      ENDIF  ! node
    ENDDO  ! loop skeleton
  ENDIF  ! ring(s) present

! -----------
! FINAL CHECK
! -----------

! check valence = 4 on each carbon, or '3' if aromatic carbon:
  DO i=1,ng
    nb=SUM(bond(i,1:ng)) 

    IF (group(i)(1:1)=='C'.AND.nb/=4 .OR. group(i)(1:1)=='c'.AND.nb/=3) THEN
      WRITE(6,'(a)') '--error-- in grbond, check valence for species:'
      WRITE(6,'(a)') TRIM(chem)
      WRITE(6,'(15X,A11,A50)') 'skeleton = ', TRIM(skelet)
      WRITE(6,'(10X,A16)') 'bond - matrix = '
      DO k=1,ng  ;  WRITE(6,'(30(i3))') (bond(k,l),l=1,ng)  ;  ENDDO
      mesg="check valence  "  ;   CALL stoperr(progname,mesg,chem)
    ENDIF
    bond(i,i) = 0
  ENDDO

! check that no 'Cd' group exists without a Cd=Cd structure
  dbflg=0
  DO i=1,ng
    IF (INDEX(group(i),'Cd').NE.0) THEN
      DO j=1,numgr
        IF (bond(i,j)==2) dbflg=1
      ENDDO
      IF (dbflg/=1) THEN
        WRITE(6,'(a)') '--error-- in routine grbond. >C=C< expected in :'
        WRITE(6,'(a)') TRIM(chem)
        WRITE(6,'(15X,A11,A50)') 'skeleton = ',TRIM(skelet)
        WRITE(6,'(10X,A16)') 'bond - matrix = '
        DO k=1,ng  ;  WRITE(6,'(30(i3))') (bond(k,l),l=1,ng)  ;  ENDDO
        mesg=">C=C< expected "  ;   CALL stoperr(progname,mesg,chem)
      ENDIF              
    ENDIF              
  ENDDO

! check that the -O- function is not a terminal position
  IF (loeter) THEN
    DO i=1,ng
      IF (group(i)(1:2)=='-O') THEN
        nok=0
        DO j=1,ng
          nok=nok+bond(i,j)
        ENDDO
        IF (nok /= 2) THEN
          WRITE(6,'(a)') '--error-- in grbond. Wrong bond # at -O- site'
          WRITE(6,'(a)') 'Expected valence is 2 at that site. '
          WRITE(6,'(a)') 'Species= ',TRIM(chem)
          WRITE(6,'(15X,A11,A50)') 'skeleton = ',TRIM(skelet)
          WRITE(6,'(10X,A16)') 'bond - matrix = '
          DO k=1,ng  ;  WRITE(6,'(30(i3))') (bond(k,l),l=1,ng)  ;  ENDDO
          mesg="Wrong bond # "  ;   CALL stoperr(progname,mesg,chem)
        ENDIF
      ENDIF
    ENDDO   

! C-O-C bonds turned to 3 (instead of 1) to be recognized later. 
    DO i=1,ng
      IF (group(i)(1:2)=='-O') THEN
        DO j=1,ng
          IF (bond(i,j)==1) THEN
            bond(i,j)=3  ;  bond(j,i)=3
          ENDIF
        ENDDO
      ENDIF
    ENDDO  
  ENDIF 

! remove the cis/trans character (/ and \) in groups and get zebond table
  IF (loze) THEN
    IF (dbflg==0) THEN
      mesg="cis-trans (/ or \) found without C=C structure" 
      CALL stoperr(progname,mesg,chem)
    ENDIF
    CALL get_zebond(chem,bond,ng,group,dbtype)
    IF (PRESENT(zebond)) zebond=dbtype 
  ENDIF

END SUBROUTINE grbond

!=======================================================================
! Purpose: 
! 1. Remove the characters (/ or \) in groups used to indicate the cis
!    or trans configuration on a C=C bond.
! 2. Create a table of cis (Z) or trans (E) C=C bond type. The element
!    zebond(i,j) is set to 1 if Cdi and Cdj are connected via a cis 
!    configuration and set 2 if Cdi and Cdj are connected via a trans 
!    configuration.
! The cis/trans configuration is kept only for RCH=CHR structures here.    
!=======================================================================
SUBROUTINE get_zebond(chem,bond,ng,group,zebond)
  USE toolbox, ONLY: stoperr
  IMPLICIT NONE

  CHARACTER(LEN=*),INTENT(IN)   :: chem        ! formula of current species
  INTEGER,INTENT(IN)            :: bond(:,:)   ! bond matrix
  INTEGER,INTENT(IN)            :: ng          ! number of nodes
  INTEGER,INTENT(OUT)           :: zebond(:,:) ! cis (1) or trans (2) [i,j] bond 
  CHARACTER(LEN=*),INTENT(INOUT):: group(:)    ! groups in chem, to be modified with tag

  INTEGER :: i,j,lengr
  INTEGER :: iaslash,islash,jaslash,jslash
  CHARACTER(LEN=12),PARAMETER   :: progname='get_zebond '
  CHARACTER(LEN=70)             :: mesg

  zebond(:,:)=0  ;  lengr=LEN(group(1))
  
! examine the various cases (//, \\, /\, \/) of Cd=Cd bonds
  grloop: DO i=1,ng-1
    IF (group(i)(1:2)=='Cd') THEN
      islash=INDEX(group(i),'/') ; iaslash=INDEX(group(i),'\')       
 
      ! start with >Cdi/= 
      IF (islash /= 0) THEN                                          
        group(i)(islash:lengr-1)=group(i)(islash+1:lengr)
        DO j=i+1,ng
          IF (bond(i,j)==2) THEN
            jslash=INDEX(group(j),'/') ; jaslash=INDEX(group(j),'\')
            IF ( (jslash+jaslash)==0) THEN                           
              mesg="missing cis-trans character on C=C structure" 
              CALL stoperr(progname,mesg,chem)
            ENDIF

            !>Cdi/=Cdj/< (trans bond)
            IF (jslash/=0) THEN                                
              group(j)(jslash:lengr-1)=group(j)(jslash+1:lengr)
              zebond(i,j)=2  ;  zebond(j,i)=2                  
	
            !>Cdi/=Cdj\< (cis bond)                            
            ELSE IF (jaslash/=0) THEN       
              group(j)(jaslash:lengr-1)=group(j)(jaslash+1:lengr)
              zebond(i,j)=1  ;  zebond(j,i)=1

            ELSE
              mesg="cis/trans structure cannot be identified" 
              CALL stoperr(progname,mesg,chem)
            ENDIF
            CYCLE grloop
          ENDIF
        ENDDO
 
      ! start with >Cdi\=
      ELSE IF (iaslash /= 0) THEN                                    
        group(i)(iaslash:lengr-1)=group(i)(iaslash+1:lengr)
        DO j=i+1,ng
          IF (bond(i,j)==2) THEN
            jslash=INDEX(group(j),'/') ; jaslash=INDEX(group(j),'\')
            IF ( (jslash+jaslash)==0) THEN                           
              mesg="missing cis-trans character on C=C structure" 
              CALL stoperr(progname,mesg,chem)
            ENDIF
            
            !>Cdi\=Cdj/< (cis case)
            IF (jslash/=0) THEN                                      
              group(j)(jslash:lengr-1)=group(j)(jslash+1:lengr)
              zebond(i,j)=1  ;  zebond(j,i)=1
            
            !>Cdi\=Cdj\< (trans case)  
            ELSE IF (jaslash/=0) THEN                                
              group(j)(jaslash:lengr-1)=group(j)(jaslash+1:lengr)
              zebond(i,j)=2  ;  zebond(j,i)=2

            ELSE
              mesg="cis/trans structure cannot be identified" 
              CALL stoperr(progname,mesg,chem)
            ENDIF
            CYCLE grloop
          ENDIF
        ENDDO
      ENDIF
 
    ENDIF 
  ENDDO grloop
 
! cis/trans only considered for RCH=CH-R. Kill other cases
  DO i=1,ng-1
    DO j=1,ng
      IF (zebond(i,j)/=0) THEN
        IF (group(i)=='CdH2')       THEN ; zebond(i,j)=0 ; zebond(j,i)=0 ; ENDIF
        IF (group(j)=='CdH2')       THEN ; zebond(i,j)=0 ; zebond(j,i)=0 ; ENDIF
        IF (INDEX(group(i),'H')==0) THEN ; zebond(i,j)=0 ; zebond(j,i)=0 ; ENDIF
        IF (INDEX(group(j),'H')==0) THEN ; zebond(i,j)=0 ; zebond(j,i)=0 ; ENDIF
      ENDIF
    ENDDO
  ENDDO

END SUBROUTINE get_zebond


END MODULE stdgrbond
