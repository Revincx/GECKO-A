MODULE namingtool
IMPLICIT NONE
CONTAINS
! SUBROUTINE naming(chem,iloc,name,fgrp)
! SUBROUTINE codefg(chem,fgrp,lfg)
!=======================================================================
! PURPOSE :  Creates a 6-character name for species given as input (chem). 
! The routine returns the index after which the "new" name must be 
! included in table of the already names (namlst). The codes of to the 
! various functional groups are returned as well, for dictionary update.                                                
!                                                                   
! The purpose is to create a 6-character name for every new initial or 
! product species, which is not already in the dictionary. According to 
! the functional groups in the species the naming system is:
!   1st char   -   most important functional group      
!   2nd char   -   second important functional group    
!   3rd char   -   number of carbon atoms (modulo 10)   
!   4-6th char -   separate species with identical names
!=======================================================================
SUBROUTINE naming(chem,iloc,name,fgrp)
  USE keyparameter, ONLY: mxnode
  USE database, ONLY: nfn,namfn,chemfn
  USE atomtool, ONLY: cnum
  USE searching, ONLY: srh5
  USE dictstackdb, ONLY: nrec,namlst
  USE toolbox, ONLY: stoperr
  IMPLICIT NONE

  CHARACTER(LEN=*),INTENT(IN)  :: chem      ! formula for which name must be given
  INTEGER,INTENT(OUT)          :: iloc      ! index after which name must be included in namlst
  CHARACTER(LEN=*),INTENT(OUT) :: name      ! name for chem
  CHARACTER(LEN=*),INTENT(OUT) :: fgrp      ! functional group code string 

  INTEGER :: lfg,il,ic
  INTEGER :: i,j,k,l
  INTEGER :: ii,jj,kk
  LOGICAL :: lockfn

  CHARACTER(LEN=7),PARAMETER :: progname='naming '
  CHARACTER(LEN=70)          :: mesg

  INTEGER,PARAMETER :: nord=20
  CHARACTER(LEN=nord),PARAMETER :: ord='1234XTEPAGHNUDVKORCZ'
  
  INTEGER,PARAMETER :: nalfa=61
  CHARACTER(LEN=nalfa),PARAMETER ::  &
   alfa='0123456789ABCDEFGHIJKLMNOPQRSTUWXYZabcdefghijklmnopqrstuvwxyz'

  ii=0 ; jj=0 ; kk=0

! check number of carbons
  ic=cnum(chem)
  IF ( (ic<2) .OR. (ic>mxnode) ) THEN
    mesg="Illegal species:"
    CALL stoperr(progname,mesg,chem)
  ENDIF

! get functional groups
  CALL codefg(chem,fgrp,lfg)
      
! ----------------------------------------
! CHECK IF THE SPECIES HAS A FIXED NAME
! ----------------------------------------

! if the species has an imposed name, find the position at which the 
! species must be inserted in the namlst table and return
  IF (nfn>0) THEN
    iloc=srh5(chem,chemfn,nfn)
    IF (iloc>0) THEN
      name=namfn(iloc)
      iloc = srh5(name,namlst,nrec)
      IF (iloc>0) THEN
        WRITE(6,*) '--error-- in naming. Species: ',TRIM(chem)
        WRITE(6,*) 'is provided with a fixed name, which was'
        WRITE(6,*) 'already used. Name = ',name
        STOP "in naming"
      ENDIF
      iloc= - iloc
      RETURN
    ENDIF 
  ENDIF 

! ---------------------------------------------------
! WRITE NAME
! ---------------------------------------------------

! first 2 characters
! ==================
! first two characters = dominant groups in fgrp (reverse order in 'ord')
  il = 0
  c2loop: DO i=1,nord
    DO j=1,lfg
      IF (fgrp(j:j)==ord(i:i)) THEN
        il = il + 1
        name(il:il) = fgrp(j:j)  
        IF (il==2) EXIT c2loop
      ENDIF
    ENDDO
  ENDDO c2loop
  
  IF (il==0) THEN
    mesg=" Cannot find 1st character"
    CALL stoperr(progname,mesg,chem)
  ENDIF

  IF( il==1) name(2:2) = '0'

! 3rd character
! ================
! write number of C's into 3rd position (only modulo 10):
  j = MOD(ic,10)  ;  WRITE(name(3:3),'(I1)') j

! character 4 to 6, counter
! ==========================
! search through the list of already used names. Find char for 4th position
  name(4:6) = '000'
  iloc = srh5(name,namlst,nrec)

! if iloc <= 0 then free slot located at iloc = - iloc + 1 and exit.
! However must check that the name does not already exist for another 
! molecule in the fixedname dict.
! The XYZ000 name may also not exist in the used list but the first 3
! letter XYZ may already be used (XYZ000 simply not used since reserved 
! for a given species). So check that 'XYZ'' is known in the used list.
  lockfn=.TRUE.
  IF (iloc<=0)  THEN
    iloc = -iloc
    IF (name(1:3)/=namlst(iloc+1)(1:3)) THEN   ! 'XYZ' not already used in namlst
      fnloop: DO l=1,nfn
        IF (name==namfn(l)) THEN               ! 'XYZ' exist in fixedname dict.
          ii=1 ; jj=1 ; kk=1
          lockfn=.FALSE.                       ! don't use XYZ000 as 1st slot
          EXIT fnloop
        ENDIF
      ENDDO fnloop
      IF (lockfn) RETURN                       ! XYZ000 is the 1st slot
    ENDIF
  ENDIF

  IF (lockfn) THEN
! loop over the species in namlst. search recorded name having the same 
! first 3 character. Find last position occupied.
    yxzloop: DO k=iloc+1,nrec+1
      IF (namlst(k)(1:3)/=name(1:3)) THEN
        iloc=k
        EXIT yxzloop
      ENDIF
    ENDDO yxzloop
    iloc=iloc-1   ! iloc is the last position with 'XYZ---'
  
! find the indexes
    IF ((namlst(iloc)(1:1)==' ') .AND. (name(4:6)=='000')) THEN
      ii=1 ; jj=1 ; kk=1
    ELSE
      ii=INDEX(alfa,namlst(iloc)(4:4))
      jj=INDEX(alfa,namlst(iloc)(5:5))
      kk=INDEX(alfa,namlst(iloc)(6:6))
    ENDIF
    
    IF ( (ii==0).OR.(jj==0).OR.(kk==0) ) THEN
       WRITE(6,*) '--error-- in naming. One of the character is not set. '
       WRITE(6,*) 'name in the list: ',namlst(iloc)
       WRITE(6,*) 'name search: ',name
       WRITE(6,*) 'chem: ',TRIM(chem)
       mesg='One of the character is not set.'
       CALL stoperr(progname,mesg,chem)
       STOP "in naming"
    ENDIF
  ENDIF

! name the species.
  ijkloop: DO
    kk=kk+1
    IF (kk>nalfa) THEN
      kk=2  ;  jj=jj+1
    ENDIF
    IF (jj>nalfa) THEN
      jj=1  ;  ii=ii+1
    ENDIF
    IF (ii>nalfa) THEN
      WRITE(6,*) '--error-- in naming. All the slot are occupied.'
      WRITE(6,*)  'name:',name,' chem: ',TRIM(chem)
      STOP "in naming"
    ENDIF
    name(4:4)=alfa(ii:ii) ; name(5:5)=alfa(jj:jj) ; name(6:6)=alfa(kk:kk)
    
! check that name does not already exist for another molecule in 
! fixedname database. If yes, then take the next in the list
    DO l=1,nfn
      IF (name==namfn(l)) CYCLE ijkloop
    ENDDO
    
    EXIT ijkloop ! if that point is reached - ii, jj, kk are set
  ENDDO ijkloop
  
END SUBROUTINE naming

!=======================================================================
! PURPOSE: returns the code string (frgp) of the species given as input.    
! Each character in fgrp is a functionnal groups borne by the molecule.   
!=======================================================================
SUBROUTINE codefg(chem,fgrp,lfg)
  USE keyparameter, ONLY: mxring
  USE rjtool
  USE atomtool, ONLY: getatoms
  USE toolbox, ONLY: stoperr
  IMPLICIT NONE
      
  CHARACTER(LEN=*),INTENT(IN) :: chem   ! formula for which fgrp is returned
  CHARACTER(LEN=*),INTENT(OUT):: fgrp   ! code string for the groups in chem   
  INTEGER,INTENT(OUT)         :: lfg    ! length of the fgrp

  INTEGER :: nc
  INTEGER :: ic, ih, in, io, ir, is, ibr, ifl, icl
  INTEGER :: i
  INTEGER :: np, na, ng, nf
  INTEGER :: rjs(mxring,2), nring
  LOGICAL :: lohc,locriegee
  CHARACTER(LEN=LEN(chem))   :: tchem
  CHARACTER(LEN=9),PARAMETER :: hcstring='CH1234().'

  CHARACTER(LEN=7),PARAMETER :: progname='codefg '
  CHARACTER(LEN=70)          :: mesg

  tchem = chem       ! chem can be modified due to removal of rings
  lfg   = 1  ;  fgrp  = ' '
  locriegee=.FALSE.

! count number of atoms and radical dots
  CALL getatoms(tchem,ic,ih,in,io,ir,is,ifl,ibr,icl)

! ----------------------------------------
! GET THE RADICAL FUNCTIONAL GROUP IN CHEM
! ----------------------------------------
  IF (ir/=0) THEN

! radical check: Assume one radical group per chemical (except criegee) 
    IF (ir>1) THEN 
      IF (INDEX(tchem,'.(OO.)')/=0)  locriegee=.TRUE.
      IF (INDEX(tchem,'.(ZOO.)')/=0) locriegee=.TRUE.
      IF (INDEX(tchem,'.(EOO.)')/=0) locriegee=.TRUE.
      IF (.NOT. locriegee) THEN
        mesg="Illegal di-radical "
        CALL stoperr(progname,mesg,chem)
      ENDIF
    ENDIF

! assume alkyl, then overwrite if needed. 
    fgrp(1:1) = '0'
    IF      (locriegee)                 THEN ; fgrp(1:1) = '4'
    ELSE IF (INDEX(tchem,'CO(OO.)')/=0) THEN ; fgrp(1:1) = '3'
    ELSE IF (INDEX(tchem,'(OO.)')/=0)   THEN ; fgrp(1:1) = '2'
    ELSE IF (INDEX(tchem,'(O.)')/=0)    THEN ; fgrp(1:1) = '1'
    ENDIF
    
! add point to fgrp and update the pointer lfg
    fgrp(2:2) = '.' ; lfg = 3
  ENDIF

! ----------------------------------------------
! GET THE (NO RADICAL) FUNCTIONAL GROUP IN CHEM
! ----------------------------------------------
       
! add flag for rings
  nring=0
  IF ((INDEX(tchem,'12')/=0).OR.(INDEX(tchem,'C2')/=0).OR.(INDEX(tchem,'-O2')/=0)) THEN
    fgrp(lfg:lfg+1) = 'TT'
    lfg = lfg+2  ;  nring=2
  ELSE IF ((INDEX(tchem,'C1')/=0).OR.(INDEX(tchem,'-O1')/=0)) THEN
   fgrp(lfg:lfg) = 'T'
   lfg = lfg+1  ;  nring=1
  ENDIF

  IF (INDEX(tchem,'c2')/=0) THEN
   fgrp(lfg:lfg+1) = 'RR'
   lfg = lfg+2  ;  nring=2
  ELSE IF (INDEX(tchem,'c1')/=0) THEN
   fgrp(lfg:lfg) = 'R'
   lfg = lfg+1  ;  nring=1
  ENDIF

! remove ring joining character to search functional groups
  IF (nring>0)  CALL rjsrm(nring,tchem,rjs)
       
! here it is possible to have more than one occurrence of the group:
! multiple occurrences are counted once for the following group.

  IF (ifl>0) THEN ; fgrp(lfg:lfg) = 'F' ; lfg = lfg+1 ; ENDIF  ! fluorine
  IF (ibr>0) THEN ; fgrp(lfg:lfg) = 'B' ; lfg = lfg+1 ; ENDIF  ! bromine
  IF (icl>0) THEN ; fgrp(lfg:lfg) = 'L' ; lfg = lfg+1 ; ENDIF  ! chlorine
!  IF (is>0) THEN  : fgrp(lfg:lfg) = 'S' ; lfg = lfg+1 ; ENDIF  ! sulfur
! amine:
!  i = INDEX(tchem,'NH')
!  IF (i>0) THEN ; fgrp(lfg:lfg) = 'M' ; lfg = lfg+1 ; ENDIF
! ketenes:
  i = INDEX(tchem,'CdO')
  IF (i>0) THEN ; fgrp(lfg:lfg) = 'X'  ; lfg = lfg+1 ; ENDIF

! for following functional groups, multiple occurrences are counted 
  nc    = INDEX(tchem,' ') - 1

! peroxy-acyl-nitrates [CO(OONO2)]
  np=0
  DO i=1,nc-8
    IF(tchem(i:i+8)=='CO(OONO2)') THEN
      fgrp(lfg:lfg) = 'P' ; lfg = lfg+1 ; np = np+1
     ENDIF
  ENDDO

! carboxylic acids [CO(OH)]
  na=0
  DO i=1,nc-5
    IF(tchem(i:i+5)=='CO(OH)') THEN
      fgrp(lfg:lfg) = 'A' ; lfg = lfg+1 ; na = na+1
    ENDIF
  ENDDO

! peroxy-acids [CO(OOH)]
  ng=0
  DO i=1,nc-6
    IF(tchem(I:I+6) == 'CO(OOH)') THEN
      fgrp(lfg:lfg) = 'G' ; lfg = lfg+1 ; ng=ng+1
    ENDIF
  ENDDO

! hydroperoxides [(OOH)]
  DO i=3,nc-4
    IF( (tchem(i:i+4)=='(OOH)').AND.(tchem(i-2:i-1)/='CO') ) THEN
       fgrp(lfg:lfg) = 'H' ; lfg = lfg + 1
    ENDIF
  ENDDO

! [(OOOH)]
  DO i=3,nc-5
    IF( (tchem(i:i+5)=='(OOOH)').AND.(tchem(i-2:i-1)/='CO') ) THEN
       fgrp(lfg:lfg) = 'Z' ; lfg = lfg + 1
    ENDIF
  ENDDO

! nitrates [(ONO2)]
  DO i=3,nc-5
    IF ( (tchem(i:i+5)=='(ONO2)').OR.  &
        ((tchem(i:i+6)=='(OONO2)').AND.(tchem(i-2:i-1)/='CO')) ) THEN
      fgrp(lfg:lfg) = 'N' ; lfg = lfg+1
    ENDIF
  ENDDO

! unsaturated 
  DO i=1,nc
    IF(tchem(i:i)=='=') THEN
      fgrp(lfg:lfg) = 'U' ; lfg = lfg+1
    ENDIF
  ENDDO

! aldehydes [CHO]
  DO i=1,nc-2
    IF(tchem(i:i+2)=='CHO') THEN
       fgrp(lfg:lfg) = 'D' ; lfg = lfg+1
    ENDIF
  ENDDO

! nitro compounds
  DO i=1,nc-4
    IF(tchem(i:i+4)=='(NO2)') THEN
      fgrp(lfg:lfg) = 'V' ; lfg = lfg + 1
    ENDIF
  ENDDO

! nitroso
!  DO i=1,nc-2
!     IF( (tchem(i:i+2)=='NO ').OR.(tchem(i:i+2)=='NO)') ) THEN
!        fgrp(lfg:lfg) = 'W' ; lfg = lfg+1
!     ENDIF
!   ENDDO 

! ketones (include RCO(ONO2), RCOBr,...). Remove pan, acid and peracid
  nf=0
  DO i=1,nc-1
    IF (tchem(i:i+1) == 'CO') THEN ; nf = nf+1 ; ENDIF
  ENDDO
  nf = nf - np - na - ng
  IF (nf>0) THEN
    DO i=1,nf
      fgrp(lfg:lfg) = 'K' ; lfg = lfg+1
    ENDDO
  ENDIF

! hydroxy groups [(OH)] (rm carboxylic acids)
  DO i=1,nc-3
    IF( (tchem(i+1:i+4)=='(OH)') .AND. (tchem(i:i+4)/='O(OH)') ) THEN
      fgrp(lfg:lfg) = 'O' ; lfg = lfg+1
    ENDIF
  ENDDO

! ether groups [-O-]
  DO i=1,nc-2
    IF(tchem(i:i+2)=='-O-') THEN
      fgrp(lfg:lfg) = 'E' ; lfg = lfg+1
    ENDIF
  ENDDO

! ---------------------------------------------------
! IF NO FUNCTIONAL GROUP THEN CHECK FOR HYDROCARBONS
! ---------------------------------------------------

! only possible characters are CH1234 - if any other then test failed
  lohc = .TRUE.
  DO i=1,nc
    IF (INDEX(hcstring,tchem(i:i))==0) lohc = .FALSE.
  ENDDO
  IF ( (lohc) .AND. (lfg == 1) ) fgrp(1:1) = 'C'

  IF (fgrp==' ') THEN
    mesg="No letter code found."
    CALL stoperr(progname,mesg,chem)
  ENDIF

  i=SCAN(fgrp,'056789FBLWSM')
  IF (i/=0) THEN
    mesg="Letter code currently not allowed."
    CALL stoperr(progname,mesg,chem)
  ENDIF
END SUBROUTINE codefg

END MODULE namingtool
