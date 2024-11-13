!****************************************************************
! PURPOSE: perform VOC+HV reactions. 
! NO UPDATE OR CLEANING PERFORMED HERE SINCE THE ROUTINE IS EXPECTED 
! TO BE REPLACED BY A NEW ROUTINE WITHIN MAGNIFY.
! BA: simply turned into f90 "free" format file, and removed the 
! include file (use keyparameter instead).
!****************************************************************
SUBROUTINE hvdiss2(idnam,chem,bond,group,nring,brch,cut_off)
      USE keyparameter, ONLY: mxnode,mxlco,mxlfo,mxlgr,mxnr,mxpd,mecu,waru
      USE references, ONLY:mxlcod
      USE keyflag, ONLY: wrtopeinfo
      USE cdtool, ONLY: alkcheck
      USE rjtool
      USE atomtool
      USE reactool
      USE ringtool, ONLY: findring, rejoin
      USE normchem
      USE bensontool, ONLY: heat
      USE database, ONLY: &
         njdat,jlabel,jchem,jprod,coprodj,jlab40,j40
      USE dictstacktool, ONLY: bratio
      USE rxwrttool, ONLY:rxwrit,rxinit
      USE radchktool, ONLY: radchk
      USE fragmenttool, ONLY: fragm
      USE stuff4hvdiss2, ONLY: openr,setchrom2,xcrieg
      USE toolbox, ONLY: stoperr,addrx,add1tonp
      IMPLICIT NONE

! input:
      CHARACTER(LEN=*),INTENT(IN) :: idnam    ! name of current species
      CHARACTER(LEN=*),INTENT(IN) :: chem     ! formula of current species
      CHARACTER(LEN=mxlgr),INTENT(IN) :: group(mxnode)
      INTEGER,INTENT(IN) ::         bond(mxnode,mxnode)
      REAL,INTENT(IN)    ::         brch
      REAL,INTENT(IN)    ::         cut_off

! internal
      CHARACTER(LEN=mxlfo) pchem(mxnr,2),tempkc,tempkc2,pchema(mxnr,2),tchem 
      CHARACTER(LEN=mxlgr) tgroup(mxnode), pold, pnew
      CHARACTER(LEN=mxlco) coprod(mxnr,2,mxnode), tprod(mxnode),temp_coprod(2)
      CHARACTER(LEN=mxlco) tempnam
      INTEGER         tbond(mxnode,mxnode), flag(mxnr)
      INTEGER         np,nr,nc,nch,nca,ich,ii,ig,i,j,k,l,m,n,kk
      INTEGER         loop,icarb,icase,ia,ia0,ib,i5,check,ia2
      INTEGER         alpha(2),icount
      INTEGER         it,hsec
      INTEGER         hgamma(mxnode,3)
      INTEGER         npos, nter, pyr, pyrpos, nring
      INTEGER         jid(mxnr),nba,nbb,carbo(2),subs(2)
      REAL            xcoeff, frct(2), yield(mxnr)
      REAL            dhfrco, dhfp1, dhfp2, dhnet(mxnr)
      REAL            brtio
      CHARACTER(LEN=mxlco) pp(mxnr)
      REAL            ss(mxnr)

      CHARACTER*1     chromtab(mxnode,4),tchromtab(mxnode,4)

      CHARACTER(LEN=mxlco) r(3), p(mxpd)
      REAL            s(mxpd),arrh(3),ar1bis

      REAL            jvcut,phot_tot,br(mxnr)

      INTEGER         idreac, nlabel
      REAL            xlabel,folow(3),fotroe(4)
 
      INTEGER         known_species
      INTEGER         gp_ia,gp_i,gp_j,subs_alk(mxnode) 
      INTEGER         kgp_ia(2),kgp_i(2),kgp_j(2),ksubs_alk(mxnode)     
      INTEGER         ic,ih,in,io,ir,is,if,ix
      INTEGER         nxc,n1,ic2,y,conj,cjpos
      CHARACTER(LEN=mxlco) copchem
      REAL            rdtcopchem 
!      INTEGER         rjg(mri,2),ring(mxnode)
      INTEGER         ring(mxnode)
      INTEGER         rngflg,opflg,dbflg,begrg,endrg
      REAL            wf
      CHARACTER(LEN=mxlfo) rdckprod(mxnode),pchem_del(mxnr,2)
      CHARACTER(LEN=mxlco) rdcktprod(mxnode,mxnode),coprod_del(mxnr,2,mxnode)
      INTEGER         nip,flag_del(mxnr,2)
      REAL            sc(mxnode),sc_del(mxnr,mxnode)
      CHARACTER(LEN=mxlcod) :: acom                             ! comment code 

      CHARACTER(LEN=10)    :: progname='*hvdiss2*'
      CHARACTER(LEN=70)    :: mesg
      
  INTEGER,PARAMETER     :: mxcom=10
  INTEGER               :: nrxref(mxnr)
  CHARACTER(LEN=mxlcod) :: rxref(mxnr,mxcom)
      
      INTEGER :: mca,mnr,mnp,lco,lcf
      
      mca=mxnode ; mnr=mxnr ; mnp=mxpd ; lco=mxlco ; lcf=mxlco+mxlfo
!      IF (wtflag.NE.0) WRITE(*,*) progname,rdct(lco+1:lcf)

!***********************************************************************
!                            INITIALIZE                               *
!***********************************************************************
  nrxref(:)=1 ; rxref(:,:)=' ' ; rxref(:,1)='OLDHVRX'  ! initialize and set flag for old reaction
      wf=0.001  ! weighting factor for the HV reactions
      brtio=brch*wf
      copchem=' '
      rdtcopchem=0. 
      nca=0
      DO i=1,mca
        IF (group(i)(1:1).NE.' ') nca=nca+1
      ENDDO

! IF RINGS EXIST remove ring-join characters from groups
!      IF (nring.gt.0) THEN
!        CALL rjgrm(nring,group,rjg)
!      ENDIF

! BA while making the f90 file : following initiation does not work.
! BA bot attempt to identify the problem since the routine is expected
! BA to be remove.
!!!  tgroup(:)=group(:) ;  tbond(:,:)=bond(:,:)
!!!  hgamma(:,:)=0      ;  chromtab(:,:)=' '  ; tprod(:)= ' '
!!!  jvcut=.1           ;  nr=0
!!!  s(:)=0.            ;   p(:)=' '
!!!  flag(:)=0          ;  pchem(:,:)=' ' 
!!!  pchema(:,:)=' '    ;  dhnet(:)=0.
!!!  yield(:)=1.        ;  flag_del(:,:)=0
!!!  pchem_del(:,:)=' ' ;  coprod(:,:,:)=' '  ; sc_del(:,:)=0
!!!  temp_coprod(:)=' '
!!!  icase=0 ;  opflg=0 ;  dbflg=0 ;  rngflg=0

      DO i=1,mca
        tgroup(i) = group(i)
        tprod(i)= ' '
        DO j=1,mca
          tbond(i,j) = bond(i,j)
        ENDDO
        DO j=1,3
          hgamma(i,j)=0
        ENDDO
        DO j=1,4
          chromtab(i,j)=' '
        ENDDO
      ENDDO

      jvcut = .1
      nr = 0

      DO i=1,mnp
         s(i)=0.
         DO j=1,lco
           p(i)(j:j)=' '
         ENDDO
      ENDDO
      DO i=1,mnr
         flag(i) = 0
         pchem(i,1) = ' '
         pchem(i,2) = ' '
         pchema(i,1) = ' '
         pchema(i,2) = ' '
         dhnet(i) = 0.
         yield(i)=1.
         flag_del(i,1) = 0
         flag_del(i,2) = 0
         pchem_del(i,1) = ' '
         pchem_del(i,2) = ' '
         DO j=1,mca
            coprod(i,1,j) = ' '
            coprod(i,2,j) = ' '
            coprod_del(i,1,j) = ' '
            coprod_del(i,2,j) = ' '
            sc_del(i,j) = 0
         ENDDO
      ENDDO
      temp_coprod(1)=' '
      temp_coprod(2)=' '
      icase=0
      opflg=0
      dbflg=0
      rngflg=0

!***********************************************************************
!  check that the species does not already exist in the input table   *
!***********************************************************************
      known_species=0
      DO i=1,njdat
!        n = lco + index(rdct(lco+1:lcf),' ')
!        tchem=rdct(lco+1:n)
        tchem=chem
        CALL stdchm(tchem)
       
!        IF (rdct(lco+1:n).EQ.jchem(i)) THEN
        IF (tchem.EQ.jchem(i)) THEN
!          CALL addrx(progname,rdct(lco+1:lcf),nr,flag)
          CALL addrx(progname,chem,nr,flag)
!          CALL addreac(nr,progname,rdct(lco+1:lcf),flag)
          Jid(nr)=jlabel(i)
          known_species = 1
          pchema(nr,1)=jprod(i,1)
          pchema(nr,2)=jprod(i,2)

! check the products
          IF (index(pchema(nr,1),'.').ne.0) THEN
            CALL radchk(pchema(nr,1),rdckprod,rdcktprod,nip,sc,nrxref(nr),rxref(nr,:))
             pchem(nr,1) = rdckprod(1)
            IF (nip.NE.1) WRITE(6,*) '2 produits hvdiss2.f 1'
            DO j=1,mca
              coprod(nr,1,j) = rdcktprod(1,j)
            ENDDO
          ELSE 
            pchem(nr,1)= pchema(nr,1)
          ENDIF
          IF (index(pchema(nr,2),'.').ne.0) THEN  
            CALL radchk(pchema(nr,2),rdckprod,rdcktprod,nip,sc,nrxref(nr),rxref(nr,:))
             pchem(nr,2) = rdckprod(1)
            IF (nip.NE.1) WRITE(6,*) '2 produits hvdiss2.f 2'
              DO j=1,mca
                coprod(nr,2,j) = rdcktprod(1,j)
              ENDDO
           ELSE 
             pchem(nr,2)= pchema(nr,2)
           ENDIF
            CALL stdchm(pchem(nr,1))
            CALL stdchm(pchem(nr,2))
        
! add third product if necessary
          IF (coprodj(i).ne.' ') THEN
            DO j=1,mca
             IF (coprod(nr,2,j)(1:2).EQ.'  ') THEN
                 coprod(nr,2,j)=coprodj(i)
                 GOTO 440
              ENDIF
            ENDDO
          ENDIF  
440      CONTINUE
                   
        ENDIF
      ENDDO
      
      IF (known_species.ne.0) GOTO 400
      

! ***********************************************************************
!                          SET THE CHROMOPHORE TABLE                  *
! **********************************************************************! 
      CALL setchrom2(chem,tbond,tgroup,chromtab)

! restore 
        tgroup = group

! copy the chromophore table
        tchromtab=chromtab

! **********************************************************************
!            DO THE REACTIONS FOR THE VARIOUS CHROMOPHORES            *
! **********************************************************************
!  101 ii=4
      DO 100 ii=1,nca
        ich = ii
        IF (chromtab(ich,1).EQ.' ') GOTO 100
! re-initialize
        tprod=' '
! ====================
! PAN like chromophore
! ====================

        IF (chromtab(ich,1).EQ.'p') THEN
          CALL addrx(progname,chem,nr,flag)
!          CALL addreac(nr,progname,rdct(lco+1:lcf),flag)

! For RCO-CO(OONO2), similar photolysis to biacetyl is assumed
! i.e. => RCO + CO + NO2
          DO i=1,mca
          IF ((tbond(i,ich).ne.0).and.(tgroup(i)(1:3).eq.'CO ')) THEN
              tbond(i,ich) = 0
              tbond(ich,i) = 0
              pold ='CO'
              pnew = 'CO.'
              CALL swap(group(i),pold,tgroup(i),pnew)    
              tgroup(ich)=' '
              CALL rebond(tbond,tgroup,tempkc,nring)
              CALL radchk(tempkc,rdckprod,rdcktprod,nip,sc,nrxref(nr),rxref(nr,:))
               pchem(nr,1) = rdckprod(1)
              IF (nip.NE.1) WRITE(6,*) '2 produits hvdiss2.f 3'
              CALL stdchm(pchem(nr,1))
              DO j=1,mca
                coprod(nr,1,j) = rdcktprod(1,j)
              ENDDO

              tempkc = 'CO.(OONO2)  '
              CALL radchk(tempkc,rdckprod,rdcktprod,nip,sc,nrxref(nr),rxref(nr,:))
              pchem(nr,2) = rdckprod(1)
              IF (nip.NE.1) WRITE(6,*) '2 produits hvdiss2.f 4'
              DO j=1,mca
                coprod(nr,2,j) = rdcktprod(1,j)
              ENDDO
              jid(nr)=5400
! reset:
              tgroup(ich) = group(ich)
              tgroup(i) = group(i)
              tbond(i,ich)=bond(i,ich)
              tbond(ich,i)=bond(ich,i)

              chromtab(ich,1)=' '
              chromtab(i,1)=' ' 
              GOTO 100
             ENDIF
          ENDDO    
             
! other case 
! change (OONO2) to (OO.)
          pold = '(OONO2)'
          pnew = '(OO.)'
          CALL swap(group(ich),pold,tgroup(ich),pnew)

! rebuild, check and rename:
          CALL rebond(tbond,tgroup,tempkc,nring)
          CALL radchk(tempkc,rdckprod,rdcktprod,nip,sc,nrxref(nr),rxref(nr,:))
           pchem(nr,1) = rdckprod(1)
          IF (nip.NE.1) WRITE(6,*) '2 produits hvdiss2.f 5'
          CALL stdchm(pchem(nr,1))

! other product are coproduct linked to pchem(nr)
          DO j=1,mca
            coprod(nr,1,j) = rdcktprod(1,j)
          ENDDO
          pchem(nr,2) = 'NO2 '

! reset:
          tgroup(ich) = group(ich)

          jid(nr)=10500
          GOTO 100
        ENDIF
        
! ==================================================
! CO(ONO2) chromophore (no data, use data for PAN) :
! ==================================================

        IF (chromtab(ich,1).EQ.'q') THEN
          CALL addrx(progname,chem,nr,flag)
!          CALL addreac(nr,progname,rdct(lco+1:lcf),flag)

! change (ONO2) to (O.)
          pold = '(ONO2)'
          pnew = '(O.)'
          CALL swap(group(ich),pold,tgroup(ich),pnew)

! rebuild, check and rename:
          CALL rebond(tbond,tgroup,tempkc,nring)
          CALL radchk(tempkc,rdckprod,rdcktprod,nip,sc,nrxref(nr),rxref(nr,:))
          pchem(nr,1) = rdckprod(1)
          IF (nip.NE.1) WRITE(6,*) '2 produits hvdiss2.f 6'
          CALL stdchm(pchem(nr,1))

! other products are coproducts linked to pchem(nr)
          DO j=1,mca
            coprod(nr,1,j) = rdcktprod(1,j)
          ENDDO
          pchem(nr,2) = 'NO2 '

! reset:

          tgroup(ich) = group(ich)

          jid(nr)=10500
          GOTO 100
        ENDIF
        
! ==================================================
! CO(OOH) chromophore (no data, use data for ROOH) :
! ==================================================

        IF (chromtab(ich,1).EQ.'g') THEN
          CALL addrx(progname,chem,nr,flag)
!          CALL addreac(nr,progname,rdct(lco+1:lcf),flag)
          
! RCO-CO(OOH) assumed to photolyse as RCO-CO(OH)
          DO i=1,mca
            IF ((tbond(ich,i).EQ.1).and.(tgroup(i)(1:3).eq.'CO '))THEN
              tbond(i,ich) = 0
              tbond(ich,i) = 0
              pold ='CO'
              pnew = 'CO(OH)'
              CALL swap(group(i),pold,tgroup(i),pnew)    
              tgroup(ich)=' '
              CALL rebond(tbond,tgroup,tempkc,nring)
              pchem(nr,1)=tempkc
              CALL stdchm(pchem(nr,1))
              DO j=1,mca
                coprod(nr,1,j) = tprod(j)
              ENDDO
              pchem(nr,2) = 'CO2  '
              jid(nr)=32100
! reset:
              tgroup(ich) = group(ich)
              tgroup(i) = group(i)
              tbond(i,ich)=bond(i,ich)
              tbond(ich,i)=bond(ich,i)

              chromtab(i,1)=' '
              chromtab(ich,1)=' '
              GOTO 100
             ENDIF
          ENDDO    

! change CO(OOH) to CO(O.)
          tgroup(ich) = 'CO(O.)'
          CALL rebond(tbond,tgroup,tempkc,nring)
          CALL radchk(tempkc,rdckprod,rdcktprod,nip,sc,nrxref(nr),rxref(nr,:))
           pchem(nr,1) = rdckprod(1)
          IF (nip.EQ.2) THEN
            sc_del(nr,1) = sc(1)
            flag_del(nr,1) = 1
            pchem_del(nr,1) = rdckprod(2)
            sc_del(nr,2) = sc(2)
            DO j=1,mca
              coprod_del(nr,1,j) = rdcktprod(2,j)
            ENDDO
          ENDIF
          CALL stdchm(pchem(nr,1))
          DO j=1,mca
            coprod(nr,1,j) = rdcktprod(1,j)
          ENDDO
          pchem(nr,2) = 'HO  '

! reset:
          tgroup(ich) = group(ich)

          jid(nr)=40100
          GOTO 100
        ENDIF
      
! =================
! -CHO CHROMOPHORE:
! =================

        IF (chromtab(ich,1).EQ.'d') THEN
 
! WARNING - WARNING - WARNING - RING cause problems
! MODIFY AS SOON AS POSSIBLE
!          IF (nring.gt.0) GOTO 100
          
! find the carbon in alpha
          DO i=1,mca
            IF (tbond(ich,i).EQ.1) THEN
              ia=i
            ENDIF          
          ENDDO
          ia0=ia

! check for "single species chromophore
! this species must be "hand written"
          IF (tchromtab(ia,1).eq.'a') THEN
            mesg = 'CHROMOPHORE FOR CO(OH)CHO must be given in the input table'
            CALL stoperr(progname,mesg,chem)
          ENDIF

          IF (tchromtab(ia,1).eq.'d') THEN
            mesg = 'CHROMOPHORE FOR CHOCHO must be given in the input table'
            CALL stoperr(progname,mesg,chem)
          ENDIF

          IF (tchromtab(ia,1).eq.'g') THEN
            mesg = 'CHROMOPHORE FOR CO(OOH)CHO must be given in the input table'
            CALL stoperr(progname,mesg,chem)
          ENDIF

!          IF (tchromtab(ia,1).eq.'p') THEN
!            WRITE(6,'(a)') '--error--'
!            WRITE(6,'(a)') 'from MASTER MECHANISM ROUTINE : hvdiss2'
!            WRITE(6,'(a)') 'CHROMOPHORE FOR CO(OONO2)CHO must be'
!            WRITE(6,'(a)') 'given in the input table'
!            WRITE(6,'(a)') rdct(lco+1:lcf)
!            WRITE(waru,*) 'hvdiss2',rdct(lco+1:lcf) !STOP
!          ENDIF

! ------------------------------------
! ALPHA DICARBONYL -CO-CHO CHROMOPHORE - see also below the same from ketone
! ------------------------------------

          IF (tchromtab(ia,1).eq.'k') THEN
                 
! FIRST PATHWAY : BREAK THE RCO-CHO BOND
! --------------------------------------

            CALL addrx(progname,chem,nr,flag)
!            CALL addreac(nr,progname,rdct(lco+1:lcf),flag)
            
! break bond (doesn't open any rings)
            tbond(ia,ich)=0
            tbond(ich,ia)=0

! add radical dots to separating groups:
            nc = INDEX(tgroup(ich),' ')
            tgroup(ich)(nc:nc) = '.'
            nc = INDEX(tgroup(ia),' ')
            tgroup(ia)(nc:nc) = '.'

! fragment and write in correct format:
  
!            IF (wtflag.NE.0) WRITE(*,*) "fragm1"
            CALL fragm(tbond,tgroup,pchem(nr,1),pchem(nr,2))
            CALL stdchm(pchem(nr,1))
            CALL stdchm(pchem(nr,2))

! check radicals, write in standard format:
            CALL radchk(pchem(nr,1),rdckprod,rdcktprod,nip,sc,nrxref(nr),rxref(nr,:))
            tempkc = rdckprod(1)
            IF (nip.NE.1) WRITE(6,*) '2 produits hvdiss2.f 8'
            DO j=1,mca
               coprod(nr,1,j) = rdcktprod(1,j)
            ENDDO
            CALL stdchm(tempkc)
            pchem(nr,1) = tempkc
            
            CALL radchk(pchem(nr,2),rdckprod,rdcktprod,nip,sc,nrxref(nr),rxref(nr,:))
            tempkc = rdckprod(1)
            IF (nip.NE.1) WRITE(6,*) '2 produits hvdiss2.f 9'
            CALL stdchm(tempkc)
            pchem(nr,2) = tempkc
            DO j=1,mca
                coprod(nr,2,j) = rdcktprod(1,j)
            ENDDO
            
            jid(nr)=21400

! reset:

            tgroup(ia)    = group(ia)
            tgroup(ich)   = group(ich)
            tbond(ia,ich) = bond(ia,ich)
            tbond(ich,ia) = bond(ich,ia)
            
! SECOND PATHWAY : ARRANGE TO RCHO+CO
! --------------------------------------

            CALL addrx(progname,chem,nr,flag)
!            CALL addreac(nr,progname,rdct(lco+1:lcf),flag)

! break bond (doesn't open any rings)
            tbond(ia,ich)=0
            tbond(ich,ia)=0

! change CO to CHO, rebuild and rename
            pold = 'CO'
            pnew = 'CHO'
            CALL swap(group(ia),pold,tgroup(ia),pnew)
            tgroup(ich)=' '

! fragment and write in correct format:
            CALL rebond(tbond,tgroup,tempkc,nring)
            pchem(nr,1)=tempkc
            CALL stdchm(pchem(nr,1))
            pchem(nr,2)='CO  '

            jid(nr)=21500

! reset
            tgroup(ia)    = group(ia)
            tgroup(ich)   = group(ich)
            tbond(ia,ich) = bond(ia,ich)
            tbond(ich,ia) = bond(ich,ia)
       
! THIRD PATHWAY : ARRANGE TO RH + 2 CO
! --------------------------------------

! Note : if R is not alkyl, then use pathway 2 above

! find carbon
            DO i=1,mca
              IF ((tbond(ia,i).GT.0) .AND. (i.NE.ich)) THEN
                ib=i
              ENDIF
            ENDDO

            IF (INDEX(tgroup(ib),'(').NE.0) GOTO 100

            CALL addrx(progname,chem,nr,flag)
!            CALL addreac(nr,progname,rdct(lco+1:lcf),flag)

! check if alkyl or "regular" >C=O 
            IF (tgroup(ib)(1:4).EQ.'CH3 ') THEN
              pchem(nr,1)='CH4 '
              coprod(nr,1,1)='CO  '
              pchem(nr,2)='CO  '
              GOTO 25
            ENDIF

            IF (tgroup(ib)(1:4).EQ.'CHO ') THEN
              pchem(nr,1)='CH2O '
              coprod(nr,1,1)='CO  '
              pchem(nr,2)='CO  '
              GOTO 25
            ENDIF

            check=0
            IF (tgroup(ib)(1:3).eq.'CH2') THEN
              pold='CH2'
              pnew='CH3'
              check=1
            ELSE IF (tgroup(ib)(1:3).eq.'CO ') THEN
              pold='CO'
              pnew='CHO'
              check=1
            ELSE IF (tgroup(ib)(1:2).eq.'CH') THEN
              pold='CH'
              pnew='CH2'
              check=1
            ELSE IF (tgroup(ib)(1:2).eq.'C ') THEN
              pold='C'
              pnew='CH'
              check=1
            ELSE IF (tgroup(ib)(1:2).eq.'C(') THEN
              pold='C('
              pnew='CH('
              check=1
            ENDIF            

            IF (check.EQ.1) THEN
              tbond(ia,ich)=0
              tbond(ich,ia)=0
              tbond(ia,ib)=0
              tbond(ib,ia)=0
              tgroup(ich)=' '
              tgroup(ia)=' '
              CALL swap(group(ib),pold,tgroup(ib),pnew)
              CALL rebond(tbond,tgroup,tempkc,nring)
              pchem(nr,1)=tempkc
              CALL stdchm(pchem(nr,1))
              coprod(nr,1,1)='CO  '
              pchem(nr,2)='CO  '
            ELSE  
              tbond(ia,ich)=0
              tbond(ich,ia)=0
              tgroup(ich)=' '
              pold = 'CO'
              pnew = 'CHO'
              CALL swap(group(ia),pold,tgroup(ia),pnew)
              CALL rebond(tbond,tgroup,tempkc,nring)
              pchem(nr,1)=tempkc
              CALL stdchm(pchem(nr,1))
              pchem(nr,2)='CO'
            ENDIF
          
25          CONTINUE

            jid(nr)=21600

! reset
            tgroup = group
            tbond = bond
            chromtab(ich,1)=' '
            chromtab(ia,1)=' '
            GOTO 100
          ENDIF
        
! ----------------------------------------------
! ALPHA CD : C=C-CHO : SEARCH THE VARIOUS CASES
! ----------------------------------------------

! case 1: -C=C-CHO 
! case 2 : -CO-C=C-CHO 
! case 3 : -CO-C=C-C=C-CHO

          IF (tgroup(ia)(1:2).eq.'Cd') THEN
            icase=0
            DO i=1,mca
              IF (tbond(ia,i).eq.2) THEN
                icase=1
                ib=i
                DO j=1,mca
                  IF ( (tbond(i,j).eq.1) .AND. (j.NE.ia) ) THEN
                    IF (tgroup(j)(1:3).EQ.'CO ') THEN
                      icase=2
                      icarb=j
                    ENDIF
                    IF (tgroup(j)(1:3).EQ.'CHO') THEN
                      icase=2
                      icarb=j
                    ENDIF
                    IF (tgroup(j)(1:2).EQ.'Cd') THEN
                      DO l=1,mca
                        IF (tbond(j,l).eq.2) THEN
                          DO k=1,mca
                            IF (tbond(l,k).eq.1) THEN
                              IF (tgroup(k)(1:3).EQ.'CO ') THEN
                                icase=3
                                i5=l
                                icarb=k
                              ENDIF
                              IF (tgroup(k)(1:3).EQ.'CHO') THEN
                                icase=3
                                icarb=k
                                i5=l
                              ENDIF
                            ENDIF
                          ENDDO
                        ENDIF
                      ENDDO
                    ENDIF
                  ENDIF
                ENDDO 
              ENDIF
            ENDDO

! check that a case was found
            IF (icase.eq.0) THEN
              mesg = 'problem 1 for a C=C-CHO structure molecule is :'
              CALL stoperr(progname,mesg,chem)
            ENDIF

! ------------------------------
! >C=C-CHO CHROMOPHORE - CASE 1
! ------------------------------

            IF (icase.eq.1) THEN

! FIRST CHANNEL : => C=C. + HCO.
! ------------------------------

              CALL addrx(progname,chem,nr,flag)
!              CALL addreac(nr,progname,rdct(lco+1:lcf),flag)

              tbond(ich,ia) = 0
              tbond(ia,ich) = 0

              nc = INDEX(tgroup(ich),' ')
              tgroup(ich)(nc:nc) = '.'
              nc = INDEX(tgroup(ia),' ')
              tgroup(ia)(nc:nc) = '.'

! fragment and write in correct format:
!              IF (wtflag.NE.0) WRITE(*,*) "fragm2"
              CALL fragm(tbond,tgroup,pchem(nr,1),pchem(nr,2))
              CALL stdchm(pchem(nr,1))
              CALL stdchm(pchem(nr,2))

! check radicals, write in standard format:
             CALL radchk(pchem(nr,1),rdckprod,rdcktprod,nip,sc,nrxref(nr),rxref(nr,:))
               tempkc = rdckprod(1)
              IF (nip.NE.1) WRITE(6,*) '2 produits hvdiss2.f 10'
              DO j=1,mca
                coprod(nr,1,j) = rdcktprod(1,j)
              ENDDO
              CALL stdchm(tempkc)
              pchem(nr,1) = tempkc
              
             CALL radchk(pchem(nr,2),rdckprod,rdcktprod,nip,sc,nrxref(nr),rxref(nr,:))
              tempkc = rdckprod(1)
              IF (nip.NE.1) WRITE(6,*) '2 produits hvdiss2.f 11'
              CALL stdchm(tempkc)
              pchem(nr,2) = tempkc
              DO j=1,mca
                coprod(nr,2,j) = rdcktprod(1,j)
              ENDDO

              jid(nr)=21100

! reset
              tgroup(ich)  = group(ich)
              tgroup(ia)   = group(ia)
              tbond(ich,ia) = bond(ich,ia)
              tbond(ia,ich) = bond(ia,ich)

! SECOND CHANNEL : => CO + >CH-C..  -> ( +O2 = CRIEGEE)
! -----------------------------------------------------
! Note : current program is not "well" done (need to be revised later).
! The criegge part MUST BE stored as first product, i.e. pchem(*,1)
              IF (nring.GT.0) CALL findring(ich,ia,nca,tbond,rngflg,ring)
              IF ((tgroup(ia)(1:7).NE.'Cd(NO2)').AND.&
                 (tgroup(ib)(1:7).NE.'Cd(NO2)') .AND.&
                 (ring(ia)==0)) THEN
              CALL addrx(progname,chem,nr,flag)
!              CALL addreac(nr,progname,rdct(lco+1:lcf),flag)

              tbond(ich,ia) = 0
              tbond(ia,ich) = 0
              tbond(ib,ia) = 1
              tbond(ia,ib) = 1

              tgroup(ich) = ' '

              pold='Cd'                            
              pnew='C'
              CALL swap(group(ia),pold,tgroup(ia),pnew)
              nc = INDEX(tgroup(ia),' ')
              tgroup(ia)(nc:nc+6) = '.(OO.)*'
!              tgroup(ia)(nc:nc+5) = criegee

              IF (tgroup(ib)(1:4).EQ.'CdH2') THEN
                pold='CdH2'
                pnew='CH3'
              ELSE IF (tgroup(ib)(1:3).EQ.'CdH') THEN
                pold='CdH'
                pnew='CH2'
              ELSE IF (tgroup(ib)(1:2).EQ.'Cd') THEN
                pold='Cd'
                pnew='CH'
              ELSE 
                mesg = 'problem 2 for a C=C-CHO structure molecule is :'
                CALL stoperr(progname,mesg,chem)
              ENDIF
              CALL swap(group(ib),pold,tgroup(ib),pnew)

              CALL rebond(tbond,tgroup,pchem(nr,1),nring)
              CALL stdchm(pchem(nr,1))
              pchem(nr,2)='CO  '

              jid(nr)=21200

! reset
              tgroup = group
              tbond = bond
              ENDIF
! THIRD CHANNEL : => >C=C-CO. + H.
! --------------------------------

              CALL addrx(progname,chem,nr,flag)
!              CALL addreac(nr,progname,rdct(lco+1:lcf),flag)

              pold = 'CHO'
              pnew = 'CO'
              CALL swap(group(ich),pold,tgroup(ich),pnew)

              nc = INDEX(tgroup(ich),' ')
              tgroup(ich)(nc:nc) = '.'
              CALL rebond(tbond,tgroup,tempkc,nring)
              CALL radchk(tempkc,rdckprod,rdcktprod,nip,sc,nrxref(nr),rxref(nr,:))
               pchem(nr,1) = rdckprod(1)
              IF (nip.NE.1) WRITE(6,*) '2 produits hvdiss2.f 12'
              DO j=1,mca
                coprod(nr,1,j) = rdcktprod(1,j)
              ENDDO
              pchem(nr,2)='HO2 '
              
              jid(nr)=21300

! reset
              tgroup = group
              tbond = bond
              
              chromtab(ich,1)=' '
              GOTO 100
            ENDIF
 
! ----------------------------------
! -CO-C=C-CHO CHROMOPHORE -  CASE 2
! ----------------------------------

            IF (icase.eq.2) THEN
!              nc = index(rdct(lco+1:lcf),' ')
             CALL getatoms(chem,ic,ih,in,io,ir,is,if,y,ix)
              IF (nring.GT.0) THEN  ! BUG correction RV 2016
                CALL findring(ib,icarb,nca,tbond,rngflg,ring)
                IF (ring(icarb).GT.0) GOTO 100
              ENDIF

! FIRST CHANNEL : MAKE FURANONE
! -----------------------------
              CALL addrx(progname,chem,nr,flag)
!              CALL addreac(nr,progname,rdct(lco+1:lcf),flag)

! FUR is furanone - species should be "hand written" in dictionary
              IF (tgroup(icarb)(1:3).EQ.'CHO') THEN
                pchem(nr,1)='-O1-COCdH=CdHC1H2'
                jid(nr)=23000
! coprod XC is introduced to keep the total carbon balance
                tempkc='xxlost'
                CALL bratio(tempkc,brtio,tempnam,nrxref(nr),rxref(nr,:))

                nxc = ic - 4
                DO i=1,nxc
                  coprod(nr,1,i) = tempnam
                ENDDO
              
              ELSE IF (tgroup(icarb)(1:3).EQ.'CO ') THEN 
                pchem(nr,1)='C1H2-O-COCdH=Cd1CH3'
                jid(nr)=23200
! coprod XC is introduced to keep the total carbon balance
                tempkc='xxlost'
                CALL bratio(tempkc,brtio,tempnam,nrxref(nr),rxref(nr,:))

                nxc = ic - 5
                DO i=1,nxc
                  coprod(nr,1,i) = tempnam
                ENDDO
               ENDIF


! SECOND CHANNEL : MAKE MALEIC ANHYDRE
! ------------------------------------
              ic2 = 0
              CALL addrx(progname,chem,nr,flag)
!              CALL addreac(nr,progname,rdct(lco+1:lcf),flag)

! MAL is maleic anhydride - species should be "hand written" in dictionary
              pchem(nr,1)='-O1-COCdH=CdHC1O'
              coprod(nr,1,1) = 'HO2'
              
              IF (tgroup(icarb)(1:3).EQ.'CHO') THEN
                jid(nr)= 23100
                pchem(nr,2)='HO2'
! coprod XC is introduced to keep the total carbon balance
                tempkc='xxlost'
                CALL bratio(tempkc,brtio,tempnam,nrxref(nr),rxref(nr,:))

                nxc = ic - 4
                DO i=2,2+(nxc-1)
                  coprod(nr,1,i) = tempnam
                ENDDO

              ELSE IF (tgroup(icarb)(1:3).EQ.'CO ') THEN
                jid(nr)=23300
                DO i=1,mca
                  IF ( (bond(icarb,i).NE.0) .AND. (i.ne.ib) ) THEN
                     tbond(icarb,i) = 0
                     tbond(i,icarb) = 0
                     nc = INDEX(tgroup(i),' ')
                     tgroup(i)(nc:nc) = '.'
                     IF (tgroup(i).EQ.'-O-.') THEN
                       DO j=1,mca
                         IF ((tbond(j,i).EQ.3).AND.(j.NE.icarb)) THEN
                           nc = INDEX(tgroup(j),' ')
                           tgroup(j)(nc:nc+3) = '(O.)'
                           tgroup(i)=' '
                           tbond(j,i) = 0
                           tbond(i,j) = 0
                           DO k=1,mca
                             IF ((tbond(j,k).EQ.3).AND.(k.NE.i).AND.(tgroup(j)=='-O-(O.)')) THEN
                               nc = INDEX(tgroup(k),' ')
                               tgroup(k)(nc:nc+4) = '(OO.)'
                               tgroup(j)=' '
                               tbond(j,k) = 0
                               tbond(k,j) = 0
                             ENDIF
                           ENDDO
                         ENDIF
                       ENDDO
                     ENDIF
! WARNING : RCO3 is used ONLY to recognise the maleic anhydride fragment
                     tgroup(icarb) = 'CO(OO.)'

!                     IF (wtflag.NE.0) WRITE(*,*) "fragm3"
                     CALL fragm(tbond,tgroup,tempkc,pchem(nr,2))
                     CALL stdchm(tempkc)
                     CALL stdchm(pchem(nr,2))
                     IF (INDEX(pchem(nr,2),'CO(OO.)').NE.0) THEN
                       pchem(nr,2)=tempkc 
                     ENDIF
                     
                     n1 = index(pchem(nr,2),' ')
                     DO j=1,n1
                     IF (index(pchem(nr,2)(j:j),'C').NE.0) ic2 = ic2 + 1
                     ENDDO
                    
! coprod XC is introduced to keep the total carbon balance
                    tempkc='xxlost'
                    CALL bratio(tempkc,brtio,tempnam,nrxref(nr),rxref(nr,:))
                    nxc = ic - 4 - ic2
                    DO k=2,2+(nxc-1)
                      coprod(nr,1,k) = tempnam
                    ENDDO
             CALL radchk(pchem(nr,2),rdckprod,rdcktprod,nip,sc,nrxref(nr),rxref(nr,:))
                     
                    tempkc = rdckprod(1)
                    IF (nip.NE.1) WRITE(6,*) '2 produits hvdiss2.f 13'
                    DO j=1,mca
                      coprod(nr,2,j) = rdcktprod(1,j)
                    ENDDO
                    CALL stdchm(tempkc)
                    pchem(nr,2) = tempkc
                  ENDIF
                ENDDO
              ENDIF


! reset
              tgroup = group
              tbond = bond
              chromtab(icarb,1)=' '
              chromtab(ich,1)=' '
              GOTO 100
             ENDIF

! --------------------------------------
! -CO-C=C-C=C-CHO CHROMOPHORE -  CASE 3
! --------------------------------------

            IF (icase.eq.3) THEN
            
! DIALDEHYDE CASE
! ===============
   
              IF (tgroup(icarb)(1:3).EQ.'CHO') THEN

! FIRST CHANNEL : BREAK THE CD-CHO BOND AT FIRST -CHO
! ---------------------------------------------------
                CALL addrx(progname,chem,nr,flag)
!                CALL addreac(nr,progname,rdct(lco+1:lcf),flag)

                tbond(ich,ia) = 0
                tbond(ia,ich) = 0

                nc = INDEX(tgroup(ich),' ')
                tgroup(ich)(nc:nc) = '.'
                nc = INDEX(tgroup(ia),' ')
                tgroup(ia)(nc:nc) = '.'

! fragment and write in correct format:
!                IF (wtflag.NE.0) WRITE(*,*) "fragm4"
                CALL fragm(tbond,tgroup,pchem(nr,1),pchem(nr,2))
                CALL stdchm(pchem(nr,1))
                CALL stdchm(pchem(nr,2))

! check radicals, write in standard format:
             CALL radchk(pchem(nr,1),rdckprod,rdcktprod,nip,sc,nrxref(nr),rxref(nr,:))
                 tempkc = rdckprod(1)
                IF (nip.NE.1) WRITE(6,*) '2 produits hvdiss2.f 14'
                DO j=1,mca
                  coprod(nr,1,j) = rdcktprod(1,j)
                ENDDO
                CALL stdchm(tempkc)
                pchem(nr,1) = tempkc
              
             CALL radchk(pchem(nr,2),rdckprod,rdcktprod,nip,sc,nrxref(nr),rxref(nr,:))
                 tempkc = rdckprod(1)
                IF (nip.NE.1) WRITE(6,*) '2 produits hvdiss2.f 15'
                CALL stdchm(tempkc)
                pchem(nr,2) = tempkc
                DO j=1,mca
                  coprod(nr,2,j) = rdcktprod(1,j)
                ENDDO

                jid(nr)=21900

! reset
                tgroup(ich)  = group(ich)
                tgroup(ia)   = group(ia)
                tbond(ich,ia) = bond(ich,ia)
                tbond(ia,ich) = bond(ia,ich)

! SECOND CHANNEL : BREAK THE CD-CHO BOND AT SECOND -CHO
! -----------------------------------------------------
                CALL addrx(progname,chem,nr,flag)
!                CALL addreac(nr,progname,rdct(lco+1:lcf),flag)

                tbond(icarb,i5) = 0
                tbond(i5,icarb) = 0

                nc = INDEX(tgroup(i5),' ')
                tgroup(i5)(nc:nc) = '.'
                nc = INDEX(tgroup(icarb),' ')
                tgroup(icarb)(nc:nc) = '.'

! fragment and write in correct format:
!                IF (wtflag.NE.0) WRITE(*,*) "fragm5"
                CALL fragm(tbond,tgroup,pchem(nr,1),pchem(nr,2))
                CALL stdchm(pchem(nr,1))
                CALL stdchm(pchem(nr,2))

! check radicals, write in standard format:
             CALL radchk(pchem(nr,1),rdckprod,rdcktprod,nip,sc,nrxref(nr),rxref(nr,:))
                 tempkc = rdckprod(1)
                IF (nip.NE.1) WRITE(6,*) '2 produits hvdiss2.f 16'
                DO j=1,mca
                  coprod(nr,1,j) = rdcktprod(1,j)
                ENDDO
                CALL stdchm(tempkc)
                pchem(nr,1) = tempkc
              
             CALL radchk(pchem(nr,2),rdckprod,rdcktprod,nip,sc,nrxref(nr),rxref(nr,:))
                 tempkc = rdckprod(1)
                IF (nip.NE.1) WRITE(6,*) '2 produits hvdiss2.f 17'
                CALL stdchm(tempkc)
                pchem(nr,2) = tempkc
                DO j=1,mca
                  coprod(nr,2,j) = rdcktprod(1,j)
                ENDDO

                jid(nr)=21900

! reset
                tgroup(icarb)  = group(icarb)
                tgroup(i5)   = group(i5)
                tbond(icarb,i5) = bond(icarb,i5)
                tbond(i5,icarb) = bond(i5,icarb)

! THIRD AND 4TH CHANNEL : => >C=C-CO. + H. (AT EACH SIDE)
! -------------------------------------------------------
                CALL addrx(progname,chem,nr,flag)
!                CALL addreac(nr,progname,rdct(lco+1:lcf),flag)

                pold = 'CHO'
                pnew = 'CO'
                CALL swap(group(ich),pold,tgroup(ich),pnew)

                nc = INDEX(tgroup(ich),' ')
                tgroup(ich)(nc:nc) = '.'
                CALL rebond(tbond,tgroup,tempkc,nring)
                CALL radchk(tempkc,rdckprod,rdcktprod,nip,sc,nrxref(nr),rxref(nr,:))
                 pchem(nr,1) = rdckprod(1)
                IF (nip.NE.1) WRITE(6,*) '2 produits hvdiss2.f 18'
                DO j=1,mca
                  coprod(nr,1,j) = rdcktprod(1,j)
                ENDDO
                pchem(nr,2)='HO2 '
                tgroup(ich) = group(ich)

                jid(nr)=21900

                CALL addrx(progname,chem,nr,flag)
!                CALL addreac(nr,progname,rdct(lco+1:lcf),flag)

                pold = 'CHO'
                pnew = 'CO'
                CALL swap(group(icarb),pold,tgroup(icarb),pnew)

                nc = INDEX(tgroup(icarb),' ')
                tgroup(icarb)(nc:nc) = '.'
                CALL rebond(tbond,tgroup,tempkc,nring)
                CALL radchk(tempkc,rdckprod,rdcktprod,nip,sc,nrxref(nr),rxref(nr,:))
                 pchem(nr,1) = rdckprod(1)
                IF (nip.NE.1) WRITE(6,*) '2 produits hvdiss2.f 19'
                 DO j=1,mca
                  coprod(nr,1,j) = rdcktprod(1,j)
                ENDDO
                pchem(nr,2)='HO2 '

                jid(nr)=21900

! reset
                tgroup = group
                tbond = bond

! KETO-ALDEHYDE
! =============   

              ELSE IF (tgroup(icarb)(1:3).EQ.'CO ') THEN

! FIRST CHANNEL : BREAK THE CD-CHO BOND AT -CHO
! -----------------------------------------------

                CALL addrx(progname,chem,nr,flag)
!                CALL addreac(nr,progname,rdct(lco+1:lcf),flag)

                tbond(ich,ia) = 0
                tbond(ia,ich) = 0

                nc = INDEX(tgroup(ich),' ')
                tgroup(ich)(nc:nc) = '.'
                nc = INDEX(tgroup(ia),' ')
                tgroup(ia)(nc:nc) = '.'

! fragment and write in correct format:
!                IF (wtflag.NE.0) WRITE(*,*) "fragm6"
                CALL fragm(tbond,tgroup,pchem(nr,1),pchem(nr,2))
                CALL stdchm(pchem(nr,1))
                CALL stdchm(pchem(nr,2))

! check radicals, write in standard format:
             CALL radchk(pchem(nr,1),rdckprod,rdcktprod,nip,sc,nrxref(nr),rxref(nr,:))
                 tempkc = rdckprod(1)
                IF (nip.NE.1) WRITE(6,*) '2 produits hvdiss2.f 20'
                DO j=1,mca
                  coprod(nr,1,j) = rdcktprod(1,j)
                ENDDO
                CALL stdchm(tempkc)
                pchem(nr,1) = tempkc
              
             CALL radchk(pchem(nr,2),rdckprod,rdcktprod,nip,sc,nrxref(nr),rxref(nr,:))
                 tempkc = rdckprod(1)
                IF (nip.NE.1) WRITE(6,*) '2 produits hvdiss2.f 21'
                CALL stdchm(tempkc)
                pchem(nr,2) = tempkc
                DO j=1,mca
                  coprod(nr,2,j) = rdcktprod(1,j)
                ENDDO

                jid(nr)=22000

! reset
                tgroup(ich)  = group(ich)
                tgroup(ia)   = group(ia)
                tbond(ich,ia) = bond(ich,ia)
                tbond(ia,ich) = bond(ia,ich)

! SECOND CHANNEL : => >C=C-CO. + H. 
! ---------------------------------

                CALL addrx(progname,chem,nr,flag)
!                CALL addreac(nr,progname,rdct(lco+1:lcf),flag)

                pold = 'CHO'
                pnew = 'CO'
                CALL swap(group(ich),pold,tgroup(ich),pnew)

                nc = INDEX(tgroup(ich),' ')
                tgroup(ich)(nc:nc) = '.'
                CALL rebond(tbond,tgroup,tempkc,nring)
                CALL radchk(tempkc,rdckprod,rdcktprod,nip,sc,nrxref(nr),rxref(nr,:))
                 pchem(nr,1) = rdckprod(1)
                IF (nip.NE.1) WRITE(6,*) '2 produits hvdiss2.f 22'
                DO j=1,mca
                  coprod(nr,1,j) = rdcktprod(1,j)
                ENDDO
                pchem(nr,2)='HO2 '
                tgroup(ich) = group(ich)

                jid(nr)=22000

! reset
                tgroup = group
                tbond = bond

              ENDIF
              chromtab(icarb,1)=' '
              chromtab(ich,1)=' '
              GOTO 100
            ENDIF
          ENDIF  

! if this point is reached then we have a "regular" R-CHO
          !print*,'regular CO'

! ---------------
! "REGULAR" R-CHO
! ---------------

! check if the C holding the -CHO is tertiary or not
          it=0
          ic=0
          DO j=1,mca
            IF (tbond(ia,j).ne.0) THEN
               ic=ic+1
            ENDIF
          ENDDO
          IF (ic.EQ.4) THEN
             it=1
          ELSE
            it=0
          ENDIF

          npos=0
! prepare Norrish type 2 reaction, search for gamma H
! no Norrish II if there is -OH or -OOH in alpha or a sp2 C between CHO
! and the gamma-H, or if sterically hindered by presence of a ring
! (2 or more of ii/alpha/beta/gamma are on ring = one on-ring bond)
! NOTE: Norrish II may be possible for ring species with complex branches
!       We assume only simple branches here, and disallow Norrish as soon
!       as a ring detected in {ii/alpha/beta/gamma}

          IF(nring.GT.0)THEN
! .. bond(ii,ia)
            CALL findring(ii,ia,nca,bond,rngflg,ring)
            IF(rngflg.GT.0) GOTO 11 ! Norrish II disallowed
! .. bond(ia,ib) 
            DO i=1,mca
              IF ((tbond(ia,i).GT.0) .AND. (i.NE.ich)) THEN
                ib=i
              ENDIF
            ENDDO
            CALL findring(ia,ib,nca,bond,rngflg,ring)
            IF(rngflg.GT.0) GOTO 11 ! Norrish II disallowed
! .. bond(ib,ig): some might be valid: correct later
            DO i=1,mca
              IF ((tbond(ia,i).GT.0) .AND. (i.NE.ich)) THEN
                ib=i
                DO j=1,mca
                  IF ((tbond(ib,j).GT.0) .AND. (j.NE.ia)) THEN
                    ig=j
                    CALL findring(ib,ig,nca,bond,rngflg,ring)
                    IF(rngflg.GT.0) GOTO 11 ! Norrish II disallowed
                  ENDIF
                ENDDO
              ENDIF
            ENDDO

          ENDIF

! find the gamma H and store the info in the hgamma table
! hgamma(j,1) <= carbon class (j being the carbon number having gamma H)
! hgamma(j,2) <= carbon in beta 
! hgamma(j,3) <= carbon in alpha 
! count number of -OH, -OOH or -ONO2 in alpha, beta or gamma position
            hsec=0
            DO j=1,mca
              hgamma(j,1) = 0
              hgamma(j,2) = 0
              hgamma(j,3) = 0
              subs_alk(j) = 0
            ENDDO

            IF (((INDEX(group(ia),'(OH)').NE.0).and.&
                  (INDEX(group(ia),'CO(OH)').EQ.0))      .OR.&
            ((INDEX(group(ia),'(OOH)').NE.0).and.&
            (INDEX(group(ia),'CO(OOH)').EQ.0))) goto 11 
            
            gp_ia = 0
            IF (INDEX(group(ia),'(ONO2)').NE.0) gp_ia = gp_ia + 1
            IF (INDEX(group(ia),'(ONO2)(ONO2)').NE.0) gp_ia = gp_ia + 1
            
            
            DO i=1,mca
              gp_i = 0
              IF ( (bond(ia,i).EQ.1) .AND. (i.ne.ich) ) THEN
                IF (group(i)(1:2).eq.'CO') goto 11
                IF (group(i)(1:2).eq.'Cd') goto 11
                IF (INDEX(group(i),'(OH)').NE.0) gp_i = gp_i + 1
                IF (INDEX(group(i),'(ONO2)').NE.0) gp_i = gp_i + 1
                IF (INDEX(group(i),'(OOH)').NE.0) gp_i = gp_i + 1
                IF (INDEX(group(i),'(OH)(OH)').NE.0) gp_i = gp_i + 1
                IF (INDEX(group(i),'(ONO2)(ONO2)').NE.0) gp_i = gp_i + 1
                IF (INDEX(group(i),'(OOH)(OOH)').NE.0) gp_i = gp_i + 1
                
                DO j=1,mca
                gp_j = 0
                  IF ( (bond(i,j).eq.1) .AND. (j.ne.ia) ) THEN
                    IF (group(j)(1:2).eq.'CO') goto 11
                    IF (group(j)(1:2).eq.'Cd') goto 11

! BA : exclude -O- group bonded to gamma node
                    DO kk=1,mca
                      IF (bond(j,kk)/=0) THEN
                        IF (group(kk)=='-O-') GOTO 11
                      ENDIF
                    ENDDO                  

                    IF (INDEX(group(j),'(OH)').NE.0) gp_j = gp_j + 1
                    IF (INDEX(group(j),'(ONO2)').NE.0) gp_j = gp_j + 1
                    IF (INDEX(group(j),'(OOH)').NE.0) gp_j = gp_j + 1
                    IF (INDEX(group(j),'(OH)(OH)').NE.0) gp_j = gp_j + 1
                    IF (INDEX(group(j),'(ONO2)(ONO2)').NE.0) gp_j=gp_j+1
                    IF (INDEX(group(j),'(OOH)(OOH)').NE.0) gp_j=gp_j+1
                    
                    
                    IF ((gp_ia.le.1).and.(gp_i.le.1).and.(gp_j.le.1))&
                              THEN
                        IF (group(j)(1:3).eq.'CH3') THEN
                          hgamma(j,1)=1
                          hgamma(j,2)=i
                          hgamma(j,3)=ia
                        ELSE IF (group(j)(1:3).eq.'CHO') THEN
                          hgamma(j,1)=0
                        ELSE IF (group(j)(1:3).eq.'CH2') THEN
                          hgamma(j,1)=2
                          hgamma(j,2)=i
                          hgamma(j,3)=ia
                          hsec=1
                        ELSE IF (group(j)(1:2).eq.'CH') THEN
                          hgamma(j,1)=3
                          hgamma(j,2)=i
                          hgamma(j,3)=ia
                          hsec=1
                        ENDIF
                        IF ((gp_ia.ge.1).or.(gp_i.ge.1).or.(gp_j.ge.1))&
                        subs_alk(j)=1                
                    ENDIF
                  ENDIF
                ENDDO
              ENDIF
            ENDDO
             
! remove from a CH3 if a CH2 or a CH exists
            IF (hsec.gt.0) THEN
              DO i=1,mca
                IF (hgamma(i,1).eq.1) hgamma(i,1)=0
              ENDDO
            ENDIF

            npos=0
            DO i=1,mca
              IF (hgamma(i,1).gt.0) npos=npos+1
            ENDDO

11     CONTINUE            

! FIRST CHANNEL : BREAK THE C-CHO BOND
! ------------------------------------

          CALL addrx(progname,chem,nr,flag)
!          CALL addreac(nr,progname,rdct(lco+1:lcf),flag)
  
          tbond(ich,ia) = 0
          tbond(ia,ich) = 0

! add radical dots to separating groups:
          nc = INDEX(tgroup(ich),' ')
          tgroup(ich)(nc:nc) = '.'
          nc = INDEX(tgroup(ia),' ')
          tgroup(ia)(nc:nc) = '.'

! fragment and write in correct format:
!          IF (wtflag.NE.0) WRITE(*,*) "fragm7"
          CALL fragm(tbond,tgroup,pchem(nr,1),pchem(nr,2))
          CALL stdchm(pchem(nr,1))
          CALL stdchm(pchem(nr,2))
                     
! check radicals, write in standard format:
          CALL radchk(pchem(nr,1),rdckprod,rdcktprod,nip,sc,nrxref(nr),rxref(nr,:))
           tempkc = rdckprod(1)
          IF (nip.EQ.2) THEN
!JMLT moved line 1 down !
            sc_del(nr,1) = sc(1)
            flag_del(nr,1) = 1
            pchem_del(nr,1) = rdckprod(2)
            CALL stdchm(pchem_del(nr,1))
            sc_del(nr,2) = sc(2)
            DO j=1,mca
              coprod_del(nr,1,j) = rdcktprod(2,j)
            ENDDO
          ENDIF
          DO j=1,mca
            coprod(nr,1,j) = rdcktprod(1,j)
          ENDDO
          CALL stdchm(tempkc)
          pchem(nr,1) = tempkc

          CALL radchk(pchem(nr,2),rdckprod,rdcktprod,nip,sc,nrxref(nr),rxref(nr,:))
           tempkc = rdckprod(1)
          IF (nip.EQ.2) THEN
!JMLT! changed from sc_del(nr,2) = sc(1)
! and moved 1 line lower
            sc_del(nr,1) = sc(1)
            flag_del(nr,2) = 1
            pchem_del(nr,2) = rdckprod(2)
            CALL stdchm(pchem_del(nr,2))
            sc_del(nr,2) = sc(2)
            DO j=1,mca
              coprod_del(nr,2,j) = rdcktprod(2,j)
            ENDDO
          ENDIF
          CALL stdchm(tempkc)
          pchem(nr,2) = tempkc
          DO j=1,mca
            coprod(nr,2,j) = rdcktprod(1,j)
          ENDDO

! atttach the right jlabel
          IF (((index(group(ia),'(OH)').NE.0).and.&
          (index(group(ia),'CO(OH)').EQ.0))   .OR.        &
          ((index(group(ia),'(OOH)').NE.0).and.&
          (index(group(ia),'CO(OOH)').EQ.0))) THEN 
            jid(nr)=21800
          ELSE  
            IF (npos.eq.0) THEN      ! with no gamma H
              IF (ic.EQ.4) THEN
                jid(nr)=20100
              ELSE IF (ic.eq.3) THEN
                jid(nr)=20200
              ELSE
                jid(nr)=20300
              ENDIF
            ELSE 
              IF (it.EQ.1) THEN      ! with gamma H
                jid(nr)=20400
              ELSE
                jid(nr)=20600
              ENDIF
            ENDIF                  
          ENDIF

! reset
          tgroup(ich)  = group(ich)
          tgroup(ia)   = group(ia)
          tbond(ich,ia) = bond(ich,ia)
          tbond(ia,ich) = bond(ia,ich)

! SECOND CHANNEL : Norrish type 2
! -------------------------------
          
          IF (npos.gt.0) THEN 
            xcoeff=1./real(npos)
            DO i=1,mca
              IF (hgamma(i,1).gt.0) THEN
                CALL addrx(progname,chem,nr,flag)
!                CALL addreac(nr,progname,rdct(lco+1:lcf),flag)
                !print*,'Norrish II'

! change gamma position (i.e. holding the H):
                IF (group(i)(1:3).EQ.'CH3') THEN
                  pold = 'CH3'
                  pnew = 'CdH2'
                ELSE IF(group(i)(1:3).EQ.'CH2') THEN
                  pold = 'CH2'
                  pnew = 'CdH'
                ELSE IF(group(i)(1:2).EQ.'CH') THEN
                  pold = 'CH'
                  pnew = 'Cd'
                ELSE      
                  WRITE(6,'(a)') '--error--'
                  WRITE(6,'(a)') 'from ROUTINE:hvdiss2 '
                  WRITE(6,'(a)') 'molecule could not be identified:'
                  WRITE(6,'(a)') chem
                  WRITE(6,'(a)') 'NORRISH2 was treated - first pb'
                  WRITE(waru,*) 'hvdiss2',chem !STOP
                ENDIF
                         
                CALL swap(group(i),pold,tgroup(i),pnew)
            
! change beta double-bond carbons:
                pold = 'C'
                pnew = 'Cd'
                ib=hgamma(i,2)
                CALL swap(group(ib),pold,tgroup(ib),pnew)

! change alpha position
                check=0
                IF (tgroup(ia)(1:3).eq.'CH2') THEN
                  pold='CH2'
                  pnew='CH3'
                  check=1
                ELSE IF (tgroup(ia)(1:3).eq.'CO ') THEN
                  pold='CO'
                  pnew='CHO'
                  check=1
                ELSE IF (tgroup(ia)(1:2).eq.'CH') THEN
                  pold='CH'
                  pnew='CH2'
                  check=1
                ELSE IF (tgroup(ia)(1:2).eq.'C ') THEN
                  pold='C'
                  pnew='CH'
                  check=1
                ELSE IF (tgroup(ia)(1:2).eq.'C(') THEN
                  pold='C('
                  pnew='CH('
                  check=1
                ENDIF            

                IF (check.EQ.1) THEN
                  CALL swap(group(ia),pold,tgroup(ia),pnew)
                ELSE  
                  WRITE(6,'(a)') '--error--'
                  WRITE(6,'(a)') 'from ROUTINE:hvdiss2 '
                  WRITE(6,'(a)') 'molecule could not be identified:'
                  WRITE(6,'(a)') chem
                  WRITE(6,'(a)') 'NORRISH2 was treated - second pb'
                  WRITE(waru,*) 'hvdiss2',chem !STOP
                ENDIF

! new bond matrix:
                tbond(i,ib) = 2
                tbond(ib,i) = 2
                         
                tbond(ia,ib) = 0
                tbond(ib,ia) = 0
 
! fragment:
!                IF (wtflag.NE.0) WRITE(*,*) "fragm8"
                CALL fragm(tbond,tgroup,pchem(nr,1),pchem(nr,2))
 
! find correct names:
                CALL stdchm(pchem(nr,1))
                CALL stdchm(pchem(nr,2))

! treat substitued alkenes                
                IF (subs_alk(i).ne.0) THEN
                  temp_coprod(1)=' '
                  temp_coprod(2)=' '
                  CALL alkcheck(pchem(nr,1),temp_coprod(1),acom)
                  CALL alkcheck(pchem(nr,2),temp_coprod(2),acom)
                ENDIF

! add coprod and check radicals, write in standard format:
                 
                 IF (INDEX(pchem(nr,1),'.').ne.0) THEN
             CALL radchk(pchem(nr,1),rdckprod,rdcktprod,nip,sc,nrxref(nr),rxref(nr,:))
                    tempkc = rdckprod(1)
                   IF (nip.NE.1) WRITE(6,*) '2 produits hvdiss2.f 25'
                   CALL stdchm(tempkc)
                   pchem(nr,1) = tempkc
                   DO j=1,mca
                     coprod(nr,1,j) = rdcktprod(1,j)
                     IF (coprod(nr,1,j).eq.' ') THEN
                       coprod(nr,1,j) = temp_coprod(1)
                       temp_coprod(1)=' '
                     ENDIF
                   ENDDO
                 ELSE 
                   coprod(nr,1,1)=temp_coprod(1)
                 ENDIF
                
                IF (INDEX(pchem(nr,2),'.').ne.0) THEN
             CALL radchk(pchem(nr,2),rdckprod,rdcktprod,nip,sc,nrxref(nr),rxref(nr,:))
                   tempkc = rdckprod(1)
                  IF (nip.NE.1) WRITE(6,*) '2 produits hvdiss2.f 26'
                  CALL stdchm(tempkc)
! next line added during the radchk modification
                  pchem(nr,2) = tempkc
                  DO j=1,mca
                     coprod(nr,2,j) = rdcktprod(1,j)
                     IF (coprod(nr,2,j).eq.' ') THEN
                       coprod(nr,2,j) = temp_coprod(2)
                       temp_coprod(2)=' '
                     ENDIF
                   ENDDO
                 ELSE 
                   coprod(nr,2,1)=temp_coprod(2)
                 ENDIF
                
       
! attach the right jlabel
                IF (it.EQ.1) THEN
                  jid(nr)=20500
                ELSE
                  jid(nr)=20700
                ENDIF
                yield(nr)=xcoeff
                

! reset
                tgroup = group
                tbond = bond
                !print*,'Norrish II products'
                !print*,pchem(nr,1)
                !print*,pchem(nr,2)
               
              ENDIF
            ENDDO
          ENDIF
          chromtab(ich,1)=' '
          GOTO 100          
        
        ENDIF

! =================
! -CO- CHROMOPHORE:
! =================

        IF (chromtab(ich,1).EQ.'k') THEN  ! corresponding endif = '@@@'

! WARNING - WARNING - WARNING - RING cause problems
! MODIFY AS SOON AS POSSIBLE
!          IF (nring.gt.0) GOTO 100

! find the 2 carbons in alpha
          alpha(1)=0
          alpha(2)=0
          ia=0
          DO i=1,mca
            IF (tbond(ich,i).EQ.1) THEN
              ia=ia+1
              alpha(ia)=i
            ENDIF
          ENDDO
          ia0=ia
          IF (ia.GT.2) THEN
            mesg = 'ia (number of C at each side of a -CO-) > 2 for molecule:'
            CALL stoperr(progname,mesg,chem)
          ENDIF

! WARNING - WARNING - WARNING - CHECK IF KETON ON RING
! IF the ketone belong to a ring => exit (this cause problems in
! the terpene chemistry). TO BE MODIFIED AS SOON AS POSSIBLE
          IF (nring.GT.0)THEN
            CALL findring(ich,alpha(1),nca,tbond,rngflg,ring)
            IF (ring(ich).eq.1) THEN
!              write(6,*) 'hv ring -', rdct(7:50)
              GOTO 100
            ENDIF
            CALL findring(ich,alpha(2),nca,tbond,rngflg,ring)
            IF (ring(ich).eq.1) GOTO 100
          ENDIF

! count the number of carbonyl or Cd in alpha 
! or nitrate (for keto-nitrate photolysis)
          icount=0
          DO i=1,2
            IF ((group(alpha(i))(1:2).EQ.'CO').or.  &
             (group(alpha(i))(1:3).EQ.'CHO').or.&
             (INDEX(group(alpha(i)),'(ONO2)').NE.0).or.&
             (group(alpha(i))(1:2).EQ.'Cd')) THEN
              icount=icount+1
! permutation required if icount=1 from second position
              IF ( (i.eq.2) .and. (icount.eq.1) ) THEN
                ia=alpha(1)
                ia0=ia
                alpha(1)=alpha(2)
                alpha(2)=ia
              ENDIF
            ENDIF
          ENDDO
          IF (icount.gt.0) THEN          ! corresponding endif 'wwwwww'

            DO 45 loop=1,icount
              ia=alpha(loop)
              ia0=ia
              ich=ii

! ---------------------
! R-CO-CO-R CHROMOPHORE
! ---------------------

              IF (tchromtab(ia,1).eq.'k') THEN

! BREAK THE RCO-COR BOND
! ----------------------
                CALL addrx(progname,chem,nr,flag)
!                CALL addreac(nr,progname,rdct(lco+1:lcf),flag)
! Find if bond on ring => ring opening follows
                IF(nring.GT.0)THEN
                  CALL findring(ich,ia,nca,tbond,rngflg,ring)
                ENDIF

! fragment or open and write in correct format:
                IF (rngflg.EQ.0) THEN
! LINEAR CHEM PHOTOLYSIS > RADICALS
! break bond & add radical dots to separating groups:
                  tbond(ia,ich)=0
                  tbond(ich,ia)=0
                  nc = INDEX(tgroup(ich),' ')
                  tgroup(ich)(nc:nc) = '.'
                  nc = INDEX(tgroup(ia),' ')
                  tgroup(ia)(nc:nc) = '.'
!                  IF (wtflag.NE.0) WRITE(*,*) "fragm9"
                  CALL fragm(tbond,tgroup,pchem(nr,1),pchem(nr,2))
                  CALL stdchm(pchem(nr,1))
                  CALL stdchm(pchem(nr,2))
! check radicals, write in standard format:
             CALL radchk(pchem(nr,1),rdckprod,rdcktprod,nip,sc,nrxref(nr),rxref(nr,:))
                   tempkc = rdckprod(1)
                  IF (nip.NE.1) WRITE(6,*) '2 produits hvdiss2.f 27'
                  CALL stdchm(tempkc)
                  pchem(nr,1) = tempkc
                  DO k=1,mca
                    coprod(nr,1,k) = rdcktprod(1,k)
                  ENDDO
             CALL radchk(pchem(nr,2),rdckprod,rdcktprod,nip,sc,nrxref(nr),rxref(nr,:))
                  tempkc2 = rdckprod(1)
                  IF (nip.NE.1) WRITE(6,*) '2 produits hvdiss2.f 28'
                  CALL stdchm(tempkc2)
                  pchem(nr,2) = tempkc2
                  DO k=1,mca
                    coprod(nr,2,k) = rdcktprod(1,k)
                  ENDDO
                ELSE
! CHEM AFTER RING-OPENING > REARRANGE DIRECT TO NON-RADICAL PRODUCTS
! break bond & add radical dots to separating groups:
!                  IF(wtflag.NE.0)&
!                  print*,'ring-opening by breaking CO-CO bond'
                  tbond(ia,ich)=0
                  tbond(ich,ia)=0
                  opflg=1
                  nring=nring-1
                  nc = INDEX(tgroup(ich),' ')
                  tgroup(ich)(nc:nc) = '.'
                  nc = INDEX(tgroup(ia),' ')
                  tgroup(ia)(nc:nc) = '.'
! .. find if ring opens at artificial break point
                  DO k=1,nca
                    IF(ring(k).NE.0) endrg=k
                  ENDDO
                  DO k=nca,1,-1
                    IF(ring(k).NE.0) begrg=k
                  ENDDO
                  IF((ich.EQ.begrg.AND.ia.EQ.endrg)&
                  .OR.(ia.EQ.begrg.AND.ich.EQ.endrg))THEN
                    CONTINUE
! .. if opens elsewhere, rearrange new linear molecule
                  ELSE
                    CALL rejoin(nca,ich,ia,m,n,tbond,tgroup)
                    ich=m
                    ia=n
                  ENDIF
! do instantaneous ring-opening chemistry
                  CALL openr(tbond,tgroup,nring,tempkc,tempkc2,tprod)
                  CALL stdchm(tempkc)
                  CALL stdchm(tempkc2)
                  pchem(nr,1) = tempkc
                  pchem(nr,2) = tempkc2
                  DO k=1,mca
                    coprod(nr,1,k) = tprod(k)
                  ENDDO
                ENDIF

                jid(nr)=31800
! reset
                IF (opflg.eq.1) THEN
                  tgroup = group
                  tbond = bond
                  ia=ia0
                  ich=ii
                  opflg=0
                  rngflg=0
                  nring=nring+1
                ELSE
                  tgroup(ia)    = group(ia)
                  tgroup(ich)   = group(ich)
                  tbond(ia,ich) = bond(ia,ich)
                  tbond(ich,ia) = bond(ich,ia)
                ENDIF
                
                chromtab(ia,1)=' '
                GOTO 45
              ENDIF

! ----------------------
! -CO-CO(OH) CHROMOPHORE
! ----------------------

              IF (tchromtab(ia,1).eq.'a') THEN

! break the RCO-COOH bond -> RCHO + CO2
                CALL addrx(progname,chem,nr,flag)
!                CALL addreac(nr,progname,rdct(lco+1:lcf),flag)

! break bond
                tbond(ia,ich)=0
                tbond(ich,ia)=0

! change CO to CHO, rebuild and rename
                pold='CO'
                pnew='CHO'
                call swap(group(ich),pold,tgroup(ich),pnew)
                tgroup(ia) = ' '
                call rebond(tbond,tgroup,tempkc,nring)
                pchem(nr,1)=tempkc
                CALL stdchm(pchem(nr,1))
                pchem(nr,2)='CO2'

                jid(nr)=32100
                
! reset
                tgroup(ia)    = group(ia)
                tgroup(ich)   = group(ich)
                tbond(ia,ich) = bond(ia,ich)
                tbond(ich,ia) = bond(ich,ia)

                chromtab(ia,1)=' '
                GOTO 45
              ENDIF

! ----------------------
! -CO-CX(ONO2) CHROMOPHORE
! ----------------------
!! keto nitrate photolysis added in the frame of oncem project
!treat keto nitrate in the nitrate block at the end of the routine
!              IF (tchromtab(ia,1).eq.'n') GOTO 100 
              IF (INDEX(group(ia),'(ONO2)').NE.0) GOTO 100
! -----------------------
! -CO-CO(OOH) CHROMOPHORE
! -----------------------

! that reaction is ignored, no data being available
              IF (tchromtab(ia,1).eq.'g') THEN
                GOTO 45
              ENDIF

! --------------------------
! -CO-CO(OONO2) CHROMOPHORE
! --------------------------

! that reaction is ignored, no data being available
              IF (tchromtab(ia,1).eq.'p') THEN
                GOTO 45
              ENDIF

! ---------------------
!  -CO-CHO CHROMOPHORE - see also above same starting from aldehyde
! ---------------------

              IF (tchromtab(ia,1).eq.'d') THEN

! FIRST PATHWAY : BREAK THE RCO-CHO BOND
! --------------------------------------

                CALL addrx(progname,chem,nr,flag)
!                CALL addreac(nr,progname,rdct(lco+1:lcf),flag)
            
! break bond
                tbond(ia,ich)=0
                tbond(ich,ia)=0

! add radical dots to separating groups:
                nc = INDEX(tgroup(ich),' ')
                tgroup(ich)(nc:nc) = '.'
                nc = INDEX(tgroup(ia),' ')
                tgroup(ia)(nc:nc) = '.'

! fragment and write in correct format:
!                IF (wtflag.NE.0) WRITE(*,*) "fragm10"
                CALL fragm(tbond,tgroup,pchem(nr,1),pchem(nr,2))
                CALL stdchm(pchem(nr,1))
                CALL stdchm(pchem(nr,2))

! check radicals, write in standard format:
             CALL radchk(pchem(nr,1),rdckprod,rdcktprod,nip,sc,nrxref(nr),rxref(nr,:))
                 tempkc = rdckprod(1)
                IF (nip.NE.1) WRITE(6,*) '2 produits hvdiss2.f 29'
                DO j=1,mca
                  coprod(nr,1,j) = rdcktprod(1,j)
                ENDDO
                CALL stdchm(tempkc)
                pchem(nr,1) = tempkc
            
             CALL radchk(pchem(nr,2),rdckprod,rdcktprod,nip,sc,nrxref(nr),rxref(nr,:))
                 tempkc = rdckprod(1)
                IF (nip.NE.1) WRITE(6,*) '2 produits hvdiss2.f 30'
                CALL stdchm(tempkc)
                pchem(nr,2) = tempkc
                DO j=1,mca
                  coprod(nr,2,j) = rdcktprod(1,j)
                ENDDO

                jid(nr)=21400

! reset
                tgroup(ia)    = group(ia)
                tgroup(ich)   = group(ich)
                tbond(ia,ich) = bond(ia,ich)
                tbond(ich,ia) = bond(ich,ia)

! SECOND PATHWAY : ARRANGE TO RCHO+CO
! ------------------------------------

               CALL addrx(progname,chem,nr,flag)
!               CALL addreac(nr,progname,rdct(lco+1:lcf),flag)

! break bond
                tbond(ia,ich)=0
                tbond(ich,ia)=0

! change CO to CHO, rebuild and rename
                pold = 'CO'
                pnew = 'CHO'
                CALL swap(group(ich),pold,tgroup(ich),pnew)
                tgroup(ia)=' '

! fragment and write in correct format:
                CALL rebond(tbond,tgroup,tempkc,nring)
                pchem(nr,1)=tempkc
                CALL stdchm(pchem(nr,1))
                pchem(nr,2)='CO  '

                jid(nr)=21500

! reset
                tgroup(ia)    = group(ia)
                tgroup(ich)   = group(ich)
                tbond(ia,ich) = bond(ia,ich)
                tbond(ich,ia) = bond(ich,ia)

! THIRD PATHWAY : arrange to RH + 2 CO 
! ------------------------------------

! if R is not alkyl, then use pathway 2 above


! find carbon
                DO i=1,mca
                  IF ((tbond(ich,i).GT.0) .AND. (i.NE.ia)) THEN
                    ib=i
                  ENDIF
                ENDDO

                IF (INDEX(tgroup(ib),'(').NE.0) GOTO 100
                CALL addrx(progname,chem,nr,flag)
!                CALL addreac(nr,progname,rdct(lco+1:lcf),flag)

! check if alkyl or "regular" >C=O 
                IF (tgroup(ib)(1:4).EQ.'CH3 ') THEN
                  pchem(nr,1)='CH4 '
                  coprod(nr,1,1)='CO  '
                  pchem(nr,2)='CO  '
                  GOTO 26
                ENDIF

                IF (tgroup(ib)(1:4).EQ.'CHO ') THEN
                  pchem(nr,1)='CH2O '
                  coprod(nr,1,1)='CO  '
                  pchem(nr,2)='CO  '
                  GOTO 26
                ENDIF
 
                check=0
                IF (tgroup(ib)(1:3).eq.'CH2') THEN
                  pold='CH2'
                  pnew='CH3'
                  check=1
                ELSE IF (tgroup(ib)(1:3).eq.'CO ') THEN
                  pold='CO'
                  pnew='CHO'
                  check=1
                ELSE IF (tgroup(ib)(1:2).eq.'CH') THEN
                  pold='CH'
                  pnew='CH2'
                  check=1
                ELSE IF (tgroup(ib)(1:2).eq.'C ') THEN
                  pold='C'
                  pnew='CH'
                  check=1
                ELSE IF (tgroup(ib)(1:2).eq.'C(') THEN
                  pold='C('
                  pnew='CH('
                  check=1
                ENDIF            
                IF (check.EQ.1) THEN
                  tbond(ia,ich)=0
                  tbond(ich,ia)=0
                  tbond(ich,ib)=0
                  tbond(ib,ich)=0
                  tgroup(ich)=' '
                  tgroup(ia)=' '
                  CALL swap(group(ib),pold,tgroup(ib),pnew)
                  CALL rebond(tbond,tgroup,tempkc,nring)
                  pchem(nr,1)=tempkc
                  CALL stdchm(pchem(nr,1))
                  coprod(nr,1,1)='CO  '
                  pchem(nr,2)='CO  '
                ELSE  
                  tbond(ia,ich)=0
                  tbond(ich,ia)=0
                  tgroup(ia)=' '
                  pold = 'CO'
                  pnew = 'CHO'
                  CALL swap(group(ich),pold,tgroup(ich),pnew)
                  CALL rebond(tbond,tgroup,tempkc,nring)
                  pchem(nr,1)=tempkc
                  CALL stdchm(pchem(nr,1))
                  pchem(nr,2)='CO  '
                ENDIF

26              CONTINUE

                jid(nr)=21600

! reset
                tgroup = group
                tbond = bond

                chromtab(ia,1)=' '
                GOTO 45
              ENDIF

! ------------------------------------------------------
! ALPHA Cd : C=C-COR  - SEARCH WHICH KIND OF CHROMOPHORE
! ------------------------------------------------------
! search for the various cases
! case 1: -C=C-COR 
! case 2 : -CO-C=C-COR or CHO-C=C-COR
! case 3 : -CO-C=C-C=C-COR or CHO-C=C-C=C-COR
              IF (tgroup(ia)(1:2).eq.'Cd') THEN
                icase=0
                DO i=1,mca
                  IF (tbond(ia,i).eq.2) THEN
                    icase=1
                    ib=i
                    DO j=1,mca
                      IF ( (tbond(i,j).eq.1) .AND. (j.NE.ia) ) THEN
                        IF (tgroup(j)(1:3).EQ.'CO ') THEN
                          icase=2
                          icarb=j
                        ENDIF
                        IF (tgroup(j)(1:3).EQ.'CHO') THEN
                          icase=2
                          icarb=j
                        ENDIF
                        IF (tgroup(j)(1:2).EQ.'Cd') THEN
                          DO l=1,mca
                            IF (tbond(j,l).eq.2) THEN
                              DO k=1,mca
                                IF (tbond(l,k).eq.1) THEN
                                  IF (tgroup(k)(1:3).EQ.'CO ') THEN
                                    icase=3
                                    i5=l
                                    icarb=k
                                  ENDIF
                                  IF (tgroup(k)(1:3).EQ.'CHO') THEN
                                    icase=3
                                    icarb=k
                                    i5=l
                                  ENDIF
                                ENDIF
                              ENDDO
                            ENDIF
                          ENDDO
                        ENDIF
                      ENDIF
                    ENDDO 
                  ENDIF
                ENDDO

! check that a case was found
                IF (icase.eq.0) THEN
                  mesg = 'problem 1a for a C=C-CHO structure molecule is:'
                  CALL stoperr(progname,mesg,chem)
                ENDIF

! --------------------------------
! >C=C-CO-R CHROMOPHORE  -  CASE 1
! CAUTION! may need to allow for ring-breaking  eventually
! --------------------------------

                IF (icase.eq.1) THEN

! FIRST CHANNEL : => C=C. + RCO.
! ------------------------------

                  CALL addrx(progname,chem,nr,flag)
!                  CALL addreac(nr,progname,rdct(lco+1:lcf),flag)

                  tbond(ich,ia) = 0
                  tbond(ia,ich) = 0

                  nc = INDEX(tgroup(ich),' ')
                  tgroup(ich)(nc:nc) = '.'
                  nc = INDEX(tgroup(ia),' ')
                  tgroup(ia)(nc:nc) = '.'

! fragment and write in correct format:
!                  IF (wtflag.NE.0) WRITE(*,*) "fragm11"
                  CALL fragm(tbond,tgroup,pchem(nr,1),pchem(nr,2))
                  CALL stdchm(pchem(nr,1))
                  !CALL stdchm(pchem(nr,2))

! check radicals, write in standard format:
             CALL radchk(pchem(nr,1),rdckprod,rdcktprod,nip,sc,nrxref(nr),rxref(nr,:))
                   tempkc = rdckprod(1)
                  IF (nip.NE.1) WRITE(6,*) '2 produits hvdiss2.f 31'
                  DO j=1,mca
                    coprod(nr,1,j) = rdcktprod(1,j)
                  ENDDO
                  CALL stdchm(tempkc)
                  pchem(nr,1) = tempkc
                
             CALL radchk(pchem(nr,2),rdckprod,rdcktprod,nip,sc,nrxref(nr),rxref(nr,:))
                   tempkc = rdckprod(1)
                  IF (nip.NE.1) WRITE(6,*) '2 produits hvdiss2.f 32'
                  CALL stdchm(tempkc)
                  pchem(nr,2) = tempkc
                  DO j=1,mca
                    coprod(nr,2,j) = rdcktprod(1,j)
                  ENDDO

                  jid(nr)=31700

! reset
                  tgroup(ich)  = group(ich)
                  tgroup(ia)   = group(ia)
                  tbond(ich,ia) = bond(ich,ia)
                  tbond(ia,ich) = bond(ia,ich)

! SECOND CHANNEL : => >C=C-R + CO 
! -------------------------------

! apply only if R is alkyl

! find carbon
                  DO i=1,mca
                    IF ((tbond(ich,i).GT.0) .AND. (i.NE.ia)) THEN
                      ib=i
                    ENDIF
                  ENDDO

! check if alkyl or "regular" >C=O 
                  check=0
                  IF (tgroup(ib)(1:3).EQ.'CHO') check=1
                  IF (tgroup(ib)(1:2).EQ.'Cd') check=1
                  IF (tgroup(ib)(1:2).EQ.'CO') check=1
 
                  IF (check.EQ.0) THEN
                    CALL addrx(progname,chem,nr,flag)
!                    CALL addreac(nr,progname,rdct(lco+1:lcf),flag)
    
                    tbond(ia,ich)=0
                    tbond(ich,ia)=0
                    tbond(ich,ib)=0
                    tbond(ib,ich)=0
                    tgroup(ich)=' '

                    tbond(ia,ib)=1
                    tbond(ib,ia)=1

                    CALL rebond(tbond,tgroup,tempkc,nring)
                    pchem(nr,1)=tempkc
                    CALL stdchm(pchem(nr,1))
                    pchem(nr,2)='CO  '
                  ENDIF
 
                  jid(nr)=31600

! reset
                  tgroup = group
                  tbond = bond

                  GOTO 45
                ENDIF

! ------------------------------------
! -CO-C=C-CO- CHROMOPHORE  -  CASE 2
! ------------------------------------

                IF (icase.eq.2) THEN
                  IF (nring.GT.0) THEN
                    CALL findring(icarb,ia,nca,tbond,rngflg,ring)
                    IF (ring(icarb).EQ.1) GOTO 33
                  ENDIF

!                  nc = index(rdct(lco+1:lcf),' ')
                  CALL getatoms(chem,ic,ih,in,io,ir,is,if&
                       ,y,ix)
                                      
                  DO i=1,mca
                    IF (bond(ia,i).EQ.2) ib =i
                  ENDDO 


! FIRST CHANNEL : MAKE FURANONE
! -----------------------------

! FUR is furanone - species should be "hand written" in dictionary
                  IF (tgroup(icarb)(1:3).EQ.'CHO') THEN
                    CALL addrx(progname,chem,nr,flag)
!                    CALL addreac(nr,progname,rdct(lco+1:lcf),flag)
                    pchem(nr,1)='C1H2-O-COCdH=Cd1CH3'
                    jid(nr)=23200

                    tempkc='xxlost'
                    CALL bratio(tempkc,brtio,tempnam,nrxref(nr),rxref(nr,:))

                    nxc = ic - 5
                    DO i=1,nxc
                      coprod(nr,1,i) = tempnam
                    ENDDO
                  ENDIF

! For di-ketones, break the CO-R bond to form the more substitued R.                  
                  IF (tgroup(icarb)(1:3).EQ.'CO ') THEN
                    nba=4
                    nbb=4
                    carbo(1)= 0
                    carbo(2)= 0
                    subs(1) = 0
                    subs(2) = 0
                    DO i=1,mca
                      IF ((tbond(ich,i).eq.1).AND.(i.ne.ia)) THEN
                        carbo(1)=ich
                        subs(1)=i
                        IF (tgroup(i)(1:2).eq.'CH')  nba=3
                        IF (tgroup(i)(1:3).eq.'CH2') nba=2
                        IF (tgroup(i)(1:3).eq.'CH3') nba=1
                      ELSE IF ((tbond(icarb,i).eq.1).AND.(i.ne.ib)) THEN
                        carbo(2)=icarb
                        subs(2)=i
                        IF (tgroup(i)(1:2).eq.'CH')  nbb=3
                        IF (tgroup(i)(1:3).eq.'CH2') nbb=2
                        IF (tgroup(i)(1:3).eq.'CH3') nbb=1
                      ELSE IF ((tbond(icarb,i).eq.3).AND.(i.ne.ib)) THEN
!***! WARNING in case of ester CO linked to -O-, need update
                        GOTO 33
                      ENDIF
                    ENDDO
                    
                    frct(1)=0.
                    frct(2)=0.
                    IF (nba.gt.nbb) THEN
                      frct(1) = 1.
                      frct(2) = 0.
                    ELSE IF (nba.eq.nbb) THEN
                      frct(1) = 0.5
                      frct(2) = 0.5
                    ELSE IF (nba.lt.nbb) THEN
                      frct(1) = 0.
                      frct(2) = 1. 
                    ENDIF 

                    DO j=1,2
                      IF (frct(j).gt.0) THEN

                        m=carbo(j)
                        n=subs(j)
                        CALL findring(n,m,nca,tbond,rngflg,ring)
                        IF (rngflg==1) CYCLE
                        tbond(m,n) = 0
                        tbond(n,m) = 0

                        CALL addrx(progname,chem,nr,flag)
!                        CALL addreac(nr,progname,rdct(lco+1:lcf),flag)

                        nc = INDEX(tgroup(m),' ')
                        tgroup(m)(nc:nc) = '.'
                        nc = INDEX(tgroup(n),' ')
                        tgroup(n)(nc:nc) = '.'

! fragment and write in correct format:
!                        IF (wtflag.NE.0) WRITE(*,*) "fragm12"
                        CALL fragm(tbond,tgroup,pchem(nr,1),pchem(nr,2))
                        CALL stdchm(pchem(nr,1))
                        CALL stdchm(pchem(nr,2))

! check radicals, write in standard format:
                        CALL radchk(pchem(nr,1),rdckprod,rdcktprod,&
                                   nip,sc,nrxref(nr),rxref(nr,:))
                         tempkc = rdckprod(1)
                        IF (nip.NE.1) WRITE(6,*) '2 prod hvdiss2.f 33'
                        DO l=1,mca
                          coprod(nr,1,l) = rdcktprod(1,l)
                        ENDDO
                        CALL stdchm(tempkc)
                        pchem(nr,1) = tempkc
              
                        CALL radchk(pchem(nr,2),rdckprod,rdcktprod,&
                                   nip,sc,nrxref(nr),rxref(nr,:))
                         tempkc = rdckprod(1)
                        IF (nip.NE.1) WRITE(6,*) '2 prod hvdiss2.f 34'
                        CALL stdchm(tempkc)
                        pchem(nr,2) = tempkc
                        DO l=1,mca
                          coprod(nr,2,l) = rdcktprod(1,l)
                        ENDDO

! assign J ID
                        jid(nr)=23400
                        yield(nr)=frct(j)

! reset
                        tgroup(m)  = group(m)
                        tgroup(n)   = group(n)
                        tbond(m,n) = bond(m,n)
                        tbond(n,m) = bond(n,m)

                     ENDIF
                    ENDDO
                  ENDIF

! SECOND CHANNEL : MAKE MALEIC ANHYDRIDE + R.
! -----------------------------------------

                  CALL addrx(progname,chem,nr,flag)
!                  CALL addreac(nr,progname,rdct(lco+1:lcf),flag)

                  ic2 = 0
                  n1 = 0
! first fragment
! ==============

! find the carbon
                  DO i=1,mca
                    IF ( (tbond(ich,i).eq.1) .AND. (i.NE.ia) ) THEN
                      ia2=i
                    ENDIF
                  ENDDO
                 
                  tbond(ich,ia2) = 0
                  tbond(ia2,ich) = 0
                  nc = INDEX(tgroup(ia2),' ')
                  tgroup(ia2)(nc:nc) = '.'
! WARNING : RCO3 is only used to recognise the maleic anhydride fragment
                  tgroup(ich) = 'CO(OO.)'

!                  IF (wtflag.NE.0) WRITE(*,*) "fragm13"
                  CALL fragm(tbond,tgroup,tempkc,pchem(nr,1))
                  IF (INDEX(pchem(nr,1),'CO(OO.)').NE.0) THEN
                    pchem(nr,1)=tempkc 
                  ENDIF
   
! determine number of C in the fragment (ic2)
                  CALL stdchm(pchem(nr,1))
                  n1 = index(pchem(nr,1),' ')
                  DO j=1,n1
                    IF (index(pchem(nr,1)(j:j),'C').NE.0) ic2 = ic2 + 1
                  ENDDO
        
             CALL radchk(pchem(nr,1),rdckprod,rdcktprod,nip,sc,nrxref(nr),rxref(nr,:))
                   tempkc = rdckprod(1)
                  IF (nip.NE.1) WRITE(6,*) '2 produits hvdiss2.f 35'
                  DO j=1,mca
                    coprod(nr,1,j) = rdcktprod(1,j)
                  ENDDO
                  CALL stdchm(tempkc)
                  pchem(nr,1) = tempkc
! reset
                  tbond(ich,ia2) = bond(ich,ia2)
                  tbond(ia2,ich) = bond(ia2,ich)
                  tgroup(ia2) = group(ia2)
                  tgroup(ich) = group(ich)

! second fragment
! ===============

                  IF (tgroup(icarb)(1:3).EQ.'CHO') THEN
                    pchem(nr,2)='HO2'
                    jid(nr)=23300

! coprod XC is introduced to keep the total carbon balance
                    tempkc='xxlost'
                    CALL bratio(tempkc,brtio,tempnam,nrxref(nr),rxref(nr,:))
                    nxc = ic - 4 - ic2
                    DO i=1,mca
                      IF (coprod(nr,1,i).EQ.' ') GOTO 5
                    ENDDO
5                    DO k=i,i+(nxc-1)
                     coprod(nr,1,k) = tempnam
                    ENDDO
                       
                    
                  ELSE IF (tgroup(icarb)(1:3).EQ.'CO ') THEN
                    jid(nr)=23500
                    
                    DO i=1,mca
                      IF ( (bond(icarb,i).eq.1) .AND. (i.ne.ib) ) THEN
                        tbond(icarb,i) = 0
                        tbond(i,icarb) = 0

                        nc = INDEX(tgroup(i),' ')
                        tgroup(i)(nc:nc) = '.'
! WARNING : RCO3 is used to recognise the maleic anhydride fragment
                        tgroup(icarb) = 'CO(OO.)'
 
!                        IF (wtflag.NE.0) WRITE(*,*) "fragm14"
                        CALL fragm(tbond,tgroup,tempkc,pchem(nr,2))
                        IF (INDEX(pchem(nr,2),'CO(OO.)').NE.0) THEN
                          pchem(nr,2)=tempkc 
                        ENDIF
                        CALL stdchm(pchem(nr,2))

! ic2 = total of C of the 2 fragments (to determine 'XC')
                        n1 = index(pchem(nr,2),' ')
                        DO j=1,n1
                     IF (index(pchem(nr,2)(j:j),'C').NE.0) ic2 = ic2 + 1
                        ENDDO
                      
                        CALL radchk(pchem(nr,2),rdckprod,rdcktprod,&
                                   nip,sc,nrxref(nr),rxref(nr,:))
                         tempkc = rdckprod(1)
                        IF (nip.EQ.2) THEN
! JMLT changed from: sc_del(nr,2) = sc(1)
! and moved down 1 line
                          sc_del(nr,1) = sc(1)
                          flag_del(nr,2) = 1
                          pchem_del(nr,2) = rdckprod(2)
                          sc_del(nr,2) = sc(2)
                          DO j=1,mca
                            coprod_del(nr,2,j) = rdcktprod(2,j)
                          ENDDO
                        ENDIF
                        DO j=1,mca
                          coprod(nr,2,j) = rdcktprod(1,j)
                        ENDDO
                        CALL stdchm(tempkc)
                        CALL stdchm(pchem_del(nr,2))
                        pchem(nr,2) = tempkc

! coprod XC is introduced to keep the total carbon balance
                        tempkc='xxlost'
                        CALL bratio(tempkc,brtio,tempnam,nrxref(nr),rxref(nr,:))
                         nxc = ic - 4 - ic2
                         DO j=1,mca
                           IF (coprod(nr,2,j).EQ.' ') GOTO 6
                         ENDDO
6                         DO k=j,j+(nxc-1)
                           coprod(nr,2,k) = tempnam
                         ENDDO
                         
                      ENDIF
                    ENDDO
                  ENDIF

! add maleic anhydride as a coproduct (short name is used here instead of
! the long name) 
                  DO j=1,mca
                    IF (coprod(nr,1,j)(1:1).eq.' ') THEN
!c                      coprod(nr,1,j)='#-O1-COCdH=CdHC1O'
                       tempkc='-O1-COCdH=CdHC1O'
                       CALL bratio(tempkc,brtio,tempnam,nrxref(nr),rxref(nr,:))
                        coprod(nr,1,j)=tempnam
                        GOTO 33
                    ENDIF
                  ENDDO
33                CONTINUE            
                  
! reset
                  tgroup = group
                  tbond = bond

                  chromtab(icarb,1)=' '
                  GOTO 45
                ENDIF
                

! -----------------------------------------
! -CO-C=C-C=C-CO- CHROMOPHORE    -  CASE 3
! -----------------------------------------

                IF (icase.eq.3) THEN
                  CALL findring(icarb,ia,nca,tbond,rngflg,ring)
! DI-KETONES
! ==========
  
                  IF (tgroup(icarb)(1:3).EQ.'CO ') THEN
                    IF ((ring(ia).EQ.1).AND.(ring(ich).EQ.1)) GOTO 44

! FIRST CHANNEL : BREAK THE Cd-COR BOND AT FIRST -COR
! ---------------------------------------------------

                    CALL addrx(progname,chem,nr,flag)
!                    CALL addreac(nr,progname,rdct(lco+1:lcf),flag)

                    tbond(ich,ia) = 0
                    tbond(ia,ich) = 0

                    nc = INDEX(tgroup(ich),' ')
                    tgroup(ich)(nc:nc) = '.'
                    nc = INDEX(tgroup(ia),' ')
                    tgroup(ia)(nc:nc) = '.'

! fragment and write in correct format:
!                    IF (wtflag.NE.0) WRITE(*,*) "fragm15"
                    CALL fragm(tbond,tgroup,pchem(nr,1),pchem(nr,2))
                    CALL stdchm(pchem(nr,1))
                    CALL stdchm(pchem(nr,2))

! check radicals, write in standard format:
             CALL radchk(pchem(nr,1),rdckprod,rdcktprod,nip,sc,nrxref(nr),rxref(nr,:))
                     tempkc = rdckprod(1)
                    IF (nip.NE.1) WRITE(6,*) '2 produits hvdiss2.f 37'
                    DO j=1,mca
                      coprod(nr,1,j) = rdcktprod(1,j)
                    ENDDO
                    CALL stdchm(tempkc)
                    pchem(nr,1) = tempkc
              
             CALL radchk(pchem(nr,2),rdckprod,rdcktprod,nip,sc,nrxref(nr),rxref(nr,:))
                     tempkc = rdckprod(1)
                    IF (nip.NE.1) WRITE(6,*) '2 produits hvdiss2.f 38'
                    CALL stdchm(tempkc)
                    pchem(nr,2) = tempkc
                    DO j=1,mca
                      coprod(nr,2,j) = rdcktprod(1,j)
                    ENDDO

                    jid(nr)=22000

! reset
                    tgroup(ich)  = group(ich)
                    tgroup(ia)   = group(ia)
                    tbond(ich,ia) = bond(ich,ia)
                    tbond(ia,ich) = bond(ia,ich)
44                  CONTINUE
! SECOND CHANNEL : BREAK THE Cd-COR BOND AT SECOND -COR
! -----------------------------------------------------
                    IF ((ring(icarb).EQ.1).AND.(ring(i5).EQ.1)) GOTO 45
                    CALL addrx(progname,chem,nr,flag)
!                    CALL addreac(nr,progname,rdct(lco+1:lcf),flag)

                    tbond(icarb,i5) = 0
                    tbond(i5,icarb) = 0

                    nc = INDEX(tgroup(i5),' ')
                    tgroup(i5)(nc:nc) = '.'
                    nc = INDEX(tgroup(icarb),' ')
                    tgroup(icarb)(nc:nc) = '.'

! fragment and write in correct format:
!                    IF (wtflag.NE.0) WRITE(*,*) "fragm16"
                    CALL fragm(tbond,tgroup,pchem(nr,1),pchem(nr,2))
                    CALL stdchm(pchem(nr,1))
                    CALL stdchm(pchem(nr,2))

! check radicals, write in standard format:
             CALL radchk(pchem(nr,1),rdckprod,rdcktprod,nip,sc,nrxref(nr),rxref(nr,:))
                     tempkc = rdckprod(1)
                    IF (nip.NE.1) WRITE(6,*) '2 produits hvdiss2.f 39'
                    DO j=1,mca
                      coprod(nr,1,j) = rdcktprod(1,j)
                    ENDDO
                    CALL stdchm(tempkc)
                    pchem(nr,1) = tempkc
              
             CALL radchk(pchem(nr,2),rdckprod,rdcktprod,nip,sc,nrxref(nr),rxref(nr,:))
                     tempkc = rdckprod(1)
                    IF (nip.NE.1) WRITE(6,*) '2 produits hvdiss2.f 40'
                    CALL stdchm(tempkc)
                    pchem(nr,2) = tempkc
                    DO j=1,mca
                      coprod(nr,2,j) = rdcktprod(1,j)
                    ENDDO

                    jid(nr)=22000

! reset
                    tgroup = group
                    tbond = bond

! KETO-ALDEHYDE (SEE ALSO SAME ABOVE FROM ALDEHYDE SIDE) 
! =============  

                  ELSE IF (tgroup(icarb)(1:3).EQ.'CHO ') THEN

! FIRST CHANNEL : BREAK THE CD-CHO BOND AT -CHO
! ---------------------------------------------

                    CALL addrx(progname,chem,nr,flag)
!                    CALL addreac(nr,progname,rdct(lco+1:lcf),flag)

                    tbond(icarb,i5) = 0
                    tbond(i5,icarb) = 0

                    nc = INDEX(tgroup(icarb),' ')
                    tgroup(icarb)(nc:nc) = '.'
                    nc = INDEX(tgroup(i5),' ')
                    tgroup(i5)(nc:nc) = '.'

! fragment and write in correct format:
!                    IF (wtflag.NE.0) WRITE(*,*) "fragm17"
                    CALL fragm(tbond,tgroup,pchem(nr,1),pchem(nr,2))
                    CALL stdchm(pchem(nr,1))
                    CALL stdchm(pchem(nr,2))

! check radicals, write in standard format:
             CALL radchk(pchem(nr,1),rdckprod,rdcktprod,nip,sc,nrxref(nr),rxref(nr,:))
                     tempkc = rdckprod(1)
                    IF (nip.NE.1) WRITE(6,*) '2 produits hvdiss2.f 41'
                    DO j=1,mca
                      coprod(nr,1,j) = rdcktprod(1,j)
                    ENDDO
                    CALL stdchm(tempkc)
                    pchem(nr,1) = tempkc
              
             CALL radchk(pchem(nr,2),rdckprod,rdcktprod,nip,sc,nrxref(nr),rxref(nr,:))
                     tempkc = rdckprod(1)
                    IF (nip.NE.1) WRITE(6,*) '2 produits hvdiss2.f 42'
                    CALL stdchm(tempkc)
                    pchem(nr,2) = tempkc
                    DO j=1,mca
                      coprod(nr,2,j) = rdcktprod(1,j)
                    ENDDO

                    jid(nr)=22000

! reset
                    tgroup(icarb)  = group(icarb)
                    tgroup(i5)   = group(i5)
                    tbond(icarb,i5) = bond(icarb,i5)
                    tbond(i5,icarb) = bond(i5,icarb)

! SECOND CHANNEL : => >C=C-CO. + H. 
! ---------------------------------

                    CALL addrx(progname,chem,nr,flag)
!                    CALL addreac(nr,progname,rdct(lco+1:lcf),flag)

                    pold = 'CHO'
                    pnew = 'CO'
                    CALL swap(group(icarb),pold,tgroup(icarb),pnew)

                    nc = INDEX(tgroup(icarb),' ')
                    tgroup(icarb)(nc:nc) = '.'
                    CALL rebond(tbond,tgroup,tempkc,nring)
                  CALL radchk(tempkc,rdckprod,rdcktprod,nip,sc,nrxref(nr),rxref(nr,:))
                     pchem(nr,1) = rdckprod(1)
                    IF (nip.NE.1) WRITE(6,*) '2 produits hvdiss2.f 43'
                    DO j=1,mca
                      coprod(nr,1,j) = rdcktprod(1,j)
                    ENDDO
                    pchem(nr,2)='HO2 '

                    jid(nr)=22000

! reset
                    tgroup = group
                    tbond = bond

                  ENDIF
                  chromtab(icarb,1)=' '
                  chromtab(ich,1)=' '
                  GOTO 45
                ENDIF
              ENDIF

45          CONTINUE

! end of conjugated -CO- chromophore
            chromtab(ich,1)=' '
            GOTO 100
          
          ELSE        ! ELSE correponding to 'wwwwwww'
! if this point is reached , then "regular" -CO-
        
! --------------------------
! REGULAR R-CO-R CHROMOPHORE 
! --------------------------
! First, detect R-CO-C-Cd=Cd< structures , to force an unity quantum
! yield for the path forming stabilised radical
            conj=0
            cjpos = 0
            DO j=1,2
              ia=alpha(j)
              ia0=ia
              DO k=1,mca
                IF ((bond(k,ia).ne.0).and.(k.ne.ich)) THEN
                  IF (group(k)(1:2).EQ.'Cd') THEN
                    conj = conj + 1
                    cjpos = j
                  ENDIF
                ENDIF
              ENDDO 
            ENDDO
            IF (conj.GT.2) THEN
              WRITE(6,*) 'from subroutine hvdiss2'
              WRITE(6,*) 'problem with : conj > 2'
              WRITE(waru,*) 'hvdiss2',chem !STOP
            ENDIF  
            
! estimate bond dissociation energy            
!            dhfrco=heat(rdct(lco+1:lcf),nbson,bsongrp,bsonval)
            dhfrco=heat(chem)
            DO j=1,2
              ia=alpha(j)
              ia0 = ia
              ich = ii

! IF ring opening, temporarily rearrange molecule rather than dissociating
              IF(nring.GT.0)THEN
                CALL findring(ich,ia,nca,tbond,rngflg,ring)
              ENDIF
             
! fragment or open and write in correct format (don't retain co-products at this stage)
! LINEAR CHEM PHOTOLYSIS > RADICALS
              IF (rngflg.EQ.0) THEN
                tbond(ia,ich)=0
                tbond(ich,ia)=0
                nc = INDEX(tgroup(ich),' ')
                tgroup(ich)(nc:nc) = '.'
                nc = INDEX(tgroup(ia),' ')
                tgroup(ia)(nc:nc) = '.'
!                IF (wtflag.NE.0) WRITE(*,*) "fragm18"
                CALL fragm(tbond,tgroup,tempkc,tempkc2)
              ELSE
! CHEM AFTER RING-OPENING > DIRECT TO NON-RADICAL PRODUCTS
! break bond & add radical dots to separating groups:
!                IF(wtflag.NE.0)&
!                print*,'ring-opening by breaking R-CO bond'
                tbond(ia,ich)=0
                tbond(ich,ia)=0
                opflg=1
                nring=nring-1
                nc = INDEX(tgroup(ich),' ')
                tgroup(ich)(nc:nc) = '.'
                nc = INDEX(tgroup(ia),' ')
                tgroup(ia)(nc:nc) = '.'
! .. find if ring opens at artificial break point
                DO k=1,nca
                  IF(ring(k).NE.0) endrg=k
                ENDDO
                DO k=nca,1,-1
                  IF(ring(k).NE.0) begrg=k
                ENDDO
                IF((ich.EQ.begrg.AND.ia.EQ.endrg)&
                .OR.(ia.EQ.begrg.AND.ich.EQ.endrg))THEN
                  CONTINUE
! .. if opens elsewhere, rearrange new linear molecule
                ELSE
                  CALL rejoin(nca,ich,ia,m,n,tbond,tgroup)
                  ich=m
                  ia=n
                ENDIF
! do instantaneous ring-opening chemistry
                CALL openr(tbond,tgroup,nring,tempkc,tempkc2,tprod)
              ENDIF

              CALL stdchm(tempkc)
              CALL stdchm(tempkc2)

! calculate formation heat of the two fragments/opened ring and reaction heat
! R-CO-C-Cd=Cd< structures (values from CH3COCH2CdH=CdH2 but doesn't matter)
              IF (conj.EQ.2) THEN
                dhnet(1) = 72.1
                dhnet(2) = 72.1
              ELSE IF (conj.EQ.1) THEN
                IF (j.EQ.cjpos) THEN
                  dhnet(j) = 72.1
                  dhnet(3-j) = 81.2
                ENDIF  
              ELSE 
! calculate formation heat of the two fragments and reaction heat
!                dhfp1=heat(tempkc,nbson,bsongrp,bsonval)
!                dhfp2=heat(tempkc2,nbson,bsongrp,bsonval)
                dhfp1=heat(tempkc)
                dhfp2=heat(tempkc2)
                dhnet(1)=dhfp1+dhfp2-dhfrco
              ENDIF
! reset
              IF (opflg.eq.1) THEN
                tgroup = group
                tbond = bond
                ia=ia0
                ich=ii
                opflg=0
                rngflg=0
                nring=nring+1
              ELSE
                tbond(ia,ich)=bond(ia,ich)
                tbond(ich,ia)=bond(ich,ia)
                tgroup(ich)=group(ich)
                tgroup(ia)=group(ia)
              ENDIF
            ENDDO

! select only the most exothermic reaction (at more or less 0.5 Kcal)
            IF (dhnet(1).gt.dhnet(2)+0.5) THEN
              frct(1)=0.
              frct(2)=1.
            ELSE IF (dhnet(2).gt.dhnet(1)+0.5) THEN
              frct(1)=1.
              frct(2)=0.
            ELSE
              frct(1)=0.5
              frct(2)=0.5
            ENDIF

! check if substituted in alpha with -OH or -OOH
            pyr=0
            pyrpos=0
            DO j=1,2
              ia=alpha(j)
              IF (INDEX(group(ia),'CO(OH)').NE.0) THEN
                 pyr=1
                 pyrpos=j
              ENDIF
            ENDDO  

! check whether carbons in alpha are tertiary or not
            nter=0
            DO j=1,2
              ia=alpha(j)
              icount=0
              DO i=1,mca
                IF (tbond(ia,i).eq.1) THEN
                  icount=icount+1
                ENDIF
              ENDDO
              IF (icount.ge.3) THEN
                nter=nter+1
              ENDIF
            ENDDO  
        
! check if Norrish type 2 can be performed
                          
! find the gamma H and store the info in the hgamma table
! hgamma(j,1) <= carbon class (j is the carbon having gamma H)
! hgamma(j,2) <= carbon in beta 
! hgamma(j,3) <= carbon in alpha 
! if a sp2 carbon is found before gamma H, then do not take into account
! if > carbons before gamma H are in the same ring, do not take into account
            hsec=0
! initialise hgamma              
            DO i=1,mca
              DO j=1,3
                hgamma(i,j)=0
              ENDDO
              ksubs_alk(i)=0
            ENDDO

            DO 12  l=1,2
              ia=alpha(l)

! no Norrish if sterically hindered by presence of a ring
! (2 or more of ii/alpha/beta/gamma are on ring)
! NOTE: Norrish II may be possible for ring species with complex branches
!       We assume only simple branches here, and disallow Norrish as soon
!       as a ring detected in {ii/alpha/beta/gamma}

              IF(nring.GT.0)THEN
! .. bond(ii,ia)
                CALL findring(ii,ia,nca,bond,rngflg,ring)
                IF(rngflg.GT.0) GOTO 12 ! Norrish II disallowed
! .. bond(ia,ib) 
                DO i=1,mca
                  IF ((tbond(ia,i).GT.0) .AND. (i.NE.ich)) THEN
                    ib=i
                    CALL findring(ia,ib,nca,bond,rngflg,ring)
                    IF(rngflg.GT.0) GOTO 12 ! Norrish II disallowed

! .. bond(ib,ig): some might be valid: correct later
                    DO k=1,mca
                      IF ((tbond(ia,k).GT.0) .AND. (k.NE.ich)) THEN
                        ib=i
                        DO j=1,mca
                          IF ((tbond(ib,j).GT.0) .AND. (j.NE.ia)) THEN
                            ig=j
                            CALL findring(ib,ig,nca,bond,rngflg,ring)
                            IF(rngflg.GT.0) GOTO 12 ! Norrish II disallowed
                          ENDIF
                        ENDDO
                      ENDIF
                    ENDDO
                  ENDIF ! correction, JLT, Apr 17
                ENDDO ! correction, JLT, Apr 17
              ENDIF

              kgp_ia(l) = 0
         IF(((INDEX(group(ia),'(OH)').NE.0).and.&
             (INDEX(group(ia),'CO(OH)').EQ.0))    .OR.&
             ((INDEX(group(ia),'(OOH)').NE.0).and.&
             (INDEX(group(ia),'CO(OOH)').EQ.0)))   GOTO 12      
         
         IF(INDEX(group(ia),'(ONO2)').NE.0) kgp_ia(l)=kgp_ia(l) + 1
         IF(INDEX(group(ia),'(ONO2)(ONO2)').NE.0) kgp_ia(l)=kgp_ia(l)+1
                
              DO i=1,mca
                kgp_i(l) = 0
                IF ( (bond(ia,i).EQ.1) .AND. (i.ne.ich) ) THEN
                  IF (group(i)(1:2).eq.'CO') goto 12
                  IF (group(i)(1:2).eq.'Cd') goto 12
        IF(INDEX(group(i),'(OH)').NE.0) kgp_i(l) = kgp_i(l) + 1
        IF(INDEX(group(i),'(ONO2)').NE.0) kgp_i(l) = kgp_i(l) + 1
        IF(INDEX(group(i),'(OOH)').NE.0) kgp_i(l) = kgp_i(l) + 1
        IF(INDEX(group(i),'(OH)(OH)').NE.0) kgp_i(l) = kgp_i(l) + 1
        IF(INDEX(group(i),'(ONO2)(ONO2)').NE.0) kgp_i(l) = kgp_i(l) + 1
        IF(INDEX(group(i),'(OOH)(OOH)').NE.0) kgp_i(l) = kgp_i(l) + 1
                DO j=1,mca
                  kgp_j(l) = 0
                  IF ( (bond(i,j).eq.1) .AND. (j.ne.ia) ) THEN
                    IF (group(j)(1:2).eq.'CO') GOTO 12
                    IF (group(j)(1:2).eq.'Cd') GOTO 12

! BA : exclude -O- group bonded to gamma node
                    DO kk=1,mca
                      IF (bond(j,kk)/=0) THEN
                        IF (group(kk)=='-O-') GOTO 12
                      ENDIF
                    ENDDO                  

        IF(INDEX(group(j),'(OH)').NE.0) kgp_j(l) = kgp_j(l) + 1
        IF(INDEX(group(j),'(ONO2)').NE.0) kgp_j(l) = kgp_j(l) + 1
        IF(INDEX(group(j),'(OOH)').NE.0) kgp_j(l) = kgp_j(l) + 1
        IF(INDEX(group(j),'(OH)(OH)').NE.0) kgp_j(l) = kgp_j(l) + 1
        IF(INDEX(group(j),'(ONO2)(ONO2)').NE.0) kgp_j(l) = kgp_j(l) + 1
        IF(INDEX(group(j),'(OOH)(OOH)').NE.0) kgp_j(l) = kgp_j(l) + 1
                    IF ((kgp_ia(l).le.1).and.(kgp_i(l).le.1).and.&
                          (kgp_j(l).le.1)) THEN
                      IF (group(j)(1:3).eq.'CH3') THEN
                        hgamma(j,1)=1
                        hgamma(j,2)=i
                        hgamma(j,3)=ia
                      ELSE IF (group(j)(1:3).eq.'CHO') THEN
                        hgamma(j,1)=0
                      ELSE IF (group(j)(1:3).eq.'CH2') THEN
                        hgamma(j,1)=2
                        hgamma(j,2)=i
                        hgamma(j,3)=ia
                        hsec=1
                      ELSE IF (group(j)(1:2).eq.'CH') THEN
                        hgamma(j,1)=3
                        hgamma(j,2)=i
                        hgamma(j,3)=ia
                        hsec=1
                      ENDIF
                      IF ((kgp_ia(l).GE.1).or.(kgp_i(l).GE.1).or.&
                         (kgp_j(l).GE.1)) ksubs_alk(j)=1
                      ENDIF
                    ENDIF
                  ENDDO
                ENDIF
              ENDDO
12          CONTINUE             
       
! remove from a CH3 if beta carbon is CH or CH2
            IF (hsec.gt.0) THEN
              DO i=1,mca
                IF (hgamma(i,1).eq.1) hgamma(i,1)=0
              ENDDO
            ENDIF

            npos=0
            DO i=1,mca
              IF (hgamma(i,1).gt.0) THEN
                npos=npos+1
              ENDIF
            ENDDO

! RCO-C(O)OH
! -------------

            IF (pyr.ne.0) THEN
              CALL addrx(progname,chem,nr,flag)
!              CALL addreac(nr,progname,rdct(lco+1:lcf),flag)
              
              ia=alpha(pyrpos)
              tbond(ich,ia)=0
              tbond(ia,ich)=0
              tgroup(ia)=' '
              
              pold='CO'
              pnew='CHO'
              
              CALL swap(group(ich),pold,tgroup(ich),pnew)
              CALL rebond(tbond,tgroup,tempkc,nring)
              pchem(nr,1)=tempkc
              CALL stdchm(pchem(nr,1))
              pchem(nr,2)='CO2'
              
              jid(nr)=32100
              
! reset
              tgroup(ich)  = group(ich)
              tgroup(ia)   = group(ia)
              tbond(ich,ia) = bond(ich,ia)
              tbond(ia,ich) = bond(ia,ich)
              GOTO 100
                
            ENDIF   

! FIRST CHANNEL - BREAK THE C-CO BONDS
! ------------------------------------

            DO j=1,2
              IF (frct(j).gt.0) THEN
                CALL addrx(progname,chem,nr,flag)
!                CALL addreac(nr,progname,rdct(lco+1:lcf),flag)

                ia=alpha(j)
                ia0=ia
                ich=ii

! IF ring opening, rearrange molecule rather than dissociating
! (first re-set tgroup, tbond to original values if second time through)
                IF(nring.GT.0)THEN
                  IF(j.EQ.2)THEN
                    DO k=1,mca
                      tgroup(k)=group(k)
                      DO l=1,mca
                        tbond(k,l)=bond(k,l)
                      ENDDO
                    ENDDO
                    opflg=0
                  ENDIF
                  CALL findring(ich,ia,nca,tbond,rngflg,ring)
                ENDIF

! fragment or open and write in correct format:
                IF (rngflg.EQ.0) THEN
! LINEAR CHEM PHOTOLYSIS > RADICALS
! break bond & add radical dots to separating groups:
                  tbond(ia,ich)=0
                  tbond(ich,ia)=0
                  nc = INDEX(tgroup(ich),' ')
                  tgroup(ich)(nc:nc) = '.'
                  nc = INDEX(tgroup(ia),' ')
                  tgroup(ia)(nc:nc) = '.'
!                  IF (wtflag.NE.0) WRITE(*,*) "fragm19"
                  CALL fragm(tbond,tgroup,pchem(nr,1),pchem(nr,2))
                  CALL stdchm(pchem(nr,1))
                  CALL stdchm(pchem(nr,2))
! check radicals, write in standard format:
             CALL radchk(pchem(nr,1),rdckprod,rdcktprod,nip,sc,nrxref(nr),rxref(nr,:))
                   tempkc = rdckprod(1)
                  CALL stdchm(tempkc)
                  pchem(nr,1) = tempkc
                  IF (nip.EQ.2) THEN
                    sc_del(nr,1) = sc(1)
                    flag_del(nr,1) = 1
                    pchem_del(nr,1) = rdckprod(2)
                    sc_del(nr,2) = sc(2)
                    DO k=1,mca
                      coprod_del(nr,1,k) = rdcktprod(2,k)
                    ENDDO
                  ENDIF
                  DO k=1,mca 
                    coprod(nr,1,k) = rdcktprod(1,k)
                  ENDDO
                  CALL radchk(pchem(nr,2),rdckprod,rdcktprod,nip,sc,nrxref(nr),rxref(nr,:))
                   tempkc2 = rdckprod(1)
                  CALL stdchm(tempkc2)
                  pchem(nr,2) = tempkc2
                  IF (nip.EQ.2) THEN
                    sc_del(nr,1) = sc(1)
                    flag_del(nr,2) = 1
                    pchem_del(nr,2) = rdckprod(2)
                    sc_del(nr,2) = sc(2)
                    DO k=1,mca
                      coprod_del(nr,1,k) = rdcktprod(2,k)
                    ENDDO
                  ENDIF
                  DO k=1,mca
                    coprod(nr,2,k) = rdcktprod(1,k)
                  ENDDO
                ELSE
! CHEM AFTER RING-OPENING > REARRANGE DIRECT TO NON-RADICAL PRODUCTS
! (uniqring should not rearrange molecule at this point)
! break bond & add radical dots to separating groups:
!                  IF(wtflag.NE.0)&
!                  print*,'ring-opening by breaking R-CO bond II'
                  tbond(ia,ich)=0
                  tbond(ich,ia)=0
                  opflg=1
                  nc = INDEX(tgroup(ich),' ')
                  tgroup(ich)(nc:nc) = '.'
                  nc = INDEX(tgroup(ia),' ')
                  tgroup(ia)(nc:nc) = '.'
                  nring=nring-1
! .. find if ring opens at artificial break point
                  DO k=1,nca
                    IF(ring(k).NE.0) endrg=k
                  ENDDO
                  DO k=nca,1,-1
                    IF(ring(k).NE.0) begrg=k
                  ENDDO
                  IF((ich.EQ.begrg.AND.ia.EQ.endrg)&
                  .OR.(ia.EQ.begrg.AND.ich.EQ.endrg))THEN
                    CONTINUE
                  ELSE
                    CALL rejoin(nca,ich,ia,m,n,tbond,tgroup)
                    ich=m
                    ia=n
                  ENDIF
! do instantaneous ring-opening chemistry
                  CALL openr(tbond,tgroup,nring,tempkc,tempkc2,tprod)
                  CALL stdchm(tempkc)
                  CALL stdchm(tempkc2)
                  pchem(nr,1) = tempkc
                  pchem(nr,2) = tempkc2
                  DO k=1,mca
                    coprod(nr,1,k) = tprod(k)
                  ENDDO
                ENDIF

! assign J ID
                IF (((index(tgroup(ia),'(OH)').NE.0).and.&
                   (index(tgroup(ia),'CO(OH)').EQ.0))   .OR.   &     
                  ((index(tgroup(ia),'(OOH)').NE.0).and.&
                   (index(tgroup(ia),'CO(OOH)').EQ.0))) THEN        
                    jid(nr)=31900
                ELSE                                
                  IF (npos.eq.0) THEN        ! no gammaH
                     IF (nter.eq.2) jid(nr)=31100
                     IF (nter.eq.1) jid(nr)=30600
                     IF (nter.eq.0) jid(nr)=30100
                  ELSE                       ! with gammaH
                     IF (nter.eq.2) jid(nr)=31200
                     IF (nter.eq.1) jid(nr)=30700
                     IF (nter.eq.0) jid(nr)=30200
                  ENDIF
                ENDIF
                yield(nr)=frct(j)

! reset
                IF (opflg.eq.1) THEN
                  tgroup = group
                  tbond = bond
                  ia=ia0
                  ich=ii
                  opflg=0
                  nring=nring+1
                ELSE
                  tgroup(ich)  = group(ich)
                  tgroup(ia)   = group(ia)
                  tbond(ich,ia) = bond(ich,ia)
                  tbond(ia,ich) = bond(ia,ich)
                ENDIF
              ENDIF
            ENDDO

              
! SECOND CHANNEL - NORRISH TYPE 2
! -------------------------------

           IF (npos.gt.0) THEN
              xcoeff=1./real(npos)
              DO i=1,mca
                IF (hgamma(i,1).gt.0) THEN
                  CALL addrx(progname,chem,nr,flag)
!                  CALL addreac(nr,progname,rdct(lco+1:lcf),flag)

! change gamma position (i.e. holding the H):
                  IF (group(i)(1:3).EQ.'CH3') THEN
                    pold = 'CH3'
                    pnew = 'CdH2'
                  ELSE IF(group(i)(1:3).EQ.'CH2') THEN
                    pold = 'CH2'
                    pnew = 'CdH'
                  ELSE IF(group(i)(1:2).EQ.'CH') THEN
                    pold = 'CH'
                    pnew = 'Cd'
                  ELSE      
                    WRITE(6,'(a)') '--error--'
                    WRITE(6,'(a)') 'from ROUTINE:hvdiss2 '
                    WRITE(6,'(a)') 'molecule couldn t be identified:'
                    WRITE(6,'(a)') chem
                    WRITE(6,'(a)') 'NORRISH2 was treated - first pb'
                    WRITE(waru,*) 'hvdiss2',chem !STOP
                  ENDIF
                         
                  CALL swap(group(i),pold,tgroup(i),pnew)
             
! change beta double-bond carbons:
                  pold = 'C'
                  pnew = 'Cd'
                  ib=hgamma(i,2)
                  CALL swap(group(ib),pold,tgroup(ib),pnew)

! change alpha position
                  check=0
                  ia=hgamma(i,3)
                  IF (tgroup(ia)(1:3).eq.'CH2') THEN
                    pold='CH2'
                    pnew='CH3'
                    check=1
                  ELSE IF (tgroup(ia)(1:3).eq.'CO ') THEN
                    pold='CO'
                    pnew='CHO'
                    check=1
                  ELSE IF (tgroup(ia)(1:2).eq.'CH') THEN
                    pold='CH'
                    pnew='CH2'
                    check=1
                  ELSE IF (tgroup(ia)(1:2).eq.'C ') THEN
                    pold='C'
                    pnew='CH'
                    check=1
                  ELSE IF (tgroup(ia)(1:2).eq.'C(') THEN
                    pold='C('
                    pnew='CH('
                    check=1
                  ENDIF            

                  IF (check.EQ.1) THEN
                    CALL swap(group(ia),pold,tgroup(ia),pnew)
                  ELSE  
                    WRITE(6,'(a)') '--error--'
                    WRITE(6,'(a)') 'from ROUTINE:hvdiss2 '
                    WRITE(6,'(a)') 'molecule couldn t be identified:'
                    WRITE(6,'(a)') chem
                    WRITE(6,'(a)') 'NORRISH2 was treated - second pb'
                    WRITE(waru,*) 'hvdiss2',chem !STOP
                  ENDIF

! new bond matrix:
                  tbond(i,ib) = 2
                  tbond(ib,i) = 2
                          
                  tbond(ia,ib) = 0
                  tbond(ib,ia) = 0
 
! fragment:
!                  IF (wtflag.NE.0) WRITE(*,*) "fragm20"
                  CALL fragm(tbond,tgroup,pchem(nr,1),pchem(nr,2))
! find correct names:
                  CALL stdchm(pchem(nr,1))
                  CALL stdchm(pchem(nr,2))
                    
! treat substitued alkenes                
                 IF (ksubs_alk(i).ne.0) THEN
                   temp_coprod(1)=' '
                   temp_coprod(2)=' '
                   CALL alkcheck(pchem(nr,1),temp_coprod(1),acom)
                   CALL alkcheck(pchem(nr,2),temp_coprod(2),acom)
                  ENDIF

! add coprod and check radicals, write in standard format:
!c                  WRITE(44,*) 'rdct=',rdct(1:50)
!c                  WRITE(44,*) 'pchem(nr,1)=',pchem(nr,1)(1:50)
!c                  WRITE(44,*) '  temp_coprod(1)=',temp_coprod(1)
!c                  WRITE(44,*) 'pchem(nr,2)=',pchem(nr,2)(1:50)
!c                  WRITE(44,*) '  temp_coprod(2)=',temp_coprod(2)
                  IF (INDEX(pchem(nr,1),'.').ne.0) THEN
!c                    CALL radchk(pchem(nr,1),tempkc,tprod)
             CALL radchk(pchem(nr,1),rdckprod,rdcktprod,nip,sc,nrxref(nr),rxref(nr,:))
                     tempkc = rdckprod(1)
                    IF (nip.NE.1) WRITE(6,*) '2 produits hvdiss2.f 46'
                    CALL stdchm(tempkc)
                    pchem(nr,1) = tempkc
                    DO j=1,mca
                      coprod(nr,1,j) = rdcktprod(1,j)
                      IF (coprod(nr,1,j).eq.' ') THEN
                        coprod(nr,1,j) = temp_coprod(1)
                        temp_coprod(1)=' '
                      ENDIF
                    ENDDO
                  ELSE 
                    coprod(nr,1,1)=temp_coprod(1)
                  ENDIF
                 
                  IF (INDEX(pchem(nr,2),'.').ne.0) THEN
             CALL radchk(pchem(nr,2),rdckprod,rdcktprod,nip,sc,nrxref(nr),rxref(nr,:))
                     tempkc = rdckprod(1)
                    IF (nip.NE.1) WRITE(6,*) '2 produits hvdiss2.f 47'
                    CALL stdchm(tempkc)
! next line added during the radchk modification
                    pchem(nr,2) = tempkc
                    DO j=1,mca
                      coprod(nr,2,j) = rdcktprod(1,j)
                      IF (coprod(nr,2,j).eq.' ') THEN
                        coprod(nr,2,j) = temp_coprod(2)
                        temp_coprod(2)=' '
                      ENDIF
                    ENDDO
                  ELSE 
                    coprod(nr,2,1)=temp_coprod(2)
                  ENDIF
       
! reset
                  tgroup = group
                  tbond = bond
! assign J ID
                  IF (nter.eq.2) jid(nr)=31300
                  IF (nter.eq.1) jid(nr)=30800
                  IF (nter.eq.0) jid(nr)=30300
                  xcoeff=1./real(npos)
      
                ENDIF
              ENDDO
            ENDIF    ! Norrish II
          chromtab(ich,1)=' '
          GOTO 100       
         
          ENDIF      ! correponding if statement'wwwwww'

! end of 'k' chromophore        
       ENDIF                    ! corresponding if statement = '@@@'

! OTHER "SIMPLE" CHROMOPHORE
! ==========================
       DO 66 l=1,4
          IF (chromtab(ich,l).EQ.' ') GOTO 66

! ===================
! C(OOH) CHROMOPHORE
! ===================

          IF (chromtab(ich,l).EQ.'h') THEN
            CALL addrx(progname,chem,nr,flag)
!            CALL addreac(nr,progname,rdct(lco+1:lcf),flag)
            
            pold='(OOH)'
            pnew='(O.)'
            CALL swap(group(ich),pold,tgroup(ich),pnew)
            CALL rebond(tbond,tgroup,tempkc,nring)
            CALL radchk(tempkc,rdckprod,rdcktprod,nip,sc,nrxref(nr),rxref(nr,:))
             pchem(nr,1) = rdckprod(1)
            IF (nip.NE.1) WRITE(waru,*) '2 produits hvdiss2.f 48'
            CALL stdchm(pchem(nr,1))
            DO j=1,mca
              coprod(nr,1,j)=rdcktprod(1,j)
            ENDDO
            
            pchem(nr,2)='HO  '
            
            jid(nr)=40100

! reset
            tgroup(ich)=group(ich)
          ENDIF


! ====================
! C(ONO2) CHROMOPHORE
! ====================

          IF (chromtab(ich,l).EQ.'n') THEN

! find the ketones in alpha
            icount=0
            DO i=1,mca
              IF ((tbond(ich,i).EQ.1).AND.(group(i)(1:3).EQ.'CO ')) THEN
                  icount=icount+1
              ENDIF
            ENDDO
              
            IF (icount.NE.0) THEN
! keto nitrate photolysis added in the frame of oncem projec
              CALL addrx(progname,chem,nr,flag)
!              CALL addreac(nr,progname,rdct(lco+1:lcf),flag)
!! change RONO2 in RO + NO2, rebuild and rename
! terminal nitrate
              IF (group(ich).EQ.'CH2(ONO2)') THEN
                pold='CH2(ONO2)'
                pnew='CH2(O.)'
                call swap(group(ich),pold,tgroup(ich),pnew)
                call rebond(tbond,tgroup,tempkc,nring)
                pchem(nr,1)=tempkc
                CALL stdchm(pchem(nr,1))
                pchem(nr,2)='NO2'

                jid(nr)=40001

! reset
                tgroup(ich)   = group(ich)

                chromtab(ich,1)=' '
              ELSE
! internal nitrate
!                IF (group(ich).EQ.'CH(ONO2)') THEN 
!                  pold='CH(ONO2)'
!                  pnew='CH(O.)'
                  pold='(ONO2)'
                  pnew='(O.)'
!                ELSE IF (group(ich).EQ.'C(ONO2)') THEN 
!                  pold='C(ONO2)'
!                  pnew='C(O.)'
!                ENDIF

                call swap(group(ich),pold,tgroup(ich),pnew)
                call rebond(tbond,tgroup,tempkc,nring)
                pchem(nr,1)=tempkc
                CALL stdchm(pchem(nr,1))
                pchem(nr,2)='NO2'

                jid(nr)=40002

! reset
                tgroup(ich)   = group(ich)

                chromtab(ich,1)=' '
              ENDIF
!!!!!!!!!! regular nitrate
            ELSE
            CALL addrx(progname,chem,nr,flag)
!            CALL addreac(nr,progname,rdct(lco+1:lcf),flag)
            
            pold='(ONO2)'
            pnew='(O.)'
            CALL swap(group(ich),pold,tgroup(ich),pnew)
            CALL rebond(tbond,tgroup,tempkc,nring)
            CALL radchk(tempkc,rdckprod,rdcktprod,nip,sc,nrxref(nr),rxref(nr,:))
             pchem(nr,1) = rdckprod(1)
            IF (nip.NE.1) WRITE(waru,*) '2 produits hvdiss2.f'
            CALL stdchm(pchem(nr,1))
            DO j=1,mca
              coprod(nr,1,j)=rdcktprod(1,j)
            ENDDO
            
            pchem(nr,2)='NO2 '
            
            IF (group(ich)(1:3).eq.'CH2') THEN
              jid(nr)=10100
            ELSE IF (group(ich)(1:2).eq.'CH') THEN
              jid(nr)=10200
            ELSE 
              jid(nr)=10300
            ENDIF  

! reset
            tgroup(ich)=group(ich)

            ENDIF
          ENDIF

66      CONTINUE
          
        
100   CONTINUE
             
400   CONTINUE

!**********************************************************************
!                            WRITE REACTIONS                         *
!**********************************************************************
       
! check for duplicate reactions and eliminate reaction
      DO 200 i=1,nr-1
         DO 210 j=i+1,nr
            IF (pchem(j,1).EQ.pchem(i,1)  .AND.&
               pchem(j,2).EQ.pchem(i,2)) THEN
               DO k=1,mca
                 IF (coprod(i,1,k).NE.coprod(j,1,k)) GO TO 210
                 IF (coprod(i,2,k).NE.coprod(j,2,k)) GO TO 210
               ENDDO
               IF (jid(i).ne.jid(j)) GO TO 210
               flag(j) = 0
               yield(i) = yield(i) + yield(j)
            ENDIF
210      CONTINUE
        
         DO 220 j=i+1,nr
            IF (pchem(j,1).EQ.pchem(i,2)  .AND.&
               pchem(j,2).EQ.pchem(i,1)) THEN
               DO k=1,mca
                  IF (coprod(i,1,k).NE.coprod(j,2,k)) GO TO 220
                  IF (coprod(i,2,k).NE.coprod(j,1,k)) GO TO 220
               ENDDO
               IF (jid(i).ne.jid(j)) GO TO 220
               flag(j) = 0
               yield(i) = yield(i) + yield(j)
            ENDIF
220      CONTINUE
200   CONTINUE
        
!*****************************
! eliminate minor reactions *
!*****************************
      IF (nr==0) RETURN  ! no reaction found
       
! calculate total photolysis rate (phot_tot)    
      phot_tot = 0.
      IF (nr.GT.0)THEN
        IF(jid(nr).GE.10000.AND.known_species.EQ.0) THEN
!          DO i = 1,mnr
          DO i = 1,nr
            DO j = 1,42
! BA: bug here, the following condition never met            
              IF (jlab40(j).EQ.jid(i)) THEN
                phot_tot= phot_tot + J40(j)
              ENDIF
            ENDDO
          ENDDO
        ENDIF
      ENDIF
     
! determine each branching ratios and delete reaction when < 5%
     
      IF (nr.GT.0)THEN
        IF(jid(nr).GE.10000.AND.known_species.EQ.0) THEN
!          DO i = 1,mnr
          DO i = 1,nr
            DO j = 1,42
! BA: bug here, the following condition never met            
              IF (jlab40(j).eq.jid(i)) THEN
                br(i) = J40(j)*yield(i) / phot_tot
                IF (br(i).lt.cut_off) flag(i) = 0
              ENDIF
            ENDDO
          ENDDO
        ENDIF
      ENDIF
! BA: following statements to be revisited once the bug above has been
! solved !!!
!      kmax=MAXVAL(br)
!      IF (kmax==0) THEN
!        WRITE(61,*) 'no HV for:',TRIM(rdct)
!        RETURN
!      ENDIF
!      IF (SUM(flag)==0) THEN
!         WHERE (br.GT.kmax/2.) flag = 1
!         sumflg=SUM(flag)
!         IF (sumflg==0) STOP 'in hvdiss2, no path'
!      ENDIF
     

! count total number of reactions:
      nch = 0
      DO i=1,mnr
         IF (flag(i).EQ.1) nch = nch + 1
      ENDDO
       
! output
! ------
      ich = 10
      DO 300 i=1,mnr
         IF (flag(i).EQ.0) GO TO 300
      
      IF (wrtopeinfo) & !Write information required for operator&
       WRITE(10,'(A15,A1,A6,A4,1X,a5,i5,a3,f4.2,a1)')'**** INIT HV + '&
      ,'G',idnam,'*****',' HV /',jid(i),'   ',yield(i),'/'
         
! initialize reaction
         CALL rxinit(r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe)
         ich = ich + 1
         r(1) = idnam
         r(2) = 'HV   '

! cut off for the reaction need to be reactivated in the future

! current program is not "well" done (need to be revised latter)
! problem is that hot criegge may be formed in the routine
! If that case occurs, then it should appears only as a first 
! product. Hot criegge are therefore checked as first product
         np=0

         IF (INDEX(pchem(i,1),'.(OO.)*').EQ.0) THEN
           CALL add1tonp(progname,chem,np)
!          CALL addprod(np,progname,rdct(lco+1:lcf))
           s(np)=1.
!c           brtio = brch
           CALL bratio(pchem(i,1),brtio,p(np),nrxref(i),rxref(i,:))

           IF (wrtopeinfo) WRITE(10,'(f5.3,2X,A1,A6)') s(np),'G', p(np)
     
           IF (flag_del(i,1).NE.0) THEN
             CALL add1tonp(progname,chem,np)
!             CALL addprod(np,progname,rdct(lco+1:lcf))
             s(np-1)= sc_del(i,1)
             s(np) = sc_del(i,2)
             CALL bratio(pchem_del(i,1),brtio,p(np),nrxref(i),rxref(i,:))
           ENDIF
           IF (rdtcopchem.GT.0.) THEN
             np=np+1 ; IF (np>mnr) STOP "in hvdiss2, np>mnr"
!             CALL addcoprod(np,progname,rdct(lco+1:lcf))
             s(np) = rdtcopchem
             p(np) = copchem
             IF (wrtopeinfo) WRITE(10,'(f5.3,2X,A1,A6)') s(np),'G', p(np)
           ENDIF

         ELSE
           CALL xcrieg(pchem(i,1),brtio,ss,pp,cut_off,nrxref(i),rxref(i,:))
           DO j=1,mnr
             IF (pp(j)(1:1).NE.' ') THEN
               CALL add1tonp(progname,chem,np)
!               CALL addprod(np,progname,rdct(lco+1:lcf))
               s(np) = ss(j)
               p(np) = pp(j)
               IF (wrtopeinfo) WRITE(10,'(f5.3,2X,A1,A6)') s(np),'G', p(np)
             ENDIF
           ENDDO
         ENDIF
         IF (pchem(i,2)(1:1).NE.' ') THEN
            CALL add1tonp(progname,chem,np)
!            CALL addprod(np,progname,rdct(lco+1:lcf))
            s(np) = 1.
            CALL bratio(pchem(i,2),brtio,p(np),nrxref(i),rxref(i,:))
            IF (wrtopeinfo) WRITE(10,'(f5.3,2X,A1,A6)') s(np),'G', p(np)
            IF (flag_del(i,2).NE.0) THEN
              CALL add1tonp(progname,chem,np)
!              CALL addprod(np,progname,rdct(lco+1:lcf))
              s(np-1)= sc_del(i,1)
              s(np) = sc_del(i,2)
              CALL bratio(pchem_del(i,2),brtio,p(np),nrxref(i),rxref(i,:))
            ENDIF
            IF (rdtcopchem.GT.0.) THEN
              np=np+1 ; IF (np>mnr) STOP "in hvdiss2, np>mnr"
!              CALL addcoprod(np,progname,rdct(lco+1:lcf))
              s(np) = rdtcopchem
              p(np) = copchem
              IF (wrtopeinfo) WRITE(10,'(f5.3,2X,A1,A6)') s(np),'G', p(np)
            ENDIF
         ENDIF

         DO j=1,mca
            IF (coprod(i,1,j)(1:1).NE.' ') THEN
              np=np+1 ; IF (np>mnr) STOP "in hvdiss2, np>mnr"
!               CALL addcoprod(np,progname,rdct(lco+1:lcf))
! following test on coef. sto. added because of criegee above
               IF (s(np).eq. 0.) s(np)=1.
               p(np) = coprod(i,1,j)
               IF (wrtopeinfo) WRITE(10,'(f5.3,2X,A1,A6)') s(np),'G', p(np)
            ENDIF
            IF (coprod(i,2,j)(1:1).NE.' ') THEN
               np=np+1 ; IF (np>mnr) STOP "in hvdiss2, np>mnr"
!                CALL addcoprod(np,progname,rdct(lco+1:lcf))
               s(np) = 1.
               p(np) = coprod(i,2,j)
               IF (wrtopeinfo) WRITE(10,'(f5.3,2X,A1,A6)')s(np),'G',p(np)
            ENDIF
            IF ((coprod_del(i,1,j)(1:1).NE.' ').AND.&
               (flag_del(i,1).NE.0)) THEN
               np=np+1 ; IF (np>mnr) STOP "in hvdiss2, np>mnr"
!               CALL addcoprod(np,progname,rdct(lco+1:lcf))
               s(np) = sc_del(i,2)
               p(np) = coprod_del(i,1,j)
               IF (wrtopeinfo) WRITE(10,'(f5.3,2X,A1,A6)')s(np),'G',p(np)
            ENDIF
            IF ((coprod_del(i,2,j)(1:1).NE.' ').AND.&
               (flag_del(i,2).NE.0)) THEN
               np=np+1 ; IF (np>mnr) STOP "in hvdiss2, np>mnr"
!               CALL addcoprod(np,progname,rdct(lco+1:lcf))
               s(np) = 1.
               p(np) = coprod_del(i,2,j)
               IF (wrtopeinfo) WRITE(10,'(f5.3,2X,A1,A6)')s(np),'G',p(np)
            ENDIF
         ENDDO

         arrh(1) = yield(i)
         
!  write out:
         ar1bis = 1.
!c         CALL rxwrit2(a1,a2,a3,a4,r,s,p,
!c     &                ar1,ar2,ar3,f298,fratio,15)

!* write jID         
!c         write(15,'(a5,i5,a3,f4.2,a1)') ' HV /',jid(i),'   ',ar1,'/'
!c         write(16,'(a5,i5,a3,f4.2,a1)') ' HV /',jid(i),'   ',ar1,'/'

!*  write out: HV reaction => idreac=1
         idreac=1
         nlabel=jid(i)
         xlabel=arrh(1)
         CALL rxwrit(mecu,r,s,p,arrh,idreac,nlabel,xlabel,folow,fotroe,rxref(nr,:))


      IF (wrtopeinfo) WRITE(10,*)'end'
300   CONTINUE
!      IF(wtflag.NE.0) WRITE(*,*) '----  end of hvdiss2'
!* end of hvdiss
      RETURN
      END SUBROUTINE hvdiss2
        
