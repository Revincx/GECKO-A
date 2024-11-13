MODULE outtool
IMPLICIT NONE
CONTAINS
!SUBROUTINE dictelement()
!SUBROUTINE diffusion_vol()
!SUBROUTINE wrt_dict()
!SUBROUTINE wrt_mxyield()
!SUBROUTINE wrt_psat(myrflg,nanflg,simflg)
!SUBROUTINE wrt_henry()
!SUBROUTINE wrt_depo()
!SUBROUTINE wrt_heatf()
!SUBROUTINE wrt_Tg()
!SUBROUTINE wrt_size()

!=======================================================================
! PURPOSE: Compute the molecular weight and elemental composition for
! the species in the dictionary. Data are are stored in dctmw & dctatom.
!=======================================================================
SUBROUTINE dictelement()
  USE keyparameter, ONLY: mxlfo
  USE dictstackdb, ONLY: nrec,dict,dctmw,dctatom  !  dictionary data
  USE atomtool, ONLY: getatoms, molweight
  IMPLICIT NONE
  
  INTEGER :: i,tc,ih,in,io,ir,is,ifl,ibr,icl
  CHARACTER(LEN=mxlfo) :: chem
  INTEGER :: noelement(SIZE(dctatom,2))
  
  dctmw(:)=1.  ; dctatom(:,:)=0  ; noelement(:)=0

  recloop: DO i=2, nrec
    chem=dict(i)(10:129)

    IF (chem(1:1)/='C'.AND.chem(1:1)/='c'.AND.chem(1:2)/='-O' ) THEN
      IF (chem(1:3)=='#mm' ) THEN
        chem(1:)=chem(4:)
      ELSE IF (chem(1:1)=='#' ) THEN
        chem(1:)=chem(2:)
      ENDIF
    ENDIF

    CALL getatoms(chem,tc,ih,in,io,ir,is,ifl,ibr,icl)
    CALL molweight(chem,dctmw(i))
    dctatom(i,1)=ir  ;  dctatom(i,2)=tc  ;  dctatom(i,3)=ih
    dctatom(i,4)=in  ;  dctatom(i,5)=io  ;  dctatom(i,6)=is
    dctatom(i,7)=ifl ;  dctatom(i,8)=ibr ;  dctatom(i,9)=icl

! overwrite some special cases
    IF (dict(i)(10:12)=='#mm' ) THEN
      CYCLE
    ELSE IF (dict(i)(10:10)=='#' ) THEN
      IF (INDEX(chem,'CATEC')/=0)  THEN
        dctmw(i)=0. ; dctatom(i,:)=noelement(:)
      ELSE IF (chem(1:5)=='CH2OO') THEN
        dctmw(i)=46. ; dctatom(i,2)=1  ! overwrite (because "CH2OOC" is misleading)
      ELSE IF (chem(1:2)=='-O') THEN ! these species are allowed
        CYCLE
      ELSE IF (chem(1:1)/='C') THEN
        dctmw(i)=1. ; dctatom(i,:)=noelement(:)
      ENDIF
    ENDIF

  ENDDO recloop
  
END SUBROUTINE dictelement


!=======================================================================
! PURPOSE: Compute diffusion volume.
! Molecular diffusion volume are estimated using the atomic diffusion 
! volume and method described in "The Properties of Gases and Liquids", 
! 5th Fifth Edition, McGraw-Hill 2004, B.E. Poling, J.M. Prausnitz, 
! J.P. O’Connell.
!=======================================================================
SUBROUTINE diffusion_vol()
  USE keyparameter, ONLY:mxnode,mxlgr,mxlfo,mxlco,tfu1,dirout
  USE stdgrbond, ONLY: grbond
  USE ringtool, ONLY: ring_data
  USE dictstackdb, ONLY: nrec,dict,dctatom  !  dictionary data

  CHARACTER(LEN=mxlgr) :: group(mxnode)
  CHARACTER(LEN=mxlfo) :: chem
  INTEGER :: bond(mxnode,mxnode)
  INTEGER :: i,k,ngr,nvd
  INTEGER :: dbflg, nring
  REAL    :: vdmol

  INTEGER,PARAMETER :: mxirg=6           ! max # of distinct rings 
  INTEGER :: ndrg                        ! # of distinct rings
  INTEGER :: trackrg(mxirg,SIZE(group))  ! (a,:)== track (node #) belonging to ring a
  LOGICAL :: lorgnod(mxirg,SIZE(group))  ! (a,b)== true if node b belongs to ring a

  ! atomic diffusion volumes
  REAL,DIMENSION(4),PARAMETER :: vdatom = &
    (/ 15.9, &     ! 1 => C
       2.31, &     ! 2 => H
       6.11, &     ! 3 => O
       4.54 /)     ! 4 => N
  ! correction for aromatic or heteroyclic ring       
  REAL,PARAMETER :: vdring = -18.3  
      
  OPEN(tfu1,FILE=TRIM(dirout)//'vdiffusion.dat')
  WRITE(tfu1,'(a)') '          ! number of records (species) in the file'
  WRITE(tfu1,'(a)') '! Diffusion volume computed using the Fuller et al. approach,'
  WRITE(tfu1,'(a)') '! as described in Poling et al. book (5th edition, 2004)'

! ----------------------------
! LOOP over the dictionary
! ----------------------------
  nvd=0
! scroll the dictionary
  recloop: DO i=2,nrec
    chem=dict(i)(10:129)
    IF (dctatom(i,2) < 2) CYCLE recloop           ! rm C1 and inorg. 
    IF (dctatom(i,1) /=0) CYCLE recloop           ! rm radicals

! check if '#' can be managed
    IF (chem(1:1)/='C'.AND.chem(1:1)/='c'.AND.chem(1:2)/='-O' ) THEN
      IF (chem(1:3)=='#mm' ) THEN
        chem(1:)=chem(4:)
      ELSE IF (chem(1:1)=='#' ) THEN
        IF (INDEX(chem,'CH2OO')/=0)  CYCLE recloop  ! rm crieggee
        chem(1:)=chem(2:)
      ELSE
        CYCLE recloop                               ! rm unexpected species 
      ENDIF
    ENDIF
    
! compute atomic diffusion volume (see reference in the header)
    nvd = nvd + 1
    vdmol = dctatom(i,2)*vdatom(1) + &  ! # of C
            dctatom(i,3)*vdatom(2) + &  ! # of H
            dctatom(i,4)*vdatom(4) + &  ! # of N
            dctatom(i,5)*vdatom(3)      ! # of O

! check for ring correction
    CALL grbond(chem,group,bond,dbflg,nring)
    IF (nring > 0) THEN

      ! aromatic ring
      IF  (INDEX(chem,'c') /= 0) THEN
        vdmol = vdmol + vdring 

      ! heterocyclic ring (i.e. -O in a ring)
      ELSE IF (INDEX(chem,'-O') /= 0) THEN
        ngr=COUNT(group/=' ')
        grloop: DO k=1,ngr
          IF (group(k)(1:2)=='-O') THEN
            CALL ring_data(k,ngr,bond,group,ndrg,lorgnod,trackrg)
            IF (ndrg > 0) THEN
              vdmol = vdmol + vdring 
              EXIT grloop
            ENDIF
          ENDIF
        ENDDO grloop
      ENDIF

    ENDIF

! write the diffusion volume
    WRITE(tfu1,'(A1,A6,1x,f9.2)') 'G',dict(i)(1:6), vdmol
  ENDDO recloop

! write record and close file
  OPEN(tfu1,POSITION='REWIND')  ;  WRITE(tfu1,'(i8$)') nvd
  OPEN(tfu1,POSITION='APPEND')  ;  WRITE(tfu1,'(a4)') "END " 
  CLOSE(tfu1)
END SUBROUTINE diffusion_vol

! ======================================================================
! PURPOSE: write the dictionary in the corresponding outputs files, i.e.
! the "gas phase species" file and the actual dictionary of species. 
! ======================================================================
SUBROUTINE wrt_dict()
  USE keyparameter, ONLY: dctu,gasu,dirout
  USE dictstackdb, ONLY: nrec,ninorg,inorglst,dict,dctmw,dctatom  !  dictionary data
  IMPLICIT NONE

  REAL :: wmass
  INTEGER :: i,j,i0,ndic

  wmass=0.  ;  i0=0

  OPEN(dctu,FILE=TRIM(dirout)//'dictionary.out')

! write number of records in the dictionary
  ndic=ninorg+nrec-1
  WRITE(dctu,'(i7,a)') ndic, "   ! number of record in the dictionary"

! write inorganic species
  DO i=1,ninorg
    WRITE(dctu,'(a, f5.1,9(i3))') inorglst(i),wmass, (i0,j=1,9)
    WRITE(gasu,'(A1,A6,10X,A4)') "G",inorglst(i),"/1./"
  ENDDO

  recloop: DO i=2, nrec

! write C1 species
    IF (dctatom(i,2) < 2) THEN
      WRITE(gasu,'(A1,A6,10X,A1,F6.1,A1)') "G",dict(i),"/",dctmw(i) ,"/" 
      WRITE(dctu,'(a, f5.1,9(i3))') dict(i),dctmw(i), (dctatom(i,j),j=1,9)
      CYCLE recloop   
    ENDIF
    
! write the species in the dictionary
    WRITE(gasu,'(A1,A6,10X,A1,F6.1,A1)') "G",dict(i),"/",dctmw(i) ,"/"
    WRITE(dctu,'(a,f5.1,9(i3))') dict(i),dctmw(i),(dctatom(i,j),j=1,9)

  ENDDO recloop

  WRITE(dctu,'(a)') 'END '    ;  CLOSE(dctu)
END SUBROUTINE wrt_dict

! ======================================================================
! PURPOSE: sort the peroxy species and write in the corresponding file
! ======================================================================
SUBROUTINE wrt_ro2()
  USE keyparameter, ONLY:mxring,mxnode,mxlgr,mxlfo,mxlco,tfu1,tfu2,tfu3,dirout
  USE keyflag, ONLY:rx_ro2_multiclass
  USE rjtool, ONLY: rjgrm
  USE stdgrbond, ONLY: grbond
  USE dictstackdb, ONLY: nrec,dict,dctatom
  IMPLICIT NONE

! local
  CHARACTER(LEN=mxlgr) :: group(mxnode)
  CHARACTER(LEN=mxlfo) :: chem
  INTEGER :: bond(mxnode,mxnode)
  INTEGER :: i,j,ip,ngr,nsub,subon,icd,ipero
  INTEGER :: dbflg, nring
  CHARACTER(LEN=512) :: filename
  INTEGER, PARAMETER :: ncpero=9
  INTEGER :: npero(ncpero)

  npero(:)=0
  
! open counter's file
  DO i=1,ncpero                                              
    WRITE(filename,*) i
    filename=TRIM(dirout)//'pero'//ADJUSTL(filename)
    filename=filename(1:LEN_TRIM(filename))//'.dat'
    OPEN(40+i,FILE=filename)
    WRITE(40+i,'(a)') '          ! number of records (RO2 species) in the file'
  ENDDO


! scroll the dictionary
  recloop: DO i=2,nrec
    chem=dict(i)(10:129)
    IF (dctatom(i,2) < 2) CYCLE recloop           ! rm C1 and inorg. 
    IF (dctatom(i,1) /=1) CYCLE recloop           ! rm radicals
    
! exclude non peroxy names
    IF ((dict(i)(1:1)/='2').AND.(dict(i)(1:1)/='3')) CYCLE recloop          

! check if '#' can be managed
    IF (chem(1:1)/='C'.AND.chem(1:1)/='c'.AND.chem(1:2)/='-O' ) THEN
      IF (chem(1:3)=='#mm')    THEN ; chem(1:)=chem(4:)
      ELSE IF (chem(1:1)=='#') THEN ; chem(1:)=chem(2:)
      ELSE                          ; CYCLE recloop     ! rm unexpected species 
      ENDIF
    ENDIF

    CALL grbond(chem,group,bond,dbflg,nring)
    
    ngr=COUNT(group/=' ')
    ip=0
    DO j=1,ngr
      IF (INDEX(group(j),'(OO.)')/=0) THEN
        ip=j ; EXIT
      ENDIF      
    ENDDO
    IF (ip==0) STOP "in wrt_ro2, no peroxy found"     

! peroxy group (primary, secondary, tertiary)
    nsub=0 ; subon=0 ; icd=0
    DO j=1,ngr
      IF (bond(ip,j)/=0) THEN
        nsub=nsub+1                                      ! peroxy group (prim, sec, ter)
        IF (INDEX(group(j),'O')/=0)       subon=subon+1  ! class ID: beta O, N
        IF (INDEX(group(j)(1:2),'Cd')/=0) icd=icd+1      ! allylic
        IF (INDEX(group(j)(1:1),'c')/=0)  icd=icd+1      ! aryl
      ENDIF
    ENDDO

! set the peroxy class (to be written in peroxy counting files)
    ipero=0
    IF (rx_ro2_multiclass) THEN
      IF (INDEX(group(ip),'CO(OO.)')/=0)  THEN  ! acyl peroxy
        ipero=9 
	  
      ELSE IF (nsub==3) THEN       ! tertiary peroxy
        IF (subon/=0) THEN  
          IF (icd==0) THEN ; ipero=3 ; ELSE ; ipero=5 ; ENDIF
        ELSE
          ipero=1
        ENDIF
      
      ELSE IF (nsub==2) THEN  ! secondary peroxy
        IF (subon/=0) THEN  
          IF (icd==0) THEN ; ipero=7 ; ELSE ; ipero=8 ; ENDIF
        ELSE
          ipero=4
          IF (chem(1:14)=='CH3CH(OO.)CH3 ') ipero=2  ! special case - overwrite 
        ENDIF
      
      ELSE                    ! primary peroxy 
        IF (subon==0) THEN
          ipero=7
          IF (chem(1:12)=='CH3CH2(OO.) ') ipero=4    ! special case - overwrite 
        ELSE 
          ipero=8
        ENDIF
      ENDIF
    ELSE
      ipero=1
    ENDIF
 
    WRITE(40+ipero,'(a1,a6)') "G",dict(i)(1:6)
    IF (ipero==0) STOP "ipero=0"
    npero(ipero)=npero(ipero)+1

  ENDDO recloop

  IF (.NOT. rx_ro2_multiclass) THEN 
    WRITE(41,'(a7)') "GCH3O2 "
    npero(1)=npero(1)+1
  ENDIF

! close files 
  DO i=1,ncpero
    OPEN(40+i,POSITION='REWIND')  ;  WRITE(40+i,'(i8$)') npero(i)
    OPEN(40+i,POSITION='APPEND')  ;  WRITE(40+i,'(a4)') "END " 
    CLOSE(40+i)
  ENDDO
  
END SUBROUTINE wrt_ro2


! ======================================================================
! PURPOSE: write the maximum yield for the species in the dictionary.
! ======================================================================
SUBROUTINE wrt_mxyield()
  USE keyparameter, ONLY: tfu1,dirout
  USE dictstackdb, ONLY: nrec,dict,dctatom,dbrch
  IMPLICIT NONE

  INTEGER :: i

  OPEN(tfu1, FILE=TRIM(dirout)//'maxyield.dat')
  recloop: DO i=2, nrec
    IF (dctatom(i,1) /= 0) CYCLE recloop   ! rm radicals
    WRITE(tfu1,'(f15.12,2x,a)') dbrch(i),dict(i)
  ENDDO recloop
  CLOSE(tfu1)
END SUBROUTINE wrt_mxyield

! ======================================================================
! PURPOSE: write the vapor pressure and the latent heat of the species
! in the dictionary. Various methods can be used simultaneously, if 
! the flag (provided as input) for the corresponding method is raised.
! ======================================================================
SUBROUTINE wrt_psat(myrflg,nanflg,simflg)
  USE keyparameter, ONLY:mxring,mxnode,mxlgr,mxlfo,mxlco,tfu1,tfu2,tfu3,dirout
  USE rjtool, ONLY: rjgrm
  USE stdgrbond, ONLY: grbond
  USE simpoltool, ONLY: simpol          ! for simpol SAR
  USE nannoolaltool, ONLY:nannoolalprop ! for nannoonlal SAR
  USE myrdaltool, ONLY: myrdalprop      ! for Myrdal & Yalkowsky SAR
  USE dictstackdb, ONLY: nrec,dict,dctatom,dctmw,dctnan,dctsim,dctmyr
  IMPLICIT NONE
  
  INTEGER,INTENT(IN) :: myrflg
  INTEGER,INTENT(IN) :: nanflg 
  INTEGER,INTENT(IN) :: simflg   
     
! local
  CHARACTER(LEN=mxlgr) :: group(mxnode)
  CHARACTER(LEN=mxlfo) :: chem
  INTEGER :: bond(mxnode,mxnode)
  INTEGER :: i
  INTEGER :: dbflg, nring
  REAL    :: molmass
  REAL    :: Tb, log10Psat, psat, latentheat
  INTEGER :: rjg(mxring,2)    ! ring-join group pairs
  REAL,PARAMETER :: temp=298. ! ref temperature for Psat and latent heat
  INTEGER :: nmyr,nnan,nsim   ! number of records

  nnan=0          ; nmyr=0           ;  nsim=0
  dctnan(:,:)=0.  ;  dctmyr(:,:)=0.  ;  dctsim(:,:)=0.

! OPEN FILES
! -----------

! write vapor pressure using myrdal SAR
  IF (myrflg==1) THEN
    OPEN(tfu1,FILE=TRIM(dirout)//'pvap.myr.dat')
    WRITE(tfu1,'(a)') '          ! number of records (species) in the file'
    WRITE(tfu1,'(a)') '! data generated using MYRDAL SAR'
    WRITE(tfu1,'(a)') '! 1st line give the T used to estimate vapor pressure' 
    WRITE(tfu1,'(a)') '! Psat (atm) @Tref (1st) and evaporation heat (kJ.mol-1) (2nd)' 
    WRITE(tfu1,'(a5,f5.1)') 'TREF ', temp 
  ENDIF
    
! write vapor pressure using NANNOOLAL SAR
  IF (nanflg==1) THEN
    OPEN(tfu2,FILE=TRIM(dirout)//'pvap.nan.dat')
    WRITE(tfu2,'(a)') '          ! number of records (species) in the file'
    WRITE(tfu2,'(a)') '! data generated using NANNOOLAL SAR'
    WRITE(tfu2,'(a)') '! 1st line give the T used to estimate vapor pressure' 
    WRITE(tfu2,'(a)') '! Psat (atm) @Tref (1st) and evaporation heat (kJ.mol-1) (2nd)' 
    WRITE(tfu2,'(a5,f5.1)') 'TREF ', temp 
  ENDIF

! write vapor pressure using SIMPOL SAR
  IF (simflg==1) THEN
    OPEN(tfu3,FILE=TRIM(dirout)//'pvap.sim.dat')
    WRITE(tfu3,'(a)') '          ! number of records (species) in the file'
    WRITE(tfu3,'(a)') '! data generated using SIMPOL SAR'
    WRITE(tfu3,'(a)') '! 1st line give the T used to estimate vapor pressure' 
    WRITE(tfu3,'(a)') '! Psat (atm) @Tref (1st) and evaporation heat (kJ.mol-1) (2nd)' 
    WRITE(tfu3,'(a5,f5.1)') 'TREF ', temp 
  ENDIF

! ----------------------------
! LOOP over the dictionary
! ----------------------------
! scroll the dictionary
  recloop: DO i=2,nrec
    chem=dict(i)(10:129)
    IF (dctatom(i,2) < 2) CYCLE recloop           ! rm C1 and inorg. 
    IF (dctatom(i,1) /=0) CYCLE recloop           ! rm radicals
    molmass=dctmw(i)

! check if '#' can be managed
    IF (chem(1:1)/='C'.AND.chem(1:1)/='c'.AND.chem(1:2)/='-O' ) THEN
      IF (chem(1:3)=='#mm' ) THEN
        chem(1:)=chem(4:)
      ELSE IF (chem(1:1)=='#' ) THEN
        IF (INDEX(chem,'CH2OO')/=0)  CYCLE recloop  ! rm crieggee
        chem(1:)=chem(2:)
      ELSE
        CYCLE recloop                               ! rm unexpected species 
      ENDIF
    ENDIF

    CALL grbond(chem,group,bond,dbflg,nring)
    CALL rjgrm(nring,group,rjg)

! vapor pressure & latent heat using Myrdal
! ----------------------------
    IF (myrflg==1) THEN
      CALL myrdalprop(chem,bond,group,nring,rjg,molmass,Tb,log10Psat,latentheat)
      psat=EXP(log10Psat*LOG(10.))
      nmyr=nmyr+1
      WRITE(tfu1,'(A1,A6,2x,1p,E10.3,0p,f8.1)') 'G',dict(i)(1:6),psat,latentheat/1000.
      dctmyr(i,1)=psat  ;  dctmyr(i,2)=latentheat/1000.
    ENDIF
        
! vapor pressure & latent heat using Nanoonal
! ----------------------------
    IF (nanflg==1) THEN
      CALL nannoolalprop(chem,bond,group,nring,rjg,Tb,log10Psat,latentheat)
      psat=EXP(log10Psat*LOG(10.))
      nnan=nnan+1
      WRITE(tfu2,'(A1,A6,2x,1p,E10.3,0p,f8.1)') 'G',dict(i)(1:6),psat,latentheat/1000.
      dctnan(i,1)=psat  ;  dctnan(i,2)=latentheat/1000.
    ENDIF

! vapor pressure & latent heat using Simpol
! ----------------------------
    IF (simflg==1) THEN
      CALL simpol(chem,bond,group,nring,log10Psat,latentheat)
      psat=EXP(log10Psat*LOG(10.))
      nsim=nsim+1
      WRITE(tfu3,'(A1,A6,2x,1p,E10.3,0p,f8.1,1X)') 'G',dict(i)(1:6),psat,latentheat/1000.
      dctsim(i,1)=psat  ;  dctsim(i,2)=latentheat/1000.
    ENDIF
  ENDDO recloop
    
! WRITE RECORD & CLOSE FILE
! ----------------------------
  IF (myrflg==1) THEN
    OPEN(tfu1,POSITION='REWIND')  ;  WRITE(tfu1,'(i8$)') nmyr
    OPEN(tfu1,POSITION='APPEND')  ;  WRITE(tfu1,'(a4)') "END " 
    CLOSE(tfu1)
  ENDIF
  IF (nanflg==1) THEN 
    OPEN(tfu2,POSITION='REWIND')  ;  WRITE(tfu2,'(i8$)') nnan
    OPEN(tfu2,POSITION='APPEND')  ;  WRITE(tfu2,'(a4)') "END " 
    CLOSE(tfu2)
  ENDIF
  IF (simflg==1) THEN 
    OPEN(tfu3,POSITION='REWIND')  ;  WRITE(tfu3,'(i8$)') nsim
    OPEN(tfu3,POSITION='APPEND')  ;  WRITE(tfu3,'(a4)') "END " 
    CLOSE(tfu3)
  ENDIF  
END SUBROUTINE wrt_psat  

! ======================================================================
! PURPOSE: compute and write the Henry's law coefficient for the 
! species in the dictionary.
! ======================================================================
SUBROUTINE wrt_henry()
  USE keyparameter, ONLY: mxlfo,tfu1,dirout
  USE gromhetool, ONLY: gromhe                  ! for gromhe SAR
  USE dictstackdb, ONLY: nrec,dict,dctatom,dcthenry
  IMPLICIT NONE

  CHARACTER(LEN=mxlfo) :: chem
  REAL    :: keff
  INTEGER :: i
  
  dcthenry(:) = 0. 
  
! open file
  OPEN(tfu1,FILE=TRIM(dirout)//'henry_gromhe.dat')
  WRITE(tfu1,'(a)') '! data generated using GROMHE SAR'

! scroll the dictionary
  recloop: DO i=2, nrec
    chem=dict(i)(10:129)
    IF (dctatom(i,1) /= 0) CYCLE recloop  ! remove radicals
    IF (dctatom(i,2) < 2) CYCLE recloop   ! rm C1
    IF (chem(1:1)/='C'.AND.chem(1:1)/='c'.AND.chem(1:2)/='-O' ) THEN
      IF (chem(1:3)=='#mm' ) THEN
        chem(1:)=chem(4:)
      ELSE IF (chem(1:1)=='#' ) THEN
        IF (INDEX(chem,'CH2OO')/=0) CYCLE recloop
        chem(1:)=chem(2:)
      ELSE
        CYCLE recloop                     ! rm unexpected species 
      ENDIF
    ENDIF

! compute the Henry's law coefficient
    CALL gromhe(chem,keff)
    Keff=10.**Keff
    Keff = MIN(Keff,1.0E18)
    dcthenry(i)=keff
    
! write the henry's law coefficient
    WRITE(tfu1,'(A1,A6,2x,1p,E10.3,0p)')  'G',dict(i)(1:6),Keff
  ENDDO recloop
      
  WRITE(tfu1,'(a4)') 'END '  ;  CLOSE(tfu1)

END SUBROUTINE wrt_henry

! ======================================================================
! PURPOSE: compute and write data for the (Wesely) deposition 
! parameterisation for the species in the dictionary. Note: the 
! subroutine next adds the parameters for C1 and inorganic species into
! the henry.dep file.
! ======================================================================
SUBROUTINE wrt_depo()
  USE keyparameter, ONLY: mxlfo,tfu1,dirout,dirgecko
  USE dictstackdb, ONLY: nrec,dict,dctmw,dcthenry
  IMPLICIT NONE

  CHARACTER(LEN=mxlfo) :: chem
  REAL    :: rf, cf
  INTEGER :: i
  CHARACTER(LEN=512)    :: mesg
 
! open file
  OPEN(tfu1, FILE=TRIM(dirout)//'henry.dat')

! scroll the dictionary
  recloop: DO i=2, nrec
    chem=dict(i)(10:129)
    IF (dcthenry(i) /= 0.) THEN
      rf=0.
      IF (INDEX(chem,'(OOH)') /= 0) rf=0.1
      IF (INDEX(chem,'CO(OONO2)') /= 0) rf=0.1
      cf= (dctmw(i)/18.01)**0.5
      WRITE(tfu1,'(1pe7.1,2x,1pe7.1,2x,1pe7.1,2x,a1,a6,a1)') &
                cf,dcthenry(i),rf,'G',dict(i)(1:6),' '
    ENDIF            
  ENDDO recloop
  WRITE(tfu1,'(a4)') "END "  ;  CLOSE(tfu1)
  mesg="cat "//TRIM(dirgecko)//"DATA/henry.dep.clean "//TRIM(dirout)//"henry.dat > "//TRIM(dirout)//"henry.dep"
  CALL SYSTEM(mesg)
END SUBROUTINE wrt_depo

! ======================================================================
! PURPOSE: compute and write the heat of formation (Benson based) for the 
! species in the dictionary.
! ======================================================================
SUBROUTINE wrt_heatf()
  USE keyparameter, ONLY: mxlfo,tfu1,dirout
  USE dictstackdb, ONLY: nrec,dict,dctatom,dctmw
  USE bensontool, ONLY: heat
  IMPLICIT NONE

  CHARACTER(LEN=mxlfo) :: chem
  INTEGER :: i
  REAL    :: dHeatf

! open file
  OPEN(tfu1, FILE=TRIM(dirout)//'dHeatf.dat')

! scroll the dictionary
  recloop: DO i=2, nrec
    IF (dctatom(i,1) /= 0) CYCLE recloop   ! rm radicals
    IF (dctatom(i,2) < 2) CYCLE recloop    ! rm C1
    IF (dctmw(i) ==0. ) CYCLE recloop      ! rm special species
    chem=dict(i)(10:129)

    IF (chem(1:1)/='C'.AND.chem(1:1)/='c'.AND.chem(1:2)/='-O' ) THEN
      IF      (chem(1:3)=='#mm') THEN ; chem(1:)=chem(4:)
      ELSE IF (chem(1:1)=='#'  ) THEN ; chem(1:)=chem(2:)
      ENDIF
    ENDIF
    dHeatf = heat(chem)
    WRITE(tfu1,'(A1,A6,2x,f10.2)') 'G',dict(i)(1:6),dHeatf
  ENDDO recloop
  WRITE(tfu1,'(a4)') 'END '  ;  CLOSE(tfu1)
END SUBROUTINE wrt_heatf


! ======================================================================
! PURPOSE: compute and write Tg for the species in the dictionary.
! ======================================================================
SUBROUTINE wrt_Tg()
  USE keyparameter, ONLY: mxlfo,tfu1,dirout
  USE dictstackdb, ONLY: nrec,dict,dctatom,dctmw
  IMPLICIT NONE
      
  CHARACTER(LEN=mxlfo) :: chem
  INTEGER :: i
  INTEGER :: tca,iha,ina,ioa
  REAL    :: Tg
  INTEGER :: nspe

! OPEN FILES
! -----------
  nspe=0
  OPEN(tfu1, FILE=TRIM(dirout)//'Tg.dat')
  WRITE(tfu1,'(a)') '          ! number of records (species) in the file'
  
! ----------------------------
! LOOP over the dictionary
! ----------------------------
  recloop: DO i=2,nrec
    IF (dctatom(i,1)/= 0) CYCLE recloop   ! rm radicals
    IF (dctatom(i,2) < 2) CYCLE recloop   ! rm C1
    IF (dctmw(i) ==0. )   CYCLE recloop   ! rm special species
    chem=dict(i)(10:129)

    IF (chem(1:1)/='C'.AND.chem(1:1)/='c'.AND.chem(1:2)/='-O' ) THEN
      IF      (chem(1:3)=='#mm') THEN ; chem(1:)=chem(4:)
      ELSE IF (chem(1:1)=='#'  ) THEN ; chem(1:)=chem(2:)
      ENDIF
    ENDIF

    tca=dctatom(i,2)  ;  iha=dctatom(i,3)
    ina=dctatom(i,4)  ;  ioa=dctatom(i,5)

    IF ((ioa==0) .AND. (ina==0) .AND. (iha/=0)) THEN
      Tg= (1.96 + log(REAL(tca)))*61.99 + log(REAL(iha))*(-113.33) +  &
          log(REAL(tca))*log(REAL(iha))*28.74

    ELSE IF ((ina==0) .AND. (iha/=0)) THEN
      Tg= (12.13 + log(REAL(tca)))*10.95 + log(REAL(iha))*(-41.82) +  &
          log(REAL(tca))*log(REAL(iha))*21.61 + log(REAL(ioa))*(118.96) + &
          log(REAL(tca))*log(REAL(ioa))*(-24.38)

    ELSE IF (iha/=0) THEN
      Tg= 551.0146 - log(REAL(tca))*61.8525 + &
          log(REAL(iha))*(-254.1697) + log(REAL(ioa))*(146.0169) + &
          log(REAL(ina))*(136.7473) + &
          log(REAL(tca))*log(REAL(iha))*82.1893  + &
          log(REAL(tca))*log(REAL(ioa))*(-57.9076) + &
          log(REAL(tca))*log(REAL(ina))*(-43.7932)

    ELSE
      Tg=0.
    ENDIF    
    
    nspe=nspe+1
    WRITE(tfu1,'(A1,A6,2x,f7.2)') 'G',dict(i)(1:6),Tg
  ENDDO recloop
   
! CLOSE FILE
! ----------------------------
  OPEN(tfu1,POSITION='REWIND')  ;  WRITE(tfu1,'(i8$)') nspe
  OPEN(tfu1,POSITION='APPEND')  ;  WRITE(tfu1,'(a4)') "END " 
  CLOSE(tfu1)
  
END SUBROUTINE wrt_Tg

! ======================================================================
! PURPOSE: Write size of the mechanism
! ======================================================================
SUBROUTINE wrt_size()
  USE keyparameter, ONLY: tfu1,dirout
  USE dictstackdb, ONLY: nrec,ninorg,nwpspe
  USE rxwrttool, ONLY: nrx_all,nrx_hv,nrx_extra,nrx_fo,nrx_tb,nrx_o2, &
                       nrx_meo2,nrx_isom,nrx_tabcf,nrx_ain,nrx_win,nrx_ro2
  IMPLICIT NONE
  INTEGER :: totspecies

  totspecies=nrec+ninorg+nwpspe
  OPEN(tfu1, FILE=TRIM(dirout)//'size.dum')
  WRITE(tfu1,'(a)') "SIZE     ! as generated + 1"

  WRITE(tfu1,'(i8,a)') totspecies+1, "  ! total number of species"        ! 1
  WRITE(tfu1,'(i8,a)') nrx_all+1,    "  ! total number of reactions"      ! 2
  WRITE(tfu1,'(i8,a)') nrx_hv+1,     "  ! HV reactions"                   ! 3
  WRITE(tfu1,'(i8,a)') nrx_tb+1,     "  ! third body M reactions"         ! 4
  WRITE(tfu1,'(i8,a)') nrx_o2+1,     "  ! O2 reactions"                   ! 5

  WRITE(tfu1,'(i8,a)') nrx_extra+1,  "  ! EXTRA reactions"                ! 6
  WRITE(tfu1,'(i8,a)') nrx_ro2+1,    "  ! RO2/RO2 reactions"              ! 7
  WRITE(tfu1,'(i8,a)') nrx_fo+1,     "  ! fall off reactions"             ! 8
  WRITE(tfu1,'(i8,a)') nrx_isom+1,   "  ! R(O.) isomerization reactions"  ! 9
  WRITE(tfu1,'(i8,a)') nrx_ain+1,    "  ! gas <-> part. equilibrium"      !10

  WRITE(tfu1,'(i8,a)') nrx_win+1,    "  ! gas <-> wall equilibrium"       !11
  WRITE(tfu1,'(i8,a)') nrx_tabcf+1,  "  ! variable stoi. coef. reactions" !13
  CLOSE(tfu1)

END SUBROUTINE wrt_size

END MODULE outtool
