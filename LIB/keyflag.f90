MODULE keyflag
  IMPLICIT NONE

  REAL              :: critvp    ! log(Pvap) below which chemistry is ignored

! various cut-offs
  REAL              :: brcut     ! rx cut off branching ratio below which a reaction pathway is ignored
  REAL              :: yldcut    ! rx cut off if: (rx branching ratio)*(yield) is below threshold
  REAL              :: rxloss    ! send to lost carbon if species yield below threshold

  INTEGER           :: maxgen    ! maximum # of generations allowed
  REAL              :: TK        ! user-set temperature

  REAL,PARAMETER    :: Tref=298. ! reference temperature
  REAL,PARAMETER    :: Mref=2.5E19 ! default 3rd body M (molec cm-3)


! -- CHEMISTRY FLAG - SAR SELECTOR
  INTEGER :: kisom_sar                     ! R(O.) selector for isomerisation reaction rate (1=Atkinson 2007; 2=Vereecken 2009)
  INTEGER :: kdiss_sar                     ! R(O.) selector for decomposition reaction rate (1=Atkinson 2007; 2=Vereecken 2009)
  INTEGER :: kohadd_sar                    ! OH addition to alkene (1=Peeters 1997; 2=Ziemann 2009; 3=Jenkin 2016)
  INTEGER :: pvap_sar                      ! select the SAR for vapor pressure: 1=JR-MY; 2=Nannoolal; 3=SIMPOL-1
						                   
! -- FLAGS                                 
  LOGICAL :: g2pfg                         ! flag to write for gas/particle partitioning "reactions"
  LOGICAL :: g2wfg                         ! flag to write for gas/wall partitioning "reactions" (chamber simulations)
  LOGICAL :: isomerfg                      ! allow isomer substitution (O=no substitution, 1=allow)
  LOGICAL :: enolflg                       ! allow to switch enols to ketones
  LOGICAL :: pamfg                         ! "PAM" simulations: H chemistry could be considered
  LOGICAL :: highnoxfg                     ! only NO chemistry for RO2 considered 
  LOGICAL :: rx_ro_no2                     ! add the RO+NO2 reaction to the RO chemistry
  REAL    :: C_NO2                         ! NO2 concentration (molec.cm-3) used in roselector (rochem) for RO+NO2 reaction
  LOGICAL :: rx_ro2_no2                    ! add the RO2+NO2 reversible reaction to the RO2 chemistry
  LOGICAL :: dhffg                         ! flag to activate DHF formation 
  LOGICAL :: chafg                         ! consider cyclic hemi acetal (CHA) species as non volatile
  LOGICAL :: bimolecrx4criegee             ! allow bimolecular reaction for stabilized criegee (keep H2O)
  LOGICAL :: rx_ro2_multiclass             ! treat RO2+RO2 with the 9 classes (true) of ro2 or only CH3O2 (false)
  LOGICAL :: rx_ro2_oh                     ! add the RO2+OH reaction in RO2 chemistry  
										   
! -- INFOS                                 
  LOGICAL :: sar_only_fg                   ! write kOH, kNO3, kO3 for SAR assessment (rate in database not longer used) - if true, run the reactions for the input species only
  LOGICAL :: wrtsarinfo                    ! write info about SAR (groups, ...) in dedicated files  
  LOGICAL, PARAMETER :: wrtopeinfo=.FALSE. ! write "operator info" in the output files - OLD STUFF

! -- OUTPUT FLAGS
! -------------
  LOGICAL :: screenfg                      ! write additional info on screen during generation
  LOGICAL :: wrtdhf                        ! write output file with formation enthalpies information
  LOGICAL :: wrtpvap                       ! write output file with the vapor pressures
  LOGICAL :: wrthenry                      ! write output file with the Henry's law coefs
  LOGICAL :: wrtref                        ! write the references for the reactions in the mechanism
  LOGICAL :: wrttg                         ! write output file with Tg data
  LOGICAL :: wrtdepo                       ! write output file with data for deposition
  LOGICAL :: wrtkivoci                     ! write output file with kivoci and kjvocj
  LOGICAL :: wrtmaxyield                   ! write output file with maximum yields

END MODULE keyflag
!-----------------------------------
! default environmental parameters for mechanism generation
! => values to use if not supplied by user
!-----------------------------------
SUBROUTINE define_defaults
  USE keyparameter, ONLY: dirgecko, dirout
  USE keyflag, ONLY: maxgen, critvp, TK, &
                     brcut, yldcut, rxloss, &
                     rx_ro_no2, rx_ro2_no2, C_NO2, &
                     rx_ro2_multiclass, &
                     g2pfg,g2wfg,isomerfg,highnoxfg,dhffg,chafg, &
                     enolflg,pamfg,bimolecrx4criegee,rx_ro2_oh,  &
                     pvap_sar,kisom_sar,kdiss_sar,kohadd_sar, &
                     wrtsarinfo,screenfg,wrtdhf,wrtpvap,wrtmaxyield, & 
                     wrthenry,wrtref,wrttg,wrtdepo,wrtkivoci,sar_only_fg

IMPLICIT NONE

! various cut-offs
    maxgen = 20              ! maximum # of generations allowed
    TK     = 298.            ! temperature for which to generate the mechanism
    critvp = -13             ! log(Pvap) below which chemistry is ignored
    brcut  = 0.05            ! rx cut off branching ratio below which a reaction pathway is ignored
    yldcut = 1E-3            ! rx cut off if: (rx branching ratio)*(yield) is below threshold
    rxloss = 1E-10           ! send to lost carbon if species yield below threshold
    C_NO2  = 2.5E10          ! NO2 concentration (molec.cm-3) used in roselector (rochem) for RO+NO2 reaction

! reaction switches
    rx_ro_no2 = .FALSE.      ! add the RO+NO2 reaction to the RO chemistry
    rx_ro2_no2= .FALSE.      ! add the RO2+NO2 reversible reaction to the RO2 chemistry
    rx_ro2_multiclass=.TRUE. ! treat RO2+RO2 with the 9 classes (true) of ro2 or only CH3O2 (false)

! SAR 
    pvap_sar   = 2           ! SAR for vapor pressure: 1=JR-MY; 2=Nannoolal; 3=SIMPOL-1
    kisom_sar  = 2           ! R(O.) selector for isomerisation reaction rate (1=Atkinson 2007; 2=Vereecken 2009)
    kdiss_sar  = 2           ! R(O.) selector for decomposition reaction rate (1=Atkinson 2007; 2=Vereecken 2009)
    kohadd_sar = 3           ! OH addition to alkene (1=Peeters 1997; 2=Ziemann 2009; 3=Jenkin 2016)
    
! flags
    g2pfg      = .TRUE.      ! flag to write for gas/particle partitioning "reactions"
    g2wfg      = .FALSE.     ! flag to write for gas/wall partitioning "reactions" (chamber simulations)
    isomerfg   = .TRUE.      ! allow isomer substitution
    highnoxfg  = .FALSE.     ! only NO chemistry for RO2 considered 
    dhffg      = .FALSE.     ! flag to activate DHF formation 
    chafg      = .FALSE.     ! consider cyclic hemi acetal (CHA) species as non volatile
    pamfg      = .FALSE.     ! "PAM" simulations: H chemistry could be considered
    bimolecrx4criegee = .FALSE. ! allow bimolecular reaction for stabilized criegee (keep H2O)
    enolflg    = .TRUE.      ! allow to switch enols to ketones -------------- WARNING --------- /!\ GECKO-A does NOT treat enol for now /!\
    rx_ro2_oh  = .FALSE.     ! add the RO2+OH reaction in RO2 chemistry  ----- WARNING --------- /!\ if true, leads to ROOOH product for which GECKO has a really simple chemistry so far /!\

! Debug/infos
    wrtsarinfo = .FALSE.     ! write info about SAR (groups, ...) in dedicated files  
    sar_only_fg= .FALSE.     ! write kOH, kNO3, kO3 for SAR assessment (rate in database not longer used) - if true, run the reactions for the input species only

! default directories
    dirgecko  = '../'
    dirout = 'OUT/'

! write output files
    screenfg   = .FALSE.     ! write additional info on screen during generation
    wrtdhf     = .FALSE.     ! write output file with formation enthalpies information
    wrtpvap    = .TRUE.      ! write output file with the vapor pressures
    wrthenry   = .TRUE.      ! write output file with the Henry's law coefs
    wrtref     = .TRUE.      ! write the references for the reactions in the mechanism
    wrttg      = .TRUE.      ! write output file with Tg data
    wrtdepo    = .FALSE.     ! write output file with data fo deposition
    wrtkivoci  = .FALSE.     ! write output file with kivoci and kjvocj
    wrtmaxyield= .FALSE.     ! write output file with maximum yields

END SUBROUTINE define_defaults
