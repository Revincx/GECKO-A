MODULE database
  USE keyparameter, ONLY: mxlfo, mxldi, mxpd, mxlco,mxnr,mxcopd
  USE references, ONLY: mxlcod
  IMPLICIT NONE

! Rate constants for VOC+oxidant
! ==============================
  INTEGER,PARAMETER :: mxkdb=1000 ! max # of rate constant in a given database

! VOC+OH database
  INTEGER,SAVE               :: nkohdb             ! # of species in the database
  CHARACTER(LEN=mxlfo),SAVE  :: kohdb_chem(mxkdb)  ! formula of the species 
  REAL,SAVE                  :: kohdb_298(mxkdb)   ! rate constant @ 298 K
  REAL,SAVE                  :: kohdb_arr(mxkdb,3) ! arrhenius parameter for the rate constant
  CHARACTER(LEN=mxlcod),SAVE :: kohdb_com(mxkdb,3) ! comment's code for the rate

! VOC+O3 database
  INTEGER,SAVE               :: nko3db             ! # of species in the database
  CHARACTER(LEN=mxlfo),SAVE  :: ko3db_chem(mxkdb)  ! formula of the species 
  REAL,SAVE                  :: ko3db_298(mxkdb)   ! rate constant @ 298 K
  REAL,SAVE                  :: ko3db_arr(mxkdb,3) ! arrhenius parameter for the rate constant
  CHARACTER(LEN=mxlcod)      :: ko3db_com(mxkdb,3) ! comment's code for the rate

! VOC+NO3 database
  INTEGER,SAVE               :: nkno3db            ! # of species in the database
  CHARACTER(LEN=mxlfo),SAVE  :: kno3db_chem(mxkdb) ! formula of the species 
  REAL,SAVE                  :: kno3db_298(mxkdb)  ! rate constant @ 298 K
  REAL,SAVE                  :: kno3db_arr(mxkdb,3)! arrhenius parameter for the rate constant
  CHARACTER(LEN=mxlcod)      :: kno3db_com(mxkdb,3)! comment's code for the rate

! Known mechanism for VOC+oxidant ("kw" data base)
! ===============================
  INTEGER,PARAMETER :: mxkwr=110  ! max # of known reaction in a given database

! VOC+OH mechanism
  INTEGER,SAVE               :: nkwoh               ! # of species in the database 
  CHARACTER(LEN=mxlfo),SAVE  :: kwoh_rct(mxkwr)     ! formula of the VOC reacting with OH         
  INTEGER,SAVE               :: nkwoh_pd(mxkwr)     ! # of channel involved in mechanism 
  REAL,SAVE                  :: kwoh_yld(mxkwr,mxnr)! yield of each channel
  CHARACTER(LEN=mxlfo),SAVE  :: kwoh_pd(mxkwr,mxnr) ! formula of the main product per channel 
  CHARACTER(LEN=mxlco),SAVE  :: kwoh_copd(mxkwr,mxnr,mxcopd) ! coproducts per channel
  CHARACTER(LEN=mxlcod),SAVE :: kwoh_com(mxkwr,3)   ! comment's code for the rate 

! VOC+O3 mechanism
  INTEGER,SAVE               :: nkwo3               ! # of species in the database 
  CHARACTER(LEN=mxlfo),SAVE  :: kwo3_rct(mxkwr)     ! formula of the VOC reacting with OH         
  INTEGER,SAVE               :: nkwo3_pd(mxkwr)     ! # of channel involved in mechanism 
  REAL,SAVE                  :: kwo3_yld(mxkwr,mxnr)! yield of each channel
  CHARACTER(LEN=mxlfo),SAVE  :: kwo3_pd(mxkwr,mxnr) ! formula of the main product per channel 3
  CHARACTER(LEN=mxlco),SAVE  :: kwo3_copd(mxkwr,mxnr,mxcopd) ! coproducts per channel
  CHARACTER(LEN=mxlcod),SAVE :: kwo3_com(mxkwr,3)   ! comment's code for the rate 

! VOC+NO3 mechanism
  INTEGER,SAVE               :: nkwno3               ! # of species in the database 
  CHARACTER(LEN=mxlfo),SAVE  :: kwno3_rct(mxkwr)     ! formula of the VOC reacting with OH         
  INTEGER,SAVE               :: nkwno3_pd(mxkwr)     ! # of channel involved in mechanism 
  REAL,SAVE                  :: kwno3_yld(mxkwr,mxnr)! yield of each channel
  CHARACTER(LEN=mxlfo),SAVE  :: kwno3_pd(mxkwr,mxnr) ! formula of the main product per channel 3
  CHARACTER(LEN=mxlco),SAVE  :: kwno3_copd(mxkwr,mxnr,mxcopd) ! coproducts per channel
  CHARACTER(LEN=mxlcod),SAVE :: kwno3_com(mxkwr,3)   ! comment's code for the rate 

! Known mechanism for radicals (RO2, RCOO2, RO and criegee)
! ===============================

! RO2 chemistry
  INTEGER,SAVE               :: nkwro2   ! # of entry (reaction) in the database
  REAL,SAVE                  :: kwro2_arrh(mxkwr,3)! arrhenius parameter for the rate constant
  REAL,SAVE                  :: kwro2_stoi(mxkwr,4)! stoi. coef. of the products
  CHARACTER(LEN=mxlfo),SAVE  :: kwro2_rct(mxkwr,2) ! formula of the reactants
  CHARACTER(LEN=mxlfo),SAVE  :: kwro2_prd(mxkwr,4) ! formula of the products
  CHARACTER(LEN=mxlcod),SAVE :: kwro2_com(mxkwr,4) ! comment's code for the ro2 reaction 

! RCOO2 chemistry
  INTEGER,SAVE               :: nkwrco3   ! # of entry (reaction) in the database
  REAL,SAVE                  :: kwrco3_arrh(mxkwr,3)! arrhenius parameter for the rate constant
  REAL,SAVE                  :: kwrco3_stoi(mxkwr,4)! stoi. coef. of the products
  CHARACTER(LEN=mxlfo),SAVE  :: kwrco3_rct(mxkwr,2) ! formula of the reactants
  CHARACTER(LEN=mxlfo),SAVE  :: kwrco3_prd(mxkwr,4) ! formula of the products
  CHARACTER(LEN=mxlcod),SAVE :: kwrco3_com(mxkwr,4) ! comment's code for the rcoo2 reaction 

! RO chemistry
  INTEGER,SAVE               :: nkwro   ! # of entry (reaction) in the database
  REAL,SAVE                  :: kwro_arrh(mxkwr,3)! arrhenius parameter for the rate constant
  REAL,SAVE                  :: kwro_stoi(mxkwr,4)! stoi. coef. of the products
  CHARACTER(LEN=mxlfo),SAVE  :: kwro_rct(mxkwr,2) ! formula of the reactants
  CHARACTER(LEN=mxlfo),SAVE  :: kwro_prd(mxkwr,4) ! formula of the products
  CHARACTER(LEN=mxlcod),SAVE :: kwro_com(mxkwr,4) ! comment's code for the ro reaction 
  
! Criegee chemistry
  INTEGER,SAVE               :: nkwcri   ! # of entry (reaction) in the database
  REAL,SAVE                  :: kwcri_arrh(mxkwr,3)! arrhenius parameter for the rate constant
  REAL,SAVE                  :: kwcri_stoi(mxkwr,4)! stoi. coef. of the products
  CHARACTER(LEN=mxlfo),SAVE  :: kwcri_rct(mxkwr,2) ! formula of the reactants
  CHARACTER(LEN=mxlfo),SAVE  :: kwcri_prd(mxkwr,4) ! formula of the products
  CHARACTER(LEN=mxlcod),SAVE :: kwcri_com(mxkwr,3) ! comment's code for the ro reaction 
  
! photolysis labels & known photolysis reactions
! ===============================
  INTEGER,SAVE               :: njdat           ! # of known photo. reactions 
  CHARACTER(LEN=mxlfo),SAVE  :: jchem(mxkwr)    ! formula of the species being photolyzed
  CHARACTER(LEN=mxlfo),SAVE  :: jprod(mxkwr,2)  ! the two main photodissociation fragments
  INTEGER,SAVE               :: jlabel(mxkwr)   ! ID number of photolytic reaction i
  CHARACTER(LEN=mxlco),SAVE  :: coprodj(mxkwr)  ! additional inorganic coproduct
  CHARACTER(LEN=mxlcod),SAVE :: j_com(mxkwr,3)  ! comment's code for the photolysis reactions
  INTEGER,SAVE               :: nj40            ! # of J40 data in database
  INTEGER,SAVE               :: jlab40(mxkwr)   ! label for which J4O is provided
  REAL,SAVE                  :: j40(mxkwr)      ! J40 values for various labels

! Benson groups database
! ==============================
  INTEGER,PARAMETER :: mxbg=500   ! max # of benson group
  INTEGER,PARAMETER :: lbg=24   ! max length of a benson group (string)
  INTEGER,SAVE      :: nbson    ! # of group in the database
  REAL,SAVE               :: bsonval(mxbg) ! values of Benson group
  CHARACTER(LEN=lbg),SAVE :: bsongrp(mxbg) ! Benson group string
  
! Simpol SAR for vapor pressure
! ==============================
  REAL,SAVE :: p0simpgrp298      ! constant term for psat @ 298K 
  REAL,SAVE :: psimpgrp298(30)   ! grp contribution for psat @ 298K
  REAL,SAVE :: h0simpgrp298      ! constant term for heat of vaporisation  @ 298K
  REAL,SAVE :: hsimpgrp298(30)   ! grp contribution for heat of vaporisation  @ 298K

! Nannoolal SAR for vapor pressure and Tb
! ==============================
  REAL,SAVE :: nanweight(219,2)     ! grp contribution. 1st column is Tb weight, 2nd is Psat  

! Henry's law coefficients (HLC)
! ==============================
  INTEGER, SAVE              :: nhlcdb   ! # number of species in the hlc database
  CHARACTER(LEN=mxlfo),SAVE  :: hlcdb_chem(mxkdb)  ! formula of the species
  REAL,SAVE                  :: hlcdb_dat(mxkdb,3) ! data
  CHARACTER(LEN=mxlcod),SAVE :: hlcdb_com(mxkdb,3) ! comment's code for the data

! Hydration constants (for effective HLC)
! ==============================
  INTEGER,SAVE               :: nkhydb    ! # number of species in the khyd database
  CHARACTER(LEN=mxlfo),SAVE  :: khydb_chem(mxkdb)  ! formula of the species
  REAL,SAVE                  :: khydb_dat(mxkdb,1)   ! data
  CHARACTER(LEN=mxlcod),SAVE :: khydb_com(mxkdb,3) ! comment's code for the data 

! special species list (#species)
! ==============================
  INTEGER,PARAMETER :: mxspsp=3000       ! max # of special species
  INTEGER,SAVE              :: nspsp          ! # of "special" species (e.g. furane) 
  LOGICAL,SAVE              :: lospsp(mxspsp) ! logical, "true" if the species has been used.
  CHARACTER(LEN=mxldi),SAVE :: dictsp(mxspsp) ! dictionary line of the special species

! special species chemistry (oge = Out GEnerator)
! ==============================
  INTEGER, PARAMETER :: mxog=10000  ! max # of entry (reaction) in the oge database     
  INTEGER,SAVE              :: noge             ! # of reactions given (oge = "out generator")
  INTEGER,SAVE              :: ogelab(mxog)     ! label for the reaction (if EXTRA or HV)
  CHARACTER(LEN=mxlfo),SAVE :: ogertve(mxog,3)   ! reactants for reaction i
  CHARACTER(LEN=mxlfo),SAVE :: ogeprod(mxog,mxpd)! products for reaction i
  REAL,SAVE                 :: ogearh(mxog,3)    ! arrehnius coefficients (A, n, Ea)
  REAL,SAVE                 :: ogestoe(mxog,mxpd)! stochiometric coefficients for product j 
  REAL,SAVE                 :: ogeaux(mxog,7)    ! aux. info. for reaction i (e.g. falloff react.)

! fixed names data (fn database)
! ==============================
  INTEGER,PARAMETER :: mxfn=3000  ! max # of species with fixed names (i.e. not set by gecko)
  INTEGER,SAVE              :: nfn          ! # of species with fixed names 
  CHARACTER(LEN=mxlco),SAVE :: namfn(mxfn)  ! names of fixed species
  CHARACTER(LEN=mxlfo),SAVE :: chemfn(mxfn) ! standardized formula of species with fixed name

END MODULE database
