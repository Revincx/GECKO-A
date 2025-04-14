! Set the list of key parameters for gecko
MODULE keyparameter
  IMPLICIT NONE

  INTEGER, PARAMETER :: mxlco=6    ! maximum length of species names (code)
  INTEGER, PARAMETER :: mxlfo=120  ! maximum length of a formula
  INTEGER, PARAMETER :: mxlfl=15   ! maximum length of functionality list (string)
  INTEGER, PARAMETER :: mxldi=146  ! string length in the dictionary
  INTEGER, PARAMETER :: mxnode=35  ! maximum number of nodes allowed
  INTEGER, PARAMETER :: mxlgr=21   ! maximum length of a string in a group
  INTEGER, PARAMETER :: mxring=4   ! maximum rings allowed
  INTEGER, PARAMETER :: mxcp=99    ! maximum # of copies of formula allowed & max # of "tracks" 
  INTEGER, PARAMETER :: mxhyd=10   ! maximum # of distinct position to add H2O (i.e. # of carbonyls)
  INTEGER, PARAMETER :: mxhiso=5000 ! maximum # of hydrate isomer a molecule can have
  INTEGER, PARAMETER :: mxpd=38    ! maximum # of products per reaction
  INTEGER, PARAMETER :: mxnr=35    ! maximum # of reactions per species
  INTEGER, PARAMETER :: mxcopd=10  ! maximum # of coproducts per generated species
  INTEGER, PARAMETER :: mxnp=4     ! maximum # of product per reaction in the output mechanism
  INTEGER, PARAMETER :: mxps=600   ! maximum # of "primary" species that can be given as input 
  INTEGER, PARAMETER :: mxtrk=5    ! maximum # of distinct Cd track in a molecule
  INTEGER, PARAMETER :: mxlcd=4    ! maximum length of a Cd track (current max is 4 for C=C-C=C)
  INTEGER, PARAMETER :: mxlest=6   ! maximum length of an ester track (current max is 6: CO-O-CO-O-CO-O)
  
! file unit
  INTEGER, PARAMETER :: dctu=7     ! dictionary file unit
  INTEGER, PARAMETER :: prmu=14    ! list of primary species used (for findname)
  INTEGER, PARAMETER :: logu=15    ! log file unit
  INTEGER, PARAMETER :: refu=16    ! mechanism with reference file unit
  INTEGER, PARAMETER :: mecu=17    ! mechanism file unit
  INTEGER, PARAMETER :: scru=18    ! file unit to redirect info from screen to file
  INTEGER, PARAMETER :: kohu=19    ! kivoci file unit (record k_OH for the various species) 
  INTEGER, PARAMETER :: kno3u=59   ! kjvocj file unit (record k_NO3 for the various species) 
  INTEGER, PARAMETER :: waru=20    ! warning file unit 
  INTEGER, PARAMETER :: gasu=21    ! gas phase species file unit
  INTEGER, PARAMETER :: prtu=22    ! particle phase species file unit
  INTEGER, PARAMETER :: walu=23    ! wall phase species file unit
  INTEGER, PARAMETER :: ohu =26    ! rate constant OH file unit (e.g. SAR assessment)
  INTEGER, PARAMETER :: o3u =27    ! rate constant O3 file unit (e.g. SAR assessment)
  INTEGER, PARAMETER :: no3u=28    ! rate constant NO3 file unit (e.g. SAR assessment)
  INTEGER, PARAMETER :: dhfu=29    ! dummy info file for dhf reactions
  INTEGER, PARAMETER :: saru=30    ! SAR info about group ...
  INTEGER, PARAMETER :: tfu1=31    ! temporary file unit 1
  INTEGER, PARAMETER :: tfu2=32    ! temporary file unit 2
  INTEGER, PARAMETER :: tfu3=33    ! temporary file unit 3
  INTEGER, PARAMETER :: tfu4=34    ! temporary file unit 4
  INTEGER, PARAMETER :: mcou=37    ! mechanism file unit with comment!
  INTEGER, PARAMETER :: pfu1=41    ! peroxy file unit 1
  INTEGER, PARAMETER :: pfu2=42    ! peroxy file unit 2
  INTEGER, PARAMETER :: pfu3=43    ! peroxy file unit 3
  INTEGER, PARAMETER :: pfu4=44    ! peroxy file unit 4
  INTEGER, PARAMETER :: pfu5=45    ! peroxy file unit 5
  INTEGER, PARAMETER :: pfu6=46    ! peroxy file unit 6
  INTEGER, PARAMETER :: pfu7=47    ! peroxy file unit 7
  INTEGER, PARAMETER :: pfu8=48    ! peroxy file unit 8
  INTEGER, PARAMETER :: pfu9=49    ! peroxy file unit 9
  INTEGER, PARAMETER :: nmlu=70    ! namelist file unit
  INTEGER, PARAMETER :: rszu=789   ! record stack size unit (information file)

! directory for the input/output files
  CHARACTER(LEN=512) :: dirgecko  
  CHARACTER(LEN=512) :: dirout  

! digit (to rate groups and ring joining char)     
  CHARACTER(LEN=1), DIMENSION(4), PARAMETER :: digit= (/ '1','2','3','4'/)                                                     

! priorities to rate/sort functionalities in groups
  CHARACTER(LEN=70), PARAMETER ::  &
     pri='CHOCHC CH(CH2 CH3O.)OO. OOH ONOF  Br2Br)Cl2Cl)NO2NO)OH)'

! prime numbers 
  INTEGER, DIMENSION(35), PARAMETER :: prim=  &  
    (/ 2,  3,  5,  7, 11, 13, 17, 19, 23, 29, &     
      31, 37, 41, 43, 47, 53, 59, 61, 67, 71, &
      73, 79, 83, 89, 97,101,103,107,109,113, &
      127,131,137,139,149/)


  
END MODULE keyparameter
