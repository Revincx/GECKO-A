MODULE dictstackdb
  USE keyparameter, ONLY: mxldi,mxlco  ! string length in the dictionary

  IMPLICIT NONE

! dictionaries 
  INTEGER,PARAMETER :: mxspe=7000000 ! max # of species allowed in the mechanism
  INTEGER,PARAMETER :: mxc1=3000     ! max # of inorganic and C1 species allowed in the mechanism

  INTEGER,SAVE               :: nrec            ! # of species recorded in dict      
  INTEGER,SAVE               :: ninorg          ! # of inorganic species
  INTEGER,SAVE               :: nwpspe          ! # of species recorded in the condensed (wall,part.) phase
  CHARACTER(LEN=mxldi),SAVE  :: dict(mxspe)     ! dictionary line (code+formula+fg)
  CHARACTER(LEN=mxlco),SAVE  :: namlst(mxspe)   ! name (lco=6 characters) of the species used
  CHARACTER(LEN=mxldi),SAVE  :: inorglst(mxc1)  ! list of inorganic species (code+formula+fg)
  REAL,SAVE                  :: dbrch(mxspe)    ! yield attach to a formula in dict
! to be introduced in a "type dict"
  REAL,SAVE                  :: dctmw(mxspe)     ! molecular weight of species in dict
  INTEGER,SAVE               :: dctatom(mxspe,9) ! atoms & radical for species in dict
  REAL,SAVE                  :: dcthenry(mxspe)  ! Henry's law coeff. for species in dict
  REAL,SAVE                  :: dctnan(mxspe,2)  ! Pvap(:,1) & heat(:,2) (Nannoonal) for species in dict
  REAL,SAVE                  :: dctsim(mxspe,2)  ! Pvap(:,1) & heat(:,2) (Simpol) for species in dict
  REAL,SAVE                  :: dctmyr(mxspe,2)  ! Pvap(:,1) & heat(:,2) (Myrdal&Yalkowski) for species in dict

! stack                                         
  INTEGER,PARAMETER          :: mxsvoc=1500000  ! max # of species in the voc stack
  INTEGER,PARAMETER          :: mxsrad=200      ! max # of species in the radical stack
  INTEGER,PARAMETER          :: lenss=132       ! length of a stack string (code+formula+i3+i3)

  INTEGER,SAVE               :: nhldvoc         ! # of (stable) VOC in the stack
  CHARACTER(LEN=lenss),SAVE  :: holdvoc(mxsvoc) ! VOCs in the stack (name[a6]+formula[a120]+stabl[i3]+level[i3])
  INTEGER,SAVE               :: nhldrad         ! # of radical in the stack
  CHARACTER(LEN=lenss),SAVE  :: holdrad(mxsrad) ! radicals in the stack(name[a6]+formula[a120]+stabl[i3]+level[i3])
  INTEGER,SAVE               :: stabl           ! # of generations needed to produce the current species
  INTEGER,SAVE               :: level           ! # of intermediates (rad + stable) needed to produce the current species
  LOGICAL,SAVE               :: lotopstack      ! if true, add species on top of VOC stack (default is false)

! isomer stack
  INTEGER,PARAMETER          :: mxiso=2000      ! max # of isomers for a given molecule 
  INTEGER,PARAMETER          :: mxcri=21        ! max # of criteria used to discriminate isomers
  INTEGER,SAVE               :: diccri(mxspe,mxcri) !

! tetrahydrofuran (CHA) stack (to manage phase partitioning)
  INTEGER,PARAMETER          :: mxcha=10000      ! max # of tetrahydrofuran (CHA) allowed 
  INTEGER,SAVE               :: ncha             ! # of species recorded in dict      
  CHARACTER(LEN=mxlco)       :: chatab(mxcha)    ! idnam of the CHA 

END MODULE dictstackdb
