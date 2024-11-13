MODULE tempflag
  IMPLICIT NONE 

! ignore reaction products (e.g. after a given # of generation)
  INTEGER,SAVE :: iflost  ! flag raised if products are not written 
  INTEGER,SAVE :: xxc     ! # of C atoms in the reactant
  INTEGER,SAVE :: xxh     ! # of H atoms in the reactant
  INTEGER,SAVE :: xxn     ! # of N atoms in the reactant
  INTEGER,SAVE :: xxo     ! # of O atoms in the reactant
  INTEGER,SAVE :: xxr     ! # of radical dots in the reactant
  INTEGER,SAVE :: xxs     ! # of S atoms in the reactant
  INTEGER,SAVE :: xxfl    ! # of F atoms in the reactant
  INTEGER,SAVE :: xxbr    ! # of Br atoms in the reactant
  INTEGER,SAVE :: xxcl    ! # of Cl atoms in the reactant

END MODULE tempflag
