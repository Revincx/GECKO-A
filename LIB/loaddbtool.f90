MODULE loaddbtool
IMPLICIT NONE
CONTAINS

SUBROUTINE loaddb()
  USE simpoltool, ONLY: load_simpoldat  ! for simpol SAR
  USE nannoolaltool, ONLY:load_nandat   ! for nannoonlal SAR
  USE bensontool, ONLY: rdbenson        ! Benson group
  USE spsptool, ONLY: rdspsp, rdoutgene ! for special species mechanism (i.e. #species)
  USE rdkratetool, ONLY: rdkoxfiles     ! to read VOC+Ox files
  USE rdmecatool, ONLY: rdmeca          ! read known mechanisms reactions
  USE rdchemprop, ONLY: rdhenry, rdhydrat
  USE keyparameter, ONLY: dirgecko
  IMPLICIT NONE
  CHARACTER(LEN=100) :: filename
  
! read fixed names for the most common species
  WRITE (6,*) '  ...reading species whose names are fixed'
  CALL rdfixnam()

! Read the rate constants (OH, NO3, O3) in the databases
  CALL rdkoxfiles()
      
! read  prescribed reactions (OH, NO3, O3, RO2, RO, RCOO2)
  CALL rdmeca()

! read dictionary and reactions for the special species (#xxx) 
  WRITE (6,*) '  ...reading special species'
  filename=TRIM(dirgecko)//'DATA/dic_ncar_furan.dat'
  CALL rdspsp(filename)
  WRITE (6,*) '  ...reading reactions for "special species" ...'
  filename=TRIM(dirgecko)//'DATA/mch_ncar_furan.dat'
  CALL rdoutgene(filename)

! read bsongrp group names and their heat of formation :
  WRITE (6,*) '  ...reading benson groups'
  CALL rdbenson()

! load data for vapor pressure SAR
  CALL load_simpoldat()
  CALL load_nandat()

! read known hydration constants & Henry's law constants
  WRITE (6,*) ' ...reading known hydration constants ...'
  CALL rdhydrat()
  WRITE (6,*) ' ...reading known Henrys law constants ...'
  CALL rdhenry()
  
END SUBROUTINE loaddb

END MODULE loaddbtool
