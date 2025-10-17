!=======================================================================
! GECKO-A C-compatible wrapper module
! 
! PURPOSE: Provide C-compatible interfaces for Python bindings
!
! This module wraps the major functions from GECKO-A to be callable
! from Python using ctypes.
!=======================================================================
MODULE gecko_wrapper
  USE iso_c_binding
  IMPLICIT NONE
  
CONTAINS

!=======================================================================
! Initialize default parameters and flags
!=======================================================================
SUBROUTINE c_define_defaults() BIND(C, name="gecko_define_defaults")
  CALL define_defaults()
END SUBROUTINE c_define_defaults

!=======================================================================
! Read namelist file
!=======================================================================
SUBROUTINE c_rd_nml() BIND(C, name="gecko_rd_nml")
  USE rdnml, ONLY: rd_nml
  CALL rd_nml()
END SUBROUTINE c_rd_nml

!=======================================================================
! Initialize dictionaries and stacks
!=======================================================================
SUBROUTINE c_initdictstack() BIND(C, name="gecko_initdictstack")
  USE maintool, ONLY: initdictstack
  CALL initdictstack()
END SUBROUTINE c_initdictstack

!=======================================================================
! Write log information
!=======================================================================
SUBROUTINE c_wrtlog() BIND(C, name="gecko_wrtlog")
  USE logtool, ONLY: wrtlog
  CALL wrtlog()
END SUBROUTINE c_wrtlog

!=======================================================================
! Read reaction info/references
!=======================================================================
SUBROUTINE c_rdreac_info() BIND(C, name="gecko_rdreac_info")
  USE reac_infotool, ONLY: rdreac_info
  CALL rdreac_info()
END SUBROUTINE c_rdreac_info

!=======================================================================
! Load C1 mechanism
!=======================================================================
SUBROUTINE c_loadc1mch() BIND(C, name="gecko_loadc1mch")
  USE loadc1tool, ONLY: loadc1mch
  CALL loadc1mch()
END SUBROUTINE c_loadc1mch

!=======================================================================
! Write koh for C1 species
!=======================================================================
SUBROUTINE c_wrt_kc1() BIND(C, name="gecko_wrt_kc1")
  CALL wrt_kc1()
END SUBROUTINE c_wrt_kc1

!=======================================================================
! Load database
!=======================================================================
SUBROUTINE c_loaddb() BIND(C, name="gecko_loaddb")
  USE loaddbtool, ONLY: loaddb
  CALL loaddb()
END SUBROUTINE c_loaddb

!=======================================================================
! Read chemical input from file
! filename: path to cheminput.dat file (null-terminated C string)
! ninp: output - number of primary species
!=======================================================================
SUBROUTINE c_rdchemin(filename, filename_len, ninp) BIND(C, name="gecko_rdchemin")
  USE loadchemin, ONLY: rdchemin
  USE keyparameter, ONLY: mxps, mxlfo
  INTEGER(C_INT), INTENT(IN), VALUE :: filename_len
  CHARACTER(KIND=C_CHAR), INTENT(IN) :: filename(filename_len)
  INTEGER(C_INT), INTENT(OUT) :: ninp
  
  CHARACTER(LEN=256) :: f_filename
  INTEGER :: ninp_f
  CHARACTER(LEN=mxlfo) :: input(mxps)
  INTEGER :: i
  
  ! Convert C string to Fortran string
  f_filename = ' '
  DO i = 1, filename_len
    IF (filename(i) == C_NULL_CHAR) EXIT
    f_filename(i:i) = filename(i)
  ENDDO
  
  CALL rdchemin(f_filename, ninp_f, input)
  ninp = ninp_f
END SUBROUTINE c_rdchemin

!=======================================================================
! Check parenthesis in chemical formula
!=======================================================================
SUBROUTINE c_check_parenthesis(chem, chem_len) BIND(C, name="gecko_check_parenthesis")
  USE maintool, ONLY: check_parenthesis
  USE keyparameter, ONLY: mxlfo
  INTEGER(C_INT), INTENT(IN), VALUE :: chem_len
  CHARACTER(KIND=C_CHAR), INTENT(INOUT) :: chem(chem_len)
  
  CHARACTER(LEN=mxlfo) :: f_chem
  INTEGER :: i
  
  ! Convert C string to Fortran string
  f_chem = ' '
  DO i = 1, MIN(chem_len, mxlfo)
    IF (chem(i) == C_NULL_CHAR) EXIT
    f_chem(i:i) = chem(i)
  ENDDO
  
  CALL check_parenthesis(f_chem)
  
  ! Convert back to C string
  DO i = 1, MIN(chem_len, mxlfo)
    chem(i) = f_chem(i:i)
  ENDDO
END SUBROUTINE c_check_parenthesis

!=======================================================================
! Process input chemical species
!=======================================================================
SUBROUTINE c_in1chm(chem, chem_len, idnam, idnam_len) BIND(C, name="gecko_in1chm")
  USE loadchemin, ONLY: in1chm
  USE keyparameter, ONLY: mxlfo, mxlco
  INTEGER(C_INT), INTENT(IN), VALUE :: chem_len, idnam_len
  CHARACTER(KIND=C_CHAR), INTENT(IN) :: chem(chem_len)
  CHARACTER(KIND=C_CHAR), INTENT(OUT) :: idnam(idnam_len)
  
  CHARACTER(LEN=mxlfo) :: f_chem
  CHARACTER(LEN=mxlco) :: f_idnam
  INTEGER :: i
  
  ! Convert C string to Fortran string
  f_chem = ' '
  DO i = 1, MIN(chem_len, mxlfo)
    IF (chem(i) == C_NULL_CHAR) EXIT
    f_chem(i:i) = chem(i)
  ENDDO
  
  CALL in1chm(f_chem, f_idnam)
  
  ! Convert Fortran string to C string
  DO i = 1, MIN(idnam_len-1, mxlco)
    idnam(i) = f_idnam(i:i)
  ENDDO
  idnam(MIN(idnam_len, mxlco+1)) = C_NULL_CHAR
END SUBROUTINE c_in1chm

!=======================================================================
! Sort name list
!=======================================================================
SUBROUTINE c_sort_namlst() BIND(C, name="gecko_sort_namlst")
  USE sortstring, ONLY: sort_string
  USE dictstackdb, ONLY: namlst, nrec
  CALL sort_string(namlst(1:nrec))
END SUBROUTINE c_sort_namlst

!=======================================================================
! Compute dictionary elements (molar mass, atoms)
!=======================================================================
SUBROUTINE c_dictelement() BIND(C, name="gecko_dictelement")
  USE outtool, ONLY: dictelement
  CALL dictelement()
END SUBROUTINE c_dictelement

!=======================================================================
! Write dictionary
!=======================================================================
SUBROUTINE c_wrt_dict() BIND(C, name="gecko_wrt_dict")
  USE outtool, ONLY: wrt_dict
  CALL wrt_dict()
END SUBROUTINE c_wrt_dict

!=======================================================================
! Write max yields
!=======================================================================
SUBROUTINE c_wrt_mxyield() BIND(C, name="gecko_wrt_mxyield")
  USE outtool, ONLY: wrt_mxyield
  CALL wrt_mxyield()
END SUBROUTINE c_wrt_mxyield

!=======================================================================
! Write RO2 species
!=======================================================================
SUBROUTINE c_wrt_ro2() BIND(C, name="gecko_wrt_ro2")
  USE outtool, ONLY: wrt_ro2
  CALL wrt_ro2()
END SUBROUTINE c_wrt_ro2

!=======================================================================
! Change phase (gas to particle/wall)
!=======================================================================
SUBROUTINE c_changephase() BIND(C, name="gecko_changephase")
  USE masstranstool, ONLY: changephase
  CALL changephase()
END SUBROUTINE c_changephase

!=======================================================================
! Write size information
!=======================================================================
SUBROUTINE c_wrt_size() BIND(C, name="gecko_wrt_size")
  USE outtool, ONLY: wrt_size
  CALL wrt_size()
END SUBROUTINE c_wrt_size

!=======================================================================
! Write vapor pressure
!=======================================================================
SUBROUTINE c_wrt_psat(nmy, nan, nsi) BIND(C, name="gecko_wrt_psat")
  USE outtool, ONLY: wrt_psat
  INTEGER(C_INT), INTENT(IN), VALUE :: nmy, nan, nsi
  CALL wrt_psat(nmy, nan, nsi)
END SUBROUTINE c_wrt_psat

!=======================================================================
! Write Henry's law coefficient
!=======================================================================
SUBROUTINE c_wrt_henry() BIND(C, name="gecko_wrt_henry")
  USE outtool, ONLY: wrt_henry
  CALL wrt_henry()
END SUBROUTINE c_wrt_henry

!=======================================================================
! Write deposition parameters
!=======================================================================
SUBROUTINE c_wrt_depo() BIND(C, name="gecko_wrt_depo")
  USE outtool, ONLY: wrt_depo
  CALL wrt_depo()
END SUBROUTINE c_wrt_depo

!=======================================================================
! Write heat of formation
!=======================================================================
SUBROUTINE c_wrt_heatf() BIND(C, name="gecko_wrt_heatf")
  USE outtool, ONLY: wrt_heatf
  CALL wrt_heatf()
END SUBROUTINE c_wrt_heatf

!=======================================================================
! Write glass transition temperature
!=======================================================================
SUBROUTINE c_wrt_Tg() BIND(C, name="gecko_wrt_Tg")
  USE outtool, ONLY: wrt_Tg
  CALL wrt_Tg()
END SUBROUTINE c_wrt_Tg

!=======================================================================
! Write diffusion volume
!=======================================================================
SUBROUTINE c_diffusion_vol() BIND(C, name="gecko_diffusion_vol")
  USE outtool, ONLY: diffusion_vol
  CALL diffusion_vol()
END SUBROUTINE c_diffusion_vol

!=======================================================================
! Get number of records in dictionary
!=======================================================================
FUNCTION c_get_nrec() BIND(C, name="gecko_get_nrec") RESULT(nrec_out)
  USE dictstackdb, ONLY: nrec
  INTEGER(C_INT) :: nrec_out
  nrec_out = nrec
END FUNCTION c_get_nrec

!=======================================================================
! Get number of VOCs in stack
!=======================================================================
FUNCTION c_get_nhldvoc() BIND(C, name="gecko_get_nhldvoc") RESULT(nhldvoc_out)
  USE dictstackdb, ONLY: nhldvoc
  INTEGER(C_INT) :: nhldvoc_out
  nhldvoc_out = nhldvoc
END FUNCTION c_get_nhldvoc

!=======================================================================
! Get number of radicals in stack
!=======================================================================
FUNCTION c_get_nhldrad() BIND(C, name="gecko_get_nhldrad") RESULT(nhldrad_out)
  USE dictstackdb, ONLY: nhldrad
  INTEGER(C_INT) :: nhldrad_out
  nhldrad_out = nhldrad
END FUNCTION c_get_nhldrad

!=======================================================================
! Get output directory
!=======================================================================
SUBROUTINE c_get_dirout(dirout_out, dirout_len) BIND(C, name="gecko_get_dirout")
  USE keyparameter, ONLY: dirout
  INTEGER(C_INT), INTENT(IN), VALUE :: dirout_len
  CHARACTER(KIND=C_CHAR), INTENT(OUT) :: dirout_out(dirout_len)
  INTEGER :: i, dlen
  
  dlen = LEN_TRIM(dirout)
  DO i = 1, MIN(dirout_len-1, dlen)
    dirout_out(i) = dirout(i:i)
  ENDDO
  dirout_out(MIN(dirout_len, dlen+1)) = C_NULL_CHAR
END SUBROUTINE c_get_dirout

END MODULE gecko_wrapper
