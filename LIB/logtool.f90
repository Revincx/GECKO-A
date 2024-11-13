MODULE logtool
IMPLICIT NONE
CONTAINS
! ======================================================================
! Purpose: write operating conditions and flags for memo to a log file.
! ======================================================================
SUBROUTINE wrtlog()
  USE keyparameter, ONLY: logu
  USE keyflag
  IMPLICIT NONE
  
! write mechanism parameters to stdout, for logging
  WRITE(logu,*) ' ---------------------------- '
  WRITE(logu,*) ' GECKO VERSION=v0.0.0 '
  WRITE(logu,*) ' ---------------------------- '
  WRITE(logu,*) ' maxgen = ', maxgen
  WRITE(logu,'(a54,f6.1)')    '  mechanism temperature (K) used for branching ratio: ', TK
  WRITE(logu,'(a29,f6.1)')    '  reference temperature (K): ', Tref
  WRITE(logu,'(a35,1pe10.2)') '  default 3rd body M (molec cm-3): ', Mref
  WRITE(logu,*) '          '

! -----

  WRITE(logu,*) ' ---------------------------- '
  WRITE(logu,*) '    ----  MECHANISM  ---- '
  WRITE(logu,*) ' ---------------------------- '

  IF   (g2pfg) THEN; WRITE(logu,*) ' gas/particle mass transfer: included'
  ELSE             ; WRITE(logu,*) ' gas/particle mass transfer: not included'
  ENDIF

  IF   (g2wfg) THEN; WRITE(logu,*) ' gas/wall mass transfer: included'
  ELSE             ; WRITE(logu,*) ' gas/wall mass transfer: not included'
  ENDIF

  IF (dhffg) THEN; WRITE(logu,*) ' DHF isomerisation: yes'
  ELSE           ; WRITE(logu,*) ' DHF isomerisation: no'
  ENDIF
  WRITE(logu,*) '          '

! ----

  WRITE(logu,*) ' ---------------------------- '
  WRITE(logu,*) '    --------- SAR -------- '
  WRITE(logu,*) ' ---------------------------- '
  IF     (pvap_sar==1) THEN; WRITE(logu,*) ' vapor pressure scheme=M&Y '
  ELSEIF (pvap_sar==2) THEN; WRITE(logu,*) ' vapor pressure scheme=Nannoolal 2008'
  ELSEIF (pvap_sar==3) THEN; WRITE(logu,*) ' vapor pressure scheme=SIMPOL-1 '
  ENDIF  

  IF     (kohadd_sar==1) THEN; WRITE(logu,*) ' OH addition SAR =Peeters 1997 '
  ELSEIF (kohadd_sar==2) THEN; WRITE(logu,*) ' OH addition SAR =Ziemann 2009 '
  ELSEIF (kohadd_sar==3) THEN; WRITE(logu,*) ' OH addition SAR =Jenkin 2018  '
  ENDIF

  IF     (kisom_sar==1) THEN; WRITE(logu,*) ' Alkoxy isomerisation: Atkinson 2007 '
  ELSEIF (kisom_sar==2) THEN; WRITE(logu,*) ' Alkoxy isomerisation: Vereecken 2009 '
  ENDIF
  IF     (kdiss_sar==1) THEN; WRITE(logu,*) ' Alkoxy dissociation: Atkinson 2007 '
  ELSEIF (kdiss_sar==2) THEN; WRITE(logu,*) ' Alkoxy dissociation: Vereecken 2009 '
  ENDIF

  WRITE(logu,*) '          '
  WRITE(logu,*) ' ---------------------------- '
  WRITE(logu,*) '   -- OPTIONAL CHEMISTRY --   '
  WRITE(logu,*) ' ---------------------------- '
  IF (rx_ro2_oh)  WRITE(logu,*) ' RO2 + OH reactions activated'
  IF (rx_ro2_oh)  WRITE(logu,*) '    >> Warning : chemistry for (OOOH) is simple here and needs to be improved'
  IF (rx_ro2_no2) WRITE(logu,*) ' RO2 + NO2 reactions activated'
  IF (rx_ro2_no2) WRITE(logu,*) '    >> Warning : species bearing (OONO2) only react by ROONO2 -> RO2 + NO2 decomposition'

  IF (rx_ro_no2) WRITE(logu,*) ' RO + NO2 -> R(ONO2) reactions activated'
  IF (rx_ro_no2) WRITE(logu,'(a44,1pe10.2,a55)') '    >> Warning : rate calculated with [NO2]=',C_NO2, &
                              ' molec.cm-3 to eliminate minor pathways of RO reactions'
  IF (bimolecrx4criegee)  WRITE(logu,*) ' Bimolecular reaction for stabilized criegee allowed (keep H2O)'
  WRITE(logu,*) '          '


! ----

  WRITE(logu,*) ' ---------------------------- '
  WRITE(logu,*) '    ----- Reductions ----- '
  WRITE(logu,*) ' ---------------------------- '
  WRITE(logu,*) ' critical vapor pressure (atm) =', critvp
  IF (isomerfg) THEN; WRITE(logu,*) ' isomerisation allowed '
  ELSE                 ; WRITE(logu,*) ' no isomerisation '
  ENDIF
  WRITE(logu,*) ' Cut off branching ratio below which a reaction pathway is ignored: ', brcut
  WRITE(logu,*) ' high-NOx flag:',highnoxfg
  WRITE(logu,*) ' All classes of RO2 in RO2+RO2: ', rx_ro2_multiclass
  WRITE(logu,*) '          '
  WRITE(logu,*) ' ---------------------------- '

END SUBROUTINE wrtlog

END MODULE logtool
