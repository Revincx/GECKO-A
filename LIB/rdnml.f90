!*****************************************************************
! PURPOSE: Read and set the parameters for gecko from the namelit
!   - dirgecko  : directory of GECKO-A
!   - dirout : directory for the output files
!   - maxgen : maximum # of generations allowed
!*****************************************************************
module rdnml

  use keyparameter, only: dirgecko, dirout, nmlu 
  use keyflag, only: critvp, brcut, yldcut, rxloss, maxgen, TK,&
                     rx_ro_no2,rx_ro2_no2,C_NO2, pvap_sar, &
                     g2pfg,g2wfg,isomerfg,highnoxfg,dhffg,chafg, &
                     rx_ro2_multiclass,bimolecrx4criegee

  implicit none
  private 

  character(len=9), parameter :: filein='gecko.nml'
  integer :: ierr

  public rd_nml

contains

! ======================================================================
  subroutine rd_nml

! variables in the namelist in DATA must keep the same sections /!\
    namelist /dir/ dirgecko, dirout
    namelist /thresholds/ critvp, brcut, yldcut, rxloss
    namelist /env_cond/ TK, C_NO2
    namelist /reductions/ maxgen,isomerfg,highnoxfg,rx_ro2_multiclass
    namelist /sar/ pvap_sar
    namelist /process/ g2pfg,g2wfg,dhffg,chafg,rx_ro2_no2, &
                       rx_ro_no2,bimolecrx4criegee

    open(nmlu, file=filein, iostat=ierr)
      if (ierr/=0) then
        print *,'--error--, problem reading namelist file'
        stop
      endif

      read(nmlu,nml=dir, iostat=ierr)
      if (ierr/=0) then
        print *,'--error--, problem reading "dir" arguments in namelist'
        stop
      endif

      read(nmlu,nml=thresholds, iostat=ierr)
      if (ierr/=0) then
        print *,'--error--, problem reading "thresholds" arguments in namelist'
        stop
      endif

      read(nmlu,nml=env_cond, iostat=ierr)
      if (ierr/=0) then
        print *,'--error--, problem reading "env_cond" arguments in namelist'
        stop
      endif

      read(nmlu,nml=reductions, iostat=ierr)
      if (ierr/=0) then
        print *,'--error--, problem reading "reductions" arguments in namelist'
        stop
      endif

      read(nmlu,nml=sar, iostat=ierr)
      if (ierr/=0) then
        print *,'--error--, problem reading "sar" arguments in namelist'
        stop
      endif

      read(nmlu,nml=process, iostat=ierr)
      if (ierr/=0) then
        print *,'--error--, problem reading "process" arguments in namelist'
        stop
      endif
    close(nmlu)

  end subroutine rd_nml

end module rdnml
