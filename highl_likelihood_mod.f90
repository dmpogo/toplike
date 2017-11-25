MODULE highl_likelihood_mod
  !
  !Program to calculate th likelihood of a topology model
  !wrt the CMB data with various cuts.
  !
  USE Topology_types
  IMPLICIT NONE

  real(DP) :: ampl    ! amplitude available to local routines

  PRIVATE ampl

  CONTAINS

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! This module should contain subroutines that set up likelhood calculation
! on high l Planck data for a given theoretical Cl's.
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

     SUBROUTINE highl_likelihood_init()
     !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     ! Initialization routine should read in necessary theory and data file
     ! and setup plik/whateven software to be ready to give likelihoods 
     !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

     END SUBROUTINE highl_likelihood_init

     FUNCTION lnL_highl(ampl_in)
     real(DP), intent(in)  :: ampl_in
     real(DP)  :: lnL_highl
     !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     ! This function should a wrapper on a call to plik/whatever likelhood
     ! evaluator for a preset theory scaled by ampl_in
     !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

     ! Take theory Cl and multiply by ampl_in (check units) and call
     ! plik/whatever library. Would be good not to have to write into files
     if ( do_highl) then
        lnL_highL = 0.0d0
     else
        lnL_highL = 0.0d0
     endif
     return
     END FUNCTION lnL_highl

END MODULE highl_likelihood_mod
