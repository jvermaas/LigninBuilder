!------------Module containing parameters for PDI generation---------
!------------Use in conjunction with pdigen.f90----------------------

MODULE PDI_PARAMS

  IMPLICIT NONE

  !--------Input parameters------------------------------------------
  INTEGER,PARAMETER::maxsteps = 100000 ! Integration steps
  REAL,PARAMETER::range = 5 ! Max value of Sigma
  INTEGER::maxiteration = 10000 ! Iteration counts
  REAL :: tol = 5 ! Tolerance value for PDI(%). Value between 1-100
  INTEGER :: nch_types = 1 ! Number of chain types

  !---------File names------------------------------------------------
  CHARACTER(len=256) :: pdi_fname, out_fname
  INTEGER, PARAMETER :: pdi_fid = 300, log_fid = 400, out_fid = 500

  !---------------------Arrays----------------------------------------
  INTEGER,ALLOCATABLE,DIMENSION(:)::avg_mw_arr !i/p Avg MW of chain types
  INTEGER,ALLOCATABLE,DIMENSION(:)::nchain_list !i/p list of num chains
  INTEGER,ALLOCATABLE,DIMENSION(:)::MolWt_arr ! Molecular wt o/p list
  REAL,ALLOCATABLE,DIMENSION(:) :: PDI_arr !i/p PDI of each type

  REAL,DIMENSION(1:maxsteps)::Prob_SZ !Probability of SZ distribution 
  REAL,DIMENSION(1:maxsteps)::NormalDist !Normalized distributation 
  REAL,DIMENSION(1:maxsteps)::IntNormalDist !Integrated area for normdist
  REAL,DIMENSION(1:maxsteps)::rat_arr !Ratio of chainlength to meanlength

END MODULE PDI_PARAMS
