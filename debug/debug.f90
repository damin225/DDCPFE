include 'Modules.f90'
include 'SetUpCrystalProps.f90'
include 'Utility.f90'
include 'ROMUtility.f90'
include 'IntegrateCrystalEqns.f90'
include 'umat.f90'
!include 'StressSolve.f90'

program debug
	
	use Numtype
	USE DataType
	USE CPVars
	use CPVars_n
	use FileIO
    USE PowerLawBoundPar
    USE Plapar
    use slipgeo
    USE  IterPar
    	
	implicit none
	
	integer        				 :: nstatv, grainid
    real(kind=8), allocatable    :: statev(:)
    REAL(KIND=8) ::  gstress_n (NVEC)
	REAL(KIND=8) ::  gestran_n (NVEC)
	REAL(KIND=8) ::  gkappa_n  (MaxNumSlip)
	REAL(KIND=8) ::  gcrot_n   (NSD, NSD)
	REAL(KIND=8) ::  grrot_n   (NSD, NSD)
	REAL(KIND=8) ::  gcrot0    (NSD, NSD)


	REAL(KIND=8) ::  geqvalues (NEQVA)
	REAL(KIND=8) ::  gmu        (NVEC)  
	REAL(KIND=8) ::  gmu_n        (NVEC)  
	REAL(KIND=8) ::  ggamdot   (MaxNumSlip)
	REAL(KIND=8) ::  grss        (MaxNumSlip)
	REAL(KIND=8) ::  gmudot(NVEC)

	REAL(KIND=8) ::  gdd_for_n   (MaxNumSlip)	
	REAL(KIND=8) ::  gdd_revp_n  (MaxNumSlip)	
	REAL(KIND=8) ::  gdd_revn_n  (MaxNumSlip)
	REAL(KIND=8) ::  gdd_stressf_n    (MaxNumSlip)
	REAL(KIND=8) ::  gdd_rev0_n    (MaxNumSlip)
	REAL(KIND=8) ::  gdd_deb_n     (MaxNumSlip)

	REAL(KIND=8) ::  gdd_for     (MaxNumSlip)	
	REAL(KIND=8) ::  gdd_revp    (MaxNumSlip)	
	REAL(KIND=8) ::  gdd_revn    (MaxNumSlip)
	REAL(KIND=8) ::  gdd_stressf    (MaxNumSlip)
	REAL(KIND=8) ::  gdd_rev0    (MaxNumSlip)
	REAL(KIND=8) ::  gdd_deb     (MaxNumSlip)
	REAL(KIND=8) ::   gstress  (NVEC)
	REAL(KIND=8) ::   gestran  (NVEC)
	REAL(KIND=8) ::   gkappa   (MaxNumSlip)
	
	REAL(KIND=8) ::  dd_for_n      (MaxNumSlip)
	REAL(KIND=8) ::  dd_revp_n     (MaxNumSlip)
	REAL(KIND=8) ::  dd_revn_n     (MaxNumSlip)
	REAL(KIND=8) ::  dd_stressf_n  (MaxNumSlip)
	REAL(KIND=8) ::  dd_rev0_n     (MaxNumSlip)
	REAL(KIND=8) ::  dd_deb_n      (MaxNumSlip)
	
	REAL(KIND=8) ::  dd_for        (MaxNumSlip)
	REAL(KIND=8) ::  dd_revp       (MaxNumSlip)
	REAL(KIND=8) ::  dd_revn       (MaxNumSlip)
	REAL(KIND=8) ::  dd_stressf    (MaxNumSlip)
	REAL(KIND=8) ::  dd_rev0       (MaxNumSlip)
	REAL(KIND=8) ::  dd_deb        (MaxNumSlip)
	
	REAL(KIND=8) ::  dd_all 
	REAL(KIND=8) ::  dd_max
	REAL(KIND=8) ::  dd_his   
	REAL(KIND=8) ::  dd_his_max
	REAL(KIND=8) ::  det_dd_max
	REAL(KIND=8) ::  ave_estran (NVEC)
	REAL(KIND=8) ::  det_dd_max_his

	REAL(KIND=8) ::  stress_n(NVEC), estran_n(NVEC), kappa_n(MaxNumSlip),mu_n(NVEC) 
	REAL(KIND=8) ::  stress(NVEC), estran(NVEC), kappa(MaxNumSlip)
	REAL(KIND=8) ::  s_ij(NSD,NSD), c_ijkl(NVEC,NVEC)
	REAL(KIND=8) ::  gamdot(MaxNumSlip), mu(NVEC)
	REAL(KIND=8) ::  rss(MaxNumSlip), crss(MaxNumSlip) 
	REAL(KIND=8) ::  InnerProductVec, SSKineticEqn
	REAL(KIND=8) ::  mudot(NVEC)
	REAL(KIND=8) ::  eqvalues(NEQVA)	
	REAL(KIND = 8) :: VecM(3,MaxNumSlip)
	REAL(KIND = 8) :: VecS(3,MaxNumSlip)
	REAL(KIND = 8) :: ZTen(3,3,MaxNumSlip),ZVec(6,MaxNumSlip)
	
	INTEGER varsPerPart, ip, dex, id, numslip
	real(kind=8)	:: dtime, argmin, argmax, time, stressf  (MaxNumSlip),stressf_n  (MaxNumSlip),dd_tot  
	INTEGER      	:: kstep, kinc, info, ierr, noel,iterState
	REAL(KIND=8) 	::  stressini(NVEC), dstrain(NVEC), cepmod(6,6)
	REAL(KIND=8) 	:: rhs(6)	
	type (xtalVars) CPV0
    type (xtalVars_n) CPV_n0
    type(POWERLOW_BOUNDPar) PLBP0
!-----------------------------------------------------------------------    
    nstatv = 471
    allocate(statev(nstatv))
    grainid = 5
    dtime = 1.d-3
    stress = 1.d0
    stressini = 0.d0
    kstep = 1
    kinc = 1
    time = 0.d0
    noel = 1
	call SetUpCrystalProps( )
	
!	call ComputeStress(rhs, stress, stressini, estran, rss, mudot,     &
!		   gamdot, CPV_n0%gkappa_n, dtime, PLBP0%argmin, PLBP0%argmax, &
!		   dstrain,SlipG%VecM0(:,:,grainid), SlipG%VecS0, SlipG%ZVEC0, iterP%MaxIterNewt, iterP%tolerNewt, kstep, kinc, ierr, grainid) 
!	call PlasticModuli(cepmod, stress,CPV_n0%gkappa_n, dtime, SlipG%VecM0(:,:,grainid), SlipG%VecS0(:,:,grainid), SlipG%ZVEC0(:,:,grainid), grainid)
 
	call SetUpStateVars(nstatv, statev, grainid) 
!	call RecoverStateVars(statev, nstatv, gstress_n, gestran_n,      &
!           gkappa_n, ggamdot, grss, gmudot,geqvalues, gmu ,gmu_n, 	   &
!           gdd_for_n, gdd_for, gdd_revp_n, gdd_revp, gdd_revn_n,       &
!           gdd_revn,gdd_stressf_n, gdd_stressf, gdd_rev0_n, gdd_rev0,  & 
!           gdd_deb_n, gdd_deb, grainid)
           
!	call FetchCrystalVariablesAtPart(gstress_n, gestran_n, gkappa_n, &
!		   geqvalues,stress_n, estran_n, kappa_n, eqvalues, mu_n, 	   &
!		   gmu_n, dd_for_n, gdd_for_n, dd_revp_n, gdd_revp_n,          &
!		   dd_revn_n, gdd_revn_n, dd_stressf_n, gdd_stressf_n,         &
!		   dd_rev0_n, gdd_rev0_n, dd_deb_n, gdd_deb_n)
!	call ResetCrystalQnts(stress, estran, kappa, statev, eqvalues,   &
!           gamdot, rss, mudot, mu, mu_n, s_ij, c_ijkl, stress_n,       &
!           estran_n, kappa_n, dtime, VecM, VecS, ZVec, dd_for_n,       &
!           dd_for, dd_revp_n, dd_revp, dd_revn_n, dd_revn,             &
!           dd_stressf_n, dd_stressf, dd_rev0_n, dd_rev0, dd_deb_n,     &
!           dd_deb, PLBP0%argmin, PLBP0%argmax, grainid   ) 
!    call SaveCrystalVariablesAtPart(stress, estran, kappa, gamdot,   &
!	       rss, mudot, eqvalues, gstress, gestran, gkappa, gmu, mu,    &
!	       ggamdot, grss, gmudot, geqvalues, dd_for, gdd_for, dd_revp, &
!	       gdd_revp, dd_revn, gdd_revn, dd_stressf, gdd_stressf,       &
!	       dd_rev0, gdd_rev0, dd_deb, gdd_deb)
!	call SaveStateVars(statev, nstatv, gstress, gestran, gkappa,     &
!           ggamdot, grss, gmudot, gmu, geqvalues, gdd_for, gdd_revp,   &
!           gdd_revn, gdd_stressf, gdd_rev0, gdd_deb, dd_all, dd_max,   &
!           dd_his, dd_his_max, det_dd_max, ave_estran, det_dd_max_his, &
!           grainid )

	call IntegrateCrystalEqns2(s_ij, stress, estran, kappa, mu,      &
		   mu_n, eqvalues, gamdot, rss, mudot, stress_n, estran_n,     &
		   kappa_n, dtime, iterP%maxIterstate, iterP%MaxIterNewt,  ierr,          &
		   iterState, iterP%tolerState, iterP%tolerNewt, cepmod, dstrain,  kinc,   &
		   kstep,  PLBP0%argmin, PLBP0%argmax, SlipG%VecM0(:,:,grainid), SlipG%VecS0(:,:,grainid), SlipG%ZVEC0(:,:,grainid), dd_for_n, dd_for,  &
		   dd_revp_n, dd_revp, dd_revn_n, dd_revn, stressf_n, stressf, &
		   dd_rev0_n, dd_rev0, dd_deb_n, dd_deb, time, dd_tot, dd_max, &
		   noel, grainid) 
		    
end program debug
