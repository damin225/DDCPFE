!***************
 module NumType
!***************

    ! This module defines the kind of integer and real numbers.
    ! Every module, subroutine or func must use this module.
    ! To change the precision from double to single,
    ! only this module needs to be changed.
    implicit none

    integer(kind(1)),parameter :: Ikind=kind(1),Rkind=kind(0.0d0)

    real(Rkind), parameter		:: pzero=0.d0,pone=1.d0,ptwo=2.d0,pthree=3.d0
    real(Rkind), parameter		:: pfour=4.d0,pfive=5.d0,psix=6.d0,pseven=7.d0
    real(Rkind), parameter		:: peight=8.d0,pnine=9.d0,pten=10.d0
    real(Rkind), parameter		:: pthird= 0.333333333333333d0, phalf = 0.5d0, ptwothrd= 0.666666666666667d0    
    real(Rkind), parameter      :: sqr23 = 0.816496580927726d0
    real(Rkind), parameter      :: sqr32 = 1.224744871391589d0
    real(Rkind), parameter      :: sqr2  = 1.414213562373095d0
    real(Rkind), parameter      :: sqr3  = 1.732050807568877d0 
    real(Rkind), parameter      :: Ident2nd(9) = (/1.d0, 0.d0, 0.d0,                &
                                                   0.d0, 1.d0, 0.d0,                &
                                                   0.d0, 0.d0, 1.d0/)
!                                              
    real(Rkind), parameter      :: Ident(6) = (/1.d0, 1.d0, 1.d0, 0.d0,  0.d0,  0.d0/)
!
    real(Rkind), parameter      :: Ident4th(6,6)= (/  1.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0,             &
                                                      0.d0, 1.d0, 0.d0, 0.d0, 0.d0, 0.d0,             &
                                                      0.d0, 0.d0, 1.d0, 0.d0, 0.d0, 0.d0,             &
                                                      0.d0, 0.d0, 0.d0, 1.d0, 0.d0, 0.d0,             &
                                                      0.d0, 0.d0, 0.d0, 0.d0, 1.d0, 0.d0,             &
                                                      0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 1.d0/)  
    real(Rkind), parameter      :: HCPmoduli(6,6)= (/0.0854d-10, -0.0363d-10, -0.0173d-10,      0.d0,      0.d0,      0.d0,&
                                                    -0.0363d-10,  0.0854d-10, -0.0173d-10,      0.d0,      0.d0,      0.d0,&
                                                    -0.0173d-10, -0.0173d-10,  0.0693d-10,      0.d0,      0.d0,      0.d0,&
                                                          0.d0,       0.d0,    	  0.d0,   0.2433d-10,      0.d0,      0.d0,&
                                                          0.d0,       0.d0,       0.d0,         0.d0,0.2062d-10,      0.d0,&
                                                          0.d0,       0.d0,       0.d0,         0.d0,	   0.d0, 0.2062d-10/) 
!    
    real(Rkind), parameter      :: factoTwo(6) = (/1.d0, 1.d0, 1.d0, 2.d0,  2.d0,  2.d0/)
    real(Rkind), parameter      :: factoHalf(6) = (/1.d0, 1.d0, 1.d0, 0.5d0, 0.5d0, 0.5d0/) 
!
    real(ikind), parameter      :: debug=0, iterprint=0     
!
    integer(ikind), parameter   :: NSTAV=5, NEQVA=6, NSD=3, Ndim=6   
    integer(ikind), parameter   :: nvec=6
    integer(ikind), parameter   :: kFCC=1, kBCC=2, kHCP=3 
    integer(ikind), parameter   :: XTAL_CONVERGED=0, XTAL_SING_JACOBIAN=1,  &
                                   XTAL_LS_FAILED=2, XTAL_MAX_ITERS_HIT=3 
!
!--ly--------------------------------------------------------------
!-- NPH: number of phases ; NPart: number of parts (here,the part number=grain number)
!-- ParPerPh: part number per each phase
    integer(ikind), parameter   :: MaxSlipType=4,MaxNumSlip=48,NPh=2, NPart=9,ParPerPh(NPh)=(/8, 1/)
    integer(ikind), parameter   :: PhSlip(NPh)=(/30,48/), SlipSysType(NPh)=(/4,1/)
!--ly--------------------------------------------------------------
    integer(ikind), parameter   :: kMISES_n=1, kSHRATE_n=2, kGAMTOT_n=3 
    integer(ikind), parameter   :: kMISES  =4, kSHRATE  =5, kGAMTOT =6
                                   
!
    integer(ikind), parameter   :: kGAMDOT=1, kdGAMdTAU=2, kdGAMdKAPP=3, kdGamDOTdTau=4, kDDGAM=5   
    integer(ikind), parameter   :: kHARD_EXPL=1, kHARD_MIDP=2, kHARD_ANAL=3, kHARD_DD=4
!
!                        
!
!--- A flag to indicate it is the first increment of an analysis of of an restarted analysis
!    
    logical :: FirstIncr                       
!
end module NumType
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
module AbaData
    ! print the results at element NO. numqpt_aba integration point at NO. numel_aba
	integer, parameter  :: numel_aba=-1, numqpt_aba=-1
! 
end
!-----
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
module PardisoVar
!
    USE Numtype
! variables needed for pardiso
!..     Internal solver memory pointer 
    INTEGER(ikind) pt(64)
!
!..     All other variables 
    INTEGER(Ikind) maxfct, mnum, mtype, phase, nrhs, errorpardiso, msglvl, nnz
    INTEGER(Ikind) iparm(64)
    INTEGER(Ikind), ALLOCATABLE :: ia( : ), ja(:) 
    REAL(Rkind)  dparm(64) 
    REAL(Rkind), ALLOCATABLE :: amatrix( : )
    INTEGER(Ikind) idum, solver, ddum

!.. predifined data

    DATA  nrhs /1/, maxfct /1/, mnum /1/
! 
end module PardisoVar
!-----
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
module WorkDir
!
!
!----Define the file path and root name
    character(len=255) :: jobname, outdir
    character(len=255), parameter :: filepre='test'
    character(len=255), parameter :: CoefTensFile = 'CoefTens'
    character(len=255), parameter :: textureFile='Texture'
    character(len=255), parameter :: PhasePart='phase' 
end module WorkDir
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
Module FILEIO
      INTEGER :: FILE_I, FILE_E, FILE_O, TXT_I, CoefTens_I, iter_O, Ph1PartFile, Ph2PartFile, OutF(2)
      parameter (FILE_I=80, FILE_O=81,FILE_E=82,                       &
                  TXT_I=83, CoefTens_I=84, iter_O=85, Ph1PartFile=86,  &
		  Ph2PartFile=87, NbDataFile=88, OutF=(/89, 90/))
END
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!
module timing
  integer c1, c2, cm, cr
  real(kind=8) :: tss, ths, tsj, longest_step,printed_step
  integer :: most_psis,printed_psis
  integer :: most_updates,printed_updates
  logical :: printing

end module
!
!-----------------------------------------------------------------------
module DataType
    use numtype
!-----------------------------------------------------------------------
    type NBG
        integer, ALLOCATABLE :: NBGID(:)      
        integer  NNG
    end type NBG
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
    type POWERLOW_BOUNDPar
        REAL(KIND=8) :: argMin
        REAL(KIND=8) :: argMax         
    end type POWERLOW_BOUNDPar  
!-----------------------------------------------------------------------
    type XtalPar
        REAL(KIND = 8) :: h0,xm,gam0,tausi,taus0,xms, gamss0, crss0
	REAL(KIND = 8) :: kbolt
	REAL(KIND = 8) :: detF		(NPh,MaxSliptype)
	REAL(KIND = 8) :: detV		(NPh,MaxSliptype)
	REAL(KIND = 8) :: ddavg		(NPh,MaxSliptype)
	REAL(KIND = 8) :: vid		(NPh,MaxSliptype)
	REAL(KIND = 8) :: burgers	(NPh,MaxSliptype)
	REAL(KIND = 8) :: sref		(NPh,MaxSliptype)
	REAL(KIND = 8) :: sita		(NPh,MaxSliptype)
	REAL(KIND = 8) :: sita_ref	(NPh,MaxSliptype)
	REAL(KIND = 8) :: s0ref		(NPh,MaxSliptype)
	REAL(KIND = 8) :: k1 		(NPh,MaxSliptype)
	REAL(KIND = 8) :: gref		(NPh,MaxSliptype)
	REAL(KIND = 8) :: Drag		(NPh,MaxSliptype)
	REAL(KIND = 8) :: Href		(NPh,MaxSliptype)
	REAL(KIND = 8) :: qsub		(NPh,MaxSliptype)
	REAL(KIND = 8) :: shmu		(NPh,MaxSliptype)
	REAL(KIND = 8) :: crss0in	(NPh,MaxSliptype)

        !REAL(KIND = 8), ALLOCATABLE :: kappa0(:,:)
        REAL(KIND = 8) :: kappa0(MaxNumSlip,NPart)    
!----CrystalID can represent the phaseID-------------------------------        
        integer(kind=ikind)   CrystalID(NPart)
    end type XtalPar  
!-----------------------------------------------------------------------
    type SlipSys
        !real(kind=rkind), allocatable ::  VecM0(:, :, :), VecS0(:, :, :)
        real(kind=rkind)    ::  VecM0(NSD, MaxNumSlip, NPart), VecS0(NSD, MaxNumSlip, NPart)
        !real(kind=rkind), allocatable ::  ZTen0(:, :, :,:)
        real(kind=rkind)    ::  ZTen0(NSD, NSD, MaxNumSlip,NPart)
        !real(kind=rkind), allocatable ::  ZVec0(:, :, :)
        real(kind=rkind)    ::  ZVec0(NVEC, MaxNumSlip, NPart)
        !real(kind=rkind), allocatable ::  ZZT0(:,:,:,:)
        real(kind=rkind)    ::  ZZT0(NVEC,NVEC,MaxNumSlip,NPart)
    end type SlipSys 
!-----------------------------------------------------------------------
    type OriData
        !real (kind=rkind), allocatable :: euler(:, :)
        real (kind=rkind)   :: euler(NSD, NPart)
        !real (kind=rkind), allocatable :: gcrot0(:, :, :)
        real (kind=rkind)   :: gcrot0(NSD, NSD, NPart)          
    end type OriData 
!-----------------------------------------------------------------------
    type xtalVars
        !real(kind=rkind), allocatable ::  gstress    (:, :)
        real(kind=rkind)    ::  gstress    (NVEC)
        !real(kind=rkind), allocatable ::  gestran    (:, :)
        real(kind=rkind)    ::  gestran    (NVEC)
        !real(kind=rkind), allocatable ::  gkappa     (:, :)
        real(kind=rkind)    ::  gkappa     (MaxNumSlip, NPart)
        !real(kind=rkind), allocatable ::  geqvalues  (:, :)
        real(kind=rkind)    ::  geqvalues  (NEQVA)
        !real(kind=rkind), allocatable ::  ggamdot    (:, :)
        real(kind=rkind)    ::  ggamdot    (NVEC)
        !real(kind=rkind), allocatable ::  gmu        (:, :) 
        real(kind=rkind)    ::  gmu        (NVEC)  

		real(kind=rkind)    ::  gdd_for     (MaxNumSlip)
		real(kind=rkind)    ::  gdd_revp    (MaxNumSlip)
		real(kind=rkind)    ::  gdd_revn    (MaxNumSlip)
		real(kind=rkind)    ::  gdd_stressf (MaxNumSlip)
		real(kind=rkind)    ::  gdd_rev0    (MaxNumSlip)
		real(kind=rkind)    ::  gdd_deb     (MaxNumSlip)

       
    end type xtalVars
!-----------------------------------------------------------------------
    type xtalVars_n
        !real(kind=rkind), allocatable ::  gstress_n    (:, :)
        real(kind=rkind)    ::  gstress_n    (NVEC)
        !real(kind=rkind), allocatable ::  gestran_n    (:, :)
        real(kind=rkind)    ::  gestran_n    (NVEC)
        !real(kind=rkind), allocatable ::  gkappa_n     (:, :)
        real(kind=rkind)    ::  gkappa_n     (MaxNumSlip, NPart)
        !real(kind=rkind), allocatable ::  gmu_n        (:, :)    
        real(kind=rkind)    ::  gmu_n        (NVEC)   

		real(kind=rkind)    ::  gdd_for_n     (MaxNumSlip)
		real(kind=rkind)    ::  gdd_revp_n    (MaxNumSlip)
		real(kind=rkind)    ::  gdd_revn_n    (MaxNumSlip)
		real(kind=rkind)    ::  gdd_stressf_n (MaxNumSlip)
		real(kind=rkind)    ::  gdd_rev0_n    (MaxNumSlip)
		real(kind=rkind)    ::  gdd_deb_n     (MaxNumSlip)
		
               
    end  type xtalVars_n
!-----------------------------------------------------------------------    
    type  IterData
        integer maxIterstate, MaxIterNewt
        real*8  tolerState, tolerNewt 
    end type  IterData
end module DataType 
!-----------------------------------------------------------------------
!
!
module IterPar
    use datatype
    type(IterData) iterP
end module IterPar
!-----------------------------------------------------------------------
module PlaPar
    use datatype
    type(XtalPar) Plap
end module PlaPar
!-----------------------------------------------------------------------
module SlipGeo
    use datatype
    type(SlipSys) SlipG
end module SlipGeo
!-----------------------------------------------------------------------
module OriPar
    use datatype
    type(OriData) OriP
end module OriPar
!-----------------------------------------------------------------------
module CPVars
!
    use datatype
    type(xtalVars) CPV
!
end module CPVars
!-----------------------------------------------------------------------
module CPVars_n
!
    use datatype
    type(xtalVars_n) CPV_n
!
end module CPVars_n
!-----------------------------------------------------------------------
module PowerLawBoundPar
!
    use datatype
    type(POWERLOW_BOUNDPar) PLBP
!
end module PowerLawBoundPar
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------

module NbData

  use numtype
  use Datatype
  type (NBG) Elset_nb(NPart)

end module NbData
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!-->damin: for direct simulation
module elasmod

	use numtype
	real (kind=8) :: Lijkl(NVec,NVec,Npart)
	
end module elasmod
!-----------------------------------------------------------------------
