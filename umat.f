! Dislocation-density based CPFE model
! Copyrights: Damin Xia
INCLUDE 'Modules.f90' 
INCLUDE 'uexternaldb.f90'
INCLUDE 'Utility.f90'
INCLUDE 'StressSolve.f90'   
INCLUDE 'ROMUtility.f90'

!***********************************************************************
SUBROUTINE UMAT (stress,  statev,  ddsdde,  sse,     spd,    &
                scd,     rpl,     ddsddt,  drplde,  drpldt,  &
                stran,  dstran, time,    dtime,   temp,      &
                dtemp,   predef,  dpred,   cmname,  ndi,     &
                nshr,    ntens,   nstatv,  props,   nprops,  &
                coords,  drot,    pnewdt,  celent,  dfgrd0,  &
                dfgrd1,  noel,    npt,     layer,   kspt,    &
                kstep,   kinc )
!***********************************************************************
    USE NumType
	USE DataType
	USE CPVars
	use CPVars_n
	use FileIO
    USE PowerLawBoundPar
    
    implicit none                        
                     
!   Variables passed into the UMAT sub

    character*80, intent(in)  :: cmname
    integer(kind=ikind), intent(in)  :: ntens, nstatv, nprops,         &
                                        ndi, nshr, noel, npt,          &
                                        kspt, kstep, kinc, layer                                       
    real(kind=rkind) :: sse, spd, scd, rpl, drpldt,                    &
                        dtime, temp, dtemp, pnewdt, celent

!   Dimension arrays passed into the UMAT sub
    real(kind=rkind)        &
    stress(ntens),          &! Cauchy stress (vector form)
    statev(nstatv),         &! State variables
    ddsdde(ntens,ntens),    &! Tangent Stiffness Matrix
    ddsddt(ntens),          &! Change in stress per change in temperature
    drplde(ntens),          &! Change in heat generation per change in strain
    stran(ntens),           &! Strain tensor (vector form)
    dstran(ntens),          &! Strain increment tensor (vector form)
    time(2),                &! Time Step and Total Time
    predef(1),              &! Predefined state vars dependent on field variables
    dpred(1),               &! Change in predefined state variables
    props(nprops),          &! Material properties
    coords(3),              &! Coordinates of Gauss pt. being evaluated
    drot(3,3),              &! Incremental rotation matrix
    dfgrd0(3,3),            &! Deformation gradient at t_n
    dfgrd1(3,3)              ! Deformation gradient at t_(n+1)
    
    integer numel_aba, numqpt_aba    

!   local variables
    type (xtalVars) CPV0
    type (xtalVars_n) CPV_n0
    type(POWERLOW_BOUNDPar) PLBP0
    
    integer     :: grainid
!-----------------------------------------------------------------------

!   initialize variables
    CPV0=CPV
    CPV_n0= CPV_n
    PLBP0=PLBP 

!   The first two Props record the total element number and number of integration points in each element. 
!   NOT USING!
	numel_aba  = nint (props(1))
	numqpt_aba = nint (props(2))    

!   step 1: determine grainID
    if ((index(cmname, 'ELSET') + index(cmname, 'elset')) <1  .or. &
        (index(cmname, '_DAMAGE')+index(cmname, '_damage') ) <1 )  &
        call RunTimeError(FILE_O, 'Make sure your material name follows the Elset_grainid_DAMAGE')
    read(cmname(6:(index(cmname, '_')-1)), *) grainid

!   Step 2: set up initial state at the first increment
	IF ( time(2)==0.d0  )  THEN 
	    CALL SetUpStateVars (nstatv, statev, grainid)
    ENDIF  

!   Step 3: solve for stress
    CALL StressSolve(  stress,  ddsdde,stran,           & 
                    dstran, statev, nstatv,             &
                    time, dtime, kinc, kstep,           &    
                    pnewdt, PLBP0%argmin, PLBP0%argmax, &
                    CPV0%geqvalues,                     &                       
                    CPV_n0%gstress_n,                   &
                    CPV_n0%gestran_n,                   &
                    CPV_n0%gkappa_n(:,grainid),         &
                    CPV0%gmu,                           &
                    CPV_n0%gmu_n,                       &
                    CPV_n0%gdd_for_n,                   &
                    CPV_n0%gdd_revp_n,                  &
                    CPV_n0%gdd_revn_n,                  &
                    CPV_n0%gdd_stressf_n,               &
                    CPV_n0%gdd_rev0_n,                  &
                    CPV_n0%gdd_deb_n,                   &
                    noel,                               &    
			        grainid                              )   

end subroutine



!*************************************************
SUBROUTINE SetUpStateVars(nstatv, statev, grainid) 
!*************************************************

    USE NumType
    USE PlaPar
    USE OriPar 
    USE FileIO
    use WorkDir
    use NbData
    
	IMPLICIT NONE

    integer         :: nstatv, grainid
    real(kind=8)    :: statev(nstatv)
    integer         :: varsPerPart, dex, numslip
!-----------------------------------------------------------------------

!   By default, the nstatv is equal to the one with more parameters
    varsPerPart =  2*NVEC          & ! stress, estran                1-12
                  + PhSlip(2)      & ! kappa                         13-60 
                  + NEQVA/2        & ! VonMises, Shearate, gamtot    61-63              
                  + PhSlip(2)      & ! gamdot                        64-111              
                  + NVEC           & ! ep                            112-117
                  + PhSlip(2)      & ! rss                           118-165
                  + NVEC           & ! epdot                         166-171
                  + 3*PhSlip(2)    & ! dd_for, dd_revp, dd_revn      172-315
                  + PhSlip(2)      & ! stressf                       316-363
		          + PhSlip(2)      & ! dd_rev0                       364-411
		          + PhSlip(2)        ! dd_deb                        412-459

    if (nstatv .ne. varsPerPart+12)   &
        call RunTimeError(FILE_O, 'check the number of state variables')
        

!  the 2+3+7 is the sequentially defined as:
!  dd in each grain: 1
!  dd_max in each part: 1
!  dd_history_max in each grain: 1
!  dd_history_max in each part: 1
!  dd_det_max_DD in each part: 1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               
!  average plastic strain in each part : 6
!  dd_history_det_max_DD in each part: 1


!   set up state variables
	dex=0

!   stress
    dex = dex + 1           
    statev(dex:(dex+NVEC-1))=pzero

    if (Plap%crystalID(grainid) .eq. kHCP) then
        numslip=PhSlip(1)
    else if (Plap%crystalID(grainid) .eq. kBCC) then
        numslip=PhSlip(2)
    else 
        call RunTimeError( FILE_O,'Error: crystalID exceed the existing lattice!')
    end if
    if (debug==1) then         
        write(FILE_E,*) '*----Initial Stress in ', grainid, '-th grain----'
        write(FILE_E,'(6F12.5)')  statev(dex:(dex+NVEC-1))
    endif
        
!   estrain
    dex = dex + NVEC                      
    statev(dex:(dex+NVEC-1))=pzero
    if (debug==1) then          
        write(FILE_E,*) '*----Initial estrain in ', grainid, '-th grain----'
        write(FILE_E,'(6F12.5)')  statev(dex:(dex+NVEC-1))  
    endif     

!   kappa
    dex = dex + NVEC                               
    statev(dex:(dex+NumSlip-1))=Plap%kappa0(1:numslip,grainid)       
    if (debug==1) then  
        write(FILE_E,*) '*----Initial strength in ', grainid, '-th grain----'
        write(FILE_E,'(6F12.5)')  statev(dex:(dex+NumSlip-1))
    endif
        
!   VonMises, Shearate, gamtot
    dex = dex + numslip                
    statev(dex:(dex+NEQVA/2-1))=pzero      
    if (debug==1) then  
        write(FILE_E,*) '*----Initial VonMises, Shearate, gamtot in ', grainid, '-th grain----'
        write(FILE_E,'(6F12.5)')  statev(dex:(dex+NEQVA/2-1))
    endif         

!   gamdot        
    dex = dex + NEQVA/2                                                           
    statev(dex:(dex+NumSlip-1))=pzero             
    if (debug==1) then  
        write(FILE_E,*) '*----Initial gamdot in ', grainid, '-th grain----'
        write(FILE_E,'(6F12.5)')  statev(dex:(dex+NumSlip-1))
    endif

!   mu 
    dex = dex + NumSlip                               
    statev(dex:(dex+NVEC-1))=pzero       
    if (debug==1) then  
        write(FILE_E,*) '*----Initial mu in ', grainid, '-th grain----'
        write(FILE_E,'(6F12.5)')  statev(dex:(dex+NVEC-1))
    endif        

!   rss
    dex = dex + NVEC                                
    statev(dex:(dex+Numslip-1))=pzero       
    if (debug==1) then  
        write(FILE_E,*) '*----Initial rss in ', grainid, '-th grain----'
        write(FILE_E,'(6F12.5)')  statev(dex:(dex+NumSlip-1))
    endif	    

!   mudot
    dex = dex + NumSlip                                
    statev(dex:(dex+NVEC-1))=pzero       
    if (debug==1) then  
        write(FILE_E,*) '*----Initial mudot in ', grainid, '-th grain----'
        write(FILE_E,'(6F12.5)')  statev(dex:(dex+NVEC-1))
    endif	    

!   dd_for
    dex = dex + NVEC                      
    statev(dex:(dex+Numslip-1))=1.0d12
    if (debug==1) then  
        write(FILE_E,*) '*----Initial dd_for_n in ', grainid, '-th grain----'
        write(FILE_E,*)   statev(dex:(dex+NumSlip-1))
    endif   
        
!   dd_revp 	
    dex = dex + NumSlip                  
    statev(dex:(dex+Numslip-1))=pzero
    if (debug==1) then  
        write(FILE_E,*) '*----Initial dd_revp_n in ', grainid, '-th grain----'
        write(FILE_E,'(6F12.5)')   statev(dex:(dex+NumSlip-1))
    endif	    

!   dd_revn 
    dex = dex + NumSlip                  
    statev(dex:(dex+Numslip-1))=pzero
    if (debug==1) then  
        write(FILE_E,*) '*----Initial dd_revn_n in ', grainid, '-th grain----'
        write(FILE_E,'(6F12.5)')   statev(dex:(dex+NumSlip-1))
    endif	

!   stressf
    dex = dex + NumSlip                 
    statev(dex:(dex+NumSlip-1))=pone
    if (debug==1) then  
        write(FILE_E,*) '*----Initial stressf in ', grainid, '-th grain----'
        write(FILE_E,'(6F12.5)')   statev(dex:(dex+NumSlip-1))
    endif	

!   dd_rev0 
    dex = dex + NumSlip                  
    statev(dex:(dex+NumSlip-1))=1.0d12	 
    if (debug==1) then  
        write(FILE_E,*) '*----Initial dd_rev0 in ', grainid, '-th grain----'
        write(FILE_E,*)  statev(dex:(dex+NumSlip-1))
    endif

!   dd_deb 
    dex = dex + NumSlip                  
    statev(dex:(dex+NumSlip-1))=1.0d10
    if (debug==1) then  
        write(FILE_E,*) '*----Initial dd_deb in ', grainid, '-th grain----'
        write(FILE_E,*)  statev(dex:(dex+NumSlip-1))
    endif

    statev(nstatv-11)=0.d0                                              ! DD in each grain
    statev(nstatv-10)=0.d0                                              ! dd_history max in each grain 
	statev(nstatv-9)=0.d0                                               ! dd_max among all parts	
	statev(nstatv-8)=0.0d0                                              ! dd_history_max in each part
	statev(nstatv-7)=0.0d0                                              ! dd_det_max_DD in each part
	statev((nstatv-6):(nstatv-1))=0.0d0                                 ! average plastic strain in each part
	statev(nstatv)=0.0d0                                                ! dd_history_det_max_DD in each part
    
	return
END

!***********************************************************************
