INCLUDE 'IntegrateCrystalEqns.f90'
!***************************************
SUBROUTINE StressSolve(                 &
			stress,  ddsdde,strain,     & 
			dstrain, statev, nstatv,    &
			time, dtime, kinc, kstep,   &    
			pnewdt, argmin, argmax,     &
			geqvalues,                  &          
			gstress_n,                  &
			gestran_n,                  &
			gkappa_n,                   &
			gmu,                        &
			gmu_n,                      &
			gdd_for_n,                  &
			gdd_revp_n,                 &
			gdd_revn_n,                 &
			gdd_stressf_n,              &
			gdd_rev0_n,                 &
			gdd_deb_n,                  &
			noel,						&
			grainid)  
!***************************************
   
	USE  PlaPar
    USE  NumType
	USE  FILEIO
	USE  SlipGeo
	USE  IterPar
	USE	 NbData
	IMPLICIT NONE  

	INTEGER      ::  nstatv ,kinc, kstep, iNG, noel
	REAL(KIND=8) ::  argmin, argmax, dtime, pnewdt, time(2)
	REAL(KIND=8) ::  stress(NVEC), ddsdde(NVEC, NVEC), PartStress(NVEC), &
	                 strain(NVEC), dstrain(NVEC), statev(nstatv)
     
	REAL(KIND=8) ::  gstress_n  (NVEC)
	REAL(KIND=8) ::  gstress    (NVEC)	
	REAL(KIND=8) ::  stress_n   (NVEC)		
	REAL(KIND=8) ::  estran     (NVEC)
	REAL(KIND=8) ::  gestran    (NVEC)			
	REAL(KIND=8) ::  estran_n   (NVEC)		
	REAL(KIND=8) ::  gestran_n  (NVEC)
	REAL(KIND=8) ::  geqvalues  (NEQVA)
	REAL(KIND=8) ::  eqvalues   (NEQVA)			
	REAL(KIND=8) ::  gmu        (NVEC)	
	REAL(KIND=8) ::  gmu_n      (NVEC)
	REAL(KIND=8) ::  mu         (NVEC)	
	REAL(KIND=8) ::  mu_n       (NVEC)		
	REAL(KIND=8) ::  gkappa_n   (MaxNumSlip)
	REAL(KIND=8) ::  kappa      (MaxNumSlip)
	REAL(KIND=8) ::  gkappa     (MaxNumSlip)		
	REAL(KIND=8) ::  kappa_n    (MaxNumSlip)	

	REAL(KIND=8) ::  gdd_for_n     (MaxNumSlip)
	REAL(KIND=8) ::  gdd_revp_n    (MaxNumSlip)
	REAL(KIND=8) ::  gdd_revn_n    (MaxNumSlip)
	REAL(KIND=8) ::  gdd_stressf_n (MaxNumSlip)
	REAL(KIND=8) ::  gdd_rev0_n    (MaxNumSlip)
	REAL(KIND=8) ::  gdd_deb_n     (MaxNumSlip)

	REAL(KIND=8) ::  gdd_for       (MaxNumSlip)
	REAL(KIND=8) ::  gdd_revp      (MaxNumSlip)
	REAL(KIND=8) ::  gdd_revn      (MaxNumSlip)
	REAL(KIND=8) ::  gdd_stressf   (MaxNumSlip)
	REAL(KIND=8) ::  gdd_rev0      (MaxNumSlip)
	REAL(KIND=8) ::  gdd_deb       (MaxNumSlip)

	REAL(KIND=8) ::  dd_for        (MaxNumSlip)
	REAL(KIND=8) ::  dd_revp       (MaxNumSlip)
	REAL(KIND=8) ::  dd_revn       (MaxNumSlip)
	REAL(KIND=8) ::  dd_stressf    (MaxNumSlip)
	REAL(KIND=8) ::  dd_rev0       (MaxNumSlip)
	REAL(KIND=8) ::  dd_deb        (MaxNumSlip)

	REAL(KIND=8) ::  dd_for_n      (MaxNumSlip)
	REAL(KIND=8) ::  dd_revp_n     (MaxNumSlip)
	REAL(KIND=8) ::  dd_revn_n     (MaxNumSlip)
	REAL(KIND=8) ::  dd_stressf_n  (MaxNumSlip)
	REAL(KIND=8) ::  dd_rev0_n     (MaxNumSlip)
	REAL(KIND=8) ::  dd_deb_n      (MaxNumSlip)


	REAL(KIND=8) ::  ggamdot    (MaxNumSlip)
	REAL(KIND=8) ::  gamdot     (MaxNumSlip)	
    REAL(KIND=8) ::  grss       (MaxNumSlip)	
	REAL(KIND=8) ::  rss        (MaxNumSlip)
	REAL(KIND=8) ::  gmudot     (NVEC)
	REAL(KIND=8) ::  mudot      (NVEC)
	REAL(KIND=8) ::  dd_all   
    REAL(KIND=8) ::  dd_max
    
!	history maximum dislocation density 
    REAL(KIND=8) ::  dd_his     
    REAL(KIND=8) ::  dd_his_max
	REAL(KIND=8) ::  det_dd, det_dd_tmp
	REAL(KIND=8) ::  det_dd_max	
	REAL(KIND=8) ::  ave_estran (NVEC)
	REAL(KIND=8) ::  det_dd_max_his		
		

	REAL(KIND=8) ::  savg_ij(NSD, NSD) 
	REAL(KIND=8) ::  cavg_ijkl(NVEC, NVEC)
	REAL(KIND=8) ::  c_ijkl(NVEC, NVEC)
	REAL(KIND=8) ::  s_ij(NSD, NSD) 		

	INTEGER statusFlag  

	REAL(KIND=8) :: VecM(NSD,MaxNumSlip)
	REAL(KIND=8) :: VecS(NSD,MaxNumSlip)
	REAL(KIND=8) :: ZVec(NVEC,MaxNumSlip)

	INTEGER      :: grainid,is, IPart1, IPart2, ierr, iterState
	CHARACTER(LEN=30) errinfo

!-----------------------------------------------------------------------      

    if (debug==1) then
		write(FILE_E,*) '====***====     Begin ROMSolve    ====***===='
        Write(FILE_E,*) 'Dstrain  and strain is:'
        Write(FILE_E,'(6e18.8)') dstrain
        Write(FILE_E,'(6e18.8)') strain
        write(FILE_E,1000) time(1),dtime,kinc      
    endif

!	fetch state variables from abaqus vector array
	call   RecoverStateVars(statev, nstatv, gstress_n, gestran_n,      &
           gkappa_n, ggamdot, grss, gmudot,geqvalues, gmu ,gmu_n, 	   &
           gdd_for_n, gdd_for, gdd_revp_n, gdd_revp, gdd_revn_n,       &
           gdd_revn,gdd_stressf_n, gdd_stressf, gdd_rev0_n, gdd_rev0,  & 
           gdd_deb_n, gdd_deb, grainid) 

!	compute state
	statusFlag = XTAL_CONVERGED
	
!	Initialize the streess of each part and system Jacobian
	savg_ij=pzero
	cavg_ijkl=pzero           

!	Initialize some local arrays
    stress=pzero
    PartStress=pzero
    estran=pzero
    kappa=pzero
    mu=pzero

	c_ijkl=pzero
	rss=pzero
	mudot=pzero
	gamdot=pzero

    stress_n=pzero
    kappa_n=pzero
    mu_n=pzero

	dd_for_n=pzero
	dd_revp_n=pzero
	dd_revn_n=pzero
	dd_stressf_n=pzero
	dd_rev0_n=pzero
	dd_deb_n=pzero

    dd_all=pzero
    dd_max=pzero
    
!	Initialize the history dd max and max of dd discrepancy
    dd_his=statev(nstatv-10)! dd_his_max in each part
    dd_his_max=statev(nstatv-8)                                 		! dd_history_max in all parts
    ave_estran(:) = statev((nstatv-6):(nstatv-1))    					! average plastic strain in each part
    det_dd_max_his=statev(nstatv)                                 		! dd_history_det_max_DD in all parts
	
!	Fetch Crystal Variables At Part  
    if (debug==1) then
        Write(FILE_E,*) 'initial estimate stress_n='
        Write(FILE_E,*) gstress_n
    end if

	call   FetchCrystalVariablesAtPart(gstress_n, gestran_n, gkappa_n, &
		   geqvalues,stress_n, estran_n, kappa_n, eqvalues, mu_n, 	   &
		   gmu_n, dd_for_n, gdd_for_n, dd_revp_n, gdd_revp_n,          &
		   dd_revn_n, gdd_revn_n, dd_stressf_n, gdd_stressf_n,         &
		   dd_rev0_n, gdd_rev0_n, dd_deb_n, gdd_deb_n)

    if (debug==1) then
        write(FILE_E,*) 'initial estimate stress_n='
        write(FILE_E,*) stress_n
    end if

	call   IntegrateCrystalEqns2(s_ij, partstress, estran, kappa, mu,      &
		   mu_n, eqvalues, gamdot, rss, mudot, stress_n, estran_n,     &
		   kappa_n, dtime, iterP%maxIterstate, iterP%MaxIterNewt,  ierr,          &
		   iterState, iterP%tolerState, iterP%tolerNewt, c_ijkl, dstrain,  kinc,   &
		   kstep,  argmin, argmax, SlipG%VecM0(:,:,grainid), SlipG%VecS0(:,:,grainid), SlipG%ZVEC0(:,:, grainid),dd_for_n, dd_for,  &
		   dd_revp_n, dd_revp, dd_revn_n, dd_revn, dd_stressf_n, dd_stressf, &
		   dd_rev0_n, dd_rev0, dd_deb_n, dd_deb, time(2), dd_all, dd_max, &
		   noel, grainid)  
	   
    if (ierr==XTAL_CONVERGED) then
        errinfo='Converged'
    elseif (ierr==XTAL_SING_JACOBIAN) then
        errinfo= 'SING_JACOBIAN'
    elseif (ierr==XTAL_MAX_ITERS_HIT) then   
        errinfo= 'MAX_ITERS_HIT'
    else
        errinfo= 'unknown error message, impossible!'  
        call RuntimeError(FILE_O, 'unknown error message, impossible!')
    endif     
        

    if (iterprint==1) write(iter_O, '(3(I6,6x), 6x, A20)'), kstep, kinc, iterstate, errinfo

    if (debug ==1) write(FILE_E,*) 'ierr= ', ierr

	IF (ierr .ne. XTAL_CONVERGED) then
		call WriteMessage(FILE_O, 'Resetting xtal quantities')

		call RecoverStateVars(statev, nstatv, gstress_n, gestran_n,      &
             gkappa_n, ggamdot, grss, gmudot,geqvalues, gmu ,gmu_n, 	 &
             gdd_for_n, gdd_for, gdd_revp_n, gdd_revp, gdd_revn_n,       &
             gdd_revn,gdd_stressf_n, gdd_stressf, gdd_rev0_n, gdd_rev0,  & 
             gdd_deb_n, gdd_deb, grainid)
           
		call FetchCrystalVariablesAtPart(gstress_n, gestran_n, gkappa_n, &
			 geqvalues,stress_n, estran_n, kappa_n, eqvalues, mu_n, 	 &
			 gmu_n, dd_for_n, gdd_for_n, dd_revp_n, gdd_revp_n,          &
			 dd_revn_n, gdd_revn_n, dd_stressf_n, gdd_stressf_n,         &
			 dd_rev0_n, gdd_rev0_n, dd_deb_n, gdd_deb_n)
	
        call ResetCrystalQnts(partstress, estran, kappa, statev, eqvalues,             &
             gamdot, rss, mudot, mu, mu_n, s_ij, c_ijkl, stress_n,                     &
             estran_n, kappa_n, dtime, SlipG%VecM0, SlipG%VecS0, SlipG%ZVEC0, dd_for_n,&
             dd_for, dd_revp_n, dd_revp, dd_revn_n, dd_revn,                           &
             dd_stressf_n, dd_stressf, dd_rev0_n, dd_rev0, dd_deb_n,                   &
             dd_deb, argmin, argmax, grainid       ) 
             
		call SaveCrystalVariablesAtPart(partstress, estran, kappa, gamdot,   &
	         rss, mudot, eqvalues, gstress, gestran, gkappa, gmu, mu,        &
	         ggamdot, grss, gmudot, geqvalues, dd_for, gdd_for, dd_revp,     &
	         gdd_revp, dd_revn, gdd_revn, dd_stressf, gdd_stressf,           &
	         dd_rev0, gdd_rev0, dd_deb, gdd_deb)

		call SaveStateVars(statev, nstatv, gstress, gestran, gkappa,     &
             ggamdot, grss, gmudot, gmu, geqvalues, gdd_for, gdd_revp,   &
             gdd_revn, gdd_stressf, gdd_rev0, gdd_deb, dd_all, dd_max,   &
             dd_his, dd_his_max, det_dd_max, ave_estran, det_dd_max_his, &
             grainid)

		write(FILE_O, *) ' ** Umat did not converged       **'
		write(FILE_O, *) ' ** re-scaling time step by 0.75 **'
		pnewdt = 0.75

		savg_ij=s_ij
	    cavg_ijkl=c_ijkl
		CALL SaveStressModuli(stress, ddsdde, savg_ij, cavg_ijkl)
		
		RETURN        
	ENDIF
		
	call SaveCrystalVariablesAtPart(partstress, estran, kappa, gamdot,   &
	     rss, mudot, eqvalues, gstress, gestran, gkappa, gmu, mu,        &
	     ggamdot, grss, gmudot, geqvalues, dd_for, gdd_for, dd_revp,     &
	     gdd_revp, dd_revn, gdd_revn, dd_stressf, gdd_stressf,           &
	     dd_rev0, gdd_rev0, dd_deb, gdd_deb)
         
    savg_ij=s_ij
	cavg_ijkl=c_ijkl
	
!***********************************************************************
!-->damin: comment this part out, incompatible with direct simulation
	
!	Output the max of dd discrepancy and local history max
!    det_dd_tmp=0.d0
!	if  (dd_all .ge. dd_his) dd_his=dd_all

!    DO iNG = 1, Elset_nb(grainid)%NNG
!		det_dd_tmp(Elset_nb(grainid)%NBGID(iNG))=abs(dd_all(grianid)-dd_all(Elset_nb(grainid)%NBGID(iNG)))
!    ENDDO
!	det_dd=maxval(det_dd_tmp)

!	det_dd_max=maxval(det_dd)

!!	Output the max inside history dd
!    IF  (dd_max .ge. dd_his_max) then
!		dd_his_max=dd_max
!	ENDIF

!	Output the average estrain in parts
	ave_estran=pzero
!	DO ip=1,NPart     
!		ave_estran(1) = ave_estran(1)+gestran(1,ip)*VolFrac(ip)
!		ave_estran(2) = ave_estran(2)+gestran(2,ip)*VolFrac(ip)
!		ave_estran(3) = ave_estran(3)+gestran(3,ip)*VolFrac(ip)
!		ave_estran(4) = ave_estran(4)+gestran(4,ip)*VolFrac(ip)
!		ave_estran(5) = ave_estran(5)+gestran(5,ip)*VolFrac(ip)
!		ave_estran(6) = ave_estran(6)+gestran(6,ip)*VolFrac(ip)
!	ENDDO

!	Output the history max of the max spatial dd discrepancy
!    IF  (det_dd_max .ge. det_dd_max_his) then
!	     det_dd_max_his=det_dd_max
!	ENDIF

!-->damin: comment this part out, incompatible with direct simulation
!***********************************************************************

	call SaveStateVars(statev, nstatv, gstress, gestran, gkappa,     &
         ggamdot, grss, gmudot, gmu, geqvalues, gdd_for, gdd_revp,   &
         gdd_revn, gdd_stressf, gdd_rev0, gdd_deb, dd_all, dd_max,   &
         dd_his, dd_his_max, det_dd_max, ave_estran, det_dd_max_his, &
         grainid) 
                    	
!	stresses and algorithmic moduli in abaqus format
	CALL SaveStressModuli(stress, ddsdde, savg_ij, cavg_ijkl)
	
	if(debug==1 .and. noel==1) write(File_e,*) 'kinc=',kinc
	if(debug==1 .and. noel==1) write(File_e,'(6E20.10)') stress
    if(debug==1 .and. noel==1) write(File_e,'(6E20.10)') strain
    
1000    FORMAT(/'*---time(1)---dtime---kinc---*'/, 7x, F16.10,F16.10,i8)  
  
	RETURN
END


!***********************************************************************
SUBROUTINE RecoverStateVars(statev, nstatv, gstress_n, gestran_n,      &
           gkappa_n, ggamdot, grss, gmudot,geqvalues, gmu ,gmu_n, 	   &
           gdd_for_n, gdd_for, gdd_revp_n, gdd_revp, gdd_revn_n,       &
           gdd_revn,gdd_stressf_n, gdd_stressf, gdd_rev0_n, gdd_rev0,  & 
           gdd_deb_n, gdd_deb, grainid)
!***********************************************************************
      
    USE NumType
	USE FILEIO
	Use Plapar
	IMPLICIT NONE

	INTEGER nstatv,kinc
	REAL(KIND=8) ::  statev(nstatv)

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

	INTEGER varsPerPart, ip, dex, id, numslip, grainid
      
	REAL(KIND=8) :: time(2)

!-----------------------------------------------------------------------

!	recover state variables from abaqus vector array
	dex=0

    if (Plap%crystalID(grainid) .eq. kHCP) then
        numslip=PhSlip(1)
    else if (Plap%crystalID(grainid) .eq. kBCC) then
        numslip=PhSlip(2)
    else 
        call RunTimeError( FILE_O,'Error: crystalID exceed the existing lattice!')
    end if

! 	stress
    dex = dex + 1   		
    gstress_n(1:NVec)=statev(dex:(dex+NVEC-1))
    if (debug==1) then   
        write(FILE_E,*) 'gstress_n:'
        write(FILE_E,'(6F12.5)')  gstress_n(:)
    endif     

!   estrain
    dex = dex + NVEC                      
    gestran_n(1:NVec)=statev(dex:(dex+NVEC-1))
    if (debug==1) then  		
        write(FILE_E,*) 'estran_n:'
        write(FILE_E,'(6F12.5)')  gestran_n(:) 
    end if
    
!	kappa
    dex = dex + NVEC                    
    gkappa_n(1:numslip)=statev(dex:(dex+NumSlip-1))
    if (debug==1) then  		
        write(FILE_E,*) 'gkappa_n:'
        write(FILE_E,'(6F12.5)') gkappa_n(:)
    endif

! 	VonMises, Shearate, gamtot
    dex = dex + numslip                 
    statev(dex:(dex+NEQVA/2-1))=pzero
    if (debug==1) then  
        write(FILE_E,*) '*----Initial VonMises, Shearate, gamtot in ', ip, '-th part----'
        write(FILE_E,'(6F12.5)')  statev(dex:(dex+NEQVA/2-1))
    endif         

!	gammadot
    dex = dex + NEQVA/2   
    ggamdot(1:numslip)=statev(dex:(dex+numslip-1))		
    if (debug==1) then  		
        write(FILE_E,*) 'ggamdot:'
        write(FILE_E,'(6F12.5)')  ggamdot(:)
    endif

!	mu
    dex = dex + NumSlip                 
    gmu_n(1:NVec)=statev(dex:(dex+NVEC-1)) 
    if (debug==1) then   
        write(FILE_E,*) 'gmu_n:'
        write(FILE_E,'(6F12.5)')  gmu_n(:)  
    endif

!	rss
    dex = dex + NVEC                  
    grss(1:numslip)=statev(dex:(dex+numslip-1))
    if (debug==1) then   
        write(FILE_E,*) 'grss:'
        write(FILE_E,'(6F12.5)')  grss(:)  
    endif		

! 	mudot
    dex = dex + NumSlip                
    gmudot(1:NVec)=statev(dex:(dex+NVEC-1))		
    if (debug==1) then   
        write(FILE_E,*) 'gmu_n:'
        write(FILE_E,'(6F12.5)')  gmudot(:)  
    endif

!	dd_for
    dex = dex + NVEC		
    gdd_for_n(1:numslip)=statev(dex:(dex+numslip-1)) 
    if (debug==1) then   
        write(FILE_E,*) 'dd_for:'
        write(FILE_E,*)  gdd_for_n(:)  
    endif

!	dd_revp
    dex = dex + NumSlip		
    gdd_revp_n(1:numslip)=statev(dex:(dex+numslip-1))
    if (debug==1) then   
        write(FILE_E,*) 'dd_revp:'
        write(FILE_E,'(6F12.5)')  gdd_revp_n(:)  
    endif
    
!	dd_revn
    dex = dex + NumSlip		
    gdd_revn_n(1:numslip)=statev(dex:(dex+numslip-1))
    if (debug==1) then   
        write(FILE_E,*) 'dd_revn:'
        write(FILE_E,'(6F12.5)')  gdd_revn_n(:)  
    endif

!	stressf
    dex = dex + NumSlip		
    gdd_stressf_n(1:numslip)=statev(dex:(dex+numslip-1)) 
    if (debug==1) then   
        write(FILE_E,*) 'stressf:'
        write(FILE_E,'(6F12.5)')  gdd_stressf_n(:)  
    endif

!	dd_rev0
    dex = dex + NumSlip		
    gdd_rev0_n(1:numslip)=statev(dex:(dex+numslip-1))
    if (debug==1) then   
        write(FILE_E,*) 'dd_rev0:'
        write(FILE_E,*)  gdd_rev0_n(:)  
    endif

!	dd_deb
    dex = dex + NumSlip		
    gdd_deb_n(1:numslip)=statev(dex:(dex+numslip-1))
    if (debug==1) then   
        write(FILE_E,*) 'dd_deb:'
        write(FILE_E,*)  gdd_deb_n(:)  
    endif

    RETURN
END


!***********************************************************************
SUBROUTINE FetchCrystalVariablesAtPart(gstress_n, gestran_n, gkappa_n, &
		   geqvalues,stress_n, estran_n, kappa_n, eqvalues, mu_n, 	   &
		   gmu_n, dd_for_n, gdd_for_n, dd_revp_n, gdd_revp_n,          &
		   dd_revn_n, gdd_revn_n, dd_stressf_n, gdd_stressf_n,         &
		   dd_rev0_n, gdd_rev0_n, dd_deb_n, gdd_deb_n)
!***********************************************************************         
         
	USE	NumType
	USE	FILEIO
	IMPLICIT NONE

	REAL(KIND=8) ::  gstress_n (NVEC)
	REAL(KIND=8) ::  gestran_n (NVEC)
	REAL(KIND=8) ::  gkappa_n  (MaxNumSlip)
	REAL(KIND=8) ::  geqvalues (NEQVA)
	REAL(KIND=8) ::  gmu_n    (NVEC) 

	REAL(KIND=8) ::  gdd_for_n     (MaxNumSlip)
	REAL(KIND=8) ::  gdd_revp_n    (MaxNumSlip)
	REAL(KIND=8) ::  gdd_revn_n    (MaxNumSlip)
	REAL(KIND=8) ::  gdd_stressf_n (MaxNumSlip)
	REAL(KIND=8) ::  gdd_rev0_n    (MaxNumSlip)
	REAL(KIND=8) ::  gdd_deb_n     (MaxNumSlip)

	REAL(KIND=8) ::  dd_for_n      (MaxNumSlip)
	REAL(KIND=8) ::  dd_revp_n     (MaxNumSlip)
	REAL(KIND=8) ::  dd_revn_n     (MaxNumSlip)
	REAL(KIND=8) ::  dd_stressf_n  (MaxNumSlip)
	REAL(KIND=8) ::  dd_rev0_n     (MaxNumSlip)
	REAL(KIND=8) ::  dd_deb_n      (MaxNumSlip)

	REAL(KIND=8) ::  stress_n(NVEC), estran_n(NVEC), kappa_n(MaxNumSlip),mu_n(NVEC) 
	REAL(KIND=8) ::  eqvalues(NEQVA)

!-----------------------------------------------------------------------

	stress_n= gstress_n
	estran_n= gestran_n
	kappa_n = gkappa_n
	mu_n=gmu_n

	dd_for_n=gdd_for_n
	dd_revp_n=gdd_revp_n
	dd_revn_n=gdd_revn_n
	dd_stressf_n=gdd_stressf_n
	dd_rev0_n=gdd_rev0_n
	dd_deb_n=gdd_deb_n
	
	eqvalues(kMISES_n)  = geqvalues(kMISES_n)
	eqvalues(kSHRATE_n) = geqvalues(kSHRATE_n)
	eqvalues(kGAMTOT_n) = geqvalues(kGAMTOT_n)		
    	
	RETURN
END

!***********************************************************************
SUBROUTINE ResetCrystalQnts(stress, estran, kappa, statev, eqvalues,   &
           gamdot, rss, mudot, mu, mu_n, s_ij, c_ijkl, stress_n,       &
           estran_n, kappa_n, dtime, VecM, VecS, ZVec, dd_for_n,       &
           dd_for, dd_revp_n, dd_revp, dd_revn_n, dd_revn,             &
           dd_stressf_n, dd_stressf, dd_rev0_n, dd_rev0, dd_deb_n,     &
           dd_deb, argmin, argmax, grainid   ) 
!***********************************************************************

	USE FILEIO
    USE NumType
	use PlaPar
	IMPLICIT NONE 
	
	INTEGER  	 ::  is, grainid, numslip
	REAL(KIND=8) ::  dtime, argmin,argmax

	REAL(KIND=8) ::  stress(NVEC), estran(NVEC), kappa(MaxNumSlip)
	REAL(KIND=8) ::  statev(NSTAV)
	REAL(KIND=8) ::  s_ij(NSD,NSD), c_ijkl(NVEC,NVEC)

	REAL(KIND=8) ::  stress_n(NVEC), estran_n(NVEC), kappa_n(MaxNumSlip)
	REAL(KIND=8) ::  gamdot(MaxNumSlip), mu(NVEC),mu_n(NVEC)
	REAL(KIND=8) ::  rss(MaxNumSlip), crss(MaxNumSlip) 
	REAL(KIND=8) ::  InnerProductVec, SSKineticEqn
	REAL(KIND=8) ::  mudot(NVEC)
	REAL(KIND=8) ::  eqvalues(NEQVA)	

	REAL(KIND=8) ::  dd_for        (MaxNumSlip)
	REAL(KIND=8) ::  dd_revp       (MaxNumSlip)
	REAL(KIND=8) ::  dd_revn       (MaxNumSlip)
	REAL(KIND=8) ::  dd_stressf    (MaxNumSlip)
	REAL(KIND=8) ::  dd_rev0       (MaxNumSlip)
	REAL(KIND=8) ::  dd_deb        (MaxNumSlip)

	REAL(KIND=8) ::  dd_for_n      (MaxNumSlip)
	REAL(KIND=8) ::  dd_revp_n     (MaxNumSlip)
	REAL(KIND=8) ::  dd_revn_n     (MaxNumSlip)
	REAL(KIND=8) ::  dd_stressf_n  (MaxNumSlip)
	REAL(KIND=8) ::  dd_rev0_n     (MaxNumSlip)
	REAL(KIND=8) ::  dd_deb_n      (MaxNumSlip)

	REAL(KIND = 8) :: VecM(3,MaxNumSlip)
	REAL(KIND = 8) :: VecS(3,MaxNumSlip)
	REAL(KIND = 8) :: ZTen(3,3,MaxNumSlip),ZVec(6,MaxNumSlip)
	
!-----------------------------------------------------------------------

	stress=stress_n
	estran=estran_n
	kappa = kappa_n

	dd_for=dd_for_n
	dd_revp=dd_revp_n
	dd_revn=dd_revn_n
	dd_stressf=dd_stressf_n
	dd_rev0=dd_rev0_n
	dd_deb=dd_deb_n
	
	eqvalues(kMISES)  = eqvalues(kMISES_n)             
	eqvalues(kSHRATE) = eqvalues(kSHRATE_n)
	eqvalues(kGAMTOT) = eqvalues(kGAMTOT_n)

!	Cauchy stress (tensor form)
	CALL Vec6x1ToMat3x3Symm(stress(1), s_ij(1,1), NSD)

	mudot(:)=pzero
	if (Plap%crystalID(grainid) .eq. kHCP) then
		numslip=PhSlip(1)
	else if (Plap%crystalID(grainid) .eq. kBCC) then
		numslip=PhSlip(2)
	else 
		call RunTimeError( FILE_O,'Error: crystalID exceed the existing lattice!')
	end if
	
	DO is = 1, numslip
!		Resolve shear stresses
	    rss(is) = InnerProductVec(stress(1), ZVec(1,is), NVEC)  
	    crss(is)= kappa(is)
!		Shear strain rate
	    gamdot(is) = SSKineticEqn(rss(is),crss(is),kDDGAM, argmin, argmax,is,grainid)
	    mudot(:)=mudot(:)+gamdot(is)*ZVec(:, is) 
	ENDDO

!	consistent tangent
	!CALL PlasticModuli(c_ijkl, stress, kappa, dtime, VecM, VecS, ZVec) 
	RETURN
END


!***********************************************************************
SUBROUTINE SaveCrystalVariablesAtPart(stress, estran, kappa, gamdot,   &
	       rss, mudot, eqvalues, gstress, gestran, gkappa, gmu, mu,    &
	       ggamdot, grss, gmudot, geqvalues, dd_for, gdd_for, dd_revp, &
	       gdd_revp, dd_revn, gdd_revn, dd_stressf, gdd_stressf,       &
	       dd_rev0, gdd_rev0, dd_deb, gdd_deb)
!*********************************************************************** 

    Use NumType
	USE FILEIO
	IMPLICIT NONE

	INTEGER      ::   igrn, iqpt, ielem,ip

	REAL(KIND=8) ::   gstress  (NVEC)
	REAL(KIND=8) ::   gestran  (NVEC)
	REAL(KIND=8) ::   gkappa   (MaxNumSlip)
	REAL(KIND=8) ::   eqvalues (NEQVA)
	REAL(KIND=8) ::   geqvalues(NEQVA)		
	REAL(KIND=8) ::   gmu      (NVEC)
	REAL(KIND=8) ::   ggamdot  (MaxNumSlip)
	REAL(KIND=8) ::   grss     (MaxNumSlip)

	REAL(KIND=8) ::   stress(NVEC), estran(NVEC), kappa(MaxNumSlip)
	REAL(KIND=8) ::   gamdot(MaxNumSlip),mu(NVEC), rss(MaxNumSlip)

	REAL(KIND=8) ::   gmudot(NVEC),mudot(NVEC)

	REAL(KIND=8) ::   gdd_for       (MaxNumSlip)
	REAL(KIND=8) ::   gdd_revp      (MaxNumSlip)
	REAL(KIND=8) ::   gdd_revn      (MaxNumSlip)
	REAL(KIND=8) ::   gdd_stressf   (MaxNumSlip)
	REAL(KIND=8) ::   gdd_rev0      (MaxNumSlip)
	REAL(KIND=8) ::   gdd_deb       (MaxNumSlip)

	REAL(KIND=8) ::   dd_for        (MaxNumSlip)
	REAL(KIND=8) ::   dd_revp       (MaxNumSlip)
	REAL(KIND=8) ::   dd_revn       (MaxNumSlip)
	REAL(KIND=8) ::   dd_stressf    (MaxNumSlip)
	REAL(KIND=8) ::   dd_rev0       (MaxNumSlip)
	REAL(KIND=8) ::   dd_deb        (MaxNumSlip)
	
!-----------------------------------------------------------------------

    gstress=stress
    gestran=estran
    gkappa =kappa
    ggamdot=gamdot
    grss   =rss
    gmudot =mudot
    gmu=mu

	gdd_for= dd_for
	gdd_revp= dd_revp
	gdd_revn= dd_revn
	gdd_stressf= dd_stressf
	gdd_rev0= dd_rev0
	gdd_deb= dd_deb
    
    geqvalues(kMISES) = eqvalues(kMISES)
    geqvalues(kSHRATE) = eqvalues(kSHRATE)
    geqvalues(kGAMTOT) = eqvalues(kGAMTOT)    

     
	RETURN 
END

!***********************************************************************
SUBROUTINE SaveStateVars(statev, nstatv, gstress, gestran, gkappa,     &
           ggamdot, grss, gmudot, gmu, geqvalues, gdd_for, gdd_revp,   &
           gdd_revn, gdd_stressf, gdd_rev0, gdd_deb, dd_all, dd_max,   &
           dd_his, dd_his_max, det_dd_max, ave_estran, det_dd_max_his, &
           grainid )
!***********************************************************************         
    USE NumType
	USE FILEIO
	Use Plapar
	IMPLICIT NONE
      
	INTEGER      ::  nstatv
	REAL(KIND=8) ::  statev(nstatv)

	REAL(KIND=8) ::  gstress   (NVEC)
	REAL(KIND=8) ::  gestran   (NVEC)
	REAL(KIND=8) ::  gkappa    (MaxNumSlip)

	REAL(KIND=8) ::  gmu       (NVEC)
	REAL(KIND=8) ::  ggamdot   (MaxNumSlip)
	REAL(KIND=8) ::  grss      (MaxNumSlip)

	REAL(KIND=8) ::  gdd_for   (MaxNumSlip)	
	REAL(KIND=8) ::  gdd_revp  (MaxNumSlip)	
	REAL(KIND=8) ::  gdd_revn  (MaxNumSlip)
	REAL(KIND=8) ::  gdd_stressf  (MaxNumSlip)
	REAL(KIND=8) ::  gdd_rev0  (MaxNumSlip)
	REAL(KIND=8) ::  gdd_deb   (MaxNumSlip)
	REAL(KIND=8) ::  dd_all 
	REAL(KIND=8) ::  dd_max
	REAL(KIND=8) ::  dd_his   
	REAL(KIND=8) ::  dd_his_max
	REAL(KIND=8) ::  det_dd_max
	REAL(KIND=8) ::  ave_estran (NVEC)
	REAL(KIND=8) ::  det_dd_max_his
	REAL(KIND=8) ::  geqvalues (NEQVA)	
	REAL(KIND=8) ::  gmudot(NVEC)

	INTEGER      ::  grainid, dex, id ,numslip
!-----------------------------------------------------------------------
       
!	save state variables in abaqus vector array
    ave_estran=0.d0		
	dex=0

	if (Plap%crystalID(grainid) .eq. kHCP) then
		numslip=PhSlip(1)
	else if (Plap%crystalID(grainid) .eq. kBCC) then
		numslip=PhSlip(2)
	else 
		call RunTimeError(FILE_O, 'Error: crystalID exceed the existing lattice!')
	end if

! 	stress         
	dex = dex + 1           
	statev(dex:(dex+NVEC-1))=gstress(1:NVEC)
	if (debug==1) then
		 Write(FILE_E,*) 'initial estimate gstress='
		 Write(FILE_E,*) gstress
	end if

! 	estrain        
	dex = dex + NVEC                    
	statev(dex:(dex+NVEC-1))=gestran(1:NVEC)

! 	kappa         
	dex = dex + NVEC                      
	statev(dex:(dex+NumSlip-1))=gkappa(1:numslip)
	
! 	VonMises, Shearate, gamtot
	dex = dex + numslip                 
	statev(dex:(dex+NEQVA/2-1))=geqvalues ((NEQVA/2+1):NEQVA)		

! 	gamdot         
	dex = dex + NEQVA/2           
	statev(dex:(dex+numslip-1))=ggamdot(1:numslip)

! 	average estran: this value is wrong    
	dex = dex + NumSlip                 
	statev(dex:(dex+NVEC-1))=gmu(1:NVEC) 
	ave_estran=ave_estran+gmu(1:NVEC)!*VolFrac(grainid)

! 	rss
	dex = dex + NVEC                  
	statev(dex:(dex+numslip-1))=grss(1:numslip)

! 	mudot
	dex = dex + NumSlip                  
	statev(dex:(dex+NVEC-1))=gmudot(1:NVEC)

! 	dd_for
	dex = dex + NVEC                     
	statev(dex:(dex+numslip-1))=gdd_for(1:numslip)

! 	dd_revp
	dex = dex + NumSlip                  
	statev(dex:(dex+numslip-1))=gdd_revp(1:numslip)

! 	dd_revn
	dex = dex + NumSlip                  
	statev(dex:(dex+numslip-1))=gdd_revn(1:numslip)

! 	stressf
	dex = dex + NumSlip                  
	statev(dex:(dex+numslip-1))=gdd_stressf(1:numslip)

! 	dd_rev0
	dex = dex + NumSlip                  
	statev(dex:(dex+numslip-1))=gdd_rev0(1:numslip)

! 	dd_deb
	dex = dex + NumSlip                  
	statev(dex:(dex+numslip-1))=gdd_deb(1:numslip)

	dex = dex + NumSlip-1    

!	dd in each part		
	statev(nstatv-11)=dd_all
	
!	dd_his max in each part
	statev(nstatv-10)=dd_his   
           
	statev(nstatv-9)=dd_max 											! dd_max among all parts
	statev(nstatv-8)=dd_his_max   										! dd_history_max among all parts
	statev(nstatv-7)=det_dd_max     									! dd_det_max_DD in each part
	statev((nstatv-6):(nstatv-1))=ave_estran(1:6) 						! average plastic strain in each part
	statev(nstatv)=det_dd_max_his										! dd_history_det_max_DD in each part

	RETURN
END

!***********************************************************************
SUBROUTINE SaveStressModuli(stress, ddsdde, savg_ij, cavg_ijkl)                                                
!***********************************************************************
   
    USE NumType
	USE FILEIO
	IMPLICIT NONE


	REAL(KIND=8) ::  stress(NVEC), ddsdde(NVEC, NVEC)
	REAL(KIND=8) ::  savg_ij(NSD, NSD), cavg_ijkl(NVEC, NVEC)

!-----------------------------------------------------------------------

!	stresses
	stress=pzero    
	stress(1) = savg_ij(1,1)
	stress(2) = savg_ij(2,2)
	stress(3) = savg_ij(3,3)
	stress(4) = savg_ij(1,2)
	stress(5) = savg_ij(1,3)
	stress(6) = savg_ij(2,3)

!	algorithmic moduli
	ddsdde=pzero
	ddsdde=ddsdde+ cavg_ijkl(:,:)

	
	RETURN
END

