!***********************************************************************
subroutine ComputeStressResidual(rhs, stress, stressini, estran, rss,  &
		   mudot, gamdot, kappa, dtime, argmin, argmax, dstrain, VecM, &
		   VecS, ZVec, grainid)      
!***********************************************************************
   
    use NumType
    use FILEIO
    use PlaPar
    use ElasMod
	implicit none

	REAL(KIND=8) 	:: stress(NVEC),stressini(NVEC)
	REAL(KIND=8) 	:: kappa(MaxNumSlip)

	REAL(KIND=8) 	:: argmin, argmax, dtime,dstrain(NVEC)
	REAL(KIND=8) 	:: rhs(NVec), PartResid(NVec)
	REAL(KIND=8) 	:: estran(NVEC)

	INTEGER 		:: is,ip, i,j, ip1,ip2, numslip, grainid         
	REAL(KIND=8) 	:: InnerProductVec, SSKineticEqn
	REAL(KIND=8) 	:: crss(MaxNumSlip),rss(MaxNumSlip), gamdot(MaxNumSlip),mudot(6)
    REAL(KIND=8) 	:: ZTen(NSD,NSD,MaxNumSlip),ZVec(NVEC,MaxNumSlip)      
	INTEGER 		:: StressSolveflat
    REAL(KIND=8) 	:: VecM(NSD,MaxNumSlip)
    REAL(KIND=8) 	:: VecS(NSD,MaxNumSlip)
!-----------------------------------------------------------------------

!	initialize eigenstrain rate
    mudot=pzero     

!	determine number of slip systems
	if (Plap%crystalID(grainid) .eq. kHCP) then
		numslip = PhSlip(1)
	else if (Plap%crystalID(grainid) .eq. kBCC) then
		numslip = PhSlip(2)
	else 
		call RunTimeError(FILE_O,'Error: crystalID exceed the existing lattice!')
	end if
		
    do is = 1, numslip
		rss(is) = InnerProductVec(stress(:), ZVec(:,is), NVEC)  
		crss(is) = kappa(is)          
        gamdot(is) = SSKineticEqn(rss(is),crss(is),KDDGAM, argmin, argmax,is,grainid)
        mudot(:) = mudot(:)+gamdot(is)*ZVec(:,is)           
    end do
                                   
    rhs = (stress-stressini)/dtime-MATMUL(Lijkl(:,:,grainid),(dstrain/dtime-mudot))

    return
END


!***********************************************************************
SUBROUTINE ComputeStressJacobian(lhs, rss, kappa, dtime, stress, VecM, &
		   VecS, Zvec, grainid)        
!***********************************************************************
		   
    use NumType
    use PlaPar
    use FILEIO
    use ElasMod
    implicit none

    real(kind=8)		:: dtime
    real(kind=8)		:: lhs(NVEC,NVEC), rss(MaxNumSlip)
    real(kind=8)  		:: kappa(MaxNumSlip)
    real(kind=8) 		:: lhsPart(6,6)

    real(kind=8) 		:: dmudotdsigma(6,6), I_P(6,6), Term1(6,6), Term2(6,6)
    real(kind=8) 		:: delta,tem
    integer 			:: is, ip, ip1, ip2, i, numslip, grainid
    real(kind=8)		:: crss(MaxnumSlip)
    real(kind=8)		:: stress(6)


    real(kind=8)		:: temp2(6,6),omega(MaxNumSlip)
    real(kind=8)		:: InnerProductVec
    real(kind=8)		:: SignOf


    real(kind=8)		:: VecM(3,MaxNumSlip)
    real(kind=8)		:: VecS(3,MaxNumSlip)
    real(kind=8)		:: ZTen(3,3,MaxNumSlip),ZVec(6,MaxNumSlip)
!-----------------------------------------------------------------------

!	initialization
    tem = 295.0
    dmudotdsigma = pzero
    omega = pzero
    
!	omega and dudotdsigma

	if (Plap%crystalID(grainid) .eq. kHCP) then
		numslip=PhSlip(1)
		do is = 1, numslip          

			if (is.le.3) then
				omega(is)=0.5*Plap%ddavg(1,1)*Plap%vid(1,1)*Plap%burgers(1,1)**2.*exp(((dabs(rss(is))&
				-kappa(is))-Plap%detF(1,1)/Plap%detV(1,1))*Plap%detV(1,1)/Plap%kbolt/tem)*Plap%detV(1,1)/Plap%kbolt/tem!*SignOf(rss(is,ip))
			else if (is.le.6) then
				omega(is)=0.5*Plap%ddavg(1,2)*Plap%vid(1,2)*Plap%burgers(1,2)**2.*exp(((dabs(rss(is))&
				-kappa(is))-Plap%detF(1,2)/Plap%detV(1,2))*Plap%detV(1,2)/Plap%kbolt/tem)*Plap%detV(1,2)/Plap%kbolt/tem!*SignOf(rss(is,ip))
			else if (is.le.12) then
				omega(is)=0.5*Plap%ddavg(1,3)*Plap%vid(1,3)*Plap%burgers(1,3)**2.*exp(((dabs(rss(is))&
				-kappa(is))-Plap%detF(1,3)/Plap%detV(1,3))*Plap%detV(1,3)/Plap%kbolt/tem)*Plap%detV(1,3)/Plap%kbolt/tem!*SignOf(rss(is,ip))
			else if (is.le.30) then
				omega(is)=0.5*Plap%ddavg(1,4)*Plap%vid(1,4)*Plap%burgers(1,4)**2.*exp(((dabs(rss(is))&
				-kappa(is))-Plap%detF(1,4)/Plap%detV(1,4))*Plap%detV(1,4)/Plap%kbolt/tem)*Plap%detV(1,4)/Plap%kbolt/tem!*SignOf(rss(is,ip))
			else
				call RunTimeError(File_o, 'ComputeJacobian: is exceed number of slip')
			end if

			call OuterProductVec(ZVec(1,is),ZVec(1,is),temp2,6) 
			 
			dmudotdsigma(:,:)=dmudotdsigma(:,:)+omega(is)*temp2           
	   enddo

	else if (Plap%crystalID(grainid) .eq. kBCC) then
		numslip=PhSlip(2)
		do is = 1, numslip  
			omega(is)=0.5*Plap%ddavg(2,1)*Plap%vid(2,1)*Plap%burgers(2,1)**2.*exp(((dabs(rss(is))&
			-kappa(is))-Plap%detF(2,1)/Plap%detV(2,1))*Plap%detV(2,1)/Plap%kbolt/tem)*Plap%detV(2,1)/Plap%kbolt/tem!*SignOf(rss(is,ip)

			call OuterProductVec(ZVec(1,is),ZVec(1,is),temp2,6)  
			dmudotdsigma(:,:)=dmudotdsigma(:,:)+omega(is)*temp2 
		enddo
	else 
		call RunTimeError(FILE_O,'Error: crystalID exceed the existing lattice!')
	end if


!	loacal jacobian
	lhs = 1/dtime*Ident4th+MATMUL(Lijkl(:,:,grainid),dmudotdsigma)
	        
	return  
end

!***********************************************************************
SUBROUTINE ComputeStress(rhs, stress, stressini, estran, rss, mudot,     &
		   gamdot, kappa, dtime, argmin, argmax, dstrain, VecM, VecS,  &
		   ZVec, iterCounterN, tolerNewt, kstep, kinc, ierr, noel, grainid) 
!***********************************************************************
         
    use NumType
    use FILEIO
	implicit none

    INTEGER      	:: iterCounterN, kstep, kinc, info, ierr, grainid
	REAL(KIND=8) 	:: stress(NVEC), stress0(NVEC), stressini(NVEC)
	REAL(KIND=8) 	:: kappa(MaxNumSlip)
	REAL(KIND=8) 	:: argmin, argmax, dtime,dstrain(NVEC), tolerNewt
	REAL(KIND=8) 	:: rhs(6), PartResid(6), lhs(6, 6)
	REAL(KIND=8) 	:: estran(NVEC)
	REAL(KIND=8) 	:: norm_stress0, norm_kappa0, norm_stress, norm_kappa
	REAL(KIND=8) 	:: rhs_norm, rhs_norm0	
	REAL(KIND=8) 	:: delstress(NVEC), delstressvec(NVEC)	

	INTEGER 		:: is,ip, i,j, ip1,ip2, iterNewt, noel
	INTEGER 		:: indx(NVEC)	
    REAL(KIND=8) 	:: InnerProductVec, SSKineticEqn, search
	REAL(KIND=8) 	:: crss(MaxNumSlip),rss(MaxNumSlip), gamdot(MaxNumSlip)
	REAL(KIND=8) 	:: delIdent4th(6,6)
	REAL(KIND=8) 	:: mudot(6), I_P(6,6), mult1(6), mult2(6),temp1(6)
		
	INTEGER 		:: StressSolveflat
	LOGICAL 		:: converged, converged2


    REAL(KIND=8) 	:: VecM(NSD,MaxNumSlip)
    REAL(KIND=8) 	:: VecS(NSD,MaxNumSlip)
    REAL(KIND=8) 	:: ZTen(NSD,NSD,MaxNumSlip),ZVec(NVEC,MaxNumSlip)
!-----------------------------------------------------------------------  
    
    if (debug==1) write(file_e,*) '-------------Begin Stress Solve-----------------'        
      
    call ComputeStressResidual(rhs, stress, stressini, estran, rss,  &
		   mudot, gamdot, kappa, dtime, argmin, argmax, dstrain, VecM, &
		   VecS, ZVec, grainid)    
     
    if (debug==1 .and. noel==1)   write(FILE_E, '(A30, /, 6e12.5)') 'Initial Residual is:', rhs
   
!	initial valuies for the norm of stress of each part        
    rhs_norm0=norm2(rhs)   
  
    iterNewt = 0   
    converged = .false. 
           
    do while (iterNewt .lt. iterCounterN  .and. .not. converged )               
        iterNewt=iterNewt+1
        if (debug==1 .and. noel==1) write(file_e, *) 'iterNewt= ', iterNewt        
        stress0=stress  

        if (debug==1) then
            write(file_E, *) 'stress= '      
            write(file_E, *) stress(:)        
        endif

!	solve for the crystal stresses
!		compute local jacobian
        call  ComputeStressJacobian(lhs, rss, kappa, dtime, stress, VecM, &
		   VecS, Zvec, grainid)     

        if (debug==1 .and. noel==1) then
			write(file_E, *) 'kinc=', kinc
            write(file_E, *) 'lhs= '      
            write(file_E, '( 6(2x, E12.5))') lhs
            write(file_E, *) 'rhs= '      
            write(file_E, '( 6(2x, E12.5))') rhs   
            write(file_E,*) 'dstrain='
            write(file_E, '( 6(2x, E12.5))') dstrain      
        endif                      
!        
!		solve for the increment of stress
        CALL DGESV( 6, 1, lhs, 6, indx, rhs, 6, INFO )  
	    if (debug==1) write(File_E,*) 'INFO= ', INFO
        if (INFO<0) then
            ierr=XTAL_SING_JACOBIAN          
            write(File_o,*) 'StressSolve: Jacobian is singular when kinc=', kinc
            return
        endif         
        
!		Please notice that here both delstress and rhs are 6*NPart by one vector
!		we just rearrange them to be in the matrix form
        delstressVec=rhs 
        delstress=reshape(delstressVec, shape(delstress)) 
        search = pone                  
        stress=stress0-search*delstress   
             
        if (debug==1 .and. noel==1) then
            write(file_E,*) 'delstress=' 
            write(file_E, '( 6(2x, E12.5))') delstress      
        endif    
               
        call ComputeStressResidual(rhs, stress, stressini, estran, rss,  &
		   mudot, gamdot, kappa, dtime, argmin, argmax, dstrain, VecM, &
		   VecS, ZVec, grainid) 

		if (debug==1)   write(FILE_E, '(A30, /, 12e12.5)') 'New Residual is:', rhs
		                
!	update stresses
        rhs_norm=norm2(rhs)
!       Add a linear search
        do while ((rhs_norm .gt. rhs_norm0) .and. (kinc .gt. 1) )        
            search = search*0.5
            if (search .lt. TINY(1.d0)) then
               call WriteWarning(File_O, 'StressSolveDeviatoric: LS Failed, search < TINY')
               ierr = XTAL_LS_FAILED
               return
            endif                  
            stress=stress0-search*delstress

            call ComputeStressResidual(rhs, stress, stressini, estran, rss,  &
		   mudot, gamdot, kappa, dtime, argmin, argmax, dstrain, VecM, &
		   VecS, ZVec, grainid) 

            rhs_norm=norm2(rhs)                                       
        enddo 
        rhs_norm0 = rhs_norm 
        converged=converged2(rhs, tolerNewt, 6)
    enddo
	
	if (noel==1 .and. debug==1) then
		write(file_e,*) 'at kinc = ', kinc, 'stress increment is:'
		write(file_e,'(6E20.10)') stressini
		write(file_E,*) 'dstrain='
        write(file_E, '( 6(2x, E12.5))') dstrain    
	end if
	
	if (iterNewt.ge. iterCounterN) then
		call WriteWarning(FILE_O, 'ComputeStress: Netwon iters > iterCounterN')
		ierr = XTAL_MAX_ITERS_HIT
		return
    endif 

!	echo data if convergence is achieved
	if (converged .and. debug==1) then
		write(File_E,*) 'stress is'
		write(File_E,'(6E15.8)') stress(:)
	end if
    return
end      

!***********************************************************************
SUBROUTINE PlasticModuli(cepmod, stress, kappa, dtime, VecM, VecS, ZVec, grainid, noel) 
!***********************************************************************

    use NumType
    use PlaPar
    USE FILEIO
    use ElasMod
    implicit none
         
    real(kind=8)		:: cepmod(NVEC,NVEC)
    real(kind=8)		:: Stress(NVEC)
    real(kind=8)		:: kappa(MaxNumSlip) 
    real(kind=8)		:: dmudotdsigma(6,6),dmudotdrss
    real(kind=8)		:: rss(MaxNumSlip), ZZT(6,6,MaxNumSlip)

    integer				:: ip,is,i,j,numslip, grainid, noel

    real(kind=8)		:: del, I_P(6,6)   
    real(kind=8)		:: InnerProductVec

    real(kind=8), ALLOCATABLE :: A(:,:), ATemp(:,:), F(:,:)
    real(kind=8)		:: phiab(6,6),etaab(6,6),fb(6,6)
    real(kind=8)		:: dtime,tem
    real(kind=8)		:: ddd
    real(kind=8)		:: SignOf
    integer 			:: Cepmodflag
    integer 			:: indx(6)

    REAL(kind=8)		:: VecM(3,MaxNumSlip)
    REAL(kind=8)		:: VecS(3,MaxNumSlip)
    REAL(kind=8)		:: ZTen(3,3,MaxNumSlip),ZVec(6,MaxNumSlip)
 	REAL(kind=8)		:: start_time, end_time
!----------------------------------------------------------------------- 

!	Caculate d(sigma_n+1^(beta))/d(epsilonbar_(n+1))
    allocate (A(6,6), ATemp(6,6), F(6,6))
    tem = 295.0d0
    Cepmodflag = 0
    A = pzero
    dmudotdsigma = pzero
    dmudotdrss = pzero 

!	Caculate the rss    
	if (Plap%crystalID(grainid) .eq. kHCP) then
		numslip=PhSlip(1)

		do is=1,NumSlip
			rss(is)=InnerProductVec(stress(:),ZVec(:,is),6)
			Call OuterProductVec(ZVec(:,is),ZVec(:,is),ZZT(:,:,is),6)

			if (is.le.3) then
				dmudotdrss=0.5*Plap%ddavg(1,1)*Plap%vid(1,1)*Plap%burgers(1,1)**2.*exp(((dabs(rss(is))&
				-kappa(is))-Plap%detF(1,1)/Plap%detV(1,1))*Plap%detV(1,1)/Plap%kbolt/tem)*Plap%detV(1,1)/Plap%kbolt/tem!*SignOf(rss(is,ip))
			else if (is.le.6) then
				dmudotdrss=0.5*Plap%ddavg(1,2)*Plap%vid(1,2)*Plap%burgers(1,2)**2.*exp(((dabs(rss(is))&
				-kappa(is))-Plap%detF(1,2)/Plap%detV(1,2))*Plap%detV(1,2)/Plap%kbolt/tem)*Plap%detV(1,2)/Plap%kbolt/tem!*SignOf(rss(is,ip))
			else if (is.le.12) then
				dmudotdrss=0.5*Plap%ddavg(1,3)*Plap%vid(1,3)*Plap%burgers(1,3)**2.*exp(((dabs(rss(is))&
				-kappa(is))-Plap%detF(1,3)/Plap%detV(1,3))*Plap%detV(1,3)/Plap%kbolt/tem)*Plap%detV(1,3)/Plap%kbolt/tem!*SignOf(rss(is,ip))
			else if (is.le.30) then
				dmudotdrss=0.5*Plap%ddavg(1,4)*Plap%vid(1,4)*Plap%burgers(1,4)**2.*exp(((dabs(rss(is))&
				-kappa(is))-Plap%detF(1,4)/Plap%detV(1,4))*Plap%detV(1,4)/Plap%kbolt/tem)*Plap%detV(1,4)/Plap%kbolt/tem!*SignOf(rss(is,ip))
			else
				 call RunTimeError(File_o, 'Plastic Moduli: is exceed number of slip')
			end if

			dmudotdsigma(:,:)=dmudotdsigma(:,:)+dmudotdrss*ZZT(:,:,is)
		enddo

	else if (Plap%crystalID(grainid) .eq. kBCC) then
		numslip=PhSlip(2)
		do is=1,NumSlip
			rss(is)=InnerProductVec(stress(:),ZVec(:,is),6)
			Call OuterProductVec(ZVec(:,is),ZVec(:,is),ZZT(:,:,is),6)

			dmudotdrss=0.5*Plap%ddavg(2,1)*Plap%vid(2,1)*Plap%burgers(2,1)**2.*exp(((dabs(rss(is))&
			-kappa(is))-Plap%detF(2,1)/Plap%detV(2,1))*Plap%detV(2,1)/Plap%kbolt/tem)*Plap%detV(2,1)/Plap%kbolt/tem!*SignOf(rss(is,ip))

			dmudotdsigma(:,:)=dmudotdsigma(:,:)+dmudotdrss*ZZT(:,:,is)
		enddo

	else 
		call RunTimeError( FILE_O, 'Error: crystalID exceed the existing lattice!')
	end if
	

!---> LY changed
    call CheckMatriceZero(dmudotdsigma(1,1),6)
!---> LY changed

	A(:,:) = Ident4th/dtime + MATMUL(Lijkl(:,:,grainid),dmudotdsigma)
	F(:,:) = Lijkl(:,:,grainid)/dtime

!	Solve for cepmod as a whole and cepmod in each part is stored in cepmod     
	CALL DGESV( 6, 6, A, 6, indx, F, 6, Cepmodflag )   

!	Split cepmodm into cepmod(:,:,ip)   
    cepmod = F

    deallocate (A, ATemp, F)

!	echo data for debugging
	if(debug==1) then
		write(FILE_E,*)  'at noel',noel,'tangential moduli is'
		write(FILE_E,'(6E20.10)') cepmod(:,:)
	end if
	
    return
END            


!***********************************************************************
SUBROUTINE IntegrateHardening(kappa, kappa_n, eqvalues, mu, mu_n,      &
		   estran, estran_n, stress, rss, mudot, gamdot, dtime, 	   &
		   kInteg_Code, VecM, VecS,ZVec, argmin, argmax, dd_for,       &
		   dd_revp, dd_revn, stressf, dd_rev0, dd_deb, time, dstrain,  &
		   dd_all, dd_max, noel, grainid) 
!***********************************************************************
                           
    use NumType
    use PlaPar
    use FileIO    
    implicit none

    real(kind=8)	:: estran(6),estran_n(6)
    real(kind=8)	:: kappa(MaxNumSlip), kappa_n(MaxNumSlip)
    real(kind=8)	:: Stress(NVEC), eqvalues(NEQVA)
    real(kind=8)	:: crss(MaxNumSlip),rss(MaxNumSlip)
    real(kind=8)	:: gamdot(MaxNumSlip), SHEARATE(MaxNumSlip)
    real(kind=8)	:: mu(NVEC), mudot(NVEC)
	real(kind=8)	:: dd_for   (MaxNumSlip)	
	real(kind=8)	:: dd_revp  (MaxNumSlip)	
	real(kind=8)	:: dd_revn  (MaxNumSlip)
	real(kind=8)	:: stressf  (MaxNumSlip)
	real(kind=8)	:: dd_rev0  (MaxNumSlip)
	real(kind=8)	:: dd_deb   (MaxNumSlip)
	real(kind=8)	:: dstrain  (NVEC)
	real(kind=8)	:: dd_deb_tot
	real(kind=8)	:: dd_all    
    real(kind=8)	:: argmin, argmax, kappa_sat
    real(kind=8)	:: c, g_n, g_s, g
    real(kind=8)	:: dkappa, gamtot_n, mu_n(NVEC), delgam, fac
    real(kind=8)	:: kTHETA
    real(kind=8)	:: dtime, time
    real(kind=8)	:: erate
    real(kind=8)	:: k2
    real(kind=8)	:: dfordsh
    real(kind=8)	:: drevpdsh
    real(kind=8)	:: drevndsh
    real(kind=8)	:: ddebdsh
    real(kind=8)	:: dkappa1
    real(kind=8)	:: dkappa2
    real(kind=8)	:: dkappa3
    real(kind=8)	:: dkappa3_s
    real(kind=8)	:: tem
    real(kind=8)	:: temrate
    real(kind=8)	:: pcoef
    real(kind=8)	:: rate_ref
    real(kind=8)	:: inter_ref
    real(kind=8)	:: kdeb
    real(kind=8)	:: dd_max
    real(kind=8)	:: dd_min
    real(kind=8)	:: dd_tot(MaxNumSlip)
    real(kind=8)	:: InnerProductVec, SSKineticEqn
    real(kind=8)	:: VecM(3,MaxNumSlip)
    real(kind=8)	:: VecS(3,MaxNumSlip)
    real(kind=8)	:: ZTen(3,3,MaxNumSlip),ZVec(6,MaxNumSlip)
    
    integer 		:: is, ip, jk, kInteg_Code, numslip, grainid, noel
    
    data    tem             /295.0d0/
    data    temrate         /0.0d0/
    data    pcoef           /0.5d0/
    data    rate_ref        /1.0d7/
    data    inter_ref       /0.9d0/
    data    kdeb            /0.86d0/
    data    kTHETA 			/1.0d0/
!-----------------------------------------------------------------------

	if (debug==1)   write(File_E,*) '----  Start IntegrateHardening!  ----'
    
	crss = pzero
    dd_deb_tot = pzero

!	dd
	mudot(:)=pzero
	
!	determine phase
	if (Plap%crystalID(grainid) .eq. kHCP) then
		numslip=PhSlip(1)
	else if (Plap%crystalID(grainid) .eq. kBCC) then
		numslip=PhSlip(2)
	else 
		call RunTimeError( FILE_O, 'Error: crystalID exceed the existing lattice!')
	end if

	dd_all=0.0d0
	do is = 1, numslip
		rss(is) = DOT_PRODUCT(stress(:), ZVec(:,is))  
		crss(is)= kappa(is)
		gamdot(is) = SSKineticEqn(rss(is),crss(is),KDDGAM, argmin, argmax,is,grainid)     
		mudot(:)=mudot(:) + gamdot(is)*ZVec(:, is)
		dd_tot(is)=dd_for(is)+dd_revp(is)+dd_revn(is)
		dd_deb_tot=dd_deb_tot+dd_deb(is)

		if ((dd_deb(is)+dd_tot(is)) .GE. dd_all) dd_all=dd_deb(is)+dd_tot(is)

		if (sign(1.0,rss(is)) .ne. stressf(is)) then
			stressf(is)=sign(1.0,rss(is))
			dd_rev0(is)=dd_tot(is)
		end if
	enddo
	if (debug==1) write(FILE_E,*) dd_all   

!	dd_max
    dd_max=0.0d0
	if (dd_all .GE. dd_max) dd_max=dd_all

!	echo some data for debugging
	if (debug==1 .and. noel==1) then 
		write(FILE_E,*) 'rss:'
		write(FILE_E,*)  rss(:)  
		write(FILE_E,*) 'dd_for:'
        write(FILE_E,*)  dd_for(:)   
        write(FILE_E,*) 'dd_revp:'
        write(FILE_E,*)  dd_revp(:)   
        write(FILE_E,*) 'dd_revn:'
        write(FILE_E,*)  dd_revn(:)    
        write(FILE_E,*) 'stressf:'
        write(FILE_E,*)  stressf(:)   
        write(FILE_E,*) 'dd_rev0:'
        write(FILE_E,*)  dd_rev0(:)    
        write(FILE_E,*) 'dd_deb:'
        write(FILE_E,*)  dd_deb(:)  
	endif

    eqvalues(kSHRATE) = pzero
  

	if (Plap%crystalID(grainid) .eq. kHCP) then
		numslip=PhSlip(1)
	else if (Plap%crystalID(grainid) .eq. kBCC) then
		numslip=PhSlip(2)
	else 
		call RunTimeError( FILE_O, 'Error: crystalID exceed the existing lattice!')
	end if

	do is = 1, numslip
		eqvalues(kSHRATE) = eqvalues(kSHRATE) + dabs(gamdot(is))
	enddo

    
!	accumulated shear strain: gamtot 
	eqvalues(kGAMTOT) = eqvalues(kGAMTOT_n) + eqvalues(kSHRATE)*dtime
	eqvalues(kMISES) = 0.5d0 *sqrt((Stress(1))**2  		   &
		                  + (Stress(2))**2                 &
		                  + (Stress(3))**2                 &
		+ 6.d0*(Stress(6)**2+Stress(5)**2+Stress(4)**2))        
   
!	inelastic strain mu and total elastic strain for each part.
    mu = mu_n
	mu(:) = mu_n(:) + mudot(:)*dtime   

    if (debug==1) then 
        write(file_e,*) 'gamdot(:): '  
        write(file_e,*) gamdot(:) 
        write(file_e,*) 'eqvalues(kshrate), eqvalues(kgamtot),  mu: '
        write(file_e,*) eqvalues(kshrate), eqvalues(kgamtot), mu(:)
    endif
    
!	integration of hardening law (one hardness/slip system)
    if (debug==1) write(FILE_E, *) 'kinteg_code =', kinteg_code

    if (kInteg_Code .eq. kHARD_EXPL) then    
!		explicit update     
		if (Plap%crystalID(grainid) .eq. kHCP) then
			numslip=PhSlip(1)
		else if (Plap%crystalID(grainid) .eq. kBCC) then
			numslip=PhSlip(2)
		else 
			call RunTimeError( FILE_O,'Error: crystalID exceed the existing lattice!')
		end if

		do is = 1, numslip
			kappa_sat = Plap%taus0 * ((eqvalues(kSHRATE) / Plap%gamss0)**Plap%xms) !saturated stress
			c = dtime*Plap%h0
			g_s = kappa_sat - Plap%tausi
			g_n = kappa_n(is) - Plap%tausi

			if ( (g_n/g_s) .le. 1.0 ) then
			   g = g_n + c*(1.0-g_n/g_s)*eqvalues(kSHRATE)             
			else
			   g = g_n
			endif
			kappa(is) = g + Plap%tausi
			
			if(debug==1) then
				write(file_e,*) 'dkappa, dkappa/dtime='
				write(file_e,*) kappa(is)-kappa_n(is), (kappa(is)-kappa_n(is))/dtime
				write(File_E,*) 'kappa(is) is:' 
				write(File_E,*) kappa(is)
			end if
		 enddo

      

	else if (kInteg_Code .eq. kHARD_MIDP) then
!		generalized mid-point rule 
		if (Plap%crystalID(grainid) .eq. kHCP) then
			numslip=PhSlip(1)
		else if (Plap%crystalID(grainid) .eq. kBCC) then
			numslip=PhSlip(2)
		else 
			call RunTimeError( FILE_O,'Error: crystalID exceed the existing lattice!')
		end if
		
		do is = 1, numslip
			kappa_sat = Plap%taus0 * ((eqvalues(kSHRATE) / Plap%gamss0)**Plap%xms)
			c = dtime*Plap%h0
			g_s = kappa_sat - Plap%tausi
			g_n = kappa_n(is) - Plap%tausi
			if ( (g_n/g_s) .le. 1.0 ) then
			   g = g_n + c*((1.0-kTHETA)*(1.0-g_n/g_s)*eqvalues(kSHRATE_n)+(kTHETA)*eqvalues(kSHRATE))
			   g = g / (1.0 + c*kTHETA*eqvalues(kSHRATE)/g_s)
			else
			   g = g_n
			endif
			kappa(is) = g + Plap%tausi
		enddo

        
	else if (kInteg_Code .eq. kHARD_ANAL) then
!		annihilation          
		gamtot_n = eqvalues(kGAMTOT_n)
		delgam   = eqvalues(kSHRATE) * dtime
		!write(*,*) 'gamtot_n, delgam=', gamtot_n, delgam
		if (Plap%crystalID(grainid) .eq. kHCP) then
			numslip=PhSlip(1)
		else if (Plap%crystalID(grainid) .eq. kBCC) then
			numslip=PhSlip(2)
		else 
			call RunTimeError( FILE_O,'Error: crystalID exceed the existing lattice!')
		end if
		
		do is= 1, numslip
			kappa_sat = Plap%taus0 * ((eqvalues(kSHRATE) / Plap%gamss0)**Plap%xms)
			g_s = kappa_sat - Plap%tausi
			fac = dabs(Plap%h0/g_s)
			dkappa = 0.0
			
			do jk = 1, numslip
				dkappa = dkappa + dabs(gamdot(jk))*dtime/delgam 
			enddo

			dkappa = dkappa * g_s * exp(-gamtot_n*fac) * (1.0 - exp(-delgam*fac))
			if (debug==1) then 
				write(file_e,*) 'kappa_sat, g_s,  fac, dkappa: ', kappa_sat, g_s,  fac, dkappa
			endif                                              
			kappa(is) = kappa_n(is) + dkappa	
			if (debug==1) write(*,*)  kappa(is)
		enddo

	else if (kInteg_Code .eq. kHARD_DD) then
	
		if (debug==1)   write(File_E,*) '-------------------  kHARD_DD begin  ------------------------------'
		
		dkappa3=pzero

		if (Plap%crystalID(grainid) .eq. kHCP) then
			numslip=PhSlip(1)
            do is= 1, numslip
				erate=dabs(gamdot(is))
				if (is.le.3) then
					k2 = Plap%k1(1,1)*Plap%burgers(1,1)*inter_ref/Plap%gref(1,1)*(1.-Plap%kbolt*tem &
						 *log(erate/rate_ref)/(Plap%Drag(1,1))/Plap%burgers(1,1)**3.0)
					dfordsh=(1.-pcoef)*Plap%k1(1,1)*sqrt(dd_tot(is))-k2*dd_for(is)

					if (stressf(is) .ge. 0.0d0) then
						drevpdsh=pcoef*Plap%k1(1,1)*sqrt(dd_tot(is))-k2*dd_revp(is)
						drevndsh=-Plap%k1(1,1)*sqrt(dd_tot(is))*(dd_revn(is)/dd_rev0(is))
					else
						drevndsh=pcoef*Plap%k1(1,1)*sqrt(dd_tot(is))-k2*dd_revn(is)
						drevpdsh=-Plap%k1(1,1)*sqrt(dd_tot(is))*(dd_revp(is)/dd_rev0(is))
					end if

					dkappa1=Plap%sref(1,1)/Plap%sita(1,1)*(1.-exp((tem-Plap%sita_ref(1,1))&
					/Plap%sita(1,1)))*temrate*dtime

					dkappa2=0.5*Plap%burgers(1,1)*inter_ref*Plap%shmu(1,1)/sqrt(dd_tot(is))&
					*(dfordsh+drevndsh+drevpdsh)*abs(gamdot(is))*dtime

					kappa(is)=kappa(is)+dkappa1+dkappa2

					dd_for(is)=dd_for(is)+dfordsh*abs(gamdot(is))*dtime

					dd_revp(is)=dd_revp(is)+drevpdsh*abs(gamdot(is))*dtime

					dd_revn(is)=dd_revn(is)+drevndsh*abs(gamdot(is))*dtime

					ddebdsh=Plap%qsub(1,1)*Plap%burgers(1,1)*sqrt(dd_deb_tot)*k2*dd_for(is)
	
					dkappa3_s=-0.5*kdeb*Plap%shmu(1,1)*Plap%burgers(1,1)/sqrt(dd_deb_tot)*&
					(log(Plap%burgers(1,1)*sqrt(dd_deb_tot))+1.0)*ddebdsh*abs(gamdot(is))*dtime

					dkappa3=dkappa3+dkappa3_s

					dd_deb(is)=dd_deb(is)+ddebdsh*abs(gamdot(is))*dtime

				else if (is.le.6) then
					k2=Plap%k1(1,2)*Plap%burgers(1,2)*inter_ref/Plap%gref(1,2)*(1.-Plap%kbolt*tem		&
					*log(erate/rate_ref)/(Plap%Drag(1,2))/Plap%burgers(1,2)**3.0)

					dfordsh=(1.-pcoef)*Plap%k1(1,2)*sqrt(dd_tot(is))-k2*dd_for(is)

					if (stressf(is) .ge. 0.0d0) then
						drevpdsh=pcoef*Plap%k1(1,2)*sqrt(dd_tot(is))-k2*dd_revp(is)
						drevndsh=-Plap%k1(1,2)*sqrt(dd_tot(is))*(dd_revn(is)/dd_rev0(is))

					else
						drevndsh=pcoef*Plap%k1(1,2)*sqrt(dd_tot(is))-k2*dd_revn(is)
						drevpdsh=-Plap%k1(1,2)*sqrt(dd_tot(is))*(dd_revp(is)/dd_rev0(is))

					end if

					dkappa1=Plap%sref(1,2)/Plap%sita(1,2)*(1.-exp((tem-Plap%sita_ref(1,2))&
					/Plap%sita(1,2)))*temrate*dtime

					dkappa2=0.5*Plap%burgers(1,2)*inter_ref*Plap%shmu(1,2)/sqrt(dd_tot(is))&
					*(dfordsh+drevndsh+drevpdsh)*abs(gamdot(is))*dtime

					kappa(is)=kappa(is)+dkappa1+dkappa2

					dd_for(is)=dd_for(is)+dfordsh*abs(gamdot(is))*dtime

					dd_revp(is)=dd_revp(is)+drevpdsh*abs(gamdot(is))*dtime

					dd_revn(is)=dd_revn(is)+drevndsh*abs(gamdot(is))*dtime

					ddebdsh=Plap%qsub(1,2)*Plap%burgers(1,2)*sqrt(dd_deb_tot)*k2*dd_for(is)
	
					dkappa3_s=-0.5*kdeb*Plap%shmu(1,2)*Plap%burgers(1,2)/sqrt(dd_deb_tot)*&
					(log(Plap%burgers(1,2)*sqrt(dd_deb_tot))+1.0)*ddebdsh*abs(gamdot(is))*dtime

					dkappa3=dkappa3+dkappa3_s

					dd_deb(is)=dd_deb(is)+ddebdsh*abs(gamdot(is))*dtime

				else if (is.le.12) then

					k2=Plap%k1(1,3)*Plap%burgers(1,3)*inter_ref/Plap%gref(1,3)*(1.-Plap%kbolt*tem		&
					*log(erate/rate_ref)/(Plap%Drag(1,3))/Plap%burgers(1,3)**3.0)

					dfordsh=(1.-pcoef)*Plap%k1(1,3)*sqrt(dd_tot(is))-k2*dd_for(is)

					if (stressf(is) .ge. 0.0d0) then
						drevpdsh=pcoef*Plap%k1(1,3)*sqrt(dd_tot(is))-k2*dd_revp(is)
						drevndsh=-Plap%k1(1,3)*sqrt(dd_tot(is))*(dd_revn(is)/dd_rev0(is))

					else
						drevndsh=pcoef*Plap%k1(1,3)*sqrt(dd_tot(is))-k2*dd_revn(is)
						drevpdsh=-Plap%k1(1,3)*sqrt(dd_tot(is))*(dd_revp(is)/dd_rev0(is))

					end if

					dkappa1=Plap%sref(1,3)/Plap%sita(1,3)*(1.-exp((tem-Plap%sita_ref(1,3))&
					/Plap%sita(1,3)))*temrate*dtime

					dkappa2=0.5*Plap%burgers(1,3)*inter_ref*Plap%shmu(1,3)/sqrt(dd_tot(is))&
					*(dfordsh+drevndsh+drevpdsh)*abs(gamdot(is))*dtime

					kappa(is)=kappa(is)+dkappa1+dkappa2

					dd_for(is)=dd_for(is)+dfordsh*abs(gamdot(is))*dtime

					dd_revp(is)=dd_revp(is)+drevpdsh*abs(gamdot(is))*dtime

					dd_revn(is)=dd_revn(is)+drevndsh*abs(gamdot(is))*dtime

					ddebdsh=Plap%qsub(1,3)*Plap%burgers(1,3)*sqrt(dd_deb_tot)*k2*dd_for(is)
	
					dkappa3_s=-0.5*kdeb*Plap%shmu(1,3)*Plap%burgers(1,3)/sqrt(dd_deb_tot)*&
					(log(Plap%burgers(1,3)*sqrt(dd_deb_tot))+1.0)*ddebdsh*abs(gamdot(is))*dtime

					dkappa3=dkappa3+dkappa3_s

					dd_deb(is)=dd_deb(is)+ddebdsh*abs(gamdot(is))*dtime

				else if (is.le.30) then
				
					k2=Plap%k1(1,4)*Plap%burgers(1,4)*inter_ref/Plap%gref(1,4)*(1.-Plap%kbolt*tem		&
					*log(erate/rate_ref)/(Plap%Drag(1,4))/Plap%burgers(1,4)**3.0)

					dfordsh=(1.-pcoef)*Plap%k1(1,4)*sqrt(dd_tot(is))-k2*dd_for(is)

					if (stressf(is) .ge. 0.0d0) then
						drevpdsh=pcoef*Plap%k1(1,4)*sqrt(dd_tot(is))-k2*dd_revp(is)
						drevndsh=-Plap%k1(1,4)*sqrt(dd_tot(is))*(dd_revn(is)/dd_rev0(is))

					else
						drevndsh=pcoef*Plap%k1(1,4)*sqrt(dd_tot(is))-k2*dd_revn(is)
						drevpdsh=-Plap%k1(1,4)*sqrt(dd_tot(is))*(dd_revp(is)/dd_rev0(is))

					end if

					dkappa1=Plap%sref(1,4)/Plap%sita(1,4)*(1.-exp((tem-Plap%sita_ref(1,4))&
					/Plap%sita(1,4)))*temrate*dtime

					dkappa2=0.5*Plap%burgers(1,4)*inter_ref*Plap%shmu(1,4)/sqrt(dd_tot(is))&
					*(dfordsh+drevndsh+drevpdsh)*abs(gamdot(is))*dtime

					kappa(is)=kappa(is)+dkappa1+dkappa2

					dd_for(is)=dd_for(is)+dfordsh*abs(gamdot(is))*dtime

					dd_revp(is)=dd_revp(is)+drevpdsh*abs(gamdot(is))*dtime

					dd_revn(is)=dd_revn(is)+drevndsh*abs(gamdot(is))*dtime

					ddebdsh=Plap%qsub(1,4)*Plap%burgers(1,4)*sqrt(dd_deb_tot)*k2*dd_for(is)
	
					dkappa3_s=-0.5*kdeb*Plap%shmu(1,4)*Plap%burgers(1,4)/sqrt(dd_deb_tot)*&
					(log(Plap%burgers(1,4)*sqrt(dd_deb_tot))+1.0)*ddebdsh*abs(gamdot(is))*dtime

					dkappa3=dkappa3+dkappa3_s

					dd_deb(is)=dd_deb(is)+ddebdsh*abs(gamdot(is))*dtime

				else
					call RunTimeError(File_o, 'Hardening rule: is exceed maxinum number')
				endif
					
		    if (debug==1) then
				write(file_e,*)  'is =', is
				write(file_e,*)  'k2 =', k2
     		    write(file_e,*) 'dfordsh, drevpdsh,drevndsh,ddebdsh='
     		    write(file_e,*) dfordsh, drevpdsh,drevndsh,ddebdsh
     		    write(file_e,*) 'dkappa1, dkappa2,dkappa3='
     		    write(file_e,*) dkappa1,dkappa2,dkappa3_s
		    endif
			enddo
		
			do is= 1, numslip
			   if (debug==1) then
				   write(file_e,*)  'is =', is
				   write(file_e,*) 'dkappa, akappa3,ratio='
				   write(file_e,*) dkappa1+dkappa2+dkappa3,dkappa3,(dkappa1+dkappa2+dkappa3)/kappa(is)
				   kappa(is)=kappa(is)+dkappa3
			   endif
			end do

		else if (Plap%crystalID(grainid) .eq. kBCC) then
			numslip=PhSlip(2)
			
			do is= 1, numslip
				erate=dabs(gamdot(is))
				k2=Plap%k1(2,1)*Plap%burgers(2,1)*inter_ref/Plap%gref(2,1)*(1.-Plap%kbolt*tem		&
				*log(erate/rate_ref)/(Plap%Drag(2,1))/Plap%burgers(2,1)**3.0)

				dfordsh=(1.-pcoef)*Plap%k1(2,1)*sqrt(dd_tot(is))-k2*dd_for(is)

				if (stressf(is) .ge. 0.0d0) then
					drevpdsh=pcoef*Plap%k1(2,1)*sqrt(dd_tot(is))-k2*dd_revp(is)
					drevndsh=-Plap%k1(2,1)*sqrt(dd_tot(is))*(dd_revn(is)/dd_rev0(is))

				else
					drevndsh=pcoef*Plap%k1(2,1)*sqrt(dd_tot(is))-k2*dd_revn(is)
					drevpdsh=-Plap%k1(2,1)*sqrt(dd_tot(is))*(dd_revp(is)/dd_rev0(is))

				end if

				dkappa1=Plap%sref(2,1)/Plap%sita(2,1)*(1.-exp((tem-Plap%sita_ref(2,1))&
				/Plap%sita(2,1)))*temrate*dtime

				dkappa2=0.5*Plap%burgers(2,1)*inter_ref*Plap%shmu(2,1)/sqrt(dd_tot(is))&
				*(dfordsh+drevndsh+drevpdsh)*abs(gamdot(is))*dtime

				kappa(is)=kappa(is)+dkappa1+dkappa2

				dd_for(is)=dd_for(is)+dfordsh*abs(gamdot(is))*dtime

				dd_revp(is)=dd_revp(is)+drevpdsh*abs(gamdot(is))*dtime

				dd_revn(is)=dd_revn(is)+drevndsh*abs(gamdot(is))*dtime

				ddebdsh=Plap%qsub(2,1)*Plap%burgers(2,1)*sqrt(dd_deb_tot)*k2*dd_for(is)
	
				dkappa3_s=-0.5*kdeb*Plap%shmu(2,1)*Plap%burgers(2,1)/sqrt(dd_deb_tot)*&
				(log(Plap%burgers(2,1)*sqrt(dd_deb_tot))+1.0)*ddebdsh*abs(gamdot(is))*dtime

				dkappa3=dkappa3+dkappa3_s

				dd_deb(is)=dd_deb(is)+ddebdsh*abs(gamdot(is))*dtime
				
				if (debug==1) then
					write(file_e,*)  'is =', is
					write(file_e,*)  'k2 =', k2
					write(file_e,*) 'dfordsh, drevpdsh,drevndsh,ddebdsh='
					write(file_e,*) dfordsh, drevpdsh,drevndsh,ddebdsh
					write(file_e,*) 'dkappa1, dkappa2,dkappa3='
					write(file_e,*) dkappa1,dkappa2,dkappa3_s
				end if
			end do

			do is= 1, numslip
				if (debug==1) then
					write(file_e,*)  'is =', is
					write(file_e,*) 'dkappa, akappa3,ratio='
					write(file_e,*) dkappa1+dkappa2+dkappa3,dkappa3,(dkappa1+dkappa2+dkappa3)/kappa(is)
				end if
			    kappa(is)=kappa(is)+dkappa3
			end do

		else 
			call RunTimeError( FILE_O, 'Error: crystalID exceed the existing lattice!')
	    end if	

    else
!	wrong code number
		call RunTimeError(FILE_O, 'IntegrateHardening: Wrong kInteg_Code!')

    endif

    return
END

!***********************************************************************
SUBROUTINE IntegrateCrystalEqns2(s_ij, stress, estran, kappa, mu,      &
		   mu_n, eqvalues, gamdot, rss, mudot, stress_n, estran_n,     &
		   kappa_n, dtime, iterCounterS, iterCounterN,  ierr,          &
		   iterState, tolerState, tolerNewt, cepmod, dstrain,  kinc,   &
		   kstep,  argmin, argmax, VecM, VecS,ZVec, dd_for_n, dd_for,  &
		   dd_revp_n, dd_revp, dd_revn_n, dd_revn, stressf_n, stressf, &
		   dd_rev0_n, dd_rev0, dd_deb_n, dd_deb, time, dd_tot, dd_max, &
		   noel, grainid)  
!***********************************************************************

    USE NumType
	USE FileIO
	USE Timing
	USE SlipGeo
	IMPLICIT NONE

	REAL(KIND=8) 	::  s_ij(NSD,NSD)
	REAL(KIND=8) 	::  Stress(NVEC),Stressini(NVEC), Stressvec(NVEC)
	REAL(KIND=8) 	::  estran(NVEC), kappa(MaxNumSlip)   
	INTEGER	     	::  iterCounterS, iterCounterN, ierr
	REAL(KIND=8) 	::  argmin, argmax, dtime, dstrain(NVEC)
	REAL(KIND=8) 	::  search, tolerState, tolerNewt
      
	REAL(KIND=8) 	::  dd_for_n   (MaxNumSlip)	
	REAL(KIND=8) 	::  dd_revp_n  (MaxNumSlip)	
	REAL(KIND=8) 	::  dd_revn_n  (MaxNumSlip)
	REAL(KIND=8) 	::  stressf_n  (MaxNumSlip)
	REAL(KIND=8) 	::  dd_rev0_n  (MaxNumSlip)
	REAL(KIND=8) 	::  dd_deb_n   (MaxNumSlip)

	REAL(KIND=8) 	::  dd_for   (MaxNumSlip)	
	REAL(KIND=8) 	::  dd_revp  (MaxNumSlip)	
	REAL(KIND=8) 	::  dd_revn  (MaxNumSlip)
	REAL(KIND=8) 	::  stressf  (MaxNumSlip)
	REAL(KIND=8) 	::  dd_rev0  (MaxNumSlip)
	REAL(KIND=8) 	::  dd_deb   (MaxNumSlip)
	REAL(KIND=8) 	::  dd_tot  
	REAL(KIND=8) 	::  R_G_all 
	REAL(KIND=8) 	::  dd_min, dd_max
	REAL(KIND=8) 	::  R_Ghosh
	REAL(KIND=8) 	::  time

	REAL(KIND=8) 	::  gamdot(MaxNumSlip), mu(NVEC), eqvalues(NEQVA)
	REAL(KIND=8) 	::  stress_n(NVEC), estran_n(NVEC), kappa_n(MaxNumSlip), mu_n(NVEC)

	REAL(KIND=8) 	::  cepmod(NVEC,NVEC)
	REAL(KIND=8) 	::  rss(MaxNumSlip),  mudot(NVEC)


	LOGICAL			::  converged, NConverged
	REAL(KIND=8) 	::  norm_stress0, norm_kappa0, norm_stress, norm_kappa, normtau, normk, normtau0, normk0
	REAL(KIND=8) 	::  vec1(NVEC), vec2(MaxNumSlip)

	REAL(KIND=8) 	::  delstress(NVEC), stress0(NVEC)
	REAL(KIND=8) 	::  delstressvec(NVEC)
    

	REAL(KIND=8) 	::  rhs(NVEC), lhs(NVEC,NVEC)
	REAL(KIND=8) 	::  rhs_norm0, rhs_norm
    
	REAL(KIND=8) 	::  ddd

	REAL(KIND=8) 	::  T_n, T_eff
	REAL(KIND=8) 	::  T_tot(NSD),T_t(NSD)

    REAL(KIND=8) 	::  VecMT(NSD)

    INTEGER         ::  c10, c20

	LOGICAL 		::  ConvergeState 
	REAL(KIND=8)    ::  InnerProductVec
    
	INTEGER 		::	is, ip, iterState, kinc, INFO, kstep, noel, grainid

    REAL(KIND=8) 	::  VecM(NSD,MaxNumSlip)
    REAL(KIND=8) 	::  VecS(NSD,MaxNumSlip)
    REAL(KIND=8) 	::  ZVec(NVEC,MaxNumSlip)
	REAL(KIND=8)  	::  start_time, end_time
    integer 		::  iflag
    CHARACTER(LEN=20) :: ErroInfo   
    
!-----------------------------------------------------------------------
      
    if (debug==1) write(FILE_E, *) '*----Begin IntegrateCrystalEqns!----*'
 
! 	initialize 
	stress=Stress_n
    kappa=kappa_n      
	dd_for=dd_for_n
	dd_revp=dd_revp_n
	dd_revn=dd_revn_n
	stressf=stressf_n
	dd_rev0=dd_rev0_n
	dd_deb=dd_deb_n
          
! 	initial stress      
    stressini = stress 

!	echo stress and slip system strength at the beginning
    if (debug==1) then
        write(FILE_E, *) 'Initial estimate of stress is:'
        write(FILE_E, '(6e24.16)') , stress(:)
        write(FILE_E, *) 'Initial estimate of kappa is:'
        write(FILE_E, '(6e24.16)') , kappa(:)      
    endif    
                   
!	echo data of dislocation density for debugging
	if (debug==1) then   
		write(FILE_E,*) 'ip, dd_for:'
        write(FILE_E,*)  ip,dd_for(:)  
        write(FILE_E,*) 'dd_revp:'
        write(FILE_E,*)  dd_revp(:)    
        write(FILE_E,*) 'dd_revn:'
        write(FILE_E,*)  dd_revn(:)  
        write(FILE_E,*) 'stressf:'
        write(FILE_E,*)  stressf(:)  
        write(FILE_E,*) 'dd_rev0:'
        write(FILE_E,*)  dd_rev0(:)    
        write(FILE_E,*) 'dd_deb:'
        write(FILE_E,*)  dd_deb(:)  
    endif
        
!	initial valuies for the norm of separation, stress and strength	
    normtau0 = norm2(stress) 
    normk0   = norm2(kappa)    

!	initialized global flag to monitor Newton/State convergence
    ierr = XTAL_CONVERGED

! 	iterate for the material state
    iterState    = 0
    converged = .false.

    do while(iterState .lt. iterCounterS  .and. .not.converged )
		iflag = 0
		           
        iterState = iterState + 1

!		solve for stress  
		call ComputeStress(rhs, stress, stressini, estran, rss, mudot,     &
		   gamdot, kappa, dtime, argmin, argmax, dstrain, VecM, VecS,  &
		   ZVec, iterCounterN, tolerNewt, kstep, kinc, ierr, noel, grainid) 

        if (ierr .ne. XTAL_CONVERGED) return
        normtau= norm2(stress)

!		solve for hardness
        call IntegrateHardening(kappa, kappa_n, eqvalues, mu, mu_n,      &
		   estran, estran_n, stress, rss, mudot, gamdot, dtime, 	   &
		   kHARD_DD, VecM, VecS,ZVec, argmin, argmax, dd_for,       &
		   dd_revp, dd_revn, stressf, dd_rev0, dd_deb, time, dstrain,  &
		   dd_tot, dd_max, noel, grainid)                                         
        normk= norm2(kappa)     

!		check convergence    
        converged= convergestate( normtau, normk, normtau0, normk0, tolerState)
        if (debug==1) write(file_E, *) 'Converged= ', converged
        if (.not.converged) then 
            normtau0 = normtau
            normk0= normk 
        endif
        
!		echo data during the iteration
		if (noel==37 .and. debug==1) then
			write(File_e,*) 'kinc = ', kinc 
			write(FILE_E, *) 'stress is:'
			write(FILE_E, '(6e24.16)') , stress(:)
			write(FILE_E, *) 'kappa is:'
			write(FILE_E, '(6e24.16)') , kappa(:)      
		endif    
		      
    enddo
	
!	echo data for debugging when convergence is achieved
    if (converged .and. debug==1) then
		write(file_e,*) 'stress is'
		write(file_e,*) stress(:)
    end if

	call Vec6x1ToMat3x3Symm(stress, s_ij, NSD)
	estran(:) = estran_n(:) + mudot*dtime ! plastic strain?   

!	keep track of state iteration and check number of iterations
    if (iterState .ge. iterCounterS) then
        Write(FILE_O,*) 'Stress and Strength Solve: iters > maxIters for element', noel
        ierr = XTAL_MAX_ITERS_HIT
        return
    endif

!	CONSISTENT TANGENT
    call system_clock(c10) 
    call PlasticModuli(cepmod, stress, kappa, dtime, VecM, VecS, ZVec, grainid, noel) 
    call system_clock(c20)    
    tsj=tsj+ (c20-c10)/REAL(cr)                 

!	Ghosh Crack nucleation Metrics
    R_Ghosh=0.0d0 
	VecMT(:)=SlipG%VecM0(:,1,grainid)
	T_tot(:)=MATMUL(VecMT(:),s_ij(:,:))
    T_n=dot_product(T_tot(:),SlipG%VecM0(:,1,grainid))
    T_t(:)=T_tot(:)-T_n*VecMT(:)
    T_eff=sqrt(T_n**2.0d0+norm2(T_t(:))**2.0d0)
   ! R_G_all=T_eff*VolFrac**(1.0d0/3.0d0)
	R_G_all=T_eff**(1.0d0/3.0d0)
    if (R_G_all .ge.R_Ghosh) then
		R_Ghosh=R_G_all
    end if

    if (debug==1) write(FILE_E,*) '----- End IntegrateCrystalEqns.f90' 
    
    return
      
END

					     
!******************************************
logical FUNCTION Converged2(res, toler, n)
!******************************************
    use NumType
    implicit none

    real(kind=8)		:: toler
    real(kind=8)		:: res(n)

    integer				:: i,n

!-----------------------------------------------------------------------
  
!	check convergence on residual
    Converged2 = ( dabs(res(1)) .lt. toler )
    do i = 2, n
        Converged2 = ( (dabs(res(i)) .lt. toler) .and. Converged2)
    enddo

    return
END

!**********************************************************************
logical FUNCTION ConvergeState(norm_a, norm_b, norm_a0, norm_b0, toler)
!**********************************************************************
      
      use FILEIO
      use NumType
      implicit none

      real(kind=8)	:: norm_a, norm_b,norm_a0, norm_b0
      real(kind=8)  :: toler
!
!---------------------------------------------------------------------72
!
      ConvergeState = (dabs(norm_a - norm_a0) .lt. toler*norm_a0) .and. & 
                      (dabs(norm_b - norm_b0) .lt. toler*norm_b0)
	  
	  if (.not. convergeState .and. debug == 1) then
	      write(file_e,*) dabs(norm_a - norm_a0)-toler*norm_a0
	      write(file_e,*) '*********************************************'
	      write(file_e,*) dabs(norm_b - norm_b0)-toler*norm_b0
	      write(file_e,*) '*********************************************'	           
	  end if
	  
      if (.not.ConvergeState) then
          norm_a0 = norm_a
          norm_b0 = norm_b
      endif

      return
END

!************************************
SUBROUTINE WriteWarning(io, message)
!************************************

	implicit none
	character*(*)		:: message
    integer 			:: io

!-----------------------------------------------------------------------

    write(io, 1000) message

1000  format(/,'***WARNING Message: '/, 3x, a)

     return
END

