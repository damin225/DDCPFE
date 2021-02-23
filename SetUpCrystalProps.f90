!********************************
SUBROUTINE SetUpCrystalProps( )
!********************************      

    implicit none
!-----------------------------------------------------------------------

!	Open files
    CALL CrystalOpenIOFiles( )

!	Read in neighboring grain info
    CALL Read_NBgrain( )

!	Read in data
    CALL CrystalModelData( )

!	Read in Orientation data
    CALL CrystalOrientationData( )     

!	Read in elastic moduli
	CALL Read_ElasMod( )
	
!	Initialize some of the arrays                                                  	
    CALL CrystalInitializeArrays( )

!	parameters for iterations
    CALL CrystalSolverData( )

    return
END 


!***************************
subroutine Read_NBgrain()
!***************************

    use NumType
    use workdir
    use fileIO
    use NbData

	implicit none

    integer(kind=IKind) :: ICTR, I, K1, K2, J, IPH1, I1, ILine, TotalLine, ip
    integer(kind=IKind) :: lenoutdir, lenCoefFile, NPartctf, NumIntctf, Reason,NumNb, Nbtmp(Npart), NumCoef!number of coefficient tensors
	character*255       :: filename
	real(kind=RKInd)	:: tempR(nvec, nsd)

!-----------------------------------------------------------------------


!	neighboring info input
    DO ip = 1, NPart 
        READ(NbDataFile, *) Nbtmp
        NumNb=minloc(Nbtmp,1)-1
        if (.not. allocated(Elset_nb(ip)%NBGID))  then
			ALLOCATE( Elset_nb(ip)%NBGID(NumNb))
			do i = 1,NumNb
				Elset_nb(ip)%NBGID(i)=Nbtmp(i)
            end do
            Elset_nb(ip)%NNG=NumNb
        endif
    END DO

!	echo neighboring info for debugging
    if (debug==1)  then
        write(FILE_E,*)   '-----Neighboring Grain info----- ' 
        DO ip = 1, NPart
            write(File_E, '(1000I10)') Elset_nb(ip)%NBGID
            write(*,*) 
       END DO
    endif

	return
end subroutine Read_NBgrain

!******************************
SUBROUTINE CrystalOpenIOFiles()
!******************************

    use NumType
    use WorkDir
    use FileIO
    implicit none

    character :: filename1*255, filename2*255, filename3*255,          &
				 filename4*255, filename5*255, filename6*255,          &
				 filename7*255, filename8*255, filename9*255

!-----------------------------------------------------------------------

    filename1 = adjustl(trim(outdir))//adjustl(trim(FilePre))//'.xtali'
    !filename1 = './'//adjustl(trim(FilePre))//'.xtali'
    open(unit=FILE_I, file=filename1, status='unknown',access='sequential')
    rewind(FILE_I)
!
    filename2 = adjustl(trim(outdir))//adjustl(trim(FilePre))//'.xtale'
    !filename2 = './'//adjustl(trim(FilePre))//'.xtale'
    open(unit=FILE_E, file=filename2, status='unknown')
    rewind(FILE_E)
    if (debug==0) write(File_E, *) 'debug=0, nothing written into this file!'
!
    filename3 = adjustl(trim(outdir))//adjustl(trim(FilePre))//'.xtalo'
    !filename3 = './'//adjustl(trim(FilePre))//'.xtalo'
    open(unit=FILE_O, file=filename3, status='unknown')
!
    filename4 = adjustl(trim(outdir))//adjustl(trim(textureFile))//'.txti'
    !filename4 = './'//adjustl(trim(textureFile))//'.txti'
    open(unit=TXT_I, file=filename4, status='unknown',access='sequential')
    rewind(TXT_I)
!    
    filename5 = adjustl(trim(outdir))//adjustl(trim(FilePre))//'.itero'
    !filename5 = './'//adjustl(trim(FilePre))//'.itero'
    open(unit=iter_O, file=filename5, status='unknown')
    rewind(iter_O) 
    if (iterprint==1) write(iter_O, '(3A12, A20)'), '----KStep----', '----KInc----', '----Iters----', '----errorinfo----'
!    
!	part number for each phase
    filename6 = adjustl(trim(outdir))//adjustl(trim(PhasePart))//'1.dat'
    !filename6 = './'//adjustl(trim(PhasePart))//'1.dat'
    open(unit=Ph1PartFile, file=filename6, status='unknown')
    rewind(Ph1PartFile) 
!
    filename7 = adjustl(trim(outdir))//adjustl(trim(PhasePart))//'2.dat'
    !filename7 = './'//adjustl(trim(PhasePart))//'2.dat'
    open(unit=Ph2PartFile, file=filename7, status='unknown')
    rewind(Ph2PartFile)   
!    
    filename8 = adjustl(trim(outdir))//'NBgrain.dat'
    !filename8 = './'//'NBgrain.dat'
    open(unit=NbDataFile, file=filename8, status='unknown')
    rewind(NbDataFile) 

!	echo file path
    if (debug==1) then 
        write(File_E, '(A15, 1x, A55)') 'xtalifile:',  filename1
        write(File_E, '(A15, 1x, A55)') 'xtalefile:',  filename2
        write(File_E, '(A15, 1x, A55)') 'xtalofile:',  filename3
        write(File_E, '(A15, 1x, A55)') 'texturefile:',filename4
        write(File_E, '(A15, 1x, A55)') 'iterofile:',  filename5        
    endif
    return
END


!*****************************
SUBROUTINE CrystalModelData()
!*****************************         
    use numtype
    use FileIO
    USE WorkDir
    use PlaPar
    use PowerLawBoundPar
      
    implicit none
    integer(kind=ikind) :: NPhin,NPartin, CrystalIDin(NPh), iPart, IPh
    integer(kind=ikind) :: iline, PartID,numslip

!-----------------------------------------------------------------------

!   check the correctness of max number of slip system and number of 
!	parts per phase defined in the module file
    if (maxnumslip .ne. maxval(PhSlip)) call RunTimeError(File_o, 'maxnumslip not equat to  max(PhSlip)')
    if (NPart .ne. sum(ParPerPh))  call RunTimeError(File_o, 'Numpart  is not equal  to sum(ParPerPh) ') 
    
!	read number of phases, parts and interface partitions
    read(FILE_I, *) NPhin, NPartin
    if (NPart .ne. NPartin)   call RunTimeError(File_o, 'Number of parts from test.xtali is not equal from that the module') 
    if (NPh .ne. NPhin)  	  call RunTimeError(File_o, 'Number of phases from test.xtali is not equal from that the module') 
	 
!	echo some data 
    if (debug==1) write(FILE_E, '(A20, 1x, I8, 1x, A20, 1x, I8)')     &
                        'Phase number:',    NPhin,                    &
                        'NPart from xtali', NPartin

!	read cyrstal ID
    read(FILE_I, *) CrystalIDin
        
!	read the partid for phase 1 
	Plap%CrystalID=0   
    do iline=1,ParPerPh(1)
        read(Ph1PartFile, *) PartID
        if(Plap%CrystalID(PartID) .ne. 0) call RunTimeError(File_o, ' Plap%PhID is not initialized to be zero! ')    
        Plap%CrystalID(PartID)= CrystalIDin(1)
    enddo       
    close(Ph1PartFile)
    
!	read the partid for phase 2   
    do iline=1,ParPerPh(2)
        read(Ph2PartFile, *) PartID
        if(Plap%CrystalID(PartID) .ne. 0) call RunTimeError(File_o, ' There is overlap between ph1.dat and ph2.dat ')          
        Plap%CrystalID(PartID)= CrystalIDin(2)
    enddo       
    close(Ph2PartFile)

!	read material parameters
    read(FILE_I, *)      Plap%h0, Plap%xm, Plap%gam0, Plap%tausi,      &
                         Plap%taus0, Plap%xms, Plap%gamss0, Plap%crss0 
    read(FILE_I, *) Plap%kbolt
    do iPh=1, NPH  
		read(FILE_I, *) Plap%detF(iPh,1:Slipsystype(iPh))
		read(FILE_I, *) Plap%detV(iPh,1:Slipsystype(iPh))
		read(FILE_I, *) Plap%ddavg(iPh,1:Slipsystype(iPh))
		read(FILE_I, *) Plap%vid(iPh,1:Slipsystype(iPh))
		read(FILE_I, *) Plap%burgers(iPh,1:Slipsystype(iPh))
		read(FILE_I, *) Plap%sref(iPh,1:Slipsystype(iPh))     !MPa-->Pa
		read(FILE_I, *) Plap%sita(iPh,1:Slipsystype(iPh))     
		read(FILE_I, *) Plap%sita_ref(iPh,1:Slipsystype(iPh))
		read(FILE_I, *) Plap%s0ref(iPh,1:Slipsystype(iPh))    !MPa-->Pa
		read(FILE_I, *) Plap%k1(iPh,1:Slipsystype(iPh))
		read(FILE_I, *) Plap%gref(iPh,1:Slipsystype(iPh))
		read(FILE_I, *) Plap%Drag(iPh,1:Slipsystype(iPh))     !MPa-->Pa
		read(FILE_I, *) Plap%Href(iPh,1:Slipsystype(iPh))
		read(FILE_I, *) Plap%qsub(iPh,1:Slipsystype(iPh))
		read(FILE_I, *)	Plap%shmu(iPh,1:Slipsystype(iPh))     !MPa-->Pa
		read(FILE_I, *) Plap%crss0in(iPh,1:Slipsystype(iPh))  !MPa-->Pa
    enddo
    
!	echo the plastic parameters
    if (debug==1) write(FILE_E, 1000)  Plap%CrystalID(1), NPhin, NPart,& 
									   Plap%h0, Plap%xm, Plap%gam0, Plap%tausi,      &
						               Plap%taus0, Plap%xms, Plap%gamss0, Plap%crss0

    if (debug==1) write(FILE_E, 1002)  Plap%CrystalID(1),Plap%detF(1,1), Plap%detV(1,1),&
									   Plap%ddavg(1,1),Plap%vid(1,1),Plap%burgers(1,1),Plap%sref(1,1),Plap%sita(1,1),Plap%sita_ref(1,1),&
									   Plap%s0ref(1,1),Plap%k1(1,1), Plap%gref(1,1),  Plap%Drag(1,1),Plap%Href(1,1),   Plap%qsub(1,1),  &
			                           Plap%shmu(1,1), Plap%crss0in(1,1)

    if (debug==1) write(FILE_E, 1002)  Plap%CrystalID(1),Plap%detF(1,2), Plap%detV(1,2),&
									   Plap%ddavg(1,2),Plap%vid(1,2),Plap%burgers(1,2),Plap%sref(1,2),Plap%sita(1,2),Plap%sita_ref(1,2),&
									   Plap%s0ref(1,2),Plap%k1(1,2), Plap%gref(1,2),  Plap%Drag(1,2),Plap%Href(1,2),   Plap%qsub(1,2),  &
									   Plap%shmu(1,2), Plap%crss0in(1,2)

    if (debug==1) write(FILE_E, 1003)  Plap%CrystalID(2),Plap%detF(2,1), Plap%detV(2,1),&
									   Plap%ddavg(2,1),Plap%vid(2,1),Plap%burgers(2,1),Plap%sref(2,1),Plap%sita(2,1),Plap%sita_ref(2,1),&
									   Plap%s0ref(2,1),Plap%k1(2,1), Plap%gref(2,1),  Plap%Drag(2,1),Plap%Href(2,1),   Plap%qsub(2,1),  &
								       Plap%shmu(2,1), Plap%crss0in(2,1)

!	bounds for power law
    CALL BoundForArgPowLaw(Plap%xm, PLBP%argmin, PLBP%argmax)
    
!	set up slip systems
    do ipart=1,NPart
        if (Plap%crystalID(iPart) .eq. kFCC) then
            call SetSlipSystemFCC( iPart )
			Plap%kappa0(:, iPart)=Plap%crss0in(2,1)
        else if (Plap%crystalID(iPart) .eq. kBCC) then
			call SetSlipSystemBCC( iPart )
			numslip=PhSlip(2)
			Plap%kappa0(1:numslip, iPart)=Plap%crss0in(2,1)
		else if (Plap%crystalID(iPart) .eq. kHCP) then
			call SetSlipSystemHCP( iPart )
			Plap%kappa0(1:3, iPart)=Plap%crss0in(1,1)
			Plap%kappa0(4:6, iPart)=Plap%crss0in(1,2)
			Plap%kappa0(7:12, iPart)=Plap%crss0in(1,3)
			Plap%kappa0(13:30, iPart)=Plap%crss0in(1,4)	
		else	    
            call RunTimeError( FILE_O,'Error: Wrong lattice code number!')
        endif
    enddo

! 	formats
1000  format(/'*-----   Crystal Plasticity parameters -----*'/,        &
              7x,'  elasticity type  = ', i5/,                         &
              7x,'  number of phase  = ', i5/,                         &
              7x,'  number of parts  = ', i5/,                         &
              7x,'  h0               = ', e12.5/,                      & 
              7x,'  xm               = ', e12.5/,                      & 
              7x,'  gam0             = ', e12.5/,                      & 
              7x,'  tausi            = ', e12.5/,                      & 
              7x,'  taus0            = ', e12.5/,                      & 
              7x,'  xms              = ', e12.5/,                      & 
              7x,'  gamss0           = ', e12.5/,                      &               
              7x,'  crss0            = ', e12.5/ )

1001  format(/'*-----   Crystal Plasticity parameters of phase(2)-----*'/,        &
              7x,'  crystal type     = ', i5/,                         &
              7x,'  boltz 		     = ', e12.5/,                      &
              7x,'  detF             = ', e12.5/,                      & 
              7x,'  detV             = ', e12.5/,                      & 
              7x,'  ddavg            = ', e12.5/,                      & 
              7x,'  vid              = ', e12.5/,                      & 
              7x,'  burger           = ', e12.5/,                      & 
              7x,'  sref             = ', e12.5/,                      & 
              7x,'  sita             = ', e12.5/,                      & 
              7x,'  sita_ref         = ', e12.5/,                      & 
              7x,'  s0ref            = ', e12.5/,                      & 
              7x,'  k1               = ', e12.5/,                      & 
              7x,'  gref             = ', e12.5/,                      & 
              7x,'  Drag             = ', e12.5/,                      & 
              7x,'  Href             = ', e12.5/,                      & 
              7x,'  qsub             = ', e12.5/,                      & 
              7x,'  shmu             = ', e12.5/,                      & 
              7x,'  crss0in          = ', e12.5/)
              
1002  format(/'*-----   Crystal Plasticity parameters of phase(1)-----*'/,        &
              7x,'  crystal type     = ', i5/,                         &
              7x,'  detF             = ', e12.5/,                      & 
              7x,'  detV             = ', e12.5/,                      & 
              7x,'  ddavg            = ', e12.5/,                      & 
              7x,'  vid              = ', e12.5/,                      & 
              7x,'  burger           = ', e12.5/,                      & 
              7x,'  sref             = ', e12.5/,                      & 
              7x,'  sita             = ', e12.5/,                      & 
              7x,'  sita_ref         = ', e12.5/,                      & 
              7x,'  s0ref            = ', e12.5/,                      & 
              7x,'  k1               = ', e12.5/,                      & 
              7x,'  gref             = ', e12.5/,                      & 
              7x,'  Drag             = ', e12.5/,                      & 
              7x,'  Href             = ', e12.5/,                      & 
              7x,'  qsub             = ', e12.5/,                      & 
              7x,'  shmu             = ', e12.5/,                      & 
              7x,'  crss0in          = ', e12.5/)
              
1003  format(/'*-----   Crystal Plasticity parameters of phase(2)-----*'/,        &
              7x,'  crystal type  	 = ', i5/,                         &
              7x,'  detF             = ', e12.5/,                      & 
              7x,'  detV             = ', e12.5/,                      & 
              7x,'  ddavg            = ', e12.5/,                      & 
              7x,'  vid              = ', e12.5/,                      & 
              7x,'  burger           = ', e12.5/,                      & 
              7x,'  sref             = ', e12.5/,                      & 
              7x,'  sita             = ', e12.5/,                      & 
              7x,'  sita_ref         = ', e12.5/,                      & 
              7x,'  s0ref            = ', e12.5/,                      & 
              7x,'  k1               = ', e12.5/,                      & 
              7x,'  gref             = ', e12.5/,                      & 
              7x,'  Drag             = ', e12.5/,                      & 
              7x,'  Href             = ', e12.5/,                      & 
              7x,'  qsub             = ', e12.5/,                      & 
              7x,'  shmu             = ', e12.5/,                      & 
              7x,'  crss0in          = ', e12.5/)

      return
END

!************************************
SUBROUTINE CrystalOrientationData( )
!************************************
    USE NumType
	USE FILEIO
	USE WorkDir
	USE OriPar
	USE SlipGeo

	IMPLICIT none
   
	INTEGER          iPart, iSlip, i
	REAL(KIND=8)  :: pi180, EulerMax
	REAL(KIND=8)  :: sps, cps, sth, cth, sph, cph     

	INTEGER    NumEulerAngle, iidr, iikc
	CHARACTER  filename*255
!-----------------------------------------------------------------------  

!	open texture file
    filename = adjustl(trim(outdir))//adjustl(trim(textureFile))//'.txti'
    !filename = './'//adjustl(trim(textureFile))//'.txti'
    open(unit=TXT_I, file=filename, status='unknown', access='sequential')
  	rewind(TXT_I)

!	read number of angles
	READ (TXT_I, *)  NumEulerAngle
	IF (NumEulerAngle .NE. NPart)  CALL RunTimeError(FILE_O, 'AssignCrystalODF: NumEulerAngle != NumPart')
	
!	read the flag for angle convention and set some constants
!	iikc = 0 : angles input in Kocks convention :  (psi,the,phi)
!          1 : angles input in Canova convention : (ph,th,om)
!          	   ph = 90 + phi; th = the; om = 90 - psi
!	iidr = 0 : angles input in degrees
!		   1 : angles input in radians
	READ(TXT_I, *) iikc, iidr
	
!	Euler angles for each partitition
	pi180 = 4.0 * datan(1.0d+00)/180.
	IF (iidr .eq. 1) pi180=1.0         
	OriP%Euler=pzero  
	EulerMax=pzero     

    IF (debug==1) WRITE(FILE_E,*) 'Euler Angles are:'       
	DO iPart=1,NPart
		READ(TXT_I, *) OriP%Euler(:,iPart)
		OriP%Euler(:,iPart)=OriP%Euler(:,iPart)*pi180
		if (debug ==1) WRITE(FILE_E,'(3F12.5)')  OriP%Euler(:,iPart)
		DO i=1,3
			IF (EulerMax .LT. abs(OriP%Euler(i,iPart))) EulerMax=abs(OriP%Euler(i,iPart)) 
		END DO
	ENDDO
	        
	IF (EulerMax .GT. 10)  &
	CALL RunTimeError(FILE_o,'Check the units of Euler angles to make sure they are in Rad') 
	CLOSE (TXT_I)
    
!	build rotation matrices C0: {x}_sm = [C0] {x}_cr    
	DO iPart=1,NPart
		CALL AnglesToRotMatrix(OriP%Euler(:,iPart),OriP%gcrot0(:,:,iPart), NSD)
!		rotate the slip system into global coordinate system
        DO ISlip=1,MaxNumSlip            
			SlipG%VecM0(:,ISlip, IPart)=MATMUL(OriP%gcrot0(:,:,IPart),SlipG%VecM0(:,ISlip, IPart))
			SlipG%VecS0(:,ISlip, IPart)=MATMUL(OriP%gcrot0(:,:,IPart),SlipG%VecS0(:,ISlip, IPart))
            CALL OuterProductVec(SlipG%VecS0(:,ISlip, IPart), SlipG%VecM0(:,ISlip, IPart), SlipG%ZTen0(:,:,ISlip,IPart), NSD)                              
            CALL ZTenToVec(SlipG%ZTen0(:,:,ISlip,IPart), SlipG%ZVec0(:, ISlip,IPart),NSD)	
        ENDDO                        
                  
		if (debug==1 .and. iPart==1) then
		    WRITE(FILE_E,'(A15, 1x, I5/, 3(5x, F12.5))') 'gcrot0-part-', IPart, OriP%gcrot0(:,:,iPart)
            write(FILE_E,*)   '-----Slip system in global coordinates----- ' 
            write(FILE_E,*)   'vecS0 of part 1 is:' 
            write(FILE_E,'(3F10.4)')  SlipG%vecS0(:,:,1)
            write(FILE_E,*)   'vecM0 of part 1 is:' 
            write(FILE_E,'(3F10.4)')  SlipG%vecM0(:,:,1)   
            write(FILE_E,*)   'ZVec0 of part 1 is:' 
            write(FILE_E,'(6E12.5)')  SlipG%ZVec0(:,:,1)      
        endif
	ENDDO 

RETURN
END

!************************************
SUBROUTINE CrystalInitializeArrays()         
!************************************         
      use numtype
      use CPVars
      use CPVars_n
      use oriPar
      use SlipGeo
      use PlaPar
    
      implicit none
      integer(kind=ikind) iPart

!-----------------------------------------------------------------------

!	initialize arrays used in XTAL constitutive integration method
	CPV%     gstress     =   pzero
    CPV%     gestran     =   pzero
    CPV%     gkappa      =   pzero
    CPV_n%   gstress_n   =   pzero
    CPV_n%   gestran_n   =   pzero
    CPV_n%   gkappa_n    =   pzero
    CPV%     gmu         =   pzero
    CPV_n%   gmu_n       =   pzero
    CPV%     ggamdot     =   pzero
    
	CPV%     gdd_for     =   1.0d12
    CPV_n%   gdd_for_n   =   1.0d12
	CPV%     gdd_revp    =   pzero
    CPV_n%   gdd_revp_n  =   pzero
	CPV%     gdd_revn    =   pzero
    CPV_n%   gdd_revn_n  =   pzero
    CPV%     gdd_stressf =   pzero
    CPV_n%   gdd_stressf_n=   pzero
    CPV%     gdd_rev0    =   1.0d12
    CPV_n%   gdd_rev0_n  =   1.0d12
    CPV%     gdd_deb     =   1.0d10
    CPV_n%   gdd_deb_n   =   1.0d10

    do iPart = 1, NPart
        CPV%gkappa(:, iPart) =Plap%kappa0(:, IPart)
        CPV_n%gkappa_n(:, iPart) =Plap%kappa0(:, IPart)
    enddo
    return
END

!******************************
SUBROUTINE CrystalSolverData()
!******************************
    use FILEIO
    use IterPar
    implicit none
    
!---------------------------------------------------------------------

!	number iterations and tolerance for state iterations
    read(FILE_I, *) iterP%maxIterState, iterP%tolerState

!	number iterations and tolerance for newton method
    read(FILE_I, *) iterP%maxIterNewt, iterP%tolerNewt

!	echo input data
    if (debug==1) then                                                     
		write(FILE_E, 1000) iterP%maxIterState, iterP%maxIterNewt, iterP%tolerState, iterP%tolerNewt
	end if

!	format
1000  format(/'*-----   Local Convergence Control Parameters  -----*'/,&
              7x,'  max iters State  = ',i5 /,                         &
              7x,'  max iters Newton = ',i5 /,                         &
              7x,'  tolerance State  = ',e12.5 /,                      &
              7x,'  tolerance Newton = ',e12.5)
    return
END
   
!***************************************   
SUBROUTINE SetSlipSystemBCC(IPart)
!***************************************
    use NumType
    use SlipGeo 
    use FileIO 
    use PlaPar
    
   	IMPLICIT NONE
 
	REAL(KIND=8) ::  sDOtm
	REAL(KIND=8) :: indexM(3,48), indexS(3,48)
	REAL(KIND=8) :: InnerProductVec
	INTEGER(kind=iKind)      ::  ISlip,IPart
      
!	Same as Marin's Order      
	data indexM /0.,  1.,  1.,  &
                 1.,  0.,  1.,  &
                 1., -1.,  0.,  &

                 0.,  1., -1.,  &
                 1.,  0.,  1.,  &
                 1.,  1.,  0.,  &
                                 
                 0.,  1.,  1.,  &
                 1.,  0., -1.,  &
                 1.,  1.,  0.,  &
                                 
                 0.,  1., -1.,  &
                 1.,  0., -1.,  &
                 1., -1.,  0.,  &

				-2.,  1.,  1.,  &
                 1., -2.,  1.,  &
                 1.,  1., -2.,  &
                                  
                 2.,  1.,  1.,  &
                -1., -2.,  1.,  &
                -1.,  1., -2.,  &
                                   
                -2., -1.,  1.,  &
                 1.,  2.,  1.,  &
                 1., -1., -2.,  &
                             
                 2., -1.,  1.,  &
                -1.,  2.,  1.,  &
                -1., -1., -2.,  &

				 3., -1., -2.,  &
                 3., -2., -1.,  &
                -1.,  3., -2.,  &
                                   
                -2.,  3., -1.,  &
                -1., -2.,  3.,  &
                -2., -1.,  3.,  &
                                    
                 3.,  1.,  2.,  &
                 3.,  2.,  1.,  &
                 1.,  3., -2.,  &
                                  
                 2.,  3., -1.,  &
                 1., -2.,  3.,  &
                 2., -1.,  3.,  &
				
				 3.,  1., -2.,  &
                 3.,  2., -1.,  &
                 1.,  3.,  2.,  &
                           
                 2.,  3.,  1.,  &
                -1.,  2.,  3.,  &
                -2.,  1.,  3.,  &
                                  
                 3., -1.,  2.,  &
                 3., -2.,  1.,  &
                -1.,  3.,  2.,  &
                                    
                -2.,  3.,  1.,  &
                 1.,  2.,  3.,  &
                 2.,  1.,  3./
       
	data indexS /1.,  1., -1.,  &
                 1.,  1., -1.,  &
                 1.,  1., -1.,  &

                 1., -1., -1.,  &
                 1., -1., -1.,  &
                 1., -1., -1.,  &
                                    
                 1., -1.,  1.,  &
                 1., -1.,  1.,  &
                 1., -1.,  1.,  &
                                  
                 1.,  1.,  1.,  &
                 1.,  1.,  1.,  &
                 1.,  1.,  1.,  &

				 1.,  1.,  1.,  &
                 1.,  1.,  1.,  &
                 1.,  1.,  1.,  &
                                   
                -1.,  1.,  1.,  &
                -1.,  1.,  1.,  &
                -1.,  1.,  1.,  &
                                    
                 1., -1.,  1.,  &
                 1., -1.,  1.,  &
                 1., -1.,  1.,  &
                                    
                -1.,  -1., 1.,  &
                -1.,  -1., 1.,  &
                -1.,  -1., 1.,  &
                                    
				 1.,  1.,  1.,  &
                 1.,  1.,  1.,  &
                 1.,  1.,  1.,  &
                                    
		         1.,  1.,  1.,  &
                 1.,  1.,  1.,  &
                 1.,  1.,  1.,  &
                                    
                -1.,  1.,  1.,  &
                -1.,  1.,  1.,  &
                -1.,  1.,  1.,  &
                                   
                -1.,  1.,  1.,  &
                -1.,  1.,  1.,  &
                -1.,  1.,  1.,  &
                                    
                 1., -1.,  1.,  &
                 1., -1.,  1.,  &
                 1., -1.,  1.,  &
                                
                 1., -1.,  1.,  &
                 1., -1.,  1.,  &
                 1., -1.,  1.,  &
                                 
                -1., -1., 1.,   &
                -1., -1., 1.,   &
                -1., -1., 1.,   &
                                   
                -1., -1., 1.,   &
                -1., -1., 1.,   &
                -1., -1., 1./ 
      
!---------------------------------------------------------------------

!	echo the slip systems
    if (debug==1) then
        write(FILE_E,*)  'indexS is:' 
        write(FILE_E,'(3F10.4)')  indexS
        write(FILE_E,*)  'indexM is:' 
        write(FILE_E,'(3F10.4)')  indexM     
    endif

!	slip system normals and slip directions: unit vectors for HCP, 
!	numslip=30 is by default, make changes if necessary
    DO ISlip = 1, PhSlip(2)		
		SlipG%VecS0(:, ISlip, IPart)=indexS(:,ISlip)/DSQRT(DOT_PRODUCT(indexS(:,ISlip), indexS(:,ISlip)))
		SlipG%VecM0(:, ISlip, IPart)=indexM(:,ISlip)/DSQRT(DOT_PRODUCT(indexM(:,ISlip), indexM(:,ISlip)))

!		check normality of vecS and vecM
		sDOtm = DOT_PRODUCT(SlipG%vecS0(:,ISlip,IPart),SlipG%vecM0(:,ISlip,IPart))
			                                       
		if (sDOtm .ge. 1.0d-3)  then
			write(file_e,*) 'ISlip,sDOTm=',ISlip,sDOTm
			CALL RunTimeError(FILE_O,'Normality of vecS and vecM is not satisfied!')
		end if

        CALL OuterProductVec(SlipG%VecS0(:,iSlip,iPart),           	   &
                             SlipG%VecM0(:,iSlip,iPart),               &
                             SlipG%zTen0(:,:,iSlip,iPart), NSD)	
        CALL ZTenToVec(SlipG%ZTen0(:,:,iSlip,iPart),                   &
                       SlipG%ZVec0(:,iSlip,iPart),NSD)		
	ENDDO  
       
!	echo the normilized slip system
      IF (debug==1) then
         write(FILE_E,*)   '-----Slip system in HCP----- ' 
         write(FILE_E,*)   'vecS0 of part 1 is:' 
         write(FILE_E,'(3F10.4)')  SlipG%vecS0(:,:,1)
         write(FILE_E,*)   'vecM0 of part 1 is:' 
         write(FILE_E,'(3F10.4)')  SlipG%vecM0(:,:,1)
     ENDIF

    RETURN
END

!**********************************
SUBROUTINE SetSlipSystemFCC(iPart)
!**********************************
    use NumType
    use SlipGeo 
    use FileIO 
    use PlaPar
    
   	IMPLICIT NONE
 
	REAL(KIND=8) :: sDOtm
	REAL(KIND=8) :: indexM(3,12), indexS(3,12)
	REAL(KIND=8) :: InnerProductVec
	INTEGER(kind=iKind)      :: iPart, ISlip
      
!	Same as Marin's Order      
    data indexM /1.,  1.,  1.,  &
                 1.,  1.,  1.,  &
                 1.,  1.,  1.,  &
                                   
                -1.,  1.,  1.,  &
                -1.,  1.,  1.,  &
                -1.,  1.,  1.,  &
                                  
                -1., -1.,  1.,  &
                -1., -1.,  1.,  &
                -1., -1.,  1.,  &
                                   
                 1., -1.,  1.,  &
                 1., -1.,  1.,  &
                 1., -1.,  1./
       
	data indexS /0.,  1., -1.,  &
                 1.,  0., -1.,  &
                 1., -1.,  0.,  &
                                   
                 0.,  1., -1.,  &
                 1.,  0.,  1.,  &
                 1.,  1.,  0.,  &
                                  
                 0.,  1.,  1.,  &
                 1.,  0.,  1.,  &
                 1., -1.,  0.,  &
                                    
                 0.,  1.,  1.,  &
                 1.,  0., -1.,  &
                 1.,  1.,  0./ 

!-----------------------------------------------------------------------

!	echo the slip systems
    if (debug==1) then
        write(FILE_E,*)  'indexS is:' 
        write(FILE_E,'(3F10.4)')  indexS
        write(FILE_E,*)  'indexM is:' 
        write(FILE_E,'(3F10.4)')  indexM     
    endif

!	slip system normals and slip directions: unit vectors for FCC, 
!	numslip=12 is by default, make changes IF necessary.

    if (debug==1) write(FILE_E,*) 'Plap%kappa0(:,IPart) is assigned!' 

    DO ISlip = 1, MaxNumSlip
		SlipG%VecS0(:, ISlip, IPart)=indexS(:,ISlip)/DSQRT         &
			            (DOT_PRODUCT(indexS(:,ISlip), indexS(:,ISlip)))
		SlipG%VecM0(:, ISlip, IPart)=indexM(:,ISlip)/DSQRT         &
			            (DOT_PRODUCT(indexM(:,ISlip), indexM(:,ISlip)))

!		check normality of vecS and vecM
		sDOtm = DOT_PRODUCT(SlipG%vecS0(:,ISlip,IPart),            &
			                SlipG%vecM0(:,ISlip,IPart))			                                       
		IF (sDOtm .ne. 0)  CALL RunTimeError(FILE_O,'Normality of vecS and vecM is not satisfied!')
		    			
		CALL OuterProductVec(SlipG%VecS0(:,iSlip,iPart),           &
                             SlipG%VecM0(:,iSlip,iPart),           &
                             SlipG%zTen0(:,:,iSlip,iPart), NSD)	
        CALL ZTenToVec(SlipG%ZTen0(:,:,iSlip,iPart),               &
                       SlipG%ZVec0(:,iSlip,iPart),NSD)		
	ENDDO
 
       
!	echo the normilized slip system
    IF (debug==1) then
		write(FILE_E,*)   '-----Slip system in FCC----- ' 
        write(FILE_E,*)   'vecS0 of part 1 is:' 
        write(FILE_E,'(3F10.4)')  SlipG%vecS0(:,:,1)
        write(FILE_E,*)   'vecM0 of part 1 is:' 
        write(FILE_E,'(3F10.4)')  SlipG%vecM0(:,:,1)
    ENDIF
          
    RETURN
END

!***********************************
SUBROUTINE SetSlipSystemHCP(IPart)
!**********************************
    use NumType
    use SlipGeo 
    use FileIO 
    use PlaPar
    
   	IMPLICIT NONE
 
	REAL(KIND=8) ::  sDOtm
	REAL(KIND=8) :: indexM(3,30), indexS(3,30)
	REAL(KIND=8) :: InnerProductVec
	INTEGER(kind=iKind)      ::  ISlip,IPart

!	Same as Marin's Order      
     data indexM   /0.,  0.,  1.,   &
                    0.,  0.,  1.,   &
                    0.,  0.,  1.,   &
				          
					0.5,   -0.866,  0.,   &
                    0.5,    0.866,  0.,   &
                    -1.,       0.,  0.,   &

                    0.439,  0.761,  0.478,   &
                    0.439, -0.761,  0.478,   &
                   -0.439,  0.761,  0.478,   &
                   -0.439, -0.761,  0.478,   &
                    0.8785,    0.,  0.478,   &		   
                   -0.8785,    0.,  0.478,   &

                    0.439,   0.761,    0.478,  &
                    0.439,   0.761,    0.478,  &
                    0.439,  -0.761,    0.478,  &
                    0.439,  -0.761,    0.478,  &
                   -0.439,   0.761,    0.478,  &
                   -0.439,   0.761,    0.478,  &
                   -0.439,  -0.761,    0.478,  &
                   -0.439,  -0.761,    0.478,  &
                   0.8785,      0.,    0.478,  &
                   0.8785,      0.,    0.478,  &
                  -0.8785,      0.,    0.478,  &
                  -0.8785,      0.,    0.478,  &

                    0.7331,  0.42325, 0.53239, &
                   -0.7331,  0.42325, 0.53239, &
                    0.7331, -0.42325, 0.53239, &
                   -0.7331, -0.42325, 0.53239, &
                        0.,   0.8465, 0.53239, &
                        0.,  -0.8465, 0.53239/
       
       data indexS /0.866,  0.5,  0.,  &
                    0.,     -1.,  0.,  &
                   -0.866,  0.5,  0.,  &

					0.866,  0.5,  0.,   &
				   -0.866,  0.5,  0.,   &
                    0.,  -1.,  0.,   &

                   -0.866,  0.5,  0.,  &
                    0.866,  0.5,  0.,  &
                   -0.866, -0.5,  0.,  &
                    0.866, -0.5,  0.,  &
                       0.,   1.,  0.,  &
                       0.,  -1.,  0.,  &

                   -0.461, -0.266,  0.8465,   &
                    0.,   -0.5324,  0.8465,   &
                    0.,    0.5324,  0.8465,   &
                   -0.461,  0.266,  0.8465,   &
                    0.461, -0.266,  0.8465,   &
                    0.,   -0.5324,  0.8465,   &
                    0.461,  0.266,  0.8465,   &
                    0.,    0.5324,  0.8465,   &
                   -0.461,  0.266,  0.8465,   &
                   -0.461, -0.266,  0.8465,   &
                    0.461, -0.266,  0.8465,   &
                    0.461,  0.266,  0.8465,   &

                 -0.46105, -0.2662,  0.8465,  &
                  0.46105, -0.2662,  0.8465,  &
                 -0.46105,  0.2662,  0.8465,  &
                  0.46105,  0.2662,  0.8465,  &
                       0., -0.5324,  0.8465,  &
                       0.,  0.5324,  0.8465/ 
      
!-----------------------------------------------------------------------

!	echo the slip systems
    if (debug==1) then
        write(FILE_E,*)  'indexS is:' 
        write(FILE_E,'(3F10.4)')  indexS
        write(FILE_E,*)  'indexM is:' 
        write(FILE_E,'(3F10.4)')  indexM     
    endif

!	slip system normals and slip directions: unit vectors for HCP, 
!	numslip=30 is by default, make changes IF necessary.

    DO ISlip = 1, PhSlip(1)
		SlipG%VecS0(:, ISlip, IPart)=indexS(:,ISlip)/DSQRT         &
			            (DOT_PRODUCT(indexS(:,ISlip), indexS(:,ISlip)))
		SlipG%VecM0(:, ISlip, IPart)=indexM(:,ISlip)/DSQRT         &
			            (DOT_PRODUCT(indexM(:,ISlip), indexM(:,ISlip)))

!		check normality of vecS and vecM
		sDOtm = DOT_PRODUCT(SlipG%vecS0(:,ISlip,IPart),            &
			                SlipG%vecM0(:,ISlip,IPart))
			                                       
		if (sDOtm .ge. 1.0d-3)  then
			write(file_e,*) 'ISlip,sDOTm=',ISlip,sDOTm
			CALL RunTimeError(FILE_O, 'Normality of vecS and vecM is not satisfied!')
		endif

        CALL OuterProductVec(SlipG%VecS0(:,iSlip,iPart),           &
                             SlipG%VecM0(:,iSlip,iPart),           &
                             SlipG%zTen0(:,:,iSlip,iPart), NSD)	
        CALL ZTenToVec(SlipG%ZTen0(:,:,iSlip,iPart),               &
                       SlipG%ZVec0(:,iSlip,iPart),NSD)		
   	    IF (debug==1) then
			write(file_e, *) 'ISlip,ipart ='
			write(file_e, *) ISlip, ipart
			write(file_e, *) 'ZVec0='
			write(file_e, *) SlipG%ZVec0(:,ISlip,iPart)
 	    ENDIF		
		
	ENDDO  

       
!	echo the normilized slip system
    IF (debug==1) then
		write(FILE_E,*)   '-----Slip system in HCP----- ' 
        write(FILE_E,*)   'vecS0 of part 1 is:' 
        write(FILE_E,'(3F10.4)')  SlipG%vecS0(:,:,1)
        write(FILE_E,*)   'vecM0 of part 1 is:' 
        write(FILE_E,'(3F10.4)')  SlipG%vecM0(:,:,1)
    ENDIF
           
    RETURN
END


!**************************
SUBROUTINE Read_ElasMod( )
!**************************
	
	use ElasMod
	use FILEIO
	use NumType
	use OriPar
	use Plapar
	use WorkDir
	
	implicit none
	character*255		:: filename
	real(kind=8)		:: elascoef1(5), elascoef2(3)
	real(kind=8)		:: ELASMOD0(6,6), ELASMOD1(6,6), ELASMOD2(6,6), elasModTmp(6,6), Lijkltemp(6,6)
		
	integer(kind=8) 	:: i, i1, i2, i3, i4, j1, j2, j3, j4, ii1, ii2, jj1, jj2, ip, mymod, switch(6), ii, jj
	real(kind=8)        :: sps, cps, sth, cth, sph, cph, crot(nsd,nsd), RotM6x6(6, 6) 
!-----------------------------------------------------------------------
	switch = (/1,2,3,6,5,4/)
	
!	open elasmod.dat file
	filename = adjustl(trim(outdir))//'elasmod.dat'
	!filename = './elasmod.dat'
    open(unit=225, file=filename, status='unknown', access='sequential')
    rewind(225)

!	read elastic coefficients: 1st line hcp; 2nd line bcc
	read(225,*) elascoef1(:)
	read(225,*) elascoef2(:)

!	construct elastic moduli
	ELASMOD1 = 0.d0
	ELASMOD2 = 0.d0

	! HCP
	DO i=1,2
		ELASMOD1(1,i+1) = elasCoef1(i+1)
		ELASMOD1(i+1,1) = elasCoef1(i+1)

		ELASMOD1(i,i) = elasCoef1(1)
		ELASMOD1(i+3,i+3) = elasCoef1(5)
	ENDDO
	ELASMOD1(2,3) = elasCoef1(3)
	ELASMOD1(3,2) = elasCoef1(3)
	ELASMOD1(6,6) = 0.5d0*(elasCoef1(1)-elasCoef1(2))
	ELASMOD1(3,3) = elasCoef1(4)
	
	!BCC
	DO i=1,3
		ELASMOD2(i,i) = elasCoef2(1)
		ELASMOD2(i+3,i+3) = elasCoef2(3)
	ENDDO
	ELASMOD2(1,2:3) = elasCoef2(2)
	ELASMOD2(2:3,1) = elasCoef2(2)
	ELASMOD2(3,2) = elasCoef2(2)
	ELASMOD2(2,3) = elasCoef2(2)

!	rotate elastic moduli: Bunge's convention

	do ip =1,NPart
		sps = dsin(OriP%Euler(1,ip))
		cps = dcos(OriP%Euler(1,ip))
		sth = dsin(OriP%Euler(2,ip))
		cth = dcos(OriP%Euler(2,ip))
		sph = dsin(OriP%Euler(3,ip))
		cph = dcos(OriP%Euler(3,ip))
		
		crot(1,1) = cph*cps - sps * sph * cth
		crot(1,2) =  -cps * sph - sps * cph * cth
		crot(1,3) =  sps * sth
		crot(2,1) =  sps * cph + cps * sph * cth
		crot(2,2) = -sps * sph + cps * cph * cth
		crot(2,3) =  -cps * sth
		crot(3,1) =  sph * sth
		crot(3,2) =  cph*sth
		crot(3,3) =  cth
		
		do JJ1=1,3
			J1=mymod(JJ1+1)
			J2=mymod(JJ1+2)
			do II1=1,3
				I1=mymod(II1+1)
				I2=mymod(II1+2)
				RotM6x6(II1,JJ1)=crot(II1,JJ1)**2
				RotM6x6(II1,JJ1+3)=2.*crot(II1,J1)*crot(II1,J2)
				RotM6x6(II1+3,JJ1)=crot(I1,JJ1)*crot(I2,JJ1)
				RotM6x6(II1+3,JJ1+3)=crot(I1,J1)*crot(I2,J2)+crot(I1,J2)*crot(I2,J1)
			end do
		end do	
		
		!	determine phase
		if (Plap%crystalID(ip) .eq. kHCP) then
			elasmod0 = elasmod1
		else if (Plap%crystalID(ip) .eq. kBCC) then
			elasmod0 = elasmod2
		else 
		    call RunTimeError( FILE_O,'Error: crystalID exceed the existing lattice!')
		end if
		
		!	rotate elastic moduli
		elasModTmp = MATMUL(RotM6x6, elasMod0) 
		Lijkltemp = 0.d0
		Lijkl(:,:,ip) = MATMUL(elasModTmp, TRANSPOSE(RotM6x6) )*1.d6 	! MPa->Pa 
		Lijkltemp = Lijkl(:,:,ip)
		do ii=1,6
			do jj=1,6
				Lijkl(ii,jj,ip) = Lijkltemp(switch(ii),switch(jj))
			end do
		end do
		
		!	echo data for debugging
		if (debug == 1) then
			write(FILE_E,*) 'Lijkl at', ip, 'part is:'
			write(FILE_E,'(6E20.10)') Lijkl(:,:,ip)
		end if
	end do

end subroutine

!*********************************
SUBROUTINE CrystalCloseIOFiles( )
!*********************************
    use FileIO
    implicit none

!---------------------------------------------------------------------72

!	close files
    close(FILE_I)
    close(FILE_E)    
    close(FILE_O)

    return
END
