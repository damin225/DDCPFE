!=====================================================================72
!
      real*8 FUNCTION SSKineticEqn(rss1, crss1, kflag, argmin, argmax, &
                                   is, ip)
      
      use NumType
      use PlaPar
      use FILEIO

      implicit none
	

      integer kflag, is, ip
      real*8  crss1, rss1, pow,arg, argmin, argmax



      real*8  CheckArgPowLaw, SignOf, Power,tem
	  
	  tem=295.0d0

!------- check argument of power law

    arg = rss1/crss1
    if (dabs(arg) .ne. 0) arg = CheckArgPowLaw(arg, SignOf(arg), argmin, argmax)

      pow = Power(dabs(arg), 1./PlaP%xm-1.)

      if (kflag .eq. kGAMDOT) then 
         SSKineticEqn = PlaP%gam0 * arg * pow
      else if (kflag .eq. kdGAMdTAU) then
         SSKineticEqn = PlaP%gam0 / (PlaP%xm*crss1) * pow
      else if (kflag .eq. kdGAMdKAPP) then
         SSKineticEqn = -PlaP%gam0 / (PlaP%xm*crss1) * arg * pow
!---- This is used to caculate the dGamDOTdTau for caculating Jacobian when sloving the first set of equations
      else if (kflag .eq. kDDGAM) then

	if (Plap%crystalID(ip) .eq. kHCP) then

		if (is.le.3) then
	 SSKineticEqn =0.5*Plap%ddavg(1,1)*Plap%vid(1,1)*Plap%burgers(1,1)**2.*exp(((dabs(rss1)-crss1)-Plap%detF(1,1)/Plap%detV(1,1))*Plap%detV(1,1)/Plap%kbolt/tem)*SignOf(rss1)
		else if (is.le.6) then
	 SSKineticEqn =0.5*Plap%ddavg(1,2)*Plap%vid(1,2)*Plap%burgers(1,2)**2.*exp(((dabs(rss1)-crss1)-Plap%detF(1,2)/Plap%detV(1,2))*Plap%detV(1,2)/Plap%kbolt/tem)*SignOf(rss1)
		else if (is.le.12) then
	 SSKineticEqn =0.5*Plap%ddavg(1,3)*Plap%vid(1,3)*Plap%burgers(1,3)**2.*exp(((dabs(rss1)-crss1)-Plap%detF(1,3)/Plap%detV(1,3))*Plap%detV(1,3)/Plap%kbolt/tem)*SignOf(rss1)
		else if (is.le.30) then
	 SSKineticEqn =0.5*Plap%ddavg(1,4)*Plap%vid(1,4)*Plap%burgers(1,4)**2.*exp(((dabs(rss1)-crss1)-Plap%detF(1,4)/Plap%detV(1,4))*Plap%detV(1,4)/Plap%kbolt/tem)*SignOf(rss1)
		else
		     call RunTimeError(File_o, 'SSKineticEqn: is exceed maxinum number')
		endif

	 else  if (Plap%crystalID(ip) .eq. kBCC) then
	 SSKineticEqn =0.5*Plap%ddavg(2,1)*Plap%vid(2,1)*Plap%burgers(2,1)**2.*exp(((dabs(rss1)-crss1)-Plap%detF(2,1)/Plap%detV(2,1))*Plap%detV(2,1)/Plap%kbolt/tem)*SignOf(rss1)

	 else 
		call RunTimeError( FILE_O,'Error: crystalID exceed the existing lattice!')
	 end if
     else
         call RunTimeError(File_o, 'SSKineticEqn: unknown kflag!')
     endif

      return
      END
!
!=====================================================================72


!=====================================================================72
!
      logical FUNCTION NConverged(res, toler, n)
      
      use NumType

      implicit none


      real*8  toler
      real*8  res(NVEC)

      integer i,n
 
!------- check convergence on residual

      NConverged = .true.
      do i = 1, n
         NConverged = ( (dabs(res(i)) .lt. toler) .and. NConverged)
      enddo

      return
      END
!     
!=====================================================================72


!=====================================================================72
!
      SUBROUTINE BoundForArgPowLaw(xm, argmin, argmax)

      use NumType      
      implicit none

      real*8  xm, xmm, argmin, argmax

      real*8  EXPON
      data EXPON /280.d0/

!---- limits (bounds) on the value of the argument for power law
!---- note: In general: xm <= 1.0 (1/xm >= 1.0)

      xmm = pone/xm - pone
      if (dabs(xm-pone) .lt. TINY(1.d0)) xmm = pone

      argMin = dexp(-EXPON/(xmm)*dlog(10.d0))
      argMax = dexp( EXPON/(xmm)*dlog(10.d0))

      return
      END
!
!=====================================================================72


!=====================================================================72
!

    real*8 FUNCTION CheckArgPowLaw(arg, sgn, argmin, argmax)
      
      use NumType
 
      implicit none
      
      real*8  arg, sgn, argmin, argmax



!---- check range of argument for power law

      if (dabs(arg) .lt. argMin) then
         CheckArgPowLaw = argMin * sgn
      else if (dabs(arg) .ge. argMax) then
         CheckArgPowLaw = argMax * sgn
      else
         CheckArgPowLaw = arg
      endif

      return
      END
!=====================================================================72

