
!=====================================================================72
!
	SUBROUTINE RunTimeError(io, message)

    implicit none
      
    character*(*) message
    integer io
      
    write(io, 1000) message
    write(*,  1000) message
    call CrystalCloseIOFiles( )

    stop

1000  format(/,'***ERROR Message: '/, 3x, a)

    END
!
!=====================================================================72


!=====================================================================72
	subroutine ZTenToVec(ZTen,ZVec,n)

	use numtype    
	implicit none

	integer n,i
	real(kind=8) ZTen(n,n), ZVec(NVEC)

	do i=1,3
		ZVec(i)=ZTen(i,i)
	enddo

	ZVec(4)=ZTen(1,2)+ZTen(2,1)
	ZVec(5)=ZTen(1,3)+ZTen(3,1)
	ZVec(6)=ZTen(2,3)+ZTen(3,2)    
		   
	return     
	end
!=====================================================================72


!=====================================================================72
	subroutine CheckMatriceZero(Matrice,n)

	USE FILEIO
	implicit none

	integer n,i,j
	real(kind=8) Matrice(n,n)
	real(kind=8) tolerance

	tolerance=tiny(1.0E0)
	do i=1,n
		do j=1,n
		if (abs(Matrice(i,j)) .le. tolerance) then
			Matrice(i,j) = 0.0
		endif
		enddo
	enddo

	return
	end
!
!=====================================================================72


!====================================================================72
!
      SUBROUTINE OuterProductVec(vecU, vecV, outer, n )
      use Numtype
      implicit none


      integer n
      real*8  vecU(n), vecV(n)
      real*8  outer(n, n)
     
      integer i, j

      do i = 1, n
         do j = 1, n
            outer(i,j) = vecU(i) * vecV(j)
         enddo
      enddo

      return
      END
!
!=====================================================================72


!=====================================================================72
!
      SUBROUTINE AnglesToRotMatrix(angle, crot, n)

      implicit none

      integer n
      real*8  angle(n)
      real*8  crot(n,n)

      real*8  sps, cps, sth, cth, sph, cph

   
!------- Construct [C] matrix from euler angles (in radians)
!-------     {a}_sm = [C] {a}_xtal     
      
      sps = dsin(angle(1))
      cps = dcos(angle(1))
      sth = dsin(angle(2))
      cth = dcos(angle(2))
      sph = dsin(angle(3))
      cph = dcos(angle(3))
!--- Marin's code      
      !crot(1,1) = -sps * sph - cps * cph * cth
      !crot(2,1) =  cps * sph - sps * cph * cth
      !crot(3,1) =  cph * sth
      !crot(1,2) =  cph * sps - sph * cps * cth
      !crot(2,2) = -cps * cph - sps * sph * cth
      !crot(3,2) =  sph * sth
      !crot(1,3) =  cps * sth
      !crot(2,3) =  sps * sth
      !crot(3,3) =  cth

!---Mine
      !crot(1,1) = -sph*sps + cps * cph * cth
      !crot(1,2) =  -sps * cph - cps * sph * cth
      !crot(1,3) =  -cps * sth
      !crot(2,1) =  cps * sph + sps * cth * cph
      !crot(2,2) = cps * cph - sps * sph * cth
      !crot(2,3) =  -sps * sth
      !crot(3,1) =  cph * sth
      !crot(3,2) =  -sph*sth
      !crot(3,3) =  cth 
      
!---Bunge
      crot(1,1) = cph*cps - sps * sph * cth
      crot(1,2) =  -cps * sph - sps * cph * cth
      crot(1,3) =  sps * sth
      crot(2,1) =  sps * cph + cps * sph * cth
      crot(2,2) = -sps * sph + cps * cph * cth
      crot(2,3) =  -cps * sth
      crot(3,1) =  sph * sth
      crot(3,2) =  cph*sth
      crot(3,3) =  cth
      
      return
      END
!
!=====================================================================72!


!=====================================================================72
!
SUBROUTINE RotateSlipGeometry( SlipG0 )

    USE NumType
    use DataType
    USE FILEIO
    USE SlipGeo
    USE OriPar
    use PlaPar

	IMPLICIT NONE
!	
    type (SlipSys) SlipG0
	INTEGER IPart,  ISlip, Numslip

	REAL(KIND=8) :: localx(3), localy(3),globalx(3),globaly(3),Tempv(3)
!---------------------------------------------------------------------------------
! It is assumed that A, M and P are already caculated in the global Coordinates.	
! So only need to rotate the slip syste from local to global directions

    if (debug==1) write(FILE_E, *) '---------Begin RotateSlipGeometry----------'
	!localx=(/ -1.d0, 0.d0,1.d0/)
	!localy=(/0.d0,1.d0,0.d0/)
	!tempv=localx/norm2(localx)
	!localx=tempv
	!tempv=localy/norm2(localy)
	!localy=tempv

	DO IPart=1,NPart

  	  if (Plap%crystalID(iPart) .eq. kHCP) then
		numslip=PhSlip(1)
	  else if (Plap%crystalID(iPart) .eq. kBCC) then
		numslip=PhSlip(2)
	  else 
		call RunTimeError( FILE_O,                                         &
                       'Error: crystalID exceed the existing lattice!')
	  end if

		DO ISlip = 1, NumSlip
		    write(*,*) IPART, ISLip
		    write(*,*) SlipG%VecS0(:,ISLIP, IPart)
		    write(*,*) SlipG%VecM0(:,ISLIP, IPart)
		    write(*,*) OriP%gcrot0(:,:,IPart)
		    localx=MATMUL(OriP%gcrot0(:,:,IPart),SlipG%VecM0(:,ISlip, IPart))
		    write(*,*) 'localx'
		    write(*,*) localx
			SlipG0%VecM0(:,ISlip, IPart)=MATMUL(OriP%gcrot0(:,:,IPart),SlipG%VecM0(:,ISlip, IPart))
			SlipG0%VecS0(:,ISlip, IPart)=MATMUL(OriP%gcrot0(:,:,IPart),SlipG%VecS0(:,ISlip, IPart))
!
            write(*,*) SlipG0%VecS0(:,ISLIP, IPart)
            write(*,*) SlipG0%VecM0(:,ISLIP, IPart)
            write(*,*) 'flag1'
            CALL OuterProductVec(SlipG0%VecS0(:,ISlip, IPart), SlipG0%VecM0(:,ISlip, IPart), &
                                                    SlipG0%ZTen0(:,:,ISlip,IPart), 3)
            write(*,*) 'flag2'                                                    
            CALL ZTenToVec(SlipG0%ZTen0(:,:,ISlip,IPart), SlipG0%ZVec0(:, ISlip,IPart),3)
		ENDDO

    if (debug==1) then
	    WRITE(FILE_E,*)  'Rotation.f90,crot0 ', IPart, '-th partition is:'
	    WRITE(FILE_E,'(3e16.8)')  OriP%gcrot0(:,:,IPart)
    endif

!--This part is used to check the correctness of the rotatension tensor as well as the two definitions of rotation tensor from two vectors or from the euler angles. 


	!globalx(:)=MATMUL(crot(:,:,IPart), localx(:) )
	!globaly(:)=MATMUL(crot(:,:,IPart), localy(:) )	

	!WRITE(FILE_E,*)  'localx=================================globalx'
	!WRITE(FILE_E,'(3(e16.8) 5x,3(e16.8))' ) localx,globalx
	!WRITE(FILE_E,*)  'localy=================================globaly'
	!WRITE(FILE_E,'(3(e16.8) 5x,3(e16.8))' ) localy,globaly	
	
	ENDDO

RETURN
END
!=====================================================================72
!
!=====================================================================72
!
      real*8 FUNCTION InnerProductVec(vecU, vecV, n)
      use NumType
      implicit none

      integer n
      real*8  vecU(n), vecV(n)

      integer i
      
!---- compute inner product of two vectors: product = {u}^T*{v}
!
      InnerProductVec = pzero
      do i = 1, n
         InnerproductVec = InnerProductVec + vecU(i) * vecV(i)
      enddo
      InnerproductVec = InnerProductVec

      return
      END
!
!=====================================================================72


!=====================================================================72
!
      real*8 FUNCTION Power(x,y)
      use NumType
      implicit none

      real*8  x, y

!---- evaluates  x^y

      if (x .eq. 0.0) then
         if (y .gt. 0.0) then
            Power = 0.d0
         elseif (y .lt. 0.0) then
            Power = 1.d+300
         else
            Power = 1.d0
         endif
      else
         Power = y * log10(dabs(x))
         if (Power .gt. 300.0) then
            Power = 1.d+300
         else
            Power = 10.d0 ** Power
         endif
         if (x .lt. 0.0) Power = -Power
      endif

      return
      END
!
!=====================================================================72


!=====================================================================72
!
      SUBROUTINE MultAxu(A, u, v, n)
        
      use NumType
      implicit none
      

      integer n
      real*8  A(n,n), u(n), v(n)

      integer i, j

      v=pzero

      do i = 1, n
         do j = 1, n
            v(i) = v(i) + A(i,j) * u(j)
         end do
      end do

      return
      END

!=====================================================================72


!=====================================================================72
!
      SUBROUTINE MultAxB(A, B, C, n)
      
      use NumType
      implicit none

     
      integer n
      real*8  A(n, n), B(n, n), C(n, n)

      integer i, j, k
      
      call SetTensor(C, pzero, n*n)

      do i = 1, n 
         do j = 1, n
            do k = 1, n
               C(i,j) = C(i,j) + A(i,k) * B(k,j)
            enddo
         enddo
      enddo
   
      return
      END
!
!=====================================================================72


!=====================================================================72
!
      SUBROUTINE MultAxB_G(A, B, C, m1, m2, m3)
      use NumType
      implicit none
      
     
      integer m1, m2, m3
      real*8  A(m1, m2), B(m2, m3), C(m1, m3)

      integer i, j, k

!---- Dimensions: A is m1 x m2; B is m2 x m3; C is m1 x m3 
!
      call SetTensor(C, pzero, m1*m3)

      do i = 1, m1 
         do j = 1, m3
            do k = 1, m2
               C(i,j) = C(i,j) + A(i,k) * B(k,j)
            enddo
         enddo
      enddo
   
      return
      END
!
!=====================================================================72


!=====================================================================72
!
      SUBROUTINE SetTensor(tensor, value, n)
      use NumType
      implicit none

      integer n
      real*8  value
      real*8  tensor(n)

      integer i

      do i = 1, n
         tensor(i) = value
      enddo

      return
      END
!
!=====================================================================72


!=====================================================================72
!
      SUBROUTINE SetToScaledTensor(fac_A, tensA, tensB, n)
      use NumType
      implicit none

      integer n
      real*8  fac_A
      real*8  tensA(n), tensB(n)

      integer i

      do i = 1, n
          tensB(i) = fac_A * tensA(i)
      enddo

      return
      END
!
!=====================================================================72


!=====================================================================72
!
      SUBROUTINE Vec6x1ToMat3x3Symm(vector, matrix, n)

      use NumType
      
      implicit none

      integer n
      real*8  matrix(n,n), vector(6)

!
      if (n .ne. 3)    & 
       write(*,*) 'Vect6x1ToMat3x3Symm: n =! 3'

      matrix(1,1) = vector(1)
      matrix(2,2) = vector(2)
      matrix(3,3) = vector(3)
      matrix(1,2) = vector(4)
      matrix(2,1) =matrix(1,2) 
      matrix(1,3) = vector(5)
      matrix(3,1) = matrix(1,3)
      matrix(2,3) = vector(6)
      matrix(3,2) =matrix(2,3) 

      return
      END
!
!=====================================================================72


!=====================================================================72
!
      SUBROUTINE WriteMessage(io, message)

      implicit none
      
      character*(*) message
      integer io
      
      write(io, 1000) message

1000  format('***Message: ', a)

      return
      END
!
!=====================================================================72


!=====================================================================72
!
      real*8 FUNCTION SignOf(value)

      implicit none

      real*8  value

!------- compute the sign of a number

      SignOf = 1.0
      if (value .lt. 0.0) SignOf = -1.0

      return
      END
!
!=====================================================================72

!--> zx
	integer function mymod(II)
		implicit none
		integer  II
	
		if (II .gt. 6) then
			write(*,*) "Input of mymod should not be greater than 6 in current case!"
			stop
		endif
	
		if (II .le. 3) then
			mymod=II
		else
			mymod=II-3
		endif
	end function mymod	
!--> zx
