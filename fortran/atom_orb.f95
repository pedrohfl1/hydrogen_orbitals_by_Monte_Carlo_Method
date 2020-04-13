! file contains code which creates the hydrogen's orbitals. 
! 
! Dependencies: 
!
! 	A SHTOOLS 4.4 tool and its dependencies are required to 
!	calculate Legendre's associated polynomials. See the link 
!	below for installation and more information.
!
!	https://shtools.oca.eu/shtools/
!
! Licenses: 
!
!	This file is subject to the terms and conditions defined in
!	file 'LICENSE', which is part of this source code.
!
	program atomOrb
		implicit none
		integer, parameter :: n = 4, l = 1, m = 1
		real, parameter    :: dx = 0.0 , dy = 60.0, dz = 60.0
		real, parameter    :: dr = 0.2, dtheta = 0.2, icr = 1.05
		integer, parameter :: points = 50000

		integer            :: rDiv, thetaDiv, i
		real :: theta, r, w, max, pi
		real :: angDF, radDF
		real :: wave2Max, psi2
		real :: rAbs, thetaVal
		real :: rVector(3)
	
		pi = 4.0 * atan(1.0)

		call init_random_seed()

		rDiv = int(rAbs((/ dx, dy, dz /)) / dr) + 1
		thetaDiv = int(2.0 * pi / dtheta) + 1
		max = wave2max((/dx,dy,dz/), (/n,l,m/), rDiv, thetaDiv)

		do while (i <= points)
		CALL RANDOM_NUMBER(w)
			w = w * max * icr
			call fillVector(rVector, dx, dy, dz)
			r = rAbs(rVector)
			theta = thetaVal(rVector)
			psi2 = radDF(n, l, r) * angDF(l, m, theta) 
			if(w < psi2) then
				write(*,*) rVector,  w
				i = i + 1
			end if
		end do
	end program atomOrb

	function wave2max(coordRange, nlm, rDiv, thetaDiv)
		implicit none
		interface
		function radDF(n, l, r)
			implicit none
			! dummy arguments        
			real :: radDF
			! local variables
			integer :: n, l, Z = 1
			integer :: k
			real :: r, r0, a0 = 0.52917721067
			real*8 :: normTermSquared, radialTerm
			real*8 :: arr(6), sum
		end function radDF 
		end interface
	
		interface
		function angDF(l, m, theta) 
		use SHTOOLS
		! function result     
		implicit none      
		! dummy arguments        
		real :: angDF 
		
		! local variables 
		integer 			  :: l, m, arrayindex
		real*8                            :: z, normRoot
		real*8, dimension ((l+1)*(l+2)/2) :: p
		real 				  :: theta, pi
		end function angDF
		end interface
	
		real :: wave2Max

		real :: coordRange(3)
		integer :: nlm(3), rDiv, thetaDiv
		real :: rmax, sum, dr, dtheta, pi, curr(2)
		real :: maxValue, r, theta, calc
		integer :: i, j
		
		pi = 4.0 * atan(1.0)
		
		do i = 1, 3, 1
			sum = sum + coordRange(i) ** 2.0
		end do 
	
		rmax = sqrt(sum)
	
		dr = rMax / (1.0 * rDiv)
		dTheta = pi / (1.0 * thetaDiv)
		maxValue = 0.0

		do i = 0, rDiv, 1
			r = i * dr
			do j = 0, thetaDiv, 1
				theta = j * dtheta
				curr(1) = angDF(nlm(2), nlm(3), theta) 
				curr(2) = radDF(nlm(1), nlm(2), r)
				calc = curr(1) * curr(2)

				if(calc > maxValue) then
					maxValue = calc
				end if
			end do
		end do

		wave2max = maxValue
	end function
	
	subroutine fillVector(r, dx, dy, dz)
	real, intent(out) :: r(3)
	real :: dx, dy, dz

	CALL RANDOM_NUMBER(r(1))
	r(1) = dx * (r(1) - 0.5) 

	CALL RANDOM_NUMBER(r(2))
	r(2) = dy * (r(2) - 0.5) 

	CALL RANDOM_NUMBER(r(3))
	r(3) = dz * (r(3) - 0.5) 
	
	end subroutine
	
	function rAbs(r)
		implicit none
		real :: rAbs
		real :: r(3)
		
		rAbs = sqrt(r(1)**2+r(2)**2+r(3)**2)

	end function

	function thetaVal(r)
		implicit none
		
		real :: thetaVal
		real :: r(3), sum
		integer :: i
				
		do i = 1, 3, 1
			sum = sum + r(i) ** 2.0
		end do
	
		thetaVal = acos(r(3) / sqrt(sum))
	
	end function thetaVal

	function radDF(n, l, r)
		implicit none
		INTERFACE 
		   FUNCTION Factorial(n)
		     INTEGER*8 :: Factorial
		     INTEGER, INTENT(IN) :: n
		   END FUNCTION Factorial
		END INTERFACE	
	
		! dummy arguments        
		real :: radDF
		! local variables
		integer :: n, l, Z = 1
		integer :: k
		real :: r, r0, a0 = 0.52917721067
		real*8 :: normTermSquared, radialTerm, arr(6), sum
		
		r0 = 2.0 * Z/(n * a0) * r
		sum = 0.0

		
		do k = 0, n - l - 1, 1
			arr(1) = (-1.0) ** (k + 1) 
			arr(2) = factorial(n + l) ** 2.0
			arr(3) = factorial(n - l - 1 - k)
			arr(4) = factorial(2 * l + 1 + k) 
			arr(5) = factorial(k)
			arr(6) = r0 ** k
			sum = sum+arr(1)*arr(2)/(arr(3)*arr(4)*arr(5))*arr(6) 
		end do
		
		arr(:) = 0.0
		arr(1) = (2.0 * Z / (n * a0)) ** 3.0
		arr(2) = factorial(n - l - 1)
		arr(3) = 2.0 * n
		arr(4) = factorial(n + l) ** 3.0
		
		normTermSquared = arr(1)*arr(2)/(arr(3)*arr(4))

		radialTerm = exp(-r0 / 2.0) * (r0 ** l);

		radDF = (sum * radialTerm) ** 2.0
		radDF = radDF * normTermSquared
	end function radDF

	function angDF (l, m, theta) 
	use SHTOOLS
	! function result     
	implicit none      
	INTERFACE 
	   FUNCTION Factorial(n)
	     INTEGER*8 :: Factorial
	     INTEGER, INTENT(IN) :: n
	   END FUNCTION Factorial
	END INTERFACE	
	
	! dummy arguments        
	real :: angDF 
	
	! local variables 
	integer 				:: l, m, arrayindex
	real*8                                  :: z, normRoot
	real*8, dimension ((l + 1) * (l + 2)/2) :: p
	real 					:: theta, pi
		
	z = cos(theta)
	call PLegendreA (p, l, z)
	arrayindex = PlmIndex (l, m)	

		
	pi = 4.0 * atan(1.0)
	normRoot = 1.0 * factorial(l - m)/ factorial(l + m)
	normRoot = normRoot * (2.0 * l + 1.0)/(4.0 * pi)

	angDF = normRoot * p(arrayIndex) ** 2.0
	
	end function angDF 

! Evaluate the factorial
! This was adapted from a post from user gsal at jun 2, 2015 to a 
! Physics Forums, and visited in March 2019. See link bellow.

! Reference https://www.physicsforums.com/threads/factorial-fortran-90.816501/

	RECURSIVE FUNCTION Factorial(n)  RESULT(Fact)
	
	IMPLICIT NONE
	INTEGER*8 :: Fact
	INTEGER, INTENT(IN) :: n
	
	IF (n == 0) THEN
	   Fact = 1
	ELSE
	   Fact = n * Factorial(n-1)
	END IF
	
	END FUNCTION Factorial

! Initialize a pseudo-random number sequence
! This subroutine is a partial copy found on the link below 
! and visited in March 2019.

! https://gcc.gnu.org/onlinedocs/gcc-4.2.3/gfortran/RANDOM_005fSEED.html

! Description:
!    Restarts or queries the state of the pseudorandom number 
!    generator used by RANDOM_NUMBER.

	subroutine init_random_seed()
	
	integer :: i, n, clock
	integer, dimension(:), allocatable :: seed
	
	call RANDOM_SEED(size = n)
	allocate(seed(n))
	
	call SYSTEM_CLOCK(COUNT=clock)
	
	seed = clock + 37 * (/ (i - 1, i = 1, n) /)

	call RANDOM_SEED(PUT = seed)
	deallocate (seed)
	
	end subroutine
