PROGRAM main

	!
	!	REMARKS:
	!	-	The indices in Mathematica, start with 1. I decided to start with 0 so that all the array entries with index i here
	!		correspond to the entries with i+1 in the mathematica script. Idices from the Mathematica script are denoted by
	!		m:i and the one from the fortran code by f:i
	!	
	!	-	The spline Interpolation of alpha(Energy) is in the range of 
	!		Energy: sbin(m:59 -> m:79) - sdbin(m:59 -> m:79)/2
	!		alpha: 0.2 -> 0.4
	!
	
	use improper_cuadratura
	use mcf_spline
	use param
	use mcf_tipos
	
	implicit none
	
	integer			:: n,m  ! used in loops
	
	real (kind=dp)	:: atau, alphaV,betaV,gammaV,deltaV,rhoV,sigmaV  ! parameters
	
	integer			:: integ_typ,order  ! ni type, ni order
	logical			:: partial  ! ni print steps yes/no ni
	
	integer			:: perr	! used by minuit
	real (kind=dp)	:: arg(10),pval(10),ierr(10),plo(10),phi(10)! used by minuit
	character*10 name(10)
	
	integer 		:: fitbin
	
	call init()

	!============STARTING MAIN CODE FROM HERE=======================================================================================
	
	call MNinit(5,6,7)
	call MNsetI('Fit of ALEPH Data to compute alpha s and the DV Parameters')
	call mnexcm(fcn, 'SET STRategy', (/ 1._dp /), 1, perr, 0) 
	call mnexcm(fcn,'SET PRINTout', (/2.0_dp/),1,perr,0)  ! set print level (2=full)
	
	do fitbin = 58, 78
		
		call MNparm(1,'alphaV',-2.14d0,0.1d0,0.0d0,0.0d0,perr)
		call MNparm(2,'betaV',4.17d0,0.1d0,0.0d0,0.0d0,perr)
		call MNparm(3,'gammaV',0.63d0,0.1d0,0.0d0,0.0d0,perr)
		call MNparm(4,'deltaV',3.48d0,0.1d0,0.0d0,0.0d0,perr)
		call MNparm(5,'rhoV',1.0d0,0.1d0,0.0d0,0.0d0,perr)
		call MNparm(6,'sigmaV',1.0d0,0.1d0,0.0d0,0.0d0,perr)
		call MNparm(7,'atau',0.29664d0,0.01d0,0.0d0,0.0d0,perr)
		
		arg(1) = 5.0_dp
		arg(2) = 6.0_dp
		arg(3) = 3.0_dp
		call mnexcm(fcn,'FIX', (/5.0_dp,6.0_dp/), 2, perr, 0)   
		call mnexcm(fcn,'MIGRAD',    (/10000.0_dp,2._dp/),2,perr,0)  ! perform MIGRAD minimization
		
		do n=1,7
			call mnpout(n,name(n),pval(n),ierr(n),plo(n),phi(n),perr)
		enddo
		
		do n=1,7
		write(*,*) name(n),pval(n),ierr(n)
		enddo
		write(*,*) '------------------'
		write(*,*) '##	fitbin = ',fitbin
		write(*,*) '##	dof = ',79-fitbin-7
		write(*,*) '------------------'
		call mnexcm(fcn, 'SHOw FCNvalue', arg, 0, perr, 0)
			
	end do
	
	!============END OF MAIN CODE===================================================================================================
	
	
	
	CONTAINS !======================================================================================================================
	
	FUNCTION VFO(a, b, c, d, e, f, a0)
	!
	!	produces vector-like array, with entries 
	!	VFO[k]=I_th(x_k)-I_exp(x_k)
	!
		real (kind=dp) :: a,b,c,d,e,f, a0
		real (kind=dp), dimension(fitbin:78) :: VFO
		
		do n = fitbin, 78
			VFO(n) = IntegralPTFO(n, a0) - ni_DVV(n, a, b, c, d, e, f) - FESR(n)
		end do
	
	END FUNCTION VFO
		
	FUNCTION f_DVV(s)
	!
	!	represents the second model of the DV term)	
	!
		real (kind=dp) :: f_DVV
		real (kind=dp), intent(in) :: s
	
		f_DVV = exp(-deltaV)*exp(-gammaV*s**rhoV)*sin(alphaV+betaV*s**sigmaV)
	
	END FUNCTION f_DVV
	
	FUNCTION ni_DVV(bin,a,b,c,d,e,f)
	!
	!	computes the integral of the DV correction
	!
		real (kind=dp) :: a,b,c,d,e,f
		real (kind=dp) :: bb, error, flag
		real (kind=dp) :: ni_DVV, s
		integer, intent(in) :: bin
		integer :: it
		
		bb = huge(1.0d0)
		s = sbin(bin)+dsbin(bin)/2._dp		

		call set_param(a,b,c,d,e,f)
		
		!ni_DVV=exp(-deltaV)*(gammaV*sin(alphaV+betaV*s)+ &
		!betaV*cos(alphaV+betaV*s))/((betaV**2+gammaV**2)*exp(gammaV*s))
!		do while (abs(f_DVV(bb)) == 0.0d0)
!			bb=bb*0.99d0
!		end do
!		call QUANC8(f_DVV,s,bb,1.0d-3,1.0d-3,ni_DVV,error,it,FLAG)
!		if (isnan(f_DVV(s))) then
!			write(*,*) 'NAN'
!			write(*,*) a,b,c,d,e,f
!			stop
!		end if
		!write(*,*) a,b,c,d,e,f
		
		call romberg_improper(f_DVV, sbin(bin)+dsbin(bin)/2._dp, bb, ni_DVV, 1,1.0d-7, .false., 9)
		
	END FUNCTION ni_DVV
	
	FUNCTION Chi2FOV(a, b, c, d, e, f, a0)
		
		real (kind=dp) :: Chi2FOV
		real (kind=dp), dimension(fitbin:78) :: FOV
		real (kind=dp) :: a,b,c,d,e,f,a0
		
		FOV = VFO(a, b, c, d, e, f, a0)
		Chi2FOV = dot_product(FOV(fitbin:78),matmul(SigInvV(fitbin:78,fitbin:78),FOV(fitbin:78)))
		
	END FUNCTION Chi2FOV
	
	SUBROUTINE set_atau(a0)
	
		real (kind=dp) :: a0
		
		atau = a0
		
	END SUBROUTINE set_atau
	
	SUBROUTINE set_param(a, b, c, d, e, f)
	
		real (kind=dp) :: a,b,c,d,e,f
		
		alphaV	=	a
		betaV	=	b
		gammaV	=	c
		deltaV	=	d
		rhoV	=	e
		sigmaV	=	f
	
	END SUBROUTINE set_param
	
	FUNCTION IntegralPTFO(bin,atau) result (iPTFO)
	!
	!	gives the integral of the FOPT expansion
	!	(see arXiv:0806.3156v1, eq. 3.2, p. 7)
	!
		real (kind=dp), intent(in) :: atau
		real (kind=dp) :: iPTFO
		integer, intent(in) :: bin
		iPTFO = &
		(sbin(bin) + dsbin(bin) / 2.0) / (4.0*pi**2) * (1.0 + c(1,1) * (alphaS(bin,atau) / pi) + &
		(c(2,1) + 2.0*c(2,2) * jj(1)) * (alphaS(bin,atau) / pi)**2 + &
		(c(3,1) + 2.0*c(3,2) * jj(1) + 3.0*c(3,3) * jj(2)) * (alphaS(bin,atau) / pi)**3 + &
		(c(4,1) + 2.0*c(4,2) * jj(1) + 3.0*c(4,3) * jj(2) + 4.0*c(4,4) * jj(3)) * (alphaS(bin,atau) / pi)**4 + &
		(c(5,1) + 2.0*c(5,2) * jj(1) + 3.0*c(5,3) * jj(2) + 4.0*c(5,4) * jj(3) + 5.0*c(5,5) * jj(4)) * (alphaS(bin,atau) / pi)**5)
		
	END FUNCTION IntegralPTFO	
		
	REAL(kind=dp) FUNCTION alphaS(bin, atau)
	!
	!	Computes alphas(m) via spline interpolation
	!	spline is computed by values from the mathematica code of santi
	!	
		integer, intent(in) :: bin
		real(kind=dp), intent(in) :: atau
		call spline(Xa,Ya(bin,0:2000),2001,Y2(bin,0:2000))
		call splint(Xa,Ya(bin,0:2000),Y2(bin,0:2000),2001,atau,alphaS)
	
	END FUNCTION alphaS
	
	SUBROUTINE fcn(npar, grad, fval, xval, iflag, futil)
	
		integer :: npar, iflag
		external :: futil
		real (kind=dp) :: grad(0:6), xval(0:6), fval
			
		fval = Chi2FOV(xval(0),xval(1),xval(2),xval(3),xval(4),xval(5),xval(6))
			
	END SUBROUTINE fcn
			
END PROGRAM main
	
