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

	integer			:: n,m, state, statenull, nbin, NPARI, NPARX, istat  ! used in loops

	real (kind=dp)	:: alpha,beta,gamma,delta,rho,sigma  ! parameters

	real (kind=dp)	:: p, q, df, bound,fval, fedm, ERRDEF, twopi=2.0d0*pi

	integer			:: integ_typ,order, status  ! ni type, ni order
	logical			:: partial  ! ni print steps yes/no ni

	integer			:: perr	! used by minuit
	real (kind=dp)	:: arg(13),pval(13),ierr(13),plo(13),phi(13)! used by minuit
	character*20 name(13)

	!============STARTING MAIN CODE FROM HERE=======================================================================================

	call MNinit(5,6,7)
	call MNsetI('Fit of ALEPH Data to compute alpha s and the DV Parameters')

	npar = 7

	npar = npar + 1

	allocate(par(npar), stat=state)

	call readConfig(par, npar)

	nbin=par(npar)

	call init(nbin)

	pval(1:npar) = par(1:npar)

		call MNparm(1,'delta',pval(1),1.0d-4,0.0d0,0.0d0,perr)
		call MNparm(2,'gamma',pval(2),1.0d-4,0.0d0,0.0d0,perr)
		call MNparm(3,'alpha',pval(3),1.0d-4,0.0d0,0.0d0,perr)
		call MNparm(4,'beta',pval(4),1.0d-4,0.0d0,0.0d0,perr)
		call MNparm(5,'rho',pval(5),1.0d-4,0.0d0,0.0d0,perr)
		call MNparm(6,'sigma',pval(6),1.0d-4,0.0d0,0.0d0,perr)
		call MNparm(7,'atau',pval(7),1.0d-4,0.0d0,0.0d0,perr)

	arg(1) = 5.0_dp
	arg(2) = 6.0_dp
	arg(3) = 3.0_dp
	!call mnexcm(fcn, 'Set LIMits', (/3._dp, -2*twopi, twopi /), 3, perr, 0)
	!call mnexcm(fcn, 'Set LIMits', (/ 4._dp, -twopi, 2*twopi /), 3, perr, 0)
	!call mnexcm(fcn, 'Set LIMits', (/2._dp, 0.0d0, 4.0d0 /), 3, perr, 0)
	!call mnexcm(fcn, 'Set LIMits', (/ 2._dp, 0.0d0, 10d0 /), 3, perr, 0)
	call mnexcm(fcn,'FIX', (/3.0d0,5.0_dp,6._dp,11.0_dp,12._dp,13._dp/), 3, perr, 0)
	!call mnexcm(fcn, 'HESse', (/ 2._dp /), 0, perr, 0)
	!call mnexcm(fcn, 'SET EPS', (/ 1.0d-7 /), 1, perr, 0)
	call mnexcm(fcn, 'SET STRategy', (/ 2._dp /), 1, perr, 0)
!	call mnexcm(fcn, 'HESse', (/ 2._dp /), 0, perr, 0)
!	call mnexcm(fcn,'SIMplex', (/100000.0_dp,0.1_dp/),2,perr,0)
!	call mnexcm(fcn, 'SCAn', arg, 0, perr, 0)
	!call mnexcm(fcn,'RELeas', (/5.0_dp,6.0_dp/), 2, perr, 0)
	!call mnexcm(fcn, 'SET ERRordef', (/4._dp/),1,perr,0)
	!call mnexcm(fcn, 'HESse', (/ 2._dp /), 0, perr, 0)
	call mnexcm(fcn,'MINimize',    (/10000.0_dp,1.0d-5/),2,perr,0)  ! perform MIGRAD minimization
	call mnexcm(fcn,'RELeas', (/5.0_dp,6.0_dp/), 2, perr, 0)
!	do n=1,7
!		call mnpout(n,name(n),pval(n),ierr(n),plo(n),phi(n),perr)
!	enddo
	call mnexcm(fcn,'MINimize',    (/10000.0_dp,1.0d-5/),2,perr,0)
!	call MNparm(5,'gamma1',4._dp*pval(3),ierr(3),0.0d0,0.0d0,perr)
!	call MNparm(6,'delta1',4._dp*pval(4),ierr(4),0.0d0,0.0d0,perr)
	!call mnexcm(fcn, 'SET STRategy', (/ 2._dp /), 1, perr, 0)
	call mnexcm(fcn,'RELeas', (/3.0_dp,6.0_dp/), 1, perr, 0)
	!call mnexcm(fcn, 'HESse', (/ 2._dp /), 0, perr, 0)
	!call mnexcm(fcn,'SIMplex',    (/10000.0_dp,100.0_dp/),2,perr,0)  ! perform MIGRAD minimization
	call mnexcm(fcn,'MINimize',    (/10000.0_dp,1.0d-5/),2,perr,0)  ! perform MIGRAD minimization
	!call mnexcm(fcn,'IMProve',    (/10000.0_dp/),1,perr,0)  ! perform MIGRAD minimization
	!call mnexcm(fcn,'SEEk',    (/10000.0_dp/),1,perr,0)  ! perform MIGRAD minimization
	!call mnexcm(fcn, 'HESse', (/ 2._dp /), 0, perr, 0)
	!call mnexcm(fcn,'MINO',    (/10000.0_dp,1.0d-5/),2,perr,0)
	!call mnexcm(fcn,'Contour',(/5._dp,7._dp,25._dp/),3,perr,0)
	do n=1,7
		call mnpout(n,name(n),pval(n),ierr(n),plo(n),phi(n),perr)
	enddo

!	call mnexcm(fcn, 'SCAn', (/ 0._dp, 100.0d0,0.0d0,0.0d0  /), 0, perr, 0)

!!	call mnexcm(fcn, 'SCAn', (/ 5._dp, 40.0d0,pval(5) - 0.1d0*pval(5),pval(5) + 0.1d0*pval(5)  /), 4, perr, 0)
!!	call mnexcm(fcn, 'SCAn', (/ 6._dp, 40.0d0,pval(6) - 0.1d0*pval(6),pval(6) + 0.1d0*pval(6)  /), 4, perr, 0)
!	!call mnexcm(fcn, 'SCAn', (/ 6._dp, 40.0d0,pval(6) - pi,pval(6) + pi  /), 4, perr, 0)
!	!call mnexcm(fcn, 'SHOw EIGenvalues', (/ 2._dp /), 0, perr, 0)
!	call mnexcm(fcn, 'MIGrad', (/ 10000.0_dp, 0.1_dp /), 2, perr, 0)
!	call mnexcm(fcn,'SIMplex', (/100000.0_dp,0.1_dp/),2,perr,0)
!	!call mnexcm(fcn,'IMProve',    (/10000.0_dp/),1,perr,0)
!	call mnexcm(fcn, 'HESse', (/ 2._dp /), 0, perr, 0)

	do n=1,7
	write(*,*) name(n),pval(n),ierr(n),kind(pval(n))
	enddo

	call MNSTAT(fval, fedm, ERRDEF, NPARI, NPARX, istat)

	call cdfchi ( 1, p, q, fval, 1._dp*(79-binmin-7), status, bound )
	write(*,*) '------------------'
	write(*,*) '##	fitbin/Energy = ',binmin,' / ',sbin(binmin)+0.5_dp*dsbin(binmin)
	write(*,*) '##	dof = ',79-binmin-7
	write(*,*) '------------------'
	call mnexcm(fcn, 'SHOw FCNvalue', arg, 0, perr, 0)
	write(*,*) '------------------'
	write(*,*) 'P-Value = ',100d0*q

	!write(95,*) 79-binmin-7,binmin,sbin(binmin)+0.5_dp*dsbin(binmin),pval(1),ierr(1),pval(2), &
	!ierr(2),pval(3),ierr(3),pval(4),ierr(4),pval(5),ierr(5),pval(6),ierr(6),pval(7),ierr(7), &
	!fval,100d0*q,istat

!	call contour(1,7,100)

	!call PltData(100,0.0d0,sbin(78))

!	write(*,*)	0.5_dp*ni_DVV(binmin, pval(1),pval(2),pval(3),pval(4),pval(5)) + &
!				0.5_dp*ni_DVV(n, pval(6),pval(7),pval(8),pval(9),pval(10))

!	write(*,*)	ni_DVV(binmin, 3.48d0, 0.63d0, -2.14d0, 4.17d0, 1.0d0)

	!write(*,*) ni_DVV(binmin,pval(1),pval(2),pval(3),pval(4),pval(5),pval(6),pval(7),pval(8),pval(9),pval(10))

!	write(*,*) SigInvV(58,58)

!write (*,*) H(4._dp)

!write(*,*) truecovv(58,58)

	call set_param(pval(1),pval(2),pval(3),pval(4),pval(5),pval(6))
	write(*,*) 'CORRECTION TO FM at 43 = ',&
	(exp(-alpha - beta * (sbin(binmin)+0.5d0*dsbin(binmin)) ** rho) * &
	sin(gamma + delta * (sbin(binmin)+0.5d0*dsbin(binmin)) ** sigma))/&
	(exp(-alpha - beta * (sbin(binmin)+0.5d0*dsbin(binmin))) * &
	sin(gamma + delta * (sbin(binmin)+0.5d0*dsbin(binmin))))-1.0d0
	write(*,*) 'CORRECTION TO FM at 58 = ',&
	(exp(-alpha - beta * (sbin(58)+0.5d0*dsbin(58)) ** rho) * &
	sin(gamma + delta * (sbin(58)+0.5d0*dsbin(58)) ** sigma))/&
	(exp(-alpha - beta * (sbin(58)+0.5d0*dsbin(58))) * &
	sin(gamma + delta * (sbin(58)+0.5d0*dsbin(58))))-1.0d0

	!============END OF MAIN CODE===================================================================================================



	CONTAINS !======================================================================================================================

	FUNCTION VFO(a0,b0,c0,d0,c1,d1, atau)
	!
	!	produces vector-like array, with entries
	!	VFO[k]=I_th(x_k)-I_exp(x_k)
	!
		real (kind=dp) :: a0,b0,c0,d0,c1,d1, atau
		real (kind=dp), dimension(binmin:78) :: VFO

		do n = binmin, 78
			VFO(n) = IntegralPTFO(n, atau) - ni_DVV(n, a0,b0,c0,d0,c1,d1) - FESR(n)
		end do

	END FUNCTION VFO

	FUNCTION f_DVV(s)
	!
	!	represents the second model of the DV term)
	!
		real (kind=dp) :: f_DVV
		real (kind=dp), intent(in) :: s

!		f_DVV = deltaV * ( (cos(alphaV * s**rhoV + betaV) - exp(-gammaV * s**rhoV) ) &
!		/ (cosh(gammaV * s**rhoV) - cos(alphaV * s**rhoV + betaV)) )

		f_DVV =	exp(-alpha - beta * s ** rho) * sin(gamma + delta * s ** sigma)

	END FUNCTION f_DVV

	FUNCTION testi(x)
		real (kind=dp)	::	testi
		real, intent(in)	::	x

		testi = 1._dp / (1._dp + 10._dp * x)**2
	END FUNCTION testi

	FUNCTION ni_DVV(bin,a,b,c,d,c1,d1)
	!
	!	computes the integral of the DV correction
	!
		real (kind=dp) :: a,b,c,d,c1,d1
		real (kind=dp) :: bb
		integer     :: fier, fneval, flast, flenw
		real (kind=dp) :: ni_DVV, s, fabserr, ss
		real (kind=dp) :: fwork(400)
		integer, dimension(200) :: fiwork
		integer, intent(in) :: bin
		!external f_DVV

		bb = huge(1.0d0)
		s = sbin(bin)+dsbin(bin)/2._dp

		call set_param(a,b,c,d,c1,d1)

		!ni_DVV=exp(-deltaV)*(gammaV*sin(alphaV+betaV*s)+ &
		!betaV*cos(alphaV+betaV*s))/((betaV**2+gammaV**2)*exp(gammaV*s))
		!ni_dvv = (d1*Cos(c1 + d1*s) + 2*b*Sin(c1 + d1*s))/((4*b**2 + d1**2)*Exp(2*(a + b*s))) + &
		!exp(-a-b*s)*(d*cos(c+d*s)+b*sin(c+d*s))/(b**2+d**2)
		!call romberg_improper(f_DVV, sbin(bin)+dsbin(bin)/2._dp, bb, ss, 1,1.0d-10, .false.,9)
		call dqagi(f_DVV,sbin(bin)+dsbin(bin)/2._dp,1,1.0d-18,1.0d-40,ni_DVV,fabserr,fneval, &
		fier,100,flenw,flast,fiwork,fwork)
		!write(*,*) fier, ni_DVV
	END FUNCTION ni_DVV

	FUNCTION ni_DV(s,a,b,c,d,c1,d1)
	!
	!	computes the integral of the DV correction
	!
		real (kind=dp) :: a,b,c,d,c1,d1
		real (kind=dp) :: bb
		integer     :: fier, fneval, flast, flenw
		real (kind=dp) :: ni_DV, s, fabserr, ss
		real (kind=dp) :: fwork(400)
		integer, dimension(200) :: fiwork

		call dqagi(f_DVV,s,1,1.0d-18,1.0d-40,ni_DV,fabserr,fneval, &
		fier,100,flenw,flast,fiwork,fwork)

	END FUNCTION ni_DV

	FUNCTION Chi2FOV(a0,b0,c0,d0,c1,d1, atau)

		real (kind=dp) :: Chi2FOV
		real (kind=dp), dimension(binmin:78) :: FOV
		real (kind=dp) :: a0,b0,c0,d0,c1,d1,atau

		FOV = VFO(a0,b0,c0,d0,c1,d1, atau)
		Chi2FOV = dot_product(FOV(binmin:78),matmul(SigInvV(binmin:78,binmin:78),FOV(binmin:78)))

	END FUNCTION Chi2FOV

	SUBROUTINE set_param(a,b,c,d,c1,d1)

		real (kind=dp) :: a,b,c,d,c1,d1

		alpha	=	a
		beta	=	b
		gamma	=	c
		delta	=	d
		rho 	=	c1
		sigma	=	d1

	END SUBROUTINE set_param

	FUNCTION IntegralPTFO(bin,atau) result (iPTFO)
	!
	!	gives the integral of the FOPT expansion
	!	(see arXiv:0806.3156v1, eq. 3.2, p. 7)
	!
		real (kind=dp), intent(in) :: atau
		real (kind=dp) :: iPTFO, alphaS, s
		integer, intent(in) :: bin

		s = sbin(bin) + dsbin(bin) / 2.0_dp

		call Find_aS(atau, s, alphaS, 1.0d-14, 20, statenull)

		iPTFO = &
		s / (4.0_dp*pi**2) * (1.0_dp + c(1,1) * (alphaS / pi) + &
		(c(2,1) + 2.0_dp*c(2,2) * jj(1)) * (alphaS / pi)**2 + &
		(c(3,1) + 2.0_dp*c(3,2) * jj(1) + 3.0_dp*c(3,3) * jj(2)) * (alphaS / pi)**3 + &
		(c(4,1) + 2.0_dp*c(4,2) * jj(1) + 3.0_dp*c(4,3) * jj(2) + 4.0_dp*c(4,4) * jj(3)) * (alphaS / pi)**4 + &
		(c(5,1) + 2.0_dp*c(5,2) * jj(1) + 3.0_dp*c(5,3) * jj(2) + 4.0_dp*c(5,4) * jj(3) + 5.0_dp*c(5,5) * jj(4)) * (alphaS / pi)**5)

		!write(*,*) kind(iPTFO)

	END FUNCTION IntegralPTFO

    FUNCTION H(t)

		real(kind=dp) :: H,t

		H = -0.2222222222222222_dp/t + 0.3730002895803152_dp * ATan(0.1957622473346858_dp - 2.777520917064214_dp*t) - &
		0.3950617283950617_dp * Log(t) + 0.2727626771780712_dp * Log(0.3539687005190283_dp + 1.0_dp * t) + &
		0.06114952560849528_dp * Log(0.13459153249825306_dp - 0.14096185280332845_dp * t + 1.0_dp * t**2)

    END FUNCTION H

    FUNCTION HPrime(x)

        real (kind=dp)  ::  HPrime,x

        HPrime = -1.036016106380334_dp/(1._dp + (0.1957622473346858_dp - 2.777520917064214_dp*x)**2) + 0.2222222222222222_dp/x**2-&
        0.3950617283950617_dp/x + 0.2727626771780712_dp/(0.3539687005190283_dp + 1._dp*x) + &
        (0.06114952560849528_dp*(-0.14096185280332845_dp + 2._dp*x))/(0.13459153249825306_dp - 0.14096185280332845_dp*x +1._dp*x**2)

    END FUNCTION HPrime

    SUBROUTINE Find_aS(atau, S, x0, eps, MaxIter, statenull)

        real (kind=dp), intent(in)	::	eps, atau
        integer , intent(in)		::	MaxIter
        integer , intent(out)		::	statenull
        real (kind=dp) , intent(out)::  x0
        real (kind=dp)				::	x1, s
        integer						::	counter

        counter = 0

        statenull = 0  !!  statenull = 0, everything is fine
        x0 = atau + 0.1_dp

        !write(*,*) x0, atau, S

        x1 = x0 - (H(x0/pi) - H(atau/pi) + Log(s/mtau**2)/2._dp) * pi / HPrime(x0/pi)

        do while ( abs(x0 - x1) >= eps * abs(x1) )

			!write(*,*) "What the Fuck", counter, MaxIter

            if ( counter == MaxIter ) then

                statenull = 1  !!  Maximum number of iterations
                return

            end if

            x0 = x1
            x1 = x0 - (H(x0/pi) - H(atau/pi) + Log(s/mtau**2)/2._dp) * pi / HPrime(x0/pi)

            counter = counter + 1

            !write(*,*) x0, abs(x0 - x1)

        end do

        !write(*,*) x1, HPrime(x1/pi)

        x0 = x1

    END SUBROUTINE Find_aS

!	REAL(kind=dp) FUNCTION alphaS(bin, atau)
!	!
!	!	Computes alphas(m) via spline interpolation
!	!	spline is computed by values from the mathematica code of santi
!	!
!		integer, intent(in) :: bin
!		real(kind=dp), intent(in) :: atau
!		call spline(Xa,Ya(bin,0:2000),2001,Y2(bin,0:2000))
!		call splint(Xa,Ya(bin,0:2000),Y2(bin,0:2000),2001,atau,alphaS)

!	END FUNCTION alphaS

	SUBROUTINE fcn(npar, grad, fval, xval, iflag, futil)

		integer :: npar, iflag
		external :: futil
		real (kind=dp) :: grad(0:6), xval(0:6), fval

		fval = Chi2FOV(xval(0),xval(1),xval(2),xval(3),xval(4),xval(5),xval(6))

	END SUBROUTINE fcn

	FUNCTION spectral(s, atau)

		real (kind=dp), intent(in) :: atau
		real (kind=dp) :: alphaS, s, spectral
!		integer, intent(in) :: bin

!		s = sbin(bin) + dsbin(bin) / 2.0_dp

		call Find_aS(atau, s, alphaS, 1.0d-14, 20, statenull)

		spectral = &
		1.0d0/(4*pi**2) * ( c(0,1) + (alphaS/pi) * c(1,1) + (alphaS/pi)**2 * c(2,1) + (alphaS/pi)**3 * (c(3,1) - c(3,3)*pi**2) +&
		(alphaS/pi)**4 * (c(4,1) - c(4,3)*pi**2) + (alphaS/pi)**5 * (c(5,1) - c(5,3)*pi**2 + c(5,5)*pi**4) )

	END FUNCTION spectral

	FUNCTION fopt(s,atau)
	!
	!	gives the integral of the FOPT expansion
	!	(see arXiv:0806.3156v1, eq. 3.2, p. 7)
	!
		real (kind=dp), intent(in) :: atau
		real (kind=dp) :: fopt, alphaS, s

		call Find_aS(atau, s, alphaS, 1.0d-14, 20, statenull)

		fopt = &
		s / (4.0_dp*pi**2) * (1.0_dp + c(1,1) * (alphaS / pi) + &
		(c(2,1) + 2.0_dp*c(2,2) * jj(1)) * (alphaS / pi)**2 + &
		(c(3,1) + 2.0_dp*c(3,2) * jj(1) + 3.0_dp*c(3,3) * jj(2)) * (alphaS / pi)**3 + &
		(c(4,1) + 2.0_dp*c(4,2) * jj(1) + 3.0_dp*c(4,3) * jj(2) + 4.0_dp*c(4,4) * jj(3)) * (alphaS / pi)**4 + &
		(c(5,1) + 2.0_dp*c(5,2) * jj(1) + 3.0_dp*c(5,3) * jj(2) + 4.0_dp*c(5,4) * jj(3) + 5.0_dp*c(5,5) * jj(4)) * (alphaS / pi)**5)

		!write(*,*) kind(iPTFO)

	END FUNCTION fopt

	SUBROUTINE PltData(pnts, a, b)

		!
		!	th Data file -> [1:s, 2:specOPE, 3:specDV, 4:specALL, 5:FESR] >> data/PltDataThGEM.dat
		!
		!	ex Data file -> [1:s, 2:spec, 3:err_spec, 4:FESR, 5:err_FESR] >> data/PltDataExGEM.dat
		!

		real (kind=dp), dimension(0:pnts-1) :: spcdata, foptdata, dvdata, thspcdata, thdata
		real (kind=dp):: a, b, s
		integer:: i, pnts

		s = a

		open(55, file='data/PltDataThGEM.dat')

		do i=0,pnts-1
			spcdata(i) = spectral(s, pval(7))
			foptdata(i) = fopt(s, pval(7))
			call set_param(pval(1),pval(2),pval(3),pval(4),pval(5),pval(6))
			dvdata(i) =  f_DVV(s)
			thspcdata(i) = 2*pi**2 * (spcdata(i) + dvdata(i))
			thdata(i) = 2*pi**2/s * (foptdata(i) - ni_dv(s,pval(1),pval(2),pval(3),pval(4),pval(5),pval(6)))
			write(55,*) s, 2*pi**2 * spcdata(i), 2*pi**2 * dvdata(i), thspcdata(i), thdata(i)
			s = a + (i+1)*(b-a)/(pnts-1)
		end do

		close(55)
		open(55, file='data/PltDataExGEM.dat')

		do i=0,78
			write(55,*) &
			sbin(i), &
			spec(i), &
			sqrt(truecovv(i,i)), &
			2*pi**2/(sbin(i)+0.5d0*dsbin(i)) * fesr(i), &
			2*pi**2/(sbin(i)+0.5d0*dsbin(i)) * sqrt(fesrcovv(i,i))
		end do

		close(55)

	END SUBROUTINE PltData

	SUBROUTINE contour(par1, par2, ngrid)

		integer :: par1, par2, ngrid, m, n, state
		real (kind=dp) :: s1ax1, s1bx1, s1ax2, s1bx2, s2ax1, s2bx1, s2ax2, s2bx2, fn, fngrid
		real (kind=dp), dimension(:), allocatable :: x1, x2
		real (kind=dp), dimension(:,:), allocatable :: z
		real (kind=dp), dimension(13) :: parr
		real (kind=dp), dimension(2) :: cntrlvl

		allocate(x1(0:ngrid), stat=state)
		allocate(x2(0:ngrid), stat=state)
		allocate(z(0:ngrid,0:ngrid), stat=state)

		fngrid = 1.0d0*ngrid

		parr = pval
		s1ax1 = pval(par1)-ierr(par1)
		s1ax2 = pval(par2)-ierr(par2)
		s1bx1 = pval(par1)+ierr(par1)
		s1bx2 = pval(par2)+ierr(par2)
		parr(par1) = s1bx1
		parr(par2) = s1bx2
		cntrlvl(1) = Chi2FOV(parr(1),parr(2),parr(3),parr(4),parr(5),parr(6),parr(7))
		write(*,*) '1st cntr at par1,par2= ',s1bx1,s1bx2,' is: ',cntrlvl(1)
		parr = pval
		s2ax1 = pval(par1)-2.0d0*ierr(par1)
		s2ax2 = pval(par2)-2.0d0*ierr(par2)
		s2bx1 = pval(par1)+2.0d0*ierr(par1)
		s2bx2 = pval(par2)+2.0d0*ierr(par2)
		parr(par1) = s2bx1
		parr(par2) = s2bx2
		cntrlvl(2) = Chi2FOV(parr(1),parr(2),parr(3),parr(4),parr(5),parr(6),parr(7))
		write(*,*) '2nd cntr at par1,par2 = ',s2bx1,s2bx2,' is: ',cntrlvl(2)

		parr = pval

		do n = 0, ngrid
			fn = 1.0d0*n
			x1(n) = s2ax1 + fn*(s2bx1-s2ax1)/fngrid
			x2(n) = s2ax2 + fn*(s2bx2-s2ax2)/fngrid
		end do

		write(*,*) 'CREATING DATASET FOR THE CONTOUR. GRID = ',ngrid

		do m = 0, ngrid
			parr(par1) = x1(m)
			write(*,*) '***',m,'/',ngrid,'***'
			do n = 0, ngrid
				parr(par2) = x2(n)
				z(m,n) = Chi2FOV(parr(1),parr(2),parr(3),parr(4),parr(5),parr(6),parr(7))
			end do
		end do

		parr = pval

		call conrec(z,0,ngrid,0,ngrid,x1,x2,2,cntrlvl)

	END SUBROUTINE

	subroutine readConfig(par, npar)

		real (kind=dp), dimension(:), allocatable, intent(inout) :: par
		integer, intent(in) :: npar
		integer :: n, m , io, iscmt=0, length=0, i=0
		character p
		character*200 buffer
		character(:), allocatable :: line
		logical fexist

		open(553, file="GEM1.in")

		do
			read(553,'(A)',iostat=io) buffer
			if (io .ne. 0) exit
			do n=1,200
				p=buffer(n:n)
				if (p .eq. '#') then
					iscmt = 1
					exit
				end if
				length = length + 1
			end do
			if (len(trim(buffer)) .eq. 0) length = 0
			if (length .gt. 0) then

				if (iscmt .eq. 1) then
					line = trim(buffer(1:length))
				else
					line = adjustl(buffer)

				endif

				length = len(line)

				do m=1,length
					p = line(m:m)
					if (p .eq. '=') then
						i = i+1
						read(line(m+1:length),*) par(i)
						exit
					end if
				end do
			end if
			if (i .eq. npar) exit
			iscmt = 0
			length=0
		end do

		if (i .lt. npar) stop(" - [ERROR] - not enough parameter defined in config file")

		close(553)

	end subroutine

END PROGRAM main
