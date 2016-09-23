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

	integer			:: n,m, nn,state, nbin, bestbin, EVierr  ! used in loops

	real (kind=dp)	:: alphaV,betaV,gammaV,deltaV, rhoV, sigmaV, dummyR, pvalue  ! parameters

	real (kind=dp)	:: p, q, df, bound

	real (kind=dp)	:: parr(5)

	integer			:: integ_typ, order, status ! ni type, ni order
	logical			:: partial  ! ni print steps yes/no ni

	integer			:: perr	! used by minuit
	real (kind=dp)	:: arg(13),pval(13),ierr(13),plo(13),phi(13)! used by minuit
	character*20 name(13)

	real (kind=dp), dimension(:), allocatable :: WR,WI
	real (kind=dp), dimension(:,:), allocatable :: Z

	allocate(WR(79-binmin), stat=statearr)
	allocate(WI(79-binmin), stat=statearr)
	allocate(Z(79-binmin,79-binmin), stat=statearr)

	!============STARTING MAIN CODE FROM HERE=======================================================================================

	call MNinit(5,6,7)
	call MNsetI('Fit of ALEPH Data to compute alpha s and the DV Parameters')

	npar = 5

	npar = npar + 1

	allocate(par(npar), stat=state)

	call readConfig(par, npar)

	nbin=par(npar)

	call init(nbin)

	pval(1:npar) = par(1:npar)

	call MNparm(1,'alpha',pval(1),0.1d-4,0.0d0,0.0d0,perr)
	call MNparm(2,'beta',pval(2),0.1d-4,0.0d0,0.0d0,perr)
	call MNparm(3,'gamma',pval(3),0.1d-4,0.0d0,0.0d0,perr)
	call MNparm(4,'delta',pval(4),0.1d-4,0.0d0,0.0d0,perr)
	call MNparm(5,'atau',pval(5),0.1d-4,0.0d0,0.0d0,perr)
	call mnexcm(fcn,'SET PRINTout', (/1.0_dp/),1,perr,0)  ! set print level (2=full)

	call mnexcm(fcn, 'HESse', (/ 2._dp /), 0, perr, 0)
	call mnexcm(fcn, 'SET STRategy', (/ 2._dp /), 1, perr, 0)
	call mnexcm(fcn, 'MIN', (/ 10000.0_dp, 0.1d-5 /), 2, perr, 0)

	do n=1,13
		call mnpout(n,name(n),pval(n),ierr(n),plo(n),phi(n),perr)
	enddo
	do nn=1,5
	write(*,*) name(nn),pval(nn),ierr(nn),kind(pval(nn))
	enddo
	call cdfchi ( 1, p, q, chi2fov(pval(1),pval(2),pval(3),pval(4),pval(5)), 1._dp*(79-binmin-5), status, bound )
	write(*,*) '------------------'
	write(*,*) '##	fitbin/Energy = ',binmin,' / ',sbin(binmin)+0.5_dp*dsbin(binmin)
	write(*,*) '##	dof = ',79-binmin-5
	write(*,*) '------------------'
	call mnexcm(fcn, 'SHOw FCNvalue', arg, 0, perr, 0)
	write(*,*) '------------------'
	write(*,*) 'P-Value = ',100d0*q



	!============END OF MAIN CODE===================================================================================================



	CONTAINS !======================================================================================================================

	FUNCTION VFO(a,b,c,d, atau)
	!
	!	produces vector-like array, with entries
	!	VFO[k]=I_th(x_k)-I_exp(x_k)
	!
		real (kind=dp) :: a,b,c,d, atau
		real (kind=dp), dimension(binmin:78) :: VFO

		do n = binmin, 78
			VFO(n) = IntegralPTFO(n, atau) - ni_DVV(n, a,b,c,d) - FESR(n)
		end do

	END FUNCTION VFO

	FUNCTION f_DVV(s)
	!
	!	represents the DV term)
	!
		real (kind=dp) :: f_DVV
		real (kind=dp), intent(in) :: s

		f_DVV =	exp(-alphaV - betaV * s) * sin(gammaV + deltaV * s)

	END FUNCTION f_DVV

	FUNCTION ni_DVV(bin,a,b,c,d)
	!
	!	computes the integral of the DV term
	!
		real (kind=dp) :: a,b,c,d,e,f
		real (kind=dp) :: ni_DVV, s
		integer, intent(in) :: bin

		s = sbin(bin)+dsbin(bin)/2._dp

		ni_DVV = exp(-a-b*s)*(d*cos(c+d*s)+b*sin(c+d*s))/(b**2+d**2)

	END FUNCTION ni_DVV

	FUNCTION Chi2FOV(a,b,c,d, atau)

		real (kind=dp) :: Chi2FOV
		real (kind=dp), dimension(binmin:78) :: FOV
		real (kind=dp) :: a,b,c,d,e,f,atau

		FOV = VFO(a,b,c,d, atau)
		Chi2FOV = dot_product(FOV(binmin:78),matmul(SigInvV(binmin:78,binmin:78),FOV(binmin:78)))

	END FUNCTION Chi2FOV

	SUBROUTINE set_param(a,b,c,d)

		real (kind=dp) :: a,b,c,d

		alphaV	=	a
		betaV	=	b
		gammaV	=	c
		deltaV	=	d

	END SUBROUTINE set_param

	FUNCTION ni_dv(s,a,b,c,d)
	!
	!	computes the integral of the DV correction
	!
		real (kind=dp) :: a,b,c,d,e,f
		integer     :: fier, fneval, flast, flenw
		real (kind=dp) :: ni_dv, s, fabserr
		real (kind=dp) :: fwork(400)
		integer, dimension(200) :: fiwork

		ni_DV = exp(-a-b*s)*(d*cos(c+d*s)+b*sin(c+d*s))/(b**2+d**2)

!		call set_param(a,b,c,d,e,f)

!		call dqagi(f_DVV,s,1,1.0d-18,1.0d-40,ni_dv,fabserr,fneval,fier,100,flenw,flast,fiwork,fwork)

	END FUNCTION ni_dv

	FUNCTION spectral(s, atau)

		real (kind=dp), intent(in) :: atau
		real (kind=dp) :: alphaS, s, spectral
!		integer, intent(in) :: bin

!		s = sbin(bin) + dsbin(bin) / 2.0_dp

		call Find_aS(atau, s, alphaS, 1.0d-14, 20, state)

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

		call Find_aS(atau, s, alphaS, 1.0d-14, 20, state)

		fopt = &
		s / (4.0_dp*pi**2) * (1.0_dp + c(1,1) * (alphaS / pi) + &
		(c(2,1) + 2.0_dp*c(2,2) * jj(1)) * (alphaS / pi)**2 + &
		(c(3,1) + 2.0_dp*c(3,2) * jj(1) + 3.0_dp*c(3,3) * jj(2)) * (alphaS / pi)**3 + &
		(c(4,1) + 2.0_dp*c(4,2) * jj(1) + 3.0_dp*c(4,3) * jj(2) + 4.0_dp*c(4,4) * jj(3)) * (alphaS / pi)**4 + &
		(c(5,1) + 2.0_dp*c(5,2) * jj(1) + 3.0_dp*c(5,3) * jj(2) + 4.0_dp*c(5,4) * jj(3) + 5.0_dp*c(5,5) * jj(4)) * (alphaS / pi)**5)

		!write(*,*) kind(iPTFO)

	END FUNCTION fopt


	FUNCTION IntegralPTFO(bin,atau) result (iPTFO)
	!
	!	gives the integral of the FOPT expansion
	!	(see arXiv:0806.3156v1, eq. 3.2, p. 7)
	!
		real (kind=dp), intent(in) :: atau
		real (kind=dp) :: iPTFO, alphaS, s
		integer, intent(in) :: bin

		s = sbin(bin) + dsbin(bin) / 2.0_dp

		call Find_aS(atau, s, alphaS, 1.0d-14, 20, state)

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

    SUBROUTINE Find_aS(atau, S, x0, eps, MaxIter, state)

        real (kind=dp), intent(in)	::	eps, atau
        integer , intent(in)		::	MaxIter
        integer , intent(out)		::	state
        real (kind=dp) , intent(out)::  x0
        real (kind=dp)				::	x1, s
        integer						::	counter

        counter = 0

        State = 0  !!  state = 0, everything is fine
        x0 = atau + 0.1_dp

        !write(*,*) x0, atau, S

        x1 = x0 - (H(x0/pi) - H(atau/pi) + Log(s/mtau**2)/2._dp) * pi / HPrime(x0/pi)

        do while ( abs(x0 - x1) >= eps * abs(x1) )

			!write(*,*) "What the Fuck", counter, MaxIter

            if ( counter == MaxIter ) then

                state = 1  !!  Maximum number of iterations
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

	SUBROUTINE fcn(npar, grad, fval, xval, iflag, futil)

		integer :: npar, iflag
		external :: futil
		real (kind=dp) :: grad(0:4), xval(0:4), fval

		fval = Chi2FOV(xval(0),xval(1),xval(2),xval(3),xval(4))

	END SUBROUTINE fcn

	SUBROUTINE PltData(pnts, a, b)

		real (kind=dp), dimension(0:pnts-1) :: spcdata, foptdata, dvdata, thspcdata, thdata
		real (kind=dp):: a, b, s
		integer:: i, pnts

		s = a

		open(55, file='data/PltDataThFirstMod.dat')

		do i=0,pnts-1
			spcdata(i) = spectral(s, parr(5))
			foptdata(i) = fopt(s, parr(5))
			call set_param(parr(1),parr(2),parr(3),parr(4))
			dvdata(i) =  f_DVV(s)
			thspcdata(i) = 2._dp*pi**2 * (spcdata(i) + dvdata(i))
			thdata(i) = 2._dp*pi**2/s * (foptdata(i) - ni_dv(s,parr(1),parr(2),parr(3),parr(4)))
			write(55,*) s, 2*pi**2 * spcdata(i), 2*pi**2 * dvdata(i), thspcdata(i), thdata(i)
			s = a + (i+1)*(b-a)/(pnts-1)
		end do

		close(55)
		open(55, file='data/PltDataExFirstMod.dat')

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

	subroutine readConfig(par, npar)

		real (kind=dp), dimension(:), allocatable, intent(inout) :: par
		integer, intent(in) :: npar
		integer :: n, m , io, iscmt=0, length=0, i=0
		character p
		character*200 buffer
		character(:), allocatable :: line
		logical fexist

		open(553, file="FM.in")

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
