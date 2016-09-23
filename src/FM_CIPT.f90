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
	use cmplx_alpha2
	use complex_root

	implicit none

	integer			:: n,m, nn,state, nbin  ! used in loops

	real (kind=dp)	:: alphaV,betaV,gammaV,deltaV, rhoV, sigmaV, sss, sssatau  ! parameters

	real (kind=dp)	:: p, q, df, bound

	integer			:: integ_typ, order, status ! ni type, ni order
	logical			:: partial  ! ni print steps yes/no ni

	integer			:: perr	! used by minuit
	real (kind=dp)	:: arg(13),pval(13),ierr(13),plo(13),phi(13)! used by minuit
	character*20 name(13)

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

	call MNparm(1,'alpha',pval(1),0.1d0,0.0d0,0.0d0,perr)
	call MNparm(2,'beta',pval(2),0.1d0,0.0d0,0.0d0,perr)
	call MNparm(3,'gamma',pval(3),0.1d0,0.0d0,0.0d0,perr)
	call MNparm(4,'delta',pval(4),0.1d0,0.0d0,0.0d0,perr)
	call MNparm(5,'atau',pval(5),0.1d0,0.0d0,0.0d0,perr)
	call mnexcm(fcn,'SET PRINTout', (/1.0_dp/),1,perr,0)  ! set print level (2=full)
	call mnexcm(fcn, 'SET STRategy', (/ 2._dp /), 1, perr, 0)
	call mnexcm(fcn, 'HESse', (/ 2._dp /), 0, perr, 0)
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
			VFO(n) = IntegralCIPT(n, atau) - ni_DVV(n, a,b,c,d) - FESR(n)
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

		real (kind=dp) :: atau, spectral
		real (kind=dp) :: iCIPTspec, iCIPTspecRe, iCIPTspecIm, s, fabserr
		integer     :: fier, fneval, flast, flenw, flimit
		real (kind=dp) :: fwork(400)
		integer :: fiwork(200)

		call set_dummy(s,atau,0._dp)

		!call dqng(CIPT,0._dp,2*pi,1.0d-18,1.0d-40,iCIPT,fabserr,fneval,fier)
		!call dqagse(CIPT,0._dp,2*pi,epsabs,epsrel,limit,iCIPT,abserr,neval,ier,alist,blist,rlist,elist,iord,last)
		call dqag(CIPTspec,-pi,pi,1.0d-8,1.0d-10,1,iCIPTspec,fabserr,fneval,fier,100,flenw,flast,fiwork,fwork)

		spectral = (1.0d0/(8.0d0*pi**3)) * iCIPTspec

	END FUNCTION spectral

	function CIPTspec(angle)

		real (kind=dp) 				:: CIPTspec, x0, y0, e, a, x, y
		real (kind=dp), intent(in)	:: angle
		complex (kind=cdp)			:: DF,ZE
		integer 					:: n, k
		complex*16 RTS(10),FNV(10)
		integer NFOUND, IFAIL, maxits, known, nmore, km, guess
		real*8 ep1,ep2
		logical realrt

		CALL set_dummy(dummy(1),dummy(2),angle)

		realrt = .false.
		ep1=1.0d-4
		ep2=1.0d-8
		maxits = 100
		known=0
		nmore=1
		km=1
		guess=1

		RTS(1) = dcmplx(dummy(2),0.0d0)

		call ROOTS2(known,nmore,km,realrt,ep1,ep2,guess,maxits,RTS,NFOUND,FNV,IFAIL)

		!write(*,*) dummy(1),dummy(2),dummy(3),angle,"wtf"

		ZE = RTS(1)

		!write(*,*) 'x0 = ', ZE
		!write(*,*) 'FCV = ', C8_A(ZE)

		DF = &
		1.0_dp + c(1,1) * (ZE / pi) + &
		c(2,1) * (ZE / pi)**2 + &
		c(3,1) * (ZE / pi)**3 + &
		c(4,1) * (ZE / pi)**4 + &
		c(5,1) * (ZE / pi)**5

		!write(*,*) 'CIPTcmplx = ', DF

		CIPTspec = real (DF)

	end function CIPTspec

	FUNCTION IntegralCIPTcont(s,atau) result (iCIPTcont)
	!
	!	gives the integral of the FOPT expansion
	!	(see arXiv:0806.3156v1, eq. 3.2, p. 7)
	!
		real (kind=dp) :: atau
		real (kind=dp) :: iCIPTcont, s, fabserr
		integer     :: fier, fneval, flast, flenw, flimit
		real (kind=dp) :: fwork(400)
		integer :: fiwork(200)

		call set_dummy(s,atau,0._dp)

		call dqag(CIPT,0.0_dp,2.0_dp*pi,1.0d-8,1.0d-10,1,iCIPTcont,fabserr,fneval,fier,100,flenw,flast,fiwork,fwork)

	END FUNCTION IntegralCIPTcont


	FUNCTION IntegralCIPT(bin,atau) result (iCIPT)
	!
	!	gives the integral of the FOPT expansion
	!	(see arXiv:0806.3156v1, eq. 3.2, p. 7)
	!
		real (kind=dp) :: atau
		real (kind=dp) :: iCIPT, s, fabserr
		integer, intent(in) :: bin
		integer     :: fier, fneval, flast, flenw, flimit
		real (kind=dp) :: fwork(400)
		integer :: fiwork(200)

		s= sbin(bin) + dsbin(bin) / 2.0_dp

		call set_dummy(s,atau,0._dp)

		call dqag(CIPT,0.0_dp,2.0_dp*pi,1.0d-8,1.0d-10,1,iCIPT,fabserr,fneval,fier,100,flenw,flast,fiwork,fwork)

	END FUNCTION IntegralCIPT

	function CIPT(angle)

		real (kind=dp) 				:: CIPT, x0, y0, e, a, x, y
		real (kind=dp), intent(in)	:: angle
		complex (kind=cdp)			:: DF,ZE
		integer 					:: n, k
		complex*16 RTS(10),FNV(10)
		integer NFOUND, IFAIL, maxits, known, nmore, km, guess
		real*8 ep1,ep2
		logical realrt

		CALL set_dummy(dummy(1),dummy(2),angle)

		realrt = .false.
		ep1=1.0d-4
		ep2=1.0d-8
		maxits = 100
		known=0
		nmore=1
		km=1
		guess=1

		RTS(1) = dcmplx(dummy(2),0.0d0)

		call ROOTS(known,nmore,km,realrt,ep1,ep2,guess,maxits,RTS,NFOUND,FNV,IFAIL)

		ZE = RTS(1)

		DF = &
		1._dp / (4.0_dp*pi**2) * (1.0_dp + c(1,1) * (ZE / pi) + &
		c(2,1) * (ZE / pi)**2 + &
		c(3,1) * (ZE / pi)**3 + &
		c(4,1) * (ZE / pi)**4 + &
		c(5,1) * (ZE / pi)**5)

		CIPT = real ((exp(dcmplx(0._dp,angle))-1._dp)*DF)*dummy(1)/(-2*pi)

	end function CIPT

	COMPLEX (kind=cdp) FUNCTION C8_H(s)

		complex (kind=cdp)	:: s

		C8_H = -0.2222222222222222_dp/s + 0.3730002895803152_dp * atan(0.1957622473346858_dp - 2.777520917064214_dp*s) - &
		0.3950617283950617_dp * log(s) + 0.2727626771780712_dp * log(0.3539687005190283_dp + 1.0_dp * s) + &
		0.06114952560849528_dp * log(0.13459153249825306_dp - 0.14096185280332845_dp * s + 1.0_dp * s**2)

	END FUNCTION C8_H

    COMPLEX (kind=cdp) FUNCTION C8_A(x0)

		complex (kind=cdp) :: x0

		C8_A = (C8_H(x0/pi) - C8_H(cmplx(dummy(2),0._dp)/pi) + log(dummy(1)*exp(dcmplx(0._dp,dummy(3)))/mtau**2)/2._dp) * pi

    END FUNCTION C8_A

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

        x1 = x0 - (H(x0/pi) - H(atau/pi) + Log(s/mtau**2)/2._dp) * pi / HPrime(x0/pi)

        do while ( abs(x0 - x1) >= eps * abs(x1) )

            if ( counter == MaxIter ) then

                state = 1  !!  Maximum number of iterations
                return

            end if

            x0 = x1
            x1 = x0 - (H(x0/pi) - H(atau/pi) + Log(s/mtau**2)/2._dp) * pi / HPrime(x0/pi)

            counter = counter + 1

        end do

        x0 = x1

    END SUBROUTINE Find_aS

	SUBROUTINE fcn(npar, grad, fval, xval, iflag, futil)

		integer :: npar, iflag
		external :: futil
		real (kind=dp) :: grad(0:4), xval(0:4), fval

		fval = Chi2FOV(xval(0),xval(1),xval(2),xval(3),xval(4))

	END SUBROUTINE fcn

	SUBROUTINE PltData(pnts, a, b)

		!
		!	th Data file -> [1:s, 2:specOPE, 3:specDV, 4:specALL, 5:FESR] >> data/PltDataThFirstModCIPT.dat
		!
		!	ex Data file -> [1:s, 2:spec, 3:err_spec, 4:FESR, 5:err_FESR] >> data/PltDataExFirstModCIPT.dat
		!


		real (kind=dp), dimension(0:pnts-1) :: spcdata, ciptdata, dvdata, thspcdata, thdata
		real (kind=dp):: a, b, s
		integer:: i, pnts

		s = a

		open(55, file='data/PltDataThFirstModCIPT.dat')

		do i=0,pnts-1
			spcdata(i) = spectral(s, pval(5))
			ciptdata(i) = IntegralCIPTcont(s, pval(5))
			call set_param(pval(1),pval(2),pval(3),pval(4))
			dvdata(i) =  f_DVV(s)
			thspcdata(i) = 2*pi**2 * (spcdata(i) + dvdata(i))
			thdata(i) = 2*pi**2/s * (ciptdata(i) - ni_dv(s,pval(1),pval(2),pval(3),pval(4)))
			write(55,*) s, 2*pi**2 * spcdata(i), 2*pi**2 * dvdata(i), thspcdata(i), thdata(i)
			s = a + (i+1)*(b-a)/(pnts-1)
		end do

		close(55)
		open(55, file='data/PltDataExFirstModCIPT.dat')

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

		open(553, file="FM_CIPT.in")

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

			if (length .ge. 0) then
				if (iscmt .eq. 1) then
					line = trim(buffer(1:length))
				else
					line = trim(buffer)
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
