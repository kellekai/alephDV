PROGRAM ni
	
	use mcf_tipos
	use improper_cuadratura
	
	implicit none
	
	integer	::	n,m,p,q,r,s
	real (kind=dp)	::	rn,rm,rp,rq,rr,rs
	
	real (kind=dp)	::	a,b,c,d,e,f						! parameters
	
	real (kind=dp), parameter	::	aa = 1.475_dp			! lower bound	
	real (kind=dp)	::	bb = huge(1.0_dp)	! upper bound
	
	real (kind=dp)	::	ss								! integration value
	
	integer :: it
	real (kind=dp) :: error, flag
	
	logical :: valid = .false.
	
	!======= MAIN SECTION ==========================================================================================================
	
	rn = 0._dp
	rm = 0._dp
	rp = 0._dp
	rq = 0._dp
	rr = 0._dp
	rs = 0._dp
	n=1
	m=1
	p=1
	q=1
	r=1
	s=1
	
	call SetParam(3.0d0,4.0d0,5.0d0,5.0d0,8.0d0,0.30d0)
	call romberg_improper(DV, aa, bb, ss, 1,1.0d-3, .false., 3)
	write(*,*) ss
	call romberg_improper(DV, aa, bb, ss, 1,1.0d-7, .false., 9)	
	write(*,*) ss
	call QUANC8(DV,AA,bb,1.0d-7,1.0d-7,ss,error,it,FLAG)
	write(*,*) ss, DV(3.0d0)==0.0d0
	
		!call SetParam(0.5_dp,0.5_dp,1._dp,0.05_dp,0.1_dp,0.1_dp)
		!call romberg_improper(DV, aa, bb, ss, 1,1.0d-14, .false., 5)	
		!write(*,*) ss
	
	DO n=1,40
		rn = -20._dp + n
	DO m=1,20
		rm = 0._dp + m
	DO p=1,10
		rp = 0._dp + p
	DO q=1,10
		rq = 0._dp + q
	DO r=1,10
		rr = 0._dp + r	
	DO s=1,10
		rs = 0._dp + s	
		
		call SetParam(rn*1._dp,rm*1._dp,rp*0.5_dp,rq*0.5_dp,rr*0.8_dp,rs*0.1_dp)
		do while (abs(DV(0.5d0*(aa+bb))) == 0.0d0)
			bb=bb*0.99d0
		end do
		!call QUANC8(DV,aa,bb,1.0d-3,1.0d-3,ss,error,it,FLAG)
		!if (isnan(ss)) stop '"ss" is a NaN'
		call romberg_improper(DV, aa, bb, ss, 4,0.1d-3, .false., 9)
		!write(*,*) rn,rm,rp,rq,rr,rs
		bb=100._dp
	!rn*0.05_dp,rm*0.05_dp,rp*0.08_dp,rq*0.05_dp,rr*0.05_dp,rs*0.08_dp
	ENDDO
	ENDDO
	write(*,*) rn,rm,rp,rq,rr,rs
	ENDDO
	write(*,*) rn,rm,rp,rq,rr,rs
	ENDDO
	write(*,*) rn,rm,rp,rq,rr,rs
	ENDDO
	write(*,*) rn,rm,rp,rq,rr,rs	
	ENDDO
	
	
	!======= END MAIN SECTION ======================================================================================================
	
	CONTAINS	!===================================================================================================================
	
	FUNCTION ExpFunc(s)	result(DV_damp)
		
		real (kind=dp)	::	DV_damp
		real (kind=dp), intent(in)	::	s
		
		DV_damp = exp( -a -b * s ** c )
		
	END FUNCTION ExpFunc
	
	FUNCTION SinFunc(s)	result(DV_osci)
	
		real (kind=dp)	::	DV_osci
		real (kind=dp), intent(in)	::	s			
		
		DV_osci = sin( -d -e * s ** f )
		
	END FUNCTION SinFunc
	
	FUNCTION DV(s)	result(DV_all)
	
		real (kind=dp)	::	DV_all
		real (kind=dp), intent(in)	::	s	
		
		DV_all = exp( -a -b * s ** c ) * cos( d + e * s ** f )
		
	END FUNCTION DV
	
	SUBROUTINE SetParam(a_p, b_p, c_p, d_p, e_p, f_p)
	
		real (kind=dp), intent(in)	::	a_p, b_p, c_p, d_p, e_p, f_p
		
		a = a_p
		b = b_p
		c = c_p
		d = d_p
		e = e_p
		f = f_p
		
	END SUBROUTINE SetParam
	
END PROGRAM ni
