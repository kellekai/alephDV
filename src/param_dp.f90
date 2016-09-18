MODULE param
	
	!
	! usage -> after including call subroutine init() to initialize experimental data arrays
	!
	use mcf_tipos
	use lu
	implicit none
	
	real (kind=dp), dimension(:), allocatable :: par
	integer npar
	
	integer :: binmin	! start fitting at...
	integer :: statearr
	
	!
	!	define constants
	!
	
	real (kind=dp), parameter :: pi = 4.0_dp * atan(1.0_dp)
	real (kind=dp), parameter   ::  mtau = 1.77684d0
	
	!
	! J function used in IIntegralPTFO defined as array. Values from python code alpha_s_dev.py
	!
	
	real (kind = dp), parameter, dimension(0:4) :: jj = &
		(/1.0d0, -1.0d0, -0.128986813369645d1, 0.386960440108935d1, 0.400340060244305d1/)
	
	!
	! c parameter used in IIntegralPTFO defined as array. Values from python code alpha_s_dev.py
	!
	
	real (kind = dp), parameter, dimension(0:5,1:5) :: c = &
		reshape( (/ 1.0d0, 1.0d0, 1.6398212048969865d0, 6.3710144831009652d0, 49.075700002948139d0, 283.0d0, &
		0.0d0, 0.0d0, -1.125d0, -5.6895977110182194d0, -33.09140661672037d0, -299.17700000000002d0, &
		0.0d0, 0.0d0, 0.0d0, 1.6875d0, 15.801594849790995d0, 129.578d0, &
		0.0d0, 0.0d0, 0.0d0, 0.0d0, -2.84765625d0, -40.616100000000003d0, &
		0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 5.1257799999999998d0 /), (/ 6,5 /) )
	
	!
	! FESR defined as array. Computed by python code alpha_s_dev
	!
	
	real (kind = dp), parameter, dimension(0:78) :: fesr = &
		(/ 4.06379172d-7, 7.56859046d-5, 1.30279781d-4, 2.13442231d-4 &
		, 3.30469793d-4, 5.02301792d-4, 6.96708699d-4, 9.26224280d-4 &
		, 1.20242714d-3, 1.55472531d-3, 2.00191260d-3, 2.61843486d-3 &
		, 3.35335545d-3, 4.27202758d-3, 5.45227928d-3, 6.98548627d-3 &
		, 8.92979936d-3, 1.13791954d-2, 1.44366779d-2, 1.78166659d-2 &
		, 2.12565637d-2, 2.44345925d-2, 2.72197438d-2, 2.95119712d-2 &
		, 3.13900976d-2, 3.29064701d-2, 3.41706011d-2, 3.52181631d-2 &
		, 3.61160078d-2, 3.68949000d-2, 3.75850849d-2, 3.81658194d-2 &
		, 3.86875106d-2, 3.91710423d-2, 3.95989649d-2, 3.99943933d-2 &
		, 4.03720009d-2, 4.07166545d-2, 4.10770160d-2, 4.14330085d-2 &
		, 4.17367497d-2, 4.20546290d-2, 4.23973398d-2, 4.26949981d-2 &
		, 4.30048739d-2, 4.33111376d-2, 4.36448568d-2, 4.39713147d-2 &
		, 4.43184133d-2, 4.46687920d-2, 4.50093897d-2, 4.53602298d-2 &
		, 4.56962737d-2, 4.60332629d-2, 4.63717227d-2, 4.67685901d-2 &
		, 4.71790947d-2, 4.75684199d-2, 4.79449926d-2, 4.83872808d-2 &
		, 4.87941829d-2, 4.92019228d-2, 4.96250111d-2, 5.00374790d-2 &
		, 5.04467282d-2, 5.13659768d-2, 5.23904530d-2, 5.33631633d-2 &
		, 5.44620820d-2, 5.56219118d-2, 5.68394122d-2, 5.98780429d-2 &
		, 6.33050338d-2, 6.68704081d-2, 7.01007191d-2, 7.67944947d-2 &
		, 8.29706685d-2, 9.08508594d-2, 1.01763446d-1 /)

	!
	! compute covariance matrix
	!
	
	real (kind = dp), dimension(:,:), allocatable :: SigInvV
	real (kind=dp), dimension(1:21,1:21) :: Test
	
	!
	! define and read in experimental data trucovV, sbin, dsbin
	!

	real (kind = dp), dimension(0:78,0:78) :: truecovv, fesrcovv
	
	real (kind = dp), dimension(0:78) :: sbin, spec

	real (kind = dp), dimension(0:78) :: dsbin
	
	real(kind=dp), dimension(0:2000) :: Xa
	
	real(kind=dp), dimension(0:78,0:2000) :: Ya, Y2
	
	CONTAINS !=========================================================
	
	SUBROUTINE init(bin) ! initialize experimental data arrays and spline data
		
		integer :: n,l,bin
		
		binmin = bin
	
		open(55, file = 'data/truecovV')
		do n = 0, 78
			read(55,*) truecovv(n,0:78)
		end do
		close(55)
		
		truecovv = (0.9987d0)**2*truecovv
	
		open(55, file = 'data/sbin')
		do n = 0, 78
			read(55,*) sbin(n)
		end do
		close(55)	
		
		open(55, file = 'data/specV')
		do n = 0, 78
			read(55,*) spec(n)
		end do
		close(55)
		
		spec = 0.9987*spec
	
		open(55, file = 'data/dsbin')
		do n = 0, 78
			read(55,*) dsbin(n)
		end do
		close(55)

		do n=0,78
			fesrcovv(n,n) = dot_product(dsbin(0:n),matmul(truecovV(0:n,0:n),dsbin(0:n)))/(2._dp*pi**2)**2
		end do
		
		open(55, file='data/spline_tics.dat')
		open(56, file='data/spline.dat')
	
		DO n=0,2000
			read(55,*) Xa(n)
		END DO
	
		DO n=0,78
			read(56,*) Ya(n,0:2000)
		END DO	
	
		close(55)
		close(56)
		
		allocate(SigInvV(binmin:78,binmin:78), stat=statearr)
		
		SigInvV(binmin:78,binmin:78) = INVM()
		!Test=T()
		
	
	END SUBROUTINE init

	FUNCTION INVM()
	!
	!	COMPUTES THE INVERSE COVARIANT MATRIX // noch nicht getestet, invertieralgorithmus finden.
	!
		integer :: i,j
		real (kind=dp), dimension(:,:), allocatable :: M, INVM
		real (kind=dp), dimension(:,:), allocatable :: A,B
		integer :: INDX(78-binmin+1),D,CODE
		
		allocate(A(1:78-binmin+1,1:78-binmin+1),stat=statearr)
		allocate(B(1:78-binmin+1,1:78-binmin+1),stat=statearr)
		allocate(M(binmin:78,binmin:78),stat=statearr)
		allocate(INVM(binmin:78,binmin:78),stat=statearr)

		do i = binmin, 78
			B(i+1-binmin,i+1-binmin) = 1._dp
			M(i,i) = dot_product(dsbin(0:i),matmul(truecovV(0:i,0:i),dsbin(0:i)))/(2._dp*pi**2)**2
			do j = i+1, 78
				B(i+1-binmin,j+1-binmin) = 0._dp
				B(j+1-binmin,i+1-binmin) = 0._dp
				M(i,j)=dot_product(dsbin(0:i),matmul(truecovV(0:i,0:j),dsbin(0:j)))/(2._dp*pi**2)**2
				M(j,i)=M(i,j)
			end do
		end do
		
		A(1:78-binmin+1,1:78-binmin+1) = M(binmin:78,binmin:78)
		
		call ludcmp(A,78-binmin+1,INDX,D,CODE)
		
		do i=1, 78-binmin+1
			call LUBKSB(A,78-binmin+1,INDX,B(i,1:78-binmin+1))
		end do
		
		INVM(binmin:78,binmin:78) = B(1:78-binmin+1,1:78-binmin+1)
	
	END FUNCTION INVM

END MODULE param
                                                                          
