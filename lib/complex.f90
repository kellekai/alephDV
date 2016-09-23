module complex_root
use mcf_tipos
use param
IMPLICIT NONE

	integer, parameter			:: cdp=kind(dcmplx(0.0d0,0.0d0))
	complex (kind=cdp)			:: s, atau, ZE, ZS
	real	(kind=dp)			:: DS, HE, DE, HS, HM, DM
	integer						:: Numint
	
!      zs = ( 0.4d+00, 0.5d+00 )
!      hs = 0.25d+00
!      hm = 0.0000001d+00
!      dm = 0.00000001d+00
      
!      phi=0.3
!      s=-1.5d0
!      atau=0.3d0
	
!	call crf( ZS, HS, HM, DM, C8_A, DS, ZE, HE, DE, N )
	
!	write(*,*) ZE
		
	CONTAINS
    
	SUBROUTINE CRF ( ZS, HS, HM, DM, FUNC, DS, ZE, HE, DE, numint )
!
!  THE SUBROUTINE DETERMINES A ROOT OF A TRANSCEN-
!  DENTAL COMPLEX EQUATION F(Z)=0 BY STEP-WISE ITE-
!  RATION, (THE DOWN HILL METHOD).
!
!  INPUT-PARAMETERS.
!
!  ZS = START VALUE OF Z. (COMPLEX)
!  HS = LENGTH OF STEP AT START.
!  HM = MINIMUM LENGTH OF STEP.
!  DM = MINIMUM DEVIATION.
!
!  SUBPROGRAM.
!
!  FUNC(Z), A COMPLEX FUNCTION SUBPROGRAM FOR THE
!  CALCULATION OF THE VALUE OF F(Z) FOR A COMPLEX
!  ARGUMENT Z.
!
!  OUTPUT-PARAMETERS.
!
!  DS = CABS(FUNC(ZS)) = DEVIATION AT START.
!  ZE = END VALUE OF Z. (COMPLEX)
!  HE = LENGTH OF STEP AT END.
!  DE = CABS(FUNC(ZE)) = DEVIATION AT END.
!  N  = NUMBER OF ITERATIONS.
!
!  RESTRICTIONS.
!
!  THE FUNCTION W=F(Z) MUST BE ANALYTICAL IN THE
!  REGION WHERE ROOTS ARE SOUGHT.
!
      IMPLICIT NONE

      COMPLEX (kind=cdp) A
      COMPLEX (kind=cdp) CW
      REAL (kind=dp) DE
      REAL (kind=dp) DM
      REAL (kind=dp) DS
      COMPLEX (kind=cdp) FUNC
      REAL (kind=dp) H
      REAL (kind=dp) HE
      REAL (kind=dp) HM
      REAL (kind=dp) HS
      INTEGER I
      INTEGER K
      INTEGER numint
      INTEGER NR
      COMPLEX (kind=cdp) U(7)
      COMPLEX (kind=cdp) V
      REAL (kind=dp) W(3)
      REAL (kind=dp) W0
      COMPLEX (kind=cdp) Z(3)
      COMPLEX (kind=cdp) Z0
      COMPLEX (kind=cdp) ZE
      COMPLEX (kind=cdp) ZS

      U(1) = (  1.0E+00, 0.0E+00 )
      U(2) = (  0.8660254E+00, 0.5000000E+00 )
      U(3) = (  0.0000000E+00, 1.0000000E+00 )
      U(4) = (  0.9659258E+00, 0.2588190E+00 )
      U(5) = (  0.7071068E+00, 0.7071068E+00 )
      U(6) = (  0.2588190E+00, 0.9659258E+00 )
      U(7) = ( -0.2588190E+00, 0.9659258E+00 )
      H = HS
      Z0 = ZS
      numint = 0

      CW = FUNC ( Z0 )
      W0 = ABS ( REAL ( CW ) ) + ABS ( AIMAG ( CW ) )
      DS = W0
      IF ( W0 - DM ) 18, 18, 1
1     K = 1
      I = 0
2     V = ( -1.0E+00, 0.0E+00 )
3     A = ( -0.5E+00, 0.866E+00 )
4     Z(1) = Z0 + H * V * A
      CW = FUNC ( Z(1) )
      W(1) = ABS ( REAL ( CW ) ) + ABS ( AIMAG ( CW ) )
      Z(2) = Z0 + H * V
      CW = FUNC ( Z(2) )
      W(2) = ABS ( REAL ( CW ) ) + ABS ( AIMAG ( CW ) )
      Z(3) = Z0 + H * CONJG ( A ) * V
      CW = FUNC ( Z(3) )
      W(3) = ABS ( REAL ( CW ) ) + ABS ( AIMAG ( CW ) )
      numint = numint + 1
      IF ( W(1) - W(3) ) 5, 5, 6
5     IF ( W(1) - W(2) ) 7, 8, 8
6     IF ( W(2) - W(3) ) 8, 8, 9
7     NR = 1
      GO TO 10
8     NR = 2
      GO TO 10
9     NR = 3
10    IF ( W0 - W(NR) ) 11, 12, 12
11    GO TO ( 13, 14, 15 ), K
12    K = 1
      I = 0
      A = ( 0.707E+00, 0.707E+00 )
      V = ( Z(NR) - Z0 ) / H
      W0 = W(NR)
      Z0 = Z(NR)
      IF ( W0 - DM ) 18, 18, 4
13    K = 2
      IF ( H .LT. HM ) GO TO 18
      H = H * 0.25E+00
      GO TO 3
14    K = 3
      H = H * 4.0E+00
      GO TO 2
15    I = I + 1
      IF ( I - 7 ) 16, 16, 17
16    V = U(I)
      GO TO 3
17    IF ( H .LT. HM ) GO TO 18
      H = H * 0.25E+00
      I = 0
      GO TO 2
18    ZE = Z0
      HE = H
      DE = W0
      RETURN
      END
      
END module complex_root

