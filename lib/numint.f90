!---------------------------------------------------------------
!Las siguientes rutinas tratan de integrar integrales impropias
!tal y como se desarrolla en:
!
!Numerical Recipes in FORTRAN. The Art of Scientific Computing
!
!Y estan basadas en las subrutinas:
!   qromo    (Metodo de Romberg con la formula del punto medio)
!   midpnt   (Integracion de punto medio. No requiere evaluacion en los
!             extremos)
!   midinf   (para integrales con limites a>0, b->infinito o
!             a->-infinito, b<0)
!   midsql   (para integrales con una singularidad del tipo
!             1/sqrt(a))
!   midsqu   (para integrales con una singularidad del tipo
!             1/sqrt(b))
!   midexp   (asume que b->infinito y que la funcion decrece
!             exponencialmente en el infinito)
!
!USO:
!             use improper_cuadratura
!             real(kind=sp o dp)  :: a,b,ss
!             integer             :: integ_typ,order
!             logical             :: partial
!
!             call romberg_improper (f,a,b,ss,integ_typ,precision,partial,order)
!
!integ_typ es opcional. Posibles valores:
!                              1 => midinf
!                              2 => midsql
!                              3 => midsqu
!                              4 => midexp
!Ausencia o cualquier otro valor => midpnt
!
!precision es opcional. Por defecto es 1.0e-07 (1.0e-14) en simple (doble) precision.
!
!partial es opcional. .true. permite ver valores intermedios
!                            de la integracion. Por defecto es .false.
!
!order es opcional. Es el orden de interpolacion. Por defecto es 5.
!
!INTERFACE DE LA FUNCION:
!             interface 
!                 function f(x) result (f_x)
!                 real(kind=sp o dp), intent(in) :: x
!                 real(kind=sp o dp)             :: f_x
!                 end function func
!             end interface
!
!Un limite infinito apropiado puede conseguirse con la funcion
!intrinseca huge.
!
module improper_cuadratura_sp

use mcf_tipos

private

public :: romberg_improper

interface romberg_improper
    module procedure qromo_sp
end interface

private :: qromo_sp
private :: midpnt_sp
private :: func_sp

contains
!
!Esta funcion realiza de manera transparente los cambios de variable
!apropiados
!
function func_sp(x,funk,aa,bb,int_typ) result(func_x)

use mcf_tipos

integer, intent(in)       :: int_typ
real(kind=sp), intent(in) :: x,aa,bb
real(kind=sp)             :: func_x

interface 
   function funk(x) result (funk_x)
   use mcf_tipos
   real(kind=sp), intent(in)    :: x
   real(kind=sp)                :: funk_x
   end function funk
end interface

select case (int_typ)
	case (1)         !midinf
       func_x=funk(1.0_sp/x)/x**2.0_sp
	case (2)         !midsql
       func_x=2.0_sp*x*funk(aa+x**2.0_sp)
	case (3)         !midsqu
       func_x=2.0_sp*x*funk(bb-x**2.0_sp)
	case (4)         !midexp
	 func_x=funk(-log(x))/x
	case default     !midpnt
	 func_x=funk(x)
end select

end function func_sp
!
!Integracion con la formula compuesta del punto medio. Sucesivas
!llamadas dividen el intervalo de integracion por 3, para aprovechar
!calculos previos
!
subroutine midpnt_sp(funk,aa,bb,s,n,integ_type)

integer, intent(in)             :: n,integ_type
real(kind=sp), intent(in)       :: aa,bb
real(kind=sp), intent(out)      :: s
integer                         :: it,j
real(kind=sp)                   :: a,b,ddel,del,sum,tnm,x

interface 
   function funk(x) result (funk_x)
   use mcf_tipos
   real(kind=sp), intent(in)    :: x
   real(kind=sp)                :: funk_x
   end function funk
end interface

select case (integ_type)
	case (1)         !midinf
	   a=1.0_sp/bb
	   b=1.0_sp/aa
	case (2,3)       !midsql, midsqu
	   a=0.0_sp
	   b=sqrt(bb-aa)
	case (4)         !midexp
	   a=0.0_sp
	   b=exp(-aa)
	case default     !midpnt
	   a=aa
	   b=bb
end select

if (n==1) then
       s=(b-a)*func_sp(0.5_sp*(a+b),funk,aa,bb,integ_type)
else
       it=3**(n-2)
       tnm=it
       del=(b-a)/(3.0_sp*tnm)
       ddel=del+del
       x=a+0.5_sp*del
       sum=0.0_sp
       do j=1,it
             sum=sum+func_sp(x,funk,aa,bb,integ_type)
             x=x+ddel
             sum=sum+func_sp(x,funk,aa,bb,integ_type)
             x=x+del
       end do
       s=(s+(b-a)*sum/tnm)/3.0_sp
endif

end subroutine midpnt_sp
!
!Integracion de romberg en un intervalo abierto
!
subroutine qromo_sp(func,a,b,ss,int_typ,precision,partial,order)

use mcf_interpoli 

real(kind=sp), intent(in), optional :: precision
logical, intent(in), optional       :: partial
integer, intent(in), optional       :: order,int_typ
real(kind=sp), intent(in)           :: a,b
real(kind=sp), intent(out)          :: ss
logical, parameter                  :: partial_default=.false.
integer, parameter                  :: jmax=20, jmaxp=jmax+1, k_default=5
integer, parameter                  :: integral_type_default=-1
real(kind=sp), parameter            :: precision_default=1.0e-07
logical                             :: see_partial
integer                             :: j,jmin,k,integral_type
real(kind=sp)                       :: dss, eps
real(kind=sp), dimension(jmaxp)     :: h,s

interface 
    function func(x) result (f_x)
    use mcf_tipos
    real(kind=sp), intent(in) :: x
    real(kind=sp)             :: f_x
    end function func
end interface

if (present(int_typ)) then
   integral_type=int_typ
else
   integral_type=integral_type_default
endif
if (present(precision) .and. (precision >= precision_default)) then
   eps = precision
else
   eps = precision_default
end if

if (present(partial)) then
   see_partial = partial
else
   see_partial = partial_default
end if

if (present(order)) then
   k = order
else
   k = k_default
endif

h(1)=1.0_sp

do j=1,jmax
        call midpnt_sp(func,a,b,s(j),j,integral_type)

        if (see_partial) then
           print "(a,i3,a,es25.15)","N=",j," estimacion de la integral= ",s(j)
        end if

        if (j>=k) then
          jmin = j - k + 1
          call polint(h(jmin:j),s(jmin:j),k,0.0_sp,ss,dss)

          if (see_partial) then
             print "(2(a,es25.15))","---->Valor interpolado= ",s(j), &
             " error estimado= ",dss
          end if

          if (abs(dss)<=eps*abs(ss)) then
              exit
          end if
        end if
        s(j+1)=s(j)
!la intepolacion es en h**2 dado el error de la regla del 
!punto medio
        h(j+1)=h(j)/9.0_sp 
end do
!Nunca se deberia llegar aqui
if (j>jmax) then
    stop "Sobrepasado el numero maximo de iteraciones en Romberg"
end if

end subroutine qromo_sp

end module improper_cuadratura_sp
!
!---------------------------------------------------------------
!
module improper_cuadratura_dp

use mcf_tipos

private

public :: romberg_improper

interface romberg_improper
    module procedure qromo_dp
end interface

private :: qromo_dp
private :: midpnt_dp
private :: func_dp

contains
!
!Esta funcion realiza de manera transparente los cambios de variable
!apropiados
!
function func_dp(x,funk,aa,bb,int_typ) result(func_x)

use mcf_tipos

integer, intent(in)       :: int_typ
real(kind=dp), intent(in) :: x,aa,bb
real(kind=dp)             :: func_x

interface 
   function funk(x) result (funk_x)
   use mcf_tipos
   real(kind=dp), intent(in)    :: x
   real(kind=dp)                :: funk_x
   end function funk
end interface

select case (int_typ)
	case (1)         !midinf
       func_x=funk(1.0_dp/x)/x**2.0_dp
	case (2)         !midsql
       func_x=2.0_dp*x*funk(aa+x**2.0_dp)
	case (3)         !midsqu
       func_x=2.0_dp*x*funk(bb-x**2.0_dp)
	case (4)         !midexp
	   func_x=funk(-log(x))/x
	case default     !midpnt
	   func_x=funk(x)
end select

end function func_dp
!
!Integracion con la formula compuesta del punto medio. Sucesivas
!llamadas dividen el intervalo de integracion por 3, para aprovechar
!calculos previos
!
subroutine midpnt_dp(funk,aa,bb,s,n,integ_type)

integer, intent(in)             :: n,integ_type
real(kind=dp), intent(in)       :: aa,bb
real(kind=dp), intent(out)      :: s
integer                         :: it,j
real(kind=dp)                   :: a,b,ddel,del,sum,tnm,x

interface
   function funk(x) result (funk_x)
   use mcf_tipos
   real(kind=dp), intent(in)    :: x
   real(kind=dp)                :: funk_x
   end function funk
end interface

select case (integ_type)
	case (1)         !midinf
	   a=1.0_dp/bb
	   b=1.0_dp/aa
	case (2,3)       !midsql, midsqu
	   a=0.0_dp
	   b=sqrt(bb-aa)
	case (4)         !midexp
	   a=0.0_dp
	   b=exp(-aa)
	case default     !midpnt
	   a=aa
	   b=bb
end select

if (n==1) then
       s=(b-a)*func_dp(0.5_dp*(a+b),funk,aa,bb,integ_type)
       !if (s==0._dp .or. isnan(s)) then
	!	write(*,*) '[WARNING], x= ', 0.5_dp*(a+b), 's= ', s, 'romberg is unable to compute the integral'
	 !  end if
else
       it=3**(n-2)
       tnm=it
       del=(b-a)/(3.0_dp*tnm)
       ddel=del+del
       x=a+0.5_dp*del
       sum=0.0_dp
       do j=1,it
             sum=sum+func_dp(x,funk,aa,bb,integ_type)
!			 if (func_dp(x,funk,aa,bb,integ_type)==0._dp .or. isnan(func_dp(x,funk,aa,bb,integ_type))) then
!				write(*,*) 'Zeile 351 in midpnt: n=', n, ', x= ', x, 'wrong function value s= ', func_dp(x,funk,aa,bb,integ_type)
!			 end if
             x=x+ddel
             sum=sum+func_dp(x,funk,aa,bb,integ_type)
!			 if (func_dp(x,funk,aa,bb,integ_type)==0._dp .or. isnan(func_dp(x,funk,aa,bb,integ_type))) then
!				write(*,*) 'Zeile 356 in midpnt: n=', j, ', x= ', x, 'wrong function value s= ', func_dp(x,funk,aa,bb,integ_type)
!			 end if
             x=x+del
       end do
       s=(s+(b-a)*sum/tnm)/3.0_dp
endif

end subroutine midpnt_dp
!
!Integracion de romberg en un intervalo abierto
!
subroutine qromo_dp(func,a,b,ss,int_typ,precision,partial,order)

use mcf_interpoli 

real(kind=dp), intent(in), optional :: precision
logical, intent(in), optional       :: partial
integer, intent(in), optional       :: order,int_typ
real(kind=dp), intent(in)           :: a,b
real(kind=dp), intent(out)          :: ss
logical, parameter                  :: partial_default=.false.
integer, parameter                  :: jmax=20, jmaxp=jmax+1, k_default=5
integer, parameter                  :: integral_type_default=-1
real(kind=dp), parameter            :: precision_default=1.0e-14
logical                             :: see_partial
integer                             :: j,jmin,k,integral_type
real(kind=dp)                       :: dss, eps
real(kind=dp), dimension(jmaxp)     :: h,s

interface 
    function func(x) result (f_x)
    use mcf_tipos
    real(kind=dp), intent(in) :: x
    real(kind=dp)             :: f_x
    end function func
end interface

if (present(int_typ)) then
   integral_type=int_typ
else
   integral_type=integral_type_default
endif

if (present(precision) .and. (precision >= precision_default)) then
   eps = precision
else
   eps = precision_default
end if

if (present(partial)) then
   see_partial = partial
else
   see_partial = partial_default
end if

if (present(order)) then
   k = order
else
   k = k_default
endif

h(1)=1.0_dp

do j=1,jmax
        call midpnt_dp(func,a,b,s(j),j,integral_type)

        if (see_partial) then
           print "(a,i3,a,es25.15)","N=",j," estimacion de la integral= ",s(j)
        end if

        if (j>=k) then
          jmin = j - k + 1
          call polint(h(jmin:j),s(jmin:j),k,0.0_dp,ss,dss)

          if (see_partial) then
             print "(2(a,es25.15))","---->Valor interpolado= ",s(j), &
             " error estimado= ",dss
          end if
          if (abs(dss)<=eps*abs(ss)) then
              exit
          end if
        end if
        s(j+1)=s(j)
!la intepolacion es en h**2 dado el error de la regla del 
!punto medio
        h(j+1)=h(j)/9.0_dp
end do
!Nunca se deberia llegar aqui
if (j>jmax) then
    stop "Sobrepasado el numero maximo de iteraciones en Romberg"
end if

end subroutine qromo_dp

end module improper_cuadratura_dp
!
!---------------------------------------------------------------
!
module improper_cuadratura

use improper_cuadratura_sp
use improper_cuadratura_dp

end module improper_cuadratura
!
!---------------------------------------------------------------
!
