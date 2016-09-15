PROGRAM OUT

integer io
character p
character*200 line
real*8 double
integer :: check=0

open(55, file="config.file")

do
	read(55,*,iostat=io) line
	if (io/=0) exit
	p = line(1:1)
	if (p.ne."#") then
		
		
end do

!read(*,*) in
!buf = LEN_TRIM(in)

!do n=1,buf
!	p = in(n:n)
!	if (p.ne.'=') then
!		check=0
!	else
!		check=1
!	end if
!	if (check.eq.1) then 
!		read(in(n+1:buf),*) double
!		exit
!	end if
!enddo

!write(*,*) double

END PROGRAM
