module test

	integer, parameter :: dp=kind(1.0d0)
	real (kind=dp), dimension(:), allocatable :: par

end module

program config
	
	use test

	implicit none

	integer :: n, npar, state

	npar = 7
	
	allocate(par(npar), stat=state)
	
	call readConfig(par, npar)

	do n=1, npar	
	write(*,*) par(n)
	end do

	contains

	subroutine readConfig(par, npar)

		real (kind=dp), dimension(:), allocatable, intent(inout) :: par
		integer, intent(in) :: npar
		integer :: n, m , io, iscmt=0, length=0, i=0
		character p
		character*200 buffer
		character(:), allocatable :: line
		
		open(55, file="config.in")
		
		do
			read(55,*,iostat=io) buffer
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
		
		if (i .le. npar) stop(" - [ERROR] - not enough parameter defined in config file")
		
		close(55)
		
	end subroutine
			
end program	
