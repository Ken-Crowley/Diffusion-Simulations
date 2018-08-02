Program Diffusion_Smoothing

use array_io
implicit none

	integer :: M, N, time, i
	real*8 :: dt, t0, tf1, h
	real*8, allocatable :: A(:,:), Astore(:,:,:)
	
	! Parameters
	h = 1	
	dt = .1*h**2
	t0 = 0
	tf1 = 2
	time = floor((tf1-t0)/dt)
	
	! Grabs size of 2D array stores into M and N
	call size2D_binary('MSU_DP.dat',M,N)

	! Allocate Domain Parameter Arrays
	allocate(A(M,N))
	allocate(Astore(M,N,time))
	
	! Stores binary file to A
	A = input2D_binary('MSU_DP.dat',M,N)
	
	! Diffusion Smoothing
	do i = 1,time
		A = A + dt*dcdt_DP(A,h)
		Astore(:,:,i) = A
	end do
	
	! Writes diffusion smoothed data to file
	call output3D_binary('DP_Smoothed_Animate.dat',Astore,M,N,time)	
	call output2D_binary('DP_Smoothed.dat',A,M,N)
	
contains
	
	function dcdt_DP(A,h) result(smooth)
	implicit none
	
	real*8, intent(in), dimension(:,:) :: A
	real*8, intent(in) :: h
	real*8, dimension(1:size(A,1),1:size(A,2)) :: smooth
	
	! Middle Section
	smooth(2:(M-1),2:(N-1)) = (A(2:(M-1),3:N) + A(2:(M-1),1:(N-2)) &
	+ A(1:(M-2),2:(N-1)) + A(3:M,2:(N-1)) -4*A(2:(M-1),2:(N-1)))/h**2
	
	! Left Section
	smooth(2:(M-1),1) = (2*A(2:(M-1),2) + A(1:(M-2),1) + A(3:M,1) - 4*A(2:(M-1),1))/h**2
	
	! Right Section
	smooth(2:(M-1),N) = (2*A(2:(M-1),N-1) + A(1:(M-2),N) + A(3:M,N) - 4*A(2:(M-1),N))/h**2
	
	! Top Section
	smooth(1,2:(N-1)) = (A(1,1:(N-2)) + A(1,3:N) + 2*A(2,2:(N-1)) - 4*A(1,2:(N-1)))/h**2
	
	! Bottom Section
	smooth(M,2:(N-1)) = (A(M,1:(N-2)) + A(M,3:N) + 2*A(M-1,2:(N-1)) - 4*A(M,2:(N-1)))/h**2
	
	! NW Corner
	smooth(1,1) = (2*A(1,2) + 2*A(2,1) - 4*A(1,1))/h**2
	
	! SW Corner
	smooth(M,1) = (2*A(M,2) + 2*A((M-1),1) - 4*A(M,1))/h**2
	
	! NE Corner
	smooth(1,N) = (2*A(1,(N-1)) + 2*A(2,N) - 4*A(1,N))/h**2
	
	! SE Corner
	smooth(M,N) = (2*A(M,(N-1)) + 2*A((M-1),N) - 4*A(M,N))/h**2
	end function dcdt_DP
end Program Diffusion_Smoothing