Program MSU_Diffusion

use array_io
implicit none

	integer :: M, N, time1, time2, i
	real*8 :: dt, t0, tf1, tf2, h
	real*8, allocatable :: A(:,:), Astore(:,:,:), ARight(:,:), ALeft(:,:), AUp(:,:), ADown(:,:)
	real*8, allocatable :: FRight(:,:), FLeft(:,:), FUp(:,:), FDown(:,:)
	real*8, allocatable :: c(:,:), cstore(:,:,:)
	
	! Parameters
	h = 1	
	dt = .1*h**2
	t0 = 0
	tf1 = .2
	time1 = floor((tf1-t0)/dt)
	
	call size2D_binary('MSUFixed.dat',M,N)
	
	! Allocate Domain Parameter Array
	allocate(A(M,N))
	! allocate(Astore(M,N,time))

	A = input2D_binary('MSUFixed.dat',M,N)
	
	do i = 1,time1
		A = A + dt*dcdt_DP(A,h)
		!Astore(:,:,i) = A
	end do
	
	!call ThreeD_Write(Astore,N,M,time)
	
	tf2 = 5
	time2 = floor((tf2-t0)/dt)
	
	A = A + .00001
	
	! Allocate Domain Parameter average arrays
	allocate(ARight(M,N))
	allocate(ALeft(M,N))
	allocate(AUp(M,N))
	allocate(ADown(M,N))

	ARight = 0
	ALeft = 0
	AUp = 0
	ADown = 0
	
	ARight(2:(M-1),2:(N-1)) = (A(2:(M-1),3:N)+A(2:(M-1),2:(N-1)))/2;
	ALeft(2:(M-1),2:(N-1))  = (A(2:(M-1),2:(N-1))+A(2:(M-1),1:(N-2)))/2;
	AUp(2:(M-1),2:(N-1))  = (A(1:(M-2),2:(N-1))+A(2:(M-1),2:(N-1)))/2;
	ADown(2:(M-1),2:(N-1))  = (A(2:(M-1),2:(N-1))+A(3:M,2:(N-1)))/2;

	allocate(FRight(M,N))
	allocate(FLeft(M,N))
	allocate(FUp(M,N))
	allocate(FDown(M,N))
	
	FRight = 0
	FLeft = 0
	FUp = 0
	FDown = 0
	
	allocate(c(M,N))
	allocate(cstore(M,N,time2))
	
	do i = 1,time2
		
		c(1,:) = 0
		c(M,:) = 0
		c(:,1) = 1
		c(:,N) = 0
		
		FRight(2:(M-1),2:(N-1)) = ((c(2:(M-1),3:N)-c(2:(M-1),2:(N-1)))/h)
		FLeft(2:(M-1),2:(N-1)) = ((c(2:(M-1),2:(N-1))-c(2:(M-1),1:(N-2)))/h)
		FUp(2:(M-1),2:(N-1)) = ((c(1:(M-2),2:(N-1))-c(2:(M-1),2:(N-1)))/h)
		FDown(2:(M-1),2:(N-1)) = ((c(2:(M-1),2:(N-1))-c(3:M,2:(N-1)))/h)
		
		c = c + (dt/A)*((FRight*ARight-FLeft*ALeft)/h + (FUp*AUp-FDown*ADown)/h)
		
		cstore(:,:,i) = c*A

	end do
	
	call output3D_binary('diffusiondata.dat',cstore,M,N,time2)
	
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
	
end Program MSU_Diffusion
