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
	
	call size2D_binary('DP_Smoothed_MATLAB.dat',M,N)
	
	! Allocate Domain Parameter Array
	allocate(A(M,N))
	! allocate(Astore(M,N,time))

	A = input2D_binary('DP_Smoothed_MATLAB.dat',M,N)
	
	!call output3D_binary(Astore,N,M,time)
	
	tf2 = 2000
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
		c(30:33,60:63) = 1
		c(:,1) = 1
		c(:,N) = 1
		
		FRight(2:(M-1),2:(N-1)) = ((c(2:(M-1),3:N)-c(2:(M-1),2:(N-1)))/h)
		FLeft(2:(M-1),2:(N-1)) = ((c(2:(M-1),2:(N-1))-c(2:(M-1),1:(N-2)))/h)
		FUp(2:(M-1),2:(N-1)) = ((c(1:(M-2),2:(N-1))-c(2:(M-1),2:(N-1)))/h)
		FDown(2:(M-1),2:(N-1)) = ((c(2:(M-1),2:(N-1))-c(3:M,2:(N-1)))/h)
		
		c = c + (dt/A)*((FRight*ARight-FLeft*ALeft)/h + (FUp*AUp-FDown*ADown)/h)
		
		cstore(:,:,i) = c*A

	end do
	
	call output3D_binary('MSU_Conc_Data.dat',cstore,M,N,time2)
	
end Program MSU_Diffusion
