program Hilliard

	Use array_io
	implicit none
	
	integer :: M, Nx, Ny, t0, tf, time, i
	real*8 :: W, e, dt, dx, dy
	real*8, allocatable :: c(:,:), cstore(:,:,:)
	
	! Parameters
	W = 1
	M = 1
	e = 0.1
	
	! Grid Size
	Nx = 100
	Ny = 100
	dx = .1
	dy = .1
	
	! Time Restriction
	t0 = 0
	tf = 5
	dt = 0.0003
	time = floor((tf-t0)/dt)
	
	! Allocates array sizes
	allocate(c(Ny,Nx))
	allocate(cstore(Ny,Nx,time))
	
	! creates 2D array of random numbers between 0 and 1
	call random_number(c)
	
	! turns the random 2D array to have an average of 0.5
	c = 0.5 + 0.1*(c-0.5)*2
	
	! Loops through up to time and calculates the new concentration gradient
	do i = 1,time
	
	c = c + M*dt*(finite_dif_i(func(c,e,W,dx,dy),dx)+finite_dif_j(func(c,e,W,dx,dy),dy))
	cstore(:,:,i) = c
	
	end do
	
	! Writes 3D array to data file
	call output3D_binary('concdata.dat',cstore,Ny,Nx,time)

		
contains


! Central Finite Difference in x direction
function finite_dif_i(c_in,dx) result(xshift)
	implicit none
	
	real*8, intent(in), dimension(:,:) :: c_in
	real*8, intent(in) :: dx
	real*8, dimension(1:size(c_in,1),1:size(c_in,2)) :: xshift
	
	xshift = (cshift(c_in,1,1) - 2*c_in + cshift(c_in,-1,1))/dx**2

end function finite_dif_i


! Central Finite Difference in y direction
function finite_dif_j(c_in,dy) result(yshift)
	implicit none

	real*8, intent(in), dimension(:,:) :: c_in
	real*8, intent(in) :: dy
	real*8, dimension(1:size(c_in,1),1:size(c_in,2)) :: yshift

	yshift = (cshift(c_in,1,2) - 2*c_in + cshift(c_in,-1,2))/dy**2

end function finite_dif_j


! Functions
function dfdc(c_in,W) result(c_out)
	implicit none
	
	real*8, intent(in), dimension(:,:) :: c_in
	real*8, intent(in) :: W
	real*8, dimension(1:size(c_in,1),1:size(c_in,2)) :: c_out

	c_out = (W/2)*c_in*(1-c_in)*(1-2*c_in)

end function dfdc

function func(c_in,e,W,dx,dy) result(c_out)
	implicit none
	
	real*8, intent(in), dimension(:,:) :: c_in
	real*8, intent(in) :: e, W, dx, dy
	real*8, dimension(1:size(c_in,1),1:size(c_in,2)) :: c_out
	
	c_out = dfdc(c_in,W) - e**2*(finite_dif_i(c_in,dx)+finite_dif_j(c_in,dx))

end function func

end program Hilliard