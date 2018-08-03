clear all;
clc;

A = ReadArray_FortranBinary('DP_Smoothed.dat',2);
A = A';
pcolor(A);

[Nx Ny] = size(A);
dx = 0.5;
dy = 0.5;

dt = .1*dx^2;
t0 = 0;
tf = 5;
time = t0:dt:tf;

A = double(A);

tf = 50;
time = t0:dt:tf;

A = A + 1e-6;

Avg_Right = zeros(Nx,Ny);
Avg_Left = zeros(Nx,Ny);
Avg_Up = zeros(Nx,Ny);
Avg_Down = zeros(Nx,Ny);

Right = zeros(Nx,Ny);
Left = zeros(Nx,Ny);
Up = zeros(Nx,Ny);
Down = zeros(Nx,Ny);

Avg_Right(2:end-1,2:end-1) = (A(2:end-1,3:end)+A(2:end-1,2:end-1))/2;
Avg_Left(2:end-1,2:end-1) = (A(2:end-1,2:end-1)+A(2:end-1,1:end-2))/2;
Avg_Up(2:end-1,2:end-1) = (A(1:end-2,2:end-1)+A(2:end-1,2:end-1))/2;
Avg_Down(2:end-1,2:end-1) = (A(2:end-1,2:end-1)+A(3:end,2:end-1))/2;

dcdt_right = @ (c,dx) ((c(2:end-1,3:end)-c(2:end-1,2:end-1))/dx);
dcdt_left = @ (c,dx) ((c(2:end-1,2:end-1)-c(2:end-1,1:end-2))/dx);
dcdt_up = @ (c,dy) ((c(1:end-2,2:end-1)-c(2:end-1,2:end-1))/dy);
dcdt_down = @ (c,dy) ((c(2:end-1,2:end-1)-c(3:end,2:end-1))/dy);

for n = 1:length(time)
    c(1,:) = 0;
    c(Nx,:) = 0;
    c(:,1) = 1;
    c(:,Ny) = 1;
    
    Right(2:end-1,2:end-1) = dcdt_right(c,dx);
    Left(2:end-1,2:end-1) = dcdt_left(c,dx);
    Up(2:end-1,2:end-1) = dcdt_up(c,dx);
    Down(2:end-1,2:end-1) = dcdt_down(c,dx);

    c = c + (dt./A).*((Right.*Avg_Right-Left.*Avg_Left)/dx + ...
    (Up.*Avg_Up-Down.*Avg_Down)/dy);
    
%     c = c + (dt./A).*((dcdt_right(c,dx).*Avg_Right-dcdt_left(c,dx).*Avg_Left)/dx + ...
%     (dcdt_up(c,dy)*Avg_Up-dcdt_down(c,dy)*Avg_Down)/dy);    
    cstore(n,:,:) = c;    
end

for i = 1:50:length(time)
    rshape = reshape(cstore(i,:,:),[Nx Ny]);
    PsiC = rshape.*A;
    pcolor(PsiC)
    hold on;
    colormap jet;
    shading flat;
    axis([0 Ny 0 Nx]);
    hold off;
    pause(0.1);
end