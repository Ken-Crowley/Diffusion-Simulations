clc
clear

% Reads in concentration dat file created in FORTRAN
A = ReadArray_FortranBinary('concdata.dat',3);

% Time parameters
t0 = 0;
tf = 5;
dt = .0003;
time = t0:dt:tf;

% Grid size
Nx = 100;
Ny = 100;

% Makes gif
gif('test.gif')

for i = 1:100:length(time)
    rshape = reshape(A(i,:,:),[Nx Ny]);
    pcolor(rshape)
    hold on;
    colormap jet;
    shading flat;
    axis equal;
    hold off;
    pause(0.001);
    gif
end

web('test.gif');

