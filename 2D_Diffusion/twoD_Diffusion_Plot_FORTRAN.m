clc
clear

A = ReadArray_FortranBinary('concdata.dat',3);

t0 = 0;
tf = 5;
dt = .0003;
time = t0:dt:tf;

Nx = 100;
Ny = 100;

for i = 1:100:length(time)
    rshape = reshape(A(i,:,:),[Nx Ny]);
    pcolor(rshape)
    hold on;
    colormap jet;
    shading flat;
    axis equal;
    hold off;
    pause(0.001);
end

