function [X,Y,Z] = billiardgame3(N,massmax,vmax,T)
%generates and plays a game of billiards with N randomly placed billiard
%balls on a table.  The balls have random mass between 0 and massmax,
%velocity between 0 and vmax, and plays for T seconds, with timestep .01.
w = 9;
l = 4.5;
h = 4.5;

%M = massmax*rand(1,N);
%DY
M = ones(1,N);

vx = 2*vmax*rand(1,N)-vmax;
vy = 2*vmax*rand(1,N)-vmax;
vz = 2*vmax*rand(1,N)-vmax;

X0 = w*rand(1,N);
Y0 = l*rand(1,N);
Z0 = h*rand(1,N);

rad = 0.2;
tstep = 0.01;

[X,Y,Z] = billiards3(w,l,h,M,vx,vy,vz,X0,Y0,Z0,T,rad,tstep);

billiardplayer3(X,Y,Z,0.01,w,l,h)
end