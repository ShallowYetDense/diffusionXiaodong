function [X,Y,Z] = billiards3(w,l,h,M,vx,vy,vz,X0,Y0,Z0,T,rad,tstep)
%simulates billiards with elastic collisions on a w x l billiards table.  M
%should be a vector recording the (positive) masses of the billiard balls
%(the function will create as many balls as the length of M).  vx, vy, X0,
%Y0 will similarly be vectors giving the initial x and y velocities of each
%billiard ball, and then initial positions.  The program runs for T seconds
%with time step tstep.  A reasonable setup is
%billiards(9,4.5,randi(10,1,9),3*rand(1,9),3*rand(1,9),.5+randi(8,1,9),.4+ran
%di(4,1,9),5,.2,.01)

X = zeros(floor(T/tstep),length(M)); %initialize the three position arrays, one column per particle
Y = zeros(floor(T/tstep),length(M));
Z = zeros(floor(T/tstep),length(M));

X(1,:) = X0; %set initial position
Y(1,:) = Y0;
Z(1,:) = Z0;

for k = 2:floor(T/tstep)
tryxpos = X(k-1,:)+tstep*vx; %here we check if any collisions will happen in the next step
tryypos = Y(k-1,:)+tstep*vy;
tryzpos = Z(k-1,:)+tstep*vz;

[A,W] = collisions3(tryxpos,tryypos,tryzpos,rad,w,l,h);
for j = 1:size(A,1)
[vx(A(j,1)),vx(A(j,2))] = collide3(M(A(j,1)),M(A(j,2)),vx(A(j,1)),vx(A(j,2))); %avoiding collisions with particles
[vy(A(j,1)),vy(A(j,2))] = collide3(M(A(j,1)),M(A(j,2)),vy(A(j,1)),vy(A(j,2)));
[vz(A(j,1)),vz(A(j,2))] = collide3(M(A(j,1)),M(A(j,2)),vz(A(j,1)),vz(A(j,2)));
end
for j = 1:length(W)
if tryxpos(W(j)) >= w || tryxpos(W(j))<=0 %avoiding collisions with walls
vx(W(j)) = -vx(W(j));
elseif tryypos(W(j))>=l || tryypos(W(j))<=0
vy(W(j)) = -vy(W(j));
elseif tryzpos(W(j)) >=h || tryzpos(W(j))<=0
vz(W(j)) = -vz(W(j));
end
end
X(k,:) = X(k-1,:)+tstep*vx;%updating the position with the ?fixed? velocity vectors
Y(k,:) = Y(k-1,:)+tstep*vy;
Z(k,:) = Z(k-1,:)+tstep*vz;
end
end