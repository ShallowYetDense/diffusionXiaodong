function [v1n,v2n] = collide3(m1,m2,v1,v2)
%solves the conservation of velocity and momentum for two balls of mass m1
%and m2 colliding at velocities v1 and v2, and returns their new
%velocities.
C1 = (m1-m2)/(m1+m2);
C2 = (2*m2)/(m1+m2);
C3 = (2*m1)/(m1+m2);
v1n = C1*v1+C2*v2;
v2n = -C1*v2+C3*v1;
end