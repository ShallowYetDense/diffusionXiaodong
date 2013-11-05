function [A,W] = collisions3(X,Y,Z,rad,w,l,h)
%takes vectors X, Y, and Z where each entry contains the x- or y-position of
%a billiard ball and returns an n x 2 matrix A containing the indices where
%there is a collision (assuming balls have radius rad),
%and a matrix W containing the indices of which balls collided with walls
%in a w x l x h rectangle
W=[];
C = X>=w; %oddly, I can?t find a way to make this code nicer. Checks collisions with walls
C = +C;
D = Y>=l;
D = +D;
E = X<=0;
E = +E;
F = Y<=0;
F = +F;
G = Z>=h;
G = +G;
H = Z<=0;
H = +H;
C = C+D+E+F+G+H;
for j = 1:length(C)
if C(j) > 0
W = [W;j];
end
end

X = ones(length(X),1)*X;
Y = ones(length(Y),1)*Y;
Z = ones(length(Z),1)*Z;

B = sqrt((X-X.').^2+(Y-Y.').^2+(Z-Z.').^2); %a matrix whose (i,j)th entry is the distance between particle i and particle j
B = B<rad;
A = [];
for j = 1:size(B,2)
for k = 1:j-1
if B(k,j)==1
A = [A;k,j];
end
end
end
end