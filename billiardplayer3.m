function billiardplayer3(X,Y,Z,ptime,w,l,h)
%a program for visualizing billiard movement in 3 dimensions.

figure(1)
axis equal
%maxwindow
for j = 2:size(X,1)
%color the first particle red
plot3(X(j,1),Y(j,1),Z(j,1),'ro','MarkerFaceColor' , 'r')
hold on
%keep the rest of the particles blue circles
plot3(X(j,2:size(X,2)),Y(j,2:size(Y,2)),Z(j,2:size(Z,2)),'o')
hold off
axis([0 w 0 l 0 h])
grid on
pause(ptime)
end
hold on
%draw the path of the first particle
plot3(X(:,1),Y(:,1),Z(:,1),'k','LineSmoothing','on')
end