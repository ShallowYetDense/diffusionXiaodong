
function diffusionXiaodong
global init_speed

  path(path,'./UniformSampling');
%parameters

n_parts = 200;     part_size = .1;     n_exits = 10;         exit_size = .5;         cell_size = 15;

nuc_size = 5;       high_viscosity = 1; low_viscosity = .2; viscosity_border = 12; %use volume ratio!

init_speed = 1;   time_limit = 10000; time_unit = 1;        wallBounceConst = low_viscosity;

rand_element_factor = 1;

viscosityVDelta = low_viscosity-high_viscosity;

draw3dSwitch = true;

%init positins
partsAll = rand_pick_sphere(n_parts,nuc_size,cell_size);%move out of wall contact areas?


%init velocity
vAll = rand_pick_sphere(n_parts,1,1).*init_speed;
%generate exit centers
exitCenters =  ParticleSampleSphere('N',n_exits).*cell_size;

theTime = 0;
partsExited = double.empty(0,3);
theFig = figure;
axis([-cell_size cell_size -cell_size cell_size -cell_size cell_size]);





%draw viscosity border
[sphX,sphY,sphZ] = sphere(100);



while theTime<time_limit
    
    if (draw3dSwitch)
        
        
        
    myplot3(exitCenters,'or');
    axis equal
    
    hold on
    myplot3(partsAll,'.g');
    myplot3v(partsAll,vAll);
    
    
    surfPtr=surf(sphX*viscosity_border,sphY*viscosity_border,sphZ*viscosity_border);
    set(surfPtr,'LineStyle','none','FaceColor',[.1,.1,.1],'FaceAlpha',.4);
    surfPtr2=surf(sphX*cell_size,sphY*cell_size,sphZ*cell_size);
    set(surfPtr2,'LineStyle','none','FaceColor',[.1,.1,.1],'FaceAlpha',.05)
    
     surfPtr2=surf(sphX*nuc_size,sphY*nuc_size,sphZ*nuc_size);
    set(surfPtr2,'LineStyle','none','FaceColor',[.1,.1,.1],'FaceAlpha',.8)
    
%    surfPtr

    myplot3(partsExited,'.r');

    end
    
    

partsTmp = partsAll;%use tmp only in exqusions below!
vTmp = vAll;

%get sph coords
[ths,phs,rs] = cart2sph(partsTmp(:,1),partsTmp(:,2),partsTmp(:,3));

projections = cell_size*(partsTmp./repmat(sum(partsTmp.^2,2).^.5,[1,3]));
%myplot3(projections,'og');
%remove exiting particles and vs
exitBool = logical(zeros(n_parts,1));

for ectr = 1:n_exits%for each exit mark inliers
    c = exitCenters(ectr,:);%exit coords
    cdists = ((projections(:,1)-c(1)).^2+(projections(:,2)-c(2)).^2+(projections(:,3)-c(3)).^2).^.5;
     Condition1 = ((cdists<exit_size));
     
     exitBool = exitBool | (Condition1);
end

exitBool = exitBool & (rs>viscosity_border);

partsExited = [partsExited;partsTmp(exitBool,:)];
%remove
partsTmp(exitBool,:)=[];
vTmp(exitBool,:)=[];
projections(exitBool,:) = [];

n_parts = size(partsTmp,1);
if n_parts ==0
    break
end

%myplot3dist(partsTmp,projections);

[ths,phs,rs] = cart2sph(partsTmp(:,1),partsTmp(:,2),partsTmp(:,3));
[vThs,vPhs,vRs] = cart2sph(vTmp(:,1),vTmp(:,2),vTmp(:,3));


[~,~,rs2] = cart2sph(partsTmp(:,1)+vTmp(:,1),partsTmp(:,2)+vTmp(:,2),partsTmp(:,3)+vTmp(:,3));
% check for collision in the NEXT STEP


dVTmp = zeros(size(vTmp));

%check outer wall collision
wallCollBool = (repmat((rs2+part_size)>cell_size,[1,3]));
wallCollBool = wallCollBool|(repmat((rs2-part_size)<nuc_size,[1,3]));

poswallVector = projections-partsTmp;%which direction?
poswallVectorNormd = poswallVector./repmat(sum(poswallVector.^2,2).^.5,[1,3]);
poswallVConst = wallBounceConst * poswallVectorNormd;

poswallVConst = poswallVConst .* wallCollBool;%only apply to parts close to the wall
 

%update v
%????????????

 vTmp = poswallVConst +vTmp;
 vTmp = vTmp./repmat(sum(vTmp.^2,2).^.5,[1,3]);%normalized

vTmp(isnan(vTmp)) = 0;
%give low velocity in high viscosity zone
%viscosityBool = repmat((rs)>viscosity_border,[1,3]);
% viscosityBool = rs>viscosity_border;
% vTmp(viscosityBool,:) = vTmp(viscosityBool,:).*low_viscosity;
% 


%?????????
%dVTmp = dVTmp+poswallVConst;



%check nucleus collision

%[posWallVectorThs,posWallVectorPhs,~]=cart2sph(poswallVector(:,1),poswallVector(:,2),poswallVector(:,3));

% vThs(wallCollBool)
% vPhs(wallCollBool)
% 
% wallCollAngle = 



%check part collisions
if true
for partctr = 1:n_parts
    c = partsTmp(partctr,:);
    e = [partctr+1:n_parts];%exlude self..AND PARTICLES ALREADY CHECKED 
    %distances to all other particles
    %cdists = ((partsTmp(e,1)-c(1)).^2+(partsTmp(e,2)-c(2)).^2+(partsTmp(e,3)-c(3)).^2).^.5;
    
    cdists = ((partsTmp(e,1)-c(1)).^2+(partsTmp(e,2)-c(2)).^2+(partsTmp(e,3)-c(3)).^2).^.5;
    
    distCondition = cdists<2.*part_size;
    collParts = partsTmp(distCondition,:);
 %disp(e);
 %disp(collParts);
 
   for ctr = e(distCondition(1:end)) %each collision with this particle
       
       
       
      [collVthis(1) collVother(1)]= mycollide(vTmp(partctr,1),vTmp(ctr,1));
      [collVthis(2) collVother(2)]= mycollide(vTmp(partctr,2),vTmp(ctr,2));
      [collVthis(3) collVother(3)]= mycollide(vTmp(partctr,3),vTmp(ctr,3));
       
       dVTmp(ctr,:)=dVTmp(ctr,:)+collVthis;
       dVTmp(partctr,:)=dVTmp(partctr,:)+collVother;
       
%        
%        dVTmp(ctr,:)= dVTmp(ctr,:) + collVthis;
%        
%        dVTmp(partctr,:)
%         
   end
        
    
end

end




%normalize dVTmp
  dVTmp = dVTmp./repmat(sum(dVTmp.^2,2).^.5,[1,3]);%normalized
  dVTmp(isnan(dVTmp)) = 0;
  vTmp=dVTmp + (dVTmp == 0) .* vTmp;


  
 randElement = 2*rand(size(vTmp))-1;
 randElement = randElement./repmat(sum(randElement.^2,2).^.5,[1,3]);
  
 
 vTmp = vTmp*(1-rand_element_factor) + randElement*rand_element_factor;

% 
%  vTmp = dVTmp +vTmp;
%  vTmp = vTmp./repmat(sum(vTmp.^2,2).^.5,[1,3]);%normalized
% 
% vTmp(isnan(vTmp)) = 0;
% %viscosityBool = repmat((rs)>viscosity_border,[1,3]);


% %give low velocity in high viscosity zone
 viscosityBool = rs>viscosity_border;
 vTmp(viscosityBool,:) = vTmp(viscosityBool,:).*low_viscosity;
% %vTmp= vTmp + viscosityVDelta*(repmat((rs)>viscosity_border,[1,3]));



vAll = vTmp;
partsTmp = vTmp+partsTmp;
partsAll = partsTmp;


theTime = theTime+time_unit;

if(draw3dSwitch)
%myplot3(partsExited,'.r');
drawnow
hold off
else
    
    hist(rs);
    drawnow;

end


disp(max(rs));

end


function [v1n,v2n]= mycollide(v1,v2)
%thanks to Colin Carroll

C1 = 0;
C2 = 1;
C3 = 1;
v1n = C1*v1+C2*v2;
v2n = -C1*v2+C3*v1;




function myplot3dist(in,in2)

for vctr = 1:size(in,1)
plot3([in(vctr,1);in2(vctr,1)],[in(vctr,2);in2(vctr,2)],[in(vctr,3);in2(vctr,3)],'b');
end


function myplot3v(in,in2)%hold must be on

d = in+in2;
myplot3dist(in,d);



function myplot3(in,in2)
plot3(in(:,1),in(:,2),in(:,3),in2);


function spheres_intersect

%posFromCenter = 





%SCRATCH
   
    %collision norms (the other is negative)
    %collNorm = (collParts - repmat(c,[1,size(collParts,1)]));%to multiply vs befor summing
    %collNorm = 
    
    
    %collVs = 
