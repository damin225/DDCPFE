clear;clc
% This program is to generate multiple point constraints (MPCs) in ABAQUS
% format. Currently only periodic mesh is supported

% Following procedures are taken:
% (1) Read in nodes and element connectivity table
% (2) Find nodes in the pair faces
% (3) Identify master-slave node pairs
% (4) Write MPC to files

% define some parameters beforehand
edgeprox = 1E-5;

%=============================
% Read nodal information
%=============================

% Open the original input file
file =  './Direct_short_crack.inp';
fid = fopen(file);

% Read each line from the input file until the *Node line is reached
tline = fgetl(fid); 
while strcmp(strtrim(tline),'*Node') == 0
    tline = fgetl(fid);
end

% Write the nodes file for each of the nodes in the structure
nodesfile = './PeriodicNodes.txt';
fidw = fopen(nodesfile,'w');

tic
% Write to PeriodicNodes.txt file first
while strncmp(tline,'*Element',8) == 0
    tline = fgetl(fid);
    if strncmp(tline,'*Element',8) == 0
        fprintf(fidw,'%s \n',tline);
    end
end

% load node information into array
load('./PeriodicNodes.txt')

%============================
% Read element information
%============================

% Write the element connectivity table file for each of the element in the 
% structure
eleconFile =  './elecon.txt';
fidw2 = fopen(eleconFile,'w');

tic
% Write to elecon.txt file first
while (strncmp(tline,'*Elset',6)==0 && strncmp(tline,'*Nset',5)== 0)
    tline = fgetl(fid);
    elenum=str2num(strtok(tline, ','));
    if (strncmp(tline,'*Elset',6)==0 && strncmp(tline,'*Nset',5)== 0)
        fprintf(fidw2,'%s \n',tline);
    end
end

% Load element connectivity table into array
load('./elecon.txt')

fprintf('Element connectivity talbe Readtime: ')
toc

%===========================================
% Identify the nodes in the pair faces
%===========================================

% Store node numbering and coordinate seperately 
Node=PeriodicNodes(:,1); X=PeriodicNodes(:,2);Y=PeriodicNodes(:,3);Z=PeriodicNodes(:,4); 
loc = [X,Y,Z];

% Define limits of the model in X,Y,Z
xmin = min(X);   xmax = max(X);   xdist=xmax-xmin;
ymin = min(Y);   ymax = max(Y);   ydist=ymax-ymin;
zmin = min(Z);   zmax = max(Z);   zdist=zmax-zmin;


% Find the middle node to set to the LockedNode
xmid = (xmin+xmax)/2;
zmid = (zmin+zmax)/2;
ymid = (ymin+ymax)/2;
Lockednode = knnsearch([X,Y,Z],[xmid,ymid,zmid]);

% Initialize the count for the number of nodes in each set, these node sets
% will include corner nodes and will be updated with nodetop ... to
% exclude corner nodes.
ntop = 1; nbot = 1; nright = 1; nleft = 1; nfront = 1; nback = 1;
top = []; bot = []; right = []; left = []; front = []; back = [];

%initialize coordinates bound of face nodes 
topbound=[zdist+1, -zdist-1, xdist+1, -xdist-1];% minium z, max z, miniumx, max x
bottombound=[zdist+1, -zdist-1, xdist+1, -xdist-1];

rightbound=[ydist+1, -ydist-1, zdist+1, -zdist-1];% minium y, max y, miniumz, max z
leftbound=[ydist+1, -ydist-1, zdist+1, -zdist-1];

frontbound=[xdist+1, -xdist-1, ydist+1, -ydist-1];% minium x, max x, minium y, max y
backbound=[xdist+1, -xdist-1, ydist+1, -ydist-1];

tic
% Identify the Nodes in the top,bottom,right,and left of the model
for i = 1:length(Node)
  % find nodes on Top surface
    if abs(loc(i,2)-ymax) <= edgeprox
        top(ntop) = Node(i);
        ntop = ntop +1;
  % update the min,max z and min, max x coordinates on top surface
        if loc(i,3)<= topbound(1)
            topbound(1)=loc(i,3);
        end
        if loc(i,3)>= topbound(2)
            topbound(2)=loc(i,3);
        end       
        if loc(i,1)<= topbound(3)
            topbound(3)=loc(i,1);
        end
      
        if loc(i,1)>= topbound(4)
            topbound(4)=loc(i,1);
        end
      
   % find nodes on Bottom surface
     elseif abs(loc(i,2)-ymin) <= edgeprox
        bot(nbot) = Node(i);
        nbot = nbot +1;
   % update the min and max z and x coordinates on bottom surface
        if loc(i,3)<= bottombound(1)
            bottombound(1)=loc(i,3);
        end
        if loc(i,3)>= bottombound(2)
            bottombound(2)=loc(i,3);
        end       
        if loc(i,1)<= bottombound(3)
            bottombound(3)=loc(i,1);
        end
      
        if loc(i,1)>= bottombound(4)
            bottombound(4)=loc(i,1);
        end
    end
    
   % find nodes on Left surface
     if abs(loc(i,1)-xmin) <= edgeprox
        left(nleft) = Node(i);
        nleft = nleft +1;
   % update the min and max y and z coordinates on left surface
        if loc(i,2)<=leftbound(1)
           leftbound(1)= loc(i,2);       
        end
        
        if loc(i,2)>=leftbound(2)
           leftbound(2)= loc(i,2);       
        end 
        
        if loc(i,3)<=leftbound(3)
           leftbound(3)= loc(i,3);       
        end
        
        if loc(i,3)>=leftbound(4)
           leftbound(4)= loc(i,3);       
        end  
        
   % find nodes on Right surface
     elseif abs(loc(i,1)-xmax) <= edgeprox
        right(nright) = Node(i);
        nright = nright +1;
   % update the min and max y and z coordinates on right surface
        if loc(i,2)<=rightbound(1)
           rightbound(1)= loc(i,2);       
        end
        
        if loc(i,2)>=rightbound(2)
           rightbound(2)= loc(i,2);       
        end 
        
        if loc(i,3)<=rightbound(3)
           rightbound(3)= loc(i,3);      
        end
        
        if loc(i,3)>=rightbound(4)
           rightbound(4)= loc(i,3);       
        end         
        
     end
    
   % find nodes on Front surface
     if abs(loc(i,3)-zmax) <= edgeprox
        front(nfront) = Node(i);
        nfront = nfront +1;
   % update min x, maxx, min y, max y on front surface
        if loc(i,1)<=frontbound(1)
            frontbound(1)=loc(i,1);
        end
        if loc(i,1)>=frontbound(2)
            frontbound(2)=loc(i,1);
        end
        if loc(i,2)<=frontbound(3)
            frontbound(3)=loc(i,2);
        end
        if loc(i,2)>=frontbound(4)
            frontbound(4)=loc(i,2);
        end
 
   % find nodes on back surface
     elseif abs(loc(i,3)-zmin) <= edgeprox
        back(nback) = Node(i);
        nback = nback +1;
    % update min x, maxx, min y, max y on back surface
        if loc(i,1)<=backbound(1)
            backbound(1)=loc(i,1);
        end
        if loc(i,1)>=backbound(2)
            backbound(2)=loc(i,1);
        end
        if loc(i,2)<=backbound(3)
            backbound(3)=loc(i,2);
        end
        if loc(i,2)>=backbound(4)
            backbound(4)=loc(i,2);
        end 
     end
end
fprintf('Time to identify facenodes: ');
toc

% Find the corner nodes
tic
% Identify corner nodes
CornerCorrd = [xmin,ymin,zmin; xmin,ymin,zmax; xmin,ymax,zmin;...
               xmin,ymax,zmax; xmax,ymin,zmin; xmax,ymin,zmax;...
               xmax,ymax,zmin; xmax,ymax,zmax];
CORNERS = transpose(knnsearch([X,Y,Z],CornerCorrd));

% Top and topnocorner are used to store the nodes of a face with or without 
% corner nodes
topwithcorner=top;
botwithcorner=bot;
rightwithcorner=right;
leftwithcorner=left;
frontwithcorner=front;
backwithcorner=back;

% Exclude the corner nodes from top
for j=1:size(CORNERS,2)
% exclude the corner nodes from top    
    flag=double(ismember(top, CORNERS(j)));
    if  (norm(flag, 1))    
        top(ind2sub(size(top,2),find(flag)))=[];
    	ntop=ntop-1;
    end
    
% exclude the corner nodes from bottom    
    flag=double(ismember(bot, CORNERS(j)));
    if  (norm(flag, 1))    
        bot(ind2sub(size(bot,2),find(flag)))=[];
    	nbot=nbot-1;
    end
    
% exclude the corner nodes from right    
    flag=double(ismember(right, CORNERS(j)));
    if  (norm(flag, 1))    
        right(ind2sub(size(right,2),find(flag)))=[];
    	nright=nright-1;
    end
    
% exclude the corner nodes from left    
    flag=double(ismember(left, CORNERS(j)));
    if  (norm(flag, 1))    
        left(ind2sub(size(left,2),find(flag)))=[];
    	nleft=nleft-1;
    end
    
% exclude the corner nodes from front    
    flag=double(ismember(front, CORNERS(j)));
    if  (norm(flag, 1))    
        front(ind2sub(size(front,2),find(flag)))=[];
    	nfront=nfront-1;
    end
    
% exclude the corner nodes from back    
    flag=double(ismember(back, CORNERS(j)));
    if  (norm(flag, 1))    
        back(ind2sub(size(back,2),find(flag)))=[];
    	nback=nback-1;
    end
    
end
fprintf('Time to identify corner nodes: ');
toc

% find edge nodes

% Convention for each edge:
% (1) two alphebet indicating two axis the edge is perpenticular to;                        
% (2) lower or upper case indicating has smallest or largerst value along that axis                        
% (3) two alphebet always follow the x-y-z order                           
% Example: xzdege indicates its one of the four deges perpenticular to X and Z axis, 
%          while nodes on this dege has smallest x and z value
% (4) edge nodes are excluded from corner nodes
% (5) master edges  are XY, XZ, and YZ. Specific arrays are defined for
%     master edges to include corner nodes

tic
xyedge= find(ismember([X,Y], [xmin,ymin], 'rows'));
commonconponents=ismember(xyedge, CORNERS);
xyedge(commonconponents)=[];
[edgenode, edgenodeorder]=sort(Z(xyedge));
xyedge=xyedge(edgenodeorder);
%
xYedge= find(ismember([X,Y],[xmin,ymax], 'rows'));
commonconponents=ismember(xYedge, CORNERS);
xYedge(commonconponents)=[];
[edgenode, edgenodeorder]=sort(Z(xYedge));
xYedge=xYedge(edgenodeorder);
%
Xyedge= find(ismember([X,Y],[xmax,ymin], 'rows'));
commonconponents=ismember(Xyedge, CORNERS);
Xyedge(commonconponents)=[];
[edgenode, edgenodeorder]=sort(Z(Xyedge));
Xyedge=Xyedge(edgenodeorder);
%
XYedge= find(ismember([X,Y],[xmax,ymax], 'rows'));
XYedgeWithCorner = XYedge;
commonconponents=ismember(XYedge, CORNERS);
XYedge(commonconponents)=[];
[edgenode, edgenodeorder]=sort(Z(XYedge));
XYedge=XYedge(edgenodeorder);
%
xzedge= find(ismember([X,Z],[xmin,zmin], 'rows'));
commonconponents=ismember(xzedge, CORNERS);
xzedge(commonconponents)=[];
[edgenode, edgenodeorder]=sort(Y(xzedge));
xzedge=xzedge(edgenodeorder);
%
xZedge= find(ismember([X,Z],[xmin,zmax], 'rows'));
commonconponents=ismember(xZedge, CORNERS);
xZedge(commonconponents)=[];
[edgenode, edgenodeorder]=sort(Y(xZedge));
xZedge=xZedge(edgenodeorder);
%
Xzedge= find(ismember([X,Z],[xmax,zmin], 'rows'));
commonconponents=ismember(Xzedge, CORNERS);
Xzedge(commonconponents)=[];
[edgenode, edgenodeorder]=sort(Y(Xzedge));
Xzedge=Xzedge(edgenodeorder);
%
XZedge= find(ismember([X,Z],[xmax,zmax], 'rows'));
XZedgeWithCorner = XZedge;
commonconponents=ismember(XZedge, CORNERS);
XZedge(commonconponents)=[];
[edgenode, edgenodeorder]=sort(Y(XZedge));
XZedge=XZedge(edgenodeorder);
%
yzedge= find(ismember([Y,Z],[ymin,zmin], 'rows'));
commonconponents=ismember(yzedge, CORNERS);
yzedge(commonconponents)=[];
[edgenode, edgenodeorder]=sort(X(yzedge));
yzedge=yzedge(edgenodeorder);
%
yZedge= find(ismember([Y,Z],[ymin,zmax], 'rows'));
commonconponents=ismember(yZedge, CORNERS);
yZedge(commonconponents)=[];
[edgenode, edgenodeorder]=sort(X(yZedge));
yZedge=yZedge(edgenodeorder);
%
Yzedge= find(ismember([Y,Z],[ymax,zmin], 'rows'));
commonconponents=ismember(Yzedge, CORNERS);
Yzedge(commonconponents)=[];
[edgenode, edgenodeorder]=sort(X(Yzedge));
Yzedge=Yzedge(edgenodeorder);
%
YZedge= find(ismember([Y,Z],[ymax,zmax], 'rows'));
YZedgeWithCorner = YZedge;
commonconponents=ismember(YZedge, CORNERS);
YZedge(commonconponents)=[];
[edgenode, edgenodeorder]=sort(X(YZedge));
YZedge=YZedge(edgenodeorder);
%
fprintf('Time to identify edges: ');
toc

%==========================
% Pair master-slave faces
%==========================
% Currently only support pairing left and right faces, other faces can be
% easily added if needed

% pair left-right faces
pairchoice = 0;
pairchoice = input('Do you want to pair left-right face? 1:Yes 0:No (default=No): ');
if (pairchoice ==1) 
    tic 
    % Right-left pair face
    RLPair=PairMasterSlave(rightwithcorner, left, PeriodicNodes, elecon, 1);
    fprintf('Time to Pair right and left face: ');
    toc
end

% modify RLPair
for i=1:size(RLPair,1)
    for j=6:9
        if (RLPair(i,j)==1)
            RL(i,1) = RLPair(i,1);
            RL(i,2) = RLPair(i,j-4);
        end
    end
end
if (size(RL,1) ~= size(RLPair,1))
    disp('please use periodic mesh')
    return
end
            
%================================
% write to file in ABAQUS format
%================================
MPCfile = fopen('./MPC.txt','w');
format1 = '*Nset, nset=master-%d, instance=PART-1-1\n';
for i=1:size(RL,1)
    fprintf(MPCfile,format1,i);
    fprintf(MPCfile,'%d,\n',RL(i,2));
end
format2 = '*Nset, nset=slave-%d, instance=PART-1-1\n';
for i=1:size(RL,1)
    fprintf(MPCfile,format2,i);
    fprintf(MPCfile,'%d,\n',RL(i,1));
end
format3 = '*Surface, type=NODE, name=salve-%d_CNS_, internal\n slave-%d,';
for i=1:size(RL,1)
    fprintf(MPCfile,format3,i,i);
    fprintf(MPCfile,'%f\n',1.);
end
format4 = '** Constraint: Constraint-%d\n*MPC\n TIE, slave-%d, master-%d\n';
for i=1:size(RL,1)
    fprintf(MPCfile,format4,i,i,i);
end





                    %=============================
                    % functions are defined below
                    %=============================

%============================================================================
function Pair=PairMasterSlave(masterface, slaveface, nodecoor, elecon, faceid)
%============================================================================

% Originally developed by Xiang Zhang to handle Tet elments
% Now modified by Damin to handle Hex elements

%-------------------------------------------
% Outputs:
% Pair: master-slave pair with coefficients
%-------------------------------------------

%---------------------------------------------------
% Inputs:
% masterface: master face nodes
% slaveface:  slave face nodes
% nodecoor:   node coordinate table
% elecon:     element connectivity table
% faceid:     1-->LR pair; 2-->BT pair; 3-->BF pair
%----------------------------------------------------

% initialize the pair and number of pairs
Pair=[];
numpair=0;

%  initialize the master elements and number of master elelments
masterele=[];
nummasterele=0;
for i=1:size(elecon, 1)
    % find the master elements
    [a1, b1, c1]=intersect(elecon(i, 2:9), masterface);
    if (size(a1, 2)==4)
        nummasterele=nummasterele+1;
        masterele(nummasterele, 1:9)=elecon(i,:);
        masterele(nummasterele, 10:17)=0;
        % flag to indicate which 4 nodes are used in a master element
        masterele(nummasterele, b1(1)+9)=1;
        masterele(nummasterele, b1(2)+9)=1;
        masterele(nummasterele, b1(3)+9)=1;
        masterele(nummasterele, b1(4)+9)=1;
    end
end




% loop over each slave node
for i=1:size(slaveface, 2)
    flag=0;
    % loop over all the face elements
    for j=1:nummasterele
        switch faceid
            case 1
                [inelement,coef]=checkinelement(nodecoor(masterface(1),2), nodecoor(slaveface(i), 3), ...
                    nodecoor(slaveface(i), 4), masterele(j,:), nodecoor, 1);
                if inelement
                    numpair=numpair+1;
                    Pair(numpair,1)=slaveface(i);
                    Pair(numpair,2:9)=coef;
                    flag=1;
                    break
                end
                
            case 2
                [inelement,coef]=checkinelement(nodecoor(slaveface(i),2),nodecoor(masterface(1),3),  ...
                    nodecoor(slaveface(i),4), masterele(j,:), nodecoor, 2);
                if inelement
                    numpair=numpair+1;
                    Pair(numpair,1)=slaveface(i);
                    Pair(numpair,2:9)=coef;
                    flag=1;
                    break
                end
                
            case 3
                [inelement,coef]=checkinelement(nodecoor(slaveface(i), 2), nodecoor(slaveface(i),3),  ...
                    nodecoor(masterface(1),4), masterele(j,:), nodecoor, 3);
                if inelement
                    numpair=numpair+1;
                    Pair(numpair,1)=slaveface(i);
                    Pair(numpair,2:9)=coef;
                    flag=1;
                    break
                end
        end
    end
end

end



%========================================================================
function [inelement, coef]=checkinelement(x, y, z, ele, nodecoor, faceid)
%========================================================================

% Originally developed by Xiang Zhang to handle Tet elments
% Now modified by Damin to handle Hex elements


%-------------------------------------------------------
% Outputs:
% inelement: 0--> not in element; 1--> in element
% coef:      coeffcients to all nodes within an element
%-------------------------------------------------------


%-------------------------------------------------------
%Inputs:
% x,y,z:    coordinate of the slave node to be checked
% ele:      master element to be checked
% nodecoor: node coordinate table including node id
% elecon:   element connectivity table
% faceid:   1-->LR pair; 2-->BT pair; 3-->BF pair
%-------------------------------------------------------


% Initilization and parameter setting
tolerance=1e-5; % tolearance to identify real slave node.
inelement=0;    % flag to indicate in element or not: 0=not
coef=zeros(1,8);
eff=zeros(1,8);
flag=find(ele(10:17)==1);

% Get master surface element. Node numbering is in CCW order
[surfnode, Idx] = getsurfele(ele,nodecoor,faceid);

% Construct vectors
pointp=[x,y,z];
v1=surfnode(1,2:4)-pointp;          % PA
u1=surfnode(1,2:4)-surfnode(4,2:4); % DA
w1=surfnode(2,2:4)-surfnode(1,2:4); % AB


v2=surfnode(2,2:4)-pointp;          % PB
u2=surfnode(2,2:4)-surfnode(1,2:4); % AB
w2=surfnode(3,2:4)-surfnode(2,2:4); % BC

v3=surfnode(3,2:4)-pointp;          % PC
u3=surfnode(3,2:4)-surfnode(2,2:4); % BC
w3=surfnode(4,2:4)-surfnode(3,2:4); % CD

v4=surfnode(4,2:4)-pointp;          % PD
u4=surfnode(4,2:4)-surfnode(3,2:4); % CD
w4=surfnode(1,2:4)-surfnode(4,2:4); % DA

sign1=sign(cross(v1,w1)*cross(u1,w1)');
sign2=sign(cross(v2,w2)*cross(u2,w2)');
sign3=sign(cross(v3,w3)*cross(u3,w3)');
sign4=sign(cross(v4,w4)*cross(u4,w4)');

% The case when P are in the same side of each side of the polygon
if(sign1==1&&sign2==1&&sign3==1&&sign4==1)
    inelement=1;
    eff=geteff(nodecoor(ele(2), 2:4), nodecoor(ele(3), 2:4), nodecoor(ele(4), 2:4), nodecoor(ele(5), 2:4), ...
               nodecoor(ele(6), 2:4), nodecoor(ele(7), 2:4), nodecoor(ele(8), 2:4), nodecoor(ele(9), 2:4), pointp);
    coef=[ele(flag(Idx(1))+1), ele(flag(Idx(2))+1), ele(flag(Idx(3))+1), ele(flag(Idx(4))+1), ...
          eff(flag(Idx(1))),eff(flag(Idx(2))),eff(flag(Idx(3))),eff(flag(Idx(4)))];
    return
end

% The case when P are on one side of the polygon
if ((sign1==0 && sign2==1&&sign3==1&&sign4==1)||(sign2==0 && sign1==1&&sign3==1&&sign4==1)||...
    (sign3==0 && sign1==1&&sign2==1&&sign4==1)||(sign4==0 && sign1==1&&sign2==1&&sign3==1))
    inelement=1;
    eff=geteff(nodecoor(ele(2), 2:4), nodecoor(ele(3), 2:4), nodecoor(ele(4), 2:4), nodecoor(ele(5), 2:4), ...
               nodecoor(ele(6), 2:4), nodecoor(ele(7), 2:4), nodecoor(ele(8), 2:4), nodecoor(ele(9), 2:4), pointp);
    coef=[ele(flag(Idx(1))+1), ele(flag(Idx(2))+1), ele(flag(Idx(3))+1), ele(flag(Idx(4))+1), ...
          eff(flag(Idx(1))),eff(flag(Idx(2))),eff(flag(Idx(3))),eff(flag(Idx(4)))];
    return
end

% The case when dummmy node coincide with one master node
if norm(v1, 1)<=tolerance
    inelement=1;
    eff(flag(Idx(1)))=1.0;
    coef=[ele(flag(Idx(1))+1), ele(flag(Idx(2))+1), ele(flag(Idx(3))+1), ele(flag(Idx(4))+1), ...
          eff(flag(Idx(1))),eff(flag(Idx(2))),eff(flag(Idx(3))),eff(flag(Idx(4)))];
    return
end

if norm(v2, 1)<=tolerance
    inelement=1;
    eff(flag(Idx(2)))=1.0;
    coef=[ele(flag(Idx(1))+1), ele(flag(Idx(2))+1), ele(flag(Idx(3))+1), ele(flag(Idx(4))+1), ...
          eff(flag(Idx(1))),eff(flag(Idx(2))),eff(flag(Idx(3))),eff(flag(Idx(4)))];
    return
end

if norm(v3, 1)<=tolerance
    inelement=1;
    eff(flag(Idx(3)))=1.0;
    coef=[ele(flag(Idx(1))+1), ele(flag(Idx(2))+1), ele(flag(Idx(3))+1), ele(flag(Idx(4))+1), ...
          eff(flag(Idx(1))),eff(flag(Idx(2))),eff(flag(Idx(3))),eff(flag(Idx(4)))];
    return
end

if norm(v4, 1)<=tolerance
    inelement=1;
    eff(flag(Idx(4)))=1.0;
    coef=[ele(flag(Idx(1))+1), ele(flag(Idx(2))+1), ele(flag(Idx(3))+1), ele(flag(Idx(4))+1), ...
          eff(flag(Idx(1))),eff(flag(Idx(2))),eff(flag(Idx(3))),eff(flag(Idx(4)))];
    return
end
end



%=================================================================
function [surfnode, Idx] = getsurfele(masterele, nodecoor, faceid)
%=================================================================

% function to get master surface element (4 nodes) and sorted them in CCW
% order. Node that the element is not indexed. 

% copyright: Damin Xia

%--------------------------------------------------------------------
% Outputs:
% surfnode: the 2-D element corresponding to the 3-D master element
% Idx     : permutation index s.t. flag(Idx) = surfnode
%--------------------------------------------------------------------


%-------------------------------------------------------
% Inputs:
% masterele: master element (in 3-D)
% nodecoor:  node coordinate table
% faceid:    1-->LR pair; 2-->BT pair; 3-->BF pair
%-------------------------------------------------------


% Find the centroid of the polygon - just the arithmetic mean of the points
node =  masterele(2:9);
coef = masterele(10:17);
surfnode = nodecoor(node(find(coef==1)),:); %(nodeid x y z)
center = mean(surfnode(:,2:4),1);           %(x,y,z)

% Calculate the angle from the centroid to each point using atan2 
switch faceid
    case 1
        for i=1:4
            angle(i) = atan2(surfnode(i,3)-center(2),surfnode(i,4)-center(3));
        end
    case 2
        for i=1:4
            angle(i) = atan2(surfnode(i,4)-center(3),surfnode(i,2)-center(1));
        end
    case 3
        for i=1:4
            angle(i) = atan2(surfnode(i,2)-center(1),surfnode(i,3)-center(2));
        end
    otherwise
        disp('wrong face id')
end

% Sort in descending angle to establish a CCW node ordering
[sorted_angle,Idx] = sort(angle,'descend');
surfnode = surfnode(Idx,:);
end




