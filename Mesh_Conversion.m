function[XYM,MAPM,neH_o]=Mesh_Conversion(XYH,MAPH)
% clear;
% clc;
% Sim_Fold='1015_S235JR';
% %%%% Folder Handlings %%%%%%%%%%%%%
% str1=pwd;                     % Get current working path
% str2='Programs';                % remove programs    
% Fold_Path=strrep(str1,str2,''); % Path of the main folder
% path_sep='\';
% PATH=[Fold_Path 'Simulations' path_sep Sim_Fold path_sep];
% PATH1=[PATH 'Process_Parameters.mat'];
% PATH2=[PATH 'Material_TM.mat'];
% PATH3=[PATH 'Material.mat'];
% % [XYH4,MAPH4] = heat_mesh_gen(XYe4,ner4,nes4);
% % [XYH,MAPH]   = mesh_combine(XYH,MAPH,XYH4,MAPH4,1.0e-6);
% PP=load(PATH1);
% XYH=PP.XYH;
% MAPH=PP.MAPH;


the=pi/2;
[XYH2]=Mesh_Transform(XYH,the);
neH_o=size(MAPH,1);
[XYH,MAPH]   = mesh_combine(XYH,MAPH,XYH2,MAPH,1.0e-6);
% draw_heat(XYH,MAPH)
[XYM,MAPM] = mech_mesh_gen(XYH,MAPH);

% ele_no=1; node_no=1; m='o';
% draw_mech_9noded_2D(XYM,MAPM,ele_no,node_no,m)

function NodeIds = node_on_pt(XY,Pt,TOL)
%  Finds all nodes on the line segment
%  Apply automatically the convetive bc.
%  Sub-function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  INPUT:                                                                    %
%    [XY]    -> Nodal coordinates                                            %
%    [StPt]  -> Line Start Point                                             %
%    [EndPt] -> Line End Point                                               %
%  OUTPUT:    Draws                                                          %
%    [NodeIds] -> Node id list                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nn = size(XY,1); % Number of elements

cnt = 0;
NodeIds = [];

for in=1:nn    % loop over nodes
  if norm(Pt-XY(in,:)) <= TOL
    cnt = cnt + 1;
    NodeIds(cnt,1) = in;
  end % if
end % for(in=1:nn)    % loop over nodes

function [XYM,MAPM] = mech_mesh_gen(XYH,MAPH)
%  Mesh Generation for the disc
%  Sub-function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  INPUT:                                                                    %
%    [XYH]  -> nodal coordinates of the mesh for heat transfer               %
%    [MAPH] -> connectivity of heat mesh                                     %
%  OUTPUT:                                                                   %
%    [XYH]  -> nodal coordinates of the mesh for mechanical calculations     %
%    [MAPH] -> connectivity of mechanical mesh                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TOL   = 1.0e-6;
ne    = size(MAPH,1);
XYnew = [];
MAPM  = [MAPH zeros(ne,5)];
nn    = size(XYH,1);
nid   = 0;
for ie=1:ne % loop over elements
  % Node 5
  Pt = 0.5*(XYH(MAPH(ie,1),:) + XYH(MAPH(ie,2),:));
  NodeId = node_on_pt(XYnew,Pt,TOL);
  if size(NodeId,1)==0
    nid=nid+1;
    XYnew = [XYnew; Pt];
    MAPM(ie,5) = nid + nn; 
  else 
    MAPM(ie,5) = NodeId + nn; 
  end

  % Node 6
  Pt = 0.5*(XYH(MAPH(ie,2),:) + XYH(MAPH(ie,3),:));
  NodeId = node_on_pt(XYnew,Pt,TOL);
  if size(NodeId,1)==0
    nid=nid+1;
    XYnew = [XYnew; Pt];
    MAPM(ie,6) = nid + nn; 
  else 
    MAPM(ie,6) = NodeId + nn; 
  end

  % Node 7
  Pt = 0.5*(XYH(MAPH(ie,3),:) + XYH(MAPH(ie,4),:));
  NodeId = node_on_pt(XYnew,Pt,TOL);
  if size(NodeId,1)==0
    nid=nid+1;
    XYnew = [XYnew; Pt];
    MAPM(ie,7) = nid + nn; 
  else 
    MAPM(ie,7) = NodeId + nn; 
  end

  % Node 8
  Pt = 0.5*(XYH(MAPH(ie,1),:) + XYH(MAPH(ie,4),:));
  NodeId = node_on_pt(XYnew,Pt,TOL);
  if size(NodeId,1)==0
    nid=nid+1;
    XYnew = [XYnew; Pt];
    MAPM(ie,8) = nid + nn; 
  else 
    MAPM(ie,8) = NodeId + nn; 
  end

  % Node 9
  Pt = 0.25*(XYH(MAPH(ie,1),:) + XYH(MAPH(ie,2),:) + ...
             XYH(MAPH(ie,3),:) + XYH(MAPH(ie,4),:));
  nid=nid+1;
  XYnew = [XYnew; Pt];
  MAPM(ie,9) = nid + nn; 
end
XYM = [XYH; XYnew];

function[XY2]=Mesh_Transform(XY,the)
Tran_Mat=[cos(the) sin(the);-sin(the) cos(the)];
XY2=XY;
for i=1:size(XY,1)
    xyi=XY(i,:);
    XY2(i,:)=(Tran_Mat*xyi')';
end

function draw_mech_9noded_2D(XYM,MAPM,ele_no,node_no,m)
%  Draw geometry
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  INPUT:
%    [XY]  : Coordinates of mechanical nodes
%    [MAP] : Mechanical element connectivity
%  OUTPUT:
%    DRAWS THE MECHANICAL MESH
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ne = size(MAPM,1); % Number of elements
figure; % create a new figure for geometry plot
map = [1 5 2 6 3 7 4 8 1]; % map for line data extraction
nodes=1:9;
tol=max(XYM(:,1))*5e-3;
yshift=max(XYM(:,2))/50;
for ie=1:ne
    % nodal points in draw order
    X  = XYM(MAPM(ie,map),1);
    Y  = XYM(MAPM(ie,map),2);
    Xc = XYM(MAPM(ie,9),1);
    Yc = XYM(MAPM(ie,9),2);
    line(X,Y,'Marker',m);   % boundary line
    line(Xc,Yc,'Marker',m); % center node
    if ele_no==1
        text(Xc,Yc+yshift,int2str(ie));
    end
    if node_no==1
        x_node =XYM(MAPM(ie,nodes),1);
        y_node =XYM(MAPM(ie,nodes),2);
        for j=1:length(nodes)
            node_number=MAPM(ie,j);
            x_coord=x_node(j);y_coord=y_node(j);
            text(x_coord+tol,y_coord,int2str(node_number));
        end
    end
end
xlabel('x-coordinate [m]');
ylabel('y-coordinate [m]');
% axis('equal');


function draw_heat(XY,MAP)
%  Draw geometry
%  Sub-function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  INPUT:                                                                    %
%    [XY]  -> Nodal coordinates                                              %
%    [MAP] -> Element connectivity                                           %
%  OUTPUT:    Draws                                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ne = size(MAP,1); % Number of elements

map = [1 2 3 4 1];

for ie=1:ne
  % nodal points in draw order
  X  = XY(MAP(ie,map),1);
  Y  = XY(MAP(ie,map),2);
  
  line(X,Y,'Marker','o');

  xm=0.25*(X(1)+X(2)+X(3)+X(4));
  ym=0.25*(Y(1)+Y(2)+Y(3)+Y(4));
  text(xm,ym,int2str(ie));
end


function [XY,MAP] = mesh_combine(XY1,MAP1,XY2,MAP2,TOL)
%  Mesh Generation for the system
%  Sub-function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  This function is called only once during inputdata to initialize some     %
%  phase transformation constants                                            %
%  INPUT:                                                                    %
%    [XY1]  -> Nodal coordinates of mesh-1                                   %
%    [MAP1] -> Element connectivity of mesh-1                                %
%    [XY2]  -> Nodal coordinates of mesh-2                                   %
%    [MAP2] -> Element connectivity of mesh-2                                %
%    [TOL]  -> merge tolerance for close nodes                               %
%  OUTPUT:                                                                   %
%    [XY]  -> Combined nodal coordinates                                     %
%    [MAP] -> Combined element connectivity                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ne1 = size(MAP1,1); % number of elements in mesh-1
ne2 = size(MAP2,1); % number of elements in mesh-2
nne = size(MAP1,2); % number of nodes per element
nn1 = size(XY1,1); % number of nodes in mesh-1
nn2 = size(XY2,1); % number of nodes in mesh-2

% find common nodes
cnc  = 0; % common node counter
for in2=1:nn2
  for in1=1:nn1
    if (abs(XY2(in2,1)-XY1(in1,1)) < TOL) && ...
       (abs(XY2(in2,2)-XY1(in1,2)) < TOL)
      cnc = cnc + 1;
      cnIds(cnc,:) = [in1 in2];
    end
  end
end

XYa = XY2; % make an copy for operations
XYa(cnIds(:,2),:) = []; % merge commen nodes
XY  = [XY1; XYa]; % combine

map = zeros(1,nn2);
map(cnIds(:,2)) = cnIds(:,1);
cnt = nn1;
for in2=1:nn2
  if map(in2)==0 
    cnt = cnt + 1;
    map(in2) = cnt;
  end
end

% Modify MAP2 and combine it MAP1 to form MAP
for ie2=1:ne2
  for ine=1:nne
    MAP2(ie2,ine) = map(MAP2(ie2,ine));
  end
end

MAP = [MAP1; MAP2];

