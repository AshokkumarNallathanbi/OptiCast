% function[]=solve_mech_2()
function[]=solve_mech_2(BEAMBC,MAT_MODEL,SimF)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%SimF
% clear;clc;
global MAT_MODEL BEAMBC 
% tic
% BEAMBC
% MAT_MODEL
% Sim_Fold
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAT_MODEL=1;
% BEAMBC=zeros(3,2);
% BEAMBC=[1 0; 1 0.126; 1 0]; % Plane strain
% Sim_Fold='1015_S235JR_2';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%% Folder Handlings %%%%%%%%%%%%%
% str1=pwd;                     % Get current working path
% str2='Programs';                % remove programs    
% Fold_Path=strrep(str1,str2,''); % Path of the main folder
% path_sep='\';
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % PATH=[Fold_Path 'Simulations' path_sep Sim_Fold path_sep 'Process_Parameters.mat'];
% PATH_sim=[Fold_Path 'Simulations' path_sep Sim_Fold path_sep];
% PATH_Mech_program=[str1 '\MECH'];

addpath(SimF.Mech_Prog);

% PATH_Material_Mech=[PATH 'Material_TM.mat'];
PATH_Material=[SimF.Sim_Heat 'Material.mat'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load(PATH_Material);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FN='MECH3';
% if isdir([PATH_sim FN])==0
%     mkdir([PATH_sim FN]);
% end
% PATH_Mech=[PATH_sim FN '\'];
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% 
% MECH_MAT=MECH_Material_Prop2(TS,TL);
% Make_Mech_Material_File(MECH_MAT,SimF.Sim_Mech);

mech_solver(SimF.Sim_Mech,SimF.Sim_Heat)

% PP=load(PATH);
% PP.XYH;
% PP.MAPH;
% PP.symmetry;
% PP
% [XYM,MAPM] = mech_mesh_gen(PP.XYH,PP.MAPH);
% ele_no=1; node_no=0; m='o';
% draw_mech_9noded_2D(XYM,MAPM,ele_no,node_no,m);
% mech_solver()
function[]=mech_solver(PATH_Mech,PATH_sim)
global XYM Lz MAPM MECHBC BEAMBC neM nnM nMBC K G CTET
global nK nG nCTET TRef MAPH  neH_o 
TOL=1e-5;
Lz=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PATH_Process=[PATH_sim 'Process_Parameters.mat'];
PATH_Material_Mech=[PATH_Mech 'Material_TM.mat'];
PATH_Material=[PATH_sim 'Material.mat'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PP=load(PATH_Process);load(PATH_Material);

XYH=PP.XYH;
MAPH=PP.MAPH;
% [XYM,MAPM] = mech_mesh_gen(XYH,MAPH);
[XYM,MAPM,neH_o]=Mesh_Conversion(XYH,MAPH);
XYM=XYM.*1e-3;
[MECHBC]=MECH_BC(XYM,TOL,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ele_no=0; node_no=0; m='o';
% draw_mech_9noded_2D(XYM,MAPM,ele_no,node_no,m)


load(PATH_Material_Mech);
% [K,G] = ENuKG(EyT,NuT);

nnH       = size(XYH,1);
nnM       = size(XYM,1);
neM        = size(MAPM,1);
nMBC        = size(BEAMBC,1);
nK        = size(K,1);
nG        = size(G,1);
nCTET        = size(CTET,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nDBs =0;
PathName = [PATH_sim 'DB' int2str(0) '.MAT'];
load(PathName);
TRef=(TL+TS)*0.5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mkdir([PATH FN]);
% PATH_Mech=[PATH FN '\'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Initialize_DB(PATH_Mech,T);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PathName = [PATH_Mech 'DB_M_' int2str(0) '.MAT'];
load(PathName);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
IPTEMP_p=IPTEMP;    % Previous Integration Point Temperature
nlts=2;              % number of local time steps
%%% Count Number of DB files %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d=dir([PATH_sim 'DB*.mat']);
str = {d.name};
ndir=size(str,2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
for i=1:ndir-2
    %%%%%%%%% Read Thermal Solution
    PathName = [PATH_sim 'DB' int2str(i) '.MAT'];
    load(PathName);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for j=1:nlts
        IPTEMPt=IPTEMP_p+(IPTEMP-IPTEMP_p)*j/nlts;
        [U,IPD]=MECH_EqmIter(U,IPD,IPTEMPt);
    end
    DB_NAME = [PATH_Mech 'DB_M_' int2str(nDBs+i) '.MAT'];
    save(DB_NAME,'t','Axial','U','IPD');
    if i==10
    [MECHBC]=MECH_BC(XYM,TOL,1);
    end
    display(['DB-' num2str(i) '- Axial -' num2str(Axial) ' m'])
end % for (it=0:dt:te)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
toc

function[]=Make_Mech_Material_File(MECH_MAT,PATH_Mech)
     EyT=MECH_MAT.EyT;
     NuT=MECH_MAT.NuT;
    CTET=MECH_MAT.CTET;
     SyT=MECH_MAT.SyT;
    HRDT=MECH_MAT.HRDT;
      K=MECH_MAT.KT;
      G=MECH_MAT.GT;
%       PATH_Mech
       save([PATH_Mech 'Material_TM.mat'],'EyT','NuT','CTET','SyT','HRDT','K','G');
       
function xvl = val_at_t(X,nX,T)
%  Interpolationfor piecewise linear data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  INPUT:                                                                 %
%    [X] : Matrix with size Nx2                                           %
%    [T] : Temperature                                                    %
%  OUTPUT:                                                                %
%    [xvl] : Linearly interpolated value                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if T<X(1,1) % Make linear extrapolation
  xvl = X(1,2)-(X(1,1)-T)*(X(2,2)-X(1,2))/(X(2,1)-X(1,1));
elseif T>X(nX,1) % Make linear extrapolation
  xvl = X(nX,2)+(T-X(nX,1))*(X(nX,2)-X(nX-1,2))/(X(nX,1)-X(nX-1,1));
else
  for jj=2:nX
    if T<=X(jj,1) % Make linear interpolation
      xvl = X(jj-1,2)+(T-X(jj-1,1))*(X(jj,2)-X(jj-1,2))/(X(jj,1)-X(jj-1,1));
      break;
    end
  end
end

function [K,G] = ENuKG(E,Nu)
%  Computation of Bulk modulus K and Shear Modulus G   
%  Sub-function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  This function is called only once after inputdata to initialize some      %
%  Bulk moduslus K and Shear modulus G                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize K and G
K = E;
G = E;
npt = size(E,2); % number of phases + 1
ip=2:npt; % loop over the phases
K(:,ip) = E(:,ip)./(1. - 2.*Nu(:,ip))/3.;
G(:,ip) = E(:,ip)./(1. + Nu(:,ip))/2.;

function[MECHBC]=MECH_BC(XYM,TOL,a)
% Mechanical boundary conditions
TOL = eps;
MECHBC = zeros(0,5); %NodeId type(0/1) Fx/Ux type(0/1) Fy/Uy
Xmax=max(XYM(:,1));
Ymax=max(XYM(:,2));
Ymin=min(XYM(:,2));

NodeIds = nodes_on_line(XYM,[0.000 Ymin],[0 Ymax],TOL);
MECHBC = [MECHBC; [NodeIds repmat([1 0e0 a 0e0],size(NodeIds))]]; 

% NodeIds = nodes_on_line(XYM,[0.000 0.000],[Xmax 0],TOL);
% MECHBC = [MECHBC; [NodeIds repmat([0 0e0 1 0e0],size(NodeIds))]]; 

NodeId = node_on_pt(XYM,[0.000 0.000],TOL);
MECHBC = [MECHBC; [NodeId,[1 0e0 1 0e0]]];

% NodeId = node_on_pt(XYM,[Xmax 0],TOL);
% MECHBC = [MECHBC; [NodeId,[1 0e0 1 0e0]]];
% 
% NodeId = node_on_pt(XYM,[0 Ymax],TOL);
% MECHBC = [MECHBC; [NodeId,[1 0e0 1 0e0]]];

function[]=Initialize_DB(PATH,T)
global  MECHBC BEAMBC neM nnM  XYM MAPM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%INITIAL DATABASE, DB0.MAT
%  T      -> Nodal initial temperature
%  UV     -> Nodal diplacements
%  Ub     -> Beamd dofs displacements
%  IPETR  -> Integration point TRIP strain      [Exx  Eyy  Exy  Ezz]'
%  IPEPL  -> Integration point Plastic strain   [Exx  Eyy  Exy  Ezz]'
%  IPSTRS -> Integration point Stress           [Sxx  Syy  Sxy  Szz]'
%  IPPH   -> Integration point f,tI,tC          [f tI tC]
%  IPLP   -> Integration point Hardening variable,
%  IPTEMP -> Integration point Temperature,
%  IPDT   -> Integration point Temperature icrement
%  IPDQ   -> Integration point Generated letant heat [epl T dT dQ]'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Global time (s)
t = 0;
% Nodal Displacements (m)
nbc = size(MECHBC,1); % number of nodes with mech. BC
dof = 2*nnM+3;
U  = zeros(dof,1);
% Apply displacement boundary conditions
for ibc=1:nbc
  if MECHBC(ibc,2) % ux
    U(2*MECHBC(ibc,1)-1,1) = MECHBC(ibc,3);
  end % x
  if MECHBC(ibc,4) % uy
    U(2*MECHBC(ibc,1),1) = MECHBC(ibc,5);
  end % y
end % for (ibc = 1:nbc)
% Apply Beam Displacement boundary conditions
if BEAMBC(1,1)
  U(dof-2) = BEAMBC(1,2);
end % Exx
if BEAMBC(2,1)
  U(dof-1) = BEAMBC(2,2);
end % Kxx
if BEAMBC(3,1) 
  U(dof) = BEAMBC(3,2);
end % Kyy

% Integration point data
% IPETR  = zeros(4,9,ne);    % TRIP strain      [Exx  Eyy  Exy  Ezz]'
% IPEPL  = zeros(4,9,ne);    % Plastic strain   [Exx  Eyy  Exy  Ezz]'
STRESS = zeros(4,9);    % Stress           [Sxx  Syy  Sxy  Szz]'
% IPLP   = zeros(ne,9);      % Hardening variable
% IPTEMP = zeros(ne,9); % Temperature
% [IPTEMP] = set_ip_temp(T,IPTEMP);
%%%%%%%%%%%%%
a=struct('STRESS',STRESS);
for i=1:neM
    IPD{i}=a;
end
Axial=0;
save([PATH 'DB_M_0.MAT'],'t','U','IPD','MAPM','XYM','Axial');
 

function [IPTEMP] = set_ip_temp(T,IPTEMP)
%  Subfunction: Sets temperature and temperature increment for each IP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  INPUT:                                                                 %
%    [T]      : Nodal temperatures                                        %
%    [IPTEMP] : IP temperature data from previous time step               %
%  OUTPUT:                                                                %
%    [IPDT]   : IP temperature increment data for current time step       %
%    [IPTEMP] : IP temperature data for current time step                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global MAPH ne;

% Integration point coordinates in local coordinate system
pt = 0.774596669241483;

% Nodal Coordinates of parent element and node numbering
ri = [+1;-1;-1;+1;+0;-1;+0;+1;+0]; % Coefficient for r
si = [+1;+1;-1;-1;+1;+0;-1;+0;+0]; % Coefficient for s

IPDT = IPTEMP;

for ie=1:ne
  Te  = T(MAPH(ie,1:4)');      % Element nodal temperatures
  for ip=1:9
    r = pt*ri(ip);
    s = pt*si(ip);
    N = 0.25*[(1+r)*(1+s);
              (1-r)*(1+s); 
              (1-r)*(1-s); 
              (1+r)*(1-s)];
    Tip            = N'*Te;  % ip current temperature
    IPDT(ie,ip)    = Tip - IPTEMP(ie,ip); % set temperature increment
    IPTEMP(ie,ip)  = Tip;                 % set temperature
  end % for (ip=1:9)
end % for (ie=1:ne)

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

function NodeIds = nodes_on_line(XY,StPt,EndPt,TOL)
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
  if abs(StPt(1)*(EndPt(2)-XY(in,2)) ...
        + EndPt(1)*(XY(in,2)-StPt(2)) ...
        + XY(in,1)*(StPt(2)-EndPt(2))) <= TOL
    cnt = cnt + 1;
    NodeIds(cnt,1) = in;
  end % if
end % for(in=1:nn)    % loop over nodes
