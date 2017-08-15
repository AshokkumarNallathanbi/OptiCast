function[Axial]=Solver_New(PathName,nDB,path_sep)
% PathName='D:\HAZELETT CASTER\SWISS STEEL SOFTWARE\TABBING Version\Final_2-11-04-2014\Simulations\1045_S275JR';
PathName=[PathName path_sep];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% nDBs   - start database number
%%%% ts     - solution start time
%%%% tf     - solution end time
%%%% N      - Number of database
%%%% EL__PL - Elastic/Elasto-plastic tag
%%%% nlts   - Number local time steps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global KHT RHT SHT nKHT nRHT nSHT 
global LATENT TS TL 
global neH nIPH TOL_T emissivity nnH
global XY MAPH MHF nMHF HTC X Z  HTC MHF Vc Tinf BC_Info2
global BC_TAG C_BC_TAG Sim_Len Sec_Len Mol_Len
global PCF nPCF  symmetry
global f1 f2 d_chamfer Wi Th faces
%%%%%%%%%%%%%%%%%%%%%%%%
load([PathName 'Material.mat']);
load([PathName 'Process_Parameters.mat']);
nPCF=size(PCF,1);

X=SEC_COOL.X;
Z=SEC_COOL.Z;
HTC=SEC_COOL.HTC;

load([PathName 'DB' int2str(nDB) '.mat']);
nIPH=size(IPPH,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Sec_Len=Nozzle_Tab(size(Nozzle_Tab,1),1)+1000;
Tinf=Tatm;
emissivity=0.8;
XY=XYH.*1e-3;
neH=size(MAPH,1);nnH=size(XY,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Vc=Grade_depend.Casting_speed/60;
LATENT=Latent;
% LATENT=250e3;
% LATENT=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Simulation control parameters
maxIt = 50;           % Maximum iterations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k_max =maxIt;order=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nKHT = size(KHT,1);
nRHT  = size(RHT,1);
nSHT  = size(SHT,1);
TOL_T = TOL;  % Tolerance in Thermal problem
a     = 1e-2; % very small number used in Phase fraction calculation in Isothermal case
% Load mesh and boundary conditions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
BC_TAG=[1 Mol_Len*1e-3;         %%%  1 -HF, 2-HTC+rad, 3-rad
        2 Sec_Len*1e-3
        3 Sim_Len*1e-3];        %%%  Secondary cooling
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Mold BC
MHF(:,1)=MHF(:,1).*1e-3;
MHF(:,2)=MHF(:,2).*1e6;
nMHF=size(MHF,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
faces=[4 1;2 3;1 2;3 4];

it1_max=3;alphmina=0;alphmaxa=25;amp=25;LSTOL=0.8;
stop_flag=0;C_BC_TAG=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    t          = t+dt;  
    Axial      =   Vc*t; %%%% distance from meniscus
    IPTEMP_p    =   IPTEMP;T_p=T;
    if Axial<BC_TAG(1,2) %%%%% MOLD
        C_BC_TAG=1;
    elseif (Axial>BC_TAG(1,2)) && (Axial<BC_TAG(2,2)) %%%% Secondary cooling
        C_BC_TAG=2;
    elseif (Axial>BC_TAG(2,2)) 
        C_BC_TAG=3;
    end
    [T,IPPH,IPDQ,NIT]= solve_solidification_problem(T,IPPH,IPDQ,dt,t);
    IPTEMP = set_ip_Temperature(T,nIPH);
    display(['DB-' num2str(nDB) '- Axial -' num2str(Axial) ' m'])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DB_NAME_T   = [PathName 'DB' int2str(nDB+1) '.mat'];
save(DB_NAME_T,'t','T','IPPH','IPTEMP','IPDT','IPDQ','Axial');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [T_h] = control_incre_temp(T_h,Th_p)
          n=length(T_h);
          for i=1:n
              T1=T_h(i);
              T0=Th_p(i);
              dT=T1-T0;
              if dT>0
                  T_h(i)=T0;
              end
          end
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
function[T,IPPH,IPDQ,N]= solve_solidification_problem(T,IPPH,IPDQ,dt,t)
        global neH nIPH nnH TOL_T maxIt LSTOL
        opts.SYM    = true;
        opts.POSDEF = true;
        DFDT=zeros(neH,nIPH);
        dT   = zeros(nnH,1);     % Initial temperature increment
        norm_check=1;
        T_pre=T;
        [K,C,F,L,Cpc] = assamble_heatPH(T,IPPH,IPDQ,DFDT,t,T_pre);
        L_pre=L;
        iter=0;
        ls_it=4;
        second_norm=1;
        
        while (norm_check>TOL_T) && (second_norm>1e-20) 
            G=[];S=[];alph=[];ICO=0;k=2;
            T1=T;                % Temp. at the beginning of iteration
            [A,B] = Calculate_AB(T1,dt,T_pre,L_pre,IPDQ,t);
            dT = linsolve(A,B,opts);
            G(1) = dT'*B;S(1)=G(1)/G(1);alph(1)=0;
            T2    = T1 + dT;
            [A,B] = Calculate_AB(T2,dt,T_pre,L_pre,IPDQ,t);
            G(2) = dT'*B;S(2)=G(2)/G(1);alph(2)=1;it=2;
            ETA=1;
%             %%%%%%%%%%%%%%%%%%%%%%%%%%%%% LS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %             if abs(S(2))>LSTOL
% %                 [ETA]=call_line_search(G,S,alph,dT,T1,dt,T_pre,L_pre,IPDQ);
% %             else
% %                 ETA=1;
% %             end
            while (abs(S(it))> 0.8) && (it<ls_it)
                alph(it+1)=alph(it-1)-S(it-1)*(alph(it)-alph(it-1))/(S(it)-S(it-1));
                T2   = T1 + alph(it+1)*dT;
                [A,B] = Calculate_AB(T2,dt,T_pre,L_pre,IPDQ,t);
                G(it+1) = dT'*B; S(it+1) = G(it+1)/G(1);
                it=it+1;
                if S(it)<0
                    S(it-1)=S(1);alph(it-1)=0;
                end
            end
            ETA=alph(it);
            if ETA<0
                ETA=0.01;
            elseif ETA>1
                ETA=1;
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            T = T1+(ETA*dT);
            [DFDT,IPPH]=calculate_DFDT(T,T_pre);
            [K,C,F,L,Cpc] = assamble_heatPH(T,IPPH,IPDQ,DFDT,t,T_pre);
            B=F-(C*((T-T_pre)/dt))-(K*T)-(L-L_pre)/dt;
            A=K+(C/dt)+(Cpc/dt);
%             [A,B]=EssentialBC(A,B);
            norm_check = (B'*B)/((K*T)'*(K*T));
            second_norm=(dT'*dT);
            iter=iter+1;
            if iter>maxIt
                warning('***TemperatureSolutionDiverged');
                retVal = 1e14;
                return;
            end
        end
        N=iter;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

function [A,B]=Calculate_AB(T1,dt,T_pre1,L_pre1,IPDQ1,t)
[DFDT,IPPH]=calculate_DFDT(T1,T_pre1);
[K,C,F,L,Cpc] = assamble_heatPH(T1,IPPH,IPDQ1,DFDT,t,T_pre1);
B=F-(C*((T1-T_pre1)/dt))-(K*T1)-(L-L_pre1)/dt;
A=K+(C/dt)+(Cpc/dt);
% [A,B]=EssentialBC(A,B);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
function [KG,CG,FG,LG,CGpc] = assamble_heatPH(T,IPPH,IPDQ,DFDT,t,T_pre)
%  Subfunction: Assables heat element matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  INPUT:                                                                 %
%    [T]      : Nodal temperatures at all nodes                           %
%    [IPPH]   : Phase fractions at integration points                     %
%    [IPDQ]   : Latent heat release at integration points                 %
%    [DFDT]   : derivative of phase fraction function w.r.t. temperature         %
%  OUTPUT:                                                                %
%    [KG] : Global heat conduction matrix                                 %
%    [CG] : Global heat capacity matrix                                   %
%    [FG] : Global heat force vector                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global XY MAPH neH nnH Vc Z HTC X MHF nMHF
global   BC_TAG C_BC_TAG  BC_Info2  faces
global ns_htc TS
% Initializations
KG = zeros(nnH);
CG = zeros(nnH);
FG = zeros(nnH,1);
CGpc = zeros(nnH);
LG = zeros(nnH,1);

for ie=1:neH
    map=MAPH(ie,:);
  % Compute element matrices and vectors
  [Ke,Ce,Fe,Le,Cpce] = heat_elementPH(XY(map,1:2),T(map'),...
      IPPH(ie,:),IPDQ(ie,:),DFDT(ie,:),t);   
  % Assamble
  KG(map,map)     = KG(map,map)     + Ke;
  CG(map,map)     = CG(map,map)     + Ce;
  CGpc(map,map)   = CGpc(map,map)   + Cpce;
  FG(map')        = FG(map')        + Fe;
  LG(map')        = LG(map')        + Le;
end 
TYPE=BC_TAG(C_BC_TAG,1);
if  TYPE==1  %%% mold cooling
    HF=-Linear_Interpolation(MHF,nMHF,Vc*t);
elseif TYPE==2   %%% secondary cooling
%     HTC1=find_htc_at_z(Z.*1e-3,HTC,Vc*t);
    HTC1=Space_Interpolation(Z(1,:).*1e-3,Vc*t,HTC);
    Shtc=[X(:,1).*1e-3 HTC1];
end
faces1=[4 1;3 2;2 1;3 4];

for ie=1:size(BC_Info2,1)
    Ke=zeros(4);Fe=zeros(4,1);
    ne=BC_Info2(ie,1);
    nf=BC_Info2(ie,2);
    map2          =   MAPH(ne,:);
    f=compute_factor(nf,XY(map2,:));
    Tf=T_pre(map2(faces(nf,:)));
    FACT=zeros(4,1);
    if Tf(1)>TS
        f(1)=1;
    end
    if Tf(2)>TS
        f(2)=1;
    end
%     FACT(faces1(nf,:))=1;
%     FACT(faces1(nf,1))=BC_Info2(ie,3);FACT(faces1(nf,2))=BC_Info2(ie,4);
    FACT(faces(nf,1))=f(1);
    FACT(faces(nf,2))=f(2);
    if (TYPE==1)
        [Fe]=APPLY_HF(XY(map2,1:2),FACT.*HF,nf);
    elseif  TYPE==2        
        [Ke,Fe]=APPLY_HTC(XY(map2,1:2),T(map2'),Shtc,nf);
    else
        [Ke,Fe]=APPLY_Radiation(XY(map2,1:2),T(map2'),nf);
    end
    KG(map2,map2) =   KG(map2,map2) + Ke;
    FG(map2')     =   FG(map2')     + Fe;
    %     pause
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function[f]=compute_factor(nf,XY)
        global f1 f2 d_chamfer Wi Th faces symmetry
        LNid=faces(nf,:);
        xL=XY(LNid,1)*1e3;
        yL=XY(LNid,2)*1e3;
        if symmetry==0
            Th1=Th/2;Wi1=Wi/2;
        elseif symmetry==3
            Th1=0;Wi1=0;
        elseif  symmetry==2
            Th1=Th/2;Wi1=0;
        elseif  symmetry==1
            Th1=0;Wi1=Wi/2;
        end
        
            if (nf==1) || (nf==2)
                zL=abs(yL-Th1);z1=Th;
            else
                zL=abs(xL-Wi1);z1=Wi;
            end
        for i=1:length(zL)
            if zL(i)<(z1/2-d_chamfer)
                f(i)=1+(2*zL(i)*(f1-1)/(z1-2*d_chamfer));
            else
                f(i)=(((f2-f1)/2/d_chamfer)*(2*zL(i)-z1+2*d_chamfer))+f1;
            end
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function[y1]=find_htc_at_z(x,y,xn)
        x1=x(1,:);
        for i=1:length(x1)-1
%            [x1(i) x1(i+1) xn ]
            if x1(i)<=xn && x1(i+1)>=xn
                n1=i;n2=i+1;
                break;
            end
            if x1>xn
                n1=i+1;n2=i+1;
            end
        end
       y1=((y(:,n2)-y(:,n1))./(x1(n1+1)-x1(n1)).*(xn-x1(n1)))+y(:,n1);     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[Ke,Fe]=APPLY_HTC(XYe,Te,Shtc,nf)
        global radial_dir nCONVBC emissivity Tinf symmetry Wi Th
        pt = 0.57735026918962576;
        ipr = [1 1;-1 -1;-pt pt;-pt pt];
        ips = [-pt pt;-pt pt;1 1;-1 -1];
        dam = [2 2 1 1]; % DA mapping         
        Ke=zeros(4);Fe=zeros(4,1);
        if symmetry==0
            fa =[0 0 0 0];
        elseif symmetry==3
            fa =[Th/2 0 Wi/2 0].*1e-3;
        elseif  symmetry==2
            fa =[0 0 Wi/2 Wi/2].*1e-3;
        elseif  symmetry==1
            fa =[Th/2 Th/2 0 0].*1e-3;
        end
        for ip=1:2 % loop over ineHgration points
            r = ipr(nf,ip);
            s = ips(nf,ip);
            [N,J] = heat_interpol_NJ(XYe,r,s);
            DA = norm(J(dam(nf),:));   % Area element         
            Temp=N'*Te;
            xyip  = N'*XYe(1:4,dam(nf));
            hf=interp1(Shtc(:,1),Shtc(:,2),xyip+fa(nf));
            %%%%%%%%%%%%%%%%
            ht=hf;
%             ht=hf/(Temp-Tinf);
            %%%%%%%%%%%%%%%%%%
            emm=0.0002*Temp+0.6274;
            ht_rad=emm*5.67e-8*((Tinf+273)^2+(Temp+273)^2)*((Tinf+273)+(Temp+273));
%             ht_rad=0;
            hte=ht+ht_rad;
            Ke = Ke + DA*hte*(N*N');     % Conductivity matrix
            Fe = Fe + DA*hte*N*Tinf;

        end 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%-------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[Ke,Fe]=APPLY_Radiation(XYe,Te,nf)
        global radial_dir nCONVBC emissivity Tinf
        pt = 0.57735026918962576;
        ipr = [1 1;-1 -1;-pt pt;-pt pt];
        ips = [-pt pt;-pt pt;1 1;-1 -1];
        dam = [2 2 1 1]; % DA mapping         
        Ke=zeros(4);Fe=zeros(4,1);
        for ip=1:2 % loop over ineHgration points
            r = ipr(nf,ip);
            s = ips(nf,ip);
            [N,J] = heat_interpol_NJ(XYe,r,s);
            DA = norm(J(dam(nf),:));   % Area element         
            Temp=N'*Te;          
            hte=0.8*5.67e-8*((Tinf+273)^2+(Temp+273)^2)*((Tinf+273)+(Temp+273));
            Ke = Ke + DA*hte*(N*N');     % Conductivity matrix
            Fe = Fe + DA*hte*N*Tinf;
        end 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%-------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[Ke,Fe]=APPLY_GENERAL_HTC_BC(XYe,Te,nf,CBC,fact,TYPE)
        global radial_dir nCONVBC emissivity
        pt = 0.57735026918962576;
        ipr = [1 1;-1 -1;-pt pt;-pt pt];
        ips = [-pt pt;-pt pt;1 1;-1 -1];
        fcm = [4 1;3 2;2 1;3 4]; % Face mapping
        dam = [2 2 1 1]; % DA mapping         
        Ke=zeros(4);Fe=zeros(4,1);
         THT = zeros(4,1);  % Heat transfer coef.
         TINF = zeros(4,1); % Environment (or Radiant) Temperature
         THT(fcm(nf,1)) = val_at_t(CBC(1:nCONVBC,[1 2]),nCONVBC,Te(fcm(nf,1)))*fact(1);
         THT(fcm(nf,2)) = val_at_t(CBC(1:nCONVBC,[1 3]),nCONVBC,Te(fcm(nf,2)))*fact(2);
         TINF(fcm(nf,1))= CBC(nCONVBC+1,2);
         TINF(fcm(nf,2))= CBC(nCONVBC+1,3);
         R=1;
        for ip=1:2 % loop over ineHgration points
            r = ipr(nf,ip);
            s = ips(nf,ip);
            [N,J] = heat_interpol_NJ(XYe,r,s);
            xyip  = N'*XYe(1:4,:);
%             if radial_dir==0
%                 R=1;
%             else
%                 R=xyip(radial_dir);
%             end
            DA = norm(J(dam(nf),:))*R;   % Area element
            ht   = N'*THT ;             % Heat transfer coef.
            Tinf  = N'*TINF; % Environment (or Radiant) Temperature           
            if TYPE==2
                hte=ht;
            elseif TYPE==4
                Temp=N'*Te;
                ht_rad=0.8*5.67e-8*((Tinf+273)^2+(Temp+273)^2)*((Tinf+273)+(Temp+273));
                hte=ht+ht_rad;
            end
%             [ht,ht_rad]
            Ke = Ke + DA*hte*(N*N');     % Conductivity matrix
            Fe = Fe + DA*hte*N*Tinf;
        end 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%-------------------------------------------
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function[Fe]=APPLY_HF(XYe,HF,nf)
        global radial_dir BCEL BCFC BCNID nCONVBC
        pt = 0.57735026918962576;
        ipr = [1 1;-1 -1;-pt pt;-pt pt];
        ips = [-pt pt;-pt pt;1 1;-1 -1];
%         fcm = [4 1;2 3;1 2;3 4]; % Face mapping
        dam = [2 2 1 1]; % DA mapping
        Fe=zeros(4,1);
        for ip=1:2 % loop over ineHgration points
            r = ipr(nf,ip);
            s = ips(nf,ip);
            [N,J] = heat_interpol_NJ(XYe,r,s);
            DA = norm(J(dam(nf),:));   % Area element
            ht   = N'*HF;             % Heat transfer coef.
            Fe = Fe + DA*ht*N;
        end
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%-------------------------------------------
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 function[Fe]=APPLY_GENERAL_HF_BC(XYe,Te,nf,CBC,fact)
        global radial_dir BCEL BCFC BCNID nCONVBC
        pt = 0.57735026918962576;
        ipr = [1 1;-1 -1;-pt pt;-pt pt];
        ips = [-pt pt;-pt pt;1 1;-1 -1];
        fcm = [4 1;3 2;2 1;3 4]; % Face mapping
        dam = [2 2 1 1]; % DA mapping
        R=1;
        Fe=zeros(4,1);
         THT = zeros(4,1);  % Heat transfer coef.
         THT(fcm(nf,1)) = val_at_t(CBC(1:nCONVBC,[1 2]),nCONVBC,Te(fcm(nf,1)))*fact(1);
         THT(fcm(nf,2)) = val_at_t(CBC(1:nCONVBC,[1 3]),nCONVBC,Te(fcm(nf,2)))*fact(2);
         for ip=1:2 % loop over ineHgration points
            r = ipr(nf,ip);
            s = ips(nf,ip);
            [N,J] = heat_interpol_NJ(XYe,r,s);
            xyip  = N'*XYe(1:4,:);
%             if radial_dir==0
%                 R=1;
%             else
%                 R=xyip(radial_dir);
%             end
            DA = norm(J(dam(nf),:))*R;   % Area element
            ht   = N'*THT;              % Heat transfer coef.
            Fe = Fe + DA*ht*N;
        end 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%-------------------------------------------
function [Ke,Ce,Fe,Le,Cpce] = heat_elementPH(XYe,Te,IPPHe,...
        IPDQe,DFDTe,t)
%  Subfunction: Computes heat transfer element matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  INPUT:                                                                 %
%    [Te]      : Element nodal temperature vector                          %
%    [IPPH]   : Element IP formed phase fractions                         %
%    [IPDQ]   : Element IP latent heat geneHratoin                         %
%    [ELBC]   : Element heat transfer BCs  /HTBC for the element          %
%    [CONVBC] : Array convection values   
%     [MIe]   : Element material index
%  OUTPUT:                                                                %
%    [Ke] : Element heat stiffneHss matrix                                 %
%    [Ce] : Element heat capacity matrix                                  %
%    [Fe] : Element heat force vector                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global KHT RHT SHT  nKHT nRHT nSHT  LATENT 
global nIPH  TS TL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%9pt%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Integration point coordinates in local coordinate system
% if nIPH==16
%     Pt1=0.861136311594052;
%     Pt2=0.339981043584856;
%     wgt1=0.121002993285602;wgt2=0.425293303010694;wgt3=0.226851851851852;
% 
%     ri=[-Pt1;-Pt2;+Pt2;+Pt1;-Pt1;-Pt2;+Pt2;+Pt1;-Pt1;-Pt2;+Pt2;+Pt1;-Pt1;-Pt2;+Pt2;+Pt1];
%     si=[-Pt1;-Pt1;-Pt1;-Pt1;-Pt2;-Pt2;-Pt2;-Pt2;+Pt2;+Pt2;+Pt2;+Pt2;+Pt1;+Pt1;+Pt1;+Pt1];
%     wg=[wgt1;wgt3;wgt3;wgt1;wgt1;wgt2;wgt2;wgt3;wgt3;wgt2;wgt2;wgt3;wgt1;wgt3;wgt3;wgt1];
% elseif nIPH==9;
   % Integration point coordinates in local coordinate system
    pt = 0.774596669241483;
    wg = [25 25 25 25 40 40 40 40 64]/81;
    % %     Nodal Coordinates of parent element and node numbering
    ri = pt*[+1;-1;-1;+1;+0;-1;+0;+1;+0]; % Coefficient for r
    si = pt*[+1;+1;-1;-1;+1;+0;-1;+0;+0]; % Coefficient for s
% elseif nIPH==4
%     pt = 0.57735026918962576;
%     ri = pt*[+1;-1;-1;+1];
%     si = pt*[+1;+1;-1;-1];
%     wg = [1 1 1 1];
% end
% Initialize Element Marrices 
Ke = zeros(4);   % Conducton Matrix
Ce = zeros(4);   % Capacity Matrix
Cpce = zeros(4);   % Latent heat Matrix
Fe = zeros(4,1); % Force vector
Le = zeros(4,1); % Latent heat vector
% VOLUME INTEGRAL
for ip=1:nIPH %%% 
    r =  ri(ip);
    s =  si(ip);
    [N,Nxy,DJAC] = heat_interpol(XYe,r,s);
%     xyip  = N'*XYe(1:4,:);
%     if radial_dir==0
%         R=1;
%     else
%         R=xyip(radial_dir);
%     end
    DVOL = wg(ip)*DJAC; % Volume element
    Temp = N'*Te;       % Integration point temperature
    f    = IPPHe(ip); % Integration point phase fractions
    

KHTm = Linear_Interpolation(KHT,nKHT,Temp);
% KHTmix=KHTm*(7-6*(1-f));
KHTmix=KHTm;
SHTmix = Linear_Interpolation(SHT,nSHT,Temp);
RHTmix = Linear_Interpolation(RHT,nRHT,Temp);


%   KHTmix = mixture_rule(KHT,nKHT,f,Temp,1);% Mixture conductivity, AM (N=1)
%   Romix  = mixture_rule(Ro,nRo,f,Temp,1);  % Mixture density, AM (N=1)
%   SHmix  = mixture_rule(SH,nSH,f,Temp,1);  % Mixure specific heat, AM (N=1)


    % Gaus-weights iwr = iws = iwt = 1 so they are not included
    Ke    =   Ke      +   DVOL*KHTmix*(Nxy*Nxy');  % Conductivity matrix
    Ce    =   Ce      +   DVOL*RHTmix*SHTmix*(N*N'); %  Capacitance matrix
    Fe    =   Fe      +   DVOL*RHTmix*IPDQe(ip)*N;  % Force Vector
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Le    =   Le      +   DVOL*RHTmix*LATENT*IPPHe(ip)*N; % LATENT heat vector
    Cpce  =   Cpce    +   DVOL*RHTmix*LATENT*DFDTe(ip)*(N*N'); % LATENT heat matrix
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
function [N,Nxy,DJAC] = heat_interpol(XY,r,s)
%  Subfunction: Isoparametric element interpolation functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  INPUT:                                                                 %
%    [XY] : Element nodal coordinates of the element                      %
%    [r]  : Local r-coordinate where the output is to be computed         %
%    [s]  : Local s-coordinate where the output is to be computed         %
%  OUTPUT:                                                                %
%    [N]    : Shape funtion vector                                        %
%    [Nxy]  : Shape funtion derivative wrt x and y                        %
%    [J]    : Jacobian matrix                                             %
%    [DJAC] : Determinant of the jacobian matrix                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Nodal point coordinates and node numbering
%ri = [+1; -1; -1; +1]; % Coefficient for r
%si = [+1; +1; -1; -1]; % Coefficient for s

% Interpolation function
N = 0.25*[(1+r)*(1+s);
          (1-r)*(1+s); 
          (1-r)*(1-s); 
          (1+r)*(1-s)];

% Derivative wrt r
Nrs = 0.25*[ 1+s  1+r; 
            -1-s  1-r; 
            -1+s -1+r; 
             1-s -1-r];

J    = Nrs'*XY;
% Jacobian matrix 
DJAC = J(1,1)*J(2,2) - J(2,1)*J(1,2); % Jacobian determinant
iJt  = [ J(2,2)/DJAC -J(2,1)/DJAC;
        -J(1,2)/DJAC  J(1,1)/DJAC];
 Nxy  = Nrs*iJt;  % Derivatives wrt global coordinates

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[DFDT,IPPH]=calculate_DFDT(T_cur,T_pre)
global TS TL MAPH a nIPH
neH=size(MAPH,1);
DFDT=zeros(neH,nIPH);
IPT_cur=set_ip_Temperature(T_cur,nIPH);
IPT_pre=set_ip_Temperature(T_pre,nIPH);
for ie=1:neH
    for ip=1:nIPH
        IPPH(ie,ip)=phasefraction(IPT_cur(ie,ip));
        DF=phasefraction(IPT_cur(ie,ip))-phasefraction(IPT_pre(ie,ip));
        DT=IPT_cur(ie,ip)-IPT_pre(ie,ip);
        if abs(DT)>0
            DFDT(ie,ip)=DF/DT;
        else
            DFDT(ie,ip)=0;
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [PH]=phasefraction(Temp)
global TS TL non_iso_limit TOL a PCF nPCF
%a=1e-10;
% Non-isothermal case
if abs(TS-TL)>0 % abs(TS-TL)>TOL
  if Temp>TL
    PH=1;
  elseif Temp<=TS
    PH=0;
  else
%     PH=(Temp-TS)/(TL-TS);
%     PH = val_at_t(PCF,nPCF,Temp);
    PH = Linear_Interpolation(PCF,nPCF,Temp);
  end
else % Isothermal case
  x=Temp-0.5*(TS+TL);
  PH=((1/pi)*atan(x/a))+0.5;
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [PH]=phasefraction(Temp)
%     global PCF nPCF
%    PH = val_at_t(PCF,nPCF,Temp);
% %     PH=interp1(PCF(:,1),PCF(:,2),Temp);
%     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  function [a,b]=EssentialBC(A,B)
  global T_BC_nodes
  a=A;b=B;
  if  (~isempty(T_BC_nodes))
    for i=1:length(T_BC_nodes)
      nid=T_BC_nodes(i);
      A(:,nid)=0;
      A(nid,:)=0;
      A(nid,nid)=1;
      B(nid,1)=0;
    end
  else
    return;
  end
  a=A;b=B;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [N,J] = heat_interpol_NJ(XY,r,s)
%  Subfunction: Isoparametric element interpolation functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  INPUT:                                                                 %
%    [XY] : Element nodal coordinates of the element                      %
%    [r]  : Local r-coordinate where the output is to be computed         %
%    [s]  : Local s-coordinate where the output is to be computed         %
%  OUTPUT:                                                                %
%    [N]    : Shape funtion vector                                        %
%    [Nxy]  : Shape funtion derivative wrt x and y                        %
%    [J]    : Jacobian matrix                                             %
%    [DJAC] : Determinant of the jacobian matrix                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Nodal point coordinates and node numbering
%ri = [+1; -1; -1; +1]; % Coefficient for r
%si = [+1; +1; -1; -1]; % Coefficient for s

% Interpolation function
N = 0.25*[(1+r)*(1+s);
          (1-r)*(1+s); 
          (1-r)*(1-s); 
          (1+r)*(1-s)];
% Derivative wrt r
Nrs = 0.25*[ 1+s  1+r; 
            -1-s  1-r; 
            -1+s -1+r; 
             1-s -1-r];
J    = Nrs'*XY; % Jacobian matrix 
%%%%%%%%%%%%%%%%%%%%%%%%usage in DFDT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[IPTEMP]=set_ip_Temperature(TEMP,nIP)
global MAPH neH
% if nIP==16
%     Pt1=0.861136311594052;
%     Pt2=0.339981043584856;
%     wgt1=0.121002993285602;wgt2=0.425293303010694;wgt3=0.226851851851852;
%     ri=[-Pt1;-Pt2;+Pt2;+Pt1;-Pt1;-Pt2;+Pt2;+Pt1;-Pt1;-Pt2;+Pt2;+Pt1;-Pt1;-Pt2;+Pt2;+Pt1];
%     si=[-Pt1;-Pt1;-Pt1;-Pt1;-Pt2;-Pt2;-Pt2;-Pt2;+Pt2;+Pt2;+Pt2;+Pt2;+Pt1;+Pt1;+Pt1;+Pt1];
%     wg=[wgt1;wgt3;wgt3;wgt1;wgt1;wgt2;wgt2;wgt3;wgt3;wgt2;wgt2;wgt3;wgt1;wgt3;wgt3;wgt1];
% elseif nIP==9;
    % Integration point coordinates in local coordinate system
    pt = 0.774596669241483;
    wg = [25 25 25 25 40 40 40 40 64]/81;
    % %     Nodal Coordinates of parent element and node numbering
    ri = pt*[+1;-1;-1;+1;+0;-1;+0;+1;+0]; % Coefficient for r
    si = pt*[+1;+1;-1;-1;+1;+0;-1;+0;+0]; % Coefficient for s
% elseif nIP==4
%     pt = 0.57735026918962576;
%     ri = pt*[+1;-1;-1;+1];
%     si = pt*[+1;+1;-1;-1];
%     wg = [1 1 1 1];
% else
%     pt=0;
%     ri=0;
%     si=0;
%     wg=2;
% end
IPTEMP=zeros(neH,nIP);
for ie=1:neH
  Te  = TEMP(MAPH(ie,1:4)');      % Element nodal temperatures
  for ip=1:nIP
    r = ri(ip);
    s = si(ip);
    N = 0.25*[(1+r)*(1+s);
              (1-r)*(1+s); 
              (1-r)*(1-s); 
              (1+r)*(1-s)];
    Tip            = N'*Te;  % ip current temperature
    IPTEMP(ie,ip)  = Tip;                 % set temperature
  end % for (ip=1:9)
end % for (ie=1:neH)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[ETA]=call_line_search(G,S,alph,dT,T1,dt,T_pre,L_pre,IPDQ,t)
        global LSTOL
       %0000000  Define Line Search Parameters 0000000000000
        it1_max=3;
        alphmina=0;
        alphmaxa=10;
        amp=10;
        LSTOL=0.8;
        kmax=10;
        Nor_fact=G(1);
        k=2;ICO=0;
         while (ICO<2) & (abs(S(k))>LSTOL)
             [alp_minus,S_minus]=Find_min_ratio(alph,S,dT,T1,dt,T_pre,L_pre,Nor_fact,IPDQ,t);
             [Smin_pre,loc]=Find_min_S_rev_dir(S);
             if S_minus>0 | (S_minus> Smin_pre)
                 S_minus=Smin_pre;
                 alp_minus=alph(loc);
             end
             if S_minus<0 %%%% when negative ratio exist
                 [alp_plus,S_plus]=Find_plus_ratio(alph,S,dT,T1,dt,T_pre,L_pre,Nor_fact,IPDQ,t);
                 [Splus_pre,loc]=Find_plus_S_rev_dir(S);
                 if S_plus<0 | (S_plus> Splus_pre)
                     S_plus=Splus_pre;
                     alp_plus=alph(loc);
                 end
                 alpk1=interpol(alph(k-1),alph(k),S(k-1),S(k),0);
                 alpk2=alp_plus+(0.2*(alp_minus-alp_plus));
                 alph(k+1)=max(alpk1,alpk2);
                 if alph(k+1)<alphmina
                     alph(k+1)=alphmina;
                     if ICO==1
                         ICO=2;
                     elseif ICO==0
                         ICO=1;
                     end
                 end
             else %%% LS part 2
                 alphmaxp=max(alph);
                 alph(k+1)=extrapolate(alph(k-1),alph(k),S(k-1),S(k),0);
                 alph3=alph(k+1);
                 if alph(k+1)<0 | alph(k+1)> amp*alphmaxp
                     alph(k+1)=amp*alphmaxp;
                 end
                 if alph(k+1)>alphmaxa & ICO==1
                     ICO=2;
                 end
                 if alph(k+1)>alphmaxa & ICO==0
                     ICO=1;
                 end
             end
             k=k+1;
             G(k)=FunctionTrial(alph(k),dT,T1,dt,T_pre,L_pre,Nor_fact,IPDQ,t);
             S(k)=G(k)/Nor_fact;
         end
ETA=alph(k);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[xx, yy]=Find_min_ratio(x,y,dT,T1,dt,T_pre,L_pre,Nor_fact,IPDQ,t)
    global it1_max
    k1=length(x);
    k=k1;
    while k<=it1_max+k1
        x(k+1)=interpol(x(k-1),x(k),y(k-1),y(k),0);
        y(k+1)=FunctionTrial(x(k+1),dT,T1,dt,T_pre,L_pre,Nor_fact,IPDQ,t);
        k=k+1;
        if y(k)<0
            break;
        end
    end
    xx=x(k);yy=y(k);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function[xx, yy]=Find_plus_ratio(x,y,dT,T1,dt,T_pre,L_pre,Nor_fact,IPPH,IPDQ,DFDT,t)
    global it1_max
    k1=length(x);
    k=k1;
    while k<=it1_max+k1
        x(k+1)=interpol(x(k-1),x(k),y(k-1),y(k),0);
        y(k+1)=FunctionTrial(x(k+1),dT,T1,dt,T_pre,L_pre,Nor_fact,IPDQ,t);
        k=k+1;
        if y(k)>0
            break;
        end
    end
    xx=x(k);yy=y(k);
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   function[y]=FunctionTrial(x,dT,T1,dt,T_pre,L_pre,Nor_fact,IPDQ,t)
   T2    = T1 + (x*dT);
   [DFDT,IPPH]=calculate_DFDT(T2,T_pre);
   [K,C,F,L,Cpc] = assamble_heatPH(T2,IPPH,IPDQ,DFDT,t); % find K C L DLDT at T
    B=F-(C*((T2-T_pre)/dt))-(K*T2)-(L-L_pre)/dt;
    A=K+(C/dt)+(Cpc/dt);
   [A,B]=EssentialBC(A,B);
   P = dT'*B; 
   y=P/Nor_fact;
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[X]=interpol(x1,x2,y1,y2,y)
dy=y2-y1;
dx=x2-x1;
if dy==0
    X=x1;
else
    X=x1+((y-y1)*dx/dy);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[X]=extrapolate(x1,x2,y1,y2,y)
dy=y2-y1;
dx=x2-x1;
if dy==0
    X=x1;
else
    X=x1+((y-y1)*dx/dy);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[X,n1]=Find_min_S_rev_dir(S)
nS=length(S); 
X=S(nS);n1=nS;
for j=1:nS
    n1=nS+1-j;
    if S(n1)<0
        X=S(n1);
        break;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[X,n1]=Find_plus_S_rev_dir(S)
nS=length(S); 
X=S(nS);n1=nS;
for j=1:nS
    n1=nS+1-j;
    if S(n1)>0
        X=S(n1);
        break;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Xmix = mixture_rule(X,nX,f,T,N)
%  Subfunction: Mixture quantity intrepolations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  INPUT:                                                                 %
%    [X] : Material Property Array                                                         %
%    [f] : Phase fractions at t                                           %
%    [T] : Temperature at t                                               %
%    [N] : Exponent (1-> Aritmatic, 0-> Geometric, -1 -> Harmonic)        %
%  OUTPUT:                                                                %
%    [Xmix] : Mixture value                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 if T<=X(1,1) % Make linear extrapolation
    XX2=X(1,3);
    XX1=X(1,2);
  elseif T>=X(nX,1) % Make linear extrapolation
    XX2=X(nX,3);
    XX1=X(nX,2);
  else
    jj=nX;
    while T<X(jj,1), jj=jj-1; end
    xdv = (T-X(jj,1))/(X(jj+1,1)-X(jj,1));
    XX2 = X(jj,3)+(X(jj+1,3)-X(jj,3))*xdv;
    XX1 = X(jj,2)+(X(jj+1,2)-X(jj,2))*xdv;
  end

  % all phase fractions
  ff2 = f;      % Liquid fraction
  ff1 = 1-ff2;  % Solid fraction

  % Mixture Slope
  if ~N % GM: Geometric mean
    Xmix  = exp(ff1*log(XX1)+ff2*log(XX2));
  elseif N==1 % AM: Arithmetic mean
    Xmix  = (ff1*XX1)+(ff2*XX2);
  elseif N==-1 % HM: Harmonic mean
    Xmix  = 1/(ff1/XX1+ff2/XX2);
  else % Other means
    Xmix  = (ff1*XX1^N+ff2*XX2^N)^(1/N);
  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%555
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xvl] = val_at_t(X,nX,T)
%  Interpolation for piecewise linear data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  INPUT:                                                                 %
%    [X] : Matrix with size Nx2                                           %
%    [T] : Temperature                                                    %
%  OUTPUT:                                                                %
%    [xvl] : Linearly interpolated value                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if T<X(1,1) % Make linear extrapolation
%   xvl = X(1,2)-(X(1,1)-T)*(X(2,2)-X(1,2))/(X(2,1)-X(1,1));
xvl=X(1,2);
elseif T>X(nX,1) % Make linear extrapolation
%   xvl = X(nX,2)+(T-X(nX,1))*(X(nX,2)-X(nX-1,2))/(X(nX,1)-X(nX-1,1));
xvl=X(nX,2);
else
  for jj=2:nX
    if T<=X(jj,1) % Make linear interpolation
       xvl = X(jj-1,2)+(T-X(jj-1,1))*(X(jj,2)-X(jj-1,2))/(X(jj,1)-X(jj-1,1));
      break;
    end
  end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Yalcin defined%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [IPDT,IPTEMP] = set_ip_temp(T,IPTEMP)
%  Sets temperature and temperature increment for each IP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  INPUT:
%    [T]      : Nodal temperatures
%    [IPTEMP] : IP temperature data from previous time step
%  OUTPUT:
%    [IPDT]   : IP temperature increment data for current time step
%    [IPTEMP] : IP temperature data for current time step
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%555
    
  
    function[xvl]=phase_interpolation(X,nX,y)
        if (X(1,1)-X(nX,1))>0 %%% descending order
            if y>=X(1,1)
                xvl=X(1,2);
            elseif y<=X(nX,1)
                xvl=X(nX,2);
            else
                for jj=2:nX
                    if y>=X(jj,1) % Make linear interpolation
                        xvl = X(jj-1,2)+(y-X(jj-1,1))*(X(jj,2)-X(jj-1,2))/(X(jj,1)-X(jj-1,1));
                        break;
                    end
                end
            end
        else %%% ascending order
            if y<=X(1,1)
                xvl=X(1,2);
            elseif y<=X(nX,1)
                xvl=X(nX,2);
            else
                for jj=2:nX
                    if y<=X(jj,1) % Make linear interpolation
                        xvl = X(jj-1,2)+(y-X(jj-1,1))*(X(jj,2)-X(jj-1,2))/(X(jj,1)-X(jj-1,1));
                        break;
                    end
                end
            end
        end
        
        function[xout]=Linear_Interpolation(X,nX,yin)
            % nX=size(X,1);
            for i=1:length(yin)
                y=yin(i);
                if (X(nX,1)-X(1,1))>0 %%% descending order
                    if y<=X(1,1)
                        xvl=X(1,2);
                    elseif y>=X(nX,1)
                        xvl=X(nX,2);
                    else
                        for jj=2:nX
                            if (y<=X(jj,1)) && (y>=X(jj-1,1))% Make linear interpolation
                                x1=X(jj-1,1);x2=X(jj,1);y1=X(jj-1,2);y2=X(jj,2);
                                xvl = y1+(y-x1)*(y2-y1)/(x2-x1);
                                break;
                            end
                        end
                    end
                end
                if (X(nX,1)-X(1,1))<0 %%% descending order
                    if y>=X(1,1)
                        xvl=X(1,2);
                    elseif y<=X(nX,1)
                        xvl=X(nX,2);
                    else
                        for jj=2:nX
                            if (y>=X(jj,1)) && (y<=X(jj-1,1))% Make linear interpolation
                                x1=X(jj-1,1);x2=X(jj,1);y1=X(jj-1,2);y2=X(jj,2);
                                xvl = y1+(y-x1)*(y2-y1)/(x2-x1);
                                break;
                            end
                        end
                    end
                end
                xout(i)=xvl;
            end


function[XY_out]=Space_Interpolation(Y,X,XY)
nY=length(Y);
nR=size(XY,1);
% XY_out=zeros(nR,length(X));
for i=1:length(X)
    x=X(i);
    if x==Y(1)
        n1=1;n2=2;fr=0;
    elseif x>Y(nY)
        n1=nY-1;n2=nY;fr=1;
    else
        for j=1:nY-1
            y1=Y(j);
            y2=Y(j+1);
            if (x>=y1) && (x<=y2)% Make linear interpolation
                n1=j;n2=j+1;
                fr=(x-y1)/(y2-y1);
                break;
            end
        end
    end
    XY_out(i,:)=XY(:,n1)+fr*(XY(:,n2)-XY(:,n1));
end
XY_out=XY_out';




