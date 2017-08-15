
function[Nozzle_Tab,No_Noz_Zo,SEC_COOL]=Secondary_Cooling(Wi,cooling_info)
[POS_NOZ,WAT_FLO,NOZ_DIS,CON_ANG,No_Noz_Zo,Spray_rad,ICFM,VEL,DIA]=create_nozzle_para(cooling_info);
Nozzle_Tab=[round(POS_NOZ.*100)./100 WAT_FLO NOZ_DIS' CON_ANG Spray_rad,ICFM,VEL',DIA'];
A=[round(POS_NOZ.*100)./100 NOZ_DIS' CON_ANG];
% save('nozzle.mat','A');

a=POS_NOZ(:,1);
[a,IA,IC] = unique(a);
b=WAT_FLO(:,2);
b=b(IA);
c=NOZ_DIS';
c=c(IA).*1e-3;
d1=CON_ANG(:,1);d2=CON_ANG(:,2);
d1=d1(IA);d2=d2(IA);
area1=pi*(c.*tan(d1.*pi/360)).^2;
imp1=b./(area1.*60);

area2=pi*(c.*tan(d2.*pi/360)).^2;
imp2=b./(area2.*60);

% [a./1000]
% b
% imp1
% imp2
% bar(a,b);
% set(gca,'LineWidth',2','FontSize',20,'Box', 'off','GridLineStyle',':');
% xlabel('Axial direction [mm]','FontSize',20);
% ylabel('Water Flow rate [Liter/min]','FontSize',20);

[Z,X,Q]=Secondary_Cooling_flow(Wi,Nozzle_Tab);
HTC=find_HTC(Z,X,Q);
SEC_COOL=struct(...
    'Z',  {Z},...
    'X',{X},...
    'Q', {Q},...
    'HTC',{HTC});
% contour(Z,X,Q,[0:1:12]);
% plot(Z(51,:),Q(51,:),'r','LineWidth',3);hold on
% set(gca,'LineWidth',2','FontSize',20,'Box', 'off','GridLineStyle',':');
% xlabel('Axial direction [mm]','FontSize',20);
% ylabel('HTC [W/m^2/K]','FontSize',20);
% ylabel('Impingement density [kg/m^2/s]','FontSize',20);

function[HTC]=find_HTC(Z,X,Q)
HTC=zeros(size(Q));
for i=1:size(Q,1)
    for j=1:size(Q,2)
        Qi=Q(i,j);
    %%% C1
%     if Qi<=4
%         HTC(i,j)=312.5*Qi;
%     else
%         HTC(i,j)=1250+75*(Qi-4);
%     end
    %%% C2
%     A=312.5;B=4*A;C=120;
%     if Qi<=4
%         HTC(i,j)=A*Qi;
%     else
%         HTC(i,j)=B+C*(Qi-4);
%     end
    HTC(i,j)=170*Qi;
% HTC(i,j)=392.5.*(Qi.^0.55).*(1-0.0075*T_water);
    end
end


function[Nozzle_Parameters,Nozzle_to_Zone,Q_vs_CA]=Collect_Nozzle_DATA()
Nozzle_Parameters=struct(...
    'CAF',{0.8,0.9,0.8},... % Cone Angle Factor
    'ICF1',{0.5,0.5,0.5},... % Inner Circle Factor 1
    'ICF2',{0.5,0.5,0.5},...  % Inner Circle Factor 2
    'STD_FR',{3.5,2.5,2.0},...
    'STD_CA',{65,60,50},...
    'NOZ_DIA',{1.9,1.6,1.4}); 
Nozzle_to_Zone=[1,1,2,3,3,3]; % 6 zones nozzle number
    
 A{1}=[
         0   24.0000
    0.5000   30.3715
    1.0000   36.4720
    1.5000   42.2805
    2.0000   47.7760
    2.5000   52.9375
    3.0000   57.7440
    3.5000   62.1745
    4.0000   66.2080
    4.5000   69.8235
    5.0000   73.0000
    5.5000   75.7165
    6.0000   77.9520
    6.5000   79.6855
    7.0000   80.8960
    ];

 A{2}=[
 0         0
    1.0000   35.0000
    1.5000   40.0000
    2.5000   57.0000
    3.5000   60.0000
    4.3000   65.0000   ];

 A{3}=[
      0         0
    0.5000   15.0000
    1.5000   36.5000
    2.0000   42.0000
    3.0000   43.0000
    3.3000   51.0000];
Q_vs_CA=A;


function[Z,X,Q]=Secondary_Cooling_flow(Wi,NozzTab)
ner    = 100;        % no.of elements in x direction
nes    = 100;        % no.of elements in y direction
nozz_posS=NozzTab(:,[1 2]);
nozz_angS=NozzTab(:,6);
nozz_disS=NozzTab(:,5);
Q_nozS=NozzTab(:,[3 4]);
spray_rad=NozzTab(:,[8 9]);
ICFM=NozzTab(:,[10 11]);
VEL=NozzTab(:,12);
DIA=NozzTab(:,13);
sec_L=nozz_posS(size(nozz_posS,1),1)+1000;
dz  = 5;
z   = 0:dz:sec_L;  % z - vector
x   = 0:Wi/ner:Wi;      % x - vector
% [Z,X,Q]=assemble_spray_matrix(z,x,nozz_posS,nozz_angS,nozz_disS,Q_nozS);
[Z,X,Q]=assemble_spray_matrix(z,x,nozz_posS,nozz_angS,nozz_disS,Q_nozS,spray_rad,ICFM,VEL,DIA);



function[Z,X,Q]=assemble_spray_matrix(z,x,nozz_posS,nozz_angS,nozz_disS,Q_nozS,spray_rad,ICFM,VEL,DIA)
[Z,X]=meshgrid(z,x);
Q=zeros(size(Z));V=zeros(size(Z));D=zeros(size(Z));
ne=length(nozz_angS);
for i=1:ne
    Qe=[];zmap=[];xmap=[];
    nozz_ang=nozz_angS(i);
    nozz_pos=nozz_posS(i,:);
    nozz_dis=nozz_disS(i);
    Q_noz=Q_nozS(i,:);
    dv=VEL(i);dd=DIA(i);
    Drop_info=[dv dd]; % droplet information
    R2   = spray_rad(i,2);R1   = spray_rad(i,1);
    ICFMi=ICFM(i,:);
    LRTB=find_square_matrix_spray(nozz_pos,R2,z,x);
    zmap=LRTB(1):LRTB(2);
    xmap=LRTB(4):LRTB(3);
    Z1=Z(xmap,zmap);
    X1=X(xmap,zmap);
%   [Qe]=Spray_dist_matrix(nozz_dis,nozz_ang,Z1,X1,nozz_pos,Q_noz,[R1 R2],ICFMi);
 

[Qe]=Spray_uniform_matrix(nozz_dis,nozz_ang,Z1,X1,nozz_pos,Q_noz,[R1 R2],ICFMi);

  
  Q(xmap,zmap)=Q(xmap,zmap)+Qe; 
end


function[LRTB]=find_square_matrix_spray(nozz_pos,R,z,x)
%%% nozzle square position in global matrix
L1=nozz_pos(1)-R; %%% left line
R1=nozz_pos(1)+R;
T1=nozz_pos(2)+R;
B1=nozz_pos(2)-R;
if T1>max(x)
    T1=max(x);
end
if B1<min(x)
    B1=min(x);
end
if L1<min(z)
    L1=min(z);
end

LRTB_coord=[L1 R1 T1 B1];
LRTB=[];cntid2=0;
for i=1:length(LRTB_coord)
    La=LRTB_coord(i);A=[];
    if i<=2
        A=z;        
    else
        A=x;
    end
    if i==1 || i==4
        pos=1;
    else
        pos=2;
    end
    cntid=0;Nid=0;a=[];b=[];
    while cntid==0
        a=A(Nid+1);
        b=A(Nid+2);
        if ((a-La)<=0) && ((b-La)>=0)
            LRTB(cntid2+1)=Nid+pos;
            cntid=cntid+1;cntid2=cntid2+1;
        end
        Nid=Nid+1;
    end
end


function[Qe]=Spray_dist_matrix(nozz_dis,nozz_ang,Z1,X1,nozz_pos,Q_noz,Radius,ICFMi)
nozz_ang=nozz_ang.*pi/180;
R   = nozz_dis*tan(nozz_ang/2);
R1=Radius(1);R2=Radius(2);
Qe=zeros(size(Z1));aa=[];
Ve=zeros(size(Z1));De=zeros(size(Z1));
f1=1.2;f2=0.8;
for row=1:size(Z1,1)
    for col=1:size(Z1,2)
        za=Z1(row,col);
        xa=X1(row,col);
        r=sqrt((nozz_pos(1)-za)^2+(nozz_pos(2)-xa)^2);
        Fract=ICFMi(1)*(Q_noz(2)/Q_noz(1))^ICFMi(2);
        C=Fract*Q_noz(2)/pi/(R1*1e-3)^2/60;
        %         pause
        if r<=R1
            Qe(row,col)=C;
        else
            alpha=C*R1*R2/(R2-R1);
            beta=-(R2+R1)/(R2*R1 );
            w=C*(R2^2-R1^2)/2+(alpha*(2*R2*R1-R1^2-R2^2)/2/R1);
            rho=log(R2/R1)-(((R2^2-R1^2)/2/R1)*(beta+1/R1))+beta*(R2-R1);
            a2=(((1-Fract)*Q_noz(2)/2/pi)-w)/rho;
            a1=alpha+a2*beta;
            Qi=C+a1*(1/r-1/R1)+a2*(1/r^2-1/R1^2);
            if Qi<0
                Qe(row,col)=0;
            else
                Qe(row,col)=C+a1*(1/r-1/R1)+a2*(1/r^2-1/R1^2);
            end
        end
       
    end
end

function[Qe]=Spray_uniform_matrix(nozz_dis,nozz_ang,Z1,X1,nozz_pos,Q_noz,Radius,ICFMi)
nozz_ang=nozz_ang.*pi/180;
R   = nozz_dis*tan(nozz_ang/2);
Qe=zeros(size(Z1));
for row=1:size(Z1,1)
    for col=1:size(Z1,2)
        za=Z1(row,col);
        xa=X1(row,col);
        r=sqrt((nozz_pos(1)-za)^2+(nozz_pos(2)-xa)^2);
        if r<=R
            Qe(row,col)=Q_noz(2)/pi/(R*1e-3)^2/60;
        else
            Qe(row,col)=0;
        end
    end
end


            
function[POS_NOZ,WAT_FLO,NOZ_DIS,CON_ANG,No_Noz_Zo,Spray_rad,ICFM,VEL,DIA]=create_nozzle_para(cooling_info)
%%% create individual nozzle information

Zones={cooling_info.Zones};
Nozzles=convert_cell_matrix({cooling_info.Nozzles});
Zone_pos=convert_cell_matrix({cooling_info.Zone_pos});
WFR=convert_cell_matrix({cooling_info.WFR});
ND=convert_cell_matrix({cooling_info.ND});
CA=convert_cell_matrix({cooling_info.CA});
NNP=find_num_nozz_per_zone(Nozzles);
WFR=WFR./NNP;
[Nozzle_Parameters,Nozzle_to_Zone,Q_vs_CA]=Collect_Nozzle_DATA;

ICFM=[];
POS_NOZ=[];% final nozzle positions
WAT_FLO=[];% Nozzle water flow rate
NOZ_DIS=[];% Nozzle distance
CON_ANG=[];% Cone angle
No_Noz_Zo=[];Spray_rad=[];
ct=0;

for i=1:length(Zones)
    nRows=Nozzles(i,1);
    if nRows>1
        nRow=nRows-1;
    else
        nRow=nRows;
    end
    nCol=Nozzles(i,2);
    axial_1=Zone_pos(i,1);
    axial_2=Zone_pos(i,2);
    lateral=[Zone_pos(i,3) Zone_pos(i,4)];
    dist=(axial_2-axial_1);
    WFRi=WFR(i);NDi=ND(i);CAi=CA(i);
    for j=1:nRows
        if lateral(2)>0
            nCol=2;
        else
            nCol=1;
        end
        for k=1:nCol
            POS_NOZ(ct+1,:)=[axial_1+(j-1)*dist/nRow, lateral(k)];
            WFR2=convert_cell_matrix({Nozzle_Parameters.STD_FR});
            WAT_FLO(ct+1,1)=WFR2(Nozzle_to_Zone(i)); % standard flow rate
            WAT_FLO(ct+1,2)=WFRi;% operating flow rate
            NOZ_DIS(ct+1)=NDi;
            
            % Jet velocity
            NOZ_DIA=convert_cell_matrix({Nozzle_Parameters.NOZ_DIA});
            dn=NOZ_DIA(Nozzle_to_Zone(i));
            VEL(ct+1)=WFRi*1e-3/(0.785398*(dn*1e-3)^2)/60;
            
            % Sauter mean diatemer
            We=1.2*VEL(ct+1)^2*dn*1e-3/72e-3;
            Re=1000*VEL(ct+1)*dn*1e-3/894e-6;
            DIA(ct+1)=dn*1e-3*(3.67*(We^0.5*Re)^-0.259);
            
            Std_WF2=convert_cell_matrix({Nozzle_Parameters.STD_FR});
            Std_WF=Std_WF2(Nozzle_to_Zone(i));
            
            CAF2=convert_cell_matrix({Nozzle_Parameters.CAF});
            CAF=CAF2(Nozzle_to_Zone(i));
            
            ICF11=convert_cell_matrix({Nozzle_Parameters.ICF1});
            ICF1=ICF11(Nozzle_to_Zone(i));
            
            ICF22=convert_cell_matrix({Nozzle_Parameters.ICF2});
            ICF2=ICF22(Nozzle_to_Zone(i));
            
            CA_MAT=Q_vs_CA{(Nozzle_to_Zone(i))};
            CA2=Linear_Interpolation(CA_MAT,WFRi);
            CON_ANG(ct+1,1)=CAi;CON_ANG(ct+1,2)=CA2;
%             [CAF ICF1 ICF2]
               
r1=NDi*tan(CAF*CA2*pi/360);
fa=1.4*(WFRi/Std_WF)^.2;
r2=NDi*tan(fa*CAi*pi/360);

%             Spray_rad(ct+1)=NDi*tan(CAi*pi/360);
            Spray_rad(ct+1,1)=r1;
             Spray_rad(ct+1,2)=r2;
             
             ICFM(ct+1,:)=[ICF1 ICF2];
            ct=ct+1;
        end
    end
    No_Noz_Zo(i)=ct;
end



function[B]=convert_cell_matrix(A)
for i=1:size(A,2)
    B(i,:)=A{i};
end

function[NNP]=find_num_nozz_per_zone(Nozzles)
for i=1:size(Nozzles,1)
    nn=Nozzles(i,:);
    NNP(i,1)=4*nn(1)*nn(2);
%     NNP(i)=4*nn(1)*nn(2);
end
NNP(4)=NNP(4)+2;


 function[xout]=Linear_Interpolation(X,yin)
nX=size(X,1);
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
