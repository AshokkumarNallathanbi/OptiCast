function varargout = Mech_Results_1(varargin)
% MECH_RESULTS_1 MATLAB code for Mech_Results_1.fig
%      MECH_RESULTS_1, by itself, creates a new MECH_RESULTS_1 or raises the existing
%      singleton*.
%
%      H = MECH_RESULTS_1 returns the handle to a new MECH_RESULTS_1 or the handle to
%      the existing singleton*.
%
%      MECH_RESULTS_1('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MECH_RESULTS_1.M with the given input arguments.
%
%      MECH_RESULTS_1('Property','Value',...) creates a new MECH_RESULTS_1 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Mech_Results_1_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Mech_Results_1_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Mech_Results_1

% Last Modified by GUIDE v2.5 13-Mar-2017 22:54:46

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Mech_Results_1_OpeningFcn, ...
                   'gui_OutputFcn',  @Mech_Results_1_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before Mech_Results_1 is made visible.
function Mech_Results_1_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Mech_Results_1 (see VARARGIN)

% Choose default command line output for Mech_Results_1
handles.output = hObject;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(hObject,'Units','pixels');
% % %Replace with the tag of panel to be postiotioned
pos_main_panel=get(hObject,'Position');
% set(handles.Panel_Default,'Visible','On');
% s_panel=get(handles.Panel_Advanced_Mold,'Position');
stretch=pos_main_panel(4)/722;
set(hObject,'Units','normalized');
pos_gui_n=get(hObject,'Position');
b=pos_gui_n(3); h=pos_gui_n(4);
if h>0.85      % 85% der Displayh?he
  h_new=0.85;
  b_new=b*h_new/h *stretch; 
  set(hObject, 'Position',[0.1 0.1 b_new h_new]); 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%% plot emblam
ah=axes('unit','normalized','position',[0.62 0.7 0.2 0.28]);
bg=imread(['Mech_Result_Points.jpg']);imagesc(bg);
set(ah,'handlevisibility','off','visible','off');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Folder's Name
% str1=pwd;                     % Get current working path
% str2='Programs';                % remove programs    
% Fold_Path=strrep(str1,str2,''); % Path of the main folder
% path_sep='\';
% default=load([Fold_Path path_sep 'Settings' path_sep 'Default_Settings.mat']);
% [Folders]=Create_Folder_Struct(Fold_Path,default,path_sep);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

V1=varargin;
V2=V1{1};
Folders=V2;
% Folders.Sub_Fold

% SimF=struct('Sim_Heat',V2{1},'Sim_Mech',V2{2},'Mech_Prog',V2{3});
% Sim_Heat=[Folders.Main_Fold.Sim  Folders.Sub_Fold.Sim_Fold '\'];
% Sim_Mech=[Sim_Heat Folders.Sub_Fold.Mec_Fold '\'];
% Mech_Prog=Folders.Main_Fold.Mec;
% SimF=struct('Sim_Heat',Sim_Heat,'Sim_Mech',Sim_Mech,'Mech_Prog',Folders.Main_Fold.Mec);

addpath(Folders.Main_Fold.Mec);
handles.Folders=Folders;
% handles.SimF=SimF;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Update handles structure
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
guidata(hObject, handles);
ResF_H=[Folders.Main_Fold.Sim Folders.Sub_Fold.Res_Fold '\'];
ResF=[Folders.Main_Fold.Sim Folders.Sub_Fold.Res_Fold '\' Folders.Sub_Fold.Mec_Fold '\'];
[DL,F2,index]=Sim_Folders(handles);
set(handles.Mech_Sim_Folders,'String',F2,'Value',index);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Contour Graph %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[AXIAL,DATA,db_value]=find_axialdis_DB(handles);
handles.AXIAL=AXIAL;
handles.DB_DATA=DATA;
handles.db_value=db_value;
set(handles.Mech_R_DataList,'Value',db_value);
guidata(hObject, handles);
MECH_R_Stress(handles);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Point Graph %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Load Initial Database
% Folders=handles.Folders;
% ResF_H=[Folders.Main_Fold.Sim Folders.Sub_Fold.Res_Fold '\'];
% ResF=[Folders.Main_Fold.Sim Folders.Sub_Fold.Res_Fold '\' Folders.Sub_Fold.Mec_Fold '\'];

PATH_Process=[ResF_H 'Process_Parameters.mat'];
load(PATH_Process);
PATH_DB=[ResF 'DB_M_'];
PATH_DB0=[PATH_DB int2str(0) '.mat'];
load(PATH_DB0);
%%%% Extract Data for mech_point_Graph
nDB=db_value;
[TP_STRESS,Axial_z]=MECH_Plot_Point_Graph(PATH_DB,ner,nes,nDB);
handles.TP_STRESS=TP_STRESS;
handles.Axial_z=Axial_z;
guidata(hObject, handles);
MECH_R_POINT_GRAPH(handles);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% --- Outputs from this function are returned to the command line.
function varargout = Mech_Results_1_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in Mech_R_DataList.
function Mech_R_DataList_Callback(hObject, eventdata, handles)
% hObject    handle to Mech_R_DataList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% get(hObject,'Value')
% Hints: contents = cellstr(get(hObject,'String')) returns Mech_R_DataList contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Mech_R_DataList
handles.db_value=get(hObject,'Value');
guidata(hObject, handles);
set(handles.Mech_R_Dist,'String',handles.AXIAL(handles.db_value,1));
MECH_R_Stress(handles)

% --- Executes during object creation, after setting all properties.
function Mech_R_DataList_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Mech_R_DataList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Mech_R_Plot_Save.
function Mech_R_Plot_Save_Callback(hObject, eventdata, handles)
% hObject    handle to Mech_R_Plot_Save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Folders=handles.Folders;
ResF=[Folders.Main_Fold.Sim Folders.Sub_Fold.Res_Fold '\' Folders.Sub_Fold.Mec_Fold '\'];

File_Sugg=[ResF 'Cont_1.fig'];
[FName, PName]=uiputfile(File_Sugg,'Save Figure');
if isequal(FName,0) || isequal(PName,0)
else
    f = figure;
    copyobj([handles.Mech_R_Plot_Cont],f);
    saveas(f,[PName FName]);
    close(gcf);
end

% --- Executes on button press in Mech_R_Evol_Save.
function Mech_R_Evol_Save_Callback(hObject, eventdata, handles)
% hObject    handle to Mech_R_Evol_Save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Folders=handles.Folders;
ResF=[Folders.Main_Fold.Sim Folders.Sub_Fold.Res_Fold '\' Folders.Sub_Fold.Mec_Fold '\'];

File_Sugg=[ResF 'STR_Evo.fig'];
[FName, PName]=uiputfile(File_Sugg,'Save Figure');
if isequal(FName,0) || isequal(PName,0)
else
    f = figure;
    copyobj([handles.Mech_R_Plot_Profiles],f);
    saveas(f,[PName FName]);
    close(gcf);
end

% --- Executes on selection change in Mech_R_Choices.
function Mech_R_Choices_Callback(hObject, eventdata, handles)
% hObject    handle to Mech_R_Choices (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Mech_R_Choices contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Mech_R_Choices
MECH_R_Stress(handles)
MECH_R_POINT_GRAPH(handles)

% --- Executes during object creation, after setting all properties.
function Mech_R_Choices_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Mech_R_Choices (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function[]=MECH_R_Stress(handles)
%%%%%%%%%%%%% Stress Contours

Folders=handles.Folders;
addpath(Folders.Main_Fold.Mec);
Mech_Prog=Folders.Main_Fold.Mec;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SimF=handles.SimF;
% SimF.Mech_Prog
% Sim_Fold='1015_S235JR_2';
% Fold_Path='G:\HAZELETT CASTER\Stress Analysis\V15\';
% PATH_sim=[Fold_Path 'Simulations\' Sim_Fold '\'];
% PATH_Mech_program=[Fold_Path 'Programs\MECH'];
addpath(Mech_Prog);
% FN='MECH3';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Folders=handles.Folders;
ResF=[Folders.Main_Fold.Sim Folders.Sub_Fold.Res_Fold '\' Folders.Sub_Fold.Mec_Fold '\'];


db_value=handles.db_value; % Required DB
DBi=handles.DB_DATA{handles.AXIAL(1,2)};
PATH_DB=[ResF DBi];
load(PATH_DB);

%nDB=51;
N=10; %%% Number of Contours
DBn=handles.DB_DATA{handles.AXIAL(db_value,2)};
PATH_DB=[ResF DBn];
load(PATH_DB);

R_Choices=get(handles.Mech_R_Choices,'Value');

[STRESS]=MECH_POST_EqnStress(IPD);
% S=STRESS.S_zz;
if R_Choices==1
    S=STRESS.S_eff;
elseif R_Choices==2
    S=STRESS.S_zz;
end

% Cont_Col=handles.Cont_Col;

[ND]=MECH_POST_IP2NODE(XYM,MAPM,S);

Cont_Col=get(handles.Mech_R_Con,'Value');

if Cont_Col==1
    MECH_R_Plot_Contour(ND./1e6,XYM,N,handles);
else
    MECH_R_Plot_Contour_color(XYM,MAPM,ND./1e6,handles);
end

function[]=MECH_R_Plot_Contour(T,XYH,N,handles)
%%%%%%%% Generate Nodal solution contour plot
zh=T;
xh=XYH(:,1);yh=XYH(:,2);
xlin=linspace(min(xh),max(xh),150);
ylin=linspace(min(yh),max(yh),150);
[Xh,Yh]=meshgrid(xlin,ylin);
% F = TriScatteredInterp(xh,yh,zh);
F =scatteredInterpolant(xh,yh,zh);
Zh=F(Xh,Yh);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A1=round(max(max(Zh)));
A2=round(min(min(Zh)));
% v=A2:(A1-A2)/N
% % v=30;
% v=[A2:50:A1];
v=[-500:25:300 300:100:1000];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cla(handles.Mech_R_Plot_Cont);
axes(handles.Mech_R_Plot_Cont);
[C,h]=contour(Xh,Yh,Zh,v,'LineWidth',2.5);hold on
text_handle=clabel(C,h,'Fontsize',12,'LabelSpacing',72*5);
set(gca,'LineWidth',2,'FontSize',12);
xlabel('Width [m]', 'FontSize',12);
ylabel('Thickness [m]', 'FontSize',12);
axis equal
axis tight
colorbar


function[AXIAL,DATA,dvalue]=find_axialdis_DB(handles)
%%%% Creates Data List
% Folder_Name=[handles.Folders.Main_Fold.Sim handles.Folders.Sub_Fold.Res_Fold];

Folders=handles.Folders;
ResF=[Folders.Main_Fold.Sim Folders.Sub_Fold.Res_Fold '\' Folders.Sub_Fold.Mec_Fold '\'];
d=dir([ResF 'DB_M_*.mat']);
DATA = {d.name};
N=size(DATA,2);
AXIAL=zeros(N,2);
a=ResF;
for i=1:N
    b=DATA{i};
    c=strcat(a,b);
    load(c);
    AXIAL(i,1)=round(Axial*100)/100;
end
[AXIAL(:,1),AXIAL(:,2)]=sort(AXIAL(:,1));

set(handles.Mech_R_DataList,'String',AXIAL(:,1));
dvalue=size(AXIAL,1);
set(handles.Mech_R_DataList,'Value',dvalue);
set(handles.Mech_R_Dist,'String',AXIAL(dvalue,1));

% guidata(hObject, handles);


%%%%%%%%%%%%%%%%%  Plot Stress Point Graphs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[EL,IP]=MECH_Termial_Points(ner,nes)
%%%% Terminal Points (1-6) in Half section
%%%% EL-Element, IP- Nearest Integration Point
EL(1)=ner*nes;IP(1)=9;%%%       1. Compression corner
EL(2)=2*ner*nes;IP(2)=9;%%      2. Tension corner
EL(3)=nes;IP(3)=9;      %%%     3. Compression midsurface
EL(4)=2*ner*nes-ner+1;IP(4)=9;%%4. Tension midsurface
EL(5)=1;IP(5)=3;        %%      5. Center core
EL(6)=ner*(nes-1)+1;IP(6)=4;%%  6. Center surface

function[TP_STRESS,Axial_z]=MECH_Plot_Point_Graph(PATH_DB,ner,nes,nDB)
%%%% Stores All database Terminal Point Stresses
%%%%% TP_STRESS - Terminal point stress structure
      %%%%% SEF - Effective Stress [ne,nIP]
      %%%%% SZZ - Axial Stress [ne,nIP]
[EL,IP]=MECH_Termial_Points(ner,nes);
for i=1:nDB
    PATH_DBi=[PATH_DB int2str(i-1) '.mat'];
    load(PATH_DBi);
    [STRESS]=MECH_POST_EqnStress(IPD);
    Szz=STRESS.S_zz;
    Sef=STRESS.S_eff;
    for j=1:length(EL)
    SZZ(i,j)=Szz(EL(j),IP(j));
    SEF(i,j)=Sef(EL(j),IP(j));
    end
    Axial_z(i)=Axial;
end
TP_STRESS=struct('SZZ',SZZ./1e6,'SEF',SEF./1e6);

%%%%%%%%%%%%%%%%%  Plot Stress Point Graphs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[]=MECH_R_POINT_GRAPH(handles)
%%%%% Creates Time evolution of Stresses
R_Choices=get(handles.Mech_R_Choices,'Value');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if R_Choices==1
    S=handles.TP_STRESS.SEF;
    A='EFFECTIVE STRESS';
elseif R_Choices==2
    S=handles.TP_STRESS.SZZ;
     A='AXIAL STRESS';
end
X=handles.Axial_z;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cla(handles.Mech_R_Plot_Profiles);
axes(handles.Mech_R_Plot_Profiles);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
L{1}='-r';L{2}='-b';L{3}=':k';L{4}=':m';L{5}='-.g';L{6}='-.c';
LE{1}='1.Com.Corner';LE{2}='2.Ten.Corner';LE{3}='3.Com.Midsur';
LE{4}='4.Ten.Midsur';LE{5}='5.Mid.Center';LE{6}='6.Mid.Surfac';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:size(S,2)
%     cla(handles.Mech_R_Plot_Profiles);
%     axes(handles.Mech_R_Plot_Profiles);
    plot(X,S(:,i),L{i},'LineWidth',2.5);hold on
end
set(gca,'LineWidth',2,'FontSize',12);
xlabel('Axial Distance [m]', 'FontSize',12);
ylabel('Stress [MPa]', 'FontSize',12);
legend(LE,'Location','NorthWest');
title(A,'FontSize',12);


function[]=MECH_R_Plot_Contour_color(XYM,MAPM,ND,handles)
ne   = size(MAPM,1); % Number of elements
map  = zeros(1,4);
mape = [1 5 9 8;5 2 6 9;9 6 3 7;8 9 7 4];
% mapx = [1 2 3 4 1];
cla(handles.Mech_R_Plot_Cont);
axes(handles.Mech_R_Plot_Cont);
for ie=1:ne
  % part1: r>0 s>0 [1 5 9 8]
  map(1:4) = MAPM(ie,mape(1,1:4));
  patch(XYM(map,1),XYM(map,2),ND(map),'EdgeColor','none');
  % part1: r<0 s>0 [5 2 6 9]
  map(1:4) = MAPM(ie,mape(2,1:4));
  patch(XYM(map,1),XYM(map,2),ND(map),'EdgeColor','none');
  % part1: r<0 s<0 [9 6 3 7]
  map(1:4) = MAPM(ie,mape(3,1:4));
  patch(XYM(map,1),XYM(map,2),ND(map),'EdgeColor','none');
  % part1: r>0 s<0 [8 9 7 4]
  map(1:4) = MAPM(ie,mape(4,1:4));
  patch(XYM(map,1),XYM(map,2),ND(map),'EdgeColor','none');  
  % draw the edge
%   line(XYM(MAPM(ie,mapx),1),XYM(MAPM(ie,mapx),2),'Color','k');  
end
set(gca,'LineWidth',2,'FontSize',12);
xlabel('Width [m]', 'FontSize',12);
ylabel('Thickness [m]', 'FontSize',12);
axis equal
axis tight
colorbar


% --- Executes when selected object is changed in Mech_R_Con_Col.
function Mech_R_Con_Col_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in Mech_R_Con_Col 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
% MECH_R_Stress(handles)



% --- Executes on button press in Mech_R_Con.
function Mech_R_Con_Callback(hObject, eventdata, handles)
% hObject    handle to Mech_R_Con (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Mech_R_Con


% --- Executes when selected object is changed in Mech_R_CONT.
function Mech_R_CONT_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in Mech_R_CONT 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
% if (hObject == handles.Mech_R_Con)
%     Cont_Col=1;
%     set(handles.Mech_R_Col,'Value',0);
% elseif (hObject == handles.Mech_R_Col)
%     Cont_Col=0;
%     set(handles.Mech_R_Con,'Value',0);
% end
% handles.Cont_Col=Cont_Col;
% guidata(hObject, handles);
MECH_R_Stress(handles);


% --- Executes on selection change in Mech_Sim_Folders.
function Mech_Sim_Folders_Callback(hObject, eventdata, handles)
% hObject    handle to Mech_Sim_Folders (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Mech_Sim_Folders contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Mech_Sim_Folders
contents = cellstr(get(hObject,'String'));
Res_Fold= contents{get(hObject,'Value')};
% Cur_Res_Folder=[handles.Folders.Main_Fold.Sim Res_Fold];
handles.Folders.Sub_Fold.Res_Fold=Res_Fold;
guidata(hObject, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Contour Graph %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[AXIAL,DATA,db_value]=find_axialdis_DB(handles);
handles.AXIAL=AXIAL;
handles.DB_DATA=DATA;
handles.db_value=db_value;
set(handles.Mech_R_DataList,'Value',db_value);
guidata(hObject, handles);
MECH_R_Stress(handles);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Point Graph %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Load Initial Database
Folders=handles.Folders;
ResF_H=[Folders.Main_Fold.Sim Folders.Sub_Fold.Res_Fold '\'];
ResF=[Folders.Main_Fold.Sim Folders.Sub_Fold.Res_Fold '\' Folders.Sub_Fold.Mec_Fold '\'];

PATH_Process=[ResF_H 'Process_Parameters.mat'];
load(PATH_Process);
PATH_DB=[ResF 'DB_M_'];
PATH_DB0=[PATH_DB int2str(0) '.mat'];
load(PATH_DB0);
%%%% Extract Data for mech_point_Graph
nDB=db_value;
[TP_STRESS,Axial_z]=MECH_Plot_Point_Graph(PATH_DB,ner,nes,nDB);
handles.TP_STRESS=TP_STRESS;
handles.Axial_z=Axial_z;
guidata(hObject, handles);
MECH_R_POINT_GRAPH(handles);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Folders=handles.Folders;
% Folders.Sub_Fold
% --- Executes during object creation, after setting all properties.
function Mech_Sim_Folders_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Mech_Sim_Folders (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function[DL,F2,index]=Sim_Folders(handles)
Folders=handles.Folders;
d=dir(Folders.Main_Fold.Sim);
str={d.name};
if length(str)>2
    folders=str(3:length(str));
end
index=1;ct=0;
for i=1:length(folders)
    str=[];
    d2=dir([Folders.Main_Fold.Sim folders{i} '\' Folders.Sub_Fold.Mec_Fold '\DB_M_*.mat']);
    str={d2.name};
    str=str(3:length(str));
    if length(str)>1 %&& ct==0
        index=i;
        F2(ct+1)=folders(i);
        ct=ct+1;
    end
    DL(i)=length(str);
end
