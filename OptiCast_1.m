function varargout = OptiCast_1(varargin)
% OptiCast_1 MATLAB code for OptiCast_1.fig
%      OptiCast_1, by itself, creates a new OptiCast_1 or raises the existing
%      singleton*.
%
%      H = OptiCast_1 returns the handle to a new OptiCast_1 or the handle to
%      the existing singleton*.
%
%      OptiCast_1('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in OptiCast_1.M with the given input arguments.
%
%      OptiCast_1('Property','Value',...) creates a new OptiCast_1 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before OptiCast_1_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to OptiCast_1_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help OptiCast_1

% Last Modified by GUIDE v2.5 25-Mar-2017 19:56:39

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @OptiCast_1_OpeningFcn, ...
                   'gui_OutputFcn',  @OptiCast_1_OutputFcn, ...
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

      
        
% --- Executes just before OptiCast_1 is made visible.
function OptiCast_1_OpeningFcn(hObject, eventdata, handles, varargin)
global  pos_main_panel s_panel ss_panel
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to CastingSimSoftware (see VARARGIN)

% Choose CastingSimSoftware command line output for CastingSimSoftware#
handles.output = hObject;

set(handles.Panel_Advanced,'Visible','Off');
set(handles.Panel_Results,'Visible','Off');

set(handles.Panel_Advanced_Mold,'Visible','Off');
set(handles.Panel_Advanced_Second,'Visible','Off');
set(handles.Panel_Advanced_CCM,'Visible','Off');
set(handles.Panel_Advanced_Simulation,'Visible','Off');

set(handles.Panel_Results_Line_Plot,'Visible','Off');
set(handles.Panel_Results_Isotherms,'Visible','Off');

set(handles.Panel_NewGrade,'Visible','Off');

if isdeployed
  try
    shell=actxserver('Shell.Application');   % Zugriff auf alle Programme
    shell.MinimiseAll;                       % Minimiere alles
  end
end


set(hObject,'Units','pixels');

% % %Replace with the tag of panel to be postiotioned
pos_main_panel=get(hObject,'Position');
set(handles.Panel_Default,'Visible','On');
s_panel=get(handles.Panel_Advanced_Mold,'Position');
stretch=pos_main_panel(4)/722;
set(hObject,'Units','normalized');
pos_gui_n=get(hObject,'Position');
b=pos_gui_n(3); h=pos_gui_n(4);

if h>0.85      % 85% der Displayh?he
  h_new=0.85;
  b_new=b*h_new/h *stretch; 
  set(hObject, 'Position',[0.1 0.1 b_new h_new]);
  
end


pos_main_panel=get(handles.Panel_Default,'Position');
s_panel=get(handles.Panel_Advanced_Mold,'Position');


set(handles.Main_Default,'Value',1);
set(handles.Main_Advanced,'Value',0);
set(handles.Main_Results,'Value',0);
set(handles.Main_NewGrade,'Value',0);
set(handles.Main_Default,'BackgroundColor',[239,239,239]/255);
set(handles.Main_Advanced,'BackgroundColor',[214,214,214]/255);
set(handles.Main_Results,'BackgroundColor',[214,214,214]/255);
set(handles.Main_NewGrade,'BackgroundColor',[214,214,214]/255);
set(handles.Advanced_Second,'BackgroundColor',[239,239,239]/255);
set(handles.Advanced_Mold,'BackgroundColor',[214,214,214]/255);
set(handles.Advanced_CCM,'BackgroundColor',[214,214,214]/255);
set(handles.Advanced_Simulation,'BackgroundColor',[214,214,214]/255);
set(handles.Results_Line_Plot,'BackgroundColor',[214,214,214]/255);
set(handles.Results_Isotherms,'BackgroundColor',[214,214,214]/255);
set(handles.Advanced_Second,'Value',0);
set(handles.Advanced_Mold,'Value',0);
set(handles.Advanced_CCM,'Value',0);
set(handles.Advanced_Simulation,'Value',0);
set(handles.Results_Line_Plot,'Value',0);
set(handles.Results_Isotherms,'Value',0);




%%%% Folder Handlings %%%%%%%%%%%%%
str1=pwd;                     % Get current working path
path_sep='\';
%%%%%%%%%%%%%
% str0='CASTSIM_1';str1=[str1 path_sep str0 path_sep];
%%%%%%%%%%%%%
str2='Programs';                % remove programs    
Fold_Path=strrep(str1,str2,''); % Path of the main folder
% Fold_Path=str1;
%path_sep=Fold_Path(length(Fold_Path));

Fold_Path='';
%%%%%%%%%%%% plot emblam
%%%%%%%%%%%% plot emblam
ah=axes('unit','normalized','position',[0.01 0 0.2 0.05]);
bg=imread(['bild_ovgu.jpg']);imagesc(bg);
set(ah,'handlevisibility','off','visible','off');


if isdeployed % Stand-alone-EXE
  [no,result]=system('path');
  root=char(regexpi(result,'Path=(.*?);','tokens','once'))
else
  root=pwd;
end
root
% ah=axes('unit','normalized','position',[0.4 0  0.15 0.05]);
% bg=imread([Fold_Path path_sep str2 path_sep 'bild-swiss-steel.jpg']);imagesc(bg);
% set(ah,'handlevisibility','off','visible','off');

% ah=axes('unit','normalized','position',[0.4 0  0.15 0.05]);
% bg=imread([str1 path_sep 'bild-swiss-steel.jpg']);imagesc(bg);
% set(ah,'handlevisibility','off','visible','off');

% % % % % % % % % % 
%%%%%%%%%%% Load  settings and store
% default=load(['Settings' path_sep 'Default_Settings.mat']);
default=load(['Default_Settings.mat']);
S1=[root path_sep default.Set_Fold_Name];
if isdir(S1)==0
    mkdir(S1);
end
copyfile('Default_Settings.mat',[S1 '\']);

% S2=[root path_sep 'UDM']; 
% if isdir(S2)==0
%     mkdir(S2);
% end
% copyfile('*.mat',[S2 '\']);
% MAT_a=dir([root path_sep '*.mat']);
% MAT_b={MAT_a.name};

% for i=1:length(MAT_b)
%     FN=MAT_b{i};
%     FN2=strrep(FN,'.mat',''); % Path of the main folder
%     FN3=str2num(FN2);
%     if length(FN3)>0
%         delete(FN);
%     end
% end

[Folders]=Create_Folder_Struct(Fold_Path,default,path_sep,root);


handles.Folders=Folders;
handles.path_sep=path_sep;

% Folders.Main_Fold
%%% Load user settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d=dir([Folders.Main_Fold.Set '*.mat']);
str = {d.name};
ndir=size(str,2);
str3=strrep(str,'.mat','');
dset=find(strcmp(str3,Folders.File.SDef));
set(handles.List_Set,'String',str3);
set(handles.List_Set,'Value',dset);
set(handles.User_Set,'String',str3{dset});
% Folders.Main_Fold.Set
% Folders.File.SUse
% Folders.Main_Fold.Set
user=load([Folders.Main_Fold.Set path_sep Folders.File.SUse '.mat']);
user.Std_cooling=default.Std_cooling; %% Standard should be defined in default
user.Std_water_share=default.Std_water_share;
% handles.Set_File=[Fold_Set str3{dset} '.mat'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%% Load Grade dependent parameter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load(['Steel_Grade_Data.mat']);
dvalue=1;
[Grade_depend]=Get_Grade_depend(Steel_reciepes_SC,default.Std_cooling,dvalue);
%%% Set G_Nr and G_Ty in user and default Grade_depend
default.Grade_depend.Grade_Number=Grade_depend.Grade_Number;
default.Grade_depend.Grade_Type=Grade_depend.Grade_Type;
user.Grade_depend.Grade_Number=Grade_depend.Grade_Number;
user.Grade_depend.Grade_Type=Grade_depend.Grade_Type;

handles.Steel_reciepes_SC=Steel_reciepes_SC;
guidata(hObject, handles);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sss.Grade_depend=Grade_depend; % when come to advanced settings, sss
%%% grade dependency is set as same as Grade_depend
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
GNr={Steel_reciepes_SC.Grade_Number};
set(handles.GNList,'String',GNr,'Value',dvalue);
set(handles.GN,'String',Grade_depend.Grade_Number);
set(handles.GT,'String',Grade_depend.Grade_Type);
% handles.Steel_reciepes_SC=Steel_reciepes_SC;
set(handles.mat_on,'Value',1);
set(handles.mat_off,'Value',0);

present=user;
% user_predef=user;
%%%% present is modified according to grade dependency

present.Grade_depend=Grade_depend;

set(handles.Z1,'String',round(present.Grade_depend.WFR_Lm3(1)*10)/10);
set(handles.Z2,'String',round(present.Grade_depend.WFR_Lm3(2)*10)/10);
set(handles.Z3,'String',round(present.Grade_depend.WFR_Lm3(3)*10)/10);
set(handles.CM,'String',present.Grade_depend.Cool_Mode);
set(handles.CS,'String',present.Grade_depend.Casting_speed);
set(handles.CT,'String',Grade_depend.Init_Temp);

% handles.present=present;
% guidata(hObject, handles);
% [SEC_COOL,Nozzle_Tab,No_Noz_Zo,cooling_info]=Updade_Cooling_Info(handles);
% present.SEC_COOL=SEC_COOL;
% present.No_Noz_Zo=No_Noz_Zo;
% present.Nozzle_Tab=Nozzle_Tab;
% present.cooling_info=cooling_info;


handles.present=present;
guidata(hObject, handles);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Sim_Fold=[num2str(Grade_depend.Grade_Number) '_' Grade_depend.Grade_Type];
set(handles.FolderName,'String',Sim_Fold);
handles.Folders.Sub_Fold.Sim_Fold=Sim_Fold;
handles.Folders.Sub_Fold.Res_Fold=Sim_Fold;
handles.Folders.Sub_Fold.Mat_Fold=Sim_Fold;

% Folder_Name=[Fold_Path  default.Sim_Main_Fold '\' Sim_Fold];
% set(handles.FolderName,'String',Sim_Fold);
% present.Sim_Fold=Sim_Fold;
% present.Folder_Name=Folder_Name;

handles.Grade_depend=Grade_depend;
handles.default=default;
handles.user=user;
% handles.present=present;

%%%%% Material data
load([default.Mat_Data_Bank '.mat']);
handles.Mat_Data=Mat_Data;
guidata(hObject, handles);


[Properties]=Load_Grade_Advanced(handles);
set(handles.Liquidus,'String',round(Properties.TL*100)/100);
set(handles.Solidus,'String',round(Properties.TS*100)/100);
set(handles.Latent,'String',round(Properties.Latent/1e3*10)/10);
handles.Properties=Properties;
handles.vprop=1;
plot_material_data(handles);


set(handles.Advanced_GN,'String',Grade_depend.Grade_Number);
set(handles.Advanced_GT,'String',Grade_depend.Grade_Type);
set(handles.Advanced_FN,'String',Sim_Fold);

p=handles.present.Grade_depend;
u=handles.user.Grade_depend;
d=handles.default.Grade_depend;
g=handles.Grade_depend;


Mold_Panel_Reset(handles);

handles.C_DB=0;handles.C_DB2=0;

h_axes_status=handles.axes_status;
set(h_axes_status,'Visible','off');


handles.present_previous=handles.present;

set(handles.Con_1,'Value',1);

% Update handles structure
guidata(hObject, handles);




% --- Outputs from this function are returned to the command line.
function varargout = OptiCast_1_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get CastingSimSoftware command line output from handles structure
varargout{1} = handles.output;


function[SEC_COOL,Nozzle_Tab,No_Noz_Zo,cooling_info]=Updade_Cooling_Info(handles,var)
if var==1
    m=handles.default;
elseif var==2
     m=handles.user;
else
     m=handles.present;
end
WFR_Lm3=handles.present.Grade_depend.WFR_Lm3;
WFR_Lm6=m.Water_share_matrix*WFR_Lm3';
WFR_Lmin6=WFR_Lm6.*m.Grade_depend.Casting_speed;
cooling_info=m.cooling_info;
cooling_info(1).WFR=WFR_Lmin6(1);cooling_info(2).WFR=WFR_Lmin6(2);
cooling_info(3).WFR=WFR_Lmin6(3);cooling_info(4).WFR=WFR_Lmin6(4);
cooling_info(5).WFR=WFR_Lmin6(5);cooling_info(6).WFR=WFR_Lmin6(6);


addpath(handles.Folders.Main_Fold.Set);
[Nozzle_Tab,No_Noz_Zo,SEC_COOL]=Secondary_Cooling(m.Wi,m.cooling_info);
rmpath(handles.Folders.Main_Fold.Set);



% --- Executes on button press in Main_Default. Window_default
function Main_Default_Callback(hObject, eventdata, handles)
% hObject    handle to Main_Default (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Main_Default
global  pos_main_panel
set(handles.Panel_Default,'Position',pos_main_panel);
set(handles.Panel_Advanced,'Visible','Off');drawnow;
set(handles.Panel_Results,'Visible','Off');drawnow;
set(handles.Panel_NewGrade,'Visible','Off');drawnow;
set(handles.Panel_Default,'Visible','On');
set(hObject,'Value',1);
set(hObject,'BackgroundColor',[239,239,239]/255);
set(handles.Main_Advanced,'BackgroundColor',[214,214,214]/255);
set(handles.Main_Results,'BackgroundColor',[214,214,214]/255);
set(handles.Main_NewGrade,'BackgroundColor',[214,214,214]/255);
set(handles.Panel_Advanced_Mold,'Visible','Off');
set(handles.Panel_Advanced_Second,'Visible','Off');
set(handles.Panel_Advanced_CCM,'Visible','Off');
set(handles.Panel_Advanced_Simulation,'Visible','Off');
set(handles.Panel_Results_Line_Plot,'Visible','Off');
set(handles.Panel_Results_Isotherms,'Visible','Off');
% set(handles.Panel_Results_Report,'Visible','Off');
set(handles.Main_Advanced,'Value',0);
set(handles.Main_Results,'Value',0);
set(handles.Main_NewGrade,'Value',0);

load([handles.Folders.Main_Fold.Grd  handles.Folders.File.Grd  '.mat']);
% handles.Steel_reciepes_SC=Steel_reciepes_SC;
GNr={Steel_reciepes_SC.Grade_Number};
Grade_Number=handles.present.Grade_depend.Grade_Number;
dvalue=find(cell2mat(GNr)==Grade_Number);
set(handles.GNList,'String',GNr,'Value',dvalue);
[Grade_depend]=Get_Grade_depend(Steel_reciepes_SC,handles.default.Std_cooling,dvalue);
mat_on=get(handles.mat_on,'Value');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
handles.present=handles.user;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
guidata(hObject, handles);
if mat_on==0
    handles.user.Grade_depend.Grade_Number=Grade_depend.Grade_Number;
    handles.user.Grade_depend.Grade_Type=Grade_depend.Grade_Type;
    handles.present.Grade_depend.Grade_Number=Grade_depend.Grade_Number;
    handles.present.Grade_depend.Grade_Type=Grade_depend.Grade_Type;
    handles.default.Grade_depend.Grade_Number=Grade_depend.Grade_Number;
    handles.default.Grade_depend.Grade_Type=Grade_depend.Grade_Type;
else
     handles.present.Grade_depend=Grade_depend;
end

handles.Grade_depend=Grade_depend;
guidata(hObject, handles);
Sim_Fold=get(handles.FolderName,'String');
handles.Folders.Sub_Fold.Sim_Fold=Sim_Fold;
handles.Folders.Sub_Fold.Res_Fold=Sim_Fold;
% a=handles.present;
set(handles.RunSim,'Value',0);

% [SEC_COOL,Nozzle_Tab,No_Noz_Zo,cooling_info]=Updade_Cooling_Info(handles);
% handles.present.SEC_COOL=SEC_COOL;
% handles.present.No_Noz_Zo=No_Noz_Zo;
% handles.present.Nozzle_Tab=Nozzle_Tab;
% handles.present.cooling_info=cooling_info;
%%% Load user settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d=dir([handles.Folders.Main_Fold.Set '*.mat']);
str = {d.name};
ndir=size(str,2);
str3=strrep(str,'.mat','');
dset=find(strcmp(str3,handles.Folders.File.SDef));
set(handles.List_Set,'String',str3);
set(handles.List_Set,'Value',dset);
set(handles.User_Set,'String',str3{dset});


if isinteger(handles.C_DB)==1
    handles.C_DB=0;
end
guidata(hObject, handles);


% --- Executes on button press in Main_Advanced.
function Main_Advanced_Callback(hObject, eventdata, handles)
% hObject    handle to Main_Advanced (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Main_Advanced
global  pos_main_panel s_panel
set(handles.Panel_Advanced,'Position',pos_main_panel);

set(handles.Panel_Default,'Visible','Off');drawnow;
set(handles.Panel_Results,'Visible','Off');drawnow;
set(handles.Panel_NewGrade,'Visible','Off');drawnow;
set(handles.Panel_Advanced,'Visible','On');
set(hObject,'Value',1);
set(hObject,'BackgroundColor',[239,239,239]/255);
set(handles.Main_Default,'BackgroundColor',[214,214,214]/255);
set(handles.Main_Results,'BackgroundColor',[214,214,214]/255);
set(handles.Main_NewGrade,'BackgroundColor',[214,214,214]/255);

set(handles.Panel_Advanced_Mold,'Position',s_panel);
set(handles.Panel_Advanced_Second,'Visible','Off');
set(handles.Panel_Advanced_CCM,'Visible','Off');
set(handles.Panel_Advanced_Simulation,'Visible','Off');
set(handles.Panel_Advanced_Mold,'Visible','On');

set(handles.Panel_Results_Line_Plot,'Visible','Off');
set(handles.Panel_Results_Isotherms,'Visible','Off');
% set(handles.Panel_Results_Report,'Visible','Off');

set(handles.Advanced_Mold,'Value',1);
set(handles.Advanced_Second,'Value',0);
set(handles.Advanced_CCM,'Value',0);
set(handles.Advanced_Simulation,'Value',0);
set(handles.Advanced_Mold,'BackgroundColor',[239,239,239]/255);
set(handles.Advanced_Second,'BackgroundColor',[214,214,214]/255);
set(handles.Advanced_CCM,'BackgroundColor',[214,214,214]/255);
set(handles.Advanced_Simulation,'BackgroundColor',[214,214,214]/255);
set(handles.Main_Default,'Value',0);
set(handles.Main_Advanced,'Value',1);
set(handles.Main_Results,'Value',0);
set(handles.Main_NewGrade,'Value',0);


handles.default.Grade_depend=handles.Grade_depend;
handles.present=handles.user;
set(handles.Advanced_GN,'String',handles.present.Grade_depend.Grade_Number);
set(handles.Advanced_GT,'String',handles.present.Grade_depend.Grade_Type);
set(handles.Advanced_FN,'String',handles.Folders.Sub_Fold.Sim_Fold);
% Set_File=get(handles.User_Set,'String')
set(handles.Advanced_Set,'String',handles.Folders.File.SUse);
% handles.Set_File=Set_File;
% Set_File=[handles.Folders.Main_Fold.Set  handles.Folders.Sub_Fold.Set '.mat']

guidata(hObject,handles);


[SEC_COOL,Nozzle_Tab,No_Noz_Zo,cooling_info]=Updade_Cooling_Info(handles,3);
handles.present.SEC_COOL=SEC_COOL;
handles.present.No_Noz_Zo=No_Noz_Zo;
handles.present.Nozzle_Tab=Nozzle_Tab;
handles.present.cooling_info=cooling_info;
[SEC_COOL,Nozzle_Tab,No_Noz_Zo,cooling_info]=Updade_Cooling_Info(handles,1);
handles.default.SEC_COOL=SEC_COOL;
handles.default.No_Noz_Zo=No_Noz_Zo;
handles.default.Nozzle_Tab=Nozzle_Tab;
handles.default.cooling_info=cooling_info;
guidata(hObject,handles);


set(handles.MCS_d,'String',handles.default.Grade_depend.Casting_speed);
set(handles.MCT_d,'String',handles.default.Grade_depend.Init_Temp);
set(handles.MCWT_d,'String',handles.default.Tw);
set(handles.MAMT_d,'String',handles.default.Tatm);
set(handles.MMWT_d,'String',handles.default.mold_wall);
set(handles.MDT_d,'String',handles.default.water_temp_inc);
set(handles.MMWF_d,'String',handles.default.Qwater);
set(handles.MML_d,'String',handles.default.Mol_Len);
%%% Modified
Mold_Panel_Reset(handles);
cla(handles.Mold_Plot);
axes(handles.Mold_Plot);

[MHF]=mold_heat_flux(handles.present);
plot(MHF(:,1),MHF(:,2),'-b','LineWidth',3);hold on
set(gca,'FontSize',12,'LineWidth',2.0,'Box','off');
ylabel(' Heat Flux [MW/m^2]','FontSize',13);
xlabel(' Distance from Meniscus [mm]','FontSize',13);
guidata(hObject,handles);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isinteger(handles.C_DB)==1
    handles.C_DB=0;
end
set(handles.RunSim,'Value',0);
set(handles.Advanced_Run,'Value',0);
if isinteger(handles.C_DB2)==1
    handles.C_DB2=0;
end
guidata(hObject, handles);


% --- Executes on button press in Main_Results.
function Main_Results_Callback(hObject, eventdata, handles)
% hObject    handle to Main_Results (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Main_Results
global  pos_main_panel s_panel
set(handles.Panel_Results,'Position',pos_main_panel);
set(handles.Panel_Default,'Visible','Off');drawnow;
set(handles.Panel_Advanced,'Visible','Off');drawnow;
set(handles.Panel_NewGrade,'Visible','Off');drawnow;
set(handles.Panel_Results,'Visible','On');
set(hObject,'Value',1);
set(hObject,'BackgroundColor',[239,239,239]/255);
set(handles.Main_Default,'BackgroundColor',[214,214,214]/255);
set(handles.Main_Advanced,'BackgroundColor',[214,214,214]/255);
set(handles.Main_NewGrade,'BackgroundColor',[214,214,214]/255);

set(handles.Panel_Advanced_Second,'Visible','Off');
set(handles.Panel_Advanced_Mold,'Visible','Off');
set(handles.Panel_Advanced_CCM,'Visible','Off');
set(handles.Panel_Advanced_Simulation,'Visible','Off');

% set(handles.Panel_Results_Line_Plot,'Position',s_panel);
% set(handles.Panel_Results_Isotherms,'Visible','Off');
% set(handles.Panel_Results_Report,'Visible','Off');

set(handles.Panel_Results_Line_Plot,'Visible','On');
set(handles.Panel_Results_Isotherms,'Visible','Off');
% set(handles.Panel_Results_Report,'Visible','Off');

set(handles.Results_Line_Plot,'Value',1);
set(handles.Results_Isotherms,'Value',0);
set(handles.Results_Report,'Value',0);

set(handles.Results_Line_Plot,'BackgroundColor',[239,239,239]/255);
set(handles.Results_Isotherms,'BackgroundColor',[214,214,214]/255);
set(handles.Results_Report,'BackgroundColor',[214,214,214]/255);

set(handles.Main_Default,'Value',0);
set(handles.Main_Advanced,'Value',0);
set(handles.Main_Results,'Value',1);
set(handles.Main_NewGrade,'Value',0);

% Folder_Name=handles.present.Folder_Name;
% Folder_Name=[handles.Fold_Path  handles.default.Sim_Main_Fold];
set(handles.RunSim,'Value',0);
set(handles.Advanced_Run,'Value',0);

Sim_Main_Fold=handles.Folders.Main_Fold.Sim; 
d=dir(Sim_Main_Fold);

str={d.name};
if length(str)>2
    folders=str(3:length(str));
end
index=1;ct=0;
for i=1:length(folders)
    str=[];
    d2=dir([Sim_Main_Fold folders{i}]);
    str={d2.name};
    str=str(3:length(str));
    if length(str)>2 && ct==0
        index=i;ct=1;
    end
    DL(i)=length(str);
end

d2=dir([Sim_Main_Fold handles.Folders.Sub_Fold.Sim_Fold]);
if length({d2.name})>=2
    dset=find(strcmp(folders,handles.Folders.Sub_Fold.Sim_Fold));
    handles.Folders.Sub_Fold.Res_Fold=handles.Folders.Sub_Fold.Sim_Fold;
else
    dset=index;
    handles.Folders.Sub_Fold.Res_Fold=folders{index};
end
set(handles.Results_Folders,'String',folders,'Value',dset);

curr_res_fold=[Sim_Main_Fold handles.Folders.Sub_Fold.Res_Fold];
Parameters=load([curr_res_fold handles.path_sep 'Process_Parameters.mat']);
handles.WSM=handles.default.Water_share_matrix;
Display_Results_Parameters(Parameters,handles);

% handles.Current_Results_Folder=Current_Results_Folder;
handles.Parameters=Parameters;
Line_Plots(curr_res_fold,handles);
set(handles.AC_Color_No,'Value',1);
guidata(hObject, handles);
% [AXIAL,DATA]=find_axialdis_DB(handles);
% handles.AXIAL=AXIAL;
% handles.DB=DATA;
% dvalue=size(AXIAL,1);
% set(handles.Results_AxialDist,'String',AXIAL(dvalue,1));
% set(handles.Results_DataList,'Value',dvalue);
% dbindex=DATA{AXIAL(dvalue,2)};handles.dbindex=dbindex;
% guidata(hObject, handles);
% Temp_contour(handles);
% handles.daxial=AXIAL(dvalue,1);
% set(handles.Results_Midplane,'Value',1);
% guidata(hObject, handles);
% 
% Plot_Axial_Contours(handles)

% --- Executes on button press in Main_NewGrade.
function Main_NewGrade_Callback(hObject, eventdata, handles)
% hObject    handle to Main_NewGrade (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Main_NewGrade
global  pos_main_panel s_panel
set(handles.Panel_NewGrade,'Position',pos_main_panel);
set(handles.Panel_Default,'Visible','Off');drawnow;
set(handles.Panel_Advanced,'Visible','Off');drawnow;
set(handles.Panel_Results,'Visible','Off');drawnow;
set(handles.Panel_NewGrade,'Visible','On');
set(hObject,'Value',1);
set(hObject,'BackgroundColor',[239,239,239]/255);
set(handles.Main_Default,'BackgroundColor',[214,214,214]/255);
set(handles.Main_Advanced,'BackgroundColor',[214,214,214]/255);
set(handles.Main_Results,'BackgroundColor',[214,214,214]/255);

set(handles.Panel_Advanced_Second,'Visible','Off');
set(handles.Panel_Advanced_Mold,'Visible','Off');
set(handles.Panel_Advanced_CCM,'Visible','Off');
set(handles.Panel_Advanced_Simulation,'Visible','Off');

% set(handles.Panel_Results_Line_Plot,'Position',s_panel);
% set(handles.Panel_Results_Isotherms,'Visible','Off');
% set(handles.Panel_Results_Report,'Visible','Off');

set(handles.Panel_Results_Line_Plot,'Visible','Off');
set(handles.Panel_Results_Isotherms,'Visible','Off');
% set(handles.Panel_Results_Report,'Visible','Off');

set(handles.Results_Line_Plot,'Value',0);
set(handles.Results_Isotherms,'Value',0);
set(handles.Results_Report,'Value',0);

set(handles.Main_Default,'Value',0);
set(handles.Main_Advanced,'Value',0);
set(handles.Main_Results,'Value',0);
set(handles.Main_NewGrade,'Value',1);

set(handles.RunSim,'Value',0);
set(handles.Advanced_Run,'Value',0);

Grade_depend=handles.Grade_depend;
set(handles.NewGrade_GN,'String',Grade_depend.Grade_Number);
set(handles.NewGrade_GT,'String',Grade_depend.Grade_Type);
NewGrade_FN=[num2str(Grade_depend.Grade_Number) '_' num2str(Grade_depend.Grade_Type)];
set(handles.NewGrade_FolderName,'String',NewGrade_FN);
% handles.Steel_reciepes_SC=Steel_reciepes_SC;
% set(handles.mat_on,'Value',1);
% set(handles.mat_off,'Value',0);

handles.New_Grade=Grade_depend;
set(handles.NewGrade_Zone_1,'String',Grade_depend.WFR_Lm3(1));
set(handles.NewGrade_Zone_2,'String',Grade_depend.WFR_Lm3(2));
set(handles.NewGrade_Zone_3,'String',Grade_depend.WFR_Lm3(3));
set(handles.NewGrade_CoolMode,'Value',Grade_depend.Cool_Mode_Value);
set(handles.NewGrade_CS,'String',Grade_depend.Casting_speed);
set(handles.NewGrade_CT,'String',Grade_depend.Init_Temp);


guidata(hObject, handles);



% --- Executes on button press in Advanced_Mold.
function Advanced_Mold_Callback(hObject, eventdata, handles)
% hObject    handle to Advanced_Mold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global  pos_main_panel s_panel
% set(handles.Panel_Advanced_Mold,'Position',s_panel);
set(handles.Panel_Advanced_Second,'Visible','Off');
set(handles.Panel_Advanced_CCM,'Visible','Off');
set(handles.Panel_Advanced_Simulation,'Visible','Off');
set(handles.Panel_Advanced_Mold,'Visible','On');

set(handles.Panel_Results_Line_Plot,'Visible','Off');
set(handles.Panel_Results_Isotherms,'Visible','Off');
% set(handles.Panel_Results_Report,'Visible','Off');

set(handles.Panel_Advanced,'Position',pos_main_panel);
set(handles.Panel_Default,'Visible','Off');drawnow;
set(handles.Panel_Results,'Visible','Off');drawnow;
set(handles.Panel_NewGrade,'Visible','Off');drawnow;
set(handles.Panel_Advanced,'Visible','On');
set(hObject,'Value',1);
set(handles.Advanced_Second,'Value',0);
set(handles.Advanced_CCM,'Value',0);
set(handles.Advanced_Simulation,'Value',0);
set(hObject,'BackgroundColor',[239,239,239]/255);
set(handles.Advanced_Second,'BackgroundColor',[214,214,214]/255);
set(handles.Advanced_CCM,'BackgroundColor',[214,214,214]/255);
set(handles.Advanced_Simulation,'BackgroundColor',[214,214,214]/255);

handles.default.Grade_depend=handles.Grade_depend;
guidata(hObject, handles);

set(handles.MCS_d,'String',handles.default.Grade_depend.Casting_speed);
set(handles.MCT_d,'String',handles.default.Grade_depend.Init_Temp);
set(handles.MCWT_d,'String',handles.default.Tw);
set(handles.MAMT_d,'String',handles.default.Tatm);
set(handles.MMWT_d,'String',handles.default.mold_wall);
set(handles.MDT_d,'String',handles.default.water_temp_inc);
set(handles.MMWF_d,'String',handles.default.Qwater);
set(handles.MML_d,'String',handles.default.Mol_Len);
%%% Modified
% handles.present.index
Mold_Panel_Reset(handles);
cla(handles.Mold_Plot);
axes(handles.Mold_Plot);
plot(handles.present.MHF(:,1),handles.present.MHF(:,2),'-b','LineWidth',3);hold on
set(gca,'FontSize',12,'LineWidth',2.0,'Box','off');
ylabel(' Heat Flux [MW/m^2]','FontSize',13);
xlabel(' Distance from Meniscus [mm]','FontSize',13);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes on button press in s2.
function Advanced_Second_Callback(hObject, eventdata, handles)
% hObject    handle to s2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of s2
global  pos_main_panel s_panel
% set(handles.Panel_Advanced_Second,'Position',s_panel);
set(handles.Panel_Advanced_Second,'Visible','On');drawnow;
set(handles.Panel_Advanced_Mold,'Visible','Off');
set(handles.Panel_Advanced_CCM,'Visible','Off');
set(handles.Panel_Advanced_Simulation,'Visible','Off');

set(handles.Panel_Results_Line_Plot,'Visible','Off');
set(handles.Panel_Results_Isotherms,'Visible','Off');
% set(handles.Panel_Results_Report,'Visible','Off');

set(handles.Panel_Advanced,'Position',pos_main_panel);
set(handles.Panel_Default,'Visible','Off');drawnow;
set(handles.Panel_Advanced,'Visible','On');drawnow;
set(handles.Panel_Results,'Visible','Off');drawnow;
set(handles.Panel_NewGrade,'Visible','Off');drawnow;
set(hObject,'Value',1);
set(handles.Advanced_Mold,'Value',0);
set(handles.Advanced_CCM,'Value',0);
set(handles.Advanced_Simulation,'Value',0);
set(hObject,'BackgroundColor',[239,239,239]/255);
set(handles.Advanced_Mold,'BackgroundColor',[214,214,214]/255);
set(handles.Advanced_CCM,'BackgroundColor',[214,214,214]/255);
set(handles.Advanced_Simulation,'BackgroundColor',[214,214,214]/255);

drawnow;
set(handles.SC_CM_d,'String',handles.Grade_depend.Cool_Mode);
set(handles.SC_CM_u,'Value',handles.present.Grade_depend.Cool_Mode_Value);
set(handles.SC_CS_d,'String',handles.default.Grade_depend.Casting_speed);
set(handles.SC_CS_u,'String',handles.present.Grade_depend.Casting_speed);

set(handles.SC_Z1_d,'String',handles.default.Grade_depend.WFR_Lm3(1));
set(handles.SC_Z2_d,'String',handles.default.Grade_depend.WFR_Lm3(2));
set(handles.SC_Z3_d,'String',handles.default.Grade_depend.WFR_Lm3(3));
set(handles.SC_Z1L_d,'String',handles.default.Grade_depend.WFR_Lm3(1)*handles.default.Grade_depend.Casting_speed);
set(handles.SC_Z2L_d,'String',handles.default.Grade_depend.WFR_Lm3(2)*handles.default.Grade_depend.Casting_speed);
set(handles.SC_Z3L_d,'String',handles.default.Grade_depend.WFR_Lm3(3)*handles.default.Grade_depend.Casting_speed);
b=handles.present.Grade_depend;
a=handles.present.Grade_depend.Cool_Mode_Value;
if a==5
    set(handles.SC_Z1_u,'Style','edit','BackgroundColor',[1 1 1]);
    set(handles.SC_Z2_u,'Style','edit','BackgroundColor',[1 1 1]);
    set(handles.SC_Z3_u,'Style','edit','BackgroundColor',[1 1 1]);
    
    set(handles.SC_Z1L_u,'Style','edit','BackgroundColor',[1 1 1]);
    set(handles.SC_Z2L_u,'Style','edit','BackgroundColor',[1 1 1]);
    set(handles.SC_Z3L_u,'Style','edit','BackgroundColor',[1 1 1]);
else
    set(handles.SC_Z1_u,'Style','text','BackgroundColor',0.941*[1 1 1]);
    set(handles.SC_Z2_u,'Style','text','BackgroundColor',0.941*[1 1 1]);
    set(handles.SC_Z3_u,'Style','text','BackgroundColor',0.941*[1 1 1]);
    
    set(handles.SC_Z1L_u,'Style','text','BackgroundColor',0.941*[1 1 1]);
    set(handles.SC_Z2L_u,'Style','text','BackgroundColor',0.941*[1 1 1]);
    set(handles.SC_Z3L_u,'Style','text','BackgroundColor',0.941*[1 1 1]);
end
set(handles.SC_Z1_u,'String',handles.present.Grade_depend.WFR_Lm3(1));
set(handles.SC_Z2_u,'String',handles.present.Grade_depend.WFR_Lm3(2));
set(handles.SC_Z3_u,'String',handles.present.Grade_depend.WFR_Lm3(3));
set(handles.SC_Z1L_u,'String',handles.present.Grade_depend.WFR_Lm3(1)*handles.present.Grade_depend.Casting_speed);
set(handles.SC_Z2L_u,'String',handles.present.Grade_depend.WFR_Lm3(2)*handles.present.Grade_depend.Casting_speed);
set(handles.SC_Z3L_u,'String',handles.present.Grade_depend.WFR_Lm3(3)*handles.present.Grade_depend.Casting_speed);


cla(handles.Secondary_plot);
axes(handles.Secondary_plot);
set(handles.Secondary_Flow_or_HTC,'Value',1);
guidata(hObject, handles);

plot_Secondary_Cooling(handles)



% --- Executes on button press in s3.
function Advanced_CCM_Callback(hObject, eventdata, handles)
% hObject    handle to s3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of s3
global  pos_main_panel s_panel ss_panel
% set(handles.spanel3,'Position',s_panel);
set(handles.Panel_Advanced_Second,'Visible','Off');
set(handles.Panel_Advanced_Mold,'Visible','Off');
set(handles.Panel_Advanced_CCM,'Visible','On');
set(handles.Panel_Advanced_Simulation,'Visible','Off');

set(handles.Panel_Results_Line_Plot,'Visible','Off');
set(handles.Panel_Results_Isotherms,'Visible','Off');
% set(handles.Panel_Results_Report,'Visible','Off');

set(handles.Panel_Advanced,'Position',pos_main_panel);
set(handles.Panel_Default,'Visible','Off');drawnow;
set(handles.Panel_Advanced,'Visible','On');
set(handles.Panel_Results,'Visible','Off');drawnow;
set(handles.Panel_NewGrade,'Visible','Off');drawnow;
set(hObject,'Value',1);
set(handles.Advanced_Mold,'Value',0);
set(handles.Advanced_Second,'Value',0);
set(handles.Advanced_Simulation,'Value',0);
set(hObject,'BackgroundColor',[239,239,239]/255);
set(handles.Advanced_Mold,'BackgroundColor',[214,214,214]/255);
set(handles.Advanced_Second,'BackgroundColor',[214,214,214]/255);
set(handles.Advanced_Simulation,'BackgroundColor',[214,214,214]/255);

[SEC_COOL,Nozzle_Tab,No_Noz_Zo,cooling_info]=Updade_Cooling_Info(handles,3);
handles.present.SEC_COOL=SEC_COOL;
handles.present.No_Noz_Zo=No_Noz_Zo;
handles.present.Nozzle_Tab=Nozzle_Tab;
handles.present.cooling_info=cooling_info;
P=handles.present.Grade_depend;
guidata(hObject, handles);

[SEC_COOL,Nozzle_Tab,No_Noz_Zo,cooling_info]=Updade_Cooling_Info(handles,1);
handles.default.SEC_COOL=SEC_COOL;
handles.default.No_Noz_Zo=No_Noz_Zo;
handles.default.Nozzle_Tab=Nozzle_Tab;
handles.default.cooling_info=cooling_info;
guidata(hObject, handles);
D=handles.default.Grade_depend;

Nozzle_Tab=handles.default.Nozzle_Tab;
POS_NOZ=Nozzle_Tab(:,[1 2]);
WAT_FLO=Nozzle_Tab(:,[4]);
NOZ_DIS=Nozzle_Tab(:,[5]);
CON_ANG=Nozzle_Tab(:,[7]);
No_Noz_Zo=handles.default.No_Noz_Zo;
Spray_rad=Nozzle_Tab(:,[9]);
[D_DATA_struct,D_DATA_table]=make_CCM_Data(POS_NOZ,WAT_FLO,NOZ_DIS,CON_ANG,No_Noz_Zo,Spray_rad);
set(handles.CCM_Default,'Data',D_DATA_table);

Nozzle_Tab=handles.present.Nozzle_Tab;
POS_NOZ=Nozzle_Tab(:,[1 2]);
WAT_FLO=Nozzle_Tab(:,[4]);
NOZ_DIS=Nozzle_Tab(:,[5]);
CON_ANG=Nozzle_Tab(:,[7]);
No_Noz_Zo=handles.present.No_Noz_Zo;
Spray_rad=Nozzle_Tab(:,[9]);
[M_DATA_struct,M_DATA_table]=make_CCM_Data(POS_NOZ,WAT_FLO,NOZ_DIS,CON_ANG,No_Noz_Zo,Spray_rad);
set(handles.CCM_Modified,'Data',M_DATA_table);
guidata(hObject, handles);


% --- Executes on button press in s4.
function Advanced_Simulation_Callback(hObject, eventdata, handles)
% hObject    handle to s4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of s4
global  pos_main_panel s_panel
set(handles.Panel_Advanced_Second,'Visible','Off');
set(handles.Panel_Advanced_Mold,'Visible','Off');
set(handles.Panel_Advanced_CCM,'Visible','Off');
set(handles.Panel_Advanced_Simulation,'Visible','On');


set(handles.Panel_Results_Line_Plot,'Visible','Off');
set(handles.Panel_Results_Isotherms,'Visible','Off');
% set(handles.Panel_Results_Report,'Visible','Off');

set(handles.Panel_Advanced,'Position',pos_main_panel);
set(handles.Panel_Default,'Visible','Off');drawnow;
set(handles.Panel_Advanced,'Visible','On');
set(handles.Panel_Results,'Visible','Off');drawnow;
set(handles.Panel_NewGrade,'Visible','Off');drawnow;
set(hObject,'Value',1);
set(handles.Advanced_Mold,'Value',0);
set(handles.Advanced_Second,'Value',0);
set(handles.Advanced_CCM,'Value',0);
set(hObject,'BackgroundColor',[239,239,239]/255);
set(handles.Advanced_Mold,'BackgroundColor',[214,214,214]/255);
set(handles.Advanced_Second,'BackgroundColor',[214,214,214]/255);
set(handles.Advanced_CCM,'BackgroundColor',[214,214,214]/255);




set(handles.SWi_d,'String',handles.default.Wi);
set(handles.STh_d,'String',handles.default.Th);
set(handles.Sner_d,'String',handles.default.ner);
set(handles.Snes_d,'String',handles.default.nes);
set(handles.Sele_arrange_d,'String',handles.default.ele_arrange);
set(handles.Sele_order_d,'String',handles.default.ele_order);
set(handles.SSim_Len_d,'String',handles.default.Sim_Len);
set(handles.Sdt_d,'String',handles.default.dt);
set(handles.STOL_d,'String',handles.default.TOL);

set(handles.SWi_u,'String',handles.present.Wi);
set(handles.STh_u,'String',handles.present.Th);
set(handles.Sner_u,'String',handles.present.ner);
set(handles.Snes_u,'String',handles.present.nes);
set(handles.pop_ele_arrange,'Value',handles.present.ele_arrange_V);
handles.ele_arrange_V_u=handles.present.ele_arrange_V;
set(handles.pop_ele_order,'Value',handles.present.ele_order_V);
set(handles.SSim_Len_u,'String',handles.present.Sim_Len);
set(handles.Sdt_u,'String',handles.present.dt);
set(handles.STOL_u,'String',handles.present.TOL);
% Update handles structure

a=[handles.present.Wi,handles.present.Th];
b=[handles.present.ner,handles.present.nes];
c=[handles.present.ele_arrange_V,handles.present.ele_order_V];

XYH=handles.present.XYH;
MAPH=handles.present.MAPH;

% [XYH,MAPH,BC_Info]=Create_mesh(user)
handles.el_no=0;
handles.no_no=0;
set(handles.nodes,'Value',handles.no_no);
set(handles.elements,'Value',handles.el_no);
set(handles.SSymmetry,'Value',handles.present.symmetry+1);
guidata(hObject, handles);
% cla(handles.mesh);
% axes(handles.mesh);
draw_heat_4n(handles);


% guidata(hObject, handles);















% --- Executes on button press in Results_Line_Plot.
function Results_Line_Plot_Callback(hObject, eventdata, handles)
% hObject    handle to Results_Line_Plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global  pos_main_panel s_panel
set(handles.Panel_Results_Line_Plot,'Visible','On');
set(handles.Panel_Results_Isotherms,'Visible','Off');
% set(handles.Panel_Results_Report,'Visible','Off');

set(handles.Panel_Results,'Position',pos_main_panel);
set(handles.Panel_Default,'Visible','Off');drawnow;
set(handles.Panel_Advanced,'Visible','Off');drawnow;
set(handles.Panel_Results,'Visible','On');
set(handles.Panel_NewGrade,'Visible','Off');drawnow;

set(handles.Panel_Advanced_Second,'Visible','Off');
set(handles.Panel_Advanced_Mold,'Visible','Off');
set(handles.Panel_Advanced_CCM,'Visible','Off');
set(handles.Panel_Advanced_Simulation,'Visible','Off');

set(hObject,'Value',1);
set(handles.Results_Isotherms,'Value',0);
set(handles.Results_Report,'Value',0);
set(hObject,'BackgroundColor',[239,239,239]/255);
set(handles.Results_Isotherms,'BackgroundColor',[214,214,214]/255);
set(handles.Results_Report,'BackgroundColor',[214,214,214]/255);
% Hint: get(hObject,'Value') returns toggle state of Results_Line_Plot

handles.Folders.Sub_Fold.Sim_Fold;

Cur_Res_Folder=[handles.Folders.Main_Fold.Sim handles.Folders.Sub_Fold.Res_Fold];
Parameters=load([Cur_Res_Folder handles.path_sep 'Process_Parameters.mat']);
handles.WSM=handles.default.Water_share_matrix;
Display_Results_Parameters(Parameters,handles);

handles.Parameters=Parameters;
Line_Plots(Cur_Res_Folder,handles)
guidata(hObject, handles);
% [AXIAL,DATA]=find_axialdis_DB(handles);
% handles.AXIAL=AXIAL;
% handles.DB=DATA;
% dvalue=size(AXIAL,1);
% set(handles.Results_AxialDist,'String',AXIAL(dvalue,1));
% set(handles.Results_DataList,'Value',dvalue);
% dbindex=DATA{AXIAL(dvalue,2)};handles.dbindex=dbindex;
% guidata(hObject, handles);
% Temp_contour(handles);
% handles.daxial=AXIAL(dvalue,1);
% set(handles.Results_Midplane,'Value',1);
% guidata(hObject, handles);
% 
% Plot_Axial_Contours(handles)


% --- Executes on button press in Results_Line_Plot.
function Results_Isotherms_Callback(hObject, eventdata, handles)
% hObject    handle to Results_Line_Plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global  pos_main_panel s_panel
set(handles.Panel_Results_Line_Plot,'Visible','Off');
set(handles.Panel_Results_Isotherms,'Visible','On');
% set(handles.Panel_Results_Report,'Visible','Off');

set(handles.Panel_Results,'Position',pos_main_panel);
set(handles.Panel_Default,'Visible','Off');drawnow;
set(handles.Panel_Advanced,'Visible','Off');drawnow;
set(handles.Panel_Results,'Visible','On');
set(handles.Panel_NewGrade,'Visible','Off');drawnow;

set(handles.Panel_Advanced_Second,'Visible','Off');
set(handles.Panel_Advanced_Mold,'Visible','Off');
set(handles.Panel_Advanced_CCM,'Visible','Off');
set(handles.Panel_Advanced_Simulation,'Visible','Off');

set(hObject,'Value',1);
set(handles.Results_Line_Plot,'Value',0);
set(handles.Results_Report,'Value',0);
set(hObject,'BackgroundColor',[239,239,239]/255);
set(handles.Results_Line_Plot,'BackgroundColor',[214,214,214]/255);
set(handles.Results_Report,'BackgroundColor',[214,214,214]/255);


Cur_Res_Folder=[handles.Folders.Main_Fold.Sim handles.Folders.Sub_Fold.Res_Fold];
% Current_Results_Folder=[handles.Sim_Main_Fold handles.present.Sim_Fold];
Parameters=load([Cur_Res_Folder handles.path_sep 'Process_Parameters.mat']);
handles.WSM=handles.default.Water_share_matrix;
Display_Results_Parameters(Parameters,handles);

% handles.Current_Results_Folder=Current_Results_Folder;
handles.Parameters=Parameters;
% Line_Plots(Current_Results_Folder,handles)
guidata(hObject, handles);

[AXIAL,DATA]=find_axialdis_DB(handles);
handles.AXIAL=AXIAL;
handles.DB=DATA;
dvalue=size(AXIAL,1);
set(handles.Results_AxialDist,'String',AXIAL(dvalue,1));
set(handles.Results_DataList,'Value',dvalue);
dbindex=DATA{AXIAL(dvalue,2)};handles.dbindex=dbindex;
guidata(hObject, handles);
if get(handles.Con_1,'Value')==1
    Temp_contour(handles);
else
    Color_plot(handles);
end

handles.daxial=AXIAL(dvalue,1);
set(handles.Results_Midplane,'Value',1);
guidata(hObject, handles);

Plot_Axial_Contours(handles)

load([Cur_Res_Folder handles.path_sep 'Material.mat']);
set(handles.Results_TS,'String',round(TS*1)/1);
set(handles.Results_TL,'String',round(TL*1)/1);

set(handles.AC_Color_No,'Value',1);



% Hint: get(hObject,'Value') returns toggle state of Results_Line_Plot
% --- Executes on button press in Results_Line_Plot.

function Results_Report_Callback(hObject, eventdata, handles)
% hObject    handle to Results_Line_Plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global  pos_main_panel s_panel
set(handles.Panel_Results_Line_Plot,'Visible','Off');
set(handles.Panel_Results_Isotherms,'Visible','Off');
% set(handles.Panel_Results_Report,'Visible','On');

set(handles.Panel_Results,'Position',pos_main_panel);
set(handles.Panel_Default,'Visible','Off');drawnow;
set(handles.Panel_Advanced,'Visible','Off');drawnow;
set(handles.Panel_Results,'Visible','On');
set(handles.Panel_NewGrade,'Visible','Off');drawnow;

set(handles.Panel_Advanced_Second,'Visible','Off');
set(handles.Panel_Advanced_Mold,'Visible','Off');
set(handles.Panel_Advanced_CCM,'Visible','Off');
set(handles.Panel_Advanced_Simulation,'Visible','Off');
% 
set(hObject,'Value',1);
set(handles.Results_Line_Plot,'Value',0);
set(handles.Results_Isotherms,'Value',0);
set(hObject,'BackgroundColor',[239,239,239]/255);
set(handles.Results_Line_Plot,'BackgroundColor',[214,214,214]/255);
set(handles.Results_Isotherms,'BackgroundColor',[214,214,214]/255);
% % Hint: get(hObject,'Value') returns toggle state of Results_Line_Plot

% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton13.
function pushbutton13_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton12.
function pushbutton12_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton11.
function pushbutton11_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton10.
function pushbutton10_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)




% --- Executes on button press in pushbutton14.
function pushbutton14_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton15.
function pushbutton15_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton19.
function pushbutton19_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton18.
function pushbutton18_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton17.
function pushbutton17_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton16.
function pushbutton16_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit20_Callback(hObject, eventdata, handles)
% hObject    handle to edit20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit20 as text
%        str2double(get(hObject,'String')) returns contents of edit20 as a double


% --- Executes during object creation, after setting all properties.
function edit20_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function GN_Callback(hObject, eventdata, handles)
% hObject    handle to GN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of GN as text
%        str2double(get(hObject,'String')) returns contents of GN as a double
GN_s=str2double(get(hObject,'String'));
set(handles.status_text,'String','');
cla(handles.axes_status);axes(handles.axes_status);
patch([0,0,0,0],[0,0,0,0],'w');
set(handles.axes_status,'Visible','off');
handles.C_DB=0;
handles.C_DB2=0;
load([handles.Folders.Main_Fold.Grd  handles.Folders.File.Grd  '.mat']);
GNr=cell2mat({Steel_reciepes_SC.Grade_Number});
Lsel=find(GNr==GN_s);
if isempty(Lsel)
    Lsel=1;
    errordlg('Grade is NOT in List','Not Exist');
end

[Grade_depend]=Get_Grade_depend(Steel_reciepes_SC,handles.default.Std_cooling,Lsel);
%%%% updating grade_nr and grade_Ty in sss and usersettings
handles.user.Grade_depend.Grade_Number=Grade_depend.Grade_Number;
handles.user.Grade_depend.Grade_Type=Grade_depend.Grade_Type;
handles.default.Grade_depend=Grade_depend;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
handles.Grade_depend=Grade_depend;
% guidata(hObject,handles);
% if get(handles.mat_off,'Value')==1
%     handles.present.Grade_depend=handles.user.Grade_depend;
% else
%     handles.present.Grade_depend=handles.Grade_depend;
% end
handles.present.Grade_depend=Grade_depend;
guidata(hObject,handles);
% set(handles.GN,'String',handles.present.Grade_depend.Grade_Number);
set(handles.GNList,'Value',Lsel);
set(handles.GT,'String',handles.present.Grade_depend.Grade_Type);
set(handles.CT,'String',handles.present.Grade_depend.Init_Temp);
set(handles.CS,'String',handles.present.Grade_depend.Casting_speed);
set(handles.CM,'String',handles.present.Grade_depend.Cool_Mode);
WFR_Lm3=handles.present.Grade_depend.WFR_Lm3;
set(handles.Z1,'String',WFR_Lm3(1));
set(handles.Z2,'String',WFR_Lm3(2));
set(handles.Z3,'String',WFR_Lm3(3));

% [SEC_COOL,Nozzle_Tab,No_Noz_Zo,cooling_info]=Updade_Cooling_Info(handles);
% handles.present.SEC_COOL=SEC_COOL;
% handles.present.No_Noz_Zo=No_Noz_Zo;
% handles.present.Nozzle_Tab=Nozzle_Tab;
% handles.present.cooling_info=cooling_info;
% 
% handles.default.SEC_COOL=SEC_COOL;
% handles.default.Nozzle_Tab=Nozzle_Tab;
% handles.default.cooling_info=cooling_info;
% handles.default.No_Noz_Zo=No_Noz_Zo;

[Properties]=Load_Grade_Advanced(handles);
set(handles.Liquidus,'String',round(Properties.TL*10)/10);
set(handles.Solidus,'String',round(Properties.TS*10)/10);
set(handles.Latent,'String',round(Properties.Latent/1e3*10)/10);
handles.Properties=Properties;
guidata(hObject,handles);
plot_material_data(handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Sim_Fold=[num2str(handles.Grade_depend.Grade_Number) '_' handles.Grade_depend.Grade_Type];
% Folder_Name=[handles.Fold_Path  handles.default.Sim_Main_Fold '\' Sim_Fold];
set(handles.FolderName,'String',Sim_Fold);
handles.Folders.Sub_Fold.Sim_Fold=Sim_Fold;
handles.Folders.Sub_Fold.Res_Fold=Sim_Fold;
handles.Folders.Sub_Fold.Mat_Fold=Sim_Fold;
% handles.present.Sim_Fold=Sim_Fold;
% handles.Folder_Name=Folder_Name;


guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function GN_CreateFcn(hObject, eventdata, handles)
% hObject    handle to GN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in RunSim.
function RunSim_Callback(hObject, eventdata, handles)
% hObject    handle to RunSim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Load_Grade(handles);
% Sim_Len=str2double(get(handles.SimLength,'String'));

%%% value=1, running, value=0, stop
FN=get(handles.FolderName,'String');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sim_Fold=[num2str(handles.present.Grade_depend.Grade_Number) '_' handles.present.Grade_depend.Grade_Type];
% Folder_Name1=[handles.Folders.Main_Fold.Sim  handles.Folders.Main_Fold.Sim];
Folder_Name=[handles.Folders.Main_Fold.Sim FN];
d=dir(handles.Folders.Main_Fold.Sim);
str={d.name};
if length(str)>2
    str=str(3:length(str));
end
Exist=strcmp(str,FN);

val=0;
if sum(Exist)>0 && get(hObject, 'Value') == 1
 % if folder exists
    val=2;
    choice = questdlg('Do you want to overwrite existing Simulations???', ...
        'Atttention', ...
        'YES','NO','Continue Simulation','NO');

    switch choice
        case 'YES'
            val=1;
        case 'NO' 
            val=2;            
        case 'Continue Simulation' 
            val=3;
        case ''
            val=2;
    end
else % if folder doesnot exist
    val=4;
    if isdir(Folder_Name)==0
        mkdir(Folder_Name);
        handles.C_DB=0;
    end
end
if val==1 || val==3 || val==4
    handles.C_DB2=0;
    set(handles.advanced_status_text,'String','');
    h_axes=handles.advanced_status;set(h_axes,'Visible','Off');
    cla(handles.advanced_status);axes(handles.advanced_status);axis off;
    patch([0,0,0,0],[0,0,0,0],'w');
    set(handles.Advanced_Run,'Value',0);
        if val==1
            delete([Folder_Name handles.path_sep 'DB*.mat']);
            handles.C_DB=0;
        elseif val==3
            pa=[Folder_Name handles.path_sep];
            delete_data_file(pa,handles.C_DB);
        end
    handles.present_previous=handles.present;    
    guidata(hObject, handles);

    %%% Material --------------------------------------------------------------
    KHT=handles.Properties.KHT;
    SHT=handles.Properties.SHT;
    RHT=handles.Properties.RHT;
    PCF=handles.Properties.PCF;
    ENTH=handles.Properties.ENTH;
    Latent=handles.Properties.Latent;
    TS=handles.Properties.TS;
    TL=handles.Properties.TL;

[SEC_COOL,Nozzle_Tab,No_Noz_Zo,cooling_info]=Updade_Cooling_Info(handles,3);
handles.present.SEC_COOL=SEC_COOL;
handles.present.No_Noz_Zo=No_Noz_Zo;
handles.present.Nozzle_Tab=Nozzle_Tab;
handles.present.cooling_info=cooling_info;
guidata(hObject, handles);
    m=handles.present;    
   save([Folder_Name handles.path_sep 'Material.mat'],'KHT','SHT','RHT','TS','TL','PCF','ENTH',...
        'Latent');
    save([Folder_Name handles.path_sep 'Process_Parameters.mat'],'-struct','m');
%     
    Init_Temp=handles.present.Grade_depend.Init_Temp;
    XYH=handles.present.XYH;
    MAPH=handles.present.MAPH;
    create_intial_DB(Init_Temp,XYH,MAPH,Folder_Name,handles);
    Vc=m.Grade_depend.Casting_speed/60;
    dx=m.dt*Vc;
    NDB=round(m.Sim_Len*1e-3/dx);
%     %     Iterative_Solver(Folder_Name);
    
    h_axes_status=handles.axes_status;set(h_axes_status,'Visible','On');
    cla(handles.axes_status);
    axes(handles.axes_status);
    C_DB=handles.C_DB;
    pr=(C_DB+1)/(NDB+1)*100;
    patch([0,1,1,0],[0,0,1,1],'w');
    patch([0,pr/100,pr/100,0],[0,0,1,1],'r');
    axis([0,1,0,1]);axis off;drawnow;

    if get(hObject, 'Value') == 1
        set(hObject, 'String','STOP Simulation');
        while C_DB<NDB+1
%             p1=handles.present;
% %             p1.Grade_depend 
% %             Axial=C_DB+1
%             break
% Folder_Name,C_DB,handles.path_sep
            Axial=Solver_new(Folder_Name,C_DB,handles.path_sep);
%             Axial=Solver_enthalpy(Folder_Name,C_DB,handles.path_sep);
            
                        
            pause(1e-3);
            pr=(C_DB+1)/(NDB+1)*100;
            patch([0,1,1,0],[0,0,1,1],'w');
            patch([0,pr/100,pr/100,0],[0,0,1,1],'r');
            axis([0,1,0,1]);axis off;drawnow;
            str=[ num2str(round(pr*1)/1) '   % Completed : ' num2str(round(Axial*100)/100) '   m reached'];
            set(handles.status_text,'String',str);
            if get(hObject, 'Value') == 0
                set(hObject, 'String','RUN Simulation');
                handles.C_DB=handles.C_DB+1;
                axis off;
                guidata(hObject, handles);
                break;
            end
            C_DB=C_DB+1;
            handles.C_DB=handles.C_DB+1;
            axis off;
            guidata(hObject, handles);
        end

        if C_DB==NDB+1
            set(hObject, 'Value',0);
            set(hObject, 'String','RUN Simulation');
        end
        
    else
        set(hObject, 'String','RUN Simulation');
    end
    axis off;
    
else
    set(hObject, 'Value',0);
    set(hObject, 'String','RUN Simulation');
end
 guidata(hObject, handles);

% --- Executes on button press in pushbutton23.
function pushbutton23_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function SimLength_Callback(hObject, eventdata, handles)
% hObject    handle to SimLength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SimLength as text
%        str2double(get(hObject,'String')) returns contents of SimLength as a double
SimLeni=str2double(get(hObject,'String'));
% handles.SimLeni=[];
% handles.SimLeni=SimLen;
% --- Executes during object creation, after setting all properties.
function SimLength_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SimLength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in GNList.
function GNList_Callback(hObject, eventdata, handles)
% hObject    handle to GNList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns GNList contents as cell array
%        contents{get(hObject,'Value')} returns selected item from GNList
Lsel=get(hObject,'Value');
set(handles.status_text,'String','');
cla(handles.axes_status);axes(handles.axes_status);
patch([0,0,0,0],[0,0,0,0],'w');
set(handles.axes_status,'Visible','off');
handles.C_DB=0;
handles.C_DB2=0;

Steel_reciepes_SC=handles.Steel_reciepes_SC;

% load([handles.Folders.Main_Fold.Grd  handles.Folders.File.Grd  '.mat']);
[Grade_depend]=Get_Grade_depend(Steel_reciepes_SC,handles.default.Std_cooling,Lsel);
%%%% updating grade_nr and grade_Ty in sss and usersettings
handles.user.Grade_depend.Grade_Number=Grade_depend.Grade_Number;
handles.user.Grade_depend.Grade_Type=Grade_depend.Grade_Type;
handles.default.Grade_depend=Grade_depend;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
handles.Grade_depend=Grade_depend;
% guidata(hObject,handles);
% if get(handles.mat_off,'Value')==1
%     handles.present.Grade_depend=handles.user.Grade_depend;
% else
%     handles.present.Grade_depend=handles.Grade_depend;
% end
handles.present.Grade_depend=Grade_depend;
guidata(hObject,handles);
% handles.present.Grade_depend;
set(handles.GN,'String',handles.present.Grade_depend.Grade_Number);
set(handles.GT,'String',handles.present.Grade_depend.Grade_Type);
set(handles.CT,'String',handles.present.Grade_depend.Init_Temp);
set(handles.CS,'String',handles.present.Grade_depend.Casting_speed);
set(handles.CM,'String',handles.present.Grade_depend.Cool_Mode);
WFR_Lm3=handles.present.Grade_depend.WFR_Lm3;
set(handles.Z1,'String',WFR_Lm3(1));
set(handles.Z2,'String',WFR_Lm3(2));
set(handles.Z3,'String',WFR_Lm3(3));




% [SEC_COOL,Nozzle_Tab,No_Noz_Zo,cooling_info]=Updade_Cooling_Info(handles);
% handles.present.SEC_COOL=SEC_COOL;
% handles.present.Nozzle_Tab=Nozzle_Tab;
% handles.present.cooling_info=cooling_info;
% HTC1=SEC_COOL.HTC;
% HTC1(50,:)*HTC1(50,:)'
% 
% handles.default.SEC_COOL=SEC_COOL;
% handles.default.Nozzle_Tab=Nozzle_Tab;
% handles.default.cooling_info=cooling_info;
% guidata(hObject,handles)



[Properties]=Load_Grade_Advanced(handles);
set(handles.Liquidus,'String',round(Properties.TL*10)/10);
set(handles.Solidus,'String',round(Properties.TS*10)/10);
set(handles.Latent,'String',round(Properties.Latent/1e3*10)/10);
handles.Properties=Properties;
guidata(hObject,handles);
plot_material_data(handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Sim_Fold=[num2str(handles.Grade_depend.Grade_Number) '_' handles.Grade_depend.Grade_Type];
% Folder_Name=[handles.Fold_Path  handles.default.Sim_Main_Fold '\' Sim_Fold];
set(handles.FolderName,'String',Sim_Fold);
handles.Folders.Sub_Fold.Sim_Fold=Sim_Fold;
handles.Folders.Sub_Fold.Res_Fold=Sim_Fold;
handles.Folders.Sub_Fold.Mat_Fold=Sim_Fold;
% handles.present.Sim_Fold=Sim_Fold;
% handles.Folder_Name=Folder_Name;


guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function GNList_CreateFcn(hObject, eventdata, handles)
% hObject    handle to GNList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function z1a_CreateFcn(hObject, eventdata, handles)
% hObject    handle to z1a (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function z1b_CreateFcn(hObject, eventdata, handles)
% hObject    handle to z1b (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function z3a_CreateFcn(hObject, eventdata, handles)
% hObject    handle to z3a (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function z3b_CreateFcn(hObject, eventdata, handles)
% hObject    handle to z3b (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function z2a_CreateFcn(hObject, eventdata, handles)
% hObject    handle to z2a (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function z2b_CreateFcn(hObject, eventdata, handles)
% hObject    handle to z2b (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called



function Z2_Callback(hObject, eventdata, handles)
% hObject    handle to Z2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Z2 as text
%        str2double(get(hObject,'String')) returns contents of Z2 as a double


% --- Executes during object creation, after setting all properties.
function Z2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Z2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Z3_Callback(hObject, eventdata, handles)
% hObject    handle to Z3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Z3 as text
%        str2double(get(hObject,'String')) returns contents of Z3 as a double


% --- Executes during object creation, after setting all properties.
function Z3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Z3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Z1_Callback(hObject, eventdata, handles)
% hObject    handle to Z1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Z1 as text
%        str2double(get(hObject,'String')) returns contents of Z1 as a double


% --- Executes during object creation, after setting all properties.
function Z1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Z1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function F1a_Callback(hObject, eventdata, handles)
% hObject    handle to F1a (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of F1a as text
%        str2double(get(hObject,'String')) returns contents of F1a as a double


% --- Executes during object creation, after setting all properties.
function F1a_CreateFcn(hObject, eventdata, handles)
% hObject    handle to F1a (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function F1b_Callback(hObject, eventdata, handles)
% hObject    handle to F1b (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of F1b as text
%        str2double(get(hObject,'String')) returns contents of F1b as a double


% --- Executes during object creation, after setting all properties.
function F1b_CreateFcn(hObject, eventdata, handles)
% hObject    handle to F1b (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function F3a_Callback(hObject, eventdata, handles)
% hObject    handle to F3a (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of F3a as text
%        str2double(get(hObject,'String')) returns contents of F3a as a double


% --- Executes during object creation, after setting all properties.
function F3a_CreateFcn(hObject, eventdata, handles)
% hObject    handle to F3a (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function F3b_Callback(hObject, eventdata, handles)
% hObject    handle to F3b (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of F3b as text
%        str2double(get(hObject,'String')) returns contents of F3b as a double


% --- Executes during object creation, after setting all properties.
function F3b_CreateFcn(hObject, eventdata, handles)
% hObject    handle to F3b (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function F2a_Callback(hObject, eventdata, handles)
% hObject    handle to F2a (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of F2a as text
%        str2double(get(hObject,'String')) returns contents of F2a as a double


% --- Executes during object creation, after setting all properties.
function F2a_CreateFcn(hObject, eventdata, handles)
% hObject    handle to F2a (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function F2b_Callback(hObject, eventdata, handles)
% hObject    handle to F2b (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of F2b as text
%        str2double(get(hObject,'String')) returns contents of F2b as a double


% --- Executes during object creation, after setting all properties.
function F2b_CreateFcn(hObject, eventdata, handles)
% hObject    handle to F2b (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in reset1.
function reset1_Callback(hObject, eventdata, handles)
% hObject    handle to reset1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in ok1.
function ok1_Callback(hObject, eventdata, handles)
% hObject    handle to ok1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function tolerance_Callback(hObject, eventdata, handles)
% hObject    handle to tolerance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tolerance as text
%        str2double(get(hObject,'String')) returns contents of tolerance as a double


% --- Executes during object creation, after setting all properties.
function tolerance_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tolerance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function stepsize_Callback(hObject, eventdata, handles)
% hObject    handle to stepsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of stepsize as text
%        str2double(get(hObject,'String')) returns contents of stepsize as a double


% --- Executes during object creation, after setting all properties.
function stepsize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to stepsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in arrangement.
function arrangement_Callback(hObject, eventdata, handles)
% hObject    handle to arrangement (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns arrangement contents as cell array
%        contents{get(hObject,'Value')} returns selected item from arrangement


% --- Executes during object creation, after setting all properties.
function arrangement_CreateFcn(hObject, eventdata, handles)
% hObject    handle to arrangement (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in order.
function order_Callback(hObject, eventdata, handles)
% hObject    handle to order (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns order contents as cell array
%        contents{get(hObject,'Value')} returns selected item from order


% --- Executes during object creation, after setting all properties.
function order_CreateFcn(hObject, eventdata, handles)
% hObject    handle to order (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in elements.
function elements_Callback(hObject, eventdata, handles)
% hObject    handle to elements (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of elements
handles.el_no=get(hObject,'Value');
% cla(handles.mesh);
% axes(handles.mesh);
guidata(hObject, handles);
draw_heat_4n(handles);


% --- Executes on button press in nodes.
function nodes_Callback(hObject, eventdata, handles)
% hObject    handle to nodes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of nodes
% [XYH,MAPH,BC_Info]=Create_mesh(user)
handles.no_no=get(hObject,'Value');
% cla(handles.mesh);
% axes(handles.mesh);
guidata(hObject, handles);
draw_heat_4n(handles);

function NrWidth_Callback(hObject, eventdata, handles)
% hObject    handle to NrWidth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of NrWidth as text
%        str2double(get(hObject,'String')) returns contents of NrWidth as a double


% --- Executes during object creation, after setting all properties.
function NrWidth_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NrWidth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function NrThick_Callback(hObject, eventdata, handles)
% hObject    handle to NrThick (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of NrThick as text
%        str2double(get(hObject,'String')) returns contents of NrThick as a double


% --- Executes during object creation, after setting all properties.
function NrThick_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NrThick (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function width_Callback(hObject, eventdata, handles)
% hObject    handle to width (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of width as text
%        str2double(get(hObject,'String')) returns contents of width as a double


% --- Executes during object creation, after setting all properties.
function width_CreateFcn(hObject, eventdata, handles)
% hObject    handle to width (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function thick_Callback(hObject, eventdata, handles)
% hObject    handle to thick (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of thick as text
%        str2double(get(hObject,'String')) returns contents of thick as a double


% --- Executes during object creation, after setting all properties.
function thick_CreateFcn(hObject, eventdata, handles)
% hObject    handle to thick (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu5.
function popupmenu5_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu5 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu5


% --- Executes during object creation, after setting all properties.
function popupmenu5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ok2.
function ok2_Callback(hObject, eventdata, handles)
% hObject    handle to ok2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in reset2.
function reset2_Callback(hObject, eventdata, handles)
% hObject    handle to reset2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function CT_Callback(hObject, eventdata, handles)
% hObject    handle to CT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of CT as text
%        str2double(get(hObject,'String')) returns contents of CT as a double


% --- Executes during object creation, after setting all properties.
function CT_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function CS_Callback(hObject, eventdata, handles)
% hObject    handle to CS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of CS as text
%        str2double(get(hObject,'String')) returns contents of CS as a double


% --- Executes during object creation, after setting all properties.
function CS_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in boundary.
function boundary_Callback(hObject, eventdata, handles)
% hObject    handle to boundary (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns boundary contents as cell array
%        contents{get(hObject,'Value')} returns selected item from boundary


% --- Executes during object creation, after setting all properties.
function boundary_CreateFcn(hObject, eventdata, handles)
% hObject    handle to boundary (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in reset3.
function reset3_Callback(hObject, eventdata, handles)
% hObject    handle to reset3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in ok3.
function ok3_Callback(hObject, eventdata, handles)
% hObject    handle to ok3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function Flow_Callback(hObject, eventdata, handles)
% hObject    handle to Flow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Flow as text
%        str2double(get(hObject,'String')) returns contents of Flow as a double


% --- Executes during object creation, after setting all properties.
function Flow_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Flow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ND_Callback(hObject, eventdata, handles)
% hObject    handle to ND (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ND as text
%        str2double(get(hObject,'String')) returns contents of ND as a double


% --- Executes during object creation, after setting all properties.
function ND_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ND (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function CA_Callback(hObject, eventdata, handles)
% hObject    handle to CA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of CA as text
%        str2double(get(hObject,'String')) returns contents of CA as a double


% --- Executes during object creation, after setting all properties.
function CA_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function strt_Callback(hObject, eventdata, handles)
% hObject    handle to strt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of strt as text
%        str2double(get(hObject,'String')) returns contents of strt as a double


% --- Executes during object creation, after setting all properties.
function strt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to strt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit133_Callback(hObject, eventdata, handles)
% hObject    handle to edit133 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit133 as text
%        str2double(get(hObject,'String')) returns contents of edit133 as a double


% --- Executes during object creation, after setting all properties.
function edit133_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit133 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ZL_Callback(hObject, eventdata, handles)
% hObject    handle to ZL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ZL as text
%        str2double(get(hObject,'String')) returns contents of ZL as a double


% --- Executes during object creation, after setting all properties.
function ZL_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ZL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function lateral_Callback(hObject, eventdata, handles)
% hObject    handle to lateral (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lateral as text
%        str2double(get(hObject,'String')) returns contents of lateral as a double


% --- Executes during object creation, after setting all properties.
function lateral_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lateral (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function nozzels_Callback(hObject, eventdata, handles)
% hObject    handle to nozzels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nozzels as text
%        str2double(get(hObject,'String')) returns contents of nozzels as a double


% --- Executes during object creation, after setting all properties.
function nozzels_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nozzels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listbox3.
function listbox3_Callback(hObject, eventdata, handles)
% hObject    handle to listbox3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox3


% --- Executes during object creation, after setting all properties.
function listbox3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function[CM,CMV]= cool_mode_sel(CoMo)
if CoMo(1)==1
    % set(handles.coolingmode, 'SelectedObject', handles.ultramild);
    CM='Ultra Mild';
    CMV=1;
elseif CoMo(2)==1
    % set(handles.coolingmode, 'SelectedObject', handles.mild);
    CM='Mild';CMV=2;
elseif CoMo(3)==1
    % set(handles.coolingmode, 'SelectedObject', handles.moderate);
    CM='Moderate';CMV=3;
elseif CoMo(4)==1
    % set(handles.coolingmode, 'SelectedObject', handles.strong);
    CM='Strong';CMV=4;
else
    CM='Other';CMV=5;
end



% function[]=Load_Grade(handles)
% % Data_FN         =   handles.FNi;
% Grade_Ty        =   handles.GTi;
% Grade_Nr        =   handles.GNi;
% Casting_speed   =   handles.CSi;
% Data_FN         =   handles.Data_FNi;
% % Sim_FN          =   handles.FNi;
% Init_Temp       =   handles.CTi;
% 
% Sim_FN=get(handles.FolderName,'String')
% 
% if isdir(Sim_FN)==0
%     mkdir(Sim_FN);
% end
% % [Sim_FN '\Material.mat']
% % isdir([Sim_FN '\Material.mat'])
% % if isdir([Sim_FN '\Material.mat'])==1
% %     delete([Sim_FN '\Material.mat']);
% % end
% FN=[Data_FN '\' num2str(Grade_Nr) '_' Grade_Ty '\']
% RT1=importdata([FN num2str(Grade_Nr) '_density.txt']);
% KT1=importdata([FN num2str(Grade_Nr) '_thermal_conduct.txt']);
% CT1=importdata([FN num2str(Grade_Nr) '_specific_heat.txt']);
% SF1=importdata([FN num2str(Grade_Nr) '_phases.txt']);
% ENTH1=importdata([FN num2str(Grade_Nr) '_enthalpy.txt']);
% a1=SF1.data;
% KT=KT1.data;
% RT=RT1.data;
% CT=CT1.data;
% ENTH=ENTH1.data;
% n=1;
% for i=1:size(a1,1)
%         if a1(i,2)>=0
%             n=i;
%         end
% end
% a=a1(1:n,[1 2]);
% SF=a;
% n1=1;n2=1;
% for i=1:size(a,1)-1
%     if (a(i+1,2)>0 && a(i,2)<=1e-2) || a(i+1,2)<=1e-2 && a(i,2)>0
%         n2=i+1;
%     end
%     
%     if (a(i+1,2)>=100 && a(i,2)<100) || (a(i+1,2)<100 && a(i,2)>=100)
%         n1=i;
%     end
% end
% 
% TLi=a(n1);
% TSi=a(n2);
% 
% set(handles.Solidus,'String',round(TSi*10)/10);
% set(handles.Liquidus,'String',round(TLi*10)/10);
% 
% f1=phase_interpolation(ENTH,round(TLi*10)/10);
% f2=phase_interpolation(ENTH,round(TSi*10)/10);
% Latent=abs(f1-f2);
% set(handles.Latent,'String',round(Latent*10)/10);
% 
% RHT=RT;RHT(:,2)=RHT(:,2).*1000;
% KHT=KT;
% SHT=CT;SHT(:,2)=SHT(:,2).*1000;
% TS=TSi;TL=TLi;
% PCF=SF;PCF(:,2)=PCF(:,2)./100;
% 
% Latent=Latent*1e3;
% save([Sim_FN '\' 'Material.mat'],'KHT','SHT','RHT','TS','TL','PCF',...
%     'Init_Temp','Grade_Nr','Grade_Ty','Latent','Casting_speed');
% 
% cool_mode=get(handles.CM,'String');
% SBC=[Sim_FN  '\cool_mode.mat'];
% save(SBC,'cool_mode');

% function[xvl]=phase_interpolation(X,y)
% nX=size(X,1);
%         if (X(1,1)-X(nX,1))>0 %%% descending order
%             if y>=X(1,1)
%                 xvl=X(1,2);
%             elseif y<=X(nX,1)
%                 xvl=X(nX,2);
%             else
%                 for jj=2:nX
%                     if y>=X(jj,1) % Make linear interpolation
%                         xvl = X(jj-1,2)+(y-X(jj-1,1))*(X(jj,2)-X(jj-1,2))/(X(jj,1)-X(jj-1,1));
%                         break;
%                     end
%                 end
%             end
%         else %%% ascending order
%             if y<=X(1,1)
%                 xvl=X(1,2);
%             elseif y<=X(nX,1)
%                 xvl=X(nX,2);
%             else
%                 for jj=2:nX
%                     if y<=X(jj,1) % Make linear interpolation
%                         xvl = X(jj-1,2)+(y-X(jj-1,1))*(X(jj,2)-X(jj-1,2))/(X(jj,1)-X(jj-1,1));
%                         break;
%                     end
%                 end
%             end
%         end

function[Properties]=Load_Grade_Advanced(handles)
%%%%% Load material data from the stored material data bank   
Mat_Data=handles.Mat_Data;
GradeList=cell2mat({Mat_Data.Grade_Nr});
Grade_Nr        =   handles.present.Grade_depend.Grade_Number;
% FN=[handles.Folders.Main_Fold.Mat num2str(Grade_Nr) '.mat'];
pos=find(GradeList==Grade_Nr);
if length(pos)==0
    warndlg('Material Data for the selected Grade does not exist. Default material properties will be used for calculations.','Materail Data Not Exist');
    Grade_Nr=1000;
    pos=find(GradeList==Grade_Nr);
    uiwait
end
% FN=[handles.Folders.Main_Fold.Mat num2str(Grade_Nr) '.mat'];
Properties=Mat_Data(pos).Property;

%%%% second old method
% Grade_Nr        =   handles.present.Grade_depend.Grade_Number;
% FN=[handles.Folders.Main_Fold.Mat num2str(Grade_Nr) '.mat'];
% if exist(FN, 'file')==0
%     warndlg('Material Data for the selected Grade does not exist. Default material properties will be used for calculations.','Materail Data Not Exist');
%     Grade_Nr=1000;
%     uiwait
% end
% FN=[handles.Folders.Main_Fold.Mat num2str(Grade_Nr) '.mat'];
% load(FN);

%%%%%%%%%%%%%%%% OLD METHOD 
% %%%%% with material database creation
% Grade_Nr        =   handles.present.Grade_depend.Grade_Number;
% Grade_Ty        =   handles.present.Grade_depend.Grade_Type;
% Property_Tag    =   handles.default.Property_Tag;
% FN=[handles.Folders.Main_Fold.Mat num2str(Grade_Nr) '_' Grade_Ty handles.path_sep];
% if isdir(FN)==0
%     warndlg('Material Data for the selected Grade does not exist. Default material properties will be used for calculations.','Materail Data Not Exist');
%     FN=[handles.Folders.Main_Fold.Mat 'Default' handles.path_sep];
%     Grade_Nr=1000;
%     uiwait
% end
% RT1=importdata([FN num2str(Grade_Nr) '_' Property_Tag.RHT]);
% KT1=importdata([FN num2str(Grade_Nr) '_' Property_Tag.KHT]);
% CT1=importdata([FN num2str(Grade_Nr) '_' Property_Tag.SHT]);
% SF1=importdata([FN num2str(Grade_Nr) '_' Property_Tag.PCF]);
% ENTH1=importdata([FN num2str(Grade_Nr) '_' Property_Tag.ENTH]);
% a1=SF1.data;
% KT=KT1.data;
% RT=RT1.data;
% CT=CT1.data;
% ENTH=ENTH1.data;
% n=1;
% for i=1:size(a1,1)
%         if a1(i,2)>=0
%             n=i;
%         end
% end
% a=a1(1:n,[1 2]);
% SF=a;
% n1=1;n2=1;
% for i=1:size(a,1)-1
%     if (a(i+1,2)>0 && a(i,2)<=1e-2) || a(i+1,2)<=1e-2 && a(i,2)>0
%         n2=i+1;
%     end
%     
%     if (a(i+1,2)>=100 && a(i,2)<100) || (a(i+1,2)<100 && a(i,2)>=100)
%         n1=i;
%     end
% end
% 
% TLi=a(n1);
% TSi=a(n2);
% f1=Linear_Interpolation(ENTH,round(TLi*100)/100);
% f2=Linear_Interpolation(ENTH,round(TSi*100)/100);
% Latent=abs(f1-f2);
% RHT=RT;RHT(:,2)=RHT(:,2).*1000;
% KHT=KT;
% SHT=CT;SHT(:,2)=SHT(:,2).*1000;
% TS=TSi;TL=TLi;
% PCF=SF;PCF(:,2)=PCF(:,2)./100;
% [a,b]= unique(PCF(:,2),'last');
% A(:,1)=PCF(b,1);A(:,2)=a;
% PCF=[];PCF=A;     
% 
% ENTH(:,2)=ENTH(:,2).*1000;
% Latent=Latent*1e3;
% 
% 
% [KHT]=find_unique_descending_order(KHT);
% [RHT]=find_unique_descending_order(RHT);
% [SHT_r]=find_unique_descending_order(SHT);
% [PCF]=find_unique_descending_order(PCF);
% 
% TS=max(PCF(find(PCF(:,2)<=0.01),1));
% TL=min(PCF(find(PCF(:,2)>=0.99),1));
% 
% [SHT,L_corr]=find_new_specific_heat(SHT,TL,TS);
% 
% Latent=Latent-L_corr;
% 
% 
% Properties=struct('KHT',{KHT},'SHT_r',{SHT_r},'SHT',{SHT},'RHT',{RHT},...
%     'TS',{TS},'TL',{TL},'Latent',{Latent}, 'PCF',{PCF},'ENTH',{ENTH});
% % handles.Properties=Properties;

    function[A,L_corr]= find_new_specific_heat(SHT,TL,TS)
        dT=50;
        CpL=Linear_Interpolation(SHT,TL+dT);
        CpS=Linear_Interpolation(SHT,TS-dT);
        L_corr=CpL*(TL-TS)-((CpL-CpS))/((TL-TS))*(TL^2+TS^2-2*TL*TS)/2;
        A=[];ct=0;
        N=size(SHT,1);
        for i=1:N-1
            if SHT(i,1)>=TL+dT && SHT(i+1,1)<TL+dT
                n1=i;
            end
            if SHT(i+1,1)<=TS-dT && SHT(i,1)>TS-dT
                n2=i+1;
            end
        end
        A(1:n1,:)=SHT(1:n1,:);
        A(n1+1,:)=[TL+dT CpL];
        A(n1+2,:)=[TL CpL];
        A(n1+3,:)=[TS CpS];
        A(n1+4,:)=[TS-dT CpS];
        A(n1+5:N-(n2-n1-5),:)=SHT(n2:N,:);
% A(n1+5,:)=[0 CpS];



    function[C]=find_unique_descending_order(A)
        [B(:,1),ix]=unique(A(:,1),'last');
        B(:,2)=A(ix,2);
        [C(:,1),ix]=sort(B(:,1),'descend');
        C(:,2)=B(ix,2);


function[]=plot_material_data(handles)
cla(handles.mat_prop);
axes(handles.mat_prop);
vprop=handles.vprop;
if vprop==1
    TData=handles.Properties.KHT;
elseif vprop==2
    TData=handles.Properties.RHT;
elseif vprop==3
    TData=handles.Properties.SHT_r;
elseif vprop==4
    TData=handles.Properties.PCF;
else
    TData=handles.Properties.ENTH;
end
switch vprop
    case 1
        %%%%% Thermal conductivity [KHT1]
%         KHT1=get(handles.KHT1,'Data')
%         KHT1=str2double(get(handles.KHT1,'Data'));
        plot(TData(:,1),TData(:,2),'r','LineWidth',3);
        ylabel('Thermal Conductivity [W/mK]','FontSize',12);
        xlabel('Temperature [degC]','FontSize',12);
    case 2
        %%%%% Specific Heat Capacity [SHT1]
%         SHT1=str2double(get(handles.SHT1,'Data'));
       plot(TData(:,1),TData(:,2)./1e3,'b','LineWidth',3);
        ylabel('Density [g/cm^3]','FontSize',12);
        xlabel('Temperature [degC]','FontSize',12);
    case 3
        %%%%% Density [RHT1]
%         RHT1=str2double(get(handles.RHT1,'Data'));
      plot(TData(:,1),TData(:,2)./1e3,'k','LineWidth',3);
       ylabel('Specific heat [J/gK]','FontSize',12);
        xlabel('Temperature [degC]','FontSize',12);
      
    case 4
              plot(TData(:,1),TData(:,2),'m','LineWidth',3);
        ylabel('Liquid fraction [-]','FontSize',12);
        xlabel('Temperature [degC]','FontSize',12);
    case 5
        plot(TData(:,1),TData(:,2)./1e3,'m','LineWidth',3);
        ylabel('Enthalpy [J/g]','FontSize',12);
        xlabel('Temperature [degC]','FontSize',12);        
end
set(gca,'FontSize',12,'LineWidth',2.0,'Box','off');

function FolderName_Callback(hObject, eventdata, handles)
% hObject    handle to FolderName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of FolderName as text
%        str2double(get(hObject,'String')) returns contents of FolderName as a double
Sim_Fold=get(hObject,'String');
handles.Folders.Sub_Fold.Sim_Fold=Sim_Fold;
handles.Folders.Sub_Fold.Res_Fold=Sim_Fold;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function FolderName_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FolderName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in simulate.
function simulate_Callback(hObject, eventdata, handles)
% hObject    handle to simulate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of simulate
Sim_Len=str2double(get(handles.SimLength,'String'));
dts=handles.dti;
Vc=handles.CSi;
Folder_Name=handles.FNi;
TOL=handles.TOLi;
% delete([Folder_Name '\DB*.mat']);
Z1=str2double(get(handles.Z1,'String'));
Z2=str2double(get(handles.Z2,'String'));
Z3=str2double(get(handles.Z3,'String'));
WFR_m=[Z1;Z2;Z3];
GN=handles.GNi;
TAB_MAIN_program(Sim_Len,WFR_m);
tf=Sim_Len/(Vc*1000/60);
tt=0:dts:tf;
nDBe=length(tt);
nDBe=150;
ct=0;
while get(hObject,'Value')
    if ct==0
    a=dir([Folder_Name handles.path_sep 'DB*.mat']);
    ndb=size(a,1);
    ct=ct+1;
    end
%     
    t=tt(ndb+1)
    TAB_Iterative_Solver2(Folder_Name,dts,TOL,t,ndb-1)
    ndb=ndb+1;
% if ndb>nDBe
%     break
%  end
end





% --- Executes on selection change in Prop_Sel.
function Prop_Sel_Callback(hObject, eventdata, handles)
% hObject    handle to Prop_Sel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Prop_Sel contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Prop_Sel
handles.vprop=get(hObject,'Value');
guidata(hObject, handles);
plot_material_data(handles);
% --- Executes during object creation, after setting all properties.
function Prop_Sel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Prop_Sel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function User_Set_Callback(hObject, eventdata, handles)
% hObject    handle to User_Set (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of User_Set as text
%        str2double(get(hObject,'String')) returns contents of User_Set as a double


% --- Executes during object creation, after setting all properties.
function User_Set_CreateFcn(hObject, eventdata, handles)
% hObject    handle to User_Set (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton63.
function pushbutton63_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton63 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in Browse.
function Browse_Callback(hObject, eventdata, handles)
% hObject    handle to Browse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in List_Set.
function List_Set_Callback(hObject, eventdata, handles)
% hObject    handle to List_Set (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns List_Set contents as cell array
%        contents{get(hObject,'Value')} returns selected item from List_Set
dval=get(hObject,'Value');
List=get(hObject,'String');
set(handles.User_Set,'String',List{dval});
handles.Folders.File.SUse=List{dval};
handles.user=load([handles.Folders.Main_Fold.Set handles.path_sep List{dval} '.mat']); % file loads the variables
guidata(hObject,handles);
handles.default.Grade_depend.Grade_Number=handles.Grade_depend.Grade_Number;
handles.default.Grade_depend.Grade_Type=handles.Grade_depend.Grade_Type;
handles.user.Grade_depend.Grade_Number=handles.Grade_depend.Grade_Number;
handles.user.Grade_depend.Grade_Type=handles.Grade_depend.Grade_Type;
guidata(hObject,handles);
handles.present=handles.user;
guidata(hObject,handles);


% if get(handles.mat_off,'Value')==0 % material dependent
%     handles.present.Grade_depend=handles.Grade_depend;
% end
% guidata(hObject,handles);
% [SEC_COOL,Nozzle_Tab,No_Noz_Zo,cooling_info]=Updade_Cooling_Info(handles);
% handles.present.SEC_COOL=SEC_COOL;
% handles.present.No_Noz_Zo=No_Noz_Zo;
% handles.present.Nozzle_Tab=Nozzle_Tab;
% handles.present.cooling_info=cooling_info;
% guidata(hObject,handles);

set(handles.Z1,'String',handles.present.Grade_depend.WFR_Lm3(1));
set(handles.Z2,'String',handles.present.Grade_depend.WFR_Lm3(2));
set(handles.Z3,'String',handles.present.Grade_depend.WFR_Lm3(3));
set(handles.CM,'String',handles.present.Grade_depend.Cool_Mode);
set(handles.CS,'String',handles.present.Grade_depend.Casting_speed);
set(handles.CT,'String',handles.present.Grade_depend.Init_Temp);
handles.C_DB=0;handles.C_DB2=0;
set(handles.RunSim,'Value',0);

guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function List_Set_CreateFcn(hObject, eventdata, handles)
% hObject    handle to List_Set (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in mat_on.
function mat_on_Callback(hObject, eventdata, handles)
% hObject    handle to mat_on (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of mat_on


% --- Executes on button press in mat_off.
function mat_off_Callback(hObject, eventdata, handles)
% hObject    handle to mat_off (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of mat_off


% --- Executes when selected object is changed in Mat_Depend.
function Mat_Depend_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in Mat_Depend 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
if (hObject == handles.mat_off)% mat_off=1, mat_on=0 % material independent
    handles.present.Grade_depend=handles.user.Grade_depend;
end
if (hObject == handles.mat_on)% mat_off=0, mat_on=1 % material dependent
    handles.present.Grade_depend=handles.Grade_depend;
end
guidata(hObject,handles);
WFR_Lm3=handles.present.Grade_depend.WFR_Lm3;
set(handles.Z1,'String',WFR_Lm3(1));
set(handles.Z2,'String',WFR_Lm3(2));
set(handles.Z3,'String',WFR_Lm3(3));
set(handles.CM,'String',handles.present.Grade_depend.Cool_Mode);
set(handles.CS,'String',handles.present.Grade_depend.Casting_speed);
set(handles.CT,'String',handles.present.Grade_depend.Init_Temp);


% [SEC_COOL,Nozzle_Tab,No_Noz_Zo,cooling_info]=Updade_Cooling_Info(handles);
% handles.present.SEC_COOL=SEC_COOL;
% handles.present.No_Noz_Zo=No_Noz_Zo;
% handles.present.Nozzle_Tab=Nozzle_Tab;
% handles.present.cooling_info=cooling_info;
% SEC_COOL.HTC(50,:)
% figure();
% plot(HTC(50,:));hold on
guidata(hObject,handles);



% function[Z,X,Q]=Secondary_Cooling_flow(Wi,NozzTab)
% ner    = 100;        % no.of elements in x direction
% nes    = 100;        % no.of elements in y direction
% nozz_posS=NozzTab(:,[1 2]);
% nozz_angS=NozzTab(:,5);
% nozz_disS=NozzTab(:,4);
% Q_nozS=NozzTab(:,3);
% sec_L=nozz_posS(size(nozz_posS,1),1)+1000;
% dz  = 5;
% z   = 0:dz:sec_L;  % z - vector
% x   = 0:Wi/ner:Wi;      % x - vector
% [Z,X,Q]=assemble_spray_matrix(z,x,nozz_posS,nozz_angS,nozz_disS,Q_nozS);


function[Z,X,Q]=Secondary_Cooling_flow(Wi,NozzTab)
ner    = 100;        % no.of elements in x direction
nes    = 100;        % no.of elements in y direction
nozz_posS=NozzTab(:,[1 2]);
nozz_angS=NozzTab(:,6);
nozz_disS=NozzTab(:,5);
Q_nozS=NozzTab(:,[3 4]);
spray_rad=NozzTab(:,[8 9]);
ICFM=NozzTab(:,[10 11]);
sec_L=nozz_posS(size(nozz_posS,1),1)+1000;
dz  = 5;
z   = 0:dz:sec_L;  % z - vector
x   = 0:Wi/ner:Wi;      % x - vector
% [Z,X,Q]=assemble_spray_matrix(z,x,nozz_posS,nozz_angS,nozz_disS,Q_nozS);
[Z,X,Q]=assemble_spray_matrix(z,x,nozz_posS,nozz_angS,nozz_disS,Q_nozS,spray_rad,ICFM);




% function[Z,X,Q]=assemble_spray_matrix(z,x,nozz_posS,nozz_angS,nozz_disS,Q_nozS)
% [Z,X]=meshgrid(z,x);
% Q=zeros(size(Z));
% ne=length(nozz_angS);
% for i=1:ne
%     Qe=[];zmap=[];xmap=[];
%     nozz_ang=nozz_angS(i);
%     nozz_pos=nozz_posS(i,:);
%     nozz_dis=nozz_disS(i);
%     Q_noz=Q_nozS(i);
%     R   = nozz_dis*tan(nozz_ang*pi/360);
%     LRTB=find_square_matrix_spray(nozz_pos,R,z,x);
%     zmap=LRTB(1):LRTB(2);
%     xmap=LRTB(4):LRTB(3);
%     Z1=Z(xmap,zmap);
%     X1=X(xmap,zmap);
%     Q_noz=Q_noz/(pi*(R*1e-3)^2);
%     [Qe]=Spray_dist_matrix(nozz_dis,nozz_ang,Z1,X1,nozz_pos,Q_noz);
%     Q(xmap,zmap)=Q(xmap,zmap)+Qe;
% end


function[Z,X,Q]=assemble_spray_matrix(z,x,nozz_posS,nozz_angS,nozz_disS,Q_nozS,spray_rad,ICFM)
[Z,X]=meshgrid(z,x);
Q=zeros(size(Z));
ne=length(nozz_angS);
for i=1:ne
    Qe=[];zmap=[];xmap=[];
    nozz_ang=nozz_angS(i);
    nozz_pos=nozz_posS(i,:);
    nozz_dis=nozz_disS(i);
    Q_noz=Q_nozS(i,:);
    %     R   = nozz_dis*tan(nozz_ang*pi/360);
    R2   = spray_rad(i,2);R1   = spray_rad(i,1);
    ICFMi=ICFM(i,:);
    LRTB=find_square_matrix_spray(nozz_pos,R2,z,x);
    zmap=LRTB(1):LRTB(2);
    xmap=LRTB(4):LRTB(3);
    Z1=Z(xmap,zmap);
    X1=X(xmap,zmap);
    %     Q_noz=Q_noz/(pi*(R*1e-3)^2);
    %     [Qe]=Spray_dist_matrix(nozz_dis,nozz_ang,Z1,X1,nozz_pos,Q_noz);
    [Qe]=Spray_dist_matrix(nozz_dis,nozz_ang,Z1,X1,nozz_pos,Q_noz,[R1 R2],ICFMi);
    Q(xmap,zmap)=Q(xmap,zmap)+Qe;
%     contour(Z1,X1,Qe)
%     pause
end

% function[LRTB]=find_square_matrix_spray(nozz_pos,R,z,x)
% %%% nozzle square position in global matrix
% L1=nozz_pos(1)-R; %%% left line
% R1=nozz_pos(1)+R;
% T1=nozz_pos(2)+R;
% B1=nozz_pos(2)-R;
% if T1>max(x)
%     T1=max(x);
% end
% if B1<min(x)
%     B1=min(x);
% end
% if L1<min(z)
%     L1=min(z);
% end
% 
% LRTB_coord=[L1 R1 T1 B1];
% LRTB=[];cntid2=0;
% for i=1:length(LRTB_coord)
%     La=LRTB_coord(i);A=[];
%     if i<=2
%         A=z;        
%     else
%         A=x;
%     end
%     if i==1 || i==4
%         pos=1;
%     else
%         pos=2;
%     end
%     cntid=0;Nid=0;a=[];b=[];
%     while cntid==0
%         a=A(Nid+1);
%         b=A(Nid+2);
%         if ((a-La)<=0) && ((b-La)>=0)
%             LRTB(cntid2+1)=Nid+pos;
%             cntid=cntid+1;cntid2=cntid2+1;
%         end
%         Nid=Nid+1;
%     end
% end

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

% function[Qe]=Spray_dist_matrix(nozz_dis,nozz_ang,Z1,X1,nozz_pos,Q_noz)
% nozz_ang=nozz_ang.*pi/180;
% R   = nozz_dis*tan(nozz_ang/2);
% Qe=zeros(size(Z1));aa=[];
% for row=1:size(Z1,1)
%     for col=1:size(Z1,2)
%         za=Z1(row,col);
%         xa=X1(row,col);
%         r=sqrt((nozz_pos(1)-za)^2+(nozz_pos(2)-xa)^2);       
%         if r<=R            
%             Qe(row,col)=Qe(row,col)+0.5*Q_noz*tan(nozz_ang/2)^2/(1-cos(nozz_ang/2))*(1/(1+r^2/nozz_dis^2)^1.5);
%         end
%     end
% end

function[Qe]=Spray_dist_matrix(nozz_dis,nozz_ang,Z1,X1,nozz_pos,Q_noz,Radius,ICFMi)
nozz_ang=nozz_ang.*pi/180;
R   = nozz_dis*tan(nozz_ang/2);
R1=Radius(1);R2=Radius(2);
Qe=zeros(size(Z1));aa=[];
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
%             Qe(row,col)=Qe(row,col)+0.5*Q_noz*tan(nozz_ang/2)^2/(1-cos(nozz_ang/2))*(1/(1+r^2/nozz_dis^2)^1.5);
        end
    end
end


% function[POS_NOZ,WAT_FLO,NOZ_DIS,CON_ANG,No_Noz_Zo,Spray_rad]=create_nozzle_para(cooling_info)
% %%% create individual nozzle information
% Zones={cooling_info.Zones};
% Nozzles=convert_cell_matrix({cooling_info.Nozzles});
% Zone_pos=convert_cell_matrix({cooling_info.Zone_pos});
% WFR=convert_cell_matrix({cooling_info.WFR});
% ND=convert_cell_matrix({cooling_info.ND});
% CA=convert_cell_matrix({cooling_info.CA});
% 
% NNP=find_num_nozz_per_zone(Nozzles);
% WFR=WFR./NNP;
% 
% POS_NOZ=[];% final nozzle positions
% WAT_FLO=[];% Nozzle water flow rate
% NOZ_DIS=[];% Nozzle distance
% CON_ANG=[];% Cone angle
% No_Noz_Zo=[];Spray_rad=[];
% ct=0;
% for i=1:length(Zones)
%     nRows=Nozzles(i,1);
%     if nRows>1
%         nRow=nRows-1;
%     else
%         nRow=nRows;
%     end
%     nCol=Nozzles(i,2);
%     axial_1=Zone_pos(i,1);
%     axial_2=Zone_pos(i,2);
%     lateral=[Zone_pos(i,3) Zone_pos(i,4)];
%     dist=(axial_2-axial_1);
%     WFRi=WFR(i);NDi=ND(i);CAi=CA(i);
%     for j=1:nRows
%         if lateral(2)>0
%             nCol=2;
%         else
%             nCol=1;
%         end
%         for k=1:nCol
%             POS_NOZ(ct+1,:)=[axial_1+(j-1)*dist/nRow, lateral(k)];
%             WAT_FLO(ct+1)=WFRi;NOZ_DIS(ct+1)=NDi;CON_ANG(ct+1)=CAi;
%             Spray_rad(ct+1)=NDi*tan(CAi*pi/360);
%             ct=ct+1;
%         end
%     end
%     No_Noz_Zo(i)=ct;
% end


function[POS_NOZ,WAT_FLO,NOZ_DIS,CON_ANG,No_Noz_Zo,Spray_rad,ICFM]=create_nozzle_para(cooling_info)
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
length(Zones)
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
%             CON_ANG(ct+1)=CAi;
            
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
% [r1 r2]
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


function[]=create_intial_DB(Init_Temp,XYH,MAPH,FN,handles)
neH=size(MAPH,1);nnH=size(XYH,1);
nIPH=9;
T = zeros(nnH,1)+ Init_Temp;        % Nodal Temperatures (?C)
IPTEMP = zeros(neH,nIPH) + Init_Temp ;
IPPH   = zeros(neH,nIPH)+1  ;     % [Liquid fraction]
IPDT   = zeros(neH,nIPH);         % Temperature icrement
IPDQ   = zeros(neH,nIPH);         % Generated letant heat
t=0;
Axial=0;
FN=[FN handles.path_sep 'DB0.mat'];
save(FN,'T','IPTEMP','IPPH','IPDT','IPDQ','t','Axial');


function[Grade_depend]=Get_Grade_depend(Steel_reciepes_SC,Std_cooling,dvalue)
GNr={Steel_reciepes_SC.Grade_Number};Grade_Number=GNr{dvalue};
GTy={Steel_reciepes_SC.Grade_Type};Grade_Type=GTy{dvalue};
T_in={Steel_reciepes_SC.T_init};Init_Temp=T_in{dvalue};
Cspeed={Steel_reciepes_SC.Vc_target};Casting_speed=Cspeed{dvalue};
CoMos={Steel_reciepes_SC.Cooling_mode};cool_mode=CoMos{dvalue};
[Cool_Mode,Cool_Mode_Value]=cool_mode_sel(cool_mode);
Water_flow_rate={Steel_reciepes_SC.Water_flow_rate};
if Cool_Mode_Value<5
    WFR_Lm3=Std_cooling(Cool_Mode_Value,:);
else
    WF_matrix_5x3=zeros(5,3);
    WF_matrix_5x3(1,1)=1;
    WF_matrix_5x3(2,2)=1;WF_matrix_5x3(3,2)=1;
    WF_matrix_5x3(4,3)=1;WF_matrix_5x3(5,3)=1;
    
    WFR_Lm3=round(Water_flow_rate{dvalue}*WF_matrix_5x3*10)/10;
end

Grade_depend=struct(...
    'Cool_Mode',      {Cool_Mode},...
    'Cool_Mode_Value',{Cool_Mode_Value},...
    'WFR_Lm3', {WFR_Lm3},...
    'Casting_speed',{Casting_speed},...
    'Init_Temp',{Init_Temp},...
    'Grade_Number',Grade_Number,...
    'Grade_Type',Grade_Type);



% --- Executes on button press in pushbutton67.
function pushbutton67_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton67 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in togglebutton38.
function togglebutton38_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton38 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebutton38


% --- Executes on button press in togglebutton39.
function togglebutton39_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton39 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebutton39


% --- Executes on button press in togglebutton40.
function togglebutton40_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton40 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebutton40


% --- Executes on button press in togglebutton41.
function togglebutton41_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton41 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebutton41


% --- Executes on button press in pushbutton68.
function pushbutton68_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton68 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in Advanced_Run.
function Advanced_Run_Callback(hObject, eventdata, handles)
% hObject    handle to Advanced_Run (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
FN=get(handles.FolderName,'String');
Folder_Name1=[handles.Folders.Main_Fold.Sim  handles.Folders.Main_Fold.Sim];
Folder_Name=[handles.Folders.Main_Fold.Sim handles.path_sep FN];
d=dir(handles.Folders.Main_Fold.Sim);
str={d.name};
if length(str)>2
    str=str(3:length(str));
end
Exist=strcmp(str,FN);
val=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if sum(Exist)>0 && get(hObject, 'Value') == 1
    %     errordlg('Steel Grade Simulation Already Exists. Change Folder Name', 'Grade Simulation Folder Already Exist');
    val=2;
    choice = questdlg('Do you want to overwrite existing Simulations???', ...
        'Atttention', ...
        'YES','NO','Continue Simulation','NO');

    switch choice
        case 'YES'
            val=1;
        case 'NO' 
            val=2;            
        case 'Continue Simulation' 
            val=3;
        case ''
            val=2;
    end
else % if folder doesnot exist
    val=4;
    if isdir(Folder_Name)==0
        mkdir(Folder_Name);
        handles.C_DB2=0;
    end
end
if val==3 &&  isequal(handles.present_previous,handles.present)==0

    choice2 = questdlg('Simulations Data moidified: simulation can not be continued: Start Overwriting?', ...
        'Atttention', ...
        'YES','NO','NO');
    switch choice2
        case 'YES'
            val=1;
            handles.C_DB2=0;guidata(hObject, handles);
        case 'NO'
            val=2;
        case ''
            val=2;
    end
end

if val==1 || val==3 || val==4
    %reset 
    handles.C_DB=0;
    set(handles.status_text,'String','');
    h_axes=handles.axes_status;set(h_axes,'Visible','Off');
    cla(handles.axes_status);axes(handles.axes_status);
    patch([0,0,0,0],[0,0,0,0],'w');axis([0,1,0,1]);axis off;
    set(handles.RunSim,'Value',0);
    
        if val==1
            delete([Folder_Name handles.path_sep 'DB*.mat']);
            handles.C_DB2=0;
        elseif val==3
            pa=[Folder_Name handles.path_sep];
            delete_data_file(pa,handles.C_DB2);
        end
    handles.present_previous=handles.present;    
   guidata(hObject, handles);

    %%% Material --------------------------------------------------------------
    KHT=handles.Properties.KHT;
    SHT=handles.Properties.SHT;
    RHT=handles.Properties.RHT;
    PCF=handles.Properties.PCF;
    ENTH=handles.Properties.ENTH;
    Latent=handles.Properties.Latent;
    TS=handles.Properties.TS;
    TL=handles.Properties.TL;
    
%     [SEC_COOL,Nozzle_Tab,No_Noz_Zo,cooling_info]=Updade_Cooling_Info(handles,3);
%     handles.present.SEC_COOL=SEC_COOL;
%     handles.present.No_Noz_Zo=No_Noz_Zo;
%     handles.present.Nozzle_Tab=Nozzle_Tab;
%     handles.present.cooling_info=cooling_info;
%     guidata(hObject, handles);
    
    m=handles.present;
    
    XYH=m.XYH;
    MAPH=m.MAPH;
%     BC_Info=m.BC_Info;
%     [BC_Info2]=find_reduction_factor(XYH,MAPH,BC_Info,m.f1,m.f2,m.Wi,m.Th);
%     m.BC_Info2=BC_Info2;
%     handles.present=m;
    
    save([Folder_Name handles.path_sep 'Material.mat'],'KHT','SHT','RHT','TS','TL','PCF','ENTH',...
        'Latent');
    save([Folder_Name handles.path_sep 'Process_Parameters.mat'],'-struct','m');
    Init_Temp=handles.present.Grade_depend.Init_Temp;

    create_intial_DB(Init_Temp,XYH,MAPH,Folder_Name,handles);
    
    Vc=m.Grade_depend.Casting_speed/60;
    dx=m.dt*Vc;
    NDB=round(m.Sim_Len*1e-3/dx);
    
    %     Iterative_Solver(Folder_Name);
    
    h_axes_status=handles.advanced_status;set(h_axes_status,'Visible','On');
    cla(handles.advanced_status);axes(handles.advanced_status);axis off;
    C_DB=handles.C_DB2;
                pr=(C_DB+1)/(NDB+1)*100;
            patch([0,1,1,0],[0,0,1,1],'w');
            patch([0,pr/100,pr/100,0],[0,0,1,1],'r');
            axis([0,1,0,1]);axis off;drawnow;
            
    if get(hObject, 'Value') == 1
        set(hObject, 'String','STOP Simulation');
        while C_DB<NDB+1
            Axial=Solver_new(Folder_Name,C_DB,handles.path_sep);
%             Axial=Solver_enthalpy(Folder_Name,C_DB,handles.path_sep);
%             Axial=C_DB
            pause(1e-3);
            pr=(C_DB+1)/(NDB+1)*100;
            patch([0,1,1,0],[0,0,1,1],'w');
            patch([0,pr/100,pr/100,0],[0,0,1,1],'r');
            axis([0,1,0,1]);axis off;drawnow;
            str=[ num2str(round(pr*1)/1) '   % Completed : ' num2str(round(Axial*100)/100) '   m reached'];
            % str
            set(handles.advanced_status_text,'String',str);
            if get(hObject, 'Value') == 0
                set(hObject, 'String','RUN Simulation');
                handles.C_DB2=handles.C_DB2+1;
                axis off;
                guidata(hObject, handles);
                break;
            end
            C_DB=C_DB+1;
            handles.C_DB2=handles.C_DB2+1;
            axis off;
            guidata(hObject, handles);
        end
        if C_DB==NDB+1
            set(hObject, 'Value',0);
            set(hObject, 'String','RUN Simulation');
        end
            
    else
        set(hObject, 'String','RUN Simulation');
    end
    axis off;
    %handles.present_previous=handles.present;
else
    set(hObject, 'Value',0);
    set(hObject, 'String','RUN Simulation');
end


 guidata(hObject, handles);


% if sum(Exist)>0
%     %     errordlg('Steel Grade Simulation Already Exists. Change Folder Name', 'Grade Simulation Folder Already Exist');
%     val=2;
%     choice = questdlg('Do you want to overwrite existing Simulations???', ...
%         'Atttention', ...
%         'YES','NO','NO');
%     switch choice
%         case 'YES'
%             val=1;
%         case 'NO'
%             val=2;
%             
%     end
% else
%     val=3;
% end
% if val==3 || val==1
%     
%     if isdir(Folder_Name)==0
%         mkdir(Folder_Name);
%     else
%         delete([Folder_Name handles.path_sep 'DB*.mat']);
%     end
%     %%% Material --------------------------------------------------------------
%     KHT=handles.Properties.KHT;
%     SHT=handles.Properties.SHT;
%     RHT=handles.Properties.RHT;
%     PCF=handles.Properties.PCF;
%     ENTH=handles.Properties.ENTH;
%     Latent=handles.Properties.Latent;
%     TS=handles.Properties.TS;
%     TL=handles.Properties.TL;
%     m=handles.present;
%     save([Folder_Name handles.path_sep 'Material.mat'],'KHT','SHT','RHT','TS','TL','PCF','ENTH',...
%         'Latent');
%     save([Folder_Name handles.path_sep 'Process_Parameters.mat'],'-struct','m');
%     
%     Init_Temp=handles.present.Grade_depend.Init_Temp;
%     XYH=handles.present.XYH;
%     MAPH=handles.present.MAPH;
%     create_intial_DB(Init_Temp,XYH,MAPH,Folder_Name,handles);
% %     Iterative_Solver(Folder_Name);
% end



% --- Executes on button press in Advanced_Save.
function Advanced_Save_Callback(hObject, eventdata, handles)
% hObject    handle to Advanced_Save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Set_File=get(handles.User_Set,'String');
Def_File_Name=handles.Folders.File.SDef;
% str2='.mat';
% Def_File_Name=strrep(Def_File_Name,str2,'');
if strcmp(Def_File_Name,Set_File)==1
  Set_File2=[handles.Folders.Main_Fold.Set  'User_Sett_1.mat'];
else
    Set_File2=[handles.Folders.Main_Fold.Set  Set_File '.mat'];
end
[FName, PName]=uiputfile(Set_File2,'Save file name');
User_Struct=handles.present;
if isequal(FName,0) || isequal(PName,0)
else
    save([PName FName],'-struct','User_Struct');
end


% --- Executes on button press in Mold_UserSet.
function Mold_UserSet_Callback(hObject, eventdata, handles)
% hObject    handle to Mold_UserSet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.present.Grade_depend.Casting_speed=handles.user.Grade_depend.Casting_speed;
handles.present.Grade_depend.Init_Temp=handles.user.Grade_depend.Init_Temp;
handles.present.Tw=handles.user.Tw;
handles.present.Tatm=handles.user.Tatm;
handles.present.mold_wall=handles.user.mold_wall;
handles.present.water_temp_inc=handles.user.water_temp_inc;
handles.present.Qwater=handles.user.Qwater;
handles.present.Mol_Len=handles.user.Mol_Len;
handles.present.index=handles.user.index;
handles.present.MHF=handles.user.MHF;


guidata(hObject,handles);
cla(handles.Mold_Plot);
axes(handles.Mold_Plot);
% [handles.present.z,handles.present.hf]=mold_heat_flux(handles.present);
plot(handles.present.MHF(:,1),handles.present.MHF(:,2),'-b','LineWidth',3);hold on
set(gca,'FontSize',12,'LineWidth',2.0,'Box','off');
ylabel(' Heat Flux [MW/m^2]','FontSize',13);
xlabel(' Distance from Meniscus [mm]','FontSize',13);
Mold_Panel_Reset(handles);
guidata(hObject,handles);

% --- Executes on button press in Mold_ApplyDefault.
function Mold_ApplyDefault_Callback(hObject, eventdata, handles)
% hObject    handle to Mold_ApplyDefault (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.present.Grade_depend.Casting_speed=handles.default.Grade_depend.Casting_speed;
handles.present.Grade_depend.Init_Temp=handles.default.Grade_depend.Init_Temp;
handles.present.Tw=handles.default.Tw;
handles.present.Tatm=handles.default.Tatm;
handles.present.mold_wall=handles.default.mold_wall;
handles.present.water_temp_inc=handles.default.water_temp_inc;
handles.present.Qwater=handles.default.Qwater;
handles.present.Mol_Len=handles.default.Mol_Len;
handles.present.index=handles.default.index;
handles.present.MHF=handles.user.MHF;


guidata(hObject,handles);
cla(handles.Mold_Plot);
axes(handles.Mold_Plot);
% [handles.present.z,handles.present.hf]=mold_heat_flux(handles.present);
plot(handles.present.MHF(:,1),handles.present.MHF(:,2),'-b','LineWidth',3);hold on
set(gca,'FontSize',12,'LineWidth',2.0,'Box','off');
ylabel(' Heat Flux [MW/m^2]','FontSize',13);
xlabel(' Distance from Meniscus [mm]','FontSize',13);
Mold_Panel_Reset(handles);
guidata(hObject,handles);




function MML_u_Callback(hObject, eventdata, handles)
% hObject    handle to MML_u (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MML_u as text
%        str2double(get(hObject,'String')) returns contents of MML_u as a double
handles.present.Mol_Len = str2double(get(hObject,'String'));
guidata(hObject,handles);
cla(handles.Mold_Plot);
axes(handles.Mold_Plot);
[MHF]=mold_heat_flux(handles.present);
plot(MHF(:,1),MHF(:,2),'-b','LineWidth',3);hold on
set(gca,'FontSize',12,'LineWidth',2.0,'Box','off');
ylabel(' Heat Flux [MW/m^2]','FontSize',13);
xlabel(' Distance from Meniscus [mm]','FontSize',13);
handles.present.MHF=MHF;
guidata(hObject,handles);
% --- Executes when selected object is changed in Advanced_Mold_Modelsel.
function Advanced_Mold_Modelsel_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in Advanced_Mold_Modelsel 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
if (hObject == handles.MB)
handles.present.index=1;
elseif (hObject == handles.MH)
handles.present.index=2;
elseif (hObject == handles.MA)
handles.present.index=3;
else
handles.present.index=4;
end
m=handles.present.index;
guidata(hObject,handles);
cla(handles.Mold_Plot);
axes(handles.Mold_Plot);
[MHF]=mold_heat_flux(handles.present);
plot(MHF(:,1),MHF(:,2),'-b','LineWidth',3);hold on
set(gca,'FontSize',12,'LineWidth',2.0,'Box','off');
ylabel(' Heat Flux [MW/m^2]','FontSize',13);
xlabel(' Distance from Meniscus [mm]','FontSize',13);
handles.present.MHF=MHF;
guidata(hObject,handles);


function[]=Mold_Panel_Reset(handles)
    set(handles.MCS_u,'String',handles.present.Grade_depend.Casting_speed);
    set(handles.MCT_u,'String',handles.present.Grade_depend.Init_Temp);
    set(handles.MCWT_u,'String',handles.present.Tw);
    set(handles.MAMT_u,'String',handles.present.Tatm);
    set(handles.MML_u,'String',handles.present.Mol_Len);
        
    set(handles.MDT_u,'Style','text','BackgroundColor',0.941*[1 1 1]);
    set(handles.MMWF_u,'Style','text','BackgroundColor',0.941*[1 1 1]);
    set(handles.MMWT_u,'Style','text','BackgroundColor',0.941*[1 1 1]);
    
    set(handles.MA,'Value',0);
    set(handles.MH,'Value',0);
    set(handles.MB,'Value',0);
    set(handles.MO,'Value',0);
    a=handles.present.index;
    if handles.present.index==1
        set(handles.MB,'Value',1);        
    elseif handles.present.index==2
        set(handles.MH,'Value',1);        
    elseif handles.present.index==3
        set(handles.MA,'Value',1);
        set(handles.MDT_u,'Style','edit','BackgroundColor',[1 1 1]);
        set(handles.MMWF_u,'Style','edit','BackgroundColor',[1 1 1]);
        set(handles.MMWT_u,'Style','edit','BackgroundColor',[1 1 1]);
    else
        set(handles.MO,'Value',1);
    end
    set(handles.MDT_u,'String',handles.present.water_temp_inc);
    set(handles.MMWF_u,'String',handles.present.Qwater);
    set(handles.MMWT_u,'String',handles.present.mold_wall);

% --- Executes during object creation, after setting all properties.
function MML_u_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MML_u (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function MAMT_u_Callback(hObject, eventdata, handles)
% hObject    handle to MAMT_u (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MAMT_u as text
%        str2double(get(hObject,'String')) returns contents of MAMT_u as a double
handles.present.Tatm=str2double(get(hObject,'String'));
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function MAMT_u_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MAMT_u (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




function MCWT_u_Callback(hObject, eventdata, handles)
% hObject    handle to MCWT_u (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MCWT_u as text
%        str2double(get(hObject,'String')) returns contents of MCWT_u as a double
handles.present.Tw=str2double(get(hObject,'String'));
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function MCWT_u_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MCWT_u (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function MCT_u_Callback(hObject, eventdata, handles)
% hObject    handle to MCT_u (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MCT_u as text
%        str2double(get(hObject,'String')) returns contents of MCT_u as a double
handles.present.Grade_depend.Init_Temp=str2double(get(hObject,'String'));
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function MCT_u_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MCT_u (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




function[MHF]=mold_heat_flux(m)
Cast_Speeds=m.Grade_depend.Casting_speed*1000/60;
z=0:30:m.Mol_Len;
Qwater=m.Qwater*1e-3/60;


if m.index==1 % Brimbacombe
    hf=2.21e6-2.67e5*sqrt(z./Cast_Speeds);
elseif m.index==2 % Ha
    t=z./Cast_Speeds;
    hf=0.07128*exp(-t)+2.2328*exp(-t./9.5)+0.698;
    hf=hf.*1e6;
elseif m.index==3 % Ali
%     hf=Alizadeh(Wi,Mol_Len,z);
    mold_alpha=1.25;
%     mold_Pm=4*m.Wi*1e-3;
%     mold_Pm=(2*(m.Wi+m.mold_wall)+2*(m.Th+m.mold_wall))*1e-3;
%     mold_Pm=2*m.Wi*(m.fact2+1)*1e-3;
%      mold_Pm=(2*(m.Wi*(m.f1+1)+2*m.d_chamfer*(m.f1-3))+4*(m.f2*(3*m.d_chamfer-m.Wi)+m.f1*(m.Wi-m.d_chamfer)))*1e-3;
    mold_Pm=((2*(m.Wi-2*m.d_chamfer)*(m.f1+1))+(4*m.d_chamfer*(m.f1+m.f2)))*1e-3;
hf=mold_alpha*1.35*m.water_temp_inc*4186*1000*Qwater*(1/mold_Pm)*(1/(1-exp(-mold_alpha*m.Mol_Len*1e-3)))*exp(-mold_alpha.*z.*1e-3);
else % ovgu
    zhf=[0   3.5
   0.02 3.4
   0.03 3.3
   0.06 3.1
   0.07 3.0
   0.08 2.93
   0.1 2.827
   0.1088 2.767
   0.118  2.7
   0.13   2.6
   0.15   2.426
   0.16  2.4
   0.19  2.2
   0.23   2.071
   0.2585 1.981
   0.282  1.893
   0.3088 1.802
   0.3312 1.734
   0.3636  1.647
   0.4078 1.55
   0.4463 1.481
   0.4854 1.426
   0.5029 1.406
   0.5576 1.354
   0.6049 1.323
   0.6737 1.29
   0.722  1.27
   0.7795 1.245
   0.8351 1.212
   0.8995 1.154
   1      0.9937];
z=zhf(:,1).*1e3;z=z';
hf=zhf(:,2).*1e6;hf=hf';
end
 MHF=[z' (hf./1e6)'];
 
 
function MCS_u_Callback(hObject, eventdata, handles)
% hObject    handle to MCS_u (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MCS_u as text
%        str2double(get(hObject,'String')) returns contents of MCS_u as a double
handles.present.Grade_depend.Casting_speed=str2double(get(hObject,'String'));
guidata(hObject,handles);
cla(handles.Mold_Plot);
axes(handles.Mold_Plot);
[MHF]=mold_heat_flux(handles.present);
plot(MHF(:,1),MHF(:,2),'-b','LineWidth',3);hold on
set(gca,'FontSize',12,'LineWidth',2.0,'Box','off');
ylabel(' Heat Flux [MW/m^2]','FontSize',13);
xlabel(' Distance from Meniscus [mm]','FontSize',13);
handles.present.MHF=MHF;
guidata(hObject,handles);



% --- Executes during object creation, after setting all properties.
function MCS_u_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MCS_u (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




function MDT_u_Callback(hObject, eventdata, handles)
% hObject    handle to MDT_u (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MDT_u as text
%        str2double(get(hObject,'String')) returns contents of MDT_u as a double
handles.present.water_temp_inc=str2double(get(hObject,'String'));
guidata(hObject,handles);
cla(handles.Mold_Plot);
axes(handles.Mold_Plot);
[MHF]=mold_heat_flux(handles.present);
plot(MHF(:,1),MHF(:,2),'-b','LineWidth',3);hold on
set(gca,'FontSize',12,'LineWidth',2.0,'Box','off');
ylabel(' Heat Flux [MW/m^2]','FontSize',13);
xlabel(' Distance from Meniscus [mm]','FontSize',13);
handles.present.MHF=MHF;
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function MDT_u_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MDT_u (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function MMWT_u_Callback(hObject, eventdata, handles)
% hObject    handle to MMWT_u (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MMWT_u as text
%        str2double(get(hObject,'String')) returns contents of MMWT_u as a double
handles.present.mold_wall=str2double(get(hObject,'String'));
guidata(hObject,handles);
cla(handles.Mold_Plot);
axes(handles.Mold_Plot);
[MHF]=mold_heat_flux(handles.present);
plot(MHF(:,1),MHF(:,2),'-b','LineWidth',3);hold on
set(gca,'FontSize',12,'LineWidth',2.0,'Box','off');
ylabel(' Heat Flux [MW/m^2]','FontSize',13);
xlabel(' Distance from Meniscus [mm]','FontSize',13);
handles.present.MHF=MHF;
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function MMWT_u_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MMWT_u (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





function MMWF_u_Callback(hObject, eventdata, handles)
% hObject    handle to MMWF_u (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MMWF_u as text
%        str2double(get(hObject,'String')) returns contents of MMWF_u as a double
handles.present.Qwater=str2double(get(hObject,'String'));
guidata(hObject,handles);
cla(handles.Mold_Plot);
axes(handles.Mold_Plot);
[MHF]=mold_heat_flux(handles.present);
plot(MHF(:,1),MHF(:,2),'-b','LineWidth',3);hold on
set(gca,'FontSize',12,'LineWidth',2.0,'Box','off');
ylabel(' Heat Flux [MW/m^2]','FontSize',13);
xlabel(' Distance from Meniscus [mm]','FontSize',13);
handles.present.MHF=MHF;
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function MMWF_u_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MMWF_u (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function Panel_Advanced_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Panel_Advanced (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on selection change in Secondary_Flow_or_HTC.
function Secondary_Flow_or_HTC_Callback(hObject, eventdata, handles)
% hObject    handle to Secondary_Flow_or_HTC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Secondary_Flow_or_HTC contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Secondary_Flow_or_HTC
cla(handles.Secondary_plot);
axes(handles.Secondary_plot);
plot_Secondary_Cooling(handles);


% --- Executes during object creation, after setting all properties.
function Secondary_Flow_or_HTC_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Secondary_Flow_or_HTC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Second_ApplyDefault.
function Secondary_ApplyDefault_Callback(hObject, eventdata, handles)
% hObject    handle to Second_ApplyDefault (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in Second_UserSet.
function Secondary_UserSet_Callback(hObject, eventdata, handles)
% hObject    handle to Second_UserSet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in SC_CM_u.
function SC_CM_u_Callback(hObject, eventdata, handles)
% hObject    handle to SC_CM_u (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns SC_CM_u contents as cell array
%        contents{get(hObject,'Value')} returns selected item from SC_CM_u
handles.present.Grade_depend.Cool_Mode_Value=get(hObject,'Value');
contents = cellstr(get(hObject,'String'));
handles.present.Grade_depend.Cool_Mode=contents{get(hObject,'Value')};
a=handles.present.Grade_depend.Cool_Mode_Value;
if a==5
    set(handles.SC_Z1_u,'Style','edit','BackgroundColor',[1 1 1]);
    set(handles.SC_Z2_u,'Style','edit','BackgroundColor',[1 1 1]);
    set(handles.SC_Z3_u,'Style','edit','BackgroundColor',[1 1 1]);
    WFR_Lm3(1)=str2double(get(handles.SC_Z1_u,'String'));
    WFR_Lm3(2)=str2double(get(handles.SC_Z2_u,'String'));
    WFR_Lm3(3)=str2double(get(handles.SC_Z3_u,'String'));
    
    set(handles.SC_Z1L_u,'Style','edit','BackgroundColor',[1 1 1]);
    set(handles.SC_Z2L_u,'Style','edit','BackgroundColor',[1 1 1]);
    set(handles.SC_Z3L_u,'Style','edit','BackgroundColor',[1 1 1]);
    
else
    set(handles.SC_Z1_u,'Style','text','BackgroundColor',0.941*[1 1 1]);
    set(handles.SC_Z2_u,'Style','text','BackgroundColor',0.941*[1 1 1]);
    set(handles.SC_Z3_u,'Style','text','BackgroundColor',0.941*[1 1 1]);
    
    set(handles.SC_Z1L_u,'Style','text','BackgroundColor',0.941*[1 1 1]);
    set(handles.SC_Z2L_u,'Style','text','BackgroundColor',0.941*[1 1 1]);
    set(handles.SC_Z3L_u,'Style','text','BackgroundColor',0.941*[1 1 1]);
    
    WFR_Lm3=handles.default.Std_cooling(a,:);
end
    handles.present.Grade_depend.WFR_Lm3=WFR_Lm3;
    set(handles.SC_Z1_u,'String',WFR_Lm3(1));
    set(handles.SC_Z2_u,'String',WFR_Lm3(2));
    set(handles.SC_Z3_u,'String',WFR_Lm3(3));    
    CS=handles.present.Grade_depend.Casting_speed;
    set(handles.SC_Z1L_u,'String',WFR_Lm3(1)*CS);
    set(handles.SC_Z2L_u,'String',WFR_Lm3(2)*CS);
    set(handles.SC_Z3L_u,'String',WFR_Lm3(3)*CS);
[cooling_info]=Update_Cooling_Info(handles); 

handles.present.cooling_info=cooling_info;
guidata(hObject,handles);
addpath(handles.Folders.Main_Fold.Set);
[SEC_COOL,Nozzle_Tab,No_Noz_Zo,cooling_info]=Updade_Cooling_Info(handles,3);
handles.present.SEC_COOL=SEC_COOL;
handles.present.No_Noz_Zo=No_Noz_Zo;
handles.present.Nozzle_Tab=Nozzle_Tab;
handles.present.cooling_info=cooling_info;

handles.user.SEC_COOL=SEC_COOL;
handles.user.No_Noz_Zo=No_Noz_Zo;
handles.user.Nozzle_Tab=Nozzle_Tab;
handles.user.cooling_info=cooling_info;

guidata(hObject,handles);

cla(handles.Secondary_plot);
axes(handles.Secondary_plot);
plot_Secondary_Cooling(handles)

guidata(hObject,handles);

    function[cooling_info]=Update_Cooling_Info(handles)
        WFR_Lm3=handles.present.Grade_depend.WFR_Lm3;
        Casting_speed=handles.present.Grade_depend.Casting_speed;
        Water_share_matrix=handles.present.Water_share_matrix;
        WFR_Lmin6=Water_share_matrix*WFR_Lm3'*Casting_speed
        cooling_info=handles.present.cooling_info;
%         cooling_info.WFR={WFR_Lmin6(1),WFR_Lmin6(2),WFR_Lmin6(3),WFR_Lmin6(4),WFR_Lmin6(5),WFR_Lmin6(6)};
for i=1:length(WFR_Lmin6)
    cooling_info(i).WFR=WFR_Lmin6(i);
end
        % --- Executes during object creation, after setting all properties.
function SC_CM_u_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SC_CM_u (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





% --- Executes during object creation, after setting all properties.
function SC_Z2_u_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SC_Z2_u (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes during object creation, after setting all properties.
function SC_Z3_u_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SC_Z3_u (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




function SC_Z1_u_Callback(hObject, eventdata, handles)
% hObject    handle to SC_Z1_u (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SC_Z1_u as text
%        str2double(get(hObject,'String')) returns contents of SC_Z1_u as a double
Z1=str2double(get(hObject,'String'));
handles.present.Grade_depend.WFR_Lm3(1)=Z1;
set(handles.SC_Z1L_u,'String',Z1*handles.present.Grade_depend.Casting_speed);
guidata(hObject,handles);
[SEC_COOL,Nozzle_Tab,No_Noz_Zo,cooling_info]=Updade_Cooling_Info(handles,3);
handles.present.SEC_COOL=SEC_COOL;
handles.present.No_Noz_Zo=No_Noz_Zo;
handles.present.Nozzle_Tab=Nozzle_Tab;
handles.present.cooling_info=cooling_info;

cla(handles.Secondary_plot);
axes(handles.Secondary_plot);
plot_Secondary_Cooling(handles);



function SC_Z2_u_Callback(hObject, eventdata, handles)
% hObject    handle to SC_Z2_u (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SC_Z2_u as text
%        str2double(get(hObject,'String')) returns contents of SC_Z2_u as a double
Z2=str2double(get(hObject,'String'));
handles.present.Grade_depend.WFR_Lm3(2)=Z2;
set(handles.SC_Z2L_u,'String',Z2*handles.present.Grade_depend.Casting_speed);
guidata(hObject,handles);

[SEC_COOL,Nozzle_Tab,No_Noz_Zo,cooling_info]=Updade_Cooling_Info(handles,3);
handles.present.SEC_COOL=SEC_COOL;
handles.present.No_Noz_Zo=No_Noz_Zo;
handles.present.Nozzle_Tab=Nozzle_Tab;
handles.present.cooling_info=cooling_info;
guidata(hObject,handles);

cla(handles.Secondary_plot);
axes(handles.Secondary_plot);
plot_Secondary_Cooling(handles);


function SC_Z3_u_Callback(hObject, eventdata, handles)
% hObject    handle to SC_Z3_u (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SC_Z3_u as text
%        str2double(get(hObject,'String')) returns contents of SC_Z3_u as a double
Z3=str2double(get(hObject,'String'));
handles.present.Grade_depend.WFR_Lm3(3)=Z3;
set(handles.SC_Z3L_u,'String',Z3*handles.present.Grade_depend.Casting_speed);
guidata(hObject,handles);
[SEC_COOL,Nozzle_Tab,No_Noz_Zo,cooling_info]=Updade_Cooling_Info(handles,3);
handles.present.SEC_COOL=SEC_COOL;
handles.present.No_Noz_Zo=No_Noz_Zo;
handles.present.Nozzle_Tab=Nozzle_Tab;
handles.present.cooling_info=cooling_info;
guidata(hObject,handles);

cla(handles.Secondary_plot);
axes(handles.Secondary_plot);
plot_Secondary_Cooling(handles);


% --- Executes during object creation, after setting all properties.
function SC_Z1_u_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SC_Z1_u (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Secondary_Reset.
function Secondary_Reset_Callback(hObject, eventdata, handles)
% hObject    handle to Secondary_Reset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
CS=handles.D_CSi;
D_WFR=handles.D_WFR;
handles.M_WFR=D_WFR;

set(handles.SC_Z1_u,'String',D_WFR(1),'Style','text','BackgroundColor',0.941*[1 1 1]);
set(handles.SC_Z2_u,'String',D_WFR(2),'Style','text','BackgroundColor',0.941*[1 1 1]);
set(handles.SC_Z3_u,'String',D_WFR(3),'Style','text','BackgroundColor',0.941*[1 1 1]);


% handles.SC_Z1_u=D_WFR(1);
% handles.SC_Z2_u=D_WFR(2);
% handles.SC_Z3_u=D_WFR(3);

set(handles.SC_Z1L_u,'String',D_WFR(1)*CS);
set(handles.SC_Z2L_u,'String',D_WFR(2)*CS);
set(handles.SC_Z3L_u,'String',D_WFR(3)*CS);

set(handles.SC_CM_u,'Value',handles.Secondary_D_CMV);
handles.M_cooling_info=handles.D_cooling_info;

guidata(hObject,handles);
[SEC_COOL,Nozzle_Tab,No_Noz_Zo,cooling_info]=Updade_Cooling_Info(handles);
handles.present.SEC_COOL=SEC_COOL;
handles.present.No_Noz_Zo=No_Noz_Zo;
handles.present.Nozzle_Tab=Nozzle_Tab;
handles.present.cooling_info=cooling_info;

cla(handles.Secondary_plot);
axes(handles.Secondary_plot);
plot_Secondary_Cooling(handles)

guidata(hObject,handles);

% --- Executes on button press in Secondary_Default.
function Secondary_Default_Callback(hObject, eventdata, handles)
% hObject    handle to Secondary_Default (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
CS=handles.D_CSi;
D_WFR=handles.D_WFR;
handles.M_WFR=D_WFR;

set(handles.SC_Z1_u,'String',D_WFR(1),'Style','text','BackgroundColor',0.941*[1 1 1]);
set(handles.SC_Z2_u,'String',D_WFR(2),'Style','text','BackgroundColor',0.941*[1 1 1]);
set(handles.SC_Z3_u,'String',D_WFR(3),'Style','text','BackgroundColor',0.941*[1 1 1]);


% handles.SC_Z1_u=D_WFR(1);
% handles.SC_Z2_u=D_WFR(2);
% handles.SC_Z3_u=D_WFR(3);

set(handles.SC_Z1L_u,'String',D_WFR(1)*CS);
set(handles.SC_Z2L_u,'String',D_WFR(2)*CS);
set(handles.SC_Z3L_u,'String',D_WFR(3)*CS);

% Secondary_M_CMV=get(handles.SC_CM_u,'Value');
% cool_mode=zeros(5,1);cool_mode(handles.Secondary_D_CMV)=1;
% [handles.Secondary_D_CM,Secondary_D_CMV]=cool_mode_sel(cool_mode);
% handles.Secondary_M_CMV=handles.Secondary_D_CMV;
% handles.Secondary_M_CM=handles.Secondary_D_CM;
Cool_Mode=handles.Secondary_D_CM;

set(handles.SC_CM_u,'Value',handles.Secondary_D_CMV);
handles.M_cooling_info=handles.D_cooling_info;

guidata(hObject,handles);
[SEC_COOL,Nozzle_Tab,No_Noz_Zo,cooling_info]=Updade_Cooling_Info(handles);
handles.present.SEC_COOL=SEC_COOL;
handles.present.No_Noz_Zo=No_Noz_Zo;
handles.present.Nozzle_Tab=Nozzle_Tab;
handles.present.cooling_info=cooling_info;
cla(handles.Secondary_plot);
axes(handles.Secondary_plot);
plot_Secondary_Cooling(handles);

Folder_Name=handles.FNi;
SBC=[Folder_Name '\Secondary.mat'];
save(SBC,'Z','X','Q','HTC','Nozzle_Tab','cooling_info','Cool_Mode');


guidata(hObject,handles);

% --- Executes on button press in Secondary_Modified.
function Secondary_Modified_Callback(hObject, eventdata, handles)
% hObject    handle to Secondary_Modified (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
CS=handles.M_CSi;
M_WFR=handles.M_WFR;


set(handles.SC_Z1_u,'String',M_WFR(1),'Style','text','BackgroundColor',0.941*[1 1 1]);
set(handles.SC_Z2_u,'String',M_WFR(2),'Style','text','BackgroundColor',0.941*[1 1 1]);
set(handles.SC_Z3_u,'String',M_WFR(3),'Style','text','BackgroundColor',0.941*[1 1 1]);


% handles.SC_Z1_u=D_WFR(1);
% handles.SC_Z2_u=D_WFR(2);
% handles.SC_Z3_u=D_WFR(3);

set(handles.SC_Z1L_u,'String',M_WFR(1)*CS);
set(handles.SC_Z2L_u,'String',M_WFR(2)*CS);
set(handles.SC_Z3L_u,'String',M_WFR(3)*CS);

Secondary_M_CMV=get(handles.SC_CM_u,'Value');
cool_mode=zeros(5,1);cool_mode(Secondary_M_CMV)=1;
[Secondary_M_CM,Secondary_M_CMV]=cool_mode_sel(cool_mode);
Cool_Mode=Secondary_M_CM;
% handles.M_cooling_info=handles.D_cooling_info;

guidata(hObject,handles);
[SEC_COOL,Nozzle_Tab,No_Noz_Zo,cooling_info]=Updade_Cooling_Info(handles);
handles.present.SEC_COOL=SEC_COOL;
handles.present.No_Noz_Zo=No_Noz_Zo;
handles.present.Nozzle_Tab=Nozzle_Tab;
handles.present.cooling_info=cooling_info;
cla(handles.Secondary_plot);
axes(handles.Secondary_plot);
plot_Secondary_Cooling(handles);

Folder_Name=handles.FNi;
SBC=[Folder_Name '\Secondary.mat'];
save(SBC,'Z','X','Q','HTC','V','D','Nozzle_Tab','cooling_info','Cool_Mode');
guidata(hObject,handles);



function[]=plot_Secondary_Cooling(handles)
m=handles.present;
Z=m.SEC_COOL.Z;
X=m.SEC_COOL.X;
Q=m.SEC_COOL.Q;
HTC=m.SEC_COOL.HTC;
Nozzle_Tab=m.Nozzle_Tab;
Mol_Len=m.Mol_Len;
Wi=m.Wi;
sel=get(handles.Secondary_Flow_or_HTC,'Value');
sec_L=Nozzle_Tab(size(Nozzle_Tab,1),1)+700;
if sel==2
    [C,h] = contour(Z,X,HTC);
    title('HTC in W/m^2K','FontSize',12);
else
    [C,h] = contour(Z,X,Q);
    title('Flow rate in Litre/min/m^2','FontSize',12);
end
cla(handles.Secondary_plot);
axes(handles.Secondary_plot);
rectangle('Position',[Mol_Len,0,sec_L-Mol_Len,Wi],'LineWidth',2,'Edgecolor','r');hold on
xlabel('Axial direction [mm]','FontSize',12);
ylabel('Lateral direction [mm]','FontSize',12);
set(gca,'LineWidth',2,'FontSize',10);
% axis equal
axis([Mol_Len-200 sec_L -30 Wi+30]);
colorbar 


% --- Executes on button press in Second_ApplyDefault.
function Second_ApplyDefault_Callback(hObject, eventdata, handles)
% hObject    handle to Second_ApplyDefault (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.present.Grade_depend.WFR_Lm3=handles.Grade_depend.WFR_Lm3;
handles.present.Grade_depend.Cool_Mode_Value=handles.Grade_depend.Cool_Mode_Value;
handles.present.Grade_depend.Cool_Mode=handles.Grade_depend.Cool_Mode;
set(handles.SC_CM_u,'Value',handles.Grade_depend.Cool_Mode_Value);
a=handles.Grade_depend.Cool_Mode_Value;
if a==5
    set(handles.SC_Z1_u,'Style','edit','BackgroundColor',[1 1 1]);
    set(handles.SC_Z2_u,'Style','edit','BackgroundColor',[1 1 1]);
    set(handles.SC_Z3_u,'Style','edit','BackgroundColor',[1 1 1]);
else
    set(handles.SC_Z1_u,'Style','text','BackgroundColor',0.941*[1 1 1]);
    set(handles.SC_Z2_u,'Style','text','BackgroundColor',0.941*[1 1 1]);
    set(handles.SC_Z3_u,'Style','text','BackgroundColor',0.941*[1 1 1]);
 end
  
guidata(hObject,handles);

set(handles.SC_Z1_u,'String',handles.Grade_depend.WFR_Lm3(1));
set(handles.SC_Z2_u,'String',handles.Grade_depend.WFR_Lm3(2));
set(handles.SC_Z3_u,'String',handles.Grade_depend.WFR_Lm3(3));
set(handles.SC_Z1L_u,'String',handles.Grade_depend.WFR_Lm3(1)*handles.present.Grade_depend.Casting_speed);
set(handles.SC_Z2L_u,'String',handles.Grade_depend.WFR_Lm3(2)*handles.present.Grade_depend.Casting_speed);
set(handles.SC_Z3L_u,'String',handles.Grade_depend.WFR_Lm3(3)*handles.present.Grade_depend.Casting_speed);
handles.present.cooling_info=handles.default.cooling_info;
guidata(hObject,handles);
[SEC_COOL,Nozzle_Tab,No_Noz_Zo,cooling_info]=Updade_Cooling_Info(handles);
handles.present.SEC_COOL=SEC_COOL;
handles.present.No_Noz_Zo=No_Noz_Zo;
handles.present.Nozzle_Tab=Nozzle_Tab;
handles.present.cooling_info=cooling_info;
guidata(hObject,handles);

cla(handles.Secondary_plot);
axes(handles.Secondary_plot);
plot_Secondary_Cooling(handles);

% --- Executes on button press in Second_UserSet.
function Second_UserSet_Callback(hObject, eventdata, handles)
% hObject    handle to Second_UserSet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.present.Grade_depend.WFR_Lm3=handles.user.Grade_depend.WFR_Lm3;
handles.present.Grade_depend.Cool_Mode_Value=handles.user.Grade_depend.Cool_Mode_Value;
handles.present.Grade_depend.Cool_Mode=handles.user.Grade_depend.Cool_Mode;
set(handles.SC_CM_u,'Value',handles.user.Grade_depend.Cool_Mode_Value);
[handles.present.cooling_info]=Update_Cooling_Info(handles);
a=handles.user.Grade_depend.Cool_Mode_Value;
if a==5
    set(handles.SC_Z1_u,'Style','edit','BackgroundColor',[1 1 1]);
    set(handles.SC_Z2_u,'Style','edit','BackgroundColor',[1 1 1]);
    set(handles.SC_Z3_u,'Style','edit','BackgroundColor',[1 1 1]);
else
    set(handles.SC_Z1_u,'Style','text','BackgroundColor',0.941*[1 1 1]);
    set(handles.SC_Z2_u,'Style','text','BackgroundColor',0.941*[1 1 1]);
    set(handles.SC_Z3_u,'Style','text','BackgroundColor',0.941*[1 1 1]);
end
  
guidata(hObject,handles);

set(handles.SC_Z1_u,'String',handles.user.Grade_depend.WFR_Lm3(1));
set(handles.SC_Z2_u,'String',handles.user.Grade_depend.WFR_Lm3(2));
set(handles.SC_Z3_u,'String',handles.user.Grade_depend.WFR_Lm3(3));
set(handles.SC_Z1L_u,'String',handles.user.Grade_depend.WFR_Lm3(1)*handles.present.Grade_depend.Casting_speed);
set(handles.SC_Z2L_u,'String',handles.user.Grade_depend.WFR_Lm3(2)*handles.present.Grade_depend.Casting_speed);
set(handles.SC_Z3L_u,'String',handles.user.Grade_depend.WFR_Lm3(3)*handles.present.Grade_depend.Casting_speed);
guidata(hObject,handles);
[SEC_COOL,Nozzle_Tab,No_Noz_Zo,cooling_info]=Updade_Cooling_Info(handles);
handles.present.SEC_COOL=SEC_COOL;
handles.present.No_Noz_Zo=No_Noz_Zo;
handles.present.Nozzle_Tab=Nozzle_Tab;
handles.present.cooling_info=cooling_info;
guidata(hObject,handles);

cla(handles.Secondary_plot);
axes(handles.Secondary_plot);
plot_Secondary_Cooling(handles);


% --- Executes on button press in CCM_ApplyDefault.
function CCM_ApplyDefault_Callback(hObject, eventdata, handles)
% hObject    handle to CCM_ApplyDefault (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Nozzle_Tab=handles.default.Nozzle_Tab;
POS_NOZ=Nozzle_Tab(:,[1 2]);
WAT_FLO=Nozzle_Tab(:,[4]);
NOZ_DIS=Nozzle_Tab(:,[5]);
CON_ANG=Nozzle_Tab(:,[7]);
No_Noz_Zo=handles.default.No_Noz_Zo;
Spray_rad=Nozzle_Tab(:,[9]);
[D_DATA_struct,D_DATA_table]=make_CCM_Data(POS_NOZ,WAT_FLO,NOZ_DIS,CON_ANG,No_Noz_Zo,Spray_rad);
set(handles.CCM_Modified,'Data',D_DATA_table);
handles.present.cooling_info=handles.default.cooling_info;
guidata(hObject,handles);
[SEC_COOL,Nozzle_Tab,No_Noz_Zo,cooling_info]=Updade_Cooling_Info(handles,3);
handles.present.SEC_COOL=SEC_COOL;
handles.present.No_Noz_Zo=No_Noz_Zo;
handles.present.Nozzle_Tab=Nozzle_Tab;
handles.present.cooling_info=cooling_info;
guidata(hObject,handles);

% --- Executes on button press in CCM_UserSet.
function CCM_UserSet_Callback(hObject, eventdata, handles)
% hObject    handle to CCM_UserSet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Nozzle_Tab=handles.user.Nozzle_Tab;
POS_NOZ=Nozzle_Tab(:,[1 2]);
WAT_FLO=Nozzle_Tab(:,[4]);
NOZ_DIS=Nozzle_Tab(:,[5]);
CON_ANG=Nozzle_Tab(:,[7]);
No_Noz_Zo=handles.user.No_Noz_Zo;
Spray_rad=Nozzle_Tab(:,[9]);
[D_DATA_struct,D_DATA_table]=make_CCM_Data(POS_NOZ,WAT_FLO,NOZ_DIS,CON_ANG,No_Noz_Zo,Spray_rad);
D_DATA_table
set(handles.CCM_Modified,'Data',D_DATA_table);
handles.present.cooling_info=handles.user.cooling_info;
guidata(hObject,handles);
[SEC_COOL,Nozzle_Tab,No_Noz_Zo,cooling_info]=Updade_Cooling_Info(handles,3);
handles.present.SEC_COOL=SEC_COOL;
handles.present.No_Noz_Zo=No_Noz_Zo;
handles.present.Nozzle_Tab=Nozzle_Tab;
handles.present.cooling_info=cooling_info;

guidata(hObject,handles);


function[DATA,DATAP]=make_CCM_Data(POS_NOZ,WAT_FLO,NOZ_DIS,CON_ANG,No_Noz_Zo,Spray_rad)
ZONE=[];
NRO=[];
ADIST=[];
PRDIA=[];
COAN=[];
DIST=[];
WFR=[];
a1=0;a3=0;
for i=1:length(No_Noz_Zo)
    a2=No_Noz_Zo(i);
    for j=1:a2-a1
        rind=a1+j;
        ind=a1+j+a3;
        if  rem(i,2)==1
            ZONE{ind}=[num2str(round(i/2)) 'a'];
        else
            ZONE{ind}=[num2str(round(i/2)) 'b'];
        end
        NRO{ind}=1;
        ADIST{ind}=round(POS_NOZ(rind)*10)/10;
%         PRDIA{ind}=round(2*Spray_rad(rind)*10)/10;
        PRDIA{ind}=round(NOZ_DIS(rind)*tan(CON_ANG(rind)*pi/360)*10)/10;
        COAN{ind}=round(CON_ANG(rind));
        DIST{ind}=NOZ_DIS(rind);
        WFR{ind}=WAT_FLO(rind);
%         if POS_NOZ
    end
    %%%%%%% 2 nozzles per row
    for j=1:a2-a1-1
        ind1=a1+j;ind2=a1+j+1;
        if abs(POS_NOZ(ind1)-POS_NOZ(ind2))<1e-3
            ZONE{ind2}=[];
            NRO{ind1}=2;NRO{ind2}=[];
            ADIST{ind2}=[];PRDIA{ind2}=[];COAN{ind2}=[];
            DIST{ind2}=[];
            WFR{ind1}=WFR{ind1}+WFR{ind2};WFR{ind2}=[];
            a3=-1;
        end
    end
    a1=a2;
end
DATA=struct('ZONE',ZONE,...
             'NRO',NRO,...
             'ADIST', ADIST,...
             'PRDIA',PRDIA,... 
             'COAN', COAN,...
             'DIST',DIST,...
             'WFR',WFR);
DATAP=[{DATA.ZONE}' {DATA.NRO}'  {DATA.ADIST}' {DATA.PRDIA}' {DATA.COAN}' DIST' WFR'];

% --- Executes when entered data in editable cell(s) in CCM_Modified.
function CCM_Modified_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to CCM_Modified (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)
handles.M_DATA_table=get(hObject,'Data');
M_DATA_table=handles.M_DATA_table;
% a={M_DATA_table{:,7}}
% data = get(tablehandle,'Data')


handles.Changed_Data(1,:)=eventdata.Indices;
handles.Changed_Data(2,:)=[eventdata.PreviousData eventdata.NewData];
Changed_Data=handles.Changed_Data;
if Changed_Data(1,2)==7
Redistribute_Nozzle_Flow(Changed_Data,M_DATA_table)
end
guidata(hObject,handles);

function[]=Redistribute_Nozzle_Flow(Changed_Data,M_DATA_table)


% --- Executes on button press in Set_ApplyDefault.
function Set_ApplyDefault_Callback(hObject, eventdata, handles)
% hObject    handle to Set_ApplyDefault (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.present=handles.default;
guidata(hObject, handles);
set(handles.SWi_u,'String',handles.present.Wi);
set(handles.STh_u,'String',handles.present.Th);
set(handles.Sner_u,'String',handles.present.ner);
set(handles.Snes_u,'String',handles.present.nes);
set(handles.Snes_u,'String',handles.present.nes);
% set(handles.Sele_arrange_u,'String',ele_arrange);
set(handles.pop_ele_arrange,'Value',handles.present.ele_arrange_V);
handles.ele_arrange_V_u=handles.present.ele_arrange_V;
% set(handles.Sele_order_u,'String',ele_order);
set(handles.pop_ele_order,'Value',handles.present.ele_order_V);
set(handles.SSim_Len_u,'String',handles.present.Sim_Len);
set(handles.Sdt_u,'String',handles.present.dt);
set(handles.STOL_u,'String',handles.present.TOL);
% cla(handles.mesh);
% axes(handles.mesh);

draw_heat_4n(handles);


function[XYH,MAPH,BC_Info]=Create_mesh(handles)
    m=handles.present;
XYe = [m.Wi, m.Th; 
       0, m.Th;
       0, 0;
       m.Wi, 0];
  BC_faces= [4 1 3 2];  
symmetry=get(handles.SSymmetry,'Value')-1;  
  if symmetry==3 %%% two way symmetry
      XYe=XYe./2;
      BC_faces=[1 3];
  elseif symmetry==1%%% x symmetry
      XYe = [m.Wi, m.Th/2; 0, m.Th/2;0, 0;m.Wi, 0];
      BC_faces= [1 3 2];
  elseif symmetry==2%%% y symmetry
      XYe = [m.Wi/2, m.Th; 0, m.Th;0, 0;m.Wi/2, 0];
      BC_faces= [4 1 3];
  end
[XYH,MAPH]=heat_mesh_gen(XYe,m.ner,m.nes,m.ele_arrange_V); %%% meshing
[BCEL,BCFC] = make_BCEL_BCFC(m.ner,m.nes,BC_faces,1,0);
BC_Info= [BCEL' BCFC'];


function[BCEL,BCFC]=make_BCEL_BCFC(r,s,face,dir,st_el_num)
%%%%%%% INPUT
%%% r    - number of elements x direction
%%% s    - number of elements y direction
%%% face - boundary faces
%%% dir  - either clockwise[0] or anticlockwise[1]
%%% st_el_num - starting element number, fist mesh-0
%%%%%%% OUTPUT
%%% BCEL - Boundary elements
%%% BCFC - Boundary faces
%%%% Assumptions
%%% only for first rectangular mesh
BCEL=[];BCFC=[];
for i=1:length(face)
    f=face(i);
    if f==1
        el = ((r-1)*s)+1:1:r*s;
        fa = zeros(size(el))+1;
    elseif f==2
        el = 1:1:s;
        fa = zeros(size(el))+2;
    elseif f==3
        el = s:s:r*s;
        fa = zeros(size(el))+3;
    else
        el = 1:s:((r-1)*s)+1;
        fa = zeros(size(el))+4;
    end
    if dir==1
        if f==3 || f==2
            el=reverse(el);
        end
    else
        if f==1 || f==4
            el=reverse(el);
        end
    end
    el=el+st_el_num;
    BCEL=[BCEL el];BCFC=[BCFC fa];
end
function[b]=reverse(a)
n=length(a);    
for i=1:n
    b(i)=a(1+n-i);
end

function [XYH,MAPH] = heat_mesh_gen(XYe,ner,nes,m)
%  Refine heat mesh
%  Sub-function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  INPUT:                                                                    %
%    [XYH]  -> Nodal coordinates                                             %
%    [MAPH] -> Element connectivity                                          %
%  OUTPUT:    Draws                                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if m==2
    a=0;b=1;x=10.^(linspace(a,b,ner+1));
    y=-2*((x-10^a)/(10^b-10^a))+1;
    for i=1:length(y)
        Rr(i)=y(length(y)+1-i);
    end
    clear y a b x
    a=0;b=1;x=10.^(linspace(a,b,nes+1));
    y=-2*((x-10^a)/(10^b-10^a))+1;
    for i=1:length(y)
        Ss(i)=y(length(y)+1-i);
    end
else
    Rr = linspace(-1,1,ner+1);
    Ss = linspace(-1,1,nes+1);
end



[RS,MAPH] = heat_mesh_grid(Rr,Ss);

XYH = RS; % initialize
nn  = size(XYH); % number of nodes

for in=1:nn % map nodes to global coordinates
  r = RS(in,1);
  s = RS(in,2);
  N = 0.25*[(1+r)*(1+s); 
            (1-r)*(1+s);
            (1-r)*(1-s);
            (1+r)*(1-s)]; 

  XYH(in,:) = N'*XYe;
end % for(in=1:nn)

function [XY,MAP] = heat_mesh_grid(X,Y)
%  Mesh Generation for the system
%  Sub-function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  This function is called only once during inputdata to initialize some     %
%  phase transformation constants                                            %
%  INPUT:                                                                    %
%    [X] -> Array of x-coordinates.                                          %
%    [Y] -> Array of y-coordinates.                                          %
%  OUTPUT:                                                                   %
%    [XY]  -> Nodal coordinates                                              %
%    [MAP] -> Element connectivity                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nx = size(X,2);
ny = size(Y,2);
nn = nx*ny;
ne = (nx-1)*(ny-1);
XY = zeros(nn,2);
MAP = zeros(ne,4);

% Node Numbering And Node Coordinates
nc = 0;
for ix=1:nx
  for iy=1:ny
    nc = nc + 1;
    XY(nc,1) = X(ix);
    XY(nc,2) = Y(iy);
  end
end

% Element Mapping
ec = 0;
xx = ny;
yy = 1;
for ix=0:(nx-2)
  for iy=0:(ny-2)
    ec = ec + 1;
    nc = xx*ix+yy*iy+1;
    MAP(ec,1) = nc+xx+yy;
    MAP(ec,2) = nc+yy;
    MAP(ec,3) = nc;
    MAP(ec,4) = nc+xx;
  end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function draw_heat_4n(handles)

    cla(handles.mesh);
    axes(handles.mesh);
    m=handles.present;
    XY =m.XYH;
    MAP =m.MAPH;
    e=handles.el_no;
    n=handles.no_no;
ne = size(MAP,1); % Number of elements
map = [1 2 3 4 1];
nodes=[1 2 3 4];
% q=0;
% qq=0;
for ie=1:ne
    % nodal points in draw order
    X  = XY(MAP(ie,map),1);
    Y  = XY(MAP(ie,map),2);
    line(X,Y,'Marker','o');
    if e>0
        %       disp('poi')
        %       set(qq, 'visible', 'off')
        xm=0.25*(X(1)+X(2)+X(3)+X(4));
        ym=0.25*(Y(1)+Y(2)+Y(3)+Y(4));
        q=text(xm,ym,int2str(ie));
    end
end
if n>0
    tol=max(X)*5e-3;
    for ie=1:ne
        x_node =XY(MAP(ie,nodes),1);
        y_node =XY(MAP(ie,nodes),2);
        for j=1:length(nodes)
            node_number=MAP(ie,j);
            x_coord=x_node(j);y_coord=y_node(j);
            qq=text(x_coord+tol,y_coord,int2str(node_number));
        end
    end
end
xmax=max(XY(:,1));
xmin=min(XY(:,1));
ymax=max(XY(:,2));
ymin=min(XY(:,2));

% xlabel('Width direction [mm]','FontSize',11);
ylabel('Thickness direction [mm]','FontSize',11);
set(gca,'LineWidth',2,'FontSize',11,'Box','on');
axis equal
axis tight
axis([xmin-5 xmax+5  ymin-5 ymax+5]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% --- Executes on button press in Set_UserSet.
function Set_UserSet_Callback(hObject, eventdata, handles)
% hObject    handle to Set_UserSet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.present=handles.user;
guidata(hObject, handles);
set(handles.SWi_u,'String',handles.present.Wi);
set(handles.STh_u,'String',handles.present.Th);
set(handles.Sner_u,'String',handles.present.ner);
set(handles.Snes_u,'String',handles.present.nes);
set(handles.Snes_u,'String',handles.present.nes);
% set(handles.Sele_arrange_u,'String',ele_arrange);
set(handles.pop_ele_arrange,'Value',handles.present.ele_arrange_V);
handles.ele_arrange_V_u=handles.present.ele_arrange_V;
% set(handles.Sele_order_u,'String',ele_order);
set(handles.pop_ele_order,'Value',handles.present.ele_order_V);
set(handles.SSim_Len_u,'String',handles.present.Sim_Len);
set(handles.Sdt_u,'String',handles.present.dt);
set(handles.STOL_u,'String',handles.present.TOL);
% cla(handles.mesh);
% axes(handles.mesh);
% [SEC_COOL,Nozzle_Tab,No_Noz_Zo,cooling_info]=Updade_Cooling_Info(handles);
% handles.present.SEC_COOL=SEC_COOL;
% handles.present.No_Noz_Zo=No_Noz_Zo;
% handles.present.Nozzle_Tab=Nozzle_Tab;
% handles.present.cooling_info=cooling_info;

guidata(hObject, handles);

draw_heat_4n(handles);

function Sdt_u_Callback(hObject, eventdata, handles)
% hObject    handle to Sdt_u (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Sdt_u as text
%        str2double(get(hObject,'String')) returns contents of Sdt_u as a double
handles.present.dt=str2double(get(hObject,'String'));
guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function Sdt_u_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Sdt_u (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function SSim_Len_u_Callback(hObject, eventdata, handles)
% hObject    handle to SSim_Len_u (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SSim_Len_u as text
%        str2double(get(hObject,'String')) returns contents of SSim_Len_u as a double
handles.present.Sim_Len=str2double(get(hObject,'String'));
guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function SSim_Len_u_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SSim_Len_u (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function STOL_u_Callback(hObject, eventdata, handles)
% hObject    handle to STOL_u (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of STOL_u as text
%        str2double(get(hObject,'String')) returns contents of STOL_u as a double
handles.present.TOL=str2double(get(hObject,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function STOL_u_CreateFcn(hObject, eventdata, handles)
% hObject    handle to STOL_u (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function SWi_u_Callback(hObject, eventdata, handles)
% hObject    handle to SWi_u (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SWi_u as text
%        str2double(get(hObject,'String')) returns contents of SWi_u as a double
handles.present.Wi=str2double(get(hObject,'String'));
guidata(hObject, handles);
[handles.present.XYH,handles.present.MAPH,...
    handles.present.BC_Info]=Create_mesh(handles);
guidata(hObject, handles);
% cla(handles.mesh);
% axes(handles.mesh);
draw_heat_4n(handles);


% --- Executes during object creation, after setting all properties.
function SWi_u_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SWi_u (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function STh_u_Callback(hObject, eventdata, handles)
% hObject    handle to STh_u (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of STh_u as text
%        str2double(get(hObject,'String')) returns contents of STh_u as a double
handles.present.Th=str2double(get(hObject,'String'));
guidata(hObject, handles);
[handles.present.XYH,handles.present.MAPH,...
    handles.present.BC_Info]=Create_mesh(handles);
guidata(hObject, handles);
% cla(handles.mesh);
% axes(handles.mesh);
draw_heat_4n(handles);

% --- Executes during object creation, after setting all properties.
function STh_u_CreateFcn(hObject, eventdata, handles)
% hObject    handle to STh_u (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in pop_ele_arrange.
function pop_ele_arrange_Callback(hObject, eventdata, handles)
% hObject    handle to pop_ele_arrange (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pop_ele_arrange contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pop_ele_arrange
handles.present.ele_arrange_V=get(hObject,'Value');
List=get(hObject,'String');
handles.present.ele_arrange=List{get(hObject,'Value')};
guidata(hObject, handles);
[handles.present.XYH,handles.present.MAPH,...
    handles.present.BC_Info]=Create_mesh(handles);
guidata(hObject, handles);
% cla(handles.mesh);
% axes(handles.mesh);
draw_heat_4n(handles);

% --- Executes during object creation, after setting all properties.
function pop_ele_arrange_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pop_ele_arrange (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in pop_ele_order.
function pop_ele_order_Callback(hObject, eventdata, handles)
% hObject    handle to pop_ele_order (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pop_ele_order contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pop_ele_order
handles.present.ele_order_V=get(hObject,'Value');
List=get(hObject,'String');
handles.present.ele_order=List{get(hObject,'Value')};
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function pop_ele_order_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pop_ele_order (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Sner_u_Callback(hObject, eventdata, handles)
% hObject    handle to Sner_u (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Sner_u as text
%        str2double(get(hObject,'String')) returns contents of Sner_u as a double
ner=str2double(get(hObject,'String'));
handles.present.ner=ner;
guidata(hObject, handles);
m=handles.present;
[XYH,MAPH,BC_Info]=Create_mesh(handles);
[BC_Info2]=find_reduction_factor(XYH,MAPH,BC_Info,m.f1,m.f2,m.Wi,m.Th);
handles.present.BC_Info2=BC_Info2;
size(BC_Info2)
handles.present.XYH=XYH;
handles.present.MAPH=MAPH;
handles.present.BC_Info=BC_Info;

guidata(hObject, handles);
% cla(handles.mesh);
% axes(handles.mesh);
draw_heat_4n(handles);

% --- Executes during object creation, after setting all properties.
function Sner_u_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Sner_u (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Snes_u_Callback(hObject, eventdata, handles)
% hObject    handle to Snes_u (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Snes_u as text
%        str2double(get(hObject,'String')) returns contents of Snes_u as a double
nes=str2double(get(hObject,'String'));
handles.present.nes=nes;
guidata(hObject, handles);
m=handles.present;
[XYH,MAPH,BC_Info]=Create_mesh(handles);
[BC_Info2]=find_reduction_factor(XYH,MAPH,BC_Info,m.f1,m.f2,m.Wi,m.Th);
handles.present.BC_Info2=BC_Info2;
size(BC_Info2)
handles.present.XYH=XYH;
handles.present.MAPH=MAPH;
handles.present.BC_Info=BC_Info;

guidata(hObject, handles);
% cla(handles.mesh);
% axes(handles.mesh);
draw_heat_4n(handles);


% --- Executes during object creation, after setting all properties.
function Snes_u_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Snes_u (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





% --- Executes on selection change in Results_Folders.
function Results_Folders_Callback(hObject, eventdata, handles)
% hObject    handle to Results_Folders (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Results_Folders contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Results_Folders
contents = cellstr(get(hObject,'String'));
Res_Fold= contents{get(hObject,'Value')};
Cur_Res_Folder=[handles.Folders.Main_Fold.Sim Res_Fold];
handles.Folders.Sub_Fold.Res_Fold=Res_Fold;
Parameters=load([Cur_Res_Folder handles.path_sep 'Process_Parameters.mat']);
Display_Results_Parameters(Parameters,handles);
% handles.Current_Results_Folder=Current_Results_Folder;
handles.Parameters=Parameters;
cla(handles.temp);axes(handles.temp);
cla(handles.shell);axes(handles.shell);
cla(handles.Axial);axes(handles.Axial);
cla(handles.Lateral);axes(handles.Lateral);
Line_Plots(Cur_Res_Folder,handles)
guidata(hObject, handles);


[AXIAL,DATA]=find_axialdis_DB(handles);
handles.AXIAL=AXIAL;
handles.DB=DATA;
dvalue=size(AXIAL,1);
set(handles.Results_AxialDist,'String',AXIAL(dvalue,1));
set(handles.Results_DataList,'Value',dvalue);
dbindex=DATA{AXIAL(dvalue,2)};handles.dbindex=dbindex;
guidata(hObject, handles);
if get(handles.Con_1,'Value')==1
    Temp_contour(handles);
else
    Color_plot(handles);
end
handles.daxial=AXIAL(dvalue,1);
set(handles.Results_Midplane,'Value',1);
guidata(hObject, handles);

Plot_Axial_Contours(handles)

load([Cur_Res_Folder handles.path_sep 'Material.mat']);
set(handles.Results_TS,'String',round(TS*10)/10);
set(handles.Results_TL,'String',round(TL*10)/10);


% --- Executes during object creation, after setting all properties.
function Results_Folders_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Results_Folders (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Results_Line_Save.
function Results_Line_Save_Callback(hObject, eventdata, handles)
% hObject    handle to Results_Line_Save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Cur_Fol=[handles.Folders.Main_Fold.Sim handles.Folders.Sub_Fold.Res_Fold];
Grade_Nr=handles.Parameters.Grade_depend.Grade_Number;
Grade_Ty=handles.Parameters.Grade_depend.Grade_Type;
File_Sugg=[Cur_Fol handles.path_sep num2str(Grade_Nr) '_' Grade_Ty '_Temp-Shell.fig'];
[FName, PName]=uiputfile(File_Sugg,'Save Figure');
if isequal(FName,0) || isequal(PName,0)
else
    f = figure;
    copyobj([handles.temp,handles.shell],f);
    saveas(f,[PName FName]);
    close(gcf);
end

function[]=Display_Results_Parameters(m,handles)
% display result parameters
set(handles.Results_GN,'String',m.Grade_depend.Grade_Number);
set(handles.Results_GT,'String',m.Grade_depend.Grade_Type);
set(handles.Results_CT,'String',m.Grade_depend.Init_Temp);
set(handles.Results_CS,'String',m.Grade_depend.Casting_speed);
set(handles.Results_CM,'String',m.Grade_depend.Cool_Mode);

WFR_Lm3=m.Grade_depend.WFR_Lm3;
WFR_Lmin6=cell2mat({m.cooling_info.WFR});
WFR_Lm6=WFR_Lm3*handles.WSM';

set(handles.Results_Ia,'String',round(WFR_Lm6(1)*10)/10);
set(handles.Results_Ib,'String',round(WFR_Lm6(2)*10)/10);
set(handles.Results_IIa,'String',round(WFR_Lm6(3)*10)/10);
set(handles.Results_IIb,'String',round(WFR_Lm6(4)*10)/10);
set(handles.Results_IIIa,'String',round(WFR_Lm6(5)*10)/10);
set(handles.Results_IIIb,'String',round(WFR_Lm6(6)*10)/10);

set(handles.Results_ILa,'String',round(WFR_Lmin6(1)*10)/10);
set(handles.Results_ILb,'String',round(WFR_Lmin6(2)*10)/10);
set(handles.Results_IILa,'String',round(WFR_Lmin6(3)*10)/10);
set(handles.Results_IILb,'String',round(WFR_Lmin6(4)*10)/10);
set(handles.Results_IIILa,'String',round(WFR_Lmin6(5)*10)/10);
set(handles.Results_IIILb,'String',round(WFR_Lmin6(6)*10)/10);

set(handles.Results_I,'String',round(WFR_Lm3(1)*10)/10);
set(handles.Results_II,'String',round(WFR_Lm3(2)*10)/10);
set(handles.Results_III,'String',round(WFR_Lm3(3)*10)/10);



function[]=Line_Plots(Folder,handles)
   % [Folder handles.path_sep 'DB*.mat']
d=dir([Folder handles.path_sep 'DB*.mat']);
nDB=size(d,1);
load([Folder handles.path_sep 'Material.mat']);
load([Folder handles.path_sep 'Process_Parameters.mat']);

if symmetry==0
    Th1=0;Wi1=0;
elseif symmetry==3
    Th1=Th/2;Wi1=Wi/2;
elseif  symmetry==2
    Th1=0;Wi1=Wi/2;
elseif  symmetry==1
    Th1=Th/2;Wi1=0;
end

P_center=[Wi/2-Wi1,Th/2-Th1];
P_corner=[Wi-Wi1,Th-Th1];
P_mid_w=[Wi/2-Wi1,Th-Th1];
P_mid_t=[Wi-Wi1,Th/2-Th1];

A=Grade_depend;

[NodeId_cen,d11] = node_at_pt(XYH,P_center);
[NodeId_cor,d22] = node_at_pt(XYH,P_corner);
[NodeId_mw,d33] = node_at_pt(XYH,P_mid_w);
[NodeId_mt,d44] = node_at_pt(XYH,P_mid_t);

Zone_pos=convert_cell_matrix({cooling_info.Zone_pos});
Zones={cooling_info.Zones};

% NodeID=str2double(get(handles.Node,'String'));
% xy=XYH(NodeID,:);
Vc=A.Casting_speed/60;

% [ind]=Contour_node_index(MAPH,[ner nes])
% if rem(ner,2)>0
%     n1=(ner+1)/2;n2=n1+1;
% else
%     n1=round((ner+1)/2); n2=n1;
% end
% if rem(nes,2)>0
%     m1=(nes+1)/2;m2=m1+1;even_add=0;
% else
%     m1=round((nes+1)/2);m2=m1;even_add=1;
% end

% if symmetry==3
%     Nid_ind=ind(:,1);
% % elseif  symmetry==2
% %     Nid_ind=[ind(:,1) ind(:,1)];
% % elseif  symmetry==1
% %     Th1=Th/2;Wi1=0;
% end
% Nid_ind

% Nid_ind=[ind(1:m2,n1) ind(1:m2,n2)];

NodeIds = nodes_on_line(XYH,P_mid_w,P_center,1e-5);
y=XYH(:,2);Y=y(NodeIds);

Tf2=interp1(PCF(:,2),PCF(:,1),0.25,'linear');
Tf=[TS Tf2 TL];
for i=1:nDB
    load([Folder handles.path_sep 'DB' int2str(i-1) '.mat']);
    axiald(i)=Axial;
    Temp(i,1)=T(NodeId_cen);Temp(i,2)=T(NodeId_cor);
    Temp(i,3)=T(NodeId_mw);Temp(i,4)=T(NodeId_mt);
    
%     Te(:,1)=T(Nid_ind(:,1));Te(:,2)=T(Nid_ind(:,2));
%     Te2=reduce_temp_vec(Te,even_add);
%     x=reduce_temp_vec(Y,even_add);
Te2=T(NodeIds);x=Y;
    d1(i)=find_pos_solidus_line(Te2,x,Tf(1));
    d2(i)=find_pos_solidus_line(Te2,x,Tf(2));
    d3(i)=find_pos_solidus_line(Te2,x,Tf(3));
end
cla(handles.temp);
axes(handles.temp);
set(handles.temp,'Units','normalized','Visible','on');
if nDB>1    
    [axiald' Temp(:,3)]
    plot(axiald,Temp(:,1),':r','LineWidth',2.5);hold on
    plot(axiald,Temp(:,2),'-b','LineWidth',2.5);hold on
    plot(axiald,Temp(:,3),'--g','LineWidth',2.5);hold on
    plot(axiald,Temp(:,4),'--c','LineWidth',2.5);hold on
    xlabel('Distance from Meniscus [m]', 'FontSize',12);
    ylabel('Temperature [^oC]', 'FontSize',12);
    set(gca,'LineWidth',2,'FontSize',12,'Box','off');
    LG{1}='Center';LG{2}='Corner';LG{3}='Width midsurface';LG{4}='Thickness midsurface';
        %%%%%%%%%%%%%% Mark Positions
    v=axis;y=v(3:4);xs=v(2);
     x1=[Mol_Len Mol_Len]./1e3;
     if x1(1)<xs
         line(x1,y,'LineWidth',1.0,'Color','k','LineStyle','-.'); hold on
         text(x1(1),mean(y),'mold end','Rotation',90,'FontSize',10);
     end
     for i=2:size(Zone_pos,1)
         x1=[Zone_pos(i,1) Zone_pos(i,1)]./1e3;
         if x1(1)<xs
             line(x1,y,'LineWidth',1.0,'Color','k','LineStyle','-.'); hold on
             text(x1(1),mean(y),[Zones{i} '  start'],'Rotation',90,'FontSize',10);
         end
         x1=[Zone_pos(i,2) Zone_pos(i,2)]./1e3;
         if x1(1)<xs
             line(x1,y,'LineWidth',1.0,'Color','k','LineStyle','-.'); hold on
             text(x1(1),mean(y),[Zones{i} '  end'],'Rotation',90,'FontSize',10);
         end
     end
%   ['Steel-' num2str(A.Grade_Number) '_' A.Grade_Type '-Vc-' num2str(A.Casting_speed) A.Cool_Mode 'Cooling']

% isstr(A.Grade_Type)
     title(['Steel - ' num2str(A.Grade_Number) '-' num2str(A.Grade_Type) ' : Vc- ' num2str(A.Casting_speed) ' : ' A.Cool_Mode ' Cooling']);
    %%%%%%%%%%%%%% Mark Positions
    hg2=legend(LG,'Location','NorthEast');
    set(hg2,'FontSize',12);
%     [0 1.05*Axial  min(Temp(:,2))-50 A.Init_Temp]
%     axis([0 1.05*Axial  min(Temp(:,2))-50 A.Init_Temp]);
%end
% da=[axiald' d2']
  
cla(handles.shell);
axes(handles.shell);
set(handles.shell,'Units','normalized','Visible','on');
%if nDB>1
    plot(axiald,d1,'-b','LineWidth',2.5);hold on
    plot(axiald,d2,'.-k','LineWidth',2.5);hold on
    plot(axiald,d3,'-r','LineWidth',2.5);hold on
    xlabel('Distance from Meniscus [m]','FontSize',12);
    ylabel('Shell Thickness [mm]','FontSize',12);
    set(gca,'LineWidth',2,'FontSize',12,'Box','off');
    LG2{1}='Solidus Line';LG2{2}='75% solid line'; LG2{3}='Liquidus Line';
    hg2=legend(LG2,'Location','NorthEast');
    set(hg2,'FontSize',12);
        %%%%%%%%%%%%%% Mark Positions
    v=axis;y=v(3:4);xs=v(2);
     x1=[Mol_Len Mol_Len]./1e3;
     if x1(1)<xs
         line(x1,y,'LineWidth',1.0,'Color','k','LineStyle','-.'); hold on
         text(x1(1),mean(y),'mold end','Rotation',90,'FontSize',10);
     end
     for i=2:size(Zone_pos,1)
         x1=[Zone_pos(i,1) Zone_pos(i,1)]./1e3;
         if x1(1)<xs
             line(x1,y,'LineWidth',1.0,'Color','k','LineStyle','-.'); hold on
             text(x1(1),mean(y),[Zones{i} '  start'],'Rotation',90,'FontSize',10);
         end
         x1=[Zone_pos(i,2) Zone_pos(i,2)]./1e3;
         if x1(1)<xs
             line(x1,y,'LineWidth',1.0,'Color','k','LineStyle','-.'); hold on
             text(x1(1),mean(y),[Zones{i} '  end'],'Rotation',90,'FontSize',10);
         end
     end
    %%%%%%%%%%%%% Mark Positions
%     [0 1.05*Axial 0 max(y)/2]
%     axis([0 1.05*Axial 0 max(Y)/2]);
end



function [NodeId,dis] = node_at_pt(XY,Pt)
%  Finds the nNodeId at the given point
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  INPUT: 
%    [XY]    : Nodal coordinates
%    [Pt]    : Point coordinates
%  OUTPUT:
%    [NodeId] : Node id
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nn = size(XY,1); % Number of elements
% NodeId = 0;      % return 0 if no node is found!
Nid=zeros(nn,1);
for in=1:nn    % loop over nodes
    Nid(in)=norm(Pt-XY(in,:));
end % if
[dis,b]=min(Nid);
NodeId=b;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[T2]=reduce_temp_vec(T,even_add)
n1=size(T,1);
T1=0.5*(T(:,1)+T(:,2));
T2=T1;
if even_add==0
    T2(n1)=[];
    T2(n1-1)=0.5*(T1(n1)+T1(n1-1));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[ind]=Contour_node_index(MAPH,RSe)
map=[3 4 
     2 1];
ner1=RSe(1,1);nes1=RSe(1,2);
ect=0;ind=zeros(nes1+1,ner1+1);
for i=1:ner1
    c(1)= i;c(2)= i+1;
    for j=1:nes1
        r(1) = j; r(2) = j+1;
        ect=ect+1;
        mapL=MAPH(ect,:);
        ind(r,c)=mapL(map);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[d] = find_pos_solidus_line(Te,Y,Ts)
%%% finds the vertical position of solidus temperature
n1=1;n2=2;R=0; % default condition
for i=1:length(Te)-1
    T1=Te(i);T2=Te(i+1);
    if T2~=T1
        r=(Ts-T1)/abs(T2-T1);
    else
        r=10;
    end
    if r>=0 && r<=1
       n1=i;n2=i+1;R=r;
    end
    if T1<=Ts && T2<=Ts
        n1=i+1;n2=i+1;R=0;
    end
end
d1=Y(n1);d2=Y(n2);
d=max(Y)-(d1+R*(d2-d1));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function NodeIds2 = nodes_on_line(XY,StPt,EndPt,TOL)
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
%%%% points inside the line starting from stpt
NodeIds2=[];ct=0;dis=[];
for j=1:length(NodeIds)
    xy=XY(NodeIds(j),:);
    d1=norm(StPt-xy);
    d2=norm(EndPt-xy);
    d=norm(StPt-EndPt);
    if (d1+d2<=d)
        NodeIds2(ct+1)=NodeIds(j);
        dis(ct+1)=d1;
        ct=ct+1;
    end
end
[a,b]=sort(dis);
NodeIds2=NodeIds2(b);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% --- Executes during object creation, after setting all properties.
function Panel_Results_Line_Plot_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Panel_Results_Line_Plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function temp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to temp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate temp


% --- Executes on button press in Results_Line_saveTemp.
function Results_Line_saveTemp_Callback(hObject, eventdata, handles)
% hObject    handle to Results_Line_saveTemp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Cur_Fol=[handles.Folders.Main_Fold.Sim handles.Folders.Sub_Fold.Res_Fold];
Grade_Nr=handles.Parameters.Grade_depend.Grade_Number;
Grade_Ty=handles.Parameters.Grade_depend.Grade_Type;
File_Sugg=[Cur_Fol handles.path_sep num2str(Grade_Nr) '_' Grade_Ty '_Temp.fig'];
[FName, PName]=uiputfile(File_Sugg,'Save Figure');
if isequal(FName,0) || isequal(PName,0)
else
    f = figure;
    copyobj([handles.temp],f);
    set(gca,'unit','normalized','position',[0.13 0.11 0.85 0.815]);
    saveas(f,[PName FName]);
    close(gcf);    
end


% --- Executes on button press in Results_Line_saveShell.
function Results_Line_saveShell_Callback(hObject, eventdata, handles)
% hObject    handle to Results_Line_saveShell (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Cur_Fol=[handles.Folders.Main_Fold.Sim handles.Folders.Sub_Fold.Res_Fold];
Grade_Nr=handles.Parameters.Grade_depend.Grade_Number;
Grade_Ty=handles.Parameters.Grade_depend.Grade_Type;
Vc=handles.Parameters.Grade_depend.Casting_speed;
Cool_Mode=handles.Parameters.Grade_depend.Cool_Mode;
File_Sugg=[Cur_Fol handles.path_sep num2str(Grade_Nr) '_' Grade_Ty '_Shell.fig'];
[FName, PName]=uiputfile(File_Sugg,'Save Figure');
if isequal(FName,0) || isequal(PName,0)
else
    f = figure;
    copyobj([handles.shell],f);
    title(['Steel - ' num2str(Grade_Nr) '-' num2str(Grade_Ty) ' : Vc- ' num2str(Vc) ' m/min : ' Cool_Mode ' Cooling']);
    set(gca,'unit','normalized','position',[0.13 0.11 0.85 0.815]);
    saveas(f,[PName FName]);
    close(gcf);    
end


% --- Executes on selection change in Results_DataList.
function Results_DataList_Callback(hObject, eventdata, handles)
% hObject    handle to Results_DataList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Results_DataList contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Results_DataList
List=get(hObject,'String');
Lsel=get(hObject,'Value');
set(handles.Results_AxialDist,'String',List(Lsel,:));
AXIAL=handles.AXIAL;
handles.dbindex=handles.DB{AXIAL(Lsel,2)};
handles.daxial=List(Lsel,:);
guidata(hObject, handles);
if get(handles.Con_1,'Value')==1
    Temp_contour(handles);
else
    Color_plot(handles);
end

% --- Executes during object creation, after setting all properties.
function Results_DataList_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Results_DataList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Results_Contour_SaveLateral.
function Results_Contour_SaveLateral_Callback(hObject, eventdata, handles)
% hObject    handle to Results_Contour_SaveLateral (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Cur_Fol=[handles.Folders.Main_Fold.Sim handles.Folders.Sub_Fold.Res_Fold];
Grade_Nr=handles.Parameters.Grade_depend.Grade_Number;
Grade_Ty=handles.Parameters.Grade_depend.Grade_Type;
Vc=handles.Parameters.Grade_depend.Casting_speed;
Cool_Mode=handles.Parameters.Grade_depend.Cool_Mode;
daxial=handles.daxial;
File_Sugg=[Cur_Fol handles.path_sep num2str(Grade_Nr) '_' Grade_Ty '_Lateral_' num2str(daxial) 'm.fig'];
[FName, PName]=uiputfile(File_Sugg,'Save Figure');
if isequal(FName,0) || isequal(PName,0)
else
    f = figure;
    copyobj([handles.Lateral],f);
    title(['Steel - ' num2str(Grade_Nr) '-' num2str(Grade_Ty) ...
        ' : Vc- ' num2str(Vc) ' m/min : ' Cool_Mode ' Cooling : at ' num2str(daxial) ' m']);
    set(gca,'unit','normalized','position',[0.13 0.11 0.80 0.80]);
    saveas(f,[PName FName]);
    close(gcf);    
end


% --- Executes on button press in Results_Contour_SaveAxial.
function Results_Contour_SaveAxial_Callback(hObject, eventdata, handles)
% hObject    handle to Results_Contour_SaveAxial (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Cur_Fol=[handles.Folders.Main_Fold.Sim handles.Folders.Sub_Fold.Res_Fold];
Grade_Nr=handles.Parameters.Grade_depend.Grade_Number;
Grade_Ty=handles.Parameters.Grade_depend.Grade_Type;
Vc=handles.Parameters.Grade_depend.Casting_speed;
Cool_Mode=handles.Parameters.Grade_depend.Cool_Mode;
daxial=handles.daxial;
File_Sugg=[Cur_Fol handles.path_sep num2str(Grade_Nr) '_' Grade_Ty '_Axial.fig'];
[FName, PName]=uiputfile(File_Sugg,'Save Figure');
if isequal(FName,0) || isequal(PName,0)
else
    f = figure;
    copyobj([handles.Axial],f);
    title(['Steel - ' num2str(Grade_Nr) '-' num2str(Grade_Ty) ...
        ' : Vc- ' num2str(Vc) ' m/min : ' Cool_Mode ' Cooling']);
    set(gca,'unit','normalized','position',[0.13 0.11 0.80 0.80]);
    saveas(f,[PName FName]);
    close(gcf);    
end


function[AXIAL,DATA]=find_axialdis_DB(handles)

Folder_Name=[handles.Folders.Main_Fold.Sim handles.Folders.Sub_Fold.Res_Fold];
d=dir([Folder_Name handles.path_sep 'DB*.mat']);
DATA = {d.name};
N=size(DATA,2);
AXIAL=zeros(N,2);
a=[Folder_Name handles.path_sep];
for i=1:N
    b=DATA{i};
    c=strcat(a,b);
    load(c);
    AXIAL(i,1)=round(Axial*100)/100;
end
[AXIAL(:,1),AXIAL(:,2)]=sort(AXIAL(:,1));

set(handles.Results_DataList,'String',AXIAL(:,1));


function[]=Temp_contour(handles)
%% temperatur contours
Folder_Name=[handles.Folders.Main_Fold.Sim handles.Folders.Sub_Fold.Res_Fold];
a=[Folder_Name handles.path_sep];
c=strcat(a,handles.dbindex);
load(c);

m=handles.Parameters;
load([Folder_Name handles.path_sep 'Material.mat']);
zh=T;
xh=m.XYH(:,1);yh=m.XYH(:,2);
xlin=linspace(min(xh),max(xh),150);
ylin=linspace(min(yh),max(yh),150);
[Xh,Yh]=meshgrid(xlin,ylin);

F = TriScatteredInterp(xh,yh,zh);
Zh=F(Xh,Yh);

% Zh=griddata(xh,yh,zh,Xh,Yh,'natural');

T25=interp1(PCF(:,2),PCF(:,1),0.25,'linear')
% T30=0.30.*(TL-TS)+TS;
vh=[700:100:1400 round([TS T25 TL])];
cla(handles.Lateral);
axes(handles.Lateral);
[C,h]=contour(Xh,Yh,Zh,vh,'LineWidth',2.5);hold on
text_handle=clabel(C,h,'Fontsize',12,'LabelSpacing',72*5);
set(gca,'LineWidth',2,'FontSize',12);
% title([num2str(Grade_Nr) ' - ' Grade_Ty ' - ' num2str(Casting_speed) ' - ' cool_mode],...
%     'FontSize', 16);
% plot_outer_surface(BCEL,BCFC,MAPH,XYH);
% axis()
xlabel('Width [m]', 'FontSize',12);
ylabel('Thickness [m]', 'FontSize',12);
axis equal
axis tight


function[]=Plot_Axial_Contours(handles)
ndb=length(get(handles.Results_DataList,'String'));
m=handles.Parameters;
Folder_Name=[handles.Folders.Main_Fold.Sim handles.Folders.Sub_Fold.Res_Fold];
load([Folder_Name handles.path_sep 'Material.mat']);
VF=1:1:m.nes+1;
VB=m.ner*(m.nes+1)+1:1:(m.ner+1)*(m.nes+1);
VC=(0.5*m.ner)*(m.nes+1)+1:1:((0.5*m.ner)+1)*(m.nes+1);
HT=(m.nes+1):(m.nes+1):(m.ner+1)*(m.nes+1);
HC=(0.5*m.ner)+1:(m.nes+1):((m.nes)+1)*(m.ner+1);
HB=1:(m.nes+1):(m.ner+1)*(m.nes+1);
SEL=get(handles.Results_Midplane,'Value');
% if SEL==1
%     Opt=6; %% vertical center
% else
%     Opt=3; %% vertical front
% end
% if Opt==1
%     Choice=HT;coord=1;
% elseif Opt==2
%     Choice=HB;coord=1;
% elseif Opt==3
%     Choice=VF;coord=2;
% elseif Opt==4
%     Choice=VB;coord=2;
% elseif Opt==5
%     Choice=HC;coord=1;
% elseif Opt==6
%     Choice=VC;coord=2;
% end

coord=2;
if SEL==1
    Opt=6; %% vertical center
    Choice=(0.5*m.ner)*(m.nes+1)+1:1:((0.5*m.ner)+1)*(m.nes+1);
    if m.symmetry==3
        Choice=1:1:m.nes+1;
    end
else
    Opt=3; %% vertical front
    Choice=m.ner*(m.nes+1)+1:1:(m.ner+1)*(m.nes+1);
end



    cla(handles.Axial);
    axes(handles.Axial);
    if ndb>1
        for i=1:ndb
            load([Folder_Name handles.path_sep 'DB' int2str(i-1) '.mat']);
            axiald(i)=Axial;
            TEMP(:,i)=T(Choice);
        end
        xh=axiald;yh=m.XYH(Choice,coord);
        xlin=linspace(min(xh),max(xh),150);
        ylin=linspace(min(yh),max(yh),150);
        
        
%         [Xh,Yh]=meshgrid(xlin,ylin);
%         zh=TEMP;
%         % Axial_contour(XX,YY.*1e3,TEMP,Folder_Name);
%         Zh=griddata(xh,yh,zh,Xh,Yh,'natural');
%         
        
        zh=TEMP';
        [X1, X2] = ndgrid(xh,yh);
        F = griddedInterpolant(X1,X2,zh, 'linear');
        [Xh, Yh] = ndgrid(xlin,ylin);
        Zh = F(Xh, Yh);
        
        
        
        T25=interp1(PCF(:,2),PCF(:,1),0.25,'linear');
        vh=[700:100:1400 round([TS T25 TL])];
        if get(handles.AC_Color_Yes,'Value')==1
            [C,h]=contourf(Xh,Yh,Zh,vh,'LineWidth',2.5);hold on
        else
            [C,h]=contour(Xh,Yh,Zh,vh,'LineWidth',2.5);hold on
        end
        text_handle=clabel(C,h,'Fontsize',12,'LabelSpacing',72*5);
        set(gca,'LineWidth',2,'FontSize',12);
        xlabel('Axial [m]', 'FontSize',12);
        ylabel('Lateral [mm]', 'FontSize',12);
%         axis equal
axis([min(xh) max(xh) min(yh) max(yh)]);
%         axis tight
    end
    
    % --- Executes on button press in Results_Midplane.
function Results_Midplane_Callback(hObject, eventdata, handles)
% hObject    handle to Results_Midplane (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Results_Midplane
% handles.Results_Midplane=get(hObject,'Value');
% guidata(hObject,handles);
Plot_Axial_Contours(handles);


% --- Executes on button press in Results_Surface.
function Results_Surface_Callback(hObject, eventdata, handles)
% hObject    handle to Results_Surface (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Results_Surface
Plot_Axial_Contours(handles);


% --- Executes on button press in save.
function save_Callback(hObject, eventdata, handles)
% hObject    handle to save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Fld=handles.Fld;
KT=handles.KT;
RT=handles.RT;
CT=handles.CT;
ENTH=handles.ENTH;
SF=handles.SF;
D=dir('Material Data Bank\');
b={D.name};
% x=length(b)
% for a=4:x
% b{a}
% if(Fld==b{a})
%     errordlg('Cannot replace the original data, Please enter different Grade number', 'error')
% end
% end
x=sum(strcmp(Fld,b));
if x==1
   errordlg('This steel grade already exists. Please try with different name', 'error');
% else
%     save(KT,)
end

function NewGrade_Load_Callback(hObject, eventdata, handles)
% hObject    handle to NewGrade_Load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Fold_Mat=handles.Folders.Main_Fold.Mat;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% a1=[Fold_Mat '1045_S275JR\' '1045_thermal_conduct.txt'];
% a2=[Fold_Mat '1045_S275JR\' '1045_density.txt'];
% a3=[Fold_Mat '1045_S275JR\' '1045_specific_heat.txt'];
% a4=[Fold_Mat '1045_S275JR\' '1045_enthalpy.txt'];
% a5=[Fold_Mat '1045_S275JR\' '1045_phases.txt'];
% set(handles.NewGrade_SelectData_Conductivity,'String',a1);
% set(handles.NewGrade_SelectData_Density,'String',a2);
% set(handles.NewGrade_SelectData_SpecificHeat,'String',a3);
% set(handles.NewGrade_SelectData_Enthalpy,'String',a4);
% set(handles.NewGrade_SelectData_Phases,'String',a5);

a1=get(handles.NewGrade_SelectData_Conductivity,'String');
a2=get(handles.NewGrade_SelectData_Density,'String');
a3=get(handles.NewGrade_SelectData_SpecificHeat,'String');
a4=get(handles.NewGrade_SelectData_Enthalpy,'String');
a5=get(handles.NewGrade_SelectData_Phases,'String');

a11=importdata(a1);KT=a11.data;
a22=importdata(a2);RT=a22.data;
a33=importdata(a3);CT=a33.data;
a44=importdata(a4);ENTH=a44.data;
a55=importdata(a5);SF=a55.data;
vprop=get(handles.NewGrade_Properties,'Value');
n=1;
for i=1:size(SF,1)
        if SF(i,2)>=0
            n=i;
        end
end
a=SF(1:n,[1 2]);

SF=a;
n1=1;n2=1;
for i=1:size(a,1)-1
    if (a(i+1,2)>0 && a(i,2)<=1e-2) || a(i+1,2)<=1e-2 && a(i,2)>0
        n2=i+1;
    end
    
    if (a(i+1,2)>=100 && a(i,2)<100) || (a(i+1,2)<100 && a(i,2)>=100)
        n1=i;
    end
end

TLi=a(n1);
TSi=a(n2);

set(handles.NewGrade_TS,'String',round(TSi*10)/10);
set(handles.NewGrade_TL,'String',round(TLi*10)/10);

f1=Linear_Interpolation(ENTH,round(TLi*10)/10);
f2=Linear_Interpolation(ENTH,round(TSi*10)/10);
Latent=abs(f1-f2);
set(handles.NewGrade_Latent,'String',round(Latent*10)/10);


% [Table_DData]=Make_Regular_Table_DATA(KHT,SHT,RHT,PCF,ENTH);
% handles.Table_DData=Table_DData;

Properties=struct('KHT',{KT},'SHT',{CT},'RHT',{RT},...
    'TS',{TSi},'TL',{TLi},'Latent',{Latent}, 'PCF',{SF},'ENTH',{ENTH});
handles.Properties=Properties;
guidata(hObject,handles);
cla(handles.NewGrade_Plot);
axes(handles.NewGrade_Plot);
Plot_NewGrade(handles);


% --- Executes on selection change in NewGrade_Properties.
function NewGrade_Properties_Callback(hObject, eventdata, handles)
% hObject    handle to NewGrade_Properties (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% vprop=get(hObject, 'Value')
% handles.mat_file='D_Material.mat';
cla(handles.NewGrade_Plot);
axes(handles.NewGrade_Plot);
Plot_NewGrade(handles);
guidata(hObject, handles);
% Hints: contents = cellstr(get(hObject,'String')) returns NewGrade_Properties contents as cell array
%        contents{get(hObject,'Value')} returns selected item from NewGrade_Properties


% --- Executes during object creation, after setting all properties.
function NewGrade_Properties_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NewGrade_Properties (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function NewGrade_SelectData_Conductivity_Callback(hObject, eventdata, handles)
% hObject    handle to NewGrade_SelectData_Conductivity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of NewGrade_SelectData_Conductivity as text
%        str2double(get(hObject,'String')) returns contents of NewGrade_SelectData_Conductivity as a double


% --- Executes during object creation, after setting all properties.
function NewGrade_SelectData_Conductivity_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NewGrade_SelectData_Conductivity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function NewGrade_SelectData_Density_Callback(hObject, eventdata, handles)
% hObject    handle to NewGrade_SelectData_Density (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of NewGrade_SelectData_Density as text
%        str2double(get(hObject,'String')) returns contents of NewGrade_SelectData_Density as a double


% --- Executes during object creation, after setting all properties.
function NewGrade_SelectData_Density_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NewGrade_SelectData_Density (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function NewGrade_SelectData_SpecificHeat_Callback(hObject, eventdata, handles)
% hObject    handle to NewGrade_SelectData_SpecificHeat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of NewGrade_SelectData_SpecificHeat as text
%        str2double(get(hObject,'String')) returns contents of NewGrade_SelectData_SpecificHeat as a double


% --- Executes during object creation, after setting all properties.
function NewGrade_SelectData_SpecificHeat_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NewGrade_SelectData_SpecificHeat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function NewGrade_SelectData_Enthalpy_Callback(hObject, eventdata, handles)
% hObject    handle to NewGrade_SelectData_Enthalpy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of NewGrade_SelectData_Enthalpy as text
%        str2double(get(hObject,'String')) returns contents of NewGrade_SelectData_Enthalpy as a double


% --- Executes during object creation, after setting all properties.
function NewGrade_SelectData_Enthalpy_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NewGrade_SelectData_Enthalpy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function NewGrade_SelectData_Phases_Callback(hObject, eventdata, handles)
% hObject    handle to NewGrade_SelectData_Phases (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of NewGrade_SelectData_Phases as text
%        str2double(get(hObject,'String')) returns contents of NewGrade_SelectData_Phases as a double


% --- Executes during object creation, after setting all properties.
function NewGrade_SelectData_Phases_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NewGrade_SelectData_Phases (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in NewGrade_Browse_Conductivity.
function NewGrade_Browse_Conductivity_Callback(hObject, eventdata, handles)
% hObject    handle to NewGrade_Browse_Conductivity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Fold_Mat=[handles.Folders.Main_Fold.Mat handles.path_sep '.txt'];
[File,Path]=uigetfile(Fold_Mat,'Load Thermal Conductivity');
if isequal(File,0) || isequal(Path,0)
else
    MAT_FILE=[Path File];
    KT1=importdata(MAT_FILE);
    set(handles.NewGrade_SelectData_Conductivity,'String',MAT_FILE);
    KT=KT1.data;
    handles.KT=KT;
end
guidata(hObject,handles);

% --- Executes on button press in NewGrade_Browse_Density.
function NewGrade_Browse_Density_Callback(hObject, eventdata, handles)
% hObject    handle to NewGrade_Browse_Density (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Fold_Mat=[handles.Folders.Main_Fold.Mat handles.path_sep '.txt'];
[File,Path]=uigetfile(Fold_Mat,'Load Density');
if isequal(File,0) || isequal(Path,0)
else
    MAT_FILE=[Path File];
    RT1=importdata(MAT_FILE);
    set(handles.NewGrade_SelectData_Density,'String',MAT_FILE);
    RT=RT1.data;
    handles.RT=RT;
end
guidata(hObject,handles);


% --- Executes on button press in NewGrade_Browse_SpecificHeat.
function NewGrade_Browse_SpecificHeat_Callback(hObject, eventdata, handles)
% hObject    handle to NewGrade_Browse_SpecificHeat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Fold_Mat=[handles.Folders.Main_Fold.Mat handles.path_sep '.txt'];
[File,Path]=uigetfile(Fold_Mat,'Load Specific Heat');
if isequal(File,0) || isequal(Path,0)
else
MAT_FILE=[Path File];
CT1=importdata(MAT_FILE);
set(handles.NewGrade_SelectData_SpecificHeat,'String',MAT_FILE);
CT=CT1.data;
handles.CT=CT;
end
guidata(hObject,handles);

% --- Executes on button press in NewGrade_Browse_Enthalpy.
function NewGrade_Browse_Enthalpy_Callback(hObject, eventdata, handles)
% hObject    handle to NewGrade_Browse_Enthalpy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Fold_Mat=[handles.Folders.Main_Fold.Mat handles.path_sep '.txt'];
[File,Path]=uigetfile(Fold_Mat,'Load Enthalpy');
if isequal(File,0) || isequal(Path,0)
else
    
    MAT_FILE=[Path File];
    ENTH1=importdata(MAT_FILE);
    set(handles.NewGrade_SelectData_Enthalpy,'String',MAT_FILE);
    ENTH=ENTH1.data;
    handles.ENTH=ENTH;
end
guidata(hObject,handles);

% --- Executes on button press in NewGrade_Browse_Phases.
function NewGrade_Browse_Phases_Callback(hObject, eventdata, handles)
% hObject    handle to NewGrade_Browse_Phases (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Fold_Mat=[handles.Folders.Main_Fold.Mat handles.path_sep '.txt'];
[File,Path]=uigetfile(Fold_Mat,'Load Phase Function');
if isequal(File,0) || isequal(Path,0)
else
    
    MAT_FILE=[Path File];
    SF1=importdata(MAT_FILE);
    set(handles.NewGrade_SelectData_Phases,'String',MAT_FILE);
    SF=SF1.data;
    handles.SF=SF;
end
guidata(hObject,handles);



function NewGrade_GN_Callback(hObject, eventdata, handles)
% hObject    handle to NewGrade_GN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
str_val=get(hObject,'String');
GN=str2double(str_val);
handles.New_Grade.Grade_Number=GN;
GT=get(handles.NewGrade_GT,'String');
Fld=[num2str(GN) '_' GT];
set(handles.NewGrade_FolderName,'String',Fld);
guidata(hObject,handles);

% Hints: get(hObject,'String') returns contents of NewGrade_GN as text
%        str2double(get(hObject,'String')) returns contents of NewGrade_GN as a double


% --- Executes during object creation, after setting all properties.
function NewGrade_GN_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NewGrade_GN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function NewGrade_GT_Callback(hObject, eventdata, handles)
% hObject    handle to NewGrade_GT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.New_Grade.Grade_Type=get(hObject,'String');
GN=str2double(get(handles.NewGrade_GN,'String'));
GT=get(handles.NewGrade_GT,'String');
Fld=[num2str(GN) '_' GT];
set(handles.NewGrade_FolderName,'String',Fld);
guidata(hObject,handles);



% --- Executes during object creation, after setting all properties.
function NewGrade_GT_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NewGrade_GT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function NewGrade_CoolMode_Callback(hObject, eventdata, handles)
% hObject    handle to NewGrade_CoolMode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Cool_Mode_Value=get(hObject,'Value');
str=get(hObject,'String');
Cool_Mode=str{get(hObject,'Value')};
Std_water_share=handles.default.Std_water_share;
Std_cooling=handles.default.Std_cooling;
% Water_share_matrix(1,1)=Std_water_share(1)+Std_water_share(1);
% Water_share_matrix(2,2)=Std_water_share(2);
% Water_share_matrix(3,2)=1-Std_water_share(2);
% Water_share_matrix(4,3)=Std_water_share(3);
% Water_share_matrix(5,3)=1-Std_water_share(3);

if Cool_Mode_Value==5
    set(handles.NewGrade_Zone_1,'Style','edit','BackgroundColor',[1 1 1]);
    set(handles.NewGrade_Zone_2,'Style','edit','BackgroundColor',[1 1 1]);
    set(handles.NewGrade_Zone_3,'Style','edit','BackgroundColor',[1 1 1]);
    WFR_Lm3=Std_cooling(4,:);
%     WFR_Lm3(1)=str2double(get(handles.NewGrade_Zone_1,'String'));
%     WFR_Lm3(2)=str2double(get(handles.NewGrade_Zone_2,'String'));
%     WFR_Lm3(3)=str2double(get(handles.NewGrade_Zone_3,'String'));    
    
else
    set(handles.NewGrade_Zone_1,'Style','text','BackgroundColor',0.941*[1 1 1]);
    set(handles.NewGrade_Zone_2,'Style','text','BackgroundColor',0.941*[1 1 1]);
    set(handles.NewGrade_Zone_3,'Style','text','BackgroundColor',0.941*[1 1 1]);
     WFR_Lm3=Std_cooling(Cool_Mode_Value,:);
     set(handles.NewGrade_Zone_1,'String',WFR_Lm3(1));
     set(handles.NewGrade_Zone_2,'String',WFR_Lm3(2));
     set(handles.NewGrade_Zone_3,'String',WFR_Lm3(3));
end
handles.New_Grade.Cool_Mode_Value=Cool_Mode_Value;
handles.New_Grade.WFR_Lm3=WFR_Lm3;
handles.New_Grade.Cool_Mode=Cool_Mode;
guidata(hObject,handles);

% Hints: contents = cellstr(get(hObject,'String')) returns NewGrade_CoolMode contents as cell array
%        contents{get(hObject,'Value')} returns selected item from NewGrade_CoolMode


% --- Executes during object creation, after setting all properties.
function NewGrade_CoolMode_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NewGrade_CoolMode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function NewGrade_CS_Callback(hObject, eventdata, handles)
% hObject    handle to NewGrade_CS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of NewGrade_CS as text
%        str2double(get(hObject,'String')) returns contents of NewGrade_CS as a double
CSS=get(hObject,'String');
CS=str2double(CSS);
handles.New_Grade.Casting_speed=CS;
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function NewGrade_CS_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NewGrade_CS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function NewGrade_CT_Callback(hObject, eventdata, handles)
% hObject    handle to NewGrade_CT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of NewGrade_CT as text
%        str2double(get(hObject,'String')) returns contents of NewGrade_CT as a double
CTS=get(hObject,'String');
CT=str2double(CTS);
handles.New_Grade.Init_Temp=CT;
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function NewGrade_CT_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NewGrade_CT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function[]=Plot_NewGrade(handles)
vprop=get(handles.NewGrade_Properties,'Value');
m=handles.Properties;
switch vprop
    case 1
        %%%%% Thermal conductivity [KHT1]
        plot(m.KHT(:,1),m.KHT(:,2),'r','LineWidth',2);
        ylabel('Thermal Conductivity [W/mK]','FontSize',12);
        xlabel('Temperature [degC]','FontSize',12);
    case 2
        %%%%% Density [RHT1]
        plot(m.RHT(:,1),m.RHT(:,2),'b','LineWidth',2);
        ylabel('Density [g/cm^3]','FontSize',12);
        xlabel('Temperature [degC]','FontSize',12);
    case 3
        
        %%%%% Specific Heat Capacity [SHT1]
        plot(m.SHT(:,1),m.SHT(:,2),'k','LineWidth',2);
        ylabel('Specific heat [J/gK]','FontSize',12);
        xlabel('Temperature [degC]','FontSize',12);
        
    case 4
        plot(m.PCF(:,1),m.PCF(:,2)./1e2,'m','LineWidth',2);
        ylabel('Liquid fraction [-]','FontSize',12);
        xlabel('Temperature [degC]','FontSize',12);
    case 5
        plot(m.ENTH(:,1),m.ENTH(:,2),'m','LineWidth',2);
        ylabel('Enthalpy [J/g]','FontSize',12);
        xlabel('Temperature [degC]','FontSize',12);
end
set(gca,'FontSize',12,'LineWidth',2.0,'Box','off');

function NewGrade_Save_Callback(hObject, eventdata, handles)
% hObject    handle to NewGrade_Save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Fold_Mat=handles.Folders.Main_Fold.Mat;
Property_Tag=handles.default.Property_Tag;
GN=handles.New_Grade.Grade_Number;
FolderName=get(handles.NewGrade_FolderName,'String');
load([handles.Folders.Main_Fold.Grd  handles.Folders.File.Grd  '.mat']);
% Steel_reciepes_SC=handles.Steel_reciepes_SC;
% handles.Properties
KHT=handles.Properties.KHT;
RHT=handles.Properties.RHT;
SHT=handles.Properties.SHT;
ENTH=handles.Properties.ENTH;
PCF=handles.Properties.PCF;
Std_water_share=handles.default.Std_water_share;

% %%%%%%%% remove finally
% % %%%%%%%%just for running purpose
% if length(Steel_reciepes_SC)>84
%     Steel_reciepes_SC(85:length(Steel_reciepes_SC))=[];
% end
% %%%%%%%%%%%%%%%%%%%%55555
% %%%%%%%%%%%%%%%%5


WF_matrix_3x5=zeros(3,5);
WF_matrix_3x5(1,1)=1;
WF_matrix_3x5(2,2)=Std_water_share(2);WF_matrix_3x5(2,3)=1-Std_water_share(2);
WF_matrix_3x5(3,4)=Std_water_share(3);WF_matrix_3x5(3,5)=1-Std_water_share(3);

d=dir(Fold_Mat);
str={d.name};
if length(str)>2
    str=str(3:length(str));
end
Exist=strcmp(str,FolderName);
if sum(Exist)>0
    errordlg('Steel Grade Already Exists. Change Grade Type or Number', 'Grade Already Exist');
else
    FOLDER =[Fold_Mat FolderName];
    if isdir(FOLDER)==0
        mkdir(FOLDER);
    end    
%     FOLDER =[Data_FN handles.path_sep num2str(Grade_Nr) '_' Grade_Ty handles.path_sep];
    FILE1=[FOLDER handles.path_sep num2str(GN) '_' Property_Tag.RHT];
    FILE2=[FOLDER handles.path_sep num2str(GN) '_' Property_Tag.ENTH];
    FILE3=[FOLDER handles.path_sep num2str(GN) '_' Property_Tag.PCF];
    FILE4=[FOLDER handles.path_sep num2str(GN) '_' Property_Tag.SHT];
    FILE5=[FOLDER handles.path_sep num2str(GN) '_' Property_Tag.KHT];
    fid1 = fopen(FILE1,'w');
    fid2 = fopen(FILE2,'w');fid3 = fopen(FILE3,'w');
    fid4 = fopen(FILE4,'w');fid5 = fopen(FILE5,'w');
    fprintf(fid1,'%s %s \r\n ','Temperature [degC]', 'Density [g/m^3]');
    fprintf(fid1,'%6.2f %6.2f \r\n ',RHT');
    fclose(fid1);
    fprintf(fid2,'%s %s \r\n ','Temperature [degC]', 'Enthalpy [J/g]');
    fprintf(fid2,'%6.2f %6.2f \r\n ',ENTH');
    fclose(fid2);
    fprintf(fid3,'%s %s \r\n ','Temperature [degC]', 'Liquid Fraction [-]');
    fprintf(fid3,'%6.2f %6.2f \r\n ',PCF');
    fclose(fid3);
    fprintf(fid4,'%s %s \r\n ','Temperature [degC]', 'Specific heat [J/gK]');
    fprintf(fid4,'%6.2f %6.2f \r\n ',SHT');
    fclose(fid4);
    fprintf(fid5,'%s %s \r\n ','Temperature [degC]', 'Thermal Conductivity[W/mK]');
    fprintf(fid5,'%6.2f %6.2f \r\n ',KHT');
    fclose(fid5);    
    
    n=length(Steel_reciepes_SC);
    Steel_reciepes_SC(n+1).Grade_Number=handles.New_Grade.Grade_Number;
    Steel_reciepes_SC(n+1).Grade_Type=handles.New_Grade.Grade_Type;
    Cooling_mode=Invers_cool_mode(handles.New_Grade.Cool_Mode_Value);
    Steel_reciepes_SC(n+1).Cooling_mode=Cooling_mode;
    Water_flow_rate=handles.New_Grade.WFR_Lm3*WF_matrix_3x5;
    Steel_reciepes_SC(n+1).Water_flow_rate=Water_flow_rate;
    Steel_reciepes_SC(n+1).Vc_target=handles.New_Grade.Casting_speed;
    Steel_reciepes_SC(n+1).T_init=handles.New_Grade.Init_Temp;
%     Steel_reciepes_SC(n+1)
    Grd_File_Name=[handles.Folders.Main_Fold.Grd handles.Folders.File.Grd];
    save(Grd_File_Name,'Steel_reciepes_SC');
    %%% Update grade list in Default window
    GNr_new={Steel_reciepes_SC.Grade_Number};
    GNr_old = cellstr(get(handles.GNList,'String'));
    GN=str2double(GNr_old{get(handles.GNList,'Value')});
    dvalue=find(cell2mat(GNr_new)==GN);
    set(handles.GNList,'String',GNr_new,'Value',dvalue);
%     handles.Steel_reciepes_SC=Steel_reciepes_SC;
    guidata(hObject,handles);
end

function[CM]= Invers_cool_mode(CMV)
if CMV==1
    CM=[1 0 0 0];
elseif CMV==2
    CM=[0 1 0 0];
elseif  CMV==3
    CM=[0 0 1 0];
elseif  CMV==4
    CM=[0 0 0 1];
else
    CM=[0 0 0 0 1];
end


function NewGrade_Zone_1_Callback(hObject, eventdata, handles)
% hObject    handle to NewGrade_Zone_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of NewGrade_Zone_1 as text
%        str2double(get(hObject,'String')) returns contents of NewGrade_Zone_1 as a double
 handles.New_Grade.WFR_Lm3(1)=str2double(get(hObject,'String'));
 guidata(hObject,handles);
 
 function NewGrade_Zone_2_Callback(hObject, eventdata, handles)
% hObject    handle to NewGrade_Zone_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of NewGrade_Zone_2 as text
%        str2double(get(hObject,'String')) returns contents of NewGrade_Zone_2 as a double
 handles.New_Grade.WFR_Lm3(2)=str2double(get(hObject,'String'));
 guidata(hObject,handles);



function NewGrade_Zone_3_Callback(hObject, eventdata, handles)
% hObject    handle to NewGrade_Zone_3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of NewGrade_Zone_3 as text
%        str2double(get(hObject,'String')) returns contents of NewGrade_Zone_3 as a double
 handles.New_Grade.WFR_Lm3(3)=str2double(get(hObject,'String'));
 guidata(hObject,handles);

 
 
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


% --- Executes when selected object is changed in Mold_Model.
% function Mold_Model_SelectionChangeFcn(hObject, eventdata, handles)
% % hObject    handle to the selected object in Mold_Model 
% % eventdata  structure with the following fields (see UIBUTTONGROUP)
% %	EventName: string 'SelectionChanged' (read only)
% %	OldValue: handle of the previously selected object or empty if none was selected
% %	NewValue: handle of the currently selected object
% % handles    structure with handles and user data (see GUIDATA)
% if (hObject == handles.MB)
%     handles.present.index=1;
% elseif (hObject == handles.MH)
%     handles.present.index=2;
% elseif (hObject == handles.MA)
%     handles.present.index=3;
% else
%     handles.present.index=4;
% end
% guidata(hObject,handles);
% 
% cla(handles.Mold_Plot);
% axes(handles.Mold_Plot);
% [z,hf]=mold_heat_flux(handles.present);
% plot(z,hf./1e6,'-b','LineWidth',3);hold on
% set(gca,'FontSize',12,'LineWidth',2.0,'Box','off');
% ylabel(' Heat Flux [MW/m^2]','FontSize',13);
% xlabel(' Distance from Meniscus [mm]','FontSize',13);
% handles.present.z=z;
% handles.present.hf=hf;
% guidata(hObject,handles);





    function[Folders]=Create_Folder_Struct(Fold_Path,default,path_sep,root)
%         Fold_Path=[Fold_Path path_sep];
%  'Mat',{[default.Mat_Fold_Name path_sep]},...
        Main_Fold=struct('Sim',{[root path_sep default.Sim_Main_Fold path_sep]},...
                 'Mat',{[root path_sep 'UDM' path_sep]},...
                 'Grd',{[default.Grd_Fold_Name path_sep]},...
                 'Set',{[root path_sep default.Set_Fold_Name path_sep]},...
                 'Pro',{[default.Pro_Fold_Name path_sep]},...
                  'Mec',{['MECH' path_sep]});
 File=struct('Grd',{default.Grd_File_Name},...
                 'SDef',{default.Def_File_Name},...
                 'SUse',{default.Def_File_Name},...
                  'MBank',{default.Mat_Data_Bank});
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%             
 Sub_Fold=struct('Sim_Fold',{''},...
                 'Res_Fold',{''},...
                 'Mat_Fold',{''},...
                  'Mec_Fold',{'MECH'});
%%%%%%%%%%%%%%% checking %%%%%%%%%%%%%%%%%%%%%%%%    
%  Sub_Fold=struct('Sim_Fold',{'1015_S235JR_2'},...
%                  'Res_Fold',{'1015_S235JR_2'},...
%                  'Mat_Fold',{''},...
%                  'Mec_Fold',{'MECH'});
 %%%%%%%%%%%%%%% checking %%%%%%%%%%%%%%%%%%%%%%%%            
  Folders=struct('Main_Fold',{Main_Fold},...
                 'Sub_Fold',{Sub_Fold},...
                 'File',{File});


% --- Executes when selected object is changed in Mold_model.
function Mold_model_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in Mold_model 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
set(handles.MDT_u,'Style','text','BackgroundColor',0.941*[1 1 1]);
set(handles.MMWF_u,'Style','text','BackgroundColor',0.941*[1 1 1]);
set(handles.MMWT_u,'Style','text','BackgroundColor',0.941*[1 1 1]);

if (hObject == handles.MB)
    handles.present.index=1;
elseif (hObject == handles.MH)
    handles.present.index=2;
elseif (hObject == handles.MA)
    handles.present.index=3;
    set(handles.MDT_u,'Style','edit','BackgroundColor',[1 1 1]);
    set(handles.MMWF_u,'Style','edit','BackgroundColor',[1 1 1]);
    set(handles.MMWT_u,'Style','edit','BackgroundColor',[1 1 1]);
else
    handles.present.index=4;
end
guidata(hObject,handles);

cla(handles.Mold_Plot);
axes(handles.Mold_Plot);
[MHF]=mold_heat_flux(handles.present);
plot(MHF(:,1),MHF(:,2),'-b','LineWidth',3);hold on
set(gca,'FontSize',12,'LineWidth',2.0,'Box','off');
ylabel(' Heat Flux [MW/m^2]','FontSize',13);
xlabel(' Distance from Meniscus [mm]','FontSize',13);
handles.present.MHF=MHF;
guidata(hObject,handles);


function[]=delete_data_file(pa,Nr)
%%% pa - folder path
%%% Nr - above this number, files will be deleted
% pa='D:\HAZELETT CASTER\SWISS STEEL SOFTWARE\TABBING Version\Final_2-11-04-2014\Simulations\1045_S275JR\'
d=dir([pa 'DB*.mat']);
str={d.name};
ct=0;file_del=[];
% Nr=11;
for i=1:length(str)
    stri=str(i);
    strc=['DB'];
    str2=strrep(stri,strc,'');
    strc=['.mat'];
    str3=strrep(str2,strc,'');
    if str2double(str3)>Nr
        index=i;
        str_del=str{index};
        file_del{ct+1}=[pa str_del];
        ct=ct+1;
    end    
end
for i=1:length(file_del)
    delete(file_del{i});
end


function[BC_Info2]=find_reduction_factor(XYH,MAPH,BC_Info,f11,f22,Wi,Th)
BC_Info2=BC_Info;
for i=1:4
    ct=0;el=[];index=[];
    for j=1:size(BC_Info,1)
        if BC_Info(j,2)==i
            el(ct+1)=BC_Info(j,1);
            index(ct+1)=j;
            ct=ct+1;
        end
    end
    ELS{i}=el;IND{i}=index;
end
N=size(ELS,2);
face=[4 1;3 2;2 1;3 4];
f1=[f11 f22];f2=[f22 f11];
for i=1:N
    els=ELS{i};
    ind=IND{i};
    if i==3 || i==4
        pos=1;d=Wi;
    else
        pos=2;d=Th;
    end
    for j=1:length(els)
        el=els(j);in=ind(j);
        nodes=MAPH(el,face(i,:));
        coord=XYH(nodes,pos);
        z=coord-d/2;
        for k=1:length(z)
            if z(k)<0
                m=1;d1=-d/2;
            else
                m=2;d1=0;
            end
            fk(k)=(2*(f2(m)-f1(m))*(z(k)-d1))/d+f1(m);
        end
        BC_Info2(in,3:4)=fk;
    end    
end


function Color_plot(handles)
%  Draw counter-lines
%  Sub-function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  INPUT:                                                                    %
%    [XY]      -> Nodal coordinates                                          %
%    [MAP]     -> Element connectivity                                       %
%    [DATA]    -> Nodal point data structure                                 %
%  OUTPUT:    Draws                                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Folder_Name=[handles.Folders.Main_Fold.Sim handles.Folders.Sub_Fold.Res_Fold];
a=[Folder_Name handles.path_sep];
c=strcat(a,handles.dbindex);
load(c);
load([Folder_Name handles.path_sep 'Process_Parameters.mat']);
ne = size(MAPH,1); % number of time steps
% figure; % create a new figure for geometry plot
cla(handles.Lateral);
axes(handles.Lateral);
for (ie=1:ne)
  patch(XYH(MAPH(ie,:),1),XYH(MAPH(ie,:),2),T(MAPH(ie,:)));
end
axis('equal');
% axis('off');
axis tight
set(gca,'LineWidth',2,'FontSize',12);
xlabel('Width [m]', 'FontSize',12);
ylabel('Thickness [m]', 'FontSize',12);


% --- Executes when selected object is changed in uipanel163.
function uipanel163_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel163 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
Plot_Axial_Contours(handles);

% --------------------------------------------------------------------
function uipanel163_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to uipanel163 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes when selected object is changed in Cont_Col_1.
function Cont_Col_1_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in Cont_Col_1 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
if (hObject == handles.Con_1)
    Temp_contour(handles)
elseif (hObject == handles.Col_1)
    Color_plot(handles)
end


% --- Executes on selection change in SSymmetry.
function SSymmetry_Callback(hObject, eventdata, handles)
% hObject    handle to SSymmetry (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns SSymmetry contents as cell array
%        contents{get(hObject,'Value')} returns selected item from SSymmetry
handles.present.symmetry=get(hObject,'Value')-1;
[handles.present.XYH,handles.present.MAPH,...
    handles.present.BC_Info]=Create_mesh(handles);
guidata(hObject, handles);
draw_heat_4n(handles);
% --- Executes during object creation, after setting all properties.
function SSymmetry_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SSymmetry (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when selected object is changed in uipanel181.
function uipanel181_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel181 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
Plot_Axial_Contours(handles)


% --- Executes during object creation, after setting all properties.
function Secondary_plot_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Secondary_plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate Secondary_plot
% plot_Secondary_Cooling(handles.present,handles)
% plot_Secondary_Cooling(handles)


function[Nozzle_Parameters,Nozzle_to_Zone,Q_vs_CA]=Collect_Nozzle_DATA()
Nozzle_Parameters=struct(...
    'CAF',{0.7,0.98,0.65},... % Cone Angle Factor
    'ICF1',{0.5,0.55,0.5},... % Inner Circle Factor 1
    'ICF2',{0.5,0.2,0.5},...  % Inner Circle Factor 2
    'STD_FR',{3.5,2.5,2.0},...
    'STD_CA',{65,60,50}); 
Nozzle_to_Zone=[1,1,2,3,3,3]; % 6 zones nozzle number
    
 A{1}=[
     0         0
    0.5800   34.5000
    1.5000   63.5000
    3.5000   65.0000
    4.0000   90.0000
    6.0000   94.0000
    10.000   98.0000];

 A{2}=[
    0         0
    0.5600   39.0000
    1.5000   48.0000
    2.0000   52.0000
    3.5000   53.0000
    4.5000   54.0000
    7.0000   58.0000];

 A{3}=[
         0         0
    0.5500   16.0000
    1.0000   38.0000
    1.5000   45.0000
    2.0000   50.0000
    3.3000   50.0000
    5.0000   52.0000];
Q_vs_CA=A;


% --------------------------------------------------------------------
function the_man_Callback(hObject, eventdata, handles)
% hObject    handle to the_man (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

try
  winopen('Theory_Manual.pdf');
catch error
  msgbox(error.message,'Fehler','error');
end

% --------------------------------------------------------------------
function use_man_Callback(hObject, eventdata, handles)
% hObject    handle to use_man (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
  winopen('User_Manual_1.pdf');
catch error
  msgbox(error.message,'Fehler','error');
end



function SC_Z1L_u_Callback(hObject, eventdata, handles)
% hObject    handle to SC_Z1L_u (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SC_Z1L_u as text
%        str2double(get(hObject,'String')) returns contents of SC_Z1L_u as a double
Z1=str2double(get(hObject,'String'));
CS=handles.present.Grade_depend.Casting_speed;
handles.present.Grade_depend.WFR_Lm3(1)=Z1/CS;
set(handles.SC_Z1_u,'String',Z1/CS);
guidata(hObject,handles);
[SEC_COOL,Nozzle_Tab,No_Noz_Zo,cooling_info]=Updade_Cooling_Info(handles,3);
handles.present.SEC_COOL=SEC_COOL;
handles.present.No_Noz_Zo=No_Noz_Zo;
handles.present.Nozzle_Tab=Nozzle_Tab;
handles.present.cooling_info=cooling_info;

cla(handles.Secondary_plot);
axes(handles.Secondary_plot);
plot_Secondary_Cooling(handles);



function SC_Z3L_u_Callback(hObject, eventdata, handles)
% hObject    handle to SC_Z3L_u (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SC_Z3L_u as text
%        str2double(get(hObject,'String')) returns contents of SC_Z3L_u as a double
Z3=str2double(get(hObject,'String'));
CS=handles.present.Grade_depend.Casting_speed;
handles.present.Grade_depend.WFR_Lm3(3)=Z3/CS;
set(handles.SC_Z3_u,'String',Z3/CS);
guidata(hObject,handles);
[SEC_COOL,Nozzle_Tab,No_Noz_Zo,cooling_info]=Updade_Cooling_Info(handles,3);
handles.present.SEC_COOL=SEC_COOL;
handles.present.No_Noz_Zo=No_Noz_Zo;
handles.present.Nozzle_Tab=Nozzle_Tab;
handles.present.cooling_info=cooling_info;

cla(handles.Secondary_plot);
axes(handles.Secondary_plot);
plot_Secondary_Cooling(handles);



function SC_Z2L_u_Callback(hObject, eventdata, handles)
% hObject    handle to SC_Z2L_u (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SC_Z2L_u as text
%        str2double(get(hObject,'String')) returns contents of SC_Z2L_u as a double
Z2=str2double(get(hObject,'String'));
CS=handles.present.Grade_depend.Casting_speed;
handles.present.Grade_depend.WFR_Lm3(2)=Z2/CS;
set(handles.SC_Z2_u,'String',Z2/CS);
guidata(hObject,handles);
[SEC_COOL,Nozzle_Tab,No_Noz_Zo,cooling_info]=Updade_Cooling_Info(handles,3);
handles.present.SEC_COOL=SEC_COOL;
handles.present.No_Noz_Zo=No_Noz_Zo;
handles.present.Nozzle_Tab=Nozzle_Tab;
handles.present.cooling_info=cooling_info;

cla(handles.Secondary_plot);
axes(handles.Secondary_plot);
plot_Secondary_Cooling(handles);



function TIME_Callback(hObject, eventdata, handles)
% hObject    handle to TIME (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of TIME as text
%        str2double(get(hObject,'String')) returns contents of TIME as a double


% --- Executes during object creation, after setting all properties.
function TIME_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TIME (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function use_man_2_Callback(hObject, eventdata, handles)
% hObject    handle to use_man_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
  winopen('User_Manual_2.pdf');
catch error
  msgbox(error.message,'Fehler','error');
end


% --- Executes on button press in Mech.
function Mech_Callback(hObject, eventdata, handles)
% hObject    handle to Mech (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Mech
Folders=handles.Folders;
Folders.Main_Fold.Sim;
Folders=handles.Folders;
Sim_path=Folders.Main_Fold.Sim;
Sim_Fold=Folders.Sub_Fold.Sim_Fold;
d=dir([Sim_path Sim_Fold '\DB*.mat']);
str = {d.name};
ndir=size(str,2);

    val=1;
    if ndir<=2
        choice = questdlg('Stress Simulation for this alloy NOT POSSIBLE before Temperature Simulation', ...
            'Atttention', ...
            'YES','NO','YES');        
        switch choice
            case 'YES'
                val=1;
            case 'NO'
                val=2;
        end
    end
    
if val==1
   Mech_1(Folders);
end
% Sim_Fold=get(handles.FolderName,'String')


% --- Executes on button press in Adv_Stress_Run.
function Adv_Stress_Run_Callback(hObject, eventdata, handles)
% hObject    handle to Adv_Stress_Run (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Adv_Stress_Run
