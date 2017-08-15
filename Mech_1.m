function varargout = Mech_1(varargin)
% MECH_1 MATLAB code for Mech_1.fig
%      MECH_1, by itself, creates a new MECH_1 or raises the existing
%      singleton*.
%
%      H = MECH_1 returns the handle to a new MECH_1 or the handle to
%      the existing singleton*.
%
%      MECH_1('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MECH_1.M with the given input arguments.
%
%      MECH_1('Property','Value',...) creates a new MECH_1 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Mech_1_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Mech_1_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Mech_1

% Last Modified by GUIDE v2.5 13-Mar-2017 22:39:20

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Mech_1_OpeningFcn, ...
                   'gui_OutputFcn',  @Mech_1_OutputFcn, ...
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


% --- Executes just before Mech_1 is made visible.
function Mech_1_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Mech_1 (see VARARGIN)

% Choose default command line output for Mech_1
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
ah=axes('unit','normalized','position',[0.57 0.2 0.38 0.56]);
bg=imread(['model_fig.jpg']);imagesc(bg);
set(ah,'handlevisibility','off','visible','off');

set(handles.E,'Value',1);
set(handles.Axial_Strain,'Value',1);
set(handles.Curvature_X,'Value',1);
set(handles.Curvature_Y,'Value',1);
% Update handles structure
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Dependent running
% Folders=varargin;
% handles.Folders=Folders;
% Sim_path=Folders.Main_Fold.Sim;
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%% Folder Handlings %%%%%%%%%%%%%
% str1=pwd;                     % Get current working path
% str2='Programs';                % remove programs    
% Fold_Path=strrep(str1,str2,''); % Path of the main folder
% % Fold_Path=str1;
% %path_sep=Fold_Path(length(Fold_Path));
% path_sep='\';
% % default=load([Fold_Path path_sep 'Settings' path_sep 'Default_Settings.mat']);
% % [Folders]=Create_Folder_Struct(Fold_Path,default,path_sep);
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% plot emblam
% % % % % % % % % % 
%%%%%%%%%%% Load  settings and store
 V1=varargin;
 Folders=V1{1};
 handles.Folders=Folders;
 

% ResF=[Folders.Main_Fold.Sim Folders.Sub_Fold.Res_Fold '\' Folders.Sub_Fold.Mec_Fold '\'];
[DL,F1,index]=Sim_Folders(handles);
set(handles.Mech_Sim_Folders,'String',F1,'Value',index);
F2=Folders.Sub_Fold.Res_Fold;
Res_H=[Folders.Main_Fold.Sim Folders.Sub_Fold.Res_Fold '\'];
set(handles.Expect_Sim_Fold,'String',F2);
if isdir(Res_H)==0
    F2=char(F1);
end
set(handles.Cur_Sim_Fold,'String',F2);
handles.CSimF=F2;
guidata(hObject, handles);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [SimF]=create_local_folders(handles);
[SimF,MECH_MAT]=create_local_folders(handles);
handles.SimF=SimF;
handles.MECH_MAT=MECH_MAT;
handles.vprop=1;
guidata(hObject, handles);
plot_material_data(handles);

% Sim_Heat=[Folders.Main_Fold.Sim  F2 '\'];
% Sim_Mech=[Sim_Heat Folders.Sub_Fold.Mec_Fold '\'];
% Mech_Prog=Folders.Main_Fold.Mec;
% SimF=struct('Sim_Heat',Sim_Heat,'Sim_Mech',Sim_Mech,'Mech_Prog',Folders.Main_Fold.Mec);
% handles.SimF=SimF;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Independent running
addpath(SimF.Mech_Prog);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PATH_Material=[SimF.Sim_Heat '\Material.mat'];
% load(PATH_Material);
% MECH_MAT=MECH_Material_Prop2(TS,TL);
% handles.MECH_MAT=MECH_MAT;
% handles.vprop=1;
% if isdir([Sim_Mech])==0
%     mkdir(Sim_Mech);
% end
% Make_Mech_Material_File(MECH_MAT,Sim_Mech);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot_material_data(handles)





% UIWAIT makes Mech_1 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Mech_1_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;




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


% --- Executes on button press in Run_Stress_Sim.
function Run_Stress_Sim_Callback(hObject, eventdata, handles)
% hObject    handle to Run_Stress_Sim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Run_Stress_Sim
SV=get(hObject,'Value');
% set(hObject, 'String','STOP Simulation');
% guidata(hObject, handles);
%     
if SV==1
    set(hObject, 'String','STOP Simulation');
    drawnow
    Run_Mech_SIM(handles)
else
    set(hObject, 'String','RUN Simulation');
    drawnow
end
% %%%%%%%% BEAMBC Making %%%%%%%%%%%%%%%%%%%%%
% a(1,1)=get(handles.Axial_Force,'Value');
% a(1,2)=get(handles.Axial_Strain,'Value');
% a(1,3)=str2double(get(handles.Axial,'String'));
% 
% a(2,1)=get(handles.Moment_X,'Value');
% a(2,2)=get(handles.Curvature_X,'Value');
% a(2,3)=str2double(get(handles.Bend_X,'String'));
% 
% a(3,1)=get(handles.Moment_Y,'Value');
% a(3,2)=get(handles.Curvature_Y,'Value');
% a(3,3)=str2double(get(handles.Bend_Y,'String'));
% BEAMBC=zeros(3,2);
% for i=1:3
%     if a(i,2)==1
%         BEAMBC(i,1)=1;
%     end
%     BEAMBC(i,2)=a(i,3);
% end
% %%%%%%%% BEAMBC Making %%%%%%%%%%%%%%%%%%%%%
% 
% %%%%%%% Material Model %%%%%%%%%%%%%%%%%%%%%
% b(1)=get(handles.E,'Value');
% b(2)=get(handles.E_PP,'Value');
% b(3)=get(handles.E_PH,'Value');
% b(4)=get(handles.E_VP,'Value');
% MAT_MODEL=find(b);
% %%%%%%% Folders %%%%%%%%%%%%%%%%%%%%%
% SimF=handles.SimF;
% d=dir([SimF.Sim_Heat '\DB*.mat']);
% str = {d.name};
% ndir=size(str,2);
% if ndir<=2
%     errordlg('Stress Simulation - NOT POSSIBLE','NOT POSSIBLE');
% else
%     solve_mech_2(BEAMBC,MAT_MODEL,SimF);
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Axial_Force_Callback(hObject, eventdata, handles)
% hObject    handle to Axial_Force (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Axial_Force as text
%        str2double(get(hObject,'String')) returns contents of Axial_Force as a double


% --- Executes during object creation, after setting all properties.
function Axial_Force_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Axial_Force (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Axial_Callback(hObject, eventdata, handles)
% hObject    handle to Axial (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Axial as text
%        str2double(get(hObject,'String')) returns contents of Axial as a double


% --- Executes during object creation, after setting all properties.
function Axial_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Axial (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Moment_x_Callback(hObject, eventdata, handles)
% hObject    handle to Moment_x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Moment_x as text
%        str2double(get(hObject,'String')) returns contents of Moment_x as a double


% --- Executes during object creation, after setting all properties.
function Moment_x_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Moment_x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Bend_X_Callback(hObject, eventdata, handles)
% hObject    handle to Bend_X (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Bend_X as text
%        str2double(get(hObject,'String')) returns contents of Bend_X as a double


% --- Executes during object creation, after setting all properties.
function Bend_X_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Bend_X (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Moment_y_Callback(hObject, eventdata, handles)
% hObject    handle to curvature_y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of curvature_y as text
%        str2double(get(hObject,'String')) returns contents of curvature_y as a double


% --- Executes during object creation, after setting all properties.
function Moment_y_CreateFcn(hObject, eventdata, handles)
% hObject    handle to curvature_y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Bend_Y_Callback(hObject, eventdata, handles)
% hObject    handle to Bend_Y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Bend_Y as text
%        str2double(get(hObject,'String')) returns contents of Bend_Y as a double


% --- Executes during object creation, after setting all properties.
function Bend_Y_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Bend_Y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes when selected object is changed in Material_Model.
function Material_Model_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in Material_Model 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
% A=get(hObject)
% set(handles.A,'Value',1);


% --- Executes on button press in Mech_Result.
function Mech_Result_Callback(hObject, eventdata, handles)
% hObject    handle to Mech_Result (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Mech_Result
% SimF=handles.SimF;
Folders=handles.Folders;
V=Folders;
% V{1}=SimF.Sim_Heat;
% V{2}=SimF.Sim_Mech;
% V{3}=SimF.Mech_Prog;
Mech_Results_1(V);

function[]=plot_material_data(handles)
cla(handles.mech_mat)
axes(handles.mech_mat)
vprop=handles.vprop;
if vprop==1
    TData=handles.MECH_MAT.EyT;
elseif vprop==2
    TData=handles.MECH_MAT.NuT;
elseif vprop==3
    TData=handles.MECH_MAT.CTET;
% elseif Prop_Sel==4
%     TData=handles.Properties.PCF;
% else
%     TData=handles.Properties.ENTH;
end
switch vprop
    case 1
        %%%%% Thermal conductivity [KHT1]
%         KHT1=get(handles.KHT1,'Data')
%         KHT1=str2double(get(handles.KHT1,'Data'));
        plot(TData(:,1),TData(:,2)./1e9,'r','LineWidth',3);
        ylabel('Elasticity Modulus [GPa]','FontSize',12);
        xlabel('Temperature [degC]','FontSize',12);
    case 2
        %%%%% Specific Heat Capacity [SHT1]
%         SHT1=str2double(get(handles.SHT1,'Data'));
       plot(TData(:,1),TData(:,2),'b','LineWidth',3);
        ylabel('Poisons Ratio [-]','FontSize',12);
        xlabel('Temperature [degC]','FontSize',12);
    case 3
        %%%%% Density [RHT1]
%         RHT1=str2double(get(handles.RHT1,'Data'));
      plot(TData(:,1),TData(:,2),'k','LineWidth',3);
       ylabel('Thermal Strain [-]','FontSize',12);
        xlabel('Temperature [degC]','FontSize',12);
      
%     case 4
%               plot(TData(:,1),TData(:,2),'m','LineWidth',3);
%         ylabel('Liquid fraction [-]','FontSize',12);
%         xlabel('Temperature [degC]','FontSize',12);
%     case 5
%         plot(TData(:,1),TData(:,2)./1e3,'m','LineWidth',3);
%         ylabel('Enthalpy [J/g]','FontSize',12);
%         xlabel('Temperature [degC]','FontSize',12);        
end
set(gca,'FontSize',12,'LineWidth',2.0,'Box','off');



function[Folders]=Create_Folder_Struct(Fold_Path,default,path_sep)
%         Fold_Path=[Fold_Path path_sep];
        Main_Fold=struct('Sim',{[Fold_Path default.Sim_Main_Fold path_sep]},...
                 'Mat',{[Fold_Path default.Mat_Fold_Name path_sep]},...
                 'Grd',{[Fold_Path default.Grd_Fold_Name path_sep]},...
                 'Set',{[Fold_Path default.Set_Fold_Name path_sep]},...
                 'Pro',{[Fold_Path default.Pro_Fold_Name path_sep]},...
                  'Mec',{[Fold_Path default.Pro_Fold_Name path_sep 'MECH' path_sep]});
 File=struct('Grd',{default.Grd_File_Name},...
                 'SDef',{default.Def_File_Name},...
                 'SUse',{default.Def_File_Name});
 Sub_Fold=struct('Sim_Fold',{'1015_S235JR_2'},...
                 'Res_Fold',{'1015_S235JR_2'},...
                 'Mat_Fold',{''},...
                 'Mec_Fold',{'MECH'});
             
  Folders=struct('Main_Fold',{Main_Fold},...
                 'Sub_Fold',{Sub_Fold},...
                 'File',{File});


% --- Executes on selection change in Mech_Sim_Folders.
function Mech_Sim_Folders_Callback(hObject, eventdata, handles)
% hObject    handle to Mech_Sim_Folders (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Mech_Sim_Folders contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Mech_Sim_Folders
contents = cellstr(get(hObject,'String'));
CSimF= contents{get(hObject,'Value')};
% Cur_Res_Folder=[handles.Folders.Main_Fold.Sim Res_Fold];
% handles.Folders.Sub_Fold.Sim_Fold=Sim_Fold;
handles.CSimF=CSimF;
guidata(hObject, handles);
[SimF,MECH_MAT]=create_local_folders(handles);
handles.SimF=SimF;
handles.MECH_MAT=MECH_MAT;
guidata(hObject, handles);
plot_material_data(handles);

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
    d2=dir([Folders.Main_Fold.Sim folders{i}]);
    str={d2.name};
    str=str(3:length(str));
    if length(str)>1 %&& ct==0
        index=i;
        F2(ct+1)=folders(i);
        ct=ct+1;
    end
    DL(i)=length(str);
end



function Curr_Sim_Fold_Callback(hObject, eventdata, handles)
% hObject    handle to Expect_Sim_Fold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Expect_Sim_Fold as text
%        str2double(get(hObject,'String')) returns contents of Expect_Sim_Fold as a double


% --- Executes during object creation, after setting all properties.
function Expect_Sim_Fold_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Expect_Sim_Fold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function[SimF,MECH_MAT]=create_local_folders(handles)
Folders=handles.Folders;
CSimF=handles.CSimF;
Sim_Heat=[Folders.Main_Fold.Sim  CSimF '\'];
Sim_Mech=[Sim_Heat Folders.Sub_Fold.Mec_Fold '\'];
Mech_Prog=Folders.Main_Fold.Mec;
SimF=struct('Sim_Heat',Sim_Heat,'Sim_Mech',Sim_Mech,'Mech_Prog',Folders.Main_Fold.Mec);
set(handles.Cur_Sim_Fold,'String',CSimF);
SimF

PATH_Material=[SimF.Sim_Heat '\Material.mat'];
load(PATH_Material);
TS
TL
Mech_Prog
rehash
addpath(Mech_Prog);
dir
MECH_MAT=MECH_Material_Prop2(TS,TL);
% handles.MECH_MAT=MECH_MAT;
% handles.vprop=1;
if isdir([Sim_Mech])==0
    mkdir(Sim_Mech);
    Make_Mech_Material_File(MECH_MAT,Sim_Mech);
end



% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over Mech_Sim_Folders.
function Mech_Sim_Folders_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to Mech_Sim_Folders (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on key press with focus on Run_Stress_Sim and none of its controls.
function Run_Stress_Sim_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to Run_Stress_Sim (see GCBO)
% eventdata  structure with the following fields (see UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)

function[]=Run_Mech_SIM(handles)
%%%%%%%% BEAMBC Making %%%%%%%%%%%%%%%%%%%%%
a(1,1)=get(handles.Axial_Force,'Value');
a(1,2)=get(handles.Axial_Strain,'Value');
a(1,3)=str2double(get(handles.Axial,'String'));

a(2,1)=get(handles.Moment_X,'Value');
a(2,2)=get(handles.Curvature_X,'Value');
a(2,3)=str2double(get(handles.Bend_X,'String'));

a(3,1)=get(handles.Moment_Y,'Value');
a(3,2)=get(handles.Curvature_Y,'Value');
a(3,3)=str2double(get(handles.Bend_Y,'String'));
BEAMBC=zeros(3,2);
for i=1:3
    if a(i,2)==1
        BEAMBC(i,1)=1;
    end
    BEAMBC(i,2)=a(i,3);
end
%%%%%%%% BEAMBC Making %%%%%%%%%%%%%%%%%%%%%

%%%%%%% Material Model %%%%%%%%%%%%%%%%%%%%%
b(1)=get(handles.E,'Value');
b(2)=get(handles.E_PP,'Value');
b(3)=get(handles.E_PH,'Value');
b(4)=get(handles.E_VP,'Value');
MAT_MODEL=find(b);
%%%%%%% Folders %%%%%%%%%%%%%%%%%%%%%
SimF=handles.SimF;
d=dir([SimF.Sim_Heat '\DB*.mat']);
str = {d.name};
ndir=size(str,2);
if ndir<=2
    errordlg('Stress Simulation - NOT POSSIBLE','NOT POSSIBLE');
else
    solve_mech_2(BEAMBC,MAT_MODEL,SimF);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
