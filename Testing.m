function[]=Testing()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Folder's Name
str1=pwd;                     % Get current working path
str2='Programs';                % remove programs    
Fold_Path=strrep(str1,str2,''); % Path of the main folder
path_sep='\';
default=load([Fold_Path path_sep 'Settings' path_sep 'Default_Settings.mat']);
[Folders]=Create_Folder_Struct(Fold_Path,default,path_sep);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d=dir(Folders.Main_Fold.Sim);
str={d.name};
if length(str)>2
    folders=str(3:length(str));
end
index=1;ct=0;
for i=1:length(folders)
    str=[];
    d2=dir([Folders.Main_Fold.Sim folders{i} path_sep Folders.Sub_Fold.Mec_Fold '\DB_M_*.mat']);
    str={d2.name};
    str=str(3:length(str));
    if length(str)>1 %&& ct==0
        index=i;
        F2(ct+1)=folders(i);
        ct=ct+1;
    end
    DL(i)=length(str);
end
DL
 F2
%  d2.name
% d2=dir([Sim_Main_Fold handles.Folders.Sub_Fold.Sim_Fold]);
% if length({d2.name})>=2
%     dset=find(strcmp(folders,handles.Folders.Sub_Fold.Sim_Fold));
%     handles.Folders.Sub_Fold.Res_Fold=handles.Folders.Sub_Fold.Sim_Fold;
% else
%     dset=index;
%     handles.Folders.Sub_Fold.Res_Fold=folders{index};
% end
% set(handles.Results_Folders,'String',folders,'Value',dset);


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
