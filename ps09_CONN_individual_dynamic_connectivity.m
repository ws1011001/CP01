%% ---------------------------
%% [script name] ps09_CONN_individual_dynamic_connectivity.m
%%
%% SCRIPT to ...
%%
%% By Shuai Wang, [date] 2021-03-31
%%
%% ---------------------------
%% Notes:
%%   
%%
%% ---------------------------

%% clean up
close all
clear
clc
%% ---------------------------

%% set environment (packages, functions, working path etc.)
% setup working path
mdir='/media/wang/BON/Projects/CP01';
ddir=fullfile(mdir,'SEEG_LectureVWFA','derivatives');  % derivatives in the BIDS structure
bdir=fullfile(ddir,'mia_SEEG_LectureVWFA');            % path to brainsotrm database
% read the subjects list
fid=fopen(fullfile(mdir,'CP01_subjects.txt'));
subjects=textscan(fid,'%s');
fclose(fid);
subjects=subjects{1};  % the subjects list
n=length(subjects);
% set processing parameters
mtg='bipolar';                          % montage
tfa='morlet';                           % time-frequency analysis
freqb = [[10, 40, 80]; [10, 80, 120]];  % frequency band of interests (FBI)
nfreqb=size(freqb,1);                   % number of FBI
time = -0.4:0.001:0.8;                  % from -400ms to 800ms
win_noedges = [-0.3 0.7];               % time-window to remove wavelet edges
conditions={'Ap','Vp'};                 % conditions
ncond=length(conditions);               % number of conditions
data_mode={'ERP','HGA'};                % data modality : ERP or high-gamma activity
nmode=length(data_mode);                % number of data modalities
max_r=15;                               % max radius to combine contacts
min_w=0.3;                              % min ICC to combine contacts
% set switches
isGetROIs=false;
isDynConn=true;
%% ---------------------------

%% extract signals of ROIs
if isGetROIs
  % extract signals for each subject
  for i=1:n
    subj=subjects{i};
    sdir=fullfile(bdir,subj);
    fprintf('Extract %s %s signals for each ROI for subject %s ......\n',mtg,tfa,subj);
    % read the table of MNI coordinates of monopolar channels
    fmni=fullfile(mdir,'SEEG_LectureVWFA',subj,'anat',sprintf('%s_space-MNI_contacts-monopolar.csv',subj));
    if ~exist(fmni,'file')
      fcoi=fullfile(mdir,'SEEG_LectureVWFA',subj,'anat/elec2atlas.mat');
      mono_mni=mni_coi2csv(fcoi,fmni);
    else
      mono_mni=readtable(fmni);
    end      
    % prepare the RS baseline data to calculate ICC
    frsb=fullfile(sdir,'RS','RS_bipolar_LFP_data.mat');
    data_bl=load(frsb,'zs'); data_bl.F=data_bl.zs;
    % calculate the signals for each ROI
    for imode=1:nmode
      dmode=data_mode{imode};
      for icond=1:ncond
        con=conditions{icond};
        ctoken=sprintf('%s_%s',dmode,con);
        cdir=fullfile(sdir,ctoken);
        % collect data file
        switch dmode
          case 'ERP'
            fprintf('Process %s data for condition %s ......\n',dmode,con);
            fdat=fullfile(cdir,sprintf('%s_%s_LFP_data.mat',ctoken,mtg));
            data=load(fdat,'labels','zs','Time'); data.F=data.zs;
            roi_all_in_one(mono_mni, data, data_bl, fdat, max_r, min_w);
            working_files.(dmode).(con)=fdat;
          case 'HGA'
            for ifreqb=1:nfreqb
              fprintf('Process %s data for condition %s ......\n',dmode,con);
              fdat=fullfile(cdir,sprintf('%s_%s_%s_data_%d_%d_%d.mat',ctoken,mtg,tfa,freqb(ifreqb,:)));
              data=load(fdat,'labels','zs','Time'); data.F=data.zs;              
              roi_all_in_one(mono_mni, data, data_bl, fdat, max_r, min_w);
              dtoken=sprintf('%s%d',dmode,ifreqb);
              working_files.(dtoken).(con)=fdat;
            end            
          otherwise
            error('Invalid data mode. QUIT!')
        end
      end
    end
    % output the information of ROIs
    load(fdat,'roi_data');
    roi_data(:,end)=[];
    roi_table=cell2table(roi_data,'VariableNames',{'ROI','mni_x','mni_y','mni_z','radius','longside','shortside','number_of_contacts','ICC'});
    ftab=fullfile(sdir,sprintf('%s_ROI-signals_information.csv',subj));
    writetable(roi_table,ftab);
    % output the working files
    fworking=fullfile(sdir,sprintf('%s_ROI-signals_working-files.mat',subj));
    save(fworking,'working_files');
  end
end
%% ---------------------------

%% estimate dynamic connectivity
if isDynConn
  for i=1:n
    subj=subjects{i};
    sdir=fullfile(bdir,subj);
    fprintf('Estimate dynamic connectivity for subject %s ......\n',subj);
    % read ROIs information (manually inspected)
    ftab=fullfile(sdir,sprintf('%s_ROI-signals_information.csv',subj));
    roi_table=readtable(ftab);
    % select ROIs
    roi_idx=logical(roi_table.inclusion_manual);  % index of ROIs that to be included
    roi_networks=roi_table.subnetworks(roi_idx);  % ROI labels of sub-networks
    roi_vwfa=ismember(roi_networks,{'VWFA'});     % ROIs to be merged as VWFA
    % do dynamic connectivity
    load(fullfile(sdir,sprintf('%s_ROI-signals_working-files.mat',subj)));  % read working files
    dtokens=fieldnames(working_files);
    for j=1:length(dtokens)
      dtoken=dtokens{j};
      for icond=1:ncond
        con=conditions{icond};
        fdat=working_files.(dtoken).(con);
        % read ROIs data
        load(fdat,'roi_data');
        signals=roi_data(roi_idx,end);
        signals=permute(cat(3, signals{:}), [3 1 2]);  % convert cell array to matrix
        % combine clusters (as the ROI of VWFA)
        signals_vwfa=mean(signals(roi_vwfa,:,:),1);  % average the VWFA signals - do not squeeze
        signals(roi_vwfa,:,:)=[];                    % remove the original VWFA signals
        signals=[signals;signals_vwfa];              % put the averaged VWFA signal at the bottom
      end
    end
  end
end
%% ---------------------------