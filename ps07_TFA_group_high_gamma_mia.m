%% ---------------------------
%% [script name] ps07_TFA_group_high_gamma_mia.m
%%
%% SCRIPT to do statistical comparisons on high-gamma activity by grouping patients.
%%
%% By Shuai Wang, [date] 2021-03-25
%%
%% ---------------------------
%% Notes: - This script is modified from ASD's script seeg_pipeline_MAIN.m (see below)
%%   
%%
%% ---------------------------
% ========================================================================
% This file is part of MIA.
% 
% MIA is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% MIA is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%  
% Copyright (C) 2016-2018 CNRS - Universite Aix-Marseille
%
% ========================================================================
% This software was developed by
%       Anne-Sophie Dubarry (CNRS Universite Aix-Marseille)
% ***********************************************************************
%
% Script that process group analysis for GAMMABLOCK : 
% 1) Get roi from group analysis of all conditions
% 2) For the found roi, perform ttest between conditions 
% ========================================================================


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
gdir=fullfile(bdir,'group');                           % group folder
% read the subjects list
fid=fopen(fullfile(mdir,'CP01_subjects.txt'));
subjects=textscan(fid,'%s');
fclose(fid);
subjects=subjects{1};  % the subjects list
n=length(subjects);
% set processing parameters
mtg='bipolar';  % montage
tfa='morlet';   % time-frequency analysis
freq_step=10;
freq_lobd=80;
freq_upbd=120;
OPTIONS.nperm = 5000;
OPTIONS.smooth = 5;  % this number x 3 < number of ROIs (for successful smoothing)
OPTIONS.threshp = 0.05;
OPTIONS.time = -0.2:0.001:0.8;  % from -200ms to 800ms
OPTIONS.win_noedges = [-0.15 0.75];
% baseline=[-0.2,-0.01];  % baseline from -200ms to -10ms
% timewindow=[-0.2,0.8];  % time-window of the epoch i.e. -200 ms to 800 ms
% nperm=10000;            % the number of randomizations in permutation tests
% fdr_duration=20;        % duration (in ms) for FDR correction in the time domain
conditions={'AAp';'AAt';'AVp';'AVt';'VAp';'VAt';'VVp';'VVt'};
ncond=length(conditions);
%% ---------------------------

%% group the contacts by ROIs
% initialize data tables
mni_varnames={'subjects','conditions','AAL3','mni_x','mni_y','mni_z','contacts_2','contacts_1'};
mni_table=array2table(zeros(0,length(mni_varnames)),'VariableNames',mni_varnames);
% extract ROI info
for icond=1:ncond
  con=conditions{icond};
  for i=1:n
    subj=subjects{i};
    sdir=fullfile(bdir,subj);
    fprintf('Combine ROI-based signals of condition %s for subject %s ......\n',con,subj);
    % load up data
    cdir=fullfile(sdir,con);
    fdat=fullfile(cdir,sprintf('%s_%s_%s_data_%d_%d_%d.mat',con,mtg,tfa,freq_step,freq_lobd,freq_upbd));
    load(fdat,'roi');
    % group individual coordinates of contacts
    n_contacts=length(roi.bipolar_contacts);
    mni_cells=[repmat({subj},n_contacts,1),repmat({con},n_contacts,1),...
               roi.AAL3,num2cell(roi.bipolar_mni),roi.bipolar_contacts];
    mni_table=[mni_table;cell2table(mni_cells,'VariableNames',mni_varnames)];
  end
end
% save group data
fmni=fullfile(gdir,sprintf('group_space-MNI_contacts-%s_ROI-AAL3',mtg));
writetable(mni_table,[fmni '.csv']);
save([fmni '.mat'],'mni_table');
%% ---------------------------

%% group signals by ROIs
rois_full=unique(mni_table.AAL3);  % get the full list of ROIs of all patients
nrois=length(rois_full);
rois_info=cell(nrois,5);
for iroi=1:nrois
  roi_name=rois_full{iroi};
  zscores_avg=cell(n,ncond);
  zscores_avg_avg=cell(n,ncond);
  zscores_contacts_avg=cell(n,ncond);
  tvals_contacts=cell(n,ncond);
  pvals_contacts=cell(n,ncond);
  threshdur_contacts=zeros(n,ncond);
  contacts_labels=cell(n,ncond);
  intersubj_corr=zeros(1,ncond);
  interchan_corr=zeros(1,ncond);
  for icond=1:ncond
    con=conditions{icond};
    for i=1:n
      subj=subjects{i};
      sdir=fullfile(bdir,subj);
      fprintf('Combine ROI-based signals of condition %s for subject %s ......\n',con,subj);   
      % load up data
      cdir=fullfile(sdir,con);
      fdat=fullfile(cdir,sprintf('%s_%s_%s_data_%d_%d_%d.mat',con,mtg,tfa,freq_step,freq_lobd,freq_upbd));
      fsta=strrep(fdat,'_data','_stats');
      load(fdat);
      % extract ROI signals
      if isfield(roi.signals,roi_name)
        % within-subject between-contacts averaged z-scores (for permutation tests)
        zscores_avg{i,icond}=roi.signals.(roi_name).avg;
        % within-subject between-contacts across-trials averaged z-scores (for calculating ISC)
        zscores_avg_avg{i,icond}=mean(roi.signals.(roi_name).avg,2);
        % contact-level across-trials averaged z-scores
        zscores_contacts_avg{i,icond}=mean(roi.signals.(roi_name).F,3);
        % contact-level T values and P values (for showing significant contacts for each ROI)
        load(fsta);
        contacts_def_idx=roi_match_bipolar_labels(labels',roi.signals.(roi_name).channels);
        tvals_contacts{i,icond}=tvals(contacts_def_idx,:);
        pvals_contacts{i,icond}=pvals(contacts_def_idx,:);
        threshdur_contacts(i,icond)=stats.threshdur;
        % subject-contact labels
        contacts_labels{i,icond}=[repmat({subj},roi.signals.(roi_name).dims(1),1) roi.signals.(roi_name).channels];
      end
    end
    % number of subjects
    subjs_idx=any(~cellfun(@isempty,contacts_labels),2);   
    nsubjs=sum(subjs_idx);
    % calculate inter-subject correlations
    if nsubjs >= 2
      isc_zscores=zscores_avg_avg(:,icond);
      isc_zscores=cat(2,isc_zscores{:});
      intersubj_corr(icond)=icc_kendall_w(isc_zscores,0);
    end
    % calculate inter-channel (across-subjects) correlations
    icc_zscores=zscores_contacts_avg(:,icond);
    icc_zscores=cat(1,icc_zscores{:});
    interchan_corr(icond)=icc_kendall_w(icc_zscores',0);
  end
  % save the information of ROIs
  rois_info{iroi,1}=roi_name;
  rois_info{iroi,2}=nsubjs;
  rois_info{iroi,3}=subjs_idx;
  rois_info{iroi,4}=mean(intersubj_corr);
  rois_info{iroi,5}=mean(interchan_corr);
  finfo=fullfile(gdir,'group_Info_ROI-AAL3.mat');
  save(finfo,'rois_info');
  % save group data for each ROI
  froi=fullfile(gdir,sprintf('group_N%d_ROI-AAL3-%s_%s_%s_data_%d_%d_%d.mat',nsubjs,roi_name,mtg,tfa,freq_step,freq_lobd,freq_upbd));
  save(froi,'roi_name','subjects','conditions','nsubjs','contacts_labels','intersubj_corr',...
            'zscores_avg','zscores_avg_avg','zscores_contacts_avg','interchan_corr',...
            'tvals_contacts','pvals_contacts','threshdur_contacts');
end
%% ---------------------------

%% between-condition comparisons for each ROI at the group level
conditions_pairs={{'AAp','AAt'},{'AVp','AVt'},{'VAp','VAt'},{'VVp','VVt'}};
% extract ROI signals for each condition
for iroi=1:nrois
  roi_name=rois_info{iroi,1};
  roi_nsub=rois_info{iroi,2};
  if roi_nsub >= 2 
    % load up ROI data
    froi=fullfile(gdir,sprintf('group_N%d_ROI-AAL3-%s_%s_%s_data_%d_%d_%d.mat',roi_nsub,roi_name,mtg,tfa,freq_step,freq_lobd,freq_upbd));
    roi_data=load(froi);
    signals_zs=roi_data.signals_zs;
    % output ROI summary
    odir=fullfile(gdir,sprintf('group_N%d_ROI-AAL3-%s_summary',roi_nsub,roi_name));
%     if ~exist(odir,'dir'); mkdir(odir); end
%     roi_plot_group_contacts(roi_data,OPTIONS,odir);
%     % read ROI-based signals of different conditions
%     fprintf('Extract signals of ROI %s from %d subjects to group data ......\n',roi_name,roi_nsub);
%     for icond=1:ncond
%       zs_avg=roi_data.zscores_avg(:,icond);
%       zs_avg=cat(2,zs_avg{:});
%       signals_zs.(conditions{icond}).(roi_name).avg=zs_avg;
%     end 
%     save(froi,'signals_zs','-append');
    % compare conditions
    for ipair=1:length(conditions_pairs)
      cond_pair=conditions_pairs{ipair};
      rdir=fullfile(odir,sprintf('group_ROI-AAL3-%s_%s-%s',roi_name,cond_pair{1},cond_pair{2}));
      if ~exist(rdir,'dir'); mkdir(rdir); end
      fperm=fullfile(rdir,sprintf('group_ROI-AAL3-%s_%s-%s_%s_%s_%d_%d_%d.mat',roi_name,cond_pair{1},cond_pair{2},mtg,tfa,freq_step,freq_lobd,freq_upbd));
      if ~exist(fperm,'file')
        roi_perm=roi_stats_permutations(signals_zs.(cond_pair{1}),signals_zs.(cond_pair{2}),OPTIONS);
        save(fperm,'roi_perm');
      else
        load(fperm,'roi_perm');
      end
      roi_plot_conditions(roi_perm,cond_pair,OPTIONS,rdir); 
    end
    % clean up working variables
    clear signals_zs roi_perm
  end
end
%% ---------------------------



