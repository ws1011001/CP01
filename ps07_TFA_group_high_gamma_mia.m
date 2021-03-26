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
OPTIONS.smooth = 9;  % this number x 3 < number of ROIs (for successful smoothing)
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
mni_varnames={'subjects','conditions','AAL3','mni_x','mni_y','mni_z','contacts_1','contacts_2'};
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
fgrp=fullfile(gdir,sprintf('group_%s_%s_data-ROI_%d_%d_%d.mat',mtg,tfa,freq_step,freq_lobd,freq_upbd));
save(fgrp,'datatable');  % it's a very big file. Don't use it if have many patients
%% ---------------------------

%% between-condition comparisons for each ROI
conditions_pairs={{'AAp','AAt'},{'AVp','AVt'},{'VAp','VAt'},{'VVp','VVt'}};
for i=1:n
  subj=subjects{i};
  sdir=fullfile(bdir,subj);
  fprintf('Perform permutation tests between conditions for each ROI for subject %s ......\n',mtg,tfa,subj);
  % read ROI-based signals of different conditions
  for icond=1:ncond
    cdir=fullfile(sdir,conditions{icond});
    fdat=fullfile(cdir,sprintf('%s_%s_%s_data_%d_%d_%d.mat',conditions{icond},mtg,tfa,freq_step,freq_lobd,freq_upbd)); 
    load(fdat, 'roi');
    signals.(conditions{icond})=roi.signals;
  end
  % compare conditions
  for ipair=1:length(conditions_pairs)
    cond_pair=conditions_pairs{ipair};
    rdir=fullfile(sdir,sprintf('%s_ROIs_%s-%s',subj,cond_pair{1},cond_pair{2}));
    if ~exist(rdir,'dir'); mkdir(rdir); end
    fperm=fullfile(rdir,sprintf('ROIs_%s-%s_%s_%s_%d_%d_%d.mat',cond_pair{1},cond_pair{2},mtg,tfa,freq_step,freq_lobd,freq_upbd));
    if ~exist(fperm,'file')
      roi_perm=roi_stats_permutations(signals.(cond_pair{1}),signals.(cond_pair{2}),OPTIONS);
      save(fperm,'roi_perm');
    else
      load(fperm,'roi_perm');
    end
    roi_plot_conditions(roi_perm,cond_pair,OPTIONS,rdir); 
  end 
end
%% ---------------------------



