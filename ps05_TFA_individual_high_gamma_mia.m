%% ---------------------------
%% [script name] ps05_TFA_individual_high_gamma_mia.m
%%
%% SCRIPT to do statistical comparisons on high-gamma activity at the individual level.
%%
%% By Shuai Wang, [date] 2021-03-19
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
conditions={'AAp','AAt','AVp','AVt','VAp','VAt','VVp','VVt'};
ncond=length(conditions);
% set switches
isGetROIs=true;
isCompCon=false;
%% ---------------------------

%% extract signals for each ROI
if isGetROIs
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
    mono_mni_xyz=[mono_mni.x mono_mni.y mono_mni.z];
%     % calculate the signals for each ROI
%     for icond=1:ncond
%       cdir=fullfile(sdir,conditions{icond});
%       fdat=fullfile(cdir,sprintf('%s_%s_%s_data_%d_%d_%d.mat',conditions{icond},mtg,tfa,freq_step,freq_lobd,freq_upbd));
%       data=load(fdat,'labels','zs','Time');
%       % get bipolar labels and coordinates
%       [bipo_mni, bipo_labels, bipo_flags]=mni_mono2bipolar(mono_mni_xyz, mono_mni.contacts, data.labels);
%       roi.bipolar_mni=bipo_mni;
%       roi.bipolar_contacts=bipo_labels;
%       roi.bipolar_flags=bipo_flags;
%       % get regions of interest
%       roi.AAL3=mni_find_regions(bipo_mni,'AAL3');
%       % get ROI-specific signals
%       roi.signals=roi_average_channels(data.zs(bipo_flags, :, :), bipo_mni, bipo_labels, roi.AAL3);
%       roi.time=data.Time;
%       % save ROI data
%       save(fdat,'roi','-append');
%     end
  end
end
%% ---------------------------

%% between-condition comparisons for each ROI
if isCompCon
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
end
%% ---------------------------



