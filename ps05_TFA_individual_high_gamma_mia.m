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
baseline=[-0.2,-0.01];  % baseline from -200ms to -10ms
timewindow=[-0.2,0.8];  % time-window of the epoch i.e. -200 ms to 800 ms
nperm=10000;            % the number of randomizations in permutation tests
fdr_duration=20;        % duration (in ms) for FDR correction in the time domain
conditions={'AAp','AAt','AVp','AVt','VAp','VAt','VVp','VVt'};
ncond=length(conditions);
%% ---------------------------

%% extract signals for each ROI
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
  % calculate the signals for each ROI
  for icond=1:ncond
    cdir=fullfile(sdir,conditions{icond});
    fdat=fullfile(cdir,sprintf('%s_%s_%s_data_%d_%d_%d.mat',conditions{icond},mtg,tfa,freq_step,freq_lobd,freq_upbd));
    data=load(fdat,'labels','zs');
    % get bipolar labels and coordinates
    [bipo_mni, bipo_labels, bipo_flags]=mni_mono2bipolar(mono_mni_xyz, mono_mni.contacts, data.labels);
    roi.bipolar_mni=bipo_mni;
    roi.bipolar_contacts=bipo_labels;
    roi.bipolar_flags=bipo_flags;
    % get regions of interest
    roi.AAL3=mni_find_regions(bipo_mni,'AAL3');
    % get ROI-specific signals
    roi.signals=roi_average_channels(data.zs(bipo_flags, :, :), bipo_mni, bipo_labels, roi.AAL3);
    % save ROI data
    save(fdat,'roi','-append');
  end
end


grpOPTIONS.mtg = 'bipolar';
grpOPTIONS.win_noedges = [-401,1601] ; % time window to analyse (in samples)
grpOPTIONS.freqid = {1,2}; % for now, we replace the freq and method by the index
% maindir : MIA working directory, ALL CONDITIONS
grpOPTIONS.maindir = '/Volumes/ASD_400GO/GAMMABLOCK_contrast/Analysis_MarsPower/Results_BLOC';

% Loads the localizations of all patients contacts 
load('/Volumes/ASD_400GO/GAMMABLOCK_contrast/Analysis_MarsPower/commons/m_table_as.mat');

% get a group strcuture (lines : patients ; colomn : studies)
[ganalysis] = s5_group_data_v1(sFiles, grpOPTIONS); 
%% ---------------------------

%% GET TABLE  OF EFFECTS 
%-----------------------------------------------------------------------------
% Filters out regions not interessting (out  blanc, etc.0)
isGood =~(strcmp(m_table_as(:,5),'out')|strcmp(m_table_as(:,5),'lesion'));
% isGood =boolean(ones(size(m_table_as,1),1));
[m_table_effect, s, smask, all_labels] = get_table_effect(m_table_as(isGood,:), ganalysis);
%% ---------------------------


%% GET ROIS
%-----------------------------------------------------------------------------
%Define OPTIONS for group anlaysis and display
getOPTIONS.freq = 2; % 1= freqid 2
getOPTIONS.nPt= 2 ;% min numb of pt by roi
getOPTIONS.signifmode =1 ; % mode to select significant activity (if 0 no constrain)
getOPTIONS.signmode = 'signed';

gan = [ganalysis{1,:}] ;
 
% Get all rois that has significant activity for this frequency band
rois = get_roi(m_table_effect,ganalysis{1,1}.t, s, smask, all_labels,{gan.freqb},getOPTIONS);

%% FILTER OUT SOME ROIS (mincorr)
%-----------------------------------------------------------------------------
fOPTIONS.mincorr = 0.3;                 % minimum correlation
% Removes blanc
froi = filter_roi(rois,fOPTIONS);

r = [froi{:}];
froi(strcmp({r.name},'blancL')|strcmp({r.name},'blancR'))=[];
%% ---------------------------


%% DISPLAYS ROIS (one window per roi)
%-----------------------------------------------------------------------------
dOPTIONS.clr = jet(numel(unique(m_table_effect(:,1)))); % distinct colors for pt
dOPTIONS.win_noedges = [-401,1601] ;
display_roi(froi,dOPTIONS);
%% ---------------------------


%% DISPLAYS SUMMARY (MEAN ROIS : ONE WINDOW)
%-----------------------------------------------------------------------------
pOPTIONS.thresh =3; % -1 for no color chronological organization
pOPTIONS.nsub =3; % number of subplot
pOPTIONS.title = '';

[labels_o, colorm] = display_summary_roi(froi,pOPTIONS) ;
%% ---------------------------

%% ASD : 2017/12/20 explo with FXA AT RA 
pOPTIONS.nsub =1; % number of subplot

%***************************************************
% Here starts the CONTRAST analysis 
%***************************************************

%% Load MIA's working directories per condition
cOPTIONS.indir1 = '/Volumes/wolf/GAMMABLOCK_contrast_DATA/Results_HOMOG'; % CONDITION 1
cOPTIONS.indir2 = '/Volumes/wolf/GAMMABLOCK_contrast_DATA/Results_HETEROG';% CONDITION 2
cOPTIONS.mtg = 'bipolar';

% rOPTIONS.mtg = 'monopolar';
cOPTIONS.freq = cell2mat(grpOPTIONS.freqid(getOPTIONS.freq));
cOPTIONS.win_noedges = [-401,1601] ;

%% Prepare the signals for statistics (Get the trials signals when needed)
tic ; [cfroi] = browse_csignal(froi,cOPTIONS);toc

%% Compute statistical thrersholds (permutation based) : LONG PROCESSING....
pOPTIONS.win_noedges = [-401,1601] ;
pOPTIONS.threshp = 0.05;
pOPTIONS.smoth = 20 ;
pOPTIONS.nperm =1000 ;
tic ; [proi] = stats_permutations_rois(cfroi,pOPTIONS);toc

save('proi','proi'); % save the result at this stage beacause it was so long to process... 


%% Compute the statistics
dOPTIONS.clr = jet(numel(unique(m_table_as(:,1)))) ; 
dOPTIONS.win_noedges = [-401,1601] ;
dOPTIONS.threshp = 0.05;
dOPTIONS.smoth = 20 ; 
[lroi] = rois2pvalues(proi,dOPTIONS);

%% Display thresholded tvalues (uncorrected, corrected duration, corrected sumtval)
dOPTIONS.clr = jet(numel(unique(m_table_as(:,1)))) ; 
dOPTIONS.win_noedges = [-401,1601] ;
dOPTIONS.threshp = 0.05;
dOPTIONS.smoth = 20 ; 
display_rois_conditions_fdr(lroi,dOPTIONS);
