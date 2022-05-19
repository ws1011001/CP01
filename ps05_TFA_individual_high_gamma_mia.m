%% ---------------------------
%% [script name] ps05_TFA_individual_high_gamma_mia.m
%%
%% SCRIPT to do statistical comparisons on high-gamma activity at the individual level.
%%
%% By Shuai Wang, [date] 2021-03-19
%%
%% ---------------------------
%% Notes: - This script is modified from ASD's script seeg_pipeline_MAIN.m
%%   
%%
%% ---------------------------


%% clean up
close all
clear
clc
%% ---------------------------

%% set environment (packages, functions, working path etc.)
% % if use HPC
% mdir = '/CP01';
% addpath(genpath('/CP01/mia'));
% setup working path
protocol = 'mia_SEEG_LectureVWFA';
mdir = '/media/wang/BON/Projects/CP01';     % the project folder
wdir = fullfile(mdir, 'SEEG_LectureVWFA');  % the main working folder 
rdir = fullfile(wdir, 'derivatives');       % path to the BIDS derivatives
bdir = fullfile(rdir, protocol);            % path to the MIA database
% read the subjects list
fid = fopen(fullfile(mdir,'CP01_subjects_lvOT.txt'));
subjects = textscan(fid,'%s');
fclose(fid);
subjects = subjects{1};  % the subjects list
n = length(subjects);
% set processing parameters
freq_step = 10;   % frequency decomposition step (Hz)
freq_lobd = 40;   % lower boundary
freq_mdbd = 80;   % the gamma frequency to separate low and high gamma
freq_upbd = 150;  % upper boundary
freq_bands = [freq_lobd, freq_mdbd; freq_mdbd, freq_upbd; freq_lobd, freq_upbd];
nbands = size(freq_bands, 1);  % number of frequency bands
% statistic parameters
stats.nboot = 5000;
stats.alpha = 0.05;
stats.nperm = 1000;
stats.smooth = 20;  % this number x 3 < number of ROIs (for successful smoothing)
stats.threshp = 0.05;
% OPTIONS.win_noedges = [-0.2 0.6];  % plot
% OPTIONS.ylab1 = 'Gamma Activity (Z-score)';  % plot

% baseline = [-0.2,-0.01];  % baseline from -200ms to -10ms
% timewindow = [-0.2,0.8];  % time-window of the epoch i.e. -200 ms to 800 ms
% nperm = 10000;            % the number of randomizations in permutation tests
% fdr_duration = 20;        % duration (in ms) for FDR correction in the time domain
conditions = {'AAp', 'AAt', 'AVp', 'AVt', 'VAp', 'VAt', 'VVp', 'VVt', ...
              'Ap', 'Vp', 'Ap_long', 'Vp_long', 'WA', 'WV', 'PA', 'PV', 'WPA', 'WPV', ...
              'SA', 'SV', 'FA', 'FV'};
ncond = length(conditions);
pairs = {{'AAp', 'AAt'}, {'AVp', 'AVt'}, {'VAp', 'VAt'}, {'VVp', 'VVt'}, ...  % repetition 
         {'AAp', 'AVp'}, {'VVp', 'VAp'}, {'AAt', 'VAt'}, {'VVt', 'AVt'}, ...  % between-prime, between-target
         {'FA', 'FV'}};  % lip-movements with audio vs. without audio 
% set switches
isPlotChn = true;
isGetROIs = false;
isCompCon = false;
%% ---------------------------

%% initialize the working file
% create protocol if not exist
if ~exist(bdir, 'dir')
  fprintf('Create a MIA protocol %s in the folder %s. \n', protocol, rdir);
  mkdir(rdir, protocol);  % create a database in the BIDS derivatives folder
end
% initialize the working file
ftmp = fullfile(bdir, 'ps05_TFA_working.mat');
% take a note
pdate = datetime('now','Format', 'yyyy-MM-dd''T''HH:mm:SS');  % the date of the present processing
%pnote = 'Run the pre-processing pipeline using MIA for 10 subjects (from sub-01 to sub-10).';
%pnote = 'Do between-condition comparisons for 10 subjects (from sub-01 to sub-10).';
pnote = 'Run pre-processing and statistical analysis for sub-03, sub-09 and sub-10 who have channels within the left-vOT.';
fprintf('%s \nStart at %s. \n\n', pnote, pdate);
if exist(ftmp, 'file')
  load(ftmp);  % read up the working file
  commits(size(commits, 1) + 1, :) = {pnote, sprintf('%s', pdate)};  % add commit
  save(ftmp, 'commits', '-append');
else
  commits = {pnote, sprintf('%s', pdate)};  % add commit
  save(ftmp, '*dir', 'subjects', 'n', 'freq*', 'nbands', 'conditions', 'ncond', 'commits')
end
%% ---------------------------

%% Import data from brainstorm database
if ~exist('miacfg', 'var')
  bst_db = fullfile(rdir, 'brainstorm_SEEG_LectureVWFA', 'data');
  for i = 1:n
    subj = subjects{i};
    sdir = fullfile(bdir, subj);
    if ~exist(sdir, 'dir'); mkdir(sdir); end
    fprintf('Import BST data into MIA to extract LFP signals for subject %s. \n', subj);
    % copy BST data to the MIA database
    for j = 1:ncond
      icond = conditions{j};
      if ~exist(fullfile(sdir, icond), 'dir'); mkdir(sdir, icond); end  % create a folder for this condition
      copyfile(fullfile(bst_db, subj, icond), fullfile(sdir, icond));
    end
    % import BST data into MIA
    config.maindir = sdir;
    config.outdir = sdir;
    miacfg.data(i).subject = subj;
    miacfg.data(i).LFP = mia_s1_extract_bst_data(config);  % edited mia_s1_extract_bst_data
  end
  miacfg.LFP = true;
  % update the working file
  save(ftmp, 'miacfg', '-append');  
end
%% ---------------------------

%% Extract frequency
if ~isfield(miacfg, 'oscillation')
  for i = 1:n
    subj = subjects{i};
    sdir = fullfile(bdir, subj);  
    for iband = 1:nbands
      freql = freq_bands(iband, 1);  % lower band
      frequ = freq_bands(iband, 2);  % upper band           
      fprintf('Extract oscillation amplitudes within frequency band %d to %d Hz for subject %s. \n', freql, frequ, subj);
      % define config
      config.outdir = sdir;
      config.removeEvoked = true;  % remove evoked responses
      config.freqs = freql:freq_step:frequ;
      config.modetf = 'Morlet';
      config.ncycles = 7;
      % TFA on monopolar signals
      config.mtg = 'monopolar';
      TFA = sprintf('TFA_%s_%s_%d_%d_%d', config.mtg, config.modetf, freq_step, freql, frequ);
      miacfg.data(i).(TFA) = mia_s4_compute_tf_bandwise(config);  % edited mia_s4_compute_tf_bandwise
      % TFA on bipolar signals
      config.mtg = 'bipolar';
      TFA = sprintf('TFA_%s_%s_%d_%d_%d', config.mtg, config.modetf, freq_step, freql, frequ);
      miacfg.data(i).(TFA) = mia_s4_compute_tf_bandwise(config);      
    end
  end
  miacfg.oscillation = true;
  % update the working file
  save(ftmp, 'miacfg', '-append');   
end
%% ---------------------------

%% single condition statistics
if ~isfield(miacfg, 'onesample')
  for i = 1:n
    subj = subjects{i};
    sdir = fullfile(bdir, subj);  
    for iband = 1:nbands
      freql = freq_bands(iband, 1);  % lower band
      frequ = freq_bands(iband, 2);  % upper band           
      for icond = 1:ncond
        cdir = fullfile(sdir, conditions{icond});
        fprintf('Estimate one-sample significance for frequency band %d to %d Hz for condition %s for subject %s. \n', freql, frequ, conditions{icond}, subj);
        % statistics on monopolar signals
        ftfa = fullfile(cdir, sprintf('%s_monopolar_morlet_data_%d_%d_%d_removeEvoked.mat', conditions{icond}, freq_step, freql, frequ));
        mia_s5_compute_stats(ftfa, stats);  % edited mia_s5_compute_stats
        % statistics on bipolar signals
        ftfa = fullfile(cdir, sprintf('%s_bipolar_morlet_data_%d_%d_%d_removeEvoked.mat', conditions{icond}, freq_step, freql, frequ));
        mia_s5_compute_stats(ftfa, stats);  % edited mia_s5_compute_stats     
      end
    end    
  end
  miacfg.onesample = true;
  % update the working file
  save(ftmp, 'miacfg', 'stats', '-append');  
end
%% ---------------------------

%% between-condition statistics (permutation tests)
if ~isfield(miacfg, 'permutation')
  for i =1:n
    subj = subjects{i};
    sdir = fullfile(bdir, subj);
    fprintf('Perform permutation tests between conditions for each channel for subject %s. \n', subj);   
    for mtg = ["monopolar", "bipolar"]  % for each montage
      % do comparisons for each frequency band
      for iband = 1:nbands
        freql = freq_bands(iband, 1);  % lower band
        frequ = freq_bands(iband, 2);  % upper band
        ptoken = sprintf('%s_morlet_data_%d_%d_%d_removeEvoked', mtg, freq_step, freql, frequ);
        % extract signals for each channel per condition
        for icond = 1:ncond
          cdir = fullfile(sdir, conditions{icond});
          ftfa = fullfile(cdir, sprintf('%s_%s.mat', conditions{icond}, ptoken));
          vtfa = who('-file', ftfa);
          if ismember('zs_chn', vtfa)
            load(ftfa, 'zs_chn');
          else
            zs_chn = chn_signals(ftfa);  % channel-wise data structure
            save(ftfa, 'zs_chn', '-append');
          end
          % group conditions
          signals.(conditions{icond}) = zs_chn.signals;
        end
        % compare conditions
        cfg = stats;
        cfg.labels      = zs_chn.labels(1, :);
        cfg.time        = zs_chn.time;
        for ipair = 1:length(pairs)
          cpair = pairs{ipair};
          cdir = fullfile(sdir, 'contrasts', sprintf('%s-%s', cpair{1}, cpair{2}));
          if ~exist(cdir, 'dir'); mkdir(cdir); end
          % two-sample permutation test            
          fperm = fullfile(cdir, sprintf('%s-%s_%s.mat', cpair{1}, cpair{2}, ptoken));
          fprintf('Perform two-sample permutation tests between %s and %s on data %s. \n', cpair{1}, cpair{2}, ptoken);
          tic  % start permutation test
          zs_chn_perm = roi_stats_permutations(signals.(cpair{1}), signals.(cpair{2}), cfg);
          save(fperm, 'zs_chn_perm');
          toc  % done this test. e.g., 277s for FA-FV 1st band sub-10
        end  
      end
    end
  end
  miacfg.permutation = true;
  % update the working file
  save(ftmp, 'miacfg', '-append');   
end
%% ---------------------------



%% extract signals for each ROI
if isGetROIs
  for i = 1:n
    subj = subjects{i};
    sdir = fullfile(bdir, subj);
    fprintf('Extract %s %s signals for each ROI for subject %s. \n', subj);
    % read the table of MNI coordinates of monopolar channels
    fmni=fullfile(mdir, 'SEEG_LectureVWFA', subj, 'anat', sprintf('%s_space-MNI_contacts-monopolar.csv', subj));
    if ~exist(fmni, 'file')
      fcoi = fullfile(mdir, 'SEEG_LectureVWFA', subj, 'anat/elec2atlas.mat');
      mono_mni = mni_coi2csv(fcoi, fmni);
    else
      mono_mni = readtable(fmni);
    end  
    mono_mni_xyz = [mono_mni.x mono_mni.y mono_mni.z];
    % calculate the signals for each ROI
    for icond = 1:ncond
      cdir = fullfile(sdir, conditions{icond});
      ftfa = fullfile(cdir, sprintf('%s_%s.mat', conditions{icond}, ptoken));
      data = load(ftfa, 'labels', 'zs', 'Time');
      % get bipolar labels and coordinates
      [bipo_mni, bipo_labels, bipo_flags] = mni_mono2bipolar(mono_mni_xyz, mono_mni.contacts, data.labels);
      roi.bipolar_mni = bipo_mni;
      roi.bipolar_contacts = bipo_labels;
      roi.bipolar_flags = bipo_flags;
      % get regions of interest
      roi.AAL3 = mni_find_regions(bipo_mni, 'AAL3');
      % get ROI-specific signals
      roi.signals = roi_average_channels(data.zs(bipo_flags, :, :), bipo_mni, bipo_labels, roi.AAL3);
      roi.time = data.Time;
      % save ROI data
      save(ftfa, 'roi', '-append');
    end
  end
end
%% ---------------------------

%% between-condition comparisons for each ROI
if isCompCon  
  for i = 1:n
    subj = subjects{i};
    sdir = fullfile(bdir,subj);
    fprintf('Perform permutation tests between conditions for each ROI for subject %s. \n', subj);
    % read ROI-based signals of different conditions
    for icond = 1:ncond
      cdir = fullfile(sdir, conditions{icond});
      ftfa = fullfile(cdir, sprintf('%s_%s.mat', conditions{icond}, ptoken)); 
      load(ftfa, 'roi');
      signals.(conditions{icond}) = roi.signals;
    end
    % compare conditions
    for ipair = 1:length(pairs)
      cpair = pairs{ipair};
      rdir = fullfile(sdir, sprintf('%s_ROIs_%s-%s', subj, cpair{1}, cpair{2}));
      if ~exist(rdir, 'dir'); mkdir(rdir); end
      fperm = fullfile(rdir, sprintf('ROIs_%s-%s_%s.mat', cpair{1}, cpair{2}, ptoken));
      if ~exist(fperm, 'file')
        roi_perm = roi_stats_permutations(signals.(cpair{1}), signals.(cpair{2}), OPTIONS);
        save(fperm, 'roi_perm');
      else
        load(fperm, 'roi_perm');
      end
      roi_plot_conditions(roi_perm, cpair, OPTIONS, rdir); 
    end 
  end
end
%% ---------------------------



