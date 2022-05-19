%% ---------------------------
%% [script name] ps06_VIS_individual_ERP_and_HGA.m
%%
%% SCRIPT to visualize statistical comparisons of ERP and HGA at the individual level.
%%
%% By Shuai Wang, [date] 2022-04-20
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
%protocol = 'mia_SEEG_LectureVWFA';
protocol = 'brainstorm_SEEG_LectureVWFA';   % the name of the database
mdir = '/media/wang/BON/Projects/CP01';     % the project folder
wdir = fullfile(mdir, 'SEEG_LectureVWFA');  % the main working folder 
rdir = fullfile(wdir, 'derivatives');       % path to the BIDS derivatives
bdir = fullfile(rdir, protocol);            % path to the MIA database
ddir = fullfile(bdir, 'data');              % path to the functional database
vdir = fullfile(mdir, 'results', 'visualization');  % path to save images
% read the subjects list
fid = fopen(fullfile(mdir,'CP01_subjects.txt'));
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
conditions = {'AAp', 'AAt', 'AVp', 'AVt', ...
              'VAp', 'VAt', 'VVp', 'VVt', ...
              'Ap', 'Vp', 'Ap_long', 'Vp_long', ...
              'WA', 'WV', 'PA', 'PV', 'WPA', 'WPV', ...
              'SA', 'SV', 'FA', 'FV'};
clabels = {'Auditory-auditory Prime', 'Auditory-auditory Target', 'Auditory-visual Prime', 'Auditory-visual Target', ...
           'Visual-auditory Prime', 'Visual-auditory Target', 'Visual-visual Prime', 'Visual-visual Target', ...
           'Auditory Prime', 'Visual Prime', 'Auditory Prime', 'Visual Prime', ...
           'Auditory Words', 'Visual Words', 'Auditory Pseudowords', 'Visual Pseudowords', 'Auditory W+P', 'Visual W+P', ...
           'Auditory Scrambled', 'Visual Consonants', 'Lip-movement Video with Audio', 'Lip-movement Video without Audio'};
ncond = length(conditions);
pairs = {{'AAp', 'AAt'}, {'AVp', 'AVt'}, {'VAp', 'VAt'}, {'VVp', 'VVt'}, ...  % repetition 
         {'AAp', 'AVp'}, {'VVp', 'VAp'}, {'AAt', 'VAt'}, {'VVt', 'AVt'}, ...  % between-prime, between-target
         {'FA', 'FV'}};  % lip-movements with audio vs. without audio 
% figure settings
cfgERP.fontsize  = 18;
cfgERP.linewidth = 2;
cfgERP.smoothing = 1;
cfgERP.xscale    = 0.1;
cfgERP.legend    = 'far-left';
cfgERP.sigbar    = true;
cfgHGA.threshp     = 0.05;
cfgHGA.plot_signif = true;  % plot only significant results
% set switches
isFindCOI = false;
isPlotERP = true;
isPlotHGA = false;
%% ---------------------------

%% Find out channels of interest according to the fMRI results
if isFindCOI
  % ROIs defined by the fMRI results
  lvot.groi = fullfile(rdir, 'masks', 'group_space-MNI152NLin2009cAsym_mask-lvOT-visual.mat');
  lvot.gbox = [-49.1 -35.1 -43 -18.5 -27.8 -17.2];  % 3dClusterize -inset group_space-MNI152NLin2009cAsym_mask-lvOT-visual.nii.gz -ithr 0 -idat 0 -NN 1 -within_range 0.5 1.5
  lvot.box0 = [-56.2 -33.4 -64.1 -29 -27.8 -12];
  lvot.box7 = [-63.2 -26.4 -71.1 -22 -34.8 -5];
  lvot.box8 = [-64.2 -25.4 -72.1 -21 -35.8 -4];
  templates = fieldnames(lvot);
  % Find out COI for each patient
  for i = 1:n
      subj = subjects{i};
      adir = fullfile(wdir, subj, 'anat');  % anatomical data
      fprintf('Find out Channels of Interest for subject %s. \n', subj);
      mni_mono = load(fullfile(adir, 'elec2atlas.mat'), 'coi');         % monopolar channels
      mni_bipo = load(fullfile(adir, 'elecbipolar2atlas.mat'), 'coi');  % bipolar channels
      for j = 1:length(templates)
        temp = templates{j};
        [~, I] = mni_find_regions(mni_mono.coi.elecpos_mni, lvot.(temp));
        mni_mono.(temp) = mni_mono.coi.label(I);
        [~, I] = mni_find_regions(mni_bipo.coi.elecpos_mni, lvot.(temp));
        mni_bipo.(temp) = mni_bipo.coi.label(I);      
      end
      % Output COIs
      fcoi = fullfile(adir, sprintf('%s_COIs_mask-lvOT.mat', subj));
      save(fcoi, 'mni_mono', 'mni_bipo');
  end
end
%% ---------------------------

%% Visualize ERPs for single conditions
if isPlotERP
  for i = 1:n
    subj = subjects{i};
    adir = fullfile(wdir, subj, 'anat');             % anatomical data
    sdir = fullfile(vdir, 'timecourse_ERPs', subj);  % result folder
    if ~exist(sdir, 'dir'); mkdir(sdir); end    
    % Read the working file and COIs for this subject
    ferp = fullfile(ddir, sprintf('%s_ps04_ERP_working.mat', subj));  % the working file   
    fcoi = fullfile(adir, sprintf('%s_COIs_mask-lvOT.mat', subj));    % COIs
    load(ferp, 'ERPs', 'contrasts');
    load(fcoi, 'mni_mono');
    COIs = unique([mni_mono.groi; mni_mono.gbox; mni_mono.box0; mni_mono.box7; mni_mono.box8]);  % all COIs
    if ~isempty(COIs)
      for icond = 1:length(conditions) 
        c = conditions{icond};
        % Extrat ERP signals        
        davg = load(fullfile(ddir, ERPs.(c).avg.FileName), 'F', 'Time');
        dstd = load(fullfile(ddir, ERPs.(c).std.FileName), 'F');
        dste = load(fullfile(ddir, ERPs.(c).ste.FileName), 'F');
        dsig = load(fullfile(ddir, ERPs.(c).test_fdr.FileName), 'F');
        % Combine ERPs average and variations
        signals.(c).n   = length(ERPs.(c).trialdata);  % number of trials
        signals.(c).avg = davg.F .* 1e6;
        signals.(c).std = dstd.F .* 1e6;
        signals.(c).ste = dste.F .* 1e6;
        signals.(c).sig = dsig.F;  % T-values with FDR correction
        signals.(c).Time = davg.Time;  % in seconds
        % Plot ERP for a single condition
        ctoken = clabels(ismember(conditions, c));
        ctitle = sprintf('%s (N=%d)', ctoken{1}, signals.(c).n);
        ftoken = fullfile(sdir, sprintf('%s_ERP_%s', subj, c));
        cfgERP.Time = signals.(c).Time;
        channels = load(fullfile(ddir, ERPs.(c).test_fdr.ChannelFile), 'Channel');
        plot_ERPs(signals.(c).avg, signals.(c).ste, signals.(c).sig, {ctitle}, channels.Channel, COIs, cfgERP, ftoken);
      end
      % Plot ERPs for each contrasts
      for icont = 1:length(contrasts)
        % Define labels
        A = contrasts{icont, 1};  % condition A
        B = contrasts{icont, 2};  % condition B
        Atoken = clabels(ismember(conditions, A));
        Btoken = clabels(ismember(conditions, B));
        Atitle = sprintf('%s (N=%d)', Atoken{1}, signals.(A).n);
        Btitle = sprintf('%s (N=%d)', Btoken{1}, signals.(B).n);
        cont = sprintf('%s_%s', A, B);
        ftoken = fullfile(sdir, sprintf('%s_ERPs_%s', subj, cont));
        % Extract data
        cont_avg = cat(3, signals.(A).avg, signals.(B).avg);
        cont_ste = cat(3, signals.(A).ste, signals.(B).ste);
        cont_sig = load(fullfile(ddir, ERPs.(cont).test_fdr.FileName), 'F');
        % Plot ERPs
        cfgERP.Time = signals.(A).Time;
        channels = load(fullfile(ddir, ERPs.(cont).test_fdr.ChannelFile), 'Channel');
        plot_ERPs(cont_avg, cont_ste, cont_sig.F, {Atitle, Btitle}, channels.Channel, COIs, cfgERP, ftoken); 
      end
    end
  end  
end
%% ---------------------------

%% Visualize channel-wise between-condition comparisons in HGA
if isPlotHGA
  for i =1 :n
    subj = subjects{i};
    sdir = fullfile(bdir, subj);
    fprintf('Plot between-condition comparisons for each channel for subject %s. \n', subj); 
    % plot comparisons for each frequency band
    for iband = 1:nbands
      freql = freq_bands(iband, 1);  % lower band
      frequ = freq_bands(iband, 2);  % upper band
      ptoken = sprintf('bipolar_morlet_data_%d_%d_%d_removeEvoked', freq_step, freql, frequ);      

      cfgHGA.ylab1       = sprintf('Gamma Activity (%d-%d Hz)', freql, frequ);
      
      for ipair = 1:length(pairs)
        cpair = pairs{ipair};
        cdir = fullfile(sdir, 'contrasts', sprintf('%s-%s', cpair{1}, cpair{2}));
        % extarc parameters
        ftfa = fullfile(sdir, cpair{1}, sprintf('%s_%s.mat', cpair{1}, ptoken));  
        load(ftfa, 'zs_chn');
        cfgHGA.labels      = zs_chn.labels(1, :);
        cfgHGA.time        = zs_chn.time;
        cfgHGA.win_noedges = [round(zs_chn.time(1), 2) + 0.1, round(zs_chn.time(end), 2) - 0.1];        
        % plot two-sample permutation test            
        fperm = fullfile(cdir, sprintf('%s-%s_%s.mat', cpair{1}, cpair{2}, ptoken));
        load(fperm);
        cfgHGA.outputdir = fullfile(cdir, ptoken);
        cfgHGA.figprefix = subj;
        fprintf('Plot results for %s. \n', fperm);
        roi_plot_conditions(zs_chn_perm, cpair, cfgHGA);  
      end  
    end
  end  
end
%% ---------------------------