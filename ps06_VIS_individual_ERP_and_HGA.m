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
conditions = {'AAp', 'AAt', 'AVp', 'AVt', 'VAp', 'VAt', 'VVp', 'VVt', ...
              'Ap', 'Vp', 'Ap_long', 'Vp_long', 'WA', 'WV', 'PA', 'PV', 'WPA', 'WPV', ...
              'SA', 'SV', 'FA', 'FV'};
ncond = length(conditions);
pairs = {{'AAp', 'AAt'}, {'AVp', 'AVt'}, {'VAp', 'VAt'}, {'VVp', 'VVt'}, ...  % repetition 
         {'AAp', 'AVp'}, {'VVp', 'VAp'}, {'AAt', 'VAt'}, {'VVt', 'AVt'}, ...  % between-prime, between-target
         {'FA', 'FV'}};  % lip-movements with audio vs. without audio 
OPTIONS.threshp = 0.05;
% set switches
isPlotChn = true;
%% ---------------------------

%% visualize channel-wise between-condition comparisons
if isPlotChn
  for i =1 :n
    subj = subjects{i};
    sdir = fullfile(bdir, subj);
    fprintf('Plot between-condition comparisons for each channel for subject %s. \n', subj); 
    % plot comparisons for each frequency band
    for iband = 1:nbands
      freql = freq_bands(iband, 1);  % lower band
      frequ = freq_bands(iband, 2);  % upper band
      ptoken = sprintf('bipolar_morlet_data_%d_%d_%d_removeEvoked', freq_step, freql, frequ);      
      % extarc parameters
      ftfa = fullfile(sdir, 'AAp', sprintf('AAp_%s.mat', ptoken));  
      load(ftfa, 'zs_chn');
      OPTIONS.ylab1       = sprintf('Gamma Activity (%d-%d Hz)', freql, frequ);
      OPTIONS.labels      = zs_chn.labels(1, :);
      OPTIONS.time        = zs_chn.time;
      OPTIONS.win_noedges = [round(zs_chn.time(1), 2) + 0.1, round(zs_chn.time(end), 2) - 0.1];
      OPTIONS.plot_signif = true;  % plot only significant results
      for ipair = 1:length(pairs)
        cpair = pairs{ipair};
        cdir = fullfile(sdir, 'contrasts', sprintf('%s-%s', cpair{1}, cpair{2}));
        % plot two-sample permutation test            
        fperm = fullfile(cdir, sprintf('%s-%s_%s.mat', cpair{1}, cpair{2}, ptoken));
        load(fperm);
        OPTIONS.outputdir = fullfile(cdir, ptoken);
        OPTIONS.figprefix = subj;
        fprintf('Plot results for %s. \n', fperm);
        roi_plot_conditions(zs_chn_perm, cpair, OPTIONS);  
      end  
    end
  end  
end
%% ---------------------------