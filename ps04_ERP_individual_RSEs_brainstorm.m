%% ---------------------------
%% [script name] ps04_STAT_individual_ERPs_brainstorm.m
%%
%% SCRIPT to do statistical comparisions between different conditions for ERPs at the individual level.
%%
%% By Shuai Wang, [date] 2021-03-17
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
% initialize brainsotrm
brainstorm nogui   % starts the interface but hides it.
%brainstorm server  % run on headless servers (computation clusters with no screen attached)
% setup working path

% ddir = fullfile(mdir, 'SEEG_LectureVWFA', 'derivatives');      % derivatives in the BIDS structure
% bdir = fullfile(ddir, 'brainstorm_SEEG_LectureVWFA', 'data');  % path to brainsotrm database
% vdir = fullfile(mdir, 'results', 'visualization');             % path to save images

% setup working path
protocol = 'brainstorm_SEEG_LectureVWFA';   % the name of the database
mdir = '/media/wang/BON/Projects/CP01';     % the project folder
wdir = fullfile(mdir, 'SEEG_LectureVWFA');  % the main working folder 
rdir = fullfile(wdir, 'derivatives');       % path to the BIDS derivatives
bdir = fullfile(rdir, protocol);            % path to the BST database
ddir = fullfile(bdir, 'data');              % path to the functional database
% read the subjects list and information
fid = fopen(fullfile(mdir, 'CP01_subjects.txt'));
subjects = textscan(fid, '%s');
fclose(fid);
subjects = subjects{1};  % the subjects list
n = length(subjects);
subjects_info = readtable(fullfile(wdir, 'participants.tsv'), 'FileType', 'text');
% set processing parameters
load(fullfile(ddir, 'ps03_PREP_working.mat'));  % load up pre-processing parameters
contrasts = [{'AAp', 'AAt'}; {'AVp', 'AVt'}; {'VAp', 'VAt'}; {'VVp', 'VVt'}; ...  % RSE contrasts
             {'AAp', 'AVp'}; {'VVp', 'VAp'}; {'AAt', 'VAt'}; {'VVt', 'AVt'}; ...  % between-prime and between-target contrasts
             {'WV', 'PV'}; {'WV', 'SV'}; {'PV', 'SV'}; ...                        % visual stimuli
             {'FA', 'FV'}; {'FA', 'WA'}];                                         % lip movements 
nperm = 5000;    % the number of randomizations in permutation tests
fdr_p = 0.05;    % alpha level
fdr_t = 20;      % duration (in ms) for FDR correction in the time domain
% set switches
isAllTests = false;  % do statistical tests for all conditions and contrasts
isPlotCond = false;  % plot time-courses for single conditions 
isPlotCont = false;  % plot ERPs for between-condition contrasts
%% ---------------------------

%% Perform permutation tests on individual ERPs
if isAllTests
  for i = 1:n
    subj = subjects{i};
    sidx = find(ismember({sFiles.trialdata.subject}, subj));
    % Initialize working file
    ftmp = fullfile(ddir, sprintf('%s_ps04_ERP_working.mat', subj));  
    % Statistics for single condition
    for ievts = 1:length(sFiles.events_epoch)
      evts = sFiles.events_epoch{ievts};               % event category, e.g., repetition
      conditions = strsplit(events.(evts), ', ');      % conditions in the event category, e.g., AAp, AAt, AVp, AVt, VAp, VAt, VVp, VVt
      trialdata = sFiles.monopolar_norm(sidx).(evts);  % normalized monopolar trials
      bl = baseline.(evts);    % baseline for the category
      tw = timewindow.(evts);  % timewindow for the category
      for icond = 1:length(conditions)
        c = conditions{icond};           % a specific condition, e.g., AAp
        c_norm = sprintf('%s_orig', c);  % condition name for normalized data
        d = trialdata(ismember({trialdata.Condition}, c_norm));
        badtrl = load(fullfile(ddir, subj, c, 'brainstormstudy.mat'), 'BadTrials');
        if ~isempty(badtrl.BadTrials)
          badidx = cell2mat(cellfun(@(x) str2num(x(end-6:end-4)), badtrl.BadTrials, 'UniformOutput', 0));
        else
          badidx = [];
        end
        d(badidx) = [];  % remove bad trials
        % Calculate average and variation Process: By trial group (folder average) i.e. 'avgtype', 5
        ERPs.(c).avg = bst_process('CallProcess', 'process_average', d, [], 'avgtype', 5, 'avg_func', 1, 'weighted', 0, 'keepevents', 0);  % Arithmetic average:  mean(x)                                          
        ERPs.(c).std = bst_process('CallProcess', 'process_average', d, [], 'avgtype', 5, 'avg_func', 4, 'weighted', 0, 'keepevents', 0);  % Standard deviation:  sqrt(var(x))   
        ERPs.(c).ste = bst_process('CallProcess', 'process_average', d, [], 'avgtype', 5, 'avg_func', 5, 'weighted', 0, 'keepevents', 0);  % Standard error:  sqrt(var(x)/N)  
        % One-sample T-test with baseline : Y = mean_trials(X); t = (Y - mean_time(Y(baseline)) / std_time(Y(baseline))); df = length(baseline) - 1
        fprintf('Conduct one sample T-test for the condition %s for subject %s. \n', c, subj);
        ERPs.(c).test = bst_process('CallProcess', 'process_test_baseline', d, [], 'baseline', bl, 'timewindow', tw, ...
                                    'sensortypes', 'SEEG', 'isabs', 0, 'avgtime', 0, 'avgrow', 0, 'test_type', 'ttest_baseline', 'tail', 'two');
        % Apply statistic threshold: alpha = 0.05 (FDR: control time)
        ERPs.(c).test_fdr = bst_process('CallProcess', 'process_extract_pthresh', ERPs.(c).test, [], 'pthresh', fdr_p, ...
                                        'durthresh', fdr_t, 'correction', 3, 'control1', 0, 'control2', 1, 'control3', 0);  
        % Update parameters
        ERPs.(c).trialdata = d;  % only good trials
        ERPs.(c).baseline = bl;
        ERPs.(c).timewindow = tw;
        ERPs.(c).fdr_params = [fdr_p, fdr_t];
      end
    end
    % Between-condition comparisons
    for icont = 1:length(contrasts)
      A = contrasts{icont, 1};  % condition A
      B = contrasts{icont, 2};  % condition B
      cont = sprintf('%s_%s', A, B);
      % Permutation test with two samples : H0:(A=B), H1:(A<>B) T = (mean(A)-mean(B)) / (Sx * sqrt(1/nA + 1/nB)), Sx = sqrt(((nA-1)*var(A) + (nB-1)*var(B)) / (nA+nB-2))                                    
      ERPs.(cont).test = bst_process('CallProcess', 'process_test_permutation2', ERPs.(A).trialdata, ERPs.(B).trialdata, ...
                                     'timewindow', ERPs.(A).timewindow, 'sensortypes', 'SEEG', 'isabs', 0, 'avgtime', 0, ...
                                     'avgrow', 0, 'iszerobad', 1, 'Comment', '', 'test_type', 'ttest_equal', 'randomizations', nperm, 'tail', 'two');   
      % Apply statistic threshold: alpha = 0.05 (FDR: control duration 20 ms)  
      ERPs.(cont).test_fdr = bst_process('CallProcess', 'process_extract_pthresh', ERPs.(cont).test, [], 'pthresh', fdr_p, ...
                                         'durthresh', fdr_t, 'correction', 3, 'control1', 0, 'control2', 1, 'control3', 0);  
      % Update parameters
      ERPs.(cont).timewindow = ERPs.(A).timewindow;
      ERPs.(cont).nperm = nperm;
      ERPs.(cont).fdr_params = [fdr_p, fdr_t];
    end
    % Save the working files
    save(ftmp, 'subj', 'sidx', 'contrasts', 'ERPs');
  end
end
%% ---------------------------












