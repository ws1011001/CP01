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
% setup working path
mdir = '/media/wang/BON/Projects/CP01';
ddir = fullfile(mdir,'SEEG_LectureVWFA','derivatives');      % derivatives in the BIDS structure
bdir = fullfile(ddir,'brainstorm_SEEG_LectureVWFA','data');  % path to brainsotrm database
% read the subjects list
fid = fopen(fullfile(mdir,'CP01_subjects.txt'));
subjects = textscan(fid,'%s');
fclose(fid);
subjects = subjects{1};  % the subjects list
n = length(subjects);
% set processing parameters
montagetype = 'bipolar_2';
conditions = {'AAp', 'AAt', 'AVp', 'AVt', 'VAp', 'VAt', 'VVp', 'VVt', 'Ap', 'Vp'};
condtype = [repmat({'repetition'}, 1, 8), repmat({'prime'}, 1, 2)];
ncond = length(conditions);
nperm = 5000;       % the number of randomizations in permutation tests
fdr_p = 0.05;       % alpha level
fdr_duration = 20;  % duration (in ms) for FDR correction in the time domain
%% ---------------------------

%% perform permutation tests on individual ERPs
for i = 1:n
  subj = subjects{i};
  % read working filenames
  ftmp = fullfile(bdir, subj, sprintf('%s_%s_working-filenames.mat', subj, ptoken));  
  load(ftmp);
  % group all conditions
  sConds = [{bFiles.primes.Condition}, {bFiles.repetitions.Condition}];
  sFiles = cellfun(@(x) fullfile(bdir, x), [{bFiles.primes.FileName}, {bFiles.repetitions.FileName}], 'UniformOutput', 0);
  % calculate the mean, SD, and SE of ERPs for each condition
  for icond = 1:ncond
    cndt = conditions{icond};
    cidx = ismember(sConds, [cndt '_' montagetype]);
    bFiles.(cndt) = sFiles(cidx);
    fprintf('Calculate the mean, SD, and SE of ERPs for the condition %s for subject %s. \n', cndt, subj);
    % calculate average and variation Process: By trial group (folder average) i.e. 'avgtype', 5
    ERPs.(cndt).favg = bst_process('CallProcess', 'process_average', bFiles.(cndt), [], 'avgtype', 5, 'avg_func', 1, ...  % Arithmetic average:  mean(x)
                                   'weighted', 0, 'keepevents', 0);      
    ERPs.(cndt).fstd = bst_process('CallProcess', 'process_average', bFiles.(cndt), [], 'avgtype', 5, 'avg_func', 4, ...  % Standard deviation:  sqrt(var(x))
                                   'weighted', 0, 'keepevents', 0);      
    ERPs.(cndt).fste = bst_process('CallProcess', 'process_average', bFiles.(cndt), [], 'avgtype', 5, 'avg_func', 5, ...  % Standard error:  sqrt(var(x)/N)
                                   'weighted', 0, 'keepevents', 0);
    % determine the baseline and time-window according to the condition
    bl = baseline.(condtype{icond});
    tw = timewindow.(condtype{icond});
    % one sample T-test with baseline : Y = mean_trials(X); t = (Y - mean_time(Y(baseline)) / std_time(Y(baseline))); df = length(baseline) - 1
    fprintf('Conduct one sample T-test for the condition %s for subject %s. \n', cndt, subj);
    ERPs.(cndt).ftest = bst_process('CallProcess', 'process_test_baseline', bFiles.(cndt), [], 'baseline', bl, 'timewindow', tw, ...
                                    'sensortypes', 'SEEG', 'isabs', 0, 'avgtime', 0, 'avgrow', 0, 'test_type', 'ttest_baseline', 'tail', 'two');
    % apply statistic threshold: alpha = 0.05 (FDR: control time)
    ERPs.(cndt).ftest_fdr = bst_process('CallProcess', 'process_extract_pthresh', ERPs.(cndt).ftest, [], 'pthresh', fdr_p, ...
                                        'durthresh', fdr_duration, 'correction', 3, 'control1', 0, 'control2', 1, 'control3', 0);
  end
  % do permutation tests with paired samples : T = mean(A-B) / std(A-B) * sqrt(n)
  fprintf('Conduct permutation tests [AAp-AAt | AVp-AVt | VAp-VAt | VVp-VVt] for subject %s. \n', subj);
  ERPs.perm.AAp2AAt = bst_process('CallProcess', 'process_test_permutation2p', bFiles.AAp, bFiles.AAt, 'timewindow', timewindow.repetition, ...
                                  'sensortypes', 'SEEG', 'isabs', 0, 'avgtime', 0, 'avgrow', 0, 'iszerobad', 1, 'Comment', '', ...
                                  'test_type', 'ttest_paired', 'randomizations', nperm, 'tail', 'two');
  ERPs.perm.AVp2AVt = bst_process('CallProcess', 'process_test_permutation2p', bFiles.AVp, bFiles.AVt, 'timewindow', timewindow.repetition, ...
                                  'sensortypes', 'SEEG', 'isabs', 0, 'avgtime', 0, 'avgrow', 0, 'iszerobad', 1, 'Comment', '', ...
                                  'test_type', 'ttest_paired', 'randomizations', nperm, 'tail', 'two');
  ERPs.perm.VAp2VAt = bst_process('CallProcess', 'process_test_permutation2p', bFiles.VAp, bFiles.VAt, 'timewindow', timewindow.repetition, ...
                                  'sensortypes', 'SEEG', 'isabs', 0, 'avgtime', 0, 'avgrow', 0, 'iszerobad', 1, 'Comment', '', ...
                                  'test_type', 'ttest_paired', 'randomizations', nperm, 'tail', 'two');
  ERPs.perm.VVp2VVt = bst_process('CallProcess', 'process_test_permutation2p', bFiles.VVp, bFiles.VVt, 'timewindow', timewindow.repetition, ...
                                  'sensortypes', 'SEEG', 'isabs', 0, 'avgtime', 0, 'avgrow', 0, 'iszerobad', 1, 'Comment', '', ...
                                  'test_type', 'ttest_paired', 'randomizations', nperm, 'tail', 'two');                                
  % apply statistic threshold: alpha = 0.05 (FDR: control time)
  contrasts = fieldnames(ERPs.perm);
  for icont = 1:length(contrasts)
    ERPs.perm.([contrasts{icont} '_fdr']) = bst_process('CallProcess', 'process_extract_pthresh', ERPs.perm.(contrasts{icont}), [], ...
                                                        'pthresh', 0.05, 'durthresh', fdr_duration, 'correction', 3, 'control1', 0, 'control2', 1, 'control3', 0);
  end
  % output the full path of working files
  save(ftmp, 'ERPs', 'bFiles', 'nperm', 'fdr_p', 'fdr_duration', '-append');
end
%% ---------------------------

% %% visualize ERPs for each contrast
% for i=1:n
%   subj=subjects{i};
%   sdir=fullfile(bdir,subj);
%   % read working path
%   fworking=fullfile(sdir,sprintf('%s_ses-01_task-RS_run-01_ieeg_ERPs_working-files-%s.mat',subj,montagetype));
%   load(fworking);
%   % combine ERPs average and variations
%   for icond=1:ncond
%     con=conditions{icond};
%     davg=load(fullfile(bdir,ERPs.(con).favg.FileName),'F');
%     dstd=load(fullfile(bdir,ERPs.(con).fstd.FileName),'F');
%     dste=load(fullfile(bdir,ERPs.(con).fste.FileName),'F');
%     ERPs.(con).avg=davg.F;
%     ERPs.(con).std=dstd.F;
%     ERPs.(con).ste=dste.F;
%   end
% end
% %% ---------------------------








