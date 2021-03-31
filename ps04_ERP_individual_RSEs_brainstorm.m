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
mdir='/media/wang/BON/Projects/CP01';
ddir=fullfile(mdir,'SEEG_LectureVWFA','derivatives');      % derivatives in the BIDS structure
bdir=fullfile(ddir,'brainstorm_SEEG_LectureVWFA','data');  % path to brainsotrm database
% read the subjects list
fid=fopen(fullfile(mdir,'CP01_subjects.txt'));
subjects=textscan(fid,'%s');
fclose(fid);
subjects=subjects{1};  % the subjects list
n=length(subjects);
% set processing parameters
montagetype='bipolar_2';
conditions={'AAp','AAt','AVp','AVt','VAp','VAt','VVp','VVt'};
ncond=length(conditions);
baseline=[-0.2,-0.01];  % baseline from -200ms to -10ms
timewindow=[-0.2,0.8];  % time-window of the epoch i.e. -200 ms to 800 ms
nperm=10000;            % the number of randomizations in permutation tests
fdr_duration=20;        % duration (in ms) for FDR correction in the time domain
%% ---------------------------

%% perform permutation tests on individual ERPs
for i=1:n
  subj=subjects{i};
  sdir=fullfile(bdir,subj);
  fworking=fullfile(sdir,sprintf('%s_ses-01_task-RS_run-01_ieeg_ERPs_working-files-%s.mat',subj,montagetype));
  if exist(fworking,'file')
    fprintf('Load previous working files for subject %s ......\n',subj);
    load(fworking)
  end
  % calculate the mean, SD, and SE of ERPs for each condition
  for icond=1:ncond
    fprintf('Calculate the mean, SD, and SE of ERPs for the condition %s for subject %s ......\n',conditions{icond},subj);
    % extract trials
    cdir=fullfile(sdir,sprintf('%s_%s',conditions{icond},montagetype));
    if ~exist(fworking,'file')
      ftrials=dir(fullfile(cdir,'data*trial*.mat'));
      sFiles=cellfun(@(x) fullfile(cdir,x),{ftrials.name},'UniformOutput',0);
      working_files.(conditions{icond}).sFiles=sFiles;  % save the full path to working files
    else
      sFiles=working_files.(conditions{icond}).sFiles;
    end
    % calculate average and variation Process: By trial group (folder average) i.e. 'avgtype', 5
    working_files.(conditions{icond}).favg=bst_process('CallProcess','process_average',...
                                                       sFiles,[],'avgtype',5,...
                                                       'avg_func',1,...  % Arithmetic average:  mean(x)
                                                       'weighted',0,'keepevents',0);      
    working_files.(conditions{icond}).fstd=bst_process('CallProcess','process_average',...
                                                       sFiles,[],'avgtype',5,...
                                                       'avg_func',4,...  % Standard deviation:  sqrt(var(x))
                                                       'weighted',0,'keepevents',0);      
    working_files.(conditions{icond}).fste=bst_process('CallProcess','process_average',...
                                                       sFiles,[],'avgtype',5,...
                                                       'avg_func',5,...  % Standard error:  sqrt(var(x)/N)
                                                       'weighted',0,'keepevents',0);
    % one sample T-test with baseline : Y = mean_trials(X); t = (Y - mean_time(Y(baseline)) / std_time(Y(baseline))); df = length(baseline) - 1
    fprintf('Conduct one sample T-test for the condition %s for subject %s ......\n',conditions{icond},subj);
    working_files.(conditions{icond}).ftest=bst_process('CallProcess','process_test_baseline',...
                                                        sFiles,[],'baseline',baseline,'timewindow',timewindow,...
                                                        'sensortypes','SEEG','isabs',0,'avgtime',0,'avgrow',0,...
                                                        'test_type','ttest_baseline','tail','two');
    % apply statistic threshold: alpha=0.05 (FDR: control time)
    working_files.(conditions{icond}).ftest_fdr=bst_process('CallProcess','process_extract_pthresh',...
                                                            working_files.(conditions{icond}).ftest,[],...
                                                            'pthresh',0.05,'durthresh',fdr_duration,...
                                                            'correction',3,'control1',0,'control2',1,'control3',0);
  end
  % do permutation tests with paired samples : T = mean(A-B) / std(A-B) * sqrt(n)
  fprintf('Conduct permutation tests [AAp-AAt | AVp-AVt | VAp-VAt | VVp-VVt] for subject %s ......\n',subj);
  working_files.perm2p.AAp2AAt=bst_process('CallProcess','process_test_permutation2p',...
                                           working_files.AAp.sFiles,working_files.AAt.sFiles,...
                                           'timewindow',timewindow,'sensortypes','SEEG','isabs',0,...
                                           'avgtime',0,'avgrow',0,'iszerobad',1,'Comment','',...
                                           'test_type','ttest_paired','randomizations',nperm,'tail','two');
  working_files.perm2p.AVp2AVt=bst_process('CallProcess','process_test_permutation2p',...
                                           working_files.AVp.sFiles,working_files.AVt.sFiles,...
                                           'timewindow',timewindow,'sensortypes','SEEG','isabs',0,...
                                           'avgtime',0,'avgrow',0,'iszerobad',1,'Comment','',...
                                           'test_type','ttest_paired','randomizations',nperm,'tail','two');
  working_files.perm2p.VAp2VAt=bst_process('CallProcess','process_test_permutation2p',...
                                           working_files.VAp.sFiles,working_files.VAt.sFiles,...
                                           'timewindow',timewindow,'sensortypes','SEEG','isabs',0,...
                                           'avgtime',0,'avgrow',0,'iszerobad',1,'Comment','',...
                                           'test_type','ttest_paired','randomizations',nperm,'tail','two');                                           
  working_files.perm2p.VVp2VVt=bst_process('CallProcess','process_test_permutation2p',...
                                           working_files.VVp.sFiles,working_files.VVt.sFiles,...
                                           'timewindow',timewindow,'sensortypes','SEEG','isabs',0,...
                                           'avgtime',0,'avgrow',0,'iszerobad',1,'Comment','',...
                                           'test_type','ttest_paired','randomizations',nperm,'tail','two'); 
  % apply statistic threshold: alpha=0.05 (FDR: control time)
  contrasts=fieldnames(working_files.perm2p);
  for icont=1:length(contrasts)
    working_files.perm2p.([contrasts{icont} '_fdr'])=bst_process('CallProcess','process_extract_pthresh',...
                                                                 working_files.perm2p.(contrasts{icont}),[],...
                                                                 'pthresh',0.05,'durthresh',fdr_duration,...
                                                                 'correction',3,'control1',0,'control2',1,'control3',0);
  end
  % output the full path of working files
  save(fworking,'working_files');
end
%% ---------------------------

%% visualize ERPs for each contrast
for i=1:n
  subj=subjects{i};
  sdir=fullfile(bdir,subj);
  % read working path
  fworking=fullfile(sdir,sprintf('%s_ses-01_task-RS_run-01_ieeg_ERPs_working-files-%s.mat',subj,montagetype));
  load(fworking);
  % combine ERPs average and variations
  for icond=1:ncond
    con=conditions{icond};
    davg=load(fullfile(bdir,working_files.(con).favg.FileName),'F');
    dstd=load(fullfile(bdir,working_files.(con).fstd.FileName),'F');
    dste=load(fullfile(bdir,working_files.(con).fste.FileName),'F');
    ERPs.(con).avg=davg.F;
    ERPs.(con).std=dstd.F;
    ERPs.(con).ste=dste.F;
  end
end
%% ---------------------------








