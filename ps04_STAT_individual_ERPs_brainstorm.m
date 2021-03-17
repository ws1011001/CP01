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
timewindow=[-0.2,0.8];  % time-window of the epoch i.e. -200 ms to 800 ms
nperm=10000;            % the number of randomizations in permutation tests
%% ---------------------------

%% perform permutation tests on individual ERPs
for i=1:n
  subj=subjects{i};
  sdir=fullfile(bdir,subj);
  fworking=fullfile(sdir,sprintf('%s_ses-01_task-RS_run-01_ieeg_ERPs_working-files-%s.mat',subj,montagetype));
  % calculate the mean, SD, and SE of ERPs for each condition
  for icond=1:ncond
    fprintf('Calculate the mean, SD, and SE of ERPs for the condition %s for subject %s ......\n',conditions{icond},subj);
    % extract trials
    cdir=fullfile(sdir,sprintf('%s_%s',conditions{icond},montagetype));
    ftrials=dir(fullfile(cdir,'data*.mat'));
    ntrials=length(ftrials);
    sFiles=cellfun(@(x) fullfile(cdir,x),{ftrials.name},'UniformOutput',0);
    working_files.(conditions{icond}).sFiles=sFiles;  % save the full path to working files
    % Process: By trial group (folder average) i.e. 'avgtype', 5
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
  end
  % do permutation tests (paired) T = mean(A-B) / std(A-B) * sqrt(n)
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
  % output the full path of working files
  save(fworking,'working_files');
end
%% ---------------------------









