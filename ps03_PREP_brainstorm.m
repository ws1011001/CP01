%% ---------------------------
%% [script name] ps04_PREP_brainstorm.m
%%
%% SCRIPT to do pre-processing by using brainstorm.
%%
%% By Shuai Wang, [date] 2021-04-18
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
ddir = fullfile(mdir, 'SEEG_LectureVWFA', 'derivatives');      % derivatives in the BIDS structure
bdir = fullfile(ddir, 'brainstorm_SEEG_LectureVWFA', 'data');  % path to brainsotrm database
% read the subjects list
fid = fopen(fullfile(mdir, 'CP01_subjects.txt'));
subjects = textscan(fid, '%s');
fclose(fid);
subjects = subjects{1};  % the subjects list
n = length(subjects);
% set processing parameters
timewindow = [-0.4, 0.8];  % epoch window
%% ---------------------------

%% import raw data
%% ---------------------------

%% basic pre-processing : (PSD), notch, band-pass
%% ---------------------------

%% manipulate events
% read filenames
ptoken = 'ses-01_task-RS_run-01_ieeg_notch_high';
sFiles = cellfun(@(x) fullfile(bdir, x, sprintf('@raw%s_%s/data_0raw_%s_%s.mat', x, ptoken, x, ptoken)), subjects, 'UniformOutput', 0);
% merge prime trials
events_merge = [{'AAp, AVp', 'Ap'}; {'VVp, VAp', 'Vp'}];
for i = 1:length(events_merge)
  bst_process('CallProcess', 'process_evt_merge', sFiles, [], 'evtnames', events_merge{i, 1}, 'newname', events_merge{i, 2});
end
%% ---------------------------

%% epoch
% prime trials
events_prime = 'Ap, Vp';
for i = 1:n
  bst_process('CallProcess', 'process_import_data_event', sFiles(i), [], 'subjectname', subjects{i}, ...
              'condition', '', 'eventname', events_prime, 'timewindow', [], 'epochtime', timewindow, ...
              'createcond', 1, 'ignoreshort', 0, 'usectfcomp', 0, 'usessp', 0, 'freq', [], 'baseline', []);  
end
% repetition trials
%% ---------------------------