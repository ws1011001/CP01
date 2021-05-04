%% ---------------------------
%% [script name] ps03_PREP_brainstorm.m
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
wdir = fullfile(mdir, 'SEEG_LectureVWFA');                     % main working folder 
ddir = fullfile(wdir, 'derivatives');                          % derivatives in the BIDS structure
bdir = fullfile(ddir, 'brainstorm_SEEG_LectureVWFA', 'data');  % path to brainsotrm database
% read the subjects list and information
fid = fopen(fullfile(mdir, 'CP01_subjects.txt'));
subjects = textscan(fid, '%s');
fclose(fid);
subjects = subjects{1};  % the subjects list
n = length(subjects);
subjects_info = readtable(fullfile(wdir, 'participants.tsv'), 'FileType', 'text');
% set processing parameters
notch_fliters = [50, 100, 150, 200, 250];
freq_highpass = 0.3;  % 0.3 Hz recommended by AnneSo
freq_lowpass = 40;    % 40 Hz low-pass for ERP analysis
timewindow.prime = [-0.4, 0.8];  % epoch window for prime trials
%% ---------------------------

%% import raw data
ptoken = 'ses-01_task-RS_run-01_ieeg';  % raw data token
for i = 1:n
  subj = subjects{i};
  fraw = fullfile(wdir, subj, 'ses-01', 'ieeg', sprintf('%s_%s.%s', subj, ptoken, subjects_info.data_format{i}));
  bst_process('CallProcess', 'process_import_data_raw', [], [], 'subjectname', subj, 'datafile', {fraw, 'SEEG-ALL'}, ...
              'channelreplace', 1, 'channelalign', 0, 'evtmode', 'value');
end
%% ---------------------------

%% basic pre-processing : (PSD), notch, band-pass
% PSD (to check bad channels)
sFiles = cellfun(@(x) fullfile(bdir, x, sprintf('@raw%s_%s/data_0raw_%s_%s.mat', x, ptoken, x, ptoken)), subjects, 'UniformOutput', 0);
bst_process('CallProcess', 'process_psd', sFiles, [], 'timewindow', [], 'win_length', 10, 'win_overlap', 50, ...
            'units', 'physical', 'sensortypes', 'SEEG', 'win_std', 0, ...
            'edit', struct('Comment', 'Power', 'TimeBands', [], 'Freqs', [], 'ClusterFuncTime', 'none', ...
            'Measure', 'power', 'Output', 'all', 'SaveKernel', 0));
% notch filters (50, 100, 150, 200, 250Hz)
sFiles.notch = bst_process('CallProcess', 'process_notch', sFiles, [], 'sensortypes', 'SEEG', 'freqlist', notch_fliters, ...
                     'cutoffW', 2, 'useold', 0, 'read_all', 0);
% high-pass filter - 0.3 Hz
sFiles.highpass = bst_process('CallProcess', 'process_bandpass', sFiles.notch, [], 'sensortypes', 'SEEG', ...
                              'highpass', freq_highpass, 'lowpass', 0, 'tranband', 0, 'attenuation', 'strict', ...  % 60dB
                              'ver', '2019', 'mirror', 0, 'read_all', 0);
% band-pass filters - 0.3 Hz to 40 Hz (for ERPs only)
sFiles.bandpass = bst_process('CallProcess', 'process_bandpass', sFiles.notch, [], 'sensortypes', 'SEEG', ...
                              'highpass', freq_highpass, 'lowpass', freq_lowpass, 'tranband', 0, 'attenuation', 'strict', ...  % 60dB
                              'ver', '2019', 'mirror', 0, 'read_all', 0);
%% ---------------------------

%% manipulate events
% IMPORT recovered/double-checked events%
% read filenames
ptoken = 'ses-01_task-RS_run-01_ieeg_notch_high';  % data after notch and high-pass
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
              'condition', '', 'eventname', events_prime, 'timewindow', [], 'epochtime', timewindow.prime, ...
              'createcond', 1, 'ignoreshort', 0, 'usectfcomp', 0, 'usessp', 0, 'freq', [], 'baseline', []);  
end
% repetition trials
%% ---------------------------