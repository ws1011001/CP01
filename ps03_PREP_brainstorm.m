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
freq_lowpass = 300;   % to avoid aliasing effects
timewindow.prime = [-0.4 0.8];       % epoch window for prime trials
timewindow.repetition = [-0.3 0.7];  % epoch window for repetition trials
baseline.prime = [-0.4 0.01];        % baseline window for prime trials
baseline.repetition = [-0.3 0.01];   % baseline window for repetition trials
events.merge_primes = [{'AAp, AVp', 'Ap'}; {'VVp, VAp', 'Vp'}];  % auditory and visual prime trials
events.primes = 'Ap, Vp';
events.repetitions = 'AAp, AAt, AVp, AVt, VAp, VAt, VVp, VVt';
%% ---------------------------

%% import raw data
ptoken = 'ses-01_task-RS_run-01_ieeg';  % raw data token
for i = 1:n
  subj = subjects{i};
  sidx = ismember(subjects_info.participant_id, subj);
  fraw = fullfile(wdir, subj, 'ses-01', 'ieeg', sprintf('%s_%s.%s', subj, ptoken, subjects_info.data_format{sidx}));
  fprintf('Import raw data file %s \n', fraw);
  bst_process('CallProcess', 'process_import_data_raw', [], [], 'subjectname', subj, 'datafile', {fraw, 'SEEG-ALL'}, ...
              'channelreplace', 1, 'channelalign', 0, 'evtmode', 'value');
end
%% ---------------------------

%% notch filters : (PSD), notch at 50, 100, 150, 200, 250Hz
% PSD (to check bad channels)
pFiles = cellfun(@(x) fullfile(bdir, x, sprintf('@raw%s_%s/data_0raw_%s_%s.mat', x, ptoken, x, ptoken)), subjects, 'UniformOutput', 0);
bst_process('CallProcess', 'process_psd', pFiles, [], 'timewindow', [], 'win_length', 10, 'win_overlap', 50, ...
            'units', 'physical', 'sensortypes', 'SEEG', 'win_std', 0, ...
            'edit', struct('Comment', 'Power', 'TimeBands', [], 'Freqs', [], 'ClusterFuncTime', 'none', ...
            'Measure', 'power', 'Output', 'all', 'SaveKernel', 0));
% notch filters (50, 100, 150, 200, 250Hz)
sFiles.notch = bst_process('CallProcess', 'process_notch', pFiles, [], 'sensortypes', 'SEEG', 'freqlist', notch_fliters, ...
                     'cutoffW', 2, 'useold', 0, 'read_all', 0);
%% ---------------------------

%% import recovered/double-checked events
for i = 1:n
  subj = subjects{i};
  fdat = fullfile(bdir, sFiles.notch(i).FileName);
  fevt = fullfile(wdir, subj, 'ses-01', 'ieeg', sprintf('%s_%s_events-final.csv', subj, ptoken));
  fprintf('Write the final version of events %s to the data %s \n', fevt, fdat);
  bst_process('CallProcess', 'process_evt_import', fdat, [], 'evtfile', {fevt, 'CSV-TIME'}, 'evtname', '', 'delete', 1);
end
%% ---------------------------

%% band-pass
fprintf('Apply band-pass filtering (%d - %d Hz) on the data \n', freq_highpass, freq_lowpass);
% band-pass filters
sFiles.bandpass = bst_process('CallProcess', 'process_bandpass', sFiles.notch, [], 'sensortypes', 'SEEG', ...
                              'highpass', freq_highpass, 'lowpass', freq_lowpass, 'tranband', 0, 'attenuation', 'strict', ...  % 60dB
                              'ver', '2019', 'mirror', 0, 'read_all', 0);
%% ---------------------------

%% manipulate events
% merge prime trials
for i = 1:length(events.merge_primes)
  bst_process('CallProcess', 'process_evt_merge', sFiles.bandpass, [], 'evtnames', events.merge_primes{i, 1}, 'newname', events.merge_primes{i, 2});
end
%% ---------------------------

%% epoch
% prime trials 
trials.primes = bst_process('CallProcess', 'process_import_data_event', sFiles.bandpass, [], 'subjectname', '', ...
                            'condition', '', 'eventname', events.primes, 'timewindow', [], 'epochtime', timewindow.prime, ...
                            'createcond', 1, 'ignoreshort', 0, 'usectfcomp', 0, 'usessp', 0, 'freq', [], 'baseline', []);              
% repetition trials
trials.repetitions = bst_process('CallProcess', 'process_import_data_event', sFiles.bandpass, [], 'subjectname', '', ...
                                 'condition', '', 'eventname', events.repetitions, 'timewindow', [], 'epochtime', timewindow.repetition, ...
                                 'createcond', 1, 'ignoreshort', 0, 'usectfcomp', 0, 'usessp', 0, 'freq', [], 'baseline', []);
%% ---------------------------

%% bipolar montage
for i =1:n
  subj = subjects{i};
  mont = sprintf('%s: SEEG (bipolar 2)[tmp]', subj);
  ftmp = fullfile(bdir, subj, sprintf('%s_%s_working-filenames.mat', subj, ptoken));
  % group all conditions
  fprintf('Apply montage %s \n', mont);
  tFiles.primes = {trials.primes(ismember({trials.primes.SubjectName}, subj)).FileName};
  tFiles.primes = cellfun(@(x) fullfile(bdir, x), tFiles.primes, 'UniformOutput', 0);
  tFiles.repetitions = {trials.repetitions(ismember({trials.repetitions.SubjectName}, subj)).FileName};
  tFiles.repetitions = cellfun(@(x) fullfile(bdir, x), tFiles.repetitions, 'UniformOutput', 0);
  % apply montage
  mFiles.primes = bst_process('CallProcess', 'process_montage_apply', tFiles.primes, [], 'montage', mont, 'createchan', 1);
  mFiles.repetitions = bst_process('CallProcess', 'process_montage_apply', tFiles.repetitions, [], 'montage', mont, 'createchan', 1);
  % save working filenames
  save(ftmp, 'tFiles', 'mFiles', 'mont', 'timewindow', 'baseline', 'events');
end 
%% ---------------------------

%% baseline normalization
for i =1:n
  subj = subjects{i};
  % read working filenames
  ftmp = fullfile(bdir, subj, sprintf('%s_%s_working-filenames.mat', subj, ptoken));
  load(ftmp, 'mFiles');
  mFiles.primes = cellfun(@(x) fullfile(bdir, x), {mFiles.primes.FileName}, 'UniformOutput', 0);
  mFiles.repetitions = cellfun(@(x) fullfile(bdir, x), {mFiles.repetitions.FileName}, 'UniformOutput', 0);
  % baseline normalization
  bFiles.primes = bst_process('CallProcess', 'process_baseline_norm', mFiles.primes, [], 'baseline', baseline.prime, ...
                              'sensortypes', 'SEEG', 'method', 'zscore', ...  % Z-score transformation: x_std = (x - mu) / sigma
                              'overwrite', 0);  
  bFiles.repetitions = bst_process('CallProcess', 'process_baseline_norm', mFiles.repetitions, [], 'baseline', baseline.repetition, ...
                                   'sensortypes', 'SEEG', 'method', 'zscore', ...  % Z-score transformation: x_std = (x - mu) / sigma
                                   'overwrite', 0);
  % save working filenames
  save(ftmp, 'bFiles', '-append');  % bFiles would be used for next analyses
end
%% ---------------------------