%% ---------------------------
%% [script name] ps03_PREP_brainstorm.m
%%
%% SCRIPT to do pre-processing by using brainstorm.
%%
%% By Shuai Wang, [date] 2021-04-18
%%
%% ---------------------------
%% Notes: - to keep the original events while merging events, process_evt_merge() was modified (see its codes).
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
ptoken = 'ses-01_task-RS_run-01_ieeg';     % raw data token
notch_filters = [50, 100, 150, 200, 250];  % notch frequencies in Hz
freq_highpass = 0.3;                       % 0.3 Hz recommended by AnneSo
freq_lowpass = 0;                          % 0 means disable
timewindow.check = [-0.5 1];               % epoch window for sanicty check in MIA requested by Agnes
timewindow.prime = [-0.4 0.8];             % epoch window for prime trials
timewindow.repetition = [-0.3 0.7];        % epoch window for repetition trials
timewindow.video = [-0.4 1.2];             % epoch window for video-based stimuli
baseline.prime = [-0.4 0];                 % baseline window for prime trials
baseline.repetition = [-0.3 0];            % baseline window for repetition trials
baseline_method = 'bl';                    % DC offset correction: x_std = x - u;
% manipulation of conditions
events.repetitions = 'AAp, AAt, AVp, AVt, VAp, VAt, VVp, VVt';   % original eight RSE conditions
events.merge_primes = [{'AAp, AVp', 'Ap'}; {'VVp, VAp', 'Vp'}];  % merge AAp and AVp as Ap (auditory prime trials), do the same for visual trials
events.primes = 'Ap, Vp';
events.merge_checks = [{'AAp, AVp', 'Ap_long'}; {'VVp, VAp', 'Vp_long'}];  % merge auditory/visual trials for sanicty check in MIA requested by Agnes
events.checks = 'Ap_long, Vp_long';
events.controls = 'WV, PA, PV, SA, SV';  % control conditions sharing the same baseline with primes: visual words, pseudowords and scrambles
events.videos = 'FA, FV, WA';  % videos and stimuli extracted from videos
%% ---------------------------

%% import raw data
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
sFiles.notch = bst_process('CallProcess', 'process_notch', pFiles, [], 'sensortypes', 'SEEG', 'freqlist', notch_filters, ...
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
% merge prime trials with long segment (-500 ~ 1000 ms) for Agnes
for i = 1:length(events.merge_checks)
  bst_process('CallProcess', 'process_evt_merge', sFiles.bandpass, [], 'evtnames', events.merge_checks{i, 1}, 'newname', events.merge_checks{i, 2});
end
%% ---------------------------

%% epoch
% check trials
trials.checks = bst_process('CallProcess', 'process_import_data_event', sFiles.bandpass, [], 'subjectname', '', ...
                            'condition', '', 'eventname', events.checks, 'timewindow', [], 'epochtime', timewindow.check, ...
                            'createcond', 1, 'ignoreshort', 0, 'usectfcomp', 0, 'usessp', 0, 'freq', [], 'baseline', []);
% prime trials 
trials.primes = bst_process('CallProcess', 'process_import_data_event', sFiles.bandpass, [], 'subjectname', '', ...
                            'condition', '', 'eventname', events.primes, 'timewindow', [], 'epochtime', timewindow.prime, ...
                            'createcond', 1, 'ignoreshort', 0, 'usectfcomp', 0, 'usessp', 0, 'freq', [], 'baseline', []);              
% repetition trials
trials.repetitions = bst_process('CallProcess', 'process_import_data_event', sFiles.bandpass, [], 'subjectname', '', ...
                                 'condition', '', 'eventname', events.repetitions, 'timewindow', [], 'epochtime', timewindow.repetition, ...
                                 'createcond', 1, 'ignoreshort', 0, 'usectfcomp', 0, 'usessp', 0, 'freq', [], 'baseline', []);
% control trials 
trials.controls = bst_process('CallProcess', 'process_import_data_event', sFiles.bandpass, [], 'subjectname', '', ...
                              'condition', '', 'eventname', events.controls, 'timewindow', [], 'epochtime', timewindow.prime, ...
                              'createcond', 1, 'ignoreshort', 0, 'usectfcomp', 0, 'usessp', 0, 'freq', [], 'baseline', []);
% video-based trials 
trials.videos = bst_process('CallProcess', 'process_import_data_event', sFiles.bandpass, [], 'subjectname', '', ...
                            'condition', '', 'eventname', events.videos, 'timewindow', [], 'epochtime', timewindow.video, ...
                            'createcond', 1, 'ignoreshort', 0, 'usectfcomp', 0, 'usessp', 0, 'freq', [], 'baseline', []);                            
% save working filepaths and parameters
ftmp = fullfile(bdir, sprintf('ps03_PREP_%s.mat', datetime('now','Format', 'yyyy-MM-dd''T''HH:mm:SS')));
save(ftmp, 'pFiles', 'sFiles', 'trials', 'subjects', 'n', 'notch_filters', 'freq_highpass', 'freq_lowpass', 'timewindow', 'baseline', 'events');
%% ---------------------------

%% monopolar and bipolar montages
for i =1:n
  subj = subjects{i};
  mon1 = sprintf('%s: SEEG (orig)[tmp]', subj);       % monopolar
  mon2 = sprintf('%s: SEEG (bipolar 2)[tmp]', subj);  % bipolar2
  % group all conditions
  fprintf('Apply montage %s and %s \n', mon1, mon2);
  tFiles.checks = {trials.checks(ismember({trials.checks.SubjectName}, subj)).FileName};
  tFiles.checks = cellfun(@(x) fullfile(bdir, x), tFiles.checks, 'UniformOutput', 0);  
  tFiles.primes = {trials.primes(ismember({trials.primes.SubjectName}, subj)).FileName};
  tFiles.primes = cellfun(@(x) fullfile(bdir, x), tFiles.primes, 'UniformOutput', 0);
  tFiles.repetitions = {trials.repetitions(ismember({trials.repetitions.SubjectName}, subj)).FileName};
  tFiles.repetitions = cellfun(@(x) fullfile(bdir, x), tFiles.repetitions, 'UniformOutput', 0);
  tFiles.controls = {trials.controls(ismember({trials.controls.SubjectName}, subj)).FileName};
  tFiles.controls = cellfun(@(x) fullfile(bdir, x), tFiles.controls, 'UniformOutput', 0); 
  tFiles.videos = {trials.videos(ismember({trials.videos.SubjectName}, subj)).FileName};
  tFiles.videos = cellfun(@(x) fullfile(bdir, x), tFiles.videos, 'UniformOutput', 0);   
  % apply monopolar montage (n)
  nFiles.checks = bst_process('CallProcess', 'process_montage_apply', tFiles.checks, [], 'montage', mon1, 'createchan', 1);
  nFiles.primes = bst_process('CallProcess', 'process_montage_apply', tFiles.primes, [], 'montage', mon1, 'createchan', 1);
  nFiles.repetitions = bst_process('CallProcess', 'process_montage_apply', tFiles.repetitions, [], 'montage', mon1, 'createchan', 1);  
  nFiles.controls = bst_process('CallProcess', 'process_montage_apply', tFiles.controls, [], 'montage', mon1, 'createchan', 1);
  nFiles.videos = bst_process('CallProcess', 'process_montage_apply', tFiles.videos, [], 'montage', mon1, 'createchan', 1);
  % apply bipolar2 montage (m)
  mFiles.checks = bst_process('CallProcess', 'process_montage_apply', tFiles.checks, [], 'montage', mon2, 'createchan', 1);
  mFiles.primes = bst_process('CallProcess', 'process_montage_apply', tFiles.primes, [], 'montage', mon2, 'createchan', 1);
  mFiles.repetitions = bst_process('CallProcess', 'process_montage_apply', tFiles.repetitions, [], 'montage', mon2, 'createchan', 1);
  mFiles.controls = bst_process('CallProcess', 'process_montage_apply', tFiles.controls, [], 'montage', mon2, 'createchan', 1);
  mFiles.videos = bst_process('CallProcess', 'process_montage_apply', tFiles.videos, [], 'montage', mon2, 'createchan', 1);
  % save working filenames
  ftmp = fullfile(bdir, subj, sprintf('%s_%s_working-filenames.mat', subj, ptoken));
  save(ftmp, 'tFiles', 'nFiles', 'mFiles', 'mon1', 'mon2', 'timewindow', 'baseline', 'events');
end 
%% ---------------------------

%% baseline normalization
for i =1:n
  subj = subjects{i};
  % read working filenames
  ftmp = fullfile(bdir, subj, sprintf('%s_%s_working-filenames.mat', subj, ptoken));
  load(ftmp, 'nFiles', 'mFiles');
  % baseline normalization for monopolar montage (o) 
  nFiles.primes = cellfun(@(x) fullfile(bdir, x), {nFiles.primes.FileName}, 'UniformOutput', 0);
  nFiles.repetitions = cellfun(@(x) fullfile(bdir, x), {nFiles.repetitions.FileName}, 'UniformOutput', 0);
  nFiles.controls = cellfun(@(x) fullfile(bdir, x), {nFiles.controls.FileName}, 'UniformOutput', 0);
  nFiles.videos = cellfun(@(x) fullfile(bdir, x), {nFiles.videos.FileName}, 'UniformOutput', 0);
  oFiles.primes = bst_process('CallProcess', 'process_baseline_norm', nFiles.primes, [], 'baseline', baseline.prime, ...
                              'sensortypes', 'SEEG', 'method', baseline_method, 'overwrite', 0);  
  oFiles.repetitions = bst_process('CallProcess', 'process_baseline_norm', nFiles.repetitions, [], 'baseline', baseline.repetition, ...
                                   'sensortypes', 'SEEG', 'method', baseline_method, 'overwrite', 0); 
  oFiles.controls = bst_process('CallProcess', 'process_baseline_norm', nFiles.controls, [], 'baseline', baseline.prime, ...
                                'sensortypes', 'SEEG', 'method', baseline_method, 'overwrite', 0);
  oFiles.videos = bst_process('CallProcess', 'process_baseline_norm', nFiles.videos, [], 'baseline', baseline.prime, ...
                              'sensortypes', 'SEEG', 'method', baseline_method, 'overwrite', 0);                            
  % baseline normalization for bipolar montage (b) 
  mFiles.primes = cellfun(@(x) fullfile(bdir, x), {mFiles.primes.FileName}, 'UniformOutput', 0);
  mFiles.repetitions = cellfun(@(x) fullfile(bdir, x), {mFiles.repetitions.FileName}, 'UniformOutput', 0);
  mFiles.controls = cellfun(@(x) fullfile(bdir, x), {mFiles.controls.FileName}, 'UniformOutput', 0);
  mFiles.videos = cellfun(@(x) fullfile(bdir, x), {mFiles.videos.FileName}, 'UniformOutput', 0);
  bFiles.primes = bst_process('CallProcess', 'process_baseline_norm', mFiles.primes, [], 'baseline', baseline.prime, ...
                              'sensortypes', 'SEEG', 'method', baseline_method, 'overwrite', 0);  
  bFiles.repetitions = bst_process('CallProcess', 'process_baseline_norm', mFiles.repetitions, [], 'baseline', baseline.repetition, ...
                                   'sensortypes', 'SEEG', 'method', baseline_method, 'overwrite', 0);
  bFiles.controls = bst_process('CallProcess', 'process_baseline_norm', mFiles.controls, [], 'baseline', baseline.prime, ...
                                'sensortypes', 'SEEG', 'method', baseline_method, 'overwrite', 0);
  bFiles.videos = bst_process('CallProcess', 'process_baseline_norm', mFiles.videos, [], 'baseline', baseline.prime, ...
                              'sensortypes', 'SEEG', 'method', baseline_method, 'overwrite', 0);                            
  % save working filenames
  save(ftmp, 'nFiles', 'oFiles', 'mFiles', 'bFiles', '-append');  % nFiles (monopolar) would be used for ERP analysis
end
%% ---------------------------