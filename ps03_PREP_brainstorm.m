%% ---------------------------
%% [script name] ps03_PREP_brainstorm.m
%%
%% SCRIPT to do pre-processing by using brainstorm.
%%
%% By Shuai Wang, [date] 2021-04-18
%%
%% ---------------------------
%% Notes: - to keep the original events while merging events, process_evt_merge() was modified (see its codes at line 166).
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
ptoken = 'ses-01_task-RS_run-01_ieeg';     % raw data token
notch_filters = [50, 100, 150, 200, 250];  % notch frequencies in Hz
freq_highpass = 0.3;                       % 0.3 Hz recommended by AnneSo
freq_lowpass = 0;                          % 0 means disable
% manipulation of events
events.repetition = 'AAp, AAt, AVp, AVt, VAp, VAt, VVp, VVt';   % original eight RSE conditions
events.prime = 'Ap, Vp';                                        % prime conditions merged from the RSE conditions
events.check = 'Ap_long, Vp_long';                              % prime conditions with a long time-window
events.control = 'WV, PA, PV, SA, SV';                          % control conditions with the same baseline as the primes: visual words, pseudowords and scrambles
events.video = 'FA, FV, WA';                                    % videos and stimuli extracted from videos
events.pword = 'WPA, WPV';                                      % combined words/pseudowords in either visual or auditory modality
events_manipulation.merge_prime = [{'AAp, AVp', 'Ap'}; {'VVp, VAp', 'Vp'}];            % merge AAp and AVp as Ap (auditory prime trials), do the same for visual trials
events_manipulation.merge_check = [{'AAp, AVp', 'Ap_long'}; {'VVp, VAp', 'Vp_long'}];  % merge auditory/visual trials for sanicty check in MIA requested by Agnes
events_manipulation.merge_pword = [{'Ap, WA, PA', 'WPA'}; {'Vp, WV, PV', 'WPV'}];      % merge word and pseudoword together for sanicty check in MIA requested by Agnes
% epoch parameters
timewindow.repetition = [-0.3 0.7];        % epoch window for repetition trials
timewindow.prime = [-0.4 0.8];             % epoch window for prime trials
timewindow.check = [-0.5 1];               % epoch window for sanicty check in MIA requested by Agnes
timewindow.control = timewindow.prime;     % epoch window for control conditions
timewindow.video = [-0.4 1.2];             % epoch window for video-based stimuli
timewindow.pword = timewindow.prime;       % epoch window for combined words/pseudowords conditions
% baseline parameters
baseline_method = 'bl';                    % DC offset correction: x_std = x - u;
baseline.repetition = [-0.3 0];            % baseline window for the repetition trials
baseline.prime = [-0.4 0];                 % baseline window for the prime condition and other conditions
baseline.check = baseline.prime;
baseline.control = baseline.prime;
baseline.video = baseline.prime;
baseline.pword = baseline.prime;
%% ---------------------------

%% initialize the working file
% create protocol if not exist
if ~exist(bdir, 'dir')
  fprintf('Create a new protocol %s in the folder %s. \n', protocol, rdir);
  gui_brainstorm('CreateProtocol', protocol, 0, 0, rdir);  % create a database in the BIDS derivatives folder
end
% initialize the working file
ftmp = fullfile(ddir, 'ps03_PREP_working.mat');
% take a note
pdate = datetime('now','Format', 'yyyy-MM-dd''T''HH:mm:SS');  % the date of the present processing
pnote = 'Run the entire pre-processing pipeline for 10 subjects (from sub-01 to sub-10).';
fprintf('%s \nStart at %s. \n\n', pnote, pdate);
if exist(ftmp, 'file')
  load(ftmp);  % read up the working file
  commits(size(commits, 1) + 1, :) = {pnote, sprintf('%s', pdate)};  % add commit
else
  commits = {pnote, sprintf('%s', pdate)};  % add commit
  save(ftmp, '*dir', 'subjects*', 'n', 'ptoken', 'notch_filters', 'freq*', 'events*', 'timewindow', 'baseline*', 'commits')
end
%% ---------------------------

%% import raw data
if ~exist('pFiles', 'var')
  for i = 1:n
    subj = subjects{i};
    sidx = ismember(subjects_info.participant_id, subj);
    fraw = fullfile(wdir, subj, 'ses-01', 'ieeg', sprintf('%s_%s.%s', subj, ptoken, subjects_info.data_format{sidx}));  % raw SEEG recordings (maybe in different formats)
    fprintf('Import raw data file %s into the BST database. \n', fraw);
    % read up raw data
    bst_process('CallProcess', 'process_import_data_raw', [], [], 'subjectname', subj, 'datafile', {fraw, 'SEEG-ALL'}, ...
                'channelreplace', 1, 'channelalign', 0, 'evtmode', 'value');
  end
  % extract BST filenames for next processes
  pFiles = cellfun(@(x) fullfile(ddir, x, sprintf('@raw%s_%s/data_0raw_%s_%s.mat', x, ptoken, x, ptoken)), subjects, 'UniformOutput', 0);
  save(ftmp, 'pFiles', '-append');
end
%% ---------------------------

%% notch filters : (PSD), notch at 50, 100, 150, 200, 250Hz
if ~exist('sFiles', 'var') || ~isfield(sFiles, 'notch')
  % PSD (to check bad channels)
  bst_process('CallProcess', 'process_psd', pFiles, [], 'timewindow', [], 'win_length', 10, 'win_overlap', 50, ...
              'units', 'physical', 'sensortypes', 'SEEG', 'win_std', 0, ...
              'edit', struct('Comment', 'Power', 'TimeBands', [], 'Freqs', [], 'ClusterFuncTime', 'none', ...
              'Measure', 'power', 'Output', 'all', 'SaveKernel', 0));
  % notch filters (50, 100, 150, 200, 250Hz)
  sFiles.notch = bst_process('CallProcess', 'process_notch', pFiles, [], 'sensortypes', 'SEEG', 'freqlist', notch_filters, ...
                             'cutoffW', 2, 'useold', 0, 'read_all', 0);
  % update the working file
  save(ftmp, 'sFiles', '-append');
end
%% ---------------------------

%% import recovered/double-checked events
if ~isfield(sFiles, 'events')
  for i = 1:n
    subj = subjects{i};
    fdat = fullfile(ddir, sFiles.notch(i).FileName);
    fevt = fullfile(wdir, subj, 'ses-01', 'ieeg', sprintf('%s_%s_events-final.csv', subj, ptoken));  % the double-checked event file that has to be CSV.
    fprintf('Write the final version of events %s into the data %s. \n', fevt, fdat);
    % write events
    bst_process('CallProcess', 'process_evt_import', fdat, [], 'evtfile', {fevt, 'CSV-TIME'}, 'evtname', '', 'delete', 1);
  end
  sFiles.events = true;  % imported events
  save(ftmp, 'sFiles', '-append');  % update the working file
end
%% ---------------------------

%% band-pass
if ~isfield(sFiles, 'bandpass')
  fprintf('Apply band-pass filtering (%d - %d Hz) on the data. \n', freq_highpass, freq_lowpass);
  % band-pass filters
  sFiles.bandpass = bst_process('CallProcess', 'process_bandpass', sFiles.notch, [], 'sensortypes', 'SEEG', ...
                                'highpass', freq_highpass, 'lowpass', freq_lowpass, 'tranband', 0, 'attenuation', 'strict', ...  % 60dB
                                'ver', '2019', 'mirror', 0, 'read_all', 0);
  % update the working file
  save(ftmp, 'sFiles', '-append');
end
%% ---------------------------

%% manipulate events
events_manip = fieldnames(events_manipulation);  % extract manipulation labels
if ~isfield(sFiles, 'events_created')
  sFiles.events_created = events_manip;
else
  events_manip = events_manip(~ismember(events_manip, sFiles.events_created));  % derive new manips
  sFiles.events_created = [sFiles.events_created; events_manip];                % update sFile
end
if ~isempty(events_manip)
  % do each manipulation
  for i = 1:length(events_manip)
    imanip = events_manip{i};
    evts_manip = events_manipulation.(imanip);
    % merge trials
    for j = 1:length(evts_manip)
      fprintf('Create a new event %s. \n', evts_manip{j, 2});
      bst_process('CallProcess', 'process_evt_merge', sFiles.bandpass, [], 'evtnames', evts_manip{j, 1}, 'newname', evts_manip{j, 2});
    end
  end
  % update the working file
  save(ftmp, 'sFiles', '-append');  
end
%% ---------------------------

%% epoch
events_epoch = fieldnames(events);
if ~isfield(sFiles, 'trials')
  sFiles.events_epoch = events_epoch;  % update sFiles
else
  events_epoch = events_epoch(~ismember(events_epoch, fieldnames(sFiles.trials)));
  sFiles.events_epoch = [sFiles.events_epoch; events_epoch];
end
if ~isempty(events_epoch)
  % do epoch for each event
  for i = 1:length(events_epoch)
    ievt = events_epoch{i};
    fprintf('Epoch the event %s. \n', ievt);
    sFiles.trials.(ievt) = bst_process('CallProcess', 'process_import_data_event', sFiles.bandpass, [], 'subjectname', '', ...
                                       'condition', '', 'eventname', events.(ievt), 'timewindow', [], 'epochtime', timewindow.(ievt), ...
                                       'createcond', 1, 'ignoreshort', 0, 'usectfcomp', 0, 'usessp', 0, 'freq', [], 'baseline', []);
  end
  % update the working file
  save(ftmp, 'sFiles', '-append');    
end                           
%% ---------------------------

%% monopolar and bipolar montages
events_montage = fieldnames(events);
if isfield(sFiles, 'trialdata')
  events_montage = events_montage(~ismember(events_montage, fieldnames(sFiles.trialdata)));
end
if ~isempty(events_montage)
  for i = 1:n
    subj = subjects{i};
    mon1 = sprintf('%s: SEEG (orig)[tmp]', subj);       % monopolar
    mon2 = sprintf('%s: SEEG (bipolar 2)[tmp]', subj);  % bipolar2
    sFiles.trialdata(i).subject = subj;
    sFiles.monopolar(i).subject = subj;
    sFiles.bipolar(i).subject = subj;
    % do montage for each event
    for j = 1:length(events_montage)
      ievt = events_montage{j};
      % group all conditions
      fprintf('Apply montage %s and %s for the event %s. \n', mon1, mon2, ievt);
      trials_evt = sFiles.trials.(ievt);  % get trials of this event
      trials_evt = {trials_evt(ismember({trials_evt.SubjectName}, subj)).FileName};
      sFiles.trialdata(i).(ievt) = cellfun(@(x) fullfile(ddir, x), trials_evt, 'UniformOutput', 0);  
      % apply monopolar montage (n)
      sFiles.monopolar(i).(ievt) = bst_process('CallProcess', 'process_montage_apply', sFiles.trialdata(i).(ievt), [], 'montage', mon1, 'createchan', 1);
      % apply bipolar2 montage (m)
      sFiles.bipolar(i).(ievt) = bst_process('CallProcess', 'process_montage_apply', sFiles.trialdata(i).(ievt), [], 'montage', mon2, 'createchan', 1);
    end
  end 
  % update the working file
  save(ftmp, 'sFiles', '-append'); 
end
%% ---------------------------

%% baseline normalization
events_norm = fieldnames(events);
if isfield(sFiles, 'monopolar_norm')
  events_norm = events_norm(~ismember(events_norm, fieldnames(sFiles.monopolar_norm)));
end
if ~isempty(events_norm)
  for i = 1:n
    subj = subjects{i};
    % read working filenames for this subject
    mono_subj = sFiles.monopolar(i);
    bipo_subj = sFiles.bipolar(i);
    % do baseline normalization for each event
    for j = 1:length(events_norm)
      ievt = events_norm{j};
      fprintf('Apply baseline normalization for the event %s. \n', ievt);
      % baseline normalization for monopolar montage; these files would be used for ERP analysis
      mono_evts = cellfun(@(x) fullfile(ddir, x), {mono_subj.(ievt).FileName}, 'UniformOutput', 0);
      sFiles.monopolar_norm(i).(ievt) = bst_process('CallProcess', 'process_baseline_norm', mono_evts, [], 'baseline', baseline.(ievt), ...
                                                    'sensortypes', 'SEEG', 'method', baseline_method, 'overwrite', 0);  
      % baseline normalization for bipolar montage
      bipo_evts = cellfun(@(x) fullfile(ddir, x), {bipo_subj.(ievt).FileName}, 'UniformOutput', 0);
      sFiles.bipolar_norm(i).(ievt) = bst_process('CallProcess', 'process_baseline_norm', bipo_evts, [], 'baseline', baseline.(ievt), ...
                                                  'sensortypes', 'SEEG', 'method', baseline_method, 'overwrite', 0); 
    end
  end
  % update the working file
  save(ftmp, 'sFiles', '-append');   
end
%% ---------------------------