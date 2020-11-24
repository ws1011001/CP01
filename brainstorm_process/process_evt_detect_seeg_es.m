function varargout = process_evt_detect_seeg_es( varargin )
% PROCESS_EVT_DETECT_SEEG_ES: Epileptic spikes detection for a group of recordings file
%
% USAGE:  OutputFiles = process_evt_detect_seeg_es('Run', sProcess, sInputs)

% @=============================================================================
% This function aims to reduce potential ictal artifacts in SEEG data according to the following study:
% Hirshorn, E. a., Li, Y., Ward, M. J., Richardson, R. M., Fiez, J. a., & Ghuman, A. S. (2016). Decoding and disrupting
% left midfusiform gyrus activity during word reading. PNAS, 113(29), 201604126. https://doi.org/10.1073/pnas.1604126113
% 
% Copyright (c)2000-2020 University of Southern California & McGill University
% This software is distributed under the terms of the GNU General Public License
% as published by the Free Software Foundation. Further details on the GPLv3
% license can be found at http://www.gnu.org/copyleft/gpl.html.
% 
% FOR RESEARCH PURPOSES ONLY. THE SOFTWARE IS PROVIDED "AS IS," AND THE
% UNIVERSITY OF SOUTHERN CALIFORNIA AND ITS COLLABORATORS DO NOT MAKE ANY
% WARRANTY, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO WARRANTIES OF
% MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE, NOR DO THEY ASSUME ANY
% LIABILITY OR RESPONSIBILITY FOR THE USE OF THIS SOFTWARE.
%
% For more information type "brainstorm license" at command prompt.
% =============================================================================@
%
% Authors: Shuai Wang, Anne-Sophie Dubarry 2020

eval(macro_method);
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
    % Description the process
    sProcess.Comment     = 'Detect potential epileptic artifacts (Hirshorn 2016)';
    sProcess.FileTag     = 'ES';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = 'Events';
    sProcess.Index       = 46;  % the same Index as PROCESS_EVT_DETECT_BADSEGMENT
    sProcess.Description = 'See SI in Hirshorn et al., 2016, PNAS.';
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'data'};
    sProcess.OutputTypes = {'data'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    % Sensor types
    sProcess.options.sensortypes.Comment = 'Sensor types or names (empty=all): ';
    sProcess.options.sensortypes.Type    = 'text';
    sProcess.options.sensortypes.Value   = 'SEEG';    
    % Threshold0 - SDs from the average amplitude (across trials)
    sProcess.options.threshold0.Comment = 'Peak X SDs away from the Mean (across trials): ';
    sProcess.options.threshold0.Type    = 'value';
    sProcess.options.threshold0.Value   = {5, 'SDs', 0};    
    % Threshold1 - SDs from the average amplitude (within trial)
    sProcess.options.threshold1.Comment = '[deprecated] Peak X SDs away from the Mean (within trial): ';
    sProcess.options.threshold1.Type    = 'value';
    sProcess.options.threshold1.Value   = {1000, 'SDs', 0};
    % Threshold2 - absolute upper limit
    sProcess.options.threshold2.Comment = 'Peak amplitude (abs.) threshold: ';
    sProcess.options.threshold2.Type    = 'value';
    sProcess.options.threshold2.Value   = {350, 'uV', []};    
    % Threshold3 - variation between consecutive time points
    sProcess.options.threshold3.Comment = 'Variation between consecutive points (abs.): ';
    sProcess.options.threshold3.Type    = 'value';
    sProcess.options.threshold3.Value   = {25, 'uV', []}; 
    % Threshold4 - bad channel threshold : If a channel led to more than X% of the trials being rejected, this channel was instead rejected.
    sProcess.options.threshold4.Comment = 'Examine channels that led to more than X% trials being rejected: ';
    sProcess.options.threshold4.Type    = 'value';
    sProcess.options.threshold4.Value   = {20, '%', []};     
    % Validation with visually detected BAD events
    % File selection options
    SelectOptions = {...
        '', ...                               % Filename
        '', ...                               % FileFormat
        'open', ...                           % Dialog type: {open,save}
        'Import events...', ...               % Window title
        'ImportData', ...                     % LastUsedDir: {ImportData,ImportChannel,ImportAnat,ExportChannel,ExportData,ExportAnat,ExportProtocol,ExportImage,ExportScript}
        'single', ...                         % Selection mode: {single,multiple}
        'files', ...                          % Selection mode: {files,dirs,files_and_dirs}
        bst_get('FileFilters', 'events'), ... % Get all the available file formats
        'EventsIn'};                          % DefaultFormats: {ChannelIn,DataIn,DipolesIn,EventsIn,MriIn,NoiseCovIn,ResultsIn,SspIn,SurfaceIn,TimefreqIn
    sProcess.options.badevtfile.Comment = 'Event file with BAD segments:';
    sProcess.options.badevtfile.Type    = 'filename';
    sProcess.options.badevtfile.Value   = SelectOptions;
    % Enable classification
    sProcess.options.ismarkbadtrials.Comment = 'Mark bad trials';
    sProcess.options.ismarkbadtrials.Type    = 'checkbox';
    sProcess.options.ismarkbadtrials.Value   = 1;

end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>
    OutputFiles = {};
    % ===== GET OPTIONS ===== 
    % define channels
    SensorTypes = sProcess.options.sensortypes.Value;
    % get thresholds
    ThresSDat   = sProcess.options.threshold0.Value{1};
    ThresSDwt   = sProcess.options.threshold1.Value{1};
    ThresPeak   = sProcess.options.threshold2.Value{1} / 10e5;  % convert microvolts to volts
    ThresVars   = sProcess.options.threshold3.Value{1} / 10e5;  % convert microvolts to volts  
    ThresChan   = sProcess.options.threshold4.Value{1};  
    % read BAD_suspected events
    EventFile  = sProcess.options.badevtfile.Value{1};  % get filenames to import
    if exist(EventFile,'file')
      %FileFormat = sProcess.options.badevtfile.Value{2};  % only support .csv for now...
      [EventType, EventTime, EventBadTime] = eventscsv(EventFile);      
    else
      EventBadTime = [];
    end    
    % if mark bad trials or only output reports
    ismarkbadtrials = sProcess.options.ismarkbadtrials.Value;

    fprintf('------------------------------------------------------------------------------------------------------\n');
    fprintf('The following parameters were used to detect artifacts for %d trials: %d SDs (across trials), %d SDs (within trial), %d uV Peak Amplitude, %d uV consecutive variation.\n',...
      length(sInputs), ThresSDat, ThresSDwt, sProcess.options.threshold2.Value{1} , sProcess.options.threshold3.Value{1});
    
    % Get current progressbar position
    progressPos = bst_progress('get');
    
    % ===== INITIALIZE VECTORS =====    
    % select good channels - all trials have the same configuration of channels
    ChannelConfig = in_bst_channel(sInputs(1).ChannelFile);
    ChannelSelected = channel_find(ChannelConfig.Channel, SensorTypes);
    ChannelNames = {ChannelConfig.Channel(:).Name}';
    ChannelMarks = in_bst_data(sInputs(1).FileName, 'ChannelFlag');
    ChannelGood = ChannelMarks.ChannelFlag == 1;
    ChannelGoodNames = ChannelNames(ChannelSelected(ChannelGood(ChannelSelected)));
    
    % initialize the trial rejection vectors
    iOk = false(1, length(sInputs));
    ChannelTrialsMeans = zeros(sum(ChannelGood), length(sInputs));
    ChannelTrialsOutSDat = zeros(sum(ChannelGood), length(sInputs));
    ChannelTrialsOutSDwt = zeros(sum(ChannelGood), length(sInputs));
    ChannelTrialsOutPeak = zeros(sum(ChannelGood), length(sInputs));
    ChannelTrialsOutVars = zeros(sum(ChannelGood), length(sInputs));
    
    % ===== For each trial =====
    for iFile = 1:length(sInputs)
        % Progress bar
        bst_progress('text', 'Reading trial to process...');
        % ===== GET DATA =====
        % Load epochs        
        DataES = in_bst_data(sInputs(iFile).FileName, 'F', 'Time', 'ChannelFlag', 'Comment');        
        DataES.ChannelGood = ChannelGoodNames;
        DataES.FSelected = DataES.F(ChannelSelected(ChannelGood(ChannelSelected)), :);
        [DataES.TrialType, DataES.TrialIndex]=commentsplit(DataES.Comment);
        % check overlap with the bad segments
        if ~isempty(EventBadTime)
          [DataES.TrialEpoch, DataES.BadOverlaps] = badvalidate(DataES.TrialType, DataES.TrialIndex, DataES.Time, EventType, EventTime, EventBadTime);
        else
          DataES.BadOverlaps = [];
        end
        
        % ===== automatic detection Hirshorn 2016 =====
        DataES.TimeN = length(DataES.Time);
        [DataES.FPeaks, DataES.FPeaksIndex] = max(abs(DataES.FSelected), [], 2);
        % calculate the average and std for each channel
        DataES.FMeans = mean(DataES.FSelected, 2);
        DataES.FStds = std(DataES.FSelected, 0, 2);
        % peak amplitude x SDs away from the Mean across the rest of the trial
        DataES.FPeaksOut = DataES.FSelected(logical(repmat(1:DataES.TimeN, length(DataES.ChannelGood), 1) ~= DataES.FPeaksIndex));
        DataES.FPeaksOut = reshape(DataES.FPeaksOut, [], DataES.TimeN - 1);
        DataES.FPeaksOutMeans = mean(DataES.FPeaksOut, 2);
        %DataES.FStds = std(DataES.FPeaksOut, 0, 2);
        DataES.checkThresSDwt = logical(DataES.FPeaks > ThresSDwt * abs(DataES.FStds) + abs(DataES.FPeaksOutMeans));
        % peak amplitude greater than the upper limit (350uV by default)
        DataES.checkThresPeak = logical(DataES.FPeaks > ThresPeak);
        % variation threshold between consecutive points (25uV by default)
        DataES.FConsecutiveVariation = diff(DataES.FSelected, 1, 2);
        DataES.checkThresVars = any(DataES.FConsecutiveVariation > ThresVars, 2);
        % count number of channels with potential ES        
        DataES.ChannelBadSDwt = DataES.ChannelGood(DataES.checkThresSDwt);
        DataES.ChannelBadPeak = DataES.ChannelGood(DataES.checkThresPeak);
        DataES.ChannelBadVars = DataES.ChannelGood(DataES.checkThresVars);
        
        % store channel-wise data
        ChannelTrialsMeans(:, iFile) = DataES.FMeans;
        ChannelTrialsOutSDwt(:, iFile) = DataES.checkThresSDwt;
        ChannelTrialsOutPeak(:, iFile) = DataES.checkThresPeak;
        ChannelTrialsOutVars(:, iFile) = DataES.checkThresVars;
        
        % ===== SAVE RESULT =====
        % Progress bar
        bst_progress('text', 'Saving results...');
        bst_progress('set', progressPos + round(3 * iFile / length(sInputs) / 3 * 100));
        % save DataES to the trial structure
        save(file_fullpath(sInputs(iFile).FileName), 'DataES', '-append');        
    end
    
    % reject trials with peak amplitude x SDs away from the mean across the rest of the trials
    ChannelTrialsStds = abs(ThresSDat * std(ChannelTrialsMeans, [], 2));  % number of channels
    fprintf('------------------------------------------------------------------------------------------------------\n');
    fprintf('======== Trial-wise Report ========\n');
    for iTrial = 1:length(sInputs)
      % read ES data
      load(file_fullpath(sInputs(iTrial).FileName), 'DataES');
      % extract the mean values for both this trial and the rest of trials
      iTrialMean = ChannelTrialsMeans(:, iTrial);
      ChannelOtherTrialsMeans = ChannelTrialsMeans;
      ChannelOtherTrialsMeans(:, iTrial) = [];
      % check this trial
      ChannelTrialsOutSDat(:, iTrial) = iTrialMean > ChannelTrialsStds + abs(mean(ChannelOtherTrialsMeans, 2));  
      DataES.checkThresSDat = ChannelTrialsOutSDat(:, iTrial);
      DataES.ChannelBadSDat = ChannelGoodNames(logical(DataES.checkThresSDat));      
      
      % report the detected event and related channels
      if any(DataES.checkThresSDat) || any(DataES.checkThresSDwt) || any(DataES.checkThresPeak) || any(DataES.checkThresVars)
        % report detected trial and corresponding channels                       
        criterion = '';
        chan2disp = '';
        if any(DataES.BadOverlaps); criterion = 'overlaps with BAD segments + '; end
        if any(DataES.checkThresSDwt); criterion = [criterion, '[wtSD]']; chan2disp = sprintf('(%s), ', DataES.ChannelBadSDwt{:}); end
        if any(DataES.checkThresSDat); criterion = [criterion, '[atSD]']; chan2disp = [chan2disp, '|', sprintf('(%s), ', DataES.ChannelBadSDat{:})]; end
        if any(DataES.checkThresPeak); criterion = [criterion, '[Peak]']; chan2disp = [chan2disp, '|', sprintf('(%s), ', DataES.ChannelBadPeak{:})]; end                
        if any(DataES.checkThresVars); criterion = [criterion, '[Vars]']; chan2disp = [chan2disp, '|', sprintf('(%s), ', DataES.ChannelBadVars{:})]; end  
        DataES.ChannelBadN = sum((DataES.checkThresSDat + DataES.checkThresSDwt + DataES.checkThresPeak + DataES.checkThresVars) ~= 0);
        fprintf('Trial %s (#%d) %s: spikes are detected in %d channels %s.\n', DataES.TrialType, DataES.TrialIndex, criterion, DataES.ChannelBadN, chan2disp)    
      elseif any(DataES.BadOverlaps)
        fprintf('Trial %s (#%d) is in BAD segments but not detected by our method.\n', DataES.TrialType, DataES.TrialIndex)
      else
        iOk(iTrial) = true;
      end
      % save ES data
      save(file_fullpath(sInputs(iTrial).FileName), 'DataES', '-append');
    end
    
    % report the percent of rejected trials
    ChannelTrialsOutAll = ChannelTrialsOutSDat + ChannelTrialsOutSDwt + ChannelTrialsOutPeak + ChannelTrialsOutVars;
    TrialsOutN = nnz(sum(ChannelTrialsOutAll));
    fprintf('------------------------------------------------------------------------------------------------------\n');
    fprintf('%d out of %d trials are rejected in total.\n', TrialsOutN, length(sInputs));
    fprintf('%d trials are rejected according to the across-trials %d SDs criteria.\n', nnz(sum(ChannelTrialsOutSDat)), ThresSDat);
    fprintf('%d trials are rejected according to the within-trial %d SDs criteria.\n', nnz(sum(ChannelTrialsOutSDwt)), ThresSDwt);
    fprintf('%d trials are rejected according to the peak amplitude criteria.\n', nnz(sum(ChannelTrialsOutPeak)));
    fprintf('%d trials are rejected according to the consecutive variation criteria.\n', nnz(sum(ChannelTrialsOutVars)));
    
    % report bad channels
    fprintf('------------------------------------------------------------------------------------------------------\n');
    fprintf('======== Channel-wise Report ========\n');
    fprintf('Out of %d rejected trials : \n', TrialsOutN);
    chan2reject = [];
    for iChannel = 1:length(ChannelGoodNames)
      iChannelTrialsOutAll = ChannelTrialsOutAll(iChannel, :);
      iChannelTrialsOutN = nnz(iChannelTrialsOutAll);
      iChannelTrialsOutRatio = iChannelTrialsOutN / TrialsOutN;
      if iChannelTrialsOutRatio > 0
        fprintf('%.2f%% trials are rejected due to the channel %s.\n', iChannelTrialsOutRatio*100, ChannelGoodNames{iChannel});
        if iChannelTrialsOutRatio >= ThresChan/100
          chan2reject = [chan2reject, iChannel];
        end
      end
    end
    ChannelRejectNames = ChannelGoodNames(chan2reject);
    ChannelTrialsOutAllUpdate = ChannelTrialsOutAll;
    ChannelTrialsOutAllUpdate(chan2reject, :) = [];  % remove bad channels
    fprintf('In summary, each of the %d channels (%s) led to at least %d%% trials being rejected.\n', length(chan2reject), sprintf('[%s]', ChannelRejectNames{:}), ThresChan);
    fprintf('After removing the channels listed above, the number of rejected trials should be %d.\n', nnz(sum(ChannelTrialsOutAllUpdate)));
    
    % Return all the input files and bad files
    if ismarkbadtrials
      OutputFiles = {sInputs(iOk).FileName};
      BadFiles = {sInputs(~iOk).FileName};
      [BadPaths, BadFiles] = cellfun(@(x) fileparts(file_fullpath(x)), BadFiles, 'UniformOutput', 0);  % extract filenames that are rejected
      BadFiles = cellfun(@(x) sprintf('%s.mat', x), BadFiles, 'UniformOutput', 0);                     % add file extension - .mat
      OutputPaths = unique(BadPaths);
      for iOutputPath = 1:length(OutputPaths)
        iBadTrials = cell2mat(cellfun(@(x) strcmp(x, OutputPaths{iOutputPath}), BadPaths, 'UniformOutput', 0));
        BadTrials = BadFiles(logical(iBadTrials))';
        save(fullfile(OutputPaths{iOutputPath}, 'brainstormstudy.mat'), 'BadTrials', '-append');
      end
    end
end

%% sub-functions
function [TrialType, TrialIndex] = commentsplit(Comment)
  % this function is specifically used to extract trial type and index from BST trial Comment
  TrialComments = strsplit(Comment);
  TrialOrder = strsplit(TrialComments{2},{'(','#',')'});
  TrialType = TrialComments{1};         % trial type or condition  
  TrialIndex = str2double(TrialOrder{2});  % within-condition index of this trial
end

function [EventType, EventTime, EventBadTime] = eventscsv(EventFile)
  % this function is specifically used to read BST events and bad segments
  % read from CSV
  fid = fopen(EventFile);
    Events = textscan(fid,'%s%f%f','Delimiter',',');
  fclose(fid);
  EventType = Events{1};  % still cell array
  EventTime = [Events{2}, Events{3}];
  % identify BAD segments
  EventBADs = cell2mat(cellfun(@(x) strcmp(x, 'BAD_suspected'), EventType, 'UniformOutput', 0));
  EventBadTime = EventTime(EventBADs, :);
  EventBadTime(:, 2) = sum(EventBadTime, 2);  % start point and end point
end

function [ThisTrial, Overlaps]=badvalidate(TrialType, TrialIndex, TrialWindow, EventType, EventTime, EventBadTime)
  % to check the overlap between a certain trial and the bad segments
  Trials = cell2mat(cellfun(@(x) strcmp(x, TrialType), EventType, 'UniformOutput', 0));
  TrialsStartTime = EventTime(Trials, 1);
  TrialsTime = [TrialsStartTime + TrialWindow(1), TrialsStartTime + TrialWindow(end)];
  ThisTrial = TrialsTime(TrialIndex, :);
  Overlaps = ~logical((EventBadTime(:,1) > ThisTrial(2)) + (EventBadTime(:,2) < ThisTrial(1)));
end







