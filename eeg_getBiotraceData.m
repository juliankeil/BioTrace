function [data, trl] = eeg_getBiotraceData(cfg,filename)
% Import Biotrace data and switch to FieldTrip Format
%
% Requires:
% cfg.trialdef.eventvalue = Vector containing the trigger labels
% cfg.trialdef.prestim = How much time before each trigger
% cfg.trialdef.poststim = Hoch much time after each trigger
%
% Example:
% cfg=[];
% cfg.trialdef.eventvalue = [1:255];
% cfg.trialdef.prestim = .2;
% cfg.trialdef.poststim = .2;
% [data,trl] = eeg_getBiotraceData(cfg,'triggertest_altesoftware.txt');
% Hacked by Julian Keil 2018

% Fixed error with german data (07.11.2018, jk)
% Avoid fucking triggers 10 and 13 in programming! they fuck up the ASCII
% files!
% Fixed more errors related to different matlab versions and text reading
% (15.04.2019, jk)
%% 1. Set Basics of data format
tic

event = [];
headerLines = 14;
footerLines = 2;

%% 2. Parsing file
disp('Parsing File...')
iLine = 0;
feof = 0;
fid = fopen(filename);
while  feof == 0
    iLine = iLine + 1;
    tline = fgetl(fid);
    if iLine == 1   %based on content in first line of file determine Biotrace language format %%MB19.3.2014
        if strcmp(tline,'Unbearbeiteter Daten-Export (tab separat)')
            lang_format = 'GE';
        elseif strcmp(tline,'RAW Data export file (tab separated)')
            lang_format = 'UK';
        else
            error('The first line doesn''t seem to be from a Biotrace text export.');
        end
    end
    
    if iLine == 9
        if strcmp(lang_format, 'GE') %%MB19.3.2014
            freq = sscanf(tline,'Ausgabegeschwindigkeit:\t%d\tSamples/sek.');
        elseif  strcmp(lang_format, 'UK')
            freq = sscanf(tline,'Output rate:\t%d\tSamples/sec.');
        end
    end
    
    if iLine == 8
        if strcmp(lang_format, 'GE') %%MB19.3.2014
            duration = sscanf(tline,'Dauer:\t%d\tSeconds.');
        elseif  strcmp(lang_format, 'UK')
            duration = sscanf(tline,'Duration:\t%d\tSeconds.');
        end
    end
    
    if iLine == 12
        labels = strsplit(tline,'\t');
        feof = 1; % Set Exit
    end
end
fclose(fid);
disp('done!');

%% 3. Get the data labels
disp('Getting Labels...')
    if size(labels,2)==1
        tmplabels = regexprep(labels{1}, '\s+', '_');
        tmplabels = strsplit(tmplabels,'_');
        tmplabels{end+1} = 'placeholder'; % Remove last Column Name
    else
        tmplabels = labels;
    end
    nSignals = length(tmplabels);
    nEEG = find(strcmp(tmplabels,'Events'))-1;
disp('done!');

%% 4. Get the number of data lines and read in data
disp('Getting Data...')
skip_samples = 1;  %avoid data errors at beginning

tmp = repmat('%f',1,nEEG);
tmp2 = repmat('%s',1,(nSignals-nEEG)+1);
tmp3 = [tmp,tmp2];
fid = fopen(filename);
M = textscan(fid, tmp3,...
    'HeaderLines',headerLines+skip_samples,...
    'CollectOutput', 1,'ReturnOnError',1);
fclose(fid);
disp('done!');

%% 5. Get Time Index
disp('Getting the Time...')
samples = (M{1}(:,1) - M{1}(1,1));
time = samples / freq;
disp('done!');

%% 6. Get Events
disp('Getting the Events...')
clear event
j = 0;
eventCol = nSignals-nEEG;
for iEvent = 1:size(M{2},1)
    tmp = cell2mat(M{2}(iEvent,eventCol));
    if ~isempty(tmp) || ~strcmp(tmp,'')
        j = j+1;
        event(j).samples = samples(iEvent);
        event(j).time = time(iEvent);
        tmp2 = strsplit(tmp,'(');
        event(j).type = str2num(tmp2{1});
        event(j).name = tmp;
    end
end

if isempty(event)
    return
end
disp('done!');

%% 7. Force data into FT-Format
disp('Building FT Data...')
data.label = tmplabels(2:nEEG);         % cell-array containing strings, Nchan*1
data.fsample = freq;            % sampling frequency in Hz, single number
data.trial{1} = M{1}(:,2:nEEG)';             % cell-array containing a data matrix for each 
                                % trial (1 X Ntrial), each data matrix is a Nchan*Nsamples matrix 
data.time{1} = time';              % cell-array containing a time axis for each 
                                % trial (1 X Ntrial), each time axis is a 1*Nsamples vector 
disp('done!');
                
%% Build Trial definition as trl matrix to be used with fieldtrip                
disp('Building trl-Matrix...')
trl=[];

for i=1:length(event);
    %%% Find the Correct Trials
    if ismember(event(i).type,cfg.trialdef.eventvalue);
      % add this to the trl definition
      begsample = event(i).samples - cfg.trialdef.prestim*freq;
      endsample = event(i).samples + cfg.trialdef.poststim*freq - 1;
      offset = -cfg.trialdef.prestim*freq;  % Pre Trigger
      
      trl(end+1, :) = round([begsample endsample offset event(i).type]); 
      
    end % if
end % event 
disp('done!');
T = toc;
fprintf('This took %i seconds\n',T)