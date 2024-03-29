function [data, trl] = eeg_getCumulusData(cfg)
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
% Added support for continuous data without trialdef-field (03.05.2019, jk)
% Fixed missing info for german data (17.05.2021, jk)

%% 1. Set Basics of data format
tic

event = [];
headerLines = 0;
footerLines = 2;
breakLines = 8; 

filename = cfg.dataset;
%% 2. Parsing file
disp('Parsing File...')
iLine = 0;
fid = fopen(filename);
while  1
    iLine = iLine + 1;
    tline = fgetl(fid);
        
    if ~ischar(tline)
        break
    end
    
    if iLine == 4
        tmp = strsplit(tline);
        freq = str2num(tmp{4});
    end
    
    if iLine == 8
        tmp = strsplit(tline);
        for l = 3:length(tmp)
            if strcmpi(tmp{l}(end),',')
                labels{l-2} = tmp{l}(1:end-1);
            else
                labels{l-2} = tmp{l}(1:end);
            end
        end
    end
    
    if strcmpi(tline(1:12),'% CLK_ADJUST')
        breakLines(end+1) = iLine;
        
    end
end

fclose(fid);
disp('done!');

nEEG = length(labels);
nIndex = 1;
nAccel = 3;
nTime = 1;
% %% 3. Get the data labels
% disp('Getting Labels...')
%     if size(labels,2)==1
%         tmplabels = regexprep(labels{1}, '\s+', '_');
%         tmplabels = strsplit(tmplabels,'_');
%         tmplabels{end+1} = 'placeholder'; % Remove last Column Name
%     else
%         tmplabels = labels;
%     end
%     nSignals = length(tmplabels);
%     if strcmpi(lang_format,'UK')
%         nEEG = find(strcmp(tmplabels,'Events'))-1;
%     elseif strcmpi(lang_format,'GE')
%         nEEG = find(strcmp(tmplabels,'Ereignisse'))-1;
%     end
% disp('done!');

%% 4. Get the number of data lines and read in data
disp('Getting Data...')
for d = 1:length(breakLines)-1
    M{d} = csvread(filename,breakLines(d)+2,0,[breakLines(d)+2,0,breakLines(d+1)-2,(nIndex + nEEG + nAccel + nTime)-1]);

    %% 5. Get Time Index
    disp('Getting the Time...')
    samples = M{d}(:,1);
    time{d} = samples / freq;

% %% 6. Get Events
% disp('Getting the Events...')
% clear event
% j = 0;
% eventCol = nSignals-nEEG;
% for iEvent = 1:size(M{2},1)
%     tmp = cell2mat(M{2}(iEvent,eventCol));
%     if ~isempty(tmp) || ~strcmp(tmp,'')
%         j = j+1;
%         event(j).samples = samples(iEvent);
%         event(j).time = time(iEvent);
%         tmp2 = strsplit(tmp,'(');
%         event(j).type = str2num(tmp2{1});
%         event(j).name = tmp;
%     end
% end
% 
% if exist('event','var') == 0
%     event = [];
% end
% 
% %if isempty(event)
% %    return
% %end
% disp('done!');

    %% 7. Force data into FT-Format
    disp('Building FT Data...')
    data.label = labels;         % cell-array containing strings, Nchan*1
    data.fsample = freq;            % sampling frequency in Hz, single number
    data.trial{d} = M{d}(:,2:nEEG+1)';             % cell-array containing a data matrix for each 
                                    % trial (1 X Ntrial), each data matrix is a Nchan*Nsamples matrix 
    data.time{d} = time{d}';              % cell-array containing a time axis for each 
                                    % trial (1 X Ntrial), each time axis is a 1*Nsamples vector 
end
disp('done!');
                
%% Build Trial definition as trl matrix to be used with fieldtrip 
if isfield(cfg,'trialdef')
    disp('Building trl-Matrix...')
    
    if exist(event,'var')
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
    elseif isfield(cfg.trialdef,'triallength')
        trl = cell(1,length(data.trial));
        clear borders begsample endsample
  
        for t = 1:length(data.trial)
            
            tmp = data.time{t}(1):cfg.trialdef.triallength:data.time{t}(end);
            for b = 1:length(tmp)
                borders(b) = nearest(data.time{t},tmp(b));
            end
                
            begsample = borders(1:end-1)';
            endsample = borders(2:end)';
            offset = zeros(length(begsample),1);
            
            trl{t} = round([begsample endsample offset]);
        end
    end
else
    trl = [];
end
disp('done!');
%%
T = toc;
fprintf('This took %i seconds\n',T)