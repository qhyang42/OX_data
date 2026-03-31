%% manual marking of inhales 

%% 
wkdir =     '/Users/qhyang/Desktop/OX_DATA/labchart'; 
subjname = 'subj3'; 
load(fullfile(wkdir, [subjname, '_events_raw.mat'])); 

evs_manual = []; 
inhale_events = []; 
%% detection of crude events 
% nruns = length(eventdata); 
% mripulse = cell(nruns,1);
% events = cell(nruns,1);
% resp = cell(nruns,1);
% 
% for i = 1: nruns
% 
%     mripulse{i} = eventdata{i}(3, :);
%     events{i} = eventdata{i}(2, :);
%     resp{i} = eventdata{i}(1, :);
% end 
% 
% event_onsets = zeros(10,nruns); 
% cue_onsets = zeros(10, nruns); 
% 
% for evidx = 1:nruns
%     currentevents = events{evidx};
%     eventbi = currentevents>0.1;
%     event_real = eventbi;
%     for i = 1:length(eventbi)-1
%         if eventbi(i+1) == eventbi(i)
%             event_real(i+1) = 0;
%         end
%     end
%     evonsets = find(event_real>0);
%     event_o = evonsets(5:2:end);
%     event_q = evonsets(4:2:end-1);
% 
%     mrionset = find(mripulse{evidx} > 1, 1); 
%     event_o = event_o - mrionset; % set MRI onset to time 0 
%     event_o = event_o/1000;
%     event_onsets(:, evidx) = event_o; 
%     
%     event_q = event_q - mrionset; % set MRI onset to time 0 
%     event_q = event_q/1000;
%     cue_onsets(:, evidx) = event_q; 
% 
% end

%%
runidx = 7; 
resp = eventdata{runidx}(1, :); 
event = eventdata{runidx}(2, :);

%% filter resp 
resp_filt = ft_preproc_lowpassfilter(resp, 1000, 2, [], 'fir'); 
resp_filt = resp_filt-mean(resp_filt); 
%% run event picker 
inhale_events{runidx} = EventPicker([resp_filt; event], []); 

%% 
evs_manual{runidx} = inhale_events{runidx}.evs(inhale_events{runidx}.selected_evs); 


