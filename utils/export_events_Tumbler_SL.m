function [events] = export_events_Tumbler_SL(Event_stream,times)
% Interpret event XDF stream for PIONEER VEP experiemtn to put them in the
% EEGLAB set
%
% Inputs:
%   - Event_stream:         Event stream loaded from the xdf
%   - times:                Complete time points vector of the recording
%                               (to compute latencies)
%
% Outputs:
%   - events:               Structure containing all events loaded

required_fields = {'type', 'latency', 'block', 'trial'};%, 'answer', 'delay'};
required_types = {'str', 'num', 'str', 'num'};%, 'str', 'num'};

events_count = length(Event_stream.time_stamps);
command = '';

for f = 1:numel(required_fields)
    if strcmp(required_fields{f}, 'duration')
        command = [command,'''',required_fields{f},''',num2cell(ones(1, events_count)),'];
    elseif strcmp(required_types{f}, 'str')
        command = [command,'''',required_fields{f},''','''','];
    elseif strcmp(required_types{f}, 'num')
        command = [command,'''',required_fields{f},''',[],'];
    end
end

% Remove the last ','
command = command(1:end-1);
eval(['events = struct(',command,');']);

events = struct('type',[],'latency',[],'block',[],'trial',[]);

t = 0;
for e = 1:events_count
    current_event = split(Event_stream.time_series(e),';');
    if ~ismember({'event:Flip'},current_event)
        t = t+1;
        for m = 1:length(current_event)
            event_marker = split(current_event{m},':');
            switch event_marker{1}
%                 case 'event'
%                     events(t).type = event_marker{2};
                case 'block'
                    events(t).block = str2num(event_marker{2})+1;
                case 'trial'
                    events(t).trial = str2num(event_marker{2})+1;
                otherwise
                    events(t).type = event_marker{1};
            end
        end
        [~,events(t).latency] = min(abs(times - Event_stream.time_stamps(e)));
    end
end

events = orderfields(events,required_fields);






end