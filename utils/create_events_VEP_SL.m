function [events] = create_events_VEP_SL(Annotations,offset,trialtypes)
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

required_fields = {'type', 'latency', 'block', 'trial', 'trialtype','duration'};%, 'answer', 'delay'};
required_types = {'str', 'num', 'str', 'num','str','num'};%, 'str', 'num'};

events_count = min(size(Annotations,1),length(offset));

events = struct('type',[],'latency',[],'block',[],'trial',[],'trialtype','','duration',[]);

t = 0;
for e = 1:events_count
    current_event = split(Annotations(e,:).Annotations,';');
    if ~ismember({'event:Flip'},current_event)
        t = t+1;
        for m = 1:length(current_event)
            event_marker = split(current_event{m},':');
            switch event_marker{1}
                case 'event'
                    events(t).type = event_marker{2};
                case 'block'
                    events(t).block = str2num(event_marker{2});
                case 'trial'
                    events(t).trial = str2num(event_marker{2});
            end
        end
        events(t).latency = offset(e);
        try
            events(t).trialtype = trialtypes([trialtypes.block]==events(t).block & [trialtypes.trial]==events(t).trial).TrialType;
        end

    end
end

events = orderfields(events,required_fields);






end