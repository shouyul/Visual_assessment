function new_events = read_trialtypes_VEP_SL(events, trialtypes_file)
% For each block, import the trial type (i.e., what stimulus was presented
% on the screen)
%
% Inputs:
%   - events:           - event structure (output of export_events_VEP)
%   - trialtypes_file   - [string] name of the trialtypes .csv file
%                           speciying trial types
%
% Output:
%   - new_events        - modified event structure with additional
%                           'trialtype'

TrialTypes = table2struct(readtable(trialtypes_file));

event_fields = fieldnames(events);
new_event_fields = [event_fields; {'trialtype'}];

events_count = numel(events);
for e = 1:events_count
    for f = 1:numel(new_event_fields)
        if f <= numel(event_fields)
            if any(strcmp({'block','trial'},new_event_fields{f}))
                new_events(e).(new_event_fields{f}) = events(e).(new_event_fields{f});
            else
                new_events(e).(new_event_fields{f}) = events(e).(new_event_fields{f});
            end
        elseif strcmp('trialtype',new_event_fields{f})
            if ~isempty(new_events(e).trial)
                new_events(e).(new_event_fields{f}) = TrialTypes(new_events(e).trial).TrialType;
            end
        end
    end
end