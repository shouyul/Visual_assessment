%% temporal decoding of object present/absent
function [result] = decoding_tumbler_patient(data,elecgroup,num_perm,run_timecourse)
%% set parameters
srate = data.srate;
subjid = data.subject;

% set time bound for the three phases: eyes closed, eyes open, respond
duration1 = 5;
duration2 = 15;

time_bd = [srate*duration1,srate*(duration1+duration2)];


%% setment data into events

% create event data segments
data_event = struct([]);
trial_length = inf;
for e = 1:length(data.event)
    if strcmp(data.event(e).type,'TrialStart')
        trial_start_idx = floor(data.event(e).latency);
        for e2 = e:length(data.event)
            if strcmp(data.event(e2).type,'TrialEnd')
                trial_end_ind = floor(data.event(e2).latency);
                break
            elseif e2 == length(data.event)
                trial_end_ind = size(data.data,2);
                break
            end
        end

        temp_event = data.event(e);
%         temp_event.answer = data.event(e2-3).answer;
%         temp_event.delay = data.event(e2-3).delay;
        temp_event.data = data.data(:,trial_start_idx:trial_end_ind);
        
        temp_event = rmfield(temp_event,{'type','duration','latency','urevent'});

        trial_length = min(trial_length,trial_end_ind-trial_start_idx+1);

%         if strcmp(data.event(e-1).type,'WithObject')
%             temp_event.trialtype = 'WithObject';
%             temp_event.trialtype_bin = 1;
%         elseif strcmp(data.event(e-1).type, 'WithoutObject')
%             temp_event.trialtype = 'WithoutObject';
%             temp_event.trialtype_bin = 0;
%         else
%             temp_event.trialtype = 'Unknown';
%             temp_event.trialtype_bin = [];
%         end
        temp_event.trialtype_bin = strcmp(temp_event.trialtype,'WithObject');

        if isempty(data_event)
            data_event = temp_event;
        else
            data_event = [data_event,temp_event];
        end

    end
end

% equate trial length
% if strcmp(subjid,'TEST')
%     trial_length = 6250;
% else
%     trial_length = 7500;
% end



for e = 1:length(data_event)
    data_event(e).data = data_event(e).data(:,1:trial_length);
end

%% zscore data
for e = 1:length(data_event)
    data_event(e).data = zscore(data_event(e).data,[],2);
end

%% select electrodes
elecs_label = elecgroup.Elecs_label;
elecs = [];

for i = 1:length(data.chanlocs)
    if ismember(data.chanlocs(i).labels,elecs_label)
        elecs = [elecs,i];
    end
end

%data_forClass = nan(size(data_event(1).data,1),2,10,size(data_event(1).data,2));

%% goggles condition
gogglecon = 'GogglesON';
data_ON = data_event(strcmp({data_event.condition},gogglecon)); 
gogglecon = 'GogglesOFF';
data_OFF = data_event(strcmp({data_event.condition},gogglecon));

% data_OFF = data_event;

%% construct data mat
% goggles ON
data_mat_ON = [];
data_label_ON = [];
for i = 1:length(data_ON)
    data_mat_ON(i,:,:) = data_ON(i).data(elecs,:);
    data_label_ON(i,2) = data_ON(i).block;
    data_label_ON(i,1) = strcmp(data_ON(i).trialtype,'WithObject');
end

data_label_ON(:,2) = 1:length(data_label_ON);

% goggles OFF
data_mat_OFF = [];
data_label_OFF = [];
for i = 1:length(data_OFF)
    data_mat_OFF(i,:,:) = data_OFF(i).data(elecs,:);
    data_label_OFF(i,2) = data_OFF(i).block;
    data_label_OFF(i,1) = strcmp(data_OFF(i).trialtype,'WithObject');
end

data_label_OFF(:,2) = 1:length(data_label_OFF);

%% run timecourse
% class_error_time_ON = [];
% class_error_time_OFF = [];
% 
% if run_timecourse
% 
%     if ~isempty(data_mat_ON)
%         % goggles ON
%         disp('running goggles on, timecourse')
%         parfor t = 1:size(data_mat_ON,3)
%             SVMModel = fitcsvm(data_mat_ON(:,:,t),data_label_ON);
%             CVSVMModel = crossval(SVMModel);
%             class_error_time_ON(t) = kfoldLoss(CVSVMModel);
%         end
%     end
% 
%     if ~isempty(data_mat_OFF)
%         % goggles OFF
%         disp('running goggles off, timecourse')
%         parfor t = 1:size(data_mat_OFF,3)
%             SVMModel = fitcsvm(data_mat_OFF(:,:,t),data_label_OFF);
%             CVSVMModel = crossval(SVMModel);
%             class_error_time_OFF(t) = kfoldLoss(CVSVMModel);
%         end
%     end
% 
% end

%% cumulative decoding

data_mat_all_OFF_close = reshape(data_mat_OFF(:,:,1:time_bd(1)),size(data_mat_OFF,1),[]);
data_mat_all_OFF_open = reshape(data_mat_OFF(:,:,time_bd(1)+1:time_bd(2)),size(data_mat_OFF,1),[]);

[train_acc_OFF_close,test_acc_OFF_close] = compute_train_test_acc(data_mat_all_OFF_close,data_label_OFF);
[train_acc_OFF_open,test_acc_OFF_open] = compute_train_test_acc(data_mat_all_OFF_open,data_label_OFF);


data_mat_all_ON_close = reshape(data_mat_ON(:,:,1:time_bd(1)),size(data_mat_ON,1),[]);
data_mat_all_ON_open = reshape(data_mat_ON(:,:,time_bd(1)+1:time_bd(2)),size(data_mat_ON,1),[]);

[train_acc_ON_close,test_acc_ON_close] = compute_train_test_acc(data_mat_all_ON_close,data_label_ON);
[train_acc_ON_open,test_acc_ON_open] = compute_train_test_acc(data_mat_all_ON_open,data_label_ON);


% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% % class_error_ON_close = nan(1,num_perm);
% % class_error_ON_open = nan(1,num_perm);
% % class_error_ON_resp = nan(1,num_perm);
% class_error_OFF_close = nan(1,num_perm);
% class_error_OFF_open = nan(1,num_perm);
% class_error_OFF_resp = nan(1,num_perm);
% 
% % predict_ON_close = nan(length(data_label_ON),num_perm);
% % predict_ON_open = nan(length(data_label_ON),num_perm);
% % predict_ON_resp = nan(length(data_label_ON),num_perm);
% predict_OFF_close = nan(length(data_label_OFF),num_perm);
% predict_OFF_open = nan(length(data_label_OFF),num_perm);
% predict_OFF_resp = nan(length(data_label_OFF),num_perm);
% 
% % goggles ON
% 
% if ~isempty(data_mat_ON)
%     % phase 1 - eyes closed
%     disp('Running goggles on, eyes closed');
%     parfor p = 1:num_perm
%         data_mat_all_ON_close = reshape(data_mat_ON(:,:,1:time_bd(1)),size(data_mat_ON,1),[]);
%         SVMModel = fitcecoc(data_mat_all_ON_close,data_label_ON);
%         CVSVMModel = crossval(SVMModel);
%         class_error_ON_close(p) = kfoldLoss(CVSVMModel);
%         predict_ON_close(:,p) = kfoldPredict(CVSVMModel);
%     end
%     disp(sprintf('Mean class_error: %2.4f',mean(class_error_ON_close)));
% 
%     % phase 2 - eyes open
%     disp('Running goggles on, eyes open')
%     parfor p = 1:num_perm
%         data_mat_all_ON_open = reshape(data_mat_ON(:,:,time_bd(1)+1:time_bd(2)),size(data_mat_ON,1),[]);
%         SVMModel = fitcecoc(data_mat_all_ON_open,data_label_ON);
%         CVSVMModel = crossval(SVMModel);
%         class_error_ON_open(p) = kfoldLoss(CVSVMModel);
%         predict_ON_open(:,p) = kfoldPredict(CVSVMModel);
%     end
%     disp(sprintf('Mean class_error: %2.4f',mean(class_error_ON_open)));
% 
%     % phase 3 - response
%     disp('Running goggles on, response')
%     parfor p = 1:num_perm
%         data_mat_all_ON_resp = reshape(data_mat_ON(:,:,time_bd(2)+1:end),size(data_mat_ON,1),[]);
%         SVMModel = fitcecoc(data_mat_all_ON_resp,data_label_ON);
%         CVSVMModel = crossval(SVMModel);
%         class_error_ON_resp(p) = kfoldLoss(CVSVMModel);
%         predict_ON_resp(:,p) = kfoldPredict(CVSVMModel);
%     end
%     disp(sprintf('Mean class_error: %2.4f',mean(class_error_ON_resp)));
% end
% 
% % goggles OFF
% if ~isempty(data_mat_OFF)
%     % phase 1 - eyes closed
%     disp('Running goggles off, eyes closed')
%     parfor p = 1:num_perm
%         data_mat_all_OFF_close = reshape(data_mat_OFF(:,:,1:time_bd(1)),size(data_mat_OFF,1),[]);
%         SVMModel = fitcecoc(data_mat_all_OFF_close,data_label_OFF);
%         CVSVMModel = crossval(SVMModel);
%         class_error_OFF_close(p) = kfoldLoss(CVSVMModel);
%         predict_OFF_close(:,p) = kfoldPredict(CVSVMModel);
%     end
%     disp(sprintf('Mean class_error: %2.4f',mean(class_error_OFF_close)));
% 
%     % phase 2 - eyes open
%     disp('Running goggles off, eyes open')
%     parfor p = 1:num_perm
%         data_mat_all_OFF_open = reshape(data_mat_OFF(:,:,time_bd(1)+1:time_bd(2)),size(data_mat_OFF,1),[]);
%         SVMModel = fitcecoc(data_mat_all_OFF_open,data_label_OFF);
%         CVSVMModel = crossval(SVMModel);
%         class_error_OFF_open(p) = kfoldLoss(CVSVMModel);
%         predict_OFF_open(:,p) = kfoldPredict(CVSVMModel);
%     end
%     disp(sprintf('Mean class_error: %2.4f',mean(class_error_OFF_open)));
% 
%     % phase 3 - response
%     disp('Running goggles off, resp')
%     parfor p = 1:num_perm
%         data_mat_all_OFF_resp = reshape(data_mat_OFF(:,:,time_bd(2)+1:end),size(data_mat_OFF,1),[]);
%         SVMModel = fitcecoc(data_mat_all_OFF_resp,data_label_OFF);
%         CVSVMModel = crossval(SVMModel);
%         class_error_OFF_resp(p) = kfoldLoss(CVSVMModel);
%         predict_OFF_resp(:,p) = kfoldPredict(CVSVMModel);
%     end
%     disp(sprintf('Mean class_error: %2.4f',mean(class_error_OFF_resp)));
% end
% 
% disp('Done.')

%% collate results

result = struct('GogglesCondition',{'GogglesON','GogglesON','GogglesOFF','GogglesOFF'},...
    'Phase',{'Close','Open','Close','Open'},...
    'TrainingAccuracy',{train_acc_ON_close,train_acc_ON_open,train_acc_OFF_close,train_acc_OFF_open},...
    'TestingAccuracy',{test_acc_ON_close,test_acc_ON_open,test_acc_OFF_close,test_acc_OFF_open});



%% save decoding results
% file_path = '/Users/shouyuling/Desktop/Postdoc/Projects/Gensight/Data/TumblerTask/Results/';
% file_name = sprintf('%s_class_error_%s',subjid,elecgroup.Location);
% 
% save([file_path,file_name],'class_error');


end






