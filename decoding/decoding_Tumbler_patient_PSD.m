function [result]  = decoding_Tumbler_TL_PSD(data_event,elecs,perm_test)

%% set parameters
% srate = 250;
% subjid = data.subject;
% 
% %% select electrodes
% elecs_label = elecgroup.Elecs_label;
% elecs = [];
% 
% for i = 1:length(data.chanlocs)
%     if ismember(data.chanlocs(i).labels,elecs_label)
%         elecs = [elecs,i];
%     end
% end

% data_event = data_event([data_event.trialtype_bin]==pair(1)|[data_event.trialtype_bin]==pair(2));

%% construct data mat
% elecs = 1:127;

gogglecon = 'GogglesON';
data_ON = data_event(strcmp({data_event.condition},gogglecon)); 
gogglecon = 'GogglesOFF';
data_OFF = data_event(strcmp({data_event.condition},gogglecon));

data_mat = [];
data_label = [];
for i = 1:length(data_event)
    data_mat_EyesOpen(i,:,:) = data_event(i).PSD_EyesOpen(elecs,:);
    data_mat_EyesClosed(i,:,:) = data_event(i).PSD_EyesClosed(elecs,:);
    data_label(i,1) = data_event(i).trialtype_bin;
    data_label(i,2) = data_event(i).block;
end


data_mat_ON = [];
data_label_ON = [];
for i = 1:length(data_ON)
    data_mat_ON_EyesOpen(i,:,:) = data_ON(i).PSD_EyesOpen(elecs,:);
    data_mat_ON_EyesClosed(i,:,:) = data_ON(i).PSD_EyesClosed(elecs,:);
    data_label_ON(i,2) = data_ON(i).block;
    data_label_ON(i,1) = strcmp(data_ON(i).trialtype,'WithObject');
end

data_label_ON(:,2) = 1:length(data_label_ON);

% goggles OFF
data_mat_OFF = [];
data_label_OFF = [];
for i = 1:length(data_OFF)
    data_mat_OFF_EyesOpen(i,:,:) = data_OFF(i).PSD_EyesOpen(elecs,:);
    data_mat_OFF_EyesClosed(i,:,:) = data_OFF(i).PSD_EyesClosed(elecs,:);
    data_label_OFF(i,2) = data_OFF(i).block;
    data_label_OFF(i,1) = strcmp(data_OFF(i).trialtype,'WithObject');
end

data_label_OFF(:,2) = 1:length(data_label_OFF);

if perm_test
    data_label_ON(:,1) = data_label_ON(randperm(length(data_label)),1);
    data_label_OFF(:,1) = data_label_OFF(randperm(length(data_label)),1);
end

%% average trials within a block
% blocks = unique(data_label(:,2));
% cats = unique(data_label(:,1));
% 
% data_mat_new = [];
% data_label_new = [];
% 
% for i = 1:length(cats)
%     for j = 1:length(blocks)
%         temp = data_mat(data_label(:,1)==cats(i)&data_label(:,2)==blocks(j),:,:);
%         temp = mean(temp,1);
% 
%         data_mat_new = cat(1,data_mat_new,temp);
%         data_label_new = cat(1,data_label_new,[cats(i),blocks(j)]);
% 
%     end
% end
% 
% 
% data_mat = data_mat_new;
% data_label = data_label_new;

%% Cumulative decoding

% entire trial

data_mat_OFF_EyesClosed_all = reshape(data_mat_OFF_EyesClosed,size(data_mat_OFF_EyesClosed,1),[]);
data_mat_OFF_EyesOpen_all = reshape(data_mat_OFF_EyesOpen,size(data_mat_OFF_EyesOpen,1),[]);

data_mat_ON_EyesClosed_all = reshape(data_mat_OFF_EyesClosed,size(data_mat_OFF_EyesClosed,1),[]);
data_mat_ON_EyesOpen_all = reshape(data_mat_OFF_EyesOpen,size(data_mat_OFF_EyesOpen,1),[]);

train_acc_ON_close = [];
train_acc_ON_open = [];
train_acc_OFF_close = [];
train_acc_OFF_open = [];

test_acc_ON_close = [];
test_acc_ON_open = [];
test_acc_OFF_close = [];
test_acc_OFF_open = [];





% cvp = cvpartition(data_label,'LeaveOut');
[train_acc_OFF_close,test_acc_OFF_close] = compute_train_test_acc(data_mat_OFF_EyesClosed_all,data_label_OFF);
[train_acc_OFF_open,test_acc_OFF_open] = compute_train_test_acc(data_mat_OFF_EyesOpen_all,data_label_OFF);

[train_acc_ON_close,test_acc_ON_close] = compute_train_test_acc(data_mat_ON_EyesClosed_all,data_label_ON);
[train_acc_ON_open,test_acc_ON_open] = compute_train_test_acc(data_mat_ON_EyesOpen_all,data_label_ON);

% 5s chunks
% data_mat_5s = [];
% data_label_5s = [];
% chunk_len = 5;
% num_chunk = trial_length/chunk_len;
% for tr = 1:size(data_mat,1)
%     temp_trial_el = [];
%     for el = 1:size(data_mat,2)
%         temp_trial = squeeze(data_mat(tr,el,:));
%         new_trial = transpose(reshape(temp_trial,srate*chunk_len,[]));
%         temp_trial_el(el,:,:) = new_trial;
%     end
%     data_label_5s = [data_label_5s;repmat(data_label(tr,:),num_chunk,1)];
%     data_mat_5s = [data_mat_5s;permute(temp_trial_el,[2,1,3])];
% end
% 
% data_mat_5s = reshape(data_mat_5s,size(data_mat_5s,1),[]);
% 
% % cvp = cvpartition(data_label_5s,'LeaveOut');
% [train_acc_5s,test_acc_5s] = compute_train_test_acc(data_mat_5s,data_label_5s);
% 
% 
% % 2.5s chunks
% data_mat_25s = [];
% data_label_25s = [];
% chunk_len = 2.5;
% num_chunk = trial_length/chunk_len;
% for tr = 1:size(data_mat,1)
%     temp_trial_el = [];
%     temp_el_label = [];
%     for el = 1:size(data_mat,2)
%         temp_trial = squeeze(data_mat(tr,el,:));
%         new_trial = transpose(reshape(temp_trial,srate*chunk_len,[]));
%         temp_trial_el(el,:,:) = new_trial;
%     end
%     data_label_25s = [data_label_25s;repmat(data_label(tr,:),num_chunk,1)];
%     data_mat_25s = [data_mat_25s;permute(temp_trial_el,[2,1,3])];
% end
% 
% data_mat_25s = reshape(data_mat_25s,size(data_mat_25s,1),[]);
% 
% % cvp = cvpartition(data_label_25s,'LeaveOut');
% [train_acc_25s,test_acc_25s] = compute_train_test_acc(data_mat_25s,data_label_25s);
% 
% 
% % 1.25s chunks
% data_mat_125s = [];
% data_label_125s = [];
% chunk_len = 1.25;
% num_chunk = trial_length/chunk_len;
% for tr = 1:size(data_mat,1)
%     temp_trial_el = [];
%     temp_el_label = [];
%     for el = 1:size(data_mat,2)
%         temp_trial = squeeze(data_mat(tr,el,:));
%         new_trial = [];
%         for c = 1:num_chunk
%             start_ind = floor(srate*chunk_len*(c-1))+1;
%             new_trial(c,:) = temp_trial(start_ind:start_ind+311);
%         end
%         temp_trial_el(el,:,:) = new_trial;
%     end
%     data_label_125s = [data_label_125s;repmat(data_label(tr,:),num_chunk,1)];
%     data_mat_125s = [data_mat_125s;permute(temp_trial_el,[2,1,3])];
% end
% 
% data_mat_125s = reshape(data_mat_125s,size(data_mat_125s,1),[]);
% 
% % cvp = cvpartition(data_label_125s,'LeaveOut');
% [train_acc_125s,test_acc_125s] = compute_train_test_acc(data_mat_125s,data_label_125s);


%% compile results

% conditions = {'Control','Disc','VerticalBar','HorizontalBar'};
% comb = sprintf('%s VS %s',conditions{pair(1)+1},conditions{pair(2)+1});
% 
% result = struct( 'Comparison',{comb;comb;comb;comb},...
%                 'ChunkType',{'Whole trial';'5s chunks';'2.5s chunks';'1.25 chunks'},...
%                 'TrainingAccuracy',{train_acc_whole;train_acc_5s;train_acc_25s;train_acc_125s},...
%                 'TestingAccuracy',{test_acc_whole;test_acc_5s;test_acc_25s;test_acc_125s});

result = struct('GogglesCondition',{'GogglesOFF','GogglesOFF','GogglesON','GogglesON'},...
    'Phase',{'Close','Open','Close','Open'},...
    'TrainingAccuracy',{train_acc_OFF_close,train_acc_OFF_open,train_acc_ON_close,train_acc_ON_open},...
    'TestingAccuracy',{test_acc_OFF_close,test_acc_OFF_open,test_acc_ON_close,test_acc_ON_open});
end
