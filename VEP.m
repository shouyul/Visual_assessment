%% set parameter

clear; clc

main_fold = '/Users/shouyuling/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Postdoc/Projects/Gensight/Data/VEPTask';
data_fold = fullfile(main_fold,'Analysis');
fig_fold = fullfile(main_fold,'Figures');
result_fold = fullfile(main_fold,'Results');


subs = {'004-4006-s1','004-4006-s2','004-4010-s3','004-4010-s4','004-4012-s1','004-4012-s2','004-4010-s5'};

trial_len = 10;
crossval_method = 'leaveout';
kfold = 10;
num_perm = 1000;

elec_labels = {'L16Z','Z16Z','R16Z',...
    'L16L','Z16L','Z16R','R16R',...
    'L17Z','Z17Z','R17Z',...
    'L17L','Z17L','Z17R','R17R',...
    'L18Z','Z18Z','R18Z',...
    'L18L','Z18L','Z18R','R18R',...
    'L19Z','Z19Z','R19Z'...
    'L19L','Z19L','Z19R','R19R',...
    'L20Z','Z20Z','R20Z'};

%%
for sub = 1:length(subs)

    subjid = subs{sub};
    disp(subjid)

    %%
    filename = fullfile(result_fold,strcat(subjid,'_results.mat'));
    if exist(filename,'file')
        load(filename);
    else
        results = struct();
        results.Subjid = subjid;
    end

    %% load data
    if ~isfield(results,'Info')
        disp('Segmenting data')
        filename = fullfile(data_fold,'2_preprocessing','bemobil','APP','Brain', ...
            strcat(subjid,'_cleaned_with_ICA.set'));
        %     filename = fullfile(data_fold,'2_preprocessing','bemobil','APP', ...
        %         strcat(subjid,'_postICA.set'));

        EEG = pop_loadset(filename);

        % segment
        trials = [];
        Info = []; % 1st - trial type, 2nd - block #, 3rd - trial #

        event = EEG.event;

        for e = 1:length(event)
            if strcmp(event(e).type,'TrialStart')
                start_ind = round(event(e).latency);
                end_ind = start_ind + EEG.srate * trial_len-1;
                current_trial = EEG.data(:,start_ind:end_ind);
                current_trial = current_trial - nanmean(current_trial,2);
                trials = cat(3,trials,current_trial);

                switch event(e).trialtype
                    case 'Control'
                        trialtype = 1;
                    case 'Disc'
                        trialtype = 2;
                    case 'VerticalBar'
                        trialtype = 3;
                    case 'HorizontalBar'
                        trialtype = 4;
                end

                Info = cat(1,Info,[trialtype, event(e).block, event(e).trial]);

            end
        end
        results.trials = trials;
        results.Info = Info;
        results.chanlocs = EEG.chanlocs;
    end

    Info = results.Info;

    %% visual inspection to reject trials
    %     filename = fullfile(result_fold,strcat(subjid,'_goodTrials.mat'));

    if ~isfield(results,'good_trials')

        good_trials = ones(size(trials,3),1);
        figure
        for i = 1:size(trials,3)
            plot(squeeze(trials(:,:,i))');
            title(strcat('Trial ', num2str(i)));
            %         answer = questdlg('Accept or reject trial','Visual Inspection','Accept','Reject','Accept');
            answer = input('Accept (1) or reject (0)? ');
            %         if strcmp(answer,'Reject')
            %             good_trials(i) = 0;
            %         end
            if ~answer
                good_trials(i) = 0;
            end
        end

        %         save(filename,'good_trials');

        good_trials = logical(good_trials);
        results.good_trials = good_trials;
    end
    good_trials = logical(results.good_trials);
    Info = Info(good_trials,:);

    %% select electrode subgroups
    if ~isfield(results,'elec_select')

        elec_select = [];
        for i = 1:length(results.chanlocs)
            if ismember(results.chanlocs(i).labels,elec_labels)
                elec_select = [elec_select,i];
            end
        end

        results.elec_select = elec_labels;
    end

    %% decoding in time domain
    %     temp_trials = trials(:,1:250,:);
    %     [results.VEP_acc_time,results.VEP_dprime_time] = decoding_VEP_v3(temp_trials,Info,crossval_method,kfold);


    %% decoding in frequency domain
    if ~isfield(results,'trials_PSD')
        trials_PSD = [];
        trials_PSD_norm = [];

        SR = EEG.srate;
        freqs = 1:40;

        % compute PSD
        for i = 1:size(trials,3)
            temp_trials = squeeze(trials(:,:,i));
            trials_PSD = cat(3,trials_PSD,pwelch(temp_trials',SR,SR/2,freqs,SR)');
        end

        trials_PSD_select = trials_PSD(elec_select,:,good_trials);

        % normalize across trials
        %     cates = unique(Info(:,1));
        %     for i = 1:length(cates)
        %         num_trials = sum(Info(:,1)==cates(i));
        %         mu = mean(trials_PSD(:,:,Info(:,1)==cates(i)),3);
        %         trials_PSD_norm(:,:,Info(:,1)==cates(i)) = trials_PSD(:,:,Info(:,1)==cates(i)) ./ repmat(mu,1,1,num_trials);
        %     end

        % normalize across chans
        %     mu = mean(trials_PSD_select,1);
        %     sig = std(trials_PSD_select,[],1);
        %     trials_PSD_norm = (trials_PSD_select - repmat(mu,size(trials_PSD_select,1),1,1)) ./ repmat(sig,size(trials_PSD_select,1),1,1);

        % convert to dB
        trials_dB = 10*log10(trials_PSD_select);

        results.trials_PSD = trials_PSD;
        results.trials_dB = trials_dB;
    end

    trials_dB = results.trials_dB;

    %% run decoding
    if ~isfield(results,'VEP_PSD_acc')
        disp('Running decoding')
        [results.VEP_PSD_acc,results.VEP_PSD_dprime] = decoding_VEP_v3(results.trials_dB,results.Info(results.good_trials,:),crossval_method,kfold);
    end

    if ~isfield(results,'VEP_PSD_alpha_acc')
        disp('Running alpha decoding')
        [results.VEP_PSD_alpha_acc,results.VEP_PSD_alpha_dprime] = ...
        decoding_VEP_v3(results.trials_dB(:,6,:),...
        results.Info(results.Info(good_trials(:,1))),crossval_method,kfold);

    end

    results.VEP_PSD_alpha_acc


    %% run permutation
    if ~isfield(results,'VEP_PSD_acc_perm')
        disp('Running permutation')
        Info = results.Info;
        VEP_PSD_acc_perm = nan(4,num_perm);
        VEP_PSD_dprime_perm = nan(4,num_perm);
        for p = 1:num_perm
            perm_Info = Info(randperm(size(Info,1)),:,:);
            [VEP_PSD_acc_perm(:,p),VEP_PSD_dprime_perm(:,p)] = decoding_VEP_v3(trials_dB,perm_Info(good_trials,1),'kfold',5);
        end

        results.VEP_PSD_acc_perm = VEP_PSD_acc_perm;
        results.VEP_PSD_dprime_perm = VEP_PSD_dprime_perm;
    end

    %% save results
    %     results.trials = trials;
    %     results.good_trials = good_trials;
    %     results.trials_PSD = trials_PSD;
    %     results.trials_PSD_norm = trials_PSD_norm;
    %     results.trials_dB = trials_dB;

    filename = fullfile(result_fold,strcat(subjid,'_results'));
    save(filename,'results');


end