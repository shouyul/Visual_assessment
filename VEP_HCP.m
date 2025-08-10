%% set parameter

clear; clc

main_fold = '/Users/shouyuling/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Postdoc/Projects/Gensight/Data/VEPTask';
data_fold = fullfile(main_fold,'Analysis');
fig_fold = fullfile(main_fold,'Figures');
result_fold = fullfile(main_fold,'Results');

% set subject info
subs = {'HCP001','HCP002','HCP003','HCP004','HCP005'};
%goggle_cond = {'GogglesOFF','GogglesOFF','GogglesOFF','GogglesOFF','GogglesOFF'};
goggle_cond = {'GogglesON/12Hz','GogglesON/12Hz','GogglesON/12Hz','GogglesON/12Hz','GogglesON/12Hz'};

% set parameters for decoding
trial_len = 10;
crossval_method = 'leaveout';
kfold = 5;
num_perm = 100;
srate = 250;
pnts = 2500;

% script action parameters
overwrite_decoding = 0;
overwrite_permutation = 0;
recompute_power = 1;
run_permutation = 0;
resegmenting = 0;
redo_visual_rej = 0;

%% 
for sub = 5:length(subs)
    subjid = subs{sub};
    disp(subjid);

    %%
    filename = fullfile(result_fold,strcat(subjid,'_results.mat')); % results2, results presented at gensight meeting
    if exist(filename,'file')
        load(filename);
    else
        results = struct();
        results.Subjid = subjid;
        results.GogglesCondition = goggle_cond{sub};
    end

    %% load data
    if ~isfield(results,'trials') || resegmenting
        disp('Segmenting data');
        filename = fullfile(data_fold,'2_preprocessing','bemobil','APP','Brain', ...
            strcat(subjid,'_cleaned_with_ICA.set'));
        EEG = pop_loadset(filename);

        % segment
        trials = [];
        Info = []; % 1st - trial type, 2nd - block #, 3rd - trial #

        EEG_chans = strcmpi({EEG.chanlocs.type},'EEG');
        data = EEG.data(EEG_chans,:);
        chanlocs = EEG.chanlocs(EEG_chans);

        event = EEG.event;
        for e = 1:length(event)
            if strcmp(event(e).type,'TrialStart')
                start_ind = round(event(e).latency);
                end_ind = start_ind + srate * trial_len-1;
                if end_ind > size(data,2)
                    break
                end
                current_trial = data(:,start_ind:end_ind);
                %                 current_trial = current_trial - nanmean(current_trial,2);
                current_trial = detrend(current_trial')';
                trials = cat(3,trials,current_trial);

                if isfield(event,'trialtype')
                    switch event(e).trialtype
                        case 'Control'
                            trialtype = 1;
                        case 'Disc'
                            trialtype = 2;
                        case {'90Bar','VerticalBar'}
                            trialtype = 3;
                        case {'0Bar','HorizontalBar'}
                            trialtype = 4;
                        otherwise
                            trialtype = 0;
                    end
                elseif isfield(event,'stimulus')
                    switch event(e).stimulus
                        case 'Control'
                            trialtype = 1;
                        case {'Disc','3Rings','5Rings'}
                            trialtype = 2;
                        case {'90Bar','-45Bar'}
                            trialtype = 3;
                        case {'0Bar','45Bar'}
                            trialtype = 4;
                        otherwise
                            trialtype = 0;
                    end
                end

                Info = cat(1,Info,[trialtype, event(e).block, event(e).trial]);

            end
        end

        results.trials = trials;
        results.Info = Info;
        results.chanlocs = chanlocs;

    end

    %% visual inspection to reject trials

    if ~isfield(results,'good_trials') || redo_visual_rej

        filename = fullfile(result_fold,strcat(results.Subjid,'_goodTrials.mat'));
        if exist(filename,'file')
            load(filename);
        else
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

            save(filename,'good_trials');
        end
        good_trials = logical(good_trials);
        results.good_trials = good_trials;
        results.trials = results.trials(:,:,good_trials);
        results.Info = results.Info(good_trials,:);
    end

    % good_trials = logical(ones(size(results.trials,3),1));
    % results.good_trials = good_trials;

    %% convert data to power
    if ~isfield(results,'trial_power') || recompute_power
        disp('Wavelet transformation')
        [results.trial_power,results.power_time,results.dB_time, results.freq] = convert_to_power_VEP(results.trials);
        % [results.trial_dB, results.freq, results.spect] = convert_to_power_VEP_v2(results.trials);
    end

    % filename = fullfile(result_fold,strcat(subjid,'_results.mat')); % results2, results presented at gensight meeting
    % save(filename,'results');
    %% plot power
    % 
    % figure; hold on
    % 
    % times = 0:1/srate:10-1/srate;
    % contourf(times,results.freq,squeeze(mean(results.dB_time,1)),40,'linecolor','none');
    % % clour_lim = 3;
    % 
    % colormap(ft_colormap('-RdBu'));
    % % set(gca,'clim',[-clour_lim clour_lim],'xlim',[0 10])
    % set(gca,'xlim',[0 10])
    % set(gca,'ylim',[1,40])


    %% run decoding
    disp('Running decoding')

    fieldname = 'all_1to40';
    if ~isfield(results,fieldname) || overwrite_decoding
        X = results.trial_power(:,:,:);
        y = results.Info(:,1);

        acc = decoding_VEP_v3(X,y,crossval_method,kfold);
        eval(['results.',fieldname,'=acc;']);
    end

    results

    filename = fullfile(result_fold,strcat(subjid,'_results.mat'));
    save(filename,'results');
end

%% collate data
Allsub_results = struct();
trial_type = {'Control','Disc','VerticalBar','HorizontalBar','All'};

cnt = 1;
for sub = 1:length(subs)
    subjid = subs{sub};
    disp(subjid);

    filename = fullfile(result_fold,strcat(subjid,'_results.mat'));
    load(filename);
    
    temp_acc = results.all_1to40;
    temp_acc(5) = mean(temp_acc);

    for a = 1:length(temp_acc)

    
    Allsub_results(cnt).Subjid = subjid;
    Allsub_results(cnt).Task = 'VEP';
    Allsub_results(cnt).ElecGroup = 'All_electrodes';
    Allsub_results(cnt).Condition = '12Hz';
    Allsub_results(cnt).FreqGroup = '1to40';
    Allsub_results(cnt).Stimuli = trial_type{a};

    Allsub_results(cnt).Accuracy = temp_acc(a)*100;
    cnt = cnt+1;
    end
end
% 
% Allsub_results(cnt).Subjid = 'TEST';
% Allsub_results(cnt).Task = 'VEP';
% Allsub_results(cnt).ElecGroup = 'All_electrodes';
% Allsub_results(cnt).Condition = '6Hz';
% Allsub_results(cnt).FreqGroup = '1to40';
% Allsub_results(cnt).Stimuli = trial_type{a};
% 
% Allsub_results(cnt).Accuracy = NaN;


filename = fullfile(result_fold,'Allsub_results_HCP.xlsx');
writetable(struct2table(Allsub_results),filename,'writemode','overwrite')



%% plot
power_time = [];
dB_time = [];
for sub = 1:length(subs)
    subjid = subs{sub};
    disp(subjid);

    filename = fullfile(result_fold,strcat(subjid,'_results.mat'));
    load(filename);

    power_time(:,:,:,sub) = results.power_time;
    dB_time(:,:,:,sub) = results.dB_time;

end

power_time_avg = mean(power_time,4);
dB_time_avg = mean(dB_time,4);


%%
figure; hold on

for i = 1:length(occi_elec_ind)
subplot(4,7,i); hold on
times = 0:1/srate:10-1/srate;
freq = results.freq;
contourf(times,freq(1:50),squeeze(mean(power_time(occi_elec_ind,1:50,:,1),1)),40,'linecolor','none');
clour_lim = 0.4;

% colormap(ft_colormap('-RdBu'));
% set(gca,'clim',[-clour_lim clour_lim])
% set(gca,'clim',[0.5,1])
set(gca,'xlim',[0 10])
set(gca,'ylim',[4,40])

plot([0,10],[6,6],'--k');


end


%%
figure; hold on
baseline = mean(power_time(:,:,:,3),3);
plot(freq,baseline)
plot([6,6],[0,7])
plot([6,6]*2,[0,7])