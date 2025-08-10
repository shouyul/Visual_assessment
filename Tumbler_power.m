%%
clear all; clc
main_fold = '/Users/shouyuling/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Postdoc/Projects/Gensight/Data/TumblerTask';
data_fold = fullfile(main_fold,'Analysis');
fig_fold = fullfile(main_fold,'Figures');
result_fold = fullfile(main_fold,'Results');

subs = {'004-4006-s1','004-4010-s1','004-4010-s2','P1001-4','P1002-2','P1002-3','P1004-2'};
% Pitt patients, block1-3 GogglesON, block4-6 GogglesOFF
% P1002-3, block1-4 GogglesOFF, block5-8, GogglesON
% Remaining Paris patients, block1-4 GogglesON, block5-8, GogglesOFF

srate = 250;
resegmenting = 0;
recompute_power = 0;
%% segment into trials
for sub = 1:length(subs)
    subjid = subs{sub};
    disp(subjid);

    %%
    filename = fullfile(result_fold,strcat(subjid,'_results.mat'));
    if exist(filename,'file')
        load(filename);
    else
        results = struct();
        results.Subjid = subjid;
        results.Task = 'Tumbler';
    end

    %% load data
    if ~isfield(results,'trials') || resegmenting
        disp('Segmenting data');
        filename = fullfile(data_fold,'2_preprocessing','bemobil','APP','Brain', ...
            strcat(subjid,'_cleaned_with_ICA.set'));
        EEG = pop_loadset(filename);

        % get EEG channales only
        EEG_chans = strcmpi({EEG.chanlocs.type},'eeg');
        data = EEG.data(EEG_chans,:);


        % get trial type
        filename = fullfile(data_fold,"0_raw-data/",subjid,[subjid,'_TrialSequence.txt']);
        trialtype = dlmread(filename);
        trialtype = 1-trialtype; % 1-object present, 0-object absent

        trials_ON_open = [];
        trials_ON_close = [];

        trials_OFF_open = [];
        trials_OFF_close = [];

        info_ON_close = [];
        info_OFF_close = [];
        info_ON_open = [];
        info_OFF_open = [];
% 
%         if sub < 4
%             info_ON = nan(3,10);
%             info_OFF = nan(3,10);
%         else
%             info_ON = nan(4,10);
%             info_OFF = nan(4,10);
%         end

        event = EEG.event;

        % EyesClosed period - 'NewTrialSoundoff' to 'OpenEyesSoundOn' - 5s
        % EyesOpen peried - 'OpenEyesSoundOff' to 'ClosedEyesSoundOn' - 15s
        close_length = 5*EEG.srate;
        open_length = 15*EEG.srate;
        for e = 1:length(event)
            start_ind = [];
            end_ind = [];

            if strcmp(event(e).type,'NewTrialSoundOff')
                start_ind = round(event(e).latency);
                end_ind = start_ind + close_length-1;
                trial_close = data(:,start_ind:end_ind);

                if sub < 4
                    switch event(e).block
                        case {1,2,3}
                            trials_ON_close = cat(3,trials_ON_close,trial_close);
                            info_ON_close = cat(1,info_ON_close,trialtype(event(e).block,event(e).trial));
                        case {4,5,6}
                            trials_OFF_close = cat(3,trials_OFF_close,trial_close);
                            info_OFF_close = cat(1,info_OFF_close,trialtype(event(e).block,event(e).trial));
                    end

                elseif sub == 6
                    switch event(e).block
                        case {1,2,3,4}
                            trials_OFF_close = cat(3,trials_OFF_close,trial_close);
                            info_OFF_close = cat(1,info_OFF_close,trialtype(event(e).block,event(e).trial));
                        case {5,6,7,8}
                            trials_ON_close = cat(3,trials_ON_close,trial_close);
                            info_ON_close = cat(1,info_ON_close,trialtype(event(e).block,event(e).trial));
                    end
                else
                    switch event(e).block
                        case {1,2,3,4}
                            trials_ON_close = cat(3,trials_ON_close,trial_close);
                            info_ON_close = cat(1,info_ON_close,trialtype(event(e).block,event(e).trial));
                        case {5,6,7,8}
                            trials_OFF_close = cat(3,trials_OFF_close,trial_close);
                            info_OFF_close = cat(1,info_OFF_close,trialtype(event(e).block,event(e).trial));
                    end
                end

            elseif strcmp(event(e).type,'OpenEyesSoundOff')
                start_ind = round(event(e).latency);
                end_ind = start_ind + open_length-1;
                trial_open = data(:,start_ind:end_ind);

                if sub < 4
                    switch event(e).block
                        case {1,2,3}
                            trials_ON_open = cat(3,trials_ON_open,trial_open);
                            info_ON_open = cat(1,info_ON_open,trialtype(event(e).block,event(e).trial));
                        case {4,5,6}
                            trials_OFF_open = cat(3,trials_OFF_open,trial_open);
                            info_OFF_open = cat(1,info_OF_open,trialtype(event(e).block,event(e).trial));
                    end
                elseif sub == 6
                    switch event(e).block
                        case {1,2,3,4}
                            trials_OFF_open = cat(3,trials_OFF_open,trial_open);
                            info_OFF_open = cat(1,info_OFF_open,trialtype(event(e).block,event(e).trial));
                        case {5,6,7,8}
                            trials_ON_open = cat(3,trials_ON_open,trial_open);
                            info_ON_open = cat(1,info_ON_open,trialtype(event(e).block,event(e).trial));
                    end
                else
                    switch event(e).block
                        case {1,2,3,4}
                            trials_ON_open = cat(3,trials_ON_open,trial_open);
                            info_ON_open = cat(1,info_ON_open,trialtype(event(e).block,event(e).trial));
                        case {5,6,7,8}
                            trials_OFF_open = cat(3,trials_OFF_open,trial_open);
                            info_OFF_open = cat(1,info_OFF_open,trialtype(event(e).block,event(e).trial));
                    end
                end
            end

        end
    end

    results.trials.GogglesON_EyesClosed = trials_ON_close;
    results.trials.GogglesON_EyesOpen = trials_ON_open;
    results.trials.GogglesOFF_EyesClosed = trials_OFF_close;
    results.trials.GogglesOFF_EyesOpen = trials_OFF_open;

    results.info.info_GogglesON_EyesClosed = info_ON_close;
    results.info.info_GogglesON_EyesOpen = info_ON_open;
    results.info.info_GogglesOFF_EyesClosed = info_OFF_close;
    results.info.info_GogglesOFF_EyesOpen = info_OFF_open;
    %% convert data to power
    if ~isfield(results,'trial_power') || recompute_power
        disp('Wavelet transformation: GogglesON - EyesClosed')
        tic
        [results.trial_power.GogglesON_EyesClosed,...
            results.dB.GogglesON_EyesClosed.ObjectPresent,...
            results.dB.GogglesON_EyesClosed.ObjectAbsent,...
            freq_close] = ...
            convert_to_power_Tumbler(results.trials.GogglesON_EyesClosed,results.info.GogglesON);
        toc

        disp('Wavelet transformation: GogglesON - EyesOpen')
        tic
        [results.trial_power.GogglesON_EyesOpen,...
            results.dB.GogglesON_EyesOpen.ObjectPresent,...
            results.dB.GogglesON_EyesOpen.ObjectAbsent,...
            freq_open] = ...
            convert_to_power_Tumbler(results.trials.GogglesON_EyesOpen,results.info.GogglesON);
        toc

        disp('Wavelet transformation: GogglesOFF - EyesClosed')
        tic
        [results.trial_power.GogglesOFF_EyesClosed,...
            results.dB.GogglesOFF_EyesClosed.ObjectPresent,...
            results.dB.GogglesOFF_EyesClosed.ObjectAbsent] = ...
            convert_to_power_Tumbler(results.trials.GogglesOFF_EyesClosed,results.info.GogglesOFF);
        toc

        disp('Wavelet transformation: GogglesOFF - EyesOpen')
        tic
        [results.trial_power.GogglesOFF_EyesOpen,...
            results.dB.GogglesOFF_EyesOpen.ObjectPresent,...
            results.dB.GogglesOFF_EyesOpen.ObjectAbsent] = ...
            convert_to_power_Tumbler(results.trials.GogglesOFF_EyesOpen,results.info.GogglesOFF);
        toc
    end

    results.freq.freq_close = freq_close;
    results.freq.freq_open = freq_open;

    filename = fullfile(result_fold,strcat(subjid,'_results.mat'));
    save(filename,"results",'-v7.3');

end

%% aggregate data across subjects

for sub = 1:length(subs)
    subjid = subs{sub};
    disp(subjid);

    %%
    filename = fullfile(result_fold,strcat(subjid,'_results.mat'));
    if exist(filename,'file')
        load(filename);
    else
        results = struct();
        results.Subjid = subjid;
        results.Task = 'Tumbler';
    end

dB_ON_close_present(:,:,sub) = squeeze(mean(results.dB.GogglesON_EyesClosed.ObjectPresent));
dB_ON_close_absent(:,:,sub)  = squeeze(mean(results.dB.GogglesON_EyesClosed.ObjectAbsent));

dB_ON_open_present(:,:,sub)  = squeeze(mean(results.dB.GogglesON_EyesOpen.ObjectPresent));
dB_ON_open_absent(:,:,sub)  = squeeze(mean(results.dB.GogglesON_EyesOpen.ObjectAbsent));

dB_OFF_close_present(:,:,sub)  = squeeze(mean(results.dB.GogglesOFF_EyesClosed.ObjectPresent));
dB_OFF_close_absent(:,:,sub)  = squeeze(mean(results.dB.GogglesOFF_EyesClosed.ObjectAbsent));

dB_OFF_open_present(:,:,sub)  = squeeze(mean(results.dB.GogglesOFF_EyesOpen.ObjectPresent));
dB_OFF_open_absent(:,:,sub)  = squeeze(mean(results.dB.GogglesOFF_EyesOpen.ObjectAbsent));

end

freq_open = results.freq.freq_open;
freq_close = results.freq.freq_close;

%%
dB_ON_close_present = mean(dB_ON_close_present,3);
dB_ON_close_absent = mean(dB_ON_close_absent,3);
dB_ON_open_present = mean(dB_ON_open_present,3);
dB_ON_open_absent = mean(dB_ON_open_absent,3);


dB_OFF_close_present = mean(dB_OFF_close_present,3);
dB_OFF_close_absent = mean(dB_OFF_close_absent,3);
dB_OFF_open_present = mean(dB_OFF_open_present,3);
dB_OFF_open_absent = mean(dB_OFF_open_absent,3);

freq_open = results.freq.freq_open;
freq_close = results.freq.freq_close;

%
dB_ON_close_present = mean(dB_ON_close_present,3);
dB_ON_close_absent = mean(dB_ON_close_absent,3);
dB_ON_open_present = mean(dB_ON_open_present,3);
dB_ON_open_absent = mean(dB_ON_open_absent,3);


dB_OFF_close_present = mean(dB_OFF_close_present,3);
dB_OFF_close_absent = mean(dB_OFF_close_absent,3);
dB_OFF_open_present = mean(dB_OFF_open_present,3);
dB_OFF_open_absent = mean(dB_OFF_open_absent,3);

%% plot
% set up parameters
srate = 250;
times_close = [0:1/srate:5-1/srate];
times_open = [0:1/srate:15-1/srate];

clour_lim = 0.5;
fig_visibility = 'on';
fig_size = [ 168   281   828   622];
font_size = 16;

% EyesClosed
figure('visible',fig_visibility); hold on
set(gcf,'Position',fig_size)

% EyesClosed Object present
subplot(2,2,1); hold on
fig = contourf(times_close,freq_close,dB_ON_close_present,40,'linecolor','none');
colormap(ft_colormap('-RdBu'));
set(gca,'clim',[-clour_lim clour_lim],'xlim',[0 5])
set(gca,'ylim',[1,40])
set(gca,'FontSize',font_size)
xlabel('Time (s)')
ylabel('Frequency (Hz)')
title('GogglesON - ObjectPresent')

% EyesClosed Object aresent
subplot(2,2,2); hold on
fig = contourf(times_close,freq_close,dB_ON_close_absent,40,'linecolor','none');
colormap(ft_colormap('-RdBu'));
set(gca,'clim',[-clour_lim clour_lim],'xlim',[0 5])
set(gca,'ylim',[1,40])
set(gca,'FontSize',font_size)
xlabel('Time (s)')
ylabel('Frequency (Hz)')
title('GogglesON - ObjectAbsent')

% EyesClosed Object present
subplot(2,2,3); hold on
fig = contourf(times_close,freq_close,dB_OFF_close_present,40,'linecolor','none');
colormap(ft_colormap('-RdBu'));
set(gca,'clim',[-clour_lim clour_lim],'xlim',[0 5])
set(gca,'ylim',[1,40])
set(gca,'FontSize',font_size)
xlabel('Time (s)')
ylabel('Frequency (Hz)')
title('GogglesOFF - ObjectPresent')

% EyesClosed Object aresent
subplot(2,2,4); hold on
fig = contourf(times_close,freq_close,dB_OFF_close_absent,40,'linecolor','none');
colormap(ft_colormap('-RdBu'));
set(gca,'clim',[-clour_lim clour_lim],'xlim',[0 5])
set(gca,'ylim',[1,40])
set(gca,'FontSize',font_size)
xlabel('Time (s)')
ylabel('Frequency (Hz)')
title('GogglesOFF - ObjectAbsent')


sgtitle('EyesClosed','Fontsize',24)



% EyesOpen
figure('visible',fig_visibility); hold on
set(gcf,'Position',fig_size)

% GogglesON Object present
subplot(2,2,1); hold on
fig = contourf(times_open,freq_open,dB_ON_open_present,40,'linecolor','none');
colormap(ft_colormap('-RdBu'));
set(gca,'clim',[-clour_lim clour_lim],'xlim',[0 15])
set(gca,'ylim',[1,40])
set(gca,'FontSize',font_size)
xlabel('Time (s)')
ylabel('Frequency (Hz)')
title('GogglesON - ObjectPresent')

% EyesOpen Object aresent
subplot(2,2,2); hold on
fig = contourf(times_open,freq_open,dB_ON_open_absent,40,'linecolor','none');
colormap(ft_colormap('-RdBu'));
set(gca,'clim',[-clour_lim clour_lim],'xlim',[0 15])
set(gca,'ylim',[1,40])
set(gca,'FontSize',font_size)
xlabel('Time (s)')
ylabel('Frequency (Hz)')
title('GogglesON - ObjectAbsent')

% EyesOpen Object present
subplot(2,2,3); hold on
fig = contourf(times_open,freq_open,dB_OFF_open_present,40,'linecolor','none');
colormap(ft_colormap('-RdBu'));
set(gca,'clim',[-clour_lim clour_lim],'xlim',[0 15])
set(gca,'ylim',[1,40])
set(gca,'FontSize',font_size)
xlabel('Time (s)')
ylabel('Frequency (Hz)')
title('GogglesOFF - ObjectPresent')

% EyesOpen Object aresent
subplot(2,2,4); hold on
fig = contourf(times_open,freq_open,dB_OFF_open_absent,40,'linecolor','none');
colormap(ft_colormap('-RdBu'));
set(gca,'clim',[-clour_lim clour_lim],'xlim',[0 15])
set(gca,'ylim',[1,40])
set(gca,'FontSize',font_size)
xlabel('Time (s)')
ylabel('Frequency (Hz)')
title('GogglesOFF - ObjectAbsent')


sgtitle('EyesOpen','Fontsize',24)


%% plot dB difference

% EyesClosed
% ObjectPresent - ObjectAbsent

diff_ON_close = mean(dB_ON_close_present - dB_ON_close_absent,3);
diff_OFF_close = mean(dB_OFF_close_present - dB_OFF_close_absent,3);

% EyesOpen
% ObjectPresent - ObjectAbsent

diff_ON_open = mean(dB_ON_open_present - dB_ON_open_absent,3);
diff_OFF_open = mean(dB_OFF_open_present - dB_OFF_open_absent,3);

% plot
% set up parameters
srate = 250;
times_close = [0:1/srate:5-1/srate];
times_open = [0:1/srate:15-1/srate];

clour_lim = 1;
fig_visibility = 'off';
fig_size = [ 168   281   828   622];
font_size = 16;

% EyesClosed
figure('visible',fig_visibility); hold on
set(gcf,'Position',fig_size)

subplot(2,2,1); hold on
fig = contourf(times_close,freq_close,diff_ON_close,40,'linecolor','none');
colormap(ft_colormap('-RdBu'));
set(gca,'clim',[-clour_lim clour_lim],'xlim',[0 5])
set(gca,'ylim',[1,40])
set(gca,'FontSize',font_size)
xlabel('Time (s)')
ylabel('Frequency (Hz)')
title('EyesClosed - GogglesON')

subplot(2,2,2); hold on
fig = contourf(times_close,freq_close,diff_OFF_close,40,'linecolor','none');
colormap(ft_colormap('-RdBu'));
set(gca,'clim',[-clour_lim clour_lim],'xlim',[0 5])
set(gca,'ylim',[1,40])
set(gca,'FontSize',font_size)
xlabel('Time (s)')
ylabel('Frequency (Hz)')
title('EyesClosed - GogglesOFF')


subplot(2,2,3); hold on
fig = contourf(times_open,freq_open,diff_ON_open,40,'linecolor','none');
colormap(ft_colormap('-RdBu'));
set(gca,'clim',[-clour_lim clour_lim],'xlim',[0 15])
set(gca,'ylim',[1,40])
set(gca,'FontSize',font_size)
xlabel('Time (s)')
ylabel('Frequency (Hz)')
title('EyesOpen - GogglesON')

subplot(2,2,4); hold on
fig = contourf(times_open,freq_open,diff_OFF_open,40,'linecolor','none');
colormap(ft_colormap('-RdBu'));
set(gca,'clim',[-clour_lim clour_lim],'xlim',[0 15])
set(gca,'ylim',[1,40])
set(gca,'FontSize',font_size)
xlabel('Time (s)')
ylabel('Frequency (Hz)')
title('EyesOpen - GogglesOFF')

sgtitle('ObjectPresent vs ObjectAbsent','Fontsize',24)

filename = fullfile(fig_fold,'dB_diff_PreVSAbs.png');
saveas(gcf,filename);

%% plot dB difference

% GogglesON - GogglesOFF
% EyesClosed

dif_close_present = mean(dB_ON_close_present - dB_OFF_close_present,3);
diff_close_absent = mean(dB_ON_close_absent - dB_OFF_close_absent,3);

% EyesOpen

diff_open_present = mean(dB_ON_open_present - dB_OFF_open_present,3);
diff_open_absent = mean(dB_ON_open_absent - dB_OFF_open_absent,3);


% plot
% set up parameters
srate = 250;
times_close = [0:1/srate:5-1/srate];
times_open = [0:1/srate:15-1/srate];

clour_lim = 0.5;
fig_visibility = 'off';
fig_size = [ 168   281   828   622];
font_size = 16;

% EyesClosed
figure('visible',fig_visibility); hold on
set(gcf,'Position',fig_size)

subplot(2,2,1); hold on
fig = contourf(times_close,freq_close,dif_close_present,40,'linecolor','none');
colormap(ft_colormap('-RdBu'));
set(gca,'clim',[-clour_lim clour_lim],'xlim',[0 5])
set(gca,'ylim',[1,40])
set(gca,'FontSize',font_size)
xlabel('Time (s)')
ylabel('Frequency (Hz)')
title('EyesClosed - ObjectPresent')

subplot(2,2,2); hold on
fig = contourf(times_close,freq_close,diff_close_absent,40,'linecolor','none');
colormap(ft_colormap('-RdBu'));
set(gca,'clim',[-clour_lim clour_lim],'xlim',[0 5])
set(gca,'ylim',[1,40])
set(gca,'FontSize',font_size)
xlabel('Time (s)')
ylabel('Frequency (Hz)')
title('EyesClosed - ObjectAbsent')


subplot(2,2,3); hold on
fig = contourf(times_open,freq_open,diff_open_present,40,'linecolor','none');
colormap(ft_colormap('-RdBu'));
set(gca,'clim',[-clour_lim clour_lim],'xlim',[0 15])
set(gca,'ylim',[1,40])
set(gca,'FontSize',font_size)
xlabel('Time (s)')
ylabel('Frequency (Hz)')
title('EyesOpen - ObjectPresent')

subplot(2,2,4); hold on
fig = contourf(times_open,freq_open,diff_open_absent,40,'linecolor','none');
colormap(ft_colormap('-RdBu'));
set(gca,'clim',[-clour_lim clour_lim],'xlim',[0 15])
set(gca,'ylim',[1,40])
set(gca,'FontSize',font_size)
xlabel('Time (s)')
ylabel('Frequency (Hz)')
title('EyesOpen - ObjectAbsent')

sgtitle('GogglesON vs GogglesOFF','Fontsize',24)

filename = fullfile(fig_fold,'dB_diff_ONvsOFF.png');
saveas(gcf,filename);