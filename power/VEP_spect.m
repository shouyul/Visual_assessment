%% set parameter

clear; clc

main_fold = '/Users/shouyuling/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Postdoc/Projects/Gensight/Data/VEPTask';
data_fold = fullfile(main_fold,'Analysis');
fig_fold = fullfile(main_fold,'Figures');
result_fold = fullfile(main_fold,'Results');
cd(result_fold)

% variables for plotting
srate = 250;
times = [0:1/srate:10-1/srate];

%% get chanlocs
filename = fullfile(result_fold,'chanlocs.xlsx');
chanlocs = readtable(filename);

%% HCP SSVEP

subs_HCP =  {'HCP001','HCP002','HCP003','HCP004','HCP005'};

for sub = 1:length(subs_HCP)
    subjid = subs_HCP{sub};
    disp(subjid)

    filename = fullfile(result_fold,strcat(subjid,'_results.mat'));
    load(filename);
    
    sub_chanlocs = results.chanlocs;
    trials = results.trials;

    front_ind = []; tempo_ind = []; parie_ind = []; occip_ind = [];
    eeg_ind = [];

    for ch = 1:length(sub_chanlocs)
        chan_ind = find(strcmpi(chanlocs.labels,sub_chanlocs(ch).labels));
        chan_loc = chanlocs(chan_ind,:);
        
        switch chan_loc.Location{1}
            case 'Frontal'
                front_ind = [front_ind,ch];
            case 'Temporal'
                tempo_ind = [tempo_ind,ch];
            case 'Parietal'
                parie_ind = [parie_ind,ch];
            case 'Occipital'
                occip_ind = [occip_ind,ch];
        end

        if strcmpi(sub_chanlocs(ch).type,'EEG')
            eeg_ind = [eeg_ind,ch];
        end

    end

    x_psd = [];
    for ch = 1:size(trials,1)
        for t = 1:size(trials,3)
            x = squeeze(trials(ch,:,t));
            x_fft = fft(x);
            x_fft = x_fft(1:length(x)/2+1);
            freq = 0:srate/length(x):srate/2;
            temp_psd = (1/(srate*length(x))) * abs(x_fft).^2;
            temp_psd(2:end-1) = 2*temp_psd(2:end-1);
            x_psd(ch,:,t) = temp_psd;
        end
    end

    x_psd_base = squeeze(mean(x_psd(eeg_ind,:,results.Info(:,1)==1),3));
    x_psd_stim = squeeze(mean(x_psd(eeg_ind,:,results.Info(:,1)~=1),3));


    x_psd_ratio = [];
    for i = 1:size(x_psd_base,1)
        x_psd_ratio(i,:) = x_psd_stim(i,:) ./ x_psd_base(i,:);
    end

    psd_allSub(:,sub) = mean(x_psd_ratio,1);
    psd_base_allSub(:,sub) = mean(x_psd_base,1);
    psd_stim_allSub(:,sub) = mean(x_psd_stim,1);
end



% average across subjects and plot
psd_HCP = mean(psd_allSub,2);
psd_HCP_base = mean(psd_base_allSub,2);
psd_HCP_stim = mean(psd_stim_allSub,2);


figure; hold on
set(gcf,'Position',[1440         475         409         763]);

subplot(2,1,1); hold on
set(gca,'FontSize',14)

plot(freq,psd_HCP_stim)
plot(freq,psd_HCP_base)
xlim([1,40])
xlabel('Frequency (Hz)')
ylabel('Power (dB)')

y_min = min(min(psd_HCP_stim(11:400)),min(psd_HCP_base(11:400))) * 0.9;
y_max = max(max(psd_HCP_stim(11:400)),max(psd_HCP_base(11:400))) * 1.1;
ylim([y_min,3.5])
legend({'Stim','Baseline'},'FontSize',12)

subplot(2,1,2); hold on
set(gca,'FontSize',14)

plot(freq,psd_HCP);
xlim([1,40])

y_min = min(psd_HCP(:)) * 0.9;
y_max = max(psd_HCP(:)) * 1.1;
ylim([y_min,y_max])
xlabel('Frequency (Hz)')
ylabel('Power ratio')

sgtitle('Sighted Control','FontSize',16)



%% BCP 6Hz SSVEP

subs_BCP_6Hz =  {'P1009','004-4012-s1','BCP003_6Hz','BCP004_6Hz','BCP005_6Hz'};

for sub = 1:length(subs_BCP_6Hz)
    subjid = subs_BCP_6Hz{sub};
    disp(subjid)

    filename = fullfile(result_fold,strcat(subjid,'_results.mat'));
    load(filename);
    
    sub_chanlocs = results.chanlocs;
    trials = results.trials;

    front_ind = []; tempo_ind = []; parie_ind = []; occip_ind = [];
    eeg_ind = [];

    for ch = 1:length(sub_chanlocs)
        chan_ind = find(strcmpi(chanlocs.labels,sub_chanlocs(ch).labels));
        chan_loc = chanlocs(chan_ind,:);
        
        switch chan_loc.Location{1}
            case 'Frontal'
                front_ind = [front_ind,ch];
            case 'Temporal'
                tempo_ind = [tempo_ind,ch];
            case 'Parietal'
                parie_ind = [parie_ind,ch];
            case 'Occipital'
                occip_ind = [occip_ind,ch];
        end

        if strcmpi(sub_chanlocs(ch).type,'EEG')
            eeg_ind = [eeg_ind,ch];
        end

    end

    x_psd = [];
    for ch = 1:size(trials,1)
        for t = 1:size(trials,3)
            x = squeeze(trials(ch,:,t));
            x_fft = fft(x);
            x_fft = x_fft(1:length(x)/2+1);
            freq = 0:srate/length(x):srate/2;
            temp_psd = (1/(srate*length(x))) * abs(x_fft).^2;
            temp_psd(2:end-1) = 2*temp_psd(2:end-1);
            x_psd(ch,:,t) = temp_psd;
        end
    end

    x_psd_base = squeeze(mean(x_psd(eeg_ind,:,results.Info(:,1)==1),3));
    x_psd_stim = squeeze(mean(x_psd(eeg_ind,:,results.Info(:,1)~=1),3));


    x_psd_ratio = [];
    for i = 1:size(x_psd_base,1)
        x_psd_ratio(i,:) = x_psd_stim(i,:) ./ x_psd_base(i,:);
    end

    psd_allSub(:,sub) = mean(x_psd_ratio,1);
    psd_base_allSub(:,sub) = mean(x_psd_base,1);
    psd_stim_allSub(:,sub) = mean(x_psd_stim,1);
end

% average across subjects and plot
psd_BCP_6Hz = mean(psd_allSub,2);
psd_BCP_6Hz_base = mean(psd_base_allSub,2);
psd_BCP_6Hz_stim = mean(psd_stim_allSub,2);


figure; hold on
set(gcf,'Position',[1440         475         409         763]);

subplot(2,1,1); hold on
set(gca,'FontSize',14)

plot(freq,psd_BCP_6Hz_stim)
plot(freq,psd_BCP_6Hz_base)
xlim([1,40])
xlabel('Frequency (Hz)')
ylabel('Power (dB)')

y_min = min(min(psd_BCP_6Hz_stim(21:400)),min(psd_BCP_6Hz_base(21:400))) * 0.9;
y_max = max(max(psd_BCP_6Hz_stim(21:400)),max(psd_BCP_6Hz_base(21:400))) * 1.1;
ylim([y_min,y_max])
legend({'Stim','Baseline'},'FontSize',12)

subplot(2,1,2); hold on
set(gca,'FontSize',14)

plot(freq,psd_BCP_6Hz);
xlim([1,40])

y_min = min(psd_BCP_6Hz(:)) * 0.9;
y_max = max(psd_BCP_6Hz(:)) * 1.1;
ylim([y_min,y_max])
xlabel('Frequency (Hz)')
ylabel('Power ratio')

sgtitle('Blind Control 6Hz','FontSize',16)


%% BCP 12Hz SSVEP

subs_BCP_12Hz =  {'004-4012-s2','BCP003_12Hz','BCP004_12Hz','BCP005_12Hz'};

for sub = 1:length(subs_BCP_12Hz)
    subjid = subs_BCP_12Hz{sub};
    disp(subjid)

    filename = fullfile(result_fold,strcat(subjid,'_results.mat'));
    load(filename);
    
    sub_chanlocs = results.chanlocs;
    trials = results.trials;

    front_ind = []; tempo_ind = []; parie_ind = []; occip_ind = [];
    eeg_ind = [];

    for ch = 1:length(sub_chanlocs)
        chan_ind = find(strcmpi(chanlocs.labels,sub_chanlocs(ch).labels));
        chan_loc = chanlocs(chan_ind,:);
        
        switch chan_loc.Location{1}
            case 'Frontal'
                front_ind = [front_ind,ch];
            case 'Temporal'
                tempo_ind = [tempo_ind,ch];
            case 'Parietal'
                parie_ind = [parie_ind,ch];
            case 'Occipital'
                occip_ind = [occip_ind,ch];
        end

        if strcmpi(sub_chanlocs(ch).type,'EEG')
            eeg_ind = [eeg_ind,ch];
        end

    end

    x_psd = [];
    for ch = 1:size(trials,1)
        for t = 1:size(trials,3)
            x = squeeze(trials(ch,:,t));
            x_fft = fft(x);
            x_fft = x_fft(1:length(x)/2+1);
            freq = 0:srate/length(x):srate/2;
            temp_psd = (1/(srate*length(x))) * abs(x_fft).^2;
            temp_psd(2:end-1) = 2*temp_psd(2:end-1);
            x_psd(ch,:,t) = temp_psd;
        end
    end

    x_psd_base = squeeze(mean(x_psd(eeg_ind,:,results.Info(:,1)==1),3));
    x_psd_stim = squeeze(mean(x_psd(eeg_ind,:,results.Info(:,1)~=1),3));


    x_psd_ratio = [];
    for i = 1:size(x_psd_base,1)
        x_psd_ratio(i,:) = x_psd_stim(i,:) ./ x_psd_base(i,:);
    end

    psd_allSub(:,sub) = mean(x_psd_ratio,1);
    psd_base_allSub(:,sub) = mean(x_psd_base,1);
    psd_stim_allSub(:,sub) = mean(x_psd_stim,1);
end

% average across subjects and plot
psd_BCP_12Hz = mean(psd_allSub,2);
psd_BCP_12Hz_base = mean(psd_base_allSub,2);
psd_BCP_12Hz_stim = mean(psd_stim_allSub,2);


figure; hold on
set(gcf,'Position',[1440         475         409         763]);

subplot(2,1,1); hold on
set(gca,'FontSize',14)

plot(freq,psd_BCP_12Hz_stim)
plot(freq,psd_BCP_12Hz_base)
xlim([1,40])
xlabel('Frequency (Hz)')
ylabel('Power (dB)')

y_min = min(min(psd_BCP_12Hz_stim(21:400)),min(psd_BCP_12Hz_base(21:400))) * 0.9;
y_max = max(max(psd_BCP_12Hz_stim(21:400)),max(psd_BCP_12Hz_base(21:400))) * 1.1;
ylim([y_min,y_max])
legend({'Stim','Baseline'},'FontSize',12)

subplot(2,1,2); hold on
set(gca,'FontSize',14)

plot(freq,psd_BCP_12Hz);
xlim([1,40])

y_min = min(psd_BCP_12Hz(:)) * 0.9;
y_max = max(psd_BCP_12Hz(:)) * 1.1;
ylim([y_min,y_max])
xlabel('Frequency (Hz)')
ylabel('Power ratio')

sgtitle('Blind Control 12Hz','FontSize',16)


%% Gensight GogglesON 6Hz SSVEP

subs = {'004-4006-s2','004-4006-s4',...
    '004-4010-s3','004-4010-s5',....
    'P1001-2','P1002','P1004-2'};
subs_GenON_6Hz = {'P1001-2','P1002','P1004-2'};

psd_allSub = [];
psd_base_allSub = [];
psd_stim_allSub = [];
for sub = 1:length(subs_GenON_6Hz)
    subjid = subs_GenON_6Hz{sub};
    disp(subjid)

    filename = fullfile(result_fold,strcat(subjid,'_results2.mat'));
    load(filename);
    
    [x_psd_base,x_psd_stim,x_psd_ratio] = compute_PSD_VEP(results);

    psd_allSub(:,sub) = mean(x_psd_ratio,1);
    psd_base_allSub(:,sub) = mean(x_psd_base,1);
    psd_stim_allSub(:,sub) = mean(x_psd_stim,1);
end
%
% average across subjects and plot
psd_GenON_6Hz = mean(psd_allSub,2);
psd_GenON_6Hz_base = mean(psd_base_allSub,2);
psd_GenON_6Hz_stim = mean(psd_stim_allSub,2);


psd_stim = psd_GenON_6Hz_stim;
psd_base = psd_GenON_6Hz_base;
psd_ratio = psd_GenON_6Hz;

figure; hold on
set(gcf,'Position',[1440         475         409         763]);

subplot(2,1,1); hold on
set(gca,'FontSize',14)

plot(freq,psd_stim)
plot(freq,psd_base)
xlim([1,40])
xlabel('Frequency (Hz)')
ylabel('Power (dB)')

y_min = min(min(psd_stim(21:400)),min(psd_base(21:400))) * 0.9;
y_max = max(max(psd_stim(21:400)),max(psd_base(21:400))) * 1.1;
ylim([y_min,y_max])
legend({'Stim','Baseline'},'FontSize',12)

subplot(2,1,2); hold on
set(gca,'FontSize',14)

plot(freq,psd_ratio);
xlim([1,40])

y_min = min(psd_ratio(21:400)) * 0.9;
y_max = max(psd_ratio(21:400)) * 1.1;
ylim([y_min,y_max])
xlabel('Frequency (Hz)')
ylabel('Power ratio')

sgtitle('GenSight Patient GogglesON 6Hz','FontSize',16)


%% Gensight GogglesON 12Hz SSVEP

subs = {'004-4006-s2','004-4006-s4',...
    '004-4010-s3','004-4010-s5',....
    'P1001-2','P1002','P1004-2'};
subs_GenON_12Hz = {'004-4006-s2','004-4006-s4',...
    '004-4010-s3','004-4010-s5'};

psd_allSub = [];
psd_base_allSub = [];
psd_stim_allSub = [];
for sub = 1:length(subs_GenON_12Hz)
    subjid = subs_GenON_12Hz{sub};
    disp(subjid)

    filename = fullfile(result_fold,strcat(subjid,'_results2.mat'));
    load(filename);
    
    [x_psd_base,x_psd_stim,x_psd_ratio] = compute_PSD_VEP(results);

    psd_allSub(:,sub) = mean(x_psd_ratio,1);
    psd_base_allSub(:,sub) = mean(x_psd_base,1);
    psd_stim_allSub(:,sub) = mean(x_psd_stim,1);
end
%
% average across subjects and plot
psd_ratio = mean(psd_allSub,2);
psd_base = mean(psd_base_allSub,2);
psd_stim = mean(psd_stim_allSub,2);


figure; hold on
set(gcf,'Position',[1440         475         409         763]);

subplot(2,1,1); hold on
set(gca,'FontSize',14)

plot(freq,psd_stim)
plot(freq,psd_base)
xlim([1,40])
xlabel('Frequency (Hz)')
ylabel('Power (dB)')

y_min = min(min(psd_stim(21:400)),min(psd_base(21:400))) * 0.9;
y_max = max(max(psd_stim(21:400)),max(psd_base(21:400))) * 1.1;
ylim([y_min,y_max])
legend({'Stim','Baseline'},'FontSize',12)

subplot(2,1,2); hold on
set(gca,'FontSize',14)

plot(freq,psd_ratio);
xlim([1,40])

y_min = min(psd_ratio(21:400)) * 0.9;
y_max = max(psd_ratio(21:400)) * 1.1;
ylim([y_min,y_max])
xlabel('Frequency (Hz)')
ylabel('Power ratio')

sgtitle('GenSight Patient GogglesON 12Hz','FontSize',16)

%% Gensight GogglesOFF 3Hz SSVEP

subs_GenOFF_3Hz = {'004-4006-s2','004-4010-s4','004-4010-s6',...
    'P1001-3','P1002-2'};

subs = subs_GenOFF_3Hz(4:5);

psd_allSub = [];
psd_base_allSub = [];
psd_stim_allSub = [];
for sub = 1:length(subs)
    subjid = subs{sub};
    disp(subjid)

    filename = fullfile(result_fold,strcat(subjid,'_results2.mat'));
    load(filename);
    
    [x_psd_base,x_psd_stim,x_psd_ratio,freq] = compute_PSD_VEP(results);

    psd_allSub(:,sub) = mean(x_psd_ratio,1);
    psd_base_allSub(:,sub) = mean(x_psd_base,1);
    psd_stim_allSub(:,sub) = mean(x_psd_stim,1);
end
%
% average across subjects and plot
psd_ratio = mean(psd_allSub,2);
psd_base = mean(psd_base_allSub,2);
psd_stim = mean(psd_stim_allSub,2);


figure; hold on
set(gcf,'Position',[1440         475         409         763]);

subplot(2,1,1); hold on
set(gca,'FontSize',14)

plot(freq,psd_stim)
plot(freq,psd_base)
xlim([1,40])
xlabel('Frequency (Hz)')
ylabel('Power (dB)')

y_min = min(min(psd_stim(21:400)),min(psd_base(21:400))) * 0.9;
y_max = max(max(psd_stim(21:400)),max(psd_base(21:400))) * 1.1;
ylim([y_min,y_max])
legend({'Stim','Baseline'},'FontSize',12)

subplot(2,1,2); hold on
set(gca,'FontSize',14)

plot(freq,psd_ratio);
xlim([1,40])

y_min = min(psd_ratio(21:400)) * 0.9;
y_max = max(psd_ratio(21:400)) * 1.1;
ylim([y_min,y_max])
xlabel('Frequency (Hz)')
ylabel('Power ratio')

sgtitle('GenSight Patient GogglesOFF 3Hz','FontSize',16)



%% Gensight GogglesOFF 6Hz SSVEP

subs_GenOFF_3Hz = {'004-4006-s2','004-4010-s4','004-4010-s6',...
    'P1001-3','P1002-2'};

subs = subs_GenOFF_3Hz(1:3);

psd_allSub = [];
psd_base_allSub = [];
psd_stim_allSub = [];
for sub = 1:length(subs)
    subjid = subs{sub};
    disp(subjid)

    filename = fullfile(result_fold,strcat(subjid,'_results2.mat'));
    load(filename);
    
    [x_psd_base,x_psd_stim,x_psd_ratio,freq] = compute_PSD_VEP(results);

    psd_allSub(:,sub) = mean(x_psd_ratio,1);
    psd_base_allSub(:,sub) = mean(x_psd_base,1);
    psd_stim_allSub(:,sub) = mean(x_psd_stim,1);
end
%
% average across subjects and plot
psd_ratio = mean(psd_allSub,2);
psd_base = mean(psd_base_allSub,2);
psd_stim = mean(psd_stim_allSub,2);


figure; hold on
set(gcf,'Position',[1440         475         409         763]);

subplot(2,1,1); hold on
set(gca,'FontSize',14)

plot(freq,psd_stim)
plot(freq,psd_base)
xlim([1,40])
xlabel('Frequency (Hz)')
ylabel('Power (dB)')

y_min = min(min(psd_stim(21:400)),min(psd_base(21:400))) * 0.9;
y_max = max(max(psd_stim(21:400)),max(psd_base(21:400))) * 1.1;
ylim([y_min,y_max])
legend({'Stim','Baseline'},'FontSize',12)

subplot(2,1,2); hold on
set(gca,'FontSize',14)

plot(freq,psd_ratio);
xlim([1,40])

y_min = min(psd_ratio(21:400)) * 0.9;
y_max = max(psd_ratio(21:400)) * 1.1;
ylim([y_min,y_max])
xlabel('Frequency (Hz)')
ylabel('Power ratio')

sgtitle('GenSight Patient GogglesOFF 6Hz','FontSize',16)

%% spectrogram
subs_HCP =  {'HCP001','HCP002','HCP003','HCP004','HCP005'};

subs_GenON = {'004-4006-s2','004-4006-s4',...
    '004-4010-s3','004-4010-s5',....
    'P1001-2','P1002','P1004-2'};
subs_GenOFF = {'004-4006-s2','004-4010-s4','004-4010-s6',...
    'P1001-3','P1002-2'};

subs_BCP =  {'004-4012-s2','BCP003_12Hz','BCP004_12Hz','BCP005_12Hz',...
            'P1009','004-4012-s1','BCP003_6Hz','BCP004_6Hz','BCP005_6Hz'};

subs = subs_BCP;

pow_diff = single([]);
pow_base = single([]);
pow_stim = single([]);
for sub = 1:length(subs)
    subjid = subs{sub};
    disp(subjid)

    if subjid(1) == '0' || subjid(1) == 'P'
        filename = fullfile(result_fold,strcat(subjid,'_results2.mat'));
    else
        filename = fullfile(result_fold,strcat(subjid,'_results.mat'));
    end
    load(filename);

    trials = results.trials;
    if size(trials,1) == 128
        trials(88,:,:) = [];
    end
    [~,freq] = cwt(trials(1,:,1),srate);

    for ch = 1:size(trials,1)
        disp([subjid,': ',num2str(ch)]);
        temp_wavelet = nan(length(freq),size(trials,2),size(trials,3));

        parfor t = 1:size(trials,3)
            temp_wavelet(:,:,t) = cwt(trials(ch,:,t),srate);
        end
    
        temppower = abs(temp_wavelet);
        temp_base = mean(temppower(:,:,results.Info(:,1)==1),3);
        temp_stim = mean(temppower(:,:,results.Info(:,1)~=1),3);

        pow_diff(ch,:,:,sub) = single(temp_stim ./ mean(temp_base,2));
        pow_base(ch,:,:,sub) = single(temp_base);
        pow_stim(ch,:,:,sub) = single(temp_stim);
    end
end

save('BCP_wavelet','pow_diff','pow_base','pow_stim');

%% plot

load('GenOFF_wavelet.mat')
pow_base_avg = squeeze(mean(mean(pow_base,4),1));
pow_stim_avg = squeeze(mean(mean(pow_stim,4),1));
pow_diff_avg = squeeze(mean(mean(10*log10(pow_diff(1:126,:,:,:)),4),1));

% clour_lim = 1;
figure; hold on
times = [0:1/srate:10-1/srate];
contourf(times,freq,pow_diff_avg,40,'linecolor','none');
%plot([0,10],[12,12],'--k')
colormap(ft_colormap('-RdBu'))
set(gca,'FontSize',12)
ylim([1,40])
set(gca,'clim',[-0.5 0.5])
xlabel('Time (s)')
ylabel('Frequency (Hz)')
title('GenSight GogglesOFF','FontSize',16);

%% 
GenON = load('GenON_wavelet.mat');
GenOFF = load('GenOFF_wavelet.mat');

GenON_diff = squeeze(mean(mean(GenON.pow_diff(1:126,:,:,:),4),1));
GenOFF_diff = squeeze(mean(mean(GenOFF.pow_diff(1:126,:,:,:),4),1));

Gen_ONvsOFF = 10*log10(GenON_diff ./ GenOFF_diff);

figure; hold on
times = [0:1/srate:10-1/srate];
contourf(times,freq,Gen_ONvsOFF,40,'linecolor','none');
%plot([0,10],[12,12],'--k')
colormap(ft_colormap('-RdBu'))
set(gca,'FontSize',12)
ylim([1,40])
set(gca,'clim',[-0.5 0.5])
xlabel('Time (s)')
ylabel('Frequency (Hz)')
title('GogglesON vs. GogglesOFF','FontSize',16);

%% compute psd
subs_HCP =  {'HCP001','HCP002','HCP003','HCP004','HCP005'};

subs_GenON = {'004-4006-s2','004-4006-s4',...
    '004-4010-s3','004-4010-s5',....
    'P1001-2','P1002','P1004-2'};
subs_GenOFF = {'004-4006-s2','004-4010-s4','004-4010-s6',...
    'P1001-3','P1002-2'};

subs_BCP =  {'004-4012-s2','BCP003_12Hz','BCP004_12Hz','BCP005_12Hz',...
            'P1009','004-4012-s1','BCP003_6Hz','BCP004_6Hz','BCP005_6Hz'};

subs = subs_GenON;

psd_allSub = [];
psd_base_allSub = [];
psd_stim_allSub = [];
for sub = 1:length(subs)
    subjid = subs{sub};
    disp(subjid)

    filename = fullfile(result_fold,strcat(subjid,'_results2.mat'));
    load(filename);
    
    [x_psd_base,x_psd_stim,x_psd_ratio,freq] = compute_PSD_VEP(results);

    psd_allSub(:,:,sub) = x_psd_ratio;
    psd_base_allSub(:,:,sub) = x_psd_base;
    psd_stim_allSub(:,:,sub) = x_psd_stim;
end

save('GenON_PSD','psd_allSub','psd_base_allSub','psd_stim_allSub','freq');

%%
% get chanlocs
filename = fullfile(data_fold,'0_raw-data','004-4006-s1','CW-04852_wavegard.elc');
chanlocs = readlocs(filename);
eeg_chans = ~strcmpi({chanlocs.labels},'VEOGR');
chanlocs = chanlocs(eeg_chans);

% load PSD
load('HCP_PSD.mat')

%
alpha_ind = freq >= 8 & freq <= 13;


psd_base = mean(sum(psd_base_allSub(:,alpha_ind,3),2),3);
psd_stim = mean(sum(psd_stim_allSub(:,alpha_ind,3),2),3);
psd_diff = mean(sum(psd_allSub(:,alpha_ind,:),2),3);

% 
% figure
% [handle,Zi,grid,Xi,Yi] = topoplot(psd_base,chanlocs,'maplimits','maxmin');
% colormap(ft_colormap('-RdBu'))
% title('Baseline')
% 
% figure
% [handle,Zi,grid,Xi,Yi] = topoplot(psd_stim,chanlocs,'maplimits','maxmin');
% colormap(ft_colormap('-RdBu'))
% title('Stimulation')

figure
[handle,Zi,grid,Xi,Yi] = topoplot(psd_diff,chanlocs,'maplimits','maxmin');
colormap(ft_colormap('-RdBu'))
title('Difference')


%%
% %% Alpha power difference
% subs_HCP =  {'HCP001','HCP002','HCP003','HCP004','HCP005'};
% subs_BCP =  {'P1009','004-4012-s1','004-4012-s2','BCP003_6Hz','BCP003_12Hz',...
%     'BCP004_6Hz','BCP004_12Hz','BCP005_6Hz','BCP005_12Hz'};
% subs_GenON = {'004-4006-s2','004-4006-s4',...
%     '004-4010-s3','004-4010-s5',....
%     'P1001-2','P1002','P1004-2'};
% 
% psd_allSub = [];
% psd_base_allSub = [];
% psd_stim_allSub = [];
% for sub = 1:length(subs_HCP)
%     subjid = subs_HCP{sub};
%     disp(subjid)
% 
%     filename = fullfile(result_fold,strcat(subjid,'_results.mat'));
%     load(filename);
% 
%     [x_psd_base,x_psd_stim,x_psd_ratio,freq] = compute_PSD_VEP(results);
% 
%     psd_allSub(:,sub) = mean(x_psd_ratio,1);
%     psd_base_allSub(:,sub) = mean(x_psd_base,1);
%     psd_stim_allSub(:,sub) = mean(x_psd_stim,1);
% end
% 
% alpha_ind = freq >= 8 & freq <= 10;
% 
% alpha_HCP_base = sum(psd_base_allSub(alpha_ind,1:4),1);
% alpha_HCP_stim = sum(psd_stim_allSub(alpha_ind,1:4),1);
% boxplot([alpha_HCP_base-alpha_HCP_stim]')
% 
% %% Alpha power difference
% subs_HCP =  {'HCP001','HCP002','HCP003','HCP004','HCP005'};
% subs_BCP =  {'P1009','004-4012-s1','004-4012-s2','BCP003_6Hz','BCP003_12Hz',...
%     'BCP004_6Hz','BCP004_12Hz','BCP005_6Hz','BCP005_12Hz'};
% subs_GenON = {'004-4006-s2','004-4006-s4',...
%     '004-4010-s3','004-4010-s5',....
%     'P1001-2','P1002','P1004-2'};
% 
% subs = subs_GenON;
% 
% psd_allSub = [];
% psd_base_allSub = [];
% psd_stim_allSub = [];
% for sub = 1:length(subs)
%     subjid = subs{sub};
%     disp(subjid)
% 
%     filename = fullfile(result_fold,strcat(subjid,'_results2.mat'));
%     load(filename);
% 
%     [x_psd_base,x_psd_stim,x_psd_ratio,freq] = compute_PSD_VEP(results);
% 
%     psd_allSub(:,sub) = mean(x_psd_ratio,1);
%     psd_base_allSub(:,sub) = mean(x_psd_base,1);
%     psd_stim_allSub(:,sub) = mean(x_psd_stim,1);
% end
% 
% alpha_ind = freq >= 8 & freq <= 12;
% 
% alpha_GenON_base = sum(psd_base_allSub(alpha_ind,1:4),1);
% alpha_GenON_stim = sum(psd_stim_allSub(alpha_ind,1:4),1);
% boxplot([alpha_GenON_base-alpha_GenON_stim]')
% [h,p] = ttest([alpha_GenON_base-alpha_GenON_stim],0)

%%
load('HCP005_results.mat')

temp = results.power_time;

figure; hold on
set(gcf,'Position',[680   228   854   650])

subplot(4,1,1); hold on
temp2 = squeeze(mean(temp(occip_ind,15:68,:),1));
contourf(times,results.freq(15:68),temp2,40,'linecolor','none');
plot([0,10],[12,12],'--k')
title('Occipital')


temp2 = squeeze(mean(temp(front_ind,15:68,:),1));
subplot(4,1,2); hold on
contourf(times,results.freq(15:68),temp2,40,'linecolor','none');
plot([0,10],[12,12],'--k')
title('Frontal');


subplot(4,1,3); hold on
temp2 = squeeze(mean(temp(tempo_ind,15:68,:),1));
contourf(times,results.freq(15:68),temp2,40,'linecolor','none');
plot([0,10],[12,12],'--k')
title('Temporal')


subplot(4,1,4); hold on
temp2 = squeeze(mean(temp(parie_ind,15:68,:),1));
contourf(times,results.freq(15:68),temp2,40,'linecolor','none');
plot([0,10],[12,12],'--k')
title('Parietal')






