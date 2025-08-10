%%
clear; clc

main_fold = '/Users/shouyuling/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Postdoc/Projects/Gensight/Data/VEPTask';
data_fold = fullfile(main_fold,'Analysis');
result_fold = fullfile(main_fold,'Results');
fig_fold = fullfile(result_fold,'Figures');


subs = {'004-4006-s1','004-4006-s2','004-4006-s3','004-4006-s4',...
    '004-4010-s3','004-4010-s4','004-4010-s5','004-4010-s6',....
    'P1001-2','P1001-3','P1002','P1004-2'};
goggle_cond = {'GogglesOFF','GogglesON','GogglesOFF','GogglesON',...
    'GogglesON','GogglesOFF','GogglesON','GogglesOFF',...
    'GogglesON','GogglesOFF','GogglesON','GogglesON',...
    };
hasOffData = {'Yes','Yes','Yes','Yes',...
    'Yes','Yes','Yes','Yes',...
    'Yes','Yes','NO','NO'};

hasOffData = strcmp(hasOffData,'Yes');
GogglesCond = strcmp(goggle_cond,'GogglesON');%[0,1,1,0,1,0,1,0,1,1,1,0,1];
% injection = [1,1,1,1,1,1,1,1,1,1,1,0,0];

elecs_incl = 1:128;
elecs_incl(88) = [];

% elecs_incl = [14 15	16	25	30	31	32	37	38	39	40	41	42	65	66	67 ...
%     68	69	77	78	79	80	94	95	96	101	102	103	104	108	109];

%% load results
filename = fullfile(result_fold,'Allsub_results.xlsx');
results_all = readtable(filename);

% overall accuracy gogglesON vs gogglesOFF, featuretype = 'all'
temp_result = results_all(strcmp(results_all.FeatureType,'all'),:);
acc_ON = temp_result{strcmp(temp_result.GogglesCondition,'GogglesON'),'Accuracy'};
acc_OFF = temp_result{strcmp(temp_result.GogglesCondition,'GogglesOFF'),'Accuracy'};

[h,p] = ttest2(acc_ON,acc_OFF,'tail','right')

[h,p] = ttest(acc_ON,25,'tail','right')
[h,p] = ttest(acc_OFF,25,'tail','right')



% stimuli = control, ON vs OFF
acc_ON_ctrl = temp_result{strcmp(temp_result.GogglesCondition,'GogglesON')...
    & strcmp(temp_result.Stimuli,'Control'),'Accuracy'};
acc_OFF_ctrl = temp_result{strcmp(temp_result.GogglesCondition,'GogglesOFF')...
    & strcmp(temp_result.Stimuli,'Control'),'Accuracy'};


[h,p] = ttest2(acc_ON_ctrl,acc_OFF_ctrl,'tail','right')

[h,p] = ttest(acc_ON_ctrl,25,'tail','right')
[h,p] = ttest(acc_OFF_ctrl,25,'tail','right')


% stimuli = disc, ON vs OFF
acc_ON_disc = temp_result{strcmp(temp_result.GogglesCondition,'GogglesON')...
    & strcmp(temp_result.Stimuli,'Disc'),'Accuracy'};
acc_OFF_disc = temp_result{strcmp(temp_result.GogglesCondition,'GogglesOFF')...
    & strcmp(temp_result.Stimuli,'Disc'),'Accuracy'};


[h,p] = ttest2(acc_ON_disc,acc_OFF_disc,'tail','right')

[h,p] = ttest(acc_ON_disc,25,'tail','right')
[h,p] = ttest(acc_OFF_disc,25,'tail','right')


% stimuli = Vertical bar, ON vs OFF
acc_ON_vert = temp_result{strcmp(temp_result.GogglesCondition,'GogglesON')...
    & strcmp(temp_result.Stimuli,'VerticalBar'),'Accuracy'};
acc_OFF_vert = temp_result{strcmp(temp_result.GogglesCondition,'GogglesOFF')...
    & strcmp(temp_result.Stimuli,'VerticalBar'),'Accuracy'};


[h,p] = ttest2(acc_ON_vert,acc_OFF_vert,'tail','right')

[h,p] = ttest(acc_ON_vert,25,'tail','right')
[h,p] = ttest(acc_OFF_vert,25,'tail','right')

% stimuli = Horizontal bar, ON vs OFF
acc_ON_hori = temp_result{strcmp(temp_result.GogglesCondition,'GogglesON')...
    & strcmp(temp_result.Stimuli,'HorizontalBar'),'Accuracy'};
acc_OFF_hori = temp_result{strcmp(temp_result.GogglesCondition,'GogglesOFF')...
    & strcmp(temp_result.Stimuli,'HorizontalBar'),'Accuracy'};


[h,p] = ttest2(acc_ON_hori,acc_OFF_hori,'tail','right')

[h,p] = ttest(acc_ON_hori,25,'tail','right')
[h,p] = ttest(acc_OFF_hori,25,'tail','right')


%% plot spectrogram

dB_allSub = [];
dB_allSub_cond = [];

for sub = 1:length(subs)

    subjid = subs{sub};
    disp(subjid);

    %%
    sub_dB = [];
    sub_dB_cond = [];
    sub_trialtype = [];

    [~,~,pow,sub_trialtype] = convert_to_power(subjid);
    sub_dB = 10*log10(squeeze(mean(pow,4)));

    for s = 1:4 % 1-control, 2-disc, 3-verticalbar, 4-horizontalbar
        sub_dB_cond(:,:,:,s) = 10*log10(mean(pow(:,:,:,sub_trialtype==s),4));
    end

    dB_allSub(:,:,:,sub) = sub_dB;
    dB_allSub_cond(:,:,:,:,sub) = sub_dB_cond;

end
load chirp
sound(y,Fs)
%
filename = fullfile(result_fold,'dB_allSub');
save(filename,'dB_allSub','dB_allSub_cond');


%%
filename = fullfile(result_fold,'dB_allSub');
load(filename)

dB_ON = dB_allSub(:,:,:,GogglesCond==1&hasOffData==1);
dB_OFF = dB_allSub(:,:,:,GogglesCond==0&hasOffData==1);

dB_ON_avg = squeeze(mean(mean(dB_ON(elecs_incl,:,:,:),4)));
dB_OFF_avg = squeeze(mean(mean(dB_OFF(elecs_incl,:,:,:),4)));

dB_diff = dB_ON_avg - dB_OFF_avg;

%% subject 4012
dB_ON_pre = dB_allSub(:,:,:,GogglesCond==1&injection==0);
dB_OFF_pre = dB_allSub(:,:,:,GogglesCond==0&injection==0);

dB_ON_pre_avg = squeeze(mean(mean(dB_ON_pre(elecs_incl,:,:,:),4)));
dB_OFF_pre_avg = squeeze(mean(mean(dB_OFF_pre(elecs_incl,:,:,:),4)));

dB_pre_diff = dB_ON_pre_avg - dB_OFF_pre_avg;

%%

srate = 50;
min_freq = 1;
max_freq = 80;
num_freq = 40;
freq = logspace(log10(min_freq),log10(max_freq),num_freq);
times = [0:1/srate:10-1/srate];

onset = [0:1/6:10-1/6];


%%
clour_lim = 0.5;
fig_visibility = 'on';
fig_size = [168         621        1080         282];

figure('Visible',fig_visibility); hold on
set(gcf,'Position',fig_size)
fig = contourf(times,freq,dB_ON_avg,40,'linecolor','none');
colormap(ft_colormap('-RdBu'))
for i = 1:length(onset)
    plot([onset(i),onset(i)],[1,40],'Color',[0,0,0,0.4],'LineStyle','--')
end
% set(gca,'clim',[-clour_lim clour_lim],'xlim',[0 10],'yscale','log','ytick',logspace(log10(min_freq),log10(max_freq),6),'yticklabel',round(logspace(log10(min_freq),log10(max_freq),6)*10)/10)
set(gca,'clim',[-clour_lim clour_lim],'xlim',[0 10])
set(gca,'ylim',[1,40])
set(gca,'FontSize',20)
xlabel('Time (s)')
ylabel('Frequency (Hz)')
title('GogglesON')

filename = fullfile(fig_fold,'dB_GogglesON.png');
saveas(gca,filename)


figure('Visible',fig_visibility); hold on
set(gcf,'Position',fig_size)
fig = contourf(times,freq,dB_OFF_avg,40,'linecolor','none');
colormap(ft_colormap('-RdBu'))
for i = 1:length(onset)
    plot([onset(i),onset(i)],[1,40],'Color',[0,0,0,0.4],'LineStyle','--')
end
% set(gca,'clim',[-clour_lim clour_lim],'xlim',[0 10],'yscale','log','ytick',logspace(log10(min_freq),log10(max_freq),6),'yticklabel',round(logspace(log10(min_freq),log10(max_freq),6)*10)/10)
set(gca,'clim',[-clour_lim clour_lim],'xlim',[0 10])
set(gca,'ylim',[1,40])
set(gca,'FontSize',20)
xlabel('Time (s)')
ylabel('Frequency (Hz)')
title('GogglesFF')

filename = fullfile(fig_fold,'dB_GogglesOFF.png');
saveas(gca,filename)

figure('Visible',fig_visibility); hold on
set(gcf,'Position',fig_size)
contourf(times,freq,dB_diff,40,'linecolor','none');
colormap(ft_colormap('-RdBu'))
for i = 1:length(onset)
    plot([onset(i),onset(i)],[1,40],'Color',[0,0,0,0.4],'LineStyle','--')
end
% set(gca,'clim',[-clour_lim clour_lim],'xlim',[0 10],'yscale','log','ytick',logspace(log10(min_freq),log10(max_freq),6),'yticklabel',round(logspace(log10(min_freq),log10(max_freq),6)*10)/10)
set(gca,'clim',[-clour_lim clour_lim],'xlim',[0 10])
set(gca,'ylim',[1,40])
set(gca,'FontSize',20)
xlabel('Time (s)')
ylabel('Frequency (Hz)')
title('GogglesON - GogglesOFF')

filename = fullfile(fig_fold,'dB_GogglesDiff.png');
saveas(gca,filename)

%% subject 4012
clour_lim = 1;

figure('Visible','off'); hold on
set(gcf,'Position',[168         621        1080         282])
fig = contourf(times,freq,dB_ON_pre_avg,40,'linecolor','none');
for i = 1:length(onset)
    plot([onset(i),onset(i)],[1,40],'Color',[0,0,0,0.4],'LineStyle','--')
end
% set(gca,'clim',[-clour_lim clour_lim],'xlim',[0 10],'yscale','log','ytick',logspace(log10(min_freq),log10(max_freq),6),'yticklabel',round(logspace(log10(min_freq),log10(max_freq),6)*10)/10)
set(gca,'clim',[-clour_lim clour_lim],'xlim',[0 10])
set(gca,'ylim',[1,40])
set(gca,'FontSize',20)
xlabel('Time (s)')
ylabel('Frequency (Hz)')
title('GogglesON')

filename = fullfile(fig_fold,'dB_pre_GogglesON.png');
saveas(gca,filename)


figure('Visible','off'); hold on
set(gcf,'Position',[168         621        1080         282])
fig = contourf(times,freq,dB_OFF_pre_avg,40,'linecolor','none');
for i = 1:length(onset)
    plot([onset(i),onset(i)],[1,40],'Color',[0,0,0,0.4],'LineStyle','--')
end
% set(gca,'clim',[-clour_lim clour_lim],'xlim',[0 10],'yscale','log','ytick',logspace(log10(min_freq),log10(max_freq),6),'yticklabel',round(logspace(log10(min_freq),log10(max_freq),6)*10)/10)
set(gca,'clim',[-clour_lim clour_lim],'xlim',[0 10])
set(gca,'ylim',[1,40])
set(gca,'FontSize',20)
xlabel('Time (s)')
ylabel('Frequency (Hz)')
title('GogglesFF')

filename = fullfile(fig_fold,'dB_pre_GogglesOFF.png');
saveas(gca,filename)

figure('Visible','off'); hold on
set(gcf,'Position',[168         621        1080         282])
fig = contourf(times,freq,dB_pre_diff,40,'linecolor','none');
for i = 1:length(onset)
    plot([onset(i),onset(i)],[1,40],'Color',[0,0,0,0.4],'LineStyle','--')
end
% set(gca,'clim',[-clour_lim clour_lim],'xlim',[0 10],'yscale','log','ytick',logspace(log10(min_freq),log10(max_freq),6),'yticklabel',round(logspace(log10(min_freq),log10(max_freq),6)*10)/10)
set(gca,'clim',[-clour_lim clour_lim],'xlim',[0 10])
set(gca,'ylim',[1,40])
set(gca,'FontSize',20)
xlabel('Time (s)')
ylabel('Frequency (Hz)')
title('GogglesON - GogglesOFF')

filename = fullfile(fig_fold,'dB_pre_GogglesDiff.png');
saveas(gca,filename)

%%
figure;
set(gcf,'Position',[168         621        1080         282])
contourf(times,freq,dB_ON_avg,40,'linecolor','none');
% set(gca,'clim',[-clour_lim clour_lim],'xlim',[0 10],'yscale','log','ytick',logspace(log10(min_freq),log10(max_freq),6),'yticklabel',round(logspace(log10(min_freq),log10(max_freq),6)*10)/10)
set(gca,'clim',[-clour_lim clour_lim],'xlim',[0 10])
set(gca,'ylim',[1,40])
title('GogglesON')


%%

figure
contourf(times,freq,temp,40,'linecolor','none')
% set(gca,'clim',[-clour_lim clour_lim],'xlim',[0 10],'yscale','log','ytick',logspace(log10(min_freq),log10(max_freq),6),'yticklabel',round(logspace(log10(min_freq),log10(max_freq),6)*10)/10)
set(gca,'clim',[-clour_lim clour_lim],'xlim',[0 10])
set(gca,'ylim',[1.5,40])

%% spectrogram for each condition (all GogglesON)
filename = fullfile(result_fold,'dB_allSub');
load(filename,'dB_allSub_cond')

dB_ON_cond = dB_allSub_cond(elecs_incl,:,:,:,injection==1&GogglesCond==1);
dB_ON_cond = squeeze(mean(mean(dB_ON_cond,5)));

%%
stimuli = {'Control','Disc','VerticalBar','HorizontalBar'};
clour_lim = 1;
for s = 1:size(dB_ON_cond,3)

    figure('Visible','off'); hold on
    set(gcf,'Position',[ 168   573   602   330])
    contourf(times,freq,dB_ON_cond(:,:,s),40,'linecolor','none');
    for i = 1:length(onset)
        plot([onset(i),onset(i)],[1,40],'Color',[0,0,0,0.4],'LineStyle','--')
    end
    % set(gca,'clim',[-clour_lim clour_lim],'xlim',[0 10],'yscale','log','ytick',logspace(log10(min_freq),log10(max_freq),6),'yticklabel',round(logspace(log10(min_freq),log10(max_freq),6)*10)/10)
    set(gca,'clim',[-clour_lim clour_lim],'xlim',[0 10])
    set(gca,'ylim',[1,40])
    set(gca,'FontSize',20)
    xlabel('Time (s)')
    ylabel('Frequency (Hz)')
    title(['GogglesON - ',stimuli{s}])

    filename = fullfile(fig_fold,['dB_ON_all',num2str(s),'.png']);
    saveas(gca,filename)

end

close all

%%
stimuli = {'Control','Disc','VerticalBar','HorizontalBar'};
clour_lim = 1;
for s = 2:size(dB_ON_cond,3)

    figure('Visible','off'); hold on
    set(gcf,'Position',[ 168   573   602   330])
    contourf(times,freq,dB_ON_cond(:,:,s)-dB_ON_cond(:,:,1),40,'linecolor','none');
    for i = 1:length(onset)
        plot([onset(i),onset(i)],[1,40],'Color',[0,0,0,0.4],'LineStyle','--')
    end
    % set(gca,'clim',[-clour_lim clour_lim],'xlim',[0 10],'yscale','log','ytick',logspace(log10(min_freq),log10(max_freq),6),'yticklabel',round(logspace(log10(min_freq),log10(max_freq),6)*10)/10)
    set(gca,'clim',[-clour_lim clour_lim],'xlim',[0 10])
    set(gca,'ylim',[1,40])
    set(gca,'FontSize',20)
    xlabel('Time (s)')
    ylabel('Frequency (Hz)')
    title(['GogglesON - ',stimuli{s}])

    filename = fullfile(fig_fold,['dB_ON_all_vsControl',num2str(s),'.png']);
    saveas(gca,filename)

end

close all

%% spectrogram for each condition (ON vs OFF)
filename = fullfile(result_fold,'dB_allSub');
load(filename,'dB_allSub_cond')

dB_ON_cond = dB_allSub_cond(elecs_incl,:,:,:,GogglesCond==1&hasOffData==1);
dB_ON_cond = squeeze(mean(mean(dB_ON_cond,5)));

dB_OFF_cond = dB_allSub_cond(elecs_incl,:,:,:,GogglesCond==0&hasOffData==1);
dB_OFF_cond = squeeze(mean(mean(dB_OFF_cond,5)));

%
stimuli = {'Control','Disc','VerticalBar','HorizontalBar'};
clour_lim = 1;
for s = 1:size(dB_ON_cond,3)

    figure('Visible','off'); hold on
    set(gcf,'Position',[305   327   480   967])
    subplot(3,1,1); hold on
    contourf(times,freq,dB_ON_cond(:,:,s),40,'linecolor','none');
    for i = 1:length(onset)
        plot([onset(i),onset(i)],[1,40],'Color',[0,0,0,0.4],'LineStyle','--')
    end
    % set(gca,'clim',[-clour_lim clour_lim],'xlim',[0 10],'yscale','log','ytick',logspace(log10(min_freq),log10(max_freq),6),'yticklabel',round(logspace(log10(min_freq),log10(max_freq),6)*10)/10)
    set(gca,'clim',[-clour_lim clour_lim],'xlim',[0 10])
    set(gca,'ylim',[1,40])
    set(gca,'FontSize',16)
    xlabel('Time (s)')
    ylabel('Frequency (Hz)')
    title(['GogglesON'])


    subplot(3,1,2); hold on
    contourf(times,freq,dB_OFF_cond(:,:,s),40,'linecolor','none');
    for i = 1:length(onset)
        plot([onset(i),onset(i)],[1,40],'Color',[0,0,0,0.4],'LineStyle','--')
    end
    % set(gca,'clim',[-clour_lim clour_lim],'xlim',[0 10],'yscale','log','ytick',logspace(log10(min_freq),log10(max_freq),6),'yticklabel',round(logspace(log10(min_freq),log10(max_freq),6)*10)/10)
    set(gca,'clim',[-clour_lim clour_lim],'xlim',[0 10])
    set(gca,'ylim',[1,40])
    set(gca,'FontSize',16)
    xlabel('Time (s)')
    ylabel('Frequency (Hz)')
    title(['GogglesOFF'])

    subplot(3,1,3); hold on
    contourf(times,freq,dB_ON_cond(:,:,s)-dB_OFF_cond(:,:,s),40,'linecolor','none');
    for i = 1:length(onset)
        plot([onset(i),onset(i)],[1,40],'Color',[0,0,0,0.4],'LineStyle','--')
    end
    % set(gca,'clim',[-clour_lim clour_lim],'xlim',[0 10],'yscale','log','ytick',logspace(log10(min_freq),log10(max_freq),6),'yticklabel',round(logspace(log10(min_freq),log10(max_freq),6)*10)/10)
    set(gca,'clim',[-clour_lim clour_lim],'xlim',[0 10])
    set(gca,'ylim',[1,40])
    set(gca,'FontSize',16)
    xlabel('Time (s)')
    ylabel('Frequency (Hz)')
    title(['GogglesON - GogglesOFF'])

    sgtitle(stimuli{s},'FontSize',20)

    filename = fullfile(fig_fold,['dB_cond',num2str(s),'.png']);
    saveas(gca,filename)

end

close all


%% alpha for each condition (ON vs OFF)
filename = fullfile(result_fold,'dB_allSub');
load(filename,'dB_allSub_cond')

%%
freq_ind = 15:28;
time_ind = 300:500;

dB_ON_cond = dB_allSub_cond(elecs_incl,freq_ind,time_ind,:,injection==1&GogglesCond==1);
dB_ON_freq = squeeze(mean(mean(mean(dB_ON_cond))));
data = dB_ON_freq(:);


% dB_ON_cond = dB_allSub_cond(elecs_incl,freq_ind,time_ind,:,injection==1&GogglesCond==1&hasOffData==1);
% dB_ON_freq = squeeze(mean(mean(mean(dB_ON_cond))));
%
% dB_OFF_cond = dB_allSub_cond(elecs_incl,freq_ind,time_ind,:,injection==1&GogglesCond==0&hasOffData==1);
% dB_OFF_freq = squeeze(mean(mean(mean(dB_OFF_cond))));
%
% dB_Diff_freq = dB_ON_freq - dB_OFF_freq;


% data = cat(1,dB_ON_freq(:,1),mean(dB_ON_freq(:,2:3),2),dB_ON_freq(:,4), ...
%     dB_OFF_freq(:,1),mean(dB_OFF_freq(:,2:3),2),dB_OFF_freq(:,4));

%data = cat(1,dB_Diff_freq(:,1),mean(dB_Diff_freq(:,2:3),2),dB_Diff_freq(:,4));

% data = cat(1,dB_ON_freq(:,1),dB_ON_freq(:,2),dB_ON_freq(:,3),dB_ON_freq(:,4), ...
%     dB_OFF_freq(:,1),dB_OFF_freq(:,2),dB_ON_freq(:,3),dB_OFF_freq(:,4));

filename = fullfile(result_fold,'bandPow_results_allON.xlsx');
banPow = readtable(filename);
banPow.Power = data;

writetable(banPow,filename)


%% spectrogram for each individual
% 4006

data = squeeze(mean(dB_allSub(elecs_incl,:,:,2),1));
clour_lim = 1;

figure;
set(gcf,'Position',[168         621        1080         282])
contourf(times,freq,data,40,'linecolor','none');
% set(gca,'clim',[-clour_lim clour_lim],'xlim',[0 10],'yscale','log','ytick',logspace(log10(min_freq),log10(max_freq),6),'yticklabel',round(logspace(log10(min_freq),log10(max_freq),6)*10)/10)
set(gca,'clim',[-clour_lim clour_lim],'xlim',[0 10])
set(gca,'ylim',[1,40])
title('GogglesON')

data = squeeze(mean(dB_allSub(elecs_incl,:,:,1),1));

figure;
set(gcf,'Position',[168         621        1080         282])
contourf(times,freq,data,40,'linecolor','none');
% set(gca,'clim',[-clour_lim clour_lim],'xlim',[0 10],'yscale','log','ytick',logspace(log10(min_freq),log10(max_freq),6),'yticklabel',round(logspace(log10(min_freq),log10(max_freq),6)*10)/10)
set(gca,'clim',[-clour_lim clour_lim],'xlim',[0 10])
set(gca,'ylim',[1,40])
title('GogglesOFF')


% 4010
data = squeeze(mean(mean(dB_allSub(elecs_incl,:,:,[3,5]),1),4));
clour_lim = 0.5;

figure;
set(gcf,'Position',[168         621        1080         282])
contourf(times,freq,data,40,'linecolor','none');
% set(gca,'clim',[-clour_lim clour_lim],'xlim',[0 10],'yscale','log','ytick',logspace(log10(min_freq),log10(max_freq),6),'yticklabel',round(logspace(log10(min_freq),log10(max_freq),6)*10)/10)
set(gca,'clim',[-clour_lim clour_lim],'xlim',[0 10])
set(gca,'ylim',[1,40])
title('GogglesON')

data = squeeze(mean(mean(dB_allSub(elecs_incl,:,:,[4,6]),1),4));

figure;
set(gcf,'Position',[168         621        1080         282])
contourf(times,freq,data,40,'linecolor','none');
% set(gca,'clim',[-clour_lim clour_lim],'xlim',[0 10],'yscale','log','ytick',logspace(log10(min_freq),log10(max_freq),6),'yticklabel',round(logspace(log10(min_freq),log10(max_freq),6)*10)/10)
set(gca,'clim',[-clour_lim clour_lim],'xlim',[0 10])
set(gca,'ylim',[1,40])
title('GogglesOFF')

%% alpha power
filename = fullfile(result_fold,'dB_allSub');
load(filename)

freq_idx = 1:20; % alpha power (2-14Hz)

target_pow = [];

for sub = 1:length(subs)
    subjid = subs{sub};
    disp(subjid);

    sub_dB = dB_allSub(:,:,:,sub);
    sub_alpha = sub_dB(elecs_incl,freq_idx,:);
    target_pow(sub) = mean(sub_alpha(:));

end

sub_acc = [19.6843434350000	28.2769556000000	29.4444444450000	...
    28.8888888900000	30.5555555550000	16.2890869300000	...
    36.6666666650000	25.0000000000000	31.8813131325000	...
    22.2979798000000	25.5555555575000];


scatter(target_pow,sub_acc)


[r,p] = corrcoef(target_pow,sub_acc)

%% prepare trial power for fisher scores

pow_freq_allSub = struct();
pow_chan_allSub = struct();

% set up parameters
srate = 50;
min_freq = 1;
max_freq = 80;
num_freq = 40;
freq = logspace(log10(min_freq),log10(max_freq),num_freq);
freq = freq(1:34);

load(fullfile(result_fold,'chanlocs.mat'));
chan_names = {chanlocs.labels};
chan_names(88) = [];

% % set up column (i.e., feature) names
% cnt = 1;
% for c = 1:length(chan_names)
%     for f = 1:length(freq)
%         feat_names{cnt} = [chan_names{c},'_',num2str(freq(f))];
%         cnt = cnt+1;
%     end
% end

cnt = 1;
for sub = 1:length(subs)

    subjid = subs{sub};
    filename = fullfile(result_fold,strcat(subjid,'_results2.mat'));
    load(filename);
    data = results.trial_power;

    if length(subjid) == 11
        subjid = subjid(5:8);
    else
        subjid = subjid(2:5);
    end

    disp(subjid);

    for trl = 1:size(data,3)
        data_trial = data(:,:,trl);
        chan_avg = mean(data_trial,2);
        freq_avg = mean(data_trial);

        for c = 1:length(chan_names)
            pow_chan_allSub(cnt).Subjid = subjid;
            pow_chan_allSub(cnt).GogglesCondition = goggle_cond{sub};
            pow_chan_allSub(cnt).(chan_names{c}) = chan_avg(c);
        end

        for f = 1:length(freq)
            pow_freq_allSub(cnt).Subjid = subjid;
            pow_freq_allSub(cnt).GogglesCondition = goggle_cond{sub};
            colname = ['Freq_index',num2str(f)];
            pow_freq_allSub(cnt).(colname) = freq_avg(f);
        end

        cnt = cnt+1;
    end

end

filename = fullfile(result_fold,'pow_chan_allSub.xlsx');
writetable(struct2table(pow_chan_allSub),filename);

filename = fullfile(result_fold,'pow_freq_allSub.xlsx');
writetable(struct2table(pow_freq_allSub),filename);



%% Fisher analysis
fisher_rank_allSub = [];
for sub = 1:length(subs)

    subjid = subs{sub};
    filename = fullfile(result_fold,strcat(subjid,'_results2.mat'));
    load(filename);
    data = results.trial_power;
    X = data(elecs_incl,:,:);
    X = reshape(permute(X,[3,1,2]),size(X,3),[]);
    y = results.Info(:,1);


    [orderedInd, orderedPower, fisher] = rankfeat_SL(X, y, 'fisher');
    %     feat_rank = reshape(orderedInd,127,40);
    %     fisher = reshape(orderedPower,127,40);
    fisher_allSub(:,:,sub) = reshape(fisher,127,40);
    %     feat_rank_allSub(:,:,sub) = feat_rank;

end

filename = fullfile(result_fold,'fisher_allSub');
save(filename,'fisher_allSub');
%% all goggles on condition
filename = fullfile(result_fold,'fisher_allSub');
load(filename)

fisher_ON_all = mean(fisher_allSub(:,:,GogglesCond==1),3);
fisher_ON_all = zscore(fisher_ON_all,[],'all');
fisher_ON_freq = mean(fisher_ON_all);
fisher_ON_chan = mean(fisher_ON_all,2);

%% Get chanlocs
filename = fullfile(result_fold,'chanlocs.xlsx');
chanlocs = readtable(filename);

chan_lobe = []; % 1-Occipital, 2-Parietal, 3-Temporal, 4-Frontal
for i = 1:height(chanlocs)
    switch chanlocs{i,'Location'}{1}
        case 'Occipital'
            chan_lobe(i) = 1;
        case 'Parietal'
            chan_lobe(i) = 2;
        case 'Temporal'
            chan_lobe(i) = 3;
        case 'Frontal'
            chan_lobe(i) = 4;
        otherwise
            chan_lobe(i) = 0;
    end

end

chan_lobe(88) = [];

[sorted_chan_lobe,sort_ind] = sort(chan_lobe);
fisher_ON_chan_sorted = fisher_ON_chan(sort_ind);
%% plot fisher
fig_visibility = 'on';
fontsize = 20;

% set up parameters
srate = 250;
min_freq = 1;
max_freq = 80;
num_freq = 40;
freq = logspace(log10(min_freq),log10(max_freq),num_freq);

load(fullfile(result_fold,'chanlocs.mat'));
chan_names = {chanlocs.labels};
chan_names(88) = [];

figure('visible',fig_visibility); hold on
set(gcf,'Position',[ 869   676   912   380])
g = contourf([1:127],freq,fisher_ON_all(sort_ind,:)',40,'linecolor','none');
set(gca,'ylim',[1,40])
plot([38,38],[1,40],'color',[0 0 0 0.5],'LineStyle','--','LineWidth',2);
plot([72,72],[1,40],'color',[0 0 0 0.5],'LineStyle','--','LineWidth',2);
plot([104,104],[1,40],'color',[0 0 0 0.5],'LineStyle','--','LineWidth',2);
%xlabel('Channel')
ylabel('Frequency')
xticklabels(chan_names(sort_ind));
xtickangle(45)
xticks(1:127)
set(gca,'FontSize',fontsize)
ax = gca;
ax.XAxis.FontSize = 8;
%xlabel('Channel','FontSize',20)
t = colorbar;
t.Location='northoutside';
%t.Label = 'Discriminating Power (z-scored)';
t.Box = 'off';

filename = fullfile(fig_fold,'Fisher.png');
saveas(gcf,filename)


figure('visible',fig_visibility); hold on
plot(fisher_ON_freq,freq,'LineWidth',2.5)
set(gca,'ylim',[1,40])
xlabel('Discriminating Power (z-scored)')
ylabel('Frequency')
set(gca,'FontSize',fontsize)


figure('visible',fig_visibility); hold on
set(gcf,'Position',[1553         789         912         315])
plot(1:127,fisher_ON_chan_sorted,'LineWidth',2.5)
plot([38,38],[-10,10],'color',[0 0 0 0.3],'LineStyle','--');
plot([72,72],[-10,10],'color',[0 0 0 0.3],'LineStyle','--');
plot([104,104],[-10,10],'color',[0 0 0 0.3],'LineStyle','--');
set(gca,'xlim',[1,127])
set(gca,'ylim',[-1.5,1.5]);
xticklabels(chan_names(sort_ind));
xtickangle(45)
xticks(1:127)
set(gca,'FontSize',fontsize)
ax = gca;
ax.XAxis.FontSize = 8;
xlabel('Channel','FontSize',fontsize)
ylabel('Discriminating Power (z-scored)','FontSize',16)


%% plot ON vs OFF
ON_color = '#BFD1DE';
OFF_color = '#537DA9';


fisher_ON = mean(fisher_allSub(:,:,GogglesCond==1&hasOffData==1),3);
fisher_ON = zscore(fisher_ON,[],'all');
fisher_ON_freq = mean(fisher_ON);
fisher_ON_chan = mean(fisher_ON,2);


fisher_OFF = mean(fisher_allSub(:,:,GogglesCond==0&hasOffData==1),3);
fisher_OFF = zscore(fisher_OFF,[],'all')
fisher_OFF_freq = mean(fisher_OFF);
feat_rank_OFF_chan = mean(fisher_OFF,2);

feat_rank_ON_chan_sorted = fisher_ON_chan(sort_ind);
feat_rank_OFF_chan_sorted = feat_rank_OFF_chan(sort_ind);



fig_visibility = 'on';

figure('visible',fig_visibility); hold on
set(gcf,'Position',[1553         789         912         315])

contourf([1:127],freq,fisher_ON(sort_ind,:)',40,'linecolor','none');
set(gca,'ylim',[1,40])
title('GogglesON');
plot([38,38],[1,40],'color',[0 0 0 0.5],'LineStyle','--','LineWidth',2);
plot([72,72],[1,40],'color',[0 0 0 0.5],'LineStyle','--','LineWidth',2);
plot([104,104],[1,40],'color',[0 0 0 0.5],'LineStyle','--','LineWidth',2);
ylabel('Frequency')
xticks([])
set(gca,'FontSize',fontsize)
ax = gca;
ax.XAxis.FontSize = 8;
%xlabel('Channel','FontSize',20)
t = colorbar;
t.Location='northoutside';
%t.Label = 'Discriminating Power (z-scored)';
t.Box = 'off';

figure('visible',fig_visibility); hold on
set(gcf,'Position',[1553         789         912         315])

contourf([1:127],freq,fisher_OFF(sort_ind,:)',40,'linecolor','none');
set(gca,'ylim',[1,40])
title('GogglesOFF');
plot([38,38],[1,40],'color',[0 0 0 0.5],'LineStyle','--','LineWidth',2);
plot([72,72],[1,40],'color',[0 0 0 0.5],'LineStyle','--','LineWidth',2);
plot([104,104],[1,40],'color',[0 0 0 0.5],'LineStyle','--','LineWidth',2);
ylabel('Frequency')
xticks([])
set(gca,'FontSize',fontsize)
ax = gca;
ax.XAxis.FontSize = 8;
%xlabel('Channel','FontSize',20)
t = colorbar;
t.Location='northoutside';
%t.Label = 'Discriminating Power (z-scored)';
t.Box = 'off';

figure('visible',fig_visibility); hold on
set(gcf,'Position',[1586         378         404         534])
plot(fisher_ON_freq,freq,'LineWidth',2,'Color',ON_color);
plot(fisher_OFF_freq,freq,'LineWidth',2,'Color',OFF_color);
set(gca,'ylim',[2,40])
% set(gca,'ylim',[-1,2])
xlabel('Discriminating Power (z-scored)')
ylabel('Frequency')
legend({'GogglesON','GogglesOFF'},'Box','on','Location','southeast')
set(gca,'FontSize',16)



figure('visible',fig_visibility); hold on
set(gcf,'Position',[1553         877         912         227])

plot(1:127,feat_rank_ON_chan_sorted,'LineWidth',2,'Color',ON_color);
plot(1:127,feat_rank_OFF_chan_sorted,'LineWidth',2,'Color',OFF_color);
plot([38,38],[-4,4],'color',[0 0 0 0.3],'LineStyle','--');
plot([72,72],[-4,4],'color',[0 0 0 0.3],'LineStyle','--');
plot([104,104],[-4,4],'color',[0 0 0 0.3],'LineStyle','--');
ylabel('Discriminating Power (z-scored)')
set(gca,'xlim',[1,127])
set(gca,'ylim',[-1.5,1.5])
xlabel('Channel')
% legend({'GogglesON','GogglesOFF'})
set(gca,'FontSize',16)
xticklabels(chan_names(sort_ind));
xtickangle(45)
xticks(1:127)
set(gca,'FontSize',fontsize)
ax = gca;
ax.XAxis.FontSize = 8;
xlabel('Channel','FontSize',fontsize)
ylabel('Discriminating Power (z-scored)','FontSize',12)

%% spectrogram for each lobe
filename = fullfile(result_fold,'chan_lobe.mat');
load(filename);

lobe_ord = {'Occipital','Parietal','Temporal','Frontal'}; % 1,2,3,4

filename = fullfile(result_fold,'dB_allSub');
load(filename)

dB_ON = dB_allSub(:,:,:,GogglesCond==1&injection==1&hasOffData==1);
dB_OFF = dB_allSub(:,:,:,GogglesCond==0&injection==1&hasOffData==1);

% dB_ON_avg = squeeze(mean(mean(dB_ON(:,:,:,:),4)));
% dB_OFF_avg = squeeze(mean(mean(dB_OFF(:,:,:,:),4)));
%
% dB_diff = dB_ON_avg - dB_OFF_avg;

%% set up parameter
clour_lim = 0.5;
fig_visibility = 'off';
fig_size = [168   124   535   779];

srate = 50;
min_freq = 1;
max_freq = 80;
num_freq = 40;
freq = logspace(log10(min_freq),log10(max_freq),num_freq);
times = [0:1/srate:10-1/srate];

onset = [0:1/6:10-1/6];
font_size = 16;
for l = 1:4

    % prepare data for plot
    dB_ON_avg = squeeze(mean(mean(dB_ON(chan_lobe==l,:,:,:),4)));
    dB_OFF_avg = squeeze(mean(mean(dB_OFF(chan_lobe==l,:,:,:),4)));

    dB_diff = dB_ON_avg - dB_OFF_avg;


    figure('visible',fig_visibility); hold on
    set(gcf,'Position',fig_size)

    subplot(3,1,1); hold on
    fig = contourf(times,freq,dB_ON_avg,40,'linecolor','none');
    colormap(ft_colormap('-RdBu'))
    for i = 1:length(onset)
        plot([onset(i),onset(i)],[1,40],'Color',[0,0,0,0.4],'LineStyle','--')
    end
    set(gca,'clim',[-clour_lim clour_lim],'xlim',[0 10])
    set(gca,'ylim',[1,40])
    set(gca,'FontSize',font_size)
    xlabel('Time (s)')
    ylabel('Frequency (Hz)')
    title('GogglesON')


    subplot(3,1,2); hold on
    fig = contourf(times,freq,dB_OFF_avg,40,'linecolor','none');
    colormap(ft_colormap('-RdBu'))
    for i = 1:length(onset)
        plot([onset(i),onset(i)],[1,40],'Color',[0,0,0,0.4],'LineStyle','--')
    end
    set(gca,'clim',[-clour_lim clour_lim],'xlim',[0 10])
    set(gca,'ylim',[1,40])
    set(gca,'FontSize',font_size)
    xlabel('Time (s)')
    ylabel('Frequency (Hz)')
    title('GogglesFF')

    subplot(3,1,3); hold on
    contourf(times,freq,dB_diff,40,'linecolor','none');
    colormap(ft_colormap('-RdBu'))
    for i = 1:length(onset)
        plot([onset(i),onset(i)],[1,40],'Color',[0,0,0,0.4],'LineStyle','--')
    end
    set(gca,'clim',[-clour_lim clour_lim],'xlim',[0 10])
    set(gca,'ylim',[1,40])
    set(gca,'FontSize',font_size)
    xlabel('Time (s)')
    ylabel('Frequency (Hz)')
    title('GogglesON - GogglesOFF')


    sgtitle(lobe_ord{l},'FontSize',24);
    filename = fullfile(fig_fold,['dB_',lobe_ord{l},'.png']);
    saveas(gca,filename)

end


%% clustering analysis
% prepare data for python
filename = fullfile(result_fold,'dB_allSub');
load(filename)


dB_ON = dB_allSub(elecs_incl,:,:,GogglesCond==1&injection==1&hasOffData==1);
dB_OFF = dB_allSub(elecs_incl,:,:,GogglesCond==0&injection==1&hasOffData==1);

dB_ON_all = dB_allSub(elecs_incl,:,:,GogglesCond==1&injection==1);

dB_ON_avg = mean(mean(dB_ON,4),3);
dB_OFF_avg = mean(mean(dB_OFF,4),3);
dB_ON_all_avg = mean(mean(dB_ON_all,4),3);

load(fullfile(result_fold,'chanlocs.mat'));
chan_names = {chanlocs.labels};
chan_names(88) = [];

dB_struct_ON = struct();
dB_struct_OFF = struct();
dB_struct_ON_all = struct();

for i = 1:length(chan_names)
    for j = 1:size(dB_ON_avg,2)
        featname = ['dB',num2str(j)];

        dB_struct_ON(i).Channel = chan_names{i};
        eval(['dB_struct_ON(i).',featname,' = dB_ON_avg(i,j);'])

        dB_struct_OFF(i).Channel = chan_names{i};
        eval(['dB_struct_OFF(i).',featname,' = dB_OFF_avg(i,j);'])

        dB_struct_ON_all(i).Channel = chan_names{i};
        eval(['dB_struct_ON_all(i).',featname,' = dB_ON_all_avg(i,j);'])
    end
end

filename = fullfile(result_fold,'dB_forClustering.xlsx');
writetable(struct2table(dB_struct_ON),filename,'FileType','spreadsheet','Sheet','dB_ON');
writetable(struct2table(dB_struct_OFF),filename,'FileType','spreadsheet','Sheet','dB_OFF');
writetable(struct2table(dB_struct_ON_all),filename,'FileType','spreadsheet','Sheet','dB_ON_all');


%% plot cluster
% cluster computed in python

filename = fullfile(result_fold,'Cluster.xlsx');
cluster = readtable(filename);

filename = fullfile(data_fold,'0_raw-data','004-4006-s1','CW-04852_wavegard.elc');
chanlocs = readlocs(filename);
chanlocs(88) = [];

clusters = cluster.Cluster;

% [handle,Zi,grid,Xi,Yi] = topoplot(clusters,chanlocs,'maplimits','maxmin');

X = [chanlocs.X];
Y = [chanlocs.Y];
Z = [chanlocs.Z];

scatter3(X,Y,Z,100,clusters,'filled');

%%
colors = {'#EDAE49','#d1495b','#00798c'};
sz = 200;

figure; hold on
set(gcf,'Position',[1440         630         722         707]);
for i = 1:length(X)
    switch clusters(i)
        case 0
            scatter(Y(i),X(i),sz,'filled','MarkerFaceColor',colors{1},'MarkerEdgeColor',colors{1});

        case 1
            scatter(Y(i),X(i),sz,'filled','MarkerFaceColor',colors{2},'MarkerEdgeColor',colors{2});


        case 2
            scatter(Y(i),X(i),sz,'filled','MarkerFaceColor',colors{3},'MarkerEdgeColor',colors{3});
    end
end
grid off
axis off


%% plot spectrogram for each group
filename = fullfile(result_fold,'dB_allSub');
load(filename)

dB_ON_all = dB_allSub(elecs_incl,:,:,GogglesCond==1&injection==1);
dB_ON_all_avg = mean(dB_ON_all,4);

freq_to_plot = [5:34];
clust1 = squeeze(mean(dB_ON_all_avg(clusters==0,freq_to_plot,:)));
clust2 = squeeze(mean(dB_ON_all_avg(clusters==1,freq_to_plot,:)));
clust3 = squeeze(mean(dB_ON_all_avg(clusters==2,freq_to_plot,:)));

%%

clour_lim = 0.5;
fig_visibility = 'on';
fig_size = [ 1440         402         741         935];


srate = 50;
min_freq = 1;
max_freq = 80;
num_freq = 40;
freq = logspace(log10(min_freq),log10(max_freq),num_freq);
times = [0:1/srate:10-1/srate];

freq = freq(freq_to_plot);

onset = [0:1/6:10-1/6];
font_size = 16;

figure('visible',fig_visibility); hold on
set(gcf,'Position',fig_size);

subplot(3,1,1); hold on
fig = contourf(times,freq,clust1,40,'linecolor','none');
colormap(ft_colormap('-RdBu'))
for i = 1:length(onset)
    plot([onset(i),onset(i)],[1,40],'Color',[0,0,0,0.4],'LineStyle','--')
end
set(gca,'clim',[-clour_lim clour_lim],'xlim',[0 10])
set(gca,'ylim',[2,40])
set(gca,'FontSize',font_size)
xlabel('Time (s)')
ylabel('Frequency (Hz)')
title('Cluster 1')


subplot(3,1,2); hold on
fig = contourf(times,freq,clust2,40,'linecolor','none');
colormap(ft_colormap('-RdBu'))
for i = 1:length(onset)
    plot([onset(i),onset(i)],[1,40],'Color',[0,0,0,0.4],'LineStyle','--')
end
set(gca,'clim',[-clour_lim clour_lim],'xlim',[0 10])
set(gca,'ylim',[2,40])
set(gca,'FontSize',font_size)
xlabel('Time (s)')
ylabel('Frequency (Hz)')
title('Cluster 2')

subplot(3,1,3); hold on
fig = contourf(times,freq,clust3,40,'linecolor','none');
colormap(ft_colormap('-RdBu'))
for i = 1:length(onset)
    plot([onset(i),onset(i)],[1,40],'Color',[0,0,0,0.4],'LineStyle','--')
end
set(gca,'clim',[-clour_lim clour_lim],'xlim',[0 10])
set(gca,'ylim',[2,40])
set(gca,'FontSize',font_size)
xlabel('Time (s)')
ylabel('Frequency (Hz)')
title('Cluster 3')

