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

%%
% get chanlocs
filename = fullfile(data_fold,'0_raw-data','004-4006-s1','CW-04852_wavegard.elc');
chanlocs = readlocs(filename);
eeg_chans = ~strcmpi({chanlocs.labels},'VEOGR');
chanlocs = chanlocs(eeg_chans);

%%
filename = ['/Users/shouyuling/Library/CloudStorage/OneDrive-Universityof' ...
    'Pittsburgh/Postdoc/Projects/Gensight/Data/TumblerTask/Analysis/0_raw-data/P1001-4/CA-208.elc'];
chanlocs = readeetraklocs(filename);
%%
% load PSD
load('BCP_PSD.mat')

%
alpha_ind = freq >= 8 & freq <= 13;


psd_base = mean(sum(psd_base_allSub(:,alpha_ind,3),2),3);
psd_stim = mean(sum(psd_stim_allSub(:,alpha_ind,3),2),3);


psd_allSub = log10(psd_allSub);
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

%
figure
[handle,Zi,grid,Xi,Yi] = topoplot(psd_diff,chanlocs,'maplimits',[-2,2.5]);
colormap(ft_colormap('-RdBu'))
title('Blind Control', FontSize=24)
% hcb=colorbar;
% hcb.Title.String = 'dB';

%%
% get chanlocs
filename = fullfile(data_fold,'0_raw-data','004-4006-s2','CW-04852_wavegard.elc');
chanlocs = readlocs(filename,'importmode','eeglab');
eeg_chans = ~strcmpi({chanlocs.labels},'VEOGR');
chanlocs = chanlocs(eeg_chans);

% load PSD
load('HCP_PSD.mat')

%
alpha_ind = freq >= 8 & freq <= 13;


psd_base = mean(sum(psd_base_allSub(:,alpha_ind,3),2),3);
psd_stim = mean(sum(psd_stim_allSub(:,alpha_ind,3),2),3);


psd_allSub = log10(psd_allSub);
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

%
figure
[handle,Zi,grid,Xi,Yi] = topoplot(psd_diff,chanlocs,'maplimits',[-1,2]);
% colormap(ft_colormap('-RdBu'))
title('Sighted Control', FontSize=24)
hcb=colorbar;
hcb.Title.String = 'dB';

%%

% load PSD
load('GenON_PSD.mat')
subjects = {'Pt 4006','Pt 4010','Pt 1001','Pt 1002','Pt 1004'};
%
alpha_ind = freq >= 8 & freq <= 13;


psd_base = mean(sum(psd_base_allSub(:,alpha_ind,3),2),3);
psd_stim = mean(sum(psd_stim_allSub(:,alpha_ind,3),2),3);

psd_allSub = cat(3,mean(psd_allSub(:,:,1:2),3),...
    mean(psd_allSub(:,:,3:4),3),psd_allSub(:,:,5:end));
psd_allSub = log10(psd_allSub);
psd_diff = mean(sum(psd_allSub(:,alpha_ind,1:4),2),3);

%%
figure
[handle,Zi,grid,Xi,Yi] = topoplot(psd_diff,chanlocs,'maplimits',[-1,2]);
% colormap(ft_colormap('-RdBu'))
title('GenSight Patient', FontSize=24)
% hcb=colorbar;
% hcb.Title.String = 'dB';


%%
figure
for i = 1:size(psd_allSub,3)
    psd_diff = mean(sum(psd_allSub(:,alpha_ind,i),2),3);

subplot(2,3,i)

[handle,Zi,grid,Xi,Yi] = topoplot(psd_diff,chanlocs,'maplimits','maxmin');
colormap(ft_colormap('-RdBu'))
title(subjects{i}, FontSize=24)
hcb=colorbar;
hcb.Title.String = 'dB';
end

%% calculate alpha power
alpha_ind = freq >= 8 & freq <= 13;
left_occi_labels = {'Z15L','Z16Z','L16L','Z16L','L17Z','L17L','Z17L','L18Z'};
right_occi_labels = {'Z15R','R16Z','Z16R','R16R','R17Z','Z17R','R17R','R18Z'};

occi_ind = cellfun(@(c)ismember({chanlocs.labels},c),{left_occi_labels,right_occi_labels},'UniformOutput',false);


load('HCP_PSD.mat')
% psd_allSub = 100*log10(psd_stim_allSub./psd_base_allSub);
psd_allSub = log10(psd_allSub);
occi_alpha_left = squeeze(mean(mean(psd_allSub(occi_ind{1},alpha_ind,:),2)));
occi_alpha_right = squeeze(mean(mean(psd_allSub(occi_ind{2},alpha_ind,:),2)));

occi_alpha_HCP = cat(2,occi_alpha_left,occi_alpha_right);
occi_alpha_HCP = mean(occi_alpha_HCP,2);


load('BCP_PSD.mat')
%psd_allSub = 100*log10(psd_stim_allSub./psd_base_allSub);
psd_allSub = log10(psd_allSub);
occi_alpha_left = squeeze(mean(mean(psd_allSub(occi_ind{1},alpha_ind,:),2)));
occi_alpha_right = squeeze(mean(mean(psd_allSub(occi_ind{2},alpha_ind,:),2)));

occi_alpha_BCP = cat(2,occi_alpha_left,occi_alpha_right);
occi_alpha_BCP = mean(occi_alpha_BCP,2);


load('GenON_PSD.mat')
%psd_allSub = 100*log10(psd_stim_allSub./psd_base_allSub);
occi_alpha_left = squeeze(mean(mean(psd_allSub(occi_ind{1},alpha_ind,:),2)));
occi_alpha_right = squeeze(mean(mean(psd_allSub(occi_ind{2},alpha_ind,:),2)));

occi_alpha_GenON = cat(2,occi_alpha_left,occi_alpha_right);
occi_alpha_GenON = mean(occi_alpha_GenON,2);


% occi_alpha_allSub = padcat(occi_alpha_GenON(:,1),occi_alpha_GenON(:,2),...
%     occi_alpha_BCP(:,1),occi_alpha_BCP(:,2),occi_alpha_HCP(:,1),occi_alpha_HCP(:,2));


occi_alpha_allSub = padcat(occi_alpha_GenON,occi_alpha_BCP,occi_alpha_HCP);

bar(nanmean(occi_alpha_allSub))

%%
alpha_ind = freq >= 8 & freq <= 13;
occi_labels = {'Z15L','Z16Z','L16L','Z16L','L17Z','L17L','Z17L','L18Z',...
    'Z15R','R16Z','Z16R','R16R','R17Z','Z17R','R17R','R18Z'};
occi_ind = ismember({chanlocs.labels},occi_labels);

load('HCP_PSD.mat')
psd_allSub = zscore(psd_allSub,[],'all');
% psd_allSub = 10*log10(psd_stim_allSub./psd_base_allSub);
% psd_allSub = psd_stim_allSub;
occi_alpha_HCP = squeeze(mean(mean(psd_allSub(occi_ind,alpha_ind,:),2),1));

load('BCP_PSD.mat')
psd_allSub = zscore(psd_allSub,[],'all');
% psd_allSub = 10*log10(psd_stim_allSub./psd_base_allSub);
% psd_allSub = psd_stim_allSub;
occi_alpha_BCP = squeeze(mean(mean(psd_allSub(occi_ind,alpha_ind,:),2),1));

load('GenON_PSD.mat')
psd_allSub = zscore(psd_allSub,[],'all');
% psd_allSub = 10*log10(psd_stim_allSub./psd_base_allSub);
% psd_allSub = psd_stim_allSub;
occi_alpha_GenON = squeeze(mean(mean(psd_allSub(occi_ind,alpha_ind,:),2)));

occi_alpha_allSub = padcat(occi_alpha_GenON,occi_alpha_BCP,occi_alpha_HCP);
% bar(nanmean(occi_alpha_allSub))
boxplot(occi_alpha_allSub)

%% write data to .csv for plot in python
temp = struct();

subs_HCP =  {'HCP001','HCP002','HCP003','HCP004','HCP005'};
subs_BCP = {'BCP001','BCP002','BCP002','BCP003','BCP003','BCP004','BCP004','BCP005','BCP005'};
subs_GenON = {'P4006','P4006','P4010','P4010','P1001','P1002','P1004'};

n = 1;
for s = 1:length(subs_HCP)
    temp(n).Subjid = subs_HCP{s};
    temp(n).SubjectType = 'Sighted Control';
    temp(n).occi_alpha = occi_alpha_HCP(s);
    n = n+1;
end

for s = 1:length(subs_BCP)
    temp(n).Subjid = subs_BCP{s};
    temp(n).SubjectType = 'Blind Control';
    temp(n).occi_alpha = occi_alpha_BCP(s);
    n = n+1;
end

for s = 1:length(subs_GenON)
    temp(n).Subjid = subs_GenON{s};
    temp(n).SubjectType = 'GenSight Patient';
    temp(n).occi_alpha = occi_alpha_GenON(s);
    n = n+1;
end

filename = fullfile(result_fold,'Allsub_results_OcciAlpha.csv');
writetable(struct2table(temp),filename)


%% alpha power over frontal electrodes
alpha_ind = freq >= 7 & freq <= 12;
front_labels = {'L3Z','Z3Z','R3Z','L5Z','Z5Z','R5Z','L6Z',...
    'Z6Z','R6Z','L7Z','Z7Z','R7Z'};
front_ind = ismember({chanlocs.labels},front_labels);

load('HCP_PSD.mat');
psd_allSub = zscore(psd_allSub,[],'all');
% psd_allSub = 10*log10(psd_stim_allSub./psd_base_allSub);
% psd_allSub = psd_stim_allSub;
front_alpha_HCP = squeeze(mean(mean(psd_allSub(front_ind,alpha_ind,:),2),1));

load('BCP_PSD.mat');
psd_allSub = zscore(psd_allSub,[],'all');
% psd_allSub = 10*log10(psd_stim_allSub./psd_base_allSub);
% psd_allSub = psd_stim_allSub;
front_alpha_BCP = squeeze(mean(mean(psd_allSub(front_ind,alpha_ind,:),2),1));

load('GenON_PSD.mat')
psd_allSub = zscore(psd_allSub,[],'all');
% psd_allSub = 10*log10(psd_stim_allSub./psd_base_allSub);
% psd_allSub = psd_stim_allSub;
front_alpha_GenON = squeeze(mean(mean(psd_allSub(front_ind,alpha_ind,:),2)));

front_alpha_allSub = padcat(front_alpha_GenON,front_alpha_BCP,front_alpha_HCP);
% bar(nanmean(front_alpha_allSub))
boxplot(front_alpha_allSub)


%% write data to .csv for plot in python
temp = struct();

subs_HCP =  {'HCP001','HCP002','HCP003','HCP004','HCP005'};
subs_BCP = {'BCP001','BCP002','BCP002','BCP003','BCP003','BCP004','BCP004','BCP005','BCP005'};
subs_GenON = {'P4006','P4006','P4010','P4010','P1001','P1002','P1004'};

n = 1;
for s = 1:length(subs_HCP)
    temp(n).Subjid = subs_HCP{s};
    temp(n).SubjectType = 'Sighted Control';
    temp(n).front_alpha = front_alpha_HCP(s);
    n = n+1;
end

for s = 1:length(subs_BCP)
    temp(n).Subjid = subs_BCP{s};
    temp(n).SubjectType = 'Blind Control';
    temp(n).front_alpha = front_alpha_BCP(s);
    n = n+1;
end

for s = 1:length(subs_GenON)
    temp(n).Subjid = subs_GenON{s};
    temp(n).SubjectType = 'GenSight Patient';
    temp(n).front_alpha = front_alpha_GenON(s);
    n = n+1;
end

filename = fullfile(result_fold,'Allsub_results_FrontAlpha.csv');
writetable(struct2table(temp),filename)