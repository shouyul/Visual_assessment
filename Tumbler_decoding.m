%%
clear all; clc
main_fold = '/Users/shouyuling/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Postdoc/Projects/Gensight/Data/TumblerTask';
data_fold = fullfile(main_fold,'Analysis');
fig_fold = fullfile(main_fold,'Figures');
result_fold = fullfile(main_fold,'Results');

subs = {'004-4006-s1','004-4006-s2','004-4010-s2','P1001-4','P1002-2','P1002-3','P1004-2','P1009','004-4012'};
% Pitt patients, block1-3 GogglesON, block4-6 GogglesOFF
% P1002-3, block1-4 GogglesOFF, block5-8, GogglesON
% Remaining Paris patients, block1-4 GogglesON, block5-8, GogglesOFF
subs = {'004-4012','P1009'};

srate = 250;
crossval_method = 'leaveout';
kfold = 5;

%% load data

for sub = 1:length(subs)
    subjid = subs{sub};
    disp(subjid);

    filename = fullfile(result_fold,strcat(subjid,'_results.mat'));
    load(filename);

    freq_idx_close = results.freq.freq_close < 1000;% & results.freq.freq_close > 0.99;
    freq_idx_open = results.freq.freq_close < 1000;% & results.freq.freq_close > 0.99;

    % GogglesON
    % EyesClosed
    X = results.trial_power.GogglesON_EyesClosed(:,freq_idx_close,:);
    y = results.info.GogglesON;

    acc = decoding_Tumbler_v3(X,y,crossval_method,kfold);

    results.decoding.GogglesON_EyesClosed = acc; clear acc

    % GogglesON
    % EyesOpen
    X = results.trial_power.GogglesON_EyesOpen(:,freq_idx_open,:);
    y = results.info.GogglesON;

    acc = decoding_Tumbler_v3(X,y,crossval_method,kfold);

    results.decoding.GogglesON_EyesOpen = acc; clear acc


    % GogglesOFF
    % EyesClosed
    X = results.trial_power.GogglesOFF_EyesClosed(:,freq_idx_close,:);
    y = results.info.GogglesOFF;

    acc = decoding_Tumbler_v3(X,y,crossval_method,kfold);
    
    results.decoding.GogglesOFF_EyesClosed = acc; clear acc

    % GogglesONFF
    % EyesOpen
    X = results.trial_power.GogglesOFF_EyesOpen(:,freq_idx_open,:);
    y = results.info.GogglesOFF;

    acc = decoding_Tumbler_v3(X,y,crossval_method,kfold);
    
    results.decoding.GogglesOFF_EyesOpen = acc; clear acc

    results.decoding

    filename = fullfile(result_fold,strcat(subjid,'_results.mat'));
    save(filename,'results')


end

%% aggregate data for plotting
Allsub_results = struct();

cnt = 1;
for sub = 1:length(subs)
    subjid = subs{sub};
    disp(subjid);

    if subjid(1)=='0'
        sub_ID = subjid(5:8);
        sess_ID = str2num(subjid(end));
    else
        sub_ID = subjid(2:5);
        if length(subjid) == 5
            sess_ID = 1;
        else
            sess_ID = str2num(subjid(end));
        end
    end

    filename = fullfile(result_fold,strcat(subjid,'_results.mat'));
    load(filename);

    GogglesCondition = {'GogglesON','GogglesOFF','GogglesON','GogglesOFF'};
    EyesCondition = {'EyesOpen','EyesOpen','EyesClosed','EyesClosed'};
    for i = 1:4
        Allsub_results(cnt).Subjid = sub_ID;
        Allsub_results(cnt).Task = 'Tumbler';
        Allsub_results(cnt).Session = sess_ID;
        Allsub_results(cnt).ElecGroup = 'All_electrodes';
        Allsub_results(cnt).FreqGroup = '1to40';

        Allsub_results(cnt).GogglesCondition = GogglesCondition{i};
        Allsub_results(cnt).EyesCondition = EyesCondition{i};
        fieldname = [GogglesCondition{i},'_',EyesCondition{i}];
        try
            eval(['temp_acc = results.decoding.',fieldname,';']);
        catch
            temp_acc = nan;
        end
        Allsub_results(cnt).Accuracy = temp_acc*100;
        cnt = cnt+1;
    end
end

filename = fullfile(result_fold,'Allsub_results_withBCP.xlsx');
writetable(struct2table(Allsub_results),filename)
