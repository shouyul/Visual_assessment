%%
%% set parameter

clear; clc

main_fold = '/Users/shouyuling/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Postdoc/Projects/Gensight/Data/VEPTask';
data_fold = fullfile(main_fold,'Analysis');
fig_fold = fullfile(main_fold,'Figures');
result_fold = fullfile(main_fold,'Results');

% set subject info
subs = {'HCP001','HCP002','HCP003'};
goggle_cond = {'GogglesOFF','GogglesOFF','GogglesOFF'};

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

for sub = 1:length(subs)
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

     %% convert data to power
    if ~isfield(results,'trial_power') || recompute_power
        % disp('Wavelet transformation')
        % [results.trial_power,results.power_time,results.dB_time, results.freq] = convert_to_power_VEP(results.trials);

        convert_to_power_SSVEP(results.trials);

        disp('Running DFT')
        



    end












end