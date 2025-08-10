%%
clear; clc

main_fold = '/Users/shouyuling/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Postdoc/Projects/Gensight/Data';

task_name = 'Tumbler';
subjid = '004-4006-s2';

data_fold = fullfile(main_fold,[task_name,'Task'],'Analysis','0_raw-data',subjid);


%% load behav data
filename = fullfile(data_fold,'TumblerBehavData.xlsx');

behav_data = readtable(filename);
behav_data = table2array(behav_data(:,2:end)); % cols: Block, Trial, Goggles, TrialType, Answer, Confidence

% goggles on
hit_ON = behav_data((behav_data(:,3)==1 & behav_data(:,4)==1 & behav_data(:,5)==1),:);
miss_ON = behav_data((behav_data(:,3)==1 & behav_data(:,4)==1 & behav_data(:,5)==0),:);

fa_ON = behav_data((behav_data(:,3)==1 & behav_data(:,4)==0 & behav_data(:,5)==1),:);
rej_ON = behav_data((behav_data(:,3)==1 & behav_data(:,4)==0 & behav_data(:,5)==0),:);

% goggles off
hit_OFF = behav_data((behav_data(:,3)==0 & behav_data(:,4)==1 & behav_data(:,5)==1),:);
miss_OFF = behav_data((behav_data(:,3)==0 & behav_data(:,4)==1 & behav_data(:,5)==0),:);

fa_OFF = behav_data((behav_data(:,3)==0 & behav_data(:,4)==0 & behav_data(:,5)==1),:);
rej_OFF = behav_data((behav_data(:,3)==0 & behav_data(:,4)==0 & behav_data(:,5)==0),:);

% goggles on overall accuracy
data_ON = behav_data(behav_data(:,3)==1,:);
acc_ON = mean(data_ON(:,4)==data_ON(:,5));

data_OFF = behav_data(behav_data(:,3)==0,:);
acc_OFF = mean(data_OFF(:,4)==data_OFF(:,5));

% confidence level
conf_ON = mean(data_ON(:,6));
conf_OFF = mean(data_OFF(:,6));

% accuracy for each block
for b = 1:6
    temp_data = behav_data(behav_data(:,1)==b,:);
    acc_blk(b) = mean(temp_data(:,4)==temp_data(:,5));
    conf_blk(b) = mean(temp_data(:,6));
end

%%
figure
set(gcf,'Position',[414         360        1054         420]);
subplot(1,2,1); hold on
bar(acc_blk*100)
plot([1,3],[acc_ON,acc_ON]*100,'LineWidth',5,'Color',[1,0,0,0.3])
plot([4,6],[acc_OFF,acc_OFF]*100,'LineWidth',5,'Color',[1,0,0,0.3])
xlabel('Block')
ylabel('Accuracy (%)')
ylim([0,80]);
xticklabels({'1-GogglesON','2-GogglesON','3-GogglesON','4-GogglesOFF','5-GogglesOFF','6-GogglesOFF'});
title('Behavioural Performance')

subplot(1,2,2)
plot(conf_blk); hold on
plot([1,3],[conf_ON,conf_ON],'LineWidth',5,'Color',[1,0,0,0.3]);
plot([4,6],[conf_OFF,conf_OFF],'LineWidth',5,'Color',[1,0,0,0.3]);
ylim([2,4])
xlabel('Block')
ylabel('Confidence Rating')
xticklabels({'1-GogglesON','2-GogglesON','3-GogglesON','4-GogglesOFF','5-GogglesOFF','6-GogglesOFF'});
title('Confidence Level')

%%
figure
set(gcf,'Position',[680         569        1049         408])
subplot(1,3,1)
histogram(behav_data(:,6))
xlim([0,6]);
xlabel('Confidence Rating');
ylabel('Count')
title('All conditions')

subplot(1,3,2)
histogram(data_ON(:,6))
xlim([0,6]);
xlabel('Confidence Rating');
ylabel('Count')
ylim([0,20])
title('GogglesON')

subplot(1,3,3)
histogram(data_OFF(:,6))
xlim([0,6]);
xlabel('Confidence Rating');
ylabel('Count')
ylim([0,20])
title('GogglesOFF')

sgtitle('Distribution of Confidence Ratings','FontSize',24)