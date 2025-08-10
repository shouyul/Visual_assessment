%% set plotting parameters
clear all
srate = 250;
times = 0:1/srate:30-1/srate;

cus_cmap = {'#80AFBF','#608595','#DFC286','#C07A92','#E2C3C9'};

font_size = 20;

%%
cd '/Users/shouyuling/Desktop/Postdoc/Projects/Gensight/Data/VEPTask/Results'
load('AllSub_class_error_BarRing.mat');

% 
% P1001 = load('P1001_class_error.mat');
% P1004 = load('P1004_class_error.mat');
% TEST = load('TEST_class_error.mat');
% 
% class_errors = 100*[mean(P1001.class_error_all),mean(P1004.class_error_all),mean(TEST.class_error_all)];
% class_errors_1sec = 100*[mean(P1001.class_error_1sec_avg),mean(P1004.class_error_1sec_avg),mean(TEST.class_error_1sec_avg)];
% class_errors_time = 100*[P1001.class_error_time;P1004.class_error_time;TEST.class_error];

%% get labels
subjids = {'P1001','P1004'};
subjlabels = {'P1001','P1004'};
elecgroups = unique({class_error_all.ElecGroup});
elecgroups = elecgroups([3,4,1,2,5]);
segcon = {'Entire trial','1-sec segments','3-sec segments'};

%% plot all subjects for each elec group

figure;hold on
set(gcf,'Position',[302   473   923   481])
for el = 1:length(elecgroups)
    elcgroup_class_error = class_error_all(strcmp({class_error_all.ElecGroup},elecgroups{el}));
    
    data_for_plot = [];
    for s = 1:length(subjids)
        data_for_plot = [data_for_plot;...
            [elcgroup_class_error(strcmp({elcgroup_class_error.Subject},subjids{s})).MeanAccuracy]];
    end
    data_for_plot = 100 - 100*data_for_plot;
    
    
    
    x = categorical({'P1001','P1004'});
    x = reordercats(x,cellstr(x)');
    b = bar(x,data_for_plot); hold on
    b(1).FaceColor = cus_cmap{1};
    b(2).FaceColor = cus_cmap{5};
    b(3).FaceColor = cus_cmap{3};
    l = plot([0.5,2.5],[1/2,1/2]*100,'--','LineWidth',2,'Color',cus_cmap{4});
    ylabel('Classification accuracy (%)');
    ylim([35,70]);
    title('VEP Decoding');
    legend(segcon,'Location','bestoutside');
    set(gca,'FontSize',font_size);

    filepath = '/Users/shouyuling/Desktop/Postdoc/Projects/Gensight/Data/VEPTask/Figures/SL/';
    filename = sprintf('VEP_decoding_BarRing_%s.png',elecgroups{el});
    print(gcf,[filepath,filename],'-dpng', '-r300', '-painters');
    clf

end

%% compare elec groups
subjids = {'P1001','P1004'};
subjlabels = {'P1001','P1004'};
font_size = 13;
figure;hold on
set(gcf,'Position',[680    91   916   886]);
sgtitle(sprintf('VEP Task'),'FontSize',20);

for s = 1:length(subjids)
    sub_class_error = class_error_all(strcmp({class_error_all.Subject},subjids{s}));
        
    data_for_plot = [];
    for el = 1:length(elecgroups)
        data_for_plot = [data_for_plot;...
            [sub_class_error(strcmp({sub_class_error.ElecGroup},elecgroups{el})).MeanAccuracy]];
    end
    data_for_plot = 100 - 100*data_for_plot;

    subplot(length(subjids),1,s);hold on
    x = categorical({'Occipital','Parietal','Frontal','Left Temporal','Right Temporal'});
    x = reordercats(x,cellstr(x)');
    b = bar(x,data_for_plot);
    b(1).FaceColor = cus_cmap{1};
    b(2).FaceColor = cus_cmap{5};
    b(3).FaceColor = cus_cmap{3};
    l = plot([0.5,5.5],[1/2,1/2]*100,'--','LineWidth',2,'Color',cus_cmap{4});
    ylabel('Classification accuracy (%)');
    ylim([35,70]);
    title(subjlabels{s});
    legend(segcon,'Location','bestoutside');
    set(gca,'FontSize',font_size);
end
WaitSecs(0.5);
filepath = '/Users/shouyuling/Desktop/Postdoc/Projects/Gensight/Data/VEPTask/Figures/SL/';
filename = sprintf('VEP_decoding_BarRing.png');
print(gcf,[filepath,filename],'-dpng', '-r300', '-painters');

%% compute confusion matrix
num_perm = size(class_error_all(1).Predictions,2);

for i = 1:length(class_error_all)
    conf_mat = [];
    if ~isempty(class_error_all(i).Predictions)
        for n = 1:num_perm
            conf_mat(:,:,n) = confusionmat(class_error_all(i).TrueLabels,class_error_all(i).Predictions(:,n));
        end
    end
    class_error_all(i).ConfMat = mean(conf_mat,3);
end

%% plot confusion matrix
order = categorical({'Bars','Rings'});
order = reordercats(order,cellstr(order)');
seg_con = {'EntireTrial','OneSecAverage','ThreeSecAverage'};

figure
for s = 1:length(subjids)
    for el = 1:length(elecgroups)
        for sg = 1:length(seg_con)

            temp_data = class_error_all(strcmp({class_error_all.Subject},subjids{s})...
                                        &strcmp({class_error_all.ElecGroup},elecgroups{el})...
                                        &strcmp({class_error_all.TrialSegCondition},seg_con{sg}));

            sgtitle(sprintf('%s - %s - %s',subjids{s},elecgroups{el},seg_con{sg}));
            confusionchart(round(temp_data.ConfMat),order,'FontSize',15);

            filepath = '/Users/shouyuling/Desktop/Postdoc/Projects/Gensight/Data/VEPTask/Figures/SL/ConfMat/';
            filename = sprintf('VEP_conf_BarRing_%s_%s_%s',subjids{s},elecgroups{el},seg_con{sg});
            print(gcf,[filepath,filename],'-dpng', '-r100', '-painters');
            WaitSecs(0.5)
            clf
            WaitSecs(0.5)
        end
    end
end




