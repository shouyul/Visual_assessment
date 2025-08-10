%% set plotting parameters
clear all
srate = 250;
times = 0:1/srate:30-1/srate;
times_TEST = 0:1/srate:25-1/srate;


font_size = 20;


%% load data
cd '/Users/shouyuling/Desktop/Postdoc/Projects/Gensight/Data/TumblerTask/Results'

load('AllSub_class_error.mat')

%% get labels
subjids = {'P1001','P1004'};
subjlabels = {'P1001','P1004'};
phases = unique({class_error_all.Phase});
phaselabels = {'Eyes Closed','Eyes Open','Response'};
elecgroups = unique({class_error_all.ElecGroup});%{{'Occipital'},{'Parietal'},{'Frontal'},{'Left_Temporal'},{'Right_Temporal'}};


%% plot individual subjects with single elec group

for s = 1:length(subjids)
    for el = 1:length(elecgroups)
        sub_class_error = class_error_all(...
            strcmp({class_error_all.Subject},subjids{s})...
            &strcmp({class_error_all.ElecGroup},elecgroups{el})...
            &~strcmp({class_error_all.Phase},'Timecourse'));

        data_for_plot = 100 - 100*reshape([sub_class_error.MeanAccuracy],3,2);

        figure;hold on
        x = categorical({'EyesClosed','EyesOpen','Response'});
        b = bar(x,data_for_plot');
        b(1).FaceColor = cus_cmap{1};
        b(2).FaceColor = cus_cmap{5};
        l = plot([0.5,3.5],[50,50],'--','LineWidth',2,'Color',cus_cmap{4});
        ylabel('Classification accuracy (%)');
        ylim([30,90]);
        title(sprintf('Tumbler Decoding - %s',subjlabels{s}));
        set(gca,'FontSize',font_size);
        legend({'Goggles ON','Goggles OFF'})

        filepath = '/Users/shouyuling/Desktop/Postdoc/Projects/Gensight/Data/TumblerTask/Figures/SL/ConfMat';
        filename = sprintf('Tumbler_decoding_%s_%s.png',subjids{s},elecgroups{el});
        print(gcf,[filepath,filename],'-dpng', '-r100', '-painters');
        close()

    end
end

%% compare across elec groups
font_size = 13;

figure;hold on
set(gcf,'Position',[680    91   916   886]);
for s = 1:length(subjids)
    
    sgtitle(sprintf('Tumbler Task - %s',subjlabels{s}),'FontSize',20);
    
    for ph = 1:length(phases)
        sub_class_error = class_error_all(...
            strcmp({class_error_all.Subject},subjids{s})...
            &strcmp({class_error_all.Phase},phases{ph}));

        data_for_plot = 100 - 100*reshape([sub_class_error.MeanAccuracy],2,5)';
        
        subplot(length(phases),1,ph);hold on
        x = categorical({'Occipital','Parietal','Frontal','Left Temporal','Right Temporal'});
        x = reordercats(x,cellstr(x)');
        b = bar(x,data_for_plot');
        b(1).FaceColor = cus_cmap{1};
        b(2).FaceColor = cus_cmap{5};
        l = plot([0.5,5.5],[50,50],'--','LineWidth',2,'Color',cus_cmap{4});
        ylabel('Classification accuracy (%)');
        ylim([25,90]);
        title(phaselabels{ph});
        set(gca,'FontSize',font_size);
        legend({'Goggles ON','Goggles OFF'},'Location','northeastoutside')
    end

    filepath = '/Users/shouyuling/Desktop/Postdoc/Projects/Gensight/Data/TumblerTask/Figures/SL/';
    filename = sprintf('Tumbler_decoding_%s.png',subjids{s});
    print(gcf,[filepath,filename],'-dpng', '-r100', '-painters');
    WaitSecs(0.5)
    clf

end


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
order = categorical({'Object Absent', 'Object Present'});
order = reordercats(order,cellstr(order)');
gogglescon = {'GogglesON','GogglesOFF'};
phase_con = {'Close','Open','Response'};

figure
for s = 2:length(subjids)
    for el = 1:length(elecgroups)
        for g = 1:length(gogglescon)
            for ph = 1:length(phase_con)

                temp_data = class_error_all(strcmp({class_error_all.Subject},subjids{s})...
                            &strcmp({class_error_all.ElecGroup},elecgroups{el})...
                            &strcmp({class_error_all.GogglesCondition},gogglescon{g})...
                            &strcmp({class_error_all.Phase},phase_con{ph}));
                
               
                sgtitle(sprintf('%s - %s - %s - %s',subjids{s},gogglescon{g},phase_con{ph},elecgroups{el}),'FontSize',20);
                confusionchart(round(temp_data.ConfMat),order,'FontSize',15);
                
                filepath = '/Users/shouyuling/Desktop/Postdoc/Projects/Gensight/Data/TumblerTask/Figures/SL/ConfMat/';
                filename = sprintf('Tumbler_conf_%s_%s_%s_%s.png',subjids{s},elecgroups{el},gogglescon{g},phase_con{ph});
                print(gcf,[filepath,filename],'-dpng', '-r100', '-painters');
                WaitSecs(0.5)
                clf
                WaitSecs(0.5)
            end
        end
    end
end

