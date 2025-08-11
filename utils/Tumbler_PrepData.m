%%
clear all; clc
main_fold = '/Users/shouyuling/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Postdoc/Projects/Gensight/Data/TumblerTask';
data_fold = fullfile(main_fold,'Analysis');
fig_fold = fullfile(main_fold,'Figures');
result_fold = fullfile(main_fold,'Results');

%% get trialsequence from log
subjid = '004-4012-s1';
filename = fullfile(data_fold,'0_raw-data',subjid,[subjid,'_log.txt']);
fileID = fopen(filename,'r');


trialseq = [];

% before block 0 started
tline = fgets(fileID);
while ischar(tline)
    tline_cell = split(strip(tline),' ');
    if strcmpi(tline_cell(end),'started.')
        break
    end
    tline = fgets(fileID);
end

while ischar(tline)
    tline_cell = split(strip(tline),' ');

    if contains(tline,'Trial sequence')
        temp_seq = tline(17:35);


        temp_seq_num = [];
        for i = 1:length(temp_seq)
            temp_seq_num = cat(2,temp_seq_num,str2num(temp_seq(i)));
        end

        trialseq = cat(1,trialseq,temp_seq_num);
    end
    tline = fgets(fileID);
end

fclose(fileID);

filename = fullfile(data_fold,'0_raw-data',subjid,[subjid,'_TrialSequence.txt']);
dlmwrite(filename,trialseq);

%% convert trial log to trial type
subjid = '004-4012-s1';
filename = fullfile(data_fold,'0_raw-data',subjid,[subjid,'_TrialSequence.txt']);
trialseq = dlmread(filename);

for i = 1:size(trialseq,1)
    filename = fullfile(data_fold,'0_raw-data',subjid,[subjid,'_block00',num2str(i),'.txt']);
    if i < 4
        goggle_cond = 'GogglesOFF';
    else
        goggle_cond = 'GogglesON';
    end

    dlmwrite(filename,goggle_cond,'delimiter','');

    for j = 1:size(trialseq,2)
        if trialseq(i,j) == 0
            trialtype = 'WithObject';
        else
            trialtype = 'WithoutObject';
        end
        dlmwrite(filename,trialtype,'-append','delimiter','');
    end

end

%% get answers from log file
subjid = '004-4012-s1';
filename = fullfile(data_fold,'0_raw-data',subjid,[subjid,'_log.txt']);
fileID = fopen(filename,'r');

answer = [];
blk = 1;
trl = 1;

tline = fgets(fileID);
while ischar(tline)
    tline_cell = split(tline,' ');

    if sum(strcmpi(tline_cell,'Answered')) && strcmpi(tline_cell(3),'object')
        if blk > 6
            tline
        end
        if strcmpi(tline_cell(4),'present')
            answer(blk,trl) = 0;
        else
            answer(blk,trl) = 1;
        end
        trl = trl+1;

        if trl > 10
            blk = blk+1;
            trl = 1;
        end
    end
    tline = fgets(fileID);
end

fclose(fileID);

filename = fullfile(data_fold,'0_raw-data',subjid,[subjid,'_Answer.txt']);
dlmwrite(filename,answer);

%% get confidence rating from log file
subjid = '004-4012-s1';
filename = fullfile(data_fold,'0_raw-data',subjid,[subjid,'_log.txt']);
fileID = fopen(filename,'r');

conf_rate = [];
blk = 1;
trl = 1;

tline = fgets(fileID);
while ischar(tline)
    tline_cell = split(tline,' ');

    if sum(strcmpi(tline_cell,'Rated'))
        conf_rate(blk,trl) = str2num(tline_cell{5}(1));
        trl = trl+1;

        if trl > 10
            blk = blk+1;
            trl = 1;
        end
    end
    tline = fgets(fileID);
end

fclose(fileID);

filename = fullfile(data_fold,'0_raw-data',subjid,[subjid,'_ConfRating.txt']);
dlmwrite(filename,conf_rate);

%% Get trialtypes from Paris patients
subs = {'P1001-4','P1002-2','P1002-3','P1004-2','P1009'};
subs = {'004-4006-s2'};


for sub = 1:length(subs)
    subjid = subs{sub};
    disp(subjid)

    trialtype = [];
    if strcmpi(subjid,'P1002-2')
        blocks = 2:9;
    else
        blocks = 1:8;
    end

    if strcmp(subjid,'004-4010-s1')
        blocks = 1:6;
    end

    for b = 1:length(blocks) % blocks
        filename = fullfile(data_fold,'0_raw-data',subjid,[subjid, ...
            '_block00',num2str(blocks(b)),'.txt']);
        fileID = fopen(filename);

        tline = fgetl(fileID);
        trl = 1;
        while ischar(tline)
            tline = fgetl(fileID);
            if strcmp(tline,'WithObject')
                trialtype(b,trl) = 0;
                trl = trl+1;
            elseif strcmp(tline,'WithoutObject')
                trialtype(b,trl) = 1;
                trl = trl+1;
            end
        end

        fclose(fileID);

    end

    filename = fullfile(data_fold,'0_raw-data',subjid,[subjid, ...
        '_TrialSequence.txt']);
    dlmwrite(filename,trialtype);
end

%% get answer from Paris patients
subs = {'P1001-4','P1002-2','P1002-3','P1009'};

for sub = 1:length(subs)
    subjid = subs{sub};
    disp(subjid)


    filename = fullfile(data_fold,'0_raw-data',subjid,[subjid,'_log.txt']);
    fileID = fopen(filename,'r');

    answer = [];
    blk = 1;
    trl = 1;

    tline = fgets(fileID);
    while ischar(tline)
        tline_cell = split(tline,' ');

        if sum(strcmpi(tline_cell,'Answered')) && strcmpi(tline_cell(3),'object')
            if strcmpi(tline_cell(4),'present')
                answer(blk,trl) = 0;
            else
                answer(blk,trl) = 1;
            end
            trl = trl+1;

            if trl > 10
                blk = blk+1;
                trl = 1;
            end
        end
        tline = fgets(fileID);
    end

    fclose(fileID);

    filename = fullfile(data_fold,'0_raw-data',subjid,[subjid,'_Answer.txt']);
    dlmwrite(filename,answer);

end


%% Get confidence rating from Paris patients
subs = {'P1001-4','P1002-2','P1002-3','P1009'};

for sub = 1:length(subs)
    subjid = subs{sub};
    disp(subjid)

    conf_rate = [];
    blk = 1;
    trl = 1;

    filename = fullfile(data_fold,'0_raw-data',subjid,[subjid,'_log.txt']);
    fileID = fopen(filename,'r');

    tline = fgets(fileID);
    while ischar(tline)
        tline_cell = split(tline,' ');

        if sum(strcmpi(tline_cell,'Rated'))
            conf_rate(blk,trl) = str2num(tline_cell{5}(1));
            trl = trl+1;

            if trl > 10
                blk = blk+1;
                trl = 1;
            end
        end
        tline = fgets(fileID);
    end

    fclose(fileID);

    filename = fullfile(data_fold,'0_raw-data',subjid,[subjid,'_ConfRating.txt']);
    dlmwrite(filename,conf_rate);
end

%% aggregate behav accuracy
subs = {'004-4006-s1','004-4010-s1','004-4010-s2','P1001-4','P1002-2','P1002-3','P1004-2','P1009'};

acc_ON = [];
acc_OFF = [];
for sub = 1:length(subs)
    subjid = subs{sub};
    disp(subjid);

    % compute accuracy for ON vs OFF
    filename = fullfile(data_fold,'0_raw-data',subjid,[subjid,'_TrialSequence.txt']);
    trialtype = dlmread(filename);

    filename = fullfile(data_fold,'0_raw-data',subjid,[subjid,'_Answer.txt']);
    answer = dlmread(filename);

    ON_blk = size(trialtype,1)/2;
    acc = mean(trialtype==answer,2);
    acc_ON(sub,1) = mean(acc(1:ON_blk));
    acc_OFF(sub,1) = mean(acc(ON_blk+1:end));
end


%% all results
subs = {'004-4006-s1','004-4010-s1','004-4010-s2','P1001-4','P1002-2','P1002-3','P1004-2','P1009','004-4006-s2','004-4012-s1'};

allresults = struct();
cnt = 1;

for sub = 1:length(subs)
    subjid = subs{sub};
    disp(subjid)

    trialseq = [];
    answer = [];
    confrate = [];

    filename = fullfile(data_fold,'0_raw-data',subjid,[subjid,'_TrialSequence.txt']);
    if exist(filename,'file')
    trialseq = dlmread(filename);
    end

    filename = fullfile(data_fold,'0_raw-data',subjid,[subjid,'_Answer.txt']);
    if exist(filename,'file')
        answer = dlmread(filename);
    end

    filename = fullfile(data_fold,'0_raw-data',subjid,[subjid,'_ConfRating.txt']);
    if exist(filename,'file')
        confrate = dlmread(filename);
    end

    num_blk = size(trialseq,1);
    num_trl = size(trialseq,2);
    for b = 1:num_blk
        for t = 1:num_trl

            if length(subjid) > 7
                allresults(cnt).Subjid = ['P',subjid(5:8)];
            else
                allresults(cnt).Subjid = ['P',subjid(2:5)];
            end

            allresults(cnt).Task = 'Tumbler';

            allresults(cnt).Block = b;
            allresults(cnt).Trial = t;

            if strcmp(subjid,'P1009') || strcmp(subjid,'004-4012-s2')
                allresults(cnt).SubjectType = 'Blind Control';
            else
                allresults(cnt).SubjectType = 'GenSight Patient';
            end

            if strcmp(subjid,'P1002-3') || strcmp(subjid,'004-4006-s2')
                if b <= num_blk/2
                    allresults(cnt).GogglesCondition = 'GogglesOFF';
                else
                    allresults(cnt).GogglesCondition = 'GogglesON';
                end
            else
                if b <= num_blk/2
                    allresults(cnt).GogglesCondition = 'GogglesON';
                else
                    allresults(cnt).GogglesCondition = 'GogglesOFF';
                end
            end

            if ~isempty(trialseq)
                if trialseq(b,t)==0
                    allresults(cnt).TrialType = 'Present';
                else
                    allresults(cnt).TrialType = 'Absent';
                end
            end

            if ~isempty(answer)
                if answer(b,t)==0
                    allresults(cnt).Answer = 'Present';
                else
                    allresults(cnt).Answer = 'Absent';
                end  
            end

            if ~isempty(answer)
                if trialseq(b,t)==answer(b,t)
                    allresults(cnt).Accuracy = 100;
                    allresults(cnt).Correct = 'Correct';
                else
                    allresults(cnt).Accuracy = 0;
                    allresults(cnt).Correct = 'Incorrect';
                end
            end

            if ~isempty(confrate)
                allresults(cnt).ConfRating = confrate(b,t);
            end

            cnt = cnt+1;
        end
    end
end

filename = fullfile(result_fold,'BehavResults_AllSub_withBCP.xlsx');
writetable(struct2table(allresults),filename);







