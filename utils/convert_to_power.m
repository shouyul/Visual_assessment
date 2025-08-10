function [trial_dB,dB_time,trial_pow_time,trialtype] = convert_to_power(subjid)

% load data
main_fold = '/Users/shouyuling/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Postdoc/Projects/Gensight/Data/VEPTask';
data_fold = fullfile(main_fold,'Analysis');
trial_len = 10;

filename = fullfile(data_fold,'2_preprocessing','bemobil','APP',...
    strcat(subjid,'_cleanedForICA.set'));
EEG = pop_loadset(filename);

% segment
trials = [];
Info = []; % 1st - trial type, 2nd - block #, 3rd - trial #

event = EEG.event;
for e = 1:length(event)
    if strcmp(event(e).type,'TrialStart')
        start_ind = round(event(e).latency);
        end_ind = start_ind + EEG.srate * trial_len-1;
        if end_ind > size(EEG.data,2)
            break
        end
        current_trial = EEG.data(:,start_ind:end_ind);
        %                 current_trial = current_trial - nanmean(current_trial,2);
        current_trial = detrend(current_trial')';
        trials = cat(3,trials,current_trial);

        if isfield(event,'trialtype')
            switch event(e).trialtype
                case 'Control'
                    trialtype = 1;
                case 'Disc'
                    trialtype = 2;
                case {'90Bar','VerticalBar'}
                    trialtype = 3;
                case {'0Bar','HorizontalBar'}
                    trialtype = 4;
                otherwise
                    trialtype = 0;
            end
        elseif isfield(event,'stimulus')
            switch event(e).stimulus
                case 'Control'
                    trialtype = 1;
                case {'Disc','3Rings','5Rings'}
                    trialtype = 2;
                case {'90Bar','-45Bar'}
                    trialtype = 3;
                case {'0Bar','45Bar'}
                    trialtype = 4;
                otherwise
                    trialtype = 0;
            end
        end

        Info = cat(1,Info,[trialtype, event(e).block, event(e).trial]);

    end
end

EEG.data = trials;
EEG.trials = size(trials,3);
EEG.times = [0:1/EEG.srate:10-1/EEG.srate];
EEG.pnts = size(trials,2);


min_freq = 1;
max_freq = 80;
num_freq = 40;
% define wavelet parameters
time = -1:1/EEG.srate:1;
freq = logspace(log10(min_freq),log10(max_freq),num_freq);
s    = logspace(log10(3),log10(10),num_freq)./(2*pi*freq);

% definte convolution parameters
n_trials = size(EEG.data,3);
n_wavelet            = length(time);
n_data               = EEG.pnts*n_trials;
n_convolution        = n_wavelet+n_data-1;
n_conv_pow2          = pow2(nextpow2(n_convolution));
half_of_wavelet_size = (n_wavelet-1)/2;
times = [0:1/EEG.srate:10-1/EEG.srate];

baseidx = dsearchn(times',[0 10]');

eegpower = zeros(size(EEG.data,1),num_freq,EEG.pnts);
trial_dB = zeros(size(EEG.data,1),num_freq,n_trials);
trial_pow_time = zeros(size(trials,1),num_freq,size(trials,2)/5,n_trials);

for ch = 1:size(EEG.data,1)

    eegfft = fft(reshape(EEG.data(ch,:,:),1,EEG.pnts*n_trials),n_conv_pow2);
    for fi=1:num_freq

        wavelet = fft( sqrt(1/(s(fi)*sqrt(pi))) * exp(2*1i*pi*freq(fi).*time) .* exp(-time.^2./(2*(s(fi)^2))) , n_conv_pow2 );

        % convolution
        eegconv = ifft(wavelet.*eegfft);
        eegconv = eegconv(1:n_convolution);
        eegconv = eegconv(half_of_wavelet_size+1:end-half_of_wavelet_size);

        %                 temppower = mean(abs(reshape(eegconv,EEG.pnts,EEG.data)).^2,2);
        temppower = abs(reshape(eegconv,EEG.pnts,n_trials)).^2;

%         temppower_avg = mean(temppower,2);
%         eegpower(ch,fi,:) = 10*log10(temppower_avg./mean(temppower_avg(baseidx(1):baseidx(2))));
%         trial_dB(ch,fi,:) =  10*log10(mean(temppower./mean(temppower_avg(baseidx(1):baseidx(2)))));
        temppower = downsample(temppower,5);
%         trial_pow_time(ch,fi,:,:) = temppower;
        chan_pow(fi,:,:) = temppower;
    end

    control_trial = Info(:,1) == 1;
%     baseline = mean(mean(mean(chan_pow(:,:,control_trial),3),2)); % average across frequency
    baseline = mean(mean(chan_pow(:,:,:),3),2);

    dB_time(ch,:,:) = 10*log10(mean(chan_pow,3)./baseline);
    trial_dB(ch,:,:) = 10*log10(squeeze(mean(chan_pow,2))./baseline);
    trial_pow_time(ch,:,:,:) = chan_pow./baseline;
end

trialtype = Info(:,1);

end