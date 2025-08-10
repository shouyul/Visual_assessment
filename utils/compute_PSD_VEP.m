function [x_psd_base,x_psd_stim,x_psd_ratio,freq] = compute_PSD_VEP(results)

srate = 250;
sub_chanlocs = results.chanlocs;
trials = results.trials;

eeg_ind = [];
for ch = 1:length(sub_chanlocs)

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
% x_psd_ratio = squeeze(mean(x_psd_ratio));

end
