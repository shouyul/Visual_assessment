function [trial_power, power_time, dB_time, freq] = convert_to_power_VEP(trials)

srate = 250;
[~,freq] = cwt(trials(1,:,1),srate);

for ch = 1:size(trials,1)
    disp(ch)
    temp_wavelet = nan(length(freq),size(trials,2),size(trials,3));
    
    parfor t = 1:size(trials,3)
        temp_wavelet(:,:,t) = cwt(trials(ch,:,t),srate);
    end

    temppower = abs(temp_wavelet);
    baseline = mean(mean(temppower,3),2);
    
    trial_power(ch,:,:) = squeeze(mean(temppower,2));

    dB_time(ch,:,:) = 10*log10(mean(temppower,3)./baseline);
    power_time(ch,:,:) = squeeze(mean(temppower,3));
end