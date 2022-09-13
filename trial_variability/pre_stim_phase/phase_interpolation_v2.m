%% Load the data and the info
clear all; close all; clc
locations = cceps_files;

savedir = fullfile(locations.results_folder,'prestim_peakfreq_EJC'); 
mkdir(savedir);

%subjects = {'211','212','213','214','216','218','219','222','224','225'};
%subjects = {'222','224','225'};

subjects = {'HUP222'}; %locations.subjects;
do_plot = false;
phase_time = 0.005; % time of phase value extraction relative to end of epoch

for s = 1:length(subjects)
    
    measure_locs = [1];
    
    for n = 1:length(measure_locs)

        subject = subjects{s};
        disp(subject);

        load(fullfile(locations.results_folder,'prestim_metrics/',[subject,'_CCEPPrestimData.mat'])); % load the data

        %load(strcat('prestim_info/results_',subject,'_CCEP.mat'))

        peak_locs = load(fullfile(locations.results_folder,'powercalc_peakfreq_AL',['peaks_',subject,'.csv']));
        
        load(fullfile(locations.results_folder,['results_',subject,'_CCEP.mat']));
        %fs = out.other.stim.fs; % Get the sampling frequency

        fs = out.other.stim.fs;
        n_electrodes = length(data.processed);

        % if fs==1024
        %     choose = 20
        % else
        %     choose = 10
        % end

        % Define the high theta filter
        % fc1 = 4;
        % fc2 = 8;
        % 
        % [b,a] = butter(2,[fc1/(fs/2) fc2/(fs/2)]);

        %% Calculate the features at the peak frequency


        % main loop
        features_cell = cell([n_electrodes 3]);

        for k = 1:n_electrodes

            fprintf('Processing Stim by Electrode: %d, %s\n',k, datestr(now,'HH:MM:SS.FFF'))
            if ~isempty(data.processed{k})

                [~,n_trials,n_electrodes2] = size(data.processed{k});
                
                avg_power_list = nan(n_electrodes,n_trials);
                phase_list = nan(n_electrodes,n_trials);
                power_list = nan(n_electrodes,n_trials);

                for i = 1:n_electrodes2
                    for j = 1:n_trials
                        x_orig = data.processed{k}(:,j,i);
                        %x_bp=bandpass(x_orig,[8 12],fs);
                        if ~any(isnan(x_orig)) && ~isnan(peak_locs(i))

                            % Define the iirfilter at the peak frequency with bw of 0.5
                            wo = peak_locs(i)/(fs/2);  
                            bw = 0.5/(fs/2);
                            [b,a] = iirpeak(wo,bw);
                            
                            x_bp = filtfilt(b ,a, x_orig);

        %                     % Replace x_bp by a sine wave
        %                                         %%Time specifications:
%                              Fs = fs;                   % samples per second
%                              dt = 1/Fs;                   % seconds per sample
%                              StopTime = length(x_bp)/fs;             % seconds
%                              t = (0:dt:StopTime)';     % seconds
%                              %%Sine wave:
%                              Fc = peak_locs(i);                     % hertz
%                              x_bp = cos(2*pi*Fc*t);


                            % Compute the phase
                            y = hilbert(x_bp);
                            sigphase = atan2(imag(y),real(y)); % get the phase

                            power = x_bp.^2; % instantaneous power. use this result with caution given the filter edge effects

                            avg_power = bandpower(x_orig, fs, [peak_locs(i)-0.25,peak_locs(i)+0.25]); % this result uses the fft directly

                            % location to extract the phase and instantaneous power
                            measure_loc = round(length(x_orig) + 1 + fs*phase_time);
                            start_loc = floor(measure_loc/2); % start extrapolating phase from halfway through epoch
                            %measure_loc = measure_locs(n);

                            % mathematically interpolate the phase
                            slope = (peak_locs(i)*2*pi)/fs;
                            % measure_loc = round(length(x_orig) + 1 + fs*0.005);
                            phase_val = sigphase(start_loc)+pi+((measure_loc-start_loc)*slope);
                            phase_val = (phase_val - floor(phase_val/(2*pi))*(2*pi))-pi;

                            % can't extrapolate from phase at t = 1
                            extrap_phase = EXTRAPOLATE_PHASE(256:512,sigphase(256),256,peak_locs(i),fs);
                            
                            avg_power_list(i,j) =  x_bp(1);
                            power_list(i,j) = power(length(x_orig));
                            phase_list(i,j) = sigphase(1);
                            
                            
                        if do_plot && i<2
                            f=figure;
                            subplot(3,1,1);
                            plot(x_orig);
                            xticklabels(get(gca,'XTick')/fs);
                            title('Original')
                            subplot(3,1,2);
                            plot(x_bp);
                            xticklabels(get(gca,'XTick')/fs);
                            title('IIR Filtered')
                            subplot(3,1,3);
                            plot(sigphase); hold on;
                            plot(256:512,extrap_phase);
                            xticklabels(get(gca,'XTick')/fs);
                            title('Phase')
                        end

                        else
                            avg_power_list(i,j) =  NaN;
                            power_list(i,j) = NaN;
                            phase_list(i,j) = NaN;
                        end
                    end
                end
            else
                avg_power_list = [];
                power_list = [];
                phase_list =[];
            end

            features_cell{k,1}=avg_power_list;
            features_cell{k,2}=power_list;
            features_cell{k,3}=phase_list;

        end


data.features = features_cell;

data.feat_name = {'avg_power','inst_power','inst_phase'};
data.band_names = {'peak freq','peak freq','peak freq'};
data.peak_locs = peak_locs;

% remove the useless fields
data = rmfield(data, 'processed');
data = rmfield(data, 'raw');
data = rmfield(data, 'reshape_fxn');

save(fullfile(savedir,[subject,'_CCEPPrestimFeatures_PeakFreq_',num2str(phase_time),'.mat']), 'data')

end
end
