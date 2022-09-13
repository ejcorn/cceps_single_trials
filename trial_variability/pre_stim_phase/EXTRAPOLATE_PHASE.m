function phase_val = EXTRAPOLATE_PHASE(measure_locs,init_phase,init_phase_idx,signal_freq,fs)

if any(measure_locs<init_phase_idx)
    error('cannot do that');
end

slope = (signal_freq*2*pi)/fs;
% measure_loc = round(length(x_orig) + 1 + fs*0.005);
phase_val = init_phase+pi+((measure_locs-init_phase_idx)*slope);
phase_val = (phase_val - floor(phase_val/(2*pi))*(2*pi))-pi;
