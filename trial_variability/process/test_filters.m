
f_w = 20;
am = 1;
srate = 1000;
t = linspace(0,1,srate); 
x = am*sin(2*pi*f_w*t) + t*10 + 100;
f=figure;

subplot(1,3,1);
plot(x);
title('20 hz sin wave w/ lin trend');
xlabel('time (ms)');
ylabel('uV');


subplot(1,3,2);
plot(detrend(x));
title('20 hz sin wave detrended');
xlabel('time (ms)');
ylabel('uV');
ylim([-2*am 2*am]);

subplot(1,3,3);
[b,a] = butter(4, 1/(srate/2), 'high');
x_f = filtfilt(b,a,x);
plot(x_f);
title('20 hz sin wave 1 hz high pass');
xlabel('time (ms)');
ylabel('uV');
ylim([-2*am 2*am]);


%%
x = trials_sch_raw(:,t,:); % from get_stim_trials_pt.m
d = x(:,1);
FFT_PLOT(d,out.other.stim.fs);
d_f = FILTER_CCEPS(d',out.other.stim.fs);
FFT_PLOT(d_f,out.other.stim.fs);