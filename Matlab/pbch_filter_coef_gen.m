function coef_pbch = pbch_filter_coef_gen(fs)

if fs == 1.92e6
    coef_pbch = fir1(46, (0.18e6*6+150e3)/fs); %freqz(coef_pbch, 1, 1024);
elseif fs == 3.84e6
    coef_pbch = fir1(126, (0.18e6*6+150e3)/fs); %freqz(coef_pbch, 1, 1024);
elseif fs == 7.68e6
    coef_pbch = fir1(126, (0.18e6*6+150e3)/fs); %freqz(coef_pbch, 1, 1024);
elseif fs == 19.2e6
    coef_pbch = fir1(126, (0.18e6*6+150e3)/fs); %freqz(coef_pbch, 1, 1024);
elseif fs == 30.72e6
    coef_pbch = fir1(126, (0.18e6*6+150e3)/fs); %freqz(coef_pbch, 1, 1024);
else
    disp('fs not supported.');
end
