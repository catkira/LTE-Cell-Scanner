% Xianjun Jiao (putaoshu@msn.com)
% CellSearch.m
% Improved LTE-Cell-Scanner (written by James Peroulas: https://github.com/Evrytania/LTE-Cell-Scanner).

% Some scripts are borrowed from:
% https://github.com/JiaoXianjun/rtl-sdr-LTE
% https://github.com/Evrytania/LTE-Cell-Scanner
% https://github.com/Evrytania/Matlab-Library
% https://github.com/JiaoXianjun/multi-rtl-sdr-calibration

function [r_pbch_sub, r_full_sub, cell_info_return] = CellSearch(r_pbch, r_full, f_search_set, fc, fs)
r_full_sub = -1;
skip_TDD = 1;

% f_search_set = 20e3:5e3:30e3; % change it wider if you don't know pre-information
%pss_fo_set = pss_fo_set_gen(td_pss, f_search_set);

num_radioframe = 8; % each radio frame length 10ms. MIB period is 4 radio frame

pbch_sampling_ratio = 16;
sampling_rate = 30.72e6;
sampling_rate_pbch = sampling_rate/pbch_sampling_ratio; % LTE spec. 30.72MHz/16.

num_subframe_per_radioframe = 10;
len_time_subframe = 1e-3; % 1ms. LTE spec
num_sample_per_radioframe = num_subframe_per_radioframe*len_time_subframe*sampling_rate_pbch;
num_sample_pbch = num_radioframe*num_sample_per_radioframe;

DS_COMB_ARM = 0; % disable combining for now
FS_LTE = 30720000;
thresh1_n_nines=12;
rx_cutoff=(6*12*15e3/2+4*15e3)/(FS_LTE/16/2);
THRESH2_N_SIGMA = 3;
ex_gain = 2;

num_try = floor(length(r_pbch)/num_sample_pbch);
%% Display overview grid
figure(3); show_time_frequency_grid_according_pss(400, 1, r_full, fs);
title("overview")
%%
for try_idx = 1 : num_try
    disp(['Try idx ' num2str(try_idx)]);

    sp = (try_idx-1)*num_sample_pbch + 1;
    ep = sp + num_sample_pbch - 1;
    capbuf_pbch = r_pbch(sp:ep).';
    
    capbuf_pbch = capbuf_pbch - mean(capbuf_pbch);
    r_pbch_sub = capbuf_pbch;

    disp(['Input averaged abs: ' num2str( mean(abs([real(capbuf_pbch) imag(capbuf_pbch)])) )]);

    disp('sampling_ppm_f_search_set_by_pss: try ... ... ');
    xc = sampling_ppm_f_search_set_by_pss(capbuf_pbch.', f_search_set);
    
    %% Display grid 
    sp_full = (sp-1)*fs/sampling_rate_pbch + 1;
    ep_full = ep*fs/sampling_rate_pbch;
    if ep_full > size(r_full,1); ep_full = size(r_full,1); end
    r_full_sub = r_full(sp_full:ep_full);
    %time_displayed = 21e-3; % 21ms
    time_displayed = size(r_full_sub,1)/fs; % show everything
    figure(2); show_time_frequency_grid_according_pss(0, 1, r_full_sub(1 : time_displayed*fs), fs); % plot 21ms
    title("current segment")

    %% xcor_pss only returns first 9600 positions, which is about 5 ms (half a frame)
    [xc_incoherent_collapsed_pow, xc_incoherent_collapsed_frq, n_comb_xc, ~, ~, ~, sp_incoherent, ~]= ...
    xcorr_pss(capbuf_pbch,f_search_set,DS_COMB_ARM,fc,xc);

    R_th1=chi2inv(1-(10.0^(-thresh1_n_nines)), 2*n_comb_xc*(2*DS_COMB_ARM+1));
    Z_th1=ex_gain.*R_th1*sp_incoherent/rx_cutoff/137/2/n_comb_xc/(2*DS_COMB_ARM+1);

    figure(1);
    subplot(3,1,1); hold off; plot(xc_incoherent_collapsed_pow(1,:)); title('n_id_2=0', 'Interpreter', 'none'); drawnow;
    subplot(3,1,2); hold off; plot(xc_incoherent_collapsed_pow(2,:)); title('n_id_2=1', 'Interpreter', 'none'); drawnow;
    subplot(3,1,3); hold off; plot(xc_incoherent_collapsed_pow(3,:)); title('n_id_2=2', 'Interpreter', 'none'); drawnow;
    peaks=peak_search(xc_incoherent_collapsed_pow,xc_incoherent_collapsed_frq,Z_th1,f_search_set,fc);

    for j = 1 : length(peaks)
        peaks(j).extra_info.num_peaks_raw = length(peaks)/2;
        subplot(3,1,peaks(j).n_id_2+1); hold on; plot(peaks(j).ind, xc_incoherent_collapsed_pow(peaks(j).n_id_2+1,peaks(j).ind),'r*'); drawnow;
    end

    % aligned pbch plot if a peak was found
    if ~isempty(peaks)
        figure(4); show_time_frequency_grid_according_pss(peaks(1).ind, peaks(1).k_factor, capbuf_pbch, 1.92e6); drawnow;
        title("pbch")
        figure(5)
        plot(abs(r_pbch(peaks(1).ind-500:peaks(1).ind+100)));
        
    end
    
    disp(['Found ' num2str(length(peaks)/2) ' peaks']);
    tdd_flags = kron(ones(1, length(peaks)/2), [0 1]); % even: tdd_flag 0; odd : tdd_flag 1
    
    detect_flag = zeros(1, length(peaks));
    figure(1);
    for i=1:length(peaks)
        disp(['try peak ' num2str(floor((i+1)/2)) ' at offset ' num2str(peaks(i).ind/((30.72e6/16)/fs))]);
        tdd_flag = tdd_flags(i);
        peak = sss_detect(peaks(i),capbuf_pbch,THRESH2_N_SIGMA,fc);
        if ~isnan( peak.n_id_1 )
            if skip_TDD == 1 && tdd_flag == 1
                continue;
            end
            peak=pss_sss_foe(peak,capbuf_pbch,fc,tdd_flag);
            [tfg, tfg_timestamp]=extract_tfg(peak,capbuf_pbch,fc, 6);
            [tfg_comp, ~, peak]=tfoec(peak,tfg,tfg_timestamp,fc, 6);
            peak=decode_mib(peak,tfg_comp);
            if isnan( peak.n_rb_dl)
                continue;
            end
            figure(1)
            subplot(3,1,peaks(j).n_id_2+1); hold on; plot(peak.ind, xc_incoherent_collapsed_pow(peaks(j).n_id_2+1,peak.ind),'g*'); drawnow;            
            if tdd_flag == 1
                disp('  Detected a TDD cell!');
            else
                disp('  Detected a FDD cell!');
            end
            disp(['    cell ID: ' num2str(peak.n_id_cell)]);
            disp(['    PSS  ID: ' num2str(peak.n_id_2+1)]);
            disp(['    RX power level: ' num2str(10*log10(peak.pow))]);
            disp(['    residual frequency offset: ' num2str(peak.freq_superfine)]);
            peaks(i) = peak;
            detect_flag(i) = 1;
        end
    end
    if sum(detect_flag)==0
        disp('No LTE cells were found...');
    else
        break;
    end
end

cell_info_return = peaks;

% show all Cells information
disp(' ');
disp('-------------------------------Cells information summary-------------------------------');

if isempty(detect_flag)
    disp('No valid PSS is found at pre-proc phase! Please try again.');
else
    if sum(detect_flag)
        hit_idx = find(detect_flag);
        for i=1:length(hit_idx);
            peak = peaks(hit_idx(i));
            
            if peak.duplex_mode == 1
                cell_mode_str = 'TDD';
            else
                cell_mode_str = 'FDD';
            end
            disp(['Cell ' num2str(i) ' information:--------------------------------------------------------']);
            disp(['               Offset: ' num2str(peak.ind/((30.72e6/16)/fs))]);
            disp(['            Cell mode: ' num2str(cell_mode_str)]);
            disp(['              Cell ID: ' num2str(peak.n_id_cell)]);
            disp(['   Num. eNB Ant ports: ' num2str(peak.n_ports)]);
            disp(['    Carrier frequency: ' num2str(fc/1e6) 'MHz']);
            disp(['Residual freq. offset: ' num2str(peak.freq_superfine/1e3) 'kHz']);
            disp(['       RX power level: ' num2str(10*log10(peak.pow))]);
            disp(['              CP type: ' peak.cp_type]);
            disp(['              Num. RB: ' num2str(peak.n_rb_dl)]);
            disp(['       PHICH duration: ' peak.phich_dur]);
            disp(['  PHICH resource type: ' num2str(peak.phich_res)]);
        end
    else
        disp('No LTE cells were found...  Please try again.');
    end
end
