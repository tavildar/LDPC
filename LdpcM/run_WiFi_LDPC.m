

% run WiFi LDPC
clear all;
%  close all;
max_runs = 100;
max_decode_iterations = 20;
ldpc_code = LDPCCode(0, 0);
min_sum = 1;
n_0 = 1/2;


block_length = 1296; % Should be one of 648, 1296, and 1944
rate = 3/4; % Should be one of 1/2, 2/3, 3/4, 5/6

constellation_name = 'ask8'; %Should be one of 'bpsk', 'ask4', 'ask8'
modulation = Constellation(constellation_name);


ebno_db_vec = (0:0.5:5) + 7;

num_block_err = zeros(length(ebno_db_vec), 1);

tic

ldpc_code.load_wifi_ldpc(block_length, rate);
info_length = ldpc_code.K;
disp(['Running LDPC with N = ', num2str(block_length), ' and rate = ' , num2str(rate), ' with constellation = ', constellation_name]);

snr_db_vec = ebno_db_vec + 10*log10(info_length/block_length) + 10*log10(modulation.n_bits);

for i_run = 1 : max_runs
    
    
    if ( mod(i_run, max_runs/10) == 1)
        disp(['Current run = ', num2str(i_run), ' percentage complete = ', num2str((i_run-1)/max_runs * 100), '%', ' time elapsed = ', num2str(toc), ' seconds']);
    end
    
    noise = sqrt(n_0) * randn(block_length/modulation.n_bits, 1);
    
    info_bits = rand(info_length, 1) < 0.5;
    
    coded_bits = ldpc_code.encode_bits(info_bits);
    
    scrambling_bits = (rand(block_length, 1) < 0.5);
    
    scrambled_bits = mod(coded_bits + scrambling_bits, 2);
    
    x = modulation.modulate(scrambled_bits);
    
    for i_snr = 1 : length(snr_db_vec)
        
        snr_db = snr_db_vec(i_snr);
        
        snr = 10^(snr_db/10);
        
        y =  x + noise/sqrt(snr);
        
        [llr, ~] = modulation.compute_llr(y, n_0/snr);
        
        llr = llr .* ( 1 - 2 * scrambling_bits);
        
        [decoded_codeword, ~] = ldpc_code.decode_llr(llr, max_decode_iterations, min_sum);
        
        if any(decoded_codeword ~= coded_bits)
            num_block_err(i_snr) = num_block_err(i_snr) + 1;
        else
            break; % Assumes the codeword will be decoded correctly for a higher SNR as well
        end
        
    end
end

bler = num_block_err/max_runs;

figure;
semilogy(ebno_db_vec, bler);
hold on; grid on;
legend(constellation_name);





