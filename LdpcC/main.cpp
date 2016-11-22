#include <iostream>
#include "LdpcCode.h"
#include "Constellation.h"
#include <random>
#include <iomanip>      // std::setprecision

using namespace std;

int main() {

    unsigned max_err = 250;
    unsigned max_runs = 10000;

    bool min_sum = false;

    unsigned code_index = 0;

    LdpcCode ldpc_code(0, 0);

    std::vector<unsigned> block_length_vec{648, 1296, 1944};

    std::vector<unsigned> constellations{1, 2, 3};

    std::vector<unsigned> rate_index_vec{0, 1, 2};

    // Used to run EbNo points of interest for the rate used
    std::vector<unsigned>  rate_offset{0, 0, 1};

    // Used to run EbNo points of interest as per the constellation used
    std::vector<double>  constellations_ebno_offset{0.0, 3.0, 7.0};

    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

    for (unsigned i_const = 0; i_const < constellations.size(); ++i_const ) {

        unsigned n_bits = constellations.at(i_const);

        Constellation modulation(n_bits);
        double ebno_log_min = 0.5 + constellations_ebno_offset.at(i_const);
        double ebno_log_max = 4.51 + constellations_ebno_offset.at(i_const);
        double ebno_log_increment = 0.25;
        std::vector<double> ebno_vec;

        for (double ebno_log = ebno_log_min; ebno_log <= ebno_log_max; ebno_log += ebno_log_increment)
            ebno_vec.push_back(ebno_log);

        for (unsigned i_block = 0; i_block < block_length_vec.size(); ++i_block) {
            unsigned block_length = block_length_vec.at(i_block);
            std::cout << "% Processing block length " << block_length << std::endl;
            for (unsigned i_rate = 0; i_rate < rate_index_vec.size(); ++i_rate) {
                unsigned rate_index = rate_index_vec.at(i_rate) + rate_offset.at(i_const);
                std::cout << "% Processing rate index " << rate_index << std::endl;
                ldpc_code.load_wifi_ldpc(block_length, rate_index);
                unsigned info_length = ldpc_code.get_info_length();

                std::vector<double> bler(ebno_vec.size(), 0);
                std::vector<double> num_err(ebno_vec.size(), 0);
                std::vector<double> num_run(ebno_vec.size(), 0);

                std::vector<double> noise(block_length / n_bits, 0);
                std::vector<double> received_signal(block_length / n_bits, 0);

                std::normal_distribution<double> gauss_dist(0.0f, 1.0);
                std::default_random_engine generator;

                for (unsigned run = 0; run < max_runs; ++run) {
                    if ((run % (max_runs / 5)) == (max_runs / 5 - 1)) {
                        std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
                        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
                        std::cout << "% Running iteration " << run << "; time elapsed = " << duration / 1000 / 1000
                                  << " seconds"
                                          "; percent complete = " << (100 * run) / max_runs << "." << std::endl;
                    }
                    for (unsigned i = 0; i < block_length / n_bits; ++i) {
                        noise.at(i) = (double) gauss_dist(generator);
                    }

                    std::vector<int> scrambling_bits(block_length, 0);

                    srand(run);

                    std::vector<uint8_t> info_bits(info_length, 0);
                    for (unsigned i_bit = 0; i_bit < info_length; ++i_bit)
                        info_bits.at(i_bit) = (uint8_t) (rand() % 2);

                    std::vector<uint8_t> coded_bits = ldpc_code.encode(info_bits);

                    std::vector<uint8_t> scrambled_bits(block_length, 0);

                    for (unsigned i_bit = 0; i_bit < block_length; ++i_bit) {
                        scrambling_bits.at(i_bit) = (rand() % 2);
                        scrambled_bits.at(i_bit) = (uint8_t) ((coded_bits.at(i_bit) + scrambling_bits.at(i_bit)) % 2);
                    }

                    std::vector<double> mod_sym = modulation.modulate(scrambled_bits);

                    std::vector<bool> prev_decoded(ebno_vec.size(), false);

                    for (unsigned i_ebno = 0; i_ebno < ebno_vec.size(); ++i_ebno) {

                        if (num_err.at(i_ebno) > max_err)
                            continue;

                        num_run.at(i_ebno)++;

                        bool run_sim = true;

                        for (unsigned i_ebno2 = 0; i_ebno2 < i_ebno; ++i_ebno2) {
                            if (prev_decoded.at(i_ebno2)) {
                                //  This is a hack to speed up simulations -- it assumes that this run will be decoded
                                // correctly since it was decoded correctly for a lower EbNo
                                run_sim = false;
                            }
                        }

                        if (!run_sim) {
                            continue;
                        }

                        double snr_linear = std::pow(10.0, ebno_vec.at(i_ebno) / 10)
                                            * ((double) info_length) / ((double) (block_length)) * ((double) n_bits);

                        double N_0 = 0.5 / snr_linear;
                        for (unsigned i = 0; i < block_length / n_bits; ++i) {
                            received_signal.at(i) = mod_sym.at(i) + std::sqrt(N_0) * noise.at(i);
                        }
                        std::vector<double> llr = modulation.llr_compute(received_signal, N_0);

                        for (unsigned i_bit = 0; i_bit < block_length; ++i_bit)
                            llr.at(i_bit) = llr.at(i_bit) * (1 - 2 * scrambling_bits.at(i_bit));

                        std::vector<uint8_t> decoded_cw = ldpc_code.decode(llr, 20, min_sum);
                        bool error = false;
                        for (unsigned i_bit = 0; i_bit < info_length; ++i_bit) {
                            if (decoded_cw.at(i_bit) != info_bits.at(i_bit)) {
                                error = true;
                                break;
                            }
                        }

                        if (!error) {
                            prev_decoded.at(i_ebno) = true;
                        } else {
                            num_err.at(i_ebno)++;
                        }

                    } // EbNo Loop End

                } // run loop end

                for (unsigned i_ebno = 0; i_ebno < ebno_vec.size(); ++i_ebno) {
                    bler.at(i_ebno) = num_err.at(i_ebno) / num_run.at(i_ebno);
                }

                code_index++;
                std::cout << "% Block length = " << block_length << "; info length = " << info_length;
                std::cout << "; bits/symbols = " << n_bits;
                std::cout << "; min-sum algorithm = " << std::boolalpha << min_sum << "." << std::endl;
                std::cout << "bler(:, :,  " << code_index << ") = [ ..." << std::endl;
                for (unsigned i_ebno = 0; i_ebno < ebno_vec.size(); ++i_ebno) {
                    std::cout << std::fixed << "\t" << std::setprecision(3) << ebno_vec.at(i_ebno) << ",\t";
                    std::cout << std::fixed << std::setprecision(6) << bler.at(i_ebno) << ";\t";
                    std::cout << std::endl;
                }
                std::cout << " ];" << std::endl;

            } // Rate loop end

        } // Block length loop end

    } // Constellation loop end

    std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
    std::cout << "% Total time taken = " << duration / 1000 / 1000 << " seconds." << std::endl;
    return 0;
}