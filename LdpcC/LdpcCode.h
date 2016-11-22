//
// Created by Saurabh on 11/4/16.
//

#ifndef LDPCC_LDPCCODE_H
#define LDPCC_LDPCCODE_H

#include <vector>

class LdpcCode {

private:

    std::vector<std::vector<uint8_t > > _H_mat;
    unsigned _N;
    unsigned _K;
    unsigned _M;
    unsigned _Z;
    std::vector<std::vector<unsigned>> _column_mat;
    std::vector<std::vector<unsigned>> _row_mat;

    void generate_compact_rep();

    void lifted_ldpc(std::vector<std::vector<int>> baseH);

public:

    void generate_gallagher_ldpc();

    bool check_codeword(std::vector<uint8_t>);

    void load_wifi_ldpc(unsigned block_length, unsigned rate_index);

    unsigned get_info_length() {return _K;};

    LdpcCode(unsigned block_length, unsigned info_length): _N(block_length), _K(info_length), _Z(0) {
        _M = _N - _K;
        _H_mat.resize(_N);
        for (unsigned j = 0; j < _N ; ++j ) {
            _H_mat.at(j).resize(_M);
        }
    };

    std::vector<uint8_t> encode(std::vector<uint8_t> info_bits);

    std::vector<uint8_t> decode(std::vector<double> llr, unsigned max_iter, bool min_sum);

};


#endif //LDPCC_LDPCCODE_H
