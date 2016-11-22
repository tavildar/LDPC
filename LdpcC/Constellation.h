//
// Created by Saurabh on 11/17/16.
//

#ifndef LDPCC_CONSTELLATION_H
#define LDPCC_CONSTELLATION_H

#include <vector>

class Constellation {
private:

    std::vector<double> _points;
    unsigned _n_bits;
    unsigned _n_syms;

    std::vector<std::vector<int>> _bit_sym_map;

public:

    Constellation(unsigned _n_bits);

    std::vector<double> modulate(std::vector<uint8_t> coded_bits);

    std::vector<double> llr_compute(std::vector<double> y, double n_0);

};

#endif //LDPCC_CONSTELLATION_H
