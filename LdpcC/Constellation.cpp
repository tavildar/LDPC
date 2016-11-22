//
// Created by Saurabh on 11/17/16.
//

#include "Constellation.h"
#include <iostream>
#include <math.h>       /* sqrt */

std::vector<double> Constellation::modulate(std::vector<uint8_t> coded_bits) {

    std::vector<double> mod_sym(coded_bits.size()/_n_bits, 0.0);

    for (unsigned c_bit = 0; c_bit < coded_bits.size(); c_bit = c_bit + _n_bits ) {
        unsigned sym = 0;
        for (unsigned i_bit = c_bit; i_bit < c_bit + _n_bits; ++i_bit)
            sym += coded_bits.at(i_bit) << (i_bit - c_bit);

        mod_sym.at(c_bit / _n_bits) = _points.at(sym);

    }
    return mod_sym;
}

std::vector<double> Constellation::llr_compute(std::vector<double> y, double n_0) {

    std::vector<double> p_0(y.size() * _n_bits, 0.0);
    std::vector<double> p_1(y.size() * _n_bits, 0.0);
    std::vector<double> llr(y.size() * _n_bits, 0.0);

    for (unsigned i_y = 0; i_y < y.size(); ++i_y ) {
        for (unsigned i_sym = 0; i_sym < _n_syms; ++i_sym ) {
            double p_sym = exp(-(y.at(i_y) - _points.at(i_sym)) * (y.at(i_y) - _points.at(i_sym))/2/n_0);
            for (unsigned i_bit = 0; i_bit < _n_bits; ++i_bit ) {
                if ( _bit_sym_map.at(i_sym).at(i_bit) == 1 )
                    p_1.at(i_y * _n_bits + i_bit) += p_sym;
                else
                    p_0.at(i_y * _n_bits + i_bit) += p_sym;
            }
        }
    }
    for (unsigned i_bit = 0; i_bit < llr.size(); ++i_bit)
        llr.at(i_bit) = log(p_0.at(i_bit)/p_1.at(i_bit));

    return llr;

}

Constellation::Constellation(unsigned n_bits): _n_bits(n_bits)  {

    _n_syms = (unsigned) (1 << _n_bits);

    _bit_sym_map.resize(_n_syms);
    for (unsigned i_sym = 0; i_sym < _n_syms; ++i_sym)
        _bit_sym_map.at(i_sym).resize(_n_bits);

    for (unsigned i_sym = 0; i_sym < _n_syms; ++i_sym ) {
        unsigned tmp_sym = i_sym;
        for (unsigned i_bit = 0; i_bit < _n_bits; ++i_bit ) {
            if ( (tmp_sym % 2) == 1 )
                _bit_sym_map.at(i_sym).at(i_bit) = 1;
            else
                _bit_sym_map.at(i_sym).at(i_bit) = 0;

            tmp_sym = ( tmp_sym - (tmp_sym % 2) ) /2;
        }
    }

    _points.resize(_n_syms);

    switch (_n_bits) {
        case 1: // BPSK
            _points = {-1.0, 1.0};
            break;
        case 2: // ASK-4 Gray mapping
            _points = {-3.0, -1.0, 3.0, 1.0};
            break;
        case 3: // ASK-8 reflected Gray
            _points = {-7.0, -5.0, -1.0, -3.0, 7.0, 5.0, 1.0, 3.0};
            break;
        default:
        std::cout << "Constellation not supported" <<std::endl;
            break;
    }

    double energy = 0;
    for (unsigned i_sym = 0; i_sym < _n_syms; ++i_sym ) {
        energy += _points.at(i_sym) * _points.at(i_sym);
    }

    energy = energy/_n_syms;

    for (unsigned i_sym = 0; i_sym < _n_syms; ++i_sym ) {
        _points.at(i_sym) =  _points.at(i_sym) / sqrt(energy);
    }

}
