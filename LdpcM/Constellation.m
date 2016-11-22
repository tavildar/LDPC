classdef Constellation < handle
    % Deals with constellation transmit and receive
    
    properties
        
        constellation_points;
        
        n_bits;
        n_sym;
        bit_sym_map;

    end
    
    properties (Constant)
        
        bpsk = [1 -1];
        
        ask4 = [-3, -1, 3, 1]/sqrt(5); % Gray mapping
        
        ask8 = [-7, -5, -1, -3, 7, 5, 1, 3]/sqrt(21); % Reflected Gray mapping
                       
    end
    
    methods
        
        function obj = Constellation ( constellion_name )
            if strcmp (constellion_name, 'bpsk')
                obj.constellation_points = obj.bpsk;
                obj.n_bits = 1;
            elseif strcmp (constellion_name, 'ask4')
                obj.constellation_points = obj.ask4;
                obj.n_bits = 2;
            elseif strcmp (constellion_name, 'ask8')
                obj.constellation_points = obj.ask8;
                obj.n_bits = 3;
            else
                disp('Unsupported constellation');
            end
            
            obj.n_sym = 2^obj.n_bits;

            obj.bit_sym_map = zeros(obj.n_sym, obj.n_bits);
            for sym_index = 0 : obj.n_sym - 1
                a = dec2bin(sym_index, obj.n_bits);
                for bit_index = 1 :obj.n_bits
                    if a(obj.n_bits + 1 - bit_index) == '1'
                        obj.bit_sym_map(sym_index + 1, bit_index) = 1;
                    end
                end
            end

            obj.constellation_points = obj.constellation_points/sqrt(mean(obj.constellation_points.^2));
            
        end
        
        function [sym] = modulate(obj, bits)
            sym = zeros(floor(length(bits)/obj.n_bits), 1);
            for i = 1 : 1 : length(bits)/obj.n_bits
                symbol = 0;
                for j = 1 : obj.n_bits
                    symbol = symbol + 2^(j-1) * bits((i-1)*obj.n_bits+j);
                end
                sym(i) = obj.constellation_points(symbol + 1);
            end
        end
        
        function [llr, p1] = compute_llr(obj, y, n_0)
            p0 = zeros(length(y) * obj.n_bits, 1);
            p1 = zeros(length(y) * obj.n_bits, 1);
            for y_index = 1 : length(y)
                for sym_index = 1  : 2^obj.n_bits
                    if length(n_0) == 1
                        p_sym = exp(-abs(y(y_index) - obj.constellation_points(sym_index))^2/2/n_0);
                    else
                        p_sym = exp(-abs(y(y_index) - obj.constellation_points(sym_index))^2/2/n_0(y_index));
                    end
                    for m_index = 1 : obj.n_bits
                        if obj.bit_sym_map(sym_index, m_index) == 0
                            p0((y_index-1) * obj.n_bits +  m_index) = p0((y_index-1) * obj.n_bits +  m_index) + p_sym;
                        else
                            p1((y_index-1) * obj.n_bits +  m_index) = p1((y_index-1) * obj.n_bits +  m_index) + p_sym;
                        end
                    end
                end
            end
            llr = log(p0./p1);
            p1 = p1./(p0 + p1);
        end
    end
end

