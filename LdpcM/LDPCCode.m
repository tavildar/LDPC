classdef LDPCCode < handle
    % Implements LDPC functionality
    %  
    
    properties
        H;
        N;
        K;
        M;
        max_col_weight;
        max_row_weight;
        col_weight_vec;
        row_weight_vec;
        col_mat;
        row_mat;
        
        Z; 
    end
    
    methods
        function obj = LDPCCode(block_length, info_length)
            obj.N = block_length;
            obj.K = info_length;
            obj.M = obj.N - obj.K;
        end
        
        function generate_Gallager_LDPC(obj)
            obj.H = zeros(obj.M, obj.N);
            w_r = 6;
            w_c = 3;
            A_0 = zeros(obj.M/w_c, obj.N);
            for i_row = 1 : obj.M/w_c
                A_0(i_row, (i_row - 1)*w_r + 1:i_row*w_r) = 1;
            end
            for i_col = 1 : w_c
                obj.H((i_col-1)*obj.M/w_c + 1 : i_col * obj.M/w_c, :) =  A_0(:,randperm(size(A_0,2)));
            end
        end
        
        function efficient_pcm(obj)
            obj.col_weight_vec = sum(obj.H, 1);
            obj.row_weight_vec = sum(obj.H, 2);
            obj.max_col_weight  = max(obj.col_weight_vec);
            obj.max_row_weight  = max(obj.row_weight_vec);
            obj.col_mat = zeros(obj.N, obj.max_col_weight); 
            obj.row_mat = zeros(obj.M, obj.max_row_weight); 
            for i_row = 1 : obj.M
                index = 1;
                for i_col = 1 : obj.N
                    if obj.H(i_row, i_col)
                        obj.row_mat(i_row, index) = i_col;
                        index = index + 1;
                    end
                end
            end
            
            for i_col = 1 : obj.N
                index = 1;
                for i_row = 1 : obj.M
                    if obj.H(i_row, i_col)
                        obj.col_mat(i_col, index) = i_row;
                        index = index + 1;
                    end
                end
            end
            
        end
        
        function [codeword] = encode_bits(obj, info_bits)
           % Does encoding by back substitution
           % Assumes a very specific structure on the partiy check matrix
           codeword = zeros(obj.N, 1);
           codeword(1:obj.K) = info_bits;
           
           parity = zeros(obj.M, 1);
           for i_row = 1 : obj.M
               for i_col = 1 : obj.max_row_weight
                   if (obj.row_mat(i_row, i_col) > 0) && (obj.row_mat(i_row, i_col) <= obj.K)
                       parity(i_row) = parity(i_row) + codeword(obj.row_mat(i_row, i_col));
                   end
               end
           end
           
           parity = mod(parity, 2);

           for i = 1 : obj.Z
               codeword(obj.K + i) = mod(sum(parity(i:obj.Z:end)), 2);
           end
           
           for i_row = 1 : obj.M
               for i_col = 1 : obj.max_row_weight
                   if (obj.row_mat(i_row, i_col) > obj.K) && (obj.row_mat(i_row, i_col) <= obj.K + obj.Z)
                       parity(i_row) = mod(parity(i_row) + codeword(obj.row_mat(i_row, i_col)), 2);
                   end
               end
           end
           
           for i_col = obj.K + obj.Z + 1 : obj.Z: obj.N
                codeword(i_col: i_col + obj.Z - 1) = parity(i_col -  obj.K - obj.Z : i_col - obj.K - 1);
                parity(i_col -  obj.K : i_col - obj.K + obj.Z - 1) = mod(parity(i_col -  obj.K - obj.Z : i_col - obj.K - 1) + ...
                    parity(i_col -  obj.K : i_col + obj.Z - obj.K - 1), 2);
           end
           
        end
        
        function [decoded_codeword, error_vec] = decode_llr(obj, input_llr_vec, max_iter, min_sum)
            eta = zeros(obj.M, obj.N);
            lasteta = zeros(obj.M, obj.N);
            updated_llr_vec = input_llr_vec;
            error_vec = zeros(max_iter, 1);
            for iter = 1 : max_iter

                for i_m = 1 : obj.M
                    for i_n1 = 1 : obj.row_weight_vec(i_m) 
                        n1 = obj.row_mat(i_m, i_n1);
                        if min_sum
                            pr = 100;
                        else
                            pr = 1;
                        end
                        for  i_n2 = 1 : obj.row_weight_vec(i_m)
                            if i_n1 == i_n2
                                continue;
                            end
                            n2 = obj.row_mat(i_m, i_n2);
                            l1 = (updated_llr_vec(n2) - lasteta(i_m, n2));
                            l1 = min(l1, 20);
                            l1 = max(l1, -20);
                            if min_sum
                                pr = sign(pr) * sign(l1) * min(abs(l1), abs(pr));
                            else
                                pr = pr * tanh(l1/2);
                            end
                        end
                        if min_sum
                             eta(i_m, n1) = pr; 
                        else
                            eta(i_m, n1) = 2 * atanh(pr);
                        end

                    end
                 end
                 
                 lasteta = eta;
                 
                 for i_n = 1 : obj.N
                     updated_llr_vec(i_n) = input_llr_vec(i_n);
                     for i_m = 1 : obj.col_weight_vec(i_n)
                         m = obj.col_mat(i_n, i_m);
                         updated_llr_vec(i_n) = updated_llr_vec(i_n) + eta(m,i_n);
                     end
                 end
                 
                 decoded_codeword = (updated_llr_vec < 0);
                 if obj.check_codeword(decoded_codeword)
                      return;
                 else
                     error_vec(iter) = 1;
                 end
            end
           
        end
        
        function [b] = check_codeword(obj, x)
           b = 1;
           for i_check = 1 : obj.M
               c = 0;
               for i_n = 1  : obj.row_weight_vec(i_check)
                   c = c + x(obj.row_mat(i_check, i_n));
               end
               if mod(c, 2) == 1
                   b = 0;
                   break;
               end                
           end
        end
                
        
        function lifted_ldpc(obj, baseH, Z)
            circ_mat_array = {};
            for i = 0 : Z - 1
                circ_mat = zeros(Z,Z);
                for j = 0 : Z - 1
                    circ_mat(j + 1, mod(j + i, Z) + 1) = 1;
                end
                circ_mat_array{i + 1} = circ_mat;
            end
            obj.N = Z *  size(baseH, 2);
            obj.M = Z *  size(baseH, 1);
            obj.K = obj.N - obj.M;
            obj.H = zeros(obj.M, obj.N);
            
            for i_row = 1 : size(baseH, 1)
                for i_col = 1 : size(baseH, 2)
                    if baseH(i_row, i_col) ~= -1
                        obj.H( (i_row-1)*Z + 1 : i_row * Z, (i_col-1)*Z + 1 : i_col*Z) = circ_mat_array{baseH(i_row, i_col) + 1};
                    end
                    
                end
            end
            
            obj.efficient_pcm();
            
        end
        
        
        function load_wifi_ldpc(obj, block_length, rate)

            H_1296_1_2 = [40 -1 -1 -1 22 -1 49 23 43 -1 -1 -1 1 0 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1;
                50 1 -1 -1 48 35 -1 -1 13 -1 30 -1 -1 0 0 -1 -1 -1 -1 -1 -1 -1 -1 -1;
                39 50 -1 -1 4 -1 2 -1 -1 -1 -1 49 -1 -1 0 0 -1 -1 -1 -1 -1 -1 -1 -1;
                33 -1 -1 38 37 -1 -1 4 1 -1 -1 -1 -1 -1 -1 0 0 -1 -1 -1 -1 -1 -1 -1;
                45 -1 -1 -1 0 22 -1 -1 20 42 -1 -1 -1 -1 -1 -1 0 0 -1 -1 -1 -1 -1 -1;
                51 -1 -1 48 35 -1 -1 -1 44 -1 18 -1 -1 -1 -1 -1 -1 0 0 -1 -1 -1 -1 -1;
                47 11 -1 -1 -1 17 -1 -1 51 -1 -1 -1 0 -1 -1 -1 -1 -1 0 0 -1 -1 -1 -1;
                5 -1 25 -1 6 -1 45 -1 13 40 -1 -1 -1 -1 -1 -1 -1 -1 -1 0 0 -1 -1 -1;
                33 -1 -1 34 24 -1 -1 -1 23 -1 -1 46 -1 -1 -1 -1 -1 -1 -1 -1 0 0 -1 -1;
                1 -1 27 -1 1 -1 -1 -1 38 -1 44 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 0 0 -1;
                -1 18 -1 -1 23 -1 -1 8 0 35 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 0 0;
                49 -1 17 -1 30 -1 -1 -1 34 -1 -1 19 1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 0];
            
            H_1296_2_3 = [39 31 22 43 -1 40 4 -1 11 -1 -1 50 -1 -1 -1 6 1 0 -1 -1 -1 -1 -1 -1;
                25 52 41 2 6 -1 14 -1 34 -1 -1 -1 24 -1 37 -1 -1 0 0 -1 -1 -1 -1 -1;
                43 31 29 0 21 -1 28 -1 -1 2 -1 -1 7 -1 17 -1 -1 -1 0 0 -1 -1 -1 -1;
                20 33 48 -1 4 13 -1 26 -1 -1 22 -1 -1 46 42 -1 -1 -1 -1 0 0 -1 -1 -1;
                45 7 18 51 12 25 -1 -1 -1 50 -1 -1 5 -1 -1 -1 0 -1 -1 -1 0 0 -1 -1;
                35 40 32 16 5 -1 -1 18 -1 -1 43 51 -1 32 -1 -1 -1 -1 -1 -1 -1 0 0 -1;
                9 24 13 22 28 -1 -1 37 -1 -1 25 -1 -1 52 -1 13 -1 -1 -1 -1 -1 -1 0 0;
                32 22 4 21 16 -1 -1 -1 27 28 -1 38 -1 -1 -1 8 1 -1 -1 -1 -1 -1 -1 0];
            
            H_1296_3_4 = [39 40 51 41 3 29 8 36 -1 14 -1 6 -1 33 -1 11 -1 4 1 0 -1 -1 -1 -1;
                48 21 47 9 48 35 51 -1 38 -1 28 -1 34 -1 50 -1 50 -1 -1 0 0 -1 -1 -1;
                30 39 28 42 50 39 5 17 -1 6 -1 18 -1 20 -1 15 -1 40 -1 -1 0 0 -1 -1;
                29 0 1 43 36 30 47 -1 49 -1 47 -1 3 -1 35 -1 34 -1 0 -1 -1 0 0 -1;
                1 32 11 23 10 44 12 7 -1 48 -1 4 -1 9 -1 17 -1 16 -1 -1 -1 -1 0 0;
                13 7 15 47 23 16 47 -1 43 -1 29 -1 52 -1 2 -1 53 -1 1 -1 -1 -1 -1 0];
            
            H_1296_5_6 = [48 29 37 52 2 16 6 14 53 31 34 5 18 42 53 31 45 -1 46 52 1 0 -1 -1;
                17 4 30 7 43 11 24 6 14 21 6 39 17 40 47 7 15 41 19 -1 -1 0 0 -1;
                7 2 51 31 46 23 16 11 53 40 10 7 46 53 33 35 -1 25 35 38 0 -1 0 0;
                19 48 41 1 10 7 36 47 5 29 52 52 31 10 26 6 3 2 -1 51 1 -1 -1 0 ];
            
            H_1944_1_2 = [57 -1 -1 -1 50 -1 11 -1 50 -1 79 -1 1 0 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1;
                3 -1 28 -1 0 -1 -1 -1 55 7 -1 -1 -1 0 0 -1 -1 -1 -1 -1 -1 -1 -1 -1;
                30 -1 -1 -1 24 37 -1 -1 56 14 -1 -1 -1 -1 0 0 -1 -1 -1 -1 -1 -1 -1 -1;
                62 53 -1 -1 53 -1 -1 3 35 -1 -1 -1 -1 -1 -1 0 0 -1 -1 -1 -1 -1 -1 -1;
                40 -1 -1 20 66 -1 -1 22 28 -1 -1 -1 -1 -1 -1 -1 0 0 -1 -1 -1 -1 -1 -1;
                0 -1 -1 -1 8 -1 42 -1 50 -1 -1 8 -1 -1 -1 -1 -1 0 0 -1 -1 -1 -1 -1;
                69 79 79 -1 -1 -1 56 -1 52 -1 -1 -1 0 -1 -1 -1 -1 -1 0 0 -1 -1 -1 -1;
                65 -1 -1 -1 38 57 -1 -1 72 -1 27 -1 -1 -1 -1 -1 -1 -1 -1 0 0 -1 -1 -1;
                64 -1 -1 -1 14 52 -1 -1 30 -1 -1 32 -1 -1 -1 -1 -1 -1 -1 -1 0 0 -1 -1;
                -1 45 -1 70 0 -1 -1 -1 77 9 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 0 0 -1;
                2 56 -1 57 35 -1 -1 -1 -1 -1 12 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 0 0;
                24 -1 61 -1 60 -1 -1 27 51 -1 -1 16 1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 0];
            
            H_1944_2_3 = [61 75 4 63 56 -1 -1 -1 -1 -1 -1 8 -1 2 17 25 1 0 -1 -1 -1 -1 -1 -1;
                56 74 77 20 -1 -1 -1 64 24 4 67 -1 7 -1 -1 -1 -1 0 0 -1 -1 -1 -1 -1;
                28 21 68 10 7 14 65 -1 -1 -1 23 -1 -1 -1 75 -1 -1 -1 0 0 -1 -1 -1 -1;
                48 38 43 78 76 -1 -1 -1 -1 5 36 -1 15 72 -1 -1 -1 -1 -1 0 0 -1 -1 -1;
                40 2 53 25 -1 52 62 -1 20 -1 -1 44 -1 -1 -1 -1 0 -1 -1 -1 0 0 -1 -1;
                69 23 64 10 22 -1 21 -1 -1 -1 -1 -1 68 23 29 -1 -1 -1 -1 -1 -1 0 0 -1;
                12 0 68 20 55 61 -1 40 -1 -1 -1 52 -1 -1 -1 44 -1 -1 -1 -1 -1 -1 0 0;
                58 8 34 64 78 -1 -1 11 78 24 -1 -1 -1 -1 -1 58 1 -1 -1 -1 -1 -1 -1 0];
            
            
            H_1944_3_4 = [48 29 28 39 9 61 -1 -1 -1 63 45 80 -1 -1 -1 37 32 22 1 0 -1 -1 -1 -1;
                4 49 42 48 11 30 -1 -1 -1 49 17 41 37 15 -1 54 -1 -1 -1 0 0 -1 -1 -1;
                35 76 78 51 37 35 21 -1 17 64 -1 -1 -1 59 7 -1 -1 32 -1 -1 0 0 -1 -1;
                9 65 44 9 54 56 73 34 42 -1 -1 -1 35 -1 -1 -1 46 39 0 -1 -1 0 0 -1;
                3 62 7 80 68 26 -1 80 55 -1 36 -1 26 -1 9 -1 72 -1 -1 -1 -1 -1 0 0;
                26 75 33 21 69 59 3 38 -1 -1 -1 35 -1 62 36 26 -1 -1 1 -1 -1 -1 -1 0];
            
            
            H_1944_5_6 = ...
              [ 13 48 80 66 4 74 7 30 76 52 37 60 -1 49 73 31 74 73 23 -1 1 0 -1 -1;
                69 63 74 56 64 77 57 65 6 16 51 -1 64 -1 68 9 48 62 54 27 -1 0 0 -1;
                51 15 0 80 24 25 42 54 44 71 71 9 67 35 -1 58 -1 29 -1 53 0 -1 0 0;
                16 29 36 41 44 56 59 37 50 24 -1 65 4 65 52 -1 4 -1 73 52 1 -1 -1 0];
            
            
            H_648_1_2 = [0 -1 -1 -1 0 0 -1 -1 0 -1 -1 0 1 0 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1;
                22 0 -1 -1 17 -1 0 0 12 -1 -1 -1 -1 0 0 -1 -1 -1 -1 -1 -1 -1 -1 -1;
                6 -1 0 -1 10 -1 -1 -1 24 -1 0 -1 -1 -1 0 0 -1 -1 -1 -1 -1 -1 -1 -1;
                2 -1 -1 0 20 -1 -1 -1 25 0 -1 -1 -1 -1 -1 0 0 -1 -1 -1 -1 -1 -1 -1;
                23 -1 -1 -1 3 -1 -1 -1 0 -1 9 11 -1 -1 -1 -1 0 0 -1 -1 -1 -1 -1 -1;
                24 -1 23 1 17 -1 3 -1 10 -1 -1 -1 -1 -1 -1 -1 -1 0 0 -1 -1 -1 -1 -1;
                25 -1 -1 -1 8 -1 -1 -1 7 18 -1 -1 0 -1 -1 -1 -1 -1 0 0 -1 -1 -1 -1;
                13 24 -1 -1 0 -1 8 -1 6 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 0 0 -1 -1 -1;
                7 20 -1 16 22 10 -1 -1 23 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 0 0 -1 -1;
                11 -1 -1 -1 19 -1 -1 -1 13 -1 3 17 -1 -1 -1 -1 -1 -1 -1 -1 -1 0 0 -1;
                25 -1 8 -1 23 18 -1 14 9 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 0 0;
                3 -1 -1 -1 16 -1 -1 2 25 5 -1 -1 1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 0];
            
            H_648_2_3 = [25 26 14 -1 20 -1 2 -1 4 -1 -1 8 -1 16 -1 18 1 0 -1 -1 -1 -1 -1 -1;
                10 9 15 11 -1 0 -1 1 -1 -1 18 -1 8 -1 10 -1 -1 0 0 -1 -1 -1 -1 -1;
                16 2 20 26 21 -1 6 -1 1 26 -1 7 -1 -1 -1 -1 -1 -1 0 0 -1 -1 -1 -1;
                10 13 5 0 -1 3 -1 7 -1 -1 26 -1 -1 13 -1 16 -1 -1 -1 0 0 -1 -1 -1;
                23 14 24 -1 12 -1 19 -1 17 -1 -1 -1 20 -1 21 -1 0 -1 -1 -1 0 0 -1 -1;
                6 22 9 20 -1 25 -1 17 -1 8 -1 14 -1 18 -1 -1 -1 -1 -1 -1 -1 0 0 -1;
                14 23 21 11 20 -1 24 -1 18 -1 19 -1 -1 -1 -1 22 -1 -1 -1 -1 -1 -1 0 0;
                17 11 11 20 -1 21 -1 26 -1 3 -1 -1 18 -1 26 -1 1 -1 -1 -1 -1 -1 -1 0];
            
            H_648_3_4 = [16 17 22 24 9 3 14 -1 4 2 7 -1 26 -1 2 -1 21 -1 1 0 -1 -1 -1 -1;
                25 12 12 3 3 26 6 21 -1 15 22 -1 15 -1 4 -1 -1 16 -1 0 0 -1 -1 -1;
                25 18 26 16 22 23 9 -1 0 -1 4 -1 4 -1 8 23 11 -1 -1 -1 0 0 -1 -1;
                9 7 0 1 17 -1 -1 7 3 -1 3 23 -1 16 -1 -1 21 -1 0 -1 -1 0 0 -1;
                24 5 26 7 1 -1 -1 15 24 15 -1 8 -1 13 -1 13 -1 11 -1 -1 -1 -1 0 0;
                2 2 19 14 24 1 15 19 -1 21 -1 2 -1 24 -1 3 -1 2 1 -1 -1 -1 -1 0];
            
            H_648_5_6 = [17 13 8 21 9 3 18 12 10 0 4 15 19 2 5 10 26 19 13 13 1 0 -1 -1;
                3 12 11 14 11 25 5 18 0 9 2 26 26 10 24 7 14 20 4 2 -1 0 0 -1;
                22 16 4 3 10 21 12 5 21 14 19 5 -1 8 5 18 11 5 5 15 0 -1 0 0;
                7 7 14 14 4 16 16 24 24 10 1 7 15 6 10 26 8 18 21 14 1 -1 -1 0];
            
            if block_length == 648
                obj.Z = 27;
            elseif block_length == 1296
                obj.Z = 54;
            elseif block_length == 1944
                obj.Z = 81;
            else
                disp('Not supported block length value');
                return;
            end
            
            H_var = ['H_', num2str(block_length), '_'];
            if rate == 1/2
                H_var = [H_var, '1_2'];
            elseif rate == 2/3
                H_var = [H_var, '2_3'];
            elseif rate == 3/4
                H_var = [H_var, '3_4'];
            elseif rate == 5/6
                H_var = [H_var, '5_6'];
            else
                disp('Not supported rate value');
                return;
            end
            
            baseH = eval(H_var);
                
            obj.lifted_ldpc(baseH, obj.Z);

        end
        
    end
    
end

