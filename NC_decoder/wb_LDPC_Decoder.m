% ***************************************************************** 
% COPYRIGHT (c) 2018 Heung-No Lee, and Woong-Bi Lee. 
% E-mail: heungno@gist.ac.kr, woongbi.lee@gmail.com
% Affiliation: INFONET Laboratory, Gwangju Institute of Science and
% Technology (GIST), Republic of Korea
% homepage: http://infonet.gist.ac.kr
% *****************************************************************  
% filename: wb_LDPC_Decoder.m
% This code is for floating-point based message passing decoding.
% *****************************************************************
% --------------------- Input --------------------- 
% H: parity-check matrix
% Q1: indices of variable nodes at each check node
% Q2: indices of check nodes at each variable node
% Num_iteration: maximum number of iterations
% LR_f: Log-likelihood ratio (LLR) of input signal
% ***************************************************************** 
% --------------------- Output --------------------- 
% LR_p: posterior LLR
% ***************************************************************** 

function [LR_p] = wb_LDPC_Decoder(H,Q1,Q2,Num_iteration,LR_f)

n = size(H,2); % codeword length
j = max(sum(H)); % Max. No. of 1's in a column
k = max(sum(H')); % Max. No. of 1's in a row
L= size(H,1);

LR_r = zeros(n,L);
LR_q = zeros(n,L);

% Iteration:
for iteration = 1:Num_iteration
    for t=1:n
        sum_LR_r = 0;
        for m=1:j
            if (Q1(m,t) ~=0) % some element of Q1 matrix is 0 due to lack of link
                sum_LR_r = LR_r(t, Q1(m,t)) + sum_LR_r;
            end
        end %m
        
        for m=1:j
            if (Q1(m,t) ~=0) % some element of Q1 matrix is 0 due to lack of link
                LR_q(t,Q1(m,t)) = LR_f(t) + sum_LR_r - LR_r(t,Q1(m,t));
            end
        end
    end
    
    for l=1:L
        kIndex = sum(Q2(:,l) ~=0); % Indicate k for LR_r
        for m=1:k
            if (Q2(m,l) ~=0) % some element of Q2 matrix is 0 due to lack of link
                temp_sgn = 1;
                temp_f = 0;
                for m_ex=1:k
                    if (m_ex ~=m)
                        if (Q2(m_ex,l) ~=0) % some element of Q2 matrix is 0 due to lack of link
                            temp_sgn = (2*(LR_q(Q2(m_ex,l), l)>0) - 1) *temp_sgn; % (+): 2*1 -1, (-):2*0 -1
                            temp_f = -log(tanh(abs(LR_q(Q2(m_ex, l), l))/2)) + temp_f;                            
                        end
                    end
                end % end m_ex
                
                if (temp_f ==0) % to prevent error occur in tanh function
%                     temp_f = 10e-100;
                    temp_f = 1e-14;
                end
%                 if temp_f < 1e-14 % to prevent error occur in tanh function
%                     temp_f = 1e-14;
%                 end
                tmp = -log(tanh(temp_f/2));
%                 if tmp > 64
%                     tmp = 64;
%                 end
                LR_r(Q2(m,l), l) = temp_sgn*tmp*(-1)^kIndex;

%                 LR_r(Q2(m,l), l) =temp_sgn*(-log(tanh(temp_f/2)))*(-1)^kIndex; %Indicate k for LR_r
                
            end % end of if (Q2(m,l) ~=0)
        end % end m
    end % end l
end % iteration

% output
sum_LR = zeros(1, n);
for t=1:n
    for m=1:j
        if (Q1(m,t) ~=0) % some element of Q1 matrix is 0 due to lack of link
            sum_LR(t) = LR_r(t, Q1(m,t)) + sum_LR(t);
        end
    end
end

LR_p = LR_f + sum_LR;
            
                