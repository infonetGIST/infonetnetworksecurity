function [LR_p] = wb_LDPC_Decoder_UR_hard2(H,Q1,Q2,Num_iteration,LR_f)

n = size(H,2); % codeword length
j = max(sum(H)); % Max. No. of 1's in a column
k = max(sum(H')); % Max. No. of 1's in a row
L= size(H,1);

LR_r = zeros(n,L);
LR_q = zeros(n,L);

% Iteration:
for iteration = 1:Num_iteration
    if iteration==1
        for t=1:n            
            for m=1:j
                if (Q1(m,t) ~=0) % some element of Q1 matrix is 0 due to lack of link
                    LR_q(t,Q1(m,t)) = LR_f(t);
                end
            end %m
        end
    else
        for t=1:n
            sum_LR_r = 0;
            for m=1:j
                if (Q1(m,t) ~=0) % some element of Q1 matrix is 0 due to lack of link
                    sum_LR_r = mod(LR_f(t) + LR_r(t, Q1(m,t)),2) + sum_LR_r;
                end
            end %m
            for m=1:j
                if (Q1(m,t) ~=0) % some element of Q1 matrix is 0 due to lack of link
                    temp_LR_r = LR_f(t) + sum_LR_r - mod(LR_r(t,Q1(m,t)),2);
                    if t<= L
                        if temp_LR_r >= k-1
                            LR_q(t,Q1(m,t)) = mod(LR_f(t)+1,2);
                        else
                            LR_q(t,Q1(m,t)) = LR_f(t);
                        end                    
                    else
                        LR_q(t,Q1(m,t))= LR_f(t);
                    end
                end
            end
        end
    end
    
    for l=1:L        
        for m=1:k
            if (Q2(m,l) ~=0) % some element of Q2 matrix is 0 due to lack of link
                temp_sgn = 0;                
                for m_ex=1:k
                    if (m_ex ~=m)
                        if (Q2(m_ex,l) ~=0)
                            temp_sgn = LR_q(Q2(m_ex,l), l) + temp_sgn;
                        end
                    end
                end % end m_ex
                
                if ( mod(temp_sgn,2)==0 )
                    LR_r(Q2(m,l), l) = 0;
                else
                    LR_r(Q2(m,l), l) = 1;
                end
            end % end of if (Q2(m,l) ~=0)
        end % end m
    end % end l    
end % iteration

sum_LR = zeros(1, n);
for t=1:n
    for m=1:j
        if (Q1(m,t) ~=0)
            sum_LR(t) = mod(LR_r(t, Q1(m,t))+LR_f(t),2)  + sum_LR(t);
        end
    end
end
    
LR_p = zeros(1,n);
for t=1:n
    if t <= n - L
        if sum_LR(t) >= k
            LR_p(t) = mod(LR_f(t)+1,2);
        else
            LR_p(t) = LR_f(t);
        end
    else
        if sum_LR(t) ==1
            LR_p(t) = mod(LR_f(t)+1,2);
        else
            LR_p(t)=LR_f(t);
        end
    end
end

% output

% LR_p = LR_f + sum_LR;