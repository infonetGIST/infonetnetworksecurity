function [H] = wb_gn_LDGC_Hs_new(H_row,H_column,H_j,H_k)
% clc; clear;
% H_row = 100; H_column = 300; H_j =10; H_k=5;

H_column_first_part = H_column - H_row;

for find_H = 1:1000000000
    clear candidate;
    
%     # of non-zero elements
    for ind_candidate=1:H_k
        candidate(ind_candidate,:) = [1:H_row];
    end
    candidate = reshape(candidate, 1, H_k*H_row);
    
    Q1 = zeros(H_j, H_column_first_part);
    for Q1_ind_i=1:H_j
        for Q1_ind_j=1:H_column_first_part
            for try_ind=1:100
%                 ind_sel = randi(1,1,[1,length(candidate)]);
                ind_sel = randi(length(candidate));
                
                if (length(find(Q1(:,Q1_ind_j)==candidate(ind_sel)))==0)
                    Q1(Q1_ind_i,Q1_ind_j)=candidate(ind_sel);
                    candidate(ind_sel)=[];
                    break;
                end
            end
        end
    end
    
    H = zeros(H_row,H_column_first_part);
    if(length(find(Q1==0))==0)
        for t = 1:H_column_first_part
            H(Q1(:,t),t)=1;
        end
    end

    if ( sum( abs(sum(H) - H_j*ones(1,H_column_first_part)) )==0 )
        if ( sum( abs(sum(H') - H_k*ones(1,H_row)) )==0 )
%             if (rank(H) >= 0.8*H_row)
            if (rank(H) >= 0.8*min(size(H)))
                break;
            end
        end
    end
    find_H
end

if (find_H < 1000000000)
    disp('random H is successfully generated !');
end

H = [H eye(H_row)];
    