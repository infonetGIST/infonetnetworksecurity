function [H_sys H_new]=wb_gn_sys_H(H_row,H_column,H_j,H_k)
for find_H=1:1000000000
    clear candidate;
    for ind_candidate=1:H_k
        candidate(ind_candidate,:)=[1:H_row];
    end
    candidate = reshape(candidate, 1, H_k*H_row);
    % dimension of candidate: [H_k x H_row]
    
    Q1=zeros(H_j,H_column);
    % dimension of Q1: [H_j, H_column]
    for Q1_ind_i=1:H_j
        for Q1_ind_j=1:H_column
            for try_ind=1:100%times to find Q1
                ind_sel=randint(1,1,[1,length(candidate)]);
                if (length(find(Q1(:,Q1_ind_j)==candidate(ind_sel)))==0)
                    Q1(Q1_ind_i,Q1_ind_j)=candidate(ind_sel);
                    candidate(ind_sel)=[];
                    break;
                end %if
            end %for try
        end% %for  Q1_ind_j
    end %for Q1_ind_i


    H=zeros(H_row, H_column);
    if (length(find(Q1==0))==0) %fin H only if Q1 is valid
       for t=1:H_column
           H(Q1(:,t),t)=1;
       end
    end
    
    %check the No. of 1's in each row &column is rightQ
    %check H is valid    
    if (sum(abs(sum(H)-H_j*ones(1,H_column)))==0)
        if (sum(abs(sum(H')-H_k*ones(1,H_row)))==0)
            if (rank(rref_mod2(H))>=H_row-1)  %the rank of mod2 is differnt than tranditional rank
                    break;
            end
        end
    end
    
    find_H %disp Find_H
end

if (find_H<1000000000)                %check H is generated before stopping criterion
    disp('random H is successfully generated !')
end

% dimension of H: [H_row, H_column]

%find the rref of rev_H matrix
%re-arrange H such that the rank of first H_column by H_column matrix is
%the Identity matrix
rev_H=H(end:-1:1, end:-1:1);%reverse and updown the H matrix 
                             %in order to use rref_mod2 function

rref_H=rref_mod2(rev_H);      %find rref form of rev_H
for index_find=1:H_row    %find the first 1 in the rref_H
    rearrange_ind(index_find)=find(rref_H(index_find,:),1);
end

%re-arrange rref_H s.t. H=[I Pt}
rest_rref_H=rref_H;   %initial setting
rest_rref_H(:,rearrange_ind)=[]; %assing the rest of H matrix
rearranged_H=[rref_H(:,rearrange_ind) rest_rref_H];

%generate corresponding new H
rest_rev_H_new=rev_H;
rest_rev_H_new(:,rearrange_ind)=[];
rev_H_new=[rev_H(:,rearrange_ind) rest_rev_H_new];

%generate systematic H
H_sys=rearranged_H(end:-1:1,end:-1:1);
H_new=rev_H_new(end:-1:1,end:-1:1);