% ***************************************************************** 
% COPYRIGHT (c) 2018 Heung-No Lee, and Woong-Bi Lee. 
% E-mail: heungno@gist.ac.kr, woongbi.lee@gmail.com
% Affiliation: INFONET Laboratory, Gwangju Institute of Science and
% Technology (GIST), Republic of Korea
% homepage: http://infonet.gist.ac.kr
% *****************************************************************  
% filename: wb_Main_gn_Initial_matrix.m
% this script generates a generator matrix and parity check matrix
% *****************************************************************
%%
clc;clear;close all;
%System Parameters
Ns = 50; %# of Senders
Nd = 100; %# of distributed relays
ICR = 5;%# of incoming connection per relay 

n=1;  %# of output coded bits in time-domain
k=1;  %# of input bits in time-domain
H_j=1;%# of 1's in a column in the H matrix
H_k=1;%# of 1's in a row i the H matrix

% function []=wb_Main_gn_Initial_matrix(Ns,Nd,ICR,n,k,H_j,H_k)
    %Generate time-domain parity-check matrix of LDPC codes
    if (n==1 & k==1)
        for IndexUser = 1:Ns
            HH=[];
            GG(:,:,IndexUser) = 1;
        end
    else
        H_row = n-k;
        H_column = n;    
        for IndexUser = 1:Ns
            [H_sys H] = wb_gn_sys_H(H_row,H_column,H_j,H_k);
            G = [eye(H_column - H_row) H_sys(:,1:H_column - H_row)'];
            HH(:,:,IndexUser) = H;
            GG(:,:,IndexUser) = G;
        end
    end
    save HH HH
    save GG GG

    %Generate spatial-domain network code/LDGM code
    %[H]=gn_LDGC_Hs(H_column,H_row,H_j,H_k), H_j/H_k: # of ones in a column/row in Pt matrix

%     H_NC = wb_gn_LDGC_Hs(Nd,(Ns+Nd),Nd*ICR/Ns,ICR);  %return Network Code/LDGC H matrix
    H_NC = wb_gn_LDGC_Hs_new(Nd,(Ns+Nd),Nd*ICR/Ns,ICR);  %return Network Code/LDGC H matrix
%     H_NC = wb_gn_LDGC_Hs_SR(Nd,(Ns+Nd),Nd*ICR/Ns,ICR);  %return Network Code/LDGC H matrix
    
    G_NC = [eye(Ns) H_NC(:,1:Ns)'];%get G=[I P]; from H_sys={Pt I]
    save H_NC H_NC
    save G_NC G_NC

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Create equivelent spatial-time mesh H matrix
    if (n==1&k==1)
        H_Mesh=H_NC;
    else
        %Upper Part
        for BlockRow=1:n
            H_Mesh_Upper([(BlockRow-1)*Nd+1:(BlockRow-1)*Nd+Nd],[(BlockRow-1)*(Ns+Nd)+1:(BlockRow-1)*(Ns+Nd)+Ns+Nd])=H_NC;
        end
        H_Mesh_Bottom=[];

        %Bottom Part
        for BlockRow=1:Ns
            H_Mesh_BottomTemp=zeros((n-k),n*(Ns+Nd));
            for Ipos=1:n
                H_Mesh_BottomTemp(:,(Ipos-1)*(Ns+Nd) + (BlockRow-1) + 1) = HH(:,Ipos,BlockRow);
            end
            H_Mesh_Bottom=[H_Mesh_Bottom;H_Mesh_BottomTemp];
        end

        %Combine Upper and Bottom Parts
        H_Mesh=cat(1,H_Mesh_Upper,H_Mesh_Bottom);
    end %end of "if loop"

    save H_Mesh H_Mesh

    %Generate Q1 and Q2 matrix for H_Mesh
    [index_i index_j]=find(H_Mesh==1); %Generate the Q1 matrix
    for mm=1:max(index_j)
        Q1(1:length(find(index_j==mm)),mm)=index_i(find(index_j==mm));
    end
    disp('Q1 is generated')
    save Q1 Q1
    for mm=1:max(index_i)
        Q2(1:length(find(index_i==mm)),mm)=index_j(find(index_i==mm));
    end
    disp('Q2 is generated')
    save Q2 Q2
% end