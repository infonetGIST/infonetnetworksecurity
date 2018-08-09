%DAS(Distributed Assistant Antenna System) Network performance simulation
% clc;clear all;close all;
% mode==1; all attack probabilities are equal
% mode==2; all attack probabilities are different
function [AveLRContradict, tmp_pos] = wb_Main_DAS_NetwrokCoding(Ns,Nd,pc,pa,EsNodB,num_avg,mode)

%System Parameters
% Ns=100; %# of Senders
% Nd=100; %# of distributed relays
% ICR=5;%# of incoming connection per relay 
% %Enter 1 for no time-domain coding
% n=1;  %# of output coded bits in time-domain
% k=1;  %# of input bits in time-domain
% H_j=1;%# of 1's in a column in the H matrix
% H_k=1;%# of 1's in a row i the H matrix
%%%%# of inner iterations in LDPC decoder
Num_iteration=100; 

AttackPercent=pc; %0.05 %# % of relay nodes are attached
Pinsist=pa;        %Insistance of the attacker to flip the binary message out of all attacked relay nodes

n_compro = Nd*AttackPercent;
n_usual = Nd*(1-AttackPercent);

load GG GG
load G_NC G_NC
load H_Mesh H_Mesh
load Q1 Q1
load Q2 Q2
%%%
% % EEsNO2=0:1:5; %SNR = 5dB
% EEsNO2=EbNo; %SNR = 5dB
% %%%

D_Prob = zeros(1,length(EsNodB));
MD_Prob = zeros(1,length(EsNodB));
FA_Prob = zeros(1,length(EsNodB));

AveLRContradict=[];

NoAttackedRelayNode = ceil(Nd*AttackPercent); % No. of Attacked Relay Nodes among total No. of Relay Nodes
% Attack position
tmp_candi = randperm(Nd);
tmp_pos = tmp_candi(1:NoAttackedRelayNode);
tmp_pos = sort(tmp_pos);

for IEsNO2=1:length(EsNodB)    
    EsNo = 10^(EsNodB/10);
    No = 1/EsNo;
    
    NumOfError = 0;
    
    noAveSamples = num_avg; %# of test samples
    
    for repeat = 1:1
        MessageMatrix = randsrc(noAveSamples, Ns, [0 1]);
        LRContradict = 0;        
        for NumOfTx = 1:noAveSamples
%             NumOfTx        
            
            Codeword = MessageMatrix(NumOfTx,:);
            MeshCodeword=[];            
            for Inc = 1:size(Codeword, 1) % Inc : 1 ~ NumOfTx
                MeshCodewordTemp = mod(Codeword(Inc,:)*G_NC ,2); % Generating Mesh Codeword                
                
                if mode ==1
                    AttackPatch = randsrc(1, NoAttackedRelayNode, [0 1; (1-Pinsist) Pinsist]); % Choose first No. of Attacked Relay Nodes
                else
%                    Pinsist = rand(1,NoAttackedRelayNode);
                   Pinsist = 1/NoAttackedRelayNode:1/NoAttackedRelayNode:1;
                   for i=1:NoAttackedRelayNode
                       AttackPatch(i) = randsrc(1, 1, [0 1; (1-Pinsist(i)) Pinsist(i)]);
                   end                    
                end

                %AttackedMeshCodeword
%                 MeshCodewordTemp(Ns+1:Ns+NoAttackedRelayNode) = mod(MeshCodewordTemp(Ns+1: Ns+NoAttackedRelayNode)+AttackPatch, 2);

                MeshCodewordTemp(Ns + tmp_pos) = mod(MeshCodewordTemp(Ns + tmp_pos) + AttackPatch, 2);
                MeshCodeword = [MeshCodeword; MeshCodewordTemp];
            end
            
            %Channel effect
            noise = sqrt(No/2)*randn(1,Ns+Nd);
            yMatrix = (2*MeshCodeword - 1) + noise;
%             yMatrix = awgn((2*MeshCodeword - 1), EsNo);
            
            %LDPC decodor
            %initial:
            y=reshape(yMatrix', 1, size(yMatrix,1)*size(yMatrix,2));
                        
            LR_f = (4*EsNo).*y;
            LR_fMatrix(NumOfTx,:) = LR_f;
            LR_f(find(LR_f>0)) = 1;
            LR_f(find(LR_f<0)) = 0;
       
%             LR_p = wb_LDPC_Decoder_UR_hard(H_Mesh, Q1,Q2, Num_iteration, LR_f);
            LR_p = wb_LDPC_Decoder_UR_hard2(H_Mesh, Q1,Q2, Num_iteration, LR_f);
%             LR_p = wb_LDPC_Decoder_MR_hard(H_Mesh, Q1,Q2, Num_iteration, LR_f);
            
            %Compare LLR from channel and extrinc LLR
            %  to find which relay nodes are attacked
            LRContradict = mod((LR_p>0)-(LR_f>0), 2) + LRContradict; %extrinc=(LR_p-LR_f)+LRContradict;
        end
        
%         AveLRContradict(IEsNO2, :) = LRContradict/(noAveSamples);
        AveLRContradict(IEsNO2, :) = LRContradict;
%         tempLR = AveLRContradict(IEsNO2, :);
        
%         %         Threshold
%         m = mean(AveLRContradict(Ns+1:Ns+Nd));
%         v = std(AveLRContradict(Ns+1:Ns+Nd));
%         thre = m + v;
%         
%         AttackRelayIndex = find(tempLR(1+Ns:end) > 0.370221508499958);
% %         AttackRelayIndex = find(tempLR(1+Ns:end) > thre)
    end
%     
%     Knon_Attack_Posi=1:1:NoAttackedRelayNode;
%     temp=0;
%     for i=1:1:length(Knon_Attack_Posi)
%         for j=1:1:length(AttackRelayIndex)
%             if Knon_Attack_Posi(i)==AttackRelayIndex(j)
%                 temp = temp +1;
%             end
%         end
%             same=temp;
%             MD = length(Knon_Attack_Posi) - temp;
%             FA = length(AttackRelayIndex) - temp;
%     end
%     D_Prob(IEsNO2) = same/n_compro;
%     MD_Prob(IEsNO2) = MD/n_compro;
%     FA_Prob(IEsNO2) = FA/n_usual;
        
end %end of EEsN02

% save AveLRContradict_5dB AveLRContradict;
% figure();
% plot(1:Nd,AveLRContradict(Ns+1:Ns+Nd));
% grid on;
% axis([0 100 0 1]);