% ***************************************************************** 
% COPYRIGHT (c) 2018 Heung-No Lee, and Woong-Bi Lee. 
% E-mail: heungno@gist.ac.kr, woongbi.lee@gmail.com
% Affiliation: INFONET Laboratory, Gwangju Institute of Science and
% Technology (GIST), Republic of Korea
% homepage: http://infonet.gist.ac.kr
% *****************************************************************  
% filename: wb_Main_DAS_NetworkCoding.m
% this script generates number of polarity contradictions given observation
% window
% *****************************************************************
% --------------------- Input --------------------- 
% Ns: number of source nodes
% Nd: number of relays
% pc: compromised rate 
% pa: attack probability at relay nodes
% EsNodB: Es/No [dB]
% num_avg: observation window
% mode: #1 all attack probabilities are equal
%       #2 all attack probabilities are different
% ***************************************************************** 
% --------------------- Output --------------------- 
% AveLRContradict: Number of polarity contradictions
% tmp_pos: True attack positions
% *****************************************************************
%%
function [AveLRContradict, tmp_pos] = wb_Main_DAS_NetwrokCoding(Ns,Nd,pc,pa,EsNodB,num_avg,mode)

Num_iteration=100; 

AttackPercent=pc; %0.05 %# % of relay nodes are attached
Pinsist=pa;        %Insistance of the attacker to flip the binary message out of all attacked relay nodes

load GG GG
load G_NC G_NC
load H_Mesh H_Mesh
load Q1 Q1
load Q2 Q2

AveLRContradict=[];

NoAttackedRelayNode = ceil(Nd*AttackPercent); % No. of Attacked Relay Nodes among total No. of Relay Nodes

% Attack position
tmp_candi = randperm(Nd);
tmp_pos = tmp_candi(1:NoAttackedRelayNode);
tmp_pos = sort(tmp_pos);

for IEsNO2=1:length(EsNodB)    
    EsNo = 10^(EsNodB/10);
    No = 1/EsNo;    
    
    noAveSamples = num_avg; %# of test samples
    
    for repeat = 1:1
        MessageMatrix = randsrc(noAveSamples, Ns, [0 1]);
        LRContradict = 0;        
        for NumOfTx = 1:noAveSamples
            
            Codeword = MessageMatrix(NumOfTx,:);
            MeshCodeword=[];            
            for Inc = 1:size(Codeword, 1) % Inc : 1 ~ NumOfTx
                MeshCodewordTemp = mod(Codeword(Inc,:)*G_NC ,2); % Generating Mesh Codeword                
                
                if mode ==1
                    AttackPatch = randsrc(1, NoAttackedRelayNode, [0 1; (1-Pinsist) Pinsist]); % Choose first No. of Attacked Relay Nodes
                else
                   Pinsist = 1/NoAttackedRelayNode:1/NoAttackedRelayNode:1;
                   for i=1:NoAttackedRelayNode
                       AttackPatch(i) = randsrc(1, 1, [0 1; (1-Pinsist(i)) Pinsist(i)]);
                   end                    
                end

                %AttackedMeshCodeword
                MeshCodewordTemp(Ns + tmp_pos) = mod(MeshCodewordTemp(Ns + tmp_pos) + AttackPatch, 2);
                MeshCodeword = [MeshCodeword; MeshCodewordTemp];
            end
            
            %Channel effect
            noise = sqrt(No/2)*randn(1,Ns+Nd);
            yMatrix = (2*MeshCodeword - 1) + noise;            

            y=reshape(yMatrix', 1, size(yMatrix,1)*size(yMatrix,2));
                        
            LR_f = (4*EsNo).*y;
            LR_fMatrix(NumOfTx,:) = LR_f;
            LR_f(find(LR_f>0)) = 1;
            LR_f(find(LR_f<0)) = 0;
       
            LR_p = wb_LDPC_Decoder_UR_hard(H_Mesh, Q1,Q2, Num_iteration, LR_f);
            
            %Compare LLR from channel and extrinc LLR
            %  to find which relay nodes are attacked
            LRContradict = mod((LR_p>0)-(LR_f>0), 2) + LRContradict; %extrinc=(LR_p-LR_f)+LRContradict;
        end        

        AveLRContradict(IEsNO2, :) = LRContradict;
    end        
end %end of EEsN02
