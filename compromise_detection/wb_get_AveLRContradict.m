%DAS(Distributed Assistant Antenna System) Network performance simulation
% clc;clear all;close all;
% mode==1; all attack probabilities are equal
% mode==2; all attack probabilities are different
function [AveLRContradict] = wb_get_AveLRContradict(rc,pa,EbNo,num_avg,mode)

%System Parameters
Ns=100; %# of Senders
Nd=100; %# of distributed relays
Num_iteration = 20; %# of inner iterations in LDPC decoder

AttackPercent = rc; %0.05 %# % of relay nodes are attached
Pinsist = pa;        %Insistance of the attacker to flip the binary message out of all attacked relay nodes

load GG GG
load G_NC G_NC
load H_Mesh H_Mesh
load Q1 Q1
load Q2 Q2
%%%
% EEsNO2=0:1:5; %SNR = 5dB
EEsNO2=EbNo; %SNR = 5dB

AveLRContradict=[];
for IEsNO2=1:length(EEsNO2)
    EsNo2_dB = EEsNO2(IEsNO2); %unit in dB
%     fprintf('Eb/No = %ddB\n\n',EsNo2_dB);
    EsNo2 = 10^(EsNo2_dB/10);
    No = 2/EsNo2;    
    
    noAveSamples = num_avg; %# of test samples
    
    for repeat = 1:1
        MessageMatrix = randsrc(noAveSamples, Ns, [0 1]);
        LRContradict = 0;
        
        for NumOfTx = 1:noAveSamples            
            Codeword = MessageMatrix(NumOfTx,:);
            
            MeshCodewordTemp = mod(Codeword*G_NC ,2); % Generating Mesh Codeword
            NoAttackedRelayNode = ceil(Nd*AttackPercent); % No. of Attacked Relay Nodes among total No. of Relay Nodes
                
            AttackPatch = 0;
            if mode ==1
                AttackPatch = randsrc(1, NoAttackedRelayNode, [0 1; (1-Pinsist) Pinsist]); % Choose first No. of Attacked Relay Nodes
            else
%                    Pinsist = rand(1,NoAttackedRelayNode);
               Pinsist = 1/NoAttackedRelayNode:1/NoAttackedRelayNode:1;                   
%                    Pinsist(1:NoAttackedRelayNode/3) = 0.3;
%                    Pinsist(NoAttackedRelayNode/3+1:NoAttackedRelayNode*2/3) = 0.5;
%                    Pinsist(NoAttackedRelayNode*2/3+1:NoAttackedRelayNode*3/3) = 1.0;
                   
               for i=1:NoAttackedRelayNode
                   AttackPatch(i) = randsrc(1, 1, [0 1; (1-Pinsist(i)) Pinsist(i)]);
               end
            end
                
            %AttackedMeshCodeword
            MeshCodeword = MeshCodewordTemp;
            MeshCodeword(Ns+1:Ns+NoAttackedRelayNode) = mod(MeshCodewordTemp(Ns+1: Ns+NoAttackedRelayNode)+AttackPatch, 2);            
            
            %Channel effect
            noise = sqrt(No/2)*randn(1,Ns+Nd);
            y = (2*MeshCodeword - 1) + noise;            
            
            %LDPC decodor                        
            LR_f = (2*EsNo2).*y;            
            
%             LR_f((LR_f > 64)) = 64; %to prevent Nan numerical error
%             LR_f((LR_f <-64)) = -64;            
%             LR_p = wb_LDPC_Decoder(H_Mesh, Q1,Q2, Num_iteration, LR_f);
%             LRContradict = mod(((LR_p-LR_f)>0)-(LR_f>0), 2) + LRContradict; %extrinc=(LR_p-LR_f)+LRContradict;
            
            LR_f(find(LR_f>0)) = 1;
            LR_f(find(LR_f<0)) = 0;
%             LR_p = wb_LDPC_Decoder_MR_hard(H_Mesh, Q1,Q2, Num_iteration, LR_f);
            LR_p = wb_LDPC_Decoder_UR_hard(H_Mesh, Q1,Q2, Num_iteration, LR_f);
            LRContradict = mod((LR_p>0)-(LR_f>0), 2) + LRContradict; %extrinc=(LR_p-LR_f)+LRContradict;
        end        
%         AveLRContradict(IEsNO2, :) = LRContradict/(noAveSamples);
        AveLRContradict = LRContradict;
    end        
end %end of EEsN02
% save AveLRContradict_5dB AveLRContradict;
% figure();
% plot(1:Nd,AveLRContradict(Ns+1:Ns+Nd));
% grid on;
% axis([0 100 0 1]);