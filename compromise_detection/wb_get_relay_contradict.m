%DAS(Distributed Assistant Antenna System) Network performance simulation
% clc;clear all;close all;
% mode==1; all attack probabilities are equal
% mode==2; all attack probabilities are different
function [relay_contradict, tmp_pos] = wb_get_relay_contradict(Ns,Nr,rc,pa,EsNodB,num_avg)
    EsNo = 10.^(EsNodB./10);
    No = 1/EsNo;
    fprintf('Eb/No = %d [dB]\n',EsNodB);
    fprintf('Noise power = %f \n',No);    
    
    Num_iteration=20; 

    AttackPercent=rc; %0.05 %# % of relay nodes are attached
    Pinsist=pa;        %Insistance of the attacker to flip the binary message out of all attacked relay nodes
    NoAttackedRelayNode = ceil(Nr*AttackPercent);
    tmp_pos = 1:NoAttackedRelayNode;    

    load GG GG
    load G_NC G_NC
    load H_Mesh H_Mesh
    load Q1 Q1
    load Q2 Q2
    
    MessageMatrix = randsrc(num_avg, Ns, [0 1]);
    AveLRContradict=[]; LRContradict = 0;
    for NumOfTx = 1:num_avg
        Codeword_origin = mod(MessageMatrix(NumOfTx,:)*G_NC,2);
        AttackPatch = randsrc(1, NoAttackedRelayNode, [0 1; (1-Pinsist) Pinsist]);
        Codeword = Codeword_origin;
        Codeword(Ns + tmp_pos) = mod(Codeword_origin(Ns + tmp_pos) + AttackPatch, 2);
                    
            %Channel effect       
        noise = sqrt(No/2)*randn(1,Ns+Nr);
        signal =  2*Codeword - 1;
        y = signal + noise;
        
        LR_f = (4*EsNo).*y;        
        
        LR_f(find(LR_f>0)) = 1;
        LR_f(find(LR_f<0)) = 0;       
        LR_p = wb_LDPC_Decoder_UR_hard2(H_Mesh, Q1,Q2, Num_iteration, LR_f);            
            
            %Compare LLR from channel and extrinc LLR to find which relay nodes are attacked
        LRContradict = mod((LR_p>0)-(LR_f>0), 2) + LRContradict; %extrinc=(LR_p-LR_f)+LRContradict;
    end
%         AveLRContradict(IEsNO2, :) = LRContradict/(noAveSamples);
    AveLRContradict= LRContradict;    
    relay_contradict = AveLRContradict(1,Ns+1:end)'; 
end
