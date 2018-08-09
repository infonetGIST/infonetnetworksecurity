clc;clear;
%%
% addpath('C:\Users\wbleelab\Dropbox\matlab_network_security\wblee_2018\H_matrix');
addpath('C:\Users\wbleelab\Dropbox\matlab_network_security\wblee_2018\LDPC_decoder');
%%
% Ns = 50; Nr = 100; dv = 10; dc = 6; num_avg = 186;
% addpath('C:\Users\wbleelab\Dropbox\matlab_network_security\wblee_2018\H_matrix\Fixed_Nr100\Ns50_icr5_dv10_dc6');
% addpath('/Users/woongbi/Dropbox/matlab_network_security/wblee_2018/H_matrix/Fixed_Nr100/Ns50_icr5_dv10_dc6');
% Ns = 100; Nr = 100; dv = 8; dc = 9; num_avg = 231;
% addpath('C:\Users\wbleelab\Dropbox\matlab_network_security\wblee_2018\H_matrix\Fixed_Nr100\Ns100_icr8_dv8_dc9');
Ns = 200; Nr = 100; dv = 6; dc = 13; num_avg = 317;
addpath('C:\Users\wbleelab\Dropbox\matlab_network_security\wblee_2018\H_matrix\Fixed_Nr100\Ns200_icr12_dv6_dc13');
code_rate = Ns / (Ns+Nr);

% Ns=100; %# of Senders
% Nd=100; %# of distributed relays
% ICR=5;%# of incoming connection per relay 

Num_iteration=20; 

AttackPercent=0.05; %0.05 %# % of relay nodes are attached
Pinsist = 1.0;        %Insistance of the attacker to flip the binary message out of all attacked relay nodes

load GG GG
load G_NC G_NC
load H_Mesh H_Mesh
load Q1 Q1
load Q2 Q2

% EEsNO2=0:1:5;
EsNodB_array = -4:1:3;
% EsNodB_array = 5;

NoAttackedRelayNode = ceil(Nr*AttackPercent); % No. of Attacked Relay Nodes among total No. of Relay Nodes
% tmp_candi = randperm(Nr);
% tmp_pos = tmp_candi(1:NoAttackedRelayNode);
% tmp_pos = sort(tmp_pos);
tmp_pos = 1:NoAttackedRelayNode;

BERRecord=[]; AveLRContradict=[];
for idx = 1:length(EsNodB_array)
    tic
%     EbNodB = EbNodB_array(Eb_idx); %unit in dB
%     EbNo = 10.^(EbNodB./10);
%     EsNo = EbNo * code_rate;
%     EsNodB = 10*log10(EsNo);    
%     No = 1/EsNo;
    
    EsNodB = EsNodB_array(idx); %unit in dB
    EsNo = 10.^(EsNodB./10);
    No = 1/EsNo;
    
    fprintf('Eb/No = %d [dB]\n',EsNodB);
    fprintf('Noise power = %f \n',No);

    NumOfError = 0;    
        
    LRContradict = 0;
    MessageMatrix = randsrc(num_avg, Ns, [0 1]);
%     Compromise Detection Part
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
    AveLRContradict(idx, :) = LRContradict;    
    relay_contradict = AveLRContradict(idx,Ns+1:end)';    
    %% K-means clustering
    ini_point = [0;num_avg];
    opts = statset('Display','off','MaxIter',100);
    [Kidx,C,sumd,D] = kmeans(relay_contradict,2,'Distance','sqeuclidean',...
       'Replicates',1,'Start',ini_point,'Options',opts);
    [a,b] = max(C);
    est_idx1 = find(Kidx==1);
    est_idx2 = find(Kidx==2);    
        
%     AttackRelayIndex = find(idx==b)
    Knon_Attack_Posi = tmp_pos;
    Knon_Attack_Posi = Knon_Attack_Posi';                
    AttackRelayIndex = find(Kidx==b);
    same = length(intersect(AttackRelayIndex,Knon_Attack_Posi));
    MD = length(setdiff(Knon_Attack_Posi,AttackRelayIndex));
    FA = length(setdiff(AttackRelayIndex,Knon_Attack_Posi)); 
    D_Prob(idx,1) = same/(NoAttackedRelayNode)*100;
    MD_Prob(idx,1) = MD/(NoAttackedRelayNode)*100;
    FA_Prob(idx,1) = FA/((Nr - NoAttackedRelayNode))*100;    
    fprintf('Detection Probability = %f \n',D_Prob(idx,1));
    
    for repeat = 1:1e6
        Message = randsrc(1, Ns, [0 1]);
        Codeword_origin = mod(Message*G_NC, 2);        
        AttackPatch = randsrc(1, NoAttackedRelayNode, [0 1; (1-Pinsist) Pinsist]);
        Codeword = Codeword_origin;
        Codeword(Ns + tmp_pos) = mod(Codeword_origin(Ns + tmp_pos) + AttackPatch, 2);
            
        %Channel effect
        noise = randn(1,Ns+Nr)*sqrt(No/2);
        signal = 2*Codeword - 1;
        y = signal + noise;
        
        LR_f = (4*EsNo).*y;        
        
        %Discard LLR from attacked relay nodes
        LR_f(Ns+AttackRelayIndex) = 0.0000001;
            
        %decoding process
        LR_f((LR_f>64))=64; %to provent Nan numerical error
        LR_f(LR_f<-64)=-64; %to provent Nan numerical error

        LR_p=wb_LDPC_Decoder(H_Mesh,Q1,Q2,Num_iteration,LR_f);            
        %Decision
        LDPCDecoderOut=(LR_p>0);
            
        %compute BER only for message bits
        NumOfError=sum(sum(mod(LDPCDecoderOut(1:Ns)+Message,2)))+NumOfError;
                
        %BER breaking criterion
        if (NumOfError>=1000) 
            break;
        end %end of if              
    end %end for repeat = 1:1e50
    
    BER=NumOfError/(repeat*Ns);
    fprintf('BER = %f\n\n',BER);
    
    BERRecord=[BERRecord BER];
    
    save BERRecord BERRecord
    save AveLRContradict AveLRContradict
    
    if (BER<=1e-4)
        break;
    end
    toc
end %end of EEsN02

figure();
semilogy(EsNodB_array,BERRecord,'b-x');
grid on;
% axis([0 8 1e-4 1e0]);
xlabel('EbNo');ylabel('BER');
%%
rmpath('C:\Users\wbleelab\Dropbox\matlab_network_security\wblee_2018\LDPC_decoder');

rmpath('C:\Users\wbleelab\Dropbox\matlab_network_security\wblee_2018\H_matrix\Fixed_Nr100\Ns50_icr5_dv10_dc6');
rmpath('C:\Users\wbleelab\Dropbox\matlab_network_security\wblee_2018\H_matrix\Fixed_Nr100\Ns100_icr8_dv8_dc9');
rmpath('C:\Users\wbleelab\Dropbox\matlab_network_security\wblee_2018\H_matrix\Fixed_Nr100\Ns200_icr12_dv6_dc13');