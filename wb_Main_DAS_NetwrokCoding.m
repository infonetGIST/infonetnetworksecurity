%DAS(Distributed Assistant Antenna System) Network performance simulation
clc;clear all;close all;

%System Parameters
Ns=100; %# of Senders
Nd=100; %# of distributed relays
ICR=5;%# of incoming connection per relay 
%Enter 1 for no time-domain coding
n=1;  %# of output coded bits in time-domain
k=1;  %# of input bits in time-domain
H_j=1;%# of 1's in a column in the H matrix
H_k=1;%# of 1's in a row i the H matrix
%%%%# of inner iterations in LDPC decoder
Num_iteration=20; 

AttackPercent=0.15; %0.05 %# % of relay nodes are attached
Pinsist=1;        %Insistance of the attacker to flip the binary message out of all attacked relay nodes

%%%
%generate the following matrix
%Main_gn_Initial_matrix(Ns,Nd,ICR,n,k,H_j,H_k);
%load HH HH
%load H_NC H_NC
load GG GG
load G_NC G_NC
load H_Mesh H_Mesh
load Q1 Q1
load Q2 Q2
%%%
EEsNO2=0:1:5;
% EEsNO2=5;
%%%


BERRecord=[];

AveLRContradict=[];

for IEsNO2=1:length(EEsNO2)
    EsNO2 = EEsNO2(IEsNO2); %unit in dB
    fprintf('Eb/No = %d dB\n',EsNO2);
    noAveSamples = 943; %# of test samples
    EsNo2 = 10^(EsNO2/10);
    No = 2/EsNo2;
    
    NumOfError = 0;
    
    p=0.5*erfc(sqrt(2*EsNo2)/2);
        
    theo_pa = 0.0; theo_pc = 0.0;
    [min_Prob_comp min_Prob_usual]=regular_ldpc_analysis_fun2(p,theo_pa,theo_pc);
    theo_pa = 1.0; theo_pc = 1.0;
    [max_Prob_comp max_Prob_usual]=regular_ldpc_analysis_fun2(p,theo_pa,theo_pc);
    
    MessageMatrix = randsrc(noAveSamples, Ns, [0 1]);
    LRContradict = 0;   
    
%     Compromise Detection Part
    for NumOfTx = 1:noAveSamples
        Codeword = MessageMatrix(NumOfTx,:);
        MeshCodeword=[];
            
        for Inc = 1:size(Codeword, 1) % Inc : 1 ~ NumOfTx
            MeshCodewordTemp = mod(Codeword(Inc,:)*G_NC ,2); % Generating Mesh Codeword
            NoAttackedRelayNode = ceil(Nd*AttackPercent); % No. of Attacked Relay Nodes among total No. of Relay Nodes
            
            Pinsist = 1/NoAttackedRelayNode:1/NoAttackedRelayNode:1;
            for i=1:NoAttackedRelayNode
                AttackPatch(i) = randsrc(1, 1, [0 1; (1-Pinsist(i)) Pinsist(i)]);
            end
                
            %AttackedMeshCodeword
            MeshCodewordTemp(Ns+1:Ns+NoAttackedRelayNode) = mod(MeshCodewordTemp(Ns+1: Ns+NoAttackedRelayNode)+AttackPatch, 2);
                
            MeshCodeword = [MeshCodeword; MeshCodewordTemp];
        end
            
            %Channel effect
            noise = sqrt(No/2)*randn(1,Ns+Nd);
            yMatrix = (2*MeshCodeword - 1) + noise;
            
            %LDPC decodor
            %initial:
        y=reshape(yMatrix', 1, size(yMatrix,1)*size(yMatrix,2));
            
        LR_f = (2*EsNo2).*y;            
        LR_fMatrix(NumOfTx,:) = LR_f;       
        LR_f(find(LR_f>0)) = 1;
        LR_f(find(LR_f<0)) = 0;
       
        LR_p = wb_LDPC_Decoder_UR_hard(H_Mesh, Q1,Q2, Num_iteration, LR_f);            
            
            %Compare LLR from channel and extrinc LLR to find which relay nodes are attacked
        LRContradict = mod((LR_p>0)-(LR_f>0), 2) + LRContradict; %extrinc=(LR_p-LR_f)+LRContradict;
    end   
%         AveLRContradict(IEsNO2, :) = LRContradict/(noAveSamples);
    AveLRContradict(IEsNO2, :) = LRContradict;    
    
    
%     Compromise Detection
    n = zeros(1,noAveSamples+1);
    for i=0:noAveSamples
        n(i+1) = histc(AveLRContradict(IEsNO2,101:200),i);
    end
    n = n/100;
    relay_contradict = AveLRContradict(IEsNO2,101:200);
           
    if sum(n(round(max_Prob_usual*noAveSamples):1:noAveSamples+1)) ~= 0
        n(round(max_Prob_usual*noAveSamples)) = sum(n(round(max_Prob_usual*noAveSamples):1:noAveSamples+1));
        n(round(max_Prob_usual*noAveSamples)+1:1:noAveSamples+1) = 0;

    %     n(round(max_Prob_usual*num_avg):1:num_avg+1) = 0;    
    %     n(round(max_Prob_usual*num_avg)) = 0.01;
        disp('1');
    else
        n(round(max_Prob_usual*noAveSamples):1:noAveSamples+1) = 0;
        n(round(max_Prob_usual*noAveSamples)) = 0.01;
        disp('2');
    end
    n = n/sum(n);

    for T = round(min_Prob_usual*noAveSamples):1:round(max_Prob_usual*noAveSamples)
    % for T = 0:1:num_avg
        m1 = 0; m2 = 0;
        for index= 0:1:T
            m1 = n(index+1)*index + m1 ;
        end
        m1 = m1 / sum(n(1:T+1));
        for index= T+1:1:noAveSamples
            m2 = n(index+1)*index + m2 ;
        end
        m2 = m2 / sum(n(T+2:noAveSamples+1));

        p1= sum(n(1:T+1));
        p2= sum(n(T+2:noAveSamples+1));

        variance(T+1) = p1*p2*(m1-m2)^2;

    end

    [maximum Optimal_T] = max(abs(variance));
    
    AttackRelayIndex = find(AveLRContradict(IEsNO2,1+Ns:end) > Optimal_T)
    
    
    for repeat = 1:1e50
        MessageMatrix = randsrc(1, Ns, [0 1]);
        LRContradict = 0;
        
        Codeword = MessageMatrix(1,:);
        MeshCodeword=[];
        
        MeshCodewordTemp = mod(Codeword(Inc,:)*G_NC ,2); % Generating Mesh Codeword
        NoAttackedRelayNode = ceil(Nd*AttackPercent); % No. of Attacked Relay Nodes among total No. of Relay Nodes
        
        Pinsist = 1/NoAttackedRelayNode:1/NoAttackedRelayNode:1;
        for i=1:NoAttackedRelayNode
            AttackPatch(i) = randsrc(1, 1, [0 1; (1-Pinsist(i)) Pinsist(i)]);
        end
        
        %AttackedMeshCodeword
        MeshCodewordTemp(Ns+1:Ns+NoAttackedRelayNode) = mod(MeshCodewordTemp(Ns+1: Ns+NoAttackedRelayNode)+AttackPatch, 2);
        MeshCodeword = [MeshCodeword; MeshCodewordTemp];        
            
        %Channel effect
        %EsN02=EEsN02(IEsN02); %unit in dB
        yMatrix = awgn((2*MeshCodeword - 1), EsNO2);
            
        %LDPC decodor
        %initial:
        y=reshape(yMatrix', 1, size(yMatrix,1)*size(yMatrix,2));
            
        EsNO2_lin = 10^(EsNO2/10);
        LR_f = (2*EsNO2_lin).*y;            
              
        
        %Discard LLR from attacked relay nodes
        LR_f(Ns+AttackRelayIndex) = 0.0000001;
            
        %decoding process
        LR_f((LR_f>64))=64; %to provent Nan numerical error
        LR_f(LR_f<-64)=-64; %to provent Nan numerical error

        LR_p=wb_LDPC_Decoder(H_Mesh,Q1,Q2,Num_iteration,LR_f);
            
        %Decision
        LDPCDecoderOut=(LR_p>0);
            
            %compute BER only for message bits
            %LDPCDecoderOutMesh=reshape(LDPCDecoderOut,size(MeshCodeword,2),size(MeshCodeword,1))';
            %NumOfError=sum(sum(mod(LDPCDecoderOutMesh(1:k,1:Ns)+MeshCodeword(1:k,1:Ns),2)))+NumOfError;
        NumOfError=sum(sum(mod(LDPCDecoderOut(1:Ns)+MessageMatrix(1,:),2)))+NumOfError;
                
        %BER breaking criterion
        if (NumOfError>=1000) 
            break;
        end %end of if
    end %end for repeat = 1:1e50
    
    BER=NumOfError/(repeat*k*Ns);
        
    BERRecord=[BERRecord BER];
    
    save BERRecord BERRecord
    save AveLRContradict AveLRContradict
    
    if (BER<=1e-4)
        break;
    end
end %end of EEsN02

figure();
semilogy(EEsNO2,BERRecord,'b-x');
grid on;
axis([0 8 1e-4 1e0]);
xlabel('EbNo');ylabel('BER');