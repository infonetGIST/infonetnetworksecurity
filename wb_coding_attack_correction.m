%DAS(Distributed Assistant Antenna System) Network performance simulation
clc;clear;close;
addpath('C:\Users\wbleelab\Dropbox\matlab_network_security\wblee_2018\H_matrix');
addpath('C:\Users\wbleelab\Dropbox\matlab_network_security\wblee_2018\LDPC_decoder');
%%
%System Parameters
Ns=100; %# of Senders
Nd=100; %# of distributed relays
ICR=5;%# of incoming connection per relay 
%Enter 1 for no time-domain coding
n=1;  %# of output coded bits in time-domain
k=1;  %# of input bits in time-domain
H_j=1;%# of 1's in a column in the H matrix
H_k=1;%# of 1's in a row i the H matrix
Num_iteration=20; %# of inner iterations in LDPC decoder

AttackPercent=0.15; %0.05 %# % of relay nodes are attached
Pinsist=1;        %Insistance of the attacker to flip the binary message out of all attacked relay nodes

load GG GG
load G_NC G_NC
load H_Mesh H_Mesh
load Q1 Q1
load Q2 Q2
%%%
EEsNO2 = 0:1:4;
%%%
BERRecord=[];
AveLRContradict=[];

for IEsNO2=1:length(EEsNO2)
    EsNO2 = EEsNO2(IEsNO2); %unit in dB
    NumOfError = 0;
    
    noAveSamples = 1000 %# of test samples
    
    for repeat = 1:1e50
        MessageMatrix = randsrc(noAveSamples, Ns, [0 1]);
        LRContradict = 0;
        
        %% Attack detection
        for NumOfTx = 1:noAveSamples
            NumOfTx
            EsNO2
            
            Codeword = MessageMatrix(NumOfTx,:);
            MeshCodeword=[];
            
            for Inc = 1:size(Codeword, 1) % Inc : 1 ~ NumOfTx
                MeshCodewordTemp = mod(Codeword(Inc,:)*G_NC ,2); % Generating Mesh Codeword
                NoAttackedRelayNode = ceil(Nd*AttackPercent); % No. of Attacked Relay Nodes among total No. of Relay Nodes
                AttackPatch = randsrc(1, NoAttackedRelayNode, [0 1; (1-Pinsist) Pinsist]); % Choose first No. of Attacked Relay Nodes
                
                %AttackedMeshCodeword
                MeshCodewordTemp(Ns+1:Ns+NoAttackedRelayNode) = mod(MeshCodewordTemp(Ns+1: Ns+NoAttackedRelayNode)+AttackPatch, 2);
                
                MeshCodeword = [MeshCodeword; MeshCodewordTemp];
            end
            
            %Channel effect
            %EsN02=EEsN02(IEsN02); %unit in dB
            yMatrix = awgn((2*MeshCodeword - 1), EsNO2);
            
            %LDPC decodor
            %initial:
            y=reshape(yMatrix', 1, size(yMatrix,1)*size(yMatrix,2));
            
            EsNO2_lin = 10^(EsNO2/10);
            LR_f = (2*EsNO2_lin).*y;            
            LR_fMatrix(NumOfTx,:) = LR_f;
            
            noise = yMatrix-(2*MeshCodeword-1);
           
            %decoding process
            LR_f((LR_f > 64)) = 64; %to prevent Nan numerical error
            LR_f((LR_f <-64)) = -64;
            
            LR_p = wb_LDPC_Decoder(H_Mesh, Q1,Q2, Num_iteration, LR_f);
            
            %Compare LLR from channel and extrinc LLR
            %  to find which relay nodes are attacked
            LRContradict = mod(((LR_p-LR_f)>0)-(LR_f>0), 2) + LRContradict; %extrinc=(LR_p-LR_f)+LRContradict;
        end
        
        AveLRContradict(IEsNO2, :) = LRContradict/(noAveSamples);
        tempLR = AveLRContradict(IEsNO2, :);        
        AttackRelayIndex = find(tempLR(1+Ns:end) > 0.5)        
        %%
        for index = 1:size(LR_fMatrix, 1)
            EsNO2_lin = 10^(EsNO2/10);
            LR_f = LR_fMatrix(index,:);

            posind1 = AttackRelayIndex(find( LR_f(Ns+AttackRelayIndex)> 0));
            negind1 = AttackRelayIndex(find( LR_f(Ns+AttackRelayIndex)< 0));
            
            LR_f(Ns+posind1) = LR_f(Ns+posind1) * (-1);
            LR_f(Ns+negind1) = LR_f(Ns+negind1) * (-1);
                                    
%             LR_f(Ns+posind1) = LR_f(Ns+posind1) - 4*EsNO2_lin;
%             LR_f(Ns+negind1) = LR_f(Ns+negind1) + 4*EsNO2_lin;
            
            %decoding process
            LR_f((LR_f>64))=64; %to provent Nan numerical error
            LR_f(LR_f<-64)=-64; %to provent Nan numerical error

            LR_p=wb_LDPC_Decoder(H_Mesh,Q1,Q2,Num_iteration,LR_f);
            
            %Decision
            LDPCDecoderOut=(LR_p>0);
            
            %compute BER only for message bits
            %LDPCDecoderOutMesh=reshape(LDPCDecoderOut,size(MeshCodeword,2),size(MeshCodeword,1))';
            %NumOfError=sum(sum(mod(LDPCDecoderOutMesh(1:k,1:Ns)+MeshCodeword(1:k,1:Ns),2)))+NumOfError;
            NumOfError=sum(sum(mod(LDPCDecoderOut(1:Ns)+MessageMatrix(index,:),2)))+NumOfError;
        end %end of index = 1:size(LR_fMatrix,1)
        
        %BER breaking criterion
        if (NumOfError>=1000) 
            break;
        end %end of if
    end %end for repeat = 1:1e50
    
    BER=NumOfError/(repeat*noAveSamples*k*Ns);
        
    BERRecord=[BERRecord BER];
    
    save BERRecord BERRecord
    save AveLRContradict AveLRContradict
    
    if (BER<=1e-4)
        break;
    end
end %end of EEsN02

% figure();
% semilogy(EEsNO2,BERRecord,'b-x');
% grid on;
% axis([0 8 1e-4 1e0]);
% xlabel('EbNo');ylabel('BER');
%%
rmpath('C:\Users\wbleelab\Dropbox\matlab_network_security\wblee_2018\H_matrix');
rmpath('C:\Users\wbleelab\Dropbox\matlab_network_security\wblee_2018\LDPC_decoder');