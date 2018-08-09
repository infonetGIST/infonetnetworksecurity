%DAS(Distributed Assistant Antenna System) Network performance simulation
clc;clear;close all;
addpath('C:\Users\wbleelab\Dropbox\matlab_network_security\wblee_2018\H_matrix');
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
%%%%# of inner iterations in LDPC decoder
Num_iteration=20; 

AttackPercent=0.15; %0.15 %# % of relay nodes are attached
Pinsist=1;        %Insistance of the attacker to flip the binary message out of all attacked relay nodes

%%
load GG GG
load G_NC G_NC
load H_Mesh H_Mesh
load Q1 Q1
load Q2 Q2

%%%
EEsNO2=0:2:20;
%%%

BERRecord=[];
AveLRContradict=[];

for IEsNO2=1:length(EEsNO2)
    EsNO2 = EEsNO2(IEsNO2); %unit in dB
    NumOfError = 0;
    EsNO2     
    
    for repeat = 1:1e50
        LRContradict = 0;        
            
        Codeword = randsrc(1, Ns, [0 1]);
        MeshCodeword=[];            
        MeshCodeword = mod(Codeword*G_NC ,2); % Generating Mesh Codeword
        
        NoAttackedRelayNode = ceil(Nd*AttackPercent); % No. of Attacked Relay Nodes among total No. of Relay Nodes
        AttackPatch = randsrc(1, NoAttackedRelayNode, [0 1; (1-Pinsist) Pinsist]); % Choose first No. of Attacked Relay Nodes
                
        %AttackedMeshCodeword
        MeshCodeword(Ns+1:Ns+NoAttackedRelayNode) = mod(MeshCodeword(Ns+1: Ns+NoAttackedRelayNode)+AttackPatch, 2);
                
        %Channel effect
        %EsN02=EEsN02(IEsN02); %unit in dB
        yMatrix = awgn((2*MeshCodeword - 1), EsNO2);
            
        %LDPC decodor
        %initial:
        y=reshape(yMatrix', 1, size(yMatrix,1)*size(yMatrix,2));
            
        EsNO2_lin = 10^(EsNO2/10);
        LR_f = (2*EsNO2_lin).*y;            
         
        %decoding process
        LR_f((LR_f > 64)) = 64; %to prevent Nan numerical error
        LR_f((LR_f <-64)) = -64;
            
        LR_p = wb_LDPC_Decoder(H_Mesh, Q1,Q2, Num_iteration, LR_f);
         
        %Decision
        LDPCDecoderOut=(LR_p>0);
            
        %compute BER only for message bits
        %LDPCDecoderOutMesh=reshape(LDPCDecoderOut,size(MeshCodeword,2),size(MeshCodeword,1))';
        %NumOfError=sum(sum(mod(LDPCDecoderOutMesh(1:k,1:Ns)+MeshCodeword(1:k,1:Ns),2)))+NumOfError;
        NumOfError=sum(sum(mod(LDPCDecoderOut(1:Ns)+Codeword(1,:),2)))+NumOfError;            

        
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
% figure();
% semilogy(EEsNO2,BERRecord,'b-x');grid on; axis([0 8 1e-4 1e0]);
% xlabel('EbNo');ylabel('BER');
%%
rmpath('C:\Users\wbleelab\Dropbox\matlab_network_security\wblee_2018\H_matrix');