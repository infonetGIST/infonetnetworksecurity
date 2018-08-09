% ***************************************************************** 
% COPYRIGHT (c) 2018 Heung-No Lee, and Woong-Bi Lee. 
% E-mail: heungno@gist.ac.kr, woongbi.lee@gmail.com
% Affiliation: INFONET Laboratory, Gwangju Institute of Science and
% Technology (GIST), Republic of Korea
% homepage: http://infonet.gist.ac.kr
% *****************************************************************  
% filename: get_message.m
% This code is for probabilistic analysis in the graph-decoding.
% *****************************************************************
% --------------------- Input --------------------- 
% pe: In a parity-check relation, probability that even number of errors
% occur at (k - 1) nodes 
% po: In a parity-check relation, probability that odd number of errors
% occur at (k - 1) nodes 
% j: number of edges per variable node
% op: #1 Unanimity Rule
%     #2 Majority Rule
% ***************************************************************** 
% --------------------- Output --------------------- 
% pe_m: probability that a message from variable to check node includes
% even number of errors.
% po_m: probability that a message from variable to check node includes
% odd number of errors.
% *****************************************************************

function[pe_m, po_m]=get_message(pe,po,j,op)
pe_m=0;
po_m=0;

% Unanimity Rule; Gallager Algorithm A
if op==1
    pe_m = pe.^(j-1);
    po_m = po.^(j-1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Majority Rule; Gallager Algorithm B
else
%     for i = (ceil((j-1)/2)):1:(j-1)
    for i = floor((j-1)/2)-1:1:(j-1)    
        pe_m = pe_m + nchoosek(j-1,i).*(pe.^i).*((1-pe).^(j-1-i));
    end
% %     for i = (ceil((j-1)/2)):1:(j-1)
%     for i = (floor((j-1)/2)):1:(j-1)
%         po_m = nchoosek(j-1,i).*(po.^i).*((1-po).^(j-1-i)) + po_m;
%     end
% end
po_m = 1 - pe_m;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end