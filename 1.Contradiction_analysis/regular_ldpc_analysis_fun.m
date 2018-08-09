% ***************************************************************** 
% COPYRIGHT (c) 2018 Heung-No Lee, and Woong-Bi Lee. 
% E-mail: heungno@gist.ac.kr, woongbi.lee@gmail.com
% Affiliation: INFONET Laboratory, Gwangju Institute of Science and
% Technology (GIST), Republic of Korea
% homepage: http://infonet.gist.ac.kr
% *****************************************************************  
% filename: regular_ldpc_analysis_fun2.m
% this script generates probabilities of even number of errors and odd
% number of errors at (k-1) source nodes at last stage of decoding
% *****************************************************************
% --------------------- Input --------------------- 
% p: error probability due to channel
% p1: error probability due to channel and pollution attack
% rc: compromised rate 
% dv: number of edges at variable node
% dc: number of edges at check node
% mode: #1 Unanimity Rule
%       #2 Majority Rule
% ***************************************************************** 
% --------------------- Output --------------------- 
% pe_out: even number of errors at (k-1) source nodes at last stage
% po_out: odd number of errors at (k-1) source nodes at last stage
% p_out: probability of errors of source nodes at last stage
% *****************************************************************
%%
function[pe_out, po_out, p_out]=regular_ldpc_analysis_fun(p,p1,rc,dv,dc,mode)

iteration=100;
p_out=p;

if rc == 0
    p2 = p;
else
    p2 = (1-rc)*p + rc*p1;
end

for i=1:1:iteration
    
    % pe: even number of errors in (k-1) digits
    % po: odd number of errors in (k-1) digits
    [pe, po]=get_pe_po2(p_out,p2,dc);    

    [pe_m, po_m]=get_message(pe,po,dv,mode);    
    
    % error probability of source node
    % {error occured, but not flipped} or {error not occured, but flipped}
    p_out=p.*pe_m + (1-p).*po_m; 
end

pe_out = (0.5*(1+(1-2.*p_out).^(dc-1)));
po_out = (0.5*(1-(1-2.*p_out).^(dc-1)));
end