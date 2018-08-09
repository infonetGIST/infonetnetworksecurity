% ***************************************************************** 
% COPYRIGHT (c) 2018 Heung-No Lee, and Woong-Bi Lee. 
% E-mail: heungno@gist.ac.kr, woongbi.lee@gmail.com
% Affiliation: INFONET Laboratory, Gwangju Institute of Science and
% Technology (GIST), Republic of Korea
% homepage: http://infonet.gist.ac.kr
% *****************************************************************  
% filename: get_pe_po2.m
% This code is for probabilistic analysis in the graph-decoding.
% *****************************************************************
% --------------------- Input --------------------- 
% p: error proabability at source nodes due to channel
% p1: error probability at relay nodes due to channel and pollution attack
% k: number of edges per check node
% ***************************************************************** 
% --------------------- Output --------------------- 
% pe: probability that (k - 1) nodes including one relay node have even
% number of errors
% po: probability that (k - 1) nodes including one relay node have odd
% number of errors
% *****************************************************************

function[pe, po]=get_pe_po2(p,p1,k)

    % even number of errors in (k-1) digits
    % error & odd number of errors + no error & even number of errors
    pe = p1.*(0.5*(1-(1-2.*p).^(k-2))) + (1-p1).*(0.5*(1+(1-2.*p).^(k-2)));
    
    % odd number of errors in (k-1) digits
    % error & even number of errors + no error & odd number of errors
    po = p1.*(0.5*(1+(1-2.*p).^(k-2)))+(1-p1).*(0.5*(1-(1-2.*p).^(k-2)));

end