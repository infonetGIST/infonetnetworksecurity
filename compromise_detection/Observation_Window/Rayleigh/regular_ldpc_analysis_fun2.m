function[Z_comp Z_usual]=regular_ldpc_analysis_fun2(p_ch,p_r,Rc)

% j=3;k=6;
j=5;k=6;

iteration=100;

p_out = p_ch;
p2 = (1-Rc)*p_ch + Rc*p_r;

for i=1:1:iteration
    
    % pe: even number of errors in (k-1) digits
    % po: odd number of errors in (k-1) digits
    [pe po]=get_pe_po2(p_out,p2,k);
    [pe_m po_m]=get_message(pe,po,j);
    p_out=p_ch.*(1 - pe_m) + (1 - p_ch).*po_m;
    
%     p_out=p.*pe_m + (1-p).*po_m;    
%     p2=p_out;
end

% [pe po]=get_pe_po(p_out,k);

pe_m = (0.5*(1+(1-2.*p_out).^(k-1)));
po_m = (0.5*(1-(1-2.*p_out).^(k-1)));

Z_comp = pe_m.*(p_r) + po_m.*(1 - p_r);
Z_usual = pe_m.*(p_ch) + po_m.*(1- p_ch);
end