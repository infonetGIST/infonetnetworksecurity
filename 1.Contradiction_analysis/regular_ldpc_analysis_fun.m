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