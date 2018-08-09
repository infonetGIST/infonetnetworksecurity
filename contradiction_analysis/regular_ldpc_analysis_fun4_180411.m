function[Z_comp, Z_usual]=regular_ldpc_analysis_fun4_180411(p,pa,rc,dv,dc,mode)

j = dv; % degree of variable nodes
k = dc; % degree of checknodes

iteration=100;
p_out=p;

p1 = p*ones(1,length(pa)) + pa - 2*p*pa;

if rc == 0
    p2 = p;
else
    p2 = (1-rc)*p + rc*p1;
end

for i=1:1:iteration

    % pe: even number of errors in (k-1) digits
    % po: odd number of errors in (k-1) digits
    [pe, po]=get_pe_po2(p_out,p2,k);
    [pe_m, po_m]=get_message(pe,po,j,mode);        
    % {error occured, but not flipped} or {error not occured, but flipped}
    p_out=p.*(1-pe_m) + (1-p).*po_m ;
    
end

% [pe po]=get_pe_po(p_out,k);

pe_m = (0.5*(1+(1-2.*p_out).^(k-1)));
po_m = (0.5*(1-(1-2.*p_out).^(k-1)));

if rc == 0
    p_comp = p*ones(1,length(p1));
else
    p_comp = p1;
end

Z_comp = pe_m.*(p_comp) + po_m.*(1-(p_comp));
Z_usual = pe_m.*(p*ones(1,length(p1))) + po_m.*(1-(p*ones(1,length(p1))));
end