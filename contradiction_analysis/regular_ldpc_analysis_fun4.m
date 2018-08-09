function[pe_m, po_m, p_out]=regular_ldpc_analysis_fun4(p,p1,rc,dv,dc,mode)


% j=3;k=6;
j = dv;
k = dc;
% j = 6; k = 5;

iteration=100;
p_out=p;
% p_out=0.3;
% p1=p+p_a-2*p.*p_a;

if rc == 0
    p2 = p;
else
    p2 = (1-rc)*p + rc*p1;
end
% rc; compromise rate
% temp = randsrc(1, iteration, [0 1;1-rc rc]);

for i=1:1:iteration

%     if temp(i)==0
%         p2=p;
%     else
%         p2=p1;
%     end
    
    % pe: even number of errors in (k-1) digits
    % po: odd number of errors in (k-1) digits
    [pe, po]=get_pe_po2(p_out,p2,k);    

    [pe_m, po_m]=get_message(pe,po,j,mode);    
    
    % {error occured, but not flipped} or {error not occured, but flipped}
%     p_out=p.*(1-pe_m) + (1-p).*po_m;
    p_out=p.*pe_m + (1-p).*po_m;
end

% [pe po]=get_pe_po(p_out,k);

pe_m = (0.5*(1+(1-2.*p_out).^(k-1)));
po_m = (0.5*(1-(1-2.*p_out).^(k-1)));




end