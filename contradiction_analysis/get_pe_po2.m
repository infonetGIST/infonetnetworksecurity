function[pe po]=get_pe_po2(p,p1,k)

    % even number of errors in (k-1) digits
    % error & odd number of errors + no error & even number of errors
    pe = p1.*(0.5*(1-(1-2.*p).^(k-2))) + (1-p1).*(0.5*(1+(1-2.*p).^(k-2)));
    
    % odd number of errors in (k-1) digits
    % error & even number of errors + no error & odd number of errors
    po = p1.*(0.5*(1+(1-2.*p).^(k-2)))+(1-p1).*(0.5*(1-(1-2.*p).^(k-2)));

end