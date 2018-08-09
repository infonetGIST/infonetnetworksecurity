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