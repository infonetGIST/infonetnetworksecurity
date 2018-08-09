function[pe_m po_m]=get_message(pe,po,j)

% 만장 일치룰; Gallager Algorithm A
pe_m = pe.^(j-1);
po_m = po.^(j-1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Majority Rule; Gallager Algorithm B
% pe_m=0;
% for i = (ceil((j-1)/2)+1):1:(j-1)
%     pe_m = nchoosek(j-1,i).*(pe.^i).*((1-pe).^(j-1-i)) + pe_m;
% end
% po_m=0;
% for i = (ceil((j-1)/2)+1):1:(j-1)
%     po_m = nchoosek(j-1,i).*(po.^i).*((1-po).^(j-1-i)) + po_m;
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end