function [Nest,idx]  = selectNormals(q,N1,N2,th)
idx = [];
A = reshape(1:size(q,2),sqrt(size(q,2)),sqrt(size(q,2)));
for i = 1 : size(q,2)
    [ic,icd] = ixneighbors(A,i);
    
    [tf, loc] = ismember(idx,icd);
    loc(loc==0) = [];
    icd(loc) = [];
    
    N1_n = N1(:,icd);
    N2_n = N2(:,icd);
    
    N1_r = repmat(N1(:,i),1,length(icd));
    N2_r = repmat(N2(:,i),1,length(icd));
    
    res_N1 = median(acosd(dot(N1_r,N1_n)));
    res_N2 = median(acosd(dot(N2_r,N2_n)));
    
    if res_N1 <= res_N2
    Nest(:,i) = N1(:,i);
    else
        Nest(:,i) = N2(:,i);
    end
    
    if min(res_N1,res_N2) > th
        idx = [idx,i];
    end
end