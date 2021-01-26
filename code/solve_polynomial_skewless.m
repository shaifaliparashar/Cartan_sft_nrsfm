function res = solve_polynomial_skewless(eq)

mpol k1
mpol k2
for i = 1: size(eq,2)
    pol = [k1^14,k1^13*k2,k1^12*k2^2,k1^11*k2^3,k1^10*k2^4,k1^9*k2^5,k1^8*k2^6,k1^7*k2^7,k1^6*k2^8,k1^5*k2^9,k1^4*k2^10,k1^3*k2^11,k1^2*k2^12,k1*k2^13,k2^14,...
           k1^13,k1^12*k2,k1^11*k2^2,k1^10*k2^3,k1^9*k2^4,k1^8*k2^5,k1^7*k2^6,k1^6*k2^7,k1^5*k2^8,k1^4*k2^9,k1^3*k2^10,k1^2*k2^11,k1*k2^12,k2^13,...
           k1^12,k1^11*k2,k1^10*k2^2,k1^9*k2^3,k1^8*k2^4,k1^7*k2^5,k1^6*k2^6,k1^5*k2^7,k1^4*k2^8,k1^3*k2^9,k1^2*k2^10,k1*k2^11,k2^12,...
           k1^11,k1^10*k2,k1^9*k2^2,k1^8*k2^3,k1^7*k2^4,k1^6*k2^5,k1^5*k2^6,k1^4*k2^7,k1^3*k2^8,k1^2*k2^9,k1*k2^10,k2^11,...
           k1^10,k1^9*k2,k1^8*k2^2,k1^7*k2^3,k1^6*k2^4,k1^5*k2^5,k1^4*k2^6,k1^3*k2^7,k1^2*k2^8,k1*k2^9,k2^10,...
           k1^9,k1^8*k2,k1^7*k2^2,k1^6*k2^3,k1^5*k2^4,k1^4*k2^5,k1^3*k2^6,k1^2*k2^7,k1*k2^8,k2^9,...
           k1^8,k1^7*k2,k1^6*k2^2,k1^5*k2^3,k1^4*k2^4,k1^3*k2^5,k1^2*k2^6,k1*k2^7,k2^8,...
           k1^7,k1^6*k2,k1^5*k2^2,k1^4*k2^3,k1^3*k2^4,k1^2*k2^5,k1*k2^6,k2^7,...
           k1^6,k1^5*k2,k1^4*k2^2,k1^3*k2^3,k1^2*k2^4,k1*k2^5,k2^6,...
           k1^5,k1^4*k2,k1^3*k2^2,k1^2*k2^3,k1*k2^4,k2^5,...
           k1^4,k1^3*k2,k1^2*k2^2,k1^1*k2^3,k2^4,...
           k1^3,k1^2*k2,k1^1*k2^2,k2^3,...
           k1^2,k1*k2,k2^2,...
           k1,k2,1]*eq(:,i);
    P = msdp(min(pol));
    [status,obj,M] = msol(P);
    if status ==1
        res(i,:) = double([k1 k2]);
    else
        res(i,:) = [0,0];
    end
    
end
