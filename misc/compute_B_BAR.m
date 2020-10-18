function [B_BAR, N_k] = compute_B_BAR(B, M, Z,N)
N_k1k2=zeros(M);
B_BAR =zeros(M,M);
for k1=1:M
    for k2=1:M
        B_BAR(k1,k2,:)=zeros(1,1);
        for i=1:N
            n_k1=sum(Z{i}(:,k1));
            n_k2=sum(Z{i}(:,k2));
            if(and(n_k1>=1,n_k2>=1))
                if not(k1==k2)
                    N_k1k2(k1,k2)=N_k1k2(k1,k2)+n_k1*n_k2;
                end                
                if k1==k2
                    N_k1k2(k1,k2)=N_k1k2(k1,k2)+n_k1;
                end
                ind1=find(Z{i}(:,k1)==1);
                ind2=find(Z{i}(:,k2)==1);
                for j1=1:n_k1
                    for j2=1:n_k2
                        B_BAR(k1,k2,:)=squish(B_BAR(k1,k2,:))+squish(B{i}(ind1(j1),ind2(j2),:));
                    end
                end
            end
        end
        if N_k1k2(k1,k2) ~= 0
            B_BAR(k1,k2,:)=B_BAR(k1,k2,:)/N_k1k2(k1,k2);
        else
            B_BAR(k1,k2,:)=0;
        end
    end
end

N_k=diag(N_k1k2);