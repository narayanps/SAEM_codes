function W = prep_whitener(E, chan_ids, rank_E, pca_flag)
         [V,D] = eig(E(chan_ids,chan_ids));
         D = diag(D);
         [D,I] = sort(D,'descend'); 
         V = V(:,I);
         if pca_flag==1
            D = 1 ./ D; 
            D(rank_E+1:end) = 0;
            W = diag(sqrt(D)) * V';
            W = W(1:rank_E,:);
         else
             D = 1 ./ D;
             W = diag(sqrt(D)) * V';
         end

end

