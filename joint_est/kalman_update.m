function [q, P, Kf] = kalman_update(q, P, Y, tmp_, XI)
Kf = P*tmp_'/(tmp_*P*tmp_'+ XI );
q = q +  Kf*(Y-tmp_*q);
P = P - Kf*tmp_*P;
end
            