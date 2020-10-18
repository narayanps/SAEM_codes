function [q, P] = kalman_predict(q, P, A, V_q )
q   = A*q;
P = A*P*A' + V_q;
end