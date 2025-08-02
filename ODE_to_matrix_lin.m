function [M,K,F] = ODE_to_matrix_lin(eq, ddq, q)
   
   for i=1:length(eq)
        eqi = eq(i);
        for j=1:length(eq)
            M(i,j) = simplify(diff(eqi, ddq(j)));
        end
   end
    C = simplify(eq - M*ddq);
    for i=1:length(eq)
        Ci = C(i);
        for j=1:length(eq)
            K(i,j) = simplify(diff(Ci, q(j)));
        end
   end
   F = simplify(C - K*q);
end