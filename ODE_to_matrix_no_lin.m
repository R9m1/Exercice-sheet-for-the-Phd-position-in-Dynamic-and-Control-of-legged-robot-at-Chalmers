function [M,C] = ODE_to_matrix_no_lin(eq, ddq)
   
   for i=1:length(eq)
        eqi = eq(i);
        for j=1:length(eq)
            M(i,j) = simplify(diff(eqi, ddq(j)));
        end
       
   end
   disp(M*ddq)
   C = simplify(eq - M*ddq);

end