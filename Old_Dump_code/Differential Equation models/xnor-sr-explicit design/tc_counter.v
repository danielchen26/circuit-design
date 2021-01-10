//xnor-SR 
module ex_counter (output a, input x, y);

    wire w1, w2, w3, s, r, f, m1, m2;
    nor (w1, x, y); 
    nor (w2, x, w1); 
    nor (w3, y, w1); 
    nor (r, w2, w3);
    not (s, r);
    nor (f, m2, s);
    nor (a, m1, r);
    assign x = f
      
endmodule




