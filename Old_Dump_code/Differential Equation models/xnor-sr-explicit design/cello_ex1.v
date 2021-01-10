//xnor 
module example (output a, input x, y);

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


// SR latch
module A(output x,y,  input in1, in2);
 always@(in1,in2)
  begin
   case({in1,in2})
    2'b00: {x,y} = 1'b00;
    2'b01: {x,y} = 1'b01;
    2'b10: {x,y} = 1'b10;
    2'b11: {x,y} = 1'b00;
   endcase
  end
endmodule
