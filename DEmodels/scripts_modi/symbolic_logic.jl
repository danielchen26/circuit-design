## Import packages
using Symbolics, SymbolicUtils


## variables
N = 3
@variables B[1:N] C[1:N-1]


##
function multi_bit_logic()
    input = 1
    mb = zeros(4)
    b1 = mb[4]; b2 = mb[3]; b3 = mb[2]; b4 = mb[1];
    for i in 1: 4
        b1 = !Bool(b1)
        @show b1
    end
end

multi_bit_logic()

exp1 = B[1] + C[1]

exp2up =
