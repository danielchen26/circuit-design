# circuit design

Folder descriptions:
1. Circuit Scaling up(bottom-up)/BESS_database : The NAND optimal network in html and csv format
2. Circuit Scaling up(bottom-up) : The BESS related NOR gates optimal boolean network generation.
3. Circuit Partition(top-down) : Related to Graph Partition.
4. Differential Equation models : On the biological sequential counter design.
5. NOR_database : The generated NOR gates optimal design for circuit scaling up problem.
6. Plots: Some statistical plots related to the 4 inputs- 1 output circuit scaling up database.
7. Plots/Statistical results : statistical data for sum and cumsum for different Max gates allowed in the cell.

# Two General Problem
## 1. Differential Equation of Xnor-SR model

- Differential Equation models/scripts/gate_para.py: Generating two libraries of gate parameters for the input-output response functions of each promoter 
- Differential Equation models/scripts : Julia and python script for DE simulations
- Differential Equation models/matlab code : Matlab code for DE simulations
- Differential Equation models/param_db : The gates parameters library need for response functions.
- Differential Equation models/xnor-sr-explicit design : Verilog file input for Cello2 (Shuyi's version currently work, but not Tim's)
## 2. Biological Circuits optimization 
### Scaling up approach
- Circuit Scaling up(bottom-up)/logic_fun.py : functions collections module
- Circuit Scaling up(bottom-up)/logic_BESS.py (**main generating functions**): Generating 4-1 boolean functions with optimal network of NOR gates from BESS NAND gates library.
- Circuit Scaling up(bottom-up)/NOR_database/NOR_df_old.csv: The generated 4input 1ouput optimal NOR gates network database (with input variable a,b,c,d)

### Top down approach (related to graph partition)
- not yet done.