# Constructing bounds for any-dimensional polynomial problems
Code producing the numerical results from the paper "Any-dimensional polynomial optimization via De Finetti theorems" by Eitan Levin and Venkat Chandrasekaran.

## Scripts:
This repo contains the following scripts, which generate the results from the paper:
1. Upper and lower bounds on the mean-field game from Example 6.1: [mean_field_game_example](https://github.com/eitangl/anyDimPolyOpt/blob/main/mean_field_game_example.m).
2. Lower bounds on symmetric functions from Example 6.4: [symmetric_func_example](https://github.com/eitangl/anyDimPolyOpt/blob/main/symmetric_func_example.m). This script requires the Matlab package [Transition Matrices between Symmetric Polynomials](https://www.mathworks.com/matlabcentral/fileexchange/136299-transition-matrices-between-symmetric-polynomials) to be available in Matlab's path.
3. Lower bounds on Ramsey multiplicity from Example 6.7: [exhaustive_search_Ramsey](https://github.com/eitangl/anyDimPolyOpt/blob/main/exhaustive_search_Ramsey.m). This script requires the list of all graphs of sizes 4-10 in graph6 format, which can be downloaded from [this website](https://users.cecs.anu.edu.au/~bdm/data/graphs.html).
4. Upper and lower bounds on the injective homomorphism numbers from Example 6.10: [graph_numbers](https://github.com/eitangl/anyDimPolyOpt/blob/main/graph_numbers.m).

All of the above scripts use [YALMIP](https://yalmip.github.io/download/) with [MOSEK](https://www.mosek.com/downloads/) (which you could replace with an SDP solver of your choice).

In case of issues or questions, please email Eitan (eitanl@caltech.edu)
