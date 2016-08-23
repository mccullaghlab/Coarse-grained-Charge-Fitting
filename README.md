# Code to perform coarse-grained charge-fitting


This code repository contains two versions of the coarse-grained charge fitting code.  These are:

## [analytic_position_minimzation](https://github.com/mccullaghlab/Coarse-grained-Charge-Fitting/tree/master/analytic_position_minimization)
This code performs the analytic coarse-grained charge-fitting and site position minimization using analytic derivatives and gradient descent minimization.  

## [grid_anlytic_charge_fit](https://github.com/mccullaghlab/Coarse-grained-Charge-Fitting/tree/master/grid_analytic_charge_fit)
This code performs three versions of coarse-grained charge-fitting.  Namely grid, grid with Lagrange multiplier to maintain charge conservation and the analytic charge optimization.  This code does not perform CG site position minimization.  

Please cite:
P. McCullagh, P. T. Lake, M. McCullagh. **Deriving Coarse-grained Charges from All-atom Systems: An Analytic Solution** *J. Chem. Theor. Comp.*, 2016.
http://pubs.acs.org/doi/abs/10.1021/acs.jctc.6b00507

