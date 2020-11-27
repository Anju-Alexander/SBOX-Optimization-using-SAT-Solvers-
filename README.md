This project tries to optimize the SBOX using SAT Solvers using CRV/PV method (https://eprint.iacr.org/2016/831.pdf)


Dependencies

The software "cnfclaimtoclaim.py" are part of the paper "Optimizing S-box Implementations for Several Criteria using SAT Solvers" by Ko Stoffelen,
published at FSE 2016, available at (https://ko.stoffelen.nl/papers/fse2016-sboxoptimization.pdf).
The link to the code "cnfclaimtoclaim.py" (https://github.com/Ko-/sboxoptimization)

cnfclaimtoclaim.py
In DIMACS CNF format, all variables are replaced by numbers. cnfclaimtoclaim.py takes care of a preprocessing step to translate a solution found by 
the SAT solver (CNF claim file) back to the original variable names. It writes to stdout.


Tool for converting ANF file to CNF
An implementation of Bard-Courtois-Jefferson for converting a sparse system of low-degree multivariate polynomials over GF(2) 
from ANF to CNF can be found here (http://www.nicolascourtois.com/software/CourtoisBardJefferson_public_distribution.zip)
  
 
Example
A typical run will go like this.  Any SAT solver that reads DIMACS CNF files can be used. Here we use SAT SOLVER cryptominisat5 (it works in linux only)

1.run code SATSolverDES.sage, the ANF file is generated (.EQS file) , press enter after the conversion (below steps) is completed

2.run the Bard-Courtois-Jefferson tool to conver the ANF file to CNF (.CNF)

3.$cryptominisat5 --verb 0 des_8.eqs.cnf > des_8.eqs.cnf.claim

4."path to cnfclaimtoclaim.py" des_8.eqs.cnf.claim > des_8.eqs.claim.txt

5.press enter in step 1.


