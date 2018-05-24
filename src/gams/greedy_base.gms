$title greedy_base.gms
$ontext
solve base case in greedy heuristic
$offtext

option solprint = off;

# solve with constraints in the objective
$ontext
VoltMagBoundSoft = 0;
PowFlowMagBoundSoft = 0;
CurrFlowMagBoundSoft = 0;
PowBalanceLinSoft = 0;
PowBalanceQuadSoft = 0;
PowBalanceDense = 1;
solve subModelDenseBase using nlp minimizing objVar;
$offtext

# step 1 subproblem solve

# fix to current values
busVoltMagVar.fx(i) = busVoltMagVar.l(i);
busVoltAngVar.fx(i) = busVoltAngVar.l(i);
genPowRealVar.fx(i,j)$GenActive(i,j) = genPowRealVar.l(i,j);
genPowImagVar.fx(i,j)$GenActive(i,j) = genPowImagVar.l(i,j);
swshAdmImagVar.fx(i)$SwshActive(i) = swshAdmImagVar.l(i);

# allow violations
VoltMagBoundSoft = 1;
PowFlowMagBoundSoft = 1;
CurrFlowMagBoundSoft = 1;
PowBalanceLinSoft = 1;
PowBalanceQuadSoft = 0;
PowBalanceDense = 0;

# solve
solve subModelBase using nlp minimizing objVar;

$ontext
MaxInfeas = subModelBase.maxinfes;

# solve with soft constraints if hard constraints are infeasible
if((subModelBase.maxinfes > 1e-4) or (subModelBase.modelstat = 5) or (subModelBase.modelstat = 6),
  VoltMagBoundSoft = 0;
  PowFlowMagBoundSoft = 0;
  CurrFlowMagBoundSoft = 0;
  PowBalanceLinSoft = 0;
  PowBalanceQuadSoft = 1;
  solve subModelBase using nlp minimizing objVar;
  MaxInfeas = subModelBase.maxinfes;
);

# something is wrong if still infeasible
if((subModelBase.maxinfes > 1e-4) or (subModelBase.modelstat = 5) or (subModelBase.modelstat = 6),
  abort 'infeasibility still - base case';
);
$offtext

$include solution_evaluation_base.gms
$include convert_solution_base.gms
$include write_solution_base.gms
