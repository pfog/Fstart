$title greedy_ctg
$ontext
solve contingencies in greedy method
$offtext

$if not set solve_ctgs $set solve_ctgs 1

# fix some complementarities
genCtgPowRealOverVar.fx(i,j,k)$genCtgActive(i,j,k) = 0;
genCtgPowRealUnderVar.fx(i,j,k)$genCtgActive(i,j,k) = 0;
#busCtgVoltMagOverVar.fx(i,k)$(Bus(i) and Ctg(k)) = 0;
#busCtgVoltMagUnderVar.fx(i,k)$(Bus(i) and Ctg(k)) = 0;
#genCtgPowRealVar.fx(i,j,k)$GenCtgPowRealFixLo(i,j,k) = GenPMin(i,j);
#genCtgPowRealVar.fx(i,j,k)$GenCtgPowRealFixUp(i,j,k) = GenPMax(i,j);
#genCtgPowImagVar.fx(i,j,k)$GenCtgPowImagFixLo(i,j,k) = GenQMin(i,j);
#genCtgPowImagVar.fx(i,j,k)$GenCtgPowImagFixUp(i,j,k) = GenQMax(i,j);

# set a start point
busCtgVoltMagVar.l(i,k)$(Bus(i) and Ctg(k)) = busVoltMagVar.l(i);
busCtgVoltAngVar.l(i,k)$(Bus(i) and Ctg(k)) = busVoltAngVar.l(i);
genCtgPowRealVar.l(i,j,k)$(Gen(i,j) and Ctg(k)) = genPowRealVar.l(i,j)$GenCtgActive(i,j,k);
genCtgPowImagVar.l(i,j,k)$(Gen(i,j) and Ctg(k)) = genPowImagVar.l(i,j)$GenCtgActive(i,j,k);
swshCtgAdmImagVar.l(i,k)$(Swsh(i) and Ctg(k)) = swshAdmImagVar.l(i)$SwshActive(i);
areaCtgPowRealChangeVar.l(i,k)$(Area(i) and Ctg(k)) = 0.0;

$ifthen %solve_ctgs%==1
# loop over contingencies to solve decomposable step 2 subproblem
loop(k0$Ctg(k0),
#loop(k0$(Ctg(k0) and sameas(k0, 'L_0012065512')),
#loop(k0$(Ctg(k0) and sameas(k0, 'L_000490HONEAPATH10-000117DONALDS0C1')),

  CtgActive(k) = no;
  CtgActive(k0) = yes;
  GenPowRealVar.fx(i,j)$GenActive(i,j) = GenPowReal(i,j);
  BusVoltMagVar.fx(i)$Bus(i) = BusVoltMag(i);

  # first solve with dense power balance in objective
$ontext
  VoltMagBoundSoft = 0;
  PowFlowMagBoundSoft = 0;
  CurrFlowMagBoundSoft = 0;
  PowBalanceLinSoft = 0;
  PowBalanceQuadSoft = 0;
  PowBalanceDense = 1;
  #subModelDenseCtg.reslim = 10;
  option solprint = off;
  #option solprint = on;
  solve subModelDenseCtg using nlp minimizing objVar;
  CtgMaxInfeas(k0) = subModelDenseCtg.maxinfes;
  CtgMaxCompViol(k0) = max(
    smax((i,j)$GenCtgParticipating(i,j,k0),
        min(genCtgPowRealOverVar.l(i,j,k0), 
            genCtgPowRealVar.l(i,j,k0) - GenPMin(i,j))),
    smax((i,j)$GenCtgParticipating(i,j,k0),
        min(genCtgPowRealUnderVar.l(i,j,k0), 
            GenPMax(i,j) - genCtgPowRealVar.l(i,j,k0))),
    smax(i$BusCtgVoltMagMaintDom(i,k0),
        min(busCtgVoltMagOverVar.l(i,k0),
	    sum((i1,j1)$GenBusCtgVoltMagMaintDom(i1,j1,i,k0),
	        genCtgPowImagVar.l(i1,j1,k0) - GenQMin(i1,j1)))),
    smax(i$BusCtgVoltMagMaintDom(i,k0),
        min(busCtgVoltMagUnderVar.l(i,k0),
	    sum((i1,j1)$GenBusCtgVoltMagMaintDom(i1,j1,i,k0),
	        GenQMax(i1,j1) - genCtgPowImagVar.l(i1,j1,k0)))));
  CtgMaxSoftViol(k0) = 0;
$offtext

#$ontext
  # first solve with hard constraints
  if(1,
    VoltMagBoundSoft = 0;
    PowFlowMagBoundSoft = 0;
    CurrFlowMagBoundSoft = 0;
    PowBalanceLinSoft = 0;
    PowBalanceQuadSoft = 0;
    subModelCtg.reslim = 10;
    option solprint = off;
    #option solprint = on;
    solve subModelCtg using nlp minimizing objVar;
  );

  # try with soft constraints if first solve fails
  if((subModelCtg.maxinfes > 1e-4) or (subModelCtg.modelstat = 5) or (subModelCtg.modelstat = 6),
  #if(1,
    VoltMagBoundSoft = 0;
    PowFlowMagBoundSoft = 0;
    CurrFlowMagBoundSoft = 0;
    PowBalanceLinSoft = 1;
    PowBalanceQuadSoft = 0;
    subModelCtg.reslim = 100; # ??
    solve subModelCtg using nlp minimizing objVar;
  );

  # if still infeasible, something is wrong
  if((subModelCtg.maxinfes > 1e-4) or (subModelCtg.modelstat = 5) or (subModelCtg.modelstat = 6),
    option solprint = on;
    subModelCtg.reslim = 300;
    solve subModelCtg using nlp minimizing objVar;
    abort 'infeasibility still';
  );

  CtgMaxInfeas(k0) = subModelCtg.maxinfes;
  CtgMaxCompViol(k0) = max(
    smax((i,j)$GenCtgParticipating(i,j,k0),
        min(genCtgPowRealOverVar.l(i,j,k0), 
            genCtgPowRealVar.l(i,j,k0) - GenPMin(i,j))),
    smax((i,j)$GenCtgParticipating(i,j,k0),
        min(genCtgPowRealUnderVar.l(i,j,k0), 
            GenPMax(i,j) - genCtgPowRealVar.l(i,j,k0))),
    smax(i$BusCtgVoltMagMaintDom(i,k0),
        min(busCtgVoltMagOverVar.l(i,k0),
	    sum((i1,j1)$GenBusCtgVoltMagMaintDom(i1,j1,i,k0),
	        genCtgPowImagVar.l(i1,j1,k0) - GenQMin(i1,j1)))),
    smax(i$BusCtgVoltMagMaintDom(i,k0),
        min(busCtgVoltMagUnderVar.l(i,k0),
	    sum((i1,j1)$GenBusCtgVoltMagMaintDom(i1,j1,i,k0),
	        GenQMax(i1,j1) - genCtgPowImagVar.l(i1,j1,k0)))));
  CtgMaxSoftViol(k0) = max(
    smax(i$Bus(i),
         busCtgPowRealBalanceViolVar.l(i,k0))$PowBalanceQuadSoft,
    smax(i$Bus(i),
         busCtgPowRealBalanceOverViolVar.l(i,k0))$PowBalanceLinSoft,
    smax(i$Bus(i),
         busCtgPowRealBalanceUnderViolVar.l(i,k0))$PowBalanceLinSoft,
    smax(i$Bus(i),
         busCtgPowImagBalanceViolVar.l(i,k0))$PowBalanceQuadSoft,
    smax(i$Bus(i),
         busCtgPowImagBalanceOverViolVar.l(i,k0))$PowBalanceLinSoft,
    smax(i$Bus(i),
         busCtgPowImagBalanceUnderViolVar.l(i,k0))$PowBalanceLinSoft,
    smax((i1,i2,j)$LineCtgActive(i1,i2,j,k0),
         lineCtgCurrMag1BoundViolVar.l(i1,i2,j,k0))$CurrFlowMagBoundSoft,
    smax((i1,i2,j)$LineCtgActive(i1,i2,j,k0),
         lineCtgCurrMag2BoundViolVar.l(i1,i2,j,k0))$CurrFlowMagBoundSoft,
    smax((i1,i2,j)$XfmrCtgActive(i1,i2,j,k0),
         xfmrCtgPowMag1BoundViolVar.l(i1,i2,j,k0))$PowFlowMagBoundSoft,
    smax((i1,i2,j)$XfmrCtgActive(i1,i2,j,k0),
         xfmrCtgPowMag2BoundViolVar.l(i1,i2,j,k0))$PowFlowMagBoundSoft);

  put LogDetailed;
  put
    k0.tl:0:0 ', '
    timeelapsed:0:0 ', '
    CtgMaxInfeas(k0):0:10 ', '
    CtgMaxCompViol(k0):0:10 ', '
    CtgMaxSoftViol(k0):0:10 /;
  putclose;

);

$ontext
# log entry
put LogFile;
put
  timeelapsed:0:0 ', '
  MaxInfeas:0:10 ', '
  max(
    smax(i$Bus(i),
         busPowRealBalanceViolVar.l(i))$PowBalanceQuadSoft,
    smax(i$Bus(i),
         busPowRealBalanceOverViolVar.l(i))$PowBalanceLinSoft,
    smax(i$Bus(i),
	 busPowRealBalanceUnderViolVar.l(i))$PowBalanceLinSoft,
    smax(i$Bus(i),
	 busPowImagBalanceViolVar.l(i))$PowBalanceQuadSoft,
    smax(i$Bus(i),
	 busPowImagBalanceOverViolVar.l(i))$PowBalanceLinSoft,
    smax(i$Bus(i),
	 busPowImagBalanceUnderViolVar.l(i))$PowBalanceLinSoft,
    smax((i1,i2,j)$LineActive(i1,i2,j),
         lineCurrMag1BoundViolVar.l(i1,i2,j))$CurrFlowMagBoundSoft,
    smax((i1,i2,j)$LineActive(i1,i2,j),
         lineCurrMag2BoundViolVar.l(i1,i2,j))$CurrFlowMagBoundSoft,
    smax((i1,i2,j)$XfmrActive(i1,i2,j),
         xfmrPowMag1BoundViolVar.l(i1,i2,j))$PowFlowMagBoundSoft,
    smax((i1,i2,j)$XfmrActive(i1,i2,j),
         xfmrPowMag2BoundViolVar.l(i1,i2,j))$PowFlowMagBoundSoft):0:10 ', '
  smax(k$Ctg(k), CtgMaxInfeas(k)):0:10 ', '
  smax(k$Ctg(k), CtgMaxCompViol(k)):0:10 ', '
  smax(k$Ctg(k), CtgMaxSoftViol(k)):0:10 /;
putclose;
$offtext

$endif

$include solution_evaluation_ctg.gms
$include convert_solution_ctg.gms
$include write_solution_ctg.gms
