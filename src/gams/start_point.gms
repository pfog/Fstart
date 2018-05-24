$title start_point
$ontext
set a start point for the optimization model
$offtext

* variables
$ontext
  totalCostVar
  genCostVar(i,j)
  genPlCoeffVar(i,j,i1)
  busVoltMagVar(i)
  busVoltAngVar(i)
  loadPowRealVar(i,j)
  loadPowImagVar(i,j)
  fxshPowRealVar(i,j)
  fxshPowImagVar(i,j)
  genPowRealVar(i,j)
  genPowImagVar(i,j)
  lineCurrReal1Var(i1,i2,j)
  lineCurrImag1Var(i1,i2,j)
  lineCurrReal2Var(i1,i2,j)
  lineCurrImag2Var(i1,i2,j)
  linePowReal1Var(i1,i2,j)
  linePowImag1Var(i1,i2,j)
  linePowReal2Var(i1,i2,j)
  linePowImag2Var(i1,i2,j)
  xfmrCurrReal1Var(i1,i2,j)
  xfmrCurrImag1Var(i1,i2,j)
  xfmrCurrReal2Var(i1,i2,j)
  xfmrCurrImag2Var(i1,i2,j)
  xfmrPowReal1Var(i1,i2,j)
  xfmrPowImag1Var(i1,i2,j)
  xfmrPowReal2Var(i1,i2,j)
  xfmrPowImag2Var(i1,i2,j)
  swshPowImagVar(i)
  swshAdmImagVar(i)
  busCtgVoltMagVar(i,k)
  busCtgVoltAngVar(i,k)
  loadCtgPowRealVar(i,j,k)
  loadCtgPowImagVar(i,j,k)
  fxshCtgPowRealVar(i,j,k)
  fxshCtgPowImagVar(i,j,k)
  genCtgPowRealVar(i,j,k)
  genCtgPowImagVar(i,j,k)
  lineCtgCurrReal1Var(i1,i2,j,k)
  lineCtgCurrImag1Var(i1,i2,j,k)
  lineCtgCurrReal2Var(i1,i2,j,k)
  lineCtgCurrImag2Var(i1,i2,j,k)
  lineCtgPowReal1Var(i1,i2,j,k)
  lineCtgPowImag1Var(i1,i2,j,k)
  lineCtgPowReal2Var(i1,i2,j,k)
  lineCtgPowImag2Var(i1,i2,j,k)
  xfmrCtgCurrReal1Var(i1,i2,j,k)
  xfmrCtgCurrImag1Var(i1,i2,j,k)
  xfmrCtgCurrReal2Var(i1,i2,j,k)
  xfmrCtgCurrImag2Var(i1,i2,j,k)
  xfmrCtgPowReal1Var(i1,i2,j,k)
  xfmrCtgPowImag1Var(i1,i2,j,k)
  xfmrCtgPowReal2Var(i1,i2,j,k)
  xfmrCtgPowImag2Var(i1,i2,j,k)
  swshCtgPowImagVar(i,k)
  swshCtgAdmImagVar(i,k)
  areaCtgPowRealChangeVar(i,k);
$offtext

# flat start
$ontext
busVoltMagVar.l(i) = 1.0;
busVoltAngVar.l(i) = 0.0;
$offtext

# random start
$ontext
$offtext

# start from data
#$ontext
busVoltMagVar.l(i) = BusVM(i);
busVoltAngVar.l(i) = BusVA(i);
genPowRealVar.l(i,j)$GenActive(i,j) = GenP(i,j);
genPowImagVar.l(i,j)$GenActive(i,j) = GenQ(i,j);
swshAdmImagVar.l(i)$SwshActive(i) = SwshBInit(i);
#$offtext
