$title equation_definitions
$ontext
define equations used by the Grid Optimization Competitions models
$offtext

objDefFull..
      objVar
  =e= costVar
   +  sum(i$Bus(i),
            PowBalanceViolPen * busPowRealBalanceOverViolVar(i)
          + PowBalanceViolPen * busPowRealBalanceUnderViolVar(i)
	  + PowBalanceViolPen * busPowImagBalanceOverViolVar(i)
          + PowBalanceViolPen * busPowImagBalanceUnderViolVar(i))$PowBalanceLinSoft
   +  PowBalanceViolPen * sum(i$Bus(i),
          sqrt(sqr(AbsSmoothing) + sqr(busPowRealBalanceViolVar(i))))$PowBalanceQuadSoft
   +  PowBalanceViolPen * sum(i$Bus(i),
          sqrt(sqr(AbsSmoothing) + sqr(busPowImagBalanceViolVar(i))))$PowBalanceQuadSoft
#   +  PowBalanceViolPen * sqrt(1e0 + sum(i$Bus(i),
#            sqr(busPowRealBalanceViolVar(i))
#          + sqr(busPowImagBalanceViolVar(i))))$PowBalanceQuadSoft
#   +  PowBalanceViolPen * sum(i$Bus(i),
#            sqr(busPowRealBalanceViolVar(i))
#          + sqr(busPowImagBalanceViolVar(i)))$PowBalanceQuadSoft
   +  sum((i1,i2,j)$LineActive(i1,i2,j),
          CurrFlowMagViolPen * lineCurrMag1BoundViolVar(i1,i2,j))$CurrFlowMagBoundSoft
   +  sum((i1,i2,j)$LineActive(i1,i2,j),
          CurrFlowMagViolPen * lineCurrMag2BoundViolVar(i1,i2,j))$CurrFlowMagBoundSoft
   +  sum((i1,i2,j)$XfmrActive(i1,i2,j),
          PowFlowMagViolPen * xfmrPowMag1BoundViolVar(i1,i2,j))$PowFlowMagBoundSoft
   +  sum((i1,i2,j)$XfmrActive(i1,i2,j),
          PowFlowMagViolPen * xfmrPowMag2BoundViolVar(i1,i2,j))$PowFlowMagBoundSoft
   +  sum((i,k)$(Bus(i) and CtgActive(k)),
            PowBalanceViolPen * busCtgPowRealBalanceOverViolVar(i,k)
          + PowBalanceViolPen * busCtgPowRealBalanceUnderViolVar(i,k)
	  + PowBalanceViolPen * busCtgPowImagBalanceOverViolVar(i,k)
          + PowBalanceViolPen * busCtgPowImagBalanceUnderViolVar(i,k))$PowBalanceLinSoft
   +  PowBalanceViolPen * sum((i,k)$(Bus(i) and CtgActive(k)),
          sqrt(sqr(AbsSmoothing) + sqr(busCtgPowRealBalanceViolVar(i,k))))$PowBalanceQuadSoft
   +  PowBalanceViolPen * sum((i,k)$(Bus(i) and CtgActive(k)),
          sqrt(sqr(AbsSmoothing) + sqr(busCtgPowImagBalanceViolVar(i,k))))$PowBalanceQuadSoft
#   +  PowBalanceViolPen * sqrt(1e0 + sum((i,k)$(Bus(i) and CtgActive(k)),
#            sqr(busCtgPowRealBalanceViolVar(i,k))
#          + sqr(busCtgPowImagBalanceViolVar(i,k))))$PowBalanceQuadSoft
#   +  PowBalanceViolPen * sum((i,k)$(Bus(i) and CtgActive(k)),
#            sqr(busCtgPowRealBalanceViolVar(i,k))
#          + sqr(busCtgPowImagBalanceViolVar(i,k)))$PowBalanceQuadSoft
   +  sum((i1,i2,j,k)$(LineCtgActive(i1,i2,j,k) and CtgActive(k)),
          CurrFlowMagViolPen * lineCtgCurrMag1BoundViolVar(i1,i2,j,k))$CurrFlowMagBoundSoft
   +  sum((i1,i2,j,k)$(LineCtgActive(i1,i2,j,k) and CtgActive(k)),
          CurrFlowMagViolPen * lineCtgCurrMag2BoundViolVar(i1,i2,j,k))$CurrFlowMagBoundSoft
   +  sum((i1,i2,j,k)$(XfmrCtgActive(i1,i2,j,k) and CtgActive(k)),
          PowFlowMagViolPen * xfmrCtgPowMag1BoundViolVar(i1,i2,j,k))$PowFlowMagBoundSoft
   +  sum((i1,i2,j,k)$(XfmrCtgActive(i1,i2,j,k) and CtgActive(k)),
          PowFlowMagViolPen * xfmrCtgPowMag2BoundViolVar(i1,i2,j,k))$PowFlowMagBoundSoft
#   +  sum((i,j,k)$(GenCtgParticipating(i,j,k) and CtgActive(k)),
#          1000 * genCtgPowRealOverVar(i,j,k))
#   +  sum((i,j,k)$(GenCtgParticipating(i,j,k) and CtgActive(k)),
#          1000 * genCtgPowRealUnderVar(i,j,k))
   +  sum((i,k)$(BusCtgVoltMagMaintDom(i,k) and CtgActive(k)),
          VoltMaintViolPen * busCtgVoltMagOverVar(i,k))$VoltMaintLinSoft
   +  sum((i,k)$(BusCtgVoltMagMaintDom(i,k) and CtgActive(k)),
          VoltMaintViolPen * busCtgVoltMagUnderVar(i,k))$VoltMaintLinSoft
   +  sum((i,k)$(BusCtgVoltMagMaintDom(i,k) and CtgActive(k)),
          VoltMaintViolPen * sqrt(1e-6 + sqr(busCtgVoltMagVar(i,k) - busVoltMagVar(i))))$VoltMaintQuadSoft;

objDefSubBase..
      objVar
  =e= costVar
   +  sum((i,j)$GenActive(i,j),
            GenPowRealCoeff1(i,j) * GenPowRealVar(i,j)
	  + GenPowRealCoeff2(i,j) * sqr(GenPowRealVar(i,j)))
   +  sum(i$Bus(i),
            BusVoltMagCoeff1(i) * BusVoltMagVar(i)
	  + BusVoltMagCoeff2(i) * sqr(BusVoltMagVar(i)))
   +  PowBalanceViolPen * sum(i$Bus(i),
            busPowRealBalanceOverViolVar(i)
          + busPowRealBalanceUnderViolVar(i)
	  + busPowImagBalanceOverViolVar(i)
          + busPowImagBalanceUnderViolVar(i))$PowBalanceLinSoft
   +  PowBalanceViolPen * sum(i$Bus(i),
            sqrt(sqr(AbsSmoothing) + sqr(busPowRealBalanceViolVar(i)))
          + sqrt(sqr(AbsSmoothing) + sqr(busPowImagBalanceViolVar(i))))$PowBalanceQuadSoft

# version with substitution and simplification for power
$ontext
   +  PowBalanceViolPen * sum(i$Bus(i), # !!!!
            sqrt(sqr(AbsSmoothing) + sqr(
                sum(j$GenActive(i,j), genPowRealVar(i,j))
              - sum(j$LoadActive(i,j), LoadP(i,j))
              - sum(j$FxshActive(i,j), FxshG(i,j) * sqr(busVoltMagVar(i)))
              - sum((i2,j)$LineActive(i,i2,j),
                      lineCurrReal1Var(i,i2,j) * busVoltMagVar(i) * cos(busVoltAngVar(i))
                    + lineCurrImag1Var(i,i2,j) * busVoltMagVar(i) * sin(busVoltAngVar(i)))
              - sum((i1,j)$LineActive(i1,i,j),
                      lineCurrReal2Var(i1,i,j) * busVoltMagVar(i) * cos(busVoltAngVar(i))
                    + lineCurrImag2Var(i1,i,j) * busVoltMagVar(i) * sin(busVoltAngVar(i)))
              - sum((i2,j)$XfmrActive(i,i2,j),
                      xfmrCurrReal1Var(i,i2,j) * busVoltMagVar(i) * cos(busVoltAngVar(i))
                    + xfmrCurrImag1Var(i,i2,j) * busVoltMagVar(i) * sin(busVoltAngVar(i)))
              - sum((i1,j)$XfmrActive(i1,i,j),
                      xfmrCurrReal2Var(i1,i,j) * busVoltMagVar(i) * cos(busVoltAngVar(i))
                    + xfmrCurrImag2Var(i1,i,j) * busVoltMagVar(i) * sin(busVoltAngVar(i)))))
          + sqrt(sqr(AbsSmoothing) + sqr(
                sum(j$GenActive(i,j), genPowImagVar(i,j))
              - sum(j$LoadActive(i,j), LoadQ(i,j))
              - sum(j$FxshActive(i,j), - FxshB(i,j) * sqr(busVoltMagVar(i)))
              - sum((i2,j)$LineActive(i,i2,j),
                      lineCurrReal1Var(i,i2,j) * busVoltMagVar(i) * sin(busVoltAngVar(i))
                    - lineCurrImag1Var(i,i2,j) * busVoltMagVar(i) * cos(busVoltAngVar(i)))
              - sum((i1,j)$LineActive(i1,i,j),
                      lineCurrReal2Var(i1,i,j) * busVoltMagVar(i) * sin(busVoltAngVar(i))
                    - lineCurrImag2Var(i1,i,j) * busVoltMagVar(i) * cos(busVoltAngVar(i)))
              - sum((i2,j)$XfmrActive(i,i2,j),
                      xfmrCurrReal1Var(i,i2,j) * busVoltMagVar(i) * sin(busVoltAngVar(i))
                    - xfmrCurrImag1Var(i,i2,j) * busVoltMagVar(i) * cos(busVoltAngVar(i)))
              - sum((i1,j)$XfmrActive(i1,i,j),
                      xfmrCurrReal2Var(i1,i,j) * busVoltMagVar(i) * sin(busVoltAngVar(i))
                    - xfmrCurrImag2Var(i1,i,j) * busVoltMagVar(i) * cos(busVoltAngVar(i)))
              + swshAdmImagVar(i) * sqr(busVoltMagVar(i)) $SwshActive(i))))$PowBalanceDense
$offtext









# version with just direct substitutione for current and voltage
   +  PowBalanceViolPen * sum(i$Bus(i), -2*AbsSmoothing # !!!!
          + sqrt(sqr(AbsSmoothing) + sqr(
                sum(j$GenActive(i,j), genPowRealVar(i,j))
              - sum(j$LoadActive(i,j), LoadP(i,j))
              - sum(j$FxshActive(i,j), FxshG(i,j) * sqr(busVoltMagVar(i)))
              - sum((i2,j)$LineActive(i,i2,j),
                      (
LineG(i,i2,j) * (busVoltMagVar(i) * cos(busVoltAngVar(i)) - busVoltMagVar(i2) * cos(busVoltAngVar(i2)))
   -  LineB(i,i2,j) * (busVoltMagVar(i) * sin(busVoltAngVar(i)) - busVoltMagVar(i2) * sin(busVoltAngVar(i2)))
   -  0.5 * LineBCh(i,i2,j) * busVoltMagVar(i) * sin(busVoltAngVar(i))
) * busVoltMagVar(i) * cos(busVoltAngVar(i))
                    + (
LineG(i,i2,j) * (busVoltMagVar(i) * sin(busVoltAngVar(i)) - busVoltMagVar(i2) * sin(busVoltAngVar(i2)))
   +  LineB(i,i2,j) * (busVoltMagVar(i) * cos(busVoltAngVar(i)) - busVoltMagVar(i2) * cos(busVoltAngVar(i2)))
   +  0.5 * LineBCh(i,i2,j) * busVoltMagVar(i) * cos(busVoltAngVar(i))
) * busVoltMagVar(i) * sin(busVoltAngVar(i)))
              - sum((i1,j)$LineActive(i1,i,j),
                      (
LineG(i1,i,j) * (busVoltMagVar(i) * cos(busVoltAngVar(i)) - busVoltMagVar(i1) * cos(busVoltAngVar(i1)))
   -  LineB(i1,i,j) * (busVoltMagVar(i) * sin(busVoltAngVar(i)) - busVoltMagVar(i1) * sin(busVoltAngVar(i1)))
   -  0.5 * LineBCh(i1,i,j) * busVoltMagVar(i) * sin(busVoltAngVar(i))
) * busVoltMagVar(i) * cos(busVoltAngVar(i))
                    + (
LineG(i1,i,j) * (busVoltMagVar(i) * sin(busVoltAngVar(i)) - busVoltMagVar(i1) * sin(busVoltAngVar(i1)))
   +  LineB(i1,i,j) * (busVoltMagVar(i) * cos(busVoltAngVar(i)) - busVoltMagVar(i1) * cos(busVoltAngVar(i1)))
   +  0.5 * LineBCh(i1,i,j) * busVoltMagVar(i) * cos(busVoltAngVar(i))
) * busVoltMagVar(i) * sin(busVoltAngVar(i)))
              - sum((i2,j)$XfmrActive(i,i2,j),
                      (

###
#xfmrCurrReal1Var(i,i2,j)
XfmrGMag(i,i2,j) * busVoltMagVar(i) * cos(busVoltAngVar(i))
   -  XfmrBMag(i,i2,j) * busVoltMagVar(i) * sin(busVoltAngVar(i))
   +    (  XfmrG(i,i2,j) * cos(XfmrAng(i,i2,j))
         - XfmrB(i,i2,j) * sin(XfmrAng(i,i2,j)))
      * (  busVoltMagVar(i) / sqr(XfmrRatio(i,i2,j)) * cos(busVoltAngVar(i) - XfmrAng(i,i2,j))
         - busVoltMagVar(i2) / XfmrRatio(i,i2,j) * cos(busVoltAngVar(i2)))
   -    (  XfmrG(i,i2,j) * sin(XfmrAng(i,i2,j))
         + XfmrB(i,i2,j) * cos(XfmrAng(i,i2,j)))
      * (  busVoltMagVar(i) / sqr(XfmrRatio(i,i2,j)) * sin(busVoltAngVar(i) - XfmrAng(i,i2,j))
         - busVoltMagVar(i2) / XfmrRatio(i,i2,j) * sin(busVoltAngVar(i2)))
) * busVoltMagVar(i) * cos(busVoltAngVar(i))
                    + (

###
#xfmrCurrImag1Var(i,i2,j)
XfmrGMag(i,i2,j) * busVoltMagVar(i) * sin(busVoltAngVar(i))
   +  XfmrBMag(i,i2,j) * busVoltMagVar(i) * cos(busVoltAngVar(i))
   +    (  XfmrG(i,i2,j) * cos(XfmrAng(i,i2,j))
         - XfmrB(i,i2,j) * sin(XfmrAng(i,i2,j)))
      * (  busVoltMagVar(i) / sqr(XfmrRatio(i,i2,j)) * sin(busVoltAngVar(i) - XfmrAng(i,i2,j))
         - busVoltMagVar(i2) / XfmrRatio(i,i2,j) * sin(busVoltAngVar(i2)))
   +    (  XfmrG(i,i2,j) * sin(XfmrAng(i,i2,j))
         + XfmrB(i,i2,j) * cos(XfmrAng(i,i2,j)))
      * (  busVoltMagVar(i) / sqr(XfmrRatio(i,i2,j)) * cos(busVoltAngVar(i) - XfmrAng(i,i2,j))
         - busVoltMagVar(i2) / XfmrRatio(i,i2,j) * cos(busVoltAngVar(i2)))
) * busVoltMagVar(i) * sin(busVoltAngVar(i)))
              - sum((i1,j)$XfmrActive(i1,i,j),
                      (

###
#xfmrCurrReal2Var(i1,i,j)
XfmrG(i1,i,j)
      * (  busVoltMagVar(i) * cos(busVoltAngVar(i))
         - busVoltMagVar(i1) / XfmrRatio(i1,i,j) * cos(busVoltAngVar(i1) - XfmrAng(i1,i,j)))
   -    XfmrB(i1,i,j)
      * (  busVoltMagVar(i) * sin(busVoltAngVar(i))
         - busVoltMagVar(i1) / XfmrRatio(i1,i,j) * sin(busVoltAngVar(i1) - XfmrAng(i1,i,j)))
) * busVoltMagVar(i) * cos(busVoltAngVar(i))
                    + (

###
#xfmrCurrImag2Var(i1,i,j)
XfmrG(i1,i,j)
      * (  busVoltMagVar(i) * sin(busVoltAngVar(i))
         - busVoltMagVar(i1) / XfmrRatio(i1,i,j) * sin(busVoltAngVar(i1) - XfmrAng(i1,i,j)))
   +    XfmrB(i1,i,j)
      * (  busVoltMagVar(i) * cos(busVoltAngVar(i))
         - busVoltMagVar(i1) / XfmrRatio(i1,i,j) * cos(busVoltAngVar(i1) - XfmrAng(i1,i,j)))
) * busVoltMagVar(i) * sin(busVoltAngVar(i)))))
          + sqrt(sqr(AbsSmoothing) + sqr(
                sum(j$GenActive(i,j), genPowImagVar(i,j))
              - sum(j$LoadActive(i,j), LoadQ(i,j))
              - sum(j$FxshActive(i,j), - FxshB(i,j) * sqr(busVoltMagVar(i)))
              - sum((i2,j)$LineActive(i,i2,j),
                      (
LineG(i,i2,j) * (busVoltMagVar(i) * cos(busVoltAngVar(i)) - busVoltMagVar(i2) * cos(busVoltAngVar(i2)))
   -  LineB(i,i2,j) * (busVoltMagVar(i) * sin(busVoltAngVar(i)) - busVoltMagVar(i2) * sin(busVoltAngVar(i2)))
   -  0.5 * LineBCh(i,i2,j) * busVoltMagVar(i) * sin(busVoltAngVar(i))
) * busVoltMagVar(i) * sin(busVoltAngVar(i))
                    - (
LineG(i,i2,j) * (busVoltMagVar(i) * sin(busVoltAngVar(i)) - busVoltMagVar(i2) * sin(busVoltAngVar(i2)))
   +  LineB(i,i2,j) * (busVoltMagVar(i) * cos(busVoltAngVar(i)) - busVoltMagVar(i2) * cos(busVoltAngVar(i2)))
   +  0.5 * LineBCh(i,i2,j) * busVoltMagVar(i) * cos(busVoltAngVar(i))
) * busVoltMagVar(i) * cos(busVoltAngVar(i)))
              - sum((i1,j)$LineActive(i1,i,j),
                      (
LineG(i1,i,j) * (busVoltMagVar(i) * cos(busVoltAngVar(i)) - busVoltMagVar(i1) * cos(busVoltAngVar(i1)))
   -  LineB(i1,i,j) * (busVoltMagVar(i) * sin(busVoltAngVar(i)) - busVoltMagVar(i1) * sin(busVoltAngVar(i1)))
   -  0.5 * LineBCh(i1,i,j) * busVoltMagVar(i) * sin(busVoltAngVar(i))
) * busVoltMagVar(i) * sin(busVoltAngVar(i))
                    - (
LineG(i1,i,j) * (busVoltMagVar(i) * sin(busVoltAngVar(i)) - busVoltMagVar(i1) * sin(busVoltAngVar(i1)))
   +  LineB(i1,i,j) * (busVoltMagVar(i) * cos(busVoltAngVar(i)) - busVoltMagVar(i1) * cos(busVoltAngVar(i1)))
   +  0.5 * LineBCh(i1,i,j) * busVoltMagVar(i) * cos(busVoltAngVar(i))
) * busVoltMagVar(i) * cos(busVoltAngVar(i)))
              - sum((i2,j)$XfmrActive(i,i2,j),
                      (

###
#xfmrCurrReal1Var(i,i2,j)
XfmrGMag(i,i2,j) * busVoltMagVar(i) * cos(busVoltAngVar(i))
   -  XfmrBMag(i,i2,j) * busVoltMagVar(i) * sin(busVoltAngVar(i))
   +    (  XfmrG(i,i2,j) * cos(XfmrAng(i,i2,j))
         - XfmrB(i,i2,j) * sin(XfmrAng(i,i2,j)))
      * (  busVoltMagVar(i) / sqr(XfmrRatio(i,i2,j)) * cos(busVoltAngVar(i) - XfmrAng(i,i2,j))
         - busVoltMagVar(i2) / XfmrRatio(i,i2,j) * cos(busVoltAngVar(i2)))
   -    (  XfmrG(i,i2,j) * sin(XfmrAng(i,i2,j))
         + XfmrB(i,i2,j) * cos(XfmrAng(i,i2,j)))
      * (  busVoltMagVar(i) / sqr(XfmrRatio(i,i2,j)) * sin(busVoltAngVar(i) - XfmrAng(i,i2,j))
         - busVoltMagVar(i2) / XfmrRatio(i,i2,j) * sin(busVoltAngVar(i2)))
) * busVoltMagVar(i) * sin(busVoltAngVar(i))
                    - (

###
#xfmrCurrImag1Var(i,i2,j)
XfmrGMag(i,i2,j) * busVoltMagVar(i) * sin(busVoltAngVar(i))
   +  XfmrBMag(i,i2,j) * busVoltMagVar(i) * cos(busVoltAngVar(i))
   +    (  XfmrG(i,i2,j) * cos(XfmrAng(i,i2,j))
         - XfmrB(i,i2,j) * sin(XfmrAng(i,i2,j)))
      * (  busVoltMagVar(i) / sqr(XfmrRatio(i,i2,j)) * sin(busVoltAngVar(i) - XfmrAng(i,i2,j))
         - busVoltMagVar(i2) / XfmrRatio(i,i2,j) * sin(busVoltAngVar(i2)))
   +    (  XfmrG(i,i2,j) * sin(XfmrAng(i,i2,j))
         + XfmrB(i,i2,j) * cos(XfmrAng(i,i2,j)))
      * (  busVoltMagVar(i) / sqr(XfmrRatio(i,i2,j)) * cos(busVoltAngVar(i) - XfmrAng(i,i2,j))
         - busVoltMagVar(i2) / XfmrRatio(i,i2,j) * cos(busVoltAngVar(i2)))
) * busVoltMagVar(i) * cos(busVoltAngVar(i)))
              - sum((i1,j)$XfmrActive(i1,i,j),
                      (

###
#xfmrCurrReal2Var(i1,i,j)
XfmrG(i1,i,j)
      * (  busVoltMagVar(i) * cos(busVoltAngVar(i))
         - busVoltMagVar(i1) / XfmrRatio(i1,i,j) * cos(busVoltAngVar(i1) - XfmrAng(i1,i,j)))
   -    XfmrB(i1,i,j)
      * (  busVoltMagVar(i) * sin(busVoltAngVar(i))
         - busVoltMagVar(i1) / XfmrRatio(i1,i,j) * sin(busVoltAngVar(i1) - XfmrAng(i1,i,j)))
) * busVoltMagVar(i) * sin(busVoltAngVar(i))
                    - (

###
#xfmrCurrImag2Var(i1,i,j)
XfmrG(i1,i,j)
      * (  busVoltMagVar(i) * sin(busVoltAngVar(i))
         - busVoltMagVar(i1) / XfmrRatio(i1,i,j) * sin(busVoltAngVar(i1) - XfmrAng(i1,i,j)))
   +    XfmrB(i1,i,j)
      * (  busVoltMagVar(i) * cos(busVoltAngVar(i))
         - busVoltMagVar(i1) / XfmrRatio(i1,i,j) * cos(busVoltAngVar(i1) - XfmrAng(i1,i,j)))
) * busVoltMagVar(i) * cos(busVoltAngVar(i)))
              + swshAdmImagVar(i) * sqr(busVoltMagVar(i)) $SwshActive(i))))$PowBalanceDense





#   +  PowBalanceViolPen * sqrt(1e0 + sum(i$Bus(i),
#            sqr(busPowRealBalanceViolVar(i))
#          + sqr(busPowImagBalanceViolVar(i))))$PowBalanceQuadSoft
#   +  PowBalanceViolPen * sum(i$Bus(i),
#            sqr(busPowRealBalanceViolVar(i))
#          + sqr(busPowImagBalanceViolVar(i)))$PowBalanceQuadSoft
   +  sum((i1,i2,j)$LineActive(i1,i2,j),
          CurrFlowMagViolPen * lineCurrMag1BoundViolVar(i1,i2,j))$CurrFlowMagBoundSoft
   +  sum((i1,i2,j)$LineActive(i1,i2,j),
          CurrFlowMagViolPen * lineCurrMag2BoundViolVar(i1,i2,j))$CurrFlowMagBoundSoft
   +  sum((i1,i2,j)$XfmrActive(i1,i2,j),
          PowFlowMagViolPen * xfmrPowMag1BoundViolVar(i1,i2,j))$PowFlowMagBoundSoft
   +  sum((i1,i2,j)$XfmrActive(i1,i2,j),
          PowFlowMagViolPen * xfmrPowMag2BoundViolVar(i1,i2,j))$PowFlowMagBoundSoft;

objDefSubCtg..
      objVar
  =e= sum((i,j)$GenActive(i,j),
            GenPowRealCoeff1(i,j) * GenPowRealVar(i,j)
	  + GenPowRealCoeff2(i,j) * sqr(GenPowRealVar(i,j)))
   +  sum(i$Bus(i),
            BusVoltMagCoeff1(i) * BusVoltMagVar(i)
	  + BusVoltMagCoeff2(i) * sqr(BusVoltMagVar(i)))
   +  PowBalanceViolPen * sum((i,k)$(Bus(i) and CtgActive(k)),
            busCtgPowRealBalanceOverViolVar(i,k)
          + busCtgPowRealBalanceUnderViolVar(i,k)
	  + busCtgPowImagBalanceOverViolVar(i,k)
          + busCtgPowImagBalanceUnderViolVar(i,k))$PowBalanceLinSoft
   +  PowBalanceViolPen * sum((i,k)$(Bus(i) and CtgActive(k)),
            sqrt(sqr(AbsSmoothing) + sqr(busCtgPowRealBalanceViolVar(i,k)))
          + sqrt(sqr(AbsSmoothing) + sqr(busCtgPowImagBalanceViolVar(i,k))))$PowBalanceQuadSoft
   +  PowBalanceViolPen * sum((i,k)$(Bus(i) and CtgActive(k)), # ????
            sqrt(sqr(AbsSmoothing) + sqr(
                sum(j$GenCtgActive(i,j,k), genCtgPowRealVar(i,j,k))
              - sum(j$LoadActive(i,j), LoadP(i,j))
              - sum(j$FxshActive(i,j), FxshG(i,j) * sqr(busCtgVoltMagVar(i,k)))
              - sum((i2,j)$LineCtgActive(i,i2,j,k), lineCtgPowReal1Var(i,i2,j,k))
              - sum((i1,j)$LineCtgActive(i1,i,j,k), lineCtgPowReal2Var(i1,i,j,k))
              - sum((i2,j)$XfmrCtgActive(i,i2,j,k), xfmrCtgPowReal1Var(i,i2,j,k))
              - sum((i1,j)$XfmrCtgActive(i1,i,j,k), xfmrCtgPowReal2Var(i1,i,j,k))))
          + sqrt(sqr(AbsSmoothing) + sqr(
                sum(j$GenCtgActive(i,j,k), genCtgPowImagVar(i,j,k))
              - sum(j$LoadActive(i,j), LoadQ(i,j))
              - sum(j$FxshActive(i,j), - FxshB(i,j) * sqr(busCtgVoltMagVar(i,k)))
              - sum((i2,j)$LineCtgActive(i,i2,j,k), lineCtgPowImag1Var(i,i2,j,k))
              - sum((i1,j)$LineCtgActive(i1,i,j,k), lineCtgPowImag2Var(i1,i,j,k))
              - sum((i2,j)$XfmrCtgActive(i,i2,j,k), xfmrCtgPowImag1Var(i,i2,j,k))
              - sum((i1,j)$XfmrCtgActive(i1,i,j,k), xfmrCtgPowImag2Var(i1,i,j,k))
              + swshCtgAdmImagVar(i,k) * sqr(busCtgVoltMagVar(i,k)) $SwshActive(i))))$PowBalanceDense
#   +  PowBalanceViolPen * sum((i,k)$(Bus(i) and CtgActive(k)),
#          sqrt(sqr(AbsSmoothing) + sqr(busCtgPowRealBalanceViolVar(i,k))))$PowBalanceQuadSoft
#   +  PowBalanceViolPen * sum((i,k)$(Bus(i) and CtgActive(k)),
#          sqrt(sqr(AbsSmoothing) + sqr(busCtgPowImagBalanceViolVar(i,k))))$PowBalanceQuadSoft
#   +  PowBalanceViolPen * sum((i,k)$(Bus(i) and CtgActive(k)),
#            sqr(busCtgPowRealBalanceViolVar(i,k))
#          + sqr(busCtgPowImagBalanceViolVar(i,k)))$PowBalanceQuadSoft
   +  sum((i1,i2,j,k)$(LineCtgActive(i1,i2,j,k) and CtgActive(k)),
          CurrFlowMagViolPen * lineCtgCurrMag1BoundViolVar(i1,i2,j,k))$CurrFlowMagBoundSoft
   +  sum((i1,i2,j,k)$(LineCtgActive(i1,i2,j,k) and CtgActive(k)),
          CurrFlowMagViolPen * lineCtgCurrMag2BoundViolVar(i1,i2,j,k))$CurrFlowMagBoundSoft
   +  sum((i1,i2,j,k)$(XfmrCtgActive(i1,i2,j,k) and CtgActive(k)),
          PowFlowMagViolPen * xfmrCtgPowMag1BoundViolVar(i1,i2,j,k))$PowFlowMagBoundSoft
   +  sum((i1,i2,j,k)$(XfmrCtgActive(i1,i2,j,k) and CtgActive(k)),
          PowFlowMagViolPen * xfmrCtgPowMag2BoundViolVar(i1,i2,j,k))$PowFlowMagBoundSoft
#   +  sum((i,j,k)$(GenCtgParticipating(i,j,k) and CtgActive(k)),
#          1000 * genCtgPowRealOverVar(i,j,k))
#   +  sum((i,j,k)$(GenCtgParticipating(i,j,k) and CtgActive(k)),
#          1000 * genCtgPowRealUnderVar(i,j,k))
   +  sum((i,k)$(BusCtgVoltMagMaintDom(i,k) and CtgActive(k)),
          VoltMaintViolPen * busCtgVoltMagOverVar(i,k))$VoltMaintLinSoft
   +  sum((i,k)$(BusCtgVoltMagMaintDom(i,k) and CtgActive(k)),
          VoltMaintViolPen * busCtgVoltMagUnderVar(i,k))$VoltMaintLinSoft
   +  sum((i,k)$(BusCtgVoltMagMaintDom(i,k) and CtgActive(k)),
          VoltMaintViolPen * sqrt(sqr(AbsSmoothing) + sqr(busCtgVoltMagVar(i,k) - busVoltMagVar(i))))$VoltMaintQuadSoft;

costDef..
      costVar
  =e= sum((i,j)$GenActive(i,j), GenCostVar(i,j));
  
genCostExpansionByPl(i,j)$GenActive(i,j)..
      genCostVar(i,j)
  =e= sum(i1$GenPl(i,j,i1), GenPlY(i,j,i1) * GenPlCoeffVar(i,j,i1));
  
genPowRealExpansionByPl(i,j)$GenActive(i,j)..
      genPowRealVar(i,j)
  =e= sum(i1$GenPl(i,j,i1), GenPlX(i,j,i1) * GenPlCoeffVar(i,j,i1));
  
genExpansionByPl(i,j)$GenActive(i,j)..
      1
  =e= sum(i1$GenPl(i,j,i1), GenPlCoeffVar(i,j,i1));

busPowRealBalance(i)$Bus(i)..
#busPowRealBalance(i)$(Bus(i) and (ord(i) ne 531) and (ord(i) ne 553))..
      sum(j$GenActive(i,j), genPowRealVar(i,j))
  =e= sum(j$LoadActive(i,j), loadPowRealVar(i,j))
   +  sum(j$FxshActive(i,j), fxshPowRealVar(i,j))
   +  sum((i2,j)$LineActive(i,i2,j), linePowReal1Var(i,i2,j))
   +  sum((i1,j)$LineActive(i1,i,j), linePowReal2Var(i1,i,j))
   +  sum((i2,j)$XfmrActive(i,i2,j), xfmrPowReal1Var(i,i2,j))
   +  sum((i1,j)$XfmrActive(i1,i,j), xfmrPowReal2Var(i1,i,j))
   +  busPowRealBalanceViolVar(i)$PowBalanceQuadSoft
   +  busPowRealBalanceOverViolVar(i)$PowBalanceLinSoft
   -  busPowRealBalanceUnderViolVar(i)$PowBalanceLinSoft;
  
busPowImagBalance(i)$Bus(i)..
#busPowImagBalance(i)$(Bus(i) and (ord(i) ne 531) and (ord(i) ne 553))..
      sum(j$GenActive(i,j), genPowImagVar(i,j))
  =e= sum(j$LoadActive(i,j), loadPowImagVar(i,j))
   +  sum(j$FxshActive(i,j), fxshPowImagVar(i,j))
   +  sum((i2,j)$LineActive(i,i2,j), linePowImag1Var(i,i2,j))
   +  sum((i1,j)$LineActive(i1,i,j), linePowImag2Var(i1,i,j))
   +  sum((i2,j)$XfmrActive(i,i2,j), xfmrPowImag1Var(i,i2,j))
   +  sum((i1,j)$XfmrActive(i1,i,j), xfmrPowImag2Var(i1,i,j))
   +  swshPowImagVar(i)$SwshActive(i)
   +  busPowImagBalanceOverViolVar(i)$PowBalanceQuadSoft
   +  busPowImagBalanceOverViolVar(i)$PowBalanceLinSoft
   -  busPowImagBalanceUnderViolVar(i)$PowBalanceLinSoft;
  
loadPowRealDef(i,j)$LoadActive(i,j)..
      loadPowRealVar(i,j)
  =e= LoadP(i,j);

loadPowImagDef(i,j)$LoadActive(i,j)..
      loadPowImagVar(i,j)
  =e= LoadQ(i,j);

fxshPowRealDef(i,j)$FxshActive(i,j)..
      fxshPowRealVar(i,j)
  =e= FxshG(i,j) * sqr(busVoltMagVar(i));

fxshPowImagDef(i,j)$FxshActive(i,j)..
      fxshPowImagVar(i,j)
  =e= -FxshB(i,j) * sqr(busVoltMagVar(i));

lineCurrMag1Bound(i1,i2,j)$LineActive(i1,i2,j)..
      sqrt(  1
           + sqr(lineCurrReal1Var(i1,i2,j))
           + sqr(lineCurrImag1Var(i1,i2,j)))
  =l= sqrt(  1
           + sqr(LineFlowMax(i1,i2,j)))
   +  lineCurrMag1BoundViolVar(i1,i2,j)$currFlowMagBoundSoft;

lineCurrMag2Bound(i1,i2,j)$LineActive(i1,i2,j)..
      sqrt(  1
           + sqr(lineCurrReal2Var(i1,i2,j))
           + sqr(lineCurrImag2Var(i1,i2,j)))
  =l= sqrt(  1
           + sqr(LineFlowMax(i1,i2,j)))
   +  lineCurrMag1BoundViolVar(i1,i2,j)$currFlowMagBoundSoft;

lineCurrMag1DenseBound(i1,i2,j)$LineActive(i1,i2,j)..
      sqrt(  1
           + sqr(
                 LineG(i1,i2,j) * (busVoltMagVar(i1) * cos(busVoltAngVar(i1)) - busVoltMagVar(i2) * cos(busVoltAngVar(i2)))
               - LineB(i1,i2,j) * (busVoltMagVar(i1) * sin(busVoltAngVar(i1)) - busVoltMagVar(i2) * sin(busVoltAngVar(i2)))
               - 0.5 * LineBCh(i1,i2,j) * busVoltMagVar(i1) * sin(busVoltAngVar(i1)))
           + sqr(
                 LineG(i1,i2,j) * (busVoltMagVar(i1) * sin(busVoltAngVar(i1)) - busVoltMagVar(i2) * sin(busVoltAngVar(i2)))
               + LineB(i1,i2,j) * (busVoltMagVar(i1) * cos(busVoltAngVar(i1)) - busVoltMagVar(i2) * cos(busVoltAngVar(i2)))
               + 0.5 * LineBCh(i1,i2,j) * busVoltMagVar(i1) * cos(busVoltAngVar(i1))))
  =l= sqrt(  1
           + sqr(LineFlowMax(i1,i2,j)))
   +  lineCurrMag1BoundViolVar(i1,i2,j)$currFlowMagBoundSoft;

lineCurrMag2DenseBound(i1,i2,j)$LineActive(i1,i2,j)..
      sqrt(  1
           + sqr(
                 LineG(i1,i2,j) * (busVoltMagVar(i2) * cos(busVoltAngVar(i2)) - busVoltMagVar(i1) * cos(busVoltAngVar(i1)))
               - LineB(i1,i2,j) * (busVoltMagVar(i2) * sin(busVoltAngVar(i2)) - busVoltMagVar(i1) * sin(busVoltAngVar(i1)))
               - 0.5 * LineBCh(i1,i2,j) * busVoltMagVar(i2) * sin(busVoltAngVar(i2)))
           + sqr(
                 LineG(i1,i2,j) * (busVoltMagVar(i2) * sin(busVoltAngVar(i2)) - busVoltMagVar(i1) * sin(busVoltAngVar(i1)))
               + LineB(i1,i2,j) * (busVoltMagVar(i2) * cos(busVoltAngVar(i2)) - busVoltMagVar(i1) * cos(busVoltAngVar(i1)))
               + 0.5 * LineBCh(i1,i2,j) * busVoltMagVar(i2) * cos(busVoltAngVar(i2))))
  =l= sqrt(  1
           + sqr(LineFlowMax(i1,i2,j)))
   +  lineCurrMag1BoundViolVar(i1,i2,j)$currFlowMagBoundSoft;

lineCurrReal1Def(i1,i2,j)$LineActive(i1,i2,j)..
      lineCurrReal1Var(i1,i2,j)
  =e= LineG(i1,i2,j) * (busVoltMagVar(i1) * cos(busVoltAngVar(i1)) - busVoltMagVar(i2) * cos(busVoltAngVar(i2)))
   -  LineB(i1,i2,j) * (busVoltMagVar(i1) * sin(busVoltAngVar(i1)) - busVoltMagVar(i2) * sin(busVoltAngVar(i2)))
   -  0.5 * LineBCh(i1,i2,j) * busVoltMagVar(i1) * sin(busVoltAngVar(i1));

lineCurrImag1Def(i1,i2,j)$LineActive(i1,i2,j)..
      lineCurrImag1Var(i1,i2,j)
  =e= LineG(i1,i2,j) * (busVoltMagVar(i1) * sin(busVoltAngVar(i1)) - busVoltMagVar(i2) * sin(busVoltAngVar(i2)))
   +  LineB(i1,i2,j) * (busVoltMagVar(i1) * cos(busVoltAngVar(i1)) - busVoltMagVar(i2) * cos(busVoltAngVar(i2)))
   +  0.5 * LineBCh(i1,i2,j) * busVoltMagVar(i1) * cos(busVoltAngVar(i1));

lineCurrReal2Def(i1,i2,j)$LineActive(i1,i2,j)..
      lineCurrReal2Var(i1,i2,j)
  =e= LineG(i1,i2,j) * (busVoltMagVar(i2) * cos(busVoltAngVar(i2)) - busVoltMagVar(i1) * cos(busVoltAngVar(i1)))
   -  LineB(i1,i2,j) * (busVoltMagVar(i2) * sin(busVoltAngVar(i2)) - busVoltMagVar(i1) * sin(busVoltAngVar(i1)))
   -  0.5 * LineBCh(i1,i2,j) * busVoltMagVar(i2) * sin(busVoltAngVar(i2));

lineCurrImag2Def(i1,i2,j)$LineActive(i1,i2,j)..
      lineCurrImag2Var(i1,i2,j)
  =e= LineG(i1,i2,j) * (busVoltMagVar(i2) * sin(busVoltAngVar(i2)) - busVoltMagVar(i1) * sin(busVoltAngVar(i1)))
   +  LineB(i1,i2,j) * (busVoltMagVar(i2) * cos(busVoltAngVar(i2)) - busVoltMagVar(i1) * cos(busVoltAngVar(i1)))
   +  0.5 * LineBCh(i1,i2,j) * busVoltMagVar(i2) * cos(busVoltAngVar(i2));

linePowReal1Def(i1,i2,j)$LineActive(i1,i2,j)..
      linePowReal1Var(i1,i2,j)
  =e= lineCurrReal1Var(i1,i2,j) * busVoltMagVar(i1) * cos(busVoltAngVar(i1))
   +  lineCurrImag1Var(i1,i2,j) * busVoltMagVar(i1) * sin(busVoltAngVar(i1));

linePowImag1Def(i1,i2,j)$LineActive(i1,i2,j)..
      linePowImag1Var(i1,i2,j)
  =e= lineCurrReal1Var(i1,i2,j) * busVoltMagVar(i1) * sin(busVoltAngVar(i1))
   -  lineCurrImag1Var(i1,i2,j) * busVoltMagVar(i1) * cos(busVoltAngVar(i1));

linePowReal2Def(i1,i2,j)$LineActive(i1,i2,j)..
      linePowReal2Var(i1,i2,j)
  =e= lineCurrReal2Var(i1,i2,j) * busVoltMagVar(i2) * cos(busVoltAngVar(i2))
   +  lineCurrImag2Var(i1,i2,j) * busVoltMagVar(i2) * sin(busVoltAngVar(i2));

linePowImag2Def(i1,i2,j)$LineActive(i1,i2,j)..
      linePowImag2Var(i1,i2,j)
  =e= lineCurrReal2Var(i1,i2,j) * busVoltMagVar(i2) * sin(busVoltAngVar(i2))
   -  lineCurrImag2Var(i1,i2,j) * busVoltMagVar(i2) * cos(busVoltAngVar(i2));

xfmrPowMag1Bound(i1,i2,j)$XfmrActive(i1,i2,j)..
      sqrt(  1
           + sqr(xfmrPowReal1Var(i1,i2,j))
           + sqr(xfmrPowImag1Var(i1,i2,j)))
  =l= sqrt(  1
           + sqr(XfmrFlowMax(i1,i2,j)))
   +  xfmrPowMag1BoundViolVar(i1,i2,j)$powFlowMagBoundSoft;

xfmrPowMag2Bound(i1,i2,j)$XfmrActive(i1,i2,j)..
      sqrt(  1
           + sqr(xfmrPowReal2Var(i1,i2,j))
           + sqr(xfmrPowImag2Var(i1,i2,j)))
  =l= sqrt(  1
           + sqr(XfmrFlowMax(i1,i2,j)))
   +  xfmrPowMag1BoundViolVar(i1,i2,j)$powFlowMagBoundSoft;

xfmrPowMag1DenseBound(i1,i2,j)$XfmrActive(i1,i2,j)..
      sqrt(  1
           + sqr(xfmrPowReal1Var(i1,i2,j))
           + sqr(xfmrPowImag1Var(i1,i2,j)))
  =l= sqrt(  1
           + sqr(XfmrFlowMax(i1,i2,j)))
   +  xfmrPowMag1BoundViolVar(i1,i2,j)$powFlowMagBoundSoft;

xfmrPowMag2DenseBound(i1,i2,j)$XfmrActive(i1,i2,j)..
      sqrt(  1
           + sqr(xfmrPowReal2Var(i1,i2,j))
           + sqr(xfmrPowImag2Var(i1,i2,j)))
  =l= sqrt(  1
           + sqr(XfmrFlowMax(i1,i2,j)))
   +  xfmrPowMag1BoundViolVar(i1,i2,j)$powFlowMagBoundSoft;

xfmrCurrReal1Def(i1,i2,j)$XfmrActive(i1,i2,j)..
      xfmrCurrReal1Var(i1,i2,j)
  =e= XfmrGMag(i1,i2,j) * busVoltMagVar(i1) * cos(busVoltAngVar(i1))
   -  XfmrBMag(i1,i2,j) * busVoltMagVar(i1) * sin(busVoltAngVar(i1))
   +    (  XfmrG(i1,i2,j) * cos(XfmrAng(i1,i2,j))
         - XfmrB(i1,i2,j) * sin(XfmrAng(i1,i2,j)))
      * (  busVoltMagVar(i1) / sqr(XfmrRatio(i1,i2,j)) * cos(busVoltAngVar(i1) - XfmrAng(i1,i2,j))
         - busVoltMagVar(i2) / XfmrRatio(i1,i2,j) * cos(busVoltAngVar(i2)))
   -    (  XfmrG(i1,i2,j) * sin(XfmrAng(i1,i2,j))
         + XfmrB(i1,i2,j) * cos(XfmrAng(i1,i2,j)))
      * (  busVoltMagVar(i1) / sqr(XfmrRatio(i1,i2,j)) * sin(busVoltAngVar(i1) - XfmrAng(i1,i2,j))
         - busVoltMagVar(i2) / XfmrRatio(i1,i2,j) * sin(busVoltAngVar(i2)));

xfmrCurrImag1Def(i1,i2,j)$XfmrActive(i1,i2,j)..
      xfmrCurrImag1Var(i1,i2,j)
  =e= XfmrGMag(i1,i2,j) * busVoltMagVar(i1) * sin(busVoltAngVar(i1))
   +  XfmrBMag(i1,i2,j) * busVoltMagVar(i1) * cos(busVoltAngVar(i1))
   +    (  XfmrG(i1,i2,j) * cos(XfmrAng(i1,i2,j))
         - XfmrB(i1,i2,j) * sin(XfmrAng(i1,i2,j)))
      * (  busVoltMagVar(i1) / sqr(XfmrRatio(i1,i2,j)) * sin(busVoltAngVar(i1) - XfmrAng(i1,i2,j))
         - busVoltMagVar(i2) / XfmrRatio(i1,i2,j) * sin(busVoltAngVar(i2)))
   +    (  XfmrG(i1,i2,j) * sin(XfmrAng(i1,i2,j))
         + XfmrB(i1,i2,j) * cos(XfmrAng(i1,i2,j)))
      * (  busVoltMagVar(i1) / sqr(XfmrRatio(i1,i2,j)) * cos(busVoltAngVar(i1) - XfmrAng(i1,i2,j))
         - busVoltMagVar(i2) / XfmrRatio(i1,i2,j) * cos(busVoltAngVar(i2)));

xfmrCurrReal2Def(i1,i2,j)$XfmrActive(i1,i2,j)..
      xfmrCurrReal2Var(i1,i2,j)
  =e=   XfmrG(i1,i2,j)
      * (  busVoltMagVar(i2) * cos(busVoltAngVar(i2))
         - busVoltMagVar(i1) / XfmrRatio(i1,i2,j) * cos(busVoltAngVar(i1) - XfmrAng(i1,i2,j)))
   -    XfmrB(i1,i2,j)
      * (  busVoltMagVar(i2) * sin(busVoltAngVar(i2))
         - busVoltMagVar(i1) / XfmrRatio(i1,i2,j) * sin(busVoltAngVar(i1) - XfmrAng(i1,i2,j)));

xfmrCurrImag2Def(i1,i2,j)$XfmrActive(i1,i2,j)..
      xfmrCurrImag2Var(i1,i2,j)
  =e=   XfmrG(i1,i2,j)
      * (  busVoltMagVar(i2) * sin(busVoltAngVar(i2))
         - busVoltMagVar(i1) / XfmrRatio(i1,i2,j) * sin(busVoltAngVar(i1) - XfmrAng(i1,i2,j)))
   +    XfmrB(i1,i2,j)
      * (  busVoltMagVar(i2) * cos(busVoltAngVar(i2))
         - busVoltMagVar(i1) / XfmrRatio(i1,i2,j) * cos(busVoltAngVar(i1) - XfmrAng(i1,i2,j)));

xfmrPowReal1Def(i1,i2,j)$XfmrActive(i1,i2,j)..
      xfmrPowReal1Var(i1,i2,j)
  =e= xfmrCurrReal1Var(i1,i2,j) * busVoltMagVar(i1) * cos(busVoltAngVar(i1))
   +  xfmrCurrImag1Var(i1,i2,j) * busVoltMagVar(i1) * sin(busVoltAngVar(i1));

xfmrPowImag1Def(i1,i2,j)$XfmrActive(i1,i2,j)..
      xfmrPowImag1Var(i1,i2,j)
  =e= xfmrCurrReal1Var(i1,i2,j) * busVoltMagVar(i1) * sin(busVoltAngVar(i1))
   -  xfmrCurrImag1Var(i1,i2,j) * busVoltMagVar(i1) * cos(busVoltAngVar(i1));

xfmrPowReal2Def(i1,i2,j)$XfmrActive(i1,i2,j)..
      xfmrPowReal2Var(i1,i2,j)
  =e= xfmrCurrReal2Var(i1,i2,j) * busVoltMagVar(i2) * cos(busVoltAngVar(i2))
   +  xfmrCurrImag2Var(i1,i2,j) * busVoltMagVar(i2) * sin(busVoltAngVar(i2));

xfmrPowImag2Def(i1,i2,j)$XfmrActive(i1,i2,j)..
      xfmrPowImag2Var(i1,i2,j)
  =e= xfmrCurrReal2Var(i1,i2,j) * busVoltMagVar(i2) * sin(busVoltAngVar(i2))
   -  xfmrCurrImag2Var(i1,i2,j) * busVoltMagVar(i2) * cos(busVoltAngVar(i2));

swshPowImagDef(i)$SwshActive(i)..
      swshPowImagVar(i)
  =e= - swshAdmImagVar(i) * sqr(busVoltMagVar(i));

busCtgPowRealBalance(i,k)$(Bus(i) and CtgActive(k))..
#busCtgPowRealBalance(i,k)$(Bus(i) and CtgActive(k) and (ord(i) ne 1206))..
#busCtgPowRealBalance(i,k)$(Bus(i) and CtgActive(k) and (ord(i) = 1206))..
      sum(j$GenCtgActive(i,j,k), genCtgPowRealVar(i,j,k))
  =e= sum(j$LoadActive(i,j), loadCtgPowRealVar(i,j,k))
   +  sum(j$FxshActive(i,j), fxshCtgPowRealVar(i,j,k))
   +  sum((i2,j)$LineCtgActive(i,i2,j,k), lineCtgPowReal1Var(i,i2,j,k))
   +  sum((i1,j)$LineCtgActive(i1,i,j,k), lineCtgPowReal2Var(i1,i,j,k))
   +  sum((i2,j)$XfmrCtgActive(i,i2,j,k), xfmrCtgPowReal1Var(i,i2,j,k))
   +  sum((i1,j)$XfmrCtgActive(i1,i,j,k), xfmrCtgPowReal2Var(i1,i,j,k))
   +  busCtgPowRealBalanceViolVar(i,k)$PowBalanceQuadSoft
   +  busCtgPowRealBalanceOverViolVar(i,k)$PowBalanceLinSoft
   -  busCtgPowRealBalanceUnderViolVar(i,k)$PowBalanceLinSoft;
  
busCtgPowImagBalance(i,k)$(Bus(i) and CtgActive(k))..
      sum(j$GenCtgActive(i,j,k), genCtgPowImagVar(i,j,k))
  =e= sum(j$LoadActive(i,j), loadCtgPowImagVar(i,j,k))
   +  sum(j$FxshActive(i,j), fxshCtgPowImagVar(i,j,k))
   +  sum((i2,j)$LineCtgActive(i,i2,j,k), lineCtgPowImag1Var(i,i2,j,k))
   +  sum((i1,j)$LineCtgActive(i1,i,j,k), lineCtgPowImag2Var(i1,i,j,k))
   +  sum((i2,j)$XfmrCtgActive(i,i2,j,k), xfmrCtgPowImag1Var(i,i2,j,k))
   +  sum((i1,j)$XfmrCtgActive(i1,i,j,k), xfmrCtgPowImag2Var(i1,i,j,k))
   +  swshCtgPowImagVar(i,k)$SwshActive(i)
   +  busCtgPowImagBalanceViolVar(i,k)$PowBalanceQuadSoft
   +  busCtgPowImagBalanceOverViolVar(i,k)$PowBalanceLinSoft
   -  busCtgPowImagBalanceUnderViolVar(i,k)$PowBalanceLinSoft;

loadCtgPowRealDef(i,j,k)$(LoadActive(i,j) and CtgActive(k))..
      loadCtgPowRealVar(i,j,k)
  =e= LoadP(i,j);

loadCtgPowImagDef(i,j,k)$(LoadActive(i,j) and CtgActive(k))..
      loadCtgPowImagVar(i,j,k)
  =e= LoadQ(i,j);

fxshCtgPowRealDef(i,j,k)$(FxshActive(i,j) and CtgActive(k))..
      fxshCtgPowRealVar(i,j,k)
  =e= FxshG(i,j) * sqr(busCtgVoltMagVar(i,k));

fxshCtgPowImagDef(i,j,k)$(FxshActive(i,j) and CtgActive(k))..
      fxshCtgPowImagVar(i,j,k)
  =e= - FxshB(i,j) * sqr(busCtgVoltMagVar(i,k));

lineCtgCurrMag1Bound(i1,i2,j,k)$(LineCtgActive(i1,i2,j,k) and CtgActive(k))..
      sqrt(  1
           + sqr(lineCtgCurrReal1Var(i1,i2,j,k))
           + sqr(lineCtgCurrImag1Var(i1,i2,j,k)))
  =l= sqrt(  1
           + sqr(LineFlowMax(i1,i2,j)))
   +  lineCtgCurrMag1BoundViolVar(i1,i2,j,k)$currFlowMagBoundSoft;

lineCtgCurrMag2Bound(i1,i2,j,k)$(LineCtgActive(i1,i2,j,k) and CtgActive(k))..
      sqrt(  1
           + sqr(lineCtgCurrReal2Var(i1,i2,j,k))
           + sqr(lineCtgCurrImag2Var(i1,i2,j,k)))
  =l= sqrt(  1
           + sqr(LineFlowMax(i1,i2,j)))
   +  lineCtgCurrMag1BoundViolVar(i1,i2,j,k)$currFlowMagBoundSoft;

lineCtgCurrMag1DenseBound(i1,i2,j,k)$(LineCtgActive(i1,i2,j,k) and CtgActive(k))..
      sqrt(  1
           + sqr(
                 LineG(i1,i2,j) * (busCtgVoltMagVar(i1,k) * cos(busCtgVoltAngVar(i1,k)) - busCtgVoltMagVar(i2,k) * cos(busCtgVoltAngVar(i2,k)))
               - LineB(i1,i2,j) * (busCtgVoltMagVar(i1,k) * sin(busCtgVoltAngVar(i1,k)) - busCtgVoltMagVar(i2,k) * sin(busCtgVoltAngVar(i2,k)))
               - 0.5 * LineBCh(i1,i2,j) * busCtgVoltMagVar(i1,k) * sin(busCtgVoltAngVar(i1,k)))
           + sqr(
                 LineG(i1,i2,j) * (busCtgVoltMagVar(i1,k) * sin(busCtgVoltAngVar(i1,k)) - busCtgVoltMagVar(i2,k) * sin(busCtgVoltAngVar(i2,k)))
               + LineB(i1,i2,j) * (busCtgVoltMagVar(i1,k) * cos(busCtgVoltAngVar(i1,k)) - busCtgVoltMagVar(i2,k) * cos(busCtgVoltAngVar(i2,k)))
               + 0.5 * LineBCh(i1,i2,j) * busCtgVoltMagVar(i1,k) * cos(busCtgVoltAngVar(i1,k))))
  =l= sqrt(  1
           + sqr(LineFlowMax(i1,i2,j)))
   +  lineCtgCurrMag1BoundViolVar(i1,i2,j,k)$currFlowMagBoundSoft;

lineCtgCurrMag2DenseBound(i1,i2,j,k)$(LineCtgActive(i1,i2,j,k) and CtgActive(k))..
      sqrt(  1
           + sqr(
                 LineG(i1,i2,j) * (busCtgVoltMagVar(i2,k) * cos(busCtgVoltAngVar(i2,k)) - busCtgVoltMagVar(i1,k) * cos(busCtgVoltAngVar(i1,k)))
               - LineB(i1,i2,j) * (busCtgVoltMagVar(i2,k) * sin(busCtgVoltAngVar(i2,k)) - busCtgVoltMagVar(i1,k) * sin(busCtgVoltAngVar(i1,k)))
               - 0.5 * LineBCh(i1,i2,j) * busCtgVoltMagVar(i2,k) * sin(busCtgVoltAngVar(i2,k)))
           + sqr(

                 LineG(i1,i2,j) * (busCtgVoltMagVar(i2,k) * sin(busCtgVoltAngVar(i2,k)) - busCtgVoltMagVar(i1,k) * sin(busCtgVoltAngVar(i1,k)))
               + LineB(i1,i2,j) * (busCtgVoltMagVar(i2,k) * cos(busCtgVoltAngVar(i2,k)) - busCtgVoltMagVar(i1,k) * cos(busCtgVoltAngVar(i1,k)))
               + 0.5 * LineBCh(i1,i2,j) * busCtgVoltMagVar(i2,k) * cos(busCtgVoltAngVar(i2,k))))
  =l= sqrt(  1
           + sqr(LineFlowMax(i1,i2,j)))
   +  lineCtgCurrMag1BoundViolVar(i1,i2,j,k)$currFlowMagBoundSoft;

lineCtgCurrReal1Def(i1,i2,j,k)$(LineCtgActive(i1,i2,j,k) and CtgActive(k))..
      lineCtgCurrReal1Var(i1,i2,j,k)
  =e= LineG(i1,i2,j) * (busCtgVoltMagVar(i1,k) * cos(busCtgVoltAngVar(i1,k)) - busCtgVoltMagVar(i2,k) * cos(busCtgVoltAngVar(i2,k)))
   -  LineB(i1,i2,j) * (busCtgVoltMagVar(i1,k) * sin(busCtgVoltAngVar(i1,k)) - busCtgVoltMagVar(i2,k) * sin(busCtgVoltAngVar(i2,k)))
   -  0.5 * LineBCh(i1,i2,j) * busCtgVoltMagVar(i1,k) * sin(busCtgVoltAngVar(i1,k));

lineCtgCurrImag1Def(i1,i2,j,k)$(LineCtgActive(i1,i2,j,k) and CtgActive(k))..
      lineCtgCurrImag1Var(i1,i2,j,k)
  =e= LineG(i1,i2,j) * (busCtgVoltMagVar(i1,k) * sin(busCtgVoltAngVar(i1,k)) - busCtgVoltMagVar(i2,k) * sin(busCtgVoltAngVar(i2,k)))
   +  LineB(i1,i2,j) * (busCtgVoltMagVar(i1,k) * cos(busCtgVoltAngVar(i1,k)) - busCtgVoltMagVar(i2,k) * cos(busCtgVoltAngVar(i2,k)))
   +  0.5 * LineBCh(i1,i2,j) * busCtgVoltMagVar(i1,k) * cos(busCtgVoltAngVar(i1,k));

lineCtgCurrReal2Def(i1,i2,j,k)$(LineCtgActive(i1,i2,j,k) and CtgActive(k))..
      lineCtgCurrReal2Var(i1,i2,j,k)
  =e= LineG(i1,i2,j) * (busCtgVoltMagVar(i2,k) * cos(busCtgVoltAngVar(i2,k)) - busCtgVoltMagVar(i1,k) * cos(busCtgVoltAngVar(i1,k)))
   -  LineB(i1,i2,j) * (busCtgVoltMagVar(i2,k) * sin(busCtgVoltAngVar(i2,k)) - busCtgVoltMagVar(i1,k) * sin(busCtgVoltAngVar(i1,k)))
   -  0.5 * LineBCh(i1,i2,j) * busCtgVoltMagVar(i2,k) * sin(busCtgVoltAngVar(i2,k));

lineCtgCurrImag2Def(i1,i2,j,k)$(LineCtgActive(i1,i2,j,k) and CtgActive(k))..
      lineCtgCurrImag2Var(i1,i2,j,k)
  =e= LineG(i1,i2,j) * (busCtgVoltMagVar(i2,k) * sin(busCtgVoltAngVar(i2,k)) - busCtgVoltMagVar(i1,k) * sin(busCtgVoltAngVar(i1,k)))
   +  LineB(i1,i2,j) * (busCtgVoltMagVar(i2,k) * cos(busCtgVoltAngVar(i2,k)) - busCtgVoltMagVar(i1,k) * cos(busCtgVoltAngVar(i1,k)))
   +  0.5 * LineBCh(i1,i2,j) * busCtgVoltMagVar(i2,k) * cos(busCtgVoltAngVar(i2,k));

lineCtgPowReal1Def(i1,i2,j,k)$(LineCtgActive(i1,i2,j,k) and CtgActive(k))..
      lineCtgPowReal1Var(i1,i2,j,k)
  =e= lineCtgCurrReal1Var(i1,i2,j,k) * busCtgVoltMagVar(i1,k) * cos(busCtgVoltAngVar(i1,k))
   +  lineCtgCurrImag1Var(i1,i2,j,k) * busCtgVoltMagVar(i1,k) * sin(busCtgVoltAngVar(i1,k));

lineCtgPowImag1Def(i1,i2,j,k)$(LineCtgActive(i1,i2,j,k) and CtgActive(k))..
      lineCtgPowImag1Var(i1,i2,j,k)
  =e= lineCtgCurrReal1Var(i1,i2,j,k) * busCtgVoltMagVar(i1,k) * sin(busCtgVoltAngVar(i1,k))
   -  lineCtgCurrImag1Var(i1,i2,j,k) * busCtgVoltMagVar(i1,k) * cos(busCtgVoltAngVar(i1,k));

lineCtgPowReal2Def(i1,i2,j,k)$(LineCtgActive(i1,i2,j,k) and CtgActive(k))..
      lineCtgPowReal2Var(i1,i2,j,k)
  =e= lineCtgCurrReal2Var(i1,i2,j,k) * busCtgVoltMagVar(i2,k) * cos(busCtgVoltAngVar(i2,k))
   +  lineCtgCurrImag2Var(i1,i2,j,k) * busCtgVoltMagVar(i2,k) * sin(busCtgVoltAngVar(i2,k));

lineCtgPowImag2Def(i1,i2,j,k)$(LineCtgActive(i1,i2,j,k) and CtgActive(k))..
      lineCtgPowImag2Var(i1,i2,j,k)
  =e= lineCtgCurrReal2Var(i1,i2,j,k) * busCtgVoltMagVar(i2,k) * sin(busCtgVoltAngVar(i2,k))
   -  lineCtgCurrImag2Var(i1,i2,j,k) * busCtgVoltMagVar(i2,k) * cos(busCtgVoltAngVar(i2,k));

xfmrCtgPowMag1Bound(i1,i2,j,k)$(XfmrCtgActive(i1,i2,j,k) and CtgActive(k))..
      sqrt(  1
           + sqr(xfmrCtgPowReal1Var(i1,i2,j,k))
           + sqr(xfmrCtgPowImag1Var(i1,i2,j,k)))
  =l= sqrt(  1
           + sqr(XfmrFlowMax(i1,i2,j)))
   +  xfmrCtgPowMag1BoundViolVar(i1,i2,j,k)$powFlowMagBoundSoft;

xfmrCtgPowMag2Bound(i1,i2,j,k)$(XfmrCtgActive(i1,i2,j,k) and CtgActive(k))..
      sqrt(  1
           + sqr(xfmrCtgPowReal2Var(i1,i2,j,k))
           + sqr(xfmrCtgPowImag2Var(i1,i2,j,k)))
  =l= sqrt(  1
           + sqr(XfmrFlowMax(i1,i2,j)))
   +  xfmrCtgPowMag2BoundViolVar(i1,i2,j,k)$powFlowMagBoundSoft;

xfmrCtgPowMag1DenseBound(i1,i2,j,k)$(XfmrCtgActive(i1,i2,j,k) and CtgActive(k))..
      sqrt(  1
           + sqr(xfmrCtgPowReal1Var(i1,i2,j,k))
           + sqr(xfmrCtgPowImag1Var(i1,i2,j,k)))
  =l= sqrt(  1
           + sqr(XfmrFlowMax(i1,i2,j)))
   +  xfmrCtgPowMag1BoundViolVar(i1,i2,j,k)$powFlowMagBoundSoft;

xfmrCtgPowMag2DenseBound(i1,i2,j,k)$(XfmrCtgActive(i1,i2,j,k) and CtgActive(k))..
      sqrt(  1
           + sqr(xfmrCtgPowReal2Var(i1,i2,j,k))
           + sqr(xfmrCtgPowImag2Var(i1,i2,j,k)))
  =l= sqrt(  1
           + sqr(XfmrFlowMax(i1,i2,j)))
   +  xfmrCtgPowMag2BoundViolVar(i1,i2,j,k)$powFlowMagBoundSoft;

xfmrCtgCurrReal1Def(i1,i2,j,k)$(XfmrCtgActive(i1,i2,j,k) and CtgActive(k))..
      xfmrCtgCurrReal1Var(i1,i2,j,k)
  =e= XfmrGMag(i1,i2,j) * busCtgVoltMagVar(i1,k) * cos(busCtgVoltAngVar(i1,k))
   -  XfmrBMag(i1,i2,j) * busCtgVoltMagVar(i1,k) * sin(busCtgVoltAngVar(i1,k))
   +    (  XfmrG(i1,i2,j) * cos(XfmrAng(i1,i2,j))
         - XfmrB(i1,i2,j) * sin(XfmrAng(i1,i2,j)))
      * (  busCtgVoltMagVar(i1,k) / sqr(XfmrRatio(i1,i2,j)) * cos(busCtgVoltAngVar(i1,k) - XfmrAng(i1,i2,j))
         - busCtgVoltMagVar(i2,k) / XfmrRatio(i1,i2,j) * cos(busCtgVoltAngVar(i2,k)))
   -    (  XfmrG(i1,i2,j) * sin(XfmrAng(i1,i2,j))
         + XfmrB(i1,i2,j) * cos(XfmrAng(i1,i2,j)))
      * (  busCtgVoltMagVar(i1,k) / sqr(XfmrRatio(i1,i2,j)) * sin(busCtgVoltAngVar(i1,k) - XfmrAng(i1,i2,j))
         - busCtgVoltMagVar(i2,k) / XfmrRatio(i1,i2,j) * sin(busCtgVoltAngVar(i2,k)));

xfmrCtgCurrImag1Def(i1,i2,j,k)$(XfmrCtgActive(i1,i2,j,k) and CtgActive(k))..
      xfmrCtgCurrImag1Var(i1,i2,j,k)
  =e= XfmrGMag(i1,i2,j) * busCtgVoltMagVar(i1,k) * sin(busCtgVoltAngVar(i1,k))
   +  XfmrBMag(i1,i2,j) * busCtgVoltMagVar(i1,k) * cos(busCtgVoltAngVar(i1,k))
   +    (  XfmrG(i1,i2,j) * cos(XfmrAng(i1,i2,j))
         - XfmrB(i1,i2,j) * sin(XfmrAng(i1,i2,j)))
      * (  busCtgVoltMagVar(i1,k) / sqr(XfmrRatio(i1,i2,j)) * sin(busCtgVoltAngVar(i1,k) - XfmrAng(i1,i2,j))
         - busCtgVoltMagVar(i2,k) / XfmrRatio(i1,i2,j) * sin(busCtgVoltAngVar(i2,k)))
   +    (  XfmrG(i1,i2,j) * sin(XfmrAng(i1,i2,j))
         + XfmrB(i1,i2,j) * cos(XfmrAng(i1,i2,j)))
      * (  busCtgVoltMagVar(i1,k) / sqr(XfmrRatio(i1,i2,j)) * cos(busCtgVoltAngVar(i1,k) - XfmrAng(i1,i2,j))
         - busCtgVoltMagVar(i2,k) / XfmrRatio(i1,i2,j) * cos(busCtgVoltAngVar(i2,k)));

xfmrCtgCurrReal2Def(i1,i2,j,k)$(XfmrCtgActive(i1,i2,j,k) and CtgActive(k))..
      xfmrCtgCurrReal2Var(i1,i2,j,k)
  =e=   XfmrG(i1,i2,j)
      * (  busCtgVoltMagVar(i2,k) * cos(busCtgVoltAngVar(i2,k))
         - busCtgVoltMagVar(i1,k) / XfmrRatio(i1,i2,j) * cos(busCtgVoltAngVar(i1,k) - XfmrAng(i1,i2,j)))
   -    XfmrB(i1,i2,j)
      * (  busCtgVoltMagVar(i2,k) * sin(busCtgVoltAngVar(i2,k))
         - busCtgVoltMagVar(i1,k) / XfmrRatio(i1,i2,j) * sin(busCtgVoltAngVar(i1,k) - XfmrAng(i1,i2,j)));

xfmrCtgCurrImag2Def(i1,i2,j,k)$(XfmrCtgActive(i1,i2,j,k) and CtgActive(k))..
      xfmrCtgCurrImag2Var(i1,i2,j,k)
  =e=   XfmrG(i1,i2,j)
      * (  busCtgVoltMagVar(i2,k) * sin(busCtgVoltAngVar(i2,k))
         - busCtgVoltMagVar(i1,k) / XfmrRatio(i1,i2,j) * sin(busCtgVoltAngVar(i1,k) - XfmrAng(i1,i2,j)))
   +    XfmrB(i1,i2,j)
      * (  busCtgVoltMagVar(i2,k) * cos(busCtgVoltAngVar(i2,k))
         - busCtgVoltMagVar(i1,k) / XfmrRatio(i1,i2,j) * cos(busCtgVoltAngVar(i1,k) - XfmrAng(i1,i2,j)));

xfmrCtgPowReal1Def(i1,i2,j,k)$(XfmrCtgActive(i1,i2,j,k) and CtgActive(k))..
      xfmrCtgPowReal1Var(i1,i2,j,k)
  =e= xfmrCtgCurrReal1Var(i1,i2,j,k) * busCtgVoltMagVar(i1,k) * cos(busCtgVoltAngVar(i1,k))
   +  xfmrCtgCurrImag1Var(i1,i2,j,k) * busCtgVoltMagVar(i1,k) * sin(busCtgVoltAngVar(i1,k));

xfmrCtgPowImag1Def(i1,i2,j,k)$(XfmrCtgActive(i1,i2,j,k) and CtgActive(k))..
      xfmrCtgPowImag1Var(i1,i2,j,k)
  =e= xfmrCtgCurrReal1Var(i1,i2,j,k) * busCtgVoltMagVar(i1,k) * sin(busCtgVoltAngVar(i1,k))
   -  xfmrCtgCurrImag1Var(i1,i2,j,k) * busCtgVoltMagVar(i1,k) * cos(busCtgVoltAngVar(i1,k));

xfmrCtgPowReal2Def(i1,i2,j,k)$(XfmrCtgActive(i1,i2,j,k) and CtgActive(k))..
      xfmrCtgPowReal2Var(i1,i2,j,k)
  =e= xfmrCtgCurrReal2Var(i1,i2,j,k) * busCtgVoltMagVar(i2,k) * cos(busCtgVoltAngVar(i2,k))
   +  xfmrCtgCurrImag2Var(i1,i2,j,k) * busCtgVoltMagVar(i2,k) * sin(busCtgVoltAngVar(i2,k));

xfmrCtgPowImag2Def(i1,i2,j,k)$(XfmrCtgActive(i1,i2,j,k) and CtgActive(k))..
      xfmrCtgPowImag2Var(i1,i2,j,k)
  =e= xfmrCtgCurrReal2Var(i1,i2,j,k) * busCtgVoltMagVar(i2,k) * sin(busCtgVoltAngVar(i2,k))
   -  xfmrCtgCurrImag2Var(i1,i2,j,k) * busCtgVoltMagVar(i2,k) * cos(busCtgVoltAngVar(i2,k));

swshCtgPowImagDef(i,k)$(SwshActive(i) and CtgActive(k))..
      swshCtgPowImagVar(i,k)
  =e= - swshCtgAdmImagVar(i,k) * sqr(busCtgVoltMagVar(i,k));

genCtgPowRealMaint(i,j,k)$(GenCtgActive(i,j,k) and not GenCtgParticipating(i,j,k) and CtgActive(k))..
      genCtgPowRealVar(i,j,k)
  =e= genPowRealVar(i,j);

genCtgPowRealDroop(i,j,k)$(GenCtgParticipating(i,j,k) and CtgActive(k))..
      genCtgPowRealVar(i,j,k)
   -  genPowRealVar(i,j)
  =e= genPartFact(i,j) * sum(i1$genArea(i,j,i1),areaCtgPowRealChangeVar(i1,k))
   +  genCtgPowRealOverVar(i,j,k)
   -  genCtgPowRealUnderVar(i,j,k);

genCtgPowRealMaxBound(i,j,k)$(GenCtgParticipating(i,j,k) and CtgActive(k))..
      GenPMax(i,j)
   -  genCtgPowRealVar(i,j,k)
  =g= 0;

genCtgPowRealMinBound(i,j,k)$(GenCtgParticipating(i,j,k) and CtgActive(k))..
      genCtgPowRealVar(i,j,k)
   -  GenPMin(i,j)
  =g= 0;

busCtgVoltMagMaint(i,k)$((not VoltMaintQuadSoft) and BusCtgVoltMagMaintDom(i,k) and CtgActive(k))..
      busCtgVoltMagVar(i,k)
   -  busVoltMagVar(i)
  =e= busCtgVoltMagOverVar(i,k)
   -  busCtgVoltMagUnderVar(i,k);

busCtgPowImagMaxBound(i,k)$(BusCtgVoltMagMaintDom(i,k) and CtgActive(k))..
      sum((i1,j1)$GenBusCtgVoltMagMaintDom(i1,j1,i,k),
      	  GenQMax(i1,j1) - genCtgPowImagVar(i1,j1,k))
  =g= 0;

busCtgPowImagMinBound(i,k)$(BusCtgVoltMagMaintDom(i,k) and CtgActive(k))..
      sum((i1,j1)$GenBusCtgVoltMagMaintDom(i1,j1,i,k),
      	  genCtgPowImagVar(i1,j1,k) - GenQMin(i1,j1))
  =g= 0;
