#include <prepare.hpp>

const Momentum momP{1,1,1};
const Momentum momM=-momP;
const Gamma P5{.iGamma=5};
const Smear smeP{.kappa=0.4,.n=80,.mom=momP};
const Smear smeM{.kappa=0.4,.n=80,.mom=momM};
const Phase phaseP{.mom=momP};
const Phase phaseM{.mom=momM};
const Phase phase2M{.mom=2*momM};
const Phase phase2P{.mom=2*momP};

const Prop prop{.kappa=0.1394267,.mass=0.00066690,.r=1,.charge=0.0,.residue=1e-20};

void tri()
{
  const Source eta(0);
  
  Run run;
  
  Contracter& dir=run.getTracer("dir");
  dir.addGammas(5,5);
  
  Line bw(phaseM*smeP*prop*smeP*phaseM*eta,"bw");
  Line fw(phaseP*smeM*prop*smeM*phaseP*eta,"fw");
  dir.tr(bw,fw);
  
  Contracter& tri=run.getTracer("tri");
  tri.addGammas(1,5);
  for(size_t t=0;t<5;t++)
    tri.tr(Line(prop*smeP*phaseM*eta,"bw2"),
	   Line(prop.dag()*smeP*phase2M*P5*smeM*DeltaT{.t=t}*prop*smeM*phaseP*eta,std::format("fw2_{}",t)));
  
  run.compile();
  
  cout<<run.printSources()<<endl;
  
  run.debugContr=true;
  cout<<run.printLines()<<endl;
  
  cout<<run.printContr()<<endl;
}

void box()
{
  Run runBox;
  
  /// Original source
  const Source eta(0);
  
  /// Defines the interpolating operator for the single pion
  const Oper OP=smeM*phase2P*P5*smeP;
  const Oper OM=OP.dag();
  
  /// Defines the forward line, which contains the _| part of the corr
  const auto fwA=prop*OM*DeltaT{.t=0}*prop.dag()*OP*eta;
  const auto fwB=prop*OP*DeltaT{.t=0}*prop.dag()*OM*eta;
  
  /// Gets the backward line for a given time
  auto getBw=
    [&](const size_t t)
    {
      return OM*DeltaT{.t=t}*prop*OP*DeltaT{.t=t}*prop*eta;
    };
  
  /// Gets the original line
  Line bw=getBw(0);
  
  doNotMergeLinComb=true;
  for(size_t t=1;t<25;t++)
    bw=Line(bw+getBw(t),std::format("bwD{}",t));
  
  auto& boxA=runBox.getTracer("boxA");
  auto& boxB=runBox.getTracer("boxB");
  
  boxA.tr(Line(bw,"bw"),
	  Line(fwA,"fwA"));
  
  boxB.tr(Line(bw,"bw"),
	  Line(fwB,"fwB"));
  
  runBox.compile();
  
  cout<<runBox.printLines()<<endl;
  
  cout<<describe(OP)<<endl;
  cout<<describe(OM)<<endl;
}

int main()
{
  box();
  
  return 0;
}
