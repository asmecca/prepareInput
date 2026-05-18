#include <prepare.hpp>

int main()
{
  const Source eta{};
  
  const Momentum momP{1,1,1};
  const Momentum momM=-momP;
  const Gamma P5{.iGamma=5};
  const Smear smeP{.kappa=0.4,.n=40,.mom=momP};
  const Smear smeM{.kappa=0.4,.n=40,.mom=momM};
  const Phase phaseP{.mom=momP};
  const Phase phaseM{.mom=momM};
  const Phase phase2M{.mom=2*momM};
  
  const Prop prop{.kappa=0.133,.mass=0.02,.r=1,.charge=0.0,.residue=1e-10};
  
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
  
  return 0;
}
