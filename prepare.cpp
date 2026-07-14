#include "prepare.hpp"

const Gamma P5{.iGamma=5};

namespace tests
{
  const Momentum momP{1,1,1};
  const Momentum momM=-momP;
  const Smear smeP{.kappa=0.4,.n=80,.mom=momP};
  const Smear smeM{.kappa=0.4,.n=80,.mom=momM};
  const Phase phaseP{.mom=momP};
  const Phase phaseM{.mom=momM};
  const Phase phase2M{.mom=2*momM};
  const Phase phase2P{.mom=2*momP};
  
  const Prop prop{.kappa=0.1394267,.mass=0.00066690,.r=1,.charge=0.0,.residue=1e-20};
}

void dir()
{
  using namespace tests;
  
  const Source eta(0);
  
  Run run;
  
  Contracter& dir=run.getTracer("dir");
  dir.addGammas(5,5);
  
  Line bw(phaseM*smeP*prop*smeP*phaseM*eta,"bw");
  Line fw(phaseP*smeM*prop*smeM*phaseP*eta,"fw");
  dir.tr(bw,fw);
  
  run.compile();
  
  cout<<run.printSources()<<endl;
  
  run.debugFree=true;
  run.debugContr=true;
  cout<<run.printLines()<<endl;
  
  cout<<run.printContr()<<endl;
}

void dir2()
{
  using namespace tests;
  
  const Source eta(0);
  
  Run run;
  
  const Oper O=smeM*phase2P*P5*smeP;
  
  Line bw(prop*eta,"bw");
  Line fw1(O*prop*O.dag()*eta,"fw1");
  Line fw2(O.dag()*prop*O.dag()*eta,"fw2");
  
  Contracter& dir=run.getTracer("dir");
  dir.addGammas(0,0);
  dir.tr(bw,fw1);
  dir.tr(bw,fw2);
  dir.tr(bw,bw);
  
  run.compile();
  
  cout<<run.printSources()<<endl;
  
  run.debugFree=true;
  run.debugContr=true;
  cout<<run.printLines()<<endl;
  
  cout<<run.printContr()<<endl;
}

void dir3()
{
  using namespace tests;
  
  const Source eta(0);
  
  Run run;
  
  const Oper OP=smeM*phaseP;
  const Oper OM=smeP*phaseM;
  
  Line bw1(OM.dag()*prop*OM*eta,"bw1");
  Line bw2(OP.dag()*prop*OM*eta,"bw2");
  Line fw1(OP.dag()*prop*OP*eta,"fw1");
  Line fw2(OM.dag()*prop*OP*eta,"fw2");
  
  Contracter& dir=run.getTracer("dir");
  dir.addGammas(5,5);
  dir.tr(bw1,fw1);
  dir.tr(bw2,fw2);
  dir.tr(bw1,bw1);
  dir.tr(fw1,fw1);
  
  run.compile();
  
  cout<<run.printSources()<<endl;
  
  run.debugFree=true;
  run.debugContr=true;
  cout<<run.printLines()<<endl;
  
  cout<<run.printContr()<<endl;
}

void tri()
{
  using namespace tests;
  
  const Source eta(0);
  
  Run run;
  
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

void print(Run &run)
{
  // run.compile();
  
  cout<<run.printSources()<<endl;
  
  cout<<run.printLines()<<endl;
  
  // cout<<describe(OP)<<endl;
  // cout<<describe(OM)<<endl;
  
  
  cout<<run.printContr()<<endl;
}

void testBox()
{
  using namespace tests;
  
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
      return OM*DeltaT{.t=t}*prop*OP*DeltaT{.t=t}*prop.dag()*eta;
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
  
  boxA.addGammas(0,0);
  
  boxB.tr(Line(bw,"bw"),
	  Line(fwB,"fwB"));
  
  boxB.addGammas(0,0);
 
  auto& boxC=runBox.getTracer("boxC");
  boxC.tr(Line(eta,"So"),
	  Line(prop*DeltaT{.t=10}*OP.dag()*prop.dag()*DeltaT{.t=10}*OM.dag()*prop*OM*DeltaT{.t=0}*prop.dag()*OP*eta,"piece"));
  boxC.addGammas(0,0);
  
  print(runBox);
}

void comboA()
{
  using namespace tests;
  
  Source u(0);
  Source d(5);
  
  const Oper OP=P5;
  
  Run runCombo;
  
  auto& combo=runCombo.getTracer("Op");
  combo.addGammas(0,0);
  combo.tr(Line(u,"u"),
	   Line(OP*prop*d,"Op_d"));
  combo.tr(Line(d,"d"),
	   Line(OP.dag()*prop.dag()*u,"Op_u"));
  
  auto& base=runCombo.getTracer("Base");
  base.addGammas(5,5);
  base.tr(Line(prop*d,"straight_d"),
	  Line(prop*d,"straight_d"));
  
  print(runCombo);
}

void comboB()
{
  using namespace tests;
  
  Source u(0);
  Source d(5);
  
  // const Phase phase0P{.mom={0,0,0}};
  const Oper OP=phase2P*P5;
  
  Run runCombo;
  
  auto& combo=runCombo.getTracer("Op");
  combo.addGammas(0,0);
  combo.tr(Line(d,"d"),
	   Line(OP*prop*u,"Op_prop_u"));
  combo.tr(Line(u,"u"),
	   Line(OP.dag()*prop.dag()*d,"Opdag_propdag_d"));
  
  auto& base=runCombo.getTracer("Base");
  base.addGammas(0,0);
  base.tr(Line(OP*prop*d,"Op_prop_d"),
	  Line(prop*OP*d,"prop_OP_d"));
  
  print(runCombo);
}

void localBox()
{
  using Entry=
    std::tuple<const char*,Momentum>;
  
  std::vector<std::vector<Entry>> a{{std::tuple
				     {"PZ",Momentum{0,0,2}},
				     {"MZ",{0,0,-2}}},
				    {{"0P11",{0,2,2}},
				     {"0M11",{0,-2,-2}},
				     {"P101",{2,0,2}},
				     {"M101",{-2,0,-2}},
				     {"0M1P1",{0,-2,2}},
				     {"0P1M1",{0,2,-2}},
				     {"M10P1",{-2,0,2}},
				     {"P10M1",{2,0,-2}}},
				    {{"P111",{2,2,2}},
				     {"M111",{-2,-2,-2}},
				     {"M1P11",{-2,2,2}},
				     {"P1M11",{2,-2,-2}},
				     {"P1M1P1",{2,-2,2}},
				     {"M1P1M1",{-2,2,-2}},
				     {"M11P1",{-2,-2,2}},
				     {"P11M1",{2,2,-2}}},
				    {{"PZ2",{0,0,4}},
				     {"MZ2",{0,0,-4}}},
				    {{"P012",{0,2,4}},
				     {"M012",{0,-2,-4}}}};
  
  const Source eta(0);
  const Prop prop{.kappa=0.1394267,.mass=0.00066690,.r=-1,.charge=0.0,.residue=1e-20};
  
  Run run;
  
  auto& box=run.getTracer("box");
  box.addGammas(5,5);
  
  auto& tri=run.getTracer("tri");
  for(int ig=1;ig<=3;ig++)
    tri.addGammas(ig,5);
  
  auto& triNaz=run.getTracer("triNaz");
  for(int ig=1;ig<=3;ig++)
    triNaz.addGammas(ig,5);
  
  doNotMergeLinComb=true;
  
  for(int iSo=0;iSo<5;iSo++)
    for(int iRotSo=0;const auto& [nSo,momSo] : a[iSo])
      {
	Phase phSo{.mom=momSo};
	const Line bwLine(prop.dag()*phSo.dag()*P5*DeltaT{.t=0}*prop*phSo*eta,std::format("bw{}",nSo));
	
	tri.tr(bwLine,
	       Line(prop*eta,"-"));
	
	for(int iSi=0;iSi<5;iSi++)
	  {
	    Phase phSi{.mom=std::get<Momentum>(a[iSi].front())};
	    auto getT=
	      [&prop,
	       &eta,
	       &phSi](const size_t t)
	      {
		return prop.dag()*DeltaT{.t=t}*phSi*P5*prop*eta;
	      };
	    
	    const size_t t0=5;
	    if(iSo==0 and iRotSo==0)
	      for(size_t t=t0;t<26;t++)
		triNaz.tr(Line(prop*phSi*eta,std::format("nazBwSi{}",iSi)), ///Correct, even if on so
			  Line(getT(t),std::format("nazFwT{}Si{}",t,iSi)));
	    
	    Line cumul=DeltaT{.t=t0}*getT(t0);
	    for(size_t t=t0+1;t<26;t++)
	      cumul=cumul+DeltaT{.t=t}*getT(t);
	    
	    box.tr(bwLine,
		   Line(phSi.dag()*cumul,std::format("fwP{}",iSi)));
	  }
	
	iRotSo++;
      }
  
  run.compile();
  
  print(run);
}

void smeBox()
{
  using Entry=
    std::tuple<const char*,Momentum>;
  
  std::vector<std::vector<Entry>> a{{std::tuple
				     {"PZ",Momentum{0,0,2}},
				     {"MZ",{0,0,-2}}},
				    {{"0P11",{0,2,2}},
				     {"0M11",{0,-2,-2}},
				     {"P101",{2,0,2}},
				     {"M101",{-2,0,-2}},
				     {"0M1P1",{0,-2,2}},
				     {"0P1M1",{0,2,-2}},
				     {"M10P1",{-2,0,2}},
				     {"P10M1",{2,0,-2}}},
				    {{"P111",{2,2,2}},
				     {"M111",{-2,-2,-2}},
				     {"M1P11",{-2,2,2}},
				     {"P1M11",{2,-2,-2}},
				     {"P1M1P1",{2,-2,2}},
				     {"M1P1M1",{-2,2,-2}},
				     {"M11P1",{-2,-2,2}},
				     {"P11M1",{2,2,-2}}},
				    {{"PZ2",{0,0,4}},
				     {"MZ2",{0,0,-4}}},
				    {{"P012",{0,2,4}},
				     {"M012",{0,-2,-4}}}};
  
  const Source eta(0);
  const Prop propSameTime{.kappa=0.1394267,.mass=0.00066690,.r=+1,.charge=0.0,.residue=1e-20};
  const Prop propDiffTime{.kappa=0.1394267,.mass=0.00066690,.r=-1,.charge=0.0,.residue=1e-20};
  
  Run run;
  
  auto& box=run.getTracer("box");
  box.addGammas(5,5);
  
  auto& tri=run.getTracer("tri");
  for(int ig=1;ig<=3;ig++)
    tri.addGammas(ig,5);
  
  auto& triNaz=run.getTracer("triNaz");
  for(int ig=1;ig<=3;ig++)
    triNaz.addGammas(ig,5);
  
  doNotMergeLinComb=true;
  
  auto getOp=
    [](const Momentum& mom)
    {
      const Phase phPi1{.mom=mom};
      
      const Smear smPi1q{.kappa=0.4,.n=80,.mom=mom/2.0};
      const Smear smPi1qBar{.kappa=0.4,.n=80,.mom=-mom/2.0};
      
      return smPi1qBar*phPi1*smPi1q;
    };
  
  for(int iSo=0;iSo<5;iSo++)
    for(int iRotSo=0;const auto& [nSo,momSo] : a[iSo])
      if(iRotSo<2)
	{
	const Oper pi1=getOp(momSo);
	const Oper pi2=pi1.dag();
	const Line bwLine(propDiffTime.dag()*pi2*P5*DeltaT{.t=0}*propSameTime.dag()*pi1*eta,std::format("bw{}",nSo));
	
	tri.tr(bwLine,
	       Line(propDiffTime*eta,"sm_prop"));
	
	for(int iSi=0;iSi<5;iSi++)
	  {
	    const Oper pi3=getOp(std::get<Momentum>(a[iSi].front()));
	    const Oper pi4=pi3.dag();
	    auto getT=
	      [&propSameTime,
	       &propDiffTime,
	       &eta,
	       &pi3](const size_t t)
	      {
		return propSameTime*DeltaT{.t=t}*pi3*P5*propDiffTime*eta;
	      };
	    
	    const size_t t0=5;
	    const size_t tf=26;
	    if(iSo==0 and iRotSo==0)
	      for(size_t t=t0;t<tf;t++)
		triNaz.tr(Line(propDiffTime*pi3*eta,std::format("nazBwSi{}",iSi)), ///Correct, even if on so
			  Line(getT(t),std::format("nazFwT{}Si{}",t,iSi)));
	    
	    Line cumul=DeltaT{.t=t0}*getT(t0);
	    for(size_t t=t0+1;t<tf;t++)
	      cumul=cumul+DeltaT{.t=t}*getT(t);
	    
	    box.tr(bwLine,
		   Line(pi4*cumul,std::format("fwP{}",iSi)));
	  }
	
	iRotSo++;
      }
  
  run.compile();
  
  // run.debugContr=true;
  // run.debugFree=true;
  print(run);
}

void smeDir()
{
  using Entry=
    std::tuple<const char*,Momentum>;
  
  std::vector<std::vector<Entry>> a{{std::tuple
				     {"PZ",Momentum{0,0,2}},
				     {"MZ",{0,0,-2}}},
				    {{"0P11",{0,2,2}},
				     {"0M11",{0,-2,-2}}},
				    {{"P111",{2,2,2}},
				     {"M111",{-2,-2,-2}}},
				    {{"PZ2",{0,0,4}},
				     {"MZ2",{0,0,-4}}},
				    {{"P012",{0,2,4}},
				     {"M012",{0,-2,-4}}}};
  
  const Source eta(0);
  const Prop prop1{.kappa=0.1394267,.mass=0.00066690,.r=+1,.charge=0.0,.residue=1e-20};
  const Prop prop0{.kappa=0.1394267,.mass=0.00066690,.r=-1,.charge=0.0,.residue=1e-20};
 
  Run run;
  
  auto& dir=run.getTracer("dir");
  dir.addGammas(5,5);

  auto& current=run.getTracer("current");
  current.addGammas(1,1);
  current.addGammas(2,2);
  current.addGammas(3,3);
    
  for(int iSo=0;iSo<5;iSo++)
    {
      const auto& [nSo, momSo] = a[iSo][0];
      for(int iSi=0;iSi<5;iSi++)
	for(int iRotSi=0;const auto& [nSi,momSi] : a[iSi])
	  {	
	    const Phase phSoM{.mom=-momSo/2.0};
	    const Phase phSoP{.mom=momSo/2.0};
	    const Smear smSoP{.kappa=0.4,.n=80,.mom=momSo/2.0};
	    const Smear smSoM{.kappa=0.4,.n=80,.mom=-momSo/2.0};
	    const Phase phSiM{.mom=-momSi/2.0};
	    const Phase phSiP{.mom=momSi/2.0};
	    const Smear smSiP{.kappa=0.4,.n=80,.mom=momSi/2.0};
	    const Smear smSiM{.kappa=0.4,.n=80,.mom=-momSi/2.0};
	    const Line bwLine(phSiM*smSiM*prop0*smSoP*phSoM*eta,std::format("bw{}_{}",nSi,nSo));
	    const Line fwLine(phSiP*smSiP*prop0*smSoM*phSoP*eta,std::format("fw{}_{}",nSi,nSo));
	    dir.tr(bwLine,fwLine);
	    iRotSi++;
	  }
    }
  const Line bwLineJ(prop1*eta,"propR1");
  const Line fwLineJ(prop0*eta,"propR0");
  current.tr(bwLineJ,fwLineJ);
  run.compile();
  
  // run.debugContr=true;
  // run.debugFree=true;
  print(run);
}

int main()
{
  // comboA();
  // comboB();
  
  //localBox();
  //smeBox();
  smeDir();
  //box();
  //dir3();
  
  return 0;
}
