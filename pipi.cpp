#include "prepare.hpp"

const Gamma P5{.iGamma=5};

const std::vector<Momentum> momList{{0,0,1},{0,1,1},{1,1,1},{0,0,2},{0,1,2}};

const int nOp=momList.size();

void print(Run &run)
{
  // run.compile();
  
  cout<<run.printSources()<<endl;
  
  cout<<run.printLines()<<endl;
  
  // cout<<describe(OP)<<endl;
  // cout<<describe(OM)<<endl;
  
  
  cout<<run.printContr()<<endl;
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
  
  for(int iSo=0;iSo<nOp;iSo++)
    for(int iRotSo=0;const auto& [nSo,momSo] : a[iSo])
      {
	Phase phSo{.mom=momSo};
	const Line bwLine(prop.dag()*phSo.dag()*P5*DeltaT{.t=0}*prop*phSo*eta,std::format("bw{}",nSo));
	
	tri.tr(bwLine,
	       Line(prop*eta,"-"));
	
	for(int iSi=0;iSi<nOp;iSi++)
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

std::pair<std::set<Momentum>,std::set<Momentum>> getAllPerms(Momentum a)
{
  std::set<Momentum> source;
  std::set<Momentum> sink;
  
  do
    {
      Momentum c{a};
      std::sort(c.begin(),c.begin()+2);
      if(c[2])
	sink.insert(c);
      
      for(int i=0;i<8;i++)
	{
	  Momentum b;
	  
	  for(int mu=0;mu<3;mu++)
	    b[mu]=a[mu]*(((i>>mu)&1)*2-1);
	  
	  if(b[2])
	    source.insert(b);
	}
    }
  while(std::next_permutation(a.begin(),a.end()));
  
  return {source,sink};
}

void smeBox()
{
  const Source eta(0);
  const Prop propSameTime{.kappa=0.1394267,.mass=0.00066690,.r=+1,.charge=0.0,.residue=1e-20};
  const Prop propDiffTime{.kappa=0.1394267,.mass=0.00066690,.r=-1,.charge=0.0,.residue=1e-20};
 
  Run run;
  
  auto& box=run.getTracer("box");
  box.addGammas(5,5);
  
  auto& tri=run.getTracer("tri");
  auto& triSme=run.getTracer("triSme");
  for(auto& t : {&tri,&triSme})
    for(int ig=1;ig<=3;ig++)
      t->addGammas(ig,5);
  
  auto& triNaz=run.getTracer("triNaz");
  for(int ig=1;ig<=3;ig++)
    triNaz.addGammas(ig,5);
  
  doNotMergeLinComb=true;
  
  const size_t t0=5;
  const size_t tf=26;
  
  size_t totalCost{1};
  
  auto getOp=
    [](const Momentum& mom)
    {
      const Phase phPi1{.mom=mom};
      
      const Smear smPi1q{.kappa=0.4,.n=80,.mom=mom/2.0};
      const Smear smPi1qBar{.kappa=0.4,.n=80,.mom=-mom/2.0};
      
      return smPi1qBar*phPi1*smPi1q;
    };
  
  std::vector<Line> bwLine;
  std::vector<Line> fwLine;
  for(int iMom=0;iMom<nOp;iMom++)
    {
      auto getSourceLine=
	[&getOp,
	 &propSameTime,
	 &eta](const Momentum& momSo) -> Line
	{
	  const Oper pi1=getOp(momSo);
	  const Oper pi2=pi1.dag();
	  return {pi2*P5*DeltaT{.t=0}*propSameTime.dag()*pi1*eta};
	};
      
      const auto [source,sink]=
	getAllPerms(momList[iMom]);
      
      for(const auto& [s,tag] : {std::make_pair(source,"source"),{sink,"sink"}})
	{
	  cout<<tag<<" "<<s.size()<<" contributions"<<endl;
	  for(const Momentum& a : s)
	    cout<<a<<endl;
	}
      const size_t cost=
	source.size()+1+
	(tf-t0)*sink.size();
      cout<<"cost: "<<cost<<endl;
      totalCost+=cost;
      
      std::optional<Line> sameTimeLine;
      for(const Momentum& signedMom : source)
	{
	  const Line curContr=
	    signedMom[2]*getSourceLine(signedMom);
	  
	  if(sameTimeLine)
	    *sameTimeLine=*sameTimeLine+curContr;
	  else
	    sameTimeLine=curContr;
	}
      
      bwLine.emplace_back(propDiffTime*(*sameTimeLine),std::format("bw{}",iMom));
      
      /////////////////////////////////////////////////////////////////
      
      bool nazIncl{};
      
      auto getTcontr=
	[&propSameTime,
	 &propDiffTime,
	 &eta,
	 &getOp,
	 &iMom,
	 &triNaz,
	 &nazIncl,
	 &t0](const Momentum& momSi)
	{
	  const Oper pi1=getOp(momSi);
	  const Oper pi2=pi1.dag();
	      
	  auto getT=
	    [&](const size_t& t)
	    {
	      return propSameTime*DeltaT{.t=t}*pi1*P5*propDiffTime*eta;
	    };
	  
	  if(not nazIncl)
	    for(size_t t=t0;t<tf;t++)
	      triNaz.tr(Line(propDiffTime*pi1*eta,std::format("nazBwSi{}",iMom)),
			Line(getT(t),std::format("nazFwT{}Si{}",t,iMom)));
	  nazIncl=true;
	  
	  Line cumul=DeltaT{.t=t0}*getT(t0);
	  for(size_t t=t0+1;t<tf;t++)
	    cumul=cumul+DeltaT{.t=t}*getT(t);
	  
	  return pi2*cumul;
	};
      
      std::optional<Line> otherTimeLine;
      const int r=source.size()/sink.size();
      for(const Momentum& s : sink)
	{
	  const Line c=getTcontr(s)*(s[2]*r);
	  if(not otherTimeLine.has_value())
	    otherTimeLine=c;
	  else
	    otherTimeLine=*otherTimeLine+c;
	}
      
      fwLine.emplace_back(*otherTimeLine,std::format("fw{}",iMom));
    }
  
  for(const Line& b : bwLine)
    for(const Line& f : fwLine)
      box.tr(b,f);
  
  for(const Line& b : bwLine)
    tri.tr(b,
	   Line(propDiffTime*eta,"fw"));
  
  cout<<"Total cost: "<<totalCost<<endl;
  run.compile();
  
  // run.debugContr=true;
  // run.debugFree=true;
  print(run);
}

/*
void smeDir()
{
  const Source eta(0);
  const Prop prop1{.kappa=0.1394267,.mass=0.00066690,.r=+1,.charge=0.0,.residue=1e-20};
  const Prop prop0{.kappa=0.1394267,.mass=0.00066690,.r=-1,.charge=0.0,.residue=1e-20};
 
  Run run;
  debugPrepare=true;
  
  auto& dir=run.getTracer("dir");
  dir.addGammas(5,5);
  
  auto& current=run.getTracer("current");
  current.addGammas(1,1);
  current.addGammas(2,2);
  current.addGammas(3,3);
  
  for(int iSo=0;iSo<5;iSo++)
    {
      const Momentum& momSo=momList[iSo];
      for(int iSi=0;iSi<5;iSi++)
	loopOnAllPerm([&momSo,
		       &prop0,
		       &dir,
		       &eta,
		       &iSo](const Momentum& momSi)
	{
	  const Phase phSoM{.mom=-momSo/2.0};
	  const Phase phSoP{.mom=momSo/2.0};
	  const Smear smSoP{.kappa=0.4,.n=80,.mom=momSo/2.0};
	  const Smear smSoM{.kappa=0.4,.n=80,.mom=-momSo/2.0};
	  const Phase phSiM{.mom=-momSi/2.0};
	  const Phase phSiP{.mom=momSi/2.0};
	  const Smear smSiP{.kappa=0.4,.n=80,.mom=momSi/2.0};
	  const Smear smSiM{.kappa=0.4,.n=80,.mom=-momSi/2.0};
	  const Line bwLine(phSiM*smSiM*prop0*smSoP*phSoM*eta,std::format("bw{}_({},{},{})",iSo,momSi[0],momSi[1],momSi[2]));
	  const Line fwLine(phSiP*smSiP*prop0*smSoM*phSoP*eta,std::format("fw{}_({},{},{})",iSo,momSi[0],momSi[1],momSi[2]));
	  dir.tr(bwLine,fwLine);
	},momList[iSi]);
    }
  const Line bwLineJ(prop1*eta,"propR1");
  const Line fwLineJ(prop0*eta,"propR0");
  current.tr(bwLineJ,fwLineJ);
  run.compile();
  
  // run.debugContr=true;
  // run.debugFree=true;
  print(run);
}
*/

int main()
{
  // comboA();
  // comboB();
  
  //localBox();
  smeBox();
  //smeDir();
  //box();
  //dir3();
  
  return 0;
}
