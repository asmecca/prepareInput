#include <algorithm>
#include <array>
#include <cstdarg>
#include <cstddef>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <iostream>
#include <execinfo.h>
#include <variant>
#include <vector>

using std::cout;
using std::cerr;
using std::endl;

/// Overload pattern
template <typename...Ts>
struct Overload :
  Ts...
{
  using Ts::operator()...;
};

/// Deducing guide for overload pattern
template <typename...Ts>
Overload(Ts...) -> Overload<Ts...>;

/// Cat a list of types to a variant
template <typename V,
	  typename...T>
struct _VariantCat;

template <typename...V,
	  typename...T>
struct _VariantCat<std::variant<V...>,T...>
{
  using type=std::variant<V...,T...>;
};

template <typename V,
	  typename...T>
using VariantCat=typename _VariantCat<V,T...>::type;

/// Write the list of called routines
void print_backtrace_list()
{
  void *callstack[128];
  int frames=backtrace(callstack,128);
  char **strs=backtrace_symbols(callstack,frames);
  
  //only master rank, but not master thread
  cerr<<"Backtracing..."<<endl;
  for(int i=0;i<frames;i++) cerr<<strs[i]<<endl;
  
  free(strs);
}

#define CRASH(...) \
  internal_crash(__LINE__,__FILE__,__PRETTY_FUNCTION__,__VA_ARGS__)

/// Crashes emitting the message
__attribute__ ((format (printf,4,5)))
void internal_crash(const int& line,
		    const char* file,
		    const char* func,
		    const char* format,...)
{
  char buffer[1024];
  va_list args;
  
  va_start(args,format);
  vsprintf(buffer,format,args);
  va_end(args);
  
  cerr<<"ERROR in function "<<func<<" at line "<<line<<" of file "<<file<<": \""<<buffer<<"\""<<endl;
  print_backtrace_list();
  exit(1);
}

/////////////////////////////////////////////////////////////////

/// Holds three components of the momentum
struct Momentum :
  std::array<double,3>
{
};

std::ostream& operator<<(std::ostream& os,
			 const Momentum& mom)
{
  return os<<"{"<<mom[0]<<","<<mom[1]<<","<<mom[2]<<"}";
}

/// Negate a momentum
Momentum operator-(const Momentum& y)
{
  return {-y[0],-y[1],-y[2]};
}

/// Product of a scalar with a momentum
Momentum operator*(const double& x,
	      const Momentum& y)
{
  return {y[0]*x,y[1]*x,y[2]*x};
}

/// Product of a momentum with a scalar
Momentum operator*(const Momentum& x,
		   const double& y)
{
  return y*x;
}

/////////////////////////////////////////////////////////////////

struct Smear;

struct Phase;

struct Gamma;

struct Oper;

/////////////////////////////////////////////////////////////////

struct Smear
{
  double kappa{};
  
  double n{};
  
  Momentum mom{};
  
  Smear dag() const
  {
    return *this;
  }
  
  std::string describe() const
  {
    return (std::ostringstream()<<"smear(kappa="<<kappa<<",n="<<n<<",mom="<<mom<<")").str();
  }
  
  constexpr std::partial_ordering operator<=>(const Smear&) const=default;
};

struct Phase
{
  Momentum mom;
  
  std::string describe() const
  {
    return (std::ostringstream()<<"phase(mom="<<mom<<")").str();
  }
  
  Phase dag() const
  {
    Phase res{*this};
    
    res.mom=-mom;
    
    return res;
  }
  
  constexpr std::partial_ordering operator<=>(const Phase&) const=default;
};

struct Gamma
{
  size_t iGamma{};
  
  std::string describe() const
  {
    return (std::ostringstream()<<"gamma(iGamma="<<iGamma<<")").str();
  }
  
  Gamma dag() const
  {
    return *this;
  }
  
  constexpr std::partial_ordering operator<=>(const Gamma&) const=default;
};

struct DeltaT
{
  size_t t;
  
  std::string describe() const
  {
    return (std::ostringstream()<<"deltaT(t="<<t<<")").str();
  }
  
  DeltaT dag() const
  {
    return *this;
  }
  
  constexpr std::partial_ordering operator<=>(const DeltaT&) const=default;
};

struct Prop
{
  double kappa{};
  
  double mass{};
  
  int r{};
  
  double charge{};
  
  Momentum mom{};
  
  double residue{};
  
  std::string describe() const
  {
    return (std::ostringstream()<<"prop(kappa="<<kappa<<",mass="<<mass<<",r="<<r<<",charge="<<charge<<",mom="<<mom<<",residue="<<residue<<")").str();
  }
  
  Prop dag() const
  {
    Prop res=*this;
    
    res.r=-r;
    
    res.mom=-mom;
    
    return *this;
  }
  
  constexpr std::partial_ordering operator<=>(const Prop&) const=default;
};

/////////////////////////////////////////////////////////////////

using Pars=
  std::variant<Smear,Phase,Gamma,DeltaT,Prop>;

Pars dag(const Pars& pars)
{
  return
    std::visit([](const auto& o) -> Pars
      {
	return o.dag();
      },pars);
}

std::string describe(const Pars& pars)
{
  return
    std::visit([](const auto& o) -> std::string
      {
	return o.describe();
      },pars);
}

Oper operator*(const Oper& a,
	       const Oper& b);

using Command=
  std::tuple<size_t,Pars,int,size_t>;

struct Executable
{
  std::vector<Command> instructions;
  
  bool contains(const size_t& lhs) const
  {
    return std::find_if(instructions.begin(),instructions.end(),
			[&lhs](const Command& command)
			{
			  return lhs==std::get<0>(command);
			})!=instructions.end();
  }
  
  bool contains(const std::vector<size_t>& rhs) const
  {
    bool are{true};
    auto it=rhs.begin();
    while(are and it!=rhs.end())
      are&=contains(*it++);
    
    return are;
  }
};

struct Oper
{
  Pars pars;
  
  std::unique_ptr<Oper> rhs;
  
  // std::variant<Pars,std::vector<std::pair<double,Oper>>> pars;
  
  // template <typename T>
  // requires std::constructible_from<decltype(pars),T>
  // Oper(const T& args) :
  //   pars{args}
  // {
  // }
  
  // Oper asList() const
  // {
  //   if(std::holds_alternative<Pars>(pars))
  //     return std::vector{std::make_pair(1.0,*this)};
  //   else
  //     return *this;
  /// Dagger of an operation
  Oper dag() const
  {
    const Pars dagPars=
      ::dag(pars);
    
    if(rhs)
      return rhs->dag()*dagPars;
    else
      return dagPars;
  }
  
  Oper(const Oper& oth) :
    pars{oth.pars},
    rhs{oth.rhs?std::make_unique<Oper>(*oth.rhs):nullptr}
  {
  }
  
  Oper(const Pars& pars,
       const Oper& rhs) :
    pars{pars},
    rhs(std::make_unique<Oper>(rhs))
  {
  }
  
  template <typename T>
  requires std::constructible_from<decltype(pars),T>
  Oper(const T& pars) :
     pars{pars}
  {
  }
  
  // std::string describe() const
  // {
  //   std::string t=::describe(pars);
    
  //   if(rhs)
  //     t+="*"+rhs->describe();
    
  //   return t;
  // }
  
  // size_t hash() const
  // {
  //   return std::hash<std::string>{}(describe());
  // }
  
  constexpr bool operator==(const Oper& other) const
  {
    return pars==other.pars and
      (((not rhs) and (not other.rhs))
       or (rhs and other.rhs and rhs==other.rhs));
  }
  
  constexpr std::partial_ordering operator<=>(const Oper& other) const
  {
    if(auto cmp=pars<=>other.pars;cmp!=0)
      return cmp;
    
    if ((not rhs) and (not other.rhs))
      return std::partial_ordering::equivalent;
    
    if(not rhs)
      return std::partial_ordering::less;
    
    if(not other.rhs)
      return std::partial_ordering::greater;
    
    return *rhs<=>*other.rhs;
  }
};

Oper operator*(const Oper& a,
	       const Oper& b)
{
  Oper c(a);
  
  Oper* d=&c;
  while(d->rhs)
    d=&*d->rhs;
  
  d->rhs=std::make_unique<Oper>(b);
  
  return c;
}

struct Source
{
  inline static size_t glbCount{};
  
  size_t count;
  
  std::string describe() const
  {
    return (std::ostringstream()<<"source(count="<<count<<")").str();
  }
  
  Source()
  {
    count=glbCount++;
  }
  
  constexpr std::partial_ordering operator<=>(const Source&) const=default;
};


struct Line;

struct Extender
{
  Pars pars;
  
  std::vector<std::pair<double,Line>> rhs;
  
  int tSelect{-1};
  
  std::string describe() const;
  
  bool isScalar() const;
  
  bool isCombiner() const ;
  
  bool isExtendingSource() const;
  
  bool isSelectingTime() const
  {
    return tSelect!=-1;
  }
  
  Line* maybeGetTrivialRhs();
  
  Extender* maybeGetNestedExtender();
  
  constexpr bool operator==(const Extender&) const;
  constexpr std::partial_ordering operator<=>(const Extender&) const;
};

struct Line
{
  std::variant<Extender,Source> data;
  
  bool isExtender() const
  {
    return std::holds_alternative<Extender>(data);
  }
  
  size_t maybeRemoveDeltaT()
  {
    return
      std::visit(Overload{[](Extender& extender) -> size_t
      {
	const DeltaT* dt=
	  std::get_if<DeltaT>(&extender.pars);
	
	size_t n=
	  dt!=nullptr;
	
	if(n)
	  {
	    extender.tSelect=dt->t;
	    extender.pars=Gamma{0};
	  }
	
	for(auto& [w,l] : extender.rhs)
	  n+=l.maybeRemoveDeltaT();
	
      return n;
      },
	    [](Source& source) -> size_t
	    {
	      return 0;
	    }},data);
  }
  
  void trySimplify()
  {
    maybeRemoveDeltaT();
    maybeRemoveUselessScalar();
  }
  
  size_t maybeRemoveUselessScalar();
  
  Line(const Pars& pars,
       const std::vector<std::pair<double,Line>>& rhs,
       const int tSelect=-1) :
    data{Extender{pars,rhs,tSelect}}
  {
    trySimplify();
  }
  
  Line(const Oper& oper,
       const Line& line) :
    data(Extender{oper.pars})
  {
    Extender& e=*std::get_if<Extender>(&data);
    
    if(oper.rhs)
      e.rhs.emplace_back(1.0,Line(*oper.rhs,line));
    else
      e.rhs.emplace_back(1.0,line);
    
    trySimplify();
  }
  
  constexpr std::partial_ordering operator<=>(const Line&) const=default;
  
  Line(const Line&)=default;
  
  template <typename T>
  requires std::constructible_from<decltype(data),T>
  Line(const T& data) :
    data{data}
  {
    trySimplify();
  }
  
  std::string describe() const
  {
    return
      std::visit([](const auto& e)
      {
	return e.describe();
      },data);
  }
};

constexpr bool Extender::operator==(const Extender&) const=default;

constexpr std::partial_ordering Extender::operator<=>(const Extender& other) const
{
  if(auto cmp=pars<=>other.pars;cmp!=0)
    return cmp;
  
  if(auto cmp=rhs<=>other.rhs;cmp!=0)
    return cmp;
  
  return tSelect<=>other.tSelect;
}

bool Extender::isExtendingSource() const
{
  return rhs.size() and not rhs.front().second.isExtender();
}

Line* Extender::maybeGetTrivialRhs()
{
  if(isCombiner())
    return nullptr;
  else
    return &rhs.front().second;
}

std::string Extender::describe() const
{
  std::ostringstream os;
  os<<::describe(pars);
  
  if(rhs.size()==0)
    CRASH("cannot describe a line extender which does not extend an existing line");
  else
    {
      if(tSelect!=-1)
	os<<"*"<<::describe(DeltaT{(size_t)tSelect});
      
      os<<"*";
      
      if(rhs.size()>1)
	os<<"(";
      
      for(size_t c{};const auto& [w,line] : rhs)
	{
	  if(w<0)
	    os<<w<<"*";
	  else
	    {
	      if(c)
		os<<"+";
	      
	      if(w!=1.0)
		os<<w<<"*";
	    }
	  
	  os<<line.describe();
	  
	  c++;
	}
      
      if(rhs.size()>1)
	os<<")";
    }
  
  return os.str();
}

std::string describe(const Line& line)
{
  return std::visit([](const auto& v)
  {
    return v.describe();
  },line.data);
}

size_t Line::maybeRemoveUselessScalar()
{
  // Check if is an extender
  Extender* e=
    std::get_if<Extender>(&this->data);
  
  // Nested remove
  int nDeep{};
  if(e)
    for(auto& [w,l] : e->rhs)
    nDeep+=l.maybeRemoveUselessScalar();
  
  // maybe remove this if 1) is scalar 2) is trivially nesting 3) extends an extender 4) has no time selection
  
  Line* nestedRhs=
    e?e->maybeGetTrivialRhs():nullptr;
  
  if(e and nestedRhs)
    {
      if(e->isScalar() and e->tSelect==-1)
	{
	  *this=*nestedRhs;
	  
	  return nDeep+1;
	}
      
      Extender* nestedE=
	std::get_if<Extender>(&nestedRhs->data);
      
      if(nestedE->isScalar() and e->tSelect==-1)
	{
	  e->tSelect=nestedE->tSelect;
	  e->rhs=std::move(nestedE->rhs);
	  
	  return nDeep+1;
	}
    }
  
  return nDeep;
}

Line operator*(const Oper& oper,
	       const Line& line)
{
  return {oper,line};
}

Line operator*(const double& c,
	       const Line& line)
{
  return std::visit(Overload{[](const Source& source)->Line
  {
    return source;
  },
	[&c](const Extender& le)->Line
	{
	  Extender res{le};
	  
	  for(auto& [w,l] : res.rhs)
	    w*=c;
	  
	  return res;
	}},line.data);
}

Line operator*(const Line& line,
	       const double& c)
{
  return c*line;
}

Line operator+(const Line& a,
	       const Line& b)
{
  if(const Extender *ae=std::get_if<Extender>(&a.data),*be=std::get_if<Extender>(&b.data);ae and be and ae->pars==be->pars)
    {
      Extender res=*ae;
      
      res.rhs.insert(res.rhs.end(),be->rhs.begin(),be->rhs.end());
      
      return res;
    }
  else
    return Extender(Gamma{0},std::vector{std::make_pair(1.0,a),{1.0,b}});
}

struct Contracter
{
  std::string name;
  
  std::set<std::pair<int,int>> gammaList;
  
  std::set<std::pair<Line,Line>> lines;
  
  Contracter(const std::string& name) :
    name(name)
  {
  }
  
  void addGammas(const int& sink,
	   const int& sour)
  {
    gammaList.emplace(sink,sour);
  }
  
  void operator()(const Line& bw,
		  const Line& fw)
  {
    lines.emplace(bw,fw);
  }
};

struct Compiler
{
  std::map<Line,std::set<std::variant<Line,std::pair<Line,Line>>>> dependencies;
  
  void insert(const Line& l,
	      const std::variant<Line,std::pair<Line,Line>>& dep)
  {
    dependencies[l].insert(dep);
    
    std::visit(Overload{[](const Source& source)
    {
    },
	  [&l,
	   this](const Extender& extender)
	  {
	    for(const auto& [w,li] : extender.rhs)
	      insert(li,l);
	  }},l.data);
  }
  
  void compile(const Contracter& contr)
  {
    for(const auto& [bw,fw] : contr.lines)
      for(const Line& l : {bw,fw})
	insert(l,std::make_pair(bw,fw));
  }
};

// void compile(UnorderedInstructions& out,
// 	     const Line& line)
// {
  
//   for(const Oper& oper : operations)
//     oper.compile(res);
  
//   return res;
// }

// struct Executable
// {
//   std::map<std::pair<size_t,size_t>,std::pair<size_t,size_t>> contractions;
  
//   void contract(const Line& a,
// 		const size_t& gSi,
// 		const Line& b,
// 		const size_t& gSo)
//   {
    
//   }
// };

/////////////////////////////////////////////////////////////////

int main()
{
  const Momentum mom{0,1,0};
  
  const Smear sme{.kappa=0.4,.n=40,.mom=mom};
  const Phase phase{.mom=2*mom};
  const DeltaT deltaT{.t=3};
  
  Prop prop{.kappa=0.133,.mass=0.02,.r=1,.charge=0.0,.residue=0.0011};
  
  const Oper op1=prop*sme*deltaT*phase*sme;
  const Oper op2=phase*sme*sme;
  
  std::vector<Oper> operations{op1,op2};
  // UnorderedInstructions unorderedInstructions=compile(operations);
  
  // const auto describe=
  //   [](const Command& command)
  //   {
  //     const auto& [lhs,pars,source]=command;
      
  //     std::ostringstream os;
  //     os<<lhs<<" = " <<::describe(pars)<<" * "<<source;
      
  //     return os.str();
  //   };
  
  // cout<<"Unordered instructions ("<<unorderedInstructions.instructions.size()<<" instr):"<<endl;
  // for(const Command& command : unorderedInstructions.instructions)
  //   cout<<" "<<describe(command)<<endl;
  // cout<<endl;
  
  Source eta{};
  Line a=op2*(op1*eta+eta);
  Line b=op2*eta;
  
  const Oper* o=&op1;
  while(o)
    {
      cout<<describe(o->pars)<<endl;
      o=&*o->rhs;
    }
  
  
  cout<<describe(a)<<endl;
  
  Contracter box("box");
  box.addGammas(5,5);
  box(a,b);
  
  Compiler compiler;
  compiler.compile(box);
  
  cout<<endl;
  for(const auto& [l,deps] : compiler.dependencies)
    {
      cout<<describe(l)<<endl;
      for(const auto& d : deps)
	cout<<" "<<std::visit(Overload{[](const Line& line)
	{
	  return describe(line);
	},
	      [](const std::pair<Line,Line>& contr)->std::string
	      {
		return "contr";
	      }},d)<<endl;
      cout<<endl;
    }
  
  Executable executable;
  
  // std::map<size_t,size_t> dependenciesMultiplicity;
  // for(const auto& [lhs,pars,sources] : unorderedInstructions.instructions)
  //   for(const size_t& source : sources)
  //    dependenciesMultiplicity[source]++;
  
  // while(executable.instructions.size()!=unorderedInstructions.instructions.size())
  //   {
  //     cout<<"NResidual instructions: "<<unorderedInstructions.instructions.size()-executable.instructions.size()<<endl;
      
  //     std::vector<std::pair<size_t,size_t>> eligibleCommands;
      
  //     for(size_t i=0;const auto& [lhs,pars,sources] : unorderedInstructions.instructions)
  // 	{
  // 	  if(not executable.contains(lhs) and executable.contains(sources))
  // 	    eligibleCommands.emplace_back(i,dependenciesMultiplicity[lhs]);
	  
  // 	  i++;
  // 	}
      
  //     std::sort(eligibleCommands.begin(),eligibleCommands.end(),[](const std::pair<size_t,size_t>& a,
  // 								   const std::pair<size_t,size_t>& b)
  //     {
  // 	return a.second<b.second;
  //     });
      
  //     cout<<"Eligible commands:"<<endl;
  //     for(const auto& [iLhs,nDep] : eligibleCommands)
  // 	cout<<" "<<describe(unorderedInstructions.instructions[iLhs])<<" holds "<<nDep<<" other instructions"<<endl;
  //     cout<<endl;
      
  //     if(eligibleCommands.empty())
  // 	CRASH("no eligible command");
      
  //     executable.instructions.push_back(unorderedInstructions.instructions[eligibleCommands.back().first]);
  //   }
  
  // cout<<"Final program:"<<endl;
  // for(const Command& command : executable.instructions)
  //   cout<<" "<<describe(command)<<endl;
  // cout<<endl;
  
  return 0;
}
