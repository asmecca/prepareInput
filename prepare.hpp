#include <algorithm>
#include <array>
#include <cmath>
#include <cstdarg>
#include <cstddef>
#include <format>
#include <limits>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <iostream>
#include <execinfo.h>
#include <unordered_set>
#include <variant>
#include <vector>

using std::cout;
using std::cerr;
using std::endl;

inline bool debugPrepare{};

inline size_t hashCombine(const size_t& h,
		   const size_t& x)
{
  return h^(x+0x9e3779b9+(h<<6)+(h>>2));
}

template <typename T>
requires requires(T t){
  t.hash();}
struct std::hash<T>
{
  size_t operator()(const T& t) const
  {
    return t.hash();
  }
};

template <typename...T>
size_t hashTuple(const std::tuple<T...>& t)
{
  return std::apply([](auto&&...arg)
  {
    std::array a{std::hash<std::decay_t<decltype(arg)>>{}(arg)...};
    
    size_t h{};
    
    for(size_t i=0;i<a.size();i++)
      h=hashCombine(h,a[i]);
    
    return h;
  },t);
}

template <typename A,
	  typename B>
struct std::hash<std::pair<A,B>>
{
  size_t operator()(const std::pair<A,B>& c) const
  {
    return hashTuple(std::tie(c.first,c.second));
  }
};

template <typename T>
requires requires(T t){
  t.begin();
  t.end();}
struct std::hash<T>
{
  size_t operator()(const T& s) const
  {
    size_t h{};
    
    for(const auto& t: s)
      h=hashCombine(h,std::hash<std::decay_t<decltype(t)>>{}(t));
    
    return h;
  }
};

template <typename T>
struct std::hash<std::vector<T>>
{
  size_t operator()(const std::vector<T>& s) const
  {
    size_t h{};
    
    for(const T& t: s)
      h=hashCombine(h,std::hash<T>{}(t));
    
    return h;
  }
};

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
inline void print_backtrace_list()
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
inline void internal_crash(const int& line,
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

struct BitSet
{
  size_t size;
  
  std::vector<uint64_t> data;
  
  void resize(const size_t& size)
  {
    this->size=size;
    data.resize((size+63)/64);
  }
  
  void set(const size_t& i)
  {
    data[i/64]|=
      (uint64_t(1)<<(i%64));
  }
  
  void orWith(const BitSet& other)
  {
    for(size_t i=0;i<data.size();i++)
      data[i]|=other.data[i];
  }
  
  size_t popCount() const
  {
    size_t c=0;
    
    for(const uint64_t& w : data)
      c+=std::popcount(w);
    
    return c;
  }
};

inline std::vector<std::string> breakStringAt(const char* str,
					      const char& d)
{
  std::vector<std::string> res;
  
  const char* p=str;
  
  const char* b=p;
  
  while(*b!='\0')
    {
      const char* e=b;
      
      while(*e!='\0' and *e!=d)
	e++;
      
      res.emplace_back(b,e);
      
      b=e+(*e!='\0');
    }
  
  return res;
}

/////////////////////////////////////////////////////////////////

template <typename T,
	  size_t N>
struct std::hash<std::array<T,N>>
{
  const size_t operator()(const std::array<T,N>& m) const
  {
    size_t h{};
    
    for(size_t i=0;i<N;i++)
      h=hashCombine(h,std::hash<T>{}(m[i]));
    
    return h;
  }
};

/// Holds three components of the momentum
struct Momentum :
  std::array<double,3>
{
};

template <>
struct std::hash<Momentum>
{
  const size_t operator()(const Momentum& m) const
  {
    return std::hash<std::array<double,3>>{}((std::array<double,3>&)m);
  }
};

inline std::ostream& operator<<(std::ostream& os,
			 const Momentum& mom)
{
  return os<<"{"<<mom[0]<<","<<mom[1]<<","<<mom[2]<<"}";
}

/// Negate a momentum
inline Momentum operator-(const Momentum& y)
{
  auto minus=
    [](const double& x)
    {
      if(x)
	return -x;
      else
	return x;
    };
  
  return {minus(y[0]),minus(y[1]),minus(y[2])};
}

/// Product of a scalar with a momentum
inline Momentum operator*(const double& x,
	      const Momentum& y)
{
  return {y[0]*x,y[1]*x,y[2]*x};
}

/// Product of a momentum with a scalar
inline Momentum operator*(const Momentum& x,
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

inline int registerOp()
{
  static int i=0;
  
  return i++;
}

#define PROVIDE_GETTER(NAME,COMMAND) \
  std::string get ## NAME () const    \
  {				     \
    return COMMAND;		     \
  }
  
template <typename D>
struct BasePar
{
  inline static int id=
    registerOp();
  
  const D& operator~() const
  {
    return *static_cast<const D*>(this);
  }
  
  PROVIDE_GETTER(Kappa,"");
  
  PROVIDE_GETTER(Mass,"");
  
  PROVIDE_GETTER(R,"");
  
  PROVIDE_GETTER(Charge,"");
  
  std::string getMom(const size_t& i) const
  {
    return "";
  }
  
  PROVIDE_GETTER(Residue,"");
  
  std::string describe() const
  {
    const D& d=~*this;
    
    std::ostringstream os;
    os<<D::name<<"(";
    
    const auto listPars=
      [&os,
       names=d.memberNames()](const auto& listPars,
			      const size_t& i,
			      const auto& head,
			      const auto&...tail)
      {
	os<<names[i]<<"="<<head;
	if constexpr(sizeof...(tail))
	  {
	    os<<",";
	    listPars(listPars,i+1,tail...);
	  }
      };
    
    std::apply([&listPars](const auto&...args)
    {
      listPars(listPars,0,args...);
    },d.getMembers());
    
    os<<")";
    
    return os.str();
  }
  
  size_t hash() const
  {
    return hashCombine(std::hash<int>{}(id),
		       hashTuple((~(*this)).getMembers()));
  }
  
  constexpr std::partial_ordering operator<=>(const BasePar&) const=default;
};

#define PROVIDE_MEMBERS(NAME,				\
			TAG,				\
			ARGS...)			\
  static constexpr const char* name{#NAME};		\
  							\
  static constexpr const char* tag{#TAG};		\
  							\
  auto getMembers() const				\
  {							\
    return std::tie(ARGS);				\
  }							\
  							\
  const std::vector<std::string> memberNames() const	\
  {							\
    return breakStringAt(#ARGS,',');			\
  }

struct Smear :
  BasePar<Smear>
{
  double kappa{};
  
  double n{};
  
  Momentum mom{};
  
  PROVIDE_MEMBERS(Smear,SM,kappa,n,mom);
  PROVIDE_GETTER(Kappa,std::format("{}",kappa));
  PROVIDE_GETTER(R,std::format("{}",n));
  
  std::string getMom(const size_t& i) const
  {
    return std::format("{}",mom[i]);
  }
  
  Smear dag() const
  {
    return *this;
  }
  
  constexpr std::partial_ordering operator<=>(const Smear&) const=default;
};

struct Phase :
  BasePar<Phase>
{
  Momentum mom;
  
  PROVIDE_MEMBERS(Phase,PH,mom);
  
  std::string getMom(const size_t& i) const
  {
    return std::format("{}",mom[i]);
  }
  
  Phase dag() const
  {
    Phase res{*this};
    
    res.mom=-mom;
    
    return res;
  }
  
  constexpr std::partial_ordering operator<=>(const Phase&) const=default;
};

struct Gamma :
  BasePar<Gamma>
{
  size_t iGamma{};
  
  PROVIDE_MEMBERS(Gamma,G,iGamma);
  
  PROVIDE_GETTER(R,std::format("{}",iGamma));
  
  PROVIDE_GETTER(Charge,"0");
  
  Gamma dag() const
  {
    return *this;
  }
  
  constexpr std::partial_ordering operator<=>(const Gamma&) const=default;
};

struct DeltaT :
  BasePar<DeltaT>
{
  size_t t;
  
  PROVIDE_MEMBERS(DeltaT,S,t);
  
  DeltaT dag() const
  {
    return *this;
  }
  
  constexpr std::partial_ordering operator<=>(const DeltaT&) const=default;
};

struct Prop :
  BasePar<Prop>
{
  double kappa{};
  
  double mass{};
  
  int r{};
  
  double charge{};
  
  Momentum mom{};
  
  double residue{};
  
  PROVIDE_MEMBERS(Prop,-,kappa,mass,r,charge,mom,residue)
  
  PROVIDE_GETTER(Kappa,std::format("{}",kappa));
  
  PROVIDE_GETTER(Mass,std::format("{}",mass));
  
  PROVIDE_GETTER(R,std::format("{}",(r==-1)?0:1));
  
  PROVIDE_GETTER(Charge,std::format("{}",charge));
  
  std::string getMom(const size_t& i) const
  {
    return std::format("{}",mom[i]);
  }
  
  PROVIDE_GETTER(Residue,std::format("{}",residue));
  
  Prop dag() const
  {
    if(r!=-1 and r!=1)
      CRASH("Impossible value for r, %d",r);
    
    Prop res=*this;
    
    res.r=-r;
    
    res.mom=-mom;
    
    return res;
  }
  
  constexpr std::partial_ordering operator<=>(const Prop&) const=default;
};

/////////////////////////////////////////////////////////////////

using Pars=
  std::variant<Smear,Phase,Gamma,DeltaT,Prop>;

inline Pars dag(const Pars& pars)
{
  return
    std::visit([](const auto& o) -> Pars
      {
	return o.dag();
      },pars);
}

inline std::string describe(const Pars& pars)
{
  return
    std::visit([](const auto& o) -> std::string
      {
	return o.describe();
      },pars);
}

template <typename...T>
std::string describe(const std::variant<T...>& t)
{
  return std::visit([](const auto& t)
  {
    return describe(t);
  },t);
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
  
  std::string describe() const
  {
    std::string res=::describe(pars);
    
    if(rhs)
      res+="*"+rhs->describe();
    
    return res;
  }
  
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

inline Oper operator*(const Oper& a,
	       const Oper& b)
{
  Oper c(a);
  
  Oper* d=&c;
  while(d->rhs)
    d=&*d->rhs;
  
  d->rhs=std::make_unique<Oper>(b);
  
  return c;
}

struct Source :
  BasePar<Source>
{
  inline static size_t glbCount{};
  
  size_t count;
  
  int tSelect;
  
  PROVIDE_MEMBERS(Source,SO,count,tSelect);
  
  size_t hash() const
  {
    return hashTuple(std::tie(count));
  }
  
  Source(const int tSelect=-1) :
    count{glbCount++},
    tSelect(tSelect)
  {
  }
  
  constexpr std::partial_ordering operator<=>(const Source&) const=default;
};

struct Line;

struct ExtendSelectingTPars
{
  Pars pars;
  
  int tSelect{-1};
  
  size_t hash() const
  {
    return hashTuple(std::tie(pars,tSelect));
  }
  
  std::string describe() const
  {
    std::string res=::describe(pars);
    if(tSelect!=-1)
      res+="deltaT("+std::format("{}",tSelect)+")";
    
    return res;
  }
  
  std::partial_ordering operator<=>(const ExtendSelectingTPars&) const=default;
};

struct Extender
{
  ExtendSelectingTPars op;
  
  std::vector<std::pair<double,Line>> rhs;
  
  std::string describe() const;
  
  bool isScalar() const;
  
  bool isCombiner() const ;
  
  bool isExtendingSource() const;
  
  bool isSelectingTime() const
  {
    return op.tSelect!=-1;
  }
  
  Line* maybeGetTrivialRhs();
  
  Extender* maybeGetNestedExtender();
  
  constexpr bool operator==(const Extender&) const;
  constexpr std::partial_ordering operator<=>(const Extender&) const;
};

struct Line
{
  std::variant<Extender,Source> data;
  
  std::string name;
  
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
	  std::get_if<DeltaT>(&extender.op.pars);
	
	size_t n=
	  dt!=nullptr;
	
	if(n)
	  {
	    extender.op.tSelect=dt->t;
	    extender.op.pars=Gamma{.iGamma=0};
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
    data{Extender{ExtendSelectingTPars{.pars=pars,.tSelect=tSelect},rhs}}
  {
    trySimplify();
  }
  
  Line(const Oper& oper,
       const Line& line) :
    data(Extender{{.pars=oper.pars}})
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
  
  Line(const Line& oth,
       const std::string& name) :
    data(oth.data),name(name)
  {
  }
  
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
  if(auto cmp=op<=>other.op;cmp!=0)
    return cmp;
  
  if(auto cmp=rhs<=>other.rhs;cmp!=0)
    return cmp;
  
  return op.tSelect<=>other.op.tSelect;
}

inline bool Extender::isExtendingSource() const
{
  return rhs.size() and not rhs.front().second.isExtender();
}

inline Line* Extender::maybeGetTrivialRhs()
{
  if(isCombiner())
    return nullptr;
  else
    return &rhs.front().second;
}

inline std::string Extender::describe() const
{
  std::ostringstream os;
  os<<op.describe();
  
  if(rhs.size()==0)
    CRASH("cannot describe a line extender which does not extend an existing line");
  else
    {
      if(op.tSelect!=-1)
	os<<"*"<<::describe(DeltaT{.t=(size_t)op.tSelect});
      
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

template <typename T>
requires requires(T t){
  t.describe();}
std::string describe(const T& t)
{
  return t.describe();
}

inline size_t Line::maybeRemoveUselessScalar()
{
  // Check if is an extender
  
  // Nested remove
  int nDeep{};
  if(Extender* e=
     std::get_if<Extender>(&this->data))
    {
      for(auto& [w,l] : e->rhs)
	nDeep+=l.maybeRemoveUselessScalar();
      
      if(Line* nestedRhs=e->maybeGetTrivialRhs())
	{
	  if(Source* nestedS=
	    std::get_if<Source>(&nestedRhs->data);
	     e->isScalar() and (e->op.tSelect==-1 or (nestedS and nestedS->tSelect==e->op.tSelect)))
	    {
	      *this=*nestedRhs;
	      
	      return nDeep+1;
	    }
	  
	  Extender* nestedE=
	    std::get_if<Extender>(&nestedRhs->data);
	  
	  if(nestedE and nestedE->isScalar() and e->op.tSelect==-1)
	    {
	      e->op.tSelect=nestedE->op.tSelect;
	      e->rhs=std::move(nestedE->rhs);
	      
	      return nDeep+1;
	    }
	}
    }
  
  return nDeep;
}

inline Line operator*(const Oper& oper,
	       const Line& line)
{
  return {oper,line};
}

inline Line operator*(const double& c,
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

inline Line operator*(const Line& line,
	       const double& c)
{
  return c*line;
}

inline bool doNotMergeLinComb{};

inline Line operator+(const Line& a,
		      const Line& b)
{
  if(const Extender *ae=std::get_if<Extender>(&a.data),*be=std::get_if<Extender>(&b.data);ae and be and ae->op==be->op and not doNotMergeLinComb)
    {
      Extender res=*ae;
      
      res.rhs.insert(res.rhs.end(),be->rhs.begin(),be->rhs.end());
      
      return res;
    }
  else
    return Extender{.op{.pars{Gamma{.iGamma=0}}},.rhs{std::vector{std::make_pair(1.0,a),{1.0,b}}}};
}

struct Contr
{
  std::string name;
  
  std::set<std::pair<int,int>> gammaList;
  
  size_t hash() const
  {
    return hashTuple(std::tie(name,gammaList));
  }
  
  std::partial_ordering operator<=>(const Contr&) const=default;
  
  std::string describe() const
  {
    std::ostringstream os;
    
    os<<"contr("<<name<<",{";
    for(const auto& [i,j] : gammaList)
      os<<i<<"-"<<j<<" ";
    os<<"})";
    
    return os.str();
  }
};

template <typename K,
	  typename V>
struct InsertionOrderedMap
{
  std::vector<std::pair<K,V>> values;
  
  std::map<K,size_t> index;
  
  size_t size() const
  {
    return values.size();
  }
  
  size_t empty() const
  {
    return values.empty();
  }
  
  auto begin() const
  {
    return values.begin();
  }
  
  auto end() const
  {
    return values.end();
  }
  
  V& operator[](const K& k)
  {
    const auto it=index.find(k);
    
    if(it==index.end())
      {
	const size_t pos=values.size();
	values.emplace_back(std::pair{k,V{}});
	index[k]=pos;
	
	return values.back().second;
      }
    else
      return values[it->second].second;
  }
  
  std::partial_ordering operator<=>(const InsertionOrderedMap&) const=default;
};

struct Node
{
  template <typename T>
  struct Hashability
  {
    size_t hash() const
    {
      const T& t=*static_cast<const T*>(this);
      
      return hashTuple(std::tie(t.task,t.deps));
    }
    
    constexpr std::partial_ordering operator<=>(const Hashability&) const=default;
  };
  
  struct Shape :
    Hashability<Shape>
  {
    using Task=
      std::variant<Contr,ExtendSelectingTPars,Source>;
    
    Task task;
    
    InsertionOrderedMap<Node*,double> deps;
    
    std::partial_ordering operator<=>(const Shape&) const=default;
    
    bool operator==(const Shape&) const=default;
    
    std::string describe() const
    {
      std::ostringstream os;
      
      os<<">>> "<<::describe(task);
      
      if(deps.size())
	{
	  os<<" on (";
	  for(int i{};const auto& [n,w] : deps)
	    {
	      if(i++)
		os<<" + ";
	      
	      if(w!=1)
		os<<w<<" * ";
	      
	      os<<::describe(*n)<<" ";
	    }
	  os<<")";
	}
      
      return os.str();
    }
  };
  
  struct Hasher
  {
    using is_transparent=void;
    
    bool operator()(const std::unique_ptr<Node>& node) const
    {
      return node->shape.hash();
    }
    
    bool operator()(const Node::Shape& key) const
    {
      return key.hash();
    }
  };
  
  struct Comparer
  {
    using is_transparent=void;
    
    bool operator()(const std::unique_ptr<Node>& a,
		    const std::unique_ptr<Node>& b) const
    {
      return a->shape==b->shape;
    }
    
    bool operator()(const Node::Shape& a,
		    const std::unique_ptr<Node>& b) const
    {
      return a==b->shape;
    }
  };
  
  Shape shape;
  
  std::string name;
  
  std::vector<Node*> users;
  
  int id;
  
  size_t nRemainingUsers;
  
  BitSet isReachable;
  
  size_t nReachable{};
  
  void computeNReachable(const size_t& n)
  {
    if(not nReachable)
      {
	isReachable.resize(n);
	isReachable.set(id);
	
	for(Node* user : users)
	  {
	    user->computeNReachable(n);
	    isReachable.orWith(user->isReachable);
	  }
	
	nReachable=isReachable.popCount();
      }
  }
  
  size_t nScheduledDeps;
  
  size_t lastUse;
  
  bool isFreeable() const
  {
    return nRemainingUsers==0;
  }
  
  bool isReadyForSchedule() const
  {
    return nScheduledDeps==shape.deps.size();
  }

  std::string getSource() const
  {
    std::ostringstream os;
    
    if(shape.deps.size()==1 and shape.deps.begin()->second==1.0)
      os<<shape.deps.begin()->first->name;
    else
      {
	os<<"LINCOMB "<<shape.deps.size()<<" ";
	for(const auto& [d,w] : shape.deps)
	  os<<d->name<<" "<<w<<" ";
      }
    
    return os.str();
  }
  
  std::string describe() const
  {
    std::ostringstream os;
    
    return shape.describe();
    
    return os.str();
  }
  
  int memoryCostIfRun() const
  {
    return std::visit(Overload{[](const Source&){return 1;},[](const ExtendSelectingTPars&){return 1;},[](const Contr&){return 0;}},shape.task);
  }
  
  int nFreedIfScheduled() const
  {
    int nFreed=0;
    for(const auto& [d,w] : shape.deps)
      if(d->nRemainingUsers==1)
	nFreed++;
    
    return nFreed;
  }
  
  int nReachableReleasedIfScheduled() const
  {
    BitSet b;
    b.resize(isReachable.size);
    
    for(const auto& [d,w] : shape.deps)
      if(d->nRemainingUsers==1)
	b.orWith(d->isReachable);
    
    return b.popCount();
  }
  
  std::array<int,3> readyness() const
  {
    return {nReachableReleasedIfScheduled(),-memoryCostIfRun(),(int)nReachable};
  }
};

struct Contracter
{
  Contr pars;
  
  std::set<std::pair<Line,Line>> traces;
  
  Contracter(const std::string& name) :
    pars{.name=name}
  {
  }
  
  void addGammas(const int& sink,
		 const int& sour)
  {
    pars.gammaList.emplace(sink,sour);
  }
  
  void tr(const Line& bw,
	  const Line& fw)
  {
    if(not traces.emplace(bw,fw).second)
      CRASH("inserting tr(%s,%s), trace already existing",bw.name.c_str(),fw.name.c_str());
  }
};

struct PrintLiner
{
  std::vector<std::string> entries;
  
  size_t nColumns{};
  
  size_t nLines{};
  
  std::vector<std::vector<std::string>> commentsPerRow{{}};
  
  size_t id(const size_t& iRow,
	    const size_t& iCol) const
  {
    return iCol+nColumns*iRow;
  }
  
  void newLine()
  {
    nLines++;
    commentsPerRow.emplace_back();
    
    if(nColumns==0)
      nColumns=entries.size();
    else
      if(entries.size()!=nColumns*nLines)
	CRASH("total entries number %zu is not equal to the product of the number of lines %zu and columns %zu",entries.size(),nLines,nColumns);
  }
  
  template <typename T>
  PrintLiner& operator<<(const T& t)
  {
    entries.push_back((std::ostringstream()<<t).str());
    
    return *this;
  }
  
  void normalize()
  {
    std::vector<size_t> maxLen(nColumns);
    for(size_t iCol=0;iCol<nColumns;iCol++)
      {
	size_t m{};
	for(size_t iRow=0;iRow<nLines;iRow++)
	  m=std::max(m,entries[id(iRow,iCol)].length());
	
	for(size_t iRow=0;iRow<nLines;iRow++)
	  for(size_t i=entries[id(iRow,iCol)].length();i<=m;i++)
	    entries[id(iRow,iCol)]+=" ";
      }
  }
  
  void comment(const std::string& comment)
  {
    commentsPerRow.back().push_back(comment);
  }
  
  std::string print() const
  {
    std::ostringstream os;
    
    for(size_t iRow=0;iRow<=nLines;iRow++)
      {
	if(const std::vector<std::string>& comments=commentsPerRow[iRow];comments.size())
	  for(const std::string& comment : comments)
	    os<<"/* "<<comment<<"*/"<<endl;
	
	if(iRow<nLines)
	  {
	    for(size_t iCol=0;iCol<nColumns;iCol++)
	      os<<entries[id(iRow,iCol)]<<" ";
	    os<<endl;
	  }
      }
    
    return os.str();
  }
};

struct Run
{
  bool debugFree{};
  
  bool debugContr{};
  
  std::vector<std::unique_ptr<Contracter>> contracters;
  
  std::unordered_set<std::unique_ptr<Node>,
		     Node::Hasher,
		     Node::Comparer> nodes;
  
  std::vector<Node*> executeList;
  
  void clearOutput()
  {
    std::apply([](auto&...c)
    {
      (c.clear(),...);
    },std::tie(nodes,executeList));
  }
  
  Node* maybeActuallyInsertNode(Node::Shape&& shape,
				const std::string& name="")
  {
    auto it=nodes.find(shape);
    
    if(it==nodes.end())
      it=nodes.insert(std::make_unique<Node>(std::move(shape),name)).first;
    
    Node& n=*it->get();
    if(n.name.empty() and not name.empty()) // might be intermediate passage
      n.name=name;
    else
      if((not name.empty()) and n.name!=name)
	CRASH("node %s has already been given name \"%s\", cannot be renamed \"%s\"",n.describe().c_str(),n.name.c_str(),name.c_str());
    
    return it->get();
  }
  
  Node* insertLine(const Line& line)
  {
    Node::Shape shape;
    
    std::visit(Overload{[this,
			 &shape](const Extender& e)
    {
      shape.task=ExtendSelectingTPars{.pars=e.op.pars,.tSelect=e.op.tSelect};
      
      for(const auto& [w,l] : e.rhs)
	shape.deps[insertLine(l)]+=w;
    },
	  [&shape](const Source& s)
	  {
	    shape.task=s;
	  }},line.data);
    
    return maybeActuallyInsertNode(std::move(shape),line.name);
  };
  
  Contracter& getTracer(const std::string& name)
  {
    contracters.push_back(std::make_unique<Contracter>(name));
    
    return *contracters.back();
  }
  
  void compile()
  {
    clearOutput();
    
    /// Avoid inserting too many times
    std::map<Line,Node*> linesToNode;
    
    for(const std::unique_ptr<Contracter>& contracter : contracters)
      for(const std::pair<Line,Line> &contr : contracter->traces)
	{
	  std::array<const Line*,2> linesOfContr{&contr.first,&contr.second};
	  
	  Node::Shape contrShape{.task=contracter->pars};
	  for(int i=0;i<2;i++)
	    {
	      const Line& l=*linesOfContr[i];
	      
	      auto it=linesToNode.find(l);
	      if(it==linesToNode.end())
		{
		  for(const auto& e : linesToNode)
		    if(e.first.name==l.name)
		      CRASH("Trying to insert with name \"%s\" a line %s but it is already associated to %s",e.first.name.c_str(),l.describe().c_str(),e.first.describe().c_str());
		  
		  it=linesToNode.insert(std::pair{l,insertLine(l)}).first;
		}
	      
	      contrShape.deps[it->second]+=1.0;
	    }
	  
	  maybeActuallyInsertNode(std::move(contrShape));
	}
    
    // Label nodes
    for(int i=0;const std::unique_ptr<Node>& n : nodes)
      n->id=i++;
    
    // Add users
    for(auto& n : nodes)
      for(auto& [d,w] : n->shape.deps)
	  d->users.push_back(&*n);
    
    // Set isReachable
    for(const std::unique_ptr<Node>& n : nodes)
      n->computeNReachable(nodes.size());
    
    // Reset the remaining user and scheduled deps
    for(auto& n : nodes)
      {
	n->nRemainingUsers=n->users.size();
	n->nScheduledDeps=0;
	n->lastUse=std::numeric_limits<size_t>::max();
      }
    
    // report
    if(debugPrepare)
      {
	cout<<"==========="<<endl;
	for(const auto& n : nodes)
	  {
	    cout<<describe(*n)<<endl;
	    cout<<" used by "<<n->users.size()<<endl;
	  }
      }
    
    std::vector<Node*> readyNodes;
    
    for(const std::unique_ptr<Node>& n : nodes)
      if(n->shape.deps.empty())
        readyNodes.push_back(n.get());
    
    while(not readyNodes.empty())
      {
	std::sort(readyNodes.begin(),
		  readyNodes.end(),
		  [](const Node* a,
		     const Node* b)
		  {
		    return a->readyness()<b->readyness();
		  });
	
	if(debugPrepare)
	  {
	    cout<<"Ready nodes:"<<endl;
	    for(const Node* n : readyNodes)
	      cout<<describe(*n)<<" nFreed if run: "<<n->nFreedIfScheduled()<<", n remaining users: "<<n->nRemainingUsers<<" nReachable: "<<n->nReachable<<" nReachableReleasedIfScheduled: "<<n->nReachableReleasedIfScheduled()<<endl;
	  }
	
	Node* toBeRun=readyNodes.back();
	readyNodes.pop_back();
	
	if(debugPrepare)
	  {
	    if(readyNodes.empty())
	      cout<<"Chosen the only possible node"<<endl;
	    else
	      cout<<"Chosen: "<<describe(*toBeRun)<<endl;
	  }
	
	executeList.push_back(toBeRun);
	
	for(const auto& [d,w] : toBeRun->shape.deps)
	  {
	    d->nRemainingUsers--;
	    if(d->isFreeable())
	      {
		if(debugPrepare)
		  cout<<"Now "<<describe(*d)<<" is freeable at step "<<executeList.size()<<endl;
		d->lastUse=executeList.size()-1;
	      }
	  }
	
	for(Node* user : toBeRun->users)
	  {
	    user->nScheduledDeps++;
	    if(user->isReadyForSchedule())
	      readyNodes.push_back(user);
	  }
	if(debugPrepare)
	  cout<<"-------"<<endl;
      }
    
    // Check
    for(const std::unique_ptr<Node>& n : nodes)
      if(not n->isReadyForSchedule())
	CRASH("Node %s not scheduled",n->describe().c_str());
    
    if(debugPrepare)
      cout<<"/////////////////////////////////////////////////////////////////"<<endl;
    
    if(debugPrepare)
      cout<<"Final schedule: "<<endl;
    int memoryPressure=0;
    
    auto freeWhatPossible=
      [this](const size_t& i)
      {
	std::vector<Node*> freeable;
	for(const std::unique_ptr<Node>& maybeFreeable : nodes)
	  if(i and maybeFreeable->lastUse==i)
	    freeable.push_back(maybeFreeable.get());
	
	size_t nReallyFreed=0;
	if(freeable.size())
	  for(const Node* f : freeable)
	    {
	      if(debugPrepare)
		cout<<" Free:("<<describe(*f)<<")"<<endl;
	      if(std::holds_alternative<ExtendSelectingTPars>(f->shape.task))
		nReallyFreed++;
	    }
	
	return nReallyFreed;
      };
    
    int maxMemoryPressure{};
    for(size_t i=0;i<executeList.size();i++)
      {
	const Node& n=*executeList[i];
	if(std::holds_alternative<ExtendSelectingTPars>(n.shape.task))
	  {
	    memoryPressure++;
	    maxMemoryPressure=std::max(maxMemoryPressure,memoryPressure);
	  }
	if(debugPrepare)
	  cout<<i<<" -- mp: "<<memoryPressure<<" -- "<<describe(n)<<endl;
	
	memoryPressure-=freeWhatPossible(i);
      }
    
    memoryPressure-=freeWhatPossible(executeList.size());
    cout<<"Final memory pressure: "<<memoryPressure<<endl;
    cout<<"Maximal memory pressure: "<<maxMemoryPressure<<endl;
    
    cout<<"/////////////////////////////////////////////////////////////////"<<endl;
    
    size_t iSource{};
    size_t iOp{};
    for(Node* n : executeList)
      if(n->name.empty())
	n->name=std::visit(Overload{[&iSource](const Source& source)
	{
	  return "SOURCE"+std::format("{}",iSource++);
	},
	      [&iOp](const ExtendSelectingTPars& pars)
	      {
		return "P"+std::format("{}",iOp++);
	      },
	      [](const Contr& contr)
	      {
		return contr.name;
	      }},n->shape.task);
  }
  
  std::string printSources() const
  {
    const size_t nSources=
      std::count_if(executeList.begin(),
		    executeList.end(),
		    [](const Node* n)
		    {
		      return std::holds_alternative<Source>(n->shape.task);
		    });
    
    cout<<"NSources "<<nSources<<endl;
    
    PrintLiner table;
    table<<"Name"<<"NoiseType"<<"Tins"<<"Store";
    table.newLine();
    
    for(const Node* n : executeList)
      if(const Source* s=std::get_if<Source>(&n->shape.task))
	{
	  table<<n->name<<"Z4"<<s->tSelect<<"0";
	  table.newLine();
	}
    
    table.normalize();
    
    return table.print();
  }
  
  std::string printLines() const
  {
    std::ostringstream os;
    
    PrintLiner table;
    table<<"Name"<<"Ins"<<"SourceName"<<"Tins"<<"Kappa"<<"Mass"<<"R"<<"Charge"<<"ThetaX"<<"ThetaY"<<"ThetaZ"<<"Residue"<<"Store";
    table.newLine();
    
    const auto includeCommentsOnFree=
      [&table,
       this](const size_t& i)
      {
	for(const std::unique_ptr<Node>& maybeFreeable : nodes)
	  if(i and maybeFreeable->lastUse==i)
	    table.comment("free "+maybeFreeable->name);
      };
    
    size_t nProps{};
    
    for(size_t i=0;i<executeList.size();i++)
      {
	const Node* n=executeList[i];
	
	if(const ExtendSelectingTPars* e=std::get_if<ExtendSelectingTPars>(&n->shape.task))
	  std::visit([&n,
		      &nProps,
		      &table,
		      &tSelect=e->tSelect](const auto& par)
	  {
	    table<<n->name;
	    
	    table<<std::remove_reference_t<decltype(par)>::tag;
	    
	    table<<n->getSource();
	    
	    table<<tSelect;
	    
	    table<<par.getKappa()<<par.getMass()<<par.getR()<<par.getCharge();
	    for(size_t i=0;i<3;i++)
	      table<<par.getMom(i);
	    
	    table<<par.getResidue();
	    
	    table<<"0";
	    
	    table.newLine();
	    nProps++;
	  },e->pars);
	
	if(debugContr)
	  if(const Contr* c=std::get_if<Contr>(&n->shape.task))
	    {
	      std::vector<std::string*> memb;
	      for(const auto& [d,w] : n->shape.deps)
		memb.push_back(&d->name);
	      
	      size_t i[2]={0,1};
	      if(memb.size()!=2)
		i[1]=0;
	      ///CRASH("contraction with %zu members, different from 2",memb.size());
	      
	      table.comment(std::format("{} <- Contr({},{})",c->name,*memb[i[0]],*memb[i[1]]));
	    }
	
	if(debugFree)
	  includeCommentsOnFree(i);
      }
    
    if(debugFree)
      includeCommentsOnFree(executeList.size());
    
    table.normalize();
    
    os<<"NProps "<<nProps<<"\n"<<endl;
    os<<table.print();
    
    return os.str();
  }
  
  std::string printContr() const
  {
    std::ostringstream os;
    os<<"NMes2PtsContrBlocks  "<<contracters.size()<<endl;
    os<<endl;
    
    for(const std::unique_ptr<Contracter>& c : contracters)
      {
	os<<"NMes2PtsContr "<<c->traces.size()<<endl;
	
	PrintLiner table;
	for(const auto& [bw,fw] : c->traces)
	  {
	    table<<c->pars.name<<" "<<bw.name<<" "<<fw.name;
	    table.newLine();
	  }
	table.normalize();
	
	os<<table.print()<<endl;
	
	const auto asGamma=
	  [](const size_t& i)->std::string
	  {
	    return (const char* []){"S0","V1","V2","V3","V0","P5","A1","A2","A3","A0","T1","T2","T3","B1","B2","B3"}[i];
	  };
	
	os<<"NGammaContr "<<c->pars.gammaList.size()<<endl;
	for(const auto& [a,b] : c->pars.gammaList)
	  os<<asGamma(a)<<asGamma(b)<<endl;
	os<<endl;
      }
    
    return os.str();
  }
};

inline bool Extender::isCombiner() const
{
  return rhs.size()>1 or (rhs.size()==1 and rhs.front().first!=1.0);
}

inline bool Extender::isScalar() const
{
  const Gamma* g=
    std::get_if<Gamma>(&op.pars);
  
  return
    g!=nullptr and g->iGamma==0;
}
