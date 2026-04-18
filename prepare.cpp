#include <array>
#include <cstdarg>
#include <cstddef>
#include <iostream>
#include <execinfo.h>
#include <variant>
#include <vector>

using std::cout;
using std::cerr;
using std::endl;

/// write the list of called routines
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

#define CRASH(...) internal_crash(__LINE__,__FILE__,__PRETTY_FUNCTION__,__VA_ARGS__)

/// crashes emitting the message
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

struct Mom :
  std::array<double,3>
{
};

Mom operator-(const Mom& y)
{
  return {-y[0],-y[1],-y[2]};
}

Mom operator*(const double& x,
	      const Mom& y)
{
  return {y[0]*x,y[1]*x,y[2]*x};
}

Mom operator*(const Mom& x,
	      const double& y)
{
  return y*x;
}


struct Smear;

struct Phase;

struct Gamma;

/////////////////////////////////////////////////////////////////

struct Smear
{
  double kappa{};
  
  double n{};
  
  Mom mom{};
  
  Smear dag() const
  {
    return *this;
  }
};

struct Phase
{
  Mom mom;
  
  Phase dag() const
  {
    Phase res{*this};
    
    res.mom=-mom;
    
    return res;
  }
};

struct Gamma
{
  size_t iGamma{};
  
  Gamma dag() const
  {
    return *this;
  }
};

struct DeltaT
{
  size_t t;
  
  DeltaT dag() const
  {
    return *this;
  }
};

using Pars=
  std::variant<Smear,Phase,Gamma,DeltaT>;

Pars dag(const Pars& pars)
{
  return
    std::visit([](const auto& o) -> Pars
      {
	return o.dag();
      },pars);
};

struct Oper;

Oper operator*(const Oper& a,
	       const Oper& b);

struct Oper
{
  Pars pars;
  
  std::vector<Oper> source;
  
  /// Dagger of an operation is meaning
  Oper dag() const
  {
    if(source.size()>1)
      CRASH("source.size()=%zu greater than 1",source.size());
    
    const Pars dagPars=
      ::dag(pars);
    
    if(source.size()==0)
      return dagPars;
    else
      return source[0].dag()*dagPars;
  }
  
  Oper(const Pars& pars,
       const std::vector<Oper>& source) :
    pars{pars},source(source)
  {
  }
  
  template <typename T>
  Oper(const T& pars) :
    pars{pars}
  {
  }
};

Oper operator*(const Oper& a,
	       const Oper& b)
{
  if(a.source.size())
    CRASH("lhs has already a source");
  
  return Oper(a.pars,{b});
}

/////////////////////////////////////////////////////////////////

int main()
{
  const Mom mom{0,1,0};
  
  const Smear sme{.kappa=0.4,.n=40,.mom=mom};
  const Phase phase{.mom=2*mom};
  
  const Oper op=sme*phase*sme.dag();
  
  return 0;
}
