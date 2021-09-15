#ifndef POTENTIAL_FUNCTORS_H
#define POTENTIAL_FUNCTORS_H


class squareWell
{
public:
  squareWell(const json_t & j) : squareWell(j["R0"].get<real_t>(),j["V0"].get<real_t>()) {}
  squareWell(real_t R0_, real_t V0_) : R0(R0_),V0(V0_) {}
  real_t operator()(real_t x) {return x<=R0 ? -V0 : 0;}
  
  static std::string name() {return "squareWell";}
  
  
private:
  
  real_t R0;
  real_t V0;  
};


class potentialBarrier
{
public:
  potentialBarrier(const json_t & j) : potentialBarrier(j["R0"].get<real_t>(),j["V0"].get<real_t>() ) {}
  potentialBarrier(real_t R0_, real_t V0_) : R0(R0_),V0(V0_) {}
  real_t operator()(real_t x) {return x<=R0 ? V0 : 0;}
  
  static std::string name() {return "potentialBarrier";}
  
  
private:
  
  real_t R0;
  real_t V0;  
};


#endif
