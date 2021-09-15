

class hardSphere
{
public:
  hardSphere(const json_t & j) : hardSphere(j["R0"].get<real_t>() ) {}
  hardSphere(real_t R0_) : _R0(R0_){}
  
  static std::string name() {return "hardSphere";}
  
  real_t R0() const  {return _R0;}
  
private:
  real_t _R0;
};


