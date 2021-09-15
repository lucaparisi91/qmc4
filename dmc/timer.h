#include <chrono>
#include <string>
#include <unordered_map>
class timer
{
public:
  timer();
  using time_t =   decltype (std::chrono::high_resolution_clock::now());
    
  void start();
  void stop();
  auto milliseconds() { return std::chrono::duration_cast<std::chrono::milliseconds>(duration).count(); }

  auto seconds() { return std::chrono::duration_cast<std::chrono::seconds>(duration).count(); }

  
protected:
  void addElapsedTime();
private:
  time_t _start;
  time_t _stop;
  size_t count;
  std::chrono::duration<double> duration;
  
};

class timers
{
private:
  static timers *singleton ;
  timers(){};

public:
  static timers* get();
  void add(std::string name);
  timer & get(const std::string & name);

  auto begin() {return _timers.begin();}
  auto end() {return _timers.end();}
  
private:
  
  std::unordered_map<std::string,timer> _timers;
};

timers & getTimers();
