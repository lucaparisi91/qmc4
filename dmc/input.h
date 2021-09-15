#ifndef INPUT_H
#define INPUT_H


#include <nlohmann/json.hpp>
#include <iostream>
#include <fstream>
#include <utility>


std::string to_string(const std::string & s);

class missingKey : public std::exception
{
public:
  missingKey(std::string message_) : message("Missing key: " + message_) {}

  missingKey(std::string message_,std::string parent) :  message("Missing key [" + message_  + "] in json object: \n" + parent ) {}
  
  
  missingKey(int index,std::string parent) : message("Missing Index [" +
						     
						     std::to_string(index)  + "] in json object: \n" + parent ) {}

  missingKey(int index) : message("Missing Index [" + std::to_string(index)  + "] " ) {}

  
  missingKey(int index,int maxLen) : message("Missing Index [" + std::to_string(index)  + "] . Max Range is : " + std::to_string(maxLen) ) {}


  
   virtual const char* what() const throw()
    {
      return message.c_str() ;
      
    }

  private:
  std::string message;
};


class invalidType : public std::exception
{
public:
	invalidType(std::string message_) : message(message_) {}
	
	 virtual const char* what() const throw()
  	{
    	return message.c_str() ;
  	}
  private:
  	std::string message;
};



class input;


input loadFromFile(std::string filename);

class input;

class inputIterator
{
public:
  using nlohmannIterator_t = decltype( std::declval<nlohmann::json>().begin() );

  inputIterator(nlohmannIterator_t it_) : it(it_) {}

  
  input  operator *();

  auto & operator++()
  {
    ++it;
    return *this;
  }

  
  inputIterator operator++(int)
  {
    it++;
    return inputIterator(it);
  }

  bool operator!=(inputIterator it2) const
  {
    return it != it2.it;
  }
  
private:

  nlohmannIterator_t it;
};



class input
{
public:

  template<class T>
  static T load(const input & j2, T value )
        {
          if (!(j2 == nullptr))
          {
            return j2.get<T>();
          }
          else
          {
            return value;
          }
          
        }

  input(nlohmann::json * j_ , input * parent , const std::string & key  ) :j(j_),ownJson(false),accessKey(key),_parent(parent) {_parent=nullptr;}
  
  
  input(nlohmann::json * j_ ) :j(j_),ownJson(false),accessKey("") {_parent=nullptr;}
  
  input(const nlohmann::json & j_) : j(new nlohmann::json(j_)),ownJson(true) {_parent=nullptr;}
  
  input(const input & jI) : input( *(jI.j) ) {_parent=nullptr;accessKey=jI.accessKey;}
  input(input * jI) :  input( jI->j ) {_parent=jI->_parent;accessKey=jI->accessKey; }
  
  size_t size() const { return j->size();}

  bool is_structured() const {return std::as_const(j)->is_structured() ;}

  template<class T>
  inputIterator find(const T & s) const {return inputIterator( std::as_const(j)->find(s));}

  
  template<class func_t>
  inputIterator find_if(func_t Func) const {

    for (auto it = begin() ; it != end() ;++it)
      {
	if ( Func(  *it ) )
		   {
		     return it;
		   }
      }

    return end();
  }


  template<class T>
  bool operator==(const T & o) const
  {
    return (*j)==o;
  }
  
  
  input() : j(new nlohmann::json()),ownJson(true),accessKey(""){_parent=nullptr;}
  
  
  ~input()
  {
    if (ownJson)
      {
	delete j;
      }
  }

  
  template<class T>
  input  operator[](const T & s) {
    
    input res( & (*j)[s] );
    
    res._parent= this;

    {
      using namespace std;
      res.accessKey=accessKey + "->" + to_string(s);
    }
    
    return res;

  }

  template<class T>
  const input  operator[](const T & s) const {
    std::string key;
    {
      using namespace std;
      key=to_string(s);
    }
    const input res( & (*j)[s], const_cast<input*>(this) , accessKey + "->" + key );
    
    return res;

  }
  
  
  template<class T>
  void operator=(const T & o)
  {
    (*j)=o;
  }

  void operator=(const input & jI)
  {
    (*j)=*(jI.j);
  }

  void push_back(const input & jI)
  {
    j->push_back(*(jI.j));
  }
  
  void operator=(const std::vector<input> & vec)
  {
    (*j)={};

    for (const input & el : vec)
      {
	j->push_back(*( el.jptr()) );
      }
  };

  
  

  
template<class T>
auto get() const  {
  T res;

  
  try {
      res=j->get<T>();
    }
  
  catch (nlohmann::detail::type_error e)
      {
	
	std::string message;
	//Inalid Type. Error in json object: \n" + j->dump()  + "\n " +
	message= std::string("") +  e.what() + "\n";
	if (accessKey != "")
	  {
	    message+="Accessed from key : " + accessKey;
	  }
	throw invalidType( message  );
      }
  return res;

  }

  
  
  
  nlohmann::json* jptr() {return j;}
  const   nlohmann::json* jptr() const {return j;}

  inputIterator begin() {return inputIterator(j->begin()); }

  inputIterator end() {return inputIterator(j->end()); }


  inputIterator begin() const  {return inputIterator( std::as_const(j)->begin()); }
  inputIterator end() const  {return inputIterator(std::as_const(j)->end()); }

  
  
  template<class T>
  bool contains(const T & o) const {return j->contains(o);}
  
  std::string dump() const
  {
    return std::as_const(j)->dump();
  }
  
  input *  parentptr() {return _parent;}

  input  parent() const {
    if ( _parent != nullptr)
      {
	return input(_parent);
      }
    {
      input nullInput;
      return nullInput;
    }
  }
private:
  nlohmann::json* j;
  input* _parent;
  
  bool ownJson;
  std::string accessKey;
  
};


std::ostream &  operator<<(std::ostream & out , const input & inputO);

std::istream &  operator>>(std::istream & is , input & inputO);

#endif
