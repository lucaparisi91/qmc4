#ifndef EXCEPTIONS_H
#define EXCEPTIONS_H

#include <exception>
#include <string>

class missingImplementation : public std::exception
{
public:
	missingImplementation(std::string message_) {message=message_;}
	
	 virtual const char* what() const throw()
  	{
      
    	return message.c_str() ;
  	}
  private:
  	std::string message;
};


class invalidInput : public std::exception
{
public:
	invalidInput(std::string message_) : message(message_) {}
	
	 virtual const char* what() const throw()
  	{
    	return message.c_str() ;
  	}
  private:
  	std::string message;
};

class invalidState : public std::exception
{
public:
  invalidState(std::string message_) : message(message_) {}
  
   virtual const char* what() const throw()
    {
      return message.c_str() ;
    }
  private:
    std::string message;
};

class factoryIdNotRecorded : public std::exception
{
public:
  factoryIdNotRecorded(std::string message_) : message(message_) {}
  
   virtual const char* what() const throw()
    {
      return message.c_str() ;
    }
  private:
    std::string message;
};






#endif
