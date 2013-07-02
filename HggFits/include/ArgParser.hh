#ifndef ArgParser_h
#define ArgParser_h

#include <vector>
#include <string>
#include <map>

class ArgParser{
public:
  ArgParser(int ac,char** av);
  void setInputs(int ac, char** av){
    argc = ac; argv=av; 
  }
  enum {noArg,reqArg};
  void addLongOption(std::string opt,bool argument,std::string desc="");
  void addShortOption(char opt,bool argument,std::string desc="");
  enum {optional,required};
  void addArgument(std::string name,bool required,std::string desc="");
  void addUnlimitedArgument(std::string,std::string desc=""); // unlimited arguments must be optional

  int process(std::string& ret);

  bool shortFlagPres(char);
  bool longFlagPres(std::string);

  std::string getLongFlag(std::string);
  std::string getShortFlag(char);
  std::string getArgument(std::string);
  std::vector<std::string>::const_iterator getUnlimitedArgument(std::string);
  std::vector<std::string>::const_iterator args_end(){return inputArgs.end();}

  int getStatus(){return status;}
  void reset();

  void printOptions(std::string);
private:
  int status;
  int argc;
  char** argv;

  //methods
  template <typename T,typename S>
  S getWithCheck(T,std::map<T,S>,S);
  int internalProcess(std::string &ret);

  //variables
  //maps between the flag name and the value
  typedef std::map<std::string,std::string> stringmap;
  stringmap longFlagMap;
  std::map<char,std::string> shortFlagMap;
  //maps between 
  std::map<std::string,bool> longFlagReqArgMap;
  std::map<char,bool> shortFlagReqArgMap;
  std::map<std::string,bool> longFlagPresMap;
  std::map<char,bool> shortFlagPresMap;
  std::vector<std::string> reqArgs;
  std::vector<std::string> optArgs;
  bool unlimitedArgMode;
  std::string unlimitedArg;
  std::vector<std::string> inputArgs;

  //descriptions
  stringmap longFlagDesc;
  std::map<char,std::string> shortFlagDesc;
  stringmap reqArgDesc;
  stringmap optArgDesc;
  std::string unlimitedArgDesc;

};

#endif
