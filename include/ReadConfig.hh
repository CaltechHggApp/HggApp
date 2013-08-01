#ifndef ReadConfig_h
#define ReadConfig_h

#include<fstream>
#include<string>
#include<vector>
#include<map>
#include<cstring>
#include <sstream>
#include <exception>
#include <stdexcept>

#include <iostream>
using namespace std;

class ReadConfig{
 public:
  enum style{kStandard, kSection};

  ReadConfig(style type=kStandard);
  ReadConfig(string,style type=kStandard);
  bool is_init(){return isInit;}
  int read(string);
  string getParameter(string s,string section="");
  vector<string> getTokens(string s,const char *tok);
  static vector<string> tokenizeString(string s, const char *tok);
  void printAll();
  vector<string> getAllParameters();
  std::string getPath(){return path;}
 private:
  std::string path;
  int parseLine(string);  // parses lines for config files in the format key=value
  int parseSectionLine(string,string &); //parses lines for section mode config files
  string stripComments(string);  //removes comments from the string
  string stripSpaces(string); //strips leading and trailing spaces
  bool isInit;
  style readerType;
  map<string,string> parameters_;
};

ReadConfig::ReadConfig(style type):
  isInit(false),
  readerType(type)
{

}

ReadConfig::ReadConfig(string s,style type):
  isInit(false),
  readerType(type)
{
  read(s);
}

int ReadConfig::read(string s){
  path = s;
  ifstream file(s.c_str());
  string line;
  if(!file.is_open()) throw runtime_error(string("cannot open file:  ")+s);
  
  int lineNo = 0;
  int ret = 0;

  string curSection=""; //monitor the current section, if in section mode
  while(file.good()){
    getline (file,line);
    
    switch(readerType){
    case kStandard:
      try{
	ret = parseLine(line);
      }catch(const runtime_error &e){
	throw runtime_error(string("ERROR in config file:  ") + s + string("\nLine Number: ") + 
			    static_cast<ostringstream*>( &(ostringstream() << lineNo) )->str() + string("\n") +string(e.what()));
      }
      break;
      
    case kSection:
      try{
	parseSectionLine(line,curSection);
      }catch(const runtime_error &e){
	throw runtime_error(string("ERROR in config file:  ") + s + string("\nLine Number: ") + 
			    static_cast<ostringstream*>( &(ostringstream() << lineNo) )->str() + string("\n") +string(e.what()));
      }
      break;
    default:
      throw runtime_error(string("INVALID CONFIG TYPE"));
      
    }
    lineNo++;
  }
  isInit = !ret; // isInit only true if return value is 0
  if(!isInit) cout << "ERROR reading cfg" << endl;
  return (isInit==true?0:10*lineNo+ret); // returns 0 or 10*linenumber+errorcode
}

int ReadConfig::parseLine(string s){
  string sc = stripComments(s);
  if(sc.find_first_not_of(' ') == string::npos) return 0; // if there was nothing other than spaces on the line

  unsigned int eqPos = sc.find_first_of('='); // location of the assignment
  //cout << eqPos << endl;
  if(eqPos == sc.find_first_not_of(' ')) throw runtime_error(string("non l-value assignment:  ")+s);; // can't start the line with =
  if(eqPos == string::npos) throw runtime_error(string("non l-value assignment:  ")+s);; // and it must exist
  
  string var = stripSpaces(sc.substr(0,eqPos)); //variable
  string val = stripSpaces(sc.substr(eqPos+1));  //value

  if(var == "" || val == "") throw runtime_error(string("error parsing line:  ")+s);

  parameters_[var] = val;
  return 0;
}

int ReadConfig::parseSectionLine(string s,string &sec){
  string sc = stripComments(s);

  if(sc.find_first_not_of(' ') == string::npos) return 0; // if there was nothing other than spaces on the line
  if(sc.find_first_not_of(' ') == string::npos) return 0; // if there was nothing other than spaces on the line
  int eqPos = sc.find_first_of('='); // location of the assignment

  if(eqPos == sc.find_first_not_of(' ')) throw runtime_error(string("non l-value assignment:  ")+s);; // can't start the line with =
  if(eqPos == string::npos){ // might be a new section, lets check.  Section format needs to be [NAME]
    int openParPos = sc.find_first_of('[');
    int closeParPos = sc.find_first_of(']');
    if(openParPos==string::npos || closeParPos==string::npos) throw runtime_error(string("non l-value assignment:  ")+s);; // need an assignment operator
    if(openParPos >= closeParPos || closeParPos==openParPos+1) throw runtime_error(string("invalid section heading:  ")+s);

    //we have a section header of the form [NAME]
    try{
      sec = stripSpaces(sc.substr(openParPos+1,closeParPos-openParPos-1));  // set the new section heading
    }catch(runtime_error &e){
      std::cout << sc << endl << openParPos << endl << closeParPos <<endl;
      throw e;
    }
    return 0;
  }
  //not a section, just an l-value assignment
  string var = stripSpaces(sc.substr(0,eqPos)); //variable
  string val = stripSpaces(sc.substr(eqPos+1));  //value

  if(var == "" || val == "") throw runtime_error(string("error parsing line:  ")+s);

  //in section mode, prepend the section name to the variable name: ___SectionName___VariableName
  var = string("___")+sec+string("___")+var;
  parameters_[var] = val;
  return 0;

}
 
string ReadConfig::stripComments(string s){  //walk down the string and remove everything after the first # _unless its in a string 
  //return s.substr(0,s.find_first_of('#'));
  bool inQuote = false;
  string out;
  for(unsigned int i=0;i<s.size();i++){
    if(s[i] == '"') inQuote = !inQuote;
    if(s[i] == '#' && !inQuote) break;
    out+=s[i];
  }
  return out;
}

string ReadConfig::stripSpaces(string s){
  return s.substr(s.find_first_not_of(' '),1+s.find_last_not_of(' ')-s.find_first_not_of(' '));
}

string ReadConfig::getParameter(string s,string section){
  if(readerType ==kSection){ // we need to prepend the section name
    s = string("___")+section+string("___")+s;
  }
  if(parameters_.find(s) == parameters_.end()) return "";
  else return parameters_.find(s)->second;
  return "";
}

vector<string> ReadConfig::getTokens(string s,const char *tok){
  return tokenizeString(getParameter(s),tok);
}


vector<string> ReadConfig::tokenizeString(string s, const char* tok){
  vector<string> out;
  char line[4000];
  strcpy(line,s.c_str());
  char *chunk = strtok(line,tok);
  while(chunk){
    out.push_back(string(chunk));
    chunk = strtok(NULL,","); 
  }
  return out;

}

void ReadConfig::printAll(){
  cout << "Config Parameters are:" <<endl;
  map<string,string>::const_iterator it;
  for(it = parameters_.begin();it !=parameters_.end(); it++){
    cout << ">> k: " << it->first << "   v: " << it->second << endl;
  }
}

vector<string> ReadConfig::getAllParameters(){
  vector<string> out;
  map<string,string>::const_iterator it;
  for(it = parameters_.begin();it !=parameters_.end(); it++) out.push_back(it->first);
  return out;
}
#endif
