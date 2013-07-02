#include "ArgParser.hh"
#include <algorithm>
#include <iostream>
#include <cstdio>

ArgParser::ArgParser(int ac,char** av){
  status=0;
  unlimitedArgMode=false;
  this->setInputs(ac,av);
}

void ArgParser::addLongOption(std::string opt,bool argument,std::string desc){
  longFlagMap[opt]="";
  longFlagReqArgMap[opt]=argument;
  longFlagPresMap[opt]=false;
  longFlagDesc[opt]=desc;
}

void ArgParser::addShortOption(char opt, bool argument,std::string desc){
  shortFlagMap[opt]="";
  shortFlagReqArgMap[opt]=argument;
  shortFlagPresMap[opt]=false;
  shortFlagDesc[opt]=desc;
}

void ArgParser::addArgument(std::string name,bool required, std::string desc){
  if(required) {
    reqArgs.push_back(name);
    reqArgDesc[name]=desc;
  }
  else {
    optArgs.push_back(name);
    optArgDesc[name]=desc;
  }
}

void ArgParser::addUnlimitedArgument(std::string name, std::string desc){
  unlimitedArg = name;
  unlimitedArgDesc=desc;
  unlimitedArgMode=true;
}

void ArgParser::printOptions(std::string appName){
  std::cout << "Usage:  " << appName << "  [options] ";
  std::vector<std::string>::const_iterator argIt;
  for(argIt = reqArgs.begin(); argIt != reqArgs.end(); argIt++){
    std::cout << *argIt << " ";
  }
  for(argIt = optArgs.begin(); argIt != optArgs.end(); argIt++){
    std::cout << "[" << *argIt << "] ";
  }
  if(unlimitedArgMode){
    std:: cout << "[" << unlimitedArg << " 1] "
	       << "[" << unlimitedArg << " 2] ...";
  }
  std::cout << std::endl
	    << "Required Arguments: " <<std::endl;
  stringmap::const_iterator smapIt;
  for(smapIt = reqArgDesc.begin(); smapIt !=reqArgDesc.end(); smapIt++){
    printf("\t%-60s%s\n",smapIt->first.c_str(),smapIt->second.c_str());
  }
  std::cout << "Optional Arguments: " << std::endl;
  for(smapIt = optArgDesc.begin(); smapIt !=optArgDesc.end(); smapIt++){
    printf("\t%-60s%s\n",smapIt->first.c_str(),smapIt->second.c_str());
  }
  if(unlimitedArgMode){
    std::cout << "Unlimited Arguemnts: " <<std::endl;
    printf("\t%-60s%s\n",unlimitedArg.c_str(),unlimitedArgDesc.c_str());
  }
  std::cout << "Options: " << std::endl;
  std::map<char,std::string>::const_iterator cmapIt;
  for(cmapIt = shortFlagDesc.begin(); cmapIt != shortFlagDesc.end(); cmapIt++){
    printf("\t-%-59c%s\n",cmapIt->first,cmapIt->second.c_str());
  }
  for(smapIt = longFlagDesc.begin(); smapIt !=longFlagDesc.end(); smapIt++){
    printf("\t--%-58s%s\n",smapIt->first.c_str(),smapIt->second.c_str());
  }
  
}

int ArgParser::process(std::string& ret){
  status = this->internalProcess(ret);
  return status;
}

int ArgParser::internalProcess(std::string& ret){
  /* Return Codes:
     0   - good
     -1  - invalid argument
     -2  - undeclared option
     -3  - option does not have required argument
     -4  - option has unexpected argument
     -6  - required argument missing
     -7  - too many arguments
  */
  for(int i=1;i<argc;i++){
    std::string thisString(argv[i]);
    
    if(thisString[0]!='-'){ //not an option, push it on the arg list
      inputArgs.push_back(thisString);
      continue;
    }
    if(thisString.length()<2) {ret=thisString; return -1;} // '-' is an invalid entry
    
    if(thisString[1]=='-'){ //check for a -- style argument
      size_t found = thisString.find_first_of('=');
      if(found!=std::string::npos){ // the argument is of the form --arg=val
	std::string arg = thisString.substr(2,found-2); //strip the '--'
	std::string val = thisString.substr(found+1,std::string::npos);
	if(longFlagMap.find(arg) == longFlagMap.end()){
	  ret =  arg;
	  return -2; 
	}
	if(val.compare("")==0 && longFlagReqArgMap[arg]){
	  ret = arg;
	  return -3;
	}
	if(val.compare("")!=0 && !longFlagReqArgMap[arg]){
	  ret = arg;
	  return -4;
	}
	longFlagMap[arg]=val;
	longFlagPresMap[arg]=true;
	continue;
      }else{ // ok, check the next entry in the arg table 
	std::string arg = thisString.substr(2,std::string::npos);
	if(longFlagMap.find(arg) == longFlagMap.end()){ ret = arg; return -2;}
	longFlagPresMap[arg]=true;
	if(!longFlagReqArgMap[arg]){ //maybe it doesn't need an argument:
	  longFlagMap[arg]="";
	  continue;
	}
	//it does need an argument
	if(i==argc-1){ //there is no next entry!
	    ret = arg;
	    return -3;
	}else{ //there is a next entry!
	  std::string nextString(argv[i+1]);
	  //std::cout << "nextString: " << nextString << std::endl;
	  //make sure it isn't a flag
	  if(nextString[0]=='-'){
	    ret=arg; return -3;
	  }
	  //its not a flag!
	  longFlagMap[arg]=nextString;
	}
	i++; //skip the next entry, since we already used it
	continue;
      }
    }  //Done with long args!
    else{ //short args.  These can only have the format -fValue or -f Value
      char flag=thisString[1];
      if(shortFlagMap.find(flag) == shortFlagMap.end()) {ret = std::string(&flag); return -2;}
      shortFlagPresMap[flag]=true;
      if(thisString.length()==2){ //its of the form -f or -f Value
	if(!shortFlagReqArgMap[flag]) continue; //doesn't need an arg, happy
	//ok, it does need an argument, get next string
	if(i==argc-1){ret=std::string(&flag); return -3;} //needs an argument, but didn't get one
	std::string nextString(argv[i+1]);
	if(nextString[0]=='-'){ret=std::string(&flag); return -3;}
	shortFlagMap[flag]=nextString;
	i++; // skip the next argument since we already used it
      }else{  //ok, its of the form -fValue
	if(!shortFlagReqArgMap[flag]) {ret=std::string(&flag); return -4;}
	std::string val = thisString.substr(2,std::string::npos);
	shortFlagMap[flag]=val;
      }
    }//done with short args!
  } //done with loop

  //check required arguments

  if(inputArgs.size() < reqArgs.size()) return -6;
  if(!unlimitedArgMode) if(inputArgs.size() > reqArgs.size()+optArgs.size()) return -7;

  return 0;
}

bool ArgParser::shortFlagPres(char f){
  return getWithCheck<char,bool>(f,shortFlagPresMap,false);
}
bool ArgParser::longFlagPres(std::string f){
  return getWithCheck<std::string,bool>(f,longFlagPresMap,false);
}

std::string ArgParser::getLongFlag(std::string arg){
  return getWithCheck<std::string,std::string>(arg,longFlagMap,"");
}
std::string ArgParser::getShortFlag(char arg){
  return getWithCheck<char,std::string>(arg,shortFlagMap,"");
}
std::string ArgParser::getArgument(std::string a){
  if(status!=0) return "";
  std::vector<std::string>::const_iterator loc = std::find(reqArgs.begin(),reqArgs.end(),a);
  if(loc != reqArgs.end())
    return inputArgs.at(loc-reqArgs.begin());
  
  loc = std::find(optArgs.begin(),optArgs.end(),a);
  if(loc != optArgs.end()){
    if(loc-optArgs.begin()+reqArgs.size() < inputArgs.size())
      return inputArgs.at(loc-optArgs.begin()+reqArgs.size());
  }
  return "";
  
}

std::vector<std::string>::const_iterator ArgParser::getUnlimitedArgument(std::string a){
  if(!unlimitedArgMode) return inputArgs.end();
  if(inputArgs.size() == reqArgs.size() + optArgs.size()) return inputArgs.end();
  return inputArgs.begin()+reqArgs.size()+optArgs.size();

}

template<typename T,typename S>
S ArgParser::getWithCheck(T key, std::map<T,S> m,S def){
  if(status!=0) return def;
  if(m.find(key) == m.end()) return def;
  return m[key];
}

void ArgParser::reset(){
  longFlagMap.clear();
  shortFlagMap.clear();
  longFlagReqArgMap.clear();
  shortFlagReqArgMap.clear();
  longFlagPresMap.clear();
  shortFlagPresMap.clear();
  reqArgs.clear();
  optArgs.clear();
  inputArgs.clear();
  status = 0;
}
