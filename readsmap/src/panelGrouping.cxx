#include "zcompress.h"

//g++ -g -O2 -march=nocona -pipe -ffast-math  -D_LARGEFILE_SOURCE -D_LARGEFILE64_SOURCE -c src/zcompress.cxx
//g++ -g -O2 -march=nocona -pipe -ffast-math  -D_LARGEFILE_SOURCE -D_LARGEFILE64_SOURCE src/panelGrouping.cxx zcompress.o -o ./panelGrouping -lz

using std::cout;
using std::cerr;

#define ABS_MAX_PANEL_ID 1900800700

void panelGroupingUsage( string arg1, bool is_short ) {
  cerr << "This indexes a csFasta reads/map file if it hasn't been already\n"
       << "indexed, and then outputs (to standard output) on each line the\n"
       << "begining and end panel for each group.  Thus there will be\n"
       << "'numSplits' groups total.\n\n";

  cerr << "Usage for " << arg1 << "\n"
       << "   panelGrouping csFasta[,csFasta2] numSplits [--minPanelNum xxx] [--maxPanelNum xxx]\n"
       << "                 [--max-number-of-panels xxx]\n"
       << "                 [--print-all-group-members] [--print-file-positions]\n"
       << "                 [--extract-random-subset xxx]\n"
       << "                 [--extractGroup]\n";
  if (is_short) exit(1);

  cerr << "\nRequired:\n"
       << "   csFasta is the reads or map file\n"
       << "       Note, if csFasta2 is specified, grouping will only be done on\n"
       << "             panels shared by both.  (Good for paired tag match files.)\n"
       << "   numSplits is the number of panel groups (files) to split into.\n"
       << "Optional:\n"
       << "   --minPanelNum xxx:\n"
       << "       min panel number [default is to have no min]\n"
       << "   --maxPanelNum xxx:\n"
       << "       max panel number [default is to have no max]\n"
       << "   --print-all-group-members:\n"
       << "       each line has all panels in that group [def. just start and end]\n"
       << "   --print-file-positions:\n"
       << "       print file position of group starts (inclusive) and ends (exclusive)\n"
       << "   --max-number-of-panels xxx:\n"
       << "       max. number of panels to use [defaults to using all]\n"
       << "   --extract-random-subset xxx:\n"
       << "       extracts a random subset of panels.\n"
       << "       Note, with this option, print-all-group-members will be forced on\n"
       << "   --extractGroup\n"
       << "       for a single csFasta file, output group that starts and ends at\n"
       << "       particular bead ids.(min/maxPanelNum required and numSplits ignored)\n"
    //  << "           --diff-random-subsets\n"
    //  << "               If multiple files are specified use different random subsets\n"
    //  << "               for each file. [defaults to using the same one]\n"
       << "\n";
  exit(1);
}

vector<string> panelGroupingParsing( int argc, char **argv ) {
  vector<string> cmdLineArgs;
  int argNum;
  for(argNum=0; argNum<argc; argNum++) {
    string curArg(argv[argNum]);
    cmdLineArgs.push_back(curArg);
  }
  return cmdLineArgs;
}

class Parameter_value {
public:
  Parameter_value() {};
  
  void set_long_parameter(string parm) {
    long_parameter = parm.substr(2,parm.size()-2); };
  void set_short_parameter(string parm) {
    short_parameter = parm.substr(1,parm.size()-1); };
  void set_value(string v_in) { value = v_in; };
  string get_value() const { return value; };
  bool operator== (const string &inParameter) const {
    return inParameter == long_parameter || inParameter == short_parameter; };
  
private:
  string long_parameter;
  string short_parameter;
  string value;
};

vector<Parameter_value> extractParameters( vector<string> &cmdLineArgsOrig) {
  vector<Parameter_value> param_values;

  vector<string> cmdLineArgs = cmdLineArgsOrig;
  cmdLineArgsOrig.clear();
  vector<string>::iterator cmdLineArgsItr = cmdLineArgs.begin();
  for(;cmdLineArgsItr!=cmdLineArgs.end();) {
    if (*cmdLineArgsItr == "--") {
      return param_values;  //empty set
    } else if (cmdLineArgsItr->substr(0,1) == "-" && cmdLineArgsItr->size() > 1) {
      Parameter_value this_param;
      this_param.set_value("1"); //default value if not set below

      if (cmdLineArgsItr->substr(0,2) == "--")
	this_param.set_long_parameter(*cmdLineArgsItr);
      else
	this_param.set_short_parameter(*cmdLineArgsItr);

      cmdLineArgsItr++;

      if (cmdLineArgsItr!=cmdLineArgs.end()) {
	if ((cmdLineArgsItr)->substr(0,1) != "-") {
	  this_param.set_value(*cmdLineArgsItr);
	  cmdLineArgsItr++;
	}
      }
      param_values.push_back(this_param);
    } else {
      cmdLineArgsOrig.push_back(*cmdLineArgsItr);
      cmdLineArgsItr++;
    }
  }
  //cout << cmdLineArgsOrig.size() << " is size.\n";
  return param_values;
}

string getParamValue( string paramName, const vector<Parameter_value> param_values ) {
  vector<Parameter_value>::const_iterator param_valuesItr = 
    find(param_values.begin(), param_values.end(),paramName);
  string paramValue;
  if(param_valuesItr!=param_values.end())
    paramValue = param_valuesItr->get_value().c_str();
  return paramValue;
}

int getIntegerParamValue( string paramName, const vector<Parameter_value> param_values,
			  int default_value) {
  int returnValue = default_value;
  string paramValue = getParamValue(paramName,param_values);
  if (paramValue.size() == 0 )
    return returnValue;
  try {
    char * last = NULL;
    if(paramValue.size() != 0 )
      returnValue = (int)strtol(paramValue.c_str(),&last,10);
    if(strcmp(last,"")) {
      throw 0;
    }
  } catch (...) {
    cerr << "Value (" << paramValue << ") given to " << paramName
	 << " was not an integer.  Exiting.\n";
    exit(-1);
  }
  return returnValue;
}

bool haveParam( string paramName, const vector<Parameter_value> param_values ) {
  //cout << paramName << param_values[0].get_value() <<" " << <<" haveParam\n";
  string printAllGroupMembersValue = getParamValue(paramName,param_values);
  return printAllGroupMembersValue.size() != 0;
}


int main(int argc, char **argv)
{

  srand ( unsigned ( time (NULL) ) );

  vector<string> cmdLineArgs = panelGroupingParsing(argc,argv);
  vector<Parameter_value> param_values = extractParameters(cmdLineArgs);
  if(haveParam("help",param_values)) panelGroupingUsage(*cmdLineArgs.begin(),false);

 
  int minPanelNum = getIntegerParamValue("minPanelNum",param_values,-1);
  int maxPanelNum = getIntegerParamValue("maxPanelNum",param_values,-1);
  if(minPanelNum<=-2 || maxPanelNum<=-2 || (maxPanelNum!=-1 && maxPanelNum<minPanelNum)) {
    cerr << "invalid min/max panelNumber" << "\n";
    exit(1);
  }

  int max_number_of_panels = getIntegerParamValue("max-number-of-panels",param_values,-1);
  bool do_printAllGroupMembers = haveParam("print-all-group-members",param_values);
  int randomSubsetSize = getIntegerParamValue("extract-random-subset",param_values,-1);
  //bool diffRandSubset = haveParam("diff-random-subsets",param_values);
  bool printFilePos = haveParam("print-file-positions",param_values);

  bool doExtractGroup = haveParam("extractGroup",param_values);
  if(doExtractGroup) {
    if(minPanelNum==-1 && maxPanelNum==-1) {
      cerr << "Either min or max (or both) panel ids are required when doing 'extractGroup'.\n";
      exit(1);
    }
  }

  if(randomSubsetSize!=-1) do_printAllGroupMembers = true;

  //cout << "min/max " << minPanelNum << " " << maxPanelNum << "\n";
  if( cmdLineArgs.size() != 3 )
    panelGroupingUsage(*cmdLineArgs.begin(),true);

  vector<string> matchFilenames = split_string(cmdLineArgs[1],",");
  if(doExtractGroup && matchFilenames.size()!=1) {
    cerr << "For 'extractGroup', can only do one csfasta at one time.\n";
    exit(1);
  }

  if (maxPanelNum==-1) maxPanelNum=ABS_MAX_PANEL_ID;

  int numSplits = atoi(cmdLineArgs[2].c_str());
  if (numSplits<=0) {
    cerr << "ERROR! '" << cmdLineArgs[2] << "' is not a valid for number of splits.\n";
    exit(-1);
  }

  vector<PanelIndexing*> panelGrpingSet;
  for(int i=0;i<(int)matchFilenames.size();i++) {
    PanelIndexing *panelGrping = new PanelIndexing();
    panelGrpingSet.push_back(panelGrping);
  }
  vector<int> rand_panel_ids;  //only select random set for first file
  //then use the same random set for all the others
  vector<string>::const_iterator matchFnItr = matchFilenames.begin();
  vector<PanelIndexing*>::iterator panelGrpingItr = panelGrpingSet.begin();
  for(; matchFnItr != matchFilenames.end(); ) {
    // check match file by openining it
    cerr << "Working with file " << *matchFnItr << ".\n";
    genFile matchFile;
    bool didOpen = matchFile.open(*matchFnItr,"r");
    if (!didOpen) {
      cerr << "Cannot open match file.\n";
      exit(1);
    }

    if(minPanelNum != -1 || maxPanelNum != -1)
      (*panelGrpingItr)->setMinMaxPanelId(minPanelNum,maxPanelNum);
    if (max_number_of_panels != -1 && matchFilenames.size()== 1)
      (*panelGrpingItr)->setMaxNumberPanels(max_number_of_panels);

    // Checks if index is there, if not makes it.
    // Also class allows for grouping of panels
    if (!doExtractGroup)  //extracting group, no need for index
      //improvement for future, can build index while extracting at the same time.
      (*panelGrpingItr)->ensureIndexFile(matchFile); 

    if (matchFnItr != matchFilenames.begin()) {  //intersection with the last one
      (*panelGrpingItr)->intersection(**(panelGrpingItr-1));
    } else if (doExtractGroup) { //Do this if you are streaming out part of a csfasta file.
      (*panelGrpingItr)->streamOutReadsFile(matchFile,minPanelNum,maxPanelNum);
      return 0;
    }
    panelGrpingItr++;
    matchFnItr++;
  }
  panelGrpingItr--;

  if (randomSubsetSize != -1 ) {
    //if (matchFnItr == matchFilenames.begin() || diffRandSubset)
    rand_panel_ids = (*panelGrpingItr)->reduceToRandomSubset(randomSubsetSize);
    //else
    //(*panelGrpingItr)->reduceToSubset(rand_panel_ids);
  }

  if (max_number_of_panels != -1 && matchFilenames.size() != 1) {
    (*panelGrpingItr)->reduceToSubset(max_number_of_panels);
  }

  cerr << "Now finding starts and ends.\n";
  (*panelGrpingItr)->determinePanelStartsEnds(numSplits);

  if (do_printAllGroupMembers)
    (*panelGrpingItr)->printAllPanelsInGroup(std::cout, printFilePos);
  else
    (*panelGrpingItr)->printStartEndGrouping(std::cout, printFilePos);
  
  //should clean up new's from above here.
  return 0;
}
