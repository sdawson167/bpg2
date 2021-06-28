#include <iostream>
#include <iomanip>
#include <utility>
#include <cstring>
#include <fstream>
#include <cmath>
#include <sstream>

#include "utils.h"


/*
 * ======================================================
 *            method to parse input files
 * ======================================================
 *
 *   method takes string, divides it into words 
 *   (separated by spaces) and returns a vector of 
 *   those words 
 *
 */

std::vector<std::string> divideWords(std::string str)
{
  using namespace std;
  string w = "";  // current word
  vector<string> strVec;

  // parse each char in the string:
  for (auto c : str)
  {
    // if we've found a space, divide the sentence here
    if (c == ' ') {
      strVec.push_back(w);
      w = "";
    // otherwise add character to current word
    } else {
      w = w + c;
    }
  }
  strVec.push_back(w);
  return strVec;
}

/* 
 * =================================================================
 * Method to parse input args and produce list of computations to do
 * =================================================================
 */

void parseInputArgs(int argc, char** argv, paramList &phasePoints, std::vector<int> &phaseIDList, bool &resetFlag, std::string &model)
{
  using namespace std;

  // choose model - LB, OK, PW
  
  // if model is not chosen as input arg, give user option to enter it:
  if (argc > 1) model = argv[1];
  else {
    cout << "choose model: LB, OK, or PW" << endl;
    getline(cin, model);
    cin.clear();
  }
                                                                                              
  // verify model choice is valid
  if ( model != "LB" && model != "OK" && model != "PW") {
    string message = "first input arg. is invalid - user must choose LB, OK, or PW";
    throw message;
  }

  // depending on chosen model paramNumber checks that we are given enough parameters
  int paramNumber = 0;
  if (model == "LB") paramNumber = 2;      // LB model: tau, gamma
  else if (model == "OK") paramNumber = 2; // OK model: tau, gamma
  else if (model == "PW") paramNumber = 4; // PW model: tau, gamma, lam1, lam2
                                                                                                                                                     
  // determine form of input arguments - input file (1), command line (2), fixed tau (3) or fixed gamma (4):
  
  // form of input should be entered as second argument:
  int inMode = 0;
  if (argc > 2) inMode = atoi(argv[2]);
  
  // if not, allow user to enter choice:
  else
  {
    string input;
    cout << "choose input mode: 1 (input file), 2 (command line input), 3 (fixed tau bdry mode), 4 (fixed gamma bdry mode)" << endl;
    getline(cin, input);
    
    stringstream stream(input);
    if (!(stream >> inMode)) {
      string message =  "second input arg is invalid - choose 1, 2, 3 or 4";
      cin.clear();
      throw message;
    }

    cin.clear();
  }

  // verify choice of input mode is valid
  if (inMode != 1 && inMode != 2 && inMode != 3 && inMode != 4) {
    string message =  "second input arg is invalid - choose 1, 2, 3 or 4";
    throw message;
  }
  
  /* initialize params from remaining input arguments
   * how we do this depends on choice of input mode:
   * 	1 (input file) - need name of file from cmd line
   * 	2 (command line) - need remaining params from cmd line
   * 	3 (bdry mode, fixed tau) - need remaining params from cmd line
   * 	4 (bdry mode, fixed gamma) - " " 
   */
  
  // use this to determine if phaseIDs were included as command line args
  // or if we should prompt the user to enter them:
  int phaseIDargNo = 0;

  // if we chose to initialize from file:
  if (inMode == 1)
  {
    // get name of file
    string inString, inFileName;
    
    // allow user to enter input file name if not given:
    if (argc < 4) {
      cout << "enter name of input file" << endl;
      getline(cin, inString);
      cin.clear();

    } else 
      inString = argv[3];
      
    inFileName = "/home/sdawson/projects/def-shi/sdawson/BPG/" + inString;
    
    // open file
    ifstream inFile;
    inFile.open(inFileName); 
    
    // verify open:
    if (!inFile) {
      string message = "could not read input file " + inFileName;
      throw message;
    }
    
    // read input params from file and store in paramList object
    string line;
    while (getline(inFile,line)) 
    {
      vector<string> brokenLine = divideWords(line);
      if ((int) brokenLine.size() < paramNumber) {
        string message = "invalid parameter file for " + model + " model, file must provide " + to_string(paramNumber) + " input params per line";
        throw message;
      }
      param otherParams;
      for (int index = 0; index < paramNumber; index++) 
        otherParams.push_back(stof(brokenLine[index]));
      
      if (model == "LB") {
        double lam1 = 0.0;
        double lam2 = 0.0;

        otherParams.push_back(lam1);
        otherParams.push_back(lam2);
      }

      if (model == "OK") {
        double lam1 = 1.0;
        double lam2 = 1.0;

        otherParams.push_back(lam1);
        otherParams.push_back(lam2);
      }

      phasePoints.push_back(otherParams);
    }

    // set reset flag value
    if (argc < 5) {
      cout << "Reset initial condition between runs? Choose 1 (yes) or 0 (no)" << endl; 
      
      string input;
      getline(cin, input);
      stringstream stream(input);
      if (!(stream >> resetFlag)) {
        cin.clear();
        throw "invalid selection for resetFlag";
      }

    } else {
      resetFlag = (bool) atoi(argv[4]);    
    }
    
    // if argc < 6 - need to request phaseIDs from user
    phaseIDargNo = 6; 

  // end inMode = 1 case 
                                                                                                                                                     
  // inMode = 2: read input params (tau, gamma, ...) and list of phaseIDs from command line:
  } else if (inMode == 2) {

    // resetFlag is always true in this case
    resetFlag = 1;
    
    // first, get input params (tau, gamma, ...):
    param otherParams;
    
    // if parameters are not provided, allow user to enter them:
    if (argc < (3 + paramNumber) ) {
      cout << "enter " << paramNumber << " input params required by " << model << " model:" << endl;
      
      double x;
      string input;
      int n = 0;
      while (getline(cin, input)) {
        
        vector<string> dividedInput = divideWords(input);
  
        for (size_t index = 0; index < dividedInput.size(); index++) {
          stringstream stream(dividedInput[index]);
    
          if (stream >> x && n < paramNumber) {
            otherParams.push_back(x);
            n++;
          } 
        } // end loop over input pts.
  
        // if we have enough params we can stop waiting for user input
        if (n == paramNumber)
          break;

      } // end input while loop
      cin.clear();
    
    // otherwise read parameters from command line args:
    } else {
      for (int index = 3; index < (3 + paramNumber); index++) 
        otherParams.push_back(atof(argv[index])); 
    }
    
    // if LB or OK models chosen - set lam1 and lam2 values accordingly
    if (model == "LB") {
      double lam1 = 0.0;
      double lam2 = 0.0;
      
      otherParams.push_back(lam1);
      otherParams.push_back(lam2);
    }
    if (model == "OK") {
      double lam1 = 1.0;
      double lam2 = 1.0;

      otherParams.push_back(lam1);
      otherParams.push_back(lam2);
    }

    phasePoints.push_back(otherParams);
    
    // if argc < 4 + paramNumber - need to prompt user to enter phaseIDs:
    phaseIDargNo = paramNumber + 4;    
  
  // end inMode = 2 case

  // inMode = 3 : bdry mode with fixed tau val. 
  } else if (inMode == 3) {
    // resetFlag is always false in this case
    resetFlag = 0;

    // get fixed tau val
    double tau;

    // if tau val not provided - allow user to enter it
    if (argc < 4) {
      cout << "choose tau value" << endl;
      string input;
      getline(cin, input);

      stringstream stream(input);
      if (!(stream >> tau)) {
        string message = "invalid choice for tau parameter";
        cin.clear();
        throw message;
      }

    } else {
      tau = atof(argv[3]);
    }

    // get range of gamma values
    double gamma1 = 0.0;
    double gamma2 = 0.0;
    
    // if parameters are not provided, allow user to enter them:
    if (argc < 6 ) {
      cout << "choose start and end vals of gamma:" << endl;
      
      double x;
      string input;
      int n = 0;
      while (getline(cin, input)) {
        
        vector<string> dividedInput = divideWords(input);
                                                                                                     
        for (size_t index = 0; index < dividedInput.size(); index++) {
          stringstream stream(dividedInput[index]);
    
          if (stream >> x && n < 2) {
            if (n == 0)
              gamma1 = x;
            if (n == 1)
              gamma2 = x;
            n++;
          } 
        } // end loop over input pts.
                                                                                                     
        // if we have enough params we can stop waiting for user input
        if (n >= 2)
          break;
                                                                                                     
      } // end input while loop
      cin.clear();
    
    // otherwise read parameters from command line args:
    } else {
      gamma1 = atof(argv[4]);
      gamma2 = atof(argv[5]);
    }

    // get lam1 and lam2 vals
    double lam1 = 0.0;
    double lam2 = 0.0;
    if (model == "OK") {
      lam1 = 1.0;
      lam2 = 1.0;
    }

    if (model == "PW") {
      // if parameters are not provided, allow user to enter them:
      if (argc < 8) {
        cout << "choose lam1 and lam2 values" << endl;
        
        double x;
        string input;
        int n = 0;
        while (getline(cin, input)) {
          
          vector<string> dividedInput = divideWords(input);
                                                                                                       
          for (size_t index = 0; index < dividedInput.size(); index++) {
            stringstream stream(dividedInput[index]);
      
            if (stream >> x && n < 2) {
              if (n == 0)
                lam1 = x;
              if (n == 1)
                lam2 = x;
              n++;
            } 
          } // end loop over input pts.
                                                                                                       
          // if we have enough params we can stop waiting for user input
          if (n >= 2)
            break;
                                                                                                       
        } // end input while loop
        cin.clear();
      
      // otherwise read parameters from command line args:
      } else {
        lam1 = atof(argv[6]);
        lam2 = atof(argv[7]);
      }
    }

    // create range of gamma values - add (tau, gamma, lam1, lam2) lists to
    // paramlist vector 
    const double dGamma = 0.01;

    int reverseFlag = gamma1 < gamma2 ? 1 : -1;
    int nGamma = (int) (abs(gamma2 - gamma1) / dGamma);
    
    for (int i = 0; i < nGamma; i++) {
      // add tau val to param list:
      param otherParams;
      otherParams.push_back(tau);

      // compute next gamma val and add: 
      double gamma = gamma1 + i * reverseFlag * dGamma;
      otherParams.push_back(gamma);

      // add lam1 and lam2 vals to param list:
      otherParams.push_back(lam1);
      otherParams.push_back(lam2);
       
      // add point to list of phasePoints
      phasePoints.push_back(otherParams);
    }
     
    // if argc < 5 + paramNumber - need to request phaseIDs from user
    phaseIDargNo = 5 + paramNumber;

  // end inMode = 3 case
  
  // inMode = 4 : bdry mode with fixed gamma  val. 
  } else if (inMode == 4) {
    // resetFlag is always false in this case
    resetFlag = 0;

    // get range of tau values
    double tau1 = 0.0;
    double tau2 = 0.0;
    
    // if parameters are not provided, allow user to enter them:
    if (argc < 5 ) {
      cout << "choose start and end vals of tau:" << endl;
      
      double x;
      string input;
      int n = 0;
      while (getline(cin, input)) {
        
        vector<string> dividedInput = divideWords(input);
                                                                                                     
        for (size_t index = 0; index < dividedInput.size(); index++) {
          stringstream stream(dividedInput[index]);
    
          if (stream >> x && n < 2) {
            if (n == 0)
              tau1 = x;
            if (n == 1)
              tau2 = x;
            n++;
          } 
        } // end loop over input pts.
                                                                                                     
        // if we have enough params we can stop waiting for user input
        if (n >= 2)
          break;
                                                                                                     
      } // end input while loop
      cin.clear();
    
    // otherwise read parameters from command line args:
    } else {
      tau1 = atof(argv[3]);
      tau2 = atof(argv[4]);
    }

    // get fixed gamma val
    double gamma;
                                                             
    // if tau val not provided - allow user to enter it
    if (argc < 6) {
      cout << "choose gamma value" << endl;
      string input;
      getline(cin, input);
                                                             
      stringstream stream(input);
      if (!(stream >> gamma)) {
        string message = "invalid choice for tau parameter";
        cin.clear();
        throw message;
      }
                                                             
    } else {
      gamma = atof(argv[5]);
    }

    // get lam1 and lam2 vals
    double lam1 = 0.0;
    double lam2 = 0.0;
    if (model == "OK") {
      lam1 = 1.0;
      lam2 = 1.0;
    }

    if (model == "PW") {
      // if parameters are not provided, allow user to enter them:
      if (argc < 8) {
        cout << "choose lam1 and lam2 values" << endl;
        
        double x;
        string input;
        int n = 0;
        while (getline(cin, input)) {
          
          vector<string> dividedInput = divideWords(input);
                                                                                                       
          for (size_t index = 0; index < dividedInput.size(); index++) {
            stringstream stream(dividedInput[index]);
      
            if (stream >> x && n < 2) {
              if (n == 0)
                lam1 = x;
              if (n == 1)
                lam2 = x;
              n++;
            } 
          } // end loop over input pts.
                                                                                                       
          // if we have enough params we can stop waiting for user input
          if (n >= 2)
            break;
                                                                                                       
        } // end input while loop
        cin.clear();
      
      // otherwise read parameters from command line args:
      } else {
        lam1 = atof(argv[6]);
        lam2 = atof(argv[7]);
      }
    }

    // create range of gamma values - add (tau, gamma, lam1, lam2) lists to
    // paramlist vector 
    const double dTau = 0.01;

    int reverseFlag = tau1 < tau2 ? 1 : -1;
    int nTau = (int) (abs(tau2 - tau1) / dTau);
    
    for (int i = 0; i < nTau; i++) {
      param otherParams;

      // compute tau value and add:
      double tau = tau1 + i * reverseFlag * dTau;
      otherParams.push_back(tau);

      // add gamma value. 
      otherParams.push_back(gamma);

      // add lam1 and lam2 vals to param list:
      otherParams.push_back(lam1);
      otherParams.push_back(lam2);
       
      // add point to list of phasePoints
      phasePoints.push_back(otherParams);
    }
     
    // if argc < 5 + paramNumber - need to request phaseIDs from user
    phaseIDargNo = 5 + paramNumber;

  } // end inMode = 4 case

 // now get phaseIDs - list of phases to optimize:
 if (argc <  phaseIDargNo) {
   cout << "choose phases (1 - 9) to optimize, enter any non-numerical character to finish:" << endl;
   
   int phaseID;
   string input;
   bool readInFlag = true;
   while (getline(cin, input)) {
 
     vector<string> dividedInput = divideWords(input);
 
     for (size_t index = 0; index < dividedInput.size(); index++) {
       stringstream stream(dividedInput[index]);
                                                                                                            
       if (stream >> phaseID) {
         if (phaseID > 0 && phaseID < 10)
           phaseIDList.push_back(phaseID);
       } else 
         readInFlag = false;
     }
     if (!readInFlag)
       break;
   }
                                                                                                             
   // make sure we have at least one phase:  
   if (phaseIDList.size() == 0) {
     string message =  "user must choose 1 or more phases to optimize";
     cin.clear();
     throw message;
   }
   cin.clear();

   cout << "starting optimization" << endl;
                                                                                                            
 } else {
   for (int p = (phaseIDargNo - 1); p < argc; p++) 
     phaseIDList.push_back(atoi(argv[p])); // end parseInputArgs method 
 }
}

double getGamma0(int phaseID) {
  double gamma0 = 0.0;

  if (phaseID == 4) gamma0 = 2.0;
 
  if (phaseID == 5) gamma0 = 2.0;

  if (phaseID == 6) gamma0 = 2.0; 

  if (phaseID == 7) gamma0 = 2.0;

  return gamma0;
}
