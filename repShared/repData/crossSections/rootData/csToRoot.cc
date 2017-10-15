// implement all (known) pbar ionisation cross sections with gas in root file 
// by ntuple, histograms and graph.
//
// Pierre Grandemange
// 29/11/2014
// last modif : 18/02/2015


//data from : "J. Phys. B: At. Mol. Opt. Phys. 44 (2011) 122001"

// standard library
#include <iostream>
#include <iomanip> 
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
// root classes
#include "TFile.h"
#include "TNtuple.h"
#include "TGraphErrors.h"
#include "TMath.h"

using namespace std;

// Prototypes ******************************************************************
void ntupleRaw3(string const& fileName);
void ntupleRaw2(string const& fileName);
void ntupleProcess3(string const& fileName, string const& fileNameMain);
void processAll();
void setCSsum(string const& gasName);

// Main ************************************************************************
int main()
{
  cout << endl;
  cout << "******************" << endl;
  cout << "Main program start" << endl;
  cout << "******************" << endl;
  cout << endl;
  
  processAll();
  setCSsum("H2");
  setCSsum("CO");
  setCSsum("CO2");
  setCSsum("CH4");
  setCSsum("N2");
  setCSsum("Ne");
  setCSsum("Ar");
  setCSsum("Kr");
  setCSsum("Xe");
  
  cout << endl;
  cout << "****************" << endl;
  cout << "Main program end" << endl;
  cout << "****************" << endl;
  cout << endl;
  return 0;
}
// *****************************************************************************


// Functions *******************************************************************
// *****************************************************************************

// processAll ******************************************************************
void processAll()
{
// _____________________________________________________________________________
// input / output files
  int nbFiles(39);  // nb of cross section files to implement
  string fName[nbFiles]; // names of the files
  // molecular H2
  fName[0]="data/CS_pbar-H2--H2plus_Hvelplund.dat";
  fName[1]="data/CS_pbar-H2--Hplus_Hvelplund.dat";
  fName[2]="data/CS_pbar-H2--H2plus_Andersen.dat";
  fName[3]="data/CS_pbar-H2--H2plus_Lhur.dat";
  // CO
  fName[4]="data/CS_pbar-CO--COplus_Knudsen.dat";
  fName[5]="data/CS_pbar-CO--Oplus_Knudsen.dat";
  fName[6]="data/CS_pbar-CO--COplusplus_Knudsen.dat";
  fName[7]="data/CS_pbar-CO--Cplus_Knudsen.dat";
  fName[8]="data/CS_pbar-CO--Oplusplus_Knudsen.dat";
  fName[9]="data/CS_pbar-CO--Cplusplus_Knudsen.dat";
  // CO2
  fName[10]="data/CS_pbar-CO2--CO2plus_Knudsen.dat";
  fName[11]="data/CS_pbar-CO2--COplus_Knudsen.dat";
  fName[12]="data/CS_pbar-CO2--CO2plusplus_Knudsen.dat";
  fName[13]="data/CS_pbar-CO2--Oplus_Knudsen.dat";
  fName[14]="data/CS_pbar-CO2--Cplus_Knudsen.dat";
  fName[15]="data/CS_pbar-CO2--Oplusplus_Knudsen.dat";
  fName[16]="data/CS_pbar-CO2--Cplusplus_Knudsen.dat";
  // CH4
  fName[17]="data/CS_pbar-CH4--CH4plus_Knudsen.dat";
  fName[18]="data/CS_pbar-CH4--CH3plus_Knudsen.dat";
  fName[19]="data/CS_pbar-CH4--CH2plus_Knudsen.dat";
  fName[20]="data/CS_pbar-CH4--CHplus_Knudsen.dat";
  fName[21]="data/CS_pbar-CH4--Cplus_Knudsen.dat";
  fName[22]="data/CS_pbar-CH4--H2plus_Knudsen.dat";
  fName[23]="data/CS_pbar-CH4--Hplus_Knudsen.dat";
  // N2
  fName[24]="data/CS_pbar-N2--N2plus_Knudsen.dat";
  fName[25]="data/CS_pbar-N2--Nplus_Knudsen.dat";
  fName[26]="data/CS_pbar-N2--Nplusplus_Knudsen.dat";
  // Ne
  fName[27]="data/CS_pbar-Ne--Neplus_Knudsen.dat";
  fName[28]="data/CS_pbar-Ne--Neplusplus_Knudsen.dat";
  fName[29]="data/CS_pbar-Ne--Neplusplusplus_Knudsen.dat";
  // Ar
  fName[30]="data/CS_pbar-Ar--Arplus_Knudsen.dat";
  fName[31]="data/CS_pbar-Ar--Arplusplus_Knudsen.dat";
  fName[32]="data/CS_pbar-Ar--Arplusplusplus_Knudsen.dat";
  // Kr
  fName[33]="data/CS_pbar-Kr--Krplus_Knudsen.dat";
  fName[34]="data/CS_pbar-Kr--Krplusplus_Knudsen.dat";
  fName[35]="data/CS_pbar-Kr--Krplusplusplus_Knudsen.dat";
  // Xe
  fName[36]="data/CS_pbar-Xe--Xeplus_Knudsen.dat";
  fName[37]="data/CS_pbar-Xe--Xeplusplus_Knudsen.dat";
  fName[38]="data/CS_pbar-Xe--Xeplusplusplus_Knudsen.dat";
  
  TFile tfile("tfPbarCrossSections.root", "RECREATE", 
  "pbar/gas cross sections");

// _____________________________________________________________________________ 
// processes
// pbar + H2
  cout << "*******************" << endl;
  cout << "pbar + H2 reactions" << endl;
  cout << "*******************" << endl;
  ntupleRaw3(fName[0]);                   // pbar + H2 -> pbar + H2+ + e-
  ntupleRaw3(fName[1]);                   // pbar + H2 -> pbar + H+ + e-
  ntupleRaw3(fName[2]);                   // pbar + H2 -> pbar + H2+ + e-
  ntupleRaw2(fName[3]);                   // pbar + H2 -> pbar + H2+ + e-

// pbar + CO
  cout << endl;
  cout << endl;
  cout << "*******************" << endl;
  cout << "pbar + CO reactions" << endl;
  cout << "*******************" << endl;
  ntupleRaw3(fName[4]);                   // pbar + CO -> pbar + CO+ + e-
  ntupleProcess3(fName[5], fName[4]);     // pbar + CO -> pbar + O+ + e-
  ntupleProcess3(fName[6], fName[4]);     // pbar + CO -> pbar + CO++ + e-
  ntupleProcess3(fName[7], fName[4]);     // pbar + CO -> pbar + C+ + e-
  ntupleProcess3(fName[8], fName[4]);     // pbar + CO -> pbar + O++ + e-
  ntupleProcess3(fName[9], fName[4]);     // pbar + CO -> pbar + C++ + e-
  
// pbar + CO2
  cout << endl;
  cout << endl;
  cout << "********************" << endl;
  cout << "pbar + CO2 reactions" << endl;
  cout << "********************" << endl;
  ntupleRaw3(fName[10]);                  // pbar + CO2 -> pbar + CO2+ + e-
  ntupleProcess3(fName[11], fName[10]);   // pbar + CO2 -> pbar + CO+ + e-
  ntupleProcess3(fName[12], fName[10]);   // pbar + CO2 -> pbar + CO2++ + e-
  ntupleProcess3(fName[13], fName[10]);   // pbar + CO2 -> pbar + O+ + e-
  ntupleProcess3(fName[14], fName[10]);   // pbar + CO2 -> pbar + C+ + e-
  ntupleProcess3(fName[15], fName[10]);   // pbar + CO2 -> pbar + O++ + e-
  ntupleProcess3(fName[16], fName[10]);   // pbar + CO2 -> pbar + C++ + e-
  
// pbar + CH4
  cout << endl;
  cout << endl;
  cout << "********************" << endl;
  cout << "pbar + CH4 reactions" << endl;
  cout << "********************" << endl;
  ntupleRaw3(fName[17]);                  // pbar + CH4 -> pbar + CH4+ + e-
  ntupleProcess3(fName[18], fName[17]);   // pbar + CH4 -> pbar + CH3+ + e-
  ntupleProcess3(fName[19], fName[17]);   // pbar + CH4 -> pbar + CH2+ + e-
  ntupleProcess3(fName[20], fName[17]);   // pbar + CH4 -> pbar + CH+ + e-
  ntupleProcess3(fName[21], fName[17]);   // pbar + CH4 -> pbar + C+ + e-
  ntupleProcess3(fName[22], fName[17]);   // pbar + CH4 -> pbar + H2+ + e-
  ntupleProcess3(fName[23], fName[17]);   // pbar + CH4 -> pbar + H+ + e-
  
// pbar + N2
  cout << endl;
  cout << endl;
  cout << "*******************" << endl;
  cout << "pbar + N2 reactions" << endl;
  cout << "*******************" << endl;
  ntupleRaw3(fName[24]);                  // pbar + N2 -> pbar + N2+ + e-
  ntupleProcess3(fName[25], fName[24]);   // pbar + N2 -> pbar + N+ + e-
                                          //            & pbar + N2++ + e-
  ntupleProcess3(fName[26], fName[24]);   // pbar + N2 -> pbar + N2++ + e-
  
// pbar + Ne
  cout << endl;
  cout << endl;
  cout << "*******************" << endl;
  cout << "pbar + Ne reactions" << endl;
  cout << "*******************" << endl;
  ntupleRaw3(fName[27]);                  // 
  ntupleRaw3(fName[28]);                  // 
  ntupleRaw3(fName[29]);                  // 
  
  // pbar + Ar
  cout << endl;
  cout << endl;
  cout << "*******************" << endl;
  cout << "pbar + Ar reactions" << endl;
  cout << "*******************" << endl;
  ntupleRaw3(fName[30]);                  // 
  ntupleRaw3(fName[31]);                  // 
  ntupleRaw3(fName[32]);                  // 
  
  // pbar + Kr
  cout << endl;
  cout << endl;
  cout << "*******************" << endl;
  cout << "pbar + Kr reactions" << endl;
  cout << "*******************" << endl;
  ntupleRaw3(fName[33]);                  // 
  ntupleRaw3(fName[34]);                  // 
  ntupleRaw3(fName[35]);                  // 
  
  // pbar + Xe
  cout << endl;
  cout << endl;
  cout << "*******************" << endl;
  cout << "pbar + Xe reactions" << endl;
  cout << "*******************" << endl;
  ntupleRaw3(fName[36]);                  // 
  ntupleRaw3(fName[37]);                  // 
  ntupleRaw3(fName[38]);                  // 
                                          
  tfile.Close();
}

// ntupleRaw3 ******************************************************************
void ntupleRaw3(string const& fileName)
{  
// declarations
  string fPath;         // entire pathname of file
  string ntName;        // name of ntuple to be crated
  string ntTitle;       // title of ntuple to be created
  string author;        // author of data
  string reaction;      // reaction involved
  string comment;       // extra comment
  std::size_t pos1(0);  // position into a string
  
  ifstream inFile;      // stream for input file
 
// intialisation 
  string dataPath="../";    // location of data directory (relative path)
  fPath=dataPath+fileName;
  
// check if file exist
  inFile.open(fPath.c_str());
  if(!inFile)
  {
    cout << endl;
    cout << "file : " << fileName << " not found!" << endl; 
    cout << "quit ntuple raw subroutine" << endl;
    cout << endl;
  }
  else
  {  
// define ntuple name, reaction, author and comment
    // NTuple name
    pos1    = fileName.find("_");
    ntName  = fileName.substr(pos1+1);      // ntuple name from file name
    // author
    pos1    = ntName.find("_");
    author  = ntName.substr(pos1+1);        // author name + comment
    // reaction
    reaction  = ntName.substr(0, pos1);
    ntName    = ntName.substr(5, ntName.size()-9); // clear "pbar" ...
    author    = author.substr(0, author.size()-4); // ... and ".dat"
    pos1      = reaction.find("-");
    reaction.replace(pos1, 1, " + ");
    pos1      = reaction.find("--");
    reaction.replace(pos1, 2, " --> ");
    // comment on reaction from file name
    pos1      = author.find("_");
    comment   = "";
    if(pos1 != string::npos)  // extract comment if there is comment!
    {
      comment = author.substr(pos1+1);
      author  = author.substr(0, pos1);
      ntTitle = reaction + " (" + author + " - " + comment + ")";
    }
    else
    {
      ntTitle = reaction + " (" + author + ")";
    }
// dump properties
    cout << "__________________________________________________________________" 
    << endl;
    cout << "file     : " << fileName << endl;
    cout << "reaction : " << reaction << endl;
    cout << "author   : " << author << endl;
    cout << "comment  : " << comment << endl;
    cout << "ntName   : " << ntName << endl;
    cout << "ntTitle  : " << ntTitle << endl;
    cout << "__________________________________________________________________" 
    << endl;
    
// write Ntuple and graph
  string varlist("energy:cross_section:sigma");
  TNtuple ntuple(ntName.c_str(), ntTitle.c_str(), varlist.c_str());
  ntuple.ReadStream(inFile);
  ntuple.Write();
  TGraphErrors tge(fPath.c_str(), "%lg %lg %lg");
  tge.SetTitle(ntTitle.c_str());
  tge.SetName(ntName.c_str());
  tge.Write();
  }
// last statements
  inFile.close();
}

// ntupleRaw2 ******************************************************************
void ntupleRaw2(string const& fileName)
{  
// declarations
  string fPath;         // entire pathname of file
  string ntName;        // name of ntuple to be crated
  string ntTitle;       // title of ntuple to be created
  string author;        // author of data
  string reaction;      // reaction involved
  string comment;       // extra comment
  std::size_t pos1(0);  // position into a string
  
  ifstream inFile;      // stream for input file
 
// intialisation 
  string dataPath="../";    // location of data directory (relative path)
  fPath=dataPath+fileName;
  
// check if file exist
  inFile.open(fPath.c_str());
  if(!inFile)
  {
    cout << endl;
    cout << "file : " << fileName << " not found!" << endl; 
    cout << "quit ntuple raw subroutine" << endl;
    cout << endl;
  }
  else
  {  
// define ntuple name, reaction, author and comment
    // NTuple name
    pos1    = fileName.find("_");
    ntName  = fileName.substr(pos1+1);      // ntuple name from file name
    // author
    pos1    = ntName.find("_");
    author  = ntName.substr(pos1+1);        // author name + comment
    // reaction
    reaction  = ntName.substr(0, pos1);
    ntName    = ntName.substr(5, ntName.size()-9); // clear "pbar" ...
    author    = author.substr(0, author.size()-4); // ... and ".dat"
    pos1      = reaction.find("-");
    reaction.replace(pos1, 1, " + ");
    pos1      = reaction.find("--");
    reaction.replace(pos1, 2, " --> ");
    // comment on reaction from file name
    pos1      = author.find("_");
    comment   = "";
    if(pos1 != string::npos)  // extract comment if there is comment!
    {
      comment = author.substr(pos1+1);
      author  = author.substr(0, pos1);
      ntTitle = reaction + " (" + author + " - " + comment + ")";
    }
    else
    {
      ntTitle = reaction + " (" + author + ")";
    }
// dump properties
    cout << "__________________________________________________________________" 
    << endl;
    cout << "file     : " << fileName << endl;
    cout << "reaction : " << reaction << endl;
    cout << "author   : " << author << endl;
    cout << "comment  : " << comment << endl;
    cout << "ntName   : " << ntName << endl;
    cout << "ntTitle  : " << ntTitle << endl;
    cout << "__________________________________________________________________" 
    << endl;
    
// get Ntuple, close file stream, write Ntuple
    string varlist("energy:cross_section");
    TNtuple ntuple(ntName.c_str(), ntTitle.c_str(), varlist.c_str());
    ntuple.ReadStream(inFile);
    ntuple.Write();
    TGraphErrors tge(fPath.c_str(), "%lg %lg");
    tge.SetTitle(ntTitle.c_str());
    tge.SetName(ntName.c_str());
    tge.Write();
  }
// last statements
  inFile.close();
}

// ntuple Process3 *************************************************************
void ntupleProcess3(string const& fileName, string const& fileNameMain)
{
// declarations for ratio file ...
  string fPath;         // entire pathname of file
  string ntName;        // name of ntuple to be crated
  string ntTitle;       // title of ntuple to be created
  string author;        // author of data
  string reaction;      // reaction involved
  string comment;       // extra comment
  std::size_t pos1(0);  // position into a string
  ifstream inFile;      // stream for input file

// ... and main ionisation process file
  string fPathMain;
  ifstream inFileMain;  // stream for main reaction input file
  
// intialisation 
  string dataPath="../";    // location of data directory (relative path)
  fPath=dataPath+fileName;
  fPathMain=dataPath+fileNameMain;

// check if files exist
  inFile.open(fPath.c_str());
  inFileMain.open(fPathMain.c_str());
  if(!inFile)
  {
    cout << endl;
    cout << "ratio file : " << fileName << " not found!" << endl; 
    cout << "quit \"ntuple process\" subroutine" << endl;
    cout << endl;
  }
  else if (!inFileMain)
  {
    cout << endl;
    cout << "main reaction file : " << fileName << " not found!" << endl; 
    cout << "quit \"ntuple process\" subroutine" << endl;
    cout << endl;
  }
  else
  {
// define ntuple name, reaction, author and comment
    // NTuple name
    pos1    = fileName.find("_");
    ntName  = fileName.substr(pos1+1);      // ntuple name from file name
    // author
    pos1    = ntName.find("_");
    author  = ntName.substr(pos1+1);        // author name + comment
    // reaction
    reaction  = ntName.substr(0, pos1);
    ntName    = ntName.substr(5, ntName.size()-9); // clear "pbar" ...
    author    = author.substr(0, author.size()-4); // ... and ".dat"
    pos1      = reaction.find("-");
    reaction.replace(pos1, 1, " + ");
    pos1      = reaction.find("--");
    reaction.replace(pos1, 2, " --> ");
    // comment on reaction from file name
    pos1      = author.find("_");
    comment   = "";
    if(pos1 != string::npos)  // extract comment if there is comment!
    {
      comment = author.substr(pos1+1);
      author  = author.substr(0, pos1);
      ntTitle = reaction + " (" + author + " - " + comment + ")";
    }
    else { ntTitle = reaction + " (" + author + ")"; }

 // dump properties
    cout << "__________________________________________________________________" 
    << endl;
    cout << "file     : " << fileName << endl;
    cout << "reaction : " << reaction << endl;
    cout << "author   : " << author << endl;
    cout << "comment  : " << comment << endl;
    cout << "ntName   : " << ntName << endl;
    cout << "ntTitle  : " << ntTitle << endl;
 
// _____________________________________________________________________________
// Get main ionisation process cs to compute sec. cs from ratio
    double energyTemp;                  // [keV] temp energy (while loop)
    double cross_sectionTemp;           // [e-20 m^2] temp cs (while loop)
    double sigmaTemp;                   // [e-20 m^2] cs sigma (while loop)
    vector<double> energyMain;          // [keV] main ionisation channel energy
    vector<double> cross_sectionMain;   // [e-20 m^2] main ionisation channel cs
    vector<double> sigmaMain;           // [e-20 m^2] main ion. channel sigma
    vector<double> energyRat;           // [keV] energy from actual channel RATIO file
    vector<double> cross_sectionRat;    // [ ] actual ionisation channel cs RATIO
    vector<double> sigmaRat;            // [ ] actual ion. channel sigma RATIO
    int dataSizeMain(0);                // size of main channel vector
    int dataSizeRat(0);                 // size of actual channel vector
  
    string ignoreHeader="#";            // to dismiss comment line in file
    string inputString;                 // string for getting line in file
  
// while loop over all lines of main ionisation process
    //cout << endl;
    //cout << "MAIN CHANNEL FILE SCAN" << endl;
    while(getline(inFileMain, inputString))
    {
      if(inputString[0]==ignoreHeader)
      { 
        //cout << inputString << endl;
        continue;
      }
      else
      {
        stringstream ss(inputString, stringstream::in);
        ss >> energyTemp;
        ss >> cross_sectionTemp;
        ss >> sigmaTemp;
        //cout << energyTemp << "   " << cross_sectionTemp 
        //<< "   " << sigmaTemp << endl;
        energyMain.push_back(energyTemp);
        cross_sectionMain.push_back(cross_sectionTemp);
        sigmaMain.push_back(sigmaTemp);
      }
    }
    dataSizeMain=energyMain.size(); // get size of vector
    inFileMain.close();             // close main ionisation channel file

// control that data are aquired
    for (int i=0; i<dataSizeMain; i++)
    {
      //cout << setw(6) << energyMain[i] << "   " << setw(6) << cross_sectionMain[i]
      //<< "   " << setw(6) << sigmaMain[i] << endl;
    }
    //cout << "NUMBER OF LINES : " << dataSizeMain << endl;
    //cout << endl;
    
// while loop over all lines of actual ionisation process
    //cout << endl;
    //cout << "RATIO CHANNEL FILE SCAN" << endl;
    while(getline(inFile, inputString))
    {
      if(inputString[0]==ignoreHeader)
      { 
        //cout << inputString << endl;
        continue;
      }
      else
      {
        stringstream ss(inputString, stringstream::in);
        ss >> energyTemp;
        ss >> cross_sectionTemp;
        ss >> sigmaTemp;
        //cout << energyTemp << "   " << cross_sectionTemp 
        //<< "   " << sigmaTemp << endl;
        energyRat.push_back(energyTemp);
        cross_sectionRat.push_back(cross_sectionTemp);
        sigmaRat.push_back(sigmaTemp);
      }
    }
    dataSizeRat=energyRat.size();   // get size of vector
    inFile.close();                 // close actual ionisation channel file

// control that data are aquired
    for (int i=0; i<dataSizeRat; i++)
    {
      //cout << setw(6) << energyRat[i] << "   " << setw(6) << cross_sectionRat[i]
      //<< "   " << setw(6) << sigmaRat[i] << endl;
    }
    //cout << "NUMBER OF LINES : " << dataSizeRat << endl;
    //cout << endl;
    
// convert actual channel ratio values from file into cs
    //cout << endl;
    //cout << "ACTUAL CHANNEL CONVERTION" << endl;
    if (dataSizeRat > dataSizeMain)
    {
      cout << "!!! BIG PROBLEM !!!" 
      << "second channel have more entries than main channel" << endl;
      cout << "No conversion will be done, please check data and code" << endl;
    }
    else
    {
      if (dataSizeRat < dataSizeMain)
      {
        cout << "!!! WARNING !!!" 
        << "second channel have less entries than main channel" << endl;
      }
      int const sizeSec(dataSizeRat);
      double energySec[sizeSec];         // [keV] actual ion. channel energy
      double  cross_sectionSec[sizeSec]; // [e-20 m^2] actual ion. channel cs
      double  sigmaSec[sizeSec];         // [e-20 m^2] actual ion. channel sigma
      
      // set Ntuple and Graph, convert data, fill, write...
      string varlist("energy:cross_section:sigma");
      TNtuple ntuple(ntName.c_str(), ntTitle.c_str(), varlist.c_str());
      TGraphErrors tge(sizeSec);
      tge.SetTitle(ntTitle.c_str());
      tge.SetName(ntName.c_str());
      for (int i=0; i<sizeSec; i++)     
      {
        energySec[i]=energyRat[i];
        cross_sectionSec[i]=cross_sectionRat[i]/1000. * cross_sectionMain[i];
        sigmaSec[i]=sigmaRat[i]/1000. * cross_sectionMain[i];
        
        //cout << setw(6) << energySec[i] << "   " << setw(8) << cross_sectionSec[i]
        //<< "   " << setw(8) << sigmaSec[i] << endl;
        
        ntuple.Fill(energySec[i], cross_sectionSec[i], sigmaSec[i]);
        tge.SetPoint(i, energySec[i], cross_sectionSec[i]);
        tge.SetPointError(i, 0, sigmaSec[i]);
      }
      ntuple.Write();
      tge.Write();
    }
// finish dumping
    cout << "__________________________________________________________________" 
    << endl;
  }
// last statements
}

// _____________________________________________________________________________
// _____________________________________________________________________________
// _____________________________________________________________________________
// _____________________________________________________________________________
// _____________________________________________________________________________
// _____________________________________________________________________________
// _____________________________________________________________________________
// _____________________________________________________________________________
// _____________________________________________________________________________
// _____________________________________________________________________________
// _____________________________________________________________________________
// _____________________________________________________________________________
// setCSsum()
void setCSsum(string const& gasName)
{
  // Open TFile with cross-sections  
  string const tfName = "tfPbarCrossSections.root";           // rootfile name
  string const tfDir = "";                                    // directory path
  string const tfPath = tfDir + tfName;                       // rootfile path
  
  TFile *tf = new TFile(tfPath.c_str(), "UPDATE");      // open root file
  if (tf->IsZombie())                                   // test if file is open
  { cout << "Error opening : " << tfPath << endl; }
  else 
  {
// *****************************************************************************
// H2    
    if (gasName=="H2")
    {
// retrieve graphs
      int const nGraphs = 4;  // number of graphs
      TGraphErrors *g1=(TGraphErrors*)tf->Get("H2--H2plus_Hvelplund");
      TGraphErrors *g2=(TGraphErrors*)tf->Get("H2--Hplus_Hvelplund");
      TGraphErrors *g3=(TGraphErrors*)tf->Get("H2--H2plus_Andersen");

// different process than for next gases
// here we sum g1 and g2, g3 comes from an other data set
// and g4 is a calculation   
// We create also a H2--H2plus_mix to have  cross section along 
// energy range of g1 and g3      
 
// gMix, g1 and g3
      int g1N = g1->GetN();   // number of points in g1
      double g1Max = TMath::MaxElement(g1N, g1->GetX());  // max element in g1
      int g3Add(0);           // g3 number of point out of range of g1
      int g3Ind(0);          // Indice of min element of g3 to add to g1
                              
      // get g3Add loop
      for (int i=0; i < g3->GetN(); ++i)
      {
        if ( (g3->GetX())[i] > g1Max )
        {
          ++g3Add;
        }
        else {g3Ind = i; }
      }
      ++g3Ind;
      
      TGraphErrors *gMix = new TGraphErrors(g1N+g3Add);
      gMix->SetName((gasName + "+" + "--Mix").c_str());
      gMix->SetTitle(("pbar/" + gasName + "+" + " cross section mix").c_str());
      double *xMix = gMix->GetX();
      double *yMix = gMix->GetY();
      double *exMix = gMix->GetEX();
      double *eyMix = gMix->GetEY();
      
      cout << endl;
      cout << "gMix" << endl;
      for (int i=0; i<g1N; ++i)
      {
        xMix[i] = (g1->GetX())[i];
        yMix[i] = (g1->GetY())[i];
        exMix[i] = (g1->GetEX())[i];
        eyMix[i] = (g1->GetEY())[i];
        
        cout << setw(10) << xMix[i] 
        << setw(10) << yMix[i] 
        << setw(10) << exMix[i] 
        << setw(10) << eyMix[i] << endl;
      }
      for (int i=0; i<g3Add; ++i)
      {
        xMix[g1N+i] = (g3->GetX())[g3Ind+i];
        yMix[g1N+i] = (g3->GetY())[g3Ind+i];
        exMix[g1N+i] = (g3->GetEX())[g3Ind+i];
        eyMix[g1N+i] = (g3->GetEY())[g3Ind+i];
        
        cout << setw(10) << xMix[g1N+i] 
        << setw(10) << yMix[g1N+i] 
        << setw(10) << exMix[g1N+i] 
        << setw(10) << eyMix[g1N+i] << endl;
      }
      gMix->Write();
      
// gSum, g1 + g2 and g3
      TGraphErrors *gSum = new TGraphErrors(g1N+g3Add);
      gSum->SetName((gasName + "--All").c_str());
      gSum->SetTitle(("pbar/" + gasName + " total cross section").c_str());
      double *xSum = gSum->GetX();
      double *ySum = gSum->GetY();
      double *exSum = gSum->GetEX();
      double *eySum = gSum->GetEY();
      
      cout << endl;
      cout << "g    " << setw(6) << "gx[i]" 
      << setw(10) << "gy[i]" << setw(10) << "gey[i]"
      << setw(10) << "xSum[i]" << setw(10) << "ySum[i]" << setw(10) 
      << "eySum[i]" << endl;

      for (int i=0; i<g1->GetN(); ++i)
      {
        cout << endl;
        
        xSum[i] = (g1->GetX())[i];
        ySum[i] = (g1->GetY())[i];
        eySum[i] = (g1->GetEY())[i];
        cout << "g1   " << setw(6) << (g1->GetX())[i] 
        << setw(10) << (g1->GetY())[i] << setw(10) << (g1->GetEY())[i]
        << setw(10) << xSum[i] << setw(10) << ySum[i] << setw(10) << eySum[i] << endl;
        
        if (i < g2->GetN())
        { 
          ySum[i] +=  (g2->GetY())[i]; 
          eySum[i] =  sqrt(pow(eySum[i], 2) + pow((g2->GetEY())[i], 2));
        }
        else { cout << "no sum" << endl; }
        cout << "g2   " << setw(6) << (g2->GetX())[i] 
        << setw(10) << (g2->GetY())[i] << setw(10) << (g2->GetEY())[i]
        << setw(10) << xSum[i] << setw(10) << ySum[i] << setw(10) << eySum[i] << endl;
      }
      for (int i=0; i<g3Add; ++i)
      {
        cout << endl;
        
        xSum[g1N+i] = (g3->GetX())[g3Ind+i];
        ySum[g1N+i] = (g3->GetY())[g3Ind+i];
        exSum[g1N+i] = (g3->GetEX())[g3Ind+i];
        eySum[g1N+i] = (g3->GetEY())[g3Ind+i];
        
        cout << "g3Add" << setw(6) << (g3->GetX())[g3Ind+i] 
        << setw(10) << (g3->GetY())[g3Ind+i] << setw(10) 
        << (g3->GetEY())[g3Ind+i] << setw(10) << xSum[g1N+i] << setw(10) 
        << ySum[g1N+i] << setw(10) << eySum[g1N+i] << endl;
      }
      gSum->Write();

    }
// *****************************************************************************
// *****************************************************************************
// CO    
    if (gasName=="CO")
    {
// retrieve graphs
      int const nGraphs = 6;  // number of graphs
      TGraphErrors *g1=(TGraphErrors*)tf->Get("CO--COplus_Knudsen");
      TGraphErrors *g2=(TGraphErrors*)tf->Get("CO--Oplus_Knudsen");
      TGraphErrors *g3=(TGraphErrors*)tf->Get("CO--COplusplus_Knudsen");
      TGraphErrors *g4=(TGraphErrors*)tf->Get("CO--Cplus_Knudsen");
      TGraphErrors *g5=(TGraphErrors*)tf->Get("CO--Oplusplus_Knudsen");
      TGraphErrors *g6=(TGraphErrors*)tf->Get("CO--Cplusplus_Knudsen");
      
      TGraphErrors *gSum = new TGraphErrors(g1->GetN());
      gSum->SetName((gasName + "--All").c_str());
      gSum->SetTitle(("pbar/" + gasName + " total cross section").c_str());
      double *xSum = gSum->GetX();
      double *ySum = gSum->GetY();
      double *exSum = gSum->GetEX();
      double *eySum = gSum->GetEY();
      
      cout << endl;
      cout << "g    " << setw(6) << "gx[i]" 
      << setw(10) << "gy[i]" << setw(10) << "gey[i]"
      << setw(10) << "xSum[i]" << setw(10) << "ySum[i]" << setw(10) 
      << "eySum[i]" << endl;
      for (int i=0; i<gSum->GetN(); ++i)
      {
        cout << endl;
        
        xSum[i] = (g1->GetX())[i];
        ySum[i] = (g1->GetY())[i];
        eySum[i] = (g1->GetEY())[i];
        cout << "g1   " << setw(6) << (g1->GetX())[i] 
        << setw(10) << (g1->GetY())[i] << setw(10) << (g1->GetEY())[i]
        << setw(10) << xSum[i] << setw(10) << ySum[i] << setw(10) << eySum[i] << endl;
        
        if (i < g2->GetN())
        { 
          ySum[i] +=  (g2->GetY())[i]; 
          eySum[i] =  sqrt(pow(eySum[i], 2) + pow((g2->GetEY())[i], 2));
        }
        else { cout << "no sum" << endl; }
        cout << "g2   " << setw(6) << (g2->GetX())[i] 
        << setw(10) << (g2->GetY())[i] << setw(10) << (g2->GetEY())[i]
        << setw(10) << xSum[i] << setw(10) << ySum[i] << setw(10) << eySum[i] << endl;
        
        if (i < g3->GetN())
        { 
          ySum[i] +=  (g3->GetY())[i]; 
          eySum[i] =  sqrt(pow(eySum[i], 2) + pow((g3->GetEY())[i], 2));
        }
        else { cout << "no sum" << endl; }
        cout << "g3   " << setw(6) << (g3->GetX())[i] 
        << setw(10) << (g3->GetY())[i] << setw(10) << (g3->GetEY())[i]
        << setw(10) << xSum[i] << setw(10) << ySum[i] << setw(10) << eySum[i] << endl;
        
        if (i < g4->GetN())
        { 
          ySum[i] +=  (g4->GetY())[i]; 
          eySum[i] =  sqrt(pow(eySum[i], 2) + pow((g4->GetEY())[i], 2));
        }
        else { cout << "no sum" << endl; }
        cout << "g4   " << setw(6) << (g4->GetX())[i] 
        << setw(10) << (g4->GetY())[i] << setw(10) << (g4->GetEY())[i]
        << setw(10) << xSum[i] << setw(10) << ySum[i] << setw(10) << eySum[i] << endl;
        
        if (i < g5->GetN())
        { 
          ySum[i] +=  (g5->GetY())[i]; 
          eySum[i] =  sqrt(pow(eySum[i], 2) + pow((g5->GetEY())[i], 2));
        }
        else { cout << "no sum" << endl; }
        cout << "g5   " << setw(6) << (g5->GetX())[i] 
        << setw(10) << (g5->GetY())[i] << setw(10) << (g5->GetEY())[i]
        << setw(10) << xSum[i] << setw(10) << ySum[i] << setw(10) << eySum[i] << endl;
        
        if (i < g6->GetN())
        { 
          ySum[i] +=  (g6->GetY())[i]; 
          eySum[i] =  sqrt(pow(eySum[i], 2) + pow((g6->GetEY())[i], 2));
        }
        else { cout << "no sum" << endl; }
        cout << "g6   " << setw(6) << (g6->GetX())[i] 
        << setw(10) << (g6->GetY())[i] << setw(10) << (g6->GetEY())[i]
        << setw(10) << xSum[i] << setw(10) << ySum[i] << setw(10) << eySum[i] << endl;
      }
      gSum->Write();
    }
// *****************************************************************************
// *****************************************************************************
// CO2    
    if (gasName=="CO2")
    {
// retrieve graphs
      int const nGraphs = 7;  // number of graphs
      TGraphErrors *g1=(TGraphErrors*)tf->Get("CO2--CO2plus_Knudsen");
      TGraphErrors *g2=(TGraphErrors*)tf->Get("CO2--COplus_Knudsen");
      TGraphErrors *g3=(TGraphErrors*)tf->Get("CO2--CO2plusplus_Knudsen");
      TGraphErrors *g4=(TGraphErrors*)tf->Get("CO2--Oplus_Knudsen");
      TGraphErrors *g5=(TGraphErrors*)tf->Get("CO2--Cplus_Knudsen");
      TGraphErrors *g6=(TGraphErrors*)tf->Get("CO2--Oplusplus_Knudsen");
      TGraphErrors *g7=(TGraphErrors*)tf->Get("CO2--Cplusplus_Knudsen");
      
      TGraphErrors *gSum = new TGraphErrors(g1->GetN());
      gSum->SetName((gasName + "--All").c_str());
      gSum->SetTitle(("pbar/" + gasName + " total cross section").c_str());
      double *xSum = gSum->GetX();
      double *ySum = gSum->GetY();
      double *exSum = gSum->GetEX();
      double *eySum = gSum->GetEY();
      
      cout << endl;
      cout << "g    " << setw(6) << "gx[i]" 
      << setw(10) << "gy[i]" << setw(10) << "gey[i]"
      << setw(10) << "xSum[i]" << setw(10) << "ySum[i]" << setw(10) 
      << "eySum[i]" << endl;
      for (int i=0; i<gSum->GetN(); ++i)
      {
        cout << endl;
        
        xSum[i] = (g1->GetX())[i];
        ySum[i] = (g1->GetY())[i];
        eySum[i] = (g1->GetEY())[i];
        cout << "g1   " << setw(6) << (g1->GetX())[i] 
        << setw(10) << (g1->GetY())[i] << setw(10) << (g1->GetEY())[i]
        << setw(10) << xSum[i] << setw(10) << ySum[i] << setw(10) << eySum[i] << endl;
        
        if (i < g2->GetN())
        { 
          ySum[i] +=  (g2->GetY())[i]; 
          eySum[i] =  sqrt(pow(eySum[i], 2) + pow((g2->GetEY())[i], 2));
        }
        else { cout << "no sum" << endl; }
        cout << "g2   " << setw(6) << (g2->GetX())[i] 
        << setw(10) << (g2->GetY())[i] << setw(10) << (g2->GetEY())[i]
        << setw(10) << xSum[i] << setw(10) << ySum[i] << setw(10) << eySum[i] << endl;
        
        if (i < g3->GetN())
        { 
          ySum[i] +=  (g3->GetY())[i]; 
          eySum[i] =  sqrt(pow(eySum[i], 2) + pow((g3->GetEY())[i], 2));
        }
        else { cout << "no sum" << endl; }
        cout << "g3   " << setw(6) << (g3->GetX())[i] 
        << setw(10) << (g3->GetY())[i] << setw(10) << (g3->GetEY())[i]
        << setw(10) << xSum[i] << setw(10) << ySum[i] << setw(10) << eySum[i] << endl;
        
        if (i < g4->GetN())
        { 
          ySum[i] +=  (g4->GetY())[i]; 
          eySum[i] =  sqrt(pow(eySum[i], 2) + pow((g4->GetEY())[i], 2));
        }
        else { cout << "no sum" << endl; }
        cout << "g4   " << setw(6) << (g4->GetX())[i] 
        << setw(10) << (g4->GetY())[i] << setw(10) << (g4->GetEY())[i]
        << setw(10) << xSum[i] << setw(10) << ySum[i] << setw(10) << eySum[i] << endl;
        
        if (i < g5->GetN())
        { 
          ySum[i] +=  (g5->GetY())[i]; 
          eySum[i] =  sqrt(pow(eySum[i], 2) + pow((g5->GetEY())[i], 2));
        }
        else { cout << "no sum" << endl; }
        cout << "g5   " << setw(6) << (g5->GetX())[i] 
        << setw(10) << (g5->GetY())[i] << setw(10) << (g5->GetEY())[i]
        << setw(10) << xSum[i] << setw(10) << ySum[i] << setw(10) << eySum[i] << endl;
        
        if (i < g6->GetN())
        { 
          ySum[i] +=  (g6->GetY())[i]; 
          eySum[i] =  sqrt(pow(eySum[i], 2) + pow((g6->GetEY())[i], 2));
        }
        else { cout << "no sum" << endl; }
        cout << "g6   " << setw(6) << (g6->GetX())[i] 
        << setw(10) << (g6->GetY())[i] << setw(10) << (g6->GetEY())[i]
        << setw(10) << xSum[i] << setw(10) << ySum[i] << setw(10) << eySum[i] << endl;
        
        if (i < g7->GetN())
        { 
          ySum[i] +=  (g7->GetY())[i]; 
          eySum[i] =  sqrt(pow(eySum[i], 2) + pow((g7->GetEY())[i], 2));
        }
        else { cout << "no sum" << endl; }
        cout << "g7   " << setw(6) << (g7->GetX())[i] 
        << setw(10) << (g7->GetY())[i] << setw(10) << (g7->GetEY())[i]
        << setw(10) << xSum[i] << setw(10) << ySum[i] << setw(10) << eySum[i] << endl;
      }
      gSum->Write();
    }
// *****************************************************************************
// *****************************************************************************
// CH4    
    if (gasName=="CH4")
    {
// retrieve graphs
      int const nGraphs = 7;  // number of graphs
      TGraphErrors *g1=(TGraphErrors*)tf->Get("CH4--CH4plus_Knudsen");
      TGraphErrors *g2=(TGraphErrors*)tf->Get("CH4--CH3plus_Knudsen");
      TGraphErrors *g3=(TGraphErrors*)tf->Get("CH4--CH2plus_Knudsen");
      TGraphErrors *g4=(TGraphErrors*)tf->Get("CH4--CHplus_Knudsen");
      TGraphErrors *g5=(TGraphErrors*)tf->Get("CH4--Cplus_Knudsen");
      TGraphErrors *g6=(TGraphErrors*)tf->Get("CH4--H2plus_Knudsen");
      TGraphErrors *g7=(TGraphErrors*)tf->Get("CH4--Hplus_Knudsen");
      
      TGraphErrors *gSum = new TGraphErrors(g1->GetN());
      gSum->SetName((gasName + "--All").c_str());
      gSum->SetTitle(("pbar/" + gasName + " total cross section").c_str());
      double *xSum = gSum->GetX();
      double *ySum = gSum->GetY();
      double *exSum = gSum->GetEX();
      double *eySum = gSum->GetEY();
      
      cout << endl;
      cout << "g    " << setw(6) << "gx[i]" 
      << setw(10) << "gy[i]" << setw(10) << "gey[i]"
      << setw(10) << "xSum[i]" << setw(10) << "ySum[i]" << setw(10) 
      << "eySum[i]" << endl;
      for (int i=0; i<gSum->GetN(); ++i)
      {
        cout << endl;
        
        xSum[i] = (g1->GetX())[i];
        ySum[i] = (g1->GetY())[i];
        eySum[i] = (g1->GetEY())[i];
        cout << "g1   " << setw(6) << (g1->GetX())[i] 
        << setw(10) << (g1->GetY())[i] << setw(10) << (g1->GetEY())[i]
        << setw(10) << xSum[i] << setw(10) << ySum[i] << setw(10) << eySum[i] << endl;
        
        if (i < g2->GetN())
        { 
          ySum[i] +=  (g2->GetY())[i]; 
          eySum[i] =  sqrt(pow(eySum[i], 2) + pow((g2->GetEY())[i], 2));
        }
        else { cout << "no sum" << endl; }
        cout << "g2   " << setw(6) << (g2->GetX())[i] 
        << setw(10) << (g2->GetY())[i] << setw(10) << (g2->GetEY())[i]
        << setw(10) << xSum[i] << setw(10) << ySum[i] << setw(10) << eySum[i] << endl;
        
        if (i < g3->GetN())
        { 
          ySum[i] +=  (g3->GetY())[i]; 
          eySum[i] =  sqrt(pow(eySum[i], 2) + pow((g3->GetEY())[i], 2));
        }
        else { cout << "no sum" << endl; }
        cout << "g3   " << setw(6) << (g3->GetX())[i] 
        << setw(10) << (g3->GetY())[i] << setw(10) << (g3->GetEY())[i]
        << setw(10) << xSum[i] << setw(10) << ySum[i] << setw(10) << eySum[i] << endl;
        
        if (i < g4->GetN())
        { 
          ySum[i] +=  (g4->GetY())[i]; 
          eySum[i] =  sqrt(pow(eySum[i], 2) + pow((g4->GetEY())[i], 2));
        }
        else { cout << "no sum" << endl; }
        cout << "g4   " << setw(6) << (g4->GetX())[i] 
        << setw(10) << (g4->GetY())[i] << setw(10) << (g4->GetEY())[i]
        << setw(10) << xSum[i] << setw(10) << ySum[i] << setw(10) << eySum[i] << endl;
        
        if (i < g5->GetN())
        { 
          ySum[i] +=  (g5->GetY())[i]; 
          eySum[i] =  sqrt(pow(eySum[i], 2) + pow((g5->GetEY())[i], 2));
        }
        else { cout << "no sum" << endl; }
        cout << "g5   " << setw(6) << (g5->GetX())[i] 
        << setw(10) << (g5->GetY())[i] << setw(10) << (g5->GetEY())[i]
        << setw(10) << xSum[i] << setw(10) << ySum[i] << setw(10) << eySum[i] << endl;
        
        if (i < g6->GetN())
        { 
          ySum[i] +=  (g6->GetY())[i]; 
          eySum[i] =  sqrt(pow(eySum[i], 2) + pow((g6->GetEY())[i], 2));
        }
        else { cout << "no sum" << endl; }
        cout << "g6   " << setw(6) << (g6->GetX())[i] 
        << setw(10) << (g6->GetY())[i] << setw(10) << (g6->GetEY())[i]
        << setw(10) << xSum[i] << setw(10) << ySum[i] << setw(10) << eySum[i] << endl;
        
        if (i < g7->GetN())
        { 
          ySum[i] +=  (g7->GetY())[i]; 
          eySum[i] =  sqrt(pow(eySum[i], 2) + pow((g7->GetEY())[i], 2));
        }
        else { cout << "no sum" << endl; }
        cout << "g7   " << setw(6) << (g7->GetX())[i] 
        << setw(10) << (g7->GetY())[i] << setw(10) << (g7->GetEY())[i]
        << setw(10) << xSum[i] << setw(10) << ySum[i] << setw(10) << eySum[i] << endl;
      }
      gSum->Write();
    }
// *****************************************************************************
// *****************************************************************************
// N2    
    if (gasName=="N2")
    {
// retrieve graphs
      int const nGraphs = 3;  // number of graphs
      TGraphErrors *g1=(TGraphErrors*)tf->Get("N2--N2plus_Knudsen");
      TGraphErrors *g2=(TGraphErrors*)tf->Get("N2--Nplus_Knudsen");
      TGraphErrors *g3=(TGraphErrors*)tf->Get("N2--Nplusplus_Knudsen");
      
      TGraphErrors *gSum = new TGraphErrors(g1->GetN());
      gSum->SetName((gasName + "--All").c_str());
      gSum->SetTitle(("pbar/" + gasName + " total cross section").c_str());
      double *xSum = gSum->GetX();
      double *ySum = gSum->GetY();
      double *exSum = gSum->GetEX();
      double *eySum = gSum->GetEY();
      
      cout << endl;
      cout << "g    " << setw(6) << "gx[i]" 
      << setw(10) << "gy[i]" << setw(10) << "gey[i]"
      << setw(10) << "xSum[i]" << setw(10) << "ySum[i]" << setw(10) 
      << "eySum[i]" << endl;
      for (int i=0; i<gSum->GetN(); ++i)
      {
        cout << endl;
        
        xSum[i] = (g1->GetX())[i];
        ySum[i] = (g1->GetY())[i];
        eySum[i] = (g1->GetEY())[i];
        cout << "g1   " << setw(6) << (g1->GetX())[i] 
        << setw(10) << (g1->GetY())[i] << setw(10) << (g1->GetEY())[i]
        << setw(10) << xSum[i] << setw(10) << ySum[i] << setw(10) << eySum[i] << endl;
        
        if (i < g2->GetN())
        { 
          ySum[i] +=  (g2->GetY())[i]; 
          eySum[i] =  sqrt(pow(eySum[i], 2) + pow((g2->GetEY())[i], 2));
        }
        else { cout << "no sum" << endl; }
        cout << "g2   " << setw(6) << (g2->GetX())[i] 
        << setw(10) << (g2->GetY())[i] << setw(10) << (g2->GetEY())[i]
        << setw(10) << xSum[i] << setw(10) << ySum[i] << setw(10) << eySum[i] << endl;
        
        if (i < g3->GetN())
        { 
          ySum[i] +=  (g3->GetY())[i]; 
          eySum[i] =  sqrt(pow(eySum[i], 2) + pow((g3->GetEY())[i], 2));
        }
        else { cout << "no sum" << endl; }
        cout << "g3   " << setw(6) << (g3->GetX())[i] 
        << setw(10) << (g3->GetY())[i] << setw(10) << (g3->GetEY())[i]
        << setw(10) << xSum[i] << setw(10) << ySum[i] << setw(10) << eySum[i] << endl;
      }
      gSum->Write();
    }
// *****************************************************************************
// *****************************************************************************
// Ne    
    if (gasName=="Ne")
    {
// retrieve graphs
      int const nGraphs = 3;  // number of graphs
      TGraphErrors *g1=(TGraphErrors*)tf->Get("Ne--Neplus_Knudsen");
      TGraphErrors *g2=(TGraphErrors*)tf->Get("Ne--Neplusplus_Knudsen");
      TGraphErrors *g3=(TGraphErrors*)tf->Get("Ne--Neplusplusplus_Knudsen");
      
      TGraphErrors *gSum = new TGraphErrors(g1->GetN());
      gSum->SetName((gasName + "--All").c_str());
      gSum->SetTitle(("pbar/" + gasName + " total cross section").c_str());
      double *xSum = gSum->GetX();
      double *ySum = gSum->GetY();
      double *exSum = gSum->GetEX();
      double *eySum = gSum->GetEY();
      
      cout << endl;
      cout << "g    " << setw(6) << "gx[i]" 
      << setw(10) << "gy[i]" << setw(10) << "gey[i]"
      << setw(10) << "xSum[i]" << setw(10) << "ySum[i]" << setw(10) 
      << "eySum[i]" << endl;
      for (int i=0; i<gSum->GetN(); ++i)
      {
        cout << endl;
        
        xSum[i] = (g1->GetX())[i];
        ySum[i] = (g1->GetY())[i];
        eySum[i] = (g1->GetEY())[i];
        cout << "g1   " << setw(6) << (g1->GetX())[i] 
        << setw(10) << (g1->GetY())[i] << setw(10) << (g1->GetEY())[i]
        << setw(10) << xSum[i] << setw(10) << ySum[i] << setw(10) << eySum[i] << endl;
        
        if (i < g2->GetN())
        { 
          ySum[i] +=  (g2->GetY())[i]; 
          eySum[i] =  sqrt(pow(eySum[i], 2) + pow((g2->GetEY())[i], 2));
        }
        else { cout << "no sum" << endl; }
        cout << "g2   " << setw(6) << (g2->GetX())[i] 
        << setw(10) << (g2->GetY())[i] << setw(10) << (g2->GetEY())[i]
        << setw(10) << xSum[i] << setw(10) << ySum[i] << setw(10) << eySum[i] << endl;
        
        if (i < g3->GetN())
        { 
          ySum[i] +=  (g3->GetY())[i]; 
          eySum[i] =  sqrt(pow(eySum[i], 2) + pow((g3->GetEY())[i], 2));
        }
        else { cout << "no sum" << endl; }
        cout << "g3   " << setw(6) << (g3->GetX())[i] 
        << setw(10) << (g3->GetY())[i] << setw(10) << (g3->GetEY())[i]
        << setw(10) << xSum[i] << setw(10) << ySum[i] << setw(10) << eySum[i] << endl;
      }
      gSum->Write();
    }
// *****************************************************************************
// *****************************************************************************
// Ar    
    if (gasName=="Ar")
    {
// retrieve graphs
      int const nGraphs = 3;  // number of graphs
      TGraphErrors *g1=(TGraphErrors*)tf->Get("Ar--Arplus_Knudsen");
      TGraphErrors *g2=(TGraphErrors*)tf->Get("Ar--Arplusplus_Knudsen");
      TGraphErrors *g3=(TGraphErrors*)tf->Get("Ar--Arplusplusplus_Knudsen");
      
      TGraphErrors *gSum = new TGraphErrors(g1->GetN());
      gSum->SetName((gasName + "--All").c_str());
      gSum->SetTitle(("pbar/" + gasName + " total cross section").c_str());
      double *xSum = gSum->GetX();
      double *ySum = gSum->GetY();
      double *exSum = gSum->GetEX();
      double *eySum = gSum->GetEY();
      
      cout << endl;
      cout << "g    " << setw(6) << "gx[i]" 
      << setw(10) << "gy[i]" << setw(10) << "gey[i]"
      << setw(10) << "xSum[i]" << setw(10) << "ySum[i]" << setw(10) 
      << "eySum[i]" << endl;
      for (int i=0; i<gSum->GetN(); ++i)
      {
        cout << endl;
        
        xSum[i] = (g1->GetX())[i];
        ySum[i] = (g1->GetY())[i];
        eySum[i] = (g1->GetEY())[i];
        cout << "g1   " << setw(6) << (g1->GetX())[i] 
        << setw(10) << (g1->GetY())[i] << setw(10) << (g1->GetEY())[i]
        << setw(10) << xSum[i] << setw(10) << ySum[i] << setw(10) << eySum[i] << endl;
        
        if (i < g2->GetN())
        { 
          ySum[i] +=  (g2->GetY())[i]; 
          eySum[i] =  sqrt(pow(eySum[i], 2) + pow((g2->GetEY())[i], 2));
        }
        else { cout << "no sum" << endl; }
        cout << "g2   " << setw(6) << (g2->GetX())[i] 
        << setw(10) << (g2->GetY())[i] << setw(10) << (g2->GetEY())[i]
        << setw(10) << xSum[i] << setw(10) << ySum[i] << setw(10) << eySum[i] << endl;
        
        if (i < g3->GetN())
        { 
          ySum[i] +=  (g3->GetY())[i]; 
          eySum[i] =  sqrt(pow(eySum[i], 2) + pow((g3->GetEY())[i], 2));
        }
        else { cout << "no sum" << endl; }
        cout << "g3   " << setw(6) << (g3->GetX())[i] 
        << setw(10) << (g3->GetY())[i] << setw(10) << (g3->GetEY())[i]
        << setw(10) << xSum[i] << setw(10) << ySum[i] << setw(10) << eySum[i] << endl;
      }
      gSum->Write();
    }
// *****************************************************************************
// *****************************************************************************
// Kr    
    if (gasName=="Kr")
    {
// retrieve graphs
      int const nGraphs = 3;  // number of graphs
      TGraphErrors *g1=(TGraphErrors*)tf->Get("Kr--Krplus_Knudsen");
      TGraphErrors *g2=(TGraphErrors*)tf->Get("Kr--Krplusplus_Knudsen");
      TGraphErrors *g3=(TGraphErrors*)tf->Get("Kr--Krplusplusplus_Knudsen");
      
      TGraphErrors *gSum = new TGraphErrors(g1->GetN());
      gSum->SetName((gasName + "--All").c_str());
      gSum->SetTitle(("pbar/" + gasName + " total cross section").c_str());
      double *xSum = gSum->GetX();
      double *ySum = gSum->GetY();
      double *exSum = gSum->GetEX();
      double *eySum = gSum->GetEY();
      
      cout << endl;
      cout << "g    " << setw(6) << "gx[i]" 
      << setw(10) << "gy[i]" << setw(10) << "gey[i]"
      << setw(10) << "xSum[i]" << setw(10) << "ySum[i]" << setw(10) 
      << "eySum[i]" << endl;
      for (int i=0; i<gSum->GetN(); ++i)
      {
        cout << endl;
        
        xSum[i] = (g1->GetX())[i];
        ySum[i] = (g1->GetY())[i];
        eySum[i] = (g1->GetEY())[i];
        cout << "g1   " << setw(6) << (g1->GetX())[i] 
        << setw(10) << (g1->GetY())[i] << setw(10) << (g1->GetEY())[i]
        << setw(10) << xSum[i] << setw(10) << ySum[i] << setw(10) << eySum[i] << endl;
        
        if (i < g2->GetN())
        { 
          ySum[i] +=  (g2->GetY())[i]; 
          eySum[i] =  sqrt(pow(eySum[i], 2) + pow((g2->GetEY())[i], 2));
        }
        else { cout << "no sum" << endl; }
        cout << "g2   " << setw(6) << (g2->GetX())[i] 
        << setw(10) << (g2->GetY())[i] << setw(10) << (g2->GetEY())[i]
        << setw(10) << xSum[i] << setw(10) << ySum[i] << setw(10) << eySum[i] << endl;
        
        if (i < g3->GetN())
        { 
          ySum[i] +=  (g3->GetY())[i]; 
          eySum[i] =  sqrt(pow(eySum[i], 2) + pow((g3->GetEY())[i], 2));
        }
        else { cout << "no sum" << endl; }
        cout << "g3   " << setw(6) << (g3->GetX())[i] 
        << setw(10) << (g3->GetY())[i] << setw(10) << (g3->GetEY())[i]
        << setw(10) << xSum[i] << setw(10) << ySum[i] << setw(10) << eySum[i] << endl;
      }
      gSum->Write();
    }
// *****************************************************************************
// *****************************************************************************
// Xe    
    if (gasName=="Xe")
    {
// retrieve graphs
      int const nGraphs = 3;  // number of graphs
      TGraphErrors *g1=(TGraphErrors*)tf->Get("Xe--Xeplus_Knudsen");
      TGraphErrors *g2=(TGraphErrors*)tf->Get("Xe--Xeplusplus_Knudsen");
      TGraphErrors *g3=(TGraphErrors*)tf->Get("Xe--Xeplusplusplus_Knudsen");
      
      TGraphErrors *gSum = new TGraphErrors(g1->GetN());
      gSum->SetName((gasName + "--All").c_str());
      gSum->SetTitle(("pbar/" + gasName + " total cross section").c_str());
      double *xSum = gSum->GetX();
      double *ySum = gSum->GetY();
      double *exSum = gSum->GetEX();
      double *eySum = gSum->GetEY();
      
      cout << endl;
      cout << "g    " << setw(6) << "gx[i]" 
      << setw(10) << "gy[i]" << setw(10) << "gey[i]"
      << setw(10) << "xSum[i]" << setw(10) << "ySum[i]" << setw(10) 
      << "eySum[i]" << endl;
      for (int i=0; i<gSum->GetN(); ++i)
      {
        cout << endl;
        
        xSum[i] = (g1->GetX())[i];
        ySum[i] = (g1->GetY())[i];
        eySum[i] = (g1->GetEY())[i];
        cout << "g1   " << setw(6) << (g1->GetX())[i] 
        << setw(10) << (g1->GetY())[i] << setw(10) << (g1->GetEY())[i]
        << setw(10) << xSum[i] << setw(10) << ySum[i] << setw(10) << eySum[i] << endl;
        
        if (i < g2->GetN())
        { 
          ySum[i] +=  (g2->GetY())[i]; 
          eySum[i] =  sqrt(pow(eySum[i], 2) + pow((g2->GetEY())[i], 2));
        }
        else { cout << "no sum" << endl; }
        cout << "g2   " << setw(6) << (g2->GetX())[i] 
        << setw(10) << (g2->GetY())[i] << setw(10) << (g2->GetEY())[i]
        << setw(10) << xSum[i] << setw(10) << ySum[i] << setw(10) << eySum[i] << endl;
        
        if (i < g3->GetN())
        { 
          ySum[i] +=  (g3->GetY())[i]; 
          eySum[i] =  sqrt(pow(eySum[i], 2) + pow((g3->GetEY())[i], 2));
        }
        else { cout << "no sum" << endl; }
        cout << "g3   " << setw(6) << (g3->GetX())[i] 
        << setw(10) << (g3->GetY())[i] << setw(10) << (g3->GetEY())[i]
        << setw(10) << xSum[i] << setw(10) << ySum[i] << setw(10) << eySum[i] << endl;
      }
      gSum->Write();
    }
// *****************************************************************************

  }
  tf->Close();
}
