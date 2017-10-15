/// elenaCool.cc
/// ----------------------------------------------------------------------------
/// to deal with cooling simulation files from Gerard Tranquille (betaCool)
/// Pierre Grandemange
/// 07.07.2017
/// ----------------------------------------------------------------------------


/// includes and namespaces
/// ----------------------------------------------------------------------------
/// standard library
#include <iostream>
#include <fstream>
#include <stdlib.h>     /// getenv
//#include <iomanip>
/// root classes
#include "TAxis.h"
#include "TGraph.h"
#include "TPad.h"
/// rep classes
/// rep namespaces
/// local functions
#include "elenaCool.h"
using namespace std;



/// ecool(string const& coolingStep)
/// ----------------------------------------------------------------------------
/// make ecooling graph from Gerard Tranquille simulation file (betaCool)
TGraph* ecool(string const& coolingStep)
{
  /// directory path
  const string home = getenv("HOME");
  const string dir("/private/programming/root/repShared/repData/coolingData/");
  string path(home+dir);
  /// file path
  string file("");
  if( coolingStep=="35MeV_c" ) {
    file="ELENA_emittance_35MeVc.cur";
  }
  else if( coolingStep=="100keV" ) {
    file="ELENA_emittance_100keV.cur";
  }
  else {
    cout << "rep WARNING - <ecool> function - do not know input: ";
    cout << coolingStep << endl;
  }
  
  /// tgrap and filestream
  TGraph* tg = new TGraph();
  ifstream stream;
  //stream.open((path + file).c_str());
  string fp = path + file;
  //string fp = file;
  stream.open(fp.c_str());
  if( !stream.is_open() ) { 
    cout << "rep WARNING - <ecool> function - opening file failed: " << endl;
    //cout << "->  " << (path+file).c_str() << endl; 
    cout << "->  " << fp << endl;
  }
  else {
    const int nVar(9);
    double var[9];
    int nPoints(0);
    while( !stream.eof() ) { /// read file till the end
      for (int i=0; i<nVar; ++i) {
        stream >> var[i];
        //cout << setw(12) << var[i] << "   " ;
      }
      tg->SetPoint(nPoints, var[0], var[1]/6.); /// see Gerard mail with data
      ++nPoints;
      //cout << endl;
    }
    stream.close();
  }
  tg->SetName("eecg");
  tg->SetTitle(file.c_str());
  tg->GetXaxis()->SetTitle("time  (s)");
  tg->GetYaxis()->SetTitle("???  (???)");
  return tg;
}







/// scool(string const& coolingStep, double xi, double yi double yf)
/// ----------------------------------------------------------------------------
/// get scaled graph from ecool - scaled ecool - scool
/// scale curve to input param zi and zf, and translate to xi
TGraph* scool(string const& coolingStep, double xi, 
              double zi, double zf)
{
  /// get graph from ecool
  TGraph* tg = ecool(coolingStep);
  /// get input graph and parameters of interest
  const int n = tg->GetN();
  double* y = tg->GetY();
  const double yi = y[0];
  const double yf = y[n-1];
  /// create a scaled graph with max=yi et min=yf
  TGraph* scg = new TGraph(n);
  double z(0);
  for(int i=0; i<n; ++i)
  {
    z = zi + (zi-zf)/(yi-yf) * (y[i]-yi);
    scg->SetPoint(i, xi+(tg->GetX())[i], z);
  }
  /// free the memory and return
  delete tg; tg=0;
  return scg;
}






/// vcool
/// ----------------------------------------------------------------------------
/// eval y(x) on electron cooling scaled curve
double vcool(TGraph* const& scool, double const& x)
{ return scool->Eval(x); }




/// ecoolplot
/// ----------------------------------------------------------------------------
/// plot ecool graph
void ecoolplot(string pathOut)
{
  string file("");
  TGraph* eco = ecool("35MeV_c");
  eco->SetMarkerStyle(4);
  eco->SetMarkerSize(0.4);
  eco->Draw("ALP");
  file = pathOut + "35MeV_c.pdf";
  gPad->Print(file.c_str(), "pdf");
  TGraph* ect = ecool("100keV");
  ect->SetMarkerStyle(4);
  ect->SetMarkerSize(0.4);
  ect->Draw("ALP");
  file = pathOut + "100keV.pdf";
  gPad->Print(file.c_str(), "pdf");
  delete eco; eco=0;
  delete ect; ect=0;
}





