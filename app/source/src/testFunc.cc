/// localFunctions.cc
/// ----------------------------------------------------------------------------
/// local functions specific to this program folder
/// Pierre Grandemange
/// 09.06.2017
/// ----------------------------------------------------------------------------


/// includes and namespaces
/// ----------------------------------------------------------------------------
/// standard library
#include <iostream>
#include <fstream>
//#include <iomanip>
/// root classes
#include "TAxis.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TMath.h"
#include "TObject.h"
#include "TPad.h"
#include "TRandom3.h"
/// rep classes
#include "RAtom.h"
#include "RBeam.h"
#include "RFlyer.h"
#include "RGas.h"
#include "RMachine.h"
#include "RMolecule.h"
#include "RObject.h"
#include "RParticle.h"
#include "RUnit.h"
/// rep namespaces
#include "RDump.h"
#include "RMath.h"
#include "RParticleList.h"
#include "RUnits.h"
/// rep data functions
#include "elenaCool.h" /// for ecooling
/// local functions
#include "testFunc.h"
using namespace std;
//using namespace RParticleList;
//using namespace RUnits;





/// couple beams, cycles and interactions
/// ****************************************************************************


/// testNewGFB() - 20170920
/// ----------------------------------------------------------------------------
/// test new gas/flyer/beam class 
void testNewGFB()
{
  /// function intro
  /// --------------------------------------------------------------------------
  /// computation time recorder start
  const clock_t chronoStart = clock();  /// chrono start
  double chrono(0);                      /// chrono initialisation
  /// function header
  cout << endl;   cout << "testNewGFB() - 20/09/2017" << endl;
  RDump::line();  cout << "test new gas/flyer/beam class" << endl;
  
  /// debug iostream parameters
  /// --------------------------------------------------------------------------
  //streamsize iop(cout.precision());
  
  /// dump iostream
  //iop = cout.precision();
  //cout << endl << "iop = " << iop << endl; 
  
  /// dump chrono
  chrono = 1000.*double(clock()-chronoStart)/CLOCKS_PER_SEC;
  cout << endl << "chrono = " << chrono << "ms" << endl; 
  
  
  
  /// simple pure H2 gas definition
  /// --------------------------------------------------------------------------
  /// gas properties
  double pre(4.e-12*RUnits::mbar);  /// overall pressure
  double tem(298.*RUnits::K);       /// overall temperature
  /// simple H2 gas
  RGas* gas = new RGas((RMolecule*)RParticleList::mH2(), pre, tem);
  
  /// dump  values
  cout << endl << endl << "define simple H2 gas " << endl;
  RDump::line(20, "-");
  gas->dump(cout);
  
  /// dump iostream
  //iop = cout.precision();
  //cout << endl << "iop = " << iop << endl; 
  
  /// dump chrono
  chrono = 1000.*double(clock()-chronoStart)/CLOCKS_PER_SEC;
  cout << endl << "chrono = " << chrono << "ms" << endl; 
  
  
  /// define machine (tbp <=> Twiss Beta Parameter)
  /// --------------------------------------------------------------------------
  string macName("elena");
  double macLen(30.4055*RUnits::m);	    /// machine length
  double macAcceptance(67.*RUnits::mm);	/// machine acceptance
  double macTwissB(3.*RUnits::m);				/// twiss beta parameter at detector
  RMachine* mac = new RMachine(macName, macLen, macAcceptance, macTwissB);
  
  /// dump  values
  cout << endl << endl << "define decelerator " << endl;
  RDump::line(18, "-");
  //cout << endl; 
  mac->dump();
  
  /// dump iostream
  //iop = cout.precision();
  //cout << endl << "iop = " << iop << endl; 
  
  /// dump chrono
  chrono = 1000.*double(clock()-chronoStart)/CLOCKS_PER_SEC;
  cout << endl << "chrono = " << chrono << "ms" << endl; 
  
  
  /// beam definition
  /// --------------------------------------------------------------------------
  /// beam name
  const string bna="testBeam";
  /// number of particles
  const int npa=1000; 
  /// particles 
  RParticle* par = (RParticle*)RParticleList::hpbar();
  /// distribution types
  string dis[6] = {"flat", "gauss", "gauss", "flat", "flat", "flat"};
  /// beam emittance
  double emx( 1.*RUnits::um * 1.*RUnits::one );  /// x emittance
  double emy( 9.*RUnits::um * 1.*RUnits::one );  /// y emittance
  double emz( 25.*RUnits::one * 1.*RUnits::MeV );  /// z emittance
  double emi[3] = { emx, emy, emz };
  /// machine parameters
  double mpx( 1.*RUnits::m );  /// twiss beta x
  double mpy( 1.*RUnits::m );  /// twiss beta y
  double mpz( 1.*RUnits::one );  /// ???
  double mpa[3] = { mpx, mpy, mpz };
  /// reference particle momentum
  double mv0 = 100.*RUnits::MeV_c;
  /// beam size
  RBeam* bea = new RBeam(bna, npa, par, dis, mv0, emi, mpa);
  
  /// dump values
  cout << endl << endl << "beam definition" << endl; RDump::line(15, "-");
  bea->dumpLine();
  /*for( int i=0; i<npa; ++i ) {
    RFlyer fly = bea->getFlyColl()[i];
    cout << setw(2) << left << i << " "; 
    fly.line6D();
  }*/
  
  /// dump chrono
  chrono = 1000.*double(clock()-chronoStart)/CLOCKS_PER_SEC;
  cout << endl << "chrono = " << chrono << "ms" << endl; 
  
  /// dump values 2! 
  cout << endl; bea->dumpFlyColl(200);
  
  /// dump iostream
  //iop = cout.precision();
  //cout << endl << "iop = " << iop << endl; 
  
  /// dump chrono
  chrono = 1000.*double(clock()-chronoStart)/CLOCKS_PER_SEC;
  cout << endl << "chrono = " << chrono << "ms" << endl; 
  
  /// dump root data folder (new RDump function using getenv)
  //cout << endl;
  //cout << getenv("repData") << endl; 
  //cout << getenv("repNamespaces") << endl; 
  //cout << getenv("repObjects") << endl; 
  //cout << RMolecule::dataFolder() << endl;
  //cout << RMolecule::crossSectionFile() << endl;
  
  
  /// beam/gas interaction rates
  /// --------------------------------------------------------------------------
  cout << endl << endl << "beam/gas interaction rates" << endl; 
  RDump::line(26, "-");
  double* rat = gas->rsRates(bea, mac);
  cout << "ion_pbarnu" << "  " << rat[0]*1./RUnits::m*1./RUnits::s
       << "+/-" << rat[0+4]*1./RUnits::m*1./RUnits::s
       << "/"+RUnits::m.s()+"/"+RUnits::s.s() << endl;
  cout << "rs_ebunu  " << "  " << rat[1]*1./RUnits::m*1./RUnits::s
       << "+/-" << rat[1+4]*1./RUnits::m*1./RUnits::s
       << "/"+RUnits::m.s()+"/"+RUnits::s.s() << endl;
  cout << "rs_lalnu  " << "  " << rat[2]*1./RUnits::m*1./RUnits::s
       << "+/-" << rat[2+4]*1./RUnits::m*1./RUnits::s
       << "/"+RUnits::m.s()+"/"+RUnits::s.s() << endl;
  cout << "rs_totnu  " << "  " << rat[3]*1./RUnits::m*1./RUnits::s
       << "+/-" << rat[3+4]*1./RUnits::m*1./RUnits::s
       << "/"+RUnits::m.s()+"/"+RUnits::s.s() << endl;
  
  /// dump chrono
  chrono = 1000.*double(clock()-chronoStart)/CLOCKS_PER_SEC;
  cout << endl << "chrono = " << chrono << "ms" << endl; 
  
  
  /// function outro
  /// --------------------------------------------------------------------------
  /// computation time recorder end
  chrono = 1000.*double(clock()-chronoStart)/CLOCKS_PER_SEC;
  /// outro header
  cout << endl; cout << "testNewGFB() function end " << endl; RDump::line();
  cout << "function computation time = " << chrono << "ms" << endl; 
}








/// testBeamDist() - 20170921
/// ----------------------------------------------------------------------------
/// test beam distributions with new flyers definition
void testBeamDist()
{
  /// function intro
  /// --------------------------------------------------------------------------
  /// computation time recorder start
  const clock_t timeFuncStart = clock();
  /// function header
  cout << endl;   cout << "testBeamDist() - 21/09/2017" << endl;
  RDump::line();  cout << "test beam distributions with new flyer." << endl;
  
  
  /// do the stuff
  /// --------------------------------------------------------------------------
  /// beam name
  const string bna="firstBeam";
  /// number of particles
  const int npa=1e4; 
  /// particles 
  RParticle* par = (RParticle*)RParticleList::hpbar();
  /// distribution types
  string dis[6] = {"flat", "gauss", "gauss", "flat", "flat", "flat"};
  /// beam size
  double dx(  1.*RUnits::mm );
  double dy(  2.*RUnits::mm );
  double dz(  3.*RUnits::mm );
  double dpx( 4.*RUnits::MeV_c );
  double dpy( 5.*RUnits::MeV_c );
  double dpz( 6.*RUnits::MeV_c );
  double siz[6] = { dx, dy, dz, dpx, dpy, dpz };
  /// beam
  ///   -> defined by beam size OR by (beam emittance AND machine parameters)
  RBeam* bea = new RBeam(bna, npa, par, dis, siz);
 
  /// dump computation time
  cout << endl << "computation time" << endl; RDump::line(16, "-");
  float t1(1000*float(clock()-timeFuncStart)/CLOCKS_PER_SEC);
  cout << "t1 = " << t1 << "ms" << endl; 
  
  /// check the stuff
  /// --------------------------------------------------------------------------
  /// beam dumpLine method
  cout << endl << "beam dumpLine method" << endl; RDump::line(20, "-");
  bea->dumpLine();
  /// dump 6D sizes 
  cout << endl << "dump beam sizes" << endl; RDump::line(15, "-");
  for( int i=0; i<6; ++i ) { 
    double size( bea->getSize()[i] );
    if( i<3 ) cout << i << "  " << size/RUnits::mm << " mm" << endl; 
    else      cout << i << "  " << size/RUnits::MeV_c << " MeV_c" << endl;  
  }
  
  /*/// dump RFlyer collection
  cout << endl << "dump RFlyer collection" << endl; RDump::line(22, "-"); 
  for(int i; i<npa; ++i) { 
    cout << i << "  "; 
    //RFlyer fly = bea->getFlyColl()[i];
    //fly.line6D();
    bea->getFlyColl()[i].line6D();
  }*/
  
  /// dump integrals
  cout << endl << "dump integrals " << endl; RDump::line(14, "-");
  for( int i=0; i<6; ++i ) {
    TF1   f( bea->getFunc()[i] );
    TH1D  h( bea->getHist()[i] );
    double xM( f.GetXmax() );
    double xm( f.GetXmin() );
    cout << "h[i].Int() "  << left << setw(6) << h.Integral()        << " ";
    cout << "h[i].Int(w) " << left << setw(7) << h.Integral("width") << " ";
    cout << "f[i].Int() "  << left << setw(8) << f.Integral(xm, xM)  << endl;
  }
  
  /// output pdf/canvas/plot parameters
  string pdfN     = "../results/20170921/testBeam.pdf";
  string pdfOpen  = pdfN + "[";
  string pdfClose = pdfN + "]";
  double cw(864); /// can width  = 0.95 * 1920/2
  double ch(486); /// can height = 0.95 * 1080/2 
  TCanvas* can = new TCanvas("can", "testBeamPDFoutput",cw, ch);
  can->SetWindowSize( cw+(cw-can->GetWw()), ch+(ch-can->GetWh()) );
  
  /// plot dist in pdf
  cout << endl << "plot functions and distributions in a pdf file" << endl;
  RDump::line(46, "-");
  can->Print(pdfOpen.c_str()); 
  for( int i=0; i<6; ++i ) { 
    /// func drawing options
    bea->getFunc()[i].Draw(); 
    /// hist drawing options
    //bea->getHist()[i].Draw("hist"); 
    bea->getHist()[i].Draw("same"); 
    //bea->getHist()[i].SetLineWidth(2); 
    /// canvas print
    can->Print(pdfN.c_str()); 
  }
  can->Print(pdfClose.c_str()); 
  
  /// test RUnit/RUnits
  //double nt=2.*RUnits::mm;
  //cout << nt << "   " << nt/RUnits::mm << RUnits::mm.getSymbol() << endl;
 
  /// dump computation time
  //cout << endl << "computation time" << endl; RDump::line(16, "-");
  //float t2(1000*float(clock()-timeFuncStart)/CLOCKS_PER_SEC);
  //cout << "t2 = " << t2 << "ms" << endl; 
  
  
  /// function outro
  /// --------------------------------------------------------------------------
  /// computation time recorder end
  float timeFuncDur(1000*float(clock()-timeFuncStart)/CLOCKS_PER_SEC);
  /// outro header
  cout << endl; cout << "testBeamDist() function end " << endl; RDump::line();
  cout << "function computation time = " << timeFuncDur << "ms" << endl; 
}













/// testRMath() - 20170921
/// ----------------------------------------------------------------------------
/// test RMath namespace new functions <mean> and <standardDev>
void testRMath()
{
  /// function intro
  /// --------------------------------------------------------------------------
  /// computation time recorder start
  const clock_t timeFuncStart = clock();
  /// function header
  cout << endl;   cout << "testRMath() - 21/09/2017" << endl;
  RDump::line();  cout << "test new RMath functions" << endl;
  
  
  /// test RMath namepace new functions
  /// --------------------------------------------------------------------------
  /// distribution function parameters
  string const dfn = "fest";        /// distribution function name 
  string const dft = "gauss";       /// distribution function type
  RUnit  const dxu = RUnits::um;    /// distribution x axis unit
  RUnit  const dx2 = RUnits::mm;
  double const dfi = 1.*dx2;      /// distribution function integral
  double const dfm = 5.*dx2;      /// distribution function mean
  double const dfs = 4.*dx2;      /// distribution function size
  double const dfp[3] = {dfi, dfm, dfs};
  /// distribution function
  TF1 fun = RMath::rf1(dfn, dft, dxu, (double*)&dfp[0]);
  cout << "debug: funxaxis = " << fun.GetXaxis()->GetTitle() << endl;
  /// double collection
  int const dcn(1.e4);                     /// distribution collection nsample
  //int const dcn(10);                     /// distribution collection nsample
  double* dco = RMath::reprand(&fun, dxu, dcn, 0); /// dist double collection
  double* cms = RMath::mpms(dco, dcn);     /// dist collec mean and sigma
  /// dist histogram parameters
  string  const dhn("hest");        /// distribution histogram name
  int     const dhb(100);           /// distribution histogram nbins
  double  const hxm(fun.GetXmin()); /// dist hist xmin
  double  const hxM(fun.GetXmax()); /// dist hist xmax
  /// dist histogram 
  TH1D his = RMath::rh1(dhn, dco, dcn, dhb, dxu, hxm, hxM); /// dist hist 
  double const hme( his.GetMean()*dxu );    /// histogram mean 
  double const hsd( his.GetStdDev()*dxu );  /// histogram stdDev
  
  /// dump double collection
  //cout << endl << "dump double collection" << endl; RDump::line(22, "-"); 
  //for(int i; i<dcn; ++i) {
  //  cout << i << "  " << dco[i]/dx2 << dx2.getSymbol() << endl; 
  //}
  
  /// dump mean and sigma
  cout << endl << "dump mean and sigma" << endl; RDump::line(19,"-");
  cout << "<f(x)> = " << "TF1 " << "method..." << endl; 
  cout << "<y(x)> = " << cms[0]/dx2 << "+/-" << cms[1]/dx2 << dx2.s() << endl;
  cout << "<h(x)> = " << hme/dxu    << "+/-" << hsd/dxu    << dxu.s() << endl;
  
  /// canvas
  /// --------------------------------------------------------------------------
  /// output pdf/canvas/plot parameters
  double cw(864); /// can width  = 0.95 * 1920/2
  double ch(486); /// can height = 0.95 * 1080/2 
  TCanvas can("can", "testRMath",cw, ch);
  can.SetWindowSize( cw+(cw-can.GetWw()), ch+(ch-can.GetWh()) );
  
  /// pdf
  /// --------------------------------------------------------------------------
  string pdfN     = "../results/20170921/testRMath.pdf";
  //string pdfOpen  = pdfN + "[";
  //string pdfClose = pdfN + "]";
  cout << endl << "pdf file writting" << endl;
  RDump::line(26, "-");
  /// function drawing options
  cout << "debug: funxaxis = " << fun.GetXaxis()->GetTitle() << endl;
  RMath::fhscale(fun, his);
  cout << "debug: funxaxis = " << fun.GetXaxis()->GetTitle() << endl;
  fun.Draw(); 
  cout << "debug: funxaxis = " << fun.GetXaxis()->GetTitle() << endl;
  /// histograms drawing options
  //his.Draw(); 
  his.Draw("same"); 
  //his.SetLineWidth(2); 
  /// canvas print
  can.Print(pdfN.c_str(), "pdf"); 
  
  
  /// function outro
  /// --------------------------------------------------------------------------
  /// computation time recorder end
  float timeFuncDur(1000*float(clock()-timeFuncStart)/CLOCKS_PER_SEC);
  /// outro header
  cout << endl; cout << "testRMath() function end " << endl; RDump::line();
  cout << "function computation time = " << timeFuncDur << "ms" << endl; 
}










/// testNewFlyers() - 20170916
/// ----------------------------------------------------------------------------
/// test the new flyer definition, with [x,y,z,px,py,pz]
void testNewFlyers()
{
  /// function header
  /// --------------------------------------------------------------------------
  cout << endl;   cout << "testNewFlyers() - 16/09/2017" << endl;
  RDump::line();  cout << "test the new flyer definition." << endl;
  
  
  /// do the stuff
  /// --------------------------------------------------------------------------
  //RFlyer* fly = new RFlyer((RParticle*)RParticleList::hrep(), 0.1);
  //RFlyer* fly = new RFlyer((RParticle*)RParticleList::hpbar(), 0.1);
  double sixD[6] = { 3.*RUnits::mm,      4.*RUnits::mm,    0.*RUnits::mm,
                     100.*RUnits::MeV_c, 0.*RUnits::MeV_c, 0.*RUnits::MeV_c };
  RFlyer* fly = new RFlyer((RParticle*)RParticleList::hpbar(), sixD);
  cout << endl; fly->dump();
  
  //cout << endl;
  //double* pos6D = fly->getPos6D();
  //for( int i=0; i<6; ++i ) {
  //  if(i<=2)  cout << pos6D[i]/RUnits::m      << " m" << endl;
  //  if(i>2)   cout << pos6D[i]/RUnits::keV_c  << " keV_c" << endl;
  //}
  
  cout << endl << "RFlyer* fly = new RFlyer(RPArticleList(), pos6D[])";
  cout << endl; fly->dump6D();
  cout << "distance " << fly->distance()/RUnits::mm     << " mm"     << endl;
  cout << "momentum " << fly->momentum()/RUnits::MeV_c  << " MeV/c"  << endl;
  cout << "beta     " << fly->beta()                    << " "      << endl;
  cout << "gamma    " << fly->gamma()                   << " "      << endl;
  cout << "e0       " << fly->e0()/RUnits::MeV          << " MeV"    << endl;
  cout << "ekin     " << fly->ekin()/RUnits::keV        << " keV"    << endl;
  cout << "velocity " << fly->velocity()/RUnits::m_us   << " m/us"   << endl;
  
  fly->setPos6D(0., 0., 0., 0., 0., 13.7*RUnits::MeV_c);
  cout << endl << "fly->setPos6D(0., 0., 0., 0., 0., 100.*RUnits::MeV_c);";
  cout << endl; fly->dump6D();
  cout << "distance " << fly->distance()/RUnits::mm     << " mm"     << endl;
  cout << "momentum " << fly->momentum()/RUnits::MeV_c  << " MeV/c"  << endl;
  cout << "beta     " << fly->beta()                    << " "      << endl;
  cout << "gamma    " << fly->gamma()                   << " "      << endl;
  cout << "e0       " << fly->e0()/RUnits::MeV          << " MeV"    << endl;
  cout << "ekin     " << fly->ekin()/RUnits::keV        << " keV"    << endl;
  cout << "velocity " << fly->velocity()/RUnits::m_us   << " m/us"   << endl;
  
  fly->ekin(5.3*RUnits::keV);
  cout << endl << "fly->ekin(100.*RUnits::keV);";
  cout << endl; fly->dump6D();
  cout << "distance " << fly->distance()/RUnits::mm     << " mm"     << endl;
  cout << "momentum " << fly->momentum()/RUnits::MeV_c  << " MeV/c"  << endl;
  cout << "beta     " << fly->beta()                    << " "      << endl;
  cout << "gamma    " << fly->gamma()                   << " "      << endl;
  cout << "e0       " << fly->e0()/RUnits::MeV          << " MeV"    << endl;
  cout << "ekin     " << fly->ekin()/RUnits::keV        << " keV"    << endl;
  cout << "velocity " << fly->velocity()/RUnits::m_us   << " m/us"   << endl;
  
  fly->momentum(13.7*RUnits::MeV_c);
  cout << endl << "fly->ekin(100.*RUnits::keV);";
  cout << endl; fly->dump6D();
  cout << "distance " << fly->distance()/RUnits::mm     << " mm"     << endl;
  cout << "momentum " << fly->momentum()/RUnits::MeV_c  << " MeV/c"  << endl;
  cout << "beta     " << fly->beta()                    << " "      << endl;
  cout << "gamma    " << fly->gamma()                   << " "      << endl;
  cout << "e0       " << fly->e0()/RUnits::MeV          << " MeV"    << endl;
  cout << "ekin     " << fly->ekin()/RUnits::keV        << " keV"    << endl;
  cout << "velocity " << fly->velocity()/RUnits::m_us   << " m/us"   << endl;
  
  
  /// end of function
  /// --------------------------------------------------------------------------
  cout << endl; cout << "testNewFlyers() function end " << endl; 
  RDump::line();
}





/// testBeamProfiles() - 20170714
/// ----------------------------------------------------------------------------
/// test beam profiles plots
/*void testBeamProfiles()
{
	/// function header
  /// --------------------------------------------------------------------------
  cout << endl;   cout << "testBeamProfiles() - 06/07/2017" << endl;
  RDump::line();  cout << "test beam profiles plots" << endl;
  
  
  /// do the stuff : declare
  /// --------------------------------------------------------------------------
  /// flyer
  RParticle* par = (RParticle*)RParticleList::hpbar();
  double beta(RFlyer::ekinToBeta(5.3*RUnits::MeV, par->getMass()));
  RFlyer* fly= new RFlyer(par, beta);
  
  /// beam
  string beaName("elenaBeam");
  double beaNum(2.e7);
  double beaEmit[3] = { 1.1*RUnits::um, 2.2*RUnits::um, 
    3.3*RUnits::keV*RUnits::ns };
  string beaDist("gauss");
  
  RBeam* bea = new RBeam(beaName, beaNum, fly);
  bea->setEmittance(beaEmit);
  bea->setDist(beaDist);
  
  /// define machine (tbp <=> Twiss Beta Parameter)
  string macName("elena");
  double macLen(30.4055*RUnits::m);	    /// machine length
  double macAcceptance(67.*RUnits::mm);	/// machine acceptance
  double macTwissB(3.*RUnits::m);				/// twiss beta parameter at detector
  RMachine* mac = new RMachine(macName, macLen, macAcceptance, macTwissB);
  
  /// define gas
  double gasPre(4.e-12*RUnits::mbar);
  double gasTem(300.);
  RGas* gas = new RGas((RMolecule*)RParticleList::mN2(), gasPre, gasTem);
  
  
  /// do the stuff : plot distributions
  /// --------------------------------------------------------------------------
  TCanvas* can = new TCanvas("can");
  double sigx(RBeam::sigmab(beaEmit[0], macTwissB));
  double darea(0);   /// check distribution area
  cout << endl << "sigx = " << sigx/RUnits::mm << "mm" << endl;
  string dnam("");  /// distribution name
  
  /// gaussian dist
  dnam = "gauss";
  TF1* tfd1 = RBeam::dist(sigx, dnam); 
  darea = tfd1->Integral(-4.*sigx, 4.*sigx);
  cout << left << setw(6) << dnam << " | " 
       << "darea = " << darea << " noUnit" << endl; 
  tfd1->Draw();
  tfd1->SetLineColor(kBlue);
  /// flat dist
  dnam = "flat";
  TF1* tfd2 = RBeam::dist(sigx, dnam); 
  darea = tfd2->Integral(-4.*sigx, 4.*sigx);
  cout << left << setw(6) << dnam << " | " 
    << "darea = " << darea << " noUnit" << endl; 
  tfd2->Draw("same"); 
  gPad->Print("dist.pdf", "pdf");
  delete tfd1; tfd1=0; delete tfd2; tfd2=0;
  
  /// do the stuff : plot beam profiles
  /// --------------------------------------------------------------------------
  int nbins(100);
  
  TH1D* thp1 = bea->x(mac, nbins); thp1->Draw();
  //gPad->Print("prof_x_gauss.pdf", "pdf"); tho=0;
  bea->setDist("flat"); 
  TH1D* thp2 = bea->x(mac, nbins); thp2->Draw("same"); 
  thp2->SetLineColor(2);
  gPad->Print("prof_x.pdf", "pdf");
  delete thp1; thp1=0; delete thp2; thp2=0; 
  
  
  /// draw ionisation;
  double ionLen(50.*RUnits::mm);
  double ionTim(100.*RUnits::ms);
  
  
  /// end of function
  /// --------------------------------------------------------------------------
  cout << endl; cout << "testBeamProfiles() function end " << endl; 
  RDump::line();
}
*/

















/// elena cycle functions
/// ****************************************************************************

/// ecplots()
/// ----------------------------------------------------------------------------
/// elena cycle plots - mv, ke, ex and ey
void ecplots()
{
  /// function header
  /// --------------------------------------------------------------------------
  cout << endl;   cout << "ecplots() - 06/07/2017" << endl;
  RDump::line();  cout << "ionistations/m/s as function of time" << endl;
  
  
  /// do the stuff
  /// --------------------------------------------------------------------------
  int const np(53);         /// number of points
  TGraph* mvg = ecmvg(np);   /// momentum cycle
  TGraph* keg = eckeg(np);   /// energy cycle
  TGraph* exg = ecxeg(np);   /// x emittance cycle
  TGraph* eyg = ecyeg(np);   /// y emittance cycle
  
  /// dump stuff
  double* mvx = mvg->GetX();
  double* mvy = mvg->GetY();
  double* key = keg->GetY();
  double* exy = exg->GetY();
  double* eyy = eyg->GetY();
  cout << endl;
  cout << left;
  cout << setw(3)   << "i"          << "   ";
  cout << setw(10)  << "time (s)"   << "   ";
  cout << setw(10)  << "mv (MeV/c)" << "   ";
  cout << setw(10)  << "ke (MeV)"   << "   ";
  cout << setw(10)  << "ex (um)"    << "   ";
  cout << setw(10)  << "ey (um)"    << "   ";
  cout << endl;
  for( int i=0; i<np; ++i ) { 
    cout << setw(3)   << i      << "   ";
    cout << setw(10)  << mvx[i] << "   ";
    cout << setw(10)  << mvy[i] << "   ";
    cout << setw(10)  << key[i] << "   ";
    cout << setw(10)  << exy[i] << "   ";
    cout << setw(10)  << eyy[i] << "   ";
    cout << endl; 
  }
  
  /// draw stuff
  mvg->SetMarkerStyle(4);
  mvg->SetMarkerSize(0.4);
  mvg->Draw("ALP");
  gPad->Print("../results/20170706/mvg.pdf", "pdf");
  keg->SetMarkerStyle(4);
  keg->SetMarkerSize(0.4);
  keg->Draw("ALP");
  gPad->Print("../results/20170706/keg.pdf", "pdf");
  exg->SetMarkerStyle(4);
  exg->SetMarkerSize(0.4);
  exg->Draw("ALP");
  gPad->Print("../results/20170706/exg.pdf", "pdf");
  eyg->SetMarkerStyle(4);
  eyg->SetMarkerSize(0.4);
  eyg->Draw("ALP");
  gPad->Print("../results/20170706/eyg.pdf", "pdf");
  
  
  /// end of function
  /// --------------------------------------------------------------------------
  cout << endl; cout << "ecplots() function end " << endl; RDump::line();
}




/// ecxeg(int const np)
/// ----------------------------------------------------------------------------
/// elena cycle x emittance graph (x=horizontal)
TGraph* ecxeg(int const np)
{ 
  /// cycle steps 
  const int cns(10);  /// cycle number of steps
  double xsa[cns];    /// x steps array
  double ysa[cns];    /// y steps array
  
  /// time steps - tc = time cycle
  const double xUnit = 1.*RUnits::s;
  const double tc_inj( 0.*xUnit); /// injection;
  const double tc_fdb( 1.*xUnit); /// first deceleration begin
  const double tc_fde( 5.*xUnit); /// first deceleration end
  const double tc_feb( 6.*xUnit); /// first electron cooling start
  const double tc_fee(16.*xUnit); /// first electron cooling end
  const double tc_sdb(17.*xUnit); /// second deceleration begin
  const double tc_sde(21.*xUnit); /// second deceleration end
  const double tc_seb(22.*xUnit); /// second electron cooling start
  const double tc_see(25.*xUnit); /// second electron cooling end
  const double tc_ext(26.*xUnit); /// extraction (rebunched);
  
  /// x emittance steps
  const double yUnit = 1.*RUnits::um;
  const double inj(0.5*yUnit);
  const double fdb(0.5*yUnit);
  const double fde(1.7*yUnit);
  const double feb(1.7*yUnit);
  const double fee(0.45*yUnit);
  const double sdb(0.45*yUnit);
  const double sde(2.2*yUnit);
  const double seb(2.2*yUnit);
  const double see(0.3*yUnit);
  const double ext(1.2*yUnit);
  
  /// attribute values
  xsa[0]=tc_inj;  xsa[1]=tc_fdb;  xsa[2]=tc_fde;  xsa[3]=tc_feb; 
  xsa[4]=tc_fee;  xsa[5]=tc_sdb;  xsa[6]=tc_sde;  xsa[7]=tc_seb; 
  xsa[8]=tc_see;  xsa[9]=tc_ext; 
  
  ysa[0]=inj;   ysa[1]=fdb;   ysa[2]=fde;   ysa[3]=feb; 
  ysa[4]=fee;   ysa[5]=sdb;   ysa[6]=sde;   ysa[7]=seb; 
  ysa[8]=see;   ysa[9]=ext;
  
  /// output
  TGraph* tg = cg(np, cns, xsa, ysa, "cooling");
  double* x = tg->GetX();
  double* y = tg->GetY();
  for( int i=0; i<tg->GetN(); ++i ) { 
    tg->SetPoint(i, x[i]/xUnit, y[i]/yUnit);
  }
  tg->SetName("ex2t");
  tg->SetTitle("elena cycle - horizontal emittance");
  tg->GetXaxis()->SetTitle("time  (s)");
  tg->GetYaxis()->SetTitle("horizontal emittance  (um)");
  return tg;
}



/// ecyeg(int const np)
/// ----------------------------------------------------------------------------
/// elena cycle y emittance graph (x=horizontal)
TGraph* ecyeg(int const np)
{ 
  /// cycle steps 
  const int cns(10);  /// cycle number of steps
  double xsa[cns];    /// x steps array
  double ysa[cns];    /// y steps array
  
  /// time steps - tc = time cycle
  const double xUnit = 1.*RUnits::s;
  const double tc_inj( 0.*xUnit); /// injection;
  const double tc_fdb( 1.*xUnit); /// first deceleration begin
  const double tc_fde( 5.*xUnit); /// first deceleration end
  const double tc_feb( 6.*xUnit); /// first electron cooling start
  const double tc_fee(16.*xUnit); /// first electron cooling end
  const double tc_sdb(17.*xUnit); /// second deceleration begin
  const double tc_sde(21.*xUnit); /// second deceleration end
  const double tc_seb(22.*xUnit); /// second electron cooling start
  const double tc_see(25.*xUnit); /// second electron cooling end
  const double tc_ext(26.*xUnit); /// extraction (rebunched);
  
  /// y emittance steps
  const double yUnit = 1.*RUnits::um;
  const double inj(0.3*yUnit);
  const double fdb(0.3*yUnit);
  const double fde(1.1*yUnit);
  const double feb(1.1*yUnit);
  const double fee(0.42*yUnit);
  const double sdb(0.42*yUnit);
  const double sde(2.5*yUnit);
  const double seb(2.5*yUnit);
  const double see(0.2*yUnit);
  const double ext(0.75*yUnit);
  
  /// attribute values
  xsa[0]=tc_inj;  xsa[1]=tc_fdb;  xsa[2]=tc_fde;  xsa[3]=tc_feb; 
  xsa[4]=tc_fee;  xsa[5]=tc_sdb;  xsa[6]=tc_sde;  xsa[7]=tc_seb; 
  xsa[8]=tc_see;  xsa[9]=tc_ext; 
  
  ysa[0]=inj;   ysa[1]=fdb;   ysa[2]=fde;   ysa[3]=feb; 
  ysa[4]=fee;   ysa[5]=sdb;   ysa[6]=sde;   ysa[7]=seb; 
  ysa[8]=see;   ysa[9]=ext;
  
  /// output
  TGraph* tg = cg(np, cns, xsa, ysa, "cooling");
  double* x = tg->GetX();
  double* y = tg->GetY();
  for( int i=0; i<tg->GetN(); ++i ) { 
    tg->SetPoint(i, x[i]/xUnit, y[i]/yUnit);
  }
  tg->SetName("ex2t");
  tg->SetTitle("elena cycle - vertical emittance");
  tg->GetXaxis()->SetTitle("time  (s)");
  tg->GetYaxis()->SetTitle("vertical emittance  (um)");
  return tg;
}



/// ecekg(int const np)
/// ----------------------------------------------------------------------------
/// elena cycle kinetic energy graph 
/// carefull: ecmvg has a unit system "hardcoded"
/// -> assuming emittance growth is proportional to pomentum
TGraph* eckeg(int const np)
{
  TGraph* keg = new TGraph(np);
  TGraph* mvg = ecmvg(np);
  double* mvx = mvg->GetX();
  double* mvy = mvg->GetY();
  const double m0 = 1.0072765*RUnits::amu;
  for( int i=0; i<np; ++i ) {
    /// warning : unit system!
    double key( RFlyer::momentumToEkin(mvy[i]*RUnits::MeV_c, m0)/RUnits::MeV );
    keg->SetPoint(i, mvx[i] , key);
  }
  /// output
  keg->SetName("ke2t");
  keg->SetTitle("elena cycle - kinetic energy");
  keg->GetXaxis()->SetTitle("time  (s)");
  keg->GetYaxis()->SetTitle("kinetic energy  (MeV)");
  return keg;
}



/// ecmvg()
/// ----------------------------------------------------------------------------
/// elena cycle momentum graph
TGraph* ecmvg(int const np)
{ 
  /// cycle steps 
  const int cns(10);  /// cycle number of steps
  double xsa[cns];    /// x steps array
  double ysa[cns];    /// y steps array
  
  /// time steps - tc = time cycle
  const double xUnit = 1.*RUnits::s;
  const double tc_inj( 0.*xUnit); /// injection;
  const double tc_fdb( 1.*xUnit); /// first deceleration begin
  const double tc_fde( 5.*xUnit); /// first deceleration end
  const double tc_feb( 6.*xUnit); /// first electron cooling start
  const double tc_fee(16.*xUnit); /// first electron cooling end
  const double tc_sdb(17.*xUnit); /// second deceleration begin
  const double tc_sde(21.*xUnit); /// second deceleration end
  const double tc_seb(22.*xUnit); /// second electron cooling start
  const double tc_see(25.*xUnit); /// second electron cooling end
  const double tc_ext(26.*xUnit); /// extraction (rebunched);
  
  /// momentum steps - mv = momentum (mass*velocity)
  /// should correspond to magnetic cycle steps:
  /// -> B*rho # momentum (# = proportional)
  /// shoud be a linear dependance of time
  const double yUnit = 1.*RUnits::MeV_c;
  const double mv_inj( 100.*yUnit );  
  const double mv_fdb( 100.*yUnit );  
  const double mv_fde(  35.*yUnit );
  const double mv_feb(  35.*yUnit );  
  const double mv_fee(  35.*yUnit );
  const double mv_sdb(  35.*yUnit );
  const double mv_sde( 13.7*yUnit );
  const double mv_seb( 13.7*yUnit );
  const double mv_see( 13.7*yUnit );
  const double mv_ext( 13.7*yUnit );
  /// attribute values
  xsa[0]=tc_inj;  xsa[1]=tc_fdb;  xsa[2]=tc_fde;  xsa[3]=tc_feb; 
  xsa[4]=tc_fee;  xsa[5]=tc_sdb;  xsa[6]=tc_sde;  xsa[7]=tc_seb; 
  xsa[8]=tc_see;  xsa[9]=tc_ext; 
  
  ysa[0] = mv_inj; ysa[1] = mv_fdb; ysa[2] = mv_fde; ysa[3] = mv_feb; 
  ysa[4] = mv_fee; ysa[5] = mv_sdb; ysa[6] = mv_sde; ysa[7] = mv_seb; 
  ysa[8] = mv_see; ysa[9] = mv_ext;
  
  /// output
  TGraph* tg = cg(np, cns, xsa, ysa, "linear");
  double* x = tg->GetX();
  double* y = tg->GetY();
  for( int i=0; i<tg->GetN(); ++i ) { 
    tg->SetPoint(i, x[i]/xUnit, y[i]/yUnit);
  }
  tg->SetName("p2t");
  tg->SetTitle("elena cycle - momentum");
  tg->GetXaxis()->SetTitle("time  (s)");
  tg->GetYaxis()->SetTitle("momentum  (MeV/c)");
  return tg;
}




/// cg(int ng, int ns, double xs[], double ys[], std::string ymodel);
/// ----------------------------------------------------------------------------
/// cycle graph from array of steps and model to interpol y(x) in between steps
/// ng: nb of graph points
/// ns: nb of steps points
/// xs: xs steps array
/// ys: ys steps array
/// ymodel: model to interpolate y(x) when x_s < x_g < x_s+1
TGraph* cg(int np, int cns, double xsa[], double ysa[], string const& ymodel)
{
  /// sampling arrays
  double gxa[np];  /// graph x array
  double gya[np];  /// graph y array
  double xstep( (xsa[cns-1]-xsa[0])/(np-1) );
  
  /// init first point
  gxa[0] = xsa[0];
  gya[0] = ysa[0];
  
  /// main loops for the different model of interpolation
  /// ---------------------------------------------------
  
  /// linear model : linear interpolation
  if( ymodel=="linear" ) {   
    for( int i=1; i<np; ++i ) { /// loop over sampling points
      /// x array
      gxa[i] = xstep*i + xsa[0];
      /// y array
      for( int j=0; j<cns-1; ++j ) {  /// loop over cycle step points
        if( gxa[i]>xsa[j] && gxa[i]<=xsa[j+1] ) { 
          gya[i]=RMath::y2x(xsa[j],xsa[j+1],gxa[i],ysa[j],ysa[j+1],"linear");
        }
      }
    } 
  }/// end of linear case
    
  /// cooling model : linear at deceleration, cooling curve at cooling
  else if( ymodel=="cooling" ) {
    /// scaled cooling graph (betaCool sim)
    int const fcs(3); /// first cooling step
    int const scs(7); /// second cooling step
    TGraph* fcg = scool("35MeV_c", xsa[fcs], ysa[fcs], ysa[fcs+1]); /// elenaCool.h
    TGraph* scg = scool("100keV",  xsa[scs], ysa[scs], ysa[scs+1]); 
    for( int i=1; i<np; ++i ) { /// loop over sampling points
      /// x array
      gxa[i] = xstep*i + xsa[0];
      /// y array
      for( int j=0; j<cns-1; ++j ) {  /// loop over cycle step points
        if( gxa[i]>xsa[j] && gxa[i]<=xsa[j+1] && j==fcs ) { /// 1st cooling
          gya[i]=vcool(fcg, gxa[i]); /// function in ecool data folder
          //gya[i]=RMath::y2x(xsa[j],xsa[j+1],gxa[i],ysa[j],ysa[j+1],"linear");
        }
        else if( gxa[i]>xsa[j] && gxa[i]<=xsa[j+1] && j==scs ) { /// second cooling
          gya[i]=vcool(scg, gxa[i]); /// function in ecool data folder
          //gya[i]=RMath::y2x(xsa[j],xsa[j+1],gxa[i],ysa[j],ysa[j+1],"linear");
        }
        else if( gxa[i]>xsa[j] && gxa[i]<=xsa[j+1] ) { /// linear
          gya[i]=RMath::y2x(xsa[j],xsa[j+1],gxa[i],ysa[j],ysa[j+1],"linear");
        }
      }
    } 
    delete fcg; fcg=0;
    delete scg; scg=0;
  } /// end of cooling case
  
  /// output!
  TGraph* tg = new TGraph(np);
  for( int i=0; i<np; ++i ) { 
    tg->SetPoint(i, gxa[i], gya[i]);
  }
  /// end of routine
  return tg;
}

















///beam gas interactions functions
/// ****************************************************************************

/// 20170918 do not work anymore with new flyers and beams
/*
/// elenaRate() - 22/06/2017
/// ----------------------------------------------------------------------------
/// compute elena rate with H2
void elenaRate()
{
  /// function header
  /// --------------------------------------------------------------------------
  cout << endl;		cout << "elenaRate() - 22.06.2017" << endl;
  RDump::line();	cout << "elena rate with beam" << endl;
  
  
  /// flyer
  double beta(0.037);
  RFlyer* fly = new RFlyer((RParticle*)RParticleList::hpbar(), beta);
  //cout << endl; fly->dump();
  
  /// beam parameters
  double npb=2.e7;			/// nb of particles pre bunch
  RBeam* beam(0);
  beam = new RBeam("elenaBeam", npb, fly);
  cout << endl; cout << "dump beam " << endl; RDump::line(9, "-");
  beam->dumpLine();
  
  
  /// machine
  double A(75.*RUnits::um);				/// acceptance
  double tBeta(3.*RUnits::m);		/// twiss beta
  double len(30.4055*RUnits::m);	/// length
  RMachine* mac = new RMachine("elena", A, tBeta, len);
  //cout << endl; mac->dump();
  
  /// gas...
  double pre(4.e-12*RUnits::mbar);
  double tem(300.);
  RGas* gas(0);
  
  cout << endl;
  cout << endl; RDump::line(80, "o"); 
  cout << "molecular Nitrogen" << endl; RDump::line(80, "o"); 
  
  /// ... = pure N2
  gas = new RGas((RMolecule*)RParticleList::mN2(), pre, tem);
  cout << endl; gas->dumpRates(beam, mac); 
  /// change energy
  fly->ekin(100.*RUnits::keV);
  cout << endl; 
  cout << endl; gas->dumpRates(beam, mac);
  /// delete gas
  delete gas; gas=0;
  
  cout << endl;
  cout << endl; RDump::line(80, "o");
  cout << "molecular hydrogen" << endl; RDump::line(80, "o"); 
  
  /// ... = pure H2
  gas = new RGas((RMolecule*)RParticleList::mH2(), pre, tem);
  
  /// change momentum
  fly->momentum(100.*RUnits::MeV_c);
  cout << endl; 
  cout << endl; gas->dumpRates(beam, mac); 
  /// change momentum
  fly->momentum(35.*RUnits::MeV_c);
  cout << endl; 
  cout << endl; gas->dumpRates(beam, mac); 
  /// change energy
  fly->ekin(100.*RUnits::keV);
  cout << endl; 
  cout << endl; gas->dumpRates(beam, mac);
  /// delete gas
  delete gas; gas=0;
  
  
  /// end of function
  /// --------------------------------------------------------------------------
  cout << endl; cout << "elenaRate() function end " << endl; RDump::line();
}
*/



/// testRutherfordScattering() - 22.06.2017
/// ----------------------------------------------------------------------------
/// test new atom/mol/gas "old style" for rutherford scattering results
void testRutherfordScattering()
{
  /// function header
  /// --------------------------------------------------------------------------
  cout << endl;cout << "testRutherfordScattering() - 22.06.2017" << endl; 
  RDump::line(); 
  cout << "test atom/mol/gas for rutherford scattering results" << endl;
  
  
  /// do stuff
  /// --------------------------------------------------------------------------
  /// define a gas
  double pre(4.e-12*RUnits::mbar);  /// overall pressure
  double tem(300.*RUnits::K);       /// overall temperature
  vector<RMolecule*> mov;        		/// molecules
  vector<double> mfv;               /// composition fract of molecules
  //mov.push_back((RMolecule*)RParticleList::mH2());  mfv.push_back(0.50);
  //mov.push_back((RMolecule*)RParticleList::mN2());  mfv.push_back(0.50);
  mov.push_back((RMolecule*)RParticleList::mN2());  mfv.push_back(1.00);
  RGas* gas = new RGas("elena test gas", pre, tem);
  gas->setCompound(mov, mfv);
  /// dump the gas
  // changed 20170929 cout << endl; gas->dump(cout, "pkr261");
  cout << endl; gas->dump(cout);
  /// define machine
  double len(30.4055*RUnits::m);	/// machine length
  double tbp(3.*RUnits::m);				/// twiss beta parameter at detector
  double tap(75.*RUnits::um);			/// twiss acceptance param = tbpMax/tbp
  RMachine* mac(0);
  mac = new RMachine("elena", tap, tbp, len);
  /// dump the machine
  cout << endl; 
  cout << endl; cout << "dump machine parameters" << endl; RDump::line(23, "-");
  mac->dumpLine();

  /// define flyers
  double beta(0.037);
  RFlyer* fly = new RFlyer((RParticle*)RParticleList::hpbar(), beta);
  /// dump the flyer
  cout << endl; 
  cout << endl; cout << "dump flyer parameters" << endl; RDump::line(21, "-");
  fly->dumpLine(); fly->dumpBeta();
  
  /// atom dumpcs test
  RAtom* ato = gas->getCompound()[0]->getComposition()[0];
  cout << endl; 
  cout << endl; ato->dumpcs(fly, mac);
  
  /// molecule dumpcs test
  cout << endl; 
  for( int i=0; i<gas->getCompound().size(); ++i ) {
    cout << endl; gas->getCompound()[i]->dumpcs(fly, mac);
  }
  
  /// dump rates
  cout << endl; 
  cout << endl; gas->dumpRates(fly, mac);
  delete gas; gas = new RGas((RMolecule*)RParticleList::mH2(), pre, tem);
  cout << endl; 
  cout << endl; gas->dumpRates(fly, mac);
  
  
  /// end of function
  /// --------------------------------------------------------------------------
  delete gas; gas=0; delete fly; fly=0; delete mac; mac=0;
  cout << endl; cout << "testRutherfordScattering() end " << endl; 
  RDump::line();
}







/// testMacFly() - 15.06.2017
/// ----------------------------------------------------------------------------
/// test Machine/Flyer class
void testMacFly()
{
  /// function header
  cout << endl; cout << "testMacFly() - 15.06.2017" << endl; RDump::line(); 
  cout << "test Machine/Flyer class" << endl;
  
  
  /// middle
  RParticle* par(0);
  RFlyer* fly(0);
  double beta(0);
  double kin(0);
  double mv(0);
  par = (RParticle*)RParticleList::hpbar();
  fly = new RFlyer(par, beta); delete par; par=0;
  
  cout << "the flyer is: ";
  fly->dumpLine();
  cout << endl; cout << "beta=" << beta << " gives: "<< endl;
  fly->dumpBeta();
  
 	kin=5.3*RUnits::MeV; fly->ekin(kin);
  cout << endl; cout << "ekin=" << kin/RUnits::MeV << "MeV" << " gives: ";
  cout << endl; fly->dumpBeta();
  
  mv=35.*RUnits::MeV_c; fly->momentum(mv);
 	cout << endl; cout << "mv=" << mv/RUnits::MeV_c << "MeV/c" << " gives: ";
 	cout << endl; fly->dumpBeta();
  
  kin=100.*RUnits::keV; fly->ekin(kin);
  cout << endl; cout << "ekin=" << kin/RUnits::MeV << "MeV" << " gives: ";
  cout << endl; fly->dumpBeta();
 
  
  /// end of function
  delete fly; fly=0;
  cout << endl; cout << "tMac() function end " << endl; RDump::line(); 
}






/// testGas() - 15.06.2017
/// ----------------------------------------------------------------------------
/// test Gas class, especially dump and dumpRates
void testGas()
{
  /// function header
  cout << endl; cout << "testGas() - 15.06.2017" << endl; RDump::line(); 
  cout << "test Gas class function start" << endl;
  
  
  /// gas properties
  double pre(4.e-12*RUnits::mbar);  /// overall pressure
  double tem(298.*RUnits::K);       /// overall temperature
  RGas* gas(0);											/// gas pointer
  /// define a gas
  vector<RMolecule*> mov(0);        /// molecules vector
  vector<double> mfv;               /// composition fract of molecules
  mov.push_back((RMolecule*)RParticleList::mH2());  mfv.push_back(0.91);
  mov.push_back((RMolecule*)RParticleList::mCO2()); mfv.push_back(0.05);
  mov.push_back((RMolecule*)RParticleList::mCO());  mfv.push_back(0.03);
  mov.push_back((RMolecule*)RParticleList::mCH4()); mfv.push_back(0.01);
  gas = new RGas("elenaVac", pre, tem);
  gas->setCompound(mov, mfv);
  
  /// dump the gas
	cout << endl << "test elena gas with molecules mixing" << endl;
  RDump::line(36, "-");
  // 20170929cout << endl; gas->dump(cout, "pkr261"); delete gas; gas = 0;
  cout << endl; gas->dump(cout); delete gas; gas = 0;
  /// simple H2 gas
  gas = new RGas((RMolecule*)RParticleList::mH2(), pre, tem);
  cout << endl << "test simple H2 gas " << endl;
  RDump::line(18, "-");
  // 20170929cout << endl; gas->dump(cout, "pkr261"); delete gas; gas = 0;
  cout << endl; gas->dump(cout); delete gas; gas = 0;
  
  /// end of function
  cout << endl; cout << "testGas() end" << endl;  RDump::line(); 
}








/// testval_csion_pbarGas() - 15.06.2017
/// ----------------------------------------------------------------------------
/// test Ionisation cross section eval/interpolation function
void testval_csion_pbarGas()
{
  /// function header
  cout << endl; cout << "testval_csion_pbarGas() - 15.06.2017" << endl; 
  RDump::line(); cout << "check pbar ionisation cs evaluation" << endl;
  
  /// do stuff 
  cout << (RParticleList::mCO2())->ion_pbarcs(500.*RUnits::keV) 
       << "m2" << endl;
  
  /// end of function
  cout << endl; cout << "testval_csion_pbarGas()" << endl; RDump::line(); 
}







/// plotIon_pbarcs() - 14.06.2017
/// ----------------------------------------------------------------------------
/// plot ion_pbarcs data graphs
void plotIon_pbarcs(string const& today)
{
  /// function header
  cout << endl; cout << "plotIon_pbarcs(string today) - 14.06.2017" << endl; 
  RDump::line(); cout << "plot ion_pbarcs data graphs" << endl;
  
  
  /// do stuff
  /// string manip
  string pPrefix("../results/");
  string path(pPrefix+today+"/");
  string fPrefix("cs_ion_pbarData_");
  string file("");
  string fSuffix(".pdf");
  /// molecule and graph
  const RMolecule* mol(0);
  TGraphErrors* tge(0);
  
  /// molecular argon
  mol = RParticleList::mAr(); tge = mol->ion_pbarcs_data(); tge->Draw(); 
  file = "mAr.pdf"; gPad->Print((path+fPrefix+file+fSuffix).c_str(), "pdf");
  delete tge; delete mol; delete gPad;
  
  /// molecular ch4
  mol = RParticleList::mCH4();  tge = mol->ion_pbarcs_data(); tge->Draw(); 
  file = "mCH4.pdf"; gPad->Print((path+fPrefix+file+fSuffix).c_str(), "pdf");
  delete tge; delete mol; delete gPad;
  
  /// molecular CO
  mol = RParticleList::mCO();   tge = mol->ion_pbarcs_data(); tge->Draw(); 
  delete tge; delete mol; delete gPad;
  
  /// molecular CO2
  mol = RParticleList::mCO2();  tge = mol->ion_pbarcs_data(); tge->Draw(); 
  file = "mCO2.pdf"; gPad->Print((path+fPrefix+file+fSuffix).c_str(), "pdf");
  delete tge; delete mol; delete gPad;
  
  /// molecular H2
  mol = RParticleList::mH2();   tge = mol->ion_pbarcs_data(); tge->Draw(); 
  file = "mH2.pdf"; gPad->Print((path+fPrefix+file+fSuffix).c_str(), "pdf");
  delete tge; delete mol; delete gPad;
  
  /// molecular krypton
  mol = RParticleList::mKr();   tge = mol->ion_pbarcs_data(); tge->Draw(); 
  file = "mKr.pdf"; gPad->Print((path+fPrefix+file+fSuffix).c_str(), "pdf");
  delete tge; delete mol; delete gPad;
  
  /// molecular nitrogen
  mol = RParticleList::mN2();   tge = mol->ion_pbarcs_data(); tge->Draw(); 
  file = "mN2.pdf"; gPad->Print((path+fPrefix+file+fSuffix).c_str(), "pdf");
  delete tge; delete mol; delete gPad;
  
  /// molecular neon
  mol = RParticleList::mNe();   tge = mol->ion_pbarcs_data(); tge->Draw(); 
  file = "mNe.pdf"; gPad->Print((path+fPrefix+file+fSuffix).c_str(), "pdf");
  delete tge; delete mol; delete gPad;
  
  /// molecular xenon
  mol = RParticleList::mXe();   tge = mol->ion_pbarcs_data(); tge->Draw(); 
  file = "mXe.pdf"; gPad->Print((path+fPrefix+file+fSuffix).c_str(), "pdf");
  delete tge; delete mol; delete gPad;
  
  
  /// end of function
  tge = 0; mol = 0; gPad = 0;
  cout << endl; cout << "plotIon_pbarcs() end" << endl; RDump::line();
}






/// testMolecule() - 14.06.2017
/// ----------------------------------------------------------------------------
/// test RMolecule class
void testMolecule()
{
  /// function header
  cout << endl << "testMolecule() - 14.06.2017" << endl; RDump::line();
  cout << "test Molecule class function start" << endl;
  
  
  /// do stuff
  RMolecule* mol(0);
  
  /// syntax 1
	vector<RAtom*> mcv;	/// molecule composition vector
  mcv.push_back((RAtom*)RParticleList::aC12());
  mcv.push_back((RAtom*)RParticleList::aO16());
  mcv.push_back((RAtom*)RParticleList::aO16());
  mol = new RMolecule("Carbon dioxide", "CO2", "C12-O16_2", mcv);
	cout << endl; cout << "syntax 1:" << endl;
  mol->dump(); delete mol; mol = 0;
  
  /// syntax 2
  mol = new RMolecule("repMol", "rM", "r30M0");
  mol->addAtom((RAtom*)RParticleList::aN14());
  mol->addAtom((RAtom*)RParticleList::aH1());
  mol->addAtom((RAtom*)RParticleList::aH1());
  mol->addAtom((RAtom*)RParticleList::aH1());
  cout << endl; cout << "syntax 2:" << endl;
  mol->dump(); delete mol; mol = 0;
  
  /// syntax 3
  mol = (RMolecule*)RParticleList::mCH4();
 	cout << endl; cout << "syntax 3:" << endl;
  mol->dump(); delete mol; mol = 0;
  
  
  /// end of function
  cout << endl; cout << "testMolecule() end" << endl; RDump::line();
}





/// listParticles() - 07.07.2017
/// ----------------------------------------------------------------------------
/// list all particles in particleList namespace
void listParticles()
{
  /// function header
  cout << endl << "listParticles() - 07.07.2017" << endl; RDump::line();
  cout << "list particles in the namespace RParticleList.cc" << endl;
  
  const RAtom* atom(0);
  cout << endl; 
  cout << "dump atoms (line)" << endl;
  RDump::line(21, "-");
  atom = RParticleList::aH1();     atom->dumpLine(); delete atom; 
  atom = RParticleList::aC12();    atom->dumpLine(); delete atom;
  atom = RParticleList::aN14();    atom->dumpLine(); delete atom;
  atom = RParticleList::aO16();    atom->dumpLine(); delete atom;
  atom = RParticleList::aNe20();   atom->dumpLine(); delete atom;
  atom = RParticleList::aAr40();   atom->dumpLine(); delete atom;
  atom = RParticleList::aKr84();   atom->dumpLine(); delete atom;
  atom = RParticleList::aXe132();  atom->dumpLine(); delete atom;
  atom = 0;
  
  const RMolecule* mol(0);
  cout << endl << endl;
  cout << "dump molecules" << endl;
  RDump::line(14, "-");
  mol = RParticleList::mAr();   mol->dumpLine(); cout << endl; delete mol;
  mol = RParticleList::mCH4();  mol->dumpLine(); cout << endl; delete mol;
  mol = RParticleList::mCO();   mol->dumpLine(); cout << endl; delete mol;
  mol = RParticleList::mCO2();  mol->dumpLine(); cout << endl; delete mol;
  mol = RParticleList::mH2();   mol->dumpLine(); cout << endl; delete mol;
  mol = RParticleList::mKr();   mol->dumpLine(); cout << endl; delete mol;
  mol = RParticleList::mN2();   mol->dumpLine(); cout << endl; delete mol;
  mol = RParticleList::mNe();   mol->dumpLine(); cout << endl; delete mol;
  mol = RParticleList::mXe();   mol->dumpLine(); cout << endl; delete mol;
  mol = 0;
  
  const RParticle* par(0);
  cout << endl << endl;
  cout << "dump hadrons" << endl;
  RDump::line(12, "-");
  par = RParticleList::hpbar();   par->dumpLine(); cout << endl; delete par;
  par = 0;
  
  
  /// end of function
  cout << endl; cout << "listParticles() function end" << endl; RDump::line();
}






/// testAtom() - 13.06.2017
/// ----------------------------------------------------------------------------
/// test RAtom class
void testAtom()
{
  /// function header
  cout << "testAtom() - 13.06.2017" << endl; RDump::line();
  cout << "test Atom class function start" << endl;
  
  
  /// test class, basics
  RAtom* ato1 = new RAtom("repAtom", "Re", 100., 200.);
  cout << endl; ato1->dump();
  
  cout << endl << "test atom get accessors" << endl;
  RDump::line(23, "-");
  cout << ato1->getName() << endl;
  cout << ato1->getSymbol() << endl;
  cout << ato1->getMass() << endl;
  cout << ato1->getCharge() << endl;
  cout << ato1->getZ() << endl;
  cout << ato1->getA() << endl;
  
  /// test polymorphism
  cout << endl << "test atom polymorphism (<< operator)" << endl;
  testPol(ato1);
  
  /// test total/geometric cross section
  cout << endl << "test atom geometric cross section" << endl;
  delete ato1; ato1=0;
  ato1 = new RAtom("Nitrogen", "N", 7, 14, 14.007*RUnits::amu);
  cout << *ato1 << endl;
  
  
  /// end of function
  cout << endl; cout << "tAtom() function end" << endl; RDump::line();
}






/// testPol - 12.06.2017
/// ----------------------------------------------------------------------------
/// test polymorphisme
void testPol(RObject* const obj)
{
  cout << *obj << endl;
}





/// testParticle() - 12.06.2017
/// ----------------------------------------------------------------------------
/// test class RParticle
void testParticle()
{
  /// function header
  cout << "testParticle() - 12.06.2017" << endl; RDump::line();
  cout << "test RParticle class function start" << endl;
  
  cout << endl;
  cout << "RObject::nb     " << RObject::nb << endl;
  cout << "RParticle::nb() " << RParticle::nb() << endl;
  
  cout << endl;
  RParticle::list();
  
  RObject* obj1 = new RObject("obj1"); 
  RParticle* par1(0); 
  RParticle* par2(0); 
  RParticle* par3(0); 
  par1 = new RParticle("proton", "p", 1.*RUnits::eV, 938.*RUnits::MeV_c2); 
  par2 = new RParticle("neutron", "n", 0., 940.*RUnits::MeV_c2); 
  par3 = new RParticle("Helium 2+", "He2+", 2.*RUnits::eV, 4.*RUnits::MeV_c2); 
  
  cout << endl; testPol(obj1);
  cout << endl; testPol(par1);
  cout << endl; testPol(par2);
  cout << endl; testPol(par3);
  
  cout << endl;
  cout << "RObject::nb     " << RObject::nb << endl;
  cout << "RParticle::nb() " << RParticle::nb() << endl;
  cout << endl;
  RParticle::list();
  
  delete par3; par3=0;
  
  cout << endl;
  cout << "RObject::nb     " << RObject::nb << endl;
  cout << "RParticle::nb() " << RParticle::nb() << endl;
  cout << endl;
  RParticle::list();
  
  
  
  /// end of function
  delete par1; par1=0;
  delete par2; par2=0;
  cout << endl; cout << "testParticle() function end" << endl; RDump::line();
}




/// testObject() - 09.06.2017
/// ----------------------------------------------------------------------------
/// test class RObject
void testObject()
{
  /// function header
  cout << "testObject() - 09.06.2017" << endl; RDump::line();
  cout << "test RObject class function start" << endl;
  
  /// declare pointers to new objects
  RObject* obj1 = new RObject();
  RObject* obj2 = new RObject();
  
  /// do stuff
  obj1->setName("name1");
  obj2->setName("name2");
  RObject obj3 = *obj1; /// call copy constructor
  RObject obj4 = obj3;  /// call copy constructor
  cout << obj1->getName() << "  " << obj2->getName() << "  " 
       << obj3.getName()  << "  " << obj4.getName()  << "  " 
       << endl;
  
  obj4 = *obj2; /// call "=" operator
  cout << obj1->getName() << "  " << obj2->getName() << "  " 
       << obj3.getName()  << "  " << obj4.getName()  << "  " 
       << endl;
  
  /// release memory
  delete obj1; obj1=0;
  delete obj2; obj2=0;
  
  bool infTest = obj3<obj4; /// test "<" operator
  cout << infTest << endl;
  cout << (obj3<obj4) << endl;
  
  /// end of function
  cout << endl; cout << "testObject() function end" << endl; RDump::line();
}




