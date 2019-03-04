// energy unit is MeV
using namespace TMath;

static double fTheta( double r , double E)
{
  double a = 1.1/115.0E3;  
  double c1 =1 + a*E;
  double c2 =1 -a*E;
  return (c2*pow(c1/c2,r)-1)/(a*E);
}

static double Acc( double r , double E)
{
  double a = 1.1/115.0E3;  
  double c1 =1 + a*E;
  double c2 =1 -a*E;
  return Log( (1.+a*E*r)/c2 )/Log(c1/c2);
}



void cosmicTest(Int_t max=100000) {

  TRandom *rand = new TRandom2();
  TFile *outfile =  new TFile("cosmicTest.root","recreate");
  TNtuple *nt= new TNtuple("test","test","r:r2:e:cos:theta:a");

  double A = 0.14;
  double b = 2.7;
  
  for(Int_t j=0; j< max; ++j)  {
    if(j%1000==0) printf("... %i \n",j);
    double r = rand->Rndm();
    double E = 1000.0*(-1/b)*Log( 1-b*r/A);
    for(Int_t i=0; i< 100; ++i)  {
      double r2 = rand->Rndm();
      //double r3 = 2.*(rand->Rndm()-0.5);
      Double_t c = fTheta(r2,E);
      double theta = ACos(c)*180/Pi();
      nt->Fill(r,r2,E,c,theta,Acc(r2,E));
    }
  }
  outfile->Write();

}


