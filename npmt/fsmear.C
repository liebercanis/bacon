using namespace TMath;

// Macro myfunc.C
Double_t ffunction(Double_t *xx, Double_t *par)
{
   Float_t x =xx[0];
   double tau = par[0];
   double sigma = par[1];
   double z = (-x/sigma+sigma/tau)/Sqrt(2);
   double a = pow(sigma/tau,2.)/2.0;
   Double_t f = 1./2./tau*Exp(-x/tau)*Exp(a)*Erfc(z);
   return f;
}
void fsmear() 
{
   TF1 *f1 = new TF1("ffunc",ffunction,-20,100,2);
   f1->SetParameters(.5,2);
   f1->SetParNames("tau","sigma");
   gStyle->SetOptFit();
   gStyle->SetOptStat();
   TCanvas *can = new TCanvas("fsmear","fsmear");
   //can->SetLogy();
   f1->Draw();
   printf(" %f \n",f1->Integral(-20,100));
}
