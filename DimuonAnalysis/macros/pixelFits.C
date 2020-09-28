


// Gaussian function
Double_t gaussian(Double_t *x, Double_t *par) {
	return par[0]*( exp(-0.5*( (x[0]-par[1])/par[2])*((x[0]-par[1])/par[2]) ) );
}

// Exponential function
Double_t exponential(Double_t *x, Double_t *par) {
	return exp(par[0]+par[1]*x[0]);
}

// Landau function
Double_t landau(Double_t *x, Double_t *par) {
	return 1.71534e+05*TMath::Landau(-x[0],2.10645e+00, 2.22117e-01);
	return par[0]*TMath::Landau(-x[0],par[1],par[2]);
}


// Quadratic background function
Double_t parabola(Double_t *x, Double_t *par) {
	return par[0] + par[1]*x[0] + par[2]*x[0]*x[0];
}

// Lorentzian Peak function
Double_t lorpeak(Double_t *x, Double_t *par) {
	return (0.5*par[0]*par[1]/TMath::Pi()) / TMath::Max(1.e-10,
	(x[0]-par[2])*(x[0]-par[2])+ .25*par[1]*par[1]);
}

// Sum of background and peak function
Double_t fitFunction(Double_t *x, Double_t *par) {

	// return landau(x,par);
	return exponential(x,par) + gaussian(x,par);

	// 
	if (*x<11.8 && *x>10){
		return landau(x,par) + gaussian(x,par);
	}
	else if (*x<5.9 && *x>7.6){
		return landau(x,par) + gaussian(x,par);
	}
	else if (*x<15 && *x>17){
		return landau(x,par) + gaussian(x,par);
	}
	else{
		return landau(x,par);
	}
}


void myfitter(){
	TFile *f = new TFile("custom_hists2018.root");
	// TFile *f = new TFile("hists2018C_radius.root");

	TH1F * h = (TH1F*)f->Get("dimuonR");

	// Layer 1
	TF1* fitFcn = new TF1("fitFcn","exp([0]+[1]*x)+[2]*( exp(-0.5*( (x-[3])/[4])*((x-[3])/[4]) ) )",1.7,3.5);
	fitFcn->SetParameters(0.9,-4.7,1e3,2,0.4);
	h->Fit(fitFcn,"R");
	// Layer 2
	TF1* fitFcn = new TF1("fitFcn","exp([0]+[1]*x)+[2]*( exp(-0.5*( (x-[3])/[4])*((x-[3])/[4]) ) )",5,8);
	fitFcn->SetParameters(0.9,-4.7,1e3,6.4,0.4);
	h->Fit(fitFcn,"R");
	// Layer 3
	TF1* fitFcn = new TF1("fitFcn","exp([0]+[1]*x)+[2]*( exp(-0.5*( (x-[3])/[4])*((x-[3])/[4]) ) )",9,14);
	fitFcn->SetParameters(5.37476e+00,-9.09412e-02,200,11,0.5);
	h->Fit(fitFcn,"R");
	// Layer 4
	TF1* fitFcn = new TF1("fitFcn","exp([0]+[1]*x)+[2]*( exp(-0.5*( (x-[3])/[4])*((x-[3])/[4]) ) )",13,18);
	fitFcn->SetParameters(5.37476e+00,-9.09412e-02,100,16,0.5);
	h->Fit(fitFcn,"R");


	// TF1* trial_exp = new TF1("trial_exp","exp([0]+[1]*x)",0,5);
	// TF1* trial_gau = new TF1("trial_gau","[0]*(    exp(-0.5*( (x-[1])/[2])*((x-[1])/[2]) )       )",0,5);
	// trial_exp->SetParameters(1.81351e+01,-4.74637e+00);
	// trial_gau->SetParameters(3.14799e+04,1.63496e+00,8.59824e-01);

	TCanvas *c1 = new TCanvas("c1","Canvas Example",200,10,600,480);
	h->Draw();
	// trial_exp->Draw("same");
	// trial_gau->Draw("same");
	c1->SetLogy();
	c1->Update();
	c1->SaveAs("dummy.png");
	gApplication->Terminate();

}