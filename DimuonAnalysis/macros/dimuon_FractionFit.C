

void dimuon_FractionFit(){
    
    // OPEN THE FILE
    TFile *f = new TFile("hists2018C_all2.root");

    // GET ALL HISTOGRAMS
    TH1F *data = (TH1F*)f->Get("dimuonL1_mass1");
    data->Add( (TH1F*)f->Get("dimuonL234_mass1") );
    data->Add( (TH1F*)f->Get("dimuonNotInLayers_mass1") );
    data->SetMarkerStyle(20);
    data->SetMarkerSize(.5);
    data->SetLineColor(2);
    data->SetTitle("Data to fit (all R-regions);m_{#mu#mu} [GeV];Events");

    TH1F *bkg1 = (TH1F*)f->Get("dimuonL234_mass1");
    bkg1->SetLineColor(3);
    bkg1->SetTitle("Template for conversions (Pixel L2-3-4);m_{#mu#mu} [GeV];Events");

    TH1F *bkg2 = (TH1F*)f->Get("dimuonNotInLayers_mass1");
    bkg2->SetLineColor(2);
    bkg2->SetTitle("Template for fakes / heavy flavor decays;m_{#mu#mu} [GeV];Events");

    TH1F *sig = (TH1F*)f->Get("dimuonL1_mass1");
    sig->SetLineColor(1);
    sig->SetTitle("Signal enriched (Pixel L1);m_{#mu#mu} [GeV];Events");


    // REBIN IF NEEDED 
    int rebin_factor=40;
    data->Rebin(rebin_factor);
    bkg1->Rebin(rebin_factor);
    bkg2->Rebin(rebin_factor);

    // DEFINE THE TEMPLATES TO USE IN THE FIT
    TObjArray *templates = new TObjArray(2);
    templates->Add(bkg1);
    templates->Add(bkg2);

    // PERFORM THE TEMPLATE FIT
    TFractionFitter* fit = new TFractionFitter(data, templates); // initialise
    fit->Constrain(0,0.0,1.0);               // constrain fraction 1 to be between 0 and 1
    fit->Constrain(1,0.0,1.0);               // constrain fraction 1 to be between 0 and 1
    Int_t status = fit->Fit();               // perform the fit
    cout << "fit status: " << status << endl;


    // DISPLAY
    gStyle->SetOptStat(0);
    TCanvas c1("c1", "FractionFitter example", 700, 700);

    c1.cd();
    gPad->SetLogy();
    TH1F *bkg1_1=(TH1F*)bkg1->Clone();
    bkg1_1->Scale(1/(bkg1->Integral()));
    bkg1_1->Draw();
    TH1F *bkg2_1=(TH1F*)bkg2->Clone();
    bkg2_1->Scale(1/(bkg2_1->Integral()));
    bkg2_1->Draw("same");
    // TH1F *sig_1=(TH1F*)sig->Clone();
    // sig_1->Scale(1/(sig_1->Integral()));
    // sig_1->Draw("same");
    gPad->BuildLegend();
    bkg1_1->SetTitle("Shape templates (normalized);m_{#mu#mu} [GeV];Events (norm.)");


    TCanvas c2("c2", "FractionFitter example", 700, 700);
    c2.cd();
    gPad->SetLogy();
    Double_t p0, p1, errP0, errP1;
    TLatex l;
    l.SetTextSize(.035);
    Char_t texte[200];

    if (/*status == 0*/1) {                       // check on fit status

        TH1F* result = (TH1F*) fit->GetPlot();
        result->SetLineColor(4);
        std::cout<<result->Integral()<<std::endl;
        std::cout<<data->Integral()<<std::endl;
        fit->GetResult( 0, p0, errP0); // p0
        fit->GetResult( 1, p1, errP1);
        
        data->GetYaxis()->SetRangeUser(1e1,1e6);
        // data->DrawClone("Ep");
        result->SetTitle("Data distribution with fitted templates");

        auto rp = new TRatioPlot(result,data);
        rp->Draw();

    }
    c2.Update();

    c1.SaveAs("/eos/user/p/pmeiring/www/CMSvDAS/DimuonScouting/TrueMuonium/dimuon_FitTemplates.png");
    c1.SaveAs("/eos/user/p/pmeiring/www/CMSvDAS/DimuonScouting/TrueMuonium/dimuon_FitTemplates.pdf");
    // c2.SaveAs("/eos/user/p/pmeiring/www/CMSvDAS/DimuonScouting/TrueMuonium/dimuon_FractionFit_L1only.png");
    // c2.SaveAs("/eos/user/p/pmeiring/www/CMSvDAS/DimuonScouting/TrueMuonium/dimuon_FractionFit_L1only.pdf");   
    c2.SaveAs("/eos/user/p/pmeiring/www/CMSvDAS/DimuonScouting/TrueMuonium/dimuon_FractionFit.png");
    c2.SaveAs("/eos/user/p/pmeiring/www/CMSvDAS/DimuonScouting/TrueMuonium/dimuon_FractionFit.pdf");   
    gApplication->Terminate();
}