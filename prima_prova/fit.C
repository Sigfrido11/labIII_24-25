void fit()
{

  //gStyle->SetOptFit(1111);

  TGraphErrors *silicio = new TGraphErrors("silicio_all.txt", "%lg %lg %lg %lg");

  TGraphErrors *germanio = new TGraphErrors("germanio.txt", "%lg %lg %lg %lg");

  TGraphErrors *calibrazione =
      new TGraphErrors("calibrazione.txt", "%lg %lg %lg %lg");

  TF1 *fitLin = new TF1("fitLinSil", "[0] + [1] * x", 0, 1);
  TF1 *fitExpSil = new TF1("fitExpSil", "[0] * exp(x/[1]) + [2]", 0, 0.51);
  TF1 *fitExpGer = new TF1("fitExpGer", "[0] * exp(x/[1])+ [2]", 0, 0.132);

  fitLin->SetParName(0, "intercetta");
  fitLin->SetParName(1, "pendenza");
  fitLin->SetParameter(0, 1);
  fitLin->SetParameter(1, 1);
  // fitLin->SetParLimits(0, 0, 1e5);
  // fitLin->SetParLimits(1, 4, 8);

  fitExpSil->SetParName(0, "I0 sil");
  fitExpSil->SetParName(1, "#eta-vt sil");
  fitExpSil->SetParameter(0, 1e-6);
  fitExpSil->SetParameter(1, 2. * 300. / 11600.);
  fitExpSil->SetParameter(2,0.02);
  // fitExpSil->SetParLimits(0, 0, 1e5);
  // fitExpSil->SetParLimits(1, 4, 8);

  fitExpGer->SetParName(0, "I0 ger");
  fitExpGer->SetParName(1, "#eta-vt ger");
  fitExpGer->SetParameter(0, 1e-3);
  fitExpGer->SetParameter(1, 300. / 11600.);
  // fitExpGer->SetParLimits(0, 0, 1e5);
  // fitExpGer->SetParLimits(1, 4, 8);

  calibrazione->Fit(fitLin, "R");
  silicio->Fit(fitExpSil, "R");
  germanio->Fit(fitExpGer, "R");

  TCanvas *c1 = new TCanvas("Silicio", "Silicio", 200, 10, 600, 400);
  c1->SetLogy();

  TCanvas *c2 = new TCanvas("Germanio", "Germanio", 200, 10, 600, 400);
  c2->SetLogy();

  TCanvas *c3 = new TCanvas("Silicio Cal", "Silicio  Cal", 200, 10, 600, 400);

  calibrazione->SetLineColor(1);
  calibrazione->SetTitle("Calibrazione O-scope DMM silicio");
  calibrazione->GetYaxis()->SetTitleOffset(1.2);
  calibrazione->GetXaxis()->SetTitleSize(0.04);
  calibrazione->GetYaxis()->SetTitleSize(0.04);
  calibrazione->GetXaxis()->SetTitle("DMM V (V)");
  calibrazione->GetYaxis()->SetTitle("O-scope V (V)");
  calibrazione->GetXaxis()->CenterTitle(true);
  calibrazione->GetXaxis()->CenterTitle(true);
  calibrazione->SetMarkerStyle(21); // Stile 21 è un cerchio
  calibrazione->SetMarkerSize(1.);  // Cambia la dimensione del marker
  calibrazione->GetXaxis()->SetRangeUser(0, 1);
  calibrazione->GetYaxis()->SetRangeUser(0, 1);

  silicio->SetLineColor(1);
  silicio->SetTitle("Caratteristica I-V silicio");
  silicio->GetYaxis()->SetTitleOffset(1.2);
  silicio->GetXaxis()->SetTitleSize(0.04);
  silicio->GetYaxis()->SetTitleSize(0.04);
  silicio->GetXaxis()->SetTitle("V (V)");
  silicio->GetYaxis()->SetTitle("I (mA)");
  silicio->GetXaxis()->CenterTitle(true);
  silicio->GetXaxis()->CenterTitle(true);
  silicio->SetMarkerStyle(21); // Stile 21 è un cerchio
  silicio->SetMarkerSize(1.);  // Cambia la dimensione del marker
  silicio->GetXaxis()->SetRangeUser(0, 0.8);

  germanio->SetLineColor(1);
  germanio->SetTitle("Caratteristica I-V germanio");
  germanio->GetYaxis()->SetTitleOffset(1.2);
  germanio->GetXaxis()->SetTitleSize(0.04);
  germanio->GetYaxis()->SetTitleSize(0.04);
  germanio->GetXaxis()->SetTitle("V (V)");
  germanio->GetYaxis()->SetTitle("I (mA)");
  germanio->GetXaxis()->CenterTitle(true);
  germanio->GetXaxis()->CenterTitle(true);
  germanio->SetMarkerStyle(21); // Stile 21 è un cerchio
  germanio->SetMarkerSize(1.);  // Cambia la dimensione del marker
  germanio->GetXaxis()->SetRangeUser(0, 0.4);

  c1->cd();
  silicio->Draw("APE");

  c2->cd();
  germanio->Draw("APE");

  c3->cd();
  calibrazione->Draw("APE");

  // finalGraph->GetYaxis()->SetLimits(1e-7, 10);
  // finalGraph->GetYaxis()->SetRangeUser(1e-7, 10);
  // auto leg1 = new TLegend(0.7, 0.1, 0.9, 0.3);
  // leg1->AddEntry(finalGraph, "Temperature 156 Mev", "");
  // leg1->AddEntry(finalGraph, "Radius 5 fm", "");
  // leg1->AddEntry(finalGraph, "|y|<0.5", "");
  // leg1->SetTextSize(0.04);
  // leg1->Draw();
}
