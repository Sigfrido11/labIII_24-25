

void sel_fit(const char *name_100, const char *name_200) {
  // Configurazione opzionale per il fit
  // gStyle->SetOptFit(1111);

  // Caricamento dei dati con TGraphErrors
  TGraphErrors *I_200 = new TGraphErrors(name_200, "%lg %lg %lg %lg");
  TGraphErrors *I_100 = new TGraphErrors(name_100, "%lg %lg %lg %lg");
  // Definizione dei fit lineari
  TF1 *fitLin_200 = new TF1("fitLin_200", "[0] + [1] * x", 0.8, 4);
  TF1 *fitLin_100 = new TF1("fitLin_100", "[0] + [1] * x", 0.5, 4);

  // Imposta parametri iniziali
  fitLin_200->SetParName(0, "intercetta 200");
  fitLin_200->SetParName(1, "pendenza 200");
  fitLin_200->SetParameter(0, 1);
  fitLin_200->SetParameter(1, 1);

  fitLin_100->SetParName(0, "intercetta 100");
  fitLin_100->SetParName(1, "pendenza 100");
  fitLin_100->SetParameter(0, 1);
  fitLin_100->SetParameter(1, 1);

  // covarianza
  TFitResultPtr res_100 = I_100->Fit("fitLin_100", "SR");

  const TMatrixDSym &cov_100 = res_100->GetCovarianceMatrix();

  double const cov_10_100 = cov_100(0, 1);

  // covarianza
  TFitResultPtr res_200 = I_200->Fit("fitLin_200", "SR");

  const TMatrixDSym &cov_200 = res_200->GetCovarianceMatrix();

  double const cov_10_200 = cov_200(0, 1);

  // Impostazione grafica per I_200
  I_200->SetLineColor(1);
  I_200->SetMarkerStyle(21);
  I_200->SetMarkerSize(1.);
  I_200->SetTitle("I_200");
  I_200->GetXaxis()->SetTitle("- V (bc) (V)");
  I_200->GetYaxis()->SetTitle("I (mA)");

  // Impostazione grafica per I_100
  I_100->SetLineColor(2);
  I_100->SetMarkerStyle(21);
  I_100->SetMarkerSize(1.);
  I_100->SetTitle("I_100");
  I_100->GetXaxis()->SetTitle("- V (bc) (V)");
  I_100->GetYaxis()->SetTitle("I (mA)");

  // Crea un multigraph
  TMultiGraph *mg = new TMultiGraph;
  mg->Add(I_200);
  mg->Add(I_100);
  mg->SetTitle(" ");
  mg->GetXaxis()->SetTitle("- V (bc) (V)");
  mg->GetYaxis()->SetTitle("I (mA)");

  const char *labels[2] = {"I=0.20 mA", "I=0.10 mA"};
  double x_100, y_100, x_200, y_200;
  I_100->GetPoint(2, x_100, y_100);
  I_200->GetPoint(2, x_200, y_200);
  TLatex *label_100 = new TLatex(x_100, y_100 * 1.8, labels[0]);
  label_100->SetTextSize(0.05);

  TLatex *label_200 = new TLatex(x_200, y_200 / 1.8, labels[1]);
  label_200->SetTextSize(0.05);

  // Calcolo di V Early aka a
  double Early_200 =
      std::abs(fitLin_200->GetParameter(0) / fitLin_200->GetParameter(1));
  double EarlyError_200 = std::sqrt(
      std::pow(fitLin_200->GetParError(0) / fitLin_200->GetParameter(1), 2) +
      std::pow(fitLin_200->GetParError(1) * fitLin_200->GetParameter(0) /
                   std::pow(fitLin_200->GetParameter(1), 2),
               2) -
      2 * fitLin_200->GetParameter(0) * cov_10_200 /
          std::pow(fitLin_200->GetParameter(1), 3));

  double Early_100 =
      std::abs(fitLin_100->GetParameter(0) / fitLin_100->GetParameter(1));
  double EarlyError_100 = std::sqrt(
      std::pow(fitLin_100->GetParError(0) / fitLin_100->GetParameter(1), 2) +
      std::pow(fitLin_100->GetParError(1) * fitLin_100->GetParameter(0) /
                   std::pow(fitLin_100->GetParameter(1), 2),
               2) -
      2 * fitLin_100->GetParameter(0) * cov_10_100 /
          std::pow(fitLin_100->GetParameter(1), 3));

  std::cout << '\n';
  std::cout << "chi_v 100 " << fitLin_100->GetChisquare() / fitLin_100->GetNDF()
            << "\n";
  std::cout << "chi_v 200 " << fitLin_200->GetChisquare() / fitLin_200->GetNDF()
            << "\n";
  std::cout << "covarianze " << cov_10_100 << " " << cov_10_200 << '\n';

  std::cout << "V Early I_100: " << Early_100 << " +/- " << EarlyError_100
            << '\n';
  std::cout << "V Early I_200: " << Early_200 << " +/- " << EarlyError_200
            << '\n';

  // calcolo resistenza in uscita aka b

  double RO_200 = std::abs(1 / fitLin_200->GetParameter(1));
  double ROError_200 = std::sqrt(std::pow(
      fitLin_200->GetParError(1) / std::pow(fitLin_200->GetParameter(1), 2),
      2));

  double RO_100 = std::abs(1 / fitLin_100->GetParameter(1));
  double ROError_100 = std::sqrt(std::pow(
      fitLin_100->GetParError(1) / std::pow(fitLin_100->GetParameter(1), 2),
      2));

  std::cout << "R Output I_100: " << RO_100 << " +/- " << ROError_100 << '\n';
  std::cout << "R Output I_200: " << RO_200 << " +/- " << ROError_200 << '\n';
  std::cout << "g I_100: " << 1 / RO_100 << " +/- " << ROError_100 / std::pow(RO_100,2) << '\n';
  std::cout << "g I_200: " << 1 / RO_200 << " +/- " << ROError_200 / std::pow(RO_200,2) << '\n';

  // Disegna i grafici
  TCanvas *c1 = new TCanvas("I_200", "I_200", 200, 10, 600, 400);
  I_200->Draw("APE");

  TCanvas *c2 = new TCanvas("I_100", "I_100", 200, 10, 600, 400);
  I_100->Draw("APE");

  TCanvas *c3 = new TCanvas("tot", "tot", 200, 10, 600, 400);
  mg->Draw("APE");
  label_100->Draw();
  label_200->Draw();

  // Calcolo del beta medio e del suo errore
  assert(I_200->GetN() == I_100->GetN());

  double beta{0};
  double betaError{0};
  double weightSum{0}; // Somma dei pesi

  for (int i{0}; i < I_200->GetN() - 10; i++) {
    double y_200, y_100;
    double err_200 = I_200->GetErrorY(i);
    double err_100 = I_100->GetErrorY(i);

    y_200 = I_200->GetPointY(i);
    y_100 = I_100->GetPointY(i);

    double beta_i = y_200 / y_100;
    double beta_i_err = beta_i * std::sqrt(std::pow(err_200 / y_200, 2) +
                                           std::pow(err_100 / y_100, 2));
    // ps è la formula delle derivate ma qui si semplifica l'espressione tirando
    // fuori il beta
    beta += beta_i / std::pow(beta_i_err, 2);   // Somma pesata
    weightSum += 1.0 / std::pow(beta_i_err, 2); // Somma dei pesi
  }

  beta /= weightSum;                      // Media pesata
  betaError = std::sqrt(1.0 / weightSum); // Errore della media pesata

  std::cout << "Beta medio: " << beta << " +/- " << betaError << '\n';
}

void fit(int sel = 1) {
  if (sel == 0) {
    sel_fit("I_100.txt", "I_100.txt");
  } else if (sel == 1) {
    sel_fit("I_100_media.txt", "I_200_media.txt");
  } else {
    std::cerr << "Errore: selezione non valida. Usa 0 o 1." << std::endl;
  }
}