

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
  mg->GetYaxis()->SetTitleOffset(0.7);
  mg->GetXaxis()->SetTitleOffset(0.85);
  mg->GetXaxis()->SetTitleSize(0.06);
  mg->GetYaxis()->SetTitleSize(0.06);
  mg->GetXaxis()->SetLabelSize(0.06);
  mg->GetYaxis()->SetLabelSize(0.06);

  const char *labels[2] = {"I=0.20 mA", "I=0.10 mA"};
  double x_100, y_100, x_200, y_200;
  I_100->GetPoint(2, x_100, y_100);
  I_200->GetPoint(2, x_200, y_200);
  TLatex *label_100 = new TLatex(x_100, y_100 * 1.8, labels[0]);
  label_100->SetTextSize(0.05);

  TLatex *label_200 = new TLatex(x_200, y_200 / 1.8, labels[1]);
  label_200->SetTextSize(0.05);

  // Calcolo di V Early aka a
  const double intercept_200 = fitLin_200->GetParameter(0);
  const double slope_200 = fitLin_200->GetParameter(1);
  const double interceptError_200 = fitLin_200->GetParError(0);
  const double slopeError_200 = fitLin_200->GetParError(1);
 
  const double Early_200 = std::abs(intercept_200 / slope_200);

  // Calcolo dell'errore di Early_200
  const double earlyErrorTerm1 = std::pow(interceptError_200 / slope_200,
                                    2); // (Errore intercetta / pendenza)^2
  const double earlyErrorTerm2 =
      std::pow(intercept_200 * slopeError_200 / std::pow(slope_200, 2),
               2); // (Errore pendenza * intercetta / pendenza^2)^2
  const double earlyErrorTerm3 = -2 * intercept_200 * cov_10_200 /
                           std::pow(slope_200, 3); // Termino di covarianza

  const double EarlyError_200 =
      std::sqrt(earlyErrorTerm1 + earlyErrorTerm2 + earlyErrorTerm3);

  // Output dei risultati
  std::cout << "V Early I_200: " << Early_200 << " +/- " << EarlyError_200
            << '\n';

  // Calcolo della tensione di Early (V Early) per I_100
  const double intercept_100 = fitLin_100->GetParameter(0);
  const double slope_100 = fitLin_100->GetParameter(1);
  const double interceptError_100 = fitLin_100->GetParError(0);
  const double slopeError_100 = fitLin_100->GetParError(1);
  const double Early_100 = std::abs(intercept_100 / slope_100);

  // Calcolo dell'errore di Early_100
  const double earlyErrorTerm1_100 = std::pow(interceptError_100 / slope_100,
                                        2); // (Errore intercetta / pendenza)^2
  const double earlyErrorTerm2_100 =
      std::pow(intercept_100 * slopeError_100 / std::pow(slope_100, 2),
               2); // (Errore pendenza * intercetta / pendenza^2)^2
  const double earlyErrorTerm3_100 = -2 * intercept_100 * cov_10_100 /
                               std::pow(slope_100, 3); // Termine di covarianza

  const double EarlyError_100 = std::sqrt(earlyErrorTerm1_100 + earlyErrorTerm2_100 +
                                    earlyErrorTerm3_100);

  std::cout << '\n';
  std::cout << "chi_v 100 " << fitLin_100->GetChisquare() / fitLin_100->GetNDF()
            << "\n";
  std::cout << "chi_v 200 " << fitLin_200->GetChisquare() / fitLin_200->GetNDF()
            << "\n";
  std::cout << "covarianze " << cov_10_100 << " " << cov_10_200 << '\n';

  std::cout << "V Early I_100: " << Early_100 << " +/- " << EarlyError_100
            << " V " << '\n';
  std::cout << "V Early I_200: " << Early_200 << " +/- " << EarlyError_200
            << " V " << '\n';

  // calcolo resistenza in uscita aka b

  const double RO_200 = std::abs(1 / fitLin_200->GetParameter(1));
  const double ROError_200 = std::sqrt(std::pow(
      fitLin_200->GetParError(1) / std::pow(fitLin_200->GetParameter(1), 2),
      2));

  const double RO_100 = std::abs(1 / fitLin_100->GetParameter(1));
  const double ROError_100 = std::sqrt(std::pow(
      fitLin_100->GetParError(1) / std::pow(fitLin_100->GetParameter(1), 2),
      2));

  std::cout << "R Output I_100: " << RO_100 << " +/- " << ROError_100 << " kohm"
            << '\n';
  std::cout << "R Output I_200: " << RO_200 << " +/- " << ROError_200 << " kohm"
            << '\n';

  std::cout << "g I_100: " << fitLin_100->GetParameter(1) << " +/- "
            << fitLin_100->GetParError(1) << " mS"
            << '\n'; // milli Simens
  std::cout << "g I_200: " << fitLin_200->GetParameter(1) << " +/- "
            << fitLin_200->GetParError(1) << " mS" << '\n';

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
  double const dCurrent{0.1};
  double const errdCurrent{(0.022+0.023)/std::sqrt(3)};
  //ho considetato la differenza come errori massimi e poi ho scalato a quelli casuali
  //così è più semplice l'espressione finale

  for (int i{0}; i < I_200->GetN() - 10; i++) {
    double y_200, y_100;
    double err_200 = I_200->GetErrorY(i);
    double err_100 = I_100->GetErrorY(i);

    y_200 = I_200->GetPointY(i);
    y_100 = I_100->GetPointY(i);

    double beta_i = (y_200 - y_100)/dCurrent;
    double beta_i_err2 = std::pow(err_200/dCurrent,2)+ std::pow(err_100/dCurrent,2)+ std::pow(beta_i* errdCurrent/dCurrent,2); 
    //è l'errore del singolo elemento beta
    beta += beta_i / beta_i_err2;   // Somma pesata
    weightSum += 1.0 / beta_i_err2; // Somma dei pesi
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