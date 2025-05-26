#include "Pythia8/Pythia.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TParameter.h"

using namespace Pythia8;

int main() {
    // Initialize Pythia for pp collisions at 13 TeV.
    Pythia pythia;
    pythia.readString("Tune:pp = 5");
    pythia.readString("Beams:idA = 2212");  // Proton
    pythia.readString("Beams:idB = 2212");  // Proton
    pythia.readString("Beams:eCM = 13000.");  // Center-of-mass energy 13 TeV

    // Turn on hard QCD processes (minimum bias).
    pythia.readString("HardQCD:all = on");
    // Set min pt
    // pythia.readString("PhaseSpace:pTHatMin = 0.5");

    // Initialize Pythia
    if (!pythia.init()) return 1;
    int bin_size = 50;
    // Create histograms
    TH1F *h_pt_distribution = new TH1F("h_pt_distribution", "pT Distribution; pT (GeV/c); Number of Charged Particles", 100, 0.0, 20.0);

    TH1F *h_dNch_deta = new TH1F("h_dNch_deta", "dNch/deta", 50, -2.5, 2.5);
    TH1F *h_dNev_dnch = new TH1F("h_dNev_dnch", "dNev/dnch", 150, -0.5, 149.5);
    TH2F *h_meanPt_vs_nch = new TH2F("h_meanPt_vs_nch", "meanPt_vs_nch", 150, -0.5, 149.5, 200, 0.0, 20.0);
    TH1F *h_meanPt_vs_nch_1D = new TH1F("h_meanPt_vs_nch_1D", "meanPt_vs_nch_1D", 150, -0.5, 149.5);
    TH2F *h_d2Nch_deta_dpT = new TH2F("h_d2Nch_deta_dpT", "d2Nch/deta/dpT", 100, 0.0, 20.0, 50, -2.5, 2.5);
    TH1F *h_finalDistribution = new TH1F("h_finalDistribution", "1/(2 $\\pi$ pT) * d$^2$Nch/deta/dpT vs pT", 100, 0.0, 20.0);

    // Counter for events with n_ch > 0
    int Nev = 0;
    double mothertau;
    int nEvents = 8820000;
    for (int iEvent = 0; iEvent < nEvents; ++iEvent) {
        if (!pythia.next()) continue;

        int nCharged = 0;
        double sumPt = 0.0;

        for (int i = 0; i < pythia.event.size(); ++i) {
            if (pythia.event[i].isFinal() && pythia.event[i].isCharged()) {
                if (pythia.event[i].pT() > 0.5 && fabs(pythia.event[i].eta()) < 2.5) {
                    double tau = pythia.event[i].tau();
                    int status = pythia.event[i].status();     // Primary particle definition
                    bool isFromDecay = (status >= 91 && status <= 99);
                    if(isFromDecay){
                        mothertau = pythia.event[pythia.event[i].mother1()].tau();
                    }
                    bool notFromDecay = ! isFromDecay;
                    bool isPrimary = ((notFromDecay && tau > 300e-12) || (isFromDecay && (tau > 300e-12) && (mothertau < 30e-12)));
                    if (isPrimary) {  // Primary particle selection
                        ++nCharged;
                        double pt = pythia.event[i].pT();
                        double eta = pythia.event[i].eta();
                        h_pt_distribution->Fill(pt);
                        sumPt += pt;

                        // Fill histograms
                        h_dNch_deta->Fill(eta);
                        h_d2Nch_deta_dpT->Fill(pt, eta);
                    }
                }
            }
        }

        if (nCharged > 0) {
            ++Nev;  // Increment Nev if n_ch > 0
            h_dNev_dnch->Fill(nCharged);
            h_meanPt_vs_nch->Fill(nCharged, sumPt / nCharged);
        }
    }

    // Normalize histograms by Nev
    if (Nev > 0) {
        h_dNch_deta->Scale(1.0 /(Nev*h_dNch_deta.getBinWidth(1)));
        h_dNev_dnch->Scale(1.0 / Nev*h_dNev_dnch.getBinWidth(1));
        h_meanPt_vs_nch->Scale(1.0 / Nev);
        //h_pt_num->Scale(1.0/Nev);
        h_pt_distribution->Scale(1.0 / Nev*h_pt_distribution.getBinWidth(1));  // Optional: normalize by the number of events if needed

        //h_pt_num->Scale(1/bin_size);

        // Compute the mean pT for each n_ch bin and fill the 1D histogram
        int nBins = h_meanPt_vs_nch->GetNbinsX();
        for (int i = 1; i <= nBins; ++i) {
            double sumMeanPt = 0.0;
            double sumEntries = 0.0;

            for (int j = 1; j <= h_meanPt_vs_nch->GetNbinsY(); ++j) {
                double binContent = h_meanPt_vs_nch->GetBinContent(i, j);
                double meanPt = h_meanPt_vs_nch->GetYaxis()->GetBinCenter(j);
                sumMeanPt += meanPt * binContent;
                sumEntries += binContent;
            }

            if (sumEntries > 0) {
                double meanPtValue = sumMeanPt / sumEntries;
                h_meanPt_vs_nch_1D->SetBinContent(i, meanPtValue);
            }
        }

        // Normalize the 1D histogram
        h_meanPt_vs_nch_1D->Scale(1.0 / h_meanPt_vs_nch_1D.getBinWidth(1));

        // Sum over all eta bins for each pt bin and fill the 1D histogram for d2Nch/deta/dpT
        for (int i = 1; i <= h_d2Nch_deta_dpT->GetNbinsX(); ++i) {
            double pt = h_d2Nch_deta_dpT->GetXaxis()->GetBinCenter(i);
            double sum = 0.0;

            for (int j = 1; j <= h_d2Nch_deta_dpT->GetNbinsY(); ++j) {
                sum += h_d2Nch_deta_dpT->GetBinContent(i, j);
            }

            h_finalDistribution->SetBinContent(i, sum / (2 * M_PI * pt * Nev*h_finalDistribution.getBinWidth(1)));
        }
    }

    // Save histograms to a ROOT file
    TFile outFile("histograms_combined_ptsp.root", "RECREATE");
    h_pt_distribution->Write();
    h_dNch_deta->Write();
    h_dNev_dnch->Write();
    h_meanPt_vs_nch->Write();
    h_meanPt_vs_nch_1D->Write();
    h_d2Nch_deta_dpT->Write();
    h_finalDistribution->Write();
    

    // Write Nev to the output file
    TParameter<int>("Nev", Nev).Write();

    outFile.Close();

    pythia.stat();

    // Clean up
    delete h_pt_distribution;
    delete h_dNch_deta;
    delete h_dNev_dnch;
    delete h_meanPt_vs_nch;
    delete h_meanPt_vs_nch_1D;
    delete h_d2Nch_deta_dpT;
    delete h_finalDistribution;
    

    return 0;
}
