#include "Pythia8/Pythia.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TParameter.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TParticle.h"
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
    // int bin_size = 50;
    // Create histograms
    TH1F *hist_pt = new TH1F("hist_pt", "pT Distribution; pT (GeV/c); Number of Charged Particles", 100, 0.0, 20.0);

    TH1F *hist_eta = new TH1F("hist_eta", "dNch/deta", 50, -2.5, 2.5);
    TH1F *hist_mult = new TH1F("hist_mult", "dnevents/dnch", 200, -0.5, 199.5);
    TH2F *hist_nchpt = new TH2F("hist_nchpt", "meanPt_vs_nch", 150, -0.5, 149.5, 200, 0.0, 20.0);
    TH1F *hist_nchpt_1D = new TH1F("hist_nchpt_1D", "meanPt_vs_nch_1D", 150, -0.5, 149.5);
    TH1F *h_finalDistribution = new TH1F("h_finalDistribution", "1/(2 $\\pi$ pT) * d$^2$Nch/deta/dpT vs pT", 100, 0.0, 20.0);

    TH1I *h_rejection = new TH1I("h_rejection", "Cut vs number of tracks", 10, 0, 10);

    // Counter for events with n_ch > 0
    int nevents = 0;
    double mothertau;
    int neventsents = 5000000;

    for (int iEvent = 0; iEvent < neventsents; ++iEvent) {
        if (!pythia.next()) continue;

        int ncharged = 0;
        double sum_pt = 0.0;

        for (int i = 0; i < pythia.event.size(); ++i) {

            h_rejection->Fill(0.5);

            // if (!pythia.event[i].isFinal()) continue;
            h_rejection->Fill(1.5);

            if (!pythia.event[i].isCharged()) continue;
            h_rejection->Fill(2.5);

            if (pythia.event[i].pT() <=  0.5) continue;
            h_rejection->Fill(3.5);

            if (pythia.event[i].eta() >= 2.5) continue;
            h_rejection->Fill(4.5);

            int status = pythia.event[i].status();
            // if(status!=0)continue;
            h_rejection->Fill(5.5);

            //fill eta
            double eta = pythia.event[i].eta();
            hist_eta->Fill(eta);

            //fill pt
            double pt = pythia.event[i].pT();
            sum_pt += pt;
            hist_pt->Fill(pt);

            ++ncharged;
        }

    if (ncharged > 0) {
            ++nevents;  // Increment nevents if n_ch > 0
            hist_mult->Fill(ncharged);
            hist_nchpt->Fill(ncharged, sum_pt/ncharged);
        }
    }


    //Norm part
    if (nevents > 0) {

        hist_eta->Scale(1.0 /(nevents*hist_eta->GetBinWidth(1)));
        hist_mult->Scale(1.0 / (nevents*hist_mult->GetBinWidth(1)));

        hist_pt->Scale(1.0 / (2*TMath::Pi()*nevents*hist_pt->GetBinWidth(1)));
        for (int ibin = 1; ibin <= hist_nchpt->GetNbinsX(); ++ibin) {
            double bincontent = hist_pt->GetBinContent(ibin);
            double bincenter = hist_pt->GetBinCenter(ibin);
            hist_pt->SetBinContent(ibin, bincontent/bincenter);
        }
    }

    TCanvas *can = new TCanvas("test","test", 800, 800);
    can->Divide(2,2);
    can->cd(1);
    hist_eta->SetMarkerStyle(20);
    hist_eta->SetMarkerSize(0.80);
    hist_eta->SetMarkerColor(kRed);
    hist_eta->Draw("hist");
    can->cd(2);
    hist_pt->SetMarkerStyle(20);
    hist_pt->SetMarkerSize(0.80);
    hist_pt->SetMarkerColor(kRed);
    hist_pt->Draw("hist");
    can->cd(3);
    h_rejection->SetMarkerStyle(20);
    //hist_mult->SetMarkerSize(0.80);
    //hist_mult->SetMarkerColor(kRed);
    h_rejection->Draw("hist");
    can->SaveAs("test.pdf");

    return 0;
}
