#include <vector>
#include <map>
#include <cmath>
#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TGraph.h>
#include <TCanvas.h>

int main() {
    // Open the ROOT file and get the tree
    TFile *inputFile = new TFile("event_data.root", "READ");
    if (!inputFile || inputFile->IsZombie()) {
        std::cerr << "Error opening file or file is corrupted." << std::endl;
        return 1;
    }

    TTree *tree = (TTree*)inputFile->Get("EventTree");
    if (!tree) {
        std::cerr << "Error: Could not find tree 'EventTree' in file." << std::endl;
        inputFile->Close();
        return 1;
    }

    // Variables to read from the tree
    float pt, eta;
    int nCharged, eventIndex;

    // Set branch addresses
    tree->SetBranchAddress("pt", &pt);
    tree->SetBranchAddress("eta", &eta);
    tree->SetBranchAddress("nCharged", &nCharged);
    tree->SetBranchAddress("eventIndex", &eventIndex);

    // Define histogram for dnch/deta
    int numBins = 50;
    TH1F *dnch_deta_hist = new TH1F("dnchdeta", "<dNch/deta>", numBins, -0.8, 0.8);

    // Variables to store event information
    std::map<int, std::vector<std::vector<double>>> multiplicityClasses; // Map to store multiplicity classes
    std::map<int, double> average_dnch_deta_per_class; // Map to store average dnch/deta for each class

    // Read data from the tree
    int nEntries = tree->GetEntries();
    int currentEventIndex = -1;
    std::vector<double> pT_values;
    int n_ch = 0;

    for (int i = 0; i < nEntries; ++i) {
        tree->GetEntry(i);

        // Detect a new event
        if (eventIndex != currentEventIndex) {
            if (n_ch > 2 && !pT_values.empty()) { // Only consider events with N_ch > 2
                multiplicityClasses[n_ch].push_back(pT_values);
            }
            // Reset for new event
            pT_values.clear();
            currentEventIndex = eventIndex;
            n_ch = 0;
        }

        // Store current particle's pT and count charged particles
        pT_values.push_back(pt);
        dnch_deta_hist->Fill(eta); // Fill histogram for dnch/deta
        n_ch = nCharged; // Update the number of charged particles for the current event
    }

    // Process the last event in the file
    if (n_ch > 2 && !pT_values.empty()) {
        multiplicityClasses[n_ch].push_back(pT_values);
    }

    // Calculate <dNch/deta> for each multiplicity class
    for (auto const& mult_class : multiplicityClasses) {
        int multiplicity = mult_class.first;
        double total_dNch_deta = 0;
        int event_count = mult_class.second.size();

        for (const auto& pt_list : mult_class.second) {
            total_dNch_deta += (double)pt_list.size() / (1.6); // Width of eta range (-0.8 to 0.8)
        }

        average_dnch_deta_per_class[multiplicity] = total_dNch_deta / event_count; // Average dnch/deta per multiplicity class
    }

    // Prepare data for graph
    std::vector<double> dnch_deta_values;
    std::vector<double> standard_skewness_values;
    std::vector<double> intensive_skewness_values;
    std::vector<double> kurtosis_values;

    // Calculate mean pT, Cm, three-point and four-point correlators for each multiplicity class
    for (auto const& mult_class : multiplicityClasses) {
        int multiplicity = mult_class.first;
        const auto& pts = mult_class.second;
        double mean_pt_m = 0;
        int total_count = 0;

        // Calculate mean pT for this class
        for (const auto& pt_list : pts) {
            for (double pt : pt_list) {
                mean_pt_m += pt;
                total_count++;
            }
        }
        if (total_count == 0) continue; // Prevent division by zero
        mean_pt_m /= total_count;

        double cm = 0;
        double three_point_corr = 0;
        double four_point_corr = 0;
        int total_pairs = 0;
        int total_triplets = 0;
        int total_quadruplets = 0;

        // Calculate Cm, three-point, and four-point correlators
        for (const auto& pt_list : pts) {
            int n = pt_list.size();
            total_pairs += n * (n - 1)/2;
            total_triplets += n * (n - 1) * (n - 2)/6;
            total_quadruplets += n * (n - 1) * (n - 2) * (n - 3)/24;

            for (size_t i = 0; i < pt_list.size(); ++i) {
                for (size_t j = i + 1; j < pt_list.size(); ++j) {
                    cm += (pt_list[i] - mean_pt_m) * (pt_list[j] - mean_pt_m);
                    for (size_t k = j + 1; k < pt_list.size(); ++k) {
                        three_point_corr += (pt_list[i] - mean_pt_m) * (pt_list[j] - mean_pt_m) * (pt_list[k] - mean_pt_m);
                        for (size_t l = k + 1; l < pt_list.size(); ++l) {
                            four_point_corr += (pt_list[i] - mean_pt_m) * (pt_list[j] - mean_pt_m) * (pt_list[k] - mean_pt_m) * (pt_list[l] - mean_pt_m);
                        }
                    }
                }
            }
        }

        if (total_pairs > 0) cm /= total_pairs; // Normalize Cm by the number of pairs
        if (total_triplets > 0) three_point_corr /= total_triplets; // Normalize three-point correlator
        if (total_quadruplets > 0) four_point_corr /= total_quadruplets; // Normalize four-point correlator

        // Calculate skewness and kurtosis
        double std_skewness = (three_point_corr > 0 && cm > 0) ? (three_point_corr / pow(cm, 1.5)) : 0;
        double intensive_skewness = (cm > 0) ? (three_point_corr * mean_pt_m) / (cm * cm) : 0;
        double kurtosis = (cm > 0) ? (four_point_corr / (cm * cm)) : 0;

        // Filter out extreme values due to insufficient statistics
        if (std::abs(std_skewness) <= 2 && std::abs(intensive_skewness) <= 20 && std::abs(kurtosis) <= 10) {
            // Store values for plotting
            double dnch_deta = average_dnch_deta_per_class[multiplicity];
            dnch_deta_values.push_back(pow(dnch_deta, 1.0/3.0));
            standard_skewness_values.push_back(std_skewness);
            intensive_skewness_values.push_back(intensive_skewness);
            kurtosis_values.push_back(kurtosis);
        }
    }

    // Create graphs for plotting
    TGraph *graph_std_skewness = new TGraph(dnch_deta_values.size(), &dnch_deta_values[0], &standard_skewness_values[0]);
    TGraph *graph_int_skewness = new TGraph(dnch_deta_values.size(), &dnch_deta_values[0], &intensive_skewness_values[0]);
    TGraph *graph_kurtosis = new TGraph(dnch_deta_values.size(), &dnch_deta_values[0], &kurtosis_values[0]);

    // Configure and plot the graphs with y-axis limits
    TCanvas *canvas1 = new TCanvas("canvas1", "Standard Skewness vs. <dNch/deta>^(1/3)", 800, 600);
    graph_std_skewness->SetTitle("Standard Skewness vs. <dNch/deta>^(1/3)");
    graph_std_skewness->GetXaxis()->SetTitle("<dNch/deta>^(1/3)");
    graph_std_skewness->GetYaxis()->SetTitle("Standard Skewness");
    graph_std_skewness->GetYaxis()->SetRangeUser(0, 2);  // Set y-axis limit
    //canvas1->SetLogx();
    graph_std_skewness->Draw("AP*");
    canvas1->SaveAs("StandardSkewnessPlotOFF.pdf");

    TCanvas *canvas2 = new TCanvas("canvas2", "Intensive Skewness vs. <dNch/deta>^(1/3)", 800, 600);
    graph_int_skewness->SetTitle("Intensive Skewness vs. <dNch/deta>^(1/3)");
    graph_int_skewness->GetXaxis()->SetTitle("<dNch/deta>^(1/3)");
    graph_int_skewness->GetYaxis()->SetTitle("Intensive Skewness");
    graph_int_skewness->GetYaxis()->SetRangeUser(0, 20);  // Set y-axis limit
    //canvas2->SetLogx();
    graph_int_skewness->Draw("AP*");
    canvas2->SaveAs("IntensiveSkewnessPlotOFF.pdf");

    TCanvas *canvas3 = new TCanvas("canvas3", "Kurtosis vs. <dNch/deta>^(1/3)", 800, 600);
    graph_kurtosis->SetTitle("Kurtosis vs. <dNch/deta>^(1/3)");
    graph_kurtosis->GetXaxis()->SetTitle("<dNch/deta>^(1/3)");
    graph_kurtosis->GetYaxis()->SetTitle("Kurtosis");
    graph_kurtosis->GetYaxis()->SetRangeUser(0, 10);  // Set y-axis limit
    //canvas3->SetLogx();
    graph_kurtosis->Draw("AP*");
    canvas3->SaveAs("KurtosisPlotOFF.pdf");

    // Close the input file
    inputFile->Close();

    return 0;
}
