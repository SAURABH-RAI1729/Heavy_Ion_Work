#include "Pythia8/Pythia.h"
#include <vector>
#include <iostream>
#include <TFile.h>
#include <TTree.h>

using namespace Pythia8;

int main() {
    // Initialize Pythia
    Pythia pythia;
    pythia.readString("Beams:eCM = 5020.");  // Set center-of-mass energy to 5.02 TeV
    pythia.readString("SoftQCD:all = on");   // Enable SoftQCD processes
    pythia.readString("ColourReconnection:reconnect = off");  // Turn CR off
    pythia.init();

    // Create a ROOT file to store the data
    TFile *outputFile = new TFile("event_data.root", "RECREATE");
    TTree *tree = new TTree("EventTree", "Tree with event data");

    // Variables to store in the tree
    float pt, eta;
    int nCharged, eventIndex;
    
    // Create branches in the tree
    tree->Branch("pt", &pt, "pt/F");
    tree->Branch("eta", &eta, "eta/F");
    tree->Branch("nCharged", &nCharged, "nCharged/I");
    tree->Branch("eventIndex", &eventIndex, "eventIndex/I");

    // Variables to store event information
    const int numEvents = 10000000; // Adjust this number based on your statistics needs

    // Event loop
    for (int iEvent = 0; iEvent < numEvents; ++iEvent) {
        if (!pythia.next()) continue;

        std::vector<double> pT_values;
        nCharged = 0;  // Reset count of selected charged particles in the event
        eventIndex = iEvent;  // Current event index

        for (int i = 0; i < pythia.event.size(); ++i) {
            if (pythia.event[i].isFinal() && pythia.event[i].isCharged()) {
                eta = pythia.event[i].eta();
                pt = pythia.event[i].pT();
                if (fabs(eta) < 0.8 && pt > 0.15 && pt < 2.0) {
                    nCharged++;  // Increment number of charged particles
                    tree->Fill();  // Fill the tree with the current particle data
                }
            }
        }
    }

    // Write the tree to the file and close the file
    tree->Write();
    outputFile->Close();

    // Print statistics
    pythia.stat();
    return 0;
}
