
#include "MyHead.h"

void obtain_histos(TH2D &evDist,TH2D &evDist_border,TH2D &acceptance,std::ofstream &log_file) {
    
    /////////////////////////////// Opening DAMPE acceptance 2D-histo
    //TFile* evDistr_file = new TFile(acceptance_final_plot.c_str());
    TFile infile(acceptance_final_plot.c_str());
    if(infile.IsZombie()) {
        std::cout << "\n\nError opening DAMPE Events Distribution TFile. Prorgram finished \n\n";
        log_file << "\n\nError opening DAMPE Events Distribution TFile. Prorgram finished \n\n";
        exit(-2);
    }
    
    new (&acceptance) TH2D(*(TH2D*)infile.Get("Acceptance"));
    new (&evDist) TH2D(*(TH2D*)infile.Get("EventsDistribution"));
    new (&evDist_border) TH2D(*(TH2D*)infile.Get("EventsDistribution"));
    
    acceptance.SetName("GAcceptance");
    evDist.SetName("evDistribution");
    evDist_border.SetName("evDistribution_border");
    
    evDist_border.Reset();
    get_evDist_border(evDist,evDist_border);
    infile.Close();
}
