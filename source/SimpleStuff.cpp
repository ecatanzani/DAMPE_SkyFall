
#include "MyHead.h"

void EvaluatePercentage(Double_t &percentage,Int_t entries,Int_t idx,std::ofstream &log_file) {
    if(((Double_t)idx/entries)>(percentage*0.01)) {
        std::cout<<"-> Tree events: [ "<<percentage<<" % ]"<<std::endl;
        log_file<<"-> Tree events: [ "<<percentage<<" % ]"<<std::endl;
        percentage++;
    }
}

void normalize_acceptance(TH2D &acceptance) {
    Int_t tot_events=get_histo_total_events(acceptance);
    Int_t tmp_events;
    
    for(Int_t bX=1; bX<=acceptance.GetNbinsX(); bX++) {
        for(Int_t bY=1; bY<=acceptance.GetNbinsY(); bY++) {
            tmp_events=acceptance.GetBinContent(bX,bY);
            acceptance.SetBinContent(bX,bY,(Double_t)tmp_events/tot_events);
        }
    }
        
}

Int_t get_histo_total_events(TH2D &histo) {
    Int_t nev=0;

    for(Int_t bX=1; bX<=histo.GetNbinsX(); bX++)
        for(Int_t bY=1; bY<=histo.GetNbinsY(); bY++)
            nev+=histo.GetBinContent(bX,bY);

    return nev;
}

void GetAcceptanceBins(TH2D &acceptance,Double_t costheta,Double_t phi,Int_t acc_bin[]) {
    
    TAxis *Yaxis = acceptance.GetYaxis();
    TAxis *Xaxis = acceptance.GetXaxis();
    
    acc_bin[0] = Xaxis->FindBin(costheta);
    acc_bin[1] = Yaxis->FindBin(phi);
    
}
