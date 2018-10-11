
#include "MyHead.h"

void obtain_orbit_rate(TH2D *Orbit,TH2D *Rate,TChain &tree,Float_t &g_lat,Float_t &g_lon,Float_t &geo_lat,Float_t &geo_lon,Float_t sat_ra[],Float_t sat_dec[],UInt_t &sec,UShort_t &n_ev,bool &good,std::ofstream &log_file) {
    
    Int_t chain_entries=tree.GetEntries();
    
    for(Int_t tree_idx=0; tree_idx<chain_entries; tree_idx++) {
     
        tree.GetEntry(tree_idx);
        if((sec%100000)==0)
            continue;
        if(!good)
            continue;
        
        ///////////// Setting correct range for galactic latitude and longitude
        if(g_lon>180)
            g_lon-=360;
        if(geo_lon>180)
            geo_lon-=360;
        
        Orbit->Fill(geo_lon,geo_lat);
        Rate->Fill(geo_lon,geo_lat,n_ev);
        
    }
    
    Rate->Divide(Orbit);
    
    std::cout << "\n\nRate and Orbit histos has been calculated \n\n";
    log_file <<  "\n\nRate and Orbit histos has been calculated \n\n";
    
}

void extract_mean_rate(TH2D *Rate,Int_t &tmp_rate,Double_t geo_lon,Double_t geo_lat) {
    
    Int_t Xbin=0,Ybin=0;
    
    TAxis *Yaxis = Rate->GetYaxis();
    TAxis *Xaxis = Rate->GetXaxis();
    Ybin=Yaxis->FindBin(geo_lat);
    Xbin=Xaxis->FindBin(geo_lon);
    //lXbin=360./Xaxis->GetNbins();
    //lYbin=180./Yaxis->GetNbins();
    
    tmp_rate = Rate->GetBinContent(Xbin,Ybin);
    
}

bool POI_discrete_plotting(Float_t sat_ra[],Float_t sat_dec[],TH2D &evDist_border,std::vector<Double_t> &peacks_costheta,std::vector<Double_t> &peacks_phi,std::vector<Double_t> &valley_costheta,std::vector<Double_t> &valley_phi,std::vector<Double_t> &Gpeack_l,std::vector<Double_t> &Gpeack_b,std::vector<Double_t> &Gvalley_l,std::vector<Double_t> &Gvalley_b,Bool_t &first_call,Int_t &n_peacks,Int_t &n_valleys) {
    
    Bool_t good_PR=true;   //Bool variable. True if the pattern recognition of the event distribution works, false in the other case.
    
    // This kind of study, at the countrary, should be done on all DAMPE acquisition period
    if(first_call) {
        //calculate_POI_stuff(sat_ra,sat_dec,evDist_border,peacks_costheta,peacks_phi,valley_costheta,valley_phi);
        get_relevant_evDist_points(evDist_border,peacks_costheta,valley_costheta,peacks_phi,valley_phi);
        n_peacks=peacks_costheta.size();
        n_valleys=valley_costheta.size();
        
        std::cout<<"\n\nPeacks: "<<n_peacks<<std::endl;
        std::cout<<"Valleys: "<<n_valleys<<std::endl;
        
        Gpeack_l.resize(n_peacks);
        Gpeack_b.resize(n_peacks);
        Gvalley_l.resize(n_valleys);
        Gvalley_b.resize(n_valleys);
        
        GetPOI(peacks_costheta,peacks_phi,valley_costheta,valley_phi,Gpeack_l,Gpeack_b,Gvalley_l,Gvalley_b,sat_ra,sat_dec,good_PR,n_peacks,n_valleys);
        first_call=false;
    }

    GetPOI(peacks_costheta,peacks_phi,valley_costheta,valley_phi,Gpeack_l,Gpeack_b,Gvalley_l,Gvalley_b,sat_ra,sat_dec,good_PR,n_peacks,n_valleys);
    
    return good_PR;
    
}

void get_evDist_border(TH2D &evDist,TH2D &evDist_border) {
    
    for(Int_t y_bin=1; y_bin<=evDist.GetNbinsY(); y_bin++)
        for(Int_t x_bin=1; x_bin<=evDist.GetNbinsX(); x_bin++)
            if(evDist.GetBinContent(x_bin,y_bin)!=0) {
                evDist_border.SetBinContent(x_bin,y_bin,evDist.GetBinContent(x_bin,y_bin));
                break; //I found the first not-empty bin regarding a such phi (or y) value. Now I have to choose a new phi bin and search for the first not-empty costheta bin
            }
}

void get_relevant_evDist_points(TH2D &evDist_border,std::vector<Double_t> &peacks_costheta,std::vector<Double_t> &valley_costheta,std::vector<Double_t> &peacks_phi,std::vector<Double_t> &valley_phi) {
    
    Double_t l_binX,l_binY,costheta,phi;
    Double_t new_costheta,new_phi;
    Int_t bunch=get_bunch_dimension(evDist_border),n_bunch=evDist_border.GetNbinsY()/bunch;
    Bool_t first_b_bin,first_b=true;
    
    Double_t tmp_min_costheta1,tmp_max_costheta1,tmp_min_phi1,tmp_max_phi1;
    Double_t tmp_min_costheta2,tmp_max_costheta2,tmp_min_phi2,tmp_max_phi2;
    Double_t tmp_min_costheta3,tmp_max_costheta3,tmp_min_phi3,tmp_max_phi3;
    
    l_binX=(Double_t).5/evDist_border.GetNbinsX();
    l_binY=(Double_t)2*TMath::Pi()/evDist_border.GetNbinsY();
    
    for(Int_t idx_b=0; idx_b<n_bunch; idx_b++) {
        first_b_bin=true;
        //for(Int_t idx_bY=(((idx_b+3)*bunch)+1); idx_bY<=((idx_b+4)*bunch); idx_bY++) {
        for(Int_t idx_bY=((idx_b*bunch)+1); idx_bY<=((idx_b+1)*bunch); idx_bY++) {
            for(Int_t idx_bX=1; idx_bX<=evDist_border.GetNbinsX(); idx_bX++) {
                if(evDist_border.GetBinContent(idx_bX,idx_bY)!=0) {
                    costheta=0.5+.5*(idx_bX*l_binX+(idx_bX-1)*l_binX);
                    phi=.5*(idx_bY*l_binY+(idx_bY-1)*l_binY);
                    if(first_b_bin) {
                        tmp_min_costheta1=tmp_max_costheta1=costheta;
                        tmp_min_phi1=tmp_max_phi1=phi;
                        first_b_bin=false;
                    }
                    else {
                        if(costheta>tmp_max_costheta1) {
                            tmp_max_costheta1=costheta;
                            tmp_max_phi1=phi;
                        }
                        if(costheta<tmp_min_costheta1) {
                            tmp_min_costheta1=costheta;
                            tmp_min_phi1=phi;
                        }
                    }
                }
            }
        }
        first_b_bin=true;
        //for(Int_t idx_bY=(((idx_b+3)*bunch)+1); idx_bY<=((idx_b+4)*bunch); idx_bY++) {
        for(Int_t idx_bY=(((idx_b+1)*bunch)+1); idx_bY<=((idx_b+2)*bunch); idx_bY++) {
            for(Int_t idx_bX=1; idx_bX<=evDist_border.GetNbinsX(); idx_bX++) {
                if(evDist_border.GetBinContent(idx_bX,idx_bY)!=0) {
                    new_costheta=0.5+.5*(idx_bX*l_binX+(idx_bX-1)*l_binX);
                    new_phi=.5*(idx_bY*l_binY+(idx_bY-1)*l_binY);
                    if(first_b_bin) {
                        tmp_min_costheta3=tmp_max_costheta3=new_costheta;
                        tmp_min_phi3=tmp_max_phi3=new_phi;
                        first_b_bin=false;
                    }
                    else {
                        if(costheta>tmp_max_costheta3) {
                            tmp_max_costheta3=new_costheta;
                            tmp_max_phi3=new_phi;
                        }
                        if(costheta<tmp_min_costheta3) {
                            tmp_min_costheta3=new_costheta;
                            tmp_min_phi3=new_phi;
                        }
                    }
                }
            }
        }
        if(first_b)
            first_b=false;
        else {
            if(tmp_min_costheta1<tmp_min_costheta2 && tmp_min_costheta1<tmp_min_costheta3) {
                peacks_costheta.push_back(tmp_min_costheta1);
                peacks_phi.push_back(tmp_min_phi1);
            }
            if(tmp_max_costheta1>tmp_max_costheta2 && tmp_max_costheta1>tmp_max_costheta3) {
                valley_costheta.push_back(tmp_max_costheta1);
                valley_phi.push_back(tmp_max_phi1);
            }
        }
        tmp_min_costheta2=tmp_min_costheta1;
        tmp_max_costheta2=tmp_max_costheta1;
        tmp_min_phi2=tmp_min_phi1;
        tmp_max_phi2=tmp_max_phi1;
    }
    
}

Int_t get_bunch_dimension(TH2D &evDist_border) {
    Int_t bin_num=20,Ybins=evDist_border.GetNbinsY();
    Bool_t found=false;
    
    while(found==false)
        if((Ybins%bin_num)==0)
            found=true;
        else
            bin_num++;
    
    return bin_num;
}

void GetPOI(std::vector<Double_t> &peacks_costheta,std::vector<Double_t> &peacks_phi,std::vector<Double_t> &valley_costheta,std::vector<Double_t> &valley_phi,std::vector<Double_t> &Gpeack_l,std::vector<Double_t> &Gpeack_b,std::vector<Double_t> &Gvalley_l,std::vector<Double_t> &Gvalley_b,Float_t sat_ra[],Float_t sat_dec[],Bool_t &good_PR,Int_t n_peacks,Int_t n_valleys) {
    
    Double_t l=0,b=0;
    
    for(Int_t idx_p=0; idx_p<n_peacks; idx_p++) {
        from_local_to_galactic(peacks_costheta.at(idx_p),peacks_phi.at(idx_p),l,b,sat_ra,sat_dec);
        if(isnan(l)) {
            good_PR=false;
            break;
        }
        if(l>180.0)
            l-=360.0;
        Gpeack_b[idx_p]=b;
        Gpeack_l[idx_p]=l;
    }
    
    if(good_PR) {
        for(Int_t idx_v=0; idx_v<n_valleys; idx_v++) {
            from_local_to_galactic(valley_costheta.at(idx_v),valley_phi.at(idx_v),l,b,sat_ra,sat_dec);
            if(isnan(l)) {
                good_PR=false;
                break;
            }
            if(l>180.0)
                l-=360.0;
            Gvalley_b[idx_v]=b;
            Gvalley_l[idx_v]=l;
        }
    }
}

Double_t obtain_min_histo(TH1D *histo) {
    Double_t min_histo=-999,bin_lenght;
    
    bin_lenght=(Double_t).5/histo->GetNbinsX();
    
    for(Int_t idx_b=1; idx_b<=histo->GetNbinsX(); idx_b++)
        if(histo->GetBinContent(idx_b)!=0) {
            min_histo=0.5*(idx_b*bin_lenght+(idx_b-1)*bin_lenght);
            break; //the first not-empty bin has been found. No need to go ahead
        }
    
    return min_histo;
}


/////////////////////////////////////////////////////////////////////////////////////////



bool points_are_inside(std::vector<Double_t> &Gpeack_l,std::vector<Double_t> &Gvalley_l,Int_t n_peacks,Int_t n_valleys) {
    bool all_points_inside=true;

    /*
    for(Int_t idx_v=0; idx_v<5; idx_v++) {
        if(idx_v==0) {
            if(valley_glon[idx_v]>peacks_glon[0]) {
                all_points_inside=false;
                break;
            }
        }
        if(idx_v>0 && idx_v<3) {
            if(peacks_glon[idx_v]>peacks_glon[0] || valley_glon[idx_v]>peacks_glon[0]) {
                all_points_inside=false;
                break;
            }
        }
        else {
            if(peacks_glon[idx_v]>peacks_glon[0]) {
                all_points_inside=false;
                break;
            }
        }
    }
     */
    
    Double_t max=Gpeack_l[0],min=max;
    
    for(Int_t idx_v=1; idx_v<n_peacks; idx_v++) {
        if(Gpeack_l[idx_v]>max)
            max=Gpeack_l[idx_v];
        if(Gpeack_l[idx_v]<min)
            min=Gpeack_l[idx_v];
    }
    
    for(Int_t idx_v=0; idx_v<n_valleys; idx_v++) {
        if(Gvalley_l[idx_v]>max)
            max=Gvalley_l[idx_v];
        if(Gvalley_l[idx_v]<min)
            min=Gvalley_l[idx_v];
    }
    
    if((max-min)>180.)
        all_points_inside=false;
        
    return all_points_inside;
}

void obtain_galactic_max_min(std::vector<Double_t> &Gpeack_l,std::vector<Double_t> &Gpeack_b,std::vector<Double_t> &Gvalley_l,std::vector<Double_t> &Gvalley_b,Int_t n_peacks,Int_t n_valleys,bool points_inside,Double_t &max_b,Double_t &min_b,Double_t &max_l,Double_t &min_l,Double_t &Nmax_b,Double_t &Nmin_b,Double_t &Pmax_b,Double_t &Pmin_b) {
    
    if(points_inside==true) {
        max_l=min_l=Gpeack_l[0];
        for(Int_t idx_v=1; idx_v<n_peacks; idx_v++) {
            if(Gpeack_l[idx_v]>max_l)
                max_l=Gpeack_l[idx_v];
            if(Gpeack_l[idx_v]<min_l)
                min_l=Gpeack_l[idx_v];
        }
        for(Int_t idx_v=0; idx_v<n_valleys; idx_v++) {
            if(Gvalley_l[idx_v]>max_l)
                max_l=Gvalley_l[idx_v];
            if(Gvalley_l[idx_v]<min_l)
                min_l=Gvalley_l[idx_v];
        }
        
        max_b=min_b=Gpeack_b[0];
        for(Int_t idx_v=1; idx_v<n_peacks; idx_v++) {
            if(Gpeack_b[idx_v]>max_b)
                max_b=Gpeack_b[idx_v];
            if(Gpeack_b[idx_v]<min_b)
                min_b=Gpeack_b[idx_v];
        }
        for(Int_t idx_v=0; idx_v<n_valleys; idx_v++) {
            if(Gvalley_b[idx_v]>max_b)
                max_b=Gpeack_b[idx_v];
            if(Gpeack_b[idx_v]<min_b)
                min_b=Gpeack_b[idx_v];
        }
    }
    else {
        max_l=-180;
        min_l=180;
        Nmax_b=Pmax_b=-90;
        Nmin_b=Pmin_b=90;
        
        for(Int_t idx_v=0; idx_v<n_peacks; idx_v++) {
            if(Gpeack_l[idx_v]<0) {
                if(Gpeack_l[idx_v]>max_l) {
                    max_l=Gpeack_l[idx_v];
                }
                if(Gpeack_b[idx_v]>Nmax_b) {
                    Nmax_b=Gpeack_b[idx_v];
                }
                if(Gpeack_b[idx_v]<Nmin_b) {
                    Nmin_b=Gpeack_b[idx_v];
                }
            }
            else {
                if(Gpeack_l[idx_v]<min_l) {
                    min_l=Gpeack_l[idx_v];
                }
                if(Gpeack_b[idx_v]>Pmax_b) {
                    Pmax_b=Gpeack_b[idx_v];
                }
                if(Gpeack_b[idx_v]<Pmin_b) {
                    Pmin_b=Gpeack_b[idx_v];
                }
            }
        }
        
        for(Int_t idx_v=0; idx_v<n_valleys; idx_v++) {
            if(Gvalley_l[idx_v]<0) {
                if(Gvalley_l[idx_v]>max_l) {
                    max_l=Gvalley_l[idx_v];
                }
                if(Gvalley_b[idx_v]>Nmax_b) {
                    Nmax_b=Gvalley_b[idx_v];
                }
                if(Gvalley_b[idx_v]<Nmin_b) {
                    Nmin_b=Gvalley_b[idx_v];
                }
            }
            else {
                if(Gvalley_l[idx_v]<min_l) {
                    min_l=Gvalley_l[idx_v];
                }
                if(Gvalley_b[idx_v]>Pmax_b) {
                    Pmax_b=Gvalley_b[idx_v];
                }
                if(Gvalley_b[idx_v]<Pmin_b) {
                    Pmin_b=Gvalley_b[idx_v];
                }
            }
        }
    }
}
