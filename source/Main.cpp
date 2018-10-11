///////////////////////////// Software description
//
//
//
//
/////////////////////////////////////////////////

#include "MyHead.h"

int main(int argc,char * argv[]) {
  
    ///////////////////////////////////////////////// VARIABLES DEFINITION
    
    ///////////// Variables for the TTrees

    Float_t g_lat = 0,g_lon = 0,geo_lat = 0,geo_lon = 0,sat_ra[3],sat_dec[3];
    UInt_t sec = 0;
    UShort_t n_ev = 0;
    bool good = true;
    Int_t chain_entries;                  // -> Variable that stores the whole number of tchain entries

    for(Int_t idx=0; idx<3; idx++) {
        sat_ra[idx] = 0;
        sat_dec[idx] = 0;
    }
    
    //////////////////////// Variable description

    /*
    
     sat_ra[3]     -> Is the array that stores the right ascension values of the satellite
     sat_dec[3]    -> Is the array that stores the celestial declination values of the satellite

     geo_lat       -> Is the geographic latitude
     geo_lon       -> Is the geographic longitude

     g_lat         -> Is the galactic latitude
     g_lon         -> Is the galactic longitude
    
                    These two are the "absolute coordinates" used to plot the final maps !!!!

   
     n_ev          -> Is the nuber of triggered events in a second
     sec           -> Is DAMPE's acquisition second number

     good          -> Is the status of the SBI

     */

    //////////////////////////////////////

    ///////////////////// Costheta flat binning variables
  
    Int_t n_bin_lon=360;                  // -> Number of bins along longitude axis

    Double_t lon_bin_min=-180.0;          // -> Set max and min for longitude binning
    Double_t lon_bin_max=180.0;

    Int_t n_bin_lat=180;                  // -> Number of bins along latitude axis
  
    Double_t lat_bin_min=-90.0;           // -> Set max and min for latitude binning
    Double_t lat_bin_max=90.0;

    Double_t* binning = nullptr;          // -> Array used to store the custom binning intervals !!!
    
    /////////////////////////////////////////////////////////////
    
    ///////////////////// Acceptance border POI vectors
  
    Double_t perc=0;                                                                             // -> Just used to store the percentage value
    
    Bool_t first_call = true;                                                                    // -> Boolean variable used to compute Events Distribution POI stuff
    Int_t n_peacks=0,n_valleys=0;                                                                // -> Number of peacks and valley into Events Distribution Histo
    
    std::vector<Double_t> peacks_costheta,valley_costheta;
    std::vector<Double_t> peacks_phi,valley_phi;

    std::vector<Double_t> Gpeack_l,Gpeack_b;
    std::vector<Double_t> Gvalley_l,Gvalley_b;
    
    ///////////////////// Stuff variables
    
    std::string log_path = output_path_creator(0),root_out_path = output_path_creator(1);        // -> Log end ROOT result file variables
    std::ofstream output_log_file(log_path);                                                     // -> log file creation !
    
    ////////////////////////////////////////////////////////////////////////////////////
    
    TFile results_file(root_out_path.c_str(),"RECREATE");                                        // -> output ROOT file
    results_file.cd();
    
    TChain tree("SBItree");                                                                      // -> Events ROOT TChain
    TRandom3 r_gen(random_seed);                                                                 // -> Set a random seed for TRandom3 objects
    
    Double_t choosen_l,choosen_b;                                                                // -> Randomly choosen galactic latitude and longitude from the Sky to Trace Back
    
    Double_t max_b=0,min_b=0,max_l=0,min_l=0;
    Double_t Nmax_b=0,Nmin_b=0,Pmax_b=0,Pmin_b=0;                                                // -> Variables to identify the Sky Zone to take events !
    
    Double_t costheta=0,phi=0;                                                                   // -> Costheta and Phi variables, final result of the trace back sky function
    Int_t tmp_rate=0;
    Int_t ev_counter=0;
    Int_t acc_bin[2];
    Double_t tmp_acc=0;
    Double_t l=-999,b=-999;
    
    /////////////////////////////////////////////////////////////////////////////////////////
    
    static TH2D evDist,evDist_border,acceptance;                                                 // -> Event Distribution and its border histos
    TH2D *Orbit = new TH2D("orbit","Exposure; Geographic Longitude (#circ);  Geographic Latitude (#circ); Exposure (s)", 180, -180, 180, 90, -90.0, 90.0);
    TH2D *Rate = new TH2D("rate","Rate; Geographic Longitude (#circ);  Geographic Latitude (#circ); Rate (hz)", 180, -180, 180, 90, -90.0, 90.0);
    
    create_and_initialize_log(output_log_file);
    read_SBI_data(g_lat,g_lon,geo_lat,geo_lon,sat_ra,sat_dec,sec,n_ev,good,tree,sbi_path,output_log_file);
  
    create_binning(n_bin_lat,lat_bin_min,lat_bin_max,binning,true);
    obtain_histos(evDist,evDist_border,acceptance,output_log_file);
    //normalize_acceptance(acceptance);
    obtain_orbit_rate(Orbit,Rate,tree,g_lat,g_lon,geo_lat,geo_lon,sat_ra,sat_dec,sec,n_ev,good,output_log_file);
    
    chain_entries = tree.GetEntries();
    results_file.cd();

    TH2D evRefMap("evRefMap","Reference Isotropic Sky Map; Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries",n_bin_lon,lon_bin_min,lon_bin_max,n_bin_lat,binning);
    
    if(evDist_back_procedure)
        gRandom->SetSeed(random_seed);
    
    /////////////////////////////////////////////////////////////////
    
    for(Int_t tree_idx=0; tree_idx<chain_entries; tree_idx++) {
        
        tree.GetEntry(tree_idx);
        if((sec%100000)==0)
            continue;                         //there was a bug in the SBI production.       !!!!!!! Ask is that problem is still present with new SBI files !!!!!!!!
        if(!good)
            continue;                         //good second
        
        ///////////// Setting correct range for galactic latitude and longitude
        if(g_lon>180)
            g_lon-=360;
        if(geo_lon>180)
            geo_lon-=360;
        
        EvaluatePercentage(perc,chain_entries,tree_idx,output_log_file);

        std::cout<<"Satellite longitude: "<<g_lon<<std::endl;
        
        ///////////// Setting arrays' values in radiants
        for(Int_t aidx=0; aidx<3; aidx++) {
            sat_ra[aidx]*=TMath::DegToRad();
            sat_dec[aidx]*=TMath::DegToRad();
        }
        
        if(evDist_back_procedure) {
            if(tree_idx!=0)
                ev_counter=0;
            extract_mean_rate(Rate,tmp_rate,geo_lon,geo_lat);
            while(ev_counter!=tmp_rate) {
                evDist.GetRandom2(costheta,phi);
            //    GetAcceptanceBins(acceptance,costheta,phi,acc_bin);
             //   tmp_acc = r_gen.Uniform();
                //std::cout<<"\nTmp acc: "<<tmp_acc<<"\t --> Acc: "<<acceptance.GetBinContent(acc_bin[0],acc_bin[1])<<std::endl;
             //   if(tmp_acc<=acceptance.GetBinContent(acc_bin[0],acc_bin[1])) {
                    from_local_to_galactic(costheta,phi,l,b,sat_ra,sat_dec);
                    if(isnan(l)) {
                        continue;
                    }
                    else {
                        if(l>180)
                            l-=360;
                    }
                    evRefMap.Fill(l,b,1);
                    ev_counter++;
            }
            
        }
        else {
            if(!POI_discrete_plotting(sat_ra,sat_dec,evDist_border,peacks_costheta,peacks_phi,valley_costheta,valley_phi,Gpeack_l,Gpeack_b,Gvalley_l,Gvalley_b,first_call,n_peacks,n_valleys)) {
                std::cout<<"\nJumped event "<<tree_idx<<" ! Event Distribution pattern ricognition error\n";
                output_log_file<<"\nJumped event "<<tree_idx<<" ! Event Distribution pattern ricognition error\n";
                continue;
            }
        
            if(points_are_inside(Gpeack_l,Gvalley_l,n_peacks,n_valleys)==true) {
                obtain_galactic_max_min(Gpeack_l,Gpeack_b,Gvalley_l,Gvalley_b,n_peacks,n_valleys,true,max_b,min_b,max_l,min_l,Nmax_b,Nmin_b,Pmax_b,Pmin_b);
                choosen_l=r_gen.Uniform(min_l,max_l);
                choosen_b=r_gen.Uniform(min_b,max_b);
            
            }
            else {
                obtain_galactic_max_min(Gpeack_l,Gpeack_b,Gvalley_l,Gvalley_b,n_peacks,n_valleys,false,max_b,min_b,max_l,min_l,Nmax_b,Nmin_b,Pmax_b,Pmin_b);
                choosen_l=r_gen.Uniform(min_l,180+(180-abs(max_l)));
                if(choosen_l>180)
                    choosen_b=r_gen.Uniform(Pmin_b,Pmax_b);
                else
                    choosen_b=r_gen.Uniform(Nmin_b,Nmax_b);
            }
            //  do {
            //    invert_map(costheta,phi,l,b,sat_ra,sat_dec,inv_vector)
        
            // Now I just have to trace back "choosen_l" and "choosen_b"
        }
        
        
    } //end loop on seconds

    ///////////////////////////////////////////////// Write output ROOT file
    
    evDist.Write();
    evDist_border.Write();
    acceptance.Write();
    
    results_file.Write();
    
    Orbit->Delete();
    Rate->Delete();

    ///////////////////////////////////////////////////////////////////////////////////
    
    std::cout<<"\n\nSimulation completed !!\n\n";
    output_log_file<<"\n\nSimulation completed !!\n\n";
  
}
