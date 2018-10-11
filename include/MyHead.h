
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <fstream>
#include <ctime>
#include <vector>
#include <string>
#include <cmath>

///// ROOT libraries 

#include "TRandom3.h"
#include "TMath.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TH2.h"
#include "TGraph.h"
#include "TStyle.h"
#include "TNamed.h"
#include "TFile.h"
#include "TChain.h"
#include "TColor.h"
#include "TLine.h"
#include "TApplication.h"
#include "TROOT.h"
#include "TMatrixD.h"
#include "TAxis.h"
#include "TTreeIndex.h"

#include "orbitStruct.h"

#define EPS 1.e-12
#define err_inv 1.e-3

///////////////////////////////////// Simulation variables

//const static Int_t sky_events = 1e+9;
const static Int_t sky_events = 1e+3;
const static UInt_t random_seed = 22;
const static time_t time_stamp = time(0);       //Setting timestamp for the out files

const static TString sbi_path = "/storage/gpfs_data/dampe/users/ecatanzani/MyRepos/DAMPE/SkyFall/SampleSBI/";
const static std::string acceptance_final_plot = "/storage/gpfs_data/dampe/users/ecatanzani/MyRepos/DAMPE/Stuff/DAMPE-GAcceptance/results/1527159147_acceptance_result.root";
const static TString sbi_subsample = "010";
const static std::string string_sbi_subsample = "010";
const static Int_t number_SBI_files = 3;       // To be precise, at the moment of sotware writing, they are 0102800000_SBI.root 0102900000_SBI.root 0103000000_SBI.root

const static std::string output_log = "/storage/gpfs_data/dampe/users/ecatanzani/MyRepos/DAMPE/SkyFall/logs/";
const static std::string output_root = "/storage/gpfs_data/dampe/users/ecatanzani/MyRepos/DAMPE/SkyFall/results/";

const static Bool_t evDist_back_procedure = true;

////////////////////////////////////////////////////////////////////////

extern void create_and_initialize_log(std::ofstream &log_file);
extern void obtain_histos(TH2D &evDist,TH2D &evDist_border,TH2D &acceptance,std::ofstream &log_file);
extern void obtain_orbit_rate(TH2D *Orbit,TH2D *Rate,TChain &tree,Float_t &g_lat,Float_t &g_lon,Float_t &geo_lat,Float_t &geo_lon,Float_t sat_ra[],Float_t sat_dec[],UInt_t &sec,UShort_t &n_ev,bool &good,std::ofstream &log_file);
extern void normalize_acceptance(TH2D &acceptance);
extern Int_t get_histo_total_events(TH2D &histo);
extern void GetAcceptanceBins(TH2D &acceptance,Double_t costheta,Double_t phi,Int_t acc_bin[]);
extern void extract_mean_rate(TH2D *Rate,Int_t &tmp_rate,Double_t geo_lon,Double_t geo_lat);
extern std::string output_path_creator(const Int_t out_choose);
extern void log_file_init(std::ofstream &out_file);
extern void read_SBI_data(Float_t &galactic_lat,Float_t &galactic_lon,Float_t &geographic_lat,Float_t &geographic_lon,Float_t sat_ra[],Float_t sat_dec[],UInt_t &sec,UShort_t &n_events,bool &good,TChain &tree,TString sbi_data_path,std::ofstream &out_file);
extern Bool_t check_sbi_loading(Float_t galactic_lat,Float_t galactic_lon,Float_t geographic_lat,Float_t geographic_lon,Float_t sat_ra[],Float_t sat_dec[],UInt_t sec,UShort_t n_events);
extern Bool_t chech_if_null_variable(Float_t in_variable);
extern void reinitialize_all_variables(Float_t &galactic_lat,Float_t &galactic_lon,Float_t &geographic_lat,Float_t &geographic_lon,Float_t sat_ra[],Float_t sat_dec[],UInt_t &sec,UShort_t &n_events);
extern void create_binning(Int_t n_bin_lat,Double_t lat_bin_min,Double_t lat_bin_max,Double_t* &binning,Bool_t cos_scale);
extern void EvaluatePercentage(Double_t &percentage,Int_t entries,Int_t idx,std::ofstream &log_file);
extern bool POI_discrete_plotting(Float_t sat_ra[],Float_t sat_dec[],TH2D &evDist_border,std::vector<Double_t> &peacks_costheta,std::vector<Double_t> &peacks_phi,std::vector<Double_t> &valley_costheta,std::vector<Double_t> &valley_phi,std::vector<Double_t> &Gpeack_l,std::vector<Double_t> &Gpeack_b,std::vector<Double_t> &Gvalley_l,std::vector<Double_t> &Gvalley_b,Bool_t &first_call,Int_t &n_peacks,Int_t &n_valleys);
// extern void calculate_POI_stuff(Float_t sat_ra[],Float_t sat_dec[],TH2D* evDist,std::vector<Double_t> &peacks_costheta,std::vector<Double_t> &peacks_phi,std::vector<Double_t> &valley_costheta,std::vector<Double_t> &valley_phi);
extern void get_evDist_border(TH2D &evDist,TH2D &evDist_border);
extern void get_relevant_evDist_points(TH2D &evDist_border,std::vector<Double_t> &peacks_costheta,std::vector<Double_t> &valley_costheta,std::vector<Double_t> &peacks_phi,std::vector<Double_t> &valley_phi);
extern Int_t get_bunch_dimension(TH2D &evDist_border);
extern void GetPOI(std::vector<Double_t> &peacks_costheta,std::vector<Double_t> &peacks_phi,std::vector<Double_t> &valley_costheta,std::vector<Double_t> &valley_phi,std::vector<Double_t> &Gpeack_l,std::vector<Double_t> &Gpeack_b,std::vector<Double_t> &Gvalley_l,std::vector<Double_t> &Gvalley_b,Float_t sat_ra[],Float_t sat_dec[],Bool_t &good_PR,Int_t n_peacks,Int_t n_valleys);
extern void from_satellite_to_celestial(Float_t ra[],Float_t dec[],double vectorin[],AtPolarVect &vector_out);
extern void AtVect_To_AtPolarVect(double in_vector[],AtPolarVect &vector_out);
extern void invert_AtPolarVect_direction(AtPolarVect vector_out,AtPolarVect &vector_out_inv);
extern void AtPolarVect_to_vector(AtPolarVect &input_polar,double out_array[]);
extern void from_celestial_to_galactic(Double_t ra,Double_t dec,Double_t &l,Double_t &b);
extern void from_local_to_galactic(Double_t costheta,Double_t phi,Double_t &l,Double_t &b,Float_t sat_ra[],Float_t sat_dec[]);
extern Double_t obtain_min_histo(TH1D *histo);
extern bool points_are_inside(std::vector<Double_t> &Gpeack_l,std::vector<Double_t> &Gvalley_l,Int_t n_peacks,Int_t n_valleys);
extern void obtain_galactic_max_min(std::vector<Double_t> &Gpeack_l,std::vector<Double_t> &Gpeack_b,std::vector<Double_t> &Gvalley_l,std::vector<Double_t> &Gvalley_b,Int_t n_peacks,Int_t n_valleys,bool points_inside,Double_t &max_b,Double_t &min_b,Double_t &max_l,Double_t &min_l,Double_t &Nmax_b,Double_t &Nmin_b,Double_t &Pmax_b,Double_t &Pmin_b);

/*
extern void project_acceptance(TFile *results,Int_t h_idx,Float_t sat_ra[],Float_t sat_dec[],TH2D* acc,TH1D *acc_costheta,TH1D* acc_phi,ofstream &out_file,Int_t n_bin_lon,Double_t lon_bin_min,Double_t lon_bin_max,Int_t n_bin_lat,Double_t* &binning);
extern void obtain_edge_coordinate(Int_t tree_idx,Float_t sat_ra[],Float_t sat_dec[],TH2D *h_acc);
extern void project_acceptance_edgeMC(TFile *results,Int_t h_idx,Float_t sat_ra[],Float_t sat_dec[],TH2D* acc,TH1D *acc_costheta,TH1D* acc_phi,ofstream &out_file,Int_t n_bin_lon,Double_t lon_bin_min,Double_t lon_bin_max,Int_t n_bin_lat,Double_t* &binning);
extern void obtain_edge(Float_t sat_ra[],Float_t sat_dec[],TH2D *h_acc,Double_t &acc_glat_min,Double_t &acc_glat_max,Double_t &acc_glon_min,Double_t &acc_glon_max);
extern void acceptance_projection_study(TFile *results_file,Int_t tree_idx,Float_t sat_ra[],Float_t sat_dec[],TH2D *acc,TH1D *acc_costheta,ofstream &output_log_file,Int_t n_bin_lon,Double_t lon_bin_min,Double_t lon_bin_max,Int_t n_bin_lat,Double_t* &binning);
extern void shift_satellite_position(Float_t new_dec,Float_t new_ra,Float_t sat_ra[],Float_t sat_dec[],TH2D* acc,Int_t n_bin_lon,Double_t lon_bin_min,Double_t lon_bin_max,Int_t n_bin_lat,Double_t* &binning);
extern Bool_t check_Us(Float_t sat_ra[],Float_t sat_dec[],const Double_t &old_ang_xy,const Double_t &old_ang_xz,const Double_t &old_ang_yz);
extern void evaluate_Us(Double_t &old_ang_xy,Double_t &old_ang_xz,Double_t &old_ang_yz,Float_t sat_ra[],Float_t sat_dec[]);
extern void obtain_full_ra_dec_arrays(Float_t sat_ra[],Float_t sat_dec[]);
extern void obtain_EvDistr_Edge(vector<Float_t> &peacks_glon,vector<Float_t> &peacks_glat,vector<Float_t> &valley_glon,vector<Float_t> &valley_glat,TH2D* h_POI_healpix,TH2D* h_POI_lb_uniform,Bool_t points_inside);
*/
