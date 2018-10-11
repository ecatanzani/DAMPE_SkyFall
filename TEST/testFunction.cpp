/*
        TEST FUNCTION
 
    This function is a test of how to read and write histos from file using external functions and not genefrating memory leacks.
        - Histos, in this case, have to be created into ROOT environment, before running the "Main" fucntion "GreenEgg()".
        - One histo il filled
        - The histo is then saved to file.
        - File is then opened again and the histo is obtained.

*/

//////////////// Include headers

#include <iostream>
#include <stdio.h>
#include <stdlib.h>

#include "TFile.h"
#include "TH1.h"
#include "TRandom3.h"
#include "TCanvas.h"

///////////////////////////////////////


void uniform_distribution_histo(TH1D& histo, size_t samples = 1000, double min_range = 0.0, double max_range = 1.0)
{
    TRandom3 rgen;
    for(size_t i=0; i != samples; ++i)  histo.Fill(rgen.Uniform(min_range,max_range));
}

void read_histo(TH1D& histo, const std::string& name, const std::string& filename = "do_not_open.root")
{
    TFile infile(filename.c_str(),"READ");
    //call the constructor of TH1D
    new (&histo) TH1D(*(TH1D*)infile.Get(name.c_str()));
    //..
    infile.Close();
}

void save_histo(TH1D& histo, const std::string& filename = "do_not_open.root")
{
    TFile outfile(filename.c_str(),"RECREATE");
    histo.Write();
    outfile.Close();
}

void GreenEgg(TH1D& histo1, TH1D& histo2)
{
    histo1 = TH1D("histo","histo",1000,0,1);
    uniform_distribution_histo(histo1);
    save_histo(histo1);
    
    read_histo(histo2,"histo");
}




