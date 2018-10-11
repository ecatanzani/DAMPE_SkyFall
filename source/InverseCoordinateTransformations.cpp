#include "MyHead.h"

///////////////////////////////////////////// Math custom functions

static Double_t acot(Double_t value) { return atan(1/value); }

/////////////////////////////////////////////////////////////////////////

void from_galactic_to_celestial(Double_t &ra,Double_t &dec,Double_t l,Double_t b) {
    
    Double_t ragc=192.85948,decgc=27.12825,lcp=122.932;
    
    Double_t ragcr=ragc*TMath::DegToRad();
    Double_t decgcr=decgc*TMath::DegToRad();
    Double_t lcpr=lcp*TMath::DegToRad();
    Double_t decr,rar;
    
    Double_t old_ra=ra;
    Double_t tmp_l=0,tmp_b=0;
    
    Double_t br,lr,sin_t,cos_t,t;
    
    ////////////////////////////////////////////////
    
    br = b/TMath::RadToDeg();
    lr = l/TMath::RadToDeg();
    
    t = lcpr-lr;
    sin_t = sin(t);
    cos_t = cos(t);
    
    //////////// Constants to easily invert the map
    
    Double_t c1 = sin(br);
    Double_t c2 = sin(decgcr);
    Double_t c3 = cos(decgcr);
    Double_t c4 = sin_t*cos(br);
    Double_t c5 = cos_t*cos(br);
    Double_t c6 = ragcr;
    
    ////////////////////////////////////////////////
    
    decr = asin(c1/c2-((c3*c4)/c2)*((c3*c1-c5*c2)/(TMath::Power(c3,2)*c4+TMath::Power(c2,2)*c4)));
    rar = acot((c3*c1-c5*c2)/(TMath::Power(c3,2)*c4+TMath::Power(c2,2)*c4))+c6;
    
    ///////////////////////////////////////////////
    
    ra = rar/TMath::DegToRad();
    dec = decr/TMath::DegToRad();
    
    //////////////////////////////// BLACK MAGIC TO SOLVE 180 DEGREES BIAS !!!
    
    from_celestial_to_galactic(ra,dec,tmp_l,tmp_b);
    if(fabs(tmp_l-l)>err_inv || fabs(tmp_b-b)>err_inv)
        ra+=180;
    from_celestial_to_galactic(ra,dec,tmp_l,tmp_b);
    if(fabs(tmp_l-l)>err_inv || fabs(tmp_b-b)>err_inv)
        ra=old_ra-180;
    
    while(ra>360)
        ra-=360;
    
    /////////////////////////////////////////////////////////////////////////////////////////
    
}

void from_celestial_to_local(AtPolarVect vector_out,Double_t vector_in[]) {
    Double_t abs_s,s;
    Double_t c,norm01;
    
    norm01 = TMath::Power(vector_out.r,2)-TMath::Power(vector_out.r*sin(vector_out.lat),2);
    c=(1-TMath::Power(tan(vector_out.lon/2.),2))/(1+TMath::Power(tan(vector_out.lon/2.),2));
    abs_s = sqrt(1-TMath::Power(c,2));
    s=((1-c)/tan(vector_out.lon/2.));
    
    if(abs_s>EPS) {
        vector_in[0]=c*sqrt(norm01);
        vector_in[1]=s*sqrt(norm01);
        vector_in[2]=vector_out.r*sin(vector_out.lat);
    }
    else {
        c=1;
        vector_in[0]=c*sqrt(norm01);
        vector_in[1]=0;
        vector_in[2]=vector_out.r*sin(vector_out.lat);
    }
    
}

void obtain_costheta_phi(Double_t &costheta,Double_t &phi,Float_t sat_ra[],Float_t sat_dec[],Double_t vector_in[]) {
    Float_t ux1[3];
    Float_t uy1[3];
    Float_t uz1[3];
    Float_t rax = sat_ra[0];
    Float_t ray = sat_ra[1];
    Float_t raz = sat_ra[2];
    Float_t decx = sat_dec[0];
    Float_t decy = sat_dec[1];
    Float_t decz = sat_dec[2];
    
    Double_t tmp_local[3],sinphi,cosphi;
    
    ux1[0] = cos(decx)*cos(rax);
    ux1[1] = cos(decx)*sin(rax);
    ux1[2] = sin(decx);
    
    uy1[0] = cos(decy)*cos(ray);
    uy1[1] = cos(decy)*sin(ray);
    uy1[2] = sin(decy);
    
    uz1[0] = cos(decz)*cos(raz);
    uz1[1] = cos(decz)*sin(raz);
    uz1[2] = sin(decz);
    
    tmp_local[2]=((ux1[0]*vector_in[1]-ux1[1]*vector_in[0])*(ux1[2]*uy1[0]-ux1[0]*uy1[2])+(ux1[0]*vector_in[2]-ux1[2]*vector_in[0])*(ux1[0]*uy1[1]-ux1[1]*uy1[0]))/((ux1[0]*uz1[2]-ux1[2]*uz1[0])*(ux1[0]*uy1[1]-ux1[1]*uy1[0])-(ux1[1]*uz1[0]-ux1[0]*uz1[1])*(ux1[2]*uy1[0]-ux1[0]*uy1[2]));
    
    tmp_local[1]=(tmp_local[2]*(ux1[1]*uz1[0]-ux1[0]*uz1[1])+ux1[0]*vector_in[1]-ux1[1]*vector_in[0])/((ux1[0]*uy1[1])-(ux1[1]*uy1[0]));
    
    tmp_local[0]=(1/ux1[0])*(vector_in[0]-tmp_local[1]*uy1[0]-tmp_local[2]*uz1[0]);
    
    costheta=tmp_local[2];
    sinphi=tmp_local[1]/sin(acos(costheta));
    cosphi=tmp_local[0]/sin(acos(costheta));
    
    if(sinphi>=0 && cosphi>=0)
        phi=acos(cosphi);
    else if(sinphi>0 && cosphi<0)
        phi=acos(cosphi);
    else if(sinphi<0 && cosphi<0) {
        phi=acos(cosphi);
        phi+=2*(TMath::Pi()-phi);
    }
    else {
        phi=acos(cosphi);
        phi=2*TMath::Pi()-phi;
    }
}

void obtain_costheta_phi_ROOTf(Double_t &costheta,Double_t &phi,Float_t sat_ra[],Float_t sat_dec[],Double_t vector_in[]) {
    Float_t ux1[3];
    Float_t uy1[3];
    Float_t uz1[3];
    Float_t rax = sat_ra[0];
    Float_t ray = sat_ra[1];
    Float_t raz = sat_ra[2];
    Float_t decx = sat_dec[0];
    Float_t decy = sat_dec[1];
    Float_t decz = sat_dec[2];
    
    Double_t tmp_local[3],sinphi,cosphi,tmp_sum;
    TMatrixD Us(3,3),Us_invert(3,3);
    
    ux1[0] = cos(decx)*cos(rax);
    ux1[1] = cos(decx)*sin(rax);
    ux1[2] = sin(decx);
    
    uy1[0] = cos(decy)*cos(ray);
    uy1[1] = cos(decy)*sin(ray);
    uy1[2] = sin(decy);
    
    uz1[0] = cos(decz)*cos(raz);
    uz1[1] = cos(decz)*sin(raz);
    uz1[2] = sin(decz);
    
    //Fill matrix with Us already calculated
    
    Us(0,0)=ux1[0];
    Us(0,1)=ux1[1];
    Us(0,2)=ux1[2];
    
    Us(1,0)=uy1[0];
    Us(1,1)=uy1[1];
    Us(1,2)=uy1[2];
    
    Us(2,0)=uz1[0];
    Us(2,1)=uz1[1];
    Us(2,2)=uz1[2];
    
    ////// Invert Us
    
    Us_invert=Us.Invert();
    
    ////// Obtain tmp_local array
    
    for(Int_t idx_c=0; idx_c<3; idx_c++) {
        tmp_sum=0;
        for(Int_t idx_r=0; idx_r<3; idx_r++)
            tmp_sum+=vector_in[idx_r]*Us_invert(idx_r,idx_c);
        tmp_local[idx_c]=tmp_sum;
    }
    
    ///// Obtain costheta and phi in local reference frame !
    
    costheta=tmp_local[2];
    sinphi=tmp_local[1]/sin(acos(costheta));
    cosphi=tmp_local[0]/sin(acos(costheta));
    
    if(sinphi>=0 && cosphi>=0)
        phi=acos(cosphi);
    else if(sinphi>0 && cosphi<0)
        phi=acos(cosphi);
    else if(sinphi<0 && cosphi<0) {
        phi=acos(cosphi);
        phi+=2*(TMath::Pi()-phi);
    }
    else {
        phi=acos(cosphi);
        phi=2*TMath::Pi()-phi;
    }
}

Bool_t outside_evDist(Double_t costheta,Double_t phi,TH2D* evDist_border) {
    Bool_t outside;
    Int_t Ybin;
    Double_t evDist_costheta=0,lXbin=0;
    
    TAxis *Yaxis = evDist_border->GetYaxis();
    TAxis *Xaxis = evDist_border->GetXaxis();
    Ybin=Yaxis->FindBin(phi);
    lXbin=.5/Xaxis->GetNbins();
    
    for(Int_t x_bin=1; x_bin<=evDist_border->GetNbinsX(); x_bin++) {
        if(evDist_border->GetBinContent(x_bin,Ybin)!=0) {
            //evDist_costheta=.5+x_bin*(.5/Xaxis->GetNbins());
            evDist_costheta = .5 + .5*((x_bin+1)*lXbin+x_bin*lXbin);
            break;
        }
    }
    
    if(evDist_costheta>costheta)
        outside=true;
    else
        outside=false;
    
    return outside;
}
