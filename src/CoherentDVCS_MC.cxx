//
//  CoherentDVCS_MC.cxx
//  
//
//  Created by Hyde, Charles E. on 18-June-2021.
//
/** @file src/CoherentDVCS_MC.cxx
 *  @brief
 *   \f$l {}^{A}Z \to l  {}^{A}Z \gamma \f$
 *   lepto-production cross sections (electron or muon)
 *
 *  Input beam kinematics allow fixed target, head-on collisions, and crossing-angle collisions.
 *
 *  Created by Hyde, Charles E. on 18-June-2021.
 *
*/

/** @file src/CoherentDVCS_MC.hpp
 *  @brief File defines global variables (constants, and event-by-event variables).
 */
#include "CoherentDVCS_MC.hpp"

/** @file src/include/Deep_event.cxx
 * @brief File generates a set of event four-vectors
 */
#include "include/Deep_event.cxx"

#include <iostream>

int Init_Kin(char *inFILE);
double GammaSmear(TLorentzVector p4gamma, double *vec4out);

/** @brief
 *  Monte-Carlo Driver
 */
int CoherentDVCS_MC(const int nEvents = 1, const double tmin=-1.0){
    char FILEin[64] = "../input/He4-DVCS.inp";
    if (Init_Kin(FILEin)!=1){
        printf("Invalid initialization \n");
        return -1;
    }
    gSystem->Load("libRflxHepMCdict");
    TH2D *h_Q2_vs_xBj = new TH2D("h_Q2_vs_xBj",
                                 "Q2 vs x_{A}; x_{A}   ;  Q^{2} (GeV^{2})   ",
                                 20,0.0,1.0, 20,0.0,Q2Max);
    const int nbin_xA = 20;
    double xA_fact = pow(10.,-1./5.);
    double xA_BinLowEdge[nbin_xA+1];
    double xxA=1.0;
    for (int ibin=0; ibin<=nbin_xA; ibin++) {
        xA_BinLowEdge[nbin_xA-ibin]=xxA;
        xxA *= xA_fact;
    }
    const int nbin_Q2 = 10;
    double Q2_BinLowEdge[nbin_Q2+1];
    double Q2_fact = pow(Q2Max/Q2Min,1.0/(double)nbin_Q2);
    double Q2val = Q2Min;
    for (int ibin=0; ibin<=nbin_Q2+1; ibin++) {
        Q2_BinLowEdge[ibin] = Q2val;
        Q2val*=Q2_fact;
    }
    h_Q2_vs_xBj->SetBins(nbin_xA,xA_BinLowEdge,nbin_Q2,Q2_BinLowEdge);
    
    TH2D *h_Px_vs_xL_Ion = new TH2D("h_Px_vs_xL_Ion",
                                       "P_{T}(A) vs x_{L} (Ion frame)",
                                       100.,-1.0,1.0,100.,0.0, sqrt(-tmin));
    TH1D *h_DeltaSq = new TH1D("h_DeltaSq","#Delta^{2} (GeV^{2})",
                               1000.,tmin,0.0);
    TH2D *h_PT_vs_t = new TH2D("h_PT_vs_t",
                               "P_{T}(A) vs #Delta^{2}; #Delta^{2} (GeV^{2}); P_{T} (GeV) ",
                               100,tmin,0.0,100.,0.0,sqrt(-tmin));
    TH2D *h_qPr_vs_eta = new TH2D("h_qPr_vs_eta",
                                  "E_{#gamma} vs #eta_{#gamma}; #eta   ; GeV    ", 40,-4.0,4.0,
                                  40,0.0,2.*k4Beam0.E());
    TH2D *h_Q2_vs_y = new TH2D("h_Q2_vs_y",
                               "Q^{2} vs q#cdot P_{A}/k#cdot P_{A}; y_{A}    ;    Q^{2} (GeV^{2})  ",
                               50, 0.0,1.0, nbin_Q2,Q2_BinLowEdge);
    TH1D *h_pT_err = new TH1D("h_pT_err",
                              "P_{T} error from #gamma resolution; #Delta P_{T} (GeV)",
                              100.,0.0,0.100);
    TH1D *h_sigma = new TH1D("h_sigma",
                              "#gamma resolution #sigma_{E}/{E}",
                              100.,0.0,0.100);

    double xL;
    TVector3 P3Scat, P3Beam, P3ScatPerp;
    P3Beam = P4Beam0.Vect();
    TLorentzVector P4perp0, P4perp1, q4Smear;
    double g4smear[4];
    double pTerr, sigma;
    // call event generation loop
    int ievt = 0, mevts=0;
    for (ievt=0; ievt<nEvents; ievt++) {
        if(Get_Event(ievt)>=0)
        {
            if ((tmin<Delta4vec.M2())&&(k4Scat.E()>0))
            {
                h_Q2_vs_xBj->Fill(xBj,Q2);
                P3Scat = P4Scat.Vect();
                xL = P3Scat.Dot(P3Beam)/P3Beam.Mag2();
                P3ScatPerp = P3Beam.Cross(P3Scat);
                h_Px_vs_xL_Ion->Fill(xL,P3ScatPerp.Mag()/P3Beam.Mag());
                h_DeltaSq->Fill(Delta4vec.M2());
                h_PT_vs_t->Fill(Delta4vec.M2(),P3ScatPerp.Mag()/P3Beam.Mag());
                h_qPr_vs_eta->Fill(q4Prime.Eta(),q4Prime.E());
                h_Q2_vs_y->Fill(yInv,Q2);
                sigma = GammaSmear(q4Prime, g4smear);
                if (sigma<0.0) continue;
                P4perp0 = k4Beam0+P4Beam0-k4Scat;
                P4perp1 = k4Beam+P4Beam-k4Scat;
                P4perp0 += (-1.0)*q4Prime;
                q4Smear.SetPxPyPzE(g4smear[1],g4smear[2],g4smear[3],g4smear[0]);
                P4perp1 += (-1.0)*q4Smear;
                P4perp0 = (-P4perp0.Dot(X4_q))*X4_q - (P4perp0.Dot(Y4_q))*Y4_q;
                P4perp1 = (-P4perp1.Dot(X4_q))*X4_q - (P4perp1.Dot(Y4_q))*Y4_q;
                P4perp1 += (-1.0)*P4perp0;
                pTerr = -P4perp1.M2();
                if (pTerr>=0.0) {
                    h_pT_err->Fill(sqrt(pTerr));
                    h_sigma->Fill(sigma);
                } else {
                    printf("-P4perp1.M2() = %10.3g \n", pTerr);
                    continue;
                }
            }
            mevts++;
        }
    }
    TCanvas *c1 = new TCanvas("c1","DVCS Kinematics",100,50,
                              900,600);
    c1->Divide(3,2);
    c1->cd(1);
    h_Q2_vs_xBj->Draw();
    gPad->SetLogy();
    gPad->SetLogx();
    c1->cd(2);
    h_Px_vs_xL_Ion->Draw();
    c1->cd(3);
    h_PT_vs_t->Draw();
    c1->cd(4);
    h_Q2_vs_y->Draw();
    gPad->SetLogy();
    c1->cd(5);
    h_pT_err->Draw();
    TF1 *f1 = new TF1("f1","x*gaus",0.0,0.2);
    f1->FixParameter(1,0.0);
//    h_pT_err->Fit(f1,"I");
//    h_DeltaSq->Draw();
    h_sigma->SetLineColor(kRed);
    h_sigma->Draw("same");
    c1->cd(6);
    h_qPr_vs_eta->Draw();
    
    return mevts;
} // CoherentDVCS_MC()

/** @brief
 *  Read input file to initialize Monte-Carlo event generation.
 */
int Init_Kin(char *inFile){
    Y4_Det.SetPxPyPzE(0.0,1.0,0.0,0.0);
    FILE  * inF;
    printf("Input control parameters from file %s \n",inFile);
    inF=fopen(inFile,"r");
    if (inF==NULL) {
        printf("Error opening input file %s \n",inFile);
        return -1;
    }
    const int nchar=132;
    double vx,vy,vz, ex,ey,ep, bx,by ;
    double ELepton, EIon;
    char line[nchar], lepton[15];
    mPion = dbPDG->GetParticle(211)->Mass();
    //  Read header line
    fgets(line,nchar,inF); printf("%s",line);
    //
    //  Read lepton header line
    fgets(line,nchar,inF); printf("%s",line);
    //  Read lepton beam descriptor
    fgets(line,nchar,inF); printf("%s",line);
    sscanf(line,"%lf %lf %lf %lf %lf %lf %lf %lf %s", &vx, &vy, &vz, &ex, &ey, &ep, &bx, &by, lepton);
    if (strncmp(lepton,"electron",7)==0) {
        mLepton = dbPDG->GetParticle(11)->Mass();
    } else if (strncmp(lepton,"muon",4)==0){
        mLepton = dbPDG->GetParticle(13)->Mass();
    } else {
        printf("invalid lepton = %s \n", lepton);
        return -1;
    }
    long double ELSq = (vx*vx) + (vy*vy) + (vz*vz) + mLepton*mLepton;
    ELepton = sqrt(ELSq);
    printf("(vx,vy,vz) = (%10.4f,%10.4f,%10.4f), mLepton = %10.3e, ELepton = %13.8lf, ELSq=%13.8Lf \n",
           vx,vy,vz, mLepton, ELepton, ELSq);
    k4Beam0.SetPxPyPzE(vx,vy,vz,ELepton);
    emitt_e[0] = ex*mLepton/ELepton; emitt_e[1] = ey*mLepton/ELepton; emitt_e[2] = ep;
    betaIP_e[0] = bx;                betaIP_e[1] = by;
    if (ex+ey+ep <= 0.0){
        eSmear = false;
    }
    //
    //  Read Ion header line
    fgets(line,nchar,inF); printf("%s",line);
    //  Read Ion Beam descripton
    fgets(line,nchar,inF); printf("%s",line);
    sscanf(line,"%lf %lf %lf %lf %lf %lf %lf %lf %d %d", &vx, &vy, &vz, &ex, &ey, &ep, &bx, &by, &ionZ, &ionA);
    MIon  = 0.0;
    // Nucleide masses from nndc.bnl.gov
    switch (ionZ) {
        case 0:
            if (ionA==1){
            // neutron
                MIon = dbPDG->GetParticle(2112)->Mass();
            }
            break;
        case 1:
            switch (ionA) {
                case 1:
                    // proton
                    MIon = dbPDG->GetParticle(2212)->Mass();
                    break;
                case 2:
                    MIon = ionA*umass+13.1357e-3; // deuteron
                    break;
                case 3:
                    MIon = ionA*umass+14.9498e-3; // triton
                default:
                    MIon  =ionA*umass;
                    break;
            }
        case 2:
            switch (ionA) {
                case 3:
                    MIon  = ionA*umass + 14.9312e-3; // 3He
                    break;
                case 4:
                    MIon  = ionA*umass + 2.4249e-3; //4He
                    break;
                default:
                    MIon  = ionA*umass;
                    break;
            }
        case 6:
            switch (ionA) {
                case 12:
                    MIon  = ionA*umass + 0.0; //12C
                    break;
                default:
                    MIon  = ionA*umass;
                    break;
            }
        default:
            MIon  = ionA*umass;
            break;
    }
    printf(" Isotope Z=%3d, A=%3d, Mass = %8.5f GeV \n",
           ionZ, ionA, MIon );
    EIon = sqrt(MIon*MIon + vx*vx + vy*vy + vz*vz);
    P4Beam0.SetPxPyPzE(vx,vy,vz, EIon);
    //  Calculate Geometrical (Transverse) Emittance
    emitt_i[0] = ex*MIon/EIon; emitt_i[1]=ey*MIon/EIon; emitt_i[2] = ep;
    betaIP_i[0] = bx;                betaIP_i[1] = by;
    if (ex+ey+ep <= 0.0){
        iSmear = false;
    }
    printf("LeptonMassSq, k4Beam0.M2(), P4Beam0.M2() = %13.6e, %13.6e, %10.6f \n",
           mLepton*mLepton, k4Beam0.M2(), P4Beam0.M2());
    
    //  Read Kinematics header line
    fgets(line,nchar,inF); printf("%s",line);
    //  Read MC Kinematics bounds
    fgets(line,nchar,inF); printf("%s",line);
    sscanf(line,"%lf %lf %lf %lf", &vx, &vy, &ex, &ey);
    Q2Min = vx;
    Q2Max = vy;
    yMin  = ex;
    yMax  = ey;
    //  Hadronic Final state bounds
    fgets(line,nchar,inF); printf("%s",line);
    fgets(line,nchar,inF); printf("%s",line);
    sscanf(line,"%lf %lf", &ex, &ey);
    csCMmin = ex;
    csCMmax = ey;
    //  Setup up nominal beam light-cone vectors and transverse coordinates
    k4Beam = k4Beam0;
    P4Beam = P4Beam0;
//    Y4_Det.SetPxPyPzE(0.0,1.0,0.0,0.0);
    EventLightCone();
    n4_e0      = n4_e;
    n4Tilde_e0 = n4Tilde_e;
    X4_e0      = X4_e;
    Y4_e0      = Y4_e;
    W2Threshold= pow(MIon+ mPion,2);
    fclose(inF);
    return 1;
}   //Init()



