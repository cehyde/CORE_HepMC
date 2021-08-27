
// dvcs_analyzer.cpp
// Charles E. Hyde, 10 July 2021
//
// Analyze the output of a Delphes simulation of DVCS event file in CORE

#include <stdio.h>

#include <TTree.h>
#include <TFile.h>
#include <TSystem.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TH2D.h>
#include <TChain.h>
#include <TF1.h>
#include <TLorentzVector.h>
#include <TRandom3.h>
#include <TStyle.h>

#ifdef __CLING__
R__LOAD_LIBRARY(libDelphes)
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootProgressBar.h"
#include "external/ExRootAnalysis/ExRootTreeBranch.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#endif


const double Crossing_angle = -0.025;  // switch future input files to positive

int dvcs_analyzer(const char *fname="output_CORE"){
    char inputfile[128];
    sprintf(inputfile,"output/%s.root",fname);
    printf("%s\n",inputfile);
    gSystem->Load("libDelphes");
    TRandom3 ran3;

    // Create chain of root trees
    TChain chain("Delphes");
    chain.Add(inputfile);

    
    // Create object of class ExRootTreeReader
    ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
    Long64_t numberOfEntries = treeReader->GetEntries();
    
    // Get pointers to branches used in this analysis
    TClonesArray *branchJet = treeReader->UseBranch("Jet");
    TClonesArray *branchElectron = treeReader->UseBranch("Electron");
    TClonesArray *branchParticle = treeReader->UseBranch("Particle");
    TClonesArray *branchPhoton = treeReader->UseBranch("Photon");
    TClonesArray *branchEFlow = treeReader->UseBranch("EFlowTrack");
    TClonesArray *branchTrack = treeReader->UseBranch("Track");

    GenParticle * thisParticle;
    Photon * thisPhoton;
    Electron * thisElectron;
    Track * thisEFlow;
    Track * thisTrack;
    // Book histograms
//    TH1 *histJetPT = new TH1F("jet_pt", "jet P_{T}", 100, 0.0, 100.0);
    TH1 *histMass = new TH1F("mass", "M_{inv}(e_{1}, e_{2})", 100, 0.0, 200.0);
    TH1 *h_MissingMassSq = new TH1F("h_MissingMassSq","M_{X}^{2} (GeV2)",100.,-5.0,5.0);
    TH1 *h_photons = new TH1F("h_photons","Photon Multiplicity",10,-0.5,9.5);
    TH1 *h_electrons = new TH1F("h_electrons","Electron Multiplicity",10,-0.5,9.5);
    TH1 *h_qprime = new TH1F("h_qprime","(Q2-Q2_{Gen})/Q2_{Gen}",250,-1.0,1.0);
    TH1 *h_PT = new TH1F("h_PT","Photon PT ; GeV ; Events", 200,0.0,20.);
    TH1 *h_E = new TH1F("h_E","Photon Energy ; GeV ; Events", 200,0.0,20.);
    TH1 *h_D_err = new TH1F("h_D_err","#Delta^{2}-#Delta^{2}_{Gen}; GeV^{2}  ; Events   ", 200.,-2.0,2.0);
    TH1 *h_Delta2gen = new TH1F("h_Delta2gen","#Delta^{2}_{Gen};GeV^{2}   ; ",
                               110,-1.0,0.10);
    TH1 *h_Delta2 = new TH1F("h_Delta2","#Delta^{2}_{Gen};GeV^{2}   ; ",
                               110,-1.0,0.10);

    TH1 *h_D_errPos = new TH1F("h_D_errPos","#Delta^{2}-#Delta^{2}_{Gen}; GeV^{2}  ; Events   ",
                            200.,-2.0,2.0);
    TH2 *h2_Q2err_y = new TH2F("h2_Q2err_y","(Q2-Q2_{Gen})/Q2_{Gen} vs y; y=(qP)(kP)  ; Q2/Q2_{Gen}-1   ",
                               50,0.0,1.0, 50,-1.0,1.0);
    TH2 *h2_Q2err_eta = new TH2F("h2_Q2err_eta","(Q2-Q2_{Gen})/Q2_{Gen} vs #eta; #eta_{e}  ; Q2/Q2_{Gen}-1   ",
                                 50,-4.0,1.0, 50,-1.0,1.0);
    char htitle[128]="#Delta^{2}-#Delta^{2}_{Gen} vs #Delta^{2}_{Gen}; #Delta^{2} GeV^{2}  ;#Delta#Delta^{2} GeV^{2}  ";
    TH2 *h2_errD_D = new TH2F("h_errD_D",htitle,50,-1.0,0.0,100,-2.0,2.0);
    TH2 *h2_qPrDiff_EvsP = new TH2F("h2_qPrDiff_EvsP","q'-q'_{Gen}; 3-vector mag; EDiff/E_{Gen}",
                                    50,0.0,1.0, 50,-0.2,0.2);
    char dtitle[128]="#Delta^{2}-#Delta^{2}_{Gen} vs q'.E-q'Gen.E; EDiff   ; #Delta#Delta^{2} GeV^{2}  ";
    TH2 *h2_errD_dE = new TH2F("h2_errD_dE",dtitle, 50,-1.0,1.0, 100,-2.0,2.0);
    TH1 *h_etaTarget = new TH1F("h_etaTarget","Recoil Ion; #eta(Ion)  ", 200,-10.,10.);
    TH1 *h_tIon = new TH1F("h_tIon","Detected vs Generated Ion; #Delta^{2}-#Delta^{2}_{Gen} (GeV2) ",
                           200,-1.0,1.0);
    TH2 *h2_P4Gen = new TH2F("h2_P4Gen","Scattered Ion; Pz ; Px ",
                             50, 100.,300., 50, 0.0,10.0);
    TH1 *h_dqPrime= new TH1F("h_dqPrime","(E(#gamma')-E(#gamma')_{Gen})/E(#gamma')_{Gen}",
                             100,-0.5,0.5);
    TH1 *h_DDelta2_naive = new TH1F("h_DDelta2_naive","Two Percent Solution; #Delta^{2}-#Delta^{2}_{Gen} (GeV2)",
                                    200,-2.0,2.0);
    TH1 *h_EMissInv = new TH1F("h_EMissInv","P4Miss.Dot(P4Beam)/P4Beam.M; GeV ",
                               200,-2.0,10.0);
    TH1 *h_MX2 = new TH1F("h_MX2","P4Miss.M2()-P4Beam.M2(); GeV^{2}  ", 200,-4.0,4.0);
    TH1 *h_MX2_0 = new TH1F("h_MX2_0","[P4Miss.M2()-P4Beam.M2()]_{Gen}; GeV^{2}  ",
                             200,-4.0,4.0);
    TH1 *h_Zero = new TH1F("h_Zero","(k+P-q'-k'-P)^{2};  GeV^{2}   ",
                           200, -4.0,4.0);
    TH1 *h_DeltaPerp = new TH1F("h_DeltaPerp","#Delta_{#perp}^{2}; GeV^{2}",
                                200,0.0,1.0);
    TH1 *h_norm = new TH1F("h_norm","LightCone Mass",200,-4.0,4.0);
    TH1 *h_DeltaPerp_err = new TH1F("h_DeltaPerp_err","-(#Delta_{#perp} -#Delta_{#perp Gen})^{2}; GeV^{2}    ",
                                    500,-0.1,0.4);
    const int nbin_xA = 20;
    double xA_fact = pow(10.,-1./5.);
    double xA_BinLowEdge[nbin_xA+1];
    double xxA=1.0;
    printf("xBj low edge" );
    for (int ibin=0; ibin<=nbin_xA; ibin++) {
        printf("%8.3g, ", xxA );
        xA_BinLowEdge[nbin_xA-ibin] = xxA;
        xxA *= xA_fact;
    }
    printf(" \n");
    const double Q2Min = 1.0, Q2Max=200.0;
    const int nbin_Q2 = 10;
    double Q2_BinLowEdge[nbin_Q2+1];
    double Q2_fact = pow(Q2Max/Q2Min,1.0/(double)nbin_Q2);
    double Q2val = Q2Min;
    for (int ibin=0; ibin<=nbin_Q2+1; ibin++) {
        Q2_BinLowEdge[ibin] = Q2val;
        Q2val*=Q2_fact;
    }

    TH2 *h2_Q2_xBj = new TH2F("h2_Q2_xBj","Q^{2} vs x_{Bj};x_{Bj}    ;Q^{2}      ",
                              nbin_xA,xA_BinLowEdge,nbin_Q2,Q2_BinLowEdge);
    //h2_Q2_xBj->SetBins(nbin_xA,xA_BinLowEdge,nbin_Q2,Q2_BinLowEdge);

    TLorentzVector k4Beam, k4Scat, P4Beam, P4Scat, q4prime, Delta, P4Miss, k4ScatCheck;
    TLorentzVector k4BeamGen, k4ScatGen, k4EFlow, P4BeamGen, P4ScatGen, q4primeGen, P4MissGen;
    TLorentzVector Delta4Gen, Zero, q4Virtual, n4_q, n4Tilde_q, Delta_perp, Delta_perpGen;
    TVector3 q3prime,q3primeGen, q3Diff;
    double DeltaSqGen, DeltaSq, Q2, Q2Gen, kScat_E, yGen, relErr, radiation;
    double aEM=0.01, bEM=0.02, esmear,gsmear;
    double q_dot_P, deltaQ, sqrt_one_d, norm;
    
    // Loop over all events
    Int_t entry=0;
    Int_t nElectron=0;
    Int_t nJet=0;
    Int_t nPhoton=0;
    // Delphes status codes, incident beam or stable produced particle
    Int_t iBeamStatus = 4;
    Int_t iStableStatus = 1;
    // PDG codes
    const int PDGelectron = 11;
    const int PDGphoton = 22;
    const int PDGproton = 2212;
    const int PDGbaryons= 2000;
    for(entry = 0; entry < numberOfEntries; ++entry)
    {
        // Load selected branches with data from specified event
        treeReader->ReadEntry(entry);
        nElectron += branchElectron->GetEntries();
        h_electrons->Fill((float)branchElectron->GetEntries());
        nJet += branchJet->GetEntries();
        nPhoton += branchPhoton->GetEntries();
        h_photons->Fill((float)branchPhoton->GetEntries());
        // Select the highest energy electron as the scattered electron
        kScat_E = 0.0;
        for (int ipart=0; ipart<branchElectron->GetEntries(); ipart++) {
            thisElectron = (Electron *) branchElectron->At(ipart);
            if (thisElectron->P4().E()>kScat_E) {
                k4Scat = thisElectron->P4();
                kScat_E = k4Beam.E();
            }
        }
        if (branchPhoton->GetEntries()>1) continue;
        //  Loop on detected photons, but only consider photons above 1 GeV
        for (Int_t iphot=0; iphot<branchPhoton->GetEntries(); iphot++) {
            thisPhoton = (Photon *) branchPhoton->At(iphot);
            q4prime = thisPhoton->P4();
            h_PT->Fill(thisPhoton->PT);
            h_E->Fill(q4prime.E());
            radiation = 2.*k4Scat.Dot(q4prime);
            if(q4prime.E()>1.0&&branchParticle->GetEntries()>4&&radiation>1.0)
            {
            
                //  branchParticle holds the generated particles
                for (int jpart=0; jpart<branchParticle->GetEntries(); jpart++){
                    thisParticle = (GenParticle *) branchParticle->At(jpart);
                    if (thisParticle->Status==iBeamStatus) {
                        if (thisParticle->PID==PDGelectron) {
                            k4BeamGen = thisParticle->P4();
                        } else if (thisParticle->PID>PDGbaryons) {
                            P4BeamGen = thisParticle->P4();
                        }
                    } else if (thisParticle->Status==iStableStatus) {
                        if (thisParticle->PID==PDGelectron) {
                            k4ScatGen = thisParticle->P4();
                        } else if (thisParticle->PID>PDGbaryons) {
                            P4ScatGen = thisParticle->P4();
                            h2_P4Gen->Fill(P4ScatGen.Pz(),P4ScatGen.Px());
                        } else if (thisParticle->PID==PDGphoton) {
                            q4primeGen = thisParticle->P4();
                        }
                    }
                }
                yGen = 1.0 - k4ScatGen.Dot(P4BeamGen)/k4BeamGen.Dot(P4BeamGen);
                //  Simplistic smearing
                esmear = 1.0+ran3.Gaus(0.0,aEM)+ran3.Gaus(0.0,bEM)/k4ScatGen.E();
                gsmear = 1.0+ran3.Gaus(0.0,aEM)+ran3.Gaus(0.0,bEM)/q4prime.E();
                k4ScatCheck.SetVectM(esmear*k4Scat.Vect(),sqrt(k4Scat.M2()));
                Delta = k4BeamGen-k4ScatCheck-gsmear*q4prime;
                DeltaSq = Delta.M2();
                //q3prime = q4prime.Vect();
                //q3primeGen=q4primeGen.Vect();
                q3Diff = q4prime.Vect() - q4primeGen.Vect();
                h2_qPrDiff_EvsP->Fill(q3Diff.Mag(),q4prime.E()/q4primeGen.E()-1.0);
                P4Miss = k4BeamGen-k4ScatGen-q4primeGen; //P4BeamGen-P4ScatGen;
                Delta4Gen = P4Miss;
                DeltaSqGen = P4Miss.M2();
                h_DDelta2_naive->Fill(DeltaSq-DeltaSqGen);
                Zero = k4BeamGen+P4BeamGen-k4ScatGen-q4primeGen-P4ScatGen;
                //h_Zero->Fill(Zero.M2());
                relErr = k4BeamGen.Dot(k4Scat)/k4BeamGen.Dot(k4ScatGen)-1.0;
                if (k4Scat.Eta()<0.0 && k4Scat.E()>1.0 )
                {
                    Delta = k4BeamGen-k4Scat-q4prime;
                    DeltaSq = Delta.M2();
                    if (q4prime.Eta()<0.0)
                    {
                        h_D_err->Fill(DeltaSq-DeltaSqGen);
                        h2_errD_D->Fill(DeltaSqGen,DeltaSq-DeltaSqGen);
                        h2_errD_dE->Fill(q4prime.E()-q4primeGen.E(),DeltaSq-DeltaSqGen);
                    } else if (q4prime.Eta()>=0.0)
                    {
                        h_D_errPos->Fill(DeltaSq-DeltaSqGen);
                        continue;
                    }
                    h_Delta2gen->Fill(DeltaSqGen);
                    h_Delta2->Fill(DeltaSq);
                    thisEFlow = (Track *) branchEFlow->At(0);
                    k4EFlow = thisEFlow->P4();
                    h_qprime->Fill(relErr);
                    // k4BeamGen.Dot(k4Scat)/k4BeamGen.Dot(k4ScatGen)-1.0);
                    Q2val = 2.*k4BeamGen.Dot(k4ScatGen);
                    xxA   = Q2val/(-2.0*(P4BeamGen.Dot(k4Scat)-P4BeamGen.Dot(k4BeamGen)));
                    h2_Q2_xBj->Fill(xxA,Q2val);
                    h2_Q2err_y->Fill(yGen,k4BeamGen.Dot(k4Scat)/k4BeamGen.Dot(k4ScatGen)-1.0);
                    h2_Q2err_eta->Fill(k4ScatGen.Eta(),
                                       k4BeamGen.Dot(k4Scat)/k4BeamGen.Dot(k4ScatGen)-1.0);
                    h_dqPrime->Fill(q4prime.E()/q4primeGen.E()-1.0);
                    //  Lorentz Invariant "Missing Energy"
                    P4Miss = P4BeamGen+Delta;
                    h_EMissInv->Fill(P4Miss.Dot(P4BeamGen)/sqrt(P4BeamGen.M2()));
                    h_MX2->Fill(P4Miss.M2()-P4BeamGen.M2());
                    h_MX2_0->Fill(2.0*P4BeamGen.Dot(Delta4Gen)+DeltaSqGen);
                    // Find Lorentz Invariant Delta_perp
                    // Detected recoil nucleus?
                    //
                    //  Make the Event Light cone vectors n4_q and n4Tilde_q
                    //
                    q4Virtual = k4BeamGen-k4Scat;
                    q_dot_P   = q4Virtual.Dot(P4BeamGen);
                    deltaQ    = -q4Virtual.M2()*P4BeamGen.M2()/(q_dot_P*q_dot_P);
                    sqrt_one_d= sqrt(1.0+deltaQ);
                    norm      = 1.0/(sqrt(q_dot_P)*sqrt_one_d);
                    n4_q  = (1.0+sqrt_one_d)* q4Virtual;
                    n4_q += ( -q4Virtual.M2()/q_dot_P ) * P4BeamGen;
                    n4_q *= norm/2.0;
                    n4Tilde_q  = (-P4BeamGen.M2()/((1.0+sqrt_one_d)*q_dot_P)) *q4Virtual ;
                    n4Tilde_q += P4BeamGen;
                    n4Tilde_q *= norm;
                    if (abs(n4Tilde_q.Dot(n4_q)-1.0)>0.001 || abs(n4Tilde_q.M2())>0.001 || abs(n4_q.M2())>0.001) {
                        n4_q.Print();
                        n4Tilde_q.Print();
                        return entry;
                    }
                    if (entry%10000==0){
                        printf("%d nqSq, nqTildeSq, nq(nqTilde) = %10.3g, %10.3g, %10.3g \n",
                               entry, n4_q.M2(), n4Tilde_q.M2(), n4Tilde_q.Dot(n4_q));
                    }
                    Delta_perp = Delta-(Delta.Dot(n4_q))*n4Tilde_q - (Delta.Dot(n4Tilde_q))*n4_q;
                    h_DeltaPerp->Fill(-Delta_perp.M2());
                    h_norm->Fill(n4Tilde_q.Dot(n4_q));
                    //h_Zero->Fill(n4_q.M2());
                    //h_Zero->Fill(n4Tilde_q.M2());
                    //
                    // Now repeat for exact kinematics
                    q4Virtual = k4BeamGen-k4ScatGen;
                    q_dot_P   = q4Virtual.Dot(P4BeamGen);
                    deltaQ    = -q4Virtual.M2()*P4BeamGen.M2()/(q_dot_P*q_dot_P);
                    sqrt_one_d= sqrt(1.0+deltaQ);
                    norm      = 1.0/(sqrt(q_dot_P)*sqrt_one_d);
                    n4_q  = (1.0+sqrt_one_d)* q4Virtual;
                    n4_q += ( -q4Virtual.M2()/q_dot_P ) * P4BeamGen;
                    n4_q *= norm/2.0;
                    n4Tilde_q  = (-P4BeamGen.M2()/((1.0+sqrt_one_d)*q_dot_P)) *q4Virtual ;
                    n4Tilde_q += P4BeamGen;
                    n4Tilde_q *= norm;
                    Delta_perpGen = Delta4Gen -(Delta4Gen.Dot(n4_q))*n4Tilde_q - (Delta4Gen.Dot(n4Tilde_q))*n4_q;
                    h_DeltaPerp_err->Fill(-Delta_perp.M2()-Delta_perpGen.M2()+2.0*Delta_perp.Dot(Delta_perpGen));
                    //
                    //
                    for (int kpart=0; kpart<branchTrack->GetEntries(); kpart++){
                        thisTrack = (Track *) branchTrack->At(kpart);
                        if(thisTrack->PID>2000){
                            h_etaTarget->Fill(thisTrack->Eta);
                            Delta = thisTrack->P4()-P4BeamGen;
                            h_tIon->Fill(Delta.M2()-DeltaSqGen);

                        }
                    }
                }
            }
        }
        // If event contains  at least 2 particles
    }
    printf("Electrons, Photons Jets found in branches = (%d, %d, %d) \n",
           nElectron,nPhoton,nJet);
   
    TFile *froot= new TFile("histOut.root","recreate");

    TCanvas *c1 = new TCanvas("c1","Delphes",75,50,900,600);
    c1->SetLogy();
    c1->Divide(3,2);
    c1->cd(1);
    h2_Q2_xBj->Draw("cont1z");
    gPad->SetLogy();  gPad->SetLogx();
    c1->cd(2);
    h_D_err->SetLineColor(kBlue);
    h_D_err->SetLineWidth(2);
    h_D_err->Draw();
    double rms = h_D_err->GetRMS();
    TF1 *f1 = new TF1("f1","gaus",-rms/2.,rms/2.);
    h_D_err->Fit("f1","","",-rms/3.,rms/3.);
    gStyle->SetOptFit(1111);
    gPad->SetLogy();
    h_D_errPos->SetLineColor(kRed);
    h_D_errPos->SetLineWidth(2);
    h_D_errPos->Draw("same");
    h_DDelta2_naive->SetLineColor(kGreen);
    h_DDelta2_naive->Draw("same");
    
    c1->cd(3);
    /*
    h_MX2->Draw();
    h_MX2_0->SetLineColor(kGreen);
    h_MX2_0->Draw("same");
    h_Zero->SetLineColor(kRed);
    h_Zero->Draw("same");
    gPad->SetLogy();
     */
    h_DeltaPerp_err->SetLineColor(kBlue);
    h_DeltaPerp_err->Draw();
    rms = h_DeltaPerp_err->GetRMS();
    double max = h_DeltaPerp_err->GetMaximum();
    TF1 *fPerp = new TF1("fPerp","[0]*exp(-x/[1])+[2]*exp(-x/[3])",0.0,0.5);
    fPerp->SetParameters(max/2.0,rms/2.0,max/2.0,rms*2.0);
    // Need at least 3 exponential terms ==>Meijer G-function
    //h_DeltaPerp_err->Fit("fPerp","I","",0.0,0.2);
    /*
    h_DeltaPerp->SetLineColor(kGreen);
    h_DeltaPerp->Draw("same");
    TF1 *ferf = new TF1("ferf","[0]*erfc((x-[1])/[2])",0.0,3.0);
    ferf->SetParameters(1500.,1.5,0.1);
    h_DeltaPerp->Fit("ferf","","",0.0,3.0);
*/

    //h_Zero->SetLineColor(kRed);
    //h_Zero->Draw("same");
    gPad->SetLogy();
    /*
    h2_errD_dE->Draw("cont1z");
    //h2_errD_D->Draw("cont1z");
    gPad->SetLogz();
    gPad->SetRightMargin(0.15);
    gPad->SetLeftMargin(0.15);
     */
    c1->cd(4);
    h_qprime->Draw();
    gPad->SetLogy();
    c1->cd(5);
    //h2_Q2err_y->Draw("cont1z");
    //h2_Q2err_eta->Draw("cont1z");
   // h2_qPrDiff_EvsP->Draw("cont1z");
    //gPad->SetLogz();
    //h_etaTarget->Draw();
    //h_tIon->Draw();
    h_dqPrime->Draw();
    gPad->SetLogy();
    /*
    h_PT->Draw();
    h_E->SetLineColor(kRed);
    h_E->Draw("same");
     */
    c1->cd(6);
    h2_P4Gen->Draw("Cont1z");
    
    /*
     h_Delta2gen->SetLineColor(kRed);
    h_Delta2gen->Draw();
    gPad->SetLogy();
    h_Delta2->SetLineColor(kBlue);
    h_Delta2->Draw("same");
     */
    
    TCanvas *c2 = new TCanvas("c2","Delta Perp",150,75,400,400);
    h_DeltaPerp_err->SetLineWidth(2);
    h_DeltaPerp_err->SetLineColor(kBlue);
    TH1 *hc = h_DeltaPerp_err->GetCumulative(kTRUE,"_hc");
    //TH1 *hc = h_DeltaPerp_err->GetMaximum()/h_DeltaPerp_err->GetIntegral();
    //hc->Draw();
    h_DeltaPerp_err_hc->Draw();
    //gPad->SetLogy();

    h_DeltaPerp_err->Write();
    hc->Write();
    h_DeltaPerp->Write();
    h_D_err->Write();
    froot->Close();
    return entry;
}
