//
//  Deep_pi_pi_MC.hpp
//  
//
//  Created by Hyde, Charles E. on 5/15/21.
//

#ifndef CoherentDVCS_MC_hpp
#define CoherentDVCS_MC_hpp

#include <stdio.h>
#include <TCanvas.h>
#include <TH2D.h>
#include <TRandom3.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TMath.h>
#include <TRandom3.h>
#include <TSystem.h>
/*
#include "HepMC3/GenEvent.h"
#include "HepMC3/WriterAscii.h"
#include "HepMC3/ReaderRoot.h"
#include "HepMC3/Print.h"
*/
#include <TDatabasePDG.h>

//using namespace HepMC3;

const double TwoPi = 2.0*TMath::Pi();
// Access physical constants, masses in pdg.lbl.gov database implemented in root
auto dbPDG = TDatabasePDG::Instance();

//  Global random number generator
TRandom3 ran3;

// HepMC Event clas constructors


// Define Global four-vectors (all in Detector Frame)
/**
 * @defgroup FourVectors Global event-by-event four-vectors
 * @{ \*
\**
 * Global (nominal) beam four-vectors
 * \f$ k_0^\mu,\quad P_0^\mu,\quad Y^\mu = [0,0,1,0] = {\rm Vertical Up} \f$
 */
TLorentzVector k4Beam0, P4Beam0;
TLorentzVector Y4det;
/** Event-by-event kinematic four-vectors
 *\f[ k^\mu,\, P^\mu,\, k^{\prime\mu},\, q^\mu = (k=k')^\mu,\, P_{\pi\pi}^\mu,\, P^{\prime\mu},
 *   \Delta^\mu = (P'-P)^\mu, p^\mu_{\pi^+}, p^\mu_{\pi^-}
 * \f] */
TLorentzVector k4Beam, P4Beam, k4Scat, q4Virt, q4Prime, P4Scat;
TLorentzVector Delta4vec, P4Tot;

/** Global (nominal) light-cone vectors
 *\f[ n_{e,0}^\mu,\, \widetilde n_{e,0}^\mu,\, X_{e,0}^\mu,\, Y_{e,0}^\mu,\, Z_{e,0}^\mu,\, T_{e,0}^\mu
 * \f] */
TLorentzVector n4_e0, n4Tilde_e0, X4_e0, Y4_e0, Z4_e0, T4_e0;

/** Event-by-event \f$eP\f$ and \f$qP\f$ light-cone four-vectors
 *\f[ n_{e}^\mu,\, \widetilde n_{e}^\mu,\, n_q^\mu \widetilde n_q^\mu
 * \f] */
TLorentzVector n4_e, n4Tilde_e, n4_q, n4Tilde_q;

/** Event-by-event cartesian four-vectors
 *\f[ X_{e}^\mu,\, Y_{e}^\mu,\, Z_{e}^\mu,\, T_{e}^\mu
 * \f] */
TLorentzVector Y4_Det;
TLorentzVector X4_e, Y4_e, Z4_e, T4_e;

/**
 *
 * cartesian four-vectors
 *\f[ X_{q}^\mu,\, Y_{q}^\mu,\, Z_{q}^\mu,\, T_{q}^\mu
 * \f] */
TLorentzVector X4_q, Y4_q, Z4_q, T4_q, X4_pipi, Y4_pipi;

/** Event-by-event, Boosted to \f$ q+P\f$ CM frame
 *\f[ X_{q,CM}^\mu,\, Y_{q,CM}^\mu,\, n_{q,CM}^\mu,\, \widetilde n_{q,CM}^\mu
 * \f] */
TLorentzVector X4_qCM, Y4_qCM, n4_qCM, n4Tilde_qCM;

/** @} */
/** Invariants, definied in routine Init()
 */
double mLepton, MIon, mPion;

//  This was a qft++ structure.  Must be redefined
//LeviCivitaTensor LeviCivita4;

//  Beam Emittance arrays
double emitt_e[3], emitt_i[3];//  Geometrical x,y emittance and sigma(p)
double betaIP_e[2], betaIP_i[2]; //  x,y beta-functions at IP
bool  eSmear(false), iSmear(false);
int nEvents;

// Monte Carlo bounds
double Q2Min, Q2Max, yMin, yMax, csCMmin, csCMmax, W2Threshold;
//  Light Cone constants
double sqrtDL, deltaL, sqrtDQ, deltaQ;

// Kinematic Lorentz Invariants
double yInv, sMinusM2, xBj, psf, phi_e, s_e, W2, Q2;
// Four-Vector contractions
double k_dot_P;

// Physical constants
const double alpha = 1.0/137.03;
const double pi = acos(-1.0);
const double umass = 0.931494; // GeV, atomic mass unit u
int ionZ, ionA;

// Detector Regions (EMCal)
const double etaMin = -4.0;
const double etaEl  = -1.7;
const double etaBe  = -1.7;
const double etaB0  =  0.0;
const double etaBi  =  1.4;
const double etaMax =  4.0;
const double re_EMCal = 1.2;
const double rBarrel  = 0.6;
const double ri_EMCal = 2.6;

#endif /* Deep_pi_pi_MC_hpp */
