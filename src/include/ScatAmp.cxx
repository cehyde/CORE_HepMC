//
//  ScatAmp.cxx
//
//
//  Created by Charles Earl Hyde on 13-June-2021.
//
//
/** @file include/ScatAmp.cxx
 *  @brief File containing  functions  to calculate scattering amplitude for  \f$e p \to e' p' \pi \pi\f$
 */

//double Diffract();
//double piN();

double LeptonSymmTensor(int lambda_q1, int lambda_q2){
    /** @brief Calculate the symmetric lepton tensor
     * Virtual photon helicities: lambda_q1, lambda_q2 \f$= \lambda_q(i) \in {-1,0,1} \f$
     */
    double L_S = 0.0;
    double KdotX = k4Beam.Dot(X4_q)+k4Scat.Dot(X4_q);
    if ((lambda_q1==0)&&(lambda_q2==0))
    {
        // neither helicity zero
        L_S = (2.0-yInv)*(2.0-yInv)/(yInv*yInv) - (1+deltaQ);
        L_S *= Q2/(1+deltaQ)
    } else if (lambda_q1*lambda_q2==0)
    {
       // one helicity zero, but not both zero
        L_S = KdotX*(2.-yInv)/yInv *sqrt(Q2/2./(1.0+deltaQ));
        L_S*= -(lambda_q1+lambda_q2);
    } else
    {
        // both non-zero
        L_S = KdotX*kDotX + Q2(1.+lambda_q1*lambda_q2);
        L_S *=lambda_q1*lambda_q2/2.0
    }
    return L_S;
}
double LeptonAntiTensor(int lambda_q1, int lambda_q2)
{
    /** @brief Calculate the anti-symmetric lepton tensor
     * Virtual photon helicities: lambda_q1, lambda_q2 \f$= \lambda_q(i) \in {-1,0,1} \f$
     */

    double L_A = 0.0;
    return L_A;
}

