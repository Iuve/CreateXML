/* 
 * File:   FermiDistribution.h
 * Author: aleksandra
 *
 * Created on 31 maj 2012, 14:57
 * 
 */
#include <vector>
#ifndef FERMIDISTRIBUTION_H
#define	FERMIDISTRIBUTION_H

class FermiDistribution{
public:
    /**
     * Calculation of beta energy distribution (works for b+ and b-)
     * Algorytm:
     * Sampling the Fermi distribution for B-dacey energy input to EGS4, W. R. Nelson, J. Liu, SLAC-TN_92-1 
     * eCharge - possible only +/- 1 value. +1 - beta plus, -1 - beta minus
     */
    FermiDistribution(int atomicNumber, double qVal, int eCharge);
    virtual ~FermiDistribution();

    //double GetRandomBetaEnergy();
    double GetAverageBetaEnergy() {return avergeBetaEnergy_;}
private:
    void CalculateBetaEnergyDistribution();
    double FindProbabilityDensityValue(double x);

    double EulerGammaFunction (double finalValue);
    int atomicNumber_;
    double qVal_;
    
    
    std::vector<double> betaEnergy_;
    std::vector<double> energProbability_; 
    
    int nrOfSteps_ ;
    

    int nrOfkeVperChan_;
    double gammaFinal_ ;//end of calculation euler gamma function
    double alpha_ ;
    double pi_ ;//2*atan(1);
    double electronMass_ ;
    
    
    int eCharge_ ;//+1 beta plus, -1 beta minus
    double v0_;
    double S_;
    double xMin_;
    double xMax_;
    double maxBetaProbability_;
    double delta0_;
    double avergeBetaEnergy_;
};

#endif	/* FERMIDISTRIBUTION_H */

