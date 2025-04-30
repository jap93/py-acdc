#pragma once

#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdlib.h>
#include <math.h>

#include "Eigen/LU"

#include "Species.h"
#include "Constants.h"

#define MAXPAR 9

struct MetalPotential
{
    int
        potentialType = -1,
        numComponents = 0,
        numPtsRep = 0,
        numPtsDens = 0,
        numPtsEmbed = 0;

    double
        startRep = 0.0,
        finishRep = 0.0,
        startDens = 0.0,
        finishDens = 0.0,
        startEmbed = 0.0,
        finishEmbed = 0.0,
        stepEmbed = 0.0,
        stepDens = 0.0,
        stepRep = 0.0,
        potCutoff = 0.0;

    double
        param1 = 0.0,                  // short- range interaction parameters for parameterised potentials
        param2 = 0.0,
        param3 = 0.0,
        param4 = 0.0,
        param5 = 0.0,
        param6 = 0.0,
        param7 = 0.0, 
        param8 = 0.0,
        param9 = 0.0;

    Eigen::VectorXd

        parRep,
        parDens,
        parEmbed,
        derivRep,
        derivDens,
        derivEmbed,
        distRep,
        distDens,
        distEmbed;         // max distance

    std::string
        atA,
        atB;

    void readEAM(std::string infile, std::ofstream& outStream);

    void readLEAM(std::string infile, std::ofstream& outStream);

    void readParameters(Eigen::VectorXd& par, int numPts, std::ifstream& inStream);

    void calculateDerivativeMesh(Eigen::VectorXd& pot, Eigen::VectorXd& deriv, double dr, int npts);

    MetalPotential();

    virtual ~MetalPotential();
    
};

