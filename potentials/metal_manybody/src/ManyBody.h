#pragma once


#include <fstream>
#include <vector>
#include <iomanip>
#include <string>

#include "Eigen/LU"

#include "Constants.h"
#include "Species.h"
#include "MetalPotential.h"

enum keymetal
{
    sutchen = 1,
    fnsc,
    xfnsc,
    gupt,
    eam,
    meam,
    metalend
};

class ManyBody
{
private:

    int 
        maxParameter = MAXPAR; // the max number of energy values per potential unless EAM where this is changed

    int
        numberPotentials = 0;

    int
        maxInteraction = 0,
        numSpec = 0,
        numPot = 0;

    std::vector<MetalPotential>
        metPots;  // stores the parameters

    int
        maxRepParameter = 0,
        maxDensParameter = 0,
        maxEmbedParameter = 0;

    Eigen::VectorXi
        meshType;                 

    double
        cutoff = 0.0;

    Eigen::MatrixXi
        meshTable;    // links atom type to look up mesh

    Eigen::MatrixXd
        
        mesh;
        
    bool
        potSutChen = false,
        potFinSin = false,
        potGupta = false,
        potEAM = false;
       

    int potenKey(std::string keyWord);

public:

    void calculateEnergyMesh(const Species& ele, double shortrangecut, std::ofstream& outStream);

    void calculateEAMEnergyMesh(const Species& ele, std::ofstream& outStream);

    double calculateManyBodyDensEnergy(int ltypea, int ltypeb, double r);

    double embed(int ltypea, int ltypeb, double rho);

    void calculateManyBodyPairEnergy(int ltypea, int ltypeb, double r, double rsq, double& vPair);

    void calculateManyBodyPairForce(int ltypea, int ltypeb, double r, double rsq, double& vPair, double& forc);

    void calculateManyBodyForce(int ltypea, int ltypeb, double r, double rhoi, double rhoj, double& fDens);

    
    /********************************************************************
    returns pair potential energy and force
    *******************************************************************/
double eamMeshEnergy(Eigen::VectorXd& energyMesh, Eigen::VectorXd& distMesh, double r, double dr, int index);

double eamPairForce(double r, double dr, int potnum, int index)
{
//force_pair: this function returns the pair potential part of the force
    double forc = 0.0;
    //double forc = -(meshRepFrc[potnum + index] + (meshRepFrc[potnum + index + 1] - meshRepFrc[potnum + index])*(r - distRepEng[potnum + index]) / dr);
    return forc;
}

    
    void checkMesh(void);

    ManyBody& operator=(const ManyBody& src);

    ManyBody(const ManyBody& src);

    ManyBody();
    virtual ~ManyBody();

    void loadPotential(int, std::ifstream&, std::ofstream&);          // dimension and load potential parameters

    void loadPotentialFit(int, std::ifstream&, std::ofstream&);

    void calculateDerivativeMesh(Eigen::VectorXd& pot, Eigen::VectorXd& deriv, int atomType);

    void printPotential(std::ofstream&); 	                    // print potentials  

    void convertPotential(void);

};

