#pragma once

#include <string>
#include <stdio.h>
#include <algorithm>

#include <Eigen/Dense>

#include "Species.h"
#include "Constants.h"

#include "Ewald.h"
#include "NbrListPBC.h"
#include "TwoBody.h"

#ifdef _OPENMP
    #include <omp.h>
#endif

class Field
{
private:

    int
        /** the type of Ewald used 0 = none, 1 = Ewald, 2 = damped and shifted */
        doEwald = 1,

        // the number of cells to be searched over when not using the image convention
        cellX = 1,
        cellY = 1,
        cellZ = 1,

        /** the maximum number of neighbours */
        maximumNbrs = 0;

       
    bool
        /** the neighbour list is used if true - faster but can have problems in ART*/
        useNbrList = false,

        /** flag to indicate that the nbr list requires updating */
        recalcNbrList = false,

        /** turns of image convention if true */
        minimumImage = true,

        /** flag to indicate first time through */
        newjob = false;

    double
        /** the accuracy of the Ewald sum */
        madelungAcc = 1.0e-6,

        /** convergence factor for shifted and damped electrostatics */
        alpha = 0.2,

        /** limit of extent of two-body interaction */
        shortRangeCut = 15.0,

        /** eta for Ewald sum */
        eta = 0.0,
        
        /** the extent of the shell around the cutoff */
        verletShell = 1.0,
        
        /** conversion of Coulomb energy to desired units */
        cFact = 1.0,

        /** conversion unit for all energies */
        energyUnit = 1.0;


    double
        delMesh,
        rMesh,
        aa, bb,
        rfld0, rfld1;

    Eigen::VectorXd
        erfc_e,
        erfc_f;

    NbrListPBC
        nbrList; 
        
    TwoBody
        vdw;

    Ewald
        coul;


    void fieldErfcGen(double alpha);

    void hfunc1(double, double&, double&);

    void hfunc2(double, double&, double&);

    void resetSimulationCellNoImage(Eigen::MatrixXd pos, Eigen::MatrixXd latVector, Eigen::MatrixXd rcpVector, int natoms);

    void resetSimulationCellMinImage(Eigen::MatrixXd pos, Eigen::MatrixXd latVector, Eigen::MatrixXd rcpVector, int natoms);

    void realSpaceForceImage(double& realenergy, double& twoenergy, const Eigen::MatrixXd &posX, Eigen::MatrixXd &force, Eigen::MatrixXd &stress, const Eigen::MatrixXd &latVector, 
                            const Eigen::MatrixXd &rcpVector, const Eigen::VectorXd &atomCharge, const Eigen::VectorXi &atmLabel, 
                            int numAtoms);

    void realSpaceForceImageNbrList(double& realenergy, double& twoenergy, const Eigen::MatrixXd &posX, Eigen::MatrixXd &force, Eigen::MatrixXd &stress, const Eigen::MatrixXd &latVector, 
                            const Eigen::MatrixXd &rcpVector, const Eigen::VectorXd &atomCharge, const Eigen::VectorXi &atmLabel, 
                            int numAtoms);

    void realSpaceForceNoImage(double& realenergy, double& twoenergy, const Eigen::MatrixXd &posX, Eigen::MatrixXd &force, Eigen::MatrixXd &stress, const Eigen::MatrixXd &latVector, 
                            const Eigen::MatrixXd &rcpVector, const Eigen::VectorXd &atomCharge, const Eigen::VectorXi &atmLabel, 
                            int numAtoms);

    //void fieldErfcGen(double alpha);
    //void interpolate_shifted_damped_potential(double, double, double, double&);
    void interpolate_shifted_damped_potential(double, double, double, double&, double&);

public:

    // constructor
    Field();

    // default destructor
    virtual ~Field();

    /**
     * return the type of periodic boundary as bool: True for minimum image convention
    */
    bool getImageConvention()
        {return minimumImage;}
        


    /**
     * sets up the environment to underatke a calculation using the desired calculator
     * @param mpi (const MPIComms&) : the mpi communicatior
     * @param latVector (Eigen::VectorXd) : lattice vectors
     * @param rcpVector (Eigen::VectorXd) : reciprocal lattice vectors
     * @param minDimension (double) : the smallest dimension of the cell
     * @param volume (double) | the cell volume
     * @param spec (Species&) : the elemental data
     * @param outStream (ofstream&) : the output stream
     */
    void setup(const Eigen::MatrixXd &latVector, const Eigen::MatrixXd &rcpVector, double minDimension, double volume, int numAtoms, const Species& spec, std::ofstream& outstream);

    /**
     * reads any information reuired to calculate the energy. If an incorrect keyword is used an invalid argument is thrown.
     * 
     */
    void readPotential(const std::string &fileName);

    //void calculateForces(Eigen::VectorXd posX, Eigen::VectorXd posY, Eigen::VectorXd posZ, Eigen::VectorXd forceX, Eigen::VectorXd forceY, Eigen::VectorXd forceZ, 
    //                     Eigen::VectorXd stress, Eigen::VectorXd latVector, Eigen::VectorXd rcpVector, Eigen::VectorXd atomcharge, int* intlabel, int* frozen, int numatoms, Energy& eng);
    Eigen::MatrixXd &calculateForces(const Eigen::MatrixXd &posX, Eigen::MatrixXd &force, Eigen::MatrixXd &stress, const Eigen::MatrixXd &latVector, 
                            const Eigen::MatrixXd &rcpVector, const Eigen::VectorXd &atomCharge, const Eigen::VectorXi &atmLabel, const Eigen::VectorXi &frozen,
                            int numAtoms);
};

