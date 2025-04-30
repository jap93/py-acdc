#pragma once

#include <fstream>
#include <iomanip>
#include <algorithm>
#include <stdlib.h>
#include <vector>
#include <complex>
#include <math.h>

#include "Constants.h"

#include <Eigen/LU>
 
class Ewald
{
private:

    double
        etaMad,
        rcpSpaceCut;                    // real space cutoff

    int
        gcellx = 1,
        gcelly = 1,
        gcellz = 1;                           // number of cells over which the real space summation

    int
        maxGVec = 0;

    int
        numGVec = 0;                         // number of reciprocal space vectors

    Eigen::VectorXd
        gx,
        gy,
        gz,
        g0,
        g1;                   

    bool newJob = true;

    Eigen::VectorXcd
        /** atom rciprocal space sums */
        rcpSum;


public:

    Ewald();
    virtual ~Ewald();

    void printNumGVec(std::ofstream&);

    void setEwaldImage(double, double, double, double, const Eigen::VectorXd &latVector, const Eigen::VectorXd &rcpVector);

    void printEwald(std::ofstream&);

    void allocateArrays(void);

    void estimateGvector(const Eigen::MatrixXd &rcpVector);

    void findGvector(double, const Eigen::MatrixXd &);

    double recipSpaceForce(const Eigen::MatrixXd &pos, Eigen::MatrixXd& force, Eigen::VectorXd , Eigen::VectorXd  atomcharge, int numatoms);


};

