#include "NbrListPBC.h"
NbrListPBC::NbrListPBC()
{
    maxNbrs = 0;

}

NbrListPBC::~NbrListPBC()
{

}

void NbrListPBC::allocNbrListArrays(int natoms, int num)
{

    maxNbrs = num;
    int mx = natoms * num;
    numNbrs.resize(natoms); 
    map.resize(mx); 

    rZero.resize(natoms,3); 
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// calculation of neighbour lists
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void NbrListPBC::calculateNbrList(const Eigen::MatrixXd &pos, const Eigen::MatrixXd &latVector, const Eigen::MatrixXd &rcpVector, 
                                  int natoms, double cutoff)
{
    int
        aatom,
        batom,
        idx,
        num;

    double
        ax, ay, az,
        bx, by, bz,
        xx, yy, zz,
        rx, ry, rz,
        radius, rsq;

    radius = cutoff * cutoff;
    
    //for (aatom = rank; aatom < natoms - 1; aatom = aatom + numProcs)
    for (aatom = 0; aatom < natoms; aatom++)
    {

        num = 0;
        idx = maxNbrs * aatom;
        //std::cout << "\n atom nbrlist " << aatom;
        //set new start for position check accummulators
        rZero(aatom,0) = 0.0;
        rZero(aatom,1) = 0.0;
        rZero(aatom,2) = 0.0;

        ax = pos(aatom,0);
        ay = pos(aatom,1);
        az = pos(aatom,2);

        for (batom = 0; batom < natoms; batom++)
        {

            if (aatom == batom)
                continue;

            bx = pos(batom,0);
            by = pos(batom,1);
            bz = pos(batom,2);

            rx = ax - bx;
            ry = ay - by;
            rz = az - bz;

            xx = rx * rcpVector(0,0) + ry * rcpVector(0,1) + rz * rcpVector(0,2);
            yy = rx * rcpVector(1,0) + ry * rcpVector(1,1) + rz * rcpVector(1,2);
            zz = rx * rcpVector(2,0) + ry * rcpVector(2,1) + rz * rcpVector(2,2);

            xx -= rint(xx);
            yy -= rint(yy);
            zz -= rint(zz);

            rx = xx * latVector(0,0) + yy * latVector(0,1) + zz * latVector(0,2);
            ry = xx * latVector(1,0) + yy * latVector(1,1) + zz * latVector(1,2);
            rz = xx * latVector(2,0) + yy * latVector(2,1) + zz * latVector(2,2);

            rsq = rx * rx + ry * ry + rz * rz;

            if (rsq <= radius)
            {

                map[idx+num] = batom;
                num++;

                if (num >= maxNbrs)
                {
                    std::cout << "\n exceeded max nbrs " << maxNbrs << " " << num << std::endl;
                    std::cout.flush();
                    exit(EXIT_FAILURE);

                }

            }

        }

        numNbrs[aatom] = num;
        //std::cout << "\n PBC nbrlist " << aatom << " " << num;

    }

}

bool NbrListPBC::checkNbrList(int i, double deltaX, double deltaY, double deltaZ, double verletShell)
{
    double
        dr,
        rMax = pow(verletShell, 2.0);

    bool
        check = false;

    rZero(i,0) += deltaX;
    rZero(i,1) += deltaY;
    rZero(i,2) += deltaZ;
    
    dr = rZero(i,0) * rZero(i,0) + rZero(i,1) * rZero(i,1) + rZero(i,2) * rZero(i,2);

    if (dr > rMax)
        check = true;

    return check;
}
