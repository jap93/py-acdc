
#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>

#include <string>
#include <stdio.h>
#include <algorithm>
#include <vector>

#include <Eigen/LU>

#include "Species.h"
#include "Constants.h"

#include "NbrListPBC.h"
#include "ManyBody.h"

#ifdef _OPENMP
    #include <omp.h>
#endif


void dcell(double*, double[10]);
std::vector<std::string> split(std::string s);

class Metal
{
private:

    int
        
        // the number of cells to be searched over when not using the image convention
        cellX = 1,
        cellY = 1,
        cellZ = 1,

        /** the size of the density work arrays */
        wrkSize = 0,

        /** the maximum number of neighbours */
        maximumNbrs = 0;


    bool
        /** the neighbour list is used if true - faster but can have problems in ART*/
        useNbrList = false,

        /** flag to indicate that the nbr list requires updating */
        recalcNbrList = false,

        /** turns of image convention if true */
        minimumImage = true,

        /** indicates first pass through */
        newJob = true;

    double

        /** limit of extent of two-body interaction */
        shortRangeCut = 15.0,

        /** the extent of the shell around the cutoff */
        verletShell = -0.5,
        
        /** conversion unit for all energies */
        energyUnit = 1.0;

        
    Eigen::VectorXd
        /** temp storage of atom density */
        tmpRho,

        /** work storage of rho */
        wrkRho;

    Eigen::MatrixXd 
        force; 
    
    NbrListPBC
        nbrList; 
        
    ManyBody
        mbdy;

    Species
        spec;

    std::ofstream
        outStream;

    std::string
        name;
  

    void resetSimulationCellNoImage(Eigen::MatrixXd& pos, const Eigen::MatrixXd& latVector, const Eigen::MatrixXd& rcpVector, int natoms);

    void resetSimulationCellMinImage(Eigen::MatrixXd& pos, const Eigen::MatrixXd& latVector, const Eigen::MatrixXd& rcpVector, int natoms);

    void realSpaceForceImage(double& manyEnergy, double& twoEnergy, const Eigen::MatrixXd &posX, Eigen::VectorXd &stress, const Eigen::MatrixXd &latVector, 
                            const Eigen::MatrixXd &rcpVector, const Eigen::VectorXi &atmLabel, 
                            int numAtoms);

    void realSpaceForceImageNbrList(double& manyEnergy, double& twoEnergy, const Eigen::MatrixXd &pos, Eigen::VectorXd &stress, const Eigen::MatrixXd &latVector, 
                            const Eigen::MatrixXd &rcpVector, const Eigen::VectorXi &atmLabel, 
                            int numAtoms);

    void realSpaceForceNoImage(double& manyEnergy, double& twoEnergy, const Eigen::MatrixXd &posX, Eigen::VectorXd &stress, const Eigen::MatrixXd &latVector, 
                            const Eigen::VectorXi &atmLabel, int numAtoms);

    double calculateAtomRho(int aatom, const Eigen::MatrixXd& pos, const Eigen::MatrixXd& latVector, const Eigen::MatrixXd& rcpVector, 
                                      const Eigen::VectorXi& atmLabel, int numAtoms);

    double calculateAtomRhoNoImage(int aatom, const Eigen::MatrixXd& pos, const Eigen::MatrixXd& latVector, 
                                      const Eigen::VectorXi& atmLabel, int numAtoms);                                  
public:

   // constructor
    //RigidIon();
    Metal(const std::string &name) : name(name) { }
    // default destructor
    //virtual ~RigidIon();

    /**
     * return the type of periodic boundary as bool: True for minimum image convention
    */
    bool getImageConvention()
        {return minimumImage;}
        
    void finalise(void)
        {outStream.flush();
         outStream.close();}

    /**
     * sets up the environment to underatke a calculation using the desired calculator
     * @param mpi (const MPIComms&) : the mpi communicatior
     * @param latVector (Eigen::VectorXd) : lattice vectors
     * @param rcpVector (Eigen::VectorXd) : reciprocal lattice vectors
     * @param minDimension (double) : the smallest dimension of the cell
     * @param volume (double) | the cell volume
     */
    void setup(const Eigen::VectorXd &cellProp, const Eigen::MatrixXd &rcpVector, const Eigen::VectorXd &rcpProp, double minDimension, double volume, int numAtoms);

    /**
     * reads any information reuired to calculate the energy. If an incorrect keyword is used an invalid argument is thrown.
     * 
     */
    void readPotential(const std::string &fileName);

    Eigen::MatrixXd &calculateForces(Eigen::MatrixXd &pos, Eigen::VectorXd &stress, const Eigen::MatrixXd &latVector, 
                            const Eigen::MatrixXd &rcpVector, const Eigen::VectorXi &atmLabel, const Eigen::VectorXi &frozen,
                            int numAtoms);

};


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// utility functions
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void Metal::setup(const Eigen::VectorXd &cellProp, const Eigen::MatrixXd &rcpVector, const Eigen::VectorXd &rcpProp, double minDimension, double volume, 
                  int numAtoms)
{


    double
        radius;

    if (shortRangeCut < 1.0)
    {
        outStream << "\n\n\n *** short range cutoff is too small : " << shortRangeCut << std::endl;
        outStream.flush();
        exit(EXIT_FAILURE);
    }

    if (useNbrList)
    {
        radius = shortRangeCut + verletShell;
    }
    else
    {
        radius = shortRangeCut;
    }

    // check the image convention is obeyed
    if (minimumImage && radius > (minDimension / 2))
    {
        outStream << "\n\n\n *** short range shortRangeCut to large for simulation box." << std::endl;
        outStream << " cutoff : " << shortRangeCut << " minimum dimension: " << minDimension << std::endl;
        outStream.flush();
        exit(EXIT_FAILURE);
    }

    if (newJob == false)
        return;

    if (useNbrList && verletShell < 1e6)
    {
        outStream << "\n\n\n *** WARNING: the verlet does not appear to have been set." << std::endl;
        //outStream.flush();
        //mpi.commsAbortWorld();
        //exit(EXIT_FAILURE);
    }

    if (!minimumImage)   // setup for non-Image convention
    {
        cellX = int(shortRangeCut / cellProp[0]) + 1;
        cellY = int(shortRangeCut / cellProp[1]) + 1;
        cellZ = int(shortRangeCut / cellProp[2]) + 1;  

        outStream << "\n\n Image convention is not being used!" << std::endl;
        outStream << "\n The number of cells used in the summation are " << cellX << " " << cellY << " " << cellZ << std::endl;       
    }
    else
    {
        outStream << "\n\n\n short range cutoffs for lattice summation (image):"
                  << "\n real space shortRangeCut       = " << shortRangeCut << " Angstroms" << std::endl;
    }
    
    //alocate the density work arrays
    outStream << "\n allocating temporary metal density arrays of size " << numAtoms << std::endl;

    tmpRho.resize(numAtoms);
    wrkRho.resize(numAtoms);

    // set up short-range potential mesh
    // set up short-range potential mesh
    if (newJob == true) mbdy.calculateEnergyMesh(spec, shortRangeCut, outStream);

    force.resize(numAtoms+1,3);

    newJob = false;

}

/******************************************************
reads in description of the species and the potential model
*******************************************************/
void Metal::readPotential(const std::string &fileName)
{
    std::string
        line,
        subWord,
        keyWord,
        dummy,
        dummy2;

    int
        num;

    std::ifstream
        inStream;

    //open inpot file for reading
    inStream.open(fileName, std::ios::in);
    //open output file for information/errors
    outStream.open("field.log", std::ofstream::out);

    // read in and print out the species and potentails
    std::vector<std::string> words;

    while (!inStream.eof())
    {                                             // or start directive

        std::getline(inStream, line);
        transform(line.begin(), line.end(), line.begin(), ::tolower);
	    words = split(line);

        keyWord = words[0];

        if (keyWord == "cutoff")
		{
			dummy = words[1];
			shortRangeCut = std::stod(dummy);
		}
        else if (keyWord == "species")
        {

            dummy = words[1];
            num = std::stoi(dummy);
            spec.loadSpecies(inStream, num);
            spec.printSpecies(outStream);
        }
		else if (keyWord == "nbrlist")
		{
			useNbrList = true;
			dummy = words[1];
			verletShell = std::stod(dummy);
		}
        else if (keyWord == "maxNbrs")
		{
			dummy = words[1];
			maximumNbrs = std::stoi(dummy);
		}
		else if (keyWord == "noimage")
		{
			minimumImage = false;
		}
        else if (keyWord == "many" || keyWord == "manybody")
        {

            dummy = words[1];
            num = std::stoi(dummy);

            // get the correct energyUnits
            dummy = words[2];

            mbdy.loadPotential(num, inStream, outStream);
            mbdy.printPotential(outStream);

        }
        else if (keyWord == "close")
        {
            break;
        }
        else
        {
            outStream << "\n unrecognised key word used in potentials file " << keyWord << "\n from line " << line << std::endl;
            outStream.flush();
            exit(EXIT_FAILURE);
            //throw std::invalid_argument("Invalid keyword!");
        }

    }

    outStream.flush();
}


/************************************************************************
 driver for the calculation of forces and energy of all ions
 ************************************************************************/
Eigen::MatrixXd &Metal::calculateForces(Eigen::MatrixXd &pos, Eigen::VectorXd &stress, const Eigen::MatrixXd &latVector, 
                            const Eigen::MatrixXd &rcpVector, const Eigen::VectorXi &atmLabel, const Eigen::VectorXi &frozen,
                            int numAtoms)
{
    double
        manyEnergy = 0.0,
        twoEnergy = 0.0;

    for (auto i = 0; i < numAtoms; i++)
    {
        //std::cout << atomCharge[i] << " " << atmLabel[i] << std::endl;
        force(i,0) = 0.0;
        force(i,1) = 0.0;
        force(i,2) = 0.0;
    }

    for (auto i = 0; i < 6; i++)
        stress(i) = 0.0;


    if (minimumImage)
    {
        resetSimulationCellMinImage(pos, latVector, rcpVector, numAtoms);

        if (useNbrList)
        {
            //std::cout << "\n should not be here";
            exit(EXIT_FAILURE);
            realSpaceForceImageNbrList(manyEnergy, twoEnergy, pos,  stress, latVector, rcpVector, atmLabel, numAtoms);
        }
        else
        {
            realSpaceForceImage(manyEnergy, twoEnergy, pos, stress, latVector, rcpVector, atmLabel, numAtoms);
        }
               
    }
    else
    {
            resetSimulationCellNoImage(pos, latVector, rcpVector, numAtoms);
            realSpaceForceNoImage(manyEnergy, twoEnergy, pos, stress, latVector, atmLabel, numAtoms);
    }
/*
    for (auto i = 0; i < numAtoms; i++)
    {
        std::cout << atomCharge[i] << " " << atmLabel[i]  << " " << 
        force(i,0)  << " " << 
        force(i,1)  << " " << 
        force(i,2) << std::endl;
    }
*/
    //tag the energies on to the force data to send back
    force(numAtoms,0) = manyEnergy;
    force(numAtoms,1) = twoEnergy;
    //std::cout << "\n energies " << rcpenergy << " " << manyEnergy << " " << twoEnergy;
    return force;

}

 
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// routines using the nearest image convention
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#if defined  _OPENMP


void Metal::realSpaceForceImage(double& manyEnergy, double& twoEnergy, const Eigen::MatrixXd &pos, Eigen::VectorXd &stress, const Eigen::MatrixXd &latVector, 
                            const Eigen::MatrixXd &rcpVector, const Eigen::VectorXi &atmLabel, int numAtoms)
{
    int
        batom,
        ltypea,
        ltypeb,
        aatom;

    double
        ax, ay, az,
        bx, by, bz,
        xx, yy, zz,
        eng_s,
        tmp_energy_s = 0.0,
        tmp_energy_m = 0.0,
        forc_m,
        forc_s,
        forceTot;

    double
        r,
        rsq,
        rx,  // dummy variables for distances
        ry,
        rz,
        radius = pow(shortRangeCut, 2.0);  // real space madelung and short range cut are same
                                           // for image convention so just squared

    double* tmpForceX = nullptr;
    double* tmpForceY = nullptr;
    double* tmpForceZ = nullptr;
    double* tmpStress = nullptr;
    double* tmpRho = nullptr;

    manyEnergy = 0.0;
    twoEnergy = 0.0;

    #pragma omp parallel default(none)  \
        reduction(+: manyEnergy, tmp_energy_s) \
        shared(atmLabel, pos, numAtoms, rcpVector, latVector, force, wrkRho, stress, radius, shortRangeCut) \
        private(ltypea, ltypeb, ax, ay, az, bx, by, bz, rx, ry, rz, r, rsq, xx, yy, zz, aatom, batom) \
        private(eng_s, forc_m, forc_s, forceTot, tmp_energy_m, tmpRho, tmpForceX, tmpForceY, tmpForceZ, tmpStress)
    {
        tmpForceX = (double*) calloc(numAtoms, sizeof(double));
        tmpForceY = (double*) calloc(numAtoms, sizeof(double));
        tmpForceZ = (double*) calloc(numAtoms, sizeof(double));
        tmpRho = (double*) calloc(numAtoms, sizeof(double));
        tmpStress = (double*) calloc(6, sizeof(double));

        #pragma omp for schedule(static,4)
        for (aatom = 0; aatom < numAtoms; aatom++)
        {
            tmpRho[aatom] = calculateAtomRho(aatom, pos, latVector, rcpVector, atmLabel, numAtoms);
        }

        #pragma omp for schedule(static,4)
        for (aatom = 0; aatom < numAtoms; aatom++)
        {
            //#pragma omp atomic
            wrkRho[aatom] = tmpRho[aatom];
        }

        #pragma omp for schedule(static,4)
        for (aatom = 0; aatom < numAtoms; aatom++)
        {
            ltypea = atmLabel[aatom];
            tmp_energy_m = mbdy.embed(ltypea, ltypea, tmpRho[aatom]);
            #pragma omp atomic
            manyEnergy += tmp_energy_m;
        }


        #pragma omp for schedule(static,4)
        for (aatom = 0; aatom < numAtoms; aatom++)
        {

            ltypea = atmLabel[aatom];

            ax = pos(aatom,0);
            ay = pos(aatom,1);
            az = pos(aatom,2);

        
            for (batom = 0; batom < numAtoms; batom++)
            {
                if (aatom == batom)
                    continue;
        
                ltypeb = atmLabel[batom];

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

                eng_s = 0.0;
                forc_m = 0.0;
                forc_s = 0.0;

                if (rsq <= radius)
                {

                    r = sqrt(rsq);

                    mbdy.calculateManyBodyPairForce(ltypea, ltypeb, r, rsq, eng_s, forc_s);
                    //std::cout << "\n pair force " << aatom << " " << batom << " " << r  << " " << eng_s << " " << forc_s;

                    mbdy.calculateManyBodyForce(ltypea, ltypeb, r, wrkRho[aatom], wrkRho[batom], forc_m);
                    //std::cout << "\n nany body " << wrkRho[aatom] << " " <<  wrkRho[batom] << " " <<  forc_m;
                    //std::cout << "\n total " << batom+1 << "  " << r << "  " << (forc_s - forc_m) / rsq;
                    //exit(-1);

                    forceTot = (forc_s - forc_m) / (rsq * 2.0);

                    tmp_energy_s += 0.5 * eng_s;

                    tmpForceX[aatom] += rx * forceTot;
                    tmpForceY[aatom] += ry * forceTot;
                    tmpForceZ[aatom] += rz * forceTot;
                    tmpForceX[batom] -= rx * forceTot;
                    tmpForceY[batom] -= ry * forceTot;
                    tmpForceZ[batom] -= rz * forceTot;

                    tmpStress[0] -= forceTot * rx * rx;
                    tmpStress[1] -= forceTot * ry * ry;
                    tmpStress[2] -= forceTot * rz * rz;
                    tmpStress[3] -= forceTot * ry * rz;
                    tmpStress[4] -= forceTot * rx * rz;
                    tmpStress[5] -= forceTot * rx * ry;

                }

            } // end of loop over batom

        } // end of loop over aatom
        for (aatom = 0; aatom < numAtoms; aatom++)
        {
            #pragma omp atomic
            force(aatom,0) += tmpForceX[aatom];
        }
        for (aatom = 0; aatom < numAtoms; aatom++)
        {
            #pragma omp atomic
            force(aatom,1) += tmpForceY[aatom];
        }
        for (aatom = 0; aatom < numAtoms; aatom++)
        {
            #pragma omp atomic
            force(aatom,2) += tmpForceZ[aatom];
        }
        for (int i = 0; i < 6; i++)
        {
            #pragma omp atomic
            stress[i] += tmpStress[i];
        }

        free(tmpForceX);
        free(tmpForceY);
        free(tmpForceZ);
        free(tmpRho);
        free(tmpStress);

    } //end of parallel


    twoEnergy = tmp_energy_s;

    /*
    for (aatom = 0; aatom < numAtoms; aatom++)
        std::cout << "frc " << aatom << " " << forceX[aatom] << " " << forceY[aatom] << " " << forceZ[aatom] << std::endl;

    std::cout << "\n two, many mody energy " << twoenergy << " " << manyEnergy << std::endl;
    exit(0);
    */
}

               

#else

void Metal::realSpaceForceImage(double& manyEnergy, double& twoEnergy, const Eigen::MatrixXd &pos, Eigen::VectorXd &stress, const Eigen::MatrixXd &latVector, 
                            const Eigen::MatrixXd &rcpVector, const Eigen::VectorXi &atmLabel, int numAtoms)
{
    int
        batom,
        ltypea,
        ltypeb,
        aatom;

    double
        ax, ay, az,
        bx, by, bz,
        xx, yy, zz,
        fxi, fyi, fzi,
        eng_s,
        tmp_energy_s = 0.0,
        forc_m,
        forc_s,
        forceTot;

    double
        r,
        rsq,
        rx,  // dummy variables for distances
        ry,
        rz,
        radius = pow(shortRangeCut, 2.0);  // real space madelung and short range cut are same
                                           // for image convention so just squared

    manyEnergy = 0.0;
    twoEnergy = 0.0;

    for (aatom = 0; aatom < numAtoms; aatom++)
    {
        ltypea = atmLabel[aatom];
        wrkRho[aatom] = calculateAtomRho(aatom, pos, latVector, rcpVector, atmLabel, numAtoms);
        manyEnergy += mbdy.embed(ltypea, ltypea, wrkRho[aatom]);
    }

    

    for (aatom = 0; aatom < numAtoms; aatom++)
    {

        ltypea = atmLabel[aatom];

        fxi = 0.0;
        fyi = 0.0;
        fzi = 0.0;

        ax = pos(aatom,0);
        ay = pos(aatom,1);
        az = pos(aatom,2);

        
        for (batom = 0; batom < numAtoms; batom++)
        {
            if (aatom == batom)
                continue;
        
            ltypeb = atmLabel[batom];

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

            eng_s = 0.0;
            forc_m = 0.0;
            forc_s = 0.0;

            if (rsq <= radius)
            {

                r = sqrt(rsq);

                mbdy.calculateManyBodyPairForce(ltypea, ltypeb, r, rsq, eng_s, forc_s);
                //std::cout << "\n pair force " << aatom << " " << batom << " " << r  << " " << eng_s << " " << forc_s;
                
                mbdy.calculateManyBodyForce(ltypea, ltypeb, r, wrkRho[aatom], wrkRho[batom], forc_m);
                //std::cout << "\n nany body " << wrkRho[aatom] << " " <<  wrkRho[batom] << " " <<  forc_m;
                //std::cout << "\n total " << batom+1 << "  " << r << "  " << (forc_s - forc_m) / rsq;
                //exit(-1);
                
                forceTot = (forc_s - forc_m) / (rsq * 2.0);
            
                tmp_energy_s += 0.5 * eng_s;

                fxi += rx * forceTot;
                fyi += ry * forceTot;
                fzi += rz * forceTot;
                force(batom,0) -= rx * forceTot;
                force(batom,1) -= ry * forceTot;
                force(batom,2) -= rz * forceTot;

                stress[0] -= forceTot * rx * rx;
                stress[1] -= forceTot * ry * ry;
                stress[2] -= forceTot * rz * rz;
                stress[3] -= forceTot * ry * rz;
                stress[4] -= forceTot * rx * rz;
                stress[5] -= forceTot * rx * ry;

            }

        } // end of loop over batom

        force(aatom,0) += fxi;
        force(aatom,1) += fyi;
        force(aatom,2) += fzi;

        //std::cout << "\n force " << aatom+1 << " " << fxi << " " << fyi << " " << fzi;
 
    } // end of loop over aatom

    twoEnergy = tmp_energy_s;
    
}
#endif

void Metal::realSpaceForceImageNbrList(double& manyEnergy, double& twoEnergy, const Eigen::MatrixXd &pos, Eigen::VectorXd &stress, const Eigen::MatrixXd &latVector, 
                            const Eigen::MatrixXd &rcpVector, const Eigen::VectorXi &atmLabel, 
                            int numAtoms)
{
    int
        batom,
        ltypea,
        ltypeb;

    int
        idx;

    double
        ax, ay, az,
        bx, by, bz,
        xx, yy, zz,
        fxi, fyi, fzi,
        eng_s,
        tmp_energy_s = 0.0,
        forc_m,
        forc_s,
        forceTot;

    double
        r,
        rsq,
        rx,  // dummy variables for distances
        ry,
        rz,
        radius = pow(shortRangeCut, 2.0);  // real space madelung and short range cut are same
                                           // for image convention so just square

    manyEnergy = 0.0;
    twoEnergy = 0.0;

    for (int aatom = 0; aatom < numAtoms; aatom++)
    {
        ltypea = atmLabel[aatom];
        wrkRho[aatom] = calculateAtomRho(aatom, pos, latVector, rcpVector, atmLabel, numAtoms);
        manyEnergy += mbdy.embed(ltypea, ltypea, wrkRho[aatom]);
    }

    for (int aatom = 0; aatom < numAtoms; aatom++)
    {

        ltypea = atmLabel[aatom];
        idx = nbrList.maxNbrs * aatom;

        ax = pos(aatom,0);
        ay = pos(aatom,1);
        az = pos(aatom,2);

        fxi = 0.0;
        fyi = 0.0;
        fzi = 0.0;
        
        for (int nbr = 0; nbr < nbrList.numNbrs[aatom]; nbr++)
        {
        
            batom = nbrList.map[idx+nbr];
            ltypeb = atmLabel[batom];

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

            eng_s = 0.0;
            forc_m = 0.0;
            forc_s = 0.0;

            if (rsq <= radius)
            {

                r = sqrt(rsq);

                mbdy.calculateManyBodyPairForce(ltypea, ltypeb, r, rsq, eng_s, forc_s);
                //std::cout << "\n pair force " << aatom << " " << batom << " " << r  << " " << eng_s << " " << forc_s;
                
                mbdy.calculateManyBodyForce(ltypea, ltypeb, r, wrkRho[aatom], wrkRho[batom], forc_m);
                //std::cout << "\n nany body " << wrkRho[aatom] << " " <<  wrkRho[batom] << " " <<  forc_m;
                //std::cout << "\n total " << batom+1 << "  " << r << "  " << (forc_s - forc_m) / rsq;
                //exit(-1);
                
                forceTot = (forc_s - forc_m) / (rsq * 2.0);
            
                tmp_energy_s += 0.5 * eng_s;

                fxi += rx * forceTot;
                fyi += ry * forceTot;
                fzi += rz * forceTot;
                force(batom,0) -= rx * forceTot;
                force(batom,1) -= ry * forceTot;
                force(batom,2) -= rz * forceTot;

                stress[0] -= forceTot * rx * rx;
                stress[1] -= forceTot * ry * ry;
                stress[2] -= forceTot * rz * rz;
                stress[3] -= forceTot * ry * rz;
                stress[4] -= forceTot * rx * rz;
                stress[5] -= forceTot * rx * ry;

            }

        } // end of loop over batom

        force(aatom,0) += fxi;
        force(aatom,1) += fyi;
        force(aatom,2) += fzi;

    } // end of loop over aatom

    twoEnergy = tmp_energy_s;

}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// functions to calculate energies using repeting cell rather than image convention
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Metal::realSpaceForceNoImage(double& manyEnergy, double& twoEnergy, const Eigen::MatrixXd &pos, Eigen::VectorXd &stress, const Eigen::MatrixXd &latVector, 
                            const Eigen::VectorXi &atmLabel, int numAtoms)
{
    int
        ltypea,
        ltypeb;

    int
        nx, ny, nz;

    double
        ax, ay, az,
        bx, by, bz,
        xx, yy, zz,
        eng_s,
        tmp_energy_s = 0.0,
        forc_m,
        forc_s,
        forceTot;

    double
        r,
        rsq,
        rx,  // dummy variables for distances
        ry,
        rz,
        radius = pow(shortRangeCut, 2.0);  // real space madelung and short range cut are same
                                           // for image convention so just square

    manyEnergy = 0.0;
    twoEnergy = 0.0;

    for (int aatom = 0; aatom < numAtoms; aatom++)
    {
        ltypea = atmLabel[aatom];
        wrkRho[aatom] = calculateAtomRhoNoImage(aatom, pos, latVector, atmLabel, numAtoms);
        manyEnergy += mbdy.embed(ltypea, ltypea, wrkRho[aatom]);
    }


    for (int aatom = 0; aatom < numAtoms; aatom++)
    {

        ltypea = atmLabel[aatom];

        ax = pos(aatom,0);
        ay = pos(aatom,1);
        az = pos(aatom,2);

        for (int batom = 0; batom < numAtoms; batom++)
        {

            ltypeb = atmLabel[batom];

            bx = pos(batom,0);
            by = pos(batom,1);
            bz = pos(batom,2);

            for(nx = -cellX; nx <= cellX; nx++)
            {
                for(ny = -cellY; ny <= cellY; ny++)
                {
                    for(nz = -cellZ; nz <= cellZ; nz++)
                    {
                        //calculate distance
                        xx = bx + nx * latVector(0,0) + ny * latVector(0,1) + nz * latVector(0,2);
                        yy = by + nx * latVector(1,0) + ny * latVector(1,1) + nz * latVector(1,2);
                        zz = bz + nx * latVector(2,0) + ny * latVector(2,1) + nz * latVector(2,2);

                        rx = ax - xx;
                        ry = ay - yy;
                        rz = az - zz;

                        rsq = rx * rx + ry * ry + rz * rz;

                        if (rsq <= radius && rsq > 1.0e-3)
                        {

                            r = sqrt(rsq);

                            mbdy.calculateManyBodyPairForce(ltypea, ltypeb, r, rsq, eng_s, forc_s);
                            //std::cout << "\n pair force " << aatom << " " << batom << " " << r  << " " << eng_s << " " << forc_s;
                
                            mbdy.calculateManyBodyForce(ltypea, ltypeb, r, wrkRho[aatom], wrkRho[batom], forc_m);
                            //std::cout << "\n nany body " << wrkRho[aatom] << " " <<  wrkRho[batom] << " " <<  forc_m;
                            //std::cout << "\n total " << batom+1 << "  " << r << "  " << (forc_s - forc_m) / rsq;
                            //exit(-1);
                
                            forceTot = (forc_s - forc_m) / (rsq * 2.0);
            
                            tmp_energy_s += 0.5 * eng_s;

                            force(aatom,0) -= rx * forceTot; // this needs to be optimised
                            force(aatom,1) -= ry * forceTot;
                            force(aatom,2) -= rz * forceTot;
                            force(batom,0) += rx * forceTot;
                            force(batom,1) += ry * forceTot;
                            force(batom,2) += rz * forceTot;

                            stress[0] -= forceTot * rx * rx;
                            stress[1] -= forceTot * ry * ry;
                            stress[2] -= forceTot * rz * rz;
                            stress[3] -= forceTot * ry * rz;
                            stress[4] -= forceTot * rx * rz;
                            stress[5] -= forceTot * rx * ry;

                        }

                    } // nz

                } // ny

            } // nz

        } // end of loop over batom

    } // end of loop over aatom

    twoEnergy = tmp_energy_s;

}


double Metal::calculateAtomRhoNoImage(int aatom, const Eigen::MatrixXd& pos, const Eigen::MatrixXd& latVector, 
                                      const Eigen::VectorXi& atmLabel, int numAtoms)

{
    int
        batom,
        ltypea,
        ltypeb,
        nx, ny, nz;

    double
        ax, ay, az,
        bx, by, bz,
        xx, yy, zz;

    double
        rho = 0.0,
        r,
        rsq,
        rx,  // dummy variables for distances
        ry,
        rz,
        radius = pow(shortRangeCut, 2.0);  // real space madelung and short range cut are same
                                           // for image convention so just square

    ltypea = atmLabel[aatom];

    ax = pos(aatom,0);
    ay = pos(aatom,1);
    az = pos(aatom,2);

    for (batom = 0; batom < numAtoms; batom++)
    {
    
        ltypeb = atmLabel[batom];

        bx = pos(batom,0);
        by = pos(batom,1);
        bz = pos(batom,2);

        for(nx = -cellX; nx <= cellX; nx++)
        {
            for(ny = -cellY; ny <= cellY; ny++)
            {
                for(nz = -cellZ; nz <= cellZ; nz++)
                {
                    //calculate distance
                    xx = bx + nx * latVector(0,0) + ny * latVector(0,1) + nz * latVector(0,2);
                    yy = by + nx * latVector(1,0) + ny * latVector(1,1) + nz * latVector(1,2);
                    zz = bz + nx * latVector(2,0) + ny * latVector(2,1) + nz * latVector(2,2);

                    rx = ax - xx;
                    ry = ay - yy;
                    rz = az - zz;

                    rsq = rx *rx + ry * ry + rz * rz;

                    if (rsq <= radius && rsq > 1.0e-3)
                    {
                        r = sqrt(rsq);

                        rho += mbdy.calculateManyBodyDensEnergy(ltypea, ltypeb, r);
                    }

                }

            }
        }

    }
    
    return rho;
}
double Metal::calculateAtomRho(int aatom, const Eigen::MatrixXd& pos, const Eigen::MatrixXd& latVector, const Eigen::MatrixXd& rcpVector, 
                                      const Eigen::VectorXi& atmLabel, int numAtoms)

{
    int
        ltypea,
        ltypeb;

    double
        ax, ay, az,
        bx, by, bz,
        xx, yy, zz;

    double
        rho = 0.0,
        r,
        rsq,
        rx,  // dummy variables for distances
        ry,
        rz,
        radius = pow(shortRangeCut, 2.0);  // real space madelung and short range cut are same
                                           // for image convention so just square

    ltypea = atmLabel[aatom];

    ax = pos(aatom,0);
    ay = pos(aatom,1);
    az = pos(aatom,2);

    for (int batom = 0; batom < numAtoms; batom++)
    {
        if (aatom == batom)
            continue;

        ltypeb = atmLabel[batom];

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
            r = sqrt(rsq);

            rho += mbdy.calculateManyBodyDensEnergy(ltypea, ltypeb, r);
        } 

    } // end of loop over aatom

    return rho;

}

void Metal::resetSimulationCellNoImage(Eigen::MatrixXd& pos, const Eigen::MatrixXd& latVector, const Eigen::MatrixXd& rcpVector, int natoms)
{

    double
        xx, yy, zz,
        rx, ry, rz;

    for (auto aatom = 0; aatom < natoms; aatom++)
    {

        rx = pos(aatom,0);
        ry = pos(aatom,1);
        rz = pos(aatom,2);

        xx = rx * rcpVector(0,0) + ry * rcpVector(0,1) + rz * rcpVector(0,2);
        yy = rx * rcpVector(1,0) + ry * rcpVector(1,1) + rz * rcpVector(1,2);
        zz = rx * rcpVector(2,0) + ry * rcpVector(2,1) + rz * rcpVector(2,2);

        //reduced coordinates should lie between 0 and 1
        if (xx < 0.0)
            xx += 1.0;
        if (yy < 0.0)
            yy += 1.0;
        if (zz < 0.0)
            zz += 1.0;
        if (xx >= 1.0)
            xx -= 1.0;
        if (yy >= 1.0)
            yy -= 1.0;
        if (zz >= 1.0)
            zz -= 1.0;

        pos(aatom,0) = xx * latVector(0,0) + yy * latVector(0,1) + zz * latVector(0,2);
        pos(aatom,1) = xx * latVector(1,0) + yy * latVector(1,1) + zz * latVector(1,2);
        pos(aatom,2) = xx * latVector(2,0) + yy * latVector(2,1) + zz * latVector(2,2);

    }
}

void Metal::resetSimulationCellMinImage(Eigen::MatrixXd& pos, const Eigen::MatrixXd& latVector, const Eigen::MatrixXd& rcpVector, int natoms)
{
    double
        rx, ry, rz,
        xx, yy, zz;

    for (int i = 0; i < natoms; i++)
    {

        rx = pos(i,0);
        ry = pos(i,1);
        rz = pos(i,2);

        xx = rx * rcpVector(0,0) + ry * rcpVector(0,1) + rz * rcpVector(0,2);
        yy = rx * rcpVector(1,0) + ry * rcpVector(1,1) + rz * rcpVector(1,2);
        zz = rx * rcpVector(2,0) + ry * rcpVector(2,1) + rz * rcpVector(2,2);

        xx -= rint(xx);
        yy -= rint(yy);
        zz -= rint(zz);

        pos(i,0) = xx * latVector(0,0) + yy * latVector(0,1) + zz * latVector(0,2);
        pos(i,1) = xx * latVector(1,0) + yy * latVector(1,1) + zz * latVector(1,2);
        pos(i,2) = xx * latVector(2,0) + yy * latVector(2,1) + zz * latVector(2,2);

    }
}


// ----------------
// Python interface
// ----------------

namespace py = pybind11;

PYBIND11_MODULE(metal, m) {
    py::class_<Metal>(m, "Metal")
        .def(py::init<const std::string &>())
        .def("readPotential", &Metal::readPotential)
        .def("setup", &Metal::setup)
        .def("finalise", &Metal::finalise)
        .def("calculateForces", &Metal::calculateForces);
}