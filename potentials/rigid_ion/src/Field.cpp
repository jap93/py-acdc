

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>

#include <string>
#include <stdio.h>
#include <algorithm>
#include <vector>

#include <Eigen/LU>

#include "Species.h"
#include "Constants.h"

#include "Ewald.h"
#include "NbrListPBC.h"
#include "TwoBody.h"

#ifdef _OPENMP
    #include <omp.h>
#endif

class RigidIon
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
        newjob = true;

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
        cFact = 14.3997584,

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

    Eigen::MatrixXd force; 

    NbrListPBC
        nbrList; 
        
    TwoBody
        vdw;

    Ewald
        coul;

    Species
        spec;

    std::ofstream
        outStream;

    std::string
        name;

    void fieldErfcGen(double alpha);

    void hfunc1(double, double&, double&);

    void hfunc2(double, double&, double&);

    std::vector<std::string> split(std::string s);

    void resetSimulationCellNoImage(Eigen::MatrixXd pos, Eigen::MatrixXd latVector, Eigen::MatrixXd rcpVector, int natoms);

    void resetSimulationCellMinImage(Eigen::MatrixXd pos, Eigen::MatrixXd latVector, Eigen::MatrixXd rcpVector, int natoms);

    void realSpaceForceImage(double& realenergy, double& twoenergy, const Eigen::MatrixXd &posX, Eigen::MatrixXd &stress, const Eigen::MatrixXd &latVector, 
                            const Eigen::MatrixXd &rcpVector, const Eigen::VectorXd &atomCharge, const Eigen::VectorXi &atmLabel, 
                            int numAtoms);

    void realSpaceForceImageNbrList(double& realenergy, double& twoenergy, const Eigen::MatrixXd &posX, Eigen::MatrixXd &stress, const Eigen::MatrixXd &latVector, 
                            const Eigen::MatrixXd &rcpVector, const Eigen::VectorXd &atomCharge, const Eigen::VectorXi &atmLabel, 
                            int numAtoms);

    void realSpaceForceNoImage(double& realenergy, double& twoenergy, const Eigen::MatrixXd &posX, Eigen::MatrixXd &stress, const Eigen::MatrixXd &latVector, 
                            const Eigen::MatrixXd &rcpVector, const Eigen::VectorXd &atomCharge, const Eigen::VectorXi &atmLabel, 
                            int numAtoms);

    //void fieldErfcGen(double alpha);
    //void interpolate_shifted_damped_potential(double, double, double, double&);
    void interpolate_shifted_damped_potential(double, double, double, double&, double&);

public:

    // constructor
    //RigidIon();
    RigidIon(const std::string &name) : name(name) { }
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

    //void calculateForces(Eigen::VectorXd posX, Eigen::VectorXd posY, Eigen::VectorXd posZ, Eigen::VectorXd forceX, Eigen::VectorXd forceY, Eigen::VectorXd forceZ, 
    //                     Eigen::VectorXd stress, Eigen::VectorXd latVector, Eigen::VectorXd rcpVector, Eigen::VectorXd atomcharge, int* intlabel, int* frozen, int numatoms, Energy& eng);
    Eigen::MatrixXd &calculateForces(const Eigen::MatrixXd &posX, Eigen::MatrixXd &stress, const Eigen::MatrixXd &latVector, 
                            const Eigen::MatrixXd &rcpVector, const Eigen::VectorXd &atomCharge, const Eigen::VectorXi &atmLabel, const Eigen::VectorXi &frozen,
                            int numAtoms);
};



//RigidIon::RigidIon()
//{
//    newjob = true;
//}

//RigidIon::~RigidIon()
//{
   
//}

/****************************************************************************************
 create error function lookup table
****************************************************************************************/
void RigidIon::fieldErfcGen(double alpha)
{

    int
        i;

    double
        tt, exp1,
        rrr, rsq;

    double
        //      epsq = job.dielec,
        alph = alpha,
        a1 = 0.254829592,
        a2 = -0.284496736,
        a3 = 1.421413741,
        a4 = -1.453152027,
        a5 = 1.061405429,
        pp = 0.3275911;

    if (alpha == 0.0)
        return;

    erfc_e.resize(MAXMESH);
    erfc_f.resize(MAXMESH);

    delMesh = shortRangeCut / (MAXMESH - 4);      // spacing between mesh points
    rMesh = (MAXMESH - 4) / shortRangeCut;

    for (i = 1; i < MAXMESH; i++)
    {

        rrr = double(i) * delMesh;
        rsq = rrr * rrr;

        tt = 1.0 / (1.0 + pp * alph * rrr);
        exp1 = exp(-pow((alph * rrr), 2.0));

        erfc_e[i] = tt * (a1 + tt * (a2 + tt * (a3 + tt * (a4 + tt * a5)))) * exp1 / rrr;
        erfc_f[i] = (erfc_e[i] + 2.0 * (alph / ROOTPI) * exp1) / rsq;
    }

    aa = erfc_f[MAXMESH - 4] * shortRangeCut;
    bb = -(erfc_e[MAXMESH - 4] + aa * shortRangeCut);

    //  b0 = 2.0 * (epsq - 1.0) / (2.0 * epsq + 1.0);
    //  rfld0 = b0 / pow(job.shortrangecut,3);
    //  rfld1 = (1.0 + 0.5 * b0) / job.shortrangecut;
    rfld0 = erfc_e[MAXMESH - 4];
    rfld1 = shortRangeCut * erfc_f[MAXMESH - 4];

    //cout << setiosflags(ios::showpoint | ios::fixed | ios::right)
    //   << setprecision(15) << "\n gen erf mesh " << job.shortrangecut << " " << rfld0  << " " << rfld1;
    //cout.flush();
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// utility functions
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void RigidIon::setup(const Eigen::VectorXd &cellProp, const Eigen::MatrixXd &rcpVector, const Eigen::VectorXd &rcpProp, double minDimension, double volume, 
                  int numAtoms)
{


    double
        eps = 0.0,
        tol = 0.0;

    if (shortRangeCut < 1.0)
    {
        outStream << "\n\n\n *** short range cutoff is too small : " << shortRangeCut << std::endl;
        outStream.flush();
        exit(EXIT_FAILURE);
    }

    // check the image convention is obeyed
    if (minimumImage && (shortRangeCut + verletShell) > (minDimension / 2))
    {
        outStream << "\n\n\n *** short range shortRangeCut to large for simulation box." << std::endl;
        outStream << " cutoff : " << shortRangeCut << " minimum dimension: " << minDimension << std::endl;
        outStream.flush();
        exit(EXIT_FAILURE);
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
        outStream << "\n\n The minimum image convention is being used!" << std::endl;
    }

    if (doEwald > 0)
    {

        if (doEwald == 2)
        {
            eta = alpha;
            fieldErfcGen(eta);
        }
        else
        {
            eps = fmin(fabs(madelungAcc), 0.5);
            tol = sqrt(fabs(log(eps * shortRangeCut)));

            eta = sqrt(fabs(log(eps * shortRangeCut * tol))) / shortRangeCut;

        }

    }
    

    if (newjob == true)
    {
        outStream << "\n cutoffs for lattice summation :"
                  << "\n real space shortRangeCut       = " << shortRangeCut << " Angstroms" << std::endl;

        outStream << "\n cutoffs for lattice summation eta      = " << eta << std::endl;
    }

    if (newjob == true && doEwald  == 1)
    {

        coul.setEwaldImage(eta, shortRangeCut, eps, tol, cellProp, rcpProp);
        coul.estimateGvector(rcpVector);
        coul.allocateArrays();
        // calculate g vectors
        coul.findGvector(volume, rcpVector);

        if (newjob == true)
        {
            //coul.printEwald(outStream);
            //coul.printNumGVec(outStream);
        }

    }

    // set up short-range potential mesh
    if (newjob == true) vdw.calculateEnergyMesh(spec, shortRangeCut, outStream);

    force.resize(numAtoms+1,3);

    newjob = false;

}


std::vector<std::string> RigidIon::split(std::string s)
{
    std::vector<std::string> words;

    std::stringstream ss(s);
    std::string word;
    
    while (ss >> word) 
    {
        words.push_back(word);
    }

    return words;
}


/******************************************************
reads in description of the species and the potential model
*******************************************************/
void RigidIon::readPotential(const std::string &fileName)
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

    #ifdef _OPENMP
    #pragma omp parallel
    {
        outStream << "\n the number of threads: " << omp_get_num_threads() << std::endl;
    }
    #endif

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
        else if (keyWord == "ewald")
		{
			subWord = words[1];

			if (subWord == "precis") 
			{
		        doEwald = 1;
				dummy = words[2];
                madelungAcc = std::stod(dummy);
			}
			else if (subWord == "damp")
			{
		        doEwald = 2;
				dummy = words[2];
                alpha = std::stod(dummy);
			}
			else if (subWord == "none") 
			{
				doEwald = 0;
			}

		}
		else if (keyWord == "noimage")
		{
			minimumImage = false;
		}
        else if (keyWord == "twobody")
        {

            dummy = words[1];
            num = std::stoi(dummy);
            
            vdw.loadPotential(inStream, outStream, num);
            vdw.printPotential(outStream);

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
Eigen::MatrixXd &RigidIon::calculateForces(const Eigen::MatrixXd &pos, Eigen::MatrixXd &stress, const Eigen::MatrixXd &latVector, 
                            const Eigen::MatrixXd &rcpVector, const Eigen::VectorXd &atomCharge, const Eigen::VectorXi &atmLabel, const Eigen::VectorXi &frozen,
                            int numAtoms)
{
    double
        realenergy = 0.0,
        rcpenergy = 0.0,
        twoenergy = 0.0;

    for (auto i = 0; i < numAtoms; i++)
    {
        //std::cout << atomCharge[i] << " " << atmLabel[i] << std::endl;
        force(i,0) = 0.0;
        force(i,1) = 0.0;
        force(i,2) = 0.0;
    }

    for (auto i = 0; i < 6; i++)
        stress(i) = 0.0;


    if (doEwald == 1)   //calculate reciprocal space part of ewald
    {

        rcpenergy = coul.recipSpaceForce(pos, force, stress, atomCharge, numAtoms);

    }
    if (minimumImage)
    {
        resetSimulationCellMinImage(pos, latVector, rcpVector, numAtoms);

        if (useNbrList)
        {
            //std::cout << "\n should not be here";
            exit(EXIT_FAILURE);
            realSpaceForceImageNbrList(realenergy, twoenergy, pos,  stress, latVector, rcpVector, 
                                    atomCharge, atmLabel, numAtoms);
        }
        else
        {
            realSpaceForceImage(realenergy, twoenergy, pos, stress, latVector, rcpVector, 
                                    atomCharge, atmLabel, numAtoms);
        }
               
    }
    else
    {
            resetSimulationCellNoImage(pos, latVector, rcpVector, numAtoms);
            realSpaceForceNoImage(realenergy, twoenergy, pos, stress, latVector, rcpVector, 
                                  atomCharge, atmLabel, numAtoms);
    }

    //tag the energies on to the force data to send back
    force(numAtoms,0) = rcpenergy;
    force(numAtoms,1) = realenergy;
    force(numAtoms,2) = twoenergy;
    //std::cout << "\n energies " << rcpenergy << " " << realenergy << " " << twoenergy;
    return force;

}
 
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// routines using the nearest image convention
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



/**************************************************************
 * calculate energy and forces using image convention- this subroutine
 * uses explicit calculation of the energy and forces
 *
 * before doing the calculation it is necessary to
 * to reset the box so that any atom that has
 * migrated out of the box can be placed back in according
 * to the image convention
 **************************************************************/
#ifdef _OPENMP

void RigidIon::realSpaceForceImage(double& realenergy, double& twoenergy, const Eigen::MatrixXd &pos, Eigen::MatrixXd &stress, const Eigen::MatrixXd &latVector, 
                            const Eigen::MatrixXd &rcpVector, const Eigen::VectorXd &atomCharge, const Eigen::VectorXi &atmLabel, int numAtoms)
{
    int
        ltypea,
        ltypeb;

    int
        batom,
        aatom;

    double
        ax, ay, az,
        bx, by, bz,
        xx, yy, zz,
        h0,
        h1,
        eng,
        eng_s,
        tmp_self = 0.0,
        tmp_energy_r = 0.0,
        tmp_energy_s = 0.0,
        forc,
        forc_s,
        chargea,
        chargeb,
        chargeprod,
        hfct0 = eta,
        hfct1 = hfct0 * hfct0 * hfct0;

    double
        r,
        rsq,
        rx,  // dummy variables for distances
        ry,
        rz,
        radius = pow(shortRangeCut, 2.0);  // real space madelung and short range cut are same
                                           // for image convention so just square

    realenergy = 0.0;
    twoenergy = 0.0;

    if (doEwald == 1)
    {

        for (int aatom = 0; aatom < numAtoms; aatom++)
        {

            hfunc2(0.0, h0, h1);

            chargea = atomCharge[aatom];
            chargeprod = 0.5 * cFact * chargea * chargea;
            tmp_self += -h0 * chargeprod * hfct0;

        }

    }

    //Eigen::VectorXd tmpForceX; //[numAtoms];
    //Eigen::VectorXd tmpForceY; //[numAtoms];
    //Eigen::VectorXd tmpForceZ; //[numAtoms];
    double* tmpForceX = nullptr;
    double* tmpForceY = nullptr;
    double* tmpForceZ = nullptr;

    #pragma omp parallel default(none)  \
        reduction(+: tmp_energy_r, tmp_energy_s) \
        shared(atmLabel, atomCharge, pos, numAtoms, rcpVector, latVector, force, radius, shortRangeCut, hfct0, hfct1, doEwald) \
        private(ltypea, ltypeb, chargea, chargeb, chargeprod, ax, ay, az, bx, by, bz, rx, ry, rz, r, rsq, xx, yy, zz, aatom, batom) \
        private(h0, h1, eng, eng_s, forc, forc_s, tmpForceX, tmpForceY, tmpForceZ)
    {

        //tmpForceX.resize(numAtoms); //calloc(numAtoms, sizeof(double));
        //tmpForceY.resize(numAtoms); //calloc(numAtoms, sizeof(double));
        //tmpForceZ.resize(numAtoms); //calloc(numAtoms, sizeof(double));
        tmpForceX = (double*) calloc(numAtoms, sizeof(double));
        tmpForceY = (double*) calloc(numAtoms, sizeof(double));
        tmpForceZ = (double*) calloc(numAtoms, sizeof(double));
        #pragma omp for schedule(static,4)
        for (aatom = 0; aatom < numAtoms; aatom++)
        {

            ltypea = atmLabel[aatom];
            chargea = cFact * atomCharge[aatom];

            ax = pos(aatom,0);
            ay = pos(aatom,1);
            az = pos(aatom,2);


            for (batom = aatom + 1; batom < numAtoms; batom++)
            {

                chargeb = atomCharge[batom];
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

                eng = 0.0;
                eng_s = 0.0;
                forc = 0.0;
                forc_s = 0.0;

                if (rsq <= radius)
                {

                    r = sqrt(rsq);

                    h0 = 0.0;
                    h1 = 0.0;

                    chargeprod = chargea * chargeb;

                    if (doEwald == 1)
                    {
                        hfunc1(hfct0 * r, h0, h1);

                        eng = h0 * chargeprod * hfct0;
                        forc = h1 * chargeprod * hfct1;
                    }
                    else if (doEwald == 2)
                    {
                        //evaluate damped and shifted
                        interpolate_shifted_damped_potential(r, chargeprod, shortRangeCut, eng, forc);

                    }

                    if (r <= shortRangeCut && r > 1.0e-4)
                    {
                        vdw.shortrangeone(ltypea, ltypeb, r, rsq, eng_s, forc_s);
                    }

                }

                tmp_energy_r += eng;
                tmp_energy_s += eng_s;

                forc -= forc_s;

                tmpForceX[aatom] -= rx * forc;
                tmpForceY[aatom] -= ry * forc;
                tmpForceZ[aatom] -= rz * forc;
                tmpForceX[batom] += rx * forc;
                tmpForceY[batom] += ry * forc;
                tmpForceZ[batom] += rz * forc;

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
            force(aatom, 1) += tmpForceY[aatom];
        }
        for (aatom = 0; aatom < numAtoms; aatom++)
        {
            #pragma omp atomic
            force(aatom, 2) += tmpForceZ[aatom];
        }
        free(tmpForceX);
        free(tmpForceY);
        free(tmpForceZ);

    }  // end of parallel

    //for (aatom = 0; aatom < numAtoms; aatom++)
    //{
    //    std::cout << "\n force rreal " << aatom << " " << forceX[aatom] << " " << forceY[aatom] << " " << forceZ[aatom];
    //}
    realenergy = tmp_energy_r + tmp_self;
    twoenergy = tmp_energy_s;

}

#else

void RigidIon::realSpaceForceImage(double& realenergy, double& twoenergy, const Eigen::MatrixXd &pos, Eigen::MatrixXd &stress, const Eigen::MatrixXd &latVector, 
                            const Eigen::MatrixXd &rcpVector, const Eigen::VectorXd &atomCharge, const Eigen::VectorXi &atmLabel, int numAtoms)
{
    int
        ltypea,
        ltypeb;

    int
        batom,
        idx,
        aatom;

    double
        ax, ay, az,
        bx, by, bz,
        xx, yy, zz,
        h0,
        h1,
        eng,
        eng_s,
        tmp_energy_r = 0.0,
        tmp_energy_s = 0.0,
        forc,
        forc_s,
        chargea,
        chargeb,
        chargeprod,
        hfct0 = eta,
        hfct1 = pow(hfct0, 3.0);

    double
        r,
        rsq,
        rx,  // dummy variables for distances
        ry,
        rz,
        radius = pow(shortRangeCut, 2.0);  // real space madelung and short range cut are same
                                           // for image convention so just square

    realenergy = 0.0;
    twoenergy = 0.0;


    if (doEwald == 1)
    {

        //for (int aatom = rank; aatom < numAtoms; aatom = aatom + numProcs)
        for (int aatom = 0; aatom < numAtoms; aatom++)
        {

            hfunc2(0.0, h0, h1);

            chargea = atomCharge[aatom];
            chargeprod = 0.5 * cFact * chargea * chargea;
            tmp_energy_r += -h0 * chargeprod * hfct0;
            //std::cout << "\n" << chargeprod << " " << h0 << " " << hfct0 << " " << tmp_energy_r;

        }

    }

    for (int aatom = 0; aatom < numAtoms; aatom++)
    {

        ltypea = atmLabel[aatom];
        chargea = cFact * atomCharge[aatom];
        idx = nbrList.maxNbrs * aatom;

        ax = pos(aatom,0);
        ay = pos(aatom,1);
        az = pos(aatom,2);

        
        for (batom = aatom + 1; batom < numAtoms; batom++)
        {
        
            chargeb = atomCharge[batom];
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

            eng = 0.0;
            eng_s = 0.0;
            forc = 0.0;
            forc_s = 0.0;

            if (rsq <= radius)
            {

                r = sqrt(rsq);

                h0 = 0.0;
                h1 = 0.0;

                chargeprod = chargea * chargeb;

                if (doEwald == 1)
                {
                    hfunc1(hfct0 * r, h0, h1);

                    eng = h0 * chargeprod * hfct0;
                    forc = h1 * chargeprod * hfct1;
                    //std::cout << "\n" << ltypea << " " << ltypeb << " " << r << " " << eng << " " << forc;
                }
                else if (doEwald == 2)
                {
                    //evaluate damped and shifted
                    interpolate_shifted_damped_potential(r, chargeprod, shortRangeCut, eng, forc);

                }

                if (r <= shortRangeCut && r > 1.0e-4)
                {
                    vdw.shortrangeone(ltypea, ltypeb, r, rsq, eng_s, forc_s);
                    //std::cout << "\n" << ltypea << " " << ltypeb << " " << r << " " << eng_s << " " << forc_s;
                }

            }

            tmp_energy_r += eng;
            tmp_energy_s += eng_s;

            forc -= forc_s;

            force(aatom,0) -= rx * forc; // this needs to be optimised
            force(aatom,1) -= ry * forc;
            force(aatom,2) -= rz * forc;
            force(batom,0) += rx * forc;
            force(batom,1) += ry * forc;
            force(batom,2) += rz * forc;

            stress(0) += forc * rx * rx;
            stress(1) += forc * ry * ry;
            stress(2) += forc * rz * rz;
            stress(3) += forc * ry * rz;
            stress(4) += forc * rx * rz;
            stress(5) += forc * rx * ry;

        } // end of loop over batom

    } // end of loop over aatom

    realenergy = tmp_energy_r;
    twoenergy = tmp_energy_s;

}

#endif

void RigidIon::realSpaceForceImageNbrList(double& realenergy, double& twoenergy, const Eigen::MatrixXd &pos, Eigen::MatrixXd &stress, const Eigen::MatrixXd &latVector, 
                            const Eigen::MatrixXd &rcpVector, const Eigen::VectorXd &atomCharge, const Eigen::VectorXi &atmLabel, int numAtoms)
{
    int
        ltypea,
        ltypeb;

    int
        batom,
        idx;

    double
        ax, ay, az,
        bx, by, bz,
        xx, yy, zz,
        h0,
        h1,
        eng,
        eng_s,
        tmp_energy_r = 0.0,
        tmp_energy_s = 0.0,
        forc,
        forc_s,
        chargea,
        chargeb,
        chargeprod,
        hfct0 = eta,
        hfct1 = pow(hfct0, 3.0);

    double
        r,
        rsq,
        rx,  // dummy variables for distances
        ry,
        rz,
        radius = pow(shortRangeCut, 2.0);  // real space madelung and short range cut are same
                                           // for image convention so just square

    realenergy = 0.0;
    twoenergy = 0.0;


    if (doEwald == 1)
    {

        //for (int aatom = rank; aatom < numAtoms; aatom = aatom + numProcs)
        for (int aatom = 0; aatom < numAtoms; aatom++)
        {

            hfunc2(0.0, h0, h1);

            chargea = atomCharge[aatom];
            chargeprod = 0.5 * cFact * chargea * chargea;
            tmp_energy_r += -h0 * chargeprod * hfct0;

        }

    }

    for (int aatom = 0; aatom < numAtoms; aatom++)
    {

        ltypea = atmLabel[aatom];
        chargea = cFact * atomCharge[aatom];
        idx = nbrList.maxNbrs * aatom;

        ax = pos(aatom,0);
        ay = pos(aatom,1);
        az = pos(aatom,2);

        
        for (int nbr = 0; nbr < nbrList.numNbrs[aatom]; nbr++)
        {
        
            batom = nbrList.map[idx+nbr];
            chargeb = atomCharge[batom];
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

            eng = 0.0;
            eng_s = 0.0;
            forc = 0.0;
            forc_s = 0.0;

            if (rsq <= radius)
            {

                r = sqrt(rsq);

                h0 = 0.0;
                h1 = 0.0;

                chargeprod = chargea * chargeb;

                if (doEwald == 1)
                {
                    hfunc1(hfct0 * r, h0, h1);

                    eng = h0 * chargeprod * hfct0;
                    forc = h1 * chargeprod * hfct1;
                }
                else if (doEwald == 2)
                {
                    //evaluate damped and shifted
                    interpolate_shifted_damped_potential(r, chargeprod, shortRangeCut, eng, forc);

                }

                if (r <= shortRangeCut && r > 1.0e-4)
                {
                    vdw.shortrangeone(ltypea, ltypeb, r, rsq, eng_s, forc_s);
                }

            }

            tmp_energy_r += eng;
            tmp_energy_s += eng_s;

            forc -= forc_s;

            force(aatom,0) -= rx * forc; // this needs to be optimised
            force(aatom,1) -= ry * forc;
            force(aatom,2) -= rz * forc;
            force(batom,0) += rx * forc;
            force(batom,1) += ry * forc;
            force(batom,2) += rz * forc;

            stress(0) += forc * rx * rx;
            stress(1) += forc * ry * ry;
            stress(2) += forc * rz * rz;
            stress(3) += forc * ry * rz;
            stress(4) += forc * rx * rz;
            stress(5) += forc * rx * ry;

        } // end of loop over batom

    } // end of loop over aatom

    realenergy = tmp_energy_r;
    twoenergy = tmp_energy_s;

}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// functions to calculate energies using repeting cell rather than image convention
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void RigidIon::realSpaceForceNoImage(double& realenergy, double& twoenergy, const Eigen::MatrixXd &pos, Eigen::MatrixXd &stress, const Eigen::MatrixXd &latVector, 
                            const Eigen::MatrixXd &rcpVector, const Eigen::VectorXd &atomCharge, const Eigen::VectorXi &atmLabel, int numAtoms)
{
    int
        ltypea,
        ltypeb;

    double
        ax, ay, az,
        bx, by, bz,
        xx, yy, zz,
        h0,
        h1,
        eng,
        eng_s,
        tmp_energy_r = 0.0,
        tmp_energy_s = 0.0,
        forc,
        forc_s,
        chargea,
        chargeb,
        chargeprod,
        hfct0 = eta,
        hfct1 = pow(hfct0, 3.0);

    double
        r,
        rsq,
        rx,  // dummy variables for distances
        ry,
        rz,
        radius = pow(shortRangeCut, 2.0);  // real space madelung and short range cut are same
                                           // for image convention so just square

    bool
        sameCell = false;

    realenergy = 0.0;
    twoenergy = 0.0;

    if (doEwald == 1)
    {

        //for (int aatom = rank; aatom < numAtoms; aatom = aatom + numProcs)
        for (auto aatom = 0; aatom < numAtoms; aatom++)
        {

            hfunc2(0.0, h0, h1);

            chargea = atomCharge[aatom];
            chargeprod = 0.5 * cFact * chargea * chargea;
            tmp_energy_r += -h0 * chargeprod * hfct0;

        }

    }

    //for (aatom = rank; aatom < numAtoms - 1; aatom = aatom + numProcs)
    for (auto aatom = 0; aatom < numAtoms; aatom++)
    {

        ltypea = atmLabel[aatom];
        chargea = cFact * atomCharge[aatom];

        ax = pos(aatom,0);
        ay = pos(aatom,1);
        az = pos(aatom,2);

        for (auto batom = aatom; batom < numAtoms; batom++)
        {

            chargeb = atomCharge[batom];
            ltypeb = atmLabel[batom];

            bx = pos(batom,0);
            by = pos(batom,1);
            bz = pos(batom,2);
            
            for(auto nx = -cellX; nx <= cellX; nx++)
            {
                for(auto ny = -cellY; ny <= cellY; ny++)
                {
                    for(auto nz = -cellZ; nz <= cellZ; nz++)
                    {
                        sameCell = false;
                        if (nx == 0 && ny == 0 && nz == 0)
                            sameCell = true;

                        if (aatom == batom && sameCell == true)
                            continue;

                        //calculate distance
                        xx = bx + nx * latVector(0,0) + ny * latVector(0,1) + nz * latVector(0,2);
                        yy = by + nx * latVector(1,0) + ny * latVector(1,1) + nz * latVector(1,2);
                        zz = bz + nx * latVector(2,0) + ny * latVector(2,1) + nz * latVector(2,2);
                        
                        rx = ax - xx;
                        ry = ay - yy;
                        rz = az - zz;

                        rsq = rx *rx + ry * ry + rz * rz;

                        eng = 0.0;
                        eng_s = 0.0;
                        forc = 0.0;
                        forc_s = 0.0;

                        if (rsq <= radius)
                        {

                            r = sqrt(rsq);

                            h0 = 0.0;
                            h1 = 0.0;

                            chargeprod = chargea * chargeb;

                            if (doEwald == 1)
                            {
                                hfunc1(hfct0 * r, h0, h1);

                                eng = h0 * chargeprod * hfct0;
                                forc = h1 * chargeprod * hfct1;
                            }
                            else if (doEwald == 2)
                            {
                                //evaluate damped and shifted
                                interpolate_shifted_damped_potential(r, chargeprod, shortRangeCut, eng, forc);

                            }

                            vdw.shortrangeone(ltypea, ltypeb, r, rsq, eng_s, forc_s);

                            tmp_energy_r += eng;
                            tmp_energy_s += eng_s;

                            forc -= forc_s;

                            force(aatom,0) -= rx * forc; // this needs to be optimised
                            force(aatom,1) -= ry * forc;
                            force(aatom,2) -= rz * forc;
                            force(batom,0) += rx * forc;
                            force(batom,1) += ry * forc;
                            force(batom,2) += rz * forc;

                            stress(0) += forc * rx * rx;
                            stress(1) += forc * ry * ry;
                            stress(2) += forc * rz * rz;
                            stress(3) += forc * ry * rz;
                            stress(4) += forc * rx * rz;
                            stress(5) += forc * rx * ry;
                        }

                    } // nz

                } // ny

            } // nz

        } // end of loop over batom

    } // end of loop over aatom

    // sum the forces over processors here - not the best way but easiest for now
    realenergy = tmp_energy_r;
    twoenergy = tmp_energy_s;

}

void inline RigidIon::interpolate_shifted_damped_potential(double r, double chg, double shortRangeCut, double& e0, double& f1)
{
    double
        tmp1,
        tmp2,
        tmp3,
        d1,
        d2,
        dlt,
        omega,
        gamma,
        rstep = 1.0 / delMesh;

    int
        point = int(r * rstep);

    dlt = r * rstep - double(point);

    e0 = f1 = 0.0;

    //interpolate energy
    tmp1 = erfc_e[point];
    tmp2 = erfc_e[point + 1];
    tmp3 = erfc_e[point + 2];

    d1 = tmp1 + (tmp2 - tmp1) * dlt;
    d2 = tmp2 + (tmp3 - tmp2) * (dlt - 1.0);

    omega = d1 + (d2 - d1) * dlt * 0.5;

    e0 = chg * (omega - rfld0 + rfld1 * (r - shortRangeCut));

    //interpolate force
    tmp1 = erfc_f[point];
    tmp2 = erfc_f[point + 1];
    tmp3 = erfc_f[point + 2];

    d1 = tmp1 + (tmp2 - tmp1) * dlt;
    d2 = tmp2 + (tmp3 - tmp2) * (dlt - 1.0);

    gamma = d1 + (d2 - d1) * dlt * 0.5;

    f1 = -chg * (gamma - rfld1 / r);

}

void RigidIon::resetSimulationCellNoImage(Eigen::MatrixXd pos, Eigen::MatrixXd latVector, Eigen::MatrixXd rcpVector, int natoms)
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

void RigidIon::resetSimulationCellMinImage(Eigen::MatrixXd pos, Eigen::MatrixXd latVector, Eigen::MatrixXd rcpVector, int natoms)
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


/*************************************************************
 * function for ewald sum
 * **********************************************************/
void inline RigidIon::hfunc1(double x, double& h0, double& h1)
{
    double
        fact = 1.128379167095512,
        expon,
        xsq;

    if (x == 0.0)
        return;

    h0 = erfc(x);
    h0 = h0 / x;

    xsq = x * x;
    expon = fact * exp(-xsq);
    h1 = (-h0 - expon) / xsq;

    return;
}

void inline RigidIon::hfunc2(double x, double& h0, double& h1)
{
    double
        fact = 1.128379167095512,
        expon,
        xsq;
    xsq = x * x;

    if (xsq < 0.01)
    {
        h0 = fact * ((0.1 * xsq - 0.33333333333) * xsq + 1.0);
        h1 = fact * ((-xsq * 0.1428571429 + 0.4) * xsq - 0.6666666666667);
    }
    else
    {
        h0 = erf(x);
        h0 = h0 / x;

        expon = fact * exp(-xsq);
        h1 = (-h0 + expon) / xsq;

    }

    return;
}

// ----------------
// Python interface
// ----------------

namespace py = pybind11;

PYBIND11_MODULE(rigidion, m) {
    py::class_<RigidIon>(m, "RigidIon")
        .def(py::init<const std::string &>())
        .def("readPotential", &RigidIon::readPotential)
        .def("setup", &RigidIon::setup)
        .def("finalise", &RigidIon::finalise)
        .def("calculateForces", &RigidIon::calculateForces);
}