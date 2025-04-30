#include "Ewald.h"

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// GPU kernals
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////
// the C++ code starts here
//////////////////////////////////////////////////////////////////////

Ewald::Ewald()
{
    rcpSpaceCut = 0.0;

    gcellx = gcelly = gcellz = 0;

}

Ewald::~Ewald()
{

    
}

void Ewald::allocateArrays(void)
{
    gx.resize(maxGVec);
    gy.resize(maxGVec);
    gz.resize(maxGVec);

    g0.resize(maxGVec);
    g1.resize(maxGVec);

    rcpSum.resize(maxGVec);
    
}

/* error function */

void Ewald::printNumGVec(std::ofstream& outStream)
{
    
    outStream << "\n the number of reciprocal lattice vectors " << numGVec << std::endl;
}



/* *************************************************************
 * calculates the ewald parameters usiming the image convention
 *
 * input: accuracy required for madelung sum
 ***************************************************************/

void Ewald::setEwaldImage(double eta, double cutoff, double eps, double tol, const Eigen::VectorXd &cellProp, 
                          const Eigen::VectorXd &rcpProp)
{
    double
        cellx,                        // cell dimensions in x y + z
        celly,
        cellz,
        tol1;

    etaMad = eta;

    /* determines the maximum number of k vectors
     * in the ewald summation */
    
    tol1 = sqrt(-log(eps * cutoff * pow(2.0 * tol * eta, 2)));

    gcellx = rint(0.25 + cellProp[0] * eta * tol1 / PI);
    gcelly = rint(0.25 + cellProp[1] * eta * tol1 / PI);
    gcellz = rint(0.25 + cellProp[2] * eta * tol1 / PI);

    cellx = double(gcellx) * rcpProp[6];
    celly = double(gcelly) * rcpProp[7];
    cellz = double(gcellz) * rcpProp[8];

    rcpSpaceCut = std::min(cellx, celly);
    rcpSpaceCut = std::min(rcpSpaceCut, cellz);

    return;

}


void Ewald::printEwald(std::ofstream& outStream)
{
    
    outStream << "\n cutoffs for lattice summation:"
              << "\n eta                     = " << etaMad
              << "\n reciprocal space cutoff = " << rcpSpaceCut << " angstroms^-1" << std::endl;

    outStream << "\n maximum k-vectors "
              << std::setw(4) << gcellx
              << std::setw(4) << gcelly
              << std::setw(4) << gcellz << std::endl;

    outStream.flush();

}

void Ewald::estimateGvector(const Eigen::MatrixXd &rcpVector)
{
    int
        nx,
        ny,
        nz;

    double
        x,
        y,
        z,
        gsq = 0.0,
        rcpcutsq = rcpSpaceCut * rcpSpaceCut;

    maxGVec = 0; // zero number of vectors

    for (nx = 0; nx <= gcellx; nx++)
    {

        for (ny = -gcelly; ny <= gcelly; ny++)
        {

            for (nz = -gcellz; nz <= gcellz; nz++)
            {
                x = nx * rcpVector(0,0) + ny * rcpVector(0,1) + nz * rcpVector(0,2);
                y = nx * rcpVector(1,0) + ny * rcpVector(1,1) + nz * rcpVector(1,2);
                z = nx * rcpVector(2,0) + ny * rcpVector(2,1) + nz * rcpVector(2,2);

                gsq = x * x + y * y + z * z;

                if (gsq < rcpcutsq && gsq > 1.0e-8)
                {
                    maxGVec++;
                }
            }
        }
    }

    maxGVec += 100; // add a bit of padding
}

/* calculates g-vectors within a given radius.
 * returns the number of g vectors.
 * nb update every time volume called */

void Ewald::findGvector(double volume, const Eigen::MatrixXd & rcpVector)
{
    int
        nx,
        ny,
        nz;

    double
        gw,
        gzz,
        etasq = etaMad * etaMad,
        gsqfct = 1.0 / (4.0 * etasq),
        gfct0 = CTOEV * PI / (etasq * volume),
        gfct1 = gsqfct * gfct0;;

    double
        x,
        y,
        z,
        xx,
        yy,
        zz,
        gsq = 0.0,
        rcpcutsq = rcpSpaceCut * rcpSpaceCut;

    numGVec = 0; // zero number of vectors

    for (nx = 0; nx <= gcellx; nx++)
    {

        for (ny = -gcelly; ny <= gcelly; ny++)
        {

            for (nz = -gcellz; nz <= gcellz; nz++)
            {

                x = nx * rcpVector(0,0) + ny * rcpVector(0,1) + nz * rcpVector(0,2);
                y = nx * rcpVector(1,0) + ny * rcpVector(1,1) + nz * rcpVector(1,2);
                z = nx * rcpVector(2,0) + ny * rcpVector(2,1) + nz * rcpVector(2,2);

                gsq = x * x + y * y + z * z;

                if (gsq < rcpcutsq && gsq > 1.0e-8)
                {

                    //add the new vector and increase counter
                    gx[numGVec] = x;
                    gy[numGVec] = y;
                    gz[numGVec] = z;

                    numGVec++;

                }
            }
        }
    }
    

    
    for (int i = 0; i < numGVec; i++)
    {

        gx[i] *= TWOPI;
        gy[i] *= TWOPI;
        gz[i] *= TWOPI;

        xx = gx[i];
        yy = gy[i];
        zz = gz[i];
        gsq = xx * xx + yy * yy + zz * zz;
        x = gsq * gsqfct;

        gzz = 2.0 * exp(-x) / x;
        gzz = exp(-x) / x;

        if (xx < 1.0e-8)
        {
            g0[i] = gzz * gfct0;
        }
        else
        {
            g0[i] = gzz * gfct0 * 2.0;
        }
        gw = -2.0 * gzz * (1.0 + 1.0 / x);
        g1[i] = gw * gfct1;

    }

}

#ifdef _OPENMP

double Ewald::recipSpaceForce(const Eigen::MatrixXd &pos, Eigen::MatrixXd& force, Eigen::VectorXd stress, Eigen::VectorXd  atomcharge, int numatoms)
{
    double
        xx, yy, zz,
        ggx, ggy, ggz, gg0,
        cossum,
        sinsum,
        phase,
        phasefact,
        energyrcp = 0.0;


    double* sinphase; //[numatoms];  // cannot dynamically allocate for OMP
    double* cosphase; //[numatoms];  // you may be able to do this within the parallel region
    double* tmpForceX; //[numatoms];
    double* tmpForceY; //[numatoms];
    double* tmpForceZ; //[numatoms];
   
    #pragma omp parallel default(none)  \
        reduction(+: energyrcp) \
        shared(gx, gy, gz, g0, pos, force, atomcharge,  numGVec, numatoms ) \
        private(cossum, sinsum, phase, cosphase, sinphase, phasefact, ggx, ggy, ggz, gg0, xx, yy, zz, tmpForceX, tmpForceY, tmpForceZ)
    { 
        sinphase = (double*) calloc(numatoms, sizeof(double));
        cosphase = (double*) calloc(numatoms, sizeof(double));
        tmpForceX = (double*) calloc(numatoms, sizeof(double));
        tmpForceY = (double*) calloc(numatoms, sizeof(double));
        tmpForceZ = (double*) calloc(numatoms, sizeof(double));

        #pragma omp for schedule(static,4)
        for (int rcp = 0; rcp < numGVec; rcp++)
        {
            cossum = 0.0;
            sinsum = 0.0;

            ggx = gx[rcp];
            ggy = gy[rcp];
            ggz = gz[rcp];
            gg0 = g0[rcp];

            for (int atom = 0; atom < numatoms; atom++)
            {

                xx = pos(atom,0);
                yy = pos(atom,1);
                zz = pos(atom,2);

                phase = ggx * xx + ggy * yy + ggz * zz;

                cosphase[atom] = atomcharge[atom] * cos(phase);
                sinphase[atom] = atomcharge[atom] * sin(phase);

                cossum = cossum + cosphase[atom];
                sinsum = sinsum + sinphase[atom];

            }

            for (int atom = 0; atom < numatoms; atom++)
            {
                energyrcp += gg0 * (cosphase[atom] * cossum + sinphase[atom] * sinsum);

                phasefact = -gg0 * (-sinphase[atom] * cossum + cosphase[atom] * sinsum);

                tmpForceX[atom] += phasefact * ggx;
                tmpForceY[atom] += phasefact * ggy;
                tmpForceZ[atom] += phasefact * ggz;

            }

        } // end of loop over reciprocal lattice vectors
        for (int atom = 0; atom < numatoms; atom++)
        {
            #pragma omp atomic
            force(atom,0) += tmpForceX[atom];
        }
        for (int atom = 0; atom < numatoms; atom++)
        {
            #pragma omp atomic
            force(atom,1) += tmpForceY[atom];
        }
        for (int atom = 0; atom < numatoms; atom++)
        {
            #pragma omp atomic
            force(atom,2) += tmpForceZ[atom];
        }

        free(tmpForceX);
        free(tmpForceY);
        free(tmpForceZ);
        free(sinphase);
        free(cosphase);

    }
    // energy, forces and stresses must be converted to eV
    energyrcp *= 0.5;

    return energyrcp;
}
#else

double Ewald::recipSpaceForce(const Eigen::MatrixXd &pos, Eigen::MatrixXd& force, Eigen::VectorXd stress, Eigen::VectorXd  atomcharge, int numatoms)
{
    double
        xx, yy, zz,
        vx, vy, vz, v0, v1,
        //dummy,
        cossum,
        sinsum,
        phase,
        phasefact,
        phSumSq,
        phaseTotal,
        factor,
        tmp_energy = 0.0,
        energyrcp = 0.0;

    Eigen::VectorXd sinphase(numatoms);
    Eigen::VectorXd cosphase(numatoms);

    for (int rcp = 0; rcp < numGVec; rcp++)
    {
        cossum = 0.0;
        sinsum = 0.0;

        vx = gx[rcp];
        vy = gy[rcp];
        vz = gz[rcp];
        v0 = g0[rcp];
        v1 = g1[rcp];

        for (int atom = 0; atom < numatoms; atom++)
        {

            xx = pos(atom,0);
            yy = pos(atom,1);
            zz = pos(atom,2);

            phase = vx * xx + vy * yy + vz * zz;

            cosphase[atom] = atomcharge[atom] * cos(phase);
            sinphase[atom] = atomcharge[atom] * sin(phase);

            cossum = cossum + cosphase[atom];
            sinsum = sinsum + sinphase[atom];

        }

        phSumSq = 0.5 * (cossum * cossum + sinsum * sinsum);

        if(phSumSq > 1.0e-16)
        {
            phaseTotal += phSumSq * v0;

            factor = 2.0 * phSumSq * v1;

            // add in stress
            stress[0] += factor * vx * vx;
            stress[1] += factor * vy * vy;
            stress[2] += factor * vz * vz;
            stress[3] += factor * vy * vz;
            stress[4] += factor * vx * vz;
            stress[5] += factor * vx * vy;

        }


        //#pragma clang loop vectorize(enable)
        //#pragma omp simd aligned(forceX:64)
        //#pragma vector aligned //all vectors aligned

        for (int atom = 0; atom < numatoms; atom++)
        {
            tmp_energy += v0 * (cosphase[atom] * cossum + sinphase[atom] * sinsum);

            phasefact = -v0 * (-sinphase[atom] * cossum + cosphase[atom] * sinsum);

            force(atom,0) += phasefact * vx;
            force(atom,1) += phasefact * vy;
            force(atom,2) += phasefact * vz;

        }

    } // end of loop over reciprocal lattice vectors

    stress[0] -= phaseTotal;
    stress[1] -= phaseTotal;
    stress[2] -= phaseTotal;

    tmp_energy *= 0.5;
    energyrcp = tmp_energy;


    return energyrcp;
}
#endif
