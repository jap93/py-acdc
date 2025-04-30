#include "MetalPotential.h"

std::vector<std::string> split(std::string s);

MetalPotential::MetalPotential()
{

}
MetalPotential::~MetalPotential()
{  
    
}


void MetalPotential::readEAM(std::string infile, std::ofstream& outStream)
{
    std::string
        dummy,
        keyWord,
        line;
    std::vector<std::string> words;
    std::ifstream inStream;

    inStream.open(infile, std::ios::in);
    if (inStream.fail())              // check to see file is there
    {
        outStream << "\n*** could not find EAM file" << std::endl;
        outStream.flush();
        outStream.close();
        exit(EXIT_FAILURE);
    }

    std::getline(inStream, line); // first lines are dummy

    std::getline(inStream, line);
	words = split(line);
    numComponents = std::stoi(words[0]);

    for (int nc = 0; nc < numComponents; nc++)
    {
                
        std::getline(inStream, line);
	    words = split(line);
        keyWord = words[0];
        transform(keyWord.begin(), keyWord.end(), keyWord.begin(), ::tolower);

        if (keyWord == "pair")
        {
            numPtsRep = std::stoi(words[3]);
            startRep = std::stoi(words[4]);
            finishRep = std::stoi(words[5]);
            stepRep = 1 + (finishRep - startRep) / numPtsRep;

            parRep.resize(numPtsRep);
            distRep.resize(numPtsRep);
            derivRep.resize(numPtsRep);

            readParameters(parRep, numPtsRep, inStream);

            for (int i = 0; i < numPtsRep; i++)
            {
                distRep[i] = startRep + i * stepRep;
            }
            calculateDerivativeMesh(parRep, derivRep, stepRep, numPtsRep);

            potCutoff = finishRep;
        }
        else if (keyWord == "dens")
        {
            numPtsDens = std::stoi(words[2]);
            startDens = std::stoi(words[3]);
            finishDens = std::stoi(words[4]);
            stepDens = (finishDens - startDens) / numPtsDens;

            parDens.resize(numPtsDens);
            distDens.resize(numPtsDens);
            derivDens.resize(numPtsDens);

            readParameters(parDens, numPtsDens, inStream);

            for (int i = 0; i < numPtsDens; i++)
            {
                distDens[i] = startDens + i * stepDens;
            }
            calculateDerivativeMesh(parDens, derivDens, stepDens, numPtsDens);

            if (potCutoff < finishDens)
                potCutoff = finishDens;
        }
        else if (keyWord == "embed")
        {
            numPtsEmbed = std::stoi(words[2]);
            startEmbed = std::stoi(words[3]);
            finishEmbed = std::stoi(words[4]);
            stepEmbed = (finishEmbed - startEmbed) / numPtsEmbed;

            parEmbed.resize(numPtsEmbed);
            distEmbed.resize(numPtsEmbed);
            derivEmbed.resize(numPtsEmbed);

            readParameters(parEmbed, numPtsEmbed, inStream);

            for (int i = 0; i < numPtsEmbed; i++)
            {
                distEmbed[i] = startEmbed + i * stepEmbed;
            }
            calculateDerivativeMesh(parEmbed, derivEmbed, stepEmbed, numPtsEmbed);
        }
        else
        {
            outStream << " wrong line in readEAM " << line << std::endl;
            outStream.flush();
            exit(EXIT_FAILURE);
        }
    }
    
}

void MetalPotential::readLEAM(std::string infile, std::ofstream& outStream)
{
    std::string
        dummy,
        keyWord,
        line;
    std::vector<std::string> words;
    std::ifstream inStream;

    inStream.open(infile, std::ios::in);
    if (inStream.fail())              // check to see file is there
    {
        outStream << "\n*** could not find LAMMPS style EAM file " << infile << std::endl;
        outStream.flush();
        outStream.close();
        exit(EXIT_FAILURE);
    }

    std::getline(inStream, line); // first 4 lines are dummy (3 appear to be for comments)
    std::getline(inStream, line); 
    std::getline(inStream, line); 
    std::getline(inStream, line); 

    std::getline(inStream, line);

	words = split(line);

    numPtsEmbed = std::stoi(words[0]);
    stepEmbed = std::stod(words[1]);
    numPtsRep = numPtsDens = std::stoi(words[2]);
    stepRep = stepDens = std::stod(words[3]);
    potCutoff = std::stod(words[4]);
    

    parEmbed.resize(numPtsEmbed);
    distEmbed.resize(numPtsEmbed);
    derivEmbed.resize(numPtsEmbed);

    parRep.resize(numPtsRep);
    distRep.resize(numPtsRep);
    derivRep.resize(numPtsRep);

    parDens.resize(numPtsDens);
    distDens.resize(numPtsDens);
    derivDens.resize(numPtsDens);

    std::getline(inStream, line); // dummy line

    //embedding or glue
    readParameters(parEmbed, numPtsEmbed, inStream);
    for (int i = 0; i < numPtsEmbed; i++)
    {
        //std::getline(inStream, line);
	    //words = split(line);
        //parEmbed[i] = std::stod(words[0]);
        distEmbed[i] = i * stepEmbed;
    }

    calculateDerivativeMesh(parEmbed, derivEmbed, stepEmbed, numPtsEmbed);
    //outStream << " finished embed " << line << std::endl;

    //density for calculating rho
    readParameters(parDens, numPtsDens, inStream);
    for (int i = 0; i < numPtsDens; i++)
    {
        //std::getline(inStream, line);
	    //words = split(line);
        //parDens[i] = std::stod(words[0]);
        distDens[i] = i * stepDens;
    }
    calculateDerivativeMesh(parDens, derivDens, stepDens, numPtsDens);
    //outStream << " finished dens " << line << std::endl;

    //pair potential
    readParameters(parRep, numPtsRep, inStream);
    for (int i = 0; i < numPtsRep; i++)
    {
        //std::getline(inStream, line);
	    //words = split(line);
        distRep[i] = i * stepRep;
        //parRep[i] = std::stod(words[0]) / distRep[i];
        parRep[i] = parRep[i] / distRep[i];
    }
    calculateDerivativeMesh(parRep, derivRep, stepRep, numPtsRep);
    //outStream << " finished rep " << line << std::endl;
}

void MetalPotential::readParameters(Eigen::VectorXd& par, int numPts, std::ifstream& inStream)
{
    int
        id = 0;

    std::string
        line;

    std::vector<std::string> words;

    while (id < numPts)
    {
        std::getline(inStream, line);
        words = split(line);

        for (unsigned int i = 0; i < words.size(); i++)
        {
            par[id] = std::stod(words[i]);
            id++;
        }

    }

}

void MetalPotential::calculateDerivativeMesh(Eigen::VectorXd& pot, Eigen::VectorXd& deriv, double dr, int npts)
{
    for(int i = 0; i < npts - 1;i++)
        deriv[i] = (pot[i+1] - pot[i]) / dr;
}

