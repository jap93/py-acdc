#pragma once

// unfortunately this is not constant as it is alterred in JobControl.cpp
// to allow for different energy units
const double BOLTZMANN = .00008617333262145; // BOLTZMANN constant in eV

const int MAXGVEC = 40000; // max number of g vectors
const int MAXMESH = 5000; // number of points along mesh

const double PI = 3.14159265358979323846264338328;
const double PISQ = PI * PI;
const double TWOPI = (2.0 * PI);
const double HALFPI = PI / 2.0;
const double ROOTPI = 1.77245385090551602729816748334;

const double TORADIANS = (PI / 180.0);

//internal conversion units
const double CTODL = 138936.0689;
const double CTOEV = 14.3997584;
const double EVTODL = 9648.530821;
const double INTERNALTOEV = 1.0 / EVTODL;
//1 Ha/Bohr = 51.42208619083232 eV/Angstrom
const double ANGTOBOHR = 1.889726134583548707935;
const double BOHRTOANG = 1.0 / 1.889726134583548707935;
const double HARTREETOEV = 27.2113846081672;
const double EVTOHARTREE = 1.0 / 27.2113846081672;
const double AUTOFS = 0.0241888468; //conversion au to femtoseconds

//status constants
const int SUCCESS = 0;
const int FAILED = -1;
const int INTERRUPTED = -2;
const int TERMINATED = -3;
const int STARTSEARCH = -4;
const int CONTINUESEARCH = -5;

//memory management
const unsigned int DOFFSET = 64;

//data transfer
const int PARENTPOSSEND = 24;   // for sending position data from parent to the workgroups
const int GRPPOSSEND = 25;      // work group send of positions to parent
