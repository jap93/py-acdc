===============
1. Introduction
===============

KMC exploits stochastic algorithms to explore rare events and coarse
grain the time evolution of the model
system. [1]_\ :sup:`,`\  [2]_\ :sup:`,`\  [3]_ In our calculations the
transition states are considered as thermally activated diffusional hops
and are governed by an Arrhenius equation

:math:`r_{D} = D_{0}exp( - \frac{Q}{\text{kT}})` [1]

where *r\ D* is the rate of an event, *D\ 0* is an exponential
pre-factor, *Q* is the activation energy of the hop, *k* is Boltzmann's
constant and *T* the temperature. We assume that the exponential
pre-factor D\ :sub:`0` is 1.0e\ :sup:`13` s\ :sup:`-1` which is the
vibrational frequency. In our simulations an event is either a cation or
an anion diffusional hop and the state of the model system is evolved by
choosing one event stochastically, according to the rate of the events
using the following equation.

:math:`\sum_{i = 1}^{m - 1}{r_{i} \leq \rho_{1}}\sum_{j = 1}^{N}r_{j} < \sum_{k = 1}^{m}r_{k}`
[2]

*where m* is the index of the chosen event and *N* is the total number
of possible events. Summation indices *i*, *j*, and *k* denote the
individual events, thus *r\ i*\ is the rate of the event *i*. *ρ\ i* is
a random number evenly distributed over the range [0, 1). This ensures
that faster events have a greater probability of being chosen than
slower events. The aKMC 500 events within each cycle allowing. Once an
event is chosen, the surface is modified to enact the diffusion event.
The simulation time is then advanced by

:math:`\mathrm{\Delta}t = \  - \frac{\ln\rho_{2}}{\sum_{i = 1}^{N}r_{i}}`
[3]

In equation 3, *t* is the elapsed time and *ρ\ 2*\ is a random number
evenly distributed over the range [0, 1).

The list of events can be pre-determined or calculated on-the-fly. A
pre-determined list requires less computation, but requires a prior
knowledge of the evolution of the simulation cell and/or material
structure and *all* the processes in operation over the entire
simulation. ACDC uses the alternative, yet computationally more
demanding, approach to calculate the activation energies on-the-fly. and
the methodology for calculating the transition states is described in
the next section.

The original C++ version has been re-written in python (the computationally expensive bits are still in C++) so that
a wider variety of forcefields can be employed, including MLIP's. In theory any method that has a ASE Calculator can
be used within the program (in reality only a few have been tested).

===============================
2. Calculation of saddle points
===============================

It is possible to do this with techniques such as molecular
dynamics, [4]_ the dimer method  [5]_ and activation relaxation
technique (ART). [6]_ At the moment the Dimer method, via the Atomic Simulation Environment (ASE), and ART method is the sole techniques
coded into ACDC. The ART method coded in closely
follows the recipe described in reference [7]_ and briefly describe
the computational method below. ART and Dimer are members of the minimum mode
following methods that employ an eigenvector-following approach to
efficiently determine transition states. 

In our ART calculations we used the FIRE and FIRE2 method [8]_ to converge the atomic
positions to the transition state. The Dimer method uses the ASE toolkits

===============
3. Installation
===============

External libraries required are Eigen3, pybind11 and ASE. In addition libraries to satisfy the ASE calculator or MLIP potentials may be necessary.

The most convenient way to run the program is using empirical potentials. The code can be built in the following way::
cd py-acdc/potentials/rigidion(or metal). There is a cmake file to build a library for either 2-body rigid ion model or many-body metal potentials.
cmake .  
cmake --build . --clean-first
cp library.so ../../python/.
cd ../../python
You will also need to copy the appropriate interface to the specific file field.py. Some have already been included (e.g. field_rigid_ion.py, field_metal.py).

If you are intending to use a MLIP and ASE then you do not need to undertake the above, just copy the appropriate field file (e.g. field_janus.py for MACE-mp).

==================
4. Parallelisation
==================

The evaluation of energies and forces are evaluated using OpenMP whilst each transition point search is task farmed using MPI. NB this has not been tested for multiple GPU's.

=========================
5. How to run the program
=========================


The program requires three input files, *control*, *potentials* and *basis* that contain the
keywords for the functionality of the program, empirical potential parameters and the atomic
positions respectively.

-----------
5.1 control
-----------

The input is broken into sections depending on the functionality. The
main input key words are

+--------------+-------------+------------------------------------------------------+
| **Key word** | **Type**    | **Functionality**                                    |
+==============+=============+======================================================+
| seed         | int         | Seed for the random number. If no number is          |
|              |             | specified the system clock is used. The workgroup    |
|              |             | number is added automatically to seed so that each   |
|              |             | workgroup will follow a different trajectory.        |
+--------------+-------------+------------------------------------------------------+
| restart      |             | Restart the calculation with new KMC iteration       |
+--------------+-------------+------------------------------------------------------+
| freeze       | int         | Freeze the following types of atoms. E.g.            |
|              |             |                                                      |
|              |             | freeze 1                                             |
|              |             | Au                                                   |
|              |             |                                                      |
|              |             | This prevents any movement of Au type atoms          |
+--------------+-------------+------------------------------------------------------+
| close        |             | Finish all input                                     |
+--------------+-------------+------------------------------------------------------+
| kmc{         |             | Starts the processing of the KMC specific key words  |
+--------------+-------------+------------------------------------------------------+
| art{         |             | Starts the processing of the ART specific key words  |
+--------------+-------------+------------------------------------------------------+
| relax{       |             | Starts the processing of the energy minimisation     |
|              |             | specific key words                                   |
+--------------+-------------+------------------------------------------------------+

Description of keywords to control the kinetic Monte Carlo functionality. This input
is started with kmc{ and ended with }. NB the enrgies are in eV and distance A.

+----------------+-------------+------------------------------------------------------+
| **Key word**   | **Type**    | **Functionality**                                    |
+================+=============+======================================================+
| kmcsteps       | int         | The number of kmc steps/cycles                       |
+----------------+-------------+------------------------------------------------------+
| mincap         | double      | The minimum activation energy.                       |
+----------------+-------------+------------------------------------------------------+
| prefactor      | double      | The prefactor used to calculate rates.               |
+----------------+-------------+------------------------------------------------------+
| window         | double      | The window for accepting kmc energies.               |
+----------------+-------------+------------------------------------------------------+
| kmctemperature | double      | The temperature to be used by the kmc simulation     |
+----------------+-------------+------------------------------------------------------+
| kmcevents      | int         | The number of events collected per cycle.            |
+----------------+-------------+------------------------------------------------------+
| basin2delta    | double      | The magnitude of the displacement away from the      |
|                |             | saddle point prior to the attempted relaxation to    |
|                |             | the second basin.                                    |
+----------------+-------------+------------------------------------------------------+
| kmcmethod      | string      | The method used for the saddle search. Only ART or   |
|                |             | dimer are available at the moment                    |
+----------------+-------------+------------------------------------------------------+
| kmcbasinradius | double      | The minimum displacement of an atom before it is     |
|                |             | considered to have entered a new basin               |
+----------------+-------------+------------------------------------------------------+
| recycle        |             | recycles saddle points. Read the code as this is     |
|                |             | experimental and may have significant impact on the  |
|                |             | results                                              |
+----------------+-------------+------------------------------------------------------+
| usegaussian    | double      | The atoms are given a random displacement weighted   |
|                |             | by a Gaussian of the specified width                 |
+----------------+-------------+------------------------------------------------------+

Description of keywords to control either the Activation Relaxation Technique. Again units are specified
in the main directives above.

+---------------------+-------------+------------------------------------------------------+
| **Key word**        | **Type**    | **Functionality**                                    |
+=====================+=============+======================================================+
| numvectors          | int         | The number of Lnanczos vectors used to obtain        |
|                     |             | eigenvalues. Default 20                              |
+---------------------+-------------+------------------------------------------------------+
| maxeigenvalue       | double      | The max eigenvalue. Once an eigenvalue falls below   |
|                     |             | this value the forces parallel to the eigenvalue are |
|                     |             | used.                                                |
+---------------------+-------------+------------------------------------------------------+
| initialdisplacement | double      | The displacement used to activate the ions at the    |
|                     |             | start of the search.                                 |
+---------------------+-------------+------------------------------------------------------+
| eigentolerence      | double      | The tolerence to converge the eigenvalues.           |
+---------------------+-------------+------------------------------------------------------+
| lanczosdisplacement | double      | The displacement of atoms used to calculate the      |
|                     |             | eigenvalues from the tri-diagonal matrix             |
+---------------------+-------------+------------------------------------------------------+
| maxstep             | int         | The number of iterations to calculate the saddle     |
|                     |             | point.                                               |
+---------------------+-------------+------------------------------------------------------+
| minmethod           | string      | The mminimisation technique to find the saddle point |
|                     |             | Only FIRE(2) is available at the moment.             |
+---------------------+-------------+------------------------------------------------------+
| timestep            | double      | The timestep used by FIRE. Typically should be       | 
|                     |             | similar to that used by MD.                          |
+---------------------+-------------+------------------------------------------------------+
| alpha               | double      | The value of alpha used in FIRE                      |
+---------------------+-------------+------------------------------------------------------+
| damp                | double      | damping factor for the parallel forces               |
+---------------------+-------------+------------------------------------------------------+
| forcetol            | double      | The convergence criterion for the minimisation       |
+---------------------+-------------+------------------------------------------------------+

Input keywords for the relaxation

+---------------------+-------------+------------------------------------------------------+
| **Key word**        | **Type**    | **Functionality**                                    |
+=====================+=============+======================================================+
| debug               |             | increases the amount of information output to files. |
+---------------------+-------------+------------------------------------------------------+
| relaxsteps          | int         | The number of iterations of the minimisaer.          |
+---------------------+-------------+------------------------------------------------------+
| initialdisplacement | double      | The displacement used to activate the ions at the    |
|                     |             | start of the search.                                 |
+---------------------+-------------+------------------------------------------------------+
| maxstep             | double      | The maximum size of the displacement in FIRE.        |
+---------------------+-------------+------------------------------------------------------+
| method              | string      | The mminimisation technique. There is a choice       |
|                     |             | between FIRE, FIRE2 and ASE the moment.              |
|                     |             | (The latter uses the minimisation technique from the |
|                     |             | librray program ASE.)                        |
+---------------------+-------------+------------------------------------------------------+
| timestep            | double      | The timestep used by FIRE. Typically should be       | 
|                     |             | similar to that used by MD.                          |
+---------------------+-------------+------------------------------------------------------+
| alpha               | double      | The value of alpha used in FIRE                      |
+---------------------+-------------+------------------------------------------------------+
| forcetol            | double      | The convergence criterion for the minimisation       |
+---------------------+-------------+------------------------------------------------------+

--------------
5.1 potentials
--------------

+---------------------+-------------+------------------------------------------------------+
| **Key word**        | **Type**    | **Functionality**                                    |
+=====================+=============+======================================================+
| cutoff              | double      | The short range cutoff for the potential and Ewald . |
+---------------------+-------------+------------------------------------------------------+
| noimage             |             | The nearest image conventin is used by default. This |
|                     |             | keyword uses a slower multiple image method that     |
|                     |             | is better suited to small simulation cells.          |
+---------------------+-------------+------------------------------------------------------+
| noewald             |             | The Ewald sum is used by default for the two-body    |
|                     |             | potential model. Thus this keyword switches it off.  |
+---------------------+-------------+------------------------------------------------------+
| ewald precision     |             | Controls the accuracy of the Ewald sum. Default      |
|                     |             | value: 1.0e-6                                        |
+---------------------+-------------+------------------------------------------------------+
| species             | int         | species keyword followed by the number of different  |
|                     |             | speccies. Each element type should be input as       |
|                     |             | follows:                                             |
|                     |             | name  mass charge atomic_number                      |
+---------------------+-------------+------------------------------------------------------+

As described in the installation section the program can be compiled with either two-body (including Ewald sum),
many-body (metal potentials) or ASE calculator. Note all parameters are in electron volts! 
Parameters compatible with the two-body are:

+---------------------+-------------+------------------------------------------------------+
| **Key word**        | **Parameters**                                                     |
+=====================+=============+======================================================+
| buck                | *A* , :math:`{\alpha}` , *C*                                       |
+---------------------+--------------------------------------------------------------------+
| morse               | *D* , r\ :sub:`eq` , *k*                                           |
+---------------------+--------------------------------------------------------------------+
| ljones              | :math:`{\eta}` , :math:`{\sigma}`                                  |
+---------------------+--------------------------------------------------------------------+
| bhm                 | *A* , :math:`{\alpha}` , *C* , *D*                                 |
+---------------------+--------------------------------------------------------------------+

Here is an example::

   cutoff 8.0
   noimage
   species 2
   Mg 24.0 2.0 12
   O 16.0 -2.0  8
   twobody 2
   buck
   Mg O  1428.5 0.2945  0.00
   buck
   O  O 22764.3 0.1490 27.879
   close

Parameters compatible with the metal potentials are:

+---------------------+-------------+------------------------------------------------------+
| **Key word**        | **Parameters**                                                     |
+=====================+=============+======================================================+
| stch                | :math:`{\eta}` , *a* , *n* , *m* , *c*                             |
+---------------------+--------------------------------------------------------------------+
| gupta               | *A* , r\ :sub:`eq` , *p*, *B*, *q*                                 |
+---------------------+--------------------------------------------------------------------+
| fnsc                | *c0* , *c1* , *c2* , *c* , *A* , *d* , :math:`{\Beta}`             |
+---------------------+--------------------------------------------------------------------+

Here is an example::

   cutoff 6.5
   species 1
   Al  25.0  0.0 13
   manybody 1 ev
   suttonchen
   Al  Al   0.033147    4.05       7.0        6.0         16.399
   close

The is also the possibility of using machine learned interatomic potentials with ASE calculators and the ASE dimer method. A potentials file is still needed
to setup the calculation:

   species 2
   B 16.0 0.0
   C 16.0 0.0
   model  CoB_v3.model
   close

By default thrdr will run on a single GPU (multiple GPU's has not been tested).

In all calculations a file search.env is required and is employed to switch on or off the ASE functionality. For rigid ion or metal calculations this should be:

   USE_ASE_RELAX = "False"
   USE_ASE_SEARCH = "False"

whilst for MLIP's and ASE dimer the False values should be changed to True:

   USE_ASE_RELAX = "True"
   USE_ASE_SEARCH = "True"

---------
5.2 basis
---------

The basis format follows that of a simplified extended xyz. The minimum format is::

   number_of_atoms
   Lattice="vectors * 9"
   name      x  y  z
   name      x  y  z

For example::

   2549
   Lattice="32.7232154959 0.0 0.0 0.0 32.7232154959 0.0 0.0 0.0 32.7232154959"
   O 1.4168603868  1.3899782202  1.3108695097
   O 4.1014301991  4.0478192159  1.4171530724
   O 3.7163970637  -1.3764087458  4.0580241564

=============
6. References
=============

.. [1]
   D.T. Gillespie, *J. Phys. Chem.*, 1997, **81**, 2340-2361.

.. [2]
   A.F. Voter, *Phys. Rev. B*, 1986, **34**, 6819-6829.

.. [3]
   C.C. Battaile, *Comput. Methods Appl. Mech. Engrg.*, 2008, **197**,
   3386-3398.

.. [4]
   D. Frenkel and B. Smit, *Understanding Molecular Simulation: From
   Algorithms to Applications*, 2002, Academic Press.

.. [5]
   G. Henkelman and H. Jónsson, *J. Chem. Phys.*, 1999, **111**,
   7010-7020.

.. [6]
   G.T. Barkema and N. Mousseau, *Computational Materials Science*,
   2001, **20**, 285-292.

.. [7]
   R.A. Olsen, G.J. Kroes, and G. Henkelman, *The Journal of Chemical Physics*, 2004,  **121(20)**, 9776–9792.

.. [8]
   E. Bitzek, P. Koskinen, F. Gähler, M. Moseler, P. Gumbsch, *Phys. Rev. Lett.*, 2006, **97** , Article 170201.
