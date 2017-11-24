#include "math/random.h"
#include "lennardjones.h"
#include "velocityverlet.h"
#include "system.h"
#include "statisticssampler.h"
#include "atom.h"
#include "io.h"
#include "unitconverter.h"
#include <iostream>
#include <iomanip>

using namespace std;

//Skriv:
//Atom a; //Lager objektet a av klassen Atom
//a. //a. gir deg naa alle funks og variabler inneholdt i Atom klassen.

int main(int numberOfArguments, char **argumentList)
{
    int numberOfUnitCells = 5;  //Want to make a cube out of several unit cells
    double initialTemperature = UnitConverter::temperatureFromSI(300.0); // measured in Kelvin
    double latticeConstant = UnitConverter::lengthFromAngstroms(5.26); // measured in angstroms

    // If a first argument is provided, it is the number of unit cells
    if(numberOfArguments > 1) numberOfUnitCells = atoi(argumentList[1]);
    // If a second argument is provided, it is the initial temperature (measured in kelvin)
    if(numberOfArguments > 2) initialTemperature = UnitConverter::temperatureFromSI(atof(argumentList[2]));
    // If a third argument is provided, it is the lattice constant determining the density (measured in angstroms)
    if(numberOfArguments > 3) latticeConstant = UnitConverter::lengthFromAngstroms(atof(argumentList[3]));

    double dt = UnitConverter::timeFromSI(1e-15); // Measured in seconds.

    cout << "One unit of length is " << UnitConverter::lengthToSI(1.0) << " meters" << endl;
    cout << "One unit of velocity is " << UnitConverter::velocityToSI(1.0) << " meters/second" << endl;
    cout << "One unit of time is " << UnitConverter::timeToSI(1.0) << " seconds" << endl;
    cout << "One unit of mass is " << UnitConverter::massToSI(1.0) << " kg" << endl;
    cout << "One unit of temperature is " << UnitConverter::temperatureToSI(1.0) << " K" << endl;

    System system; //Make an object "system" of the System class

    //Func that creates lattice and places the atoms in random positions in it:
    system.createFCCLattice(numberOfUnitCells, latticeConstant, initialTemperature); //latticeconstant=b is the size of the unit cells

    system.potential().setEpsilon(1.0); //Sets episilon equal to 1.0?
    system.potential().setSigma(1.0); //Sets sigma equal to 1.0?

    //Func that make total momentum of system equal to zero
    system.removeTotalMomentum();
    //std::cout << "Totalt moment:" << totalMomentum << std::endl;

    StatisticsSampler statisticsSampler;
    IO movie("movie.xyz"); // To write the state to file

    cout << setw(20) << "Timestep" <<
            setw(20) << "Time" <<
            setw(20) << "Temperature" <<
            setw(20) << "KineticEnergy" <<
            setw(20) << "PotentialEnergy" <<
            setw(20) << "TotalEnergy" << endl;
    for(int timestep=0; timestep<1000; timestep++) {
        system.step(dt); //The step function integrates with velocity verlet?
        statisticsSampler.sample(system);
        //Prints out parameters for every 100th timestep
        if( timestep % 100 == 0 ) {
            // Print the timestep every 100 timesteps
            cout << setw(20) << system.steps() <<
                    setw(20) << system.time() <<
                    setw(20) << statisticsSampler.temperature() <<
                    setw(20) << statisticsSampler.kineticEnergy() <<
                    setw(20) << statisticsSampler.potentialEnergy() <<
                    setw(20) << statisticsSampler.totalEnergy() << endl;
        }
        movie.saveState(system);
    }

    movie.close();

    return 0;
}
