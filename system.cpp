#include "system.h"
#include "velocityverlet.h"
#include "lennardjones.h"
#include "statisticssampler.h"
#include "unitconverter.h"
#include "math/random.h"

System::System()
{

}

System::~System()
{
    for(Atom *atom : m_atoms) {
        delete atom;
    }
    m_atoms.clear();
}

void System::applyPeriodicBoundaryConditions() {

    int nrOfAtoms = m_atoms.size();

    for (int i = 0; i < nrOfAtoms; i++){
        Atom*a = m_atoms[i];

        // 10 is latticesize, needs to be changed.
        // Positions should not larger than 5 or less than -5
        if (a->position(0) < 0 ){
            a->position(0) += systemSize().x();
        }
        else if (a->position(0) >= systemSize().x()){
            a->position(0) -= systemSize().x();
        }

        if (a->position(1) < 0){
            a->position(1) += systemSize().y();
        }
        else if (a->position(1) >= systemSize().y()){
            a->position(1) -= systemSize().y();
        }

        if (a->position(2) < 0){
            a->position(2) += systemSize().z();
        }
        else if (a->position(2) >= systemSize().z()){
            a->position(2) -= systemSize().z();
        }

    }

    // Read here: http://en.wikipedia.org/wiki/Periodic_boundary_conditions#Practical_implementation:_continuity_and_the_minimum_image_convention
}

//Function removeTotalMomentum:
void System::removeTotalMomentum() {

    // Find the total momentum and remove momentum equally on each atom so the total momentum becomes zero.
    totalMomentum.zeros();
    averageVelocity.zeros();
    totalVelocity.zeros();

    numberOfAtoms = m_atoms.size();

    for (int i=0; i < numberOfAtoms; i++){ //Loop through all atoms
        Atom * atom = m_atoms[i];
        totalVelocity += atom->velocity; //Adding up the velocities of the atoms.
        //std::cout << atom->velocity << std::endl;
    }

    //std::cout << totalMomentum << std::endl;
    averageVelocity = totalVelocity/numberOfAtoms; //Find average velocity of the atoms

    for (int i=0; i<numberOfAtoms; i++){
        Atom * atom = m_atoms[i];
        atom->velocity -= averageVelocity; //Subtract the average velocity from all atoms to get total momentum=0.
        totalMomentum += atom->mass() * atom->velocity;//Adding up the momentum of all atoms to get the total
    }
    std::cout << "Totalt moment:" << totalMomentum << std::endl;
} //End func removeTotalMomentum

//Function FCCLattice creates the lattice consisting of N unit cells in each dimension, and places out atoms with initial pos and vel.
void System::createFCCLattice(int numberOfUnitCellsEachDimension, double latticeConstant, double temperature) {//A lattice is built up by unit cells
    // You should implement this function properly. Right now, 100 atoms are created uniformly placed in the system of size (10, 10, 10).
    b = latticeConstant;
    N = numberOfUnitCellsEachDimension;
    setSystemSize(vec3(b*N, b*N, b*N)); //Evt slik
    //std::cout << m_atoms.size() << std::endl;

    /*
    for(int i=0; i<100; i++) {
        Atom *atom = new Atom(UnitConverter::massFromSI(6.63352088e-26));
        double x = Random::nextDouble(0, 10); // random number in the interval [0,10]
        double y = Random::nextDouble(0, 10);
        double z = Random::nextDouble(0, 10);
        atom->position.set(x,y,z);
        atom->resetVelocityMaxwellian(temperature);
        m_atoms.push_back(atom);
    }
    setSystemSize(vec3(10, 10, 10)); // Remember to set the correct system size!
    */

    //Oppg c) :
    for(int i=0; i<N; i++){ //x-direction
        for(int j=0; j<N; j++){ //y-direction
            for(int k=0; k<N; k++){ //z-dirextion
                // R_{i,j,k} = (i*b, j*b, k*b) is the origin of unit cell (i,j,k)
                Atom *atom1 = new Atom(UnitConverter::massFromSI(6.63352088e-26)); //Creates four new atoms per cell
                Atom *atom2 = new Atom(UnitConverter::massFromSI(6.63352088e-26));
                Atom *atom3 = new Atom(UnitConverter::massFromSI(6.63352088e-26));
                Atom *atom4 = new Atom(UnitConverter::massFromSI(6.63352088e-26));

                //Set initial positions of all atoms
                atom1->position.set( i*b, j*b, k*b);
                atom2->position.set( (0.5+i)*b, (0.5+j)*b, k*b );
                atom3->position.set( i*b, b*(0.5+j), b*(0.5+k) );
                atom4->position.set( b*(0.5+i), j*b, b*(0.5+k) );

                //Set initial velocity of all atoms. Given by a Maxwellian velocity distribution
                atom1->resetVelocityMaxwellian(temperature);
                atom2->resetVelocityMaxwellian(temperature);
                atom3->resetVelocityMaxwellian(temperature);
                atom4->resetVelocityMaxwellian(temperature);

                //'append' atoms to the 'm_atoms' array
                m_atoms.push_back(atom1);
                m_atoms.push_back(atom2);
                m_atoms.push_back(atom3);
                m_atoms.push_back(atom4);

            } //End k-loop
        } //End j-loop
    } //End i-loop

    //Set systemsize:
    //setSystemSize(vec3(systemSize().x(), systemSize().y(), systemSize().z())); // Remember to set the correct system size!
}//End create FCCLattice func

void System::calculateForces() {
    for(Atom *atom : m_atoms) {
        atom->resetForce(); //Resets force to zero for all atoms
    }
    m_potential.calculateForces(*this); // this is a pointer, *this is a reference to this object
}//End calculateForces func

//step function:
void System::step(double dt) {
    m_integrator.integrate(*this, dt); //Integrates with velocityVerlet for each timestep
    m_steps++;
    m_time += dt;
}//End step function

