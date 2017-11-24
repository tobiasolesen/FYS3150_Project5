#include "lennardjones.h"
#include "system.h"
#include "math.h"

double LennardJones::potentialEnergy() const
{
    return m_potentialEnergy;
}

double LennardJones::sigma() const
{
    return m_sigma;
}

void LennardJones::setSigma(double sigma)
{
    m_sigma = sigma;
}

double LennardJones::epsilon() const
{
    return m_epsilon;
}

void LennardJones::setEpsilon(double epsilon)
{
    m_epsilon = epsilon;
}

void LennardJones::calculateForces(System &system)
{
    //Constants:
    double sigma6 = pow( sigma(), 6 );
    double sigma12 = pow( sigma(), 12);
    double F_x;
    double F_y;
    double F_z;

    m_potentialEnergy = 0; // Remember to compute this in the loop
    //Loops over all atoms except i=j to calculate forces between two atoms i and j
    for(int i=0; i < system.atoms().size(); i++){
        for(int j=i+1; j< system.atoms().size(); j++){
            Atom * atom_i = system.atoms()[i];
            Atom * atom_j = system.atoms()[j];
            vec3 r_ij = atom_i->position - atom_j->position; //Vector from atom i to atom j
            double len_r_ij = r_ij.length(); //Length of vector r_ij
            F_x = -24*epsilon()*( sigma6/pow(len_r_ij, 7) - 2*sigma12/pow(len_r_ij, 13) ) * atom_i->position(0)/len_r_ij;
            F_y = -24*epsilon()*( sigma6/pow(len_r_ij, 7) - 2*sigma12/pow(len_r_ij, 13) ) * atom_i->position(1)/len_r_ij;
            F_z = -24*epsilon()*( sigma6/pow(len_r_ij, 7) - 2*sigma12/pow(len_r_ij, 13) ) * atom_i->position(2)/len_r_ij;
            //atom_i->force += (F_x, F_y, F_z);
            //atom_j->force += -(F_x, F_y, F_z);
            vec3 force_ij(F_x, F_y, F_z);
            vec3 force_ji(-F_x, -F_y, -F_z);

            atom_i->force += force_ij;
            atom_j->force += force_ji;

            //atom_i->force += ( -24*epsilon()*( (sigma6/pow(len_r_ij, 7)) - (2*sigma12/pow(len_r_ij, 13)) ) );




        }

    }

}//End function LennardJones
