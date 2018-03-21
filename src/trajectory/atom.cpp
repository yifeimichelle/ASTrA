#include "atom.h"
#include <array>

Atom::Atom() {};

Atom::~Atom() {};

Atom::Atom(array<double, DIM> a_coords)
{
    for (int i = 0; i < DIM; i++)
    {
        m_position[i] = a_coords[i];
    }
};

const array<double, DIM>& Atom::getPosition() const
{
    return m_position;
}

void Atom::setPosition(double a_x, double a_y, double a_z)
{
    m_position[0] = a_x;
    m_position[1] = a_y;
    m_position[2] = a_z;
}
