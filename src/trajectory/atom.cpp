#include "atom.h"
#include <array>
#include <string.h>

Atom::Atom() {};

Atom::~Atom() {};

Atom::Atom(array<double, DIM> a_coords)
{
  for (int i = 0; i < DIM; i++)
  {
    m_position[i] = a_coords[i];
  }
  m_charge = 0;
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

void Atom::setCharge(double a_q)
{
  m_charge = a_q;
}

const double Atom::getCharge() const
{
  return m_charge;
}

/// Set name of atom (a string)
void Atom::setName(const char* a_name)
{
  strcpy(m_name, a_name);
}
/// Set type of atom (an int)
void Atom::setType(int a_type)
{
  m_type = a_type;
}
/// Constant access to name.
const char* Atom::getName() const
{
  return m_name;
}
/// Constant access to type.
const int Atom::getType() const
{
  return m_type;
}
