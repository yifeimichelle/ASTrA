#ifndef _ATOM_H_
#define _ATOM_H_

#include <array>
using namespace std;

/// A class describing a single atom
/**
 * This class stores atom positions
 */

class Atom
{
  public:
    Atom();
    /// Constructor.
    Atom(array<double, DIM> a_coords);
    virtual ~Atom(); // destructor
    /// Constant access to atom position.
    const array<double, DIM>& getPosition() const;
    /// Set the atomic position.
    void setPosition(double a_x, double a_y, double a_z);
    /// Set the atomic charge.
    void setCharge(double a_q);
    /// Set name of atom (a string)
    void setName(const char* a_name);
    /// Set type of atom (an int)
    void setType(int a_type);
    /// Constant access to charge.
    const double getCharge() const;
    /// Constant access to name.
    const char* getName() const;
    /// Constant access to type.
    const int getType() const;
  private:
    array<double, DIM> m_position;
    double m_charge;
    char m_name[10];
    unsigned int m_type;
};
#endif
