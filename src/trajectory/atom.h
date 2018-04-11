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
        /// Constant access to atom position
        const array<double, DIM>& getPosition() const;
	/// Set the atomic position.
        void setPosition(double a_x, double a_y, double a_z);
    private:
        array<double, DIM> m_position;
};
#endif
