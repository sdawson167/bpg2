#ifndef _PHASE_GUARD
#define _PHASE_GUARD

#include <stdexcept>

//Represents a phase
struct Phase
{
private:
    int m_index;

    Phase(int index) : m_index{ index } { }

public:
    //Creates a disordered phase
    Phase() : m_index{ 0 } { } //Default phase is the disordered one

    static const Phase dis;         // disordered
    static const int disId = 0;
    static const Phase lam;         // lamellar
    static const int lamId = 1;
    static const Phase gyr;         // gyroid
    static const int gyrId = 2;
    static const Phase hex;         // cylindrical
    static const int hexId = 3;
    static const Phase bcc;         // BCC
    static const int bccId = 4;
    static const Phase fcc;         // FCC
    static const int fccId = 5;
    static const Phase a15;         // A15
    static const int a15Id = 6;
    static const Phase sig;         // sigma
    static const int sigId = 7;
    static const Phase c14;         // C14
    static const int c14Id = 8;
    static const Phase c15;         // C15
    static const int c15Id = 9;


    //Converts this phase to its integer representation
    int toInt() const
    {
      return m_index;
    }

    bool operator ==(const Phase &rhs)
    {
      return rhs.m_index == m_index;
    }

    bool operator !=(const Phase &rhs)
    {
      return rhs.m_index != m_index;
    }

    //Generates a phase from the given integer representation
    static Phase fromInt(int id)
    {
      if (id > 10)
        throw std::invalid_argument("'id' is not within the correct range of values");

      return Phase(id);
    }
};
#endif
