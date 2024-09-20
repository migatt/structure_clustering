#ifndef __ATOM__
#define __ATOM__

#include <iostream>

class Position {
    double _x, _y, _z;

public:
    Position(const double x, const double y, const double z) : _x(x), _y(y), _z(z) {}
    double x() const { return _x; }
    double y() const { return _y; }
    double z() const { return _z; }
};

class Atom {
    int _atomicNumber;
    Position _position;

public:
    Atom(const int &atomType, const Position &position)
        : _atomicNumber(atomType), _position(position) {}
    Atom(const int &atomType, const double x, const double y, const double z)
        : _atomicNumber(atomType), _position(Position(x, y, z)) {}

    int atomicNumber() const { return _atomicNumber; }
    Position position() const { return _position; }
};

#endif
