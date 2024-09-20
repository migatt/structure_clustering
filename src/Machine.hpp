#ifndef __MACHINE__
#define __MACHINE__

#include <iostream>
#include <vector>

#include "Result.hpp"
#include "Structure.hpp"

class Structure;
class Result;

class Machine {
    bool _onlyConnectedGraphs;
    std::vector<double> _covalentRadii; // index corresponds to (atomic number - 1)
    std::map<std::pair<int, int>, double> _pairDistances;

public:
    Machine();

    void setOnlyConnectedGraphs(bool onlyConnectedGraphs);
    void setCovalentRadius(int atomicNumber, double radius);
    void addPairDistance(int atomicNumberA, int atomicNumberB, double maxDistance);

    bool isOnlyConnectedGraphs() const;
    double getCovalentRadius(int atomicNumber) const;
    double getMaxPairDistance(int atomicNumberA, int atomicNumberB) const;

    Result cluster(std::vector<Structure> &structures);
};

#endif
