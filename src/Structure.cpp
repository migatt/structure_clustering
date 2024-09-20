#include <iostream>
#include <string>
#include <vector>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/graph/isomorphism.hpp>
#include <boost/graph/vf2_sub_graph_iso.hpp>

#include "Atom.hpp"
#include "Machine.hpp"
#include "Structure.hpp"

double distance(const Atom &atom_a, const Atom &atom_b) {
    return sqrt(pow(atom_b.position().x() - atom_a.position().x(), 2) +
                pow(atom_b.position().y() - atom_a.position().y(), 2) +
                pow(atom_b.position().z() - atom_a.position().z(), 2));
}

bool isConnected(const double &r_a, const double &r_b, const double &distance_ab) {
    return .8 < distance_ab && distance_ab < r_a + r_b + .4;
}

Structure::Structure(const int id) : _id(id), _graph({}) {}

void Structure::addAtom(const Atom &atom) { _atoms.push_back(atom); }

int Structure::numAtoms() const { return _atoms.size(); }
const Atom &Structure::getAtom(int index) const { return _atoms[index]; }
const int Structure::getNumConnections() const { return boost::num_edges(_graph); }
const Graph &Structure::getGraph() const { return _graph; }
const bool Structure::isGraphFullyConnected() const {
    std::vector<int> component(boost::num_vertices(_graph));
    int num = boost::connected_components(_graph, &component[0]);
    return num == 1;
}

void Structure::constructGraph(const Machine &machine) {
    // add graph vertices
    auto vertex_name_map = get(boost::vertex_name, _graph);
    std::vector<Graph::vertex_descriptor> vertex_descriptors(this->numAtoms());

    for (std::size_t i = 0; i < this->numAtoms(); ++i) {
        auto vd = boost::add_vertex(_graph);
        vertex_descriptors[i] = vd;
        vertex_name_map[vd] = this->getAtom(i).atomicNumber();
    }

    // check connectivity and add graph edges
    for (int i = 0; i < this->numAtoms(); i++) {
        for (int j = i + 1; j < this->numAtoms(); j++) {
            const Atom &atom_a = this->getAtom(i);
            const Atom &atom_b = this->getAtom(j);

            const int &atomicNumberA = atom_a.atomicNumber();
            const int &atomicNumberB = atom_b.atomicNumber();

            double distance_ab = distance(atom_a, atom_b);
            double r_a = machine.getCovalentRadius(atomicNumberA);
            double r_b = machine.getCovalentRadius(atomicNumberB);
            double maxPairDistance = machine.getMaxPairDistance(atomicNumberA, atomicNumberB);

            if (isConnected(r_a, r_b, distance_ab) || (distance_ab <= maxPairDistance)) {

                boost::add_edge(vertex_descriptors[i], vertex_descriptors[j], _graph);
            }
        }
    }
};
