#ifndef __STRUCTURE__
#define __STRUCTURE__

#include <iostream>
#include <string>
#include <vector>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/graph/isomorphism.hpp>
#include <boost/graph/vf2_sub_graph_iso.hpp>

#include "Atom.hpp"
#include "Machine.hpp"

class Machine; // forward declaration

using Graph = boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS,
                                    boost::property<boost::vertex_name_t, std::string>>;

namespace {
    struct VertexInvariant {
        using Map = std::map<std::string, size_t>;
        Graph const &_graph;
        Map &_mappings;

        using result_type = size_t;
        using argument_type = Graph::vertex_descriptor;

        size_t operator()(argument_type u) const {
            return _mappings.at(boost::get(boost::vertex_name, _graph, u));
        }
        size_t max() const { return _mappings.size(); }

        void collect_names() {
            for (auto vd : boost::make_iterator_range(boost::vertices(_graph))) {
                size_t next_id = _mappings.size();
                _mappings.insert({boost::get(boost::vertex_name, _graph, vd), next_id});
            }
        }
    };
} // namespace

class Structure {
    int _id;
    std::vector<Atom> _atoms;
    mutable std::optional<VertexInvariant> _vertexInvariant;
    mutable std::optional<std::vector<int>> _degreeSequence;
    Graph _graph;


public:
    Structure(const int id);

    void addAtom(const Atom &atom);

    int numAtoms() const;
    const Atom &getAtom(int index) const;
    const int getNumConnections() const;
    const Graph &getGraph() const;
    const bool isGraphFullyConnected() const;
    const VertexInvariant &getVertexInvariant(VertexInvariant::Map &shared_names) const;
    const std::vector<int> &getDegreeSequence() const;
    void constructGraph(const Machine &machine);
};

#endif
