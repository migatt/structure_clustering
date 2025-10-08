#ifndef VERTEXINVARIANT_HPP
#define VERTEXINVARIANT_HPP

#include <map>
#include <string>
#include <boost/graph/adjacency_list.hpp>

struct VertexInvariant {
    using Map = std::map<std::string, size_t>;
    const Graph &_graph;
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

    VertexInvariant(const Graph &graph, Map &mappings) : _graph(graph), _mappings(mappings) {}
};

#endif
