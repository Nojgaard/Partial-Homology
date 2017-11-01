#ifndef GRAPH_HPP
#define GRAPH_HPP
#include <boost/graph/subgraph.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/strong_components.hpp>
#include <boost/graph/topological_sort.hpp>
#include <iostream>

namespace boost {
using NoProp = no_property;
struct LabelProp { std::string label; };

using UndirectedGraph = adjacency_list<vecS, vecS,
	  undirectedS, LabelProp, property<edge_index_t, int>, LabelProp, vecS>;

using DirectedGraph = adjacency_list<vecS, vecS, directedS,
	  LabelProp, property<edge_index_t, int>, LabelProp, vecS>;

} /* boost */ 

namespace homology {
	
using UndirectedGraph = boost::UndirectedGraph;
using DirectedGraph = boost::DirectedGraph;
template <typename G> using SubGraph = boost::subgraph<G>;

#define self boost::graph_bundle
/* using self = boost::graph_bundle; */

template<typename G> auto vertices(const G& g) {
	return boost::make_iterator_range(boost::vertices(g));
}

template<typename G> auto edges(const G& g) {
	return boost::make_iterator_range(boost::edges(g));
}

template<typename G, typename V> auto adjacent_vertices(V v, const G& g) {
	return boost::make_iterator_range(boost::adjacent_vertices(v, g));
}

template<typename G, typename V> auto out_edges(V v, const G& g) {
	return boost::make_iterator_range(boost::out_edges(v, g));
}

template<typename G, typename V> auto in_edges(V v, const G& g) {
	return boost::make_iterator_range(boost::in_edges(v, g));
}

template<typename G>
std::pair<size_t, std::vector<int>> connected_components(const G& g) {
	std::vector<int> components(num_vertices(g));
	size_t num_comp = boost::connected_components(g, &components[0]);
	return std::make_pair(num_comp, components);
}

template<typename G>
std::pair<size_t, std::vector<int>> strong_connected_components(const G& g) {
	std::vector<int> components(num_vertices(g));
	size_t num_comp = boost::strong_components(g, &components[0]);
	return std::make_pair(num_comp, components);
}

//is sorted in reversed
template<typename G>
std::vector<size_t> topo_sort(const G& g) {
	std::vector<size_t> out;
	boost::topological_sort(g, std::back_inserter(out));
	return out;
}

template<typename GI, typename GO>
void merge(const GI& g, GO& out_g) {
	for (auto e : edges(g)) {
		size_t src = source(e, g), tar = target(e, g);
		assert(src < num_vertices(out_g) && tar < num_vertices(out_g));
		if (!boost::edge(src, tar, out_g).second) { add_edge(src, tar, out_g); }
	}
}

template<typename G, typename GO>
inline void merge_as_bidir(const G& g, GO& out_g) {
	for (auto e : edges(g)) {
		size_t src = source(e, g), tar = target(e, g);
		assert(src < num_vertices(out_g) && tar < num_vertices(out_g));
		if (!boost::edge(src, tar, out_g).second) { add_edge(src, tar, out_g); }
		if (!boost::edge(tar, src, out_g).second) { add_edge(tar, src, out_g); }
	}
}

} /* homology */ 

#endif /* ifndef GRAPH_HPP */
