#ifndef TREE_HPP
#define TREE_HPP

#include "graph.hpp"
#include "relation_set.hpp"
#include <iostream>
#include <array>

namespace boost {
	struct TreeProp { size_t root; };
	using Tree = adjacency_list<vecS, vecS, bidirectionalS, 
		  LabelProp, NoProp, TreeProp, vecS>;
}

namespace homology {
using Tree = boost::Tree;
using Node = Tree::vertex_descriptor;

inline Node root(const Tree& t) { return t[self].root; }

inline bool is_root(Node n, const Tree& t) {
	return t[self].root == n;
}

inline Node parent(Node n, const Tree& t) {
	assert(!is_root(n, t));
	return source(*boost::in_edges(n, t).first, t);
}

inline size_t height(Node n, const Tree& t) {
	size_t h = 0;
	while (!is_root(n, t)) { 
		++h;
		n = parent(n, t);
	}
	return h;
}

inline size_t num_children(Node n, const Tree& t) { return out_degree(n, t); }
inline auto children(Node n, const Tree& t) { return adjacent_vertices(n, t); }

inline void ordered_leaves(std::vector<Node>& out, Node n, const Tree& t) {
	if (num_children(n, t) == 0) {
		out.push_back(n);
		return;
	}
	for (auto e : out_edges(n, t)) {
		Node c = target(e, t);
		ordered_leaves(out, c, t);
	}
}

inline std::vector<Node> ordered_leaves(const Tree& t) {
	Node r = root(t);
	std::vector<Node> out;
	ordered_leaves(out, r, t);
	return out;
}

inline Node lca(Node n1, Node n2, const Tree& t) {
	assert(n1 != n2);
	assert(n1 < num_vertices(t) && n2 < num_vertices(t));
	size_t h1 = height(n1, t), h2 = height(n2, t);
	if (h1 < h2) {
		std::swap(h1,h2);
		std::swap(n1,n2);
	}
	while (h1 > h2) { n1 = parent(n1, t); --h1; }
	while (n1 != n2) {
		n1 = parent(n1, t);
		n2 = parent(n2, t);
	}
	return n1;
}

Tree collapse_tree(const Tree& t);
std::pair<Tree, bool> build_cotree(const RelationSet& rs);
std::pair<Tree, bool> build_cotree_dist(const RelationSet& rs, const std::array<double,3>& probs);

inline void print_tree_dot(const Tree& t, std::ostream& os = std::cout) {
	os << "graph {" << std::endl;
	for (Node n : vertices(t)) {
		os << n << "[label=" << t[n].label << "];" << std::endl;
	}
	for (auto e : edges(t)) {
		size_t src = source(e, t), tar = target(e, t);
		os << src << "--" << tar << ";" << std::endl;
	}
	os << "}" << std::endl;
}
	
} /* homology */ 

#endif /* ifndef TREE_HPP */
