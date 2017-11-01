#include "relation_set.hpp"
#include <iostream>
#include <math.h>
#include <random>

namespace homology {
	

RelationSet::RelationSet(size_t num_vertices): num_vertices(num_vertices) {
	for (size_t i = 0; i < num_vertices; ++i) {
		vertex_labels_.push_back(std::to_string(i));
	}
}

RelationSet::RelationSet(const Matrix<int>& ortho_matrix, double unknown_prop):
		num_vertices(ortho_matrix.size1()) {
	std::random_device rd;
    std::mt19937 e2(rd());
    std::uniform_real_distribution<> dist(0, 1);
	for (size_t i = 0; i < num_vertices; ++i) {
		vertex_labels_.push_back(std::to_string(i));
	}
	std::cout << "UNKNOWN_PROP: " << unknown_prop << std::endl;

	auto O = add_sym_relation("S").first;
	auto P = add_sym_relation("D").first;
	auto X = add_asym_relation("L").first;
	for (size_t i = 0; i < num_vertices; ++i) {
	for (size_t j = i + 1; j < num_vertices; ++j) {
		double p = dist(e2);
		if (p > unknown_prop) { 
			continue; 
		}
		if (ortho_matrix(i, j) == 1 && ortho_matrix(j, i) == 1) {
			assign(i, j, O);
		} else if (ortho_matrix(i, j) == 0 && ortho_matrix(j, i) == 0) {
			assign(i, j, P);
		} else if (ortho_matrix(i, j) == 1) {
			assign(i, j, X);
		} else {
			assign(j, i, X);
		}
	}
	}
}

RelationSet::RelationSet(const std::vector<std::string>& vertex_labels)
		: num_vertices(vertex_labels.size()), vertex_labels_(vertex_labels) { }

std::pair<SymRelId, bool> RelationSet::add_sym_relation(std::string name) {
	if (used_labels_.count(name)) { return std::make_pair(SymRelId{0}, false); }
	SymRelId id{symmetric_relations_.size()};
	symmetric_relations_.push_back(UndirectedGraph(num_vertices));
	symmetric_relations_.back()[self].label = name;
	forbidden_relations_.push_back(UndirectedGraph(num_vertices));
	forbidden_relations_.back()[self].label = name;
	return std::make_pair(id, true);
}

std::pair<AsymRelId, bool> RelationSet::add_asym_relation(std::string name) {
	if (used_labels_.count(name)) { return std::make_pair(AsymRelId{0}, false); }
	AsymRelId id{asymmetric_relations_.size()};
	asymmetric_relations_.push_back(DirectedGraph(num_vertices));
	asymmetric_relations_.back()[self].label = name;
	forbidden_relations_.push_back(UndirectedGraph(num_vertices));
	forbidden_relations_.back()[self].label = name;
	return std::make_pair(id, true);
}

bool RelationSet::add_edge(size_t src, size_t tar) {
	assert(src < num_vertices && tar < num_vertices);
	auto e = std::make_pair(src, tar), re = std::make_pair(tar, src);
	if (used_edges_.count(e) || used_edges_.count(re)) { return false; }

	used_edges_.insert(e);
	used_edges_.insert(re);
	return true;
}

bool RelationSet::is_edge_used(size_t src, size_t tar) const {
	assert(src < num_vertices && tar < num_vertices);
	auto e = std::make_pair(src, tar), re = std::make_pair(tar, src);
	return used_edges_.count(e) || used_edges_.count(re);
}

bool RelationSet::assign(size_t src, size_t tar, SymRelId rel) {
	if (add_edge(src, tar)) {
		boost::add_edge(src, tar, symmetric_relations_[rel.id]);
		return true;
	} else {
		return false;
	}
}

bool RelationSet::assign(size_t src, size_t tar, AsymRelId rel) {
	if (add_edge(src, tar)) {
		boost::add_edge(src, tar, asymmetric_relations_[rel.id]);
		return true;
	} else {
		return false;
	}
}

bool RelationSet::forbid(size_t src, size_t tar, size_t frel_id) {
	assert(frel_id < forbidden_relations_.size());
	assert(src < num_vertices && tar < num_vertices);
	auto e = std::make_pair(src, tar), re = std::make_pair(tar, src);
	if (used_edges_.count(e) || used_edges_.count(re)) { return false; }
	auto pair_edge = boost::add_edge(src, tar, forbidden_relations_[frel_id]);
	return pair_edge.second;
}

size_t RelationSet::num_relations() const { 
	return num_sym_relations() + num_asym_relations(); 
}

size_t RelationSet::num_sym_relations() const { 
	return symmetric_relations_.size(); 
}

size_t RelationSet::num_asym_relations() const { 
	return asymmetric_relations_.size(); 
}

const std::vector<UndirectedGraph>& RelationSet::forbidden_relations() const {
	return forbidden_relations_;
}

const std::vector<UndirectedGraph>& RelationSet::sym_relations() const {
	return symmetric_relations_;
}

const std::vector<DirectedGraph>& RelationSet::asym_relations() const {
	return asymmetric_relations_;
}

const std::string& RelationSet::vertex_label(size_t v) const {
	assert(v < vertex_labels_.size());
	return vertex_labels_[v];
}

void RelationSet::print_dot(std::ostream& os) const {
	os << "digraph Relations {" << std::endl;
	int n = ceil(sqrt(num_vertices));
	int i = 0;
	for (const auto& g : sym_relations()) {
		size_t offset = i * num_vertices;
		os << "subgraph cluster_" << i++ << "{" << std::endl;
		os << "edge [dir=none]" << std::endl;
		print_dot_relation(g, offset, n, os);
		os << "}" << std::endl;
	}

	for (const auto& g : asym_relations()) {
		size_t offset = i * num_vertices;
		os << "subgraph cluster_" << i++ << "{" << std::endl;
		print_dot_relation(g, offset, n, os);
		os << "}" << std::endl;
	}
	os << "}" << std::endl;
}	

} /* homology */ 
