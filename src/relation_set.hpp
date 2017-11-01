#ifndef RELATION_SET
#define RELATION_SET

#include "graph.hpp"
#include <vector>
#include <set>
#include <iostream>
/* #include "matrix.hpp" */
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp> 

namespace homology {
template<typename T> using Matrix = boost::numeric::ublas::matrix<T>;
	

struct SymRelId { size_t id; };
struct AsymRelId { size_t id; };

class RelationSet {
public:
	RelationSet(size_t num_vertices);
	RelationSet(const Matrix<int>& ortho_matrix, double unknown_prop);
	RelationSet(const std::vector<std::string>& vertex_labels);
	std::pair<SymRelId, bool> add_sym_relation(std::string name);
	std::pair<AsymRelId, bool> add_asym_relation(std::string name);
	bool assign(size_t src, size_t tar, SymRelId rel);
	bool assign(size_t src, size_t tar, AsymRelId rel);
	bool forbid(size_t src, size_t tar, size_t frel_id);
	const std::string& vertex_label(size_t v) const;
	size_t num_relations() const;
	size_t num_sym_relations() const;
	size_t num_asym_relations() const;
	bool is_edge_used(size_t src, size_t tar) const;

	const std::vector<UndirectedGraph>& sym_relations() const;
	const std::vector<UndirectedGraph>& forbidden_relations() const;
	const std::vector<DirectedGraph>& asym_relations() const;

	void print_dot(std::ostream& os = std::cout) const;
	size_t num_vertices;

private:
	template<typename R>
	void print_dot_relation(const R& g, size_t offset, int n, std::ostream& os) const {
		os << "label=" << g[self].label << std::endl;
		/* os << "node [shape=point]" << std::endl; */
		int row = 0, col = 0;
		for (auto v : vertices(g)) {
			os << (v + offset) << "[shape=point xlabel=" << v;
			os << " pos=\"" << (row/2.) << "," << (col/2.) << "!\"";
			col = col + 1;
			if (col == n) {
				row += 1;
				col = 0;
			}
			os << "];" << std::endl;
		}
		for (auto e : edges(g)) {
			size_t src = source(e, g) + offset, tar = target(e, g) + offset;
			os << src << " -> " << tar << ";" << std::endl;
		}
	}
	bool add_edge(size_t src, size_t tar);

	std::vector<UndirectedGraph> symmetric_relations_;
	std::vector<DirectedGraph> asymmetric_relations_;
	std::vector<UndirectedGraph> forbidden_relations_;
	std::vector<std::string> vertex_labels_;
	std::set<std::pair<size_t, size_t>> used_edges_;
	std::set<std::string> used_labels_;
};

} /* homology */ 

#endif /* ifndef RELATION_SET */

