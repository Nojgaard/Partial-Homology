#ifndef PARSE_ORTHO_HPP
#define PARSE_ORTHO_HPP

#include "../graph.hpp"
#include "../relation_set.hpp"
#include "../tree.hpp"
#include <iostream>
#include "../matrix.hpp"

namespace homology {
namespace parse {
namespace ortho {

struct Data {
	std::vector<RelationSet> relations;
	std::vector<std::vector<std::string>> matrices;
};

std::tuple<RelationSet, Tree, std::vector<std::string>, bool> read(std::istream& os);
std::pair<Matrix<int>, bool> read_org(std::istream& os);
void write_unknown(std::istream& os, std::ostream& oos, double up);
void read_data(std::istream& is, Data& out, int max_iter = 0);
	
} /* ortho */ 
} /* parse */ 
} /* homology */ 

#endif /* ifndef PARSE_ORTHO_HPP */
