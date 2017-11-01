#ifndef PARSE_NEWICK
#define PARSE_NEWICK 

#include "../tree.hpp"
#include <iostream>

namespace homology {
namespace parse {
namespace newick {

Tree read(const std::string& line);
void write(const Tree& t, bool strip_labels = false, std::ostream& os = std::cout);
	
} /* newick */ 
} /* parse */ 
} /* homology */ 
#endif /* ifndef PARSE_NEWICK */
