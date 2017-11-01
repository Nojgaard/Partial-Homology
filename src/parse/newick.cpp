#include "../tree.hpp"

namespace homology {
namespace parse {
namespace newick {

std::string parse_label(size_t &i, const std::string& line) {
	std::string lbl;
	while (line[i] != ',' && line[i] != ')' && line[i] != ';') {
		lbl += line[i];
		++i;
	}
	return lbl;
}

Tree read(const std::string& line) {
	Tree t;
	t[self].root = 0;

	Node cur = 0;
	size_t i = 0;
	while (i < line.size() && line[i] != '=') { ++i; }
	++i;
	while (i < line.size()) {
		char c = line[i];
		if (c == '(') {
			Node n = add_vertex(t);
			if (n != 0) {
				add_edge(cur, n, t);
			}
			cur = n;
			++i;
		} else if (c == ')') {
			t[cur].label = line[i + 1];
			if (cur != 0) {
				cur = parent(cur, t);
			}
			i += 2;
		} else if (c == ',') {
			++i;
		} else if (c == ';') {
			break;
		} else {
			Node n = add_vertex(t);
			t[n].label = parse_label(i, line);
			if (n != 0) {
				add_edge(cur, n, t);
			}
		}
	}
	return t;
}

void write(Node n, const Tree& t, bool strip_labels, std::ostream& os) {
	if (num_children(n, t) == 0) {
		os << t[n].label;
		return;
	}
	os << "(";
	size_t i = 0;
	for (Node c : children(n, t)) {
		write(c, t, strip_labels, os);
		if (i < num_children(n, t)-1) { os << ","; }
		i += 1;
	}
	os << ")";
	if (!strip_labels) {
		os << t[n].label;
	}
}

void write(const Tree& t, bool strip_labels, std::ostream& os) {
	write(root(t), t, strip_labels, os);
	os << ";";
}

} /* newick */ 
} /* parse */ 
} /* homology */ 
