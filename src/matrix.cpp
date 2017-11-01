#include "matrix.hpp"
#include <vector>
#include <map>

namespace homology {

std::map<std::string, size_t> build_ortho_map(const Tree & t) {
	std::vector<Node> ol = ordered_leaves(t);
	std::vector<size_t> leaf_map;
	for (size_t i = 0; i < ol.size(); ++i) { leaf_map.push_back(i); }
	std::sort(leaf_map.begin(), leaf_map.end(), 
			[&](size_t lhs, size_t rhs) -> bool {
				return t[ol[lhs]].label < t[ol[rhs]].label;
			});
	std::map<std::string, size_t> out;
	for (size_t i = 0; i < ol.size(); ++i) {
		out[t[ol[i]].label] = leaf_map[i];
	}
	return out;
}

std::vector<std::string> build_matrix(const Tree& t) {
	std::vector<Node> ol = ordered_leaves(t);
	std::vector<size_t> leaf_map;
	for (size_t i = 0; i < ol.size(); ++i) { leaf_map.push_back(i); }
	std::sort(leaf_map.begin(), leaf_map.end(), 
			[&](size_t lhs, size_t rhs) -> bool {
				return t[ol[lhs]].label < t[ol[rhs]].label;
			});
	std::vector<std::string> m(ol.size(), std::string(ol.size(), '0'));
	for (size_t i = 0; i < ol.size(); ++i) {
		for (size_t j = i+1; j < ol.size(); ++j) {
			size_t x = leaf_map[i], y = leaf_map[j];
			/* if (x < y || i == j) { continue; } */
			auto label = t[lca(ol[x], ol[y], t)].label;
			if (label == "S") {
				m[i][j] = '1'; m[j][i] = '1';
			} else if (label == "D") {
				m[i][j] = '0'; m[j][i] = '0';
			} else {
				if (x < y) { m[i][j] = '1'; m[j][i] = '0'; }
				else { m[i][j] = '0'; m[j][i] = '1'; }
			}
		}
	}
	return m;
}

Matrix<int> build_ortho_matrix(const Tree& t) {
	std::vector<Node> ol = ordered_leaves(t);
	std::vector<size_t> leaf_map;
	for (size_t i = 0; i < ol.size(); ++i) { leaf_map.push_back(i); }
	std::sort(leaf_map.begin(), leaf_map.end(), 
			[&](size_t lhs, size_t rhs) -> bool {
				return t[ol[lhs]].label < t[ol[rhs]].label;
			});
	Matrix<int> m(ol.size(), ol.size());
	for (size_t i = 0; i < ol.size(); ++i) { m(i, i) = 0; }
	for (size_t i = 0; i < ol.size(); ++i) {
		for (size_t j = i+1; j < ol.size(); ++j) {
			size_t x = leaf_map[i], y = leaf_map[j];
			/* if (x < y || i == j) { continue; } */
			auto label = t[lca(ol[x], ol[y], t)].label;
			if (label == "S") {
				m(x, y) = 1; m(y, x) = 1;
			} else if (label == "D") {
				m(x, y) = 0; m(y, x) = 0;
			} else {
				if (x < y) { m(x, y) = 1; m(y, x) = 0; }
				else { m(x, y) = 0; m(y, x) = 1; }
			}
			/* auto label = t[lca(ol[i], ol[j], t)].label; */
			/* if (label == "S") { */
			/* 	m(i, j) = 1; m(j, i) = 1; */
			/* } else if (label == "D") { */
			/* 	m(i, j) = 0; m(j, i) = 0; */
			/* } else { */
			/* 	m(i, j) = 1; m(j, i) = 0; */
			/* } */
		}
	}
	return m;
}

int pairwise_diff(const Matrix<int>& m1, const Matrix<int>& m2) {
	assert(m1.size1() == m2.size1() && m1.size2() == m2.size2());
	if (m1.size1() == 1) { return 0; }
	int diff = 0;
	size_t n = m1.size1(), m = m1.size2();
	for (size_t i = 0; i < n; ++i) {
		for (size_t j = i + 1; j < m; ++j) {
			if (m1(i,j) != m2(i,j) || m1(j,i) != m2(j,i)) {
				diff += 1;
			}
		}
	}

	return diff;
}

std::pair<int,int> compute_tpfn(const Matrix<int>& m1, const Matrix<int>& m2, std::pair<int,int> pat) {
	int tp = 0, fn = 0;
	size_t n = m1.size1(), m = m1.size2();
	for (size_t i = 0; i < n; ++i) {
		for (size_t j = i + 1; j < m; ++j) {
			int x = std::max(m1(i,j), m1(j,i)), y = std::min(m1(i,j), m1(j,i));
			if (x != pat.first || y != pat.second) { continue; }

			if (m1(i,j) == m2(i,j) && m1(j,i) == m2(j,i)) {
				++tp;
			} else {
				++fn;
			}
		}
	}
	return std::make_pair(tp,fn);
}

std::array<std::pair<double,double>, 3> tpfn_rates(const Matrix<int>& m1, const Matrix<int>& m2) {
	std::array<std::pair<double,double>, 3> out;
	out[0] = std::make_pair(0,0); out[1] = std::make_pair(0,0); out[2] = std::make_pair(0,0);
	if (m1.size1() == 1) { 
		return out;; 
	}
	auto [stp, sfn] = compute_tpfn(m1, m2, std::make_pair(1, 1));
	auto [dtp, dfn] = compute_tpfn(m1, m2, std::make_pair(0, 0));
	auto [htp, hfn] = compute_tpfn(m1, m2, std::make_pair(1, 0));
	if (stp + sfn > 0) {
		out[0] = std::make_pair(((double)stp)/(stp+sfn), ((double)sfn)/(stp+sfn));
	}
	if (dtp + dfn > 0) {
		out[1] = std::make_pair(((double)dtp)/(dtp+dfn), ((double)dfn)/(dtp+dfn));
	}
	if (htp + hfn > 0) {
		out[2] = std::make_pair(((double)htp)/(htp+hfn), ((double)hfn)/(htp+hfn));
	}
	return out;
}

double relative_diff(const Matrix<int>& m1, const Matrix<int>& m2) {
	if (m1.size1() == 1) { return 0; }
	int diff = pairwise_diff(m1, m2);
	size_t n = m1.size1();
	return diff / ((n*n-n)/2.);
}

int pairwise_diff(const std::vector<std::string>& m1, 
		const std::vector<std::string>& m2) {
	assert(m1.size() == m2.size());
	if (m1.size() == 1) { return 0; }
	int diff = 0;
	size_t n = m1.size();
	for (size_t i = 0; i < n; ++i) {
		for (size_t j = i + 1; j < n; ++j) {
			if (m1[i][j] != m2[i][j] || m1[j][i] != m2[j][i]) {
				diff += 1;
			}
		}
	}
	return diff;
}
double relative_diff(const std::vector<std::string>& m1, 
		const std::vector<std::string>& m2) {
	assert(m1.size() == m2.size());
	if (m1.size() == 1) { return 0; }
	int diff = pairwise_diff(m1, m2);
	size_t n = m1.size();
	return diff / ((n*n-n)/2.);
}



double relative_unknown_score(const std::vector<std::string>& om,
		const Matrix<int>& m1, const Matrix<int>& m2) {
	assert(m1.size1() == m2.size1() && m1.size2() == m2.size2());
	if (m1.size1() == 1) { return 0; }
	int diff = 0, num_unknowns = 0;
	size_t n = m1.size1(), m = m1.size2();
	for (size_t i = 0; i < n; ++i) {
		for (size_t j = i + 1; j < m; ++j) {
			if (m1(i,j) != m2(i,j) || m1(j,i) != m2(j,i)) {
				assert(om[i][j] == '-');
				diff += 1;
			}
			if (om[i][j] == '-') { num_unknowns += 1; }
		}
	}

	double res = ((double)diff) / num_unknowns;
	return (num_unknowns > 0)?res:0;
}
	
} /* homology */ 
