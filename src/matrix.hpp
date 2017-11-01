#ifndef HOMOLOGY_MATRIX_HPP
#define HOMOLOGY_MATRIX_HPP

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp> 
#include "tree.hpp"
#include <iostream>

namespace homology {

template<typename T> using Matrix = boost::numeric::ublas::matrix<T>;

std::vector<std::string> build_matrix(const Tree& t);
Matrix<int> build_ortho_matrix(const Tree& t);
std::map<std::string, size_t> build_ortho_map(const Tree & t);
int pairwise_diff(const Matrix<int>& m1, const Matrix<int>& m2);
double relative_diff(const Matrix<int>& m1, const Matrix<int>& m2);
double relative_unknown_score(const std::vector<std::string>& om, 
		const Matrix<int>& m1, const Matrix<int>& m2);
std::array<std::pair<double,double>, 3> tpfn_rates(const Matrix<int>& m1, const Matrix<int>& m2);
double relative_diff(const std::vector<std::string>& m1, 
		const std::vector<std::string>& m2);
template<typename T>
void print(const Matrix<T>& m, std::ostream& os = std::cout) {
	os << "[" << m.size1() << "," << m.size2() << "]" << std::endl;
	for (size_t i = 0; i < m.size1(); ++i) {
		os << "(";
		for (size_t j = 0; j < m.size2(); ++j) {
			os << m(i,j);
			if (j < m.size2() - 1) { os << " "; }
		}
		os << ")" << std::endl;
	}
}

inline void print(const std::vector<std::string>& m, std::ostream& os = std::cout) {
	size_t n = m.size();
	os << "[" << n << "]" << std::endl;
	for (size_t i = 0; i < n; ++i) {
		os << "(";
		for (size_t j = 0; j < n; ++j) {
			os << m[i][j];
			if (j < n - 1) { os << " "; }
		}
		os << ")" << std::endl;
	}
}
	
} /* homology */ 

#endif /* ifndef HOMOLOGY_MATRIX_HPP */
