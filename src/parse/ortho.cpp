#include "ortho.hpp"
#include "newick.hpp"
#include "../matrix.hpp"
#include <sstream>
#include <assert.h>
#include <string>
#include <random>

namespace homology {
namespace parse {
namespace ortho {

void write_unknown(std::istream& os, std::ostream& oos, double up) {
	std::string line, tok;
	std::vector<std::string> matrix;
	bool success = false;
	while (std::getline(os, line)) {
		std::stringstream ss(line); ss >> tok;
		oos << line << std::endl;
		if (tok != "#GeneFamily") { 
			continue; 
		}
		success = true;
		std::vector<std::string> vertex_labels;
		std::getline(os, line); oos << line << std::endl;
		std::getline(os, line); oos << line << std::endl;
		/* std::cout << line << std::endl; */
		ss = std::stringstream(line);
		while ((ss >> tok)) { vertex_labels.emplace_back(tok); }
		std::getline(os, line); oos << line << std::endl;
		/* std::cout << line << std::endl; */
		/* t = newick::read(line); */
		std::getline(os, line); oos << line << std::endl;
		while (line != "#OrthologyMatrix") { std::getline(os, line); oos << line << std::endl; }
		/* std::vector<std::string> matrix(vertex_labels.size()); */
		matrix = std::vector<std::string>(vertex_labels.size());
		for (size_t i = 0; i < vertex_labels.size(); ++i) {
			std::getline(os, line);
			ss = std::stringstream(line);
			while ((ss >> tok)) { matrix[i] += tok; }
		}

		std::random_device rd;
		std::mt19937 e2(rd());
		std::uniform_real_distribution<> dist(0, 1);
		for (size_t i = 0; i < vertex_labels.size(); ++i) {
			for (size_t j = i + 1; j < vertex_labels.size(); ++j) {
				double p = dist(e2);
				if (p < up) {
					matrix[i][j] = '-';
					matrix[j][i] = '-';
				}
			}
		}
		for (size_t i = 0; i < vertex_labels.size(); ++i) {
			for (size_t j = 0; j < vertex_labels.size(); ++j) {
				oos << matrix[i][j] << " ";
			}
			oos << std::endl;
		}
		oos << std::endl;
		/* break; */
	}
}

std::pair<Matrix<int>, bool> read_org(std::istream& os) {
	std::string line, tok;
	std::vector<std::string> matrix;
	bool success = false;
	while (std::getline(os, line)) {
		std::stringstream ss(line); ss >> tok;
		if (tok != "#GeneFamily") { continue; }
		success = true;
		std::vector<std::string> vertex_labels;
		std::getline(os, line);
		std::getline(os, line);
		/* std::cout << line << std::endl; */
		ss = std::stringstream(line);
		while ((ss >> tok)) { vertex_labels.emplace_back(tok); }
		std::getline(os, line);
		/* std::cout << line << std::endl; */
		/* t = newick::read(line); */
		std::getline(os, line);
		while (line != "#OrthologyMatrix") { std::getline(os, line); }
		/* std::vector<std::string> matrix(vertex_labels.size()); */
		matrix = std::vector<std::string>(vertex_labels.size());
		for (size_t i = 0; i < vertex_labels.size(); ++i) {
			std::getline(os, line);
			ss = std::stringstream(line);
			while ((ss >> tok)) { matrix[i] += tok; }
		}
		Matrix<int> ortho_matrix(vertex_labels.size(), vertex_labels.size());
		for (size_t i = 0; i < vertex_labels.size(); ++i) { ortho_matrix(i, i) = 0; }
		for (size_t i = 0; i < vertex_labels.size(); ++i) {
			for (size_t j = i + 1; j < vertex_labels.size(); ++j) {
				if (matrix[i][j] == '-') { 
					assert("has unknown" && false);
				}
				if (matrix[i][j] == '1' && matrix[j][i] == '1') {
					ortho_matrix(i,j) = 1;
					ortho_matrix(j,i) = 1;
				} else if (matrix[i][j] == '0' && matrix[j][i] == '0') {
					ortho_matrix(i,j) = 0;
					ortho_matrix(j,i) = 0;
				} else if(matrix[i][j] == '1') {
					ortho_matrix(i,j) = 1;
					ortho_matrix(j,i) = 0;
				} else {
					ortho_matrix(i,j) = 0;
					ortho_matrix(j,i) = 1;
				}
			}
		}
		/* for (const std::string& s : matrix) { */
		/* 	std::cout << "("; */
		/* 	for (size_t i = 0; i < s.size(); ++i) { */
		/* 		std::cout << s[i]; */
		/* 		if (i < s.size() - 1) { std::cout << " "; } */
		/* 	} */
		/* 	std::cout << ")" << std::endl; */
		/* 	/1* std::cout << s << std::endl; *1/ */
		/* } */
		return std::make_pair(ortho_matrix, success);
		break;
	}

	return std::make_pair(Matrix<int>(0,0), success);
}

void read_data(std::istream& is, Data& out, int max_iter) {
	std::string line, tok;
	while (std::getline(is, line)) {
		std::stringstream ss(line); ss >> tok;
		if (tok != "#GeneFamily") { continue; }
		std::vector<std::string> vertex_labels;
		std::getline(is, line);
		std::getline(is, line);
		/* std::cout << line << std::endl; */
		ss = std::stringstream(line);
		while ((ss >> tok)) { vertex_labels.emplace_back(tok); }
		std::getline(is, line);
		/* std::cout << line << std::endl; */
		auto& rs = out.relations.emplace_back(vertex_labels);
		auto O = rs.add_sym_relation("S").first;
		auto P = rs.add_sym_relation("D").first;
		auto X = rs.add_asym_relation("L").first;
		while (line != "#OrthologyMatrix") { std::getline(is, line); }
		/* std::vector<std::string> matrix(vertex_labels.size()); */
		auto& matrix = out.matrices.emplace_back(vertex_labels.size());
		for (size_t i = 0; i < vertex_labels.size(); ++i) {
			std::getline(is, line);
			ss = std::stringstream(line);
			while ((ss >> tok)) { matrix[i] += tok; }
		}
		for (size_t i = 0; i < vertex_labels.size(); ++i) {
			for (size_t j = i + 1; j < vertex_labels.size(); ++j) {
				if (matrix[i][j] == '-') { continue; }
				if (matrix[i][j] == '1' && matrix[j][i] == '1') {
					rs.assign(i, j, O);
				} else if (matrix[i][j] == '0' && matrix[j][i] == '0') {
					rs.assign(i, j, P);
				} else if(matrix[i][j] == '1') {
					rs.assign(i, j, X);
				} else {
					rs.assign(j, i, X);
				}
			}
		}
		/* for (const std::string& s : matrix) { */
		/* 	std::cout << "("; */
		/* 	for (size_t i = 0; i < s.size(); ++i) { */
		/* 		std::cout << s[i]; */
		/* 		if (i < s.size() - 1) { std::cout << " "; } */
		/* 	} */
		/* 	std::cout << ")" << std::endl; */
		/* 	/1* std::cout << s << std::endl; *1/ */
		/* } */
		if (out.relations.size() == (size_t)max_iter) {
			break;
		}
	}
}

std::tuple<RelationSet, Tree, std::vector<std::string>, bool> read(std::istream& os) {
	RelationSet rs(0);
	std::string line, tok;
	Tree t;
	std::vector<std::string> matrix;
	bool success = false;
	while (std::getline(os, line)) {
		std::stringstream ss(line); ss >> tok;
		if (tok != "#GeneFamily") { continue; }
		success = true;
		std::vector<std::string> vertex_labels;
		std::getline(os, line);
		std::getline(os, line);
		/* std::cout << line << std::endl; */
		ss = std::stringstream(line);
		while ((ss >> tok)) { vertex_labels.emplace_back(tok); }
		std::getline(os, line);
		/* std::cout << line << std::endl; */
		t = newick::read(line);
		rs = RelationSet(vertex_labels);
		auto O = rs.add_sym_relation("S").first;
		auto P = rs.add_sym_relation("D").first;
		auto X = rs.add_asym_relation("L").first;
		while (line != "#OrthologyMatrix") { std::getline(os, line); }
		/* std::vector<std::string> matrix(vertex_labels.size()); */
		matrix = std::vector<std::string>(vertex_labels.size());
		for (size_t i = 0; i < vertex_labels.size(); ++i) {
			std::getline(os, line);
			ss = std::stringstream(line);
			while ((ss >> tok)) { matrix[i] += tok; }
		}
		for (size_t i = 0; i < vertex_labels.size(); ++i) {
			for (size_t j = i + 1; j < vertex_labels.size(); ++j) {
				if (matrix[i][j] == '-') { continue; }
				if (matrix[i][j] == '1' && matrix[j][i] == '1') {
					rs.assign(i, j, O);
				} else if (matrix[i][j] == '0' && matrix[j][i] == '0') {
					rs.assign(i, j, P);
				} else if(matrix[i][j] == '1') {
					rs.assign(i, j, X);
				} else {
					rs.assign(j, i, X);
				}
			}
		}
		/* for (const std::string& s : matrix) { */
		/* 	std::cout << "("; */
		/* 	for (size_t i = 0; i < s.size(); ++i) { */
		/* 		std::cout << s[i]; */
		/* 		if (i < s.size() - 1) { std::cout << " "; } */
		/* 	} */
		/* 	std::cout << ")" << std::endl; */
		/* 	/1* std::cout << s << std::endl; *1/ */
		/* } */
		break;
	}

	assert(!success || (rs.sym_relations().size() == 2 && rs.asym_relations().size() == 1));
	return std::make_tuple(rs, t, matrix, success);
}
	
} /* ortho */ 
} /* parse */ 
} /* homology */ 
