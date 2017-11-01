#include <iostream>
#include "graph.hpp"
#include "relation_set.hpp"
#include "tree.hpp"
#include "parse/ortho.hpp"
#include "matrix.hpp"
#include <fstream>
#include <cxxopts.hpp>
#include <sstream>
#include "parse/newick.hpp"
#include <array>
#include <random>

using namespace homology;

std::vector<std::string> split(const std::string& s, char delim) {
	std::stringstream ss(s);
	std::string item;
	std::vector<std::string> out;
	while (std::getline(ss, item, delim)) { out.push_back(item); }
	return out;
}

enum class Score { matrix_diff, unknown_diff };

void print_trees(std::string p, bool strip_labels = false) {
	std::ifstream input(p);
	RelationSet rs(0);
	Tree org_tree;
	bool success;
	std::vector<std::string> om;
	std::tie(rs, org_tree, om, success) = parse::ortho::read(input);
	while (success) {
		auto [t, has_tree] = build_cotree(rs);
		assert(has_tree);

		/* parse::newick::write(org_tree, strip_labels); */
		/* print_tree_dot(org_tree); */
		org_tree = collapse_tree(org_tree);
		/* print_tree_dot(org_tree); */
		parse::newick::write(org_tree, strip_labels);
		parse::newick::write(t, strip_labels);
		std::cout << std::endl;
		std::tie(rs, org_tree, om, success) = parse::ortho::read(input);
	}
}

double compute_diff(std::string p, Score s, std::array<double, 3>& probs) {
	std::ifstream input(p);
	RelationSet rs(0);
	Tree org_tree;
	bool success;
	std::vector<std::string> om;
	std::tie(rs, org_tree, om, success) = parse::ortho::read(input);
	size_t i = 0;
	double res = 0;
	while (success) {
		i += 1;
		auto [t, has_tree] = build_cotree_dist(rs, probs);
		/* auto [t, has_tree] = build_cotree(rs); */
		assert(has_tree);
		/* std::cout << num_vertices(org_tree) << std::endl; */

		auto m1 = build_ortho_matrix(org_tree);
		/* print(m1); */
		auto m2 = build_ortho_matrix(t);
		/* print(m2); */
		/* std::cout << relative_diff(m1, m2) << std::endl; */
		if (s == Score::matrix_diff) {
			res += relative_diff(m1, m2);
		} else if (s == Score::unknown_diff) {
			res += relative_unknown_score(om, m1, m2);
		}
		/* double d = relative_diff(m1, m2); */
		/* double ud = relative_unknown_score(om, m1, m2); */
		/* res += d; */
		/* res_unknown += ud; */
		/* std::cout << i << " " << res << std::endl; */
		std::tie(rs, org_tree, om, success) = parse::ortho::read(input);
	}
	/* std::cout << "AVG REL DIFF: " << res << " " << (i-1) << " " << (res/(i-1)) << " " << (res_unknown/(i-1)) <<  std::endl; */
	return res/i;
}

void restrict_relations(double forbid_p, const Tree& t, const
		std::vector<std::string>& strm,
		RelationSet& rs) {
	std::random_device rd;
    std::mt19937 e2(rd());
    std::uniform_real_distribution<> dist(0, 1);
	auto om = build_ortho_matrix(t);
	auto omap = build_ortho_map(t);
	assert(om.size1() == rs.num_vertices);

	for (size_t i = 0; i < rs.num_vertices; ++i) {
		for (size_t j = i + 1; j < rs.num_vertices; ++j) {
			if (rs.is_edge_used(i, j)) { continue; }
			/* size_t x = omap[rs.vertex_label(i)], y = omap[rs.vertex_label(j)]; */
			size_t x = i, y = j;
			double p = dist(e2);
			assert(strm[i][j] == '-');
			if (p < forbid_p) {
				p = dist(e2);
				size_t pick = 0;
				/* std::cout << x << " " << y << std::endl; */
				if (om(x,y) == 1 && om(y,x) == 1) {
					if (p < 0.5) { pick = 1; } else { pick = 2; }
				} else if (om(x,y) == 0 && om(y,x) == 0) {
					if (p < 0.5) { pick = 0; } else { pick = 2; }
				} else {
					if (p < 0.5) { pick = 0; } else { pick = 1; }
				}
				/* std::cout << i << ":" << j << " "<< om(x,y) << "," << om(y,x) << " " << pick << std::endl; */
				/* std::cout << "<<" << x << "," << y << std::endl; */
				rs.forbid(i, j, pick);
			}
		}
	}
}

double compute_diff(std::string p, Score s, double forbid_p, std::array<double, 3>& probs) {
	std::ifstream input(p);
	RelationSet rs(0);
	Tree org_tree;
	bool success;
	std::vector<std::string> om;
	std::tie(rs, org_tree, om, success) = parse::ortho::read(input);
	size_t i = 0;
	double res = 0;
	while (success) {
		i += 1;
		auto m1 = build_ortho_matrix(org_tree);
		restrict_relations(forbid_p, org_tree, om, rs);
		/* print(m1); */
		/* for (auto& str : om) { */
		/* 	std::cout << str << std::endl; */
		/* } */

		auto [t, has_tree] = build_cotree_dist(rs, probs);
		assert(has_tree);
		/* std::cout << ((has_tree)?"TRUE":"FALSE") << std::endl; */

		/* print(m1); */
		auto m2 = build_ortho_matrix(t);
		/* std::cout << m1.size1() << " " << m2.size1()<< std::endl; */
		/* print(m2); */
		/* std::cout << relative_diff(m1, m2) << std::endl; */
		if (s == Score::matrix_diff) {
			res += relative_diff(m1, m2);
		} else if (s == Score::unknown_diff) {
			res += relative_unknown_score(om, m1, m2);
		}
		std::tie(rs, org_tree, om, success) = parse::ortho::read(input);
	}
	return res/i;
}

double compute_diff(std::string p, Score s, size_t max_iter = 0) {
	std::ifstream input(p);
	RelationSet rs(0);
	Tree org_tree;
	bool success;
	std::vector<std::string> om;
	std::tie(rs, org_tree, om, success) = parse::ortho::read(input);
	size_t i = 0;
	double res = 0;
	while (success) {
		i += 1;
		/* if (num_vertices(org_tree) <= 1) { */ 
		/* 	std::tie(rs, org_tree, om, success) = parse::ortho::read(input); */
		/* 	continue; */ 
		/* } */
		auto [t, has_tree] = build_cotree(rs);
		assert(has_tree);
		/* std::cout << num_vertices(org_tree) << std::endl; */

		auto m1 = build_ortho_matrix(org_tree);
		/* print(m1); */
		auto m2 = build_ortho_matrix(t);
		/* print(m2); */
		/* std::cout << relative_diff(m1, m2) << std::endl; */
		if (s == Score::matrix_diff) {
			res += relative_diff(m1, m2);
		} else if (s == Score::unknown_diff) {
			res += relative_unknown_score(om, m1, m2);
		}
		/* double d = relative_diff(m1, m2); */
		/* double ud = relative_unknown_score(om, m1, m2); */
		/* res += d; */
		/* res_unknown += ud; */
		/* std::cout << i << " " << res << std::endl; */
		if (i == max_iter) { break; }
		std::tie(rs, org_tree, om, success) = parse::ortho::read(input);
	}
	/* std::cout << "AVG REL DIFF: " << res << " " << (i-1) << " " << (res/(i-1)) << " " << (res_unknown/(i-1)) <<  std::endl; */
	return res/i;
}

std::array<double, 3> compute_event_distribution(std::string p) {
	std::ifstream input(p);
	RelationSet rs(0);
	Tree org_tree;
	bool success;
	std::vector<std::string> om;
	std::tie(rs, org_tree, om, success) = parse::ortho::read(input);
	size_t i = 0;
	std::array<double, 3> freq{{0,0,0}};
	while (success) {
		i += 1;
		freq[0] += num_edges(rs.sym_relations()[0]);
		freq[1] += num_edges(rs.sym_relations()[1]);
		freq[2] += num_edges(rs.asym_relations()[0]);

		std::tie(rs, org_tree, om, success) = parse::ortho::read(input);
	}
	size_t total = freq[0] + freq[1] + freq[2];
	freq[0] = freq[0] / total;
	freq[1] = freq[1] / total;
	freq[2] = freq[2] / total;
	return freq;
}



int main(int argc, char *argv[]) {
	std::string u, hl, dl;
	int num_sims;
	bool pt, strip_labels,use_forbid;
	std::string use_dist;
	Score score_type = Score::matrix_diff;
	size_t max_iter;
	try {
		using namespace cxxopts;

	cxxopts::Options options("Homology", "");
	options.add_options()
		("u,unknown", "unknown probability", cxxopts::value<std::string>()->default_value("0.2"))
		("d,dl", "dupl loss", cxxopts::value<std::string>()->default_value("50"))
		("h,hl", "hgt loss", cxxopts::value<std::string>()->default_value("50"))
		("s,score", "score type [md,ud]", value<std::string>()->default_value("md"))
		("t,print_trees", "print trees")
		("l,strip", "strip inner vertex labels")
		("m,dist", "use computed distribution", value<std::string>()->default_value(""))
		("f,forbid", "use forbidden")
		("i,sims", "number of simulations", cxxopts::value<int>()->default_value("1"))
		("max-iter", "max number of genes", cxxopts::value<int>()->default_value("0"))
		;
	options.parse(argc, argv);
	u = options["u"].as<std::string>();
	hl = options["hl"].as<std::string>();
	dl = options["dl"].as<std::string>();
	num_sims = options["sims"].as<int>();
	max_iter = options["max-iter"].as<int>();
	pt = options["t"].count();
	use_forbid = options["forbid"].count();
	/* use_dist = options["m"].count(); */
	use_dist = options["dist"].as<std::string>();
	strip_labels = options["l"].count();
	if (options["score"].as<std::string>() == "md") {
		score_type = Score::matrix_diff;
	} else if (options["score"].as<std::string>() == "ud") {
		score_type = Score::unknown_diff;
	}
	} catch (cxxopts::argument_incorrect_type e) {
		std::cout << e.what() << std::endl;
		return 0;
	}

	auto dls = split(dl, ',');
	auto hls = split(hl, ',');
	if (!pt) { 
		std::cout << "prop\thl\tdl\tdiff";
		std::cout << "\torder";
		std::cout << std::endl; 
	}
	for (const auto& d : dls) {
	for (const auto& h : hls) {
	std::vector<double> res(18, 0);
	for (int i = 1; i <= num_sims; ++i) {
	std::string path =
		"../data/DirectedSimulations/data/DL"+d+"/H"+h+"/sim"+std::to_string(i)+".org.paraphylo";
	std::ifstream input(path);
	auto  [ortho_matrix, success] = parse::ortho::read_org(input);
	/* std::cout << ortho_matrix << std::endl; */
	assert(success);
	size_t res_pos = 0;
	for (double u = 0.05; u < 0.98; u += 0.05) {
		RelationSet rs(ortho_matrix, u);
		std::string ustr = std::to_string(u).substr(0, 4);
		if (use_forbid) {
			ustr = std::to_string(0.5).substr(0,4);
		}
		if (use_dist != "") {
		} else {
			auto [t, has_tree] = build_cotree(rs);
			auto pm = build_ortho_matrix(t);
			assert(has_tree);
			res[res_pos] += relative_diff(ortho_matrix, pm);
			std::cout << "res_ps: " << res[res_pos] << std::endl;
		}
		++res_pos;
	}
	}
	size_t res_pos = 0;
	for (double u = 0.05; u < 0.98; u += 0.05) {
		res[res_pos] /= num_sims;
		std::cout << std::to_string(u).substr(0,4) << "\t" << h << "\t" << d << "\t" << res[res_pos];
		if (use_dist != "") {
			std::cout << "\t" << use_dist[0] << "/" << use_dist[1] << "/" << use_dist[2];
		}
	std::cout << std::endl;
		++res_pos;
	}
	}
	}
	/* std::string path = "../data/DirectedSimulations/data/DL50/H50/sim1.u0.20.partial.paraphylo"; */
	/* read_file(path); */

	/* std::cout << m2 << std::endl; */
	/* print_tree_dot(ot); */
	/* print_tree_dot(t); */
	return 0;
}
