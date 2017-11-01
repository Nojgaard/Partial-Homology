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
	std::cout << d << "/" << h << std::endl;
	for (int i = 1; i <= num_sims; ++i) {
	/* std::string path = */
	/* 	"../data/DirectedSimulations/data/DL"+d+"/H"+h+"/sim"+std::to_string(i)+".org.paraphylo"; */
		std::string ns = "10";
	std::string path =
		"../data/yule_trees/sp"+ns+"_D"+d+"_H"+h+"/org.paraphylo";
	/* for (double u = 0.05; u < 0.98; u += 0.05) { */
	for (double u = 0.0; u < 1.0001; u += 0.1) {
	std::string ustr = std::to_string(u).substr(0, 4);
	std::string out_path =
		"../data/yule_trees/sp"+ns+"_D"+d+"_H"+h+"/u"+ustr+".partial.paraphylo";
	std::cout << out_path << std::endl;
			/* "../data/DirectedSimulations/data/DL"+d+"/H"+h+"/sim"+std::to_string(i)+".u"+ustr+".partial.paraphylo"; */
	std::ifstream input(path);
	std::ofstream out(out_path);
	parse::ortho::write_unknown(input, out, u);
	/* std::cout << ortho_matrix << std::endl; */
	}
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
