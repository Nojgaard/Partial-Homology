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

struct MetaData {
	std::array<double, 3> avg_dist;
	size_t avg_leaves;
	size_t min_leaves = 100000;
	size_t max_leaves = 0;
};

MetaData compute_meta_data(std::string p) {
	std::ifstream input(p);
	RelationSet rs(0);
	Tree org_tree;
	bool success;
	std::vector<std::string> om;
	std::tie(rs, org_tree, om, success) = parse::ortho::read(input);
	size_t i = 0;
	std::array<double, 3> freq{{0,0,0}};
	size_t total_verts = 0;
	MetaData md;
	md.avg_dist[0] = 0;
	md.avg_dist[1] = 0;
	md.avg_dist[2] = 0;
	while (success) {
		i += 1;
		/* freq[0] = num_edges(rs.sym_relations()[0]); */
		/* freq[1] = num_edges(rs.sym_relations()[1]); */
		/* freq[2] = num_edges(rs.asym_relations()[0]); */
		freq[0] = 0;
		freq[1] = 0;
		freq[2] = 0;
		for (auto v : vertices(org_tree)) {
			if (org_tree[v].label == "S") {
				++freq[0];
			} else if (org_tree[v].label == "D") {
				++freq[1];
			} else if (org_tree[v].label == "L") {
				++freq[2];
			}
		}
		size_t total = freq[0] + freq[1] + freq[2];
		/* if (total > 0) { */
		/* md.avg_dist[0] += freq[0] / total; */
		/* md.avg_dist[1] += freq[1] / total; */
		/* md.avg_dist[2] += freq[2] / total; */
		md.avg_dist[0] += freq[0];
		md.avg_dist[1] += freq[1];
		md.avg_dist[2] += freq[2];
		/* } */
		md.min_leaves = std::min(md.min_leaves, rs.num_vertices);
		md.max_leaves = std::max(md.max_leaves, rs.num_vertices);
		//std::cout << md.avg_dist[0] << std::endl;
		total_verts += rs.num_vertices;

		std::tie(rs, org_tree, om, success) = parse::ortho::read(input);
	}
		size_t total = md.avg_dist[0] + md.avg_dist[1] + md.avg_dist[2];
	md.avg_leaves = total_verts / i;
	md.avg_dist[0] /= total;
	md.avg_dist[1] /= total;
	md.avg_dist[2] /= total;
	/* md.avg_dist[0] /= i; */
	/* md.avg_dist[1] /= i; */
	/* md.avg_dist[2] /= i; */
	return md;
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
		std::cout << "dl\thl\tnverts\tsdist\tddist\tldist\tminl\tmaxl";
		std::cout << std::endl; 
	}
	for (const auto& d : dls) {
	for (const auto& h : hls) {
		MetaData meta_data;
		meta_data.avg_leaves = 0;
		meta_data.avg_dist[0] = 0;
		meta_data.avg_dist[1] = 0;
		meta_data.avg_dist[2] = 0;
	
	for (int i = 1; i <= num_sims; ++i) {
	std::string path =
		"../data/DirectedSimulations/data/DL"+d+"/H"+h+"/sim"+std::to_string(i)+".org.paraphylo";
		auto tmp_meta_data = compute_meta_data(path);
		meta_data.avg_leaves += tmp_meta_data.avg_leaves;
		meta_data.avg_dist[0] += tmp_meta_data.avg_dist[0];
		meta_data.avg_dist[1] += tmp_meta_data.avg_dist[1];
		meta_data.avg_dist[2] += tmp_meta_data.avg_dist[2];
		meta_data.min_leaves = std::min(meta_data.min_leaves, tmp_meta_data.min_leaves);;
		meta_data.max_leaves = std::max(meta_data.max_leaves, tmp_meta_data.max_leaves);;
	}
		meta_data.avg_leaves /= num_sims;
		meta_data.avg_dist[0] /= num_sims;
		meta_data.avg_dist[1] /= num_sims;
		meta_data.avg_dist[2] /= num_sims;
		std::cout << d << "\t";
		std::cout << h << "\t";
		std::cout << meta_data.avg_leaves << "\t";
		std::cout << meta_data.avg_dist[0] << "\t";
		std::cout << meta_data.avg_dist[1] << "\t";
		std::cout << meta_data.avg_dist[2] << "\t";
		std::cout << meta_data.min_leaves << "\t";
		std::cout << meta_data.max_leaves << "\t";
		std::cout << std::endl;

	}
	}
	return 0;
}
