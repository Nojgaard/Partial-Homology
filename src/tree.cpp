#include "tree.hpp"
#include <vector>
#include <boost/graph/subgraph.hpp>
#include <random>

namespace homology {

std::vector<SubGraph<UndirectedGraph>> build_sym_constraints(const RelationSet& rs) {
	std::vector<SubGraph<UndirectedGraph>> out;
	for (size_t i = 0; i < rs.num_sym_relations(); ++i) {
		out.emplace_back(num_vertices(rs.sym_relations()[i]));
		auto& constraint_graph = out.back();
		constraint_graph[self].label = rs.sym_relations()[i][self].label;
		for (size_t j = 0; j < rs.num_sym_relations(); ++j) {
			if (i == j) { continue; }
			merge(rs.sym_relations()[j], constraint_graph);
		}
		for (const auto& g : rs.asym_relations()) {
			merge(g, constraint_graph);
		}
		merge(rs.forbidden_relations()[i], constraint_graph);
	}
	return out;
}

std::vector<SubGraph<DirectedGraph>> build_asym_constraints(const RelationSet& rs) {
	std::vector<SubGraph<DirectedGraph>> out;
	for (size_t i = 0; i < rs.num_asym_relations(); ++i) {
		out.emplace_back(num_vertices(rs.asym_relations()[i]));
		auto& cg = out.back();
		cg[self].label = rs.asym_relations()[i][self].label;
		for (size_t j = 0; j < rs.num_asym_relations(); ++j) {
			if (j == i) { continue; }
			merge_as_bidir(rs.asym_relations()[j], cg);
		}
		merge(rs.asym_relations()[i], cg);
		for (const auto& g: rs.sym_relations()) {
			merge_as_bidir(g, cg);
		}
		merge_as_bidir(rs.forbidden_relations()[2], cg);
	}
	return out;
}

struct ConstraintSet {
	/* SubGraph<UndirectedGraph>& vs; */
	std::vector<size_t> vs;
	std::vector<SubGraph<UndirectedGraph>*> scg;
	std::vector<SubGraph<DirectedGraph>*> acg;
};

std::vector<ConstraintSet> partition(ConstraintSet& cs, size_t n, 
		const std::vector<int>& comps) {
	std::vector<ConstraintSet> out;
	for (size_t i = 0; i < n; ++i) {
		out.push_back(ConstraintSet{
				/* cs.vs.create_subgraph(), */
				std::vector<size_t>(),
				std::vector<SubGraph<UndirectedGraph>*>(),
				std::vector<SubGraph<DirectedGraph>*>()
				});
		for (SubGraph<UndirectedGraph>* g : cs.scg) {
			out.back().scg.push_back(&g->create_subgraph());
		}
		for (SubGraph<DirectedGraph>* g : cs.acg) {
			out.back().acg.push_back(&g->create_subgraph());
		}
	}
	for (size_t v = 0; v < comps.size(); ++v) {
		ConstraintSet& ocs = out[comps[v]];
		size_t global_v = cs.vs[v];
		ocs.vs.push_back(global_v);
		for (SubGraph<UndirectedGraph>* g : ocs.scg) { add_vertex(global_v, *g); }
		for (SubGraph<DirectedGraph>* g : ocs.acg) { add_vertex(global_v, *g); }
	}
	return out;
}

DirectedGraph build_quotient_graph(const SubGraph<DirectedGraph>& g, 
		size_t n, std::vector<int>& comp) {
	DirectedGraph og(n);

	for (auto e : homology::edges(g)) {
		size_t src = comp[source(e, g)], tar = comp[target(e, g)];
		if (src == tar || edge(src, tar, og).second) { continue; }
		add_edge(src, tar, og);
	}
	return og;
}

std::pair<Node, bool> build_cotree(const RelationSet& rs, ConstraintSet& cs, Tree& t){
	/* std::cout << std::endl; */
	if (cs.vs.size() == 1) {
		auto v = add_vertex(t);
		/* t[v].label = std::to_string(cs.vs[0]); */
		t[v].label = rs.vertex_label(cs.vs[0]);
		return std::make_pair(v, true);
	}

	for (const auto& g : cs.scg) {
		auto [num_comp, components] = connected_components(*g);
		if (num_comp == 1) { continue; }
		/* std::cout << g->root()[self].label << std::endl; */

		std::vector<ConstraintSet> ccs = partition(cs, num_comp, components);
		Node n = add_vertex(t);
		t[n].label = g->root()[self].label;
		for (ConstraintSet& comp_cs : ccs) {
			auto [r, success] = build_cotree(rs, comp_cs, t);
			if (!success) { return std::make_pair(0, false); }
			add_edge(n, r, t);
		}
		return std::make_pair(n, true);
	}

	for (const SubGraph<DirectedGraph>* g : cs.acg) {
		auto [num_comp, components] = strong_connected_components(*g);
		if (num_comp == 1) { continue; }
		/* std::cout << g->root()[self].label << std::endl; */
		std::vector<ConstraintSet> ccs = partition(cs, num_comp, components);
		auto qg = build_quotient_graph(*g, num_comp, components);
		auto sorted_qg = topo_sort(qg);
		Node n = add_vertex(t);
		t[n].label = g->root()[self].label;
		for (auto i = sorted_qg.rbegin(); i != sorted_qg.rend(); ++i) {
			size_t v = *i;
			ConstraintSet& comp_cs = ccs[v];
			auto [r, success] = build_cotree(rs, comp_cs, t);
			if (!success) { return std::make_pair(0, false); }
			add_edge(n, r, t);
		}
		return std::make_pair(n, true);
	}
	return std::make_pair(0, false);
}

std::pair<Tree, bool> build_cotree(const RelationSet& rs) {
	assert(rs.sym_relations().size() > 0);
	Tree cotree;
	auto asym_constraints = build_asym_constraints(rs);
	auto sym_constraints = build_sym_constraints(rs);
	std::vector<size_t> vertex_set;
	for (size_t i = 0; i < num_vertices(rs.sym_relations()[0]); ++i) {
		vertex_set.push_back(i);
	}
	ConstraintSet cs;
	cs.vs = vertex_set;
	for (auto& g : sym_constraints) { cs.scg.push_back(&g); }
	for (auto& g : asym_constraints) { cs.acg.push_back(&g); }
	auto [root, success] = build_cotree(rs, cs, cotree);
	cotree[self].root = root;
	/* std::cout << ((success) ? "TRUE" : "FALSE") << std::endl; */

	return std::make_pair(cotree, true);
}

struct ApplicableRule {
	std::string label;
	size_t num_comp;
	std::vector<int> comps;
	int type;
};

std::pair<Node, bool> build_cotree_dist(const RelationSet& rs, ConstraintSet& cs, Tree& t, const std::array<double, 3>& probs){
	/* std::cout << std::endl; */
	/* std::cout << "["; */
	/* for (auto v : cs.vs) { */
	/* 	std::cout << v << " "; */
	/* } */
	/* std::cout << "]" << std::endl; */
	if (cs.vs.size() == 1) {
		auto v = add_vertex(t);
		t[v].label = rs.vertex_label(cs.vs[0]);
		return std::make_pair(v, true);
	}

	std::vector<ApplicableRule> valid_rules;

	int type = -1;
	for (const SubGraph<UndirectedGraph>* g : cs.scg) {
		auto [num_comp, components] = connected_components(*g);
		type += 1;
		if (num_comp == 1) { continue; }
		valid_rules.push_back(ApplicableRule{g->root()[self].label, num_comp, components, type});
	}

	for (const SubGraph<DirectedGraph>* g : cs.acg) {
		auto [num_comp, components] = strong_connected_components(*g);
		type += 1;
		if (num_comp == 1) { continue; }
		valid_rules.push_back(ApplicableRule{g->root()[self].label, num_comp, components, type});
	}
	/* std::cout << "VALID_RULES: " << valid_rules.size(); */
	if (valid_rules.size() == 0) {
		assert(false);
		return std::make_pair(0, false);
	}

	size_t pick = 0;
	/* std::cout << "VALID_RULES: " << valid_rules.size() << " " << valid_rules[pick].type << std::endl; */
	/* if (valid_rules.size() > 1) { */
	/* 	for (const auto& a : valid_rules) { */
	/* 		std::cout << a.type << " "; */
	/* 	} */
	/* 	std::cout << std::endl; */
	/* } */
	if (probs[0] < 0) {
		std::random_device rd;
		std::mt19937 e2(rd());
		std::uniform_int_distribution<> dist(0, valid_rules.size()-1);
		pick = dist(e2);
	} else {
		for (size_t i = 1; i < valid_rules.size(); ++i) {
			if (probs[valid_rules[i].type] > probs[valid_rules[pick].type]) { 
				pick = i; 
			}
		}
	}
	/* std::random_device rd; */
    /* std::mt19937 e2(rd()); */
    /* std::uniform_real_distribution<> dist(0, 1); */
	/* double p = dist(e2); */
	/* if (valid_rules.size() == 3) { */
	/* 	double total_p = 0; */
	/* 	for (const auto& a : valid_rules) { */
	/* 		total_p += probs[a.type]; */
	/* 		if (p <= total_p) { break; } */
	/* 		++pick; */
	/* 	} */
	/* } else if (valid_rules.size() < 3 && valid_rules.size() > 1) { */
	/* 	int missing_type = 3 - (valid_rules[0].type + valid_rules[1].type); */
	/* 	assert(missing_type >= 0 && missing_type <= 2); */
	/* 	if (probs[valid_rules[0].type] / (1 - probs[missing_type]) <= p) { */
	/* 		pick = 0; */
	/* 	} else { */
	/* 		pick = 1; */
	/* 	} */
	/* } */
	assert(pick < 3);

	std::vector<ConstraintSet> ccs = partition(cs, valid_rules[pick].num_comp, 
			valid_rules[pick].comps);
	Node n = add_vertex(t);
	t[n].label = valid_rules[pick].label;
	/* std::cout << t[n].label << " " << valid_rules[pick].num_comp << " " << valid_rules[pick].type << " css: " << ccs.size() <<  std::endl; */
	if (valid_rules[pick].label != "L") {
		for (ConstraintSet& comp_cs : ccs) {
			auto [r, success] = build_cotree_dist(rs, comp_cs, t,probs);
			if (!success) { return std::make_pair(0, false); }
			add_edge(n, r, t);
		}
	} else {
		auto qg = build_quotient_graph(*cs.acg[0], 
				valid_rules[pick].num_comp, valid_rules[pick].comps);
		auto sorted_qg = topo_sort(qg);
		for (auto i = sorted_qg.rbegin(); i != sorted_qg.rend(); ++i) {
			size_t v = *i;
			ConstraintSet& comp_cs = ccs[v];
			auto [r, success] = build_cotree_dist(rs, comp_cs, t, probs);
			if (!success) { return std::make_pair(0, false); }
			add_edge(n, r, t);
		}
	}
	return std::make_pair(n, true);
}

std::pair<Tree, bool> build_cotree_dist(const RelationSet& rs, const std::array<double,3>& probs) {
	size_t test = 0;
	assert(rs.sym_relations().size() > 0);
	Tree cotree;
	auto asym_constraints = build_asym_constraints(rs);
	auto sym_constraints = build_sym_constraints(rs);
	std::vector<size_t> vertex_set;
	for (size_t i = 0; i < num_vertices(rs.sym_relations()[0]); ++i) {
		vertex_set.push_back(i);
	}
	ConstraintSet cs;
	cs.vs = vertex_set;
	for (auto& g : sym_constraints) { cs.scg.push_back(&g); }
	for (auto& g : asym_constraints) { cs.acg.push_back(&g); }
	auto [root, success] = build_cotree_dist(rs, cs, cotree, probs);
	cotree[self].root = root;
	/* std::cout << ((success) ? "TRUE" : "FALSE") << std::endl; */

	return std::make_pair(cotree, true);
}

void collapse_tree(Node in_n, Node out_n, const Tree& in_tree, Tree& out_tree) {
	for (Node in_c : children(in_n, in_tree)) {
		Node out_c = out_n;
		if (in_tree[in_n].label != in_tree[in_c].label) {
			out_c = add_vertex(out_tree);
			out_tree[out_c].label = in_tree[in_c].label;
			add_edge(out_n, out_c, out_tree);
		}
		collapse_tree(in_c, out_c, in_tree, out_tree);
	}
}

Tree collapse_tree(const Tree& t) {
	Tree out;
	if (num_vertices(t) == 0) { return out; }
	out[self].root = add_vertex(out);
	out[root(out)].label = t[root(t)].label;
	collapse_tree(root(t), root(out), t, out);
	return out;
}
	
} /* homology */ 
