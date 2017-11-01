#include "../relation_set.hpp"
#include <catch/catch.hpp>

using namespace homology;

TEST_CASE("RelationSet", "[RelationSet]") {
	RelationSet rl(5);
	rl.add_sym_relation("h");
	auto id = rl.add_sym_relation("t").first;
	auto aid = rl.add_asym_relation("s").first;
	REQUIRE(rl.num_relations() == 3);
	REQUIRE(rl.num_sym_relations() == 2);
	REQUIRE(rl.num_asym_relations() == 1);
	
	bool success = rl.assign(0, 3, id);
	REQUIRE(success == true);
	success = rl.assign(3, 0, aid);
	REQUIRE(success == false);
}
