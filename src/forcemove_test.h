#include "forcemove.h"
#include "units.h"
#include "externalpotential.h"

namespace Faunus {
namespace Move {

using doctest::Approx;

TEST_SUITE_BEGIN("ForceMove");

TEST_CASE("[Faunus] Thermostat") {

    std::shared_ptr<ThermostatBase> andersen_ptr;

    CHECK_THROWS_WITH(ThermostatBase::makeThermostat(R"({"name": "andersen", "nu": 1.0})"_json), R"(bad syntax for thermostat)");
    CHECK_THROWS_WITH(ThermostatBase::makeThermostat(R"({"xxx": {}})"_json), R"(unknown thermostat)");
    CHECK_THROWS_WITH(ThermostatBase::makeThermostat(R"({"andersen": {}})"_json), R"([json.exception.out_of_range.403] key 'nu' not found)");

    REQUIRE_NOTHROW(andersen_ptr = ThermostatBase::makeThermostat(R"({"andersen": {"nu": 2.0}})"_json));

    json j_andersen = *andersen_ptr;
    CHECK_EQ(j_andersen["andersen"]["nu"], 2.0);
}

TEST_CASE("[Faunus] Integerator") {
    class DummyEnergy : public Energy::Energybase {
        double energy(Change &) override {return 0.0; }
    };
    Space spc;
    DummyEnergy energy;

    auto ld = LangevinVelocityVerlet(spc, energy, R"({"friction": 8.0, "time_step": 2.0 })"_json);
    json j_ld = ld;
    CHECK_EQ(j_ld["friction_coefficient"], 8.0);
    CHECK_EQ(j_ld["time_step"], 2.0);
}

TEST_CASE("[Faunus] LangevinDynamics") {
    class DummyEnergy : public Energy::Energybase {
        double energy(Change &) override {return 0.0; }
    };
    Space spc;
    DummyEnergy energy;

    SUBCASE("[Faunus] JSON init") {
      json j_in =  R"({"nsteps": 100, "integrator": {"time_step": 0.001, "friction": 2.0}})"_json;
      LangevinDynamics ld(spc, energy);
      ld.from_json(j_in);
      json j_out = ld;

      //CHECK_EQ(j_out, j_in);
    }
}
TEST_SUITE_END();
} // namespace Move
} // namespace Faunus
