#include "forcemove.h"
#include "random.h"
namespace Faunus {

namespace Move {

/**
 * @brief Compute a single dimension contribution to the mean square thermal speed of a particle, i.e., compute 〈v²_x〉.
 *
 * This is equivalent to a one third of a total mean square thermal speed of a particle in three dimensions,
 * i.e., 〈v²〉 = 3 × 〈v²_x〉.
 *
 * @param mass  particle's mass in g/mol
 * @return mean square thermal speed in (Å/ps)²
 */
static inline double meanSquareSpeedComponent(T mass) {
    return (pc::kT() / mass * 1.0_kg) * ((1.0_m * 1.0_m) / (1.0_s * 1.0_s));
}

/**
 * @brief Generate standard normal random vector
 */
class NormalRandomVector {
    std::normal_distribution<double> normal_distribution;

  public:
    NormalRandomVector(double mean = 0.0, double stddev = 1.0) : normal_distribution(mean, stddev){};
    inline Point operator()() {
        return Point(normal_distribution(random.engine),
                     normal_distribution(random.engine),
                     normal_distribution(random.engine));
    };
};

NormalRandomVector random_vector; //!< standard normal random vector generator

// =============== ThermostatBase ===============

std::shared_ptr<ThermostatBase> ThermostatBase::makeThermostat(ThermostatType type) {
    std::shared_ptr<ThermostatBase> thermostat_ptr;
    switch (type) {
    case ANDERSEN:
        thermostat_ptr = std::make_shared<AndersenThermostat>();
        break;
    default:
        throw std::invalid_argument("unknown thermostat");
    }
    return thermostat_ptr;
}

std::shared_ptr<ThermostatBase> ThermostatBase::makeThermostat(const json &j) {
    if (! j.is_object() || j.size() != 1) {
        throw std::runtime_error("bad syntax for thermostat");
    }
    const json j_name = j.begin().key();
    auto thermostat_ptr = makeThermostat(j_name.get<ThermostatBase::ThermostatType>());
    thermostat_ptr->from_json(j.begin().value());
    return thermostat_ptr;
}

ThermostatBase::ThermostatBase(ThermostatType type) : type(type) {}

void from_json(const json &j, ThermostatBase &i) { i.from_json(j); }
void to_json(json &j, const ThermostatBase &i) { i.to_json(j); }

// =============== AndersenThermostat ===============

AndersenThermostat::AndersenThermostat(double nu) : ThermostatBase(ThermostatBase::ANDERSEN), nu(nu) {}

void AndersenThermostat::from_json(const json &j) {
    nu = j.at("nu").get<double>();  // fixme units
}

void AndersenThermostat::to_json(json &j) const {
    j = {{type, {{"nu", nu}}}};  // fixme units
}

void AndersenThermostat::apply(ParticleVector &particles, PointVector &velocities) {
    std::uniform_int_distribution<std::size_t> uniform_int_distribution(0, velocities.size());
    const auto collisions_cnt = uniform_int_distribution(random.engine); // number of particles to undergo collision

    std::set<std::size_t> collisions; // particle indices undergoing collision
    while (collisions.size() < collisions_cnt) {
        auto i = random.range(0, velocities.size() - 1);
        if (collisions.insert(i).second) { // index did not exist in the set before
            velocities[i] = random_vector() * sqrt(pc::kT() / Faunus::atoms[particles[i].id].mw);
        }
    }
}

// =============== IntegratorBase  ===============

IntegratorBase::IntegratorBase(Space &spc, Energy::Energybase &energy) : spc(spc), energy(energy) {}

void from_json(const json &j, IntegratorBase &i) { i.from_json(j); }
void to_json(json &j, const IntegratorBase &i) { i.to_json(j); }

// =============== LangevinVelocityVerlet ===============

LangevinVelocityVerlet::LangevinVelocityVerlet(Space &spc, Energy::Energybase &energy)
    : IntegratorBase(spc, energy) {}

LangevinVelocityVerlet::LangevinVelocityVerlet(Space &spc, Energy::Energybase &energy, double time_step,
                                               double friction_coefficient)
    : IntegratorBase(spc, energy), time_step(time_step), friction_coefficient(friction_coefficient) {}

LangevinVelocityVerlet::LangevinVelocityVerlet(Space &spc, Energy::Energybase &energy, const json &j)
    : LangevinVelocityVerlet::LangevinVelocityVerlet(spc, energy) {
    from_json(j);
}

inline Point LangevinVelocityVerlet::positionIncrement(const Point& velocity) {
    return 0.5 * time_step * velocity;
}

inline Point LangevinVelocityVerlet::velocityIncrement(const Point& force, const double mass) {
    // As forces are in kT per ångström units (a hybrid between reduced energy units and absolute units), we use
    // the mean square speed to compute acceleration from the force and as a conversion factor.
    // Dimension analysis: (ps * 1 / Å) * (Å^2 / ps^2) = Å / ps.
    return 0.5 * time_step * force * meanSquareSpeedComponent(mass);
}

inline Point LangevinVelocityVerlet::velocityFluctuationDissipation(const Point& velocity, const double mass) {
    const double prefactor = exp(-friction_coefficient * time_step); // Ornstein-Uhlenbeck process prefactor
    return (prefactor * velocity) +
                    std::sqrt((1 - std::pow(prefactor, 2)) * meanSquareSpeedComponent(mass)) * noise();
}

void LangevinVelocityVerlet::step(PointVector &velocities, PointVector &forces) {
    assert(spc.p.size() == forces.size());     // Check that the position and force vectors are of the same size.
    assert(spc.p.size() == velocities.size()); // Check that the position and velocity vectors are of the same size.

    for (std::size_t i = 0; i < spc.p.size(); ++i) {
        const auto mass = atoms[spc.p[i].id].mw;
        auto &position = spc.p[i].pos;
        auto &velocity = velocities[i];
        velocity += velocityIncrement(forces[i], mass);            // B step
        position += positionIncrement(velocity);                   // A step
        velocity = velocityFluctuationDissipation(velocity, mass); // O step
        position += positionIncrement(velocity);                   // A step
        spc.geo.boundary(position);
    }
    // Update forces. The force vector must be initialized with zeros as the resulting force is computed
    // additively (+=).
    std::fill(forces.begin(), forces.end(), Point(0.0, 0.0, 0.0));
    energy.force(forces);

    for (std::size_t i = 0; i < spc.p.size(); ++i) {
        const auto mass = atoms[spc.p[i].id].mw;
        velocities[i] += velocityIncrement(forces[i], mass);       // B step
    }
}

void LangevinVelocityVerlet::from_json(const json &j) {
    time_step = j.at("time_step").get<double>() * 1.0_ps;
    friction_coefficient = j.at("friction").get<double>() / 1.0_ps;
}

void LangevinVelocityVerlet::to_json(json &j) const {
    j = {{"time_step", time_step / 1.0_ps}, {"friction", friction_coefficient * 1.0_ps}};
}

// =============== ForceMoveBase ===============

ForceMoveBase::ForceMoveBase(Space &spc, std::shared_ptr<IntegratorBase> integrator, unsigned int nsteps)
    : nsteps(nsteps), spc(spc), integrator(integrator) {
    forces.resize(spc.p.size());
    velocities.resize(spc.p.size());
    repeat = 1;
}

void ForceMoveBase::_move(Change &change) {
    change.clear();
    change.all = true;
    for (unsigned int step = 0; step < nsteps; ++step) {
        integrator->step(velocities, forces);
    }
}

void ForceMoveBase::_to_json(json &j) const {
    j = {{"nsteps", nsteps}};
    j["integrator"] = *integrator;
}

void ForceMoveBase::_from_json(const json &j) {
    nsteps = j.at("nsteps").get<unsigned int>();
    integrator->from_json(j["integrator"]);
    generateVelocities();
}

double ForceMoveBase::bias(Change &, double, double) {
    return pc::neg_infty; // always accept the move
}

void ForceMoveBase::setVelocities(const PointVector &velocities) {
    assert(spc.p.size() == velocities.size()); // Check that the position and velocity vectors are of the same size.
    this->velocities = velocities;
}

void ForceMoveBase::generateVelocities() {
    for (std::size_t i = 0; i < spc.p.size(); ++i) {
        const auto mass = atoms[spc.p[i].id].mw;
        velocities[i] = random_vector() * std::sqrt(meanSquareSpeedComponent(mass));
        forces[i] = {0, 0, 0};
    }
}

// =============== LangevinMove ===============

LangevinDynamics::LangevinDynamics(Space &spc, std::shared_ptr<IntegratorBase> integrator, unsigned int nsteps)
    : ForceMoveBase::ForceMoveBase(spc, integrator, nsteps) {
    name = "langevin_dynamics";
}

LangevinDynamics::LangevinDynamics(Space &spc, Energy::Energybase &energy)
    : LangevinDynamics::LangevinDynamics(spc, std::make_shared<LangevinVelocityVerlet>(spc, energy), 0) {}

LangevinDynamics::LangevinDynamics(Space &spc, Energy::Energybase &energy, const json &j)
    : LangevinDynamics::LangevinDynamics(spc, energy) {
    from_json(j);
}

void LangevinDynamics::_to_json(json &j) const {
    ForceMoveBase::_to_json(j);

}
void LangevinDynamics::_from_json(const json &j) {
    ForceMoveBase::_from_json(j);
}

} // end of namespace Move
} // end of namespace Faunus
