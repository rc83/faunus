#pragma once
#include "move.h"
#include "energy.h"
#include "units.h"

namespace Faunus {

namespace Move {
/** @section Force moves
 *
 * Force moves are pseudo-Monte-Carlo moves that run a short trajectory of Langevin or similar dynamics.
 * The move is always accepted regardless of the total potential energy of the final system's configuration.
 *
 * ForceMoveBase provides an interface for all dynamics moves. It orchestrates integrators, thermostats, etc.
 * to work in accord.
 * IntegratorBase provides an interface for all integrators that are responsible for development
 * of the simulation in time steps. The integrators update positions and velocities of particles as needed.
 */


/**
 * @brief Generate a random 3d vector from the normal distribution
 * @note A helper class only to be used within this module.
 */
class NormalRandomVector {
    std::normal_distribution<double> normal_distribution;

  public:
    NormalRandomVector(double mean = 0.0, double stddev = 1.0) : normal_distribution(mean, stddev){};
    inline Point operator()() {
        return {normal_distribution(random.engine), normal_distribution(random.engine),
                normal_distribution(random.engine)};
    };
};


/**
 * @brief Base class for MD thermostats
 *
 * To add a new derived thermostat include it into ThermostatType as well as JSON mapping and extend makeThermostat
 * static factory method.
 */
class ThermostatBase {
  public:
    enum ThermostatType { UNDEFINED, ANDERSEN, LANGEVIN, NVE }; // CANONICALVSCALE, NOSEHOVERCHAINS
  protected:
    ThermostatType type;
    ThermostatBase(ThermostatType type);

  public:
    static std::shared_ptr<ThermostatBase> makeThermostat(ThermostatType type);
    static std::shared_ptr<ThermostatBase> makeThermostat(const json &j);

    virtual void apply(ParticleVector &particles, PointVector &velocities) = 0;
    virtual void from_json(const json &j) = 0;
    virtual void to_json(json &j) const = 0;
    virtual ~ThermostatBase() = default;
};

NLOHMANN_JSON_SERIALIZE_ENUM(ThermostatBase::ThermostatType, {
    {ThermostatBase::UNDEFINED, "undefined"},
    {ThermostatBase::ANDERSEN, "andersen"},
    {ThermostatBase::LANGEVIN, "langevin"},
    {ThermostatBase::NVE, "nve"},
})

void from_json(const json &j, ThermostatBase &i);
void to_json(json &j, const ThermostatBase &i);

// todo
struct NVEThermostat : public ThermostatBase {
    NVEThermostat() : ThermostatBase(ThermostatBase::NVE){};
    void from_json(const json &) override{};
    void to_json(json &) const override{};
    void apply(Space::Tpvec &, PointVector &) override{};
};

// todo
struct LangevinThermostat : public ThermostatBase {
    LangevinThermostat() : ThermostatBase(ThermostatBase::LANGEVIN){};
    void from_json(const json &) override{};
    void to_json(json &) const override{};
    void apply(Space::Tpvec &, PointVector &) override{};
};

class AndersenThermostat : public ThermostatBase {
    double nu; //!< collision frequency
  public:
    AndersenThermostat(double nu = 0.0);
    void apply(ParticleVector &particles, PointVector &velocities) override;
    void from_json(const json &) override;
    void to_json(json &) const override;
};

/**
 * @brief Base class for dynamics integrators
 *
 * Integrators progress the simulation in time steps. Positions and velocities of the particles as well as forces
 * acting on them are updated in every time step as prescribed by the given integrator scheme.
 */
class IntegratorBase {
  protected:
    Space &spc;
    Energy::Energybase &energy;
    IntegratorBase(Space &, Energy::Energybase &);
    virtual ~IntegratorBase() = default;

  public:
    /**
     * @brief Move particles by one time step and update their positions, velocities and acting forces
     * @param velocities  a vector of particles' velocities
     * @param forces  a vector of particles' forces
     */
    virtual void step(PointVector &velocities, PointVector &forces) = 0; // todo shall we return real time progress?
    virtual void from_json(const json &j) = 0;
    virtual void to_json(json &j) const = 0;
};

void from_json(const json &j, IntegratorBase &i);
void to_json(json &j, const IntegratorBase &i);

/**
 * @brief Symmetric Langevin velocity-Verlet method (BAOAB)
 *
 * The integrator can conduct normal velocity-Verlet, when friction_coefficient = 0.
 */
class LangevinVelocityVerlet : public IntegratorBase {
    double time_step;            //!< integration time step (picoseconds)
    double friction_coefficient; //!< friction coefficient (inverse picoseconds)
    NormalRandomVector noise;    //!< generator of random 3d vectors from the normal distribution

    //! increment particle's position in an A semi-step
    inline Point positionIncrement(const Point& velocity);
    //! increment particle's velocity in a B semi-step
    inline Point velocityIncrement(const Point& force, const double mass);
    //! apply fluctuation-dissipation to particle's velocity by solving of Ornstein-Uhlenbeck process in an O-step
    inline Point velocityFluctuationDissipation(const Point& velocity, const double mass);
  public:
    LangevinVelocityVerlet(Space &spc, Energy::Energybase &energy);
    LangevinVelocityVerlet(Space &spc, Energy::Energybase &energy, double time_step, double friction_coefficient = 0.0);
    LangevinVelocityVerlet(Space &spc, Energy::Energybase &energy, const json &j);
    void step(PointVector &velocities, PointVector &forces) override;
    void from_json(const json &j) override;
    void to_json(json &j) const override;
};

void from_json(const json &j, IntegratorBase &i);
void to_json(json &j, const IntegratorBase &i);

/**
 * @brief Base class for force moves, e.g., molecular dynamics or Langevin dynamics.
 *
 * Orchestrate execution of integrators, thermostats, etc. Store vectors for velocities and forces, which are not part
 * of the particle vector in the space.
 */
class ForceMoveBase : public Movebase {
  protected:
    unsigned int nsteps; //!< number of integration steps to perform during a single MC move
    Space &spc;          //!< space with particles and their positions
    std::shared_ptr<IntegratorBase> integrator;
    NormalRandomVector random_vector; //!< generator of random 3d vectors from the normal distribution
    PointVector velocities, forces;   //!< velocities and forces; the vector index corresponds to the particle index

    void _to_json(json &j) const override;
    void _from_json(const json &j) override;
    void _move(Change &change) override;
    void setVelocities(const PointVector &velocities); //!< set velocities
    void generateVelocities(); //!< set the initial conditions for the dynamics: generate velocities and zero forces
    ForceMoveBase(Space &, std::shared_ptr<IntegratorBase> integrator, unsigned int nsteps);
    virtual ~ForceMoveBase() = default;

  public:
    double bias(Change &, double, double) override;
};

/**
 * @brief Langevin dynamics move using Langevin equation of motion
 */
class LangevinDynamics : public ForceMoveBase {
  public:
    LangevinDynamics(Space &, std::shared_ptr<IntegratorBase> integrator, unsigned int nsteps);
    LangevinDynamics(Space &, Energy::Energybase &, const json &);
    LangevinDynamics(Space &, Energy::Energybase &);
    void _to_json(json &j) const override;
    void _from_json(const json &j) override;
};

/**
 * @brief Molecular dynamics move using Newton equation of motion
 */
//class NewtonMove : public ForceMoveBase {
//  protected:
//    std::shared_ptr<ThermostatBase> thermostat;
//
//  public:
//    NewtonMove(Space &spc, std::shared_ptr<IntegratorBase> integrator, std::shared_ptr<ThermostatBase> thermostat,
//               unsigned int nsteps = 0);
//    NewtonMove(Space &spc, Energy::Energybase &energy, const json &j);
//};

} // end of namespace Move
} // end of namespace Faunus
