/* Conservative, opt-in reference-surface contact for the solid gel control. */
#ifndef PERI_GEL_SURFACE_CONTACT_H
#define PERI_GEL_SURFACE_CONTACT_H

#include "pd_surface_quadrature.h"
#include "pd_checkpoint.h"
#include "pd_contact.h"
#include <fstream>
#include <iomanip>

namespace gel_compression
{
class SurfaceContact
{
  public:
    SurfaceContact(SPH::SPHBody &body, SPH::Real radius, SPH::Real height, SPH::Real support,
                   SPH::Real penalty_length, int rings, int azimuths, SPH::Real mu, SPH::Real veps,
                   SPH::Real half_width)
        : particles_(body.getBaseParticles()), count_(particles_.TotalRealParticles()),
          x_(particles_.getVariableDataByName<SPH::Vecd>("Position")),
          v_(particles_.getVariableDataByName<SPH::Vecd>("Velocity")),
          mass_(particles_.getVariableDataByName<SPH::Real>("Mass")),
          acceleration_(particles_.getVariableDataByName<SPH::Vecd>("PDAcceleration")),
          force_(particles_.getVariableDataByName<SPH::Vecd>("ContactRepulsionFromWall")),
          friction_(particles_.getVariableDataByName<SPH::Vecd>("ContactFrictionForce")),
          energy_(particles_.getVariableDataByName<SPH::Real>("AnalyticContactEnergy")),
          power_(particles_.getVariableDataByName<SPH::Real>("ContactBoundaryPower")),
          friction_power_(particles_.getVariableDataByName<SPH::Real>("ContactFrictionBoundaryPower")),
          lower_(particles_.getVariableDataByName<SPH::Vecd>("LowerPlateReaction")),
          upper_(particles_.getVariableDataByName<SPH::Vecd>("UpperPlateReaction")),
          modulus_(dynamic_cast<SPH::HughesWingetSolid &>(body.getMatterMaterial()).ContactStiffness()),
          stiffness_(modulus_ / penalty_length), mu_(mu), veps_(veps), half_width_(half_width),
          rule_(SPH::Vecd(0, height, 0), radius, rings, azimuths),
          transfer_(x_, count_, rule_, support), active_(rule_.area.size(), 0),
          last_x_(rule_.area.size(), SPH::Vecd::Zero()), last_v_(last_x_),
          last_lower_(last_x_), last_upper_(last_x_)
    {}

    void setUpper(SPH::Real face, const SPH::Vecd &velocity) { face_ = face; wall_velocity_ = velocity; }

    SPH::StdVec<SPH::solid_dynamics::RigidBoxContactObservation> observe() const
    {
        SPH::StdVec<SPH::solid_dynamics::RigidBoxContactObservation> result(2);
        for (size_t s = 0; s < transfer_.areas.size(); ++s)
        {
            const auto x = transfer_.gather(s, x_);
            if (x[1] >= 0 && x[1] <= face_) continue;
            checkPlateSupport(x);
            const auto v = transfer_.gather(s, v_);
            for (size_t wall = 0; wall < 2; ++wall)
            {
                const auto response = evaluate(s, wall, x, v);
                auto &total = result[wall];
                total.reaction -= response.force;
                total.potential_energy += response.energy;
                total.boundary_power += response.boundary_power;
                total.max_penetration = std::max(total.max_penetration, response.penetration);
                total.contacting_particles += response.penetration > 0; // Surface samples, recorded explicitly in config.
            }
        }
        return result;
    }

    void apply(const SPH::Vecd *sampling_velocity, SPH::Real sample_time)
    {
        for (size_t i = 0; i < count_; ++i)
        {
            force_[i].setZero(); friction_[i].setZero(); lower_[i].setZero(); upper_[i].setZero();
            energy_[i] = power_[i] = friction_power_[i] = 0;
        }
        for (size_t s = 0; s < transfer_.areas.size(); ++s)
        {
            last_x_[s] = transfer_.gather(s, x_);
            last_v_[s] = transfer_.gather(s, sampling_velocity);
            last_lower_[s].setZero(); last_upper_[s].setZero(); active_[s] = 0;
            if (last_x_[s][1] >= 0 && last_x_[s][1] <= face_) continue;
            checkPlateSupport(last_x_[s]);
            SPH::Real square_sum = 0;
            for (size_t j = transfer_.offsets[s]; j < transfer_.offsets[s + 1]; ++j)
                square_sum += transfer_.weights[j] * transfer_.weights[j];
            for (size_t wall = 0; wall < 2; ++wall)
            {
                const auto response = evaluate(s, wall, last_x_[s], last_v_[s]);
                if (response.penetration == 0) continue;
                active_[s] |= size_t(1) << wall;
                (wall == 0 ? last_lower_ : last_upper_)[s] = -response.force;
                transfer_.scatter(s, response.force, force_);
                transfer_.scatter(s, response.friction, friction_);
                transfer_.scatter(s, -response.force, wall == 0 ? lower_ : upper_);
                for (size_t j = transfer_.offsets[s]; j < transfer_.offsets[s + 1]; ++j)
                {
                    const auto i = transfer_.targets[j]; const auto phi = transfer_.weights[j];
                    // Nonnegative bookkeeping allocation; force/power still use the exact transpose.
                    energy_[i] += phi * phi / square_sum * response.energy;
                    power_[i] += phi * response.boundary_power;
                    friction_power_[i] += phi * response.friction_boundary_power;
                }
            }
        }
        for (size_t i = 0; i < count_; ++i) acceleration_[i] += force_[i] / mass_[i];
        last_time_ = sample_time;
    }

    /** Gershgorin bounds of the active mapped penalty and regularized friction.
     * Predict the actual first half drift without writing positions or histories. */
    SPH::Real stableStep(SPH::Real trial, SPH::Real upper_midface, SPH::Real cfl) const
    {
        SPH::StdVec<SPH::Real> stiffness(count_, 0), damping(count_, 0);
        for (size_t s = 0; s < transfer_.areas.size(); ++s)
        {
            SPH::Real y = 0, sum_abs = 0;
            for (size_t j = transfer_.offsets[s]; j < transfer_.offsets[s + 1]; ++j)
            {
                const auto i = transfer_.targets[j]; const auto phi = transfer_.weights[j];
                y += phi * (x_[i][1] + 0.5 * trial * v_[i][1]); sum_abs += std::abs(phi);
            }
            const auto penetration = std::max({SPH::Real(0), -y, y - upper_midface});
            if (penetration == 0) continue;
            const SPH::Real k = stiffness_ * transfer_.areas[s] * (1 + mu_);
            const SPH::Real d = mu_ * stiffness_ * transfer_.areas[s] * penetration / veps_;
            for (size_t j = transfer_.offsets[s]; j < transfer_.offsets[s + 1]; ++j)
            {
                const auto i = transfer_.targets[j]; const auto row = std::abs(transfer_.weights[j]) * sum_abs;
                stiffness[i] += k * row; damping[i] += d * row;
            }
        }
        SPH::Real rate = 0;
        for (size_t i = 0; i < count_; ++i)
            rate = std::max({rate, std::sqrt(stiffness[i] / mass_[i]), damping[i] / mass_[i]});
        return rate == 0 ? trial : std::min(trial, cfl / rate);
    }

    void registerCheckpoint(SPH::peridynamics::DiskCheckpoint &checkpoint)
    {
        transfer_.registerCheckpoint(checkpoint);
        checkpoint.addArray("SurfaceActiveWallMasks", active_);
        checkpoint.addArray("SurfaceLastPositions", last_x_);
        checkpoint.addArray("SurfaceLastVelocities", last_v_);
        checkpoint.addArray("SurfaceLastLowerReaction", last_lower_);
        checkpoint.addArray("SurfaceLastUpperReaction", last_upper_);
        checkpoint.addScalar("SurfaceAppliedTime", last_time_);
    }

    void referenceFile(bool resume) const
    {
        std::ostringstream content;
        content << std::setprecision(17) << "sample,x,y,z,area,particle,weight\n";
        for (size_t s = 0; s < transfer_.areas.size(); ++s)
            for (size_t j = transfer_.offsets[s]; j < transfer_.offsets[s + 1]; ++j)
                content << s << ',' << transfer_.reference[s][0] << ',' << transfer_.reference[s][1] << ','
                        << transfer_.reference[s][2] << ',' << transfer_.areas[s] << ','
                        << transfer_.targets[j] << ',' << transfer_.weights[j] << '\n';
        const std::filesystem::path path = "output/surface_reference.csv";
        if (resume)
        {
            std::ifstream stream(path);
            const std::string stored((std::istreambuf_iterator<char>(stream)), std::istreambuf_iterator<char>());
            if (!stream || stored != content.str()) throw std::runtime_error("surface reference differs from restored checkpoint");
        }
        else
        {
            if (std::filesystem::exists(path)) throw std::runtime_error("surface reference already exists");
            std::ofstream stream(path); stream << content.str(); stream.flush();
            if (!stream) throw std::runtime_error("cannot preserve reference surface cache");
        }
    }

  private:
    SPH::peridynamics::SurfacePlaneResponse evaluate(size_t s, size_t wall, const SPH::Vecd &x, const SPH::Vecd &v) const
    {
        return SPH::peridynamics::surfacePlaneResponse(x, v, SPH::Vecd(0, wall == 0 ? 1 : -1, 0),
            wall == 0 ? 0 : -face_, wall == 0 ? SPH::Vecd::Zero().eval() : wall_velocity_,
            transfer_.areas[s], stiffness_, mu_, veps_);
    }
    void checkPlateSupport(const SPH::Vecd &x) const
    {
        if (std::abs(x[0]) >= half_width_ || std::abs(x[2]) >= half_width_)
            throw std::runtime_error("surface contact exceeded finite plate support");
    }
    SPH::BaseParticles &particles_;
    size_t count_;
    SPH::Vecd *x_, *v_;
    SPH::Real *mass_;
    SPH::Vecd *acceleration_, *force_, *friction_;
    SPH::Real *energy_, *power_, *friction_power_;
    SPH::Vecd *lower_, *upper_;
    SPH::Real modulus_, stiffness_, mu_, veps_, half_width_, face_ = 0, last_time_ = 0;
    SPH::Vecd wall_velocity_ = SPH::Vecd::Zero();
    SPH::peridynamics::SphereSurfaceRule rule_;
    SPH::peridynamics::ReferenceSurfaceTransfer transfer_;
    SPH::StdVec<size_t> active_;
    SPH::StdVec<SPH::Vecd> last_x_, last_v_, last_lower_, last_upper_;
};
} // namespace gel_compression
#endif
