/**
 * @file PD_cantilever.cpp
 * @brief 3D cantilever beam with NOSB-PD (Hughes-Winget incremental
 * constitutive update), the first validation case of the peridynamics
 * module ported from Icebreaking_PD-SPH.
 * @ref 	doi.org/10.1016/j.jcp.2013.12.012
 */
#include "sphinxsys.h"
/** Name space. */
using namespace SPH;
/** Geometry parameters. */
Real PL = 6.0;
Real PH = 1.0;
Real PW = 1.0;
Real SL = 0.5;
Real global_resolution = PH / 12.0; /**< Initial particle spacing. */
Real BW = global_resolution * 4;    /**< Boundary width. */
Vecd halfsize_cantilever(0.5 * (PL + SL), 0.5 * PH, 0.5 * PW);
Vecd translation_cantilever(0.5 * (PL - SL), 0.5 * PH, 0.5 * PW);
Vecd halfsize_holder(0.5 * SL, 0.5 * PH, 0.5 * PW);
Vecd translation_holder(-0.5 * SL, 0.5 * PH, 0.5 * PW);
/** Domain bounds of the system. */
BoundingBoxd system_domain_bounds(Vecd(-SL - BW, -BW, -BW),
                                  Vecd(PL + BW, PH + BW, PH + BW));
// Observer location
StdVec<Vecd> observation_location = {Vecd(PL, PH, PW)};
/** For material properties of the solid. */
Real rho0_s = 1100.0;
Real poisson = 0.45;
Real Youngs_modulus = 1.7e7;
Real gravity_g = 0.0;

/** Define the cantilever shape. */
class Cantilever : public ComplexShape
{
  public:
    explicit Cantilever(const std::string &shape_name) : ComplexShape(shape_name)
    {
        add<GeometricShapeBox>(Transform(translation_cantilever), halfsize_cantilever);
        add<GeometricShapeBox>(Transform(translation_holder), halfsize_holder);
    }
};
/**
 * application dependent initial condition
 */
class CantileverInitialCondition
    : public solid_dynamics::ElasticDynamicsInitialCondition
{
  public:
    explicit CantileverInitialCondition(SPHBody &sph_body)
        : solid_dynamics::ElasticDynamicsInitialCondition(sph_body) {};

    void update(size_t index_i, Real dt)
    {
        if (pos_[index_i][0] > 0.0)
        {
            vel_[index_i][1] = 5.0 * sqrt(3.0);
            vel_[index_i][2] = 5.0;
        }
    };
};
/**
 *  The main program
 */
int main(int ac, char *av[])
{
    /** Setup the system. Please the make sure the global domain bounds are correctly defined. */
    SPHSystem sph_system(system_domain_bounds, global_resolution);
    sph_system.handleCommandlineOptions(ac, av);
    /** create a Cantilever PD body, corresponding material and particles. */
    PDBody cantilever_body(sph_system, makeShared<Cantilever>("PDBody"));
    cantilever_body.defineMatterMaterial<HughesWingetSolid>(rho0_s, Youngs_modulus, poisson);
    cantilever_body.generateParticles<BaseParticles, Lattice>();
    cantilever_body.getBaseParticles().addVariableToWrite<Real>("Damage");
    cantilever_body.getBaseParticles().addVariableToWrite<Matd>("PDStress");
    /** Define Observer. */
    ObserverBody cantilever_observer(sph_system, "CantileverObserver");
    cantilever_observer.generateParticles<ObserverParticles>(observation_location);

    /** topology */
    InnerRelation cantilever_body_inner(cantilever_body);
    ContactRelation cantilever_observer_contact(cantilever_observer, {&cantilever_body});

    /**
     * This section define all numerical methods will be used in this case.
     */
    SimpleDynamics<CantileverInitialCondition> initialization(cantilever_body);
    /** prior acceleration by (zero) gravity. */
    SimpleDynamics<solid_dynamics::PDTimeStepInitialization>
        initialize_a_pd_step(cantilever_body, makeShared<Gravity>(Vecd(0.0, 0.0, -gravity_g)));
    /** calculate the shape matrix once in the reference configuration. */
    InteractionWithUpdate<solid_dynamics::NosbPDShapeMatrix> cantilever_shape_matrix(cantilever_body_inner);
    /** Time step size calculation. */
    ReduceDynamics<solid_dynamics::PDTimeStep> computing_time_step_size(cantilever_body, 0.085);
    /** NOSB-PD position-Verlet steps with hourglass control and numerical damping. */
    SimpleDynamics<solid_dynamics::NosbPDFirstStep> nosbpd_first_step(cantilever_body);
    InteractionWithUpdate<solid_dynamics::NosbPDSecondStep> nosbpd_second_step(cantilever_body_inner);
    InteractionDynamics<solid_dynamics::LittleWoodHourGlassControl> hourglass_control(cantilever_body_inner);
    InteractionDynamics<solid_dynamics::PairNumericalDampingforPD> numerical_damping(cantilever_body_inner);
    InteractionDynamics<solid_dynamics::NosbPDThirdStep> nosbpd_third_step(cantilever_body_inner);
    SimpleDynamics<solid_dynamics::NosbPDFourthStep> nosbpd_fourth_step(cantilever_body);
    /** Constrain the holder. */
    GeometricShapeBox holder_shape(Transform(translation_holder), halfsize_holder, "Holder");
    BodyRegionByParticle holder(cantilever_body, holder_shape);
    SimpleDynamics<FixBodyPartConstraint> constraint_holder(holder);
    /** Output */
    BodyStatesRecordingToVtp write_states(sph_system);
    ObservedQuantityRecording<Vecd> write_displacement("Position", cantilever_observer_contact);
    /**
     * From here the time stepping begins.
     */
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    /** apply initial condition and compute the shape matrix once. */
    initialization.exec();
    cantilever_shape_matrix.exec();
    write_states.writeToFile(0);
    write_displacement.writeToFile(0);
    /** Setup physical parameters. */
    Real &physical_time = *sph_system.getSystemVariableDataByName<Real>("PhysicalTime");
    int ite = 0;
    Real end_time = 2.0;
    Real output_period = end_time / 200.0;
    Real dt = 0.0;
    /** Statistics for computing time. */
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    /**
     * Main loop.
     * NOTE: the peridynamic inner configuration is built once above and
     * is intentionally NOT updated during the simulation, so that the
     * reference bond geometry and bond states persist.
     */
    while (physical_time < end_time)
    {
        Real integration_time = 0.0;
        while (integration_time < output_period && physical_time < end_time)
        {
            if (ite % 100 == 0)
            {
                std::cout << "N=" << ite << " Time: "
                          << physical_time << "	dt: "
                          << dt << "\n";
            }
            dt = peridynamics::boundedTimeStep(computing_time_step_size.exec(), {output_period - integration_time, end_time - physical_time});
            initialize_a_pd_step.exec(dt);
            nosbpd_first_step.exec(dt);
            nosbpd_second_step.exec(dt);
            hourglass_control.exec(dt);
            numerical_damping.exec(dt);
            nosbpd_third_step.exec(dt);
            nosbpd_fourth_step.exec(dt);
            constraint_holder.exec(dt);

            ite++;
            integration_time = dt == output_period - integration_time ? output_period : integration_time + dt;
            physical_time = dt == end_time - physical_time ? end_time : physical_time + dt;
        }
        write_displacement.writeToFile(ite);
        TickCount t2 = TickCount::now();
        write_states.writeToFile();
        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }
    TickCount t4 = TickCount::now();

    TimeInterval tt;
    tt = t4 - t1 - interval;
    std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;

    return 0;
}
