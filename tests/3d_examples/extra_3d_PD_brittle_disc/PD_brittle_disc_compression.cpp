/**
 * @file PD_brittle_disc_compression.cpp
 * @brief Diametral compression of a brittle solid disc (Brazilian-disc type
 * splitting): a generic benchmark of the NOSB-PD module for brittle fracture.
 * @details A solid elastic disc (HughesWingetSolid) is compressed between two
 * rigid plates with smooth-Coulomb contact friction until a tensile split
 * initiates at the disc center and propagates along the loading diameter.
 * Demonstrates, in a pure-PD setup (no fluid):
 *   - PDBody with the elastic Hughes-Winget material route (NosbPDSecondStep);
 *   - bond breaking by combined principal stresses (BondBreakBySigma1andSigma3);
 *   - Littlewood hourglass control and pair-wise numerical damping;
 *   - PD-wall penalty contact with friction (ContactRepulsionFromWallforPD)
 *     and reaction recording on the driven plate.
 * The case is resolution-tagged: refine resolution_s for production runs.
 */
#include "sphinxsys.h"
using namespace SPH;

//----------------------------------------------------------------------
//	Geometry and resolution.
//----------------------------------------------------------------------
Real resolution_s = 5.0e-5;          ///< PD particle spacing (coarse: 5.0e-5 for a quick smoke run) [m]
Real disc_radius = 0.5e-3;           ///< disc radius [m]
Real disc_half_thickness = 0.125e-3; ///< disc half thickness (axis along z) [m]
Real plate_radius = 1.0e-3;          ///< half-size of the square rigid plates [m]
Real plate_thickness = 0.25e-3;      ///< rigid plate thickness [m]
Real initial_gap = 0.5 * resolution_s; ///< initial plate-disc gap (small: contact engages at once) [m]
//----------------------------------------------------------------------
//	Layout along the compression (y) axis.
//----------------------------------------------------------------------
Real lower_plate_bottom = 0.0;
Real lower_plate_top = lower_plate_bottom + plate_thickness;
Real disc_center_y = lower_plate_top + initial_gap + disc_radius;
Real upper_plate_bottom = disc_center_y + disc_radius + initial_gap;
Real upper_plate_top = upper_plate_bottom + plate_thickness;
//----------------------------------------------------------------------
//	Generic brittle solid (rock-like): elastic up to fracture.
//----------------------------------------------------------------------
Real rho0_s = 2500.0;              ///< density [kg/m^3]
Real youngs_s = 5.0e9;             ///< Young's modulus [Pa]
Real poisson_s = 0.25;             ///< Poisson ratio
Real gravity_g = 9.8;              ///< gravity (pins the disc on the lower plate;
                                   ///< PD free-surface drift otherwise lets it surf ahead of the punch) [m/s^2]
Real max_tension_stress = 2.0e6;   ///< bond-breaking tensile strength sigma1 [Pa]
Real max_shear_stress = 6.0e6;     ///< bond-breaking shear strength |sigma1 - sigma3|/2 [Pa]
//----------------------------------------------------------------------
//	Plate contact friction (smooth Coulomb; 0 = frictionless).
//----------------------------------------------------------------------
Real plate_friction_mu = 0.3;      ///< friction coefficient [-]
Real plate_friction_veps = 5.0e-3; ///< regularization velocity [m/s]
//----------------------------------------------------------------------
//	Compression schedule.
//----------------------------------------------------------------------
Real compression_speed = 0.2;                           ///< plate speed (demo rate; fracture shows earlier) [m/s]
Real compression_distance = 0.3e-3;                     ///< total stroke: ~0.2 mm is absorbed by lattice-quantized
                                                        ///< plate-disc gaps, leaving ~10% real compression [m]
Real end_time = compression_distance / compression_speed; ///< 0.0015 s
Real BW = resolution_s * 4.0;                           ///< boundary width

//----------------------------------------------------------------------
//	Disc shape: cylinder with axis along z.
//----------------------------------------------------------------------
class DiscShape : public ComplexShape
{
  public:
    explicit DiscShape(const std::string &shape_name) : ComplexShape(shape_name)
    {
        add<TriangleMeshShapeCylinder>(Vec3d(0.0, 0.0, 1.0), disc_radius,
                                       disc_half_thickness, 20, Vecd(0.0, disc_center_y, 0.0));
    }
};
//----------------------------------------------------------------------
//	Rigid plate: brick centered at the plate mid-plane.
//----------------------------------------------------------------------
class PlateShape : public ComplexShape
{
  public:
    PlateShape(const std::string &shape_name, Real mid_plane_y) : ComplexShape(shape_name)
    {
        Vecd halfsize(plate_radius, 0.5 * plate_thickness, plate_radius);
        add<GeometricShapeBox>(Transform(Vecd(0.0, mid_plane_y, 0.0)), halfsize);
    }
};

// the main program with commandline options
int main(int ac, char *av[])
{
    //----------------------------------------------------------------------
    //	Build up an SPHSystem.
    //----------------------------------------------------------------------
    BoundingBoxd system_domain_bounds(
        Vecd(-plate_radius - BW, -BW, -plate_radius - BW),
        Vecd(plate_radius + BW, upper_plate_top + BW, plate_radius + BW));
    SPHSystem sph_system(system_domain_bounds, resolution_s);
    sph_system.handleCommandlineOptions(ac, av);
    //----------------------------------------------------------------------
    //	Creating bodies with corresponding materials and particles.
    //----------------------------------------------------------------------
    PDBody disc(sph_system, makeShared<DiscShape>("BrittleDisc"));
    disc.defineMatterMaterial<HughesWingetSolid>(rho0_s, youngs_s, poisson_s);
    disc.generateParticles<BaseParticles, Lattice>();
    disc.getBaseParticles().addVariableToWrite<Real>("Damage");
    std::cout << "disc generated" << std::endl;

    SolidBody lower_plate(sph_system, makeShared<PlateShape>("LowerPlate", 0.5 * (lower_plate_bottom + lower_plate_top)));
    /** plates at the disc resolution: coarsened walls leave the contact surface
     * sparsely sampled (the disc then "surfs" ahead of the plate without loading). */
    lower_plate.defineMatterMaterial<Solid>(rho0_s, youngs_s);
    lower_plate.generateParticles<BaseParticles, Lattice>();

    SolidBody upper_plate(sph_system, makeShared<PlateShape>("UpperPlate", upper_plate_bottom + 0.5 * plate_thickness));
    upper_plate.defineMatterMaterial<Solid>(rho0_s, youngs_s);
    upper_plate.generateParticles<BaseParticles, Lattice>();
    std::cout << "plates generated" << std::endl;
    //----------------------------------------------------------------------
    //	Define body relation map.
    //----------------------------------------------------------------------
    InnerRelation disc_inner(disc);
    ContactRelation disc_plates_contact(disc, {&lower_plate, &upper_plate});
    ContactRelation lower_plate_disc_contact(lower_plate, {&disc});
    ContactRelation upper_plate_disc_contact(upper_plate, {&disc});
    std::cout << "relations done" << std::endl;
    //----------------------------------------------------------------------
    //	Algorithms of NOSB-PD dynamics (elastic / brittle route).
    //----------------------------------------------------------------------
    /** gravity baseline: holds the disc against the lower plate (disc weight
     * ~5 uN is negligible against the mN-scale contact forces). */
    SimpleDynamics<solid_dynamics::PDTimeStepInitialization> initialize_a_pd_step(disc, makeShared<Gravity>(Vecd(0.0, -gravity_g, 0.0)));
    InteractionWithUpdate<solid_dynamics::NosbPDShapeMatrix> disc_shape_matrix(disc_inner);
    SimpleDynamics<solid_dynamics::NosbPDFirstStep> nosbpd_first_step(disc);
    InteractionDynamics<solid_dynamics::BondBreakBySigma1andSigma3> check_bond_live(disc_inner, max_tension_stress, max_shear_stress);
    InteractionWithUpdate<solid_dynamics::NosbPDSecondStep> nosbpd_second_step(disc_inner);
    InteractionDynamics<solid_dynamics::LittleWoodHourGlassControl> hourglass_control(disc_inner);
    InteractionDynamics<solid_dynamics::PairNumericalDampingforPD> numerical_damping(disc_inner);
    InteractionDynamics<solid_dynamics::NosbPDThirdStep> nosbpd_third_step(disc_inner);
    SimpleDynamics<solid_dynamics::NosbPDFourthStep> nosbpd_fourth_step(disc);
    std::cout << "pd methods done" << std::endl;
    //----------------------------------------------------------------------
    //	Plate contact with friction and reaction recording.
    //----------------------------------------------------------------------
    /** the repulsion factor must be constructed (and executed) before the
     * contact force classes reading "RepulsionFactor". */
    InteractionDynamics<solid_dynamics::ContactRepulsionFactorforPD> disc_contact_factor(disc_plates_contact);
    InteractionDynamics<solid_dynamics::ContactRepulsionFromWallforPD> plate_contact_on_disc(disc_plates_contact, plate_friction_mu, plate_friction_veps);
    InteractionDynamics<solid_dynamics::ContactRepulsionToWallforPD> disc_reaction_on_lower(lower_plate_disc_contact);
    InteractionDynamics<solid_dynamics::ContactRepulsionToWallforPD> disc_reaction_on_upper(upper_plate_disc_contact);
    //----------------------------------------------------------------------
    //	Plate kinematics: upper plate driven, lower plate fixed.
    //----------------------------------------------------------------------
    Vecd upper_plate_center_end(0.0, upper_plate_bottom + 0.5 * plate_thickness - compression_distance, 0.0);
    SimpleDynamics<solid_dynamics::PositionSolidBody> drive_upper_plate(upper_plate, 0.0, end_time, upper_plate_center_end);
    //----------------------------------------------------------------------
    //	Output.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp write_body_states(sph_system);
    ReducedQuantityRecording<QuantitySummation<Vecd>> write_upper_plate_reaction(upper_plate, "ContactRepulsionToWall");
    //----------------------------------------------------------------------
    //	Prepare the simulation.
    //----------------------------------------------------------------------
    std::cout << "particle generation done\n"
              << " disc: " << disc.getBaseParticles().TotalRealParticles()
              << " lower: " << lower_plate.getBaseParticles().TotalRealParticles()
              << " upper: " << upper_plate.getBaseParticles().TotalRealParticles() << std::endl;
    sph_system.initializeSystemCellLinkedLists();
    std::cout << "cell linked lists done" << std::endl;
    sph_system.initializeSystemConfigurations();
    std::cout << "configurations done" << std::endl;
    disc_shape_matrix.exec();
    std::cout << "shape matrix done" << std::endl;
    //----------------------------------------------------------------------
    //	Setup for time-stepping control.
    //----------------------------------------------------------------------
    Real &physical_time = *sph_system.getSystemVariableDataByName<Real>("PhysicalTime");
    int number_of_iterations = 0;
    int screen_output_interval = 10000;
    Real output_interval = end_time / 60.0;
    /** fixed PD step derived from the wave-speed limit with a conservative
     * margin (ice-jet practice): at E = 5 GPa the Littlewood hourglass term
     * is stiff, and a margin of ~0.1 x dx/c is required for a clean start. */
    Real dt_s = 0.085 * resolution_s / sqrt(youngs_s / rho0_s);
    //----------------------------------------------------------------------
    //	First output before the main loop.
    //----------------------------------------------------------------------
    write_body_states.writeToFile(0);
    //----------------------------------------------------------------------
    //	Main time-stepping loop (pure PD, fixed step).
    //----------------------------------------------------------------------
    while (physical_time < end_time)
    {
        Real integration_time = 0.0;
        while (integration_time < output_interval)
        {
            /** external loads (plate contact) accumulate onto PDAccelerationPrior. */
            initialize_a_pd_step.exec(dt_s);
            disc_contact_factor.exec(dt_s);
            plate_contact_on_disc.exec(dt_s);
            /** NOSB-PD integration chain. */
            nosbpd_first_step.exec(dt_s);
            /** bond breaking checked sparsely (eigen solve per bond is expensive). */
            if (number_of_iterations % 10 == 0) check_bond_live.exec(dt_s);
            nosbpd_second_step.exec(dt_s);
            hourglass_control.exec(dt_s);
            numerical_damping.exec(dt_s);
            nosbpd_third_step.exec(dt_s);
            nosbpd_fourth_step.exec(dt_s);
            /** plate kinematics (dt-dependent) and contact reactions. */
            drive_upper_plate.exec(dt_s);
            disc_reaction_on_lower.exec(dt_s);
            disc_reaction_on_upper.exec(dt_s);

            physical_time += dt_s;
            integration_time += dt_s;
            number_of_iterations++;

            if (number_of_iterations % screen_output_interval == 0)
            {
                std::cout << std::scientific << "N=" << number_of_iterations
                          << " Time: " << physical_time << "	dt_s = " << dt_s << "\n";
            }
            if (number_of_iterations % 25 == 0)
            {
                disc.updateCellLinkedList();
                upper_plate.updateCellLinkedList();
                disc_plates_contact.updateConfiguration();
                lower_plate_disc_contact.updateConfiguration();
                upper_plate_disc_contact.updateConfiguration();
            }
        }
        write_upper_plate_reaction.writeToFile(number_of_iterations);
        write_body_states.writeToFile();
    }

    std::cout << "Total wall time for computation done." << std::endl;
    return 0;
}
