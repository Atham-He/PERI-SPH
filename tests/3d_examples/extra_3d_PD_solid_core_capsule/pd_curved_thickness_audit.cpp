#define main capsule_main_unused
#include "pd_solid_capsule.cpp"
#undef main

/** One manufactured spherical thickness mode. This does not evolve a
 * compression history and does not invoke or change material parameters. */
int main(int argc,char**argv)
{
    try {
        if(argc!=2)throw std::invalid_argument("Output JSON required");
        const Real rc=500e-6,t=10e-6,dp=25e-6,mid=rc+t/2,a=5000;
        SPHSystem system(BoundingBoxd(Vecd::Constant(-.002),Vecd::Constant(.002)),dp,4);
        PDBody body(system,makeShared<DefaultShape>("RadialThicknessPatch"));
        body.defineAdaptation<PDAdaptation>(1.5075,1);
        body.defineMatterMaterial<HughesWingetSolid>(1200,2.151e9,.3);
        body.generateParticles<BaseParticles,LayeredPDCapsuleGenerator>(rc,t,dp,3);
        InnerRelation inner(body);SimpleDynamics<solid_dynamics::PDTimeStepInitialization> init(body,makeShared<Gravity>(Vecd::Zero()));
        SimpleDynamics<solid_dynamics::NosbPDFirstStep> first(body);
        InteractionWithUpdate<solid_dynamics::NosbPDShapeMatrix> shape(inner);
        size_t ncol=body.getBaseParticles().TotalRealParticles()/3;auto*xref=body.getBaseParticles().getVariableDataByName<Vecd>("Position");
        StdVec<StdVec<size_t>> cols(ncol);StdVec<Vecd>normal(ncol);for(size_t c=0;c<ncol;++c){for(size_t l=0;l<3;++l)cols[c].push_back(l*ncol+c);normal[c]=xref[c].normalized();}
        ThicknessQuadrature rule(3);peridynamics::LayeredParticleGeometry geometry(body.getBaseParticles(),cols,normal,rule.xi);
        solid_dynamics::PDLayeredQuadraticReference reference(inner,geometry);
        system.initializeSystemCellLinkedLists();system.initializeSystemConfigurations();shape.exec();reference.build();State s(body);
        peridynamics::SolidForcePort port(body,peridynamics::SolidForcePort::Storage::PDAccelerationPrior);
        StdVec<Vecd>inner_points,outer_points;for(size_t c=0;c<ncol;++c){Vecd n=s.x0[c].normalized();inner_points.push_back(rc*n);outer_points.push_back((rc+t)*n);}
        peridynamics::LayeredSolidSurface inner_trace(port,inner_points,geometry,-1,dp),outer_trace(port,outer_points,geometry,1,dp);
        peridynamics::ReferenceSolidSurface inner_volume(port,inner_points,dp,2.5,2),outer_volume(port,outer_points,dp,2.5,2);

        for(size_t i=0;i<s.n;++i){Real r=s.x0[i].norm();Vecd n=s.x0[i]/r;s.x[i]=s.x0[i]+a*(r-mid)*(r-mid)*n;}
        Real maximum_normal_error=0,maximum_F_error=0,maximum_reference_error=0;
        std::array<Real,3> normal_error={0,0,0},computed_normal_min={1e10,1e10,1e10},computed_normal_max={-1e10,-1e10,-1e10};
        const size_t columns=s.n/3;
        for(size_t i=0;i<s.n;++i)
        {
            Real r=s.x0[i].norm(),u=a*(r-mid)*(r-mid),normal=1+2*a*(r-mid);Vecd n=s.x0[i]/r;
            Matd exact=normal*n*n.transpose()+(1+u/r)*(Matd::Identity()-n*n.transpose());
            Matd F=Matd::Zero(),Fref=Matd::Zero();
            for(auto&e:reference.rows[i]){F+=(s.x[e.j]-s.x[i])*e.gradient.transpose();Fref+=(s.x0[e.j]-s.x0[i])*e.gradient.transpose();}
            size_t l=i/columns;Real value=n.dot(F*n),err=std::abs(value-normal);
            normal_error[l]=std::max(normal_error[l],err);computed_normal_min[l]=std::min(computed_normal_min[l],value);computed_normal_max[l]=std::max(computed_normal_max[l],value);
            maximum_normal_error=std::max(maximum_normal_error,err);maximum_F_error=std::max(maximum_F_error,(F-exact).norm());maximum_reference_error=std::max(maximum_reference_error,(Fref-Matd::Identity()).norm());
        }
        Real layered_surface_error=0,volume_surface_error=0;
        for(size_t c=0;c<ncol;++c){Vecd n=s.x0[c].normalized(),inner_exact=inner_points[c]+a*t*t/4*n,outer_exact=outer_points[c]+a*t*t/4*n;
            layered_surface_error=std::max({layered_surface_error,(inner_trace.position(c)-inner_exact).norm(),(outer_trace.position(c)-outer_exact).norm()});
            volume_surface_error=std::max({volume_surface_error,(inner_volume.position(c)-inner_exact).norm(),(outer_volume.position(c)-outer_exact).norm()});}
        if(layered_surface_error>1e-12)throw std::runtime_error("Curved normal polynomial surface trace failed");
        if(maximum_normal_error>1e-8||maximum_F_error>1e-8||maximum_reference_error>1e-9)throw std::runtime_error("Curved thickness polynomial gradient failed");
        std::ofstream f(argv[1]);f<<std::setprecision(17)<<"{\"geometry\":\"actual spherical three-Gauss-layer capsule shell\",\"operator\":\"layered_normal_Cartesian_P2\",\"particles\":"<<s.n
            <<",\"thickness_m\":"<<t<<",\"spacing_m\":"<<dp<<",\"mode\":\"u_r=5000*(r-r_mid)^2\",\"exact_radial_strains\":["<<-a*t*std::sqrt(.6)<<",0,"<<a*t*std::sqrt(.6)<<']'
            <<",\"maximum_layered_surface_position_error_m\":"<<layered_surface_error<<",\"maximum_volume_surface_position_error_m\":"<<volume_surface_error<<",\"maximum_radial_gradient_error\":"<<maximum_normal_error<<",\"maximum_F_error\":"<<maximum_F_error<<",\"initial_affine_error\":"<<maximum_reference_error<<",\"layers\":[";
        for(int l=0;l<3;++l){if(l)f<<',';f<<"{\"layer\":"<<l<<",\"maximum_normal_error\":"<<normal_error[l]<<",\"computed_normal_min\":"<<computed_normal_min[l]<<",\"computed_normal_max\":"<<computed_normal_max[l]<<'}';}
        f<<"],\"new_compression_histories\":0,\"material_parameters_changed\":false}\n";
        std::cout<<"SPHERICAL_THICKNESS_MODE error="<<maximum_normal_error<<" reference="<<maximum_reference_error<<std::endl;
        return 0;
    }catch(const std::exception&e){std::cerr<<e.what()<<std::endl;return 1;}
}
