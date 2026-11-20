/** Captured hardening-driven shear/cap branch change; no fitted correction. */
#include "dpc_snsni_adapter.h"
#include "pd_execution.h"
#include <string_view>
#include <iostream>
using namespace SPH::dpc;
int main(int argc, char **argv)
{
    try
    {
        using namespace SPH::peridynamics;
        PeriExecutionMode mode=PeriExecutionMode::CPU;
        if (argc>1)
        {
            const std::string_view option(argv[1]);
            if (argc!=2 || option.substr(0,12)!="--execution=")
                throw std::invalid_argument("usage: dpc_moving_junction_check [--execution=cpu|gpu]");
            mode=parsePeriExecution(option.substr(12));
        }
        const PeriExecution execution(mode);
        const double E=2950288297.4204993,d=10600161.839936104,beta=1.0821041362364843,R=.5300000002997748;
        Model model({{.765,E,.13,d,beta,R,69929248.90063852},
                     {.790,E,.13,d,beta,R,72371837.73722216},
                     {.815,E,.13,d,beta,R,75168990.63679394}},.79,.1,160e6);
        State old;old.kappa=-0.00052340178242876278;old.alpha=0.0086326538989501544;
        old.inverse_Cp<<.99653049471791832,.00061699583439536873,-.011495617369410408,
            .00061699583439536873,.99254294707704183,-.002474922739529011,
            -.011495617369410408,-.002474922739529011,1.0101021678096256;
        Tensor F;F<<.99856426671999998,-.00064334308,-.0018425635600000001,
            -.00086636836000000004,1.0047445775999999,-.00034274784000000002,
            .032795676480000001,.0071206471600000001,.97324670223999998;
        auto batch=execution.makeDpcPoints(model);
        std::vector<PointInput> inputs;
        std::vector<Result> references;
        auto append=[&](const Tensor &deformation, const State &history, const Result &result)
        {
            PointInput point;
            for(int row=0;row<3;++row)
                for(int col=0;col<3;++col)
                {
                    point.deformation[row][col]=deformation(row,col);
                    point.history.inverse_Cp[row][col]=history.inverse_Cp(row,col);
                }
            point.history.kappa=history.kappa;
            point.history.alpha=history.alpha;
            inputs.push_back(point);
            references.push_back(result);
        };
        for(double perturb:{-1e-8,0.,1e-8})
        {
            Tensor test=F;test(2,2)+=perturb;auto result=model.update(test,old);
            if(result.yield_residual>2e-8||result.plastic_stress_work< -1e-6||result.state.alpha<old.alpha)
                throw std::runtime_error("Invalid moving-junction return");
            if(perturb==0&&(result.surface!=2||std::abs(result.p-33043345.708660252)>.5||std::abs(result.q-74128630.8473845)>.5))
                throw std::runtime_error("Return differs from independently solved cap root");
            append(test,old,result);
        }
        // Isotropic stretches exercise the repeated-eigenvalue path,
        // including hydrostatic cap compaction and the tensile apex.
        const State initial;
        const auto append_hydrostatic=[&](double stretch, int expected_branch)
        {
            const Tensor deformation=stretch*Tensor::Identity();
            const auto result=model.update(deformation,initial);
            if (result.surface!=expected_branch)
                throw std::runtime_error("Unexpected branch in hydrostatic DPC reference");
            append(deformation,initial,result);
        };
        append_hydrostatic(1.,0);
        append_hydrostatic(.98,2);
        append_hydrostatic(1.005,3);
        double maximum_error=0, trial_error=0, acoustic_error=0;
        if (batch)
        {
            // One invalid point must report failure without hiding adjacent
            // successes or changing their accepted histories.
            inputs.push_back(inputs.back());
            inputs.back().deformation[2][2]=-1;
            std::vector<PointOutput> outputs;
            batch->evaluate(inputs,outputs);
            if (outputs.size()!=inputs.size() || outputs.back().failure!=PointFailure::Deformation)
                throw std::runtime_error("Device DPC failed to report an invalid point");
            for (size_t i=0;i<references.size();++i)
            {
                const auto &expected=references[i];
                const auto &computed=outputs[i];
                if (computed.failure!=PointFailure::None)
                    throw std::runtime_error("Device DPC point "+std::to_string(i)+" failed with status "+
                                             std::to_string(static_cast<int>(computed.failure)));
                if (computed.surface!=expected.surface)
                    throw std::runtime_error("Device DPC point "+std::to_string(i)+" selected branch "+
                                             std::to_string(computed.surface)+", expected "+std::to_string(expected.surface));
                Tensor sigma,inverse_Cp;
                for (int row=0;row<3;++row)
                    for (int col=0;col<3;++col)
                    {
                        sigma(row,col)=computed.sigma[row][col];
                        inverse_Cp(row,col)=computed.history.inverse_Cp[row][col];
                    }
                maximum_error=std::max({maximum_error,(sigma-expected.sigma).norm()/std::max(1.,expected.sigma.norm()),
                    (inverse_Cp-expected.state.inverse_Cp).norm(),
                    std::abs(computed.history.kappa-expected.state.kappa),
                    std::abs(computed.history.alpha-expected.state.alpha)});
                if (computed.yield_residual>2e-8 || computed.plastic_stress_work < -1e-6 || maximum_error>1e-9)
                    throw std::runtime_error("Device DPC failed stress, history or yield consistency at the moving junction");
            }

            const TensorGradient grad_F{.12*Tensor::Identity(),-.27*Tensor::Identity(),.09*Tensor::Identity()};
            TensorGradient grad_Cp{Tensor::Zero(),Tensor::Zero(),Tensor::Zero()};
            grad_Cp[1](0,2)=grad_Cp[1](2,0)=.15;
            const Eigen::Vector3d grad_kappa(.4,-.3,.17);
            std::vector<SnsniPointInput> trial_inputs(2);
            // Check a nonuniform history and the repeated spectrum at rest.
            trial_inputs[0].material=inputs[1];
            trial_inputs[1].material=inputs[3];
            for (auto &point : trial_inputs)
                for (int k=0;k<3;++k)
                {
                    point.kappa_gradient[k]=grad_kappa[k];
                    for (int row=0;row<3;++row)
                        for (int col=0;col<3;++col)
                        {
                            point.deformation_gradient[k][row][col]=grad_F[k](row,col);
                            point.inverse_Cp_gradient[k][row][col]=grad_Cp[k](row,col);
                        }
                }

            std::vector<SnsniPointOutput> trial_results;
            constexpr double density=1229.24, cfl_length=5e-6;
            batch->evaluateSnsni(trial_inputs,density,cfl_length,trial_results);
            if (trial_results.size()!=trial_inputs.size())
                throw std::runtime_error("Device S-NSNI trial output size mismatch");
            for (size_t i=0;i<trial_results.size();++i)
            {
                const auto &computed=trial_results[i];
                if (computed.material.failure!=PointFailure::None)
                    throw std::runtime_error("Device S-NSNI trial update did not converge");
                const Tensor deformation=i==0 ? F : Tensor(Tensor::Identity());
                const auto &history=i==0 ? old : initial;
                const auto &accepted=references[i==0 ? 1 : 3];
                const auto expected=henckyTrialField(model,deformation,history,grad_F,grad_Cp,grad_kappa);
                const auto parameters=model.parameters(accepted.relative_density,accepted.state.alpha);
                const auto expected_dt=acousticStep(deformation,deformation.inverse(),parameters,density,cfl_length);
                acoustic_error=std::max(acoustic_error,std::abs(computed.acoustic_dt/expected_dt-1));
                Tensor stress;
                TensorGradient gradient;
                for (int row=0;row<3;++row)
                    for (int col=0;col<3;++col)
                    {
                        stress(row,col)=computed.trial_sigma[row][col];
                        for (int k=0;k<3;++k)
                            gradient[k](row,col)=computed.trial_sigma_gradient[k][row][col];
                    }
                trial_error=std::max(trial_error,(stress-expected.sigma).norm()/std::max(1.,expected.sigma.norm()));
                for (int k=0;k<3;++k)
                    trial_error=std::max(trial_error,(gradient[k]-expected.gradient[k]).norm()/expected.gradient[k].norm());
            }
            if (!(acoustic_error<1e-10 && trial_error<1e-9))
                throw std::runtime_error("Device S-NSNI trial stress, gradient or acoustic bound differs from CPU");
        }
        std::cout<<"DPC moving shear/cap junction regression passed ("<<execution.backendName()
                 <<", maximum relative stress or absolute history difference "<<maximum_error
                 <<", trial stress/gradient "<<trial_error<<", acoustic bound "<<acoustic_error<<")"<<std::endl;
        return 0;
    }
    catch(const std::exception &e){std::cerr<<e.what()<<std::endl;return 1;}
}
