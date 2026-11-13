/** Captured hardening-driven shear/cap branch change; no fitted correction. */
#include "drucker_prager_cap_kernel.h"
#include <iostream>
using namespace SPH::dpc;
int main()
{
    try
    {
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
        for(double perturb:{-1e-8,0.,1e-8})
        {
            Tensor test=F;test(2,2)+=perturb;auto result=model.update(test,old);
            if(result.yield_residual>2e-8||result.plastic_stress_work< -1e-6||result.state.alpha<old.alpha)
                throw std::runtime_error("Invalid moving-junction return");
            if(perturb==0&&(result.surface!=2||std::abs(result.p-33043345.708660252)>.5||std::abs(result.q-74128630.8473845)>.5))
                throw std::runtime_error("Return differs from independently solved cap root");
        }
        std::cout<<"DPC moving shear/cap junction regression passed"<<std::endl;return 0;
    }
    catch(const std::exception &e){std::cerr<<e.what()<<std::endl;return 1;}
}
