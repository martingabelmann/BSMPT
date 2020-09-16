#include <BSMPT/baryo_calculation/Fluid_Type/BA_template.h>
#include <BSMPT/utility.h>

namespace BSMPT
{
    namespace Baryo
    {

        void BA_template::operator()(const state_type &omega, state_type &domega, const double z)
        {
            //TODO: Implement transport equations as differential equation system
            (void) omega;
            (void) domega;
            (void) z;
        }
        double BA_template::Calc_nL(double z_start, double z_end) const
        {
            //TODO: Use boost to solve transport equations
            const double C_AbsErr = 1e-9;
            const double C_RelErr = 1e-5;
            double stepsize_initial;
            if (z_start < z_end)
                stepsize_initial = 1e-8;
            if (z_start > z_end)
                stepsize_initial = -1e-8;
            return 0;
        }

    } // namespace Baryo
} // namespace BSMPT