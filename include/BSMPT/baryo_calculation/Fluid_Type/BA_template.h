#ifndef BA_TEMPLATE_H
#define BA_TEMPLATE_H

#include<iostream>
#include <BSMPT/models/IncludeAllModels.h>//Needed to call models from BSMPT
#include <boost/numeric/odeint.hpp>
#include <boost/math/interpolators/cubic_b_spline.hpp>
#include <BSMPT/baryo_calculation/Fluid_Type/gen_func_fluid.h>

namespace BSMPT{
namespace Baryo{

using namespace boost::numeric::odeint;
typedef std::vector<double> state_type; // Common type def for the boos implementation


class BA_template : public gen_fluid 
{
    private:
    public:
        /**
         * @brief Differential system of equations for the transport equations --> Second order possible
         * 
         * @param omega Vector containing mu and its first derivatives
         * @param domega Vector containing muprime and muprimeprime
         * @param z Wall distance z 
         */
        void operator()(const state_type & omega,state_type &domega , const double z);
        /**
         * @brief Evaluates the left-handed fermion density in front of the bubble wall --> Needed to create the grid
         * 
         * @param z_start Boundary condition for the bubble wall distance where the chemical potentials are assumed to vanish
         * @param z_end Bubble wall distance where nL is evaluated
         * @return double Returns the left-handed fermion density in front of the bubble wall evaluated at z_end
         */
        double Calc_nL(double z_start, double z_end) const;
};


}//namespace end: Baryo
}//namespace: BSMPT


#endif
