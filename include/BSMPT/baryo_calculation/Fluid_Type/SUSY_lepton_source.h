#ifndef SUSY_LEPTON_SOURCE_H
#define SUSY_LEPTON_SOURCE_H

#include <iostream>
#include <BSMPT/models/IncludeAllModels.h> //Needed to call models from BSMPT
#include <boost/numeric/odeint.hpp>
#include <boost/math/interpolators/cubic_b_spline.hpp>
#include <BSMPT/baryo_calculation/Fluid_Type/gen_func_fluid.h>

namespace BSMPT
{
    namespace Baryo
    {

        using namespace boost::numeric::odeint;
        typedef std::vector<double> state_type; // Common type def for the boos implementation

        class SUSY_lepton_source : public gen_fluid
        {
        private:
        public:
            using gen_fluid::set_class;
            /**
             * @brief Set the class object Overloaded set_class function to additionally setting kappaQL and LambdaQL
             * 
             * @param bottom_mass_inp  Defined as in gen_fluid::set_class
             * @param container Defined as in gen_fluid::set_class
             * @param Calc_Gam_inp Defined as in gen_fluid::set_class
             * @param Calc_Scp_inp Defined as in gen_fluid::set_class
             * @param Calc_kappa_inp Defined as in gen_fluid::set_class
             * @param kappaQL_inp kappaQL O(1) prefactor for the interaction rate GammaQL
             * @param LambdaQL_inp Energy scale of the dim6 operator
             */
            void set_class(
                const int bottom_mass_inp,
                struct GSL_integration_mubl &container,
                const Calc_Gam_M &Calc_Gam_inp,
                const Calc_Scp &Calc_Scp_inp,
                const Calc_kappa_t &Calc_kappa_inp,
                const double &LambdaQL_inp,
                const double &mchi_inp,
                const int &source_flag);
            double mchi;
            double LambdaQL;
            int source_flag;


            /**
         * @brief Differential system of equations for the transport equations --> Second order possible
         * 
         * @param omega Vector containing mu and its first derivatives
         * @param domega Vector containing muprime and muprimeprime
         * @param z Wall distance z 
         */
            void operator()(const state_type &omega, state_type &domega, const double z);
            /**
         * @brief Evaluates the left-handed fermion density in front of the bubble wall --> Needed to create the grid
         * 
         * @param z_start Boundary condition for the bubble wall distance where the chemical potentials are assumed to vanish
         * @param z_end Bubble wall distance where nL is evaluated
         * @return double Returns the left-handed fermion density in front of the bubble wall evaluated at z_end
         */
            double Calc_nL(double z_start, double z_end) const;
        };

    } // namespace Baryo
} // namespace BSMPT

#endif
