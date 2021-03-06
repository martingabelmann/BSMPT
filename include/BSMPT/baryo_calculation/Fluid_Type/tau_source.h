#ifndef TAU_SOURCE_H
#define TAU_SOURCE_H

/*
 * tau_source.h
 *
 *  Copyright (C) 2020  Philipp Basler, Margarete Mühlleitner and Jonas Müller

		This program is free software: you can redistribute it and/or modify
		it under the terms of the GNU General Public License as published by
		the Free Software Foundation, either version 3 of the License, or
		(at your option) any later version.

		This program is distributed in the hope that it will be useful,
		but WITHOUT ANY WARRANTY; without even the implied warranty of
		MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
		GNU General Public License for more details.

		You should have received a copy of the GNU General Public License
		along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include<iostream>
#include <BSMPT/models/IncludeAllModels.h>
#include <boost/numeric/odeint.hpp>
#include <boost/math/interpolators/cubic_b_spline.hpp>
#include <BSMPT/baryo_calculation/Fluid_Type/gen_func_fluid.h>

/**
 * @file
 */


namespace BSMPT{
namespace Baryo{

using namespace boost::numeric::odeint;
typedef std::vector<double> state_type;

/**
 * @brief The tau_source class Class instance for the transport equation system and numerical calculation of the left-handed fermion density in front of the bubble wall including the bot,top quark and the tau lepton.
 */
class tau_source : public gen_fluid
{
    private:
    public: 
    /**
         * @brief operator () Needed for the numerical solution via boost.
         * @param omega Vector of all included (rescaled) chemical potentials (top,bot,tau,massles leptons/quarks, neutrinos and bot higgs bosons).
         * @param domega Vector of all derivatives of the (rescaled) chemical potentials in respect to the wall distance z.
         * @param z Current wall distance.
         */
        void operator()(const state_type &omega, state_type &domega , const double z);
        /**
         * @brief Calc_nL Calculates the sum of all left-handed fermions in front of the bubble wall contributing to the SU(2) sphaleron transition.
         * @param z_start Boundary condition for the bubble wall distance where the chemical potentials are assumed to vanish.
         * @param z_end Bubble wall distance where nL is evaluated.
         * @return Returns the left-handed fermion density in front of the bubble wall evaluated at z_end.
         */
        double Calc_nL(double z_start , double z_end) const ;

};

}//namespace: Baryo
}//namespace: BSMPT

#endif
