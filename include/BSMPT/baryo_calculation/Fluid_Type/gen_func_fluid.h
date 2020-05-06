/*
 * gen_func_fluid.cpp
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


#ifndef GEN_FUNC_FLUID_H
#define GEN_FUNC_FLUID_H
/**
 * @file
 */
#include<iostream>
#include <BSMPT/models/IncludeAllModels.h>
#include <BSMPT/baryo_calculation/transport_equations.h>
#include <boost/numeric/odeint.hpp>
#include <boost/math/interpolators/cubic_b_spline.hpp>
#include <BSMPT/models/SMparam.h>

namespace BSMPT{
namespace Baryo{
using namespace boost::numeric::odeint;
typedef std::vector<double> state_type;

/**
 * @brief C_smallcut For numerical stability
 */
const double C_smallcut = 1e-20;

/**
 * @brief The Calc_Gam_M class Class instance for the numerical evaluation of the relaxation rate at given temperature and quark mass.
 */
class Calc_Gam_M
{
    private:
     /**
         * @brief hf Calculates the h_f function as in arXiv:1710.04061
         * @param w = sqrt(m^2+k^2)
         * @return hf(w)
         */
        double hf(double w);
    /**
         * @brief hf_prime Calculates the derivative of h_f function as in arXiv:1710.04061
         * @param w = sqrt(m^2+k^2)
         * @return hf'(w)
         */
        double hf_prime(double w);
    public:
        double m_full,beta,T,vw,GamT,YukCoupling;
        double msqrt_thermal,dmsqrt_thermal;
        /**
         * @brief set_class Defines all needed parameter for the evaluation of Gam_M
         * @param T Temperature
         * @param vw Wall velocity
         * @param mt Quark mass
         */
        void set_class(double T , double vw, double mt , double m_thermal,double dm_thermal);
        /**
         * @brief operator () Needed for the boost interface.
         * @param Gam Current state of Gam_M.
         * @param dGam Current derivative of Gam_M.
         * @param k Integration variable/Momentum.
         */
        void operator()(const state_type &Gam, state_type &dGam,const double k);
};
/**
 * @brief Nintegrate_GamM The numerical evaluation of the relaxation rates at given temperature and quark mass.
 * @param C_Gam Class reference to pass the parameters.
 * @return Relaxation rate at given temperature and quark mass.
 */
double Nintegrate_GamM(Calc_Gam_M& C_Gam);
/**
 * @brief The Calc_Scp class Class instance for the numerical evaluation of the CP-violating source terms at given temperature and quark mass.
 */
class Calc_Scp
{
    private:
    /**
         * @brief nf Distribution function of fermions at given temperature and mass.
         * @param w = sqrt(m^2+k^2)
         * @return Distribution function of a fermion.
         */
        double nf(double w);
    /**
         * @brief nf Derivative of the distribution function of fermions at given temperature and mass.
         * @param w = sqrt(m^2+k^2)
         * @return Derivative of the distribution function of a fermion.
         */
        double nf_prime(double w);
    public:
        double m_full,T,vw,theta_prime;
        double msqrt_thermal,dmsqrt_thermal;
        int YukType;
        /**
         * @brief operator () Needed for the boost interface.
         * @param Scp Current state of Scp.
         * @param dScp Current derivative of Scp.
         * @param k Integration variable/Momentum.
         */
        void operator()(const state_type &Scp , state_type &dScp , const double k);
        /**
         * @brief set_class Defines all needed parameter of the class Calc_Scp
         * @param T Temperature
         * @param vw Wall velocity
         * @param mt Quark mass
         * @param theta_prime Derivative of the quark phase in respect of the wall distance z
         */
        void set_class(double T , double vw , double mt , double theta_prime , double msqrt_thermal_in,double dmsqrt_thermal_in);
};
/**
 * @brief Nintegrate_Scp Numerical evaluation of the CP-violating source terms at given temperature and quark mass.
 * @param C_Scp Class reference to pass the parameters.
 * @return Value of Scp at given temperature and quark mass.
 */
double Nintegrate_Scp(Calc_Scp& C_Scp);

/**
 * @brief The Calc_kappa_t class Class instance for the numerical calculation of the statistical factor kappa for a given quark mass and temperature.
 */
class Calc_kappa_t
{
    private:
    public:
        double mt,Temp;
        /**
         * @brief set_class Set up of the class parameters.
         * @param Temp_in Temperature
         * @param mt_in Quark mass
         */
        void set_class(double Temp_in , double mt_in);
        /**
         * @brief operator () Needed for the boost interface.
         */
        void operator() (const state_type& , state_type& , const double );
};
/**
 * @brief NIntegrate_kappa Numerical Numerical evaluation of the statistical factor
 * @param C_kap Class reference to pass the needed parameters
 * @return Value of the statistical factor at given temperature and mass.
 */
double NIntegrate_kappa(const Calc_kappa_t& C_kap);
/**
 * Numerical evaluation of the eta value; 
 * Checked with 1710.04061 
 * GammaWS = 3*GammaWS_Huber
 */
/**
 * @brief The Calc_eta class Class instance for the numerical evaluation of the eta value.
 */
class Calc_eta
{
    private:
        double Temp,vw;
        /**
         * @brief nL_cub Boost spline for the left-handed fermion density as function of the bubble wall distance z.
         */
        boost::math::cubic_b_spline<double> nL_cub;
        double exponent_prefactor;
        std::vector<double> array_z,array_nL;
    public:
        double prefactor;
        void set_class(std::vector<double> array_z,std::vector<double> array_nL , double Temp,double vw); 
        /**
         * @brief set_class Set function to set all class parameters.
         * @param arr Grid for the cubic spline of the left-handed quark density as function of the bubble wall distance z.
         * arr.first--> z
         * arr.second--> n_L
         * @param Temp Temperature
         * @param vw Bubble wall velocity
         */
        void set_class(std::pair<std::vector<double> , std::vector<double>> arr , double Temp , double vw);
        /**
         * @brief operator () Needed for the boost interface.
         * @param eta Current state of eta.
         * @param deta Current state of the derivative of eta.
         * @param z Bubble wall distance.
         */
        void operator()(const state_type& eta, state_type& deta, const double z);
};
/**
 * @brief Nintegrate_eta Numerical evaluation of the eta value at a given distance to the bubble wall.
 * @param C_eta Class reference to pass all needed parameters.
 * @param z_start Distance to the wall where the chemical potentials are assumend to vanish.
 * @param z_end Distance at which eta is evaluated.
 * @return Numerical value of eta.
 */
double Nintegrate_eta(const Calc_eta& C_eta, const double& z_start, const double& z_end);

/**
 * @brief The gen_fluid class Class instance overhead for all transport classes. Including all common functions and parameters.
 */
class gen_fluid
{
    private:
    public: 
        //TODO: Implement the full track for the flags
        int  bot_mass_flag;
        int  tau_mass_flag = 1 ;
        /**
         * @brief modelPointer BSMPT model pointer to get a interface to the BSMPT framework.
         */
        std::shared_ptr<Class_Potential_Origin> modelPointer;
        /**
         * @brief par Input parameters.
         */
        std::vector<double> par;
        /**
         * @brief parCT Counterterms of the input parameters.
         */
        std::vector<double> parCT;
        /**
         * @brief vcritical Critical VEV configuration at critical temperature TC.
         */
        std::vector<double> vcritical;
        /**
         * @brief gen_vcritical Critical VEV configuration in the full field basis.
         */
        std::vector<double>gen_vcritical;
        /**
         * @brief vsym VEV configuration in the symmetric potential.
         */
        std::vector<double>vsym;
        /**
         * @brief tanbeta TanBeta
         */
        double tanbeta;
        /**
         * @brief Yuk_Type Type of the C2HDM (Type 1 and 2 are implemented)
         */
        int Yuk_Type;
        /**
         * @brief symmetric_CP_violating_phase The CP-violating phase in the symmetric vacuum (TOP quark).
         */
        double symmetric_CP_violating_phase;
        /**
         * @brief broken_CP_violating_phase The CP-violating phase in the broken vacuum (TOP quark).
         */
        double broken_CP_violating_phase;
        /**
         * @brief TOP_symmetric_CP_violating_phase The CP-violating phase of the top quark in the symmetric vacuum.
         */
        double TOP_symmetric_CP_violating_phase;
        /**
         * @brief TOP_broken_CP_violating_phase The CP-violating phase of the top quark in the broken vacuum.
         */
        double TOP_broken_CP_violating_phase;
        /**
         * @brief BOT_symmetric_CP_violating_phase The CP-violating phase of the bot quark in the broken vacuum.
         */
        double BOT_symmetric_CP_violating_phase;
        /**
         * @brief BOT_broken_CP_violating_phase The CP-violating phase of the bot quark in the symmetric vacuum.
         */
        double BOT_broken_CP_violating_phase;
        /**
         * @brief TAU_symmetric_CP_violating_phase The CP-violating phase of the tau lepton in the symmetric vacuum.
         */
        double TAU_symmetric_CP_violating_phase;
        /**
         * @brief TAU_broken_CP_violating_phase The CP-violating phase of the tau lepton in the broken vacuum.
         */
        double TAU_broken_CP_violating_phase;
        /**
         * @brief Temp Temperature
         */
        double Temp;
        /**
         * @brief vw Bubble wall velocity.
         */
        double vw;
        /**
         * @brief LW Bubble wall thickness in the Kink-profile ansatz.
         */
        double LW; 

        double gS       =   1.23;
        double g        =   0.65;
        double gprime   =   0.36;
        double alphaS   =   1./7;
        double alphaW   =   1./30;
        double mtop_0   =   C_MassTop;
        double mbot_0   =   C_MassBottom;
        double mtau_0   =   C_MassTau;

        /**
         * @brief Dq Diffusion constant of quarks.
         */
        double Dq;
        /**
         * @brief Dt Diffusion constant of the top quark.
         */
        double Dt;
        /**
         * @brief Dh Diffusion constant of a Higgs boson.
         */
        double Dh;
        /**
         * @brief Dlep Diffusion constant of leptons.
         */
        double Dlep;
        /**
         * @brief Dtau Diffusion constant of the tau lepton.
         */
        double Dtau;
        double yuk_q;
        /**
         * @brief Gam_SS Decay width of the strong sphaleron transition.
         */
        double Gam_SS;
        /**
         * @brief Gam_Y Decay width of Yukawa interactions.
         */
        double Gam_Y;
        /**
         * @brief Gam_Y_t Decay width of the top quark.
         */
        double Gam_Y_t; 
        /**
         * @brief Gam_Y_b Decay width of the bot quark.
         */
        double Gam_Y_b;
        /**
         * @brief Gam_Y_tau Decay width of the tau quark.
         */
        double Gam_Y_tau;

        double kappa_QL_0   =   6;
        double kappa_H_0    =   4;
        double kappa_QR_0   =   3;

        

        Calc_Scp Calc_Scp_obj; //Class reference to the numerical evaulation of the CP-violating source term
        Calc_kappa_t Calc_kappa_obj;//Class reference to the numerical evaluation of the statistical factor kappa.
        Calc_Gam_M Calc_Gam_obj; // Class reference to the numerical evaluation of the relaxation rates.

        /**
         * @brief set_class Set-up of all class parameters with given input values.
         * @param bottom_mass_inp Flag how to treat the bot mass.
         * @param container Struct containing all needed parameters.
         * @param calc_Gam_obj Class reference for the numerical evaluation of the relaxation rate.
         * @param Calc_Scp_obj Class reference for the numerical evaluation of the CP-violating sources.
         * @param Calc_kappa_inp Class reference for the numerical evaluation of the statistical factor.
         */
        void set_class
            (
            const int bottom_mass_inp,
            struct GSL_integration_mubl& container,
            const Calc_Gam_M& calc_Gam_obj,
            const Calc_Scp& Calc_Scp_obj,
            const Calc_kappa_t& Calc_kappa_inp
            );
        double msqrt_thermal_top,dmsqrt_thermal_top;
        double msqrt_thermal_bot,dmsqrt_thermal_bot;
        double msqrt_thermal_tau,dmsqrt_thermal_tau;
        /**
         * @brief Calc_ThermalMass_q Calculates the right-handed thermal mass of the quark and the difference of right- and left-handed quark thermal mass
         * @param YukCoupling_in Yukawa coupling of the quark
         * @param T_in Critical temperature
         * @return Pair of right-handed thermal mass and the mass difference of right- and left-handed thermal mass
         */
        std::pair<double,double> Calc_ThermalMass_q(double& YukCoupling_in , double& T_in);
        /**
         * @brief Calc_ThermalMass_l Calculates the right-handed thermal mass of the lepton and the difference of right- and left-handed quark thermal mass
         * @param YukCoupling_in Yukawa coupling of the lepton
         * @param T_in Critical temperature
         * @return Pair of right-handed thermal mass and the mass difference of right- and left-handed thermal mass
         */
        std::pair<double,double> Calc_ThermalMass_l(double& YukCoupling_in , double& T_in);

        //TODO: RETURN values for top_func and tau_func.
        /**
         * @brief top_func Calculation of the top and bot mass at a given distance of the bubble wall and temperature
         * @param z Wall distance.
         * @param m_quark m_quark[0]--> top quark mass; m_quark[1]-->bos_mass
         * @param m_quark_prime m_quark_prime[0/1]-->Derivative of the top/bot quark mass in respect of the wall distance z.
         */
        void top_func(double z, std::vector<double>& m_quark, std::vector<double>& m_quark_prime);
        /**
         * @brief tau_func Calculation of the tau mass and derivative.
         * @param z Wall distance.
         * @param m_lep Returns the tau lepton mass.
         * @param m_lep_prime Return the derivative of the tau lepton mass in respect of the bubble wall distance z.
         */
        void tau_func(double z, std::vector<double>& m_lep, std::vector<double>& m_lep_prime);
        std::vector<double> omegaprime(double z);
        double atan2_mod(double Im , double Re);
        /**
         * @brief Calc_theta Caclulates the phase profile over the bubble wall by assuming the Kink-profile ansatz.
         * @param z Distance of the bubble wall.
         * @param CP_sym CP-violating phase in the symmetric vacuum.
         * @param CP_brk CP-violating phase in the broken vacuum.
         * @return The value of the CP-violating phase at given distance of the bubble wall z and the boundary values CP_sym and CP_brk
         */
        std::vector<double> Calc_theta(double z , double CP_sym , double CP_brk);
};

}
}

#endif