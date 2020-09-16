

#include <BSMPT/baryo_calculation/Fluid_Type/SUSY_lepton_source.h>
#include <BSMPT/utility.h>

namespace BSMPT
{
    namespace Baryo
    {

        void SUSY_lepton_source::operator()(const state_type &omega, state_type &domega, const double z)
        {
            bool debug = false;
            if (debug)
                std::cout << "Begin of debug in " << __func__ << std::endl;

            /*
        Definition of all transport coefficients
    */
            std::vector<double> quark_mass;
            std::vector<double> quark_mass_prime;
            top_func(z, quark_mass, quark_mass_prime);
            double mt = quark_mass[0];
            double sym_phase = gen_fluid::symmetric_CP_violating_phase;
            double brk_phase = gen_fluid::broken_CP_violating_phase;
            auto theta_vec = Calc_theta(z, sym_phase, brk_phase);

            double theta_prime = theta_vec[1];
            //Calculation of kappa_t
            Calc_kappa_obj.set_class(Temp, mt);
            double num_int = NIntegrate_kappa(Calc_kappa_obj);
            double kappa_q = kappa_QL_0 * num_int;
            double kappa_tR = kappa_QR_0 * num_int;

            // Calculation for the new lepton chi
            std::vector<double> tau_mass, tau_mass_prime;
            tau_func(z, tau_mass, tau_mass_prime);
            double mtau_z = tau_mass[0];
            double mchi_z = (mtau_z / C_MassTau) * mchi;
            auto theta_vec_chi = Calc_theta(z, gen_fluid::TAU_symmetric_CP_violating_phase, gen_fluid::TAU_broken_CP_violating_phase);
            double theta_prime_chi = theta_vec_chi[1];
            //Chi Statistical factor
            Calc_kappa_obj.set_class(Temp, mchi_z);
            num_int = NIntegrate_kappa(Calc_kappa_obj);
            double kappa_chiL = kappa_LL_0 * num_int;
            double kappa_chiR = kappa_RL_0 * num_int;
            double yuk_chi;

            if (Yuk_Type == 1)
                yuk_chi = std::sqrt(2) * mchi / (C_vev0 * std::sin(std::atan(tanbeta)));
            if (Yuk_Type == 2)
                yuk_chi = std::sqrt(2) * mchi / (C_vev0 * std::cos(std::atan(tanbeta)));
            if ((Yuk_Type != 1) and (Yuk_Type != 2))
                throw(std::runtime_error("No valid yuktype in SUSY_lepton_source"));
            auto thermal_mass_chi = Calc_ThermalMass_l(yuk_chi, Temp);
            double msqrt_thermal_chi = thermal_mass_chi.first;
            double dmsqrt_thermal_chi = thermal_mass_chi.second;
            //Rescaled chemical potential like in 1811.11104
            /*
                omega[0] -> q
                omega[1] -> t
                omega[2] -> h1
                omega[3] -> h2
                omega[4] -> l
                omega[5] -> chi
                omega[6] -> nu
                omega[7] -> q_prime
                omega[8] -> t_prime
                omega[9] -> h1_prime
                omega[10] - >h2_prime
                omega[11] -> l_prime
                omega[12] -> chi_prime
                omega[13] -> nu_prime
            */
            double mu_M = (omega[1] / kappa_tR - omega[0] / kappa_q);
            double bR = -(omega[1] + omega[0]);                                   //local baryon number conservation yields t+ b + q = 0 --> b = - (q+t)
            double mu_SS = 2 * omega[0] / kappa_q - omega[1] / kappa_tR - 3 * bR; //Used q1 = -2 b and mb = 0
            double mu_Y = omega[1] / kappa_tR - omega[0] / kappa_q - omega[2] / kappa_H_0 - omega[3] / kappa_H_0;

            //Chi parts
            double mu_M_chi = omega[5] / kappa_chiR - omega[4] / kappa_chiL;
            double mu_Y_chi = omega[5] / kappa_chiR - omega[4] / kappa_chiL + omega[2] / kappa_H_0 + omega[3] / kappa_H_0;
            Calc_Gam_obj.set_class(Temp, vw, mchi_z, msqrt_thermal_chi, dmsqrt_thermal_chi);
            double Gam_M_chi = Nintegrate_GamM(Calc_Gam_obj);
            Calc_Scp_obj.set_class(Temp, vw, mchi_z, theta_prime_chi, msqrt_thermal_chi, dmsqrt_thermal_chi);
            double Scp_chi = Nintegrate_Scp(Calc_Scp_obj);
            double Gam_Y_chi = 0.28 * alphaW * yuk_chi * yuk_chi * Temp;
            double mu_QL = omega[4] / kappa_chiL - omega[5] / kappa_chiR - omega[0] / kappa_q + omega[1] / kappa_tR;
            double Gamma_QL = 0;
            if (LambdaQL > 1)
                Gamma_QL = 1 * std::pow(Temp, 5) / std::pow(LambdaQL, 4);
            else if (LambdaQL == 0)
                Gamma_QL = 0;
            else if (LambdaQL < 0)
                throw std::runtime_error("Negative energy scale LambdaQL");

            Calc_Gam_obj.set_class(Temp, vw, mt, msqrt_thermal_top, dmsqrt_thermal_top);
            double Gam_M = Nintegrate_GamM(Calc_Gam_obj);
            Calc_Scp_obj.set_class(Temp, vw, mt, theta_prime, msqrt_thermal_top, dmsqrt_thermal_top);
            double Scp = Nintegrate_Scp(Calc_Scp_obj);
            if (debug)
            {
                std::cout << "GammaQL = " << sep << Gamma_QL << std::endl;
                std::cout << "LambdaQL = " << sep << LambdaQL << std::endl;
                std::cout << "mChi_inp = " << sep << mchi << std::endl;
                std::cout << "mchi(z) = " << sep << mchi_z << std::endl;
            }
            double Chi_source = 0;
            double stretch_factor;
            /////////////////////////////////////////////////
            stretch_factor = 1;
            /////////////////////////////////////////////////
            if (source_flag == -1)
                throw std::runtime_error("No valid source flag set!");
            if (source_flag == 0)
                Chi_source = 0;
            else
            {
                //Calculation of Scp(z=0) to be the amplitude of the new source term
                std::vector<double> tau_mass0, tau_mass_prime0;
                tau_func(0, tau_mass0, tau_mass_prime0);
                double mtau0 = tau_mass0[0];
                double mchi_0 = (mtau0 / C_MassTau) * mchi;
                auto theta_vec_chi_0 = Calc_theta(0, gen_fluid::TAU_symmetric_CP_violating_phase, gen_fluid::TAU_broken_CP_violating_phase);
                double theta_prime_chi_0 = theta_vec_chi_0[1];
                Calc_Scp_obj.set_class(Temp, vw, mchi_0, theta_prime_chi_0, msqrt_thermal_chi, dmsqrt_thermal_chi);
                double Lambda0 = Nintegrate_Scp(Calc_Scp_obj);
                // std::cout<< Lambda0 << std::endl;
                if (source_flag == 1)
                    Chi_source = Lambda0;
                if (source_flag == 2)
                    Chi_source = Lambda0 * (1 - std::tanh(stretch_factor * z / LW));
                if (source_flag == 3)
                    Chi_source = Lambda0 * std::tanh(stretch_factor * z / LW);
            }
            // if(gen_fluid::)
            /*
                omega[0] -> q
                omega[1] -> t
                omega[2] -> h1
                omega[3] -> h2
                omega[4] -> l
                omega[5] -> chi
                omega[6] -> nu
                omega[7] -> q_prime
                omega[8] -> t_prime
                omega[9] -> h1_prime
                omega[10] - >h2_prime
                omega[11] -> l_prime
                omega[12] -> chi_prime
                omega[13] -> nu_prime
            */
            domega[0] = omega[7];
            domega[1] = omega[8];
            domega[2] = omega[9];
            domega[3] = omega[10];
            domega[4] = omega[11];
            domega[5] = omega[12];
            domega[6] = omega[13];
            domega[7] = (vw * omega[7] - (Gam_M * mu_M + Gam_Y * mu_Y - 2 * Gam_SS * mu_SS - Scp + Gamma_QL * mu_QL)) / Dq;
            domega[8] = (vw * omega[8] - (-Gam_M * mu_M - Gam_Y * mu_Y + Gam_SS * mu_SS + Scp - Gamma_QL * mu_QL)) / Dq;
            domega[9] = (vw * omega[9] - (Gam_Y * mu_Y - Gam_Y_chi * mu_Y_chi)) / Dh;
            domega[10] = (vw * omega[10] - (Gam_Y * mu_Y - Gam_Y_chi * mu_Y_chi)) / Dh;
            domega[11] = (vw * omega[11] - (Gam_M_chi * mu_M_chi + Gam_Y_chi * mu_M_chi - Scp_chi - Gamma_QL * mu_QL - Chi_source)) / Dlep;  //4 Point interaction + additional source term
            domega[12] = (vw * omega[12] - (-Gam_M_chi * mu_M_chi - Gam_Y_chi * mu_M_chi + Scp_chi + Gamma_QL * mu_QL + Chi_source)) / Dlep; //4 Point interaction + additional source term
            domega[13] = (vw * omega[13] - (0)) / Dlep;                                                                                      //Neutrinos are decoupled from the thermal bath --> no mass/no Higgs interactions
        }
        double SUSY_lepton_source::Calc_nL(double z_start, double z_end) const
        {
            bool debug = false;
            if (debug)
                std::cout << "start of debug in " << __func__ << std::endl;
            if (debug)
                std::cout << "mchi=" << mchi << sep << "LambdaQL=" << LambdaQL << std::endl;
            if (debug)
                std::cout << "top class called" << std::endl;
            if (debug)
                std::cout << "source_flag = " << source_flag << std::endl;

            /*
                omega[0] -> q
                omega[1] -> t
                omega[2] -> h1
                omega[3] -> h2
                omega[4] -> l
                omega[5] -> chi
                omega[6] -> nu
                omega[7] -> q_prime
                omega[8] -> t_prime
                omega[9] -> h1_prime
                omega[10] - >h2_prime
                omega[11] -> l_prime
                omega[12] -> chi_prime
                omega[13] -> nu_prime
            */
            state_type mu(14);
            mu = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
            if (debug)
                std::cout << "Before ODE Calc:" << std::endl;
            if (debug)
                for (size_t i = 0; i < mu.size(); i++)
                    std::cout << "\tmu[" << i << "] = " << mu[i] << std::endl;
            const double C_AbsErr = 1e-9;
            const double C_RelErr = 1e-5;
            double stepsize_initial;
            if (z_start < z_end)
                stepsize_initial = 1e-8;
            if (z_start > z_end)
                stepsize_initial = -1e-8;
            double abs_err = C_AbsErr;
            double rel_err = C_RelErr;
            integrate_adaptive(make_controlled(abs_err, rel_err, error_stepper_type()), *this, mu, z_start, z_end, stepsize_initial);
            if (debug)
                std::cout << "After ODE Calc:" << std::endl;
            if (debug)
                for (size_t i = 0; i < mu.size(); i++)
                    std::cout << "\tmu[" << i << "] = " << mu[i] << std::endl;

            return 5 * mu[0] + 4 * mu[1] + mu[4];
        }
        void SUSY_lepton_source::set_class(const int bottom_mass_inp,
                                           struct GSL_integration_mubl &container,
                                           const Calc_Gam_M &Calc_Gam_inp,
                                           const Calc_Scp &Calc_Scp_inp,
                                           const Calc_kappa_t &Calc_kappa_inp,
                                           const double &LambdaQL_inp,
                                           const double &mchi_inp,
                                           const int &source_flag_inp)
        {
            LambdaQL = LambdaQL_inp;
            mchi = mchi_inp;
            source_flag = source_flag_inp;
            set_class(bottom_mass_inp, container, Calc_Gam_inp, Calc_Scp_inp, Calc_kappa_inp);
        }

    } // namespace Baryo
} // namespace BSMPT
