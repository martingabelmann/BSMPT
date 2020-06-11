/*
 * Class_PNGDM.h
 *
 *  Copyright (C) 2018  Philipp Basler and Margarete Mühlleitner

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

/**
  * @file
  */

#pragma once

#include <string> // for string
#include <vector> // for vector

#include <BSMPT/models/ClassPotentialOrigin.h>
namespace BSMPT
{
  namespace Models
  {

    /**
 * @brief The Class_Template class
 * Template for implementing a new model
 */
    class Class_PNGDM : public Class_Potential_Origin
    {
    public:
      Class_PNGDM();
      virtual ~Class_PNGDM();

      double muHsqr, LamH, muSsqr, LamS, LamHS, mXsq;
      double CTmuHsq, CTLamH, CTmuSsqr, CTLamS, CTLamHS, CTmXsqr, dTh, dTs;
      double g1 = C_gs;
      double g2 = C_g;
      double vh, vs;

      void ReadAndSet(const std::string &linestr, std::vector<double> &par) override;
      std::vector<std::string> addLegendCT() const override;
      std::vector<std::string> addLegendTemp() const override;
      std::vector<std::string> addLegendTripleCouplings() const override;
      std::vector<std::string> addLegendVEV() const override;


      /**
       * @brief Set the gen object
       * 
       * @param par[0] = muHsq
       * @param par[1] = LamH
       * @param par[2] = muSsqr
       * @param par[3] = LamS
       * @param par[4] = LamHS
       * @param par[5] = mXsq
       * @param par[6] = vH
       * @param par[7] = vS
       */
      void set_gen(const std::vector<double> &par) override;
      void set_CT_Pot_Par(const std::vector<double> &par) override;
      void write() const override;

      void TripleHiggsCouplings() override;
      std::vector<double> calc_CT() const override;

      void SetCurvatureArrays() override;
      bool CalculateDebyeSimplified() override;
      bool CalculateDebyeGaugeSimplified() override;
      double VTreeSimplified(const std::vector<double> &v) const override;
      double VCounterSimplified(const std::vector<double> &v) const override;
      void Debugging(const std::vector<double> &input, std::vector<double> &output) const override;
    };

  } // namespace Models
} // namespace BSMPT
