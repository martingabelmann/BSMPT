/*
 * PNGDM.cpp
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

#include <ext/alloc_traits.h>	  // for __alloc_traits<>::value_type
#include <stddef.h>				  // for size_t
#include <algorithm>			  // for max, copy
#include <iostream>				  // for operator<<, endl, basic_o...
#include <memory>				  // for allocator_traits<>::value...
#include <BSMPT/models/SMparam.h> // for C_vev0, C_MassTop, C_g
#include "Eigen/Dense"
#include "Eigen/Eigenvalues"
#include "Eigen/IterativeLinearSolvers"

#include <BSMPT/models/ClassPotentialPNGDM.h>
#include <BSMPT/models/IncludeAllModels.h>
#include <BSMPT/utility.h>
using namespace Eigen;

/**
 * @file
 * Template for adding a new model class
 */

namespace BSMPT
{
	namespace Models
	{
		/**
 * Here you have to adjust NNeutralHiggs, NChargedHiggs, nPar (number of Lagrangian parameters AFTER
 *  using the tadpole conditions),
 * nParCT (number of counterterms) as well as nVEV (number of VEVs for minimization)
 */
		Class_PNGDM::Class_PNGDM()
		{
			Model = ModelID::ModelIDs::PNGDM; // global int constant which will be used to tell the program which model is called
			NNeutralHiggs = 4;				  // number of neutral Higgs bosons at T = 0
			NChargedHiggs = 2;				  // number of charged Higgs bosons  at T = 0 (all d.o.f.)

			nPar = 8;	// number of parameters in the tree-Level Lagrangian
			nParCT = 8; // number of parameters in the counterterm potential

			nVEV = 2; // number of VEVs to minimize the potential

			NHiggs = NNeutralHiggs + NChargedHiggs;

			VevOrder.resize(nVEV);
			/*field basis:
	- hc1		0
	- hc2 		1
	- sigmaH 	2
	- chi 		3
	- phiH 		4
	- phiS 		5
  double = { {hc1,hc2} , { phiH , sigmaH}}
  singlet = { phiS , chi} 
  */
			VevOrder[0] = 4;
			VevOrder[1] = 5;

			// Set UseVTreeSimplified to use the tree-level potential defined in VTreeSimplified
			UseVTreeSimplified = false;

			// Set UseVCounterSimplified to use the counterterm potential defined in VCounterSimplified
			UseVCounterSimplified = false;
		}

		Class_PNGDM::~Class_PNGDM()
		{
		}

		/**
 * returns a string which tells the user the chronological order of the counterterms. Use this to
 * complement the legend of the given input file
 */
		std::vector<std::string> Class_PNGDM::addLegendCT() const
		{
			std::vector<std::string> labels;
			labels.push_back("CTmuHsqr");
			labels.push_back("CTLamH");
			labels.push_back("CTmuSsqr");
			labels.push_back("CTLamS");
			labels.push_back("CTLamHS");
			labels.push_back("CTmXsqr");
			labels.push_back("dTh");
			labels.push_back("dTs");
			return labels;
		}

		/**
 * returns a string which tells the user the chronological order of the VEVs and the critical temperature. Use this to
 * complement the legend of the given input file
 */
		std::vector<std::string> Class_PNGDM::addLegendTemp() const
		{
			std::vector<std::string> labels;
			labels.push_back("T_c");
			labels.push_back("v_c");
			labels.push_back("xi_c");
			labels.push_back("omega_h");
			labels.push_back("omega_s");
			return labels;
		}

		std::vector<std::string> Class_PNGDM::addLegendTripleCouplings() const
		{
			std::vector<std::string> labels;
			std::vector<std::string> particles;
			particles.resize(NHiggs);
			//here you have to define the particle names in the vector particles

			particles[0] = "H";

			std::string out = "Tree_";
			for (size_t i = 0; i < NHiggs; i++)
			{
				for (size_t j = i; j < NHiggs; j++)
				{
					for (size_t k = j; k < NHiggs; k++)
					{
						labels.push_back("Tree_" + particles.at(i) + particles.at(j) + particles.at(k));
						labels.push_back("CT_" + particles.at(i) + particles.at(j) + particles.at(k));
						labels.push_back("CW_" + particles.at(i) + particles.at(j) + particles.at(k));
					}
				}
			}

			return labels;
		}

		/**
 * returns a string which tells the user the chronological order of the VEVs. Use this to
 * complement the legend of the given input file
 */
		std::vector<std::string> Class_PNGDM::addLegendVEV() const
		{
			std::vector<std::string> labels;
			labels.push_back("omega_1");
			labels.push_back("omega_2");
			return labels;
		}

		/**
 * Reads the string linestr and sets the parameter point
 */
		void Class_PNGDM::ReadAndSet(const std::string &linestr, std::vector<double> &par)
		{
			std::stringstream ss(linestr);
			double tmp;

			if (UseIndexCol)
			{
				ss >> tmp;
			}

			for (int k = 1; k <= 12; k++)
			{
				// to Tizian's convention: {lambda, d2, delta2} --> {2lH, 2lS, 2lHS}
				// msq = -muH^2 in Tizian's convention
				// b2 = -muS^2 in Tizian's convention
				ss >> tmp;
				if (k == 3)
					par[6] = tmp; //v
				else if (k == 4)
					par[7] = tmp; //vS
				else if (k == 5)
					par[5] = tmp * tmp; //mx^2 --> Input is mX...
				else if (k == 8)
					par[1] = tmp / 2.; //lambda converted to LamH
				else if (k == 9)
					par[3] = tmp / 2.; // d2 converted to LamS
				else if (k == 10)
					par[4] = tmp / 2.; // delta2 converted to LamHS
				else if (k == 11)
					par[0] = -tmp; //msq converted to muHsq
				else if (k == 12)
					par[2] = -tmp; //b2 converted to muSsqr
			}
			set_gen(par);
			return;
		}

		/**
 * Set Class Object as well as the VEV configuration
 */
		void Class_PNGDM::set_gen(const std::vector<double> &par)
		{

			muHsqr = par[0];
			LamH = par[1];
			muSsqr = par[2];
			LamS = par[3];
			LamHS = par[4];
			mXsq = par[5];
			vh = par[6];
			vs = par[7];

			scale = C_vev0;

			vevTreeMin.resize(nVEV);
			vevTree.resize(NHiggs);

			// Here you have to set the vector vevTreeMin. The vector vevTree will then be set by the function MinimizeOrderVEV
			vevTreeMin[0] = vh;
			vevTreeMin[1] = vs;

			vevTree = MinimizeOrderVEV(vevTreeMin);
			if (!SetCurvatureDone)
				SetCurvatureArrays();
		}

		/**
 * set your counterterm parameters from the entries of par as well as the entries of Curvature_Higgs_CT_L1 to
 * Curvature_Higgs_CT_L4.
 */
		void Class_PNGDM::set_CT_Pot_Par(const std::vector<double> &par)
		{

			/*
		 	parCT[0]=CTmXsqr 
			parCT[1]=CTmuSsqr
			parCT[2]=CTLamS
			parCT[3]=CTmuHsqr
			parCT[4]=CTLamHS
			parCT[5]=CTLamH
			parCT[6]=dTh
			parCT[7]=dTs
			*/
			CTmXsqr = par[0];
			CTmuSsqr = par[1];
			CTLamS = par[2];
			CTmuHsq = par[3];
			CTLamHS = par[4];
			CTLamH = par[5];
			dTh = par[6];
			dTs = par[7];

			Curvature_Higgs_CT_L1[4] = dTh;
			Curvature_Higgs_CT_L1[5] = dTs;

			Curvature_Higgs_CT_L2[0][0] = -CTmuHsq / 2.;
			Curvature_Higgs_CT_L2[1][1] = -CTmuHsq / 2.;
			Curvature_Higgs_CT_L2[2][2] = -CTmuHsq / 2.;
			Curvature_Higgs_CT_L2[3][3] = -CTmuSsqr / 2. + CTmXsqr / 2.;
			Curvature_Higgs_CT_L2[4][4] = -CTmuHsq / 2.;
			Curvature_Higgs_CT_L2[5][5] = -CTmuSsqr / 2. - CTmXsqr / 2.;

			Curvature_Higgs_CT_L4[0][0][0][0] = 3 * CTLamH;
			Curvature_Higgs_CT_L4[1][1][0][0] = CTLamH;
			Curvature_Higgs_CT_L4[2][2][0][0] = CTLamH;
			Curvature_Higgs_CT_L4[3][3][0][0] = CTLamHS;
			Curvature_Higgs_CT_L4[4][4][0][0] = CTLamH;
			Curvature_Higgs_CT_L4[5][5][0][0] = CTLamHS;
			Curvature_Higgs_CT_L4[1][0][1][0] = CTLamH;
			Curvature_Higgs_CT_L4[0][1][1][0] = CTLamH;
			Curvature_Higgs_CT_L4[2][0][2][0] = CTLamH;
			Curvature_Higgs_CT_L4[0][2][2][0] = CTLamH;
			Curvature_Higgs_CT_L4[3][0][3][0] = CTLamHS;
			Curvature_Higgs_CT_L4[0][3][3][0] = CTLamHS;
			Curvature_Higgs_CT_L4[4][0][4][0] = CTLamH;
			Curvature_Higgs_CT_L4[0][4][4][0] = CTLamH;
			Curvature_Higgs_CT_L4[5][0][5][0] = CTLamHS;
			Curvature_Higgs_CT_L4[0][5][5][0] = CTLamHS;
			Curvature_Higgs_CT_L4[1][0][0][1] = CTLamH;
			Curvature_Higgs_CT_L4[0][1][0][1] = CTLamH;
			Curvature_Higgs_CT_L4[0][0][1][1] = CTLamH;
			Curvature_Higgs_CT_L4[1][1][1][1] = 3 * CTLamH;
			Curvature_Higgs_CT_L4[2][2][1][1] = CTLamH;
			Curvature_Higgs_CT_L4[3][3][1][1] = CTLamHS;
			Curvature_Higgs_CT_L4[4][4][1][1] = CTLamH;
			Curvature_Higgs_CT_L4[5][5][1][1] = CTLamHS;
			Curvature_Higgs_CT_L4[2][1][2][1] = CTLamH;
			Curvature_Higgs_CT_L4[1][2][2][1] = CTLamH;
			Curvature_Higgs_CT_L4[3][1][3][1] = CTLamHS;
			Curvature_Higgs_CT_L4[1][3][3][1] = CTLamHS;
			Curvature_Higgs_CT_L4[4][1][4][1] = CTLamH;
			Curvature_Higgs_CT_L4[1][4][4][1] = CTLamH;
			Curvature_Higgs_CT_L4[5][1][5][1] = CTLamHS;
			Curvature_Higgs_CT_L4[1][5][5][1] = CTLamHS;
			Curvature_Higgs_CT_L4[2][0][0][2] = CTLamH;
			Curvature_Higgs_CT_L4[0][2][0][2] = CTLamH;
			Curvature_Higgs_CT_L4[2][1][1][2] = CTLamH;
			Curvature_Higgs_CT_L4[1][2][1][2] = CTLamH;
			Curvature_Higgs_CT_L4[0][0][2][2] = CTLamH;
			Curvature_Higgs_CT_L4[1][1][2][2] = CTLamH;
			Curvature_Higgs_CT_L4[2][2][2][2] = 3 * CTLamH;
			Curvature_Higgs_CT_L4[3][3][2][2] = CTLamHS;
			Curvature_Higgs_CT_L4[4][4][2][2] = CTLamH;
			Curvature_Higgs_CT_L4[5][5][2][2] = CTLamHS;
			Curvature_Higgs_CT_L4[3][2][3][2] = CTLamHS;
			Curvature_Higgs_CT_L4[2][3][3][2] = CTLamHS;
			Curvature_Higgs_CT_L4[4][2][4][2] = CTLamH;
			Curvature_Higgs_CT_L4[2][4][4][2] = CTLamH;
			Curvature_Higgs_CT_L4[5][2][5][2] = CTLamHS;
			Curvature_Higgs_CT_L4[2][5][5][2] = CTLamHS;
			Curvature_Higgs_CT_L4[3][0][0][3] = CTLamHS;
			Curvature_Higgs_CT_L4[0][3][0][3] = CTLamHS;
			Curvature_Higgs_CT_L4[3][1][1][3] = CTLamHS;
			Curvature_Higgs_CT_L4[1][3][1][3] = CTLamHS;
			Curvature_Higgs_CT_L4[3][2][2][3] = CTLamHS;
			Curvature_Higgs_CT_L4[2][3][2][3] = CTLamHS;
			Curvature_Higgs_CT_L4[0][0][3][3] = CTLamHS;
			Curvature_Higgs_CT_L4[1][1][3][3] = CTLamHS;
			Curvature_Higgs_CT_L4[2][2][3][3] = CTLamHS;
			Curvature_Higgs_CT_L4[3][3][3][3] = 3 * CTLamS;
			Curvature_Higgs_CT_L4[4][4][3][3] = CTLamHS;
			Curvature_Higgs_CT_L4[5][5][3][3] = CTLamS;
			Curvature_Higgs_CT_L4[4][3][4][3] = CTLamHS;
			Curvature_Higgs_CT_L4[3][4][4][3] = CTLamHS;
			Curvature_Higgs_CT_L4[5][3][5][3] = CTLamS;
			Curvature_Higgs_CT_L4[3][5][5][3] = CTLamS;
			Curvature_Higgs_CT_L4[4][0][0][4] = CTLamH;
			Curvature_Higgs_CT_L4[0][4][0][4] = CTLamH;
			Curvature_Higgs_CT_L4[4][1][1][4] = CTLamH;
			Curvature_Higgs_CT_L4[1][4][1][4] = CTLamH;
			Curvature_Higgs_CT_L4[4][2][2][4] = CTLamH;
			Curvature_Higgs_CT_L4[2][4][2][4] = CTLamH;
			Curvature_Higgs_CT_L4[4][3][3][4] = CTLamHS;
			Curvature_Higgs_CT_L4[3][4][3][4] = CTLamHS;
			Curvature_Higgs_CT_L4[0][0][4][4] = CTLamH;
			Curvature_Higgs_CT_L4[1][1][4][4] = CTLamH;
			Curvature_Higgs_CT_L4[2][2][4][4] = CTLamH;
			Curvature_Higgs_CT_L4[3][3][4][4] = CTLamHS;
			Curvature_Higgs_CT_L4[4][4][4][4] = 3 * CTLamH;
			Curvature_Higgs_CT_L4[5][5][4][4] = CTLamHS;
			Curvature_Higgs_CT_L4[5][4][5][4] = CTLamHS;
			Curvature_Higgs_CT_L4[4][5][5][4] = CTLamHS;
			Curvature_Higgs_CT_L4[5][0][0][5] = CTLamHS;
			Curvature_Higgs_CT_L4[0][5][0][5] = CTLamHS;
			Curvature_Higgs_CT_L4[5][1][1][5] = CTLamHS;
			Curvature_Higgs_CT_L4[1][5][1][5] = CTLamHS;
			Curvature_Higgs_CT_L4[5][2][2][5] = CTLamHS;
			Curvature_Higgs_CT_L4[2][5][2][5] = CTLamHS;
			Curvature_Higgs_CT_L4[5][3][3][5] = CTLamS;
			Curvature_Higgs_CT_L4[3][5][3][5] = CTLamS;
			Curvature_Higgs_CT_L4[5][4][4][5] = CTLamHS;
			Curvature_Higgs_CT_L4[4][5][4][5] = CTLamHS;
			Curvature_Higgs_CT_L4[0][0][5][5] = CTLamHS;
			Curvature_Higgs_CT_L4[1][1][5][5] = CTLamHS;
			Curvature_Higgs_CT_L4[2][2][5][5] = CTLamHS;
			Curvature_Higgs_CT_L4[3][3][5][5] = CTLamS;
			Curvature_Higgs_CT_L4[4][4][5][5] = CTLamHS;
			Curvature_Higgs_CT_L4[5][5][5][5] = 3 * CTLamS;
		}

		/**
 * console output of all Parameters
 */
		void Class_PNGDM::write() const
		{

			std::cout << "Model = " << Model << std::endl;

			std::cout << "The parameters are : " << std::endl;
			std::cout << sep << "muHsqr = " << sep << muHsqr << std::endl;
			std::cout << sep << "LamH = " << sep << LamH << std::endl;
			std::cout << sep << "muSsqr = " << sep << muSsqr << std::endl;
			std::cout << sep << "LamS = " << sep << LamS << std::endl;
			std::cout << sep << "LamHS = " << sep << LamHS << std::endl;
			std::cout << sep << "mXsqr = " << sep << mXsq << std::endl;
			std::cout << sep << "vh = " << sep << vh << std::endl;
			std::cout << sep << "vs = " << sep << vs << std::endl;

			std::cout << "The counterterm parameters are : " << std::endl;
			std::cout << sep << "CTmxsqr = " << sep << CTmXsqr << std::endl;
			std::cout << sep << "CTmuSsqr = " << sep << CTmuSsqr << std::endl;
			std::cout << sep << "CTLamS = " << sep << CTLamS << std::endl;
			std::cout << sep << "CTmuHsqr = " << sep << CTmuHsq << std::endl;
			std::cout << sep << "CTLamHS = " << sep << CTLamHS << std::endl;
			std::cout << sep << "CTLamH = " << sep << CTLamH << std::endl;
			std::cout << sep << "dTh = " << sep << dTh << std::endl;
			std::cout << sep << "dTs = " << sep << dTs << std::endl;

			std::cout << "The scale is given by mu = " << scale << " GeV " << std::endl;
		}

		/**
 * Calculates the counterterms. Here you need to work out the scheme and implement the formulas.
 */
		std::vector<double> Class_PNGDM::calc_CT() const
		{
			std::vector<double> parCT;

			if (!SetCurvatureDone)
			{
				std::string retmes = __func__;
				retmes += " was called before SetCurvatureArrays()!\n";
				throw std::runtime_error(retmes);
			}
			if (!CalcCouplingsdone)
			{
				std::string retmes = __func__;
				retmes += " was called before CalculatePhysicalCouplings()!\n";
				throw std::runtime_error(retmes);
			}

			std::vector<double> WeinbergNabla, WeinbergHesse;
			WeinbergNabla = WeinbergFirstDerivative();
			WeinbergHesse = WeinbergSecondDerivative();

			VectorXd NablaWeinberg(NHiggs);
			MatrixXd HesseWeinberg(NHiggs, NHiggs), HiggsRot(NHiggs, NHiggs);
			for (size_t i = 0; i < NHiggs; i++)
			{
				NablaWeinberg[i] = WeinbergNabla[i];
				for (size_t j = 0; j < NHiggs; j++)
					HesseWeinberg(i, j) = WeinbergHesse.at(j * NHiggs + i);
			}
			/*
		 	parCT[0]=CTmXsqr 
			parCT[1]=CTmuSsqr
			parCT[2]=CTLamS
			parCT[3]=CTmuHsqr
			parCT[4]=CTLamHS
			parCT[5]=CTLamH
			parCT[6]=dTh
			parCT[7]=dTs
			*/
			// parCT.push_back(HesseWeinberg(3, 3) - NablaWeinberg[5] / vs);
			// parCT.push_back((vh * HesseWeinberg(4, 5) + vs * (-HesseWeinberg(3, 3) + HesseWeinberg(5, 5)) - 2 * NablaWeinberg[5]) / vs);
			// parCT.push_back((vs * HesseWeinberg(5, 5) - NablaWeinberg(5)) / std::pow(vs, 3));
			// parCT.push_back(-3 * HesseWeinberg(2, 2) + HesseWeinberg(4, 4) + (vs * HesseWeinberg(4, 5)) / vh);
			// parCT.push_back(HesseWeinberg(4, 5) / (vh * vs));
			// parCT.push_back((-HesseWeinberg(2, 2) + HesseWeinberg(4, 4)) / std::pow(vh, 2));
			// parCT.push_back(-(vh * HesseWeinberg(2, 2)) + NablaWeinberg[4]);
			// parCT.push_back(0);
			parCT.push_back(-HesseWeinberg(3, 3) + NablaWeinberg(5) / vs);
			parCT.push_back((vs * HesseWeinberg(3, 3) - vh * HesseWeinberg(4, 5) - vs * HesseWeinberg(5, 5) + 2 * NablaWeinberg(5)) / vs);
			parCT.push_back((-(vs * HesseWeinberg(5, 5)) + NablaWeinberg(5)) / std::pow(vs, 3));
			parCT.push_back(3 * HesseWeinberg(2, 2) - HesseWeinberg(4, 4) - (vs * HesseWeinberg(4, 5)) / vh);
			parCT.push_back(-(HesseWeinberg(4, 5) / (vh * vs)));
			parCT.push_back((HesseWeinberg(2, 2) - HesseWeinberg(4, 4)) / std::pow(vh, 2));
			parCT.push_back(vh * HesseWeinberg(2, 2) - NablaWeinberg(4));
			parCT.push_back(0);
			return parCT;
		}

		void Class_PNGDM::TripleHiggsCouplings()
		{
			if (!SetCurvatureDone)
				SetCurvatureArrays();
			if (!CalcCouplingsdone)
				CalculatePhysicalCouplings();

			std::vector<double> HiggsOrder(NHiggs);
			// Here you have to set the vector HiggsOrder. By telling e.g. HiggsOrder[0] = 5 you always want your 6th lightest
			// particle to be the first particle in the vector (which has the index 5 because they are sorted by mass)

			// example for keeping the mass order
			for (size_t i = 0; i < NHiggs; i++)
			{
				HiggsOrder[i] = i;
			}

			std::vector<double> TripleDeriv;
			TripleDeriv = WeinbergThirdDerivative();
			std::vector<std::vector<std::vector<double>>> GaugeBasis(NHiggs, std::vector<std::vector<double>>(NHiggs,
																											  std::vector<double>(NHiggs)));
			for (size_t i = 0; i < NHiggs; i++)
			{
				for (size_t j = 0; j < NHiggs; j++)
				{
					for (size_t k = 0; k < NHiggs; k++)
					{
						GaugeBasis[i][j][k] = TripleDeriv.at(i + j * NHiggs + k * NHiggs * NHiggs);
					}
				}
			}

			MatrixXd HiggsRot(NHiggs, NHiggs);
			for (size_t i = 0; i < NHiggs; i++)
			{
				for (size_t j = 0; j < NHiggs; j++)
				{
					HiggsRot(i, j) = HiggsRotationMatrix[i][j];
				}
			}

			MatrixXd HiggsRotSort(NHiggs, NHiggs);

			for (size_t i = 0; i < NHiggs; i++)
			{
				HiggsRotSort.row(i) = HiggsRot.row(HiggsOrder[i]);
			}

			TripleHiggsCorrectionsCWPhysical.resize(NHiggs);
			TripleHiggsCorrectionsTreePhysical.resize(NHiggs);
			TripleHiggsCorrectionsCTPhysical.resize(NHiggs);
			for (size_t i = 0; i < NHiggs; i++)
			{
				TripleHiggsCorrectionsTreePhysical[i].resize(NHiggs);
				TripleHiggsCorrectionsCWPhysical[i].resize(NHiggs);
				TripleHiggsCorrectionsCTPhysical[i].resize(NHiggs);
				for (size_t j = 0; j < NHiggs; j++)
				{
					TripleHiggsCorrectionsCWPhysical[i][j].resize(NHiggs);
					TripleHiggsCorrectionsTreePhysical[i][j].resize(NHiggs);
					TripleHiggsCorrectionsCTPhysical[i][j].resize(NHiggs);
				}
			}

			for (size_t i = 0; i < NHiggs; i++)
			{
				for (size_t j = 0; j < NHiggs; j++)
				{
					for (size_t k = 0; k < NHiggs; k++)
					{
						TripleHiggsCorrectionsCWPhysical[i][j][k] = 0;
						TripleHiggsCorrectionsTreePhysical[i][j][k] = 0;
						TripleHiggsCorrectionsCTPhysical[i][j][k] = 0;
						for (size_t l = 0; l < NHiggs; l++)
						{
							for (size_t m = 0; m < NHiggs; m++)
							{
								for (size_t n = 0; n < NHiggs; n++)
								{
									(void)NHiggs;
									// double RotFac = HiggsRotSort(i, l) * HiggsRotSort(j, m) * HiggsRotSort(k, n);
									// TripleHiggsCorrectionsCWPhysical[i][j][k] += RotFac * GaugeBasis[l][m][n];
									// TripleHiggsCorrectionsTreePhysical[i][j][k] += RotFac * LamHiggs_3[l][m][n];
									// TripleHiggsCorrectionsCTPhysical[i][j][k] += RotFac * LamHiggs_3_CT[l][m][n];
								}
							}
						}
					}
				}
			}
		}

		void Class_PNGDM::SetCurvatureArrays()
		{
			/*
   *  Here you have to set the vectors
   *  Curvature_Higgs_L1,Curvature_Higgs_L2,Curvature_Higgs_L3,Curvature_Higgs_L4
   *  Curvature_Gauge_G2H2
   *  Curvature_Quark_F2H1, Curvature_Lepton_F2H1
   *  as described in the potential in the paper.
   */
			initVectors();
			for (size_t i = 0; i < NHiggs; i++)
				HiggsVev[i] = vevTree[i];

			Curvature_Higgs_L2[0][0] = -muHsqr / 2.;
			Curvature_Higgs_L2[1][1] = -muHsqr / 2.;
			Curvature_Higgs_L2[2][2] = -muHsqr / 2.;
			Curvature_Higgs_L2[3][3] = -muSsqr / 2. + mXsq / 2.;
			Curvature_Higgs_L2[4][4] = -muHsqr / 2.;
			Curvature_Higgs_L2[5][5] = -muSsqr / 2. - mXsq / 2.;

			Curvature_Higgs_L4[0][0][0][0] = 3 * LamH;
			Curvature_Higgs_L4[0][0][1][1] = LamH;
			Curvature_Higgs_L4[0][0][2][2] = LamH;
			Curvature_Higgs_L4[0][0][3][3] = LamHS;
			Curvature_Higgs_L4[0][0][4][4] = LamH;
			Curvature_Higgs_L4[0][0][5][5] = LamHS;
			Curvature_Higgs_L4[0][1][0][1] = LamH;
			Curvature_Higgs_L4[0][1][1][0] = LamH;
			Curvature_Higgs_L4[0][2][0][2] = LamH;
			Curvature_Higgs_L4[0][2][2][0] = LamH;
			Curvature_Higgs_L4[0][3][0][3] = LamHS;
			Curvature_Higgs_L4[0][3][3][0] = LamHS;
			Curvature_Higgs_L4[0][4][0][4] = LamH;
			Curvature_Higgs_L4[0][4][4][0] = LamH;
			Curvature_Higgs_L4[0][5][0][5] = LamHS;
			Curvature_Higgs_L4[0][5][5][0] = LamHS;
			Curvature_Higgs_L4[1][0][0][1] = LamH;
			Curvature_Higgs_L4[1][0][1][0] = LamH;
			Curvature_Higgs_L4[1][1][0][0] = LamH;
			Curvature_Higgs_L4[1][1][1][1] = 3 * LamH;
			Curvature_Higgs_L4[1][1][2][2] = LamH;
			Curvature_Higgs_L4[1][1][3][3] = LamHS;
			Curvature_Higgs_L4[1][1][4][4] = LamH;
			Curvature_Higgs_L4[1][1][5][5] = LamHS;
			Curvature_Higgs_L4[1][2][1][2] = LamH;
			Curvature_Higgs_L4[1][2][2][1] = LamH;
			Curvature_Higgs_L4[1][3][1][3] = LamHS;
			Curvature_Higgs_L4[1][3][3][1] = LamHS;
			Curvature_Higgs_L4[1][4][1][4] = LamH;
			Curvature_Higgs_L4[1][4][4][1] = LamH;
			Curvature_Higgs_L4[1][5][1][5] = LamHS;
			Curvature_Higgs_L4[1][5][5][1] = LamHS;
			Curvature_Higgs_L4[2][0][0][2] = LamH;
			Curvature_Higgs_L4[2][0][2][0] = LamH;
			Curvature_Higgs_L4[2][1][1][2] = LamH;
			Curvature_Higgs_L4[2][1][2][1] = LamH;
			Curvature_Higgs_L4[2][2][0][0] = LamH;
			Curvature_Higgs_L4[2][2][1][1] = LamH;
			Curvature_Higgs_L4[2][2][2][2] = 3 * LamH;
			Curvature_Higgs_L4[2][2][3][3] = LamHS;
			Curvature_Higgs_L4[2][2][4][4] = LamH;
			Curvature_Higgs_L4[2][2][5][5] = LamHS;
			Curvature_Higgs_L4[2][3][2][3] = LamHS;
			Curvature_Higgs_L4[2][3][3][2] = LamHS;
			Curvature_Higgs_L4[2][4][2][4] = LamH;
			Curvature_Higgs_L4[2][4][4][2] = LamH;
			Curvature_Higgs_L4[2][5][2][5] = LamHS;
			Curvature_Higgs_L4[2][5][5][2] = LamHS;
			Curvature_Higgs_L4[3][0][0][3] = LamHS;
			Curvature_Higgs_L4[3][0][3][0] = LamHS;
			Curvature_Higgs_L4[3][1][1][3] = LamHS;
			Curvature_Higgs_L4[3][1][3][1] = LamHS;
			Curvature_Higgs_L4[3][2][2][3] = LamHS;
			Curvature_Higgs_L4[3][2][3][2] = LamHS;
			Curvature_Higgs_L4[3][3][0][0] = LamHS;
			Curvature_Higgs_L4[3][3][1][1] = LamHS;
			Curvature_Higgs_L4[3][3][2][2] = LamHS;
			Curvature_Higgs_L4[3][3][3][3] = 3 * LamS;
			Curvature_Higgs_L4[3][3][4][4] = LamHS;
			Curvature_Higgs_L4[3][3][5][5] = LamS;
			Curvature_Higgs_L4[3][4][3][4] = LamHS;
			Curvature_Higgs_L4[3][4][4][3] = LamHS;
			Curvature_Higgs_L4[3][5][3][5] = LamS;
			Curvature_Higgs_L4[3][5][5][3] = LamS;
			Curvature_Higgs_L4[4][0][0][4] = LamH;
			Curvature_Higgs_L4[4][0][4][0] = LamH;
			Curvature_Higgs_L4[4][1][1][4] = LamH;
			Curvature_Higgs_L4[4][1][4][1] = LamH;
			Curvature_Higgs_L4[4][2][2][4] = LamH;
			Curvature_Higgs_L4[4][2][4][2] = LamH;
			Curvature_Higgs_L4[4][3][3][4] = LamHS;
			Curvature_Higgs_L4[4][3][4][3] = LamHS;
			Curvature_Higgs_L4[4][4][0][0] = LamH;
			Curvature_Higgs_L4[4][4][1][1] = LamH;
			Curvature_Higgs_L4[4][4][2][2] = LamH;
			Curvature_Higgs_L4[4][4][3][3] = LamHS;
			Curvature_Higgs_L4[4][4][4][4] = 3 * LamH;
			Curvature_Higgs_L4[4][4][5][5] = LamHS;
			Curvature_Higgs_L4[4][5][4][5] = LamHS;
			Curvature_Higgs_L4[4][5][5][4] = LamHS;
			Curvature_Higgs_L4[5][0][0][5] = LamHS;
			Curvature_Higgs_L4[5][0][5][0] = LamHS;
			Curvature_Higgs_L4[5][1][1][5] = LamHS;
			Curvature_Higgs_L4[5][1][5][1] = LamHS;
			Curvature_Higgs_L4[5][2][2][5] = LamHS;
			Curvature_Higgs_L4[5][2][5][2] = LamHS;
			Curvature_Higgs_L4[5][3][3][5] = LamS;
			Curvature_Higgs_L4[5][3][5][3] = LamS;
			Curvature_Higgs_L4[5][4][4][5] = LamHS;
			Curvature_Higgs_L4[5][4][5][4] = LamHS;
			Curvature_Higgs_L4[5][5][0][0] = LamHS;
			Curvature_Higgs_L4[5][5][1][1] = LamHS;
			Curvature_Higgs_L4[5][5][2][2] = LamHS;
			Curvature_Higgs_L4[5][5][3][3] = LamS;
			Curvature_Higgs_L4[5][5][4][4] = LamHS;
			Curvature_Higgs_L4[5][5][5][5] = 3 * LamS;
			/*
				* Copied from CxSM Model
			*/
			Curvature_Gauge_G2H2[0][0][0][0] = std::pow(C_g, 2) / 2.;
			Curvature_Gauge_G2H2[0][0][1][1] = std::pow(C_g, 2) / 2.;
			Curvature_Gauge_G2H2[0][0][2][2] = std::pow(C_g, 2) / 2.;
			Curvature_Gauge_G2H2[0][0][4][4] = std::pow(C_g, 2) / 2.;
			Curvature_Gauge_G2H2[0][3][0][4] = (C_gs * C_g) / 2.;
			Curvature_Gauge_G2H2[0][3][1][2] = (C_gs * C_g) / 2.;
			Curvature_Gauge_G2H2[0][3][2][1] = (C_gs * C_g) / 2.;
			Curvature_Gauge_G2H2[0][3][4][0] = (C_gs * C_g) / 2.;
			Curvature_Gauge_G2H2[1][1][0][0] = std::pow(C_g, 2) / 2.;
			Curvature_Gauge_G2H2[1][1][1][1] = std::pow(C_g, 2) / 2.;
			Curvature_Gauge_G2H2[1][1][2][2] = std::pow(C_g, 2) / 2.;
			Curvature_Gauge_G2H2[1][1][4][4] = std::pow(C_g, 2) / 2.;
			Curvature_Gauge_G2H2[1][3][0][2] = (C_gs * C_g) / 2.;
			Curvature_Gauge_G2H2[1][3][1][4] = -(C_gs * C_g) / 2.;
			Curvature_Gauge_G2H2[1][3][2][0] = (C_gs * C_g) / 2.;
			Curvature_Gauge_G2H2[1][3][4][1] = -(C_gs * C_g) / 2.;
			Curvature_Gauge_G2H2[2][2][0][0] = std::pow(C_g, 2) / 2.;
			Curvature_Gauge_G2H2[2][2][1][1] = std::pow(C_g, 2) / 2.;
			Curvature_Gauge_G2H2[2][2][2][2] = std::pow(C_g, 2) / 2.;
			Curvature_Gauge_G2H2[2][2][4][4] = std::pow(C_g, 2) / 2.;
			Curvature_Gauge_G2H2[2][3][0][0] = (C_gs * C_g) / 2.;
			Curvature_Gauge_G2H2[2][3][1][1] = (C_gs * C_g) / 2.;
			Curvature_Gauge_G2H2[2][3][2][2] = -(C_gs * C_g) / 2.;
			Curvature_Gauge_G2H2[2][3][4][4] = -(C_gs * C_g) / 2.;
			Curvature_Gauge_G2H2[3][0][0][4] = (C_gs * C_g) / 2.;
			Curvature_Gauge_G2H2[3][0][1][2] = (C_gs * C_g) / 2.;
			Curvature_Gauge_G2H2[3][0][2][1] = (C_gs * C_g) / 2.;
			Curvature_Gauge_G2H2[3][0][4][0] = (C_gs * C_g) / 2.;
			Curvature_Gauge_G2H2[3][1][0][2] = (C_gs * C_g) / 2.;
			Curvature_Gauge_G2H2[3][1][1][4] = -(C_gs * C_g) / 2.;
			Curvature_Gauge_G2H2[3][1][2][0] = (C_gs * C_g) / 2.;
			Curvature_Gauge_G2H2[3][1][4][1] = -(C_gs * C_g) / 2.;
			Curvature_Gauge_G2H2[3][2][0][0] = (C_gs * C_g) / 2.;
			Curvature_Gauge_G2H2[3][2][1][1] = (C_gs * C_g) / 2.;
			Curvature_Gauge_G2H2[3][2][2][2] = -(C_gs * C_g) / 2.;
			Curvature_Gauge_G2H2[3][2][4][4] = -(C_gs * C_g) / 2.;
			Curvature_Gauge_G2H2[3][3][0][0] = std::pow(C_gs, 2) / 2.;
			Curvature_Gauge_G2H2[3][3][1][1] = std::pow(C_gs, 2) / 2.;
			Curvature_Gauge_G2H2[3][3][2][2] = std::pow(C_gs, 2) / 2.;
			Curvature_Gauge_G2H2[3][3][4][4] = std::pow(C_gs, 2) / 2.;

			MatrixXcd YIJQc1(NQuarks, NQuarks), YIJQc2(NQuarks, NQuarks), YIJQc2OI(NQuarks, NQuarks), YIJQg0(NQuarks, NQuarks),
				YIJQg0OI(NQuarks, NQuarks),
				YIJQh1(NQuarks, NQuarks), YIJQh2(NQuarks, NQuarks), YIJQh3(NQuarks, NQuarks);

			std::complex<double> V11, V12, V13, V21, V22, V23, V31, V32, V33;
			V11 = C_Vud;
			V12 = C_Vus;
			V13 = C_Vub;
			V21 = C_Vcd;
			V22 = C_Vcs;
			V23 = C_Vcb;
			V31 = C_Vtd;
			V32 = C_Vts;
			V33 = C_Vtb;

			std::complex<double> VC11, VC12, VC13, VC21, VC22, VC23, VC31, VC32, VC33;
			VC11 = std::conj(C_Vud);
			VC12 = std::conj(C_Vus);
			VC13 = std::conj(C_Vub);
			VC21 = std::conj(C_Vcd);
			VC22 = std::conj(C_Vcs);
			VC23 = std::conj(C_Vcb);
			VC31 = std::conj(C_Vtd);
			VC32 = std::conj(C_Vts);
			VC33 = std::conj(C_Vtb);

			std::complex<double> II(0, 1);

			Curvature_Lepton_F2H1[0][1][3] = II / vh * C_MassElectron;
			Curvature_Lepton_F2H1[0][1][4] = 0.1e1 / vh * C_MassElectron;
			Curvature_Lepton_F2H1[1][0][2] = II / vh * C_MassElectron;
			Curvature_Lepton_F2H1[1][0][4] = 0.1e1 / vh * C_MassElectron;
			Curvature_Lepton_F2H1[1][6][0] = 0.1e1 / vh * C_MassElectron;
			Curvature_Lepton_F2H1[1][6][1] = II / vh * C_MassElectron;
			Curvature_Lepton_F2H1[2][3][2] = II / vh * C_MassMu;
			Curvature_Lepton_F2H1[2][3][4] = 0.1e1 / vh * C_MassMu;
			Curvature_Lepton_F2H1[3][2][2] = II / vh * C_MassMu;
			Curvature_Lepton_F2H1[3][2][4] = 0.1e1 / vh * C_MassMu;
			Curvature_Lepton_F2H1[3][7][0] = 0.1e1 / vh * C_MassMu;
			Curvature_Lepton_F2H1[3][7][1] = II / vh * C_MassMu;
			Curvature_Lepton_F2H1[4][5][2] = II / vh * C_MassTau;
			Curvature_Lepton_F2H1[4][5][4] = 0.1e1 / vh * C_MassTau;
			Curvature_Lepton_F2H1[5][4][2] = II / vh * C_MassTau;
			Curvature_Lepton_F2H1[5][4][4] = 0.1e1 / vh * C_MassTau;
			Curvature_Lepton_F2H1[5][8][0] = 0.1e1 / vh * C_MassTau;
			Curvature_Lepton_F2H1[5][8][1] = II / vh * C_MassTau;
			Curvature_Lepton_F2H1[6][1][0] = 0.1e1 / vh * C_MassElectron;
			Curvature_Lepton_F2H1[6][1][1] = II / vh * C_MassElectron;
			Curvature_Lepton_F2H1[7][3][0] = 0.1e1 / vh * C_MassMu;
			Curvature_Lepton_F2H1[7][3][1] = II / vh * C_MassMu;
			Curvature_Lepton_F2H1[8][5][0] = 0.1e1 / vh * C_MassTau;
			Curvature_Lepton_F2H1[8][5][1] = II / vh * C_MassTau;

			Curvature_Quark_F2H1[0][6][2] = -II / vh * C_MassUp;
			Curvature_Quark_F2H1[0][6][4] = 0.1e1 / vh * C_MassUp;
			Curvature_Quark_F2H1[0][9][0] = -0.1e1 / vh * C_MassUp * conj(V11);
			Curvature_Quark_F2H1[0][9][1] = II / vh * C_MassUp * conj(V11);
			Curvature_Quark_F2H1[0][10][0] = -0.1e1 / vh * C_MassUp * conj(V12);
			Curvature_Quark_F2H1[0][10][1] = II / vh * C_MassUp * conj(V12);
			Curvature_Quark_F2H1[0][11][0] = -0.1e1 / vh * C_MassUp * conj(V13);
			Curvature_Quark_F2H1[0][11][1] = II / vh * C_MassUp * conj(V13);
			Curvature_Quark_F2H1[1][7][2] = -II / vh * C_MassCharm;
			Curvature_Quark_F2H1[1][7][4] = 0.1e1 / vh * C_MassCharm;
			Curvature_Quark_F2H1[1][9][0] = -0.1e1 / vh * C_MassCharm * conj(V21);
			Curvature_Quark_F2H1[1][9][1] = II / vh * C_MassCharm * conj(V21);
			Curvature_Quark_F2H1[1][10][0] = -0.1e1 / vh * C_MassCharm * conj(V22);
			Curvature_Quark_F2H1[1][10][1] = II / vh * C_MassCharm * conj(V22);
			Curvature_Quark_F2H1[1][11][0] = -0.1e1 / vh * C_MassCharm * conj(V23);
			Curvature_Quark_F2H1[1][11][1] = II / vh * C_MassCharm * conj(V23);
			Curvature_Quark_F2H1[2][8][2] = -II / vh * C_MassTop;
			Curvature_Quark_F2H1[2][8][4] = 0.1e1 / vh * C_MassTop;
			Curvature_Quark_F2H1[2][9][0] = -0.1e1 / vh * C_MassTop * conj(V31);
			Curvature_Quark_F2H1[2][9][1] = II / vh * C_MassTop * conj(V31);
			Curvature_Quark_F2H1[2][10][0] = -0.1e1 / vh * C_MassTop * conj(V32);
			Curvature_Quark_F2H1[2][10][1] = II / vh * C_MassTop * conj(V32);
			Curvature_Quark_F2H1[2][11][0] = -0.1e1 / vh * C_MassTop * conj(V33);
			Curvature_Quark_F2H1[2][11][1] = II / vh * C_MassTop * conj(V33);
			Curvature_Quark_F2H1[3][6][0] = 0.1e1 / vh * C_MassDown * V11;
			Curvature_Quark_F2H1[3][6][1] = II / vh * C_MassDown * V11;
			Curvature_Quark_F2H1[3][7][0] = V21 / vh * C_MassDown;
			Curvature_Quark_F2H1[3][7][1] = II * V21 / vh * C_MassDown;
			Curvature_Quark_F2H1[3][8][0] = 0.1e1 / vh * C_MassDown * V31;
			Curvature_Quark_F2H1[3][8][1] = II / vh * C_MassDown * V31;
			Curvature_Quark_F2H1[3][9][2] = II / vh * C_MassDown;
			Curvature_Quark_F2H1[3][9][4] = 0.1e1 / vh * C_MassDown;
			Curvature_Quark_F2H1[4][6][0] = 0.1e1 / vh * C_MassStrange * V12;
			Curvature_Quark_F2H1[4][6][1] = II / vh * C_MassStrange * V12;
			Curvature_Quark_F2H1[4][7][0] = V22 / vh * C_MassStrange;
			Curvature_Quark_F2H1[4][7][1] = II * V22 / vh * C_MassStrange;
			Curvature_Quark_F2H1[4][8][0] = 0.1e1 / vh * C_MassStrange * V32;
			Curvature_Quark_F2H1[4][8][1] = II / vh * C_MassStrange * V32;
			Curvature_Quark_F2H1[4][10][2] = II / vh * C_MassStrange;
			Curvature_Quark_F2H1[4][10][4] = 0.1e1 / vh * C_MassStrange;
			Curvature_Quark_F2H1[5][6][0] = V13 / vh * C_MassBottom;
			Curvature_Quark_F2H1[5][6][1] = II / vh * C_MassBottom * V13;
			Curvature_Quark_F2H1[5][7][0] = V23 / vh * C_MassBottom;
			Curvature_Quark_F2H1[5][7][1] = II / vh * C_MassBottom * V23;
			Curvature_Quark_F2H1[5][8][0] = V33 / vh * C_MassBottom;
			Curvature_Quark_F2H1[5][8][1] = II / vh * C_MassBottom * V33;
			Curvature_Quark_F2H1[5][11][2] = II / vh * C_MassBottom;
			Curvature_Quark_F2H1[5][11][4] = 0.1e1 / vh * C_MassBottom;
			Curvature_Quark_F2H1[6][0][2] = -II / vh * C_MassUp;
			Curvature_Quark_F2H1[6][0][4] = 0.1e1 / vh * C_MassUp;
			Curvature_Quark_F2H1[6][3][0] = 0.1e1 / vh * C_MassDown * V11;
			Curvature_Quark_F2H1[6][3][1] = II / vh * C_MassDown * V11;
			Curvature_Quark_F2H1[6][4][0] = 0.1e1 / vh * C_MassStrange * V12;
			Curvature_Quark_F2H1[6][4][1] = II / vh * C_MassStrange * V12;
			Curvature_Quark_F2H1[6][5][0] = V13 / vh * C_MassBottom;
			Curvature_Quark_F2H1[6][5][1] = II / vh * C_MassBottom * V13;
			Curvature_Quark_F2H1[7][1][2] = -II / vh * C_MassCharm;
			Curvature_Quark_F2H1[7][1][4] = 0.1e1 / vh * C_MassCharm;
			Curvature_Quark_F2H1[7][3][0] = V21 / vh * C_MassDown;
			Curvature_Quark_F2H1[7][3][1] = II * V21 / vh * C_MassDown;
			Curvature_Quark_F2H1[7][4][0] = V22 / vh * C_MassStrange;
			Curvature_Quark_F2H1[7][4][1] = II * V22 / vh * C_MassStrange;
			Curvature_Quark_F2H1[7][5][0] = V23 / vh * C_MassBottom;
			Curvature_Quark_F2H1[7][5][1] = II / vh * C_MassBottom * V23;
			Curvature_Quark_F2H1[8][2][2] = -II / vh * C_MassTop;
			Curvature_Quark_F2H1[8][2][4] = 0.1e1 / vh * C_MassTop;
			Curvature_Quark_F2H1[8][3][0] = 0.1e1 / vh * C_MassDown * V31;
			Curvature_Quark_F2H1[8][3][1] = II / vh * C_MassDown * V31;
			Curvature_Quark_F2H1[8][4][0] = 0.1e1 / vh * C_MassStrange * V32;
			Curvature_Quark_F2H1[8][4][1] = II / vh * C_MassStrange * V32;
			Curvature_Quark_F2H1[8][5][0] = V33 / vh * C_MassBottom;
			Curvature_Quark_F2H1[8][5][1] = II / vh * C_MassBottom * V33;
			Curvature_Quark_F2H1[9][0][0] = -0.1e1 / vh * C_MassUp * conj(V11);
			Curvature_Quark_F2H1[9][0][1] = II / vh * C_MassUp * conj(V11);
			Curvature_Quark_F2H1[9][1][0] = -0.1e1 / vh * C_MassCharm * conj(V21);
			Curvature_Quark_F2H1[9][1][1] = II / vh * C_MassCharm * conj(V21);
			Curvature_Quark_F2H1[9][2][0] = -0.1e1 / vh * C_MassTop * conj(V31);
			Curvature_Quark_F2H1[9][2][1] = II / vh * C_MassTop * conj(V31);
			Curvature_Quark_F2H1[9][3][2] = II / vh * C_MassDown;
			Curvature_Quark_F2H1[9][3][4] = 0.1e1 / vh * C_MassDown;
			Curvature_Quark_F2H1[10][0][0] = -0.1e1 / vh * C_MassUp * conj(V12);
			Curvature_Quark_F2H1[10][0][1] = II / vh * C_MassUp * conj(V12);
			Curvature_Quark_F2H1[10][1][0] = -0.1e1 / vh * C_MassCharm * conj(V22);
			Curvature_Quark_F2H1[10][1][1] = II / vh * C_MassCharm * conj(V22);
			Curvature_Quark_F2H1[10][2][0] = -0.1e1 / vh * C_MassTop * conj(V32);
			Curvature_Quark_F2H1[10][2][1] = II / vh * C_MassTop * conj(V32);
			Curvature_Quark_F2H1[10][4][2] = II / vh * C_MassStrange;
			Curvature_Quark_F2H1[10][4][4] = 0.1e1 / vh * C_MassStrange;
			Curvature_Quark_F2H1[11][0][0] = -0.1e1 / vh * C_MassUp * conj(V13);
			Curvature_Quark_F2H1[11][0][1] = II / vh * C_MassUp * conj(V13);
			Curvature_Quark_F2H1[11][1][0] = -0.1e1 / vh * C_MassCharm * conj(V23);
			Curvature_Quark_F2H1[11][1][1] = II / vh * C_MassCharm * conj(V23);
			Curvature_Quark_F2H1[11][2][0] = -0.1e1 / vh * C_MassTop * conj(V33);
			Curvature_Quark_F2H1[11][2][1] = II / vh * C_MassTop * conj(V33);
			Curvature_Quark_F2H1[11][5][2] = II / vh * C_MassBottom;
			Curvature_Quark_F2H1[11][5][4] = 0.1e1 / vh * C_MassBottom;

			SetCurvatureDone = true;
		}

		bool Class_PNGDM::CalculateDebyeSimplified()
		{
			return false;
			/*
   * Use this function if you calculated the Debye corrections to the Higgs mass matrix and implement
   * your formula here and return true. The vector is given by DebyeHiggs[NHiggs][NHiggs]
   */
		}

		bool Class_PNGDM::CalculateDebyeGaugeSimplified()
		{
			/*
     * Use this function if you calculated the Debye corrections to the gauge mass matrix and implement
     * your formula here and return true. The vector is given by DebyeGauge[NGauge][NGauge]
     */

			return false;
		}
		double Class_PNGDM::VTreeSimplified(const std::vector<double> &v) const
		{
			(void)v;
			if (not UseVTreeSimplified)
				return 0;
			double res = 0;
			return res;
		}

		double Class_PNGDM::VCounterSimplified(const std::vector<double> &v) const
		{
			(void)v;
			if (not UseVCounterSimplified)
				return 0;
			double res = 0;
			return res;
		}

		void Class_PNGDM::Debugging(const std::vector<double> &input, std::vector<double> &output) const
		{
			(void)input;
			(void)output;
		}

	} // namespace Models
} // namespace BSMPT
