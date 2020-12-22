/*
 * ClassTemplate.cpp
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

#include <ext/alloc_traits.h>               // for __alloc_traits<>::value_type
#include <stddef.h>                         // for std::size_t
#include <algorithm>                        // for max, copy
#include <iostream>                         // for operator<<, endl, basic_o...
#include <memory>                           // for allocator_traits<>::value...
#include <BSMPT/models/SMparam.h>           // for C_vev0, C_MassTop, C_g
#include "Eigen/Dense"
#include "Eigen/Eigenvalues"
#include "Eigen/IterativeLinearSolvers"

#include <BSMPT/models/future/ClassPotentialSM.h>

#include <BSMPT/models/IncludeAllModels.h>

#include <BSMPT/utility.h>
using namespace Eigen;

/**
 * @file
 * Template for adding a new model class
 */

namespace BSMPT {
namespace Models {


/**
 * Here you have to adjust NNeutralHiggs, NChargedHiggs, nPar (number of Lagrangian parameters AFTER
 *  using the tadpole conditions),
 * nParCT (number of counterterms) as well as nVEV (number of VEVs for minimization)
 */
Class_SM::Class_SM ()
{
  Model = ModelID::ModelIDs::SM; // global int constant which will be used to tell the program which model is called
  NNeutralHiggs = 4; // number of neutral Higgs bosons at T = 0
  NChargedHiggs=0; // number of charged Higgs bosons  at T = 0 (all d.o.f.)

  nPar = 2; // number of parameters in the tree-Level Lagrangian
  nParCT = 2+4; // number of parameters in the counterterm potential

  nVEV=1; // number of VEVs to minimize the potential

  NHiggs = NNeutralHiggs+NChargedHiggs;

  VevOrder.resize(nVEV);
  // Here you have to tell which scalar field gets which VEV.
  VevOrder[0] = 3;

  // Set UseVTreeSimplified to use the tree-level potential defined in VTreeSimplified
  UseVTreeSimplified = false;

  // Set UseVCounterSimplified to use the counterterm potential defined in VCounterSimplified
  UseVCounterSimplified = false;

}

Class_SM::~Class_SM ()
{
}

/**
 * returns a string which tells the user the chronological order of the counterterms. Use this to
 * complement the legend of the given input file
 */
std::vector<std::string> Class_SM::addLegendCT() const
{
    std::vector<std::string> labels;
    labels.push_back("dlambda");
    labels.push_back("dmsquared");
    labels.push_back("dT1");
    labels.push_back("dT2");
    labels.push_back("dT3");
    labels.push_back("dT4");
    return labels;
}

/**
 * returns a string which tells the user the chronological order of the VEVs and the critical temperature. Use this to
 * complement the legend of the given input file
 */
std::vector<std::string> Class_SM::addLegendTemp() const
{
    std::vector<std::string> labels;
    labels.push_back("T_c");
    labels.push_back("v_c");
    labels.push_back("omega_c/T_c");
    labels.push_back("omega");
    return  labels;
}

/**
 * returns a string which tells the user the chronological order of the Triple Higgs couplings. Use this to
 * complement the legend of the given input file
 *
 */
std::vector<std::string> Class_SM::addLegendTripleCouplings() const
{
    std::vector<std::string> labels;
	std::vector<std::string> particles;
	particles.resize(NHiggs);
	//here you have to define the particle names in the vector particles

	particles[0]="GP";
	particles[1]="GM";
	particles[2]="H";
	particles[3]="G0";

    for(size_t i=0;i<NHiggs;i++)
    {
        for(size_t j=i;j<NHiggs;j++)
        {
            for(size_t k=j;k<NHiggs;k++)
            {
                labels.push_back("Tree_"+particles.at(i)+particles.at(j)+particles.at(k));
                labels.push_back("CT_"+particles.at(i)+particles.at(j)+particles.at(k));
                labels.push_back("CW_"+particles.at(i)+particles.at(j)+particles.at(k));
            }
        }
    }

    return labels;
}

/**
 * returns a string which tells the user the chronological order of the VEVs. Use this to
 * complement the legend of the given input file
 */
std::vector<std::string> Class_SM::addLegendVEV() const
{
    std::vector<std::string> labels;
    labels.push_back("omega");
    return  labels;
}

/**
 * Reads the string linestr and sets the parameter point
 */
void Class_SM::ReadAndSet(const std::string& linestr, std::vector<double>& par )
{
	std::stringstream ss(linestr);
	double tmp;
	if (UseIndexCol){
		ss >> tmp;
	}

	double mh = 0;
	ss >> mh;

	lambda = 2.0*std::pow(mh/C_vev0,2);
	mus = -mh*mh;

	par[0] = mus;
	par[1] = lambda;


	set_gen(par); // This you have to call so that everything will be set
	return ;
}


/**
 * Set Class Object as well as the VEV configuration
 */
void Class_SM::set_gen(const std::vector<double>& par) {


	mus = par[0];
	lambda = par[1];

	scale = C_vev0;

	vevTreeMin.resize(nVEV);
	vevTree.resize(NHiggs);

	// Here you have to set the vector vevTreeMin. The vector vevTree will then be set by the function MinimizeOrderVEV
	vevTreeMin[0] = C_vev0;

    vevTree=MinimizeOrderVEV(vevTreeMin);
	if(!SetCurvatureDone) SetCurvatureArrays();
}

/**
 * set your counterterm parameters from the entries of par as well as the entries of Curvature_Higgs_CT_L1 to
 * Curvature_Higgs_CT_L4.
 */
void Class_SM::set_CT_Pot_Par(const std::vector<double>& par){

	dlambda = par[0];
	dmus = par[1];
	dT1 = par[2];
	dT2 = par[3];
	dT3 = par[4];
	dT4 = par[5];

	Curvature_Higgs_CT_L1[0] = dT1;
	Curvature_Higgs_CT_L1[1] = dT2;
	Curvature_Higgs_CT_L1[2] = dT3;
	Curvature_Higgs_CT_L1[3] = dT4;

	Curvature_Higgs_CT_L2[0][0] = dmus / 0.2e1;
	Curvature_Higgs_CT_L2[1][1] = dmus / 0.2e1;
	Curvature_Higgs_CT_L2[2][2] = dmus / 0.2e1;
	Curvature_Higgs_CT_L2[3][3] = dmus / 0.2e1;

	Curvature_Higgs_CT_L4[0][0][0][0] = 0.3e1 / 0.2e1 * dlambda;
	Curvature_Higgs_CT_L4[0][0][1][1] = dlambda / 0.2e1;
	Curvature_Higgs_CT_L4[0][0][2][2] = dlambda / 0.2e1;
	Curvature_Higgs_CT_L4[0][0][3][3] = dlambda / 0.2e1;
	Curvature_Higgs_CT_L4[0][1][0][1] = dlambda / 0.2e1;
	Curvature_Higgs_CT_L4[0][1][1][0] = dlambda / 0.2e1;
	Curvature_Higgs_CT_L4[0][2][0][2] = dlambda / 0.2e1;
	Curvature_Higgs_CT_L4[0][2][2][0] = dlambda / 0.2e1;
	Curvature_Higgs_CT_L4[0][3][0][3] = dlambda / 0.2e1;
	Curvature_Higgs_CT_L4[0][3][3][0] = dlambda / 0.2e1;
	Curvature_Higgs_CT_L4[1][0][0][1] = dlambda / 0.2e1;
	Curvature_Higgs_CT_L4[1][0][1][0] = dlambda / 0.2e1;
	Curvature_Higgs_CT_L4[1][1][0][0] = dlambda / 0.2e1;
	Curvature_Higgs_CT_L4[1][1][1][1] = 0.3e1 / 0.2e1 * dlambda;
	Curvature_Higgs_CT_L4[1][1][2][2] = dlambda / 0.2e1;
	Curvature_Higgs_CT_L4[1][1][3][3] = dlambda / 0.2e1;
	Curvature_Higgs_CT_L4[1][2][1][2] = dlambda / 0.2e1;
	Curvature_Higgs_CT_L4[1][2][2][1] = dlambda / 0.2e1;
	Curvature_Higgs_CT_L4[1][3][1][3] = dlambda / 0.2e1;
	Curvature_Higgs_CT_L4[1][3][3][1] = dlambda / 0.2e1;
	Curvature_Higgs_CT_L4[2][0][0][2] = dlambda / 0.2e1;
	Curvature_Higgs_CT_L4[2][0][2][0] = dlambda / 0.2e1;
	Curvature_Higgs_CT_L4[2][1][1][2] = dlambda / 0.2e1;
	Curvature_Higgs_CT_L4[2][1][2][1] = dlambda / 0.2e1;
	Curvature_Higgs_CT_L4[2][2][0][0] = dlambda / 0.2e1;
	Curvature_Higgs_CT_L4[2][2][1][1] = dlambda / 0.2e1;
	Curvature_Higgs_CT_L4[2][2][2][2] = 0.3e1 / 0.2e1 * dlambda;
	Curvature_Higgs_CT_L4[2][2][3][3] = dlambda / 0.2e1;
	Curvature_Higgs_CT_L4[2][3][2][3] = dlambda / 0.2e1;
	Curvature_Higgs_CT_L4[2][3][3][2] = dlambda / 0.2e1;
	Curvature_Higgs_CT_L4[3][0][0][3] = dlambda / 0.2e1;
	Curvature_Higgs_CT_L4[3][0][3][0] = dlambda / 0.2e1;
	Curvature_Higgs_CT_L4[3][1][1][3] = dlambda / 0.2e1;
	Curvature_Higgs_CT_L4[3][1][3][1] = dlambda / 0.2e1;
	Curvature_Higgs_CT_L4[3][2][2][3] = dlambda / 0.2e1;
	Curvature_Higgs_CT_L4[3][2][3][2] = dlambda / 0.2e1;
	Curvature_Higgs_CT_L4[3][3][0][0] = dlambda / 0.2e1;
	Curvature_Higgs_CT_L4[3][3][1][1] = dlambda / 0.2e1;
	Curvature_Higgs_CT_L4[3][3][2][2] = dlambda / 0.2e1;
	Curvature_Higgs_CT_L4[3][3][3][3] = 0.3e1 / 0.2e1 * dlambda;

	return ;
}


/**
 * console output of all Parameters
 */
void Class_SM::write() const {

	std::cout << "The parameters are : " << std::endl;
	std::cout << "lambda = " << lambda << std::endl
			<< "\tm^2 = " << mus << std::endl;

	std::cout << "The counterterm parameters are : " << std::endl;
	std::cout << "dT1 = "<< dT1 << std::endl
			<< "dT2 = "<< dT2 << std::endl
			<< "dT3 = "<< dT3 << std::endl
			<< "dT4 = "<< dT4 << std::endl
			<< "dlambda = " << dlambda << std::endl
			<< "dm^2 = "<< dmus << std::endl;

	std::cout << "The scale is given by mu = " << scale << " GeV " << std::endl;

	std::vector<double> Higgsmasses;
    Higgsmasses=HiggsMassesSquared(vevTree,0,0);
	std::cout << "The Higgs masses are given by " << std::endl;
	for(auto x: Higgsmasses) std::cout << "m^2 = " << x << "\tm = " << std::sqrt(std::abs(x)) << std::endl;

}


/**
 * Calculates the counterterms. Here you need to work out the scheme and implement the formulas.
 */
std::vector<double> Class_SM::calc_CT() const
{
    std::vector<double> parCT;
	bool Debug=false;
	if(Debug) std::cout << "Start" << __func__ << std::endl;

    if(!SetCurvatureDone){
        std::string retmes = __func__;
        retmes += " was called before SetCurvatureArrays()!\n";
        throw std::runtime_error(retmes);
    }
    if(!CalcCouplingsdone){
        std::string retmes = __func__;
        retmes += " was called before CalculatePhysicalCouplings()!\n";
        throw std::runtime_error(retmes);
    }
	if(Debug) {
	std::cout << "Couplings done " << std::endl;
	}
	std::vector<double> WeinbergNabla,WeinbergHesse;
    WeinbergNabla=WeinbergFirstDerivative();
    WeinbergHesse=WeinbergSecondDerivative();

	if(Debug) std::cout << "Finished Derivatives " << std::endl;

	VectorXd NablaWeinberg(NHiggs);
	MatrixXd HesseWeinberg(NHiggs,NHiggs),HiggsRot(NHiggs,NHiggs);
	for(size_t i=0;i<NHiggs;i++)
	{
		NablaWeinberg[i] = WeinbergNabla[i];
		for(size_t j=0;j<NHiggs;j++) HesseWeinberg(i,j) = WeinbergHesse.at(j*NHiggs+i);
	}

	if(Debug){
		std::cout << "Nabla VCW = " << NablaWeinberg.transpose() << std::endl;
		std::cout << "Hesse VCW = \n" << HesseWeinberg << std::endl;
	}

	// Here you have to use your formulas for the counterterm scheme


    //dlambda
    parCT.push_back((2 * HesseWeinberg(2, 2) - 2 * HesseWeinberg(3, 3)) * pow(C_vev0, -2));
    //dmus
    parCT.push_back(-3 * HesseWeinberg(2, 2) + HesseWeinberg(3, 3));
    //dT1
    parCT.push_back(-NablaWeinberg(0));
    //dT2
    parCT.push_back(-NablaWeinberg(1));
    //dT3
    parCT.push_back(-NablaWeinberg(2));
    //dT4
    parCT.push_back(HesseWeinberg(2, 2) * C_vev0 - NablaWeinberg(3));

    return parCT;

}




void Class_SM::TripleHiggsCouplings()
{
	bool Debug=false;
	if(Debug) std::cout << "Debug turned on in " << __func__ << std::endl;

	if(!SetCurvatureDone)SetCurvatureArrays();
	if(!CalcCouplingsdone)CalculatePhysicalCouplings();


	std::vector<double> HiggsOrder(NHiggs);
	// Here you have to set the vector HiggsOrder. By telling e.g. HiggsOrder[0] = 5 you always want your 6th lightest
	// particle to be the first particle in the vector (which has the index 5 because they are sorted by mass)

	// example for keeping the mass order
	for(size_t i=0;i<NHiggs;i++) {
		HiggsOrder[i]=i;
	}

	if(Debug) std::cout << "Calculate Derivative" << std::endl;
	std::vector<double> TripleDeriv;
    TripleDeriv=WeinbergThirdDerivative();
	if(Debug) std::cout << "Finished calculating derivatives " << std::endl;
	std::vector<std::vector<std::vector<double>>> GaugeBasis(NHiggs, std::vector<std::vector<double>>(NHiggs,
				std::vector<double>(NHiggs)));
	for(size_t i=0;i<NHiggs;i++)
	  {
		for(size_t j=0;j<NHiggs;j++)
		{
		  for(size_t k=0;k<NHiggs;k++)
			{
			  GaugeBasis[i][j][k] = TripleDeriv.at(i+j*NHiggs+k*NHiggs*NHiggs);
			}
		}
	  }

	MatrixXd HiggsRot(NHiggs,NHiggs);
	for(size_t i=0;i<NHiggs;i++)
	{
		for(size_t j=0;j<NHiggs;j++)
		{
			HiggsRot(i,j) = HiggsRotationMatrix[i][j];
		}
	}


	MatrixXd HiggsRotSort(NHiggs,NHiggs);






	for(size_t i=0;i<NHiggs;i++)
	{
		HiggsRotSort.row(i) = HiggsRot.row(HiggsOrder[i]);
	}

	TripleHiggsCorrectionsCWPhysical.resize(NHiggs);
	TripleHiggsCorrectionsTreePhysical.resize(NHiggs);
	TripleHiggsCorrectionsCTPhysical.resize(NHiggs);
	for(size_t i=0;i<NHiggs;i++) {
		TripleHiggsCorrectionsTreePhysical[i].resize(NHiggs);
		TripleHiggsCorrectionsCWPhysical[i].resize(NHiggs);
		TripleHiggsCorrectionsCTPhysical[i].resize(NHiggs);
		for(size_t j=0;j<NHiggs;j++) {
			TripleHiggsCorrectionsCWPhysical[i][j].resize(NHiggs);
			TripleHiggsCorrectionsTreePhysical[i][j].resize(NHiggs);
			TripleHiggsCorrectionsCTPhysical[i][j].resize(NHiggs);
		}
	}

	if(Debug) std::cout << "Setup done " << std::endl;


	for(size_t i=0;i<NHiggs;i++)
	  {
		for(size_t j=0;j<NHiggs;j++)
		{
			for(size_t k=0;k<NHiggs;k++)
			{
			  TripleHiggsCorrectionsCWPhysical[i][j][k] = 0;
			  TripleHiggsCorrectionsTreePhysical[i][j][k] = 0;
			  TripleHiggsCorrectionsCTPhysical[i][j][k] = 0;
			  for(size_t l=0;l<NHiggs;l++)
			  {
				  for(size_t m=0;m<NHiggs;m++)
				  {
					  for(size_t n=0;n<NHiggs;n++)
					  {
						  double RotFac = HiggsRotSort(i,l)*HiggsRotSort(j,m)*HiggsRotSort(k,n);
						  TripleHiggsCorrectionsCWPhysical[i][j][k] += RotFac*GaugeBasis[l][m][n];
						  TripleHiggsCorrectionsTreePhysical[i][j][k] += RotFac*LambdaHiggs_3[l][m][n];
						  TripleHiggsCorrectionsCTPhysical[i][j][k] += RotFac*LambdaHiggs_3_CT[l][m][n];

					  }
				  }
			  }
			}
		}
	  }



}

void Class_SM::SetCurvatureArrays(){
  /*
   *  Here you have to set the vectors
   *  Curvature_Higgs_L1,Curvature_Higgs_L2,Curvature_Higgs_L3,Curvature_Higgs_L4
   *  Curvature_Gauge_G2H2
   *  Curvature_Quark_F2H1, Curvature_Lepton_F2H1
   *  as described in the potential in the paper.
   */

	initVectors();
	SetCurvatureDone=true;
	for(size_t i=0;i<NHiggs;i++) HiggsVev[i] = vevTree[i];


	std::complex<double> II(0,1);

	using std::conj;


	Curvature_Higgs_L2[0][0] = mus / 0.2e1;
	Curvature_Higgs_L2[0][1] = 0;
	Curvature_Higgs_L2[0][2] = 0;
	Curvature_Higgs_L2[0][3] = 0;
	Curvature_Higgs_L2[1][0] = 0;
	Curvature_Higgs_L2[1][1] = mus / 0.2e1;
	Curvature_Higgs_L2[1][2] = 0;
	Curvature_Higgs_L2[1][3] = 0;
	Curvature_Higgs_L2[2][0] = 0;
	Curvature_Higgs_L2[2][1] = 0;
	Curvature_Higgs_L2[2][2] = mus / 0.2e1;
	Curvature_Higgs_L2[2][3] = 0;
	Curvature_Higgs_L2[3][0] = 0;
	Curvature_Higgs_L2[3][1] = 0;
	Curvature_Higgs_L2[3][2] = 0;
	Curvature_Higgs_L2[3][3] = mus / 0.2e1;

	Curvature_Higgs_L4[0][0][0][0] = 0.3e1 / 0.2e1 * lambda;
	Curvature_Higgs_L4[0][0][1][1] = lambda / 0.2e1;
	Curvature_Higgs_L4[0][0][2][2] = lambda / 0.2e1;
	Curvature_Higgs_L4[0][0][3][3] = lambda / 0.2e1;
	Curvature_Higgs_L4[0][1][0][1] = lambda / 0.2e1;
	Curvature_Higgs_L4[0][1][1][0] = lambda / 0.2e1;
	Curvature_Higgs_L4[0][2][0][2] = lambda / 0.2e1;
	Curvature_Higgs_L4[0][2][2][0] = lambda / 0.2e1;
	Curvature_Higgs_L4[0][3][0][3] = lambda / 0.2e1;
	Curvature_Higgs_L4[0][3][3][0] = lambda / 0.2e1;
	Curvature_Higgs_L4[1][0][0][1] = lambda / 0.2e1;
	Curvature_Higgs_L4[1][0][1][0] = lambda / 0.2e1;
	Curvature_Higgs_L4[1][1][0][0] = lambda / 0.2e1;
	Curvature_Higgs_L4[1][1][1][1] = 0.3e1 / 0.2e1 * lambda;
	Curvature_Higgs_L4[1][1][2][2] = lambda / 0.2e1;
	Curvature_Higgs_L4[1][1][3][3] = lambda / 0.2e1;
	Curvature_Higgs_L4[1][2][1][2] = lambda / 0.2e1;
	Curvature_Higgs_L4[1][2][2][1] = lambda / 0.2e1;
	Curvature_Higgs_L4[1][3][1][3] = lambda / 0.2e1;
	Curvature_Higgs_L4[1][3][3][1] = lambda / 0.2e1;
	Curvature_Higgs_L4[2][0][0][2] = lambda / 0.2e1;
	Curvature_Higgs_L4[2][0][2][0] = lambda / 0.2e1;
	Curvature_Higgs_L4[2][1][1][2] = lambda / 0.2e1;
	Curvature_Higgs_L4[2][1][2][1] = lambda / 0.2e1;
	Curvature_Higgs_L4[2][2][0][0] = lambda / 0.2e1;
	Curvature_Higgs_L4[2][2][1][1] = lambda / 0.2e1;
	Curvature_Higgs_L4[2][2][2][2] = 0.3e1 / 0.2e1 * lambda;
	Curvature_Higgs_L4[2][2][3][3] = lambda / 0.2e1;
	Curvature_Higgs_L4[2][3][2][3] = lambda / 0.2e1;
	Curvature_Higgs_L4[2][3][3][2] = lambda / 0.2e1;
	Curvature_Higgs_L4[3][0][0][3] = lambda / 0.2e1;
	Curvature_Higgs_L4[3][0][3][0] = lambda / 0.2e1;
	Curvature_Higgs_L4[3][1][1][3] = lambda / 0.2e1;
	Curvature_Higgs_L4[3][1][3][1] = lambda / 0.2e1;
	Curvature_Higgs_L4[3][2][2][3] = lambda / 0.2e1;
	Curvature_Higgs_L4[3][2][3][2] = lambda / 0.2e1;
	Curvature_Higgs_L4[3][3][0][0] = lambda / 0.2e1;
	Curvature_Higgs_L4[3][3][1][1] = lambda / 0.2e1;
	Curvature_Higgs_L4[3][3][2][2] = lambda / 0.2e1;
	Curvature_Higgs_L4[3][3][3][3] = 0.3e1 / 0.2e1 * lambda;


	Curvature_Gauge_G2H2[0][0][0][0] = C_g * C_g / 0.2e1;
	Curvature_Gauge_G2H2[0][0][1][1] = C_g * C_g / 0.2e1;
	Curvature_Gauge_G2H2[0][0][2][2] = C_g * C_g / 0.2e1;
	Curvature_Gauge_G2H2[0][0][3][3] = C_g * C_g / 0.2e1;
	Curvature_Gauge_G2H2[0][3][0][3] = C_gs * C_g / 0.2e1;
	Curvature_Gauge_G2H2[0][3][1][2] = C_gs * C_g / 0.2e1;
	Curvature_Gauge_G2H2[0][3][2][1] = C_gs * C_g / 0.2e1;
	Curvature_Gauge_G2H2[0][3][3][0] = C_gs * C_g / 0.2e1;
	Curvature_Gauge_G2H2[1][1][0][0] = C_g * C_g / 0.2e1;
	Curvature_Gauge_G2H2[1][1][1][1] = C_g * C_g / 0.2e1;
	Curvature_Gauge_G2H2[1][1][2][2] = C_g * C_g / 0.2e1;
	Curvature_Gauge_G2H2[1][1][3][3] = C_g * C_g / 0.2e1;
	Curvature_Gauge_G2H2[1][3][0][2] = C_gs * C_g / 0.2e1;
	Curvature_Gauge_G2H2[1][3][1][3] = -C_gs * C_g / 0.2e1;
	Curvature_Gauge_G2H2[1][3][2][0] = C_gs * C_g / 0.2e1;
	Curvature_Gauge_G2H2[1][3][3][1] = -C_gs * C_g / 0.2e1;
	Curvature_Gauge_G2H2[2][2][0][0] = C_g * C_g / 0.2e1;
	Curvature_Gauge_G2H2[2][2][1][1] = C_g * C_g / 0.2e1;
	Curvature_Gauge_G2H2[2][2][2][2] = C_g * C_g / 0.2e1;
	Curvature_Gauge_G2H2[2][2][3][3] = C_g * C_g / 0.2e1;
	Curvature_Gauge_G2H2[2][3][0][0] = C_gs * C_g / 0.2e1;
	Curvature_Gauge_G2H2[2][3][1][1] = C_gs * C_g / 0.2e1;
	Curvature_Gauge_G2H2[2][3][2][2] = -C_gs * C_g / 0.2e1;
	Curvature_Gauge_G2H2[2][3][3][3] = -C_gs * C_g / 0.2e1;
	Curvature_Gauge_G2H2[3][0][0][3] = C_gs * C_g / 0.2e1;
	Curvature_Gauge_G2H2[3][0][1][2] = C_gs * C_g / 0.2e1;
	Curvature_Gauge_G2H2[3][0][2][1] = C_gs * C_g / 0.2e1;
	Curvature_Gauge_G2H2[3][0][3][0] = C_gs * C_g / 0.2e1;
	Curvature_Gauge_G2H2[3][1][0][2] = C_gs * C_g / 0.2e1;
	Curvature_Gauge_G2H2[3][1][1][3] = -C_gs * C_g / 0.2e1;
	Curvature_Gauge_G2H2[3][1][2][0] = C_gs * C_g / 0.2e1;
	Curvature_Gauge_G2H2[3][1][3][1] = -C_gs * C_g / 0.2e1;
	Curvature_Gauge_G2H2[3][2][0][0] = C_gs * C_g / 0.2e1;
	Curvature_Gauge_G2H2[3][2][1][1] = C_gs * C_g / 0.2e1;
	Curvature_Gauge_G2H2[3][2][2][2] = -C_gs * C_g / 0.2e1;
	Curvature_Gauge_G2H2[3][2][3][3] = -C_gs * C_g / 0.2e1;
	Curvature_Gauge_G2H2[3][3][0][0] = C_gs * C_gs / 0.2e1;
	Curvature_Gauge_G2H2[3][3][1][1] = C_gs * C_gs / 0.2e1;
	Curvature_Gauge_G2H2[3][3][2][2] = C_gs * C_gs / 0.2e1;
	Curvature_Gauge_G2H2[3][3][3][3] = C_gs * C_gs / 0.2e1;

	double v = C_vev0;
	std::complex<double> V11,V12,V13,V21,V22,V23,V31,V32,V33;
	V11 = C_Vud;
	V12 = C_Vus;
	V13 = C_Vub;
	V21 = C_Vcd;
	V22 = C_Vcs;
	V23 = C_Vcb;
	V31 = C_Vtd;
	V32 = C_Vts;
	V33 = C_Vtb;

	Curvature_Lepton_F2H1[0][1][2] = II / v * C_MassElectron;
	Curvature_Lepton_F2H1[0][1][3] = 0.1e1 / v * C_MassElectron;
	Curvature_Lepton_F2H1[1][0][2] = II / v * C_MassElectron;
	Curvature_Lepton_F2H1[1][0][3] = 0.1e1 / v * C_MassElectron;
	Curvature_Lepton_F2H1[1][6][0] = 0.1e1 / v * C_MassElectron;
	Curvature_Lepton_F2H1[1][6][1] = II / v * C_MassElectron;
	Curvature_Lepton_F2H1[2][3][2] = II / v * C_MassMu;
	Curvature_Lepton_F2H1[2][3][3] = 0.1e1 / v * C_MassMu;
	Curvature_Lepton_F2H1[3][2][2] = II / v * C_MassMu;
	Curvature_Lepton_F2H1[3][2][3] = 0.1e1 / v * C_MassMu;
	Curvature_Lepton_F2H1[3][7][0] = 0.1e1 / v * C_MassMu;
	Curvature_Lepton_F2H1[3][7][1] = II / v * C_MassMu;
	Curvature_Lepton_F2H1[4][5][2] = II / v * C_MassTau;
	Curvature_Lepton_F2H1[4][5][3] = 0.1e1 / v * C_MassTau;
	Curvature_Lepton_F2H1[5][4][2] = II / v * C_MassTau;
	Curvature_Lepton_F2H1[5][4][3] = 0.1e1 / v * C_MassTau;
	Curvature_Lepton_F2H1[5][8][0] = 0.1e1 / v * C_MassTau;
	Curvature_Lepton_F2H1[5][8][1] = II / v * C_MassTau;
	Curvature_Lepton_F2H1[6][1][0] = 0.1e1 / v * C_MassElectron;
	Curvature_Lepton_F2H1[6][1][1] = II / v * C_MassElectron;
	Curvature_Lepton_F2H1[7][3][0] = 0.1e1 / v * C_MassMu;
	Curvature_Lepton_F2H1[7][3][1] = II / v * C_MassMu;
	Curvature_Lepton_F2H1[8][5][0] = 0.1e1 / v * C_MassTau;
	Curvature_Lepton_F2H1[8][5][1] = II / v * C_MassTau;

	Curvature_Quark_F2H1[0][6][2] = -II / v * C_MassUp;
	Curvature_Quark_F2H1[0][6][3] = 0.1e1 / v * C_MassUp;
	Curvature_Quark_F2H1[0][9][0] = -0.1e1 / v * C_MassUp * conj(V11);
	Curvature_Quark_F2H1[0][9][1] = II / v * C_MassUp * conj(V11);
	Curvature_Quark_F2H1[0][10][0] = -0.1e1 / v * C_MassUp * conj(V12);
	Curvature_Quark_F2H1[0][10][1] = II / v * C_MassUp * conj(V12);
	Curvature_Quark_F2H1[0][11][0] = -0.1e1 / v * C_MassUp * conj(V13);
	Curvature_Quark_F2H1[0][11][1] = II / v * C_MassUp * conj(V13);
	Curvature_Quark_F2H1[1][7][2] = -II / v * C_MassCharm;
	Curvature_Quark_F2H1[1][7][3] = 0.1e1 / v * C_MassCharm;
	Curvature_Quark_F2H1[1][9][0] = -0.1e1 / v * C_MassCharm * conj(V21);
	Curvature_Quark_F2H1[1][9][1] = II / v * C_MassCharm * conj(V21);
	Curvature_Quark_F2H1[1][10][0] = -0.1e1 / v * C_MassCharm * conj(V22);
	Curvature_Quark_F2H1[1][10][1] = II / v * C_MassCharm * conj(V22);
	Curvature_Quark_F2H1[1][11][0] = -0.1e1 / v * C_MassCharm * conj(V23);
	Curvature_Quark_F2H1[1][11][1] = II / v * C_MassCharm * conj(V23);
	Curvature_Quark_F2H1[2][8][2] = -II / v * C_MassTop;
	Curvature_Quark_F2H1[2][8][3] = 0.1e1 / v * C_MassTop;
	Curvature_Quark_F2H1[2][9][0] = -0.1e1 / v * C_MassTop * conj(V31);
	Curvature_Quark_F2H1[2][9][1] = II / v * C_MassTop * conj(V31);
	Curvature_Quark_F2H1[2][10][0] = -0.1e1 / v * C_MassTop * conj(V32);
	Curvature_Quark_F2H1[2][10][1] = II / v * C_MassTop * conj(V32);
	Curvature_Quark_F2H1[2][11][0] = -0.1e1 / v * C_MassTop * conj(V33);
	Curvature_Quark_F2H1[2][11][1] = II / v * C_MassTop * conj(V33);
	Curvature_Quark_F2H1[3][6][0] = 0.1e1 / v * C_MassDown * V11;
	Curvature_Quark_F2H1[3][6][1] = II / v * C_MassDown * V11;
	Curvature_Quark_F2H1[3][7][0] = V21 / v * C_MassDown;
	Curvature_Quark_F2H1[3][7][1] = II * V21 / v * C_MassDown;
	Curvature_Quark_F2H1[3][8][0] = 0.1e1 / v * C_MassDown * V31;
	Curvature_Quark_F2H1[3][8][1] = II / v * C_MassDown * V31;
	Curvature_Quark_F2H1[3][9][2] = II / v * C_MassDown;
	Curvature_Quark_F2H1[3][9][3] = 0.1e1 / v * C_MassDown;
	Curvature_Quark_F2H1[4][6][0] = 0.1e1 / v * C_MassStrange * V12;
	Curvature_Quark_F2H1[4][6][1] = II * V12 / v * C_MassStrange;
	Curvature_Quark_F2H1[4][7][0] = V22 / v * C_MassStrange;
	Curvature_Quark_F2H1[4][7][1] = II * V22 / v * C_MassStrange;
	Curvature_Quark_F2H1[4][8][0] = 0.1e1 / v * C_MassStrange * V32;
	Curvature_Quark_F2H1[4][8][1] = II * V32 / v * C_MassStrange;
	Curvature_Quark_F2H1[4][10][2] = II / v * C_MassStrange;
	Curvature_Quark_F2H1[4][10][3] = 0.1e1 / v * C_MassStrange;
	Curvature_Quark_F2H1[5][6][0] = V13 / v * C_MassBottom;
	Curvature_Quark_F2H1[5][6][1] = II / v * C_MassBottom * V13;
	Curvature_Quark_F2H1[5][7][0] = V23 / v * C_MassBottom;
	Curvature_Quark_F2H1[5][7][1] = II / v * C_MassBottom * V23;
	Curvature_Quark_F2H1[5][8][0] = V33 / v * C_MassBottom;
	Curvature_Quark_F2H1[5][8][1] = II / v * C_MassBottom * V33;
	Curvature_Quark_F2H1[5][11][2] = II / v * C_MassBottom;
	Curvature_Quark_F2H1[5][11][3] = 0.1e1 / v * C_MassBottom;
	Curvature_Quark_F2H1[6][0][2] = -II / v * C_MassUp;
	Curvature_Quark_F2H1[6][0][3] = 0.1e1 / v * C_MassUp;
	Curvature_Quark_F2H1[6][3][0] = 0.1e1 / v * C_MassDown * V11;
	Curvature_Quark_F2H1[6][3][1] = II / v * C_MassDown * V11;
	Curvature_Quark_F2H1[6][4][0] = 0.1e1 / v * C_MassStrange * V12;
	Curvature_Quark_F2H1[6][4][1] = II * V12 / v * C_MassStrange;
	Curvature_Quark_F2H1[6][5][0] = V13 / v * C_MassBottom;
	Curvature_Quark_F2H1[6][5][1] = II / v * C_MassBottom * V13;
	Curvature_Quark_F2H1[7][1][2] = -II / v * C_MassCharm;
	Curvature_Quark_F2H1[7][1][3] = 0.1e1 / v * C_MassCharm;
	Curvature_Quark_F2H1[7][3][0] = V21 / v * C_MassDown;
	Curvature_Quark_F2H1[7][3][1] = II * V21 / v * C_MassDown;
	Curvature_Quark_F2H1[7][4][0] = V22 / v * C_MassStrange;
	Curvature_Quark_F2H1[7][4][1] = II * V22 / v * C_MassStrange;
	Curvature_Quark_F2H1[7][5][0] = V23 / v * C_MassBottom;
	Curvature_Quark_F2H1[7][5][1] = II / v * C_MassBottom * V23;
	Curvature_Quark_F2H1[8][2][2] = -II / v * C_MassTop;
	Curvature_Quark_F2H1[8][2][3] = 0.1e1 / v * C_MassTop;
	Curvature_Quark_F2H1[8][3][0] = 0.1e1 / v * C_MassDown * V31;
	Curvature_Quark_F2H1[8][3][1] = II / v * C_MassDown * V31;
	Curvature_Quark_F2H1[8][4][0] = 0.1e1 / v * C_MassStrange * V32;
	Curvature_Quark_F2H1[8][4][1] = II * V32 / v * C_MassStrange;
	Curvature_Quark_F2H1[8][5][0] = V33 / v * C_MassBottom;
	Curvature_Quark_F2H1[8][5][1] = II / v * C_MassBottom * V33;
	Curvature_Quark_F2H1[9][0][0] = -0.1e1 / v * C_MassUp * conj(V11);
	Curvature_Quark_F2H1[9][0][1] = II / v * C_MassUp * conj(V11);
	Curvature_Quark_F2H1[9][1][0] = -0.1e1 / v * C_MassCharm * conj(V21);
	Curvature_Quark_F2H1[9][1][1] = II / v * C_MassCharm * conj(V21);
	Curvature_Quark_F2H1[9][2][0] = -0.1e1 / v * C_MassTop * conj(V31);
	Curvature_Quark_F2H1[9][2][1] = II / v * C_MassTop * conj(V31);
	Curvature_Quark_F2H1[9][3][2] = II / v * C_MassDown;
	Curvature_Quark_F2H1[9][3][3] = 0.1e1 / v * C_MassDown;
	Curvature_Quark_F2H1[10][0][0] = -0.1e1 / v * C_MassUp * conj(V12);
	Curvature_Quark_F2H1[10][0][1] = II / v * C_MassUp * conj(V12);
	Curvature_Quark_F2H1[10][1][0] = -0.1e1 / v * C_MassCharm * conj(V22);
	Curvature_Quark_F2H1[10][1][1] = II / v * C_MassCharm * conj(V22);
	Curvature_Quark_F2H1[10][2][0] = -0.1e1 / v * C_MassTop * conj(V32);
	Curvature_Quark_F2H1[10][2][1] = II / v * C_MassTop * conj(V32);
	Curvature_Quark_F2H1[10][4][2] = II / v * C_MassStrange;
	Curvature_Quark_F2H1[10][4][3] = 0.1e1 / v * C_MassStrange;
	Curvature_Quark_F2H1[11][0][0] = -0.1e1 / v * C_MassUp * conj(V13);
	Curvature_Quark_F2H1[11][0][1] = II / v * C_MassUp * conj(V13);
	Curvature_Quark_F2H1[11][1][0] = -0.1e1 / v * C_MassCharm * conj(V23);
	Curvature_Quark_F2H1[11][1][1] = II / v * C_MassCharm * conj(V23);
	Curvature_Quark_F2H1[11][2][0] = -0.1e1 / v * C_MassTop * conj(V33);
	Curvature_Quark_F2H1[11][2][1] = II / v * C_MassTop * conj(V33);
	Curvature_Quark_F2H1[11][5][2] = II / v * C_MassBottom;
	Curvature_Quark_F2H1[11][5][3] = 0.1e1 / v * C_MassBottom;




}


bool Class_SM::CalculateDebyeSimplified(){
  return false;
  /*
   * Use this function if you calculated the Debye corrections to the Higgs mass matrix and implement
   * your formula here and return true. The vector is given by DebyeHiggs[NHiggs][NHiggs]
   */
}

bool Class_SM::CalculateDebyeGaugeSimplified()
{
	bool Debug = false;
  if(Debug) std::cout << "Debug turned on in Class_SM :: " << __func__ << std::endl;

  /*
     * Use this function if you calculated the Debye corrections to the gauge mass matrix and implement
     * your formula here and return true. The vector is given by DebyeGauge[NGauge][NGauge]
     */


  return false;
}
double Class_SM::VTreeSimplified(const std::vector<double>& v) const
{
	(void) v;
	if(not UseVTreeSimplified) return 0;
	double res = 0;

	return res;
}

double Class_SM::VCounterSimplified(const std::vector<double>& v) const
{
	(void) v;
	if(not UseVCounterSimplified) return 0;
	double res = 0;
	return res;
}

void Class_SM::Debugging(const std::vector<double>& input, std::vector<double>& output) const
{
	(void) input;
	(void) output;

}
}
}
