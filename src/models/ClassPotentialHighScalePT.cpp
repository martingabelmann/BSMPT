/*
 * Class_HighScalePT.cpp
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

#include <ext/alloc_traits.h>               // for __alloc_traits<>::value_type
#include <stddef.h>                         // for std::size_t
#include <algorithm>                        // for max, copy
#include <iostream>                         // for operator<<, endl, basic_o...
#include <memory>                           // for allocator_traits<>::value...
#include <BSMPT/models/SMparam.h>           // for C_vev0, C_MassTop, C_g
#include "Eigen/Dense"
#include "Eigen/Eigenvalues"
#include "Eigen/IterativeLinearSolvers"

#include <BSMPT/models/ClassPotentialHighScalePT.h>

#include <BSMPT/models/IncludeAllModels.h>

#include <BSMPT/utility.h>
using namespace Eigen;

/**
 * SM + EWinos
 *
 */

namespace BSMPT{
namespace Models{
/**
 * Here you have to adjust NNeutralHiggs, NChargedHiggs, nPar (number of Lagrangian parameters AFTER
 *  using the tadpole conditions),
 * nParCT (number of counterterms) as well as nVEV (number of VEVs for minimization)
 */
Class_HighScalePT::Class_HighScalePT ()
{
  Model = ModelID::ModelIDs::HSPT; // global int constant which will be used to tell the program which model is called
  NNeutralHiggs = 1; // number of neutral Higgs bosons at T = 0
  NChargedHiggs=0; // number of charged Higgs bosons  at T = 0 (all d.o.f.)

  nPar = 9; // number of parameters in the tree-Level Lagrangian
  nParCT = 3; // number of parameters in the counterterm potential

  nVEV=1; // number of VEVs to minimize the potential

  NHiggs = NNeutralHiggs+NChargedHiggs;
  NNonSMFermion=6; //Number of Non-SM like fermions

  VevOrder.resize(nVEV);
  // Here you have to tell which scalar field gets which VEV.
  VevOrder[0] = 0;

  // Set UseVTreeSimplified to use the tree-level potential defined in VTreeSimplified
  UseVTreeSimplified = false;

  // Set UseVCounterSimplified to use the counterterm potential defined in VCounterSimplified
  UseVCounterSimplified = false;

}

Class_HighScalePT::~Class_HighScalePT ()
{
}

/**
 * returns a string which tells the user the chronological order of the counterterms. Use this to
 * complement the legend of the given input file
 */
std::vector<std::string> Class_HighScalePT::addLegendCT() const{
    std::vector<std::string> labels;
    labels.push_back("dT");
    labels.push_back("dlambda");
    labels.push_back("dmsquared");
    return labels;
}

/**
 * returns a string which tells the user the chronological order of the VEVs and the critical temperature. Use this to
 * complement the legend of the given input file
 */
std::vector<std::string> Class_HighScalePT::addLegendTemp() const{
    std::vector<std::string> labels;
    labels.push_back("T_c"); // Label for the critical temperature
    labels.push_back("v_c"); // Label for the critical vev
    labels.push_back("v_c/T_c"); // Label for v_c/T_c, you could use xi_c also for example
    //out += "Your VEV order"; // Now you have to put the label for your vevs
    labels.push_back("omega");
    return labels;
}

/**
 * returns a string which tells the user the chronological order of the Triple Higgs couplings. Use this to
 * complement the legend of the given input file
 *
 */
std::vector<std::string> Class_HighScalePT::addLegendTripleCouplings() const{
    std::vector<std::string> labels;
	std::vector<std::string> particles;
	particles.resize(NHiggs);
	//here you have to define the particle names in the vector particles

	particles[0]="H";

    std::string out = "Tree_";
    for(std::size_t i=0;i<NHiggs;i++)
    {
        for(std::size_t j=i;j<NHiggs;j++)
        {
            for(std::size_t k=j;k<NHiggs;k++)
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
std::vector<std::string> Class_HighScalePT::addLegendVEV() const{
    std::vector<std::string> labels;
    //out = "Your VEV order";
    labels.push_back("omega");
    return labels;
}

/**
 * Reads the string linestr and sets the parameter point
 */
void Class_HighScalePT::ReadAndSet(const std::string& linestr, std::vector<double>& par )
{
	std::stringstream ss(linestr);
	double tmp;

	double lms,llambda;
	double lmuchi,lmb,lmw,lg1u,lg2u,lg1d,lg2d;

    if (UseIndexCol){
        ss >> tmp;
    }


	for(std::size_t k=1;k<=nPar;k++)
	{
	      ss>>tmp;
	      if(k==1) lms = tmp;
	      else if(k==2) llambda = tmp;
		  else if(k==3) lmuchi =tmp;
		  else if(k==4) lmw = tmp;
		  else if(k==5) lmb = tmp;
		  else if(k==6) lg1u = tmp;
		  else if(k==7) lg2u = tmp;
		  else if(k==8) lg1d = tmp;
		  else if(k==9) lg2d = tmp;
	}
	par[0] = lms;
	par[1] = llambda;
	par[2] = lmuchi;
	par[3] = lmw;
	par[4] = lmb;
	par[5] = lg1u;
	par[6] = lg2u;
	par[7] = lg1d;
	par[8] = lg2d;


	set_gen(par); // This you have to call so that everything will be set
	return ;
}


/**
 * Set Class Object as well as the VEV configuration
 */
void Class_HighScalePT::set_gen(const std::vector<double>& par) {
    ms = par[0]; // Class member is set accordingly to the input parameters
    lambda = par[1]; // Class member is set accordingly to the input parameters

	//nonSM Fermions
	muchi = par[2];
	mw    = par[3];
	mb    = par[4];
	g1u   = par[5];
	g2u   = par[6];
	g1d   = par[7];
	g2d   = par[8];
    
	g=C_g; // SM SU (2) gauge coupling --> SMparam .h
    yt = std::sqrt(2)/C_vev0 * C_MassTop; // Top Yukawa coupling --> SMparam .h
    scale = C_vev0; // Renormalisation scale is set to the SM VEV
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
void Class_HighScalePT::set_CT_Pot_Par(const std::vector<double>& par){

	dT = par[0];
	dlambda = par[1];
	dms = par[2];

	Curvature_Higgs_CT_L1[0] = dT;
	Curvature_Higgs_CT_L2[0][0] = dms;
	Curvature_Higgs_CT_L4[0][0][0][0] = dlambda;
}


/**
 * console output of all Parameters
 */
void Class_HighScalePT::write() const {

    std::cout << "Model = " << Model << std::endl;

	std::cout << "The parameters are : " << std::endl;
	std::cout << "lambda = " << lambda << std::endl
            << "m^2 = " << ms << std::endl
			<< "mu_chi = "<< muchi << std::endl
			<< "M_W = "<< muchi << std::endl
			<< "M_B = "<< muchi << std::endl
			<< "g_1u =" << g1u << std::endl
			<< "g_2u =" << g2u << std::endl
			<< "g_1d =" << g1d << std::endl
			<< "g_2d =" << g2d << std::endl;

	std::cout << "The counterterm parameters are : " << std::endl;
	std::cout << "dT = "<< dT << std::endl
			<< "dlambda = " << dlambda << std::endl
			<< "dm^2 = "<< dms << std::endl;

	std::cout << "The scale is given by mu = " << scale << " GeV " << std::endl;

}


/**
 * Calculates the counterterms. Here you need to work out the scheme and implement the formulas.
 */
std::vector<double> Class_HighScalePT::calc_CT() const {

    std::vector<double> parCT;

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

	std::vector<double> WeinbergNabla,WeinbergHesse;
    WeinbergNabla = WeinbergFirstDerivative();
    WeinbergHesse = WeinbergSecondDerivative();

	VectorXd NablaWeinberg(NHiggs);
	MatrixXd HesseWeinberg(NHiggs,NHiggs),HiggsRot(NHiggs,NHiggs);
	for(std::size_t i=0;i<NHiggs;i++)
	{
		NablaWeinberg[i] = WeinbergNabla[i];
		for(std::size_t j=0;j<NHiggs;j++) HesseWeinberg(i,j) = WeinbergHesse.at(j*NHiggs+i);
	}

    // Here you have to use your formulae for the counterterm scheme
	double t = 0;
    parCT.push_back(t); // dT
    parCT.push_back(3.0*t/std::pow(C_vev0,3) + 3.0/std::pow(C_vev0,3) * NablaWeinberg(0) -3.0/std::pow(C_vev0,2) *HesseWeinberg(0,0)); // dlambda
    parCT.push_back(-3.0/(2*std::pow(C_vev0,2)) *NablaWeinberg(0) + 1.0/2.0 *HesseWeinberg(0,0) -3.0*t/(2*C_vev0)); // dms


    return parCT;

}




void Class_HighScalePT::TripleHiggsCouplings()
{
	if(!SetCurvatureDone)SetCurvatureArrays();
	if(!CalcCouplingsdone)CalculatePhysicalCouplings();


	std::vector<double> HiggsOrder(NHiggs);
	// Here you have to set the vector HiggsOrder. By telling e.g. HiggsOrder[0] = 5 you always want your 6th lightest
	// particle to be the first particle in the vector (which has the index 5 because they are sorted by mass)

	// example for keeping the mass order
	for(std::size_t i=0;i<NHiggs;i++) {
		HiggsOrder[i]=i;
	}

	std::vector<double> TripleDeriv;
    TripleDeriv=WeinbergThirdDerivative();
	std::vector<std::vector<std::vector<double>>> GaugeBasis(NHiggs, std::vector<std::vector<double>>(NHiggs,
				std::vector<double>(NHiggs)));
	for(std::size_t i=0;i<NHiggs;i++)
	  {
		for(std::size_t j=0;j<NHiggs;j++)
		{
		  for(std::size_t k=0;k<NHiggs;k++)
			{
			  GaugeBasis[i][j][k] = TripleDeriv.at(i+j*NHiggs+k*NHiggs*NHiggs);
			}
		}
	  }

	MatrixXd HiggsRot(NHiggs,NHiggs);
	for(std::size_t i=0;i<NHiggs;i++)
	{
		for(std::size_t j=0;j<NHiggs;j++)
		{
			HiggsRot(i,j) = HiggsRotationMatrix[i][j];
		}
	}


	MatrixXd HiggsRotSort(NHiggs,NHiggs);




	for(std::size_t i=0;i<NHiggs;i++)
	{
		HiggsRotSort.row(i) = HiggsRot.row(HiggsOrder[i]);
	}

	TripleHiggsCorrectionsCWPhysical.resize(NHiggs);
	TripleHiggsCorrectionsTreePhysical.resize(NHiggs);
	TripleHiggsCorrectionsCTPhysical.resize(NHiggs);
	for(std::size_t i=0;i<NHiggs;i++) {
		TripleHiggsCorrectionsTreePhysical[i].resize(NHiggs);
		TripleHiggsCorrectionsCWPhysical[i].resize(NHiggs);
		TripleHiggsCorrectionsCTPhysical[i].resize(NHiggs);
		for(std::size_t j=0;j<NHiggs;j++) {
			TripleHiggsCorrectionsCWPhysical[i][j].resize(NHiggs);
			TripleHiggsCorrectionsTreePhysical[i][j].resize(NHiggs);
			TripleHiggsCorrectionsCTPhysical[i][j].resize(NHiggs);
		}
	}


	for(std::size_t i=0;i<NHiggs;i++)
	  {
		for(std::size_t j=0;j<NHiggs;j++)
		{
			for(std::size_t k=0;k<NHiggs;k++)
			{
			  TripleHiggsCorrectionsCWPhysical[i][j][k] = 0;
			  TripleHiggsCorrectionsTreePhysical[i][j][k] = 0;
			  TripleHiggsCorrectionsCTPhysical[i][j][k] = 0;
			  for(std::size_t l=0;l<NHiggs;l++)
			  {
				  for(std::size_t m=0;m<NHiggs;m++)
				  {
					  for(std::size_t n=0;n<NHiggs;n++)
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

void Class_HighScalePT::SetCurvatureArrays(){
  /*
   *  Here you have to set the vectors
   *  Curvature_Higgs_L1,Curvature_Higgs_L2,Curvature_Higgs_L3,Curvature_Higgs_L4
   *  Curvature_Gauge_G2H2
   *  Curvature_Quark_F2H1, Curvature_Lepton_F2H1
   *  as described in the potential in the paper.
   */

	initVectors();
	SetCurvatureDone=true;
	for(std::size_t i=0;i<NHiggs;i++) HiggsVev[i] = vevTree[i];


	Curvature_Higgs_L2[0][0] = ms;
	Curvature_Higgs_L4[0][0][0][0] = lambda;

	Curvature_Gauge_G2H2[0][0][0][0] = 4*std::pow(g,2);

	Curvature_Quark_F2H1[1][0][0] = yt;
	Curvature_Quark_F2H1[0][1][0] = yt;

    /* neutralino masses */
	Curvature_NonSMFermion_F2[0][0] = mb;
	Curvature_NonSMFermion_F2[0][1] = 0;
	Curvature_NonSMFermion_F2[0][2] = -g1d*C_vev0/2;
	Curvature_NonSMFermion_F2[0][3] = -g1u*C_vev0/2;
	Curvature_NonSMFermion_F2[1][1] = mw;
	Curvature_NonSMFermion_F2[1][2] = g2d*C_vev0/2;
	Curvature_NonSMFermion_F2[1][3] = -g2u*C_vev0/2;
	Curvature_NonSMFermion_F2[2][2] = 0;
	Curvature_NonSMFermion_F2[2][3] = -muchi;
	Curvature_NonSMFermion_F2[3][3] = 0;
	Curvature_NonSMFermion_F2[3][2] = Curvature_NonSMFermion_F2[2][3];
	Curvature_NonSMFermion_F2[3][1] = Curvature_NonSMFermion_F2[1][3];
	Curvature_NonSMFermion_F2[2][1] = Curvature_NonSMFermion_F2[1][2];
	Curvature_NonSMFermion_F2[3][0] = Curvature_NonSMFermion_F2[0][3];
	Curvature_NonSMFermion_F2[1][0] = Curvature_NonSMFermion_F2[0][1];

    /* chargino masses */
	Curvature_NonSMFermion_F2[4][4] = mw;
	Curvature_NonSMFermion_F2[4][5] = g2u*C_vev0/2;
	Curvature_NonSMFermion_F2[5][4] = g2d*C_vev0/2;
	Curvature_NonSMFermion_F2[5][5] = muchi;

    /* neutralino yukawas */
	Curvature_NonSMFermion_F2H1[0][2][0] = g1d/2;
	Curvature_NonSMFermion_F2H1[0][3][0] = -g1u/2;
	Curvature_NonSMFermion_F2H1[1][2][0] = -g2d/2;
	Curvature_NonSMFermion_F2H1[1][3][0] = g2u/2;
	
    Curvature_NonSMFermion_F2H1[2][0][0] = Curvature_NonSMFermion_F2H1[0][2][0];
    Curvature_NonSMFermion_F2H1[3][0][0] = Curvature_NonSMFermion_F2H1[0][3][0];
    Curvature_NonSMFermion_F2H1[2][1][0] = Curvature_NonSMFermion_F2H1[1][2][0];
    Curvature_NonSMFermion_F2H1[3][1][0] = Curvature_NonSMFermion_F2H1[1][3][0];

	Curvature_NonSMFermion_F2H1[0][0][0]=0;
	Curvature_NonSMFermion_F2H1[1][0][0]=0;
	Curvature_NonSMFermion_F2H1[0][1][0]=0;
	Curvature_NonSMFermion_F2H1[1][1][0]=0;
	Curvature_NonSMFermion_F2H1[2][2][0]=0;
	Curvature_NonSMFermion_F2H1[2][3][0]=0;
	Curvature_NonSMFermion_F2H1[3][2][0]=0;
	Curvature_NonSMFermion_F2H1[3][3][0]=0;
	
    /* chargino yukawas */
	Curvature_NonSMFermion_F2H1[4][4][0] = 0;
	Curvature_NonSMFermion_F2H1[4][5][0] = -g2u/std::sqrt(2);
	Curvature_NonSMFermion_F2H1[5][4][0] = -g2d/std::sqrt(2);
	Curvature_NonSMFermion_F2H1[5][5][0] = 0;

}


bool Class_HighScalePT::CalculateDebyeSimplified(){
  return false;
  /*
   * Use this function if you calculated the Debye corrections to the Higgs mass matrix and implement
   * your formula here and return true. The vector is given by DebyeHiggs[NHiggs][NHiggs]
   */
}

bool Class_HighScalePT::CalculateDebyeGaugeSimplified()
{
  /*
     * Use this function if you calculated the Debye corrections to the gauge mass matrix and implement
     * your formula here and return true. The vector is given by DebyeGauge[NGauge][NGauge]
     */


  return false;
}
double Class_HighScalePT::VTreeSimplified(const std::vector<double>& v) const {
    if(not UseVTreeSimplified) return  0;
	double res = 0;

	double vIn = v[0];
	res = 0.5*ms*std::pow(vIn,2) + 1.0/24.0 * lambda * std::pow(vIn,4);

	return res;
}

double Class_HighScalePT::VCounterSimplified(const std::vector<double>& v) const
{
    if(not UseVCounterSimplified) return 0;
	double res = 0;
	double vIn = v[0];
		res = 0.5*dms*std::pow(vIn,2) + 1.0/24.0 * dlambda * std::pow(vIn,4) + dT * vIn;
	return res;
}

void Class_HighScalePT::Debugging(const std::vector<double>& input, std::vector<double>& output) const
{
	(void) input;
	(void) output;

}

}
}
