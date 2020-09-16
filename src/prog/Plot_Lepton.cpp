/*
 * PlotEWBG_vw.cpp
 *
 *
 *      Copyright (C) 2020  Philipp Basler, Margarete Mühlleitner and Jonas Müller

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
 * This program calculates the EWBG eta as a function of vw and varies vw over a given array.
 */

#include <bits/exception.h>					   // for exception
#include <stdlib.h>							   // for atof, EXI...
#include <algorithm>						   // for copy, max
#include <memory>							   // for shared_ptr
#include <string>							   // for operator<<
#include <vector>							   // for vector
#include <BSMPT/models/ClassPotentialOrigin.h> // for Class_Pot...
#include <BSMPT/models/IncludeAllModels.h>
#include <BSMPT/minimizer/Minimizer.h>
#include <BSMPT/baryo_calculation/CalculateEtaInterface.h>
#include <BSMPT/utility.h>
#include <iostream>
#include <fstream>
using namespace std;
using namespace BSMPT;

int main(int argc, char *argv[])
try
{
	if ((argc != 11) and (argc != 12))
	{
		std::cerr << "./PlotEWBG_vw Model Inputfile Outputfile Line mchi_start mchi_stepsize mchi_end LambdaQL source_flag EWBGConfigFile TerminalOutput(y/n)\n";
		ShowInputError();
		return EXIT_FAILURE;
	}
	char *in_file;
	char *out_file;
	in_file = argv[2];
	out_file = argv[3];
	double LineNumb;

	auto Model = ModelID::getModel(argv[1]);
	if (Model == ModelID::ModelIDs::NotSet)
	{
		std::cerr << "Your Model parameter does not match with the implemented Models." << std::endl;
		ShowInputError();
		return EXIT_FAILURE;
	}

	LineNumb = atoi(argv[4]);
	double mchi_Start, mchi_End, mchi_Stepsize, LambdaQL;
	double source_flag;
	mchi_Start = atof(argv[5]);
	mchi_Stepsize = atof(argv[6]);
	mchi_End = atof(argv[7]);
	LambdaQL = atof(argv[8]);
	source_flag = atof(argv[9]);
	if (LineNumb < 1)
	{
		std::cerr << "Start line counting with 1" << std::endl;
		return EXIT_FAILURE;
	}
	bool TerminalOutput = false;
	if (argc == 12)
	{
		std::string s7 = argv[11];
		std::cout << "Terminal Output:" << s7 << std::endl;
		TerminalOutput = ("y" == s7);
	}
	//Set up of BSMPT/Baryo Classes
	Baryo::CalculateEtaInterface EtaInterface(argv[10] /* = Config File */);
	std::shared_ptr<Class_Potential_Origin> modelPointer = ModelID::FChoose(Model);

	std::vector<double> start, solPot;

	std::ifstream infile(in_file);
	if (!infile.good())
	{
		std::cout << "Input file not found " << std::endl;
		return EXIT_FAILURE;
	}

	std::ofstream outfile(out_file);
	if (!outfile.good())
	{
		std::cout << "Can not create file " << out_file << std::endl;
		return EXIT_FAILURE;
	}

	std::string linestr;
	int linecounter = 1;
	bool found = false;

	while (true)
	{
		if (infile.eof())
			break;
		std::getline(infile, linestr);
		if (linecounter == 1)
		{
			modelPointer->setUseIndexCol(linestr);
			//outfile: LEGEND
			outfile << linestr << sep;
			outfile << "T_c_var" << sep << "omega_c_var" << sep << "vw_var" << sep << "LW_var" << sep << "mChi" << sep << "LambdaQL" << sep << "source_flag";
			auto legend = EtaInterface.legend();
			for (const auto &x : legend)
				outfile << sep << x + "_var";
			outfile << std::endl;
		}
		else if (linecounter == LineNumb)
		{
			modelPointer->initModel(linestr);
			modelPointer->FindSignSymmetries();
			found = true;
			break;
		}
		else if (linecounter > LineNumb)
			break;
		linecounter++;
		if (infile.eof())
			break;
	}
	infile.close();
	if (!found)
	{
		std::cout << "Line not found !\n";
		return -1;
	}

	if (TerminalOutput)
		modelPointer->write();
	//CALL: BSMPT-->Phasetransition
	if (TerminalOutput)
		std::cout << "PTFinder called..." << std::endl;
	auto EWPT = Minimizer::PTFinder_gen_all(modelPointer, 0, 300);
	//SFOEWPT FOUND
	if (EWPT.StatusFlag == Minimizer::MinimizerStatus::SUCCESS and C_PT * EWPT.Tc < EWPT.vc)
	{
		if (TerminalOutput)
			std::cout << "SFOEWPT found..." << std::endl;
		std::vector<double> vcritical, vbarrier;
		vcritical = EWPT.EWMinimum;
		double TC = EWPT.Tc;
		double vc = EWPT.vc;
		std::vector<double> MinimumPlane;
		//Find the minimum in the symmetric phase. For this minimise at T = Tc + 1
		std::vector<double> vevsymmetricSolution, checksym, startpoint;
		for (size_t i = 0; i < modelPointer->get_nVEV(); i++)
			startpoint.push_back(0.5 * vcritical.at(i));
		vevsymmetricSolution = Minimizer::Minimize_gen_all(modelPointer, TC + 1, checksym, startpoint);
		double vw = 0.1;
		if (TerminalOutput)
			std::cout << "Currently calculating:" << std::endl;
		for (double mchi = mchi_Start; mchi <= mchi_End; mchi += mchi_Stepsize)
		{
			if (TerminalOutput)
			{

				std::cout << "\rmchi = " << mchi << "\n";
				std::cout << "\rsource_flag = " << source_flag << "\n";
				std::cout << "\rLambdaQL = " << LambdaQL << "\n";
			}
			auto eta = EtaInterface.CalcEta(vw, vcritical, vevsymmetricSolution, TC, modelPointer, LambdaQL, mchi, source_flag);
			outfile << linestr << sep;
			outfile << TC << sep << vc << sep << vw << sep << EtaInterface.getLW() << sep << mchi << sep << LambdaQL << sep << source_flag;
			for (auto x : eta)
				outfile << sep << x;
			outfile << std::endl;
		} //END: vw loop
	}	  //END: SFOEWPT FOUND
	else
	{
		outfile << -1 << -1 << -1 << -1 << -1 << std::endl;
	} //NO SFOEWPT
	outfile.close();

	return EXIT_SUCCESS;
} //END: Try
catch (exception &e)
{
	std::cerr << e.what() << std::endl;
	return EXIT_FAILURE;
}
