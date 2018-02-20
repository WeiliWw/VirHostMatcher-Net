/***************************************************************
*  Copyright (C) 2016 Yang Lu <ylu465@usc.edu>
*  Computational and Molecular Biology, Department of Biological Science
*  University of Southern California, LA, CA 90089, USA
*
*  Related publication:
*  TBA
***************************************************************/
#ifndef _OUTPUT_H
#define	_OUTPUT_H

#include "utils.h"
#include "SimpleMatrix.h"

enum OUTPUT_TYPE { PLAIN, PHYLIP, CYTOSCAPE, MDS };

class OutputWriter
{
public:
	static OutputWriter *getInstance();
	void writeToFile(OUTPUT_TYPE arg_output_type, smat::Matrix<double>* arg_distMat, std::vector<std::string>* arg_nameVec1, std::vector<std::string>* arg_nameVec2, std::string str_arg_outputFileURL);

private:
	OutputWriter(){}
	static OutputWriter* instance;
};

#endif