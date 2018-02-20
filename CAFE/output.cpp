#include "output.h"

OutputWriter* OutputWriter::instance = 0;

OutputWriter* OutputWriter::getInstance()
{
	if (!instance) instance = new OutputWriter();
	return instance;
}

void OutputWriter::writeToFile(OUTPUT_TYPE arg_output_type, smat::Matrix<double>* arg_distMat, std::vector<std::string>* arg_nameVec1, std::vector<std::string>* arg_nameVec2, std::string str_arg_outputFileURL)
{
	std::ofstream tmp_ofsPipe(str_arg_outputFileURL.c_str(), std::ofstream::out);

	// tmp_ofsPipe << ",";
	for (int j = 0; j<arg_nameVec2->size(); ++j)
		tmp_ofsPipe << "," << arg_nameVec2->at(j);
	tmp_ofsPipe << std::endl;

	for (int i = 0; i < arg_nameVec1->size(); ++i)
	{
		tmp_ofsPipe << arg_nameVec1->at(i);
		for (int j = 0; j<arg_nameVec2->size(); ++j)
			tmp_ofsPipe << "," << arg_distMat->get(i, j);
		tmp_ofsPipe << std::endl;
	}

	tmp_ofsPipe.close();
}

