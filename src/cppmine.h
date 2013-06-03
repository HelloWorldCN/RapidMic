#ifndef _CPPMINE_H
#define _CPPMINE_H

#include "mine.h"
#include <string>
#include <cstring>
#include <vector>

class MINE {
    
private:
    
    /* The C mine_parameter structure defined in "mine.h" */
    mine_parameter param;
    
public:
    
    /*
     MINE constructor.
     
     alpha is the exponent in B(n) = n^alpha and must be in (0,1], and
     c determines how many more clumps there will be than columns in
     every partition. c = 15 meaning that when trying to draw Gx grid
     lines on the x-axis, the algorithm will start with at most 15*Gx
     clumps. c must be > 0.
     
     The constructor throws an exception when the parameters are invalid.
     */
    MINE(double alpha=0.6, double c=15);
    
    /*
     MINE destructor.
     */
    ~MINE();
    
    
    
    /* Returns the Maximal Information Coefficient (MIC). */
    double get_mic();
    
    /* Returns the Maximum Asymmetry Score (MAS). */
    double get_mas();
    
    /* Returns the Maximum Edge Value (MEV). */
    double get_mev();
    
    /* Returns the Minimum Cell Number (MCN). */
    double get_mcn();
    /*
     Computes the maximum normalized mutual information scores between
     the variables x and y of length n.
     */
	int OnePairsAnalysis(double *x, double *y,int n);
	int AllPairsAnalysis(double **inData,int m,int n);
	int TwoSetsAnalysis(double **inDataSet,int m,int n,int betweenid);
	int MasterAnalysis(double **inData,int m,int n,int masterid );
	int run(int argc, char **argv);
    void printResult();
    bool exportResult();
    
    
private:
	ANALYSISSTYLES m_AnalysisStyle;
	
	mine_result_score *m_Results;

    int m_ResultsArrayLen;

	int m_inputDataRowNum;
	int m_inputDataColNum;
    int m_inputLabel;
    int m_outputLabel;
	int m_masterId;
	int m_betweenId;
	int m_onepair1,m_onepair2;
	std::string m_inputFileName;
	
	std::string m_outputFileName;
	
    std::vector<std::string> m_varRowNames;
    std::vector<std::string> m_varColNames;
	double **m_inputDataMatrix;

    
public:
	std::string InputFileName() const { return m_inputFileName; }
	void InputFileName(std::string val) { m_inputFileName = val; }
	std::string OutputFileName() const { return m_outputFileName; }
	void OutputFileName(std::string val) { m_outputFileName = val; }

private:
	bool readCSV(bool rowlabel,bool collabel);
	void releaseInputMatrix();
	bool parserArgs(int argc, char **argv);
};


#endif /* _CPPMINE_H */
