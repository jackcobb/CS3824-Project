/* 
 * File:   ScoreEngine.h
 * Author: Eric
 *
 * Created on November 11, 2014, 2:18 PM
 */

#ifndef SCOREENGINE_H
#define	SCOREENGINE_H
#include <vector>
#include "Nucleotide.h"
#include "IScoreEngine.h"
#include "IDnaRepository.h"

using namespace std;

class ScoreEngine : public IScoreEngine {
public:
    ScoreEngine(IDnaRepository& repo);
    void SetRepo(IDnaRepository& repo);
    double Score(vector<Nucleotide_t>& motif, vector<int>& starting_loci);
    virtual ~ScoreEngine();
private:
    IDnaRepository& DNA;
    //Indexed [Nucleotide][Position]
    vector<double> ProbabilityMatrix;
    void ResizeProbabilityMatrix(int newSize);
    bool ValidateStartingLoci(int motifSize, vector<int>& starting_loci);
    void UpdateProbabilityMatrix(vector<Nucleotide_t>& motif, vector<int>& starting_loci);
    double LogProductProbMatrix(vector<Nucleotide_t>& motif);
};

#endif	/* SCOREENGINE_H */

