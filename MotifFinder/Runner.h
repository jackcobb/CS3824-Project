/* 
 * File:   Runner.h
 * Author: Eric
 *
 * Created on November 16, 2014, 8:10 PM
 */

#ifndef RUNNER_H
#define	RUNNER_H

#include <string>
#include <vector>
#include "IDnaRepository.h"
#include "ScoreEngine.h"
#include "ISearchEngine.h"


class Runner {
public:
    Runner(IDnaRepository& dna, ISearchEngine& searchEngine);
    void Run(double runTime);
    void SetSearchEngine(ISearchEngine& newEngine);
    void SetRepo(IDnaRepository& newRepo);
private:
    ScoreEngine scoreEngine;
    IDnaRepository& repo;
    ISearchEngine& searchEngine;
    string MotifToString(vector<Nucleotide_t>& input);
};

#endif	/* RUNNER_H */

