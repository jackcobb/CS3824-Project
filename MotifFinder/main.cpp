/* 
 * File:   main.cpp
 * Author: Eric
 *
 * Created on October 28, 2014, 3:17 PM
 * David Was Here
 */

#include <cstdlib>
#include<iostream>
#include <fstream>  
#include "RandomEnumerator.h"
#include "DnaRepository.h"
#include "FastaParser.h"
#include "ScoreEngine.h"
#include "IScoreEngine.h"
#include "RandomizedSearchEngine.h"
using namespace std;

/*
 * 
 */
int main(int argc, char** argv) 
{
    srand (time(NULL));
    DnaRepository repo(512);
    FastaParser parser = FastaParser();
    ifstream stream("input.fasta", std::ifstream::in);
    
    if (stream)
    {
        parser.Parse(stream, repo);
    }
    
    RandomizedSearchEngine engine = RandomizedSearchEngine(repo, 2, 0);
    engine.Search(60);
    vector<Nucleotide_t> motif = engine.GetMotif();
    vector<int> loci = engine.GetStartingLoci();
    
    ScoreEngine scorer = ScoreEngine(repo);
    double bestScore = scorer.Score(motif, loci);
    
    cout << bestScore << "\n";
    for(int i = 0; i < motif.size(); i++)
    {
        cout << motif[i] ;
    }
    cout << "\n";
    
    for(int i = 0; i < loci.size(); i++)
    {
        cout << loci[i] << "\n";
    }
}

