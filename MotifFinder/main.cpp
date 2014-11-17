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
using namespace std;

/*
 * 
 */
int main(int argc, char** argv) 
{
    DnaRepository repo(512);
    FastaParser parser = FastaParser();
    ifstream stream("input.fasta", std::ifstream::in);
    
    if (stream)
    {
        parser.Parse(stream, repo);
    }
    
    ScoreEngine scorer = ScoreEngine(repo);
    vector<Nucleotide_t> motif = vector<Nucleotide_t>();
    motif.push_back(A);
    motif.push_back(DC);
    motif.push_back(A);
    vector<int> loci = vector<int>();
    loci.push_back(2);
    loci.push_back(2);
    
    double i = scorer.Score(motif, loci);
    return 0;
}

