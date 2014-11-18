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
    DnaRepository repo(512);
    FastaParser parser = FastaParser();
    ifstream stream("input.fasta", std::ifstream::in);
    
    if (stream)
    {
        parser.Parse(stream, repo);
    }
    
    RandomizedSearchEngine engine = RandomizedSearchEngine(repo, 1, 0);
    engine.Search(2);
    return 0;
}

