/* 
 * File:   main.cpp
 * Author: Eric
 *
 * Created on October 28, 2014, 3:17 PM
 * David Was Here
 */

#include <cstdlib>
#include <cstring>
#include<iostream>
#include <fstream>  
#include "RandomEnumerator.h"
#include "DnaRepository.h"
#include "FastaParser.h"
#include "ScoreEngine.h"
#include "IScoreEngine.h"
#include "RandomizedSearchEngine.h"
#include "Runner.h"
#include "GibbsSearchEngine.h"
using namespace std;

void parseArgs(int argc, char** argv, int& motifLength, int& dontCare, double& time, string& filePath)
{
    for(int i = 1; i < argc; i++)
    {
        if(strcmp(argv[i], "-k") == 0)
        {
            i++;
            motifLength = atoi(argv[i]);
        }
        else if(strcmp(argv[i], "-d") == 0)
        {
            i++;
            dontCare = atoi(argv[i]);
        }
        else if(strcmp(argv[i], "-t") == 0)
        {
            i++;
            time = atof(argv[i]);
        }
        else
        {
            filePath = argv[i];
        }
    }
}

/*
 * 
 */
int main(int argc, char** argv) 
{
    int motifLength = 6;
    int dontCare = 0;
    double runTime = 2;
    string filePath = "notProvided";
    parseArgs(argc, argv, motifLength, dontCare, runTime, filePath);
    if(filePath.compare("notProvided") == 0)
    {
        cout << "A file path must be provided as an arguement\n";
        return -1;
    }
    
    if(dontCare != 0 && dontCare > motifLength - 2)
    {
        cout << "Invalid number of dont cares provided.\n";
        cout << "Must be <= motifLength - 2 which is " << motifLength - 2 << ".\n";
        return -1;
    }
    srand (time(NULL));
    DnaRepository repo(512);
    FastaParser parser = FastaParser();
    ifstream stream(filePath.c_str(), std::ifstream::in);
    if (stream)
    {
        parser.Parse(stream, repo);
    }
    else
    {
        cout << "An error occured reading file at " << filePath << ".\nCheck that the path is valid.\n";
        return -1;
    }
    for(int i = 0; i < repo.Size(); i++)
    {
        if(motifLength > repo.Size(i))
        {
            cout << "Dna sequence " << i+1 << " is not long enough for a motif of length " << motifLength << ".\n";
            cout << "-k must be provided < " << repo.Size(i);
        }
    }    
    
    GibbsSearchEngine engine = GibbsSearchEngine(repo, motifLength, dontCare);
    Runner runner = Runner(repo, engine);
    runner.Run(runTime);
}

