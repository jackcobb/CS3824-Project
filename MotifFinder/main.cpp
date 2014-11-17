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
    
    
    
    const double arr[] = { .5, .25, .125, .125 };
    const int size = sizeof( arr ) / sizeof ( *arr );
    std::vector<double> dist(arr, arr+size);
    RandomEnumerator enumerator = RandomEnumerator();
    enumerator.InitializeDistribution(dist);
    double count[] = {0, 0, 0, 0};
    for(int i = 0; i < 100000; i++)
    {
        count[enumerator.EnumerateRandomVar()]++;
    }
    cout << count[0] << " " << count[1] << " " << count[2] << " " << count[3];
}

