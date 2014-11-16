/* 
 * File:   FastaParser.cpp
 * Author: Eric
 * 
 * Created on November 3, 2014, 5:48 PM
 */

#include "FastaParser.h"
#include "DnaSequenceVector.h"

#include <iostream>
#include <vector>
#include <string>
using namespace std;

FastaParser::FastaParser() {
}

FastaParser::FastaParser(const FastaParser& orig) {
}

FastaParser::~FastaParser() {
}

void FastaParser::Parse(std::ifstream& file, IDnaRepository& repo)
{
    string line;
    
    file.open("input.fasta");
    
    if (file)
    {
        while (getline(file, line))
        {
            if (line[0] == '<')
            {
                continue;
            }
            else
            {
                DnaSequenceVector vector(512);
                Nucleotide_t a = A;
                Nucleotide_t c = C;
                Nucleotide_t g = G;
                Nucleotide_t t = T;
                
                for (int i = 0; i < line.length(); i++)
                {
                    if (line[i] == T)
                    {
                        vector.Push(t);
                    }
                    else if (line[i] == A)
                    {
                        vector.Push(a);
                    }
                    else if (line[i] == G)
                    {
                        vector.Push(g);
                    }
                    else if (line[i] == C)
                    {
                        vector.Push(c);
                    }
                }
                
                repo.Add(vector);
            }
            
        }
    }
    
    else
    {
        //file is bad
    }
}

