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
#include <sstream>

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
    Nucleotide_t a = A;
    Nucleotide_t c = C;
    Nucleotide_t g = G;
    Nucleotide_t t = T;
    DnaSequenceVector vector;
    //getline(file, line);
    while (getline(file, line))
    {
        if (line[0] == '>')
        {
            istringstream iss(line);       
            string size;
            iss >> size >> size >> size;
            printf("%s\n", size.c_str());
            //repo.Add(vector);
            vector = DnaSequenceVector(stoi(size));
            continue;
        }
        else
        {
            for (int i = 0; i < line.length(); i++)
            {
                if (line[i] == 'T')
                {
                    vector.Set(i, t);
                }
                else if (line[i] == 'A')
                {
                    vector.Set(i, a);
                }
                else if (line[i] == 'G')
                {
                    vector.Set(i, g);
                }
                else if (line[i] == 'C')
                {
                    vector.Set(i, c);
                }
            }
            repo.Add(vector);
        }
    }
            //repo.Add(vector);
}

