/* 
 * File:   ScoreEngine.cpp
 * Author: Eric
 * 
 * Created on November 11, 2014, 2:18 PM
 */
#include <stdexcept>
#include <math.h>
#include "ScoreEngine.h"


ScoreEngine::ScoreEngine() {
    DNA = 0;
    ProbabilityMatrix = vector<double>();
}

ScoreEngine::ScoreEngine(IDnaRepository* repo) {
    SetRepo(repo);
    ProbabilityMatrix = vector<double>();
}

void ScoreEngine::SetRepo(IDnaRepository* repo){
    DNA = repo;
}

double ScoreEngine::Score(vector<Nucleotide_t>& motif, vector<int>& starting_loci) {
    int motifSize = motif.size();
        if(!ValidateStartingLoci(motifSize, starting_loci))
            throw std::invalid_argument("Starting loci are invalid against the motif provided.");
    int notNullMotifSize = motifSize;
    for(int i = 0; i < motif.size(); i++)
        if(motif[i] == DC)
            notNullMotifSize--;    
    ResizeProbabilityMatrix(notNullMotifSize);
    UpdateProbabilityMatrix(motif, starting_loci);
    return LogProductProbMatrix(motif);
}

void ScoreEngine::ResizeProbabilityMatrix(int newSize) 
{
    if(ProbabilityMatrix.size() < newSize)
    {
        ProbabilityMatrix.resize(newSize, 0);
    }
}

bool ScoreEngine::ValidateStartingLoci(int motifSize, vector<int>& starting_loci) {
    int loci_size = starting_loci.size();
    int dna_size = DNA->Size();
    if(loci_size != dna_size)
        return false;
    for(int i = 0; i < loci_size; i++)
    {
        if(starting_loci[i] > dna_size - motifSize)
            return false;
    }
    return true;
}

void ScoreEngine::UpdateProbabilityMatrix(vector<Nucleotide_t>& motif, vector<int>& starting_loci) {
    Nucleotide_t comparer = DC;
    int columnCount = 0;
    for(int motifPosition = 0, matrixPosition = 0; motifPosition < motif.size(); motifPosition++)
    {
        if(motif[motifPosition] != DC)
        {
            comparer = motif[motifPosition];
            columnCount = 0;
            for(int j = 0; j < DNA->Size(); j++)
            {
                if(DNA->Get(j, motifPosition + starting_loci[j]) == comparer)
                    columnCount++;
            }
            ProbabilityMatrix[matrixPosition] = columnCount / motif.size();
            matrixPosition++;        
        }
    }
    
}

double ScoreEngine::LogProductProbMatrix(vector<Nucleotide_t>& motif) {
    double glob_prob[]=  {
        DNA->Count(A) / DNA->Size(),
        DNA->Count(T) / DNA->Size(),
        DNA->Count(G) / DNA->Size(),
        DNA->Count(C) / DNA->Size()};
    double runningSum = 0;
    for(int motifIndex = 0, matrixIndex = 0; motifIndex < motif.size(); motifIndex++)
    {
        if(motif[motifIndex] != DC)
        {
            runningSum += log2(ProbabilityMatrix[matrixIndex] / glob_prob[motif[motifIndex]]);
        }
    }
    return runningSum;
}


ScoreEngine::~ScoreEngine(){
    delete &ProbabilityMatrix;
}


