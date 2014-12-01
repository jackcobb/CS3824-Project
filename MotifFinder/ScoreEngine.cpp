/* 
 * File:   ScoreEngine.cpp
 * Author: Eric
 * 
 * Created on November 11, 2014, 2:18 PM
 */
#include <stdexcept>
#include <math.h>
#include "ScoreEngine.h"

ScoreEngine::ScoreEngine(IDnaRepository& repo) :DNA(repo) {
    ProbabilityMatrix = vector<double>();
}
void ScoreEngine::SetRepo(IDnaRepository& repo){
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
    int temp = DNA.Size();
    if(loci_size != temp)
        return false;
    for(int i = 0; i < loci_size; i++)
    {
        if(starting_loci[i] > DNA.Size(i) - motifSize)
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
            for(int j = 0; j < DNA.Size(); j++)
            {
                if(DNA.Get(j, motifPosition + starting_loci[j]) == comparer)
                    columnCount++;
            }
            ProbabilityMatrix[matrixPosition] = columnCount / (double) DNA.Size();
            matrixPosition++;        
        }
    }
    
}

vector<Nucleotide_t> ScoreEngine::optimizeDontCaresInMotif(vector<Nucleotide_t>& motif, int dontCares, vector<int>& loci) {
    if(dontCares <= 0)
    {
        return motif;
    }
    
    vector<vector<double> > indexToScore (motif.size() - 2, vector<double>(2, 0));
    double score;
    ResizeProbabilityMatrix(motif.size());
    UpdateProbabilityMatrix(motif, loci);
    int i;
    for(i = 1; i < motif.size() - 1; i++){
        score = scoreAtIndex(motif[i], i);
        indexToScore[i - 1][0] = i;
        indexToScore[i - 1][1] = score;
    }
    quickSort(indexToScore, 0, indexToScore.size() - 1);
    vector<Nucleotide_t> motifWithDc = motif;
    for(i = 0; i < dontCares; i++)
    {
        motifWithDc[indexToScore[i][0]] = DC;
    }
    return motifWithDc;
}

double ScoreEngine::scoreAtIndex(Nucleotide_t toScore, int index) {
    double a = DNA.Count(A);
    double t = DNA.Count(T);
    double g = DNA.Count(G);
    double c = DNA.Count(C);
    double count = DNA.Count();
    double glob_prob[]=  {
        a / count, t / count, g / count, c / count};
    return log2(ProbabilityMatrix[index] / glob_prob[toScore]);
}

void ScoreEngine::quickSort(vector<vector<double> >& elements, int lower, int upper) {
    if(lower < upper)
    {
        int p = partition(elements, lower, upper);
        quickSort(elements, lower, p - 1);
        quickSort(elements, p + 1, upper);        
    }
}


void ScoreEngine::swap(vector<vector<double> >& elements, int i, int j) {
    double tempIndex = elements[i][0];
    double tempScore = elements[i][1];
    elements[i][0] = elements[j][0];
    elements[i][1] = elements[j][1];
    elements[j][0] = tempIndex;
    elements[j][1] = tempScore;    
}

int ScoreEngine::partition(vector<vector<double> >& elements, int lower, int upper) {
    int pivotIndex = (lower + upper) / 2;
    double pivotValueIndex = elements[pivotIndex][0];
    double pivotValueScore = elements[pivotIndex][1];
    swap(elements, pivotIndex, upper);
    int storeIndex = lower;
    for(int i = lower; i < upper; i++)
    {
        if(elements[i][1] < pivotValueScore)
        {
            swap(elements, i, storeIndex);
            storeIndex++;
        }
    }
    swap(elements, storeIndex, upper);
    return storeIndex;
}



double ScoreEngine::LogProductProbMatrix(vector<Nucleotide_t>& motif) {
    double a = DNA.Count(A);
    double t = DNA.Count(T);
    double g = DNA.Count(G);
    double c = DNA.Count(C);
    double count = DNA.Count();
    double glob_prob[]=  {
        a / count, t / count, g / count, c / count};
    double runningSum = 0;
    for(int motifIndex = 0, matrixIndex = 0; motifIndex < motif.size(); motifIndex++)
    {
        if(motif[motifIndex] != DC)
        {
            runningSum += log2(ProbabilityMatrix[matrixIndex] /(double) glob_prob[motif[motifIndex]]);
            matrixIndex++;
        }
    }
    return runningSum;
}


ScoreEngine::~ScoreEngine(){
}


