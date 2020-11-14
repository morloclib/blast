#ifndef __CORE_HPP__
#define __CORE_HPP__

#define BUFFER_SIZE 8192

#include <cstdlib>
#include <vector>
#include <tuple>
#include <string>
#include <iostream>
#include <fstream>
#include <unistd.h>


std::string classify(std::vector<std::tuple<std::string,std::string>> xs){
    std::array<char, 256> counts;
    counts.fill(0);
    for(size_t i = 0; i < xs.size(); i++){
        for(size_t j = 0; j < std::get<1>(xs[i]).size(); j++){
            counts[std::get<1>(xs[i])[j]]++;
        }
    }

    size_t N = 0;
    for(size_t i = 0; i < counts.size(); i++){
        N += counts[i];
    }

    bool dnaish = ( counts['a'] + counts['A']
                  + counts['t'] + counts['T']
                  + counts['g'] + counts['G']
                  + counts['c'] + counts['C']
                  + counts['n'] + counts['N']) / N > 0.80;

    return (dnaish ? "nucl" : "prot");
}


void writeFasta(std::vector<std::tuple<std::string,std::string>> xs, std::string path){
    std::ofstream fastafile;
    myfile.open (path);
    for(size_t i = 0; i < xs.size(); i++){
        fastafile << ">" << std::get<0>(xs[i]) << '\n';
        fastafile << std::get<1>(xs[i]) << '\n';
    }
    myfile.close();
}


void mlc_makeblastdb(std::vector<std::tuple<std::string,std::string>> xs, std::string path) {
    std::string dbtype = classify(xs);
    writeFasta(xs, path);
    std::ostringstream cmd;
    cmd << "makeblastdb -dbtype " << dbtype << "-in " << path << "\n";
    std::system(cmd.str());
    return path;
}

template <class R>
R _blast(
  std::string program
  std::string path,
  std::vector<std::string> args,
  std::vector<std::tuple<std::string,std::string>> xs)
{
    FILE *blast_fp;
    int chars_read;
    // + 1 so that the string is always NULL terminated
    char buffer[BUFFER_SIZE + 1];
    memset(buffer, '\0', sizeof(buffer));

    std::ostringstream cmd;
    cmd << program;
    for(size_t i = 0; i < args.size(); i++){
        cmd << args[i] << " ";
    }
    cmd << "-outfmt 6";

    blast_fp = popen(cmd.str(), "w");
}

  // { qaccver  :: Str
  // , saccver  :: Str
  // , pident   :: Num
  // , length   :: Int
  // , mismatch :: Int
  // , gapopen  :: Int
  // , qstart   :: Int
  // , qend     :: Int
  // , sstart   :: Int
  // , send     :: Int
  // , evalue   :: Num
  // , bitscore :: Num

template <class R>
R mlc_blastp(
  std::string path,
  std::vector<std::string> args,
  std::vector<std::tuple<std::string,std::string>> xs)
{
    return _blast("blastp", path, args, xs);
}


// blastp :: BlastDB AaSeq  -> [Str] -> Bioseq AaSeq  -> BlastOutput6

// blastn :: BlastDB DnaSeq -> [Str] -> Bioseq DnaSeq -> BlastOutput6

// blastx :: BlastDB AaSeq  -> [Str] -> Bioseq DnaSeq -> BlastOutput6

// tblastn :: BlastDB DnaSeq -> [Str] -> Bioseq AaSeq  -> BlastOutput6

#endif
