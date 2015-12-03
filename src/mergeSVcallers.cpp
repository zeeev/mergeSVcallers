/*

This program was created at:  Tue Dec  1 11:51:06 2015
This program was created by:  Zev N. Kronenberg

Contact: zev.kronenber@gmail.com

Organization: University of washington\n Genome Sciences

The MIT License (MIT)

Copyright (c) <2015> <Zev N. Kronenberg>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.


*/

#include <string>
#include <iostream>
#include <math.h>
#include <cmath>
#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <unistd.h>
#include <algorithm>
#include "split.h"
#include "Fasta.h"
#include "../vcflib/include/Variant.h"

// sorts vcflib::Variant by start;
bool sortStart(vcflib::Variant * L, vcflib::Variant * R){

  if(L->sequenceName != R->sequenceName){
    std::cerr << "FATAL: trying to sort positions on different chromosomes." << std::endl;
    exit(1);
  }

  return (L->position < R->position);
}


struct options{
  std::vector<std::string>             files;
  std::vector<std::string>              tags;
  std::map<std::string, int>          seqids;
  std::map<std::string, int>            skip;
  std::vector<vcflib::VariantCallFile*> vcfs;
  std::string                         fastaF;
  FastaReference                      fastaH;

}globalOpts;

static const char *optString = "hf:t:a:";

//-------------------------------   OPTIONS   --------------------------------
int parseOpts(int argc, char** argv)
{
  int opt = 0;
  opt = getopt(argc, argv, optString);
  while(opt != -1){
    switch(opt){
    case 'a':
      {
	globalOpts.fastaF = optarg;
	globalOpts.fastaH.open(globalOpts.fastaF);
	std::cerr << "INFO: ref: " << globalOpts.fastaF << std::endl;
	
	for(vector<std::string>::iterator it = 
	      globalOpts.fastaH.index->sequenceNames.begin(); 
	    it != globalOpts.fastaH.index->sequenceNames.end(); it++){
	  
	  globalOpts.seqids[*it] = 1;
	}
	break;
      }
    case 'h':
      {
	break;
      }
    case 'f':
      {
	globalOpts.files = split(optarg, ",");
	for(std::vector<std::string>::iterator it = globalOpts.files.begin();
            it != globalOpts.files.end(); it++){
	  std::cerr << "INFO: file: " << *it << std::endl;
	  vcflib::VariantCallFile * vcf;
	  vcf =  new vcflib::VariantCallFile;
	  vcf->open(*it);
	  globalOpts.vcfs.push_back(vcf);
	}
	break;
      }
    case 't':
      {
	globalOpts.tags = split(optarg, ",");
	for(std::vector<std::string>::iterator it = globalOpts.tags.begin();
	    it != globalOpts.tags.end(); it++){
	  std::cerr << "INFO: tag: " << *it << std::endl; 
	}
	break;
      }
    case '?':
      {
	break;
      }
    }
    
    opt = getopt( argc, argv, optString ); 
  }
  return 1;
}

//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  : a vector of var objects

 Function does   : processes a chunk of the files

 Function returns:

*/

void manageLoopOverVar(std::vector<vcflib::Variant *> & data){

  std::vector<vcflib::Variant *> svBuffer;

  long int lastPos = data.front()->position;
  
  for(std::vector<vcflib::Variant *>::iterator it = data.begin(); 
      it != data.end(); it++){
    
    if(svBuffer.empty()){
      svBuffer.push_back(*it);
    }


  }


}

//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  : a seqid

 Function does   : processes a chunk of the files

 Function returns:

*/

void processChunk(std::string  seqid){
  
  std::vector<vcflib::Variant *> data;

  int index = 0;

  // loading the data into the vector
  for(std::vector<vcflib::VariantCallFile*>::iterator it = globalOpts.vcfs.begin(); 
      it != globalOpts.vcfs.end(); it++){  
   


    bool getNext = true;

    vcflib::Variant var(**it);

    if(!(*it)->setRegion(seqid)){
      std::cerr << "WARNING: could not set region: seqid: " << seqid << " file:  " << globalOpts.files[index] << std::endl;
      std::cerr << "INFO: Seqid might not be in file" << std::endl;
    }

    // looping over read and loading it into new vcflib pointer

    while (getNext) {

      if(!(*it)->getNextVariant(var)){
	getNext = false;
      }
      else{
	
	if(var.info["SVTYPE"].front() == "BND"){
	  continue;
	}

	vcflib::Variant * v  = new vcflib::Variant;
	*v = var;
	v->info["SVCALLER"].push_back(globalOpts.tags[index]);
	data.push_back(v);
      }
    }
    index += 1;
  }

  std::cerr << "INFO: sorting: seqid: " << seqid << std::endl;
  std::sort(data.begin(), data.end(), sortStart);
  
  manageLoopOverVar(data);

  // deleting objects
  for(std::vector<vcflib::Variant *>::iterator it = data.begin(); 
      it != data.end(); it++){
    delete *it;
  }
  




}


//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  : reads from globalOpts

 Function does   : validates that the VCF lines are okay;

 Function returns:

*/
void validate (void)
{

  int index = 0;


  for(std::vector<vcflib::VariantCallFile*>::iterator it = globalOpts.vcfs.begin(); it != globalOpts.vcfs.end(); it++){
  
    std::cerr << "INFO: validating: " << globalOpts.files[index] << std::endl;

    if (!(*it)->is_open()) {
      std::cerr << "FATAL: could not open one of the VCF files" << std::endl;
      exit(1);
    }
    vcflib::Variant var(**it);

    int nvar = 0;

    while ( (*it)->getNextVariant(var) ) {

      nvar += 1;
      
      if(var.info.find("SVTYPE") == var.info.end()){
	std::cerr << "FATAL: missing SVTYPE in: " << globalOpts.files[index] << " : " << var.originalLine
		  << std::endl;
        exit(1);
      }
      
      if(var.info["SVTYPE"].front() ==  "BND"){
	std::cerr << "WARNING: BND events are skipped: " << globalOpts.files[index] << std::endl;
	continue;
      }

      if(var.info.find("END") == var.info.end()){
	std::cerr << "FATAL: missing END in: " << globalOpts.files[index] << " : " << var.originalLine << std::endl;
	exit(1);
      }
      if(var.info.find("CIPOS") == var.info.end()){
	std::cerr << "FATAL: missing CIPOS in: " << globalOpts.files[index] << std::endl;
        exit(1);
      }
      if(var.info.find("SVTYPE") == var.info.end()){
	std::cerr << "FATAL: missing SVTYPE in: " << globalOpts.files[index] << std::endl;
        exit(1);
      }

      if(var.info.find("CIEND") == var.info.end()){
	std::cerr << "FATAL: missing CIEND in: " << globalOpts.files[index] << std::endl;
        exit(1);
      }

      if(nvar == 1000){
	break;
      }
    }

    std::cerr << "INFO: " << globalOpts.files[index] << " - VCF INFO: END   [x]" << std::endl;
    std::cerr << "INFO: " << globalOpts.files[index] << " - VCF INFO: CIPOS [x]" << std::endl;
    std::cerr << "INFO: " << globalOpts.files[index] << " - VCF INFO: CIEND [x]" << std::endl;

    index += 1;

  }
}
//-------------------------------    MAIN     --------------------------------
/*
 Comments:
*/

int main( int argc, char** argv)
{
  
  parseOpts(argc, argv);
  
  if(globalOpts.fastaF.empty()){
    std::cerr << "FATAL: required: -a - fasta file" << std::endl;
    exit(1);
  }    
  if(globalOpts.tags.empty()){
    std::cerr << "FATAL: required: -t - tags for input files" << std::endl;
    exit(1);
  }
  if(globalOpts.tags.size() != globalOpts.files.size()){
    std::cerr << "FATAL: The number of tags and files are different" << std::endl;
    exit(1);
  }

  validate();
  std::cerr << "INFO: All VCF files validated" << std::endl;

  for(std::map<std::string, int>::iterator it = globalOpts.seqids.begin(); 
      it != globalOpts.seqids.end(); it++){
    
    processChunk(it->first);
  }


  return 0;
}
