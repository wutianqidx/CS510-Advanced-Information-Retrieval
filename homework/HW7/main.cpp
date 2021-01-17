#include "hmm.hpp"
#include <iostream>
#include <cstring>
#include <cstdlib>
#include <fstream>

using namespace std;


int main(int argc, char *argv[])
{
  Model *m;
  char modFile[300];
  char seqFile[300];
  strcpy(modFile,"");
  strcpy(seqFile,"");
  int action=0; 
  int N =2;
  for (int i=1; i<argc; i++) {
    if (!strncmp("-train",argv[i],2)) {
      action = 1;
    } else if (!strncmp("-decode",argv[i],2)) {
      action = 2;
    } else if (!strncmp("-count",argv[i],2)) {
      action = 3;
    } else if (!strncmp("-m",argv[i],2)) {
      if (i+1 >= argc) {
	cerr << " a file name must follow option -m \n";
	exit(1);
      } else {
	strcpy(modFile, argv[i+1]);
	i++;
      }
    } else if (!strncmp("-s",argv[i],2)) {
      if (i+1 >= argc) {
	cerr << " a file name must follow option -s \n";
	exit(1);
      } else {
	strcpy(seqFile, argv[i+1]);
	i++;
      }
    } else if (!strncmp("-n",argv[i],2)) {
      if (i+1 >= argc) {
	cerr << " a number must follow option -n \n";
	exit(1);
      } else {
	N = atoi(argv[i+1]);
	i++;
      }
    }
  }
  if (action==1) {
    m = new Model(N);
    m->Train(seqFile);
    m->Save(modFile);
  } else if (action==2) {
    m = new Model(modFile);
    m->Decode(seqFile);

  } else if (action==3) { // count sequence
    m = new Model(N);
    m->ResetCounter();
    m->CountSequence(seqFile);
    m->UpdateParameter();
    m->Save(modFile);
  } else {
    cerr << "no action specified\n";
  }
  return 0;
}
