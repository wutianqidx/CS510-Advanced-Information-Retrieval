#include "hmm.hpp"
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <cstring>
#include <vector>

using namespace std;
Model::Model(int stateCount) 
  : N(stateCount) 
{
  AllocateMemory();
}

Model::Model(char *modelFile) 
{
  ifstream ifs(modelFile);
  if (ifs.fail()) {
    cerr << "Failed to open the model file:"<< modelFile << endl;
    exit(1);
  }
  ifs >> N; // the first is the number of states
  AllocateMemory();

  char str[1000];
  int entryCount;

  // load initial state prob
  ifs >> str >> entryCount;
  if (strcmp(str, "InitPr")) {
    cerr << "Error: InitPr expected in model file\n";
  }
  int i;
  for (i=0; i<N; i++) {
    I[i]=0;
  }
  int s;
  double pr;
  for (i=0; i<entryCount; i++) {
    ifs >> s >> pr;
    I[s]=pr;
  }

  // load output prob
  ifs >> str >> entryCount;
  if (strcmp(str, "OutputPr")) {
    cerr << "Error: OutputPr expected in model file\n";
  }
  int j;
  for (i=0; i<N; i++) {
    for (j=0; j<SYMNUM; j++) {
      B[j][i] = 0;
    }
  }
  char sym;
  for (i=0; i<entryCount; i++) {
    ifs >> s >> sym >> pr;
    B[sym-baseChar][s]=pr;
  }


  // load state transition prob
  ifs >> str >> entryCount;
  if (strcmp(str, "TransPr")) {
    cerr << "Error: TransPr expected in model file\n";
  }
  for (i=0; i<N; i++) {
    for (j=0; j<N; j++) {
      A[i][j] = 0;
    }
  }
  int s1;
  for (i=0; i<entryCount; i++) {
    ifs >> s >> s1 >> pr;
    A[s][s1]=pr;
  }
  
}


void Model::AllocateMemory()
{
  I = new double[N];
  A = new double*[N];
  B = new double*[SYMNUM];
  int i;
  for (i=0; i< N; i++) {
    A[i] = new double[N];
  }
  for (i=0; i< SYMNUM; i++) {
    B[i] = new double[N];
  }

  ICounter = new double[N];
  ACounter = new double*[N];
  BCounter = new double*[SYMNUM];
  for (i=0; i< N; i++) {
    ACounter[i] = new double[N];
  }
  for (i=0; i< SYMNUM; i++) {
    BCounter[i] = new double[N];
  }
  INorm = 0;
  ANorm = new double[N];
  BNorm = new double[N];
}


void Model::ResetCounter()
{
  int i, j;
  for (i=0; i< N; i++) {
    ICounter[i] = 0;
    
    for (j=0; j<N; j++ ) {
      ACounter[i][j] = 0;
    }
    ANorm[i] =BNorm[i] = 0;
    for (j=0; j<SYMNUM; j++) {
      BCounter[j][i] = 0;
    }
  }
  INorm = 0;
}


// Using the counts stored in the counters to estimate
// the HMM parameters
void Model::UpdateParameter()
{
  double smoothConstant = 0.00001; // smoothing constant
  int i, j;
  for (i=0; i< N; i++) {
    I[i] = (smoothConstant+ICounter[i])/(N*smoothConstant+INorm);
    
    for (j=0; j<N; j++ ) {
      A[i][j] = (smoothConstant+ACounter[i][j])/(N*smoothConstant+ANorm[i]);
    }
    for (j=0; j<SYMNUM; j++) {
      B[j][i] = (smoothConstant+BCounter[j][i])/(SYMNUM*smoothConstant+BNorm[i]);
    }
  }
}


Model::~Model() 
{ 
  int i;
  for (i=0; i< N; i++) {
    delete [] A[i];
    delete [] ACounter[i];
  }
  for (i=0; i< SYMNUM; i++) {
      delete [] B[i];
      delete [] BCounter[i];
  }
  delete [] A;
  delete [] B;
  delete [] I;
  delete [] ACounter;
  delete [] BCounter;
  delete [] ANorm;
  delete [] BNorm;
}


// Supervised training: You need to finish six assignment statements corresponding to
//           counting three different events.
void Model::CountSequence(char *seqFile)
{

  ifstream ifs(seqFile);
  if (ifs.fail()) {
    cerr << "can't open sequence file: "<< seqFile <<endl;
    exit(1);
  }

  // ======== You will need to 
  // - count how many times it starts with state i
  // - count how many times a particular transition happens
  // - count how many times a particular symbol would be generated from a particular state
  // Increase the corresponding counters as well as the normalizers
  // See below for directions on where to add code

  char c;  // to store the symbol
  int s; // to store the tagged state
  int prevState= -1;
  bool first=true;
  while (ifs >> c >> s) {
    int cIndex = c-baseChar; // convert the character to an index that can be used for BCounter.
    if (first) { // this is the very first observed symbol

      //=========== uncomment and finish the following two assignment statements to count 
      // the initial state, i.e., change ICounter and INorm appropriately to reflect that 
      // you've seen once that state s is a starting state.

       ICounter[s] += 1;
       INorm += 1;

      first = false;
    }

    // ============ Uncomment and finish the following two assignment statements to count 
    // the event that character c has been generated from state s. 
    // I.e., change BCounter and BNorm approriately. 

        BCounter[cIndex][s] += 1;
        BNorm[s] += 1;

    if (prevState>=0) { // the current symbol is not the first one; prevState has the previous state.

      // =============== Uncomment and finish the following two assignment statements to count
      // the state transition event, 
      // I.e., change ACounter and ANorm appropriately to reflect that 
      // you've seen a transition from prevState to s. 

       ACounter[prevState][s] += 1;
       ANorm[prevState] += 1;

    }
    prevState =s;
  }
    
}

void Model::RandomInit(int *sym, int seqLen)
{
  /* seed the random number generator */
  srand48(time(NULL));

  int i,j;
  double sum,sumI;
  sumI = 0;

  for (i = 0; i < N; i++) {

    /* initialize the I(state_t) vector */
    I[i] = 0.99/N + 0.01*drand48();
    sumI += I[i];
    /* initialize the A(state_t, state_t+1) matrix */
    sum = 0.0;
    for (j = 0; j < N; j++) {
      A[i][j] = 0.99 / N + 0.01 * drand48();
      sum += A[i][j];
    }
    for (j = 0; j < N; j++)
      A[i][j] /= sum;
    
    /* initialize the B(output,state) matrix */
    sum = 0.0;
    for (j=0; j<seqLen;j++) {
      B[sym[j]][i] = 0.99 / N + 0.01 * drand48();
    }
    for (j = 0; j<SYMNUM; j++) {
      // B[j][i] = 0.99 / N + 0.01 * drand48();
      sum += B[j][i];
    }
    for (j = 0; j<SYMNUM; j++)
      B[j][i] /= sum;
  }
  for (i = 0; i < N; i++) {
    I[i]  /= sumI;
  }
}


void Model::UniformInit(int *sym, int seqLen)
{

  int i,j;
  double sum,sumI;
  sumI = 0;

  for (i = 0; i < N; i++) {

    /* initialize the I(state_t) vector */
    I[i] = 1;
    sumI += I[i];
    /* initialize the A(state_t, state_t+1) matrix */
    sum = 0.0;
    for (j = 0; j < N; j++) {
      A[i][j] = 1;
      sum += A[i][j];
    }
    for (j = 0; j < N; j++)
      A[i][j] /= sum;
    
    /* initialize the B(output,state) matrix */
    sum = 0.0;

    for(j=0;j<seqLen; j++) {
      
      B[sym[j]][i] = 1;
    }

    for (j = 0; j<SYMNUM; j++) {
      sum += B[j][i];
    }
    for (j = 0; j<SYMNUM; j++)
      B[j][i] /= sum;
  }
  for (i = 0; i < N; i++) {
    I[i]  /= sumI;
  }
}


void Model::Save(char *modelFile)
{
  ofstream ofs(modelFile);
  ofs << N << endl;
  int i;
  ofs << "InitPr "<< N << endl;
  for (i=0; i<N; i++) {
    ofs << i << " "<< I[i] << endl;
  }

  int count;
  count = 0;
  int j;
  for (i=0; i<N; i++) {
    for (j=0; j<SYMNUM; j++) {
      if (B[j][i]>0) {
	count++;
      }
    }
  }
  ofs << "OutputPr "<< count << endl;
  for (i=0; i<N; i++) {
    for (j=0; j<SYMNUM; j++) {
      if (B[j][i]>0) {
	ofs << i << " " << (char)(j+baseChar) << " "<< B[j][i] << endl;
      }
    }
  }


  count = 0;
  for (i=0; i<N; i++) {
    for (j=0; j<N; j++) {
      if (A[i][j]>0) {
	count++;
      }
    }
  }
  ofs << "TransPr "<< count << endl;
  for (i=0; i<N; i++) {
    for (j=0; j<N; j++) {
      if (A[i][j]>0) {
	ofs << i << " " << j << " "<< A[i][j] << endl;
      }
    }
  }
  
  ofs.close();
}


// Decoding function: You need to complete two assignment statements

void Model::Decode(char *seqFile)
{
  // a node in the path
  struct PathNode {
    int s; // the state 
    PathNode *prev; // the previous state in transition
  };

  // information related to a state needed for Viterbi iteration
  struct TrellisNode {
    PathNode *path; // the most likely path that has led to this state
    double pr; // the probability of the "path"
    double newPr; // a buffer to store the probability for the new extended path
  };
  
  PathNode **curBest = new PathNode*[N]; 
  // A buffer for constructing the best path to each state

  TrellisNode *trellis = new TrellisNode[N]; // Current trellis information

  int i;
  ifstream seq(seqFile);
  if (seq.fail()) {
    cerr << "can't open sequence file : "<< seqFile << endl;
    exit(1);
  }

  char c;  // store the symbol
  seq >> c;
  int cIndex = c- baseChar;
  // use cIndex to access matrix B, e.g., B[cIndex][i] gives you
  // the probability of generating symbol c from state i
  for (i=0; i<N; i++) {
    //

    // ===== Uncomment and complete the following assignment statement to compute the 
    // probability of starting at state i, i.e., p(i) * p(c|s_i)
    // working on the logorithm domain to avoid overflow.
    
     trellis[i].pr = log(I[i]*B[cIndex][i]);


    trellis[i].path = new PathNode();  // construct a path node for this state
    trellis[i].path->s = i;  // record the state
    trellis[i].path->prev = NULL; // prev=NULL means that it's the very first node
  }


  int j;
  while (seq >> c) {
    int cIndex = c- baseChar;
    // use cIndex to access matrix B, e.g., B[cIndex][i] gives you
    // the probability of generating symbol c from state i


    // =========================================================
    // The code in this loop is to grow 
    // the best path for each state according to the Viterbi algorithm
    // Suppose the sequence already processed is <c1, c2, ..., ck>.
    // For the next character c, we will update each of
    // the trellis node i.e., trellis[i], so that it will hold 
    // the state transition path ending at state i that most likley 
    // "generates" the
    // sequence <c1, ..., ck, c>. Remember trellis[i].path was supposed 
    // to hold the most likely path for <c1, ..., ck> that ends at state i
    //  and trellis[i].pr was supposed to hold the log of probability of data given
    // the path.

    for (i=0; i<N; i++) { // We now consider a path ending at each of the N states
      int bestFromState; // temporary variable to store the best "from-state", i.e.,
      // the best previous path's ending state.
      double bestPr;  // temporary variable to help select the best path
      for (j=0; j<N; j++) { // consider all previous states
	double thisPr;

	// ================== uncomment and complete the following assignment statement so that
	// the variable thisPr would hold the log of probability of 
	// <c1,...,ck, c> given that we follow the path stored in
	// trellis[j].path for generating <c1,...,ck> and we go from
	// state j to state i to generate c.
	// Note that we are working on the logorithm domain to avoid overflow.

    thisPr = trellis[j].pr+log(A[j][i]*B[cIndex][i]);

	if (j==0 || thisPr > bestPr) { // keep track of the best "from state".
	  bestFromState = j;
	  bestPr = thisPr;
	}
      }


       // Now bestFromState can tell us "previous path" works the best with our current
      // state i, and bestPr tells us the probability of the data given this best path. 
      trellis[i].newPr = bestPr;
    // Note that the newPr field is needed for temporarily storing the 
    // probability for the new best path you constructed. (Other states may
    // still need the probability for the old path (stored in the field pr)
    // in order to compute the new probabilities for their new best paths.
      curBest[i] = new PathNode();
      curBest[i]->prev = trellis[bestFromState].path;
      curBest[i]->s = i;
   // Note that the auxilary array curBest is necessary for temporarily storing
    // the constructed best path, so that it won't mess up with the old best
    // path which may be needed for constructing new paths for other states.
    }
    // now update pr
    for (i=0; i<N; i++) {
      trellis[i].pr = trellis[i].newPr;  // update trellis[i].pr so that it
      // now has p(<c1,....,ck,c>|path_new); before it had p(<c1,...,ck>|path_old)
      trellis[i].path = curBest[i]; // record the new, extended path

      // The following code, if uncommented, would print out the path stored in trellis[i].path
      // This may be useful for debugging...
       cerr << "stop at "<< i << " :";
      PathNode *p = trellis[i].path;
      while (p) {
	cerr << " "<< p->s;
	p = p->prev;
      }
      cerr << " =>"<< trellis[i].pr << " "<< endl;
      
    }
    
  }


  // Now, we only need to take the best trellis, i.e., the one with the highest
  // probability (trellis[i].pr), and print out the path according to trellis[i].path.
  // ========= the following code in this function is complete ============
  seq.clear();
  seq.seekg(ios::beg);
  vector<int> bestPath;
  bestPath.clear();
  int bestTrellis;
  double bestTrellisPr = 0;
  for (int i = 0; i<N; i++ ){
    if (i==0 || (trellis[i].pr > bestTrellisPr)) {
      bestTrellis =i;
      bestTrellisPr = trellis[i].pr;
    }
  }
  PathNode *p = trellis[bestTrellis].path;
  while (p) {
    bestPath.push_back(p->s);
    p = p->prev;
  }
  vector<int>::reverse_iterator it;
  it=bestPath.rbegin(); 
  while (seq >> c)  {
    if (it!=bestPath.rend()) {
      cout << c << " "<< *it << endl;
      it++;
    } else {
      cerr << "mismatched path and sequence\n";
    }
  }
  
  delete [] trellis;
}

// Computing alphas: You need to complete an assignment statement 
void Model::ComputeAlpha(int *seq, int seqLength)
{

  // alpha[t][i] = prob. of generating observations up to time t and being in state i at time t
  // seq has the sequence of observed symbols, and is an array of 
  // length seqLength. It is 0-indexed, and seq[t] is the index of 
  // the symbols at time t, which can be used directly to access matrix B,
  // e.g., B[seq[t]][i].

  // Store the normalized alpha values in the alpha array and the normalizers 
  // in the eta array.

  int t=0;
  int i,j;
  eta[t]=0;

  // the following code computes alpha[0][i] for i =0, 1,...,N-1
  for (i=0; i<N; i++) {
    alpha[t][i] = I[i]*B[seq[t]][i];
    eta[t]+= alpha[t][i];  // compute the normalizer
  }
  eta[t] = 1.0/eta[t];
  for (i=0; i<N; i++) {
    alpha[t][i] *= eta[t]; //normalize all alpha values
  }
  
  t++;
  // the following code will compute alpha[t][i] for t=1,..., seqLength-1.
  while (t<seqLength) {
    eta[t]=0;
    for (i= 0; i<N; i++) {
      alpha[t][i] = 0;

      for (j=0; j<N; j++) {


	alpha[t][i] += B[seq[t]][i] * alpha[t-1][j] * A[j][i];

      }
      
      eta[t] += alpha[t][i]; // compute the normalizer for alpha[t][i], i=0,1,...,N-1.
    }
    eta[t] = 1.0/eta[t];
    for (i=0; i<N; i++) {
      alpha[t][i] *= eta[t]; // normalize alphas.
    }
    t++;
  }

}

// Computing betas: You need to complete an assignment statement
void Model::ComputeBeta(int *seq, int seqLength)
{

  // beta[t][i] = prob. of generating observations after time t given
  // being in state i at time t
  // seq has the sequence of observed symbols, and is an array of 
  // length seqLength. It is 0-indexed, and seq[t] is the index of 
  // the symbols at time t, which can be used directly to access matrix B,
  // e.g., B[seq[t]][i].

  // the eta array has the normalizers
  // Store the normalized beta values in the beta array.

  int i,j;
  int t=seqLength-1;
  // compute beta[seqLength-1][i] for i=0,1,...,N-1.
  for (i=0; i<N; i++) {
    beta[t][i] = 1;
  }

  t--;
  // the following code will compute beta[t][i] for t = seqLength-2,...,0.
  while (t>=0) {
    for (i= 0; i<N; i++) {
      beta[t][i] = 0;
      for (j=0; j<N; j++) {


        
        beta[t][i] += beta[t+1][j] * A[i][j] * B[seq[t+1]][j]; 
      }
    }

    for (i=0; i<N; i++) {
      beta[t][i] *= eta[t+1];  // normalize betas
      //      cerr << "beta: "<< t << " state: "<< i << " ="<<beta[t][i] << endl;

    }
    t--;
  }

}

// Accumulating counts for parameter updating
void Model::AccumulateCounts(int *seq, int seqLength)
{

  double **gamma = new double*[seqLength];
  int t;
  for (t=0; t< seqLength; t++ ) {
    gamma[t] = new double[N];
  }
  
  // The following code computes the gamma's based on 
  // the alpha's and beta's

  int i,j;
  double norm;
  for (t=0; t<seqLength; t++) {
    norm = 0;
    for (j=0; j<N; j++ ) {

      gamma[t][j] = alpha[t][j] * beta[t][j];

      norm += gamma[t][j]; // compute the denominator in the gamma updating formula.
    }
    norm = 1.0/norm;
    for (j=0; j<N; j++ ) {
      gamma[t][j] *= norm; // normalizing to get the true gamma values.
    }
  }

  
  double countInc;
  for (i=0; i<N; i++) {


    countInc = gamma[0][i];
    ICounter[i] += countInc; 
    INorm += countInc;
  }



  for (t=0; t<seqLength; t++) {
    for (j=0; j<N; j++) {


      countInc = gamma[t][j]; 
      BCounter[seq[t]][j] += countInc; 
      BNorm[j] += countInc;
    }
  }

  for (t=0; t<seqLength-1; t++) {
    for (i=0; i<N; i++) {
      for (j=0; j<N; j++) {


	countInc = (gamma[t][i] * A[i][j] * B[seq[t+1]][j] * eta[t+1] * beta[t+1][j]) / beta[t][i];

	ACounter[i][j] += countInc;
	ANorm[i] += countInc;
      }
    }
  }
  for (i=0; i< seqLength; i++ ) {
    delete [] gamma[i];
  } 
  delete [] gamma;


}



void Model::PrintMatrix(double **mat, int seqLen)
{
  int i, j;
  cout << "   ";
  for (i=0; i<seqLen; i++ ){
    cout << " "<< i;
  }
  cout << endl;
  for (i=0; i<N; i++) {
    cout << i << " ";
    for (j=0; j<seqLen; j++) {
      cout << mat[j][i] << " ";
    }
    cout << endl;
  }
}


void Model::PrintOutputMatrix()
{
  int i, j;
  cout << "   ";
  for (i=0; i<N; i++ ){
    cout << " "<< i;
  }
  cout << endl;
  for (i=0; i<SYMNUM; i++) {
    bool valid = false;
    for (j=0; j<N;j++) {
      if (B[i][j]>0) {
	valid=true;
      }
    }
    if (valid) {
      char c=i+baseChar;
      cout << c << " ";
      for (j=0; j<N; j++) {
	cout << B[i][j] << " ";
      }
      cout << endl;
    }
  }
}


double Model::UpdateHMM(int *data, int seqLength)
{
  ResetCounter();
  alpha = new double*[seqLength];
  beta = new double*[seqLength];

  int i;
  for (i=0; i< seqLength; i++ ) {
    alpha[i] = new double[N];
    beta[i] = new double[N];
  }
  eta = new double[seqLength];

  cout << "Transition matrix\n";
  PrintMatrix(A, N);
  cout << "Output matrix\n";
  PrintOutputMatrix();

  ComputeAlpha(data, seqLength);
  // cout << "Alpha\n";
  //  PrintMatrix(alpha, seqLength);

  ComputeBeta(data, seqLength);
  //  cout << "Beta\n"; 
  // PrintMatrix(beta, seqLength);

  // compute data likelihood
  double prData =0;
  int t;
  // now scale back
  for (t=0; t<seqLength; t++) {
    prData += -log(eta[t])/log(1.01);
  }  

  ResetCounter();  // Initialize all the counters and the normalizers

  AccumulateCounts(data, seqLength);

  UpdateParameter();

  for (i=0; i< seqLength; i++ ) {
    delete [] alpha[i];
    delete [] beta[i];

  }
  delete [] eta;

  delete [] alpha;
  delete [] beta;

  return (prData);
}


void Model::Train(char *seqFile)
{
  vector<char> buffer(0);
  ifstream ifs(seqFile);
  if (ifs.fail()) {
    cerr << "can't open sequence file for training: "<< seqFile << endl;
  }
  char c;
  while (ifs >> c) {
    buffer.push_back(c);
  }

  /* 
     We train for many epochs. As we train, we keep track of how well
     the current HMM models the data ( P(data | hmm) ). Additionally,
     we keep track of a history of how the model has fit the data in
     previous epochs (meanFit). As the HMM settles into a stable
     state, currentFit will asymptote to a particular value. The
     time-averaged meanFit will asymptote to this value also (though
     more slowly). When the currentFit becomes very similar to the
     meanFit, the model isn't changing, so we stop training.
  */

  int maxT = buffer.size();
  int *sym = new int[maxT];
  vector<char>::iterator it;
  int k=0;
  for (it=buffer.begin(); it!=buffer.end(); it++) {
    sym[k++]= (*it)- baseChar;
    //    cerr << "sym: "<< sym[k-1]<<endl;
  }
  
  // UniformInit(sym,maxT);
   RandomInit(sym,maxT);

  double currentFit, meanFit = -1e-80;
   for(int epoch = 0; fabs( (meanFit - currentFit) / meanFit ) > 1e-6; epoch++) {
  
    /* Train for many epochs */
    cerr << "\nEpoch " <<  epoch << endl;
    currentFit = UpdateHMM(sym, maxT);
    cerr << " log-likelihood = "<< currentFit << endl;
    if (fabs(currentFit)< 0.1) {
      break;
    }
    meanFit = 0.2 * currentFit + 0.8 * meanFit;
  }
  cerr << "############## Final likelihood: " << currentFit << endl;
}





