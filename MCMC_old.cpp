#include <string>
#include <algorithm>
#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

List leapAndShift(rowvec rho, const int& L);
rowvec augPair(rowvec r, const mat& tc, const vec& constr);
double dist(const mat& R, const rowvec& rho, const std::string& metric);
double vecdist(const rowvec& R, const rowvec& rho, const std::string& metric);
double logZn(const int& n, const double& alfa, const std::string& metric, const vec& fit, const vec& Cdn, const vec& Di);
int calculate(int num);
int binomialCoeff(int n, int k);

// [[Rcpp::export]]
List MCMCfunction(mat R, List pars) {

  RNGScope scope;

  // Read in all the parameters
  int nmc = as<int>(pars["nmc"]); // number of MCMC iterations
  int thin = as<int>(pars["thin"]); // thinning parameter (for saving rho and R, if needed)
  int alphaJump = as<int>(pars["alphaJump"]); // thinning parameter (for updating alpha)
  double lambda = as<double>(pars["lambda"]);
  int L = as<int>(pars["L"]);
  double sdAlpha = as<double>(pars["sdAlpha"]);
  std::string metric = as<std::string>(pars["dist"]);
  vec fit = as<vec>(pars["fit"]);
  bool aug = as<bool>(pars["aug"]);
  vec Rmiss = as<vec>(pars["Rmiss"]);
  List Cd = as<List>(pars["Cd"]);


  // Import functions from R
  Function sample("sample");
  Function seq("seq");
  Function getElement("getElement");

  // prepare relevant variables
  int N = R.n_rows; // number of assessors
  int n = R.n_cols; // number of items
  int niter = nmc/thin; // how many iterations to save in the end (for rho and R)
  mat rho=zeros<mat>(niter,n); // latent ranks
  rowvec rhoold = as<rowvec>(pars["rho0"]); // initialize rho
  vec alpha = zeros<vec>(std::floor(nmc/alphaJump));// alpha
  alpha(0) = as<double>(pars["alpha0"]);// initialize alpha
  cube RaugSAVE=zeros(niter, N, n);      // latent ranks
  double accRho=0, accAlpha=0, accRaug=0; // acceptance rate

  // other stuff
  List leapShiftValue; // values returned by leap and shift
  uvec LSindices; // indices of the elements of rho that have changed in the leap-and-shift
  double u, ratio, alphaRatio, logAlpha, alphaProp, probForw, probBack;
  rowvec rhoProp;
  vec test, tmp;
  mat Raug=R, propRaug=R;
  uvec miss;
  int iterindex=0;
  double prevAlphaInd, currAlphaInd;

  // Initialize the ranks R (data augmentation, if needed)
  if(aug){
    propRaug = R;
    for(int i=0; i<N; i++) {
      miss = find(R.row(i) < 0); // which elements to augment
      tmp = trans(R.row(i));
      tmp(miss) = shuffle(linspace(Rmiss(i), n, n-Rmiss(i)+1));
      Raug.row(i) = trans(tmp);
    }
  }


  //create Di and Cdn for normalizing constant
  vec Di;
  vec Cdn;
  if((metric=="footrule")&(n<51)){
    if((n % 2) == 0){
      int m=n/2;
      Di=2*linspace(0,pow(m,2),pow(m,2)+1);
      Cdn=as<vec>(getElement(Cd,n));
    }else{
      int m=(n-1)/2;
      Di=2*linspace(0,pow(m,2)+m,pow(m,2)+m+1);
      Cdn=as<vec>(getElement(Cd,n));
    }}else if ((metric=="spearman")&(n<15)){
      int m= binomialCoeff(n+1,3);
      Di=2*linspace(0,m,m+1);
      Cdn=as<vec>(getElement(Cd,n));
    }else{
      Di=zeros(n);
      Cdn=zeros(n);
    }




  // MCMC
  for(int t=1; t<nmc; t++){
    prevAlphaInd = std::floor((t-1)/alphaJump);
    currAlphaInd = std::floor(t/alphaJump);

    // Update rho
    leapShiftValue = leapAndShift(rhoold, L);
    rhoProp = as<rowvec>(leapShiftValue["prop"]);   // Sample a rank proposal
    probForw = leapShiftValue["probF"];
    probBack = leapShiftValue["probB"];
    LSindices = as<uvec>(leapShiftValue["indices"]);
    ratio = -alpha(prevAlphaInd)/n*(dist(Raug.cols(LSindices), rhoProp.cols(LSindices),metric) - dist(Raug.cols(LSindices), rhoold.cols(LSindices), metric)) + log(probBack) - log(probForw);
    u = log(as_scalar(randu(1,1)));
    if(ratio > u){
      rhoold = rhoProp;
      accRho += 1;
    }

    // Update alpha
    if(t % alphaJump == 0) { // Sample an alpha proposal every alphaJump'th time:
      logAlpha = as_scalar(randn(1,1))*sdAlpha + log(alpha(prevAlphaInd)); // Sample an alpha proposal (normal on the log scale)
      alphaProp = exp(logAlpha);
      alphaRatio = (alpha(prevAlphaInd)-alphaProp)/n*dist(Raug,rhoold,metric) + lambda*(alpha(prevAlphaInd)-alphaProp) + N*(logZn(n, alpha(prevAlphaInd), metric,fit,Cdn,Di)- logZn(n, alphaProp, metric,fit,Cdn,Di)) + log(alphaProp) - log(alpha(prevAlphaInd));
      u = log(as_scalar(randu(1,1)));
      if(alphaRatio > u){
      alpha(currAlphaInd) = alphaProp;
      accAlpha += 1;
      } else {
      alpha(currAlphaInd) = alpha(prevAlphaInd);
      }
    }

      // Update the ranks R (data augmentation)
      if(aug){
      for(int i=0; i<N; i++) {
      miss = find(R.row(i) < 0); // which elements to augment
      tmp = trans(R.row(i));
      tmp(miss) = shuffle(linspace(Rmiss(i), n, n-Rmiss(i)+1));
      propRaug.row(i) = trans(tmp);
      u = log(as_scalar(randu(1,1)));
      ratio = -alpha(currAlphaInd)/n*(dist(propRaug.row(i), rhoold,metric) -dist(Raug.row(i), rhoold,metric));
      if(ratio > u){
      Raug.row(i) = propRaug.row(i);
      accRaug += 1;
      }
      }
      }

      // SAVE
      if((t+1) % thin == 0){
      rho.row(iterindex) = rhoold;// save ranks
      for(int i=0; i<n; i++){
      RaugSAVE.slice(i).row(iterindex) = trans(Raug.col(i));
      }
      iterindex++;
      }


  }


      return List::create(Named("rho")=rho,
      Named("alpha")=alpha,
      Named("Raug")=RaugSAVE,
      Named("accRho")=accRho/(nmc-1.0),
      Named("accAlpha")=accAlpha/(std::floor(nmc/alphaJump)-1.0),
      Named("accRaug")=accRaug/(nmc-1.0)/N);
}

    // [[Rcpp::export]]
    List MCMCpairwise(List ranking, List draws, List constr, mat R, List pars) {

    RNGScope scope;

    // Read in all the parameters
    int nmc = as<int>(pars["nmc"]); // number of MCMC iterations
    int thin = as<int>(pars["thin"]); // thinning parameter (for saving rho and R, if needed)
    int whenTies = as<int>(pars["whenTies"]); // thinning parameter (for randomizing ties)
    int alphaJump = as<int>(pars["alphaJump"]); // thinning parameter (for updating alpha)
    double lambda = as<double>(pars["lambda"]);
    int L = as<int>(pars["L"]);
    double sdAlpha = as<double>(pars["sdAlpha"]);
    std::string metric = as<std::string>(pars["dist"]);
    vec fit = as<vec>(pars["fit"]);
    List Cd = as<List>(pars["Cd"]);

    // Import functions from R
    Function generateTC("generateTC");
    Function generateR("generateR");
    Function getElement("getElement");

    // prepare relevant variables
    int n = R.n_cols;
    int N = R.n_rows;
    int niter = nmc/thin; // how many iterations to save in the end (for rho and R)
    mat rho=zeros<mat>(niter,n); // rho
    rowvec rhoold = as<rowvec>(pars["rho0"]);
    vec alpha = zeros<vec>(std::floor(nmc/alphaJump)); // alpha
    alpha(0) = as<double>(pars["alpha0"]); // initialize alpha
    cube RaugSAVE=zeros(niter, N, n);      // augmented latent ranks
    double accRho=0, accAlpha=0, accRaug=0; // acceptance probabilities

    // other stuff
    List leapShiftValue; // values returned from leap and shift
    uvec LSindices; // indices of the elements of rho that have changed in the leap-and-shift
    double probForw, probBack;
    double u, ratio, alphaRatio, logAlpha, alphaProp;
    rowvec rhoProp, augProp;
    mat Rprop = R; // set it equal to the initial value
    List tc;
    tc = generateTC(ranking, draws, n); // generate the first transitive closure
    double prevAlphaInd, currAlphaInd;
    int iterindex=0;

    //create Di and Cdn for normalizing constant
    vec Di;
    vec Cdn;
    if((metric=="footrule")&(n<51)){
    if((n % 2) == 0){
    //int m=2*pow(n/2,2);
    // Di = regspace( 0, 2, m );
    int m=n/2;
    Di=2*linspace(0,pow(m,2),pow(m,2)+1);
    Cdn=as<vec>(getElement(Cd,n));
    }else{
    //int m=2*pow((n-1)/2,2)+2*(n-1)/2;
    //Di = regspace( 0, 2, m );
    int m=(n-1)/2;
    Di=2*linspace(0,pow(m,2)+m,pow(m,2)+m+1);
    Cdn=as<vec>(getElement(Cd,n));
    }}else if ((metric=="spearman")&(n<15)){
    int m= binomialCoeff(n+1,3);
    Di=2*linspace(0,m,m+1);
    Cdn=as<vec>(getElement(Cd,n));
    }else{
    Di=zeros(n);
    Cdn=zeros(n);

    }


    // MCMC
    for(int t=1; t<nmc; t++){

    prevAlphaInd = std::floor((t-1)/alphaJump);
    currAlphaInd = std::floor(t/alphaJump);

    // Update rho
    leapShiftValue = leapAndShift(rhoold, L);
    rhoProp = as<rowvec>(leapShiftValue["prop"]);     // Sample a rank proposal
    probForw = leapShiftValue["probF"];
    probBack = leapShiftValue["probB"];
    u = log(as_scalar(randu(1,1)));
    LSindices = as<uvec>(leapShiftValue["indices"]);
    uvec tmp3(1); tmp3.fill(t-1);
    ratio = -alpha(prevAlphaInd)/n*(dist(R.cols(LSindices), rhoProp.cols(LSindices), metric) - dist(R.cols(LSindices), rhoold.cols(LSindices), metric)) + log(probBack) - log(probForw);
    if(ratio > u){
    rhoold = rhoProp;
    accRho += 1;
    }

    // Update alpha
    if(t % alphaJump == 0) {
    logAlpha = as_scalar(randn(1,1))*sdAlpha + log(alpha(prevAlphaInd)); // Sample an alpha proposal (normal on the log scale)
    alphaProp = exp(logAlpha);
    u = log(as_scalar(randu(1,1)));
    alphaRatio = (alpha(prevAlphaInd)-alphaProp)/n*dist(R,rhoold,metric) + lambda*(alpha(prevAlphaInd)-alphaProp) + N*(logZn(n, alpha(prevAlphaInd), metric,fit,Cdn,Di)- logZn(n, alphaProp, metric,fit,Cdn,Di))+ log(alphaProp) - log(alpha(prevAlphaInd));
    if(alphaRatio > u){
    alpha(currAlphaInd) = alphaProp;
    accAlpha += 1;
    } else {
    alpha(currAlphaInd) = alpha(prevAlphaInd);
    }

    }

    // Update the ranks R (data augmentation)
    for(int i=0; i<N; i++){
    Rprop.row(i) = augPair(R.row(i), tc[i], constr[i]);
    ratio = -alpha(currAlphaInd)/n*(dist(Rprop.row(i), rhoold,metric) -dist(R.row(i), rhoold,metric));
    u = log(as_scalar(randu(1,1)));
    if(ratio > u){
    R.row(i) = Rprop.row(i);
    accRaug += 1;
    }
    }

    // SAVE
    if((t+1) % thin == 0){
    rho.row(iterindex) = rhoold;// save ranks
    for(int i=0; i<n; i++){
    RaugSAVE.slice(i).row(iterindex) = trans(R.col(i));
    }
    iterindex++;
    }

    // randomize the ties
    if((t+1) % whenTies == 0){
    tc = generateTC(ranking, draws, n);
    R  = as<mat>(generateR(ranking, n));
    }

    }

    return List::create(Named("rho")=rho,
    Named("alpha")=alpha,
    Named("Raug")=RaugSAVE,
    Named("accRho")=accRho/(nmc-1.0),
    Named("accAlpha")=accAlpha/(std::floor(nmc/alphaJump)-1.0),
    Named("accRaug")=accRaug/(nmc-1.0)/N);

    }

    // [[Rcpp::export]]
    List MCMCclustering(mat R, List pars) {

    RNGScope scope;

    // Read in all the parameters
    std::string metric = as<std::string>(pars["dist"]);
    int nmc = as<int>(pars["nmc"]); // number of MCMC iterations
    int thin = as<int>(pars["thin"]); // thinning parameter (for saving rho and R, if needed)
    int alphaJump = as<int>(pars["alphaJump"]); // thinning parameter (for updating alpha)
    double lambda = as<double>(pars["lambda"]);
    int K = as<int>(pars["K"]);
    int L = as<int>(pars["L"]);
    double sdAlpha = as<double>(pars["sdAlpha"]);
    bool aug = as<bool>(pars["aug"]);
    vec Rmiss = as<vec>(pars["Rmiss"]);
    vec fit = as<vec>(pars["fit"]);
    List Cd = as<List>(pars["Cd"]);

    // Import functions from R
    Function sample("sample");
    Function rdirichlet("rdirichlet");
    Function rbeta("rbeta");
    Function seq("seq");
    Function getElement("getElement");

    // prepare relevant variables
    int N = R.n_rows; // number of assessors
    int n = R.n_cols; // number of items
    int niter = nmc/thin; // how many iterations to save in the end
    cube rho=zeros(niter,n,K); // rho
    mat rhoold(K,n);
    for(int i=0; i<K; i++){// Initialize rho randomly
    rhoold.row(i) = as<rowvec>(wrap(sample(n)));
    }
    mat alpha(std::floor(nmc/alphaJump), K); // alpha
    alpha.row(0) = sort(trans( randu(K) ));// initialize alpha
    mat prob(niter, K);      // cluster probabilities
    rowvec probold = 1.0/K*trans(ones(K));
    mat zeta(niter, N);      // cluster labels
    rowvec zetaold = as<rowvec>(wrap(sample(K,N,1))); // '1' means 'replace=TRUE';
    cube probAss=zeros(niter, N, K);      // assessor-specific cluster probabilities
    cube RaugSAVE=zeros(niter, N, n);      // latent ranks
    vec accRho=zeros(K); // rank acceptance probability
    vec accAlpha=zeros(K); // alpha acceptance probability

    // Initialize the ranks R (data augmentation, if needed)
    double accRaug=0;
    mat Raug=R, propRaug =R;
    uvec miss;
    vec mytmp;
    if(aug){
    for(int i=0; i<N; i++) {
    miss = find(R.row(i) < 0); // which elements to augment
    mytmp = trans(R.row(i));
    mytmp(miss) = shuffle(linspace(Rmiss(i), n, n-Rmiss(i)+1));
    Raug.row(i) = trans(mytmp);
    }
    }

    // other stuff
    int psi = N/K; // hyperparameter for probabilities: keep it high so that we avoid being stuck in local means
    mat Rclust; // samples in a given cluster
    List leapShiftValue; // values returned from leap and shift
    uvec LSindices; // indices of the elements of rho that have changed in the leap-and-shift
    mat rhotmp, newmat;
    mat updatezeta = zeros(K,N);
    rowvec rhoProp(n);
    vec probProp(K);
    vec contTab(K);
    uvec tmp;
    vec tauk, uu, wold = ones(K), wnew = ones(K), aa;
    uvec clustsel;
    double probForw, probBack, prevAlphaInd, currAlphaInd;
    double u, ratio, alphaRatio, logAlpha, alphaProp, sampleU;
    int n1, k, i, zetai, quale, somma, iterindex=0;

    //create Di and Cdn for normalizing constant
    vec Di;
    vec Cdn;
    if((metric=="footrule")&(n<51)){
    if((n % 2) == 0){
    //int m=2*pow(n/2,2);
    // Di = regspace( 0, 2, m );
    int m=n/2;
    Di=2*linspace(0,pow(m,2),pow(m,2)+1);
    Cdn=as<vec>(getElement(Cd,n));
    }else{
    //int m=2*pow((n-1)/2,2)+2*(n-1)/2;
    //Di = regspace( 0, 2, m );
    int m=(n-1)/2;
    Di=2*linspace(0,pow(m,2)+m,pow(m,2)+m+1);
    Cdn=as<vec>(getElement(Cd,n));
    }}else if ((metric=="spearman")&(n<15)){
    int m= binomialCoeff(n+1,3);
    Di=2*linspace(0,m,m+1);
    Cdn=as<vec>(getElement(Cd,n));
    }else{
    Di=zeros(n);
    Cdn=zeros(n);

    }


    // MCMC
    for(int t=1; t<nmc; t++){

    // Select current alpha indexes
    prevAlphaInd = std::floor((t-1)/alphaJump);
    currAlphaInd = std::floor(t/alphaJump);

    // Update probabilities (conjugate model)
    for(k=0; k<K; k++){
    tmp = find(zetaold == (k+1));
    contTab(k) = tmp.n_rows;
    }
    tauk = contTab + psi*ones(K);
    probold = trans(as<vec>(wrap(rdirichlet(1, tauk))));

    // Cycle over clusters
    for(k=0; k<K; k++){

    // Select cluster elements
    clustsel = find(zetaold == (k+1));
    n1 = clustsel.n_elem;

    if(n1 > 0){
    Rclust = Raug.rows(trans(clustsel)); // Find the submatrix corresponding to this cluster

    // Update rho
    leapShiftValue = leapAndShift(rhoold.row(k), L);
    rhoProp = as<rowvec>(leapShiftValue["prop"]);     // Sample a rank proposal
    probForw = leapShiftValue["probF"];
    probBack = leapShiftValue["probB"];
    u = log(as_scalar(randu(1,1)));
    LSindices = as<uvec>(leapShiftValue["indices"]);
    uvec tmp3(1); tmp3.fill(k);
    ratio = -alpha(prevAlphaInd,k)/n*(dist(Rclust.cols(LSindices), rhoProp.cols(LSindices), metric) -
    dist(Rclust.cols(LSindices), rhoold.submat(tmp3,LSindices), metric)) + log(probBack) - log(probForw);
    if(ratio > u){        // Accept or reject the rank
    rhoold.row(k) = rhoProp;
    accRho(k)++;
    }


    // Update alpha
    if(t % alphaJump == 0) {
    logAlpha = as_scalar(randn(1,1))*sdAlpha + log(alpha(prevAlphaInd,k)); // Sample an alpha proposal (normal on the log scale)
    alphaProp = exp(logAlpha);
    u = log(as_scalar(randu(1,1)));
    alphaRatio = (alpha(prevAlphaInd,k)-alphaProp)/n*dist(Rclust,rhoold.row(k),metric) + lambda*(alpha(prevAlphaInd,k)-alphaProp) + n1*(logZn(n, alpha(prevAlphaInd), metric,fit,Cdn,Di)- logZn(n, alphaProp, metric,fit,Cdn,Di))+ log(alphaProp) - log(alpha(prevAlphaInd,k));
    if(alphaRatio > u){
    alpha(currAlphaInd,k) = alphaProp;
    accAlpha(k)++;
    } else {
    alpha(currAlphaInd,k) = alpha(prevAlphaInd,k);
    }

    }

    } else { // if the cluster is empty...

    // Update rho
    leapShiftValue = leapAndShift(rhoold.row(k), L);
    rhoProp = as<rowvec>(leapShiftValue["prop"]);     // Sample a rank proposal
    probForw = leapShiftValue["probF"];
    probBack = leapShiftValue["probB"];
    u = log(as_scalar(randu(1,1)));
    ratio = log(probBack) - log(probForw); // no data, ratio depends on the prior only (not symmetrical...)
    if(ratio > u){        // Accept or reject the rank
    rhoold.row(k) = rhoProp;
    accRho(k)++;
    }

    // Update alpha
    if(t % alphaJump == 0) {
    logAlpha = as_scalar(randn(1,1))*sdAlpha + log(alpha(prevAlphaInd,k)); // Sample an alpha proposal (normal on the log scale)
    alphaProp = exp(logAlpha);
    u = log(as_scalar(randu(1,1)));
    alphaRatio = lambda*(alpha(prevAlphaInd,k)-alphaProp) + log(alphaProp) - log(alpha(prevAlphaInd,k));
    if(alphaRatio > u){
    alpha(currAlphaInd,k) = alphaProp;
    accAlpha(k)++;
    } else {
    alpha(currAlphaInd,k) = alpha(prevAlphaInd,k);
    }
    }

    }

    }

    // Update the cluster labels
    for(k=0; k<K; k++){
    for(i=0; i<N; i++){
    updatezeta(k,i) = probold(k)*exp(-alpha(currAlphaInd,k)/n*dist(Raug.row(i), rhoold.row(k),metric))/exp(logZn(n, alpha(currAlphaInd,k), metric,fit,Cdn,Di));
    }
    }
    for(i=0; i<N; i++){
    aa = updatezeta.col(i)/sum(updatezeta.col(i));
    sampleU = as_scalar(randu(1,1));
    quale = 0;
    k = 0;
    somma = 0;
    while(quale == 0){
    if( ( sampleU < as_scalar(aa(k))/(1-somma) ) || ( k == K-1 ) ) {
    quale = 1;
    }
    somma += as_scalar(aa(k));
    k++;
    }
    zetaold(i) = k;
    }

    // Update the ranks R (data augmentation)
    if(aug){
    for(int i=0; i<N; i++) {
    miss = find(R.row(i) < 0); // which elements to augment
    mytmp = trans(R.row(i));
    mytmp(miss) = shuffle(linspace(Rmiss(i), n, n-Rmiss(i)+1));
    propRaug.row(i) = trans(mytmp);
    zetai = zetaold(i) - 1;
    u = log(as_scalar(randu(1,1)));
    ratio = -alpha(currAlphaInd,zetai)/n*(dist(propRaug.row(i), rhoold.row(zetai),metric) -dist(Raug.row(i), rhoold.row(zetai),metric));
    if(ratio > u){
    Raug.row(i) = propRaug.row(i);
    accRaug += 1;
    }
    }
    }

    // SAVE
    if((t+1) % thin == 0){
    for(int i=0; i<K; i++){
    rho.slice(i).row(iterindex) = rhoold.row(i);// save ranks
    probAss.slice(i).row(iterindex) = updatezeta.row(i);// save assessor-specific probabilities
    }
    for(int i=0; i<n; i++){
    RaugSAVE.slice(i).row(iterindex) = trans(Raug.col(i));
    }
    prob.row(iterindex) = probold; // save probabilities
    zeta.row(iterindex) = zetaold; // save labels
    iterindex++;
    }

    }


    return List::create(Named("rho")=rho,
    Named("alpha")=alpha,
    Named("prob")=prob,
    Named("zeta")=zeta,
    Named("probAss")=probAss,
    Named("Raug")=RaugSAVE,
    Named("accAlpha")=accAlpha/(nmc-1.0),
    Named("accRho")=accRho/(nmc-1.0),
    Named("accRaug")=accRaug/(nmc-1.0)/N);
    }

    // [[Rcpp::export]]
    List MCMCclusteringPair(List ranking, List draws, List constr, mat R, List pars) {
    RNGScope scope;

    // Read in all the parameters
    std::string metric = as<std::string>(pars["dist"]);
    int nmc = as<int>(pars["nmc"]);
    int thin = as<int>(pars["thin"]);   // thinning parameter for saving
    int whenTies = as<int>(pars["whenTies"]);  // thinning parameter for randomizing ties to preferences
    int alphaJump = as<int>(pars["alphaJump"]); // thinning parameter (for updating alpha)
    double lambda = as<double>(pars["lambda"]);
    int K = as<int>(pars["K"]);
    int L = as<int>(pars["L"]);
    double sdAlpha = as<double>(pars["sdAlpha"]);
    vec fit = as<vec>(pars["fit"]);
    bool performPred = as<bool>(pars["performPred"]); // does prediction of missing preferences has to be performed?
    List Cd = as<List>(pars["Cd"]);

    // Import functions from R
    Function generateTC("generateTC");
    Function generateR("generateR");
    Function sample("sample");
    Function rdirichlet("rdirichlet");
    Function rbeta("rbeta");
    Function getElement("getElement");

    // prepare relevant parameters / variables
    int N = R.n_rows; // number of assessors
    int n = R.n_cols; // number of items
    int niter = nmc/thin; // how many iterations to save in the end
    cube rho=zeros(niter,n,K); // rho
    mat rhoold(K,n);
    for(int i=0; i<K; i++){// Initialize the first ranks randomly
    rhoold.row(i) = as<rowvec>(wrap(sample(n)));
    }
    mat alpha(std::floor(nmc/alphaJump), K); // alpha
    alpha.row(0) = sort(trans( randu(K) ));// initialize alpha
    mat prob(niter, K);      // cluster probabilities
    rowvec probold = 1.0/K*trans(ones(K));
    mat zeta(niter, N);      // cluster labels
    rowvec zetaold = as<rowvec>(wrap(sample(K,N,1))); // '1' means 'replace=TRUE';
    cube probAss=zeros(niter, N, K);      // matrix of assessor-specific cluster probabilities
    cube RaugSAVE=zeros(niter, N, n);      // version of Raug for saving purposes
    vec accRho=zeros(K); // rank acceptance ratio
    vec accAlpha=zeros(K); // alpha acceptance ratio

    // prepare augmented data
    double accRaug=0;
    mat Raug=R, propRaug=R;
    // generate the first transitive closure
    List tc;
    tc = generateTC(ranking, draws, n);

    // if we have to perform prediction, we need the discarded preferences:
    mat preferences = as<mat>(pars["preferences"]);    // matrix with the preference removed for each assessor
    vec PairTest = zeros(N);

    // other stuff
    List leapShiftValue; // values returned from leap and shift
    uvec LSindices; // indices of the elements of rho that have changed in the leap-and-shift
    int psi = N/K; // hyperparameter for probabilities: keep it high so that we avoid being stuck in local means
    mat Rclust; // samples in a given cluster
    mat rhotmp, newmat;
    mat updatezeta = zeros(K,N);
    rowvec rhoProp(n);
    vec probProp(K);
    vec contTab(K);
    uvec tmp;
    vec tauk, uu;
    uvec clustsel;
    double probForw, probBack, prevAlphaInd, currAlphaInd;
    double u, ratio, alphaRatio, logAlpha, alphaProp, sampleU;
    int n1,k,i,zetai, ind1, ind2;
    vec aa, wold = ones(K), wnew = ones(K);
    int quale, somma, iterindex=0;

    //create Di and Cdn for normalizing constant
    vec Di;
    vec Cdn;
    if((metric=="footrule")&(n<51)){
    if((n % 2) == 0){
    //int m=2*pow(n/2,2);
    // Di = regspace( 0, 2, m );
    int m=n/2;
    Di=2*linspace(0,pow(m,2),pow(m,2)+1);
    Cdn=as<vec>(getElement(Cd,n));
    }else{
    //int m=2*pow((n-1)/2,2)+2*(n-1)/2;
    //Di = regspace( 0, 2, m );
    int m=(n-1)/2;
    Di=2*linspace(0,pow(m,2)+m,pow(m,2)+m+1);
    Cdn=as<vec>(getElement(Cd,n));
    }}else if ((metric=="spearman")&(n<15)){
    int m= binomialCoeff(n+1,3);
    Di=2*linspace(0,m,m+1);
    Cdn=as<vec>(getElement(Cd,n));
    }else{
    Di=zeros(n);
    Cdn=zeros(n);

    }


    // MCMC
    for(int t=1; t<nmc; t++){

    // Select current alpha indexes
    prevAlphaInd = std::floor((t-1)/alphaJump);
    currAlphaInd = std::floor(t/alphaJump);

    // Update probabilities (conjugate model)
    for(k=0; k<K; k++){
    tmp = find(zetaold == (k+1));
    contTab(k) = tmp.n_rows;
    }
    tauk = contTab + psi*ones(K);
    probold = trans(as<vec>(wrap(rdirichlet(1, tauk))));

    // Cycle over clusters
    for(k=0; k<K; k++){

    // Select cluster elements
    clustsel = find(zetaold == (k+1));
    n1 = clustsel.n_elem;

    if(n1 > 0){

    Rclust = Raug.rows(trans(clustsel)); // Find the submatrix corresponding to this cluster

    // Update rho
    leapShiftValue = leapAndShift(rhoold.row(k), L);
    rhoProp = as<rowvec>(leapShiftValue["prop"]);     // Sample a rank proposal
    probForw = leapShiftValue["probF"];
    probBack = leapShiftValue["probB"];
    u = log(as_scalar(randu(1,1)));
    LSindices = as<uvec>(leapShiftValue["indices"]);
    uvec tmp3(1); tmp3.fill(k);
    ratio = -alpha(prevAlphaInd,k)/n*(dist(Rclust.cols(LSindices), rhoProp.cols(LSindices), metric) -
    dist(Rclust.cols(LSindices), rhoold.submat(tmp3,LSindices), metric)) + log(probBack) - log(probForw);
    if(ratio > u){        // Accept or reject the rank
    rhoold.row(k) = rhoProp;
    accRho(k)++;
    }

    // Update alpha
    if(t % alphaJump == 0) {
    logAlpha = as_scalar(randn(1,1))*sdAlpha + log(alpha(prevAlphaInd,k)); // Sample an alpha proposal (normal on the log scale)
    alphaProp = exp(logAlpha);
    u = log(as_scalar(randu(1,1)));
    alphaRatio = (alpha(prevAlphaInd,k)-alphaProp)/n*dist(Rclust,rhoold.row(k),metric) + lambda*(alpha(prevAlphaInd,k)-alphaProp) + n1*(logZn(n,alpha(prevAlphaInd,k), metric,fit,Cdn,Di)- logZn(n,alphaProp, metric,fit,Cdn,Di))+ log(alphaProp) - log(alpha(prevAlphaInd,k));
    if(alphaRatio > u){
    alpha(currAlphaInd,k) = alphaProp;
    accAlpha(k)++;
    } else {
    alpha(currAlphaInd,k) = alpha(prevAlphaInd,k);
    }

    }

    } else { // if the cluster is empty...

    // Update rho
    leapShiftValue = leapAndShift(rhoold.row(k), L);
    rhoProp = as<rowvec>(leapShiftValue["prop"]);     // Sample a rank proposal
    probForw = leapShiftValue["probF"];
    probBack = leapShiftValue["probB"];
    u = log(as_scalar(randu(1,1)));
    ratio = log(probBack) - log(probForw); // no data, ratio depends on the prior only (not symmetrical...)
    if(ratio > u){        // Accept or reject the rank
    rhoold.row(k) = rhoProp;
    accRho(k)++;
    }


    // Update alpha
    if(t % alphaJump == 0) {
    logAlpha = as_scalar(randn(1,1))*sdAlpha + log(alpha(prevAlphaInd,k)); // Sample an alpha proposal (normal on the log scale)
    alphaProp = exp(logAlpha);
    u = log(as_scalar(randu(1,1)));
    alphaRatio = lambda*(alpha(prevAlphaInd,k)-alphaProp) + log(alphaProp) - log(alpha(prevAlphaInd,k));
    if(alphaRatio > u){
    alpha(currAlphaInd,k) = alphaProp;
    accAlpha(k)++;
    } else {
    alpha(currAlphaInd,k) = alpha(prevAlphaInd,k);
    }
    }

    }

    }

    // Update the cluster labels
    for(k=0; k<K; k++){
    for(i=0; i<N; i++){
    updatezeta(k,i) = probold(k)*exp(-alpha(currAlphaInd,k)/n*dist(Raug.row(i), rhoold.row(k),metric))/exp(logZn(n,alpha(currAlphaInd,k), metric,fit,Cdn,Di));
    }
    }
    for(i=0; i<N; i++){
    aa = updatezeta.col(i)/sum(updatezeta.col(i));
    sampleU = as_scalar(randu(1,1));
    quale = 0;
    k = 0;
    somma = 0;
    while(quale == 0){
    if( ( sampleU < as_scalar(aa(k))/(1-somma) ) || ( k == K-1 ) ) {
    quale = 1;
    }
    somma += as_scalar(aa(k));
    k++;
    }
    zetaold(i) = k;
    }

    // Update the ranks R (data augmentation)
    for(i=0; i<N; i++){
    propRaug.row(i) = augPair(Raug.row(i), tc[i], constr[i]);
    zetai = zetaold(i) - 1;
    u = log(as_scalar(randu(1,1)));
    ratio = -alpha(currAlphaInd,zetai)/n*(dist(propRaug.row(i), rhoold.row(zetai),metric) -dist(Raug.row(i), rhoold.row(zetai),metric));
    if(ratio > u){
    Raug.row(i) = propRaug.row(i);
    accRaug += 1;
    }
    }

    // IF PREDICTION HAS TO BE PERFORMED: check the preference which has been discarded
    if(performPred){
    for(i=0; i<N; i++){
    ind1 = preferences(i,0) - 1; // preferred item in the removed preferences
    ind2 = preferences(i,1) - 1; // dispreferred item in the removed preferences
    if( Raug(i,ind1) <  Raug(i,ind2) ){ // if the items ranks are coherent with the removed preference: successful prediction
    PairTest(i) += 1;
    }
    }
    }

    // Randomize the ties
    if((t+1) % whenTies == 0){
    tc = generateTC(ranking, draws, n);
    Raug  = as<mat>(generateR(tc, n));
    }

    // SAVE
    if((t+1) % thin == 0){
    for(i=0; i<K; i++){
    rho.slice(i).row(iterindex) = rhoold.row(i);// save ranks
    probAss.slice(i).row(iterindex) = updatezeta.row(i);// save assessor-specific probabilities
    }
    for(i=0; i<n; i++){
    RaugSAVE.slice(i).row(iterindex) = trans(Raug.col(i));
    }
    prob.row(iterindex) = probold; // save probabilities
    zeta.row(iterindex) = zetaold; // save labels
    iterindex++;
    }

    }


    return List::create(Named("rho")=rho,
    Named("alpha")=alpha,
    Named("prob")=prob,
    Named("zeta")=zeta,
    Named("probAss")=probAss,
    Named("PairTest")=PairTest,
    Named("Raug")=RaugSAVE,
    Named("accAlpha")=accAlpha/(nmc-1.0),
    Named("accRaug")=accRaug/(nmc-1.0)/N,
    Named("accRho")=accRho/(nmc-1.0));
    }



    rowvec augPair(rowvec r, const mat& tc, const vec& constr){

    int n = r.n_elem, newRank, u;
    uvec tmpLeft, tmpRight;
    vec support;
    int lbound, rbound;
    IntegerVector tmp;
    rowvec rProp = r;

    u = as_scalar(randi(1, distr_param(1,n)));
    if(any(constr == u)){
    tmpLeft = find(tc.col(u-1)==1);
    if(tmpLeft.n_elem == 0){
    lbound = 0;
    } else {
    lbound = max(r(tmpLeft));
    }
    tmpRight = find(tc.row(u-1)==1);
    if(tmpRight.n_elem == 0){
    rbound = n+1;
    } else {
    rbound = min(r(tmpRight));
    }
    if((rbound -1) > (lbound +1) ){// this means that I can sample at least 1 new rank
    newRank = as_scalar(randi(1, distr_param(lbound+1, rbound-1)));
    } else {// this means that lbound+1=rbound-1, i.e., I have no ranks to sample.
    // Then, the new rank equals the old rank
    return r;
    }
    } else {
    // Item u is not constrained
    newRank = as_scalar(randi(1, distr_param(1, n)));
    }
    // Complete the leap step by assigning a new value to element u
    rProp(u-1) = newRank;
    if(!any(r==newRank)) {
    std::cout << "Error" << std::endl;
    }

    // Now do the shift step
    int deltaR = rProp(u-1) - r(u-1);
    int index;
    uvec check;
    if(deltaR > 0){
    for(int k=1; k<=deltaR; k++){
    check = find(r==r(u-1)+k);
    index = as<int>(wrap(check));
    rProp(index) -= 1;
    }
    } else if(deltaR < 0) {
    for(int k=-1; k>=deltaR; k--){
    check = find(r == r(u-1) + k);
    index = as<int>(wrap(check));
    rProp(index) += 1;
    }
    }
    return rProp;
    }

    double logZn(const int& n, const double& alfa, const std::string& metric, const vec& fit, const vec& Cdn, const vec& Di){
    double res=0;
    if(metric == "kendall"){
    for(int j=1; j<(n+1); j++){res += log((1-exp(-j*alfa/n))/(1-exp(-alfa/n)));
    }
    }else if(metric=="hamming"){
    for(int j=0; j<(n+1); j++){res += calculate(n)*exp(-alfa)*pow((exp(alfa/n)-1),j)/calculate(j);
    }
    res=log(res);
    }else if(metric == "cayley"){
    for(int j=1; j<n; j++){
    res += log(1+j*exp(-alfa/n));
    }
    }else if(metric == "footrule"){
    if (n<51){
    res = log(sum(trans(Cdn)*exp(-alfa*Di/n)));
    }else{
    int nterms = fit.n_elem;
    for(int i=0; i<nterms; i++){
    res += pow(alfa,i)*fit(i);
    }
    }
    }else if(metric == "spearman"){
    if(n<15){
    res = log(sum(trans(Cdn)*exp(-alfa*Di/n)));
    }else{
    int nterms = fit.n_elem;
    for(int i=0; i<nterms; i++){
    res += pow(alfa,i)*fit(i);
    }
    }
    }
    return res;
    }

    List leapAndShift(rowvec rho, const int& L){
    RNGScope scope;

    rowvec prop=rho;
    IntegerVector tmp;
    vec support, indices;
    int n = rho.n_elem, u, index;
    double deltaR, probF, probB, supportNew;

    // Leap step:
    // 1. sample u
    u = as_scalar(randi(1, distr_param(1,n)));
    // 2. compute the S set for sampling the new rank
    if((rho(u-1) > 1) & (rho(u-1) < n)){
    support = join_cols(linspace(std::max(1.0, rho(u-1)-L), rho(u-1)-1,std::min(rho(u-1)-1,static_cast<double>(L))),
    linspace(rho(u-1)+1,std::min(static_cast<double>(n), rho(u-1)+L),std::min(n-rho(u-1),static_cast<double>(L))));
    } else if(rho(u-1) == 1){
    support = linspace(rho(u-1)+1,std::min(static_cast<double>(n), rho(u-1)+L),std::min(n-rho(u-1),static_cast<double>(L)));
    } else if(rho(u-1) == n){
    support = linspace(std::max(1.0, rho(u-1)-L), rho(u-1)-1,std::min(rho(u-1)-1,static_cast<double>(L)));
    }
    // 3. assign a random element of the support set, this completes the leap step
    index = as_scalar(randi(1, distr_param(0, support.n_elem-1)));
    prop(u-1) = support(index);

    // Compute the associated transition probabilities (BEFORE THE SHIFT STEP, WHICH IS DETERMINISTIC --> EASIER)
    if(std::abs(prop(u-1)-rho(u-1))==1){
    // in this case the transition probabilities coincide! (and in fact for L = 1 the L&S is symmetric)
    supportNew = std::min(prop(u-1)-1,static_cast<double>(L)) + std::min(n-prop(u-1),static_cast<double>(L));
    probF = 1.0/(n*support.n_elem) + 1.0/(n*supportNew);
    probB = probF;
    }else{
    // P(proposed|current)
    probF = 1.0/(n*support.n_elem);
    // P(current|proposed)
    supportNew = std::min(prop(u-1)-1,static_cast<double>(L)) + std::min(n-prop(u-1),static_cast<double>(L));
    probB = 1.0/(n*supportNew);
    }

    // Shift step:
    deltaR = prop(u-1) - rho(u-1);
    indices = zeros(std::abs(deltaR)+1);
    indices[0] = u-1;
    if(deltaR > 0){
    for(int k=1; k<=deltaR; k++){
    index = as_scalar(find(rho == rho(u-1) + k));
    prop(index) -= 1;
    indices[k] = index;
    }
    } else if(deltaR < 0) {
    for(int k=-1; k>=deltaR; k--){
    index = as_scalar(find(rho == rho(u-1) + k));
    prop(index) += 1;
    indices[-(k)] = index;
    }
    }


    return List::create(Named("prop")=prop,
    Named("indices")=indices,
    Named("deltaR")=deltaR,
    Named("probF")=probF,
    Named("probB")=probB);
    }

    double dist(const mat& R, const rowvec& rho, const std::string& metric){
    int N = R.n_rows;
    int n = rho.n_elem;
    double totDist=0;
    if(metric == "footrule"){
    for(int i=0; i<N; i++){
    totDist += norm(R.row(i) - rho, 1);
    }
    } else if(metric == "kendall"){
    for(int i=0; i<N; i++){
    for(int u=0; u<n; u++){
    for(int t=0; t<u; t++){
    if((R(i,t) > R(i,u)) & (rho(t) < rho(u) ) || ((R(i,t) < R(i,u) ) & (rho(t) > rho(u)))) totDist +=1;
    }
    }
    }
    } else if(metric == "spearman") {
    for(int i=0; i<N; i++){
    totDist += pow(norm(R.row(i) - rho, 2),2.0);
    }
    }
    return totDist;
    }

    double vecdist(const rowvec& R, const rowvec& rho, const std::string& metric){
    int n = rho.n_elem;
    double totDist=0;
    if(metric == "footrule"){
    totDist = norm(R - rho, 1);
    } else if(metric == "kendall"){
    for(int u=0; u<n; u++){
    for(int t=0; t<u; t++){
    if((R(t) > R(u)) & (rho(t) < rho(u) ) || ((R(t) < R(u) ) & (rho(t) > rho(u)))) totDist +=1;
    }
    }
    } else if(metric == "spearman") {
    totDist += pow(norm(R - rho, 2),2.0);
    }
    return totDist;
    }

    int calculate (int num){
    int fac=1,k;
    for(int m=1; m<=num; num--)
    fac = fac * num;
    k = fac;
    return k;
    }//factorial

    int binomialCoeff(int n, int k){
    int res = 1;

    // Since C(n, k) = C(n, n-k)
    if ( k > n - k )
    k = n - k;

    // Calculate value of [n * (n-1) *---* (n-k+1)] / [k * (k-1) *----* 1]
    for (int i = 0; i < k; ++i)
    {
    res *= (n - i);
    res /= (i + 1);
    }

    return res;
    }//binomial
