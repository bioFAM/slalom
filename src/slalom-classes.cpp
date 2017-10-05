// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
//   Learn more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
using namespace Rcpp;

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(BH)]]

// one include file from Boost to access the digamma function
#include <boost/math/special_functions/digamma.hpp>


/*----------------------------------------------------------------------------*/
// SlalomModel class definition and module //

//' @title
//' SlalomModel C++ class
//' @description
//' A C++ class for SlalomModel models.
//'
//' @param Y_init matrix of expression values
//' @param pi_init G x K matrix with each entry being the prior
//' probability for a gene g being active for factor k.
//' @param X_init matrix of initial factor states (N x K)
//' @param W_init G x K matrix of initial weights
//' @param prior_alpha numeric vector of length two giving prior values
//' for the gamma hyperparameters of the precisions
//' @param prior_epsilon numeric vector of length two giving prior values
//' for the gamma hyperparameters of the residual variances
//'
//' @return
//' an object of the SlalomModel class
//'
//' @useDynLib slalom
//' @importFrom Rcpp evalCpp
//' @exportClass Rcpp_SlalomModel
class SlalomModel {
// base class for slalom holding models
 public:
    // declare necessary variables for alpha
    double alpha_pa;
    double alpha_pb;
    arma::vec alpha_a;
    arma::vec alpha_b;
    arma::vec alpha_E1;
    arma::vec alpha_lnE1;
    // declare necessary variables for epsilon
    double epsilon_pa;
    double epsilon_pb;
    arma::vec epsilon_a;
    arma::vec epsilon_b;
    arma::vec epsilon_E1;
    arma::vec epsilon_lnE1;
    arma::vec epsilon_diagSigmaS;
    // declare necessary variables for X
    arma::mat X_E1;
    arma::mat X_diagSigmaS;
    arma::mat X_init;
    // declare necessary variables for W
    arma::mat W_E1;
    arma::mat W_sigma2;
    arma::mat W_E2diag;
    arma::mat W_gamma0;
    arma::mat W_gamma1;
    // declare necessary variables for Z
    arma::mat Z_init;
    arma::mat Z_E1;
    // declare other variables
    arma::mat Pi_a;
    arma::mat Pi_b;
    arma::mat Pi_pa;
    arma::mat Pi_E1;
    arma::mat Y;
    arma::mat pseudo_Y;
    arma::vec YY;
    arma::mat Known;
    int K;
    int N;
    int G;
    int nAnnotated;
    int nHidden;
    int nKnown;
    int nIterations;
    int minIterations;
    int iterationCount;
    double tolerance;
    bool forceIterations;
    bool shuffle;
    double nScale;
    char noiseModel;
    double onF;
    arma::vec nOn;
    arma::vec iUnannotatedDense;
    arma::vec iUnannotatedSparse;
    arma::uvec doUpdate;
    bool dropFactors;
    bool learnPi;
    // for testing
    arma::vec SmTSk;
    arma::uvec setMinus;
    arma::mat tmp1;
    arma::mat tmp2;
    arma::mat tmp3;
    arma::mat tmp4;
    arma::mat tmp5;
    // methods
    void train(void);
    void update(void);
    void updateEpsilon(void);
    void updateAlpha(const int);
    void updateX(const int);
    void updateW(const int);
    void updatePi(const int);
    // constructor
    SlalomModel() {}
    SlalomModel(arma::mat Y_init, arma::mat pi_init, arma::mat X_init,
            arma::mat W_init, arma::vec prior_alpha, arma::vec prior_epsilon) :
        Y(Y_init), Pi_E1(pi_init), W_gamma0(pi_init), X_E1(X_init),
         W_E1(W_init) {
        K = pi_init.n_cols;
        N = Y_init.n_rows;
        G = Y_init.n_cols;
        nScale = 100.0;
        // initialise squared column sums
        YY = arma::zeros(G);
        for (int g=0; g < G; g++) {
            YY(g) = arma::accu(Y.col(g) % Y.col(g));
        }
        // initialise alpha variables
        alpha_pa = prior_alpha[0];
        alpha_pb = prior_alpha[1];
        alpha_a = arma::ones(K) * alpha_pa;
        alpha_b = arma::ones(K) * alpha_pb;
        alpha_E1 = alpha_b / alpha_a;
        alpha_lnE1 = arma::ones(K);
        for (int i = 0; i < K; i++) {
            alpha_lnE1(i) = boost::math::digamma(alpha_a(i));
        }
        alpha_lnE1 = alpha_lnE1 - arma::log(alpha_b);
        // initialise epsilon variables
        epsilon_pa = prior_epsilon[0];
        epsilon_pb = prior_epsilon[1];
        epsilon_a = arma::ones(G) * epsilon_pa;
        epsilon_b = arma::ones(G) * epsilon_pb;
        epsilon_E1 = epsilon_b / epsilon_a;
        epsilon_lnE1 = arma::ones(G);
        for (int i = 0; i < G; i++) {
            epsilon_lnE1(i) = boost::math::digamma(epsilon_a(i));
        }
        epsilon_lnE1 = epsilon_lnE1 - arma::log(epsilon_b);
        epsilon_diagSigmaS = arma::zeros(G);
        // initialise X variables
        X_diagSigmaS = arma::ones(N, K);
        // initialise W variables
        W_gamma1 = W_gamma0;
        W_sigma2 = arma::ones(G, K);
        W_E2diag = arma::zeros(G, K);
        // initialise iterations
        iterationCount = 0;
        tolerance = 1e-08;
        learnPi = true;
        doUpdate = arma::ones<arma::uvec>(K);
        minIterations = 1;
        nIterations = 1;
    }
};


/*----------------------------------------------------------------------------*/
// SlalomModel class train method

// train method
void SlalomModel::train(void) {
    /*
     Iterate updates of weights (with spike-and-slab prior), ARD parameters, factors, annd noise parameters.

     No arguments, but utilises the following parameters defined in the object:
     nIterations          (int): Number of iterations.
     forceIterations     (bool): Force the model to update `nIteration` times.
     tolerance          (float): Tolerance to monitor convergence of reconstruction error
     minIterations        (int): Minimum number of iterations the model should perform.
     */
    Rprintf( "     ....trying to train....     \n");
    arma::umat Ion = (this->W_gamma0 > .5);
    Rprintf( "have defined Ion\n");
    arma::mat Zr = this->X_E1 * (this->W_E1.t() % Ion.t());
    Rprintf( "have defined Zr\n");
    arma::mat Zd = this->Z_E1 - Zr;
    Rprintf( "have defined Zd\n");
    arma::rowvec error(1);
    arma::rowvec error_old(1);
    error.fill(arma::mean(arma::mean((arma::abs(Zd)))));
    double meanerr = arma::mean(error);
    Rprintf("Mean error : %f\n", meanerr);
    // mean absolute error
    error_old.fill(100);
    bool converged = false;
    Rprintf( "have defined error\n");
    // iterate model to train it
    this->iterationCount = 0;
    for (int iter = 0; iter <= this->nIterations; iter++) {
        // t = time.time();
        this->update();
        this->iterationCount++;
        if (iter % 100 == 0) {
            Rcout << "iteration " << iter << std::endl;
            // Rprintf( "iteration %d\\n", iter);
        }
        if (iter % 1 == 0) {
            Rcout << "iteration " << iter << std::endl;
            error_old = error;
            double meanW = arma::mean(arma::mean(this->W_E1));
            Rprintf("Mean W_E1 : %f\n", meanW);
            double meanX = arma::mean(arma::mean(this->X_E1));
            Rprintf("Mean X_E1 : %f\n", meanX);
            double meanAlpha = arma::mean(arma::mean(this->alpha_E1));
            Rprintf("Mean Alpha_E1 : %f\n", meanAlpha);
            double meanPi = arma::mean(arma::mean(this->Pi_E1));
            Rprintf("Mean Pi_E1 : %f\n", meanPi);
            double meanEpsilon = arma::mean(arma::mean(this->epsilon_E1));
            Rprintf("Mean Epsilon_E1 : %f\n", meanEpsilon);
            Zr = this->X_E1 * this->W_E1.t();
            Zd = this->Z_E1 - Zr;
            error.fill(arma::mean(arma::mean((arma::abs(Zd)))));
            // mean absolute error
            converged = arma::approx_equal(error_old, error, "absdiff",
                                           this->tolerance);
            double meanerr = arma::mean(error);
            Rprintf("Mean error : %f\n", meanerr);
            // converged = (arma::abs(error_old - error) < this->tolerance);
        }
        if ( converged && !(this->forceIterations) &&
             (iter > this->minIterations) ) {
            Rcout << "Model converged after " << iter << " iterations." << std::endl;
            // Rprintf( "Converged after %d iterations\n", iter);
            break;
        }
        this->Z_E1 = this->X_E1 * this->W_E1.t();
    }
    if (!converged) {
        Rcout << "Model not converged after " <<  this->nIterations << " iterations." << std::endl;
    }
}


/*----------------------------------------------------------------------------*/
// SlalomModel class update methods //

// wrapper around R's RNG such that we get a uniform distribution over
// [0,n) as required by the STL algorithm
inline int randWrapper(const int n) { return floor(unif_rand() * n); }

void SlalomModel::update(void) {
    /*
    Do one update of weights (with spike-and-slab prior), ARD parameters, factors, annd noise parameters.

    */
    this->epsilon_diagSigmaS = arma::zeros(this->K);
    arma::umat Ion = (this->W_gamma0 > .5);
    // check above: set parameter to zero with every update?
    arma::uvec kRange = arma::regspace<arma::uvec>(0,  (this->K - 1));
    if (this->shuffle == true && this->iterationCount > 0) {
        arma::uvec kunfix = arma::regspace<arma::uvec>(this->nKnown,
        (this->K - 1));
        std::random_shuffle(kunfix.begin(), kunfix.end(), randWrapper);
        for (int i = 0; i < (this->K - this->nKnown); ++i) {
            int ii = i + this->nKnown;
            kRange(ii) = kunfix(i);
        }
    }
    // switch factors off such that they can't be turned back on
    // Better syntax, use auto! automatically gets the right iterator type (C++11)
    // Regular iterator, non-C++11
    for (int m = 0; m < this->K; m++) {
        int k = kRange[m];
        if (this->doUpdate[k]) {
            if (this->dropFactors == false || this->iterationCount < 10 ||
            (this->alpha_E1[k] / arma::var(this->X_E1.col(k))) < 1e10) {
                this->updateW(k);
                if (this->learnPi) {
                    Rprintf("Update Pi");
                    if (arma::any(this->iUnannotatedSparse == k)) {
                        this->updatePi(k);
                    }
                }
                this->updateAlpha(k);
                this->updateX(k);
            } else {
                this->doUpdate[k] = 0;
                Rcout << "Switched off factor " << k << std::endl;
            }
        }
    }
    // const char noise_gauss = "gauss";
    // if (std::strcmp(this->noiseModel, noise_gauss) == 0) {
    //     // this->updateEpsilon();
    // };
    Rprintf("Update Epsilon");
    this->updateEpsilon();
}


void SlalomModel::updateAlpha(const int k) {
    // pdate alpha in a SlalomModel class object - precisions
    double Ewdwd = arma::accu(this->W_gamma0.col(k) % this->W_E2diag.col(k));
    // elementwise mult.
    this->alpha_a(k) = this->alpha_pa + 0.5 * Ewdwd;
    this->alpha_b(k) = this->alpha_pb + arma::accu(this->W_gamma0.col(k)) / 2.0;
    this->alpha_E1(k) = this->alpha_b(k) / this->alpha_a(k);
}


void SlalomModel::updatePi(const int k) {
    // update Pi in a SlalomModel class object
    this->Pi_a.col(k) = this->Pi_pa.col(k) + arma::accu(this->W_gamma0.col(k));
    this->Pi_b.col(k) = (this->Pi_pa.col(k) + this->G -
        arma::accu(this->W_gamma0.col(k)));
    this->Pi_E1.col(k) = this->Pi_a.col(k) / this->Pi_a.col(k) +
        this->Pi_b.col(k);
}


void SlalomModel::updateEpsilon(void) {
    // update Epsilon (vectorised) - noise parameters
    arma::uvec update_cols = this->doUpdate;
    arma::mat SW_sigma = (this->W_gamma0.cols(update_cols) %
                              this->W_E1.cols(update_cols));  // elementwise mult.; GxK mat
    arma::mat SW2_sigma = (this->W_gamma0.cols(update_cols) %
                               this->W_E2diag.cols(update_cols));  // elementwise mult.; GxK mat
    arma::mat muSTmuS = (this->X_E1.cols(update_cols).t() *
        this->X_E1.cols(update_cols));     // KxK matrix

    arma::mat newmat = this->Z_E1.t() * this->X_E1.cols(update_cols);
    arma::vec t1 = arma::sum(SW_sigma % newmat, 1);  // K length vec
    arma::vec t2 = arma::sum(SW2_sigma % arma::repmat(muSTmuS.diag().t() +
        this->epsilon_diagSigmaS(update_cols).t(), this->G, 1), 1);
    // set diagonals to zeros in muSTmuS for next calculation
    muSTmuS.diag().zeros();
    arma::vec t3 = arma::sum((SW_sigma * muSTmuS) % SW_sigma, 1);
    this->epsilon_E1 = (1.0 / (0.5 * (this->YY  + (-2 * t1  + t2 + t3))) /
        (0.5 * this->N));
    this->epsilon_a.fill(0.5 * this->N + this->epsilon_pa);
    this->epsilon_b = this->epsilon_pb + 0.5 * (this->YY + (-2 * t1  + t2 + t3));
    for (int i = 0; i < this->epsilon_E1.n_elem; i++) {
        if (this->epsilon_E1(i) > 1e6) {
            this->epsilon_E1(i) = 1e6;
        }
    }
}


void SlalomModel::updateX(const int k) {
    // update the factor states
    arma::uvec set1 = arma::regspace<arma::uvec>(0, 1, k - 1);   // zero-indexing
    arma::uvec set2 = arma::regspace<arma::uvec>(k + 1, 1, this->K - 1);
    arma::uvec setMinus = arma::join_cols(set1, set2);
    setMinus = setMinus(this->doUpdate(setMinus));
    arma::vec SW_sigma;
    arma::vec SW2_sigma;
    SW2_sigma = ((this->W_gamma0.col(k) % this->W_E2diag.col(k)) %
                     this->epsilon_E1);
    double alphaSm = arma::sum(SW2_sigma);
    for (int i = 0; i < this->N; i++) {
        this->X_diagSigmaS(i, k) = 1.0 / (1.0 + alphaSm);
    }
    if (k >= this->nKnown) {
        SW_sigma = (this->W_gamma0.col(k) %
                        this->W_E1.col(k)) % this->epsilon_E1;     // Gx1 vec
        arma::mat b0 = (this->X_E1.cols(setMinus) *
            (this->W_gamma0.cols(setMinus) % this->W_E1.cols(setMinus)).t());  // NxG mat.
        arma::vec b = b0 * SW_sigma;    // NxG x Gx1 -> Nx1 vec
        arma::vec barmuX = (this->Z_E1 * SW_sigma) - b;

        // update X
        this->X_E1.col(k) = barmuX / (1.0 + alphaSm);

        // keep diagSigmaS
        this->epsilon_diagSigmaS(k) = arma::accu(this->X_diagSigmaS.col(k));
    }
}


void SlalomModel::updateW(const int k) {
    // update the factor weights
    // define logPi values
    Rprintf("Factor (column) : %d\n", k);
    arma::vec logPi;
    int Muse = arma::accu(this->doUpdate);
    Rprintf("Muse : %d\n", Muse);
    if (k < this->nKnown || arma::any(this->iUnannotatedSparse == k) ||
        arma::any(this->iUnannotatedDense == k)) {
        logPi = arma::log(this->Pi_E1.col(k) / (1.0 - this->Pi_E1.col(k)));
    } else if (this->nScale > 0 && this->nScale < this->N) {
        logPi = arma::log(this->Pi_E1.col(k) / (1.0 - this->Pi_E1.col(k)));
        arma::uvec isOFF_ = arma::find(this->Pi_E1.col(k) < 0.5);
        arma::uvec kvec(1);
        kvec(0) = k;
        logPi(isOFF_) = ((this->N / this->nScale) *
                             arma::log(this->Pi_E1(isOFF_, kvec) /
                             (1 - this->Pi_E1(isOFF_, kvec))));

        arma::uvec isON_ = arma::find(this->Pi_E1.col(k) > 0.5);
        if (this->onF > 1.0) {
            logPi(isON_) = this->onF * arma::log(
                this->Pi_E1(isON_, kvec) / (1 - this->Pi_E1(isON_, kvec)));
        }
    } else {
        logPi = arma::log(this->Pi_E1.col(k) / (1 - this->Pi_E1.col(k)));
    }
    Rprintf("Mean logPi : %f\n", arma::mean(logPi));
    arma::vec sigma2Sigmaw;
    sigma2Sigmaw = (1.0 / this->epsilon_E1) * this->alpha_E1(k);
    this->tmp5 = sigma2Sigmaw;
    Rprintf("Mean sigma2Sigmaw : %f\n", arma::mean(sigma2Sigmaw));
    Rprintf("N sigma2Sigmaw < 0 : %d\n", arma::accu(sigma2Sigmaw < 0));
    arma::uvec set1 = arma::regspace<arma::uvec>(0, 1, k - 1);    // zero-indexing
    arma::uvec set2 = arma::regspace<arma::uvec>(k + 1, 1, this->K - 1);
    arma::uvec setMinus = arma::join_cols(set1, set2);
    setMinus = setMinus(this->doUpdate(setMinus));
    this->setMinus = setMinus;
    this->tmp1 = arma::repmat(this->X_E1.col(k), 1, Muse - 1);
    this->tmp2 = this->X_E1.cols(setMinus);
    arma::vec SmTSk = arma::sum(
        (arma::repmat(this->X_E1.col(k), 1, Muse - 1) %
             this->X_E1.cols(setMinus)), 0).t();
    this->SmTSk = SmTSk;
    // Rprintf("Mean SmTSk : %f\n", arma::mean(SmTSk));
    double tmp = arma::as_scalar(this->X_E1.col(k).t() * this->X_E1.col(k));
    double SmTSm = (tmp + arma::sum(this->X_diagSigmaS.col(k)));
    // Rprintf("SmTSm : %f\n", SmTSm);

    arma::vec b;
    arma::vec diff;
    b = (this->W_gamma0.cols(setMinus) % this->W_E1.cols(setMinus)) * SmTSk;
    this->tmp3 = b;
    diff = (this->X_E1.col(k).t() * this->Z_E1).t() - b;
    this->tmp4 = diff;
    // Rprintf("Mean diff : %f\n", arma::mean(diff));
    arma::vec diff2 = diff % diff;

    arma::vec SmTSmSig = SmTSm + sigma2Sigmaw;
    // Rprintf("Mean SmTSmSig : %f\n", arma::mean(SmTSmSig));
    // update gamma and W
    arma::vec u_qm = (logPi + 0.5 * arma::log(sigma2Sigmaw) - 0.5 *
                          arma::log(SmTSmSig) + (0.5 * this->epsilon_E1) %
                          (diff2 / SmTSmSig));
    this->W_gamma0.col(k) = 1.0 / (1 + arma::exp(-u_qm));
    this->W_gamma1.col(k) = 1.0 - this->W_gamma0.col(k);
    this->W_E1.col(k) = (diff / SmTSmSig);
    this->W_sigma2.col(k) = (1.0 / this->epsilon_E1) / SmTSmSig;
    // // check how to element-wise square a vector with Armadillo
    this->W_E2diag.col(k) = this->W_E1.col(k) % this->W_E1.col(k) + this->W_sigma2.col(k);
}


/*----------------------------------------------------------------------------*/
// SlalomModel Rcpp module

// Define a module to make the C++ class available to R
RCPP_MODULE(SlalomModel) {
    using namespace Rcpp;

    class_<SlalomModel>("SlalomModel")
        // expose the default constructor
        .constructor()
        .constructor<arma::mat, arma::mat, arma::mat, arma::mat, arma::vec,
        arma::vec>()
        //
        // fields
        .field("K", &SlalomModel::K)
        .field("N", &SlalomModel::N)
        .field("G", &SlalomModel::G)
        .field("nScale", &SlalomModel::nScale)
        .field("nAnnotated", &SlalomModel::nAnnotated)
        .field("nHidden", &SlalomModel::nHidden)
        .field("nKnown", &SlalomModel::nKnown)
        .field("nIterations", &SlalomModel::nIterations)
        .field("minIterations", &SlalomModel::minIterations)
        .field("iterationCount", &SlalomModel::iterationCount)
        .field("forceIterations", &SlalomModel::forceIterations)
        .field("tolerance", &SlalomModel::tolerance)
        .field("shuffle", &SlalomModel::shuffle)
        .field("noiseModel", &SlalomModel::noiseModel)
        .field("onF", &SlalomModel::onF)
        // alpha
        .field("alpha_pa", &SlalomModel::alpha_pa)
        .field("alpha_pb", &SlalomModel::alpha_pb)
        .field("alpha_a", &SlalomModel::alpha_a)
        .field("alpha_b", &SlalomModel::alpha_b)
        .field("alpha_E1", &SlalomModel::alpha_E1)
        .field("alpha_lnE1", &SlalomModel::alpha_lnE1)
        // epsilon
        .field("epsilon_pa", &SlalomModel::epsilon_pa)
        .field("epsilon_pb", &SlalomModel::epsilon_pb)
        .field("epsilon_a", &SlalomModel::epsilon_a)
        .field("epsilon_b", &SlalomModel::epsilon_b)
        .field("epsilon_E1", &SlalomModel::epsilon_E1)
        .field("epsilon_lnE1", &SlalomModel::epsilon_lnE1)
        .field("epsilon_diagSigmaS", &SlalomModel::epsilon_diagSigmaS)
        // X
        .field("X_E1", &SlalomModel::X_E1)
        .field("X_diagSigmaS", &SlalomModel::X_diagSigmaS)
        .field("X_init", &SlalomModel::X_init)
        // W
        .field("W_E1", &SlalomModel::W_E1)
        .field("W_sigma2", &SlalomModel::W_sigma2)
        .field("W_E2diag", &SlalomModel::W_E2diag)
        .field("W_gamma0", &SlalomModel::W_gamma0)
        .field("W_gamma1", &SlalomModel::W_gamma1)
        // Z
        .field("Z_E1", &SlalomModel::Z_E1)
        .field("Z_init", &SlalomModel::Z_init)
        // Pi
        .field("Pi_a", &SlalomModel::Pi_a)
        .field("Pi_pa", &SlalomModel::Pi_pa)
        .field("Pi_b", &SlalomModel::Pi_b)
        .field("Pi_E1", &SlalomModel::Pi_E1)
        // other variables
        .field("Known", &SlalomModel::Known)
        .field("Y", &SlalomModel::Y)
        .field("pseudo_Y", &SlalomModel::pseudo_Y)
        .field("YY", &SlalomModel::YY)
        .field("iUnannotatedDense", &SlalomModel::iUnannotatedDense)
        .field("iUnannotatedSparse", &SlalomModel::iUnannotatedSparse)
        .field("nOn", &SlalomModel::nOn)
        .field("doUpdate", &SlalomModel::doUpdate)
        .field("learnPi", &SlalomModel::learnPi)
        .field("dropFactors", &SlalomModel::dropFactors)
        // for testing
        .field("SmTSk", &SlalomModel::SmTSk)
        .field("setMinus", &SlalomModel::setMinus)
        .field("tmp1", &SlalomModel::tmp1)
        .field("tmp2", &SlalomModel::tmp2)
        .field("tmp3", &SlalomModel::tmp3)
        .field("tmp4", &SlalomModel::tmp4)
        .field("tmp4", &SlalomModel::tmp5)
        // methods
        .method("train", &SlalomModel::train , "Train the SlalomModel")
        .method("update", &SlalomModel::update , "Update the SlalomModel")
        .method("updateW", &SlalomModel::updateW , "Update W")
        .method("updateX", &SlalomModel::updateW , "Update X")
        .method("updatePi", &SlalomModel::updateW , "Update Pi")
        .method("updateEpsilon", &SlalomModel::updateW , "Update Epsilon")
        .method("updateAlpha", &SlalomModel::updateAlpha , "Update alpha")
        ;
}



