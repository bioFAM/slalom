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

// ' @title
// ' SlalomModel C++ class
// ' @description
// ' A C++ class for SlalomModel models.
// '
// ' @param Y_init matrix of expression values
// ' @param pi_init G x K matrix with each entry being the prior
// ' probability for a gene g being active for factor k.
// ' @param X_init matrix of initial factor states (N x K)
// ' @param W_init G x K matrix of initial weights
// ' @param prior_alpha numeric vector of length two giving prior values for the gamma hyperparameters of the precisions
// ' @param prior_epsilon numeric vector of length two giving prior values for the gamma hyperparameters of the residual variances
// '
// ' @return
// ' an object of the SlalomModel class
// '
// ' @useDynLib slalom
// ' @importFrom Rcpp evalCpp
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
    arma::mat Known;
    arma::mat pi;
    arma::mat Y;
    arma::mat pseudo_Y;
    arma::vec YY;
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
    // methods
    void train(void);
    void update(void);
    void updateEpsilon(void);
    void updateAlpha(const int);
    void updateX(const int);
    void updateW(const int);
    // constructor
    SlalomModel () {};
    SlalomModel (arma::mat Y_init, arma::mat pi_init, arma::mat X_init,
            arma::mat W_init, arma::vec prior_alpha, arma::vec prior_epsilon) :
        Y(Y_init), pi(pi_init), W_gamma0(pi_init), X_E1(X_init), W_E1(W_init) {
        K = pi_init.n_cols;
        N = Y_init.n_rows;
        G = Y_init.n_cols;
        nScale = 100.0;
        // initialise squared column sums
        YY = arma::zeros(G);
        for (int g=0; g < G; g++) {
            YY(g) = arma::accu(Y.col(g) % Y.col(g));
        };
        // initialise alpha variables
        alpha_pa = prior_alpha[0];
        alpha_pb = prior_alpha[1];
        alpha_a = arma::ones(K) * alpha_pa;
        alpha_b = arma::ones(K) * alpha_pb;
        alpha_E1 = alpha_b / alpha_a;
        alpha_lnE1 = arma::ones(K);
        for (int i = 0; i < K; i++) {
            alpha_lnE1(i) = boost::math::digamma(alpha_a(i));
        };
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
        };
        epsilon_lnE1 = epsilon_lnE1 - arma::log(epsilon_b);
        epsilon_diagSigmaS = arma::zeros(G);
        // initialise X variables
        X_diagSigmaS = arma::ones(N, K);
        // initialise W variables
        W_gamma1 = W_gamma0;
        W_sigma2 = arma::ones(G, K);
        W_E2diag = arma::zeros(G, K);
    };

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
    arma::umat Ion = (this->W_gamma0 > .5);
    arma::mat Zr = this->X_E1 * (this->W_E1.t() % Ion.t());
    arma::mat Zd = this->Z_E1 - Zr;
    arma::rowvec error(1);
    arma::rowvec error_old(1);
    error.fill(arma::mean(arma::mean((arma::abs(Zd))))); // mean absolute error
    error_old.fill(100);
    bool converged = false;
    // iterate model to train it
    this->iterationCount = 0;
    for (int iter = 0; iter <= this->nIterations; iter++) {
        // t = time.time();
        this->update();
        this->iterationCount++;
        if (iter % 100 == 0) {
            Rcout << "iteration " << iter << std::endl;
            // Rprintf( "iteration %d\\n", iter);
        };
        if (iter % 50 == 0) {
            error_old = error;
            Zr = this->X_E1 * this->W_E1.t();
            Zd = this->Z_E1 - Zr;
            error.fill(arma::mean(arma::mean((arma::abs(Zd))))); // mean absolute error
            converged = arma::approx_equal(error_old, error, "absdiff",
                                           this->tolerance);
            // converged = (arma::abs(error_old - error) < this->tolerance);
        }
        if ( converged && !(this->forceIterations) &&
             (iter > this->minIterations) ) {
            Rcout << "Model converged after " << iter << " iterations." << std::endl;
            // Rprintf( "Converged after %d iterations\n", iter);
            break;
        };
        this->Z_E1 = this->X_E1 * this->W_E1.t();
    };
    if (!converged) {
        Rcout << "Model not converged after " <<  this->nIterations << " iterations." << std::endl;
    };
};


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
        arma::uvec kunfix = arma::regspace<arma::uvec>(this->nKnown, (this->K - 1));
        std::random_shuffle(kunfix.begin(), kunfix.end(), randWrapper);
        for(int i = 0; i < (this->K - this->nKnown); ++i) {
            int ii = i + this->nKnown;
            kRange(ii) = kunfix(i);
        };
    };
    // Better syntax, use auto! automatically gets the right iterator type (C++11)
    // Regular iterator, non-C++11
    for(int i = 0; i < this->K; i++) {
        int k = kRange(i); // * gets the factor idx out of the list
        this->updateW(k);
        this->updateAlpha(k);
        this->updateX(k);
    }
    // const char noise_gauss = "gauss";
    // if (std::strcmp(this->noiseModel, noise_gauss) == 0) {
    //     // this->updateEpsilon();
    // };
    this->updateEpsilon();
};

void SlalomModel::updateAlpha(const int k) {
    //update alpha in a SlalomModel class object - precisions
    double Ewdwd = arma::accu(this->W_gamma0.col(k) % this->W_E2diag.col(k)); // elementwise mult.
    this->alpha_a(k) = this->alpha_pa + 0.5 * Ewdwd;
    this->alpha_b(k) = this->alpha_pb + arma::accu(this->W_gamma0.col(k)) / 2.0;
    this->alpha_E1(k) = this->alpha_b(k) / this->alpha_a(k);
};


void SlalomModel::updateEpsilon(void) {
    //update Epsilon (vectorised) - noise parameters
    arma::mat SW_sigma = this->W_gamma0 % this->W_E1; // elementwise mult.; GxK mat
    arma::mat SW2_sigma = this->W_gamma0 % this->W_E2diag; // elementwise mult.; GxK mat

    arma::mat muSTmuS = this->X_E1.t() * this->X_E1; // KxK matrix

    arma::vec t1 = arma::sum(SW_sigma % (this->Z_E1.t() * this->X_E1), 1); // K length vec
    arma::vec t2 = arma::sum(SW2_sigma % arma::repmat(muSTmuS.diag().t() +
        this->epsilon_diagSigmaS.t(), this->G, 1), 1);
    // set diagonals to zeros in muSTmuS for next calculation
    muSTmuS.diag().zeros();
    arma::vec t3 = arma::sum((SW_sigma * muSTmuS) % SW_sigma, 1);

    this->epsilon_E1 = 1.0 / ((this->YY  + (-2 * t1  + t2 + t3)) / this->N);
    for (int i = 0; i < this->epsilon_E1.n_elem; i++) {
        if (this->epsilon_E1(i) > 1e6) {
            this->epsilon_E1(i) = 1e6;
        };
    };
};


void SlalomModel::updateX(const int k) {
    // update the factor states
    arma::uvec set1 = arma::regspace<arma::uvec>(0, 1, k - 1); // zero-indexing
    arma::uvec set2 = arma::regspace<arma::uvec>(k + 1, 1, this->K - 1);
    arma::uvec setMinus = arma::join_cols(set1, set2);
    arma::vec SW_sigma;
    arma::vec SW2_sigma;
    SW2_sigma = ((this->W_gamma0.col(k) % this->W_E2diag.col(k)) % this->epsilon_E1);
    double alphaSm = arma::sum(SW2_sigma);
    for (int i = 0; i < this->N; i++) {
        this->X_diagSigmaS(i,k) = 1.0 / (1.0 + alphaSm);
    }
    if (k >= this->nKnown) {
        SW_sigma = (this->W_gamma0.col(k) % this->W_E1.col(k)) % this->epsilon_E1; // Gx1 vec

        arma::mat b0 = (this->X_E1.cols(setMinus) *
            (this->W_gamma0.cols(setMinus) % this->W_E1.cols(setMinus)).t()); // NxG mat.
        arma::vec b = b0 * SW_sigma; // NxG x Gx1 -> Nx1 vec
        arma::vec barmuX = (this->Z_E1 * SW_sigma) - b;

        // update X
        this->X_E1.col(k) = barmuX / (1.0 + alphaSm);

        // keep diagSigmaS
        this->epsilon_diagSigmaS(k) = arma::accu(this->X_diagSigmaS.col(k));
    };
};


void SlalomModel::updateW(const int k) {
    // update the factor weights
    // define logPi values
    arma::vec logPi;
    if (arma::any(this->iUnannotatedSparse == k) || arma::any(this->iUnannotatedDense == k)) {
        logPi = arma::log(this->pi.col(k) / (1.0 - this->pi.col(k)));
    } else if (this->nScale > 0 && this->nScale < this->N) {
        logPi = arma::log(this->pi.col(k) / (1.0 - this->pi.col(k)));
        arma::uvec isOFF_ = arma::find(this->pi.col(k) < 0.5);
        arma::uvec kvec(1);
        kvec(0) = k;
        logPi(isOFF_) = ((this->N / this->nScale) *
                             arma::log(this->pi(isOFF_, kvec) /
                             (1 - this->pi(isOFF_, kvec))));

        arma::uvec isON_ = arma::find(this->pi.col(k) > 0.5);
        if (this->onF > 1.0) {
            logPi(isON_) = this->onF * arma::log(
                this->pi(isON_, kvec) / (1 - this->pi(isON_, kvec)));
        };
    } else {
        logPi = arma::log(this->pi.col(k) / (1 - this->pi.col(k)));
    }

    arma::vec sigma2Sigmaw;
    sigma2Sigmaw = (1.0 / this->epsilon_E1) * this->alpha_E1(k);

    arma::uvec set1 = arma::regspace<arma::uvec>(0, 1, k - 1); // zero-indexing
    arma::uvec set2 = arma::regspace<arma::uvec>(k + 1, 1, this->K - 1);
    arma::uvec setMinus = arma::join_cols(set1, set2);
    arma::vec SmTSk = arma::sum(
        (arma::repmat(this->X_E1.col(k), 1, this->K - 1) %
             this->X_E1.cols(setMinus)), 0).t();
    arma::mat tmp = (this->X_E1.col(k).t() * this->X_E1.col(k));
    double SmTSm = (tmp(0, 0) + arma::sum(this->X_diagSigmaS.col(k)));

    arma::vec b;
    arma::vec diff;
    b = (this->W_gamma0.cols(setMinus) % this->W_E1.cols(setMinus)) * SmTSk;
    diff = (this->X_E1.col(k).t() * this->Y).t() - b;
    arma::vec diff2 = diff % diff;

    arma::vec SmTSmSig = SmTSm + sigma2Sigmaw;

    //update gamma and W
    arma::vec u_qm = (logPi + 0.5 * arma::log(sigma2Sigmaw) - 0.5 *
                          arma::log(SmTSmSig) + (0.5 * this->epsilon_E1) %
                          ( diff2 / SmTSmSig));
    this->W_gamma0.col(k) = 1.0 / (1 + arma::exp(-u_qm));
    this->W_gamma1.col(k) = 1.0 - this->W_gamma0.col(k);
    this->W_E1.col(k) = (diff / SmTSmSig);
    this->W_sigma2.col(k) = (1.0 / this->epsilon_E1) / SmTSmSig;
    // // check how to element-wise square a vector with Armadillo
    this->W_E2diag.col(k) = this->W_E1.col(k) % this->W_E1.col(k) + this->W_sigma2.col(k);
};


/*----------------------------------------------------------------------------*/
// SlalomModel Rcpp module

// Define a module to make the C++ class available to R
RCPP_MODULE(SlalomModel) {
    using namespace Rcpp;

    class_<SlalomModel>("SlalomModel")
        // expose the default constructor
        .constructor()
        .constructor<arma::mat,arma::mat,arma::mat,arma::mat,arma::vec,arma::vec>()
        //
        // fields
        .field( "K", &SlalomModel::K )
        .field( "N", &SlalomModel::N )
        .field( "G", &SlalomModel::G )
        .field( "nScale", &SlalomModel::nScale)
        .field( "nAnnotated", &SlalomModel::nAnnotated )
        .field( "nHidden", &SlalomModel::nHidden )
        .field( "nKnown", &SlalomModel::nKnown )
        .field( "nIterations", &SlalomModel::nIterations )
        .field( "minIterations", &SlalomModel::nIterations )
        .field( "iterationCount", &SlalomModel::iterationCount )
        .field( "forceIterations", &SlalomModel::forceIterations )
        .field( "tolerance", &SlalomModel::tolerance )
        .field( "shuffle", &SlalomModel::shuffle )
        .field( "noiseModel", &SlalomModel::noiseModel )
        .field( "onF", &SlalomModel::onF )
        // alpha
        .field( "alpha_pa", &SlalomModel::alpha_pa )
        .field( "alpha_pb", &SlalomModel::alpha_pb )
        .field( "alpha_a", &SlalomModel::alpha_a )
        .field( "alpha_b", &SlalomModel::alpha_b )
        .field( "alpha_E1", &SlalomModel::alpha_E1 )
        .field( "alpha_lnE1", &SlalomModel::alpha_lnE1 )
        // epsilon
        .field( "epsilon_pa", &SlalomModel::epsilon_pa )
        .field( "epsilon_pb", &SlalomModel::epsilon_pb )
        .field( "epsilon_a", &SlalomModel::epsilon_a )
        .field( "epsilon_b", &SlalomModel::epsilon_b )
        .field( "epsilon_E1", &SlalomModel::epsilon_E1 )
        .field( "epsilon_lnE1", &SlalomModel::epsilon_lnE1 )
        .field( "epsilon_diagSigmaS", &SlalomModel::epsilon_diagSigmaS )
        // X
        .field( "X_E1", &SlalomModel::X_E1 )
        .field( "X_diagSigmaS", &SlalomModel::X_diagSigmaS )
        .field( "X_init", &SlalomModel::X_init )
        // W
        .field( "W_E1", &SlalomModel::W_E1 )
        .field( "W_sigma2", &SlalomModel::W_sigma2 )
        .field( "W_E2diag", &SlalomModel::W_E2diag )
        .field( "W_gamma0", &SlalomModel::W_gamma0 )
        .field( "W_gamma1", &SlalomModel::W_gamma1 )
        // Z
        .field( "Z_E1", &SlalomModel::Z_E1 )
        .field( "Z_init", &SlalomModel::Z_init )
        // other variables
        .field( "Known", &SlalomModel::Known )
        .field( "pi", &SlalomModel::pi )
        .field( "Y", &SlalomModel::Y )
        .field( "pseudo_Y", &SlalomModel::pseudo_Y )
        .field( "YY", &SlalomModel::YY )
        .field( "iUnannotatedDense", &SlalomModel::iUnannotatedDense )
        .field( "iUnannotatedSparse", &SlalomModel::iUnannotatedSparse )
        .field( "nOn", &SlalomModel::nOn )

        .method("train", &SlalomModel::train , "Train the SlalomModel")
        .method("update", &SlalomModel::update , "Update the SlalomModel")
        // .method("updateAlpha", &SlalomModel::updateAlpha , "Update the SlalomModel")
        ;
}



/*----------------------------------------------------------------------------*/
// SlalomModel node classes

// class Node_alpha {
// public:
//     double pa;
//     double pb;
//     arma::vec a;
//     arma::vec b;
//     arma::vec E1;
//     arma::vec lnE1;
//     // constructor
//     Node_alpha () {};
//     Node_alpha (arma::vec prior, int K) : pa(prior[0]), pb(prior[1]) {
//         a = arma::ones(K) * pa;
//         b = arma::ones(K) * pb;
//         E1 = b / a;
//         lnE1 = arma::ones(K);
//         for (int i = 0; i < K; i++) {
//             lnE1(i) = boost::math::digamma(a(i));
//         };
//         lnE1 = lnE1 - arma::log(b); // to be assigned to: digamma(a) - ln(b)
//         };
// };
//
//
// class Node_epsilon {
// public:
//     double pa;
//     double pb;
//     arma::vec a;
//     arma::vec b;
//     arma::vec E1;
//     arma::vec lnE1;
//     arma::vec diagSigmaS;
//     // constructor
//     Node_epsilon () {};
//     Node_epsilon (arma::vec prior, int G, int K) : pa(prior[0]), pb(prior[1]) {
//         a = arma::ones(G) * pa;
//         b = arma::ones(G) * pb;
//         E1 = b / a;
//         lnE1 = arma::ones(K);
//         for (int i = 0; i < K; i++) {
//             lnE1(i) = boost::math::digamma(a(i));
//         };
//         lnE1 = lnE1 - arma::log(b); // to be assigned to: digamma(a) - ln(b)
//         diagSigmaS = arma::zeros(K);
//         };
// };
//

// class Node_pi {
// public:
//     double pa;
//     double pb;
//     arma::vec a;
//     arma::vec b;
//     arma::vec E1;
//     arma::vec lnE1;
//     // constructor
//     Node_pi (arma::vec prior, int dim, int K) : pa(prior[0]), pb(prior[1])) {
//         a = arma::ones(dim) * pa;
//         b = arma::ones(dim) * pb;
//         E1 = b / a;
//         lnE1 = void; // to be assigned to: digamma(a) - ln(b)
//         };
// };

// class Node_X {
//     // factor states
// public:
//     arma::mat E1;
//     arma::mat diagSigmaS;
//     // constructor
//     Node_X () {};
//     Node_X (int N, int K) {
//         E1 = arma::randn(N, K);
//         diagSigmaS = arma::ones(N, K);
//     };
//     Node_X (arma::mat E1_in, int N, int K) :  E1(E1_in) {
//         diagSigmaS = arma::ones(N, K);
//     };
//
// };
//
// class Node_W {
//     // factor weights
// public:
//     arma::mat E1;
//     arma::mat sigma2;
//     arma::mat E2diag;
//     arma::mat gamma;
//     // constructor
//     Node_W () {};
//     Node_W (arma::mat gamma_in, int G, int K) : gamma(gamma_in) {
//         E1 = arma::randn(G, K);
//         sigma2 = arma::ones(G, K);
//         E2diag = arma::zeros(G, K);
//     };
//     Node_W (arma::mat gamma_in, arma::mat E1_in, int G, int K) :
//         gamma(gamma_in), E1(E1_in) {
//         sigma2 = arma::ones(G, K);
//         E2diag = arma::zeros(G, K);
//     };
// };

