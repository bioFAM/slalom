// #include <RcppArmadillo.h>
// //[[Rcpp::depends(RcppArmadillo)]]
// using namespace Rcpp;
// using namespace std;
//
// // This is a simple example of exporting a C++ function to R. You can
// // source this function into an R session using the Rcpp::sourceCpp
// // function (or via the Source button on the editor toolbar). Learn
// // more about Rcpp at:
// //
// //   http://www.rcpp.org/
// //   http://adv-r.had.co.nz/Rcpp.html
// //   http://gallery.rcpp.org/
// //
//
// //' @title
// //' Initialise a Slalom model
// //' @description
// //' Initialise a Slalom model by providing the necessary data objects.
// //'
// //' @param Y_init matrix of expression values
// //' @param pi_init G x K matrix with each entry being the prior
// //' probability for a gene g being active for factor k.
// //' @param X_init matrix of initial factor states (N x K)
// //' @param W_init G x K matrix of initial weights
// //' @param prior_alpha numeric vector of length two giving prior values for the gamma hyperparameters of the precisions
// //' @param prior_epsilon numeric vector of length two giving prior values for the gamma hyperparameters of the residual variances
// //'
// //' @return
// //' an object of the Slalom class
// //'
// //' @useDynLib slalom
// //' @importFrom Rcpp evalCpp
// //' @export
// // [[Rcpp::export]]
// Slalom newSlalom(
//     arma::mat Y_init,
//     arma::mat pi_init,
//     arma::mat X_init,
//     arma::mat W_init,
//     arma::vec prior_alpha,
//     arma::vec prior_epsilon
//     ) {
//     Slalom sla = Slalom(Y_init, pi_init, X_init, W_init, prior_alpha,
//                         prior_epsilon);
//     return sla;
// }
//
//
// class Slalom {
// // base class for slalom holding models
// public:
//     Node_alpha alpha;
//     Node_epsilon epsilon;
//     //Node_pi pi;
//     arma::mat pi;
//     Node_X X;
//     Node_W W;
//     arma::mat Y;
//     arma::mat pseudo_Y;
//     arma::vec YY;
//     int K;
//     int N;
//     int G;
//     double nScale;
//     arma::vec iUnannotatedDense;
//     arma::vec iUnannotatedSparse;
//     // constructor
//     Slalom (arma::mat Y_init, arma::mat pi_init, arma::mat X_init, arma::mat W_init, arma::vec prior_alpha, arma::vec prior_epsilon) :
//         Y(Y_init) {
//         K = pi_init.n_cols();
//         N = Y_init.n_rows();
//         G = Y_init.n_cols();
//         nScale = 100;
//         alpha = Node_alpha(prior_alpha, K);
//         epsilon = Node_epsilon(prior_epsilon, G, K);
//         //pi = Node_pi();
//         pi = pi_init;
//         X = Node_X(X_init, N, K);
//         W = Node_W(pi_init, W_init, G, K);
//         YY = arma::zeros(G);
//         for (int g=0; g < G; g++) {
//             YY(g) = arma::accu(Y.col(g) % Y.col(g))
//         };
//     };
// };
//
//
// // Slalom class initialiation method //
// void Slalom::init() {
//
// };
//
//
// // Slalom class update methods //
//
// void Slalom::updateAlpha(m) {
//     //update alpha in a Slalom class object - precisions
//     Ewdwd = arma::accu(this->W.gamma.col(m) * this->W.E2diag.col(m));
//     this->alpha.a(m) = this->alpha.pa + 0.5 * Ewdwd;
//     this->alpha.b(m) = this->alpha.pb + arma::accu(this->W.gamma.col(m)) / 2.0;
//     this->alpha.E1(m) = this->alpha.b(m) / this->alpha.a(m);
// };
//
//
// void Slalom::updateEpsilon(m) {
//     //update Epsilon (vectorised) - noise parameters
//     SW_sigma = this->W.gamma * this->W.E1;
//     SW2_sigma = this->W.gamma * this->W.E2diag;
//
//     muSTmuS = this->X.E1.t() * this->X.E1;
//     muSTmuS.diag().zeros();
//
//     t1 = arma::sum(SW_sigma % (this->Y.t() * this->X.E1), 1);
//     t2 = arma::sum(SW2_sigma % arma::repmat(muSTmuS.diag().t() + this->epsilon.diag_sigma_S, this->G, 1), 1);
//     t3 = arma::sum((SW_sigma  * muSTmuS0) % SW_sigma, 1);
//
//     this->epsilon.E1 = 1.0 / ((this->YY + (-2 * t1 + t2 + t3)) / this->N);
//     this->epsilon.E1(this->epsilon.E1 > 1e6) = 1e6;
// };
//
//
// void Slalom::updateX(m) {
//     // update the factor states
//     arma::vec set1 = arma::regspace(0, 1, m - 1);
//     arma::vec set2 = arma::regspace(m + 1, 1, this->K - 1);
//     arma::vec setMinus = arma::join_cols(set1, set2);
//
//     // update X
//     arma::mat SW_sigma;
//     arma::mat SW2_sigma;
//     SW_sigma = (this->W.gamma.col(m) % this->W.E1.col(m)) % this->epsilon.E1;
//     SW2_sigma = (this->W.gamma.col(m) % this->W.E2diag.col(m)) % this->epsilon.E1;
//
//     b0 = this->X.E1.col(setMinus) * (this->W.gamma.col(setMinus) % this->.W.E1[:, setMinus]).t();
//     b = b0 * SW_sigma;
//
//     alphaSm = arma::sum(SW2_sigma, 0);
//     barmuS = (this->Y * SW_sigma) - b;
//     this->X.diag_sigma_S.col(m) = 1.0 / (1.0 + alphaSm);
//     this->X.E1.col(m) = barmuS / (1.0 + alphaSm);
//     // keep diag_sigma_S
//     this->epsilon.diag_sigma_S(m) = arma::accu(this->X.diag_sigma_S.col(m));
// };
//
//
// void Slalom::updateW(m) {
//     // update the factor weights
//     arma::vec logPi;
//     if (arma::any(this->iUnannotatedSparse == m) | arma::any(this->iUnannotatedDense == m)) {
//         logPi = arma::log(this->pi.col(m) / (1 - this->pi.col(m));
//     } else if(this->nScale > 0 & this->nScale < this->N) {
//         logPi = arma::log(this->pi.col(m) / (1 - this->pi.col(m));
//         isOFF_ = this->pi.col(m) < .5;
//         logPi(isOFF_) = (this->N / self.nScale) % arma::log(this->pi.col(m).row(isOFF_) / (1 - this->pi.col(m).row(isOFF_)));
//     } else {
//         logPi = arma::log(this->pi.col(m) / (1 - this->pi.col(m));
//     }
//
//     arma::vec sigma2Sigmaw;
//     sigma2Sigmaw = (1.0 / this->epsilon.E1) % this->alpha.E1(m);
//
//     arma::vec set1 = arma::regspace(0, 1, m - 1);
//     arma::vec set2 = arma::regspace(m + 1, 1, this->K - 1);
//     arma::vec setMinus = arma::join_cols(set1, set2);
//
//     SmTSk = arma::sum( arma::repmat(this->X.E1.col(m), 1, this->K - 1) % this->X.E1.col(setMinus), 0);
//     SmTSm = (this->X.E1.col(m).t() * this->X.E1.col(m)) + arma::accu(this->X.diag_sigma_S.cols(m));
//
//     arma::vec b;
//     arma::vec diff;
//     b = (this->W.gamma.col(setMinus) % this->W.E1.col(setMinus)) * SmTSk.t());
//     diff = (this->X.E1.col(m).t() * this->Y) - b;
//     arma::vec diff2 = diff % diff;
//
//     SmTSmSig = SmTSm + sigma2Sigmaw;
//
//     //update C and W
//     arma::vec
//     u_qm = logPi + 0.5 % arma::log(sigma2Sigmaw) - 0.5 % arma::log(SmTSmSig) + (0.5 % this->epsilon.E1) % ( diff2 / SmTSmSig);
//     this->W.gamma.col(m) = 1.0 / (1 + arma::exp(-u_qm));
//     // what do do about this line?? first time we see self.W.C[:,m,1]
//     self.W.C[:,m,1] = 1-self.W.C[:,m,0];
//     // check line above
//     this->.W.E1.col(m) = (diff / SmTSmSig);
//     this->W.sigma2.col(m) = (1.0 / this->epsilon.E1) / SmTSmSig;
//     // check how to element-wise square a vector with Armadillo
//     this->W.E2diag.col(m) = this->W.E1.col(m) % this->W.E1.col(m) + this->W.sigma2.col(m);
// };
//
//
// class Node_alpha {
// public:
//     double pa;
//     double pb;
//     arma::vec a;
//     arma::vec b;
//     arma::vec E1;
//     arma::vec lnE1;
//     // constructor
//     Node_alpha (arma::vec prior, int K) : pa(prior[0]), pb(prior[1])) {
//         a = arma::ones(K) * pa;
//         b = arma::ones(K) * pb;
//         E1 = b / a;
//         lnE1 = void; // to be assigned to: digamma(a) - ln(b)
//     };
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
//     arma::vec diag_sigma_S;
//     // constructor
//     Node_epsilon (arma::vec prior, int G, int K) : pa(prior[0]), pb(prior[1])) {
//         a = arma::ones(G) * pa;
//         b = arma::ones(G) * pb;
//         E1 = b / a;
//         lnE1 = void; // to be assigned to: digamma(a) - ln(b)
//         diag_sigma_S = arma::zeros(K);
//         };
// };
//
//
// // class Node_pi {
// // public:
// //     double pa;
// //     double pb;
// //     arma::vec a;
// //     arma::vec b;
// //     arma::vec E1;
// //     arma::vec lnE1;
// //     // constructor
// //     Node_pi (arma::vec prior, int dim, int K) : pa(prior[0]), pb(prior[1])) {
// //         a = arma::ones(dim) * pa;
// //         b = arma::ones(dim) * pb;
// //         E1 = b / a;
// //         lnE1 = void; // to be assigned to: digamma(a) - ln(b)
// //         };
// // };
//
// class Node_X {
// // factor states
// public:
//     arma::mat E1;
//     arma::mat diag_sigma_S;
//     // constructor
//     Node_X (int N, int K) :  {
//         E1 = arma::randn(N, K);
//         diag_sigma_S = arma::ones(N, K);
//     };
//     Node_X (arma::mat E1_in, int N, int K) :  E1(E1_in) {
//         diag_sigma_S = arma::ones(N, K);
//     };
//
// };
//
// class Node_W {
// // factor weights
// public:
//     arma::mat E1;
//     arma::mat sigma2;
//     arma::mat E2diag;
//     arma::mat gamma;
//     // constructor
//     Node_W (arma::mat gamma_in, int G, int K) : gamma(gamma_in) {
//         E1 = arma::randn(G, K);
//         sigma2 = arma::ones(G, K);
//         E2diag = arma::zeros(G, K);
//     };
//     Node_W (arma::mat gamma_in, arma::mat E1_in, int G, int K, ) :
//         gamma(gamma_in), E1(E1_in) {
//         sigma2 = arma::ones(G, K);
//         E2diag = arma::zeros(G, K);
//     };
// };
