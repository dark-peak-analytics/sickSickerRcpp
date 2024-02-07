#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::plugins(cpp11)]]

std::vector<double> custom_Pow( const std::vector<double> base,
                                const std::vector<double> expo )
{
  // This function is a vectorised version of std::pow
  // Each element on the base vector is raised to the corresponding element in the expo vector
  // This function requires cpp11 plugin
  std::vector<double> res(expo.size());
  std::transform(base.begin(), base.end(), expo.begin(), res.begin(),
                 [&](double lhs, double rhs) -> double {
                   return (1/ std::pow(lhs, rhs));
                 });
  return res;
}

// [[Rcpp::export]]
arma::mat ProbsV_Cpp( arma::rowvec v_S_t,
                      int& n_I,
                      int& n_S,
                      NumericVector& t_P )
{
  // v_S_t: numeric vector containing the health states occupied by the individuals at cycle t
  // n_I: number of simulated individuals
  // n_S: number of health states
  // t_P: vector containing transition probabilities.
  
  // declare other R functions:
  Function cat("cat") ;
  Function printR("print") ;
  Function roundR("round") ;
  
  // create a matrix for the state transition probabilities (m_P_t):
  arma::mat m_P_t(n_I, n_S) ;
  
  // define probabilities:
  arma::mat m_H  = arma::repmat( arma::rowvec {1 - t_P["p.HS1"] - t_P["p.HD"], t_P["p.HS1"], 0, t_P["p.HD"]}, n_I, 1 ) ;
  arma::mat m_S1 = arma::repmat( arma::rowvec {t_P["p.S1H"], 1 - t_P["p.S1H"] - t_P["p.S1S2"] - t_P["p.S1D"], t_P["p.S1S2"], t_P["p.S1D"]}, n_I, 1 ) ;
  arma::mat m_S2 = arma::repmat( arma::rowvec {0, 0, 1 - t_P["p.S2D"], t_P["p.S2D"]}, n_I, 1 ) ;
  arma::mat m_D  = arma::repmat( arma::rowvec {0, 0, 0, 1}, n_I, 1 ) ;
  
  // update v.p.it with the appropriate probabilities:
  m_P_t.rows ( arma::find( v_S_t == 1 ) ) = m_H.rows  ( arma::find( v_S_t == 1 ) ) ;
  m_P_t.rows ( arma::find( v_S_t == 2 ) ) = m_S1.rows ( arma::find( v_S_t == 2 ) ) ;
  m_P_t.rows ( arma::find( v_S_t == 3 ) ) = m_S2.rows ( arma::find( v_S_t == 3 ) ) ;
  m_P_t.rows ( arma::find( v_S_t == 4 ) ) = m_D.rows  ( arma::find( v_S_t == 4 ) ) ;
  
  // check if vector of probabilities sum to 1:
  // need to round up to the 1e-6, otherwise it breaks, e.g, once it was 1.00000001
  colvec sums = as< Rcpp::NumericVector >(Rcpp::wrap(roundR(arma::sum(m_P_t, 1), 6))) ;
  bool notSumToOne = any(sums != 1.0000) ;
  
  if(notSumToOne) {
    stop("Probabilities do not sum to 1!") ;
  }
  else
    return(m_P_t) ;
}

// [[Rcpp::export]]
arma::mat SampleV_Cpp( arma::mat m_P_t,
                       int& n_I,
                       int& n_S,
                       int m = 1)
{
  // m_P_t: numeric matrix containing the transition probabilities for individuals at cycle t
  // n_I: number of simulated individuals.
  // n_S: number of health states.
  // m: number of health states to sample.
  
  // declare some assistive functions from R:
  Function printR("print") ;
  Function roundR("round") ;
  
  // create a matrix for sampled health states (m_M_t):
  arma::mat m_M_t(n_I, m, arma::fill::ones) ;
  
  // create a matrix m_CumProb_t:
  arma::mat m_CumProb_t = arma::cumsum(m_P_t, 1) ;
  
  // recheck the probabilities:
  // need to round up to the 1e-6, otherwise it breaks, e.g, once it was 1.00000001
  colvec v_CumProb_t = as< Rcpp::NumericVector >(Rcpp::wrap(
    roundR(m_CumProb_t.col(m_CumProb_t.n_cols - 1), 6))) ;
  // the cumulative probability in the last column is expected to be equal to 1
  if(any(v_CumProb_t != 1.0000)) {
    stop("error in multinom: probabilities do not sum to 1") ;
  }
  
  // sample from a uniform U~(0, 1) distribution:
  arma::mat m_U(n_I, n_S, fill::ones) ; // a matrix to save values sampled from the U~(0, 1)
  for(int i = 0; i < m; i++) {
    // in each row, sample one random value for n_I individuals and repeat that value n_S times
    m_U = arma::repmat( randu<colvec>(n_I), 1, n_S ) ;
    // identify the first individual-specific health-state with a cumulative probability higher than the their corresponding sampled value
    // using a logical (true/false or 1/0), matrix get the value to be added to 1 (the starting health-state)
    // one plus the sum of each row of the logical matrix gives the health-state for the corresponding individuals at time t + 1
    m_M_t.col(i) = m_M_t.col(i) + arma::sum( (m_U > m_CumProb_t), 1 ) ;
  }
  
  return(m_M_t) ;
}

// [[Rcpp::export]]
arma::colvec CostsV_Cpp( arma::colvec v_S_t,
                         int& n_I,
                         int& n_S,
                         NumericVector& v_Costs,
                         bool b_Trt = false )
{
  // v_S_t: vector of health states occupied by individuals at cycle t
  // n_I: number of simulated individuals.
  // n_S: number of health states.
  // v_Costs: a vector containing cost parameters.
  // b_Trt: a bool indicating whether treatment costs (default is false).
  
  // declaring costs for the individuals at time t and set it to 0:
  arma::mat m_C_t (n_I, 1, fill::zeros) ;
  
  // define state costs:
  arma::mat m_H  = arma::repmat( arma::rowvec {(v_Costs["c.H"])}, n_I, 1 ) ;
  arma::mat m_S1 = arma::repmat( arma::rowvec {(v_Costs["c.S1"] + (v_Costs["c.Trt"] * b_Trt))}, n_I, 1 ) ;
  arma::mat m_S2 = arma::repmat( arma::rowvec {(v_Costs["c.S2"] + (v_Costs["c.Trt"] * b_Trt))}, n_I, 1 ) ;
  
  // update m_C_t with the appropriate costs:
  m_C_t.rows ( arma::find( v_S_t == 1 ) ) = m_H.rows  ( arma::find( v_S_t == 1 ) ) ;
  m_C_t.rows ( arma::find( v_S_t == 2 ) ) = m_S1.rows ( arma::find( v_S_t == 2 ) ) ;
  m_C_t.rows ( arma::find( v_S_t == 3 ) ) = m_S2.rows ( arma::find( v_S_t == 3 ) ) ;
  
  return(m_C_t.col(0)) ;
}

// [[Rcpp::export]]
arma::colvec EffsV_Cpp( arma::colvec v_S_t,
                        int& n_I,
                        int& n_S,
                        NumericVector& v_Utilities,
                        bool b_Trt = false,
                        int cl = 1 )
{
  // v_S_t: vector of health states occupied by individuals at cycle t
  // n_I: number of simulated individuals.
  // n_S: number of health states.
  // v_Utilities: a vector containing utilities values for each health state
  // b_Trt: bool indicating whether treatment costs (default is false).
  // cl: integer variable indicating cycle length (cl) - default is 1
  
  // declaring effects matrix for the individuals at time t and set it to 0:
  arma::mat m_E_t (n_I, 1, fill::zeros) ;
  
  // define state QALYs:
  arma::mat m_H  = arma::repmat( arma::rowvec {(v_Utilities["u.H"] * cl)}, n_I, 1 ) ;
  arma::mat m_S1 = arma::repmat( arma::rowvec {((b_Trt * v_Utilities["u.Trt"] + (1 - b_Trt) * v_Utilities["u.S1"]) * cl)}, n_I, 1 ) ;
  arma::mat m_S2 = arma::repmat( arma::rowvec {(v_Utilities["u.S2"] * cl)}, n_I, 1 ) ;
  
  // update m_E_t with the appropriate health effects:
  m_E_t.rows ( arma::find( v_S_t == 1 ) ) = m_H.rows  ( arma::find( v_S_t == 1 ) ) ;
  m_E_t.rows ( arma::find( v_S_t == 2 ) ) = m_S1.rows ( arma::find( v_S_t == 2 ) ) ;
  m_E_t.rows ( arma::find( v_S_t == 3 ) ) = m_S2.rows ( arma::find( v_S_t == 3 ) ) ;
  
  return(m_E_t.col(0)) ;
}

/************* Micro-simulation ***************/

// [[Rcpp::export]]
List MicroSimV_Cpp( arma::colvec& v_S_t,
                    NumericVector t_P,
                    NumericVector v_C,
                    NumericVector v_U,
                    int n_I,
                    int n_S = 4,
                    int n_T = 30,
                    int n_Cl = 1,
                    double d_dC = 0.03,
                    double d_dE = 0.03,
                    bool b_Trt = false,
                    int n_Seed = 1 ) {
  // Arguments:
  // v_S_t:  numeric vector containing the health states occupied by the individuals at cycle t
  // t_P:    vector containing transition probabilities
  // v_C:    vector containing costs
  // v_U:    vector containing utilities
  // n_I:    number of individuals
  // n_S:    number of health-states
  // n_T:    total number of cycles to run the model
  // n_Cl:   length of model cycle (default is 1 year)
  // d_dC:   discount rate for costs
  // d_dE:   discount rate for health outcome (QALYs)
  // b_Trt:  are the n.i individuals receiving treatment? (scalar with a Boolean value, default is FALSE)
  // n_Seed: starting seed number for random number generator (default is 1)
  
  // Makes use of:
  // ProbsV_Cpp:  function for the estimation of transition probabilities
  // SampleV_Cpp: function for sampling the health states that will be occupied by individuals at (t+1)
  // CostsV_Cpp:  function for the estimation of costs associated with state at (t+1)
  // EffsV_Cpp:   function for the estimation of QALYs associated with state at (t+1)
  
  // declare R's set.seed function to control setting the seed number:
  Function set_seed("set.seed") ;
  
  // declare other R functions:
  Function cat("cat") ;
  Function paste("paste") ;
  
  // declare discounting vectors and estimate them:
  arma::vec v_dwc = custom_Pow(
    as< std::vector<double> >(Rcpp::NumericVector {wrap(Rcpp::rep((1 + d_dC), (n_T + 1)))} ),
    as< std::vector<double> >(Rcpp::NumericVector {wrap(Rcpp::seq(0, n_T))} )) ;
  
  arma::vec v_dwe = custom_Pow(
    as< std::vector<double> >(Rcpp::NumericVector {wrap(Rcpp::rep((1 + d_dE), (n_T + 1)))} ),
    as< std::vector<double> >(Rcpp::NumericVector {wrap(Rcpp::seq(0, n_T))} )) ;
  
  // declare transitions, costs, and effects matrices:
  arma::mat m_M(n_I, (n_T + 1), fill::zeros) ;
  arma::mat m_C(n_I, (n_T + 1), fill::zeros) ;
  arma::mat m_E(n_I, (n_T + 1), fill::zeros) ;
  
  // indicate the initial health state:
  m_M.col(0) = v_S_t ;
  
  // set the seed for every individual for the random number generator
  // this is done in R as WM still does not know how to do the same in R
  set_seed(Named("seed", (n_Seed))) ;
  
  // estimate costs for all individuals at the initial health state conditional on treatment:
  m_C.col(0) = CostsV_Cpp(m_M.col(0), n_I, n_S, v_C, b_Trt) ;
  
  // estimate QALYs for all individuals at the initial health state conditional on treatment:
  m_E.col(0) = EffsV_Cpp (m_M.col(0), n_I, n_S, v_U, b_Trt, n_Cl) ;
  
  for(int t = 0; t < n_T; t++) {
    
    // calculate the transition probabilities for cycle (t + 1):
    arma::mat m_P  = ProbsV_Cpp (arma::trans(m_M.col(t)), n_I, n_S, t_P) ;
    // sample the health state at (t+1):
    m_M.col(t + 1) = SampleV_Cpp(m_P, n_I, n_S, 1) ;
    // estimate costs for all individuals occupying health states at (t+1) conditional on treatment:
    m_C.col(t + 1) = CostsV_Cpp (m_M.col(t + 1), n_I, n_S, v_C, b_Trt) ;
    // estimate QALYs for all individuals occupying health states at (t+1) conditional on treatment:
    m_E.col(t + 1) = EffsV_Cpp  (m_M.col(t + 1), n_I, n_S, v_U, b_Trt, n_Cl) ;
    
    // print the progress of the simulation:
    cat(Named("...", '\r'),
        Named("...", paste(Named("...", (round(((double)(t + 1)/n_T) * 100 * 10) / 10)),
                        Named("...", "% done "), Named("sep", " ")))) ;
  } // end of loops through cycles:
  
  // compute summary stats:
  // total discounted costs and QALYs (per individual):
  arma::vec tc = m_C * v_dwc ;
  arma::vec te = m_E * v_dwe ;
  
  // average discounted costs and average discounted QALYs:
  double tc_hat = arma::mean(tc) ;
  double te_hat = arma::mean(te) ;
  
  // save the results in a list:
  List results = List::create(Named("m.M") = m_M, _["m.C"] = m_C, _["m.E"] = m_E,
                              _["tc"] = tc, _["te"] = te, _["tc_hat"] = tc_hat,
                                _["te_hat"] = te_hat) ;
  
  return(results) ;
}

// [[Rcpp::export]]
List SickSickerMicroSim_Cpp( arma::colvec& v_S_t,
                             NumericVector t_P,
                             NumericVector v_C,
                             NumericVector v_U,
                             int n_I,
                             int n_S = 4,
                             int n_T = 30,
                             int n_Cl = 1,
                             double d_dC = 0.03,
                             double d_dE = 0.03,
                             int n_Seed = 1,
                             bool b_Trt = false ) {
  // Arguments:
  // v_S_t:  numeric vector containing the health states occupied by the individuals at cycle t
  // t_P:    vector containing transition probabilities
  // v_C:    vector containing costs
  // v_U:    vector containing utilities
  // n_I:    number of individuals
  // n_S:    number of health-states
  // n_T:    total number of cycles to run the model
  // n_Cl:   length of model cycle (default is 1 year)
  // d_dC:   discount rate for costs
  // d_dE:   discount rate for health outcome (QALYs)
  // b_Trt:  are the n.i individuals receiving treatment? (scalar with a Boolean value, default is FALSE)
  // n_Seed: starting seed number for random number generator (default is 1)
  
  // Makes use of:
  // ProbsV_Cpp:  function for the estimation of transition probabilities
  // SampleV_Cpp: function for sampling the health states that will be occupied by individuals at (t+1)
  // CostsV_Cpp:  function for the estimation of costs associated with state at (t+1)
  // EffsV_Cpp:   function for the estimation of QALYs associated with state at (t+1)
  
  // declare R's set.seed function to control setting the seed number:
  Function set_seed("set.seed") ;
  
  // declare othre R functions:
  Function cat("cat") ;
  Function paste("paste") ;
  
  // declare discounting vectors and estimate them:
  arma::vec v_dwc = custom_Pow(
    as< std::vector<double> >(Rcpp::NumericVector {wrap(Rcpp::rep((1 + d_dC), (n_T + 1)))} ),
    as< std::vector<double> >(Rcpp::NumericVector {wrap(Rcpp::seq(0, n_T))} )) ;
  
  arma::vec v_dwe = custom_Pow(
    as< std::vector<double> >(Rcpp::NumericVector {wrap(Rcpp::rep((1 + d_dE), (n_T + 1)))} ),
    as< std::vector<double> >(Rcpp::NumericVector {wrap(Rcpp::seq(0, n_T))} )) ;
  
  // declare transitions, costs, and effects matrices:
  arma::mat m_M(n_I, (n_T + 1), fill::zeros) ;
  arma::mat m_C(n_I, (n_T + 1), fill::zeros) ;
  arma::mat m_E(n_I, (n_T + 1), fill::zeros) ;
  
  // indicate the initial health state:
  m_M.col(0) = v_S_t ;
  
  // set the seed for every individual for the random number generator
  // this is done in R as WM still does not know how to do the same in R
  set_seed(Named("seed", (n_Seed))) ;
  
  // estimate costs for all individuals at the initial health state conditional on treatment:
  m_C.col(0) = CostsV_Cpp(m_M.col(0), n_I, n_S, v_C, b_Trt) ;
  
  // estimate QALYs for all individuals at the initial health state conditional on treatment:
  m_E.col(0) = EffsV_Cpp (m_M.col(0), n_I, n_S, v_U, b_Trt, n_Cl) ;
  
  for(int t = 0; t < n_T; t++) {
    
    // calculate the transition probabilities for cycle (t + 1):
    arma::mat m_P  = ProbsV_Cpp (arma::trans(m_M.col(t)), n_I, n_S, t_P) ;
    // sample the health state at (t+1):
    m_M.col(t + 1) = SampleV_Cpp(m_P, n_I, n_S, 1) ;
    // estimate costs for all individuals occupying health states at (t+1) conditional on treatment:
    m_C.col(t + 1) = CostsV_Cpp (m_M.col(t + 1), n_I, n_S, v_C, b_Trt) ;
    // estimate QALYs for all individuals occupying health states at (t+1) conditional on treatment:
    m_E.col(t + 1) = EffsV_Cpp  (m_M.col(t + 1), n_I, n_S, v_U, b_Trt, n_Cl) ;
    
    // print the progress of the simulation:
    cat(Named("...", '\r'),
        Named("...", paste(Named("...", (round(((double)(t + 1)/n_T) * 100 * 10) / 10)),
                        Named("...", "% done "), Named("sep", " ")))) ;
  } // end of loops through cycles:
  
  // compute summary stats:
  // total discounted costs and QALYs (per individual):
  arma::vec tc = m_C * v_dwc ;
  arma::vec te = m_E * v_dwe ;
  
  // average discounted costs and average discounted QALYs:
  double tc_hat = arma::mean(tc) ;
  double te_hat = arma::mean(te) ;
  
  // save the results in a list:
  List results = List::create(Named("m.M") = m_M,
                              _["tc_hat"] = tc_hat,
                              _["te_hat"] = te_hat) ;
  
  return(results) ;
}
