#include <Rcpp.h>

// [[Rcpp::export]]
Rcpp::List simulate_birth_death_cpp(double lambda, double mu, int n, double max_time) {
  double current_time = 0;
  double previous_time = 0;
  int Z = n;
  std::vector<double> times(1, 0);
  std::vector<int> population_sizes(1, n);
  int extinction = 0;
  int B_t = 0;
  int D_t = 0;
  double S_t = 0;
  
  double birth_prob = lambda / (lambda + mu);
  
  while (current_time < max_time) {
    double rate = Z * (lambda + mu);
    Rcpp::NumericVector uniformSample = Rcpp::runif(2);
    double wait_time = -log(1-uniformSample[0])*1/rate; // exponential distribution with rate Z * (lambda + mu).
    previous_time = current_time;
    current_time += wait_time;
    
    S_t += Z * (current_time - previous_time);

    
    bool is_birth = uniformSample[1] < birth_prob;
    
    if (is_birth) {
      Z ++;
      B_t ++;
    } else {
      Z--;
      D_t++;
      if (Z == 0) {
        extinction = 1;
        break;
      }
    }
    
    times.push_back(current_time);
    population_sizes.push_back(Z);
  }
  
  return Rcpp::List::create(Rcpp::Named("times") = times,
                      Rcpp::Named("Z") = population_sizes,
                      Rcpp::Named("extinction") = extinction,
                      Rcpp::Named("S_t") = S_t,
                      Rcpp::Named("B_t") = B_t,
                      Rcpp::Named("D_t") = D_t);
}