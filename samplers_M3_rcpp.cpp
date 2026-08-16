// =============================================================================
// sample_xi_E_cpp.cpp
//
// C++ / RcppArmadillo translation of sample_xi_E() (R function).
// E-level cluster-assignment updating for a three-DP Dirichlet Process model.
//
//   Outputs:
//   list with:
//     - xi_E    : updated E-level cluster assignments (relabelled 1,...,K)
//     - beta_k_E: updated list of E-level cluster coefficients
// =============================================================================

// RcppArmadillo for interface between R and cpp, and lineal algebra (vectors and matrices)
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <map>       // Dictionary (key -> value)
#include <set>       // unique sorted elements
#include <vector>    // Dynamic array 
#include <cmath>     // math functions (log, sqrt, etc)

// [[Rcpp::export]]
Rcpp::List sample_xi_E_cpp(         // Output is a list
    const arma::vec&  y,            // (n)        response variable
    const arma::mat&  X_E_full,     // (n x p_E)  individual-level covariate matrix
    arma::ivec        xi_E,         // (n)        current E-level cluster assignments (1-based)
    Rcpp::List        beta_k_E,     // (K_E)      E-level cluster coefficients (consecutive 1-based)
    double            sigma2_E_base,//            base measure variance sigma^2_E
    double            alpha_E,      //            E-level DP concentration parameter
    const arma::vec&  mu_E,         // (p_E)      base measure mean
    const arma::mat&  Z_full,       // (n x p_M)  municipality-level covariate matrix
    const arma::mat&  W_full,       // (n x p_D)  department-level covariate matrix
    const arma::ivec& xi_M,         // (n)        current M-level cluster assignments (1-based)
    const arma::ivec& xi_D,         // (n)        current D-level cluster assignments (1-based)
    Rcpp::List        beta_k_M,     // (K_M)      M-level cluster coefficients (consecutive 1-based)
    Rcpp::List        beta_k_D,     // (K_D)      D-level cluster coefficients (consecutive 1-based)
    const arma::vec&  mun_kappa,    // (M_mun)    kappa^2_{j,q} per municipality, 0-based in C++
    const arma::ivec& obs_to_mun,   // (n)        municipality index per obs, 1-based (R convention)
    const arma::vec&  xTx_vec,      // (n)        pre-computed x_i^T x_i per observation
    const arma::vec&  sqrt_kappa,   // (M_mun)    pre-computed sqrt(kappa^2_{j,q}), 0-based in C++
    double            beta_int      //            global intercept beta
) {
  // -------------------------------------------------------------------------
  // 0. Dimensions and global constants
  // -------------------------------------------------------------------------
  const int    n           = static_cast<int>(y.n_elem);         // .n_elem: number of elements in Armadillo vector y
  const int    p_E         = static_cast<int>(X_E_full.n_cols);  // n_cols num of cols. static_cast<int> converts the size type to an integer.
  
  // Pre-compute constants used in new cluster beta sampling
  const double sigma2_inv  = 1.0 / sigma2_E_base;
 
  // sigma2_inv * I_{p_E} -- used every time a new cluster is born
  const arma::mat sigma2_diag = sigma2_inv * arma::eye<arma::mat>(p_E, p_E);

  // -------------------------------------------------------------------------
  // 1. Load beta_k_E into a map indexed by 1-based cluster ID
  //    Input is consecutive [1..K_E_init]; but during the sweep IDs can gap
  //     as clusters may disappear.
  // -------------------------------------------------------------------------
  // Some clusters may disappear during the algorithm, maping let us store coefficients indexed by their
  //explicit cluster IDs. Introduces a mapping for efficient cluster ID handling.
  std::map<int, arma::vec> beta_map_E;  //Declares a C++ map
  {
    const int K_E_init = static_cast<int>(beta_k_E.size()); // number (integer) of initial cluster coefficient vectors (initial k_E)
    for (int k = 1; k <= K_E_init; k++) {                   // From k to initial K_E. Initialised in 1 because R code initialise 1-based
      beta_map_E[k] = Rcpp::as<arma::vec>(beta_k_E[k - 1]); // 0-based R list -> 1-based map. Rcpp::as<arma::vec> Convert the list into Armadillo vector
    }                                                       // beta_k_E[k_E - 1] transform to 0-based C++ indexing
  }
  
  // -------------------------------------------------------------------------
  // 2. Load beta_k_M and beta_k_D into maps for fast contrib lookup
  //    These are fixed during the E-level sweep (not updated here)
  // -------------------------------------------------------------------------
  // same structure as for E level 
  std::map<int, arma::vec> beta_map_M;
  {
    const int K_M = static_cast<int>(beta_k_M.size());       // Const for non updating
    for (int k = 1; k <= K_M; k++) {
      beta_map_M[k] = Rcpp::as<arma::vec>(beta_k_M[k - 1]);
    }
  }
  std::map<int, arma::vec> beta_map_D;
  {
    const int K_D = static_cast<int>(beta_k_D.size());
    for (int k = 1; k <= K_D; k++) {
      beta_map_D[k] = Rcpp::as<arma::vec>(beta_k_D[k - 1]);
    }
  }
  // -------------------------------------------------------------------------
  // 3. Working copy of xi_E as std::vector<int> (1-based throughout)
  // -------------------------------------------------------------------------
  std::vector<int> xi_vec(n);                         // declares a vector of length n (integer array) std::vector
  for (int i = 0; i < n; i++) xi_vec[i] = xi_E[i];    // Copies the current labels from xi into xi_vec
 
  // -------------------------------------------------------------------------
  // 4. Running cluster-count map over all n observations
  //    Decrement/re-add observation i each iteration -> reflects xi_E_{-i}
  //    Replaces O(n) table(xi_E[-i]) calls with O(1) updates
  // -------------------------------------------------------------------------
  // Build a map from cluster ID to number of observation currently assigned to that cluster
  std::map<int, int> cnt;                         // stores ordered key–value pairs. cluster_id (1-based) ; cnt = count
  for (int i = 0; i < n; i++) cnt[xi_vec[i]]++;   // increments the count for the cluster label of observation i.

  // =========================================================================
  // 5. Main sweep (sequential, i = 0 .. n-1, xi_E updated in place)
  // =========================================================================
  for (int i = 0; i < n; i++) {         // Main loop Iterate over i = 1,...,n
    
    // --- 5a. Pre-computed values for observation i ---
    const int    mun0        = obs_to_mun[i] - 1;    // 0-based municipality index (-1 to 1-based to 0-based)
    const double kappa2_i    = mun_kappa[mun0];
    const double sqrt_kap_i  = sqrt_kappa[mun0];
    const double kappa2_inv  = 1.0 / kappa2_i;
    const arma::vec x_i      = X_E_full.row(i).t();  // Extracts covariates i row and transforms to column vector (p_E x 1) (Armadillo is column oriented)
    const double xTx_i       = xTx_vec[i];           // Precomputed xTx for i
    
    // --- 5b. Compute M-level and D-level contributions for observation i ---
    // These are fixed while we sample xi_E[i]
    // contrib_M = z_{j,q}^T beta_k_M  (municipality contribution)
    const arma::vec z_i      = Z_full.row(i).t();                       // (p_M x 1)
    const int    km_i        = xi_M[i];                                 // M-level cluster of obs i (1-based)
    const double contrib_M   = arma::dot(z_i, beta_map_M.at(km_i));     // Product z_i^T beta_k_M
    
    // contrib_D = w_q^T beta_k_D  (department contribution)
    const arma::vec w_i      = W_full.row(i).t();    // (p_D x 1)
    const int    kd_i        = xi_D[i];               // D-level cluster of obs i (1-based)
    const double contrib_D   = arma::dot(w_i, beta_map_D.at(kd_i));     // Product w_i^T beta_k_D
    
    // --- 5c. Residual r_i^(E) = y_i - (beta + contrib_M + contrib_D) ---
    // Used for: (1) new cluster marginal likelihood
    //           (2) new beta_k_E posterior sampling
    const double r_i_E = y[i] - (beta_int + contrib_M + contrib_D);      // Residual
    
    // --- 5d. Temporarily remove obs i from count map -> cnt = xi_E_{-i} ---
    const int old_cluster  = xi_vec[i];                    // Retrieves current cluster of observation i 
    cnt[old_cluster]--;                                    // Decrements the count in the mapping. -- subtracts 1
    const bool old_empty   = (cnt[old_cluster] == 0);      // boolean that checks if the count is 0
    if (old_empty) cnt.erase(old_cluster);                 // if count is zero (empty cluster), ereases it from cnt
    
    // --- 5e. Unique clusters and counts from xi_E_{-i} ---
    const int K_minus = static_cast<int>(cnt.size());      // Counts the number of unique clusters
    std::vector<int> uniq_cl;  uniq_cl.reserve(K_minus);   // Creates an array uniq_cl with k_minus elements
    std::vector<int> uniq_cnt; uniq_cnt.reserve(K_minus);  // Creates an array uniq_cnt with k_minus elements
    for (const auto& kv : cnt) {                           // iterates over every entry in cnt map
      uniq_cl.push_back(kv.first);   // cluster ID         // push_back adds elements to the end of the vector: i cluster ID 
      uniq_cnt.push_back(kv.second); // count n_{-i, k}    // push_back adds elements to the end of the vector: i cluster count
    }
    
    // --- 5f. Log-probability vector: slots 0..K_minus-1 = existing, K_minus = new ---
    //Allocate log probability vector
    // Creates a vector with one slot per existing cluster plus one for new cluster
    arma::vec log_probs(K_minus + 1);
    
    // - Probabilities for EXISTING clusters -
    // p(xi_i^E = k | -) propto n_{-i,k} * N(y_i | beta + x^T beta_k_E + contrib_M + contrib_D, kappa2)
    for (int k = 0; k < K_minus; k++) {                         // Iterates over clusters
      //Get beta for this cluster    
      const arma::vec& bk  = beta_map_E.at(uniq_cl[k]);         // Looks-up the coefficient vector for cluster uniq_cl[k] in beta_k(mapped)
      // Mean: beta + x_i^T beta_k_E + contrib_M + contrib_D
      const double mean_k  = beta_int + arma::dot(x_i, bk) + contrib_M + contrib_D; // performs dot product (covariates and betas by leve)
      // Log probability: log(n_{-i,k}) + log N(y_i | mean_k, kappa2_i)
      log_probs[k] = std::log(static_cast<double>(uniq_cnt[k])) // Prior probability 
        + R::dnorm(y[i], mean_k, sqrt_kap_i, /*log=*/1);        // likelihood in log
    }
    
    // - Probability for NEW cluster -
    // p(xi_i^E = K+1 | -) propto alpha_E * N(r_i^E | x_i^T mu_E, kappa2 + sigma2_E * x_i^T x_i)
    const double var_new  = kappa2_i + (sigma2_E_base * xTx_i); // kappa2 + sigma2_E * x_i^T x_i
    const double mean_new = arma::dot(x_i, mu_E);               // R: sum(x_i * mu_E)  x_i^T mu_E (evaluated on r_i^E)
    
    // compute the log of alpha_E * N(r_i^E | x_i^T mu_E, kappa2 + sigma2_E * x_i^T x_i)
    log_probs[K_minus] = std::log(alpha_E)
      + R::dnorm(r_i_E, mean_new, std::sqrt(var_new), /*log=*/1);
  
    // --- 5g. Sample new cluster assignment ---
    //# Normalize probabilities (log-sum-exp trick for numerical stability)
    const double max_lp = log_probs.max();                  // Finds the max value in log_probs
    arma::vec probs     = arma::exp(log_probs - max_lp);    // subtracts max value for stability
    probs              /= arma::sum(probs);                 // arma::sum(probs) and /= divides by the total sum
  
    // Inverse-CDF sampling using R's RNG (equivalent to R's sample.int)
    //        Index 0..K_minus-1 = existing clusters
    //        Index K_minus      = new cluster
    //        Equivalent to R: sample.int(K_new+1, ...) shifted by -1.
    const double u = R::unif_rand();           // random uniform betwen 0 and 1
    int sampled    = K_minus;                  // default: new cluster
    {  // manually implement inverse-CDF sampling from a uniform draw.
      double cum = 0.0;                        // Initialises a cummulative sum
      for (int k = 0; k <= K_minus; k++) {     // Iterates over k (including k_new)
        cum += probs[k];                       // Accumulates probability
        if (u <= cum) { sampled = k; break; }  // Once where the number landed is found, stops
      }
    }
  
    // --- 5h. Update assignment ---
    int new_cluster;                     // Declares the variable with chosen cluster ID
  
    if (sampled < K_minus) {
      // Assign to existing cluster (use original cluster ID)
      new_cluster = uniq_cl[sampled];   // If sampler choses a already existing cluster, look up for its real ID in uniq_c
      
    } else {
      // ---- New cluster: sample beta_{K+1}^E from its posterior ----
      
      //#Find next available cluster ID
      // New-ID search replicates R's:
      //   used_ids <- unique(xi)   # FULL xi, so xi[i] = old_cluster
      //
      // cnt currently holds xi_{-i}.  Adding {old_cluster} back gives
      // the same set as unique(full xi), because xi[i] = old_cluster.
      //
      // Replicate used_ids <- unique(xi)
      std::set<int> used_ids;        // keeps numbers sorted and unique. We want to create a list of take IDs
      for (const auto& kv : cnt) used_ids.insert(kv.first); // Collects IDs in cnt map and  // clusters in xi_E_{-i}
      used_ids.insert(old_cluster);  // add back old_cluster (xi_E[i] still at old value)
      
      const int max_id = *used_ids.rbegin();       // largest used ID
      //# Find first gap or use max+1
      int new_id       = max_id + 1;               // fallback: one past current max
      // Find first gap or use max+1 (replicates R's setdiff logic)
      for (int id = 1; id <= max_id; id++) {       // Prepares to use cluster label max_id + 1 (new). fallback: one past current max
        if (used_ids.find(id) == used_ids.end()) { // If gap exists,
          new_id = id;                             //uses the smallest missing label, (if a number isn't in the set, takes this as new ID)
          break;                                  //otherwise use max_id + 1
        }
      }
      
      // ---- Sample beta_new from posterior ----
      //
      // Precision: A = (1/kappa^2) x_i x_i^T + (1/sigma^2_E) I_{p_E}
      // Location:  b = (1/kappa^2) r_i^E x_i  + (1/sigma^2_E) mu_E
      // Posterior mean: M = A^{-1} b
      // Sample:    beta_new = M + R^{-1} z,  z ~ N(0, I_{p_E})
      //            where A = R^T R  (upper Cholesky)
      //
      // Note: r_i^E is used instead of t_i = y[i] - beta_int because the
      // residual already absorbs the M and D contributions
      
      // Precision matrix: A = (1/kappa^2) x_i x_i^T + (1/sigma^2_E) I_{p_E}
      const arma::mat A      = kappa2_inv * (x_i * x_i.t()) + sigma2_diag;  // R:  kappa2_inv * tcrossprod(x_i) + sigma2_diag
      
      // Upper Cholesky factor: A = R^T R  (matches R's chol())
      // Use cholesky decomposition instead of solve() for A
      const arma::mat R_chol = arma::chol(A, "upper");
      
      //  Mean vector: M = A^{-1} b
      // b = (1/kappa^2) r_i^E x_i + (1/sigma^2_E) mu_E
      // Note: only observation i is in this new cluster, so no sum needed
      // (Note: notice that the sum of kappa2_inv does not appear because we are sampling
      // based on observation i, the only observation in that cluster)
      const arma::vec b = kappa2_inv * r_i_E * x_i + sigma2_inv * mu_E;
      
      // M = A^{-1} b via two triangular solves
      // R's: backsolve(R, backsolve(R, b, transpose=TRUE))
      // Step 1: solve R^T w = b  (forward substitution, lower triangular)
      const arma::vec w  = arma::solve(arma::trimatl(R_chol.t()), b);
      // Step 2: solve R M = w    (back substitution, upper triangular)
      const arma::vec M_post = arma::solve(arma::trimatu(R_chol), w);
      
      // z ~ N(0, I_{p_E}) via R's RNG
      arma::vec z(p_E);                                         // Creates the p_E x 1 vector
      for (int j = 0; j < p_E; j++) z[j] = R::rnorm(0.0, 1.0);  // Draws p numbers from N(0,1)
      
      // beta_new = M + R^{-1} z  (R's: M + backsolve(R, z))
      // To get a sample from N(M, A^-1), we take the mean M and add the white noise z scaled by
      // the square root of the variance
      const arma::vec beta_new = M_post + arma::solve(arma::trimatu(R_chol), z);  // Gives a draw from N(M,A^-1)
      
      // Store new beta at new cluster ID
      beta_map_E[new_id] = beta_new;   // Saves the new coefficient vector under the new cluster ID
      new_cluster        = new_id;     // Assign observation to new cluster
    }
    
    // --- 5i. Commit: update count map and xi_vec ---
    cnt[new_cluster]++;            // Increments the count for the cluster that observation i is now assigned to
    xi_vec[i] = new_cluster;       // Updates the label vector xi
    
  } // end main sweep
  
  // =========================================================================
  // 6. Final cleanup: relabel clusters consecutively 1, 2, ..., K_E_final
  //    Done ONCE at the end, not inside the loop
  //    Mirrors R's as.integer(factor(xi_E, levels=used_clusters, ...))
  // =========================================================================
  
  // Collect used labels (sorted, unique)
  const std::set<int>    used_set(xi_vec.begin(), xi_vec.end());      // R equivalent sort(unique(xi_E))
  const std::vector<int> used_vec(used_set.begin(), used_set.end());  // Copies from the set to a vector, already sorted
  // Compute the final number of clusters
  const int K_final = static_cast<int>(used_vec.size());              // R equivalent length(used_clusters)
  
  // Remapping from old labels to new labels.
  // Old cluster ID  ->  new consecutive 1-based ID
  // i.e used_vec[0] -> new label 1  
  // mapping old cluster D to new consecutive cluster ID; # Relabel xi to be 1, 2, 3, ..., K_E
  std::map<int, int> remap;                 // R equivalent: factor(xi, levels = used_clusters, labels = seq_along(used_clusters))   
  for (int idx = 0; idx < K_final; idx++) { // Iterates over indexes 1,...,k_final (K)
    remap[used_vec[idx]] = idx + 1;         // 1-based output (R convention)
  }
  
  // beta_k_E_clean: Rcpp::List of length K_final (consecutive, 1-based in R)
  Rcpp::List beta_k_E_clean(K_final);                                // creates a list of length k_final (K_E) (R: vector("list", K_final))
  for (int idx = 0; idx < K_final; idx++) {                          // Iterates over k_E
    const arma::vec& bk = beta_map_E.at(used_vec[idx]);              // Retrieves the coefficient vector associated to k (R: beta_k[[used_clusters[idx]]])
    beta_k_E_clean[idx] = Rcpp::NumericVector(bk.begin(), bk.end()); // Converts an Armadillo vector to a numeric vector
  }
  
  // xi_E_final: relabelled, 1-based
  // # Relabel xi to be 1, 2, 3, ..., K
  // Equivalent to R: xi_final <- as.integer(factor(xi, levels = used_clusters, labels = seq_along(used_clusters)))
  
  Rcpp::IntegerVector xi_E_final(n);     // Creates an integer vector of length n
  for (int i = 0; i < n; i++) {          // Iterates over obervations
    xi_E_final[i] = remap.at(xi_vec[i]); // For each observation replaces old cluster label with the new consecutive label
  }
  
  // R equivalent: return(list(xi_E = xi_E_final, beta_k_E = beta_k_E_clean))
  return Rcpp::List::create(                  // Returns a list
    Rcpp::Named("xi_E")     = xi_E_final,     // xi vector 
    Rcpp::Named("beta_k_E") = beta_k_E_clean  // beta_k list
  );
}





// =============================================================================
// sample_xi_M_cpp.cpp
//
// C++ / RcppArmadillo translation of sample_xi_M() (R function).
// M-level cluster-assignment updating for a three-DP Dirichlet Process model.
//
//   Outputs:
//   list with:
//     - xi_M    : updated M-level cluster assignments (relabelled 1,...,K_M)
//     - beta_k_M: updated list of M-level cluster coefficients
// =============================================================================


// [[Rcpp::export]]
Rcpp::List sample_xi_M_cpp(
    const arma::vec&  y,            // (n)        response variable
    const arma::mat&  Z_full,       // (n x p_M)  municipality-level covariate matrix
    arma::ivec        xi_M,         // (n)        current M-level cluster assignments (1-based)
    Rcpp::List        beta_k_M,     // (K_M)      M-level cluster coefficients (consecutive 1-based)
    double            sigma2_M_base,//            base measure variance sigma^2_M
    double            alpha_M,      //            M-level DP concentration parameter
    const arma::vec&  mu_M,         // (p_M)      base measure mean
    const arma::mat&  X_E_full,     // (n x p_E)  individual-level covariate matrix
    const arma::mat&  W_full,       // (n x p_D)  department-level covariate matrix
    const arma::ivec& xi_E,         // (n)        current E-level cluster assignments (1-based)
    const arma::ivec& xi_D,         // (n)        current D-level cluster assignments (1-based)
    Rcpp::List        beta_k_E,     // (K_E)      E-level cluster coefficients (consecutive 1-based)
    Rcpp::List        beta_k_D,     // (K_D)      D-level cluster coefficients (consecutive 1-based)
    const arma::vec&  mun_kappa,    // (M_mun)    kappa^2_{j,q} per municipality, 0-based in C++
    const arma::ivec& obs_to_mun,   // (n)        municipality index per obs, 1-based (R convention)
    const arma::vec&  zTz_vec,      // (n)        pre-computed z_i^T z_i per observation
    const arma::vec&  sqrt_kappa,   // (M_mun)    pre-computed sqrt(kappa^2_{j,q}), 0-based in C++
    double            beta_int      //            global intercept beta
) {
  // -------------------------------------------------------------------------
  // 0. Dimensions and global constants
  // -------------------------------------------------------------------------
  const int    n           = static_cast<int>(y.n_elem);
  const int    p_M         = static_cast<int>(Z_full.n_cols);
  
  // Pre-compute constants used in new cluster beta sampling
  const double sigma2_inv  = 1.0 / sigma2_M_base;

  // sigma2_inv * I_{p_M} -- used every time a new cluster is born
  const arma::mat sigma2_diag = sigma2_inv * arma::eye<arma::mat>(p_M, p_M);

  // -------------------------------------------------------------------------
  // 1. Load beta_k_M into a map indexed by 1-based cluster ID
  //    Input is consecutive [1..K_M_init]; gaps can develop during the sweep
  // -------------------------------------------------------------------------
  std::map<int, arma::vec> beta_map_M;
  {
    const int K_M_init = static_cast<int>(beta_k_M.size());
    for (int k = 1; k <= K_M_init; k++) {
      beta_map_M[k] = Rcpp::as<arma::vec>(beta_k_M[k - 1]); // 0-based R list -> 1-based map
    }
  }

  // -------------------------------------------------------------------------
  // 2. Load beta_k_E and beta_k_D into maps for fast contrib lookup
  //    These are fixed during the M-level sweep (not updated here)
  // -------------------------------------------------------------------------
  std::map<int, arma::vec> beta_map_E;
  {
    const int K_E = static_cast<int>(beta_k_E.size());
    for (int k = 1; k <= K_E; k++) {
      beta_map_E[k] = Rcpp::as<arma::vec>(beta_k_E[k - 1]);
    }
  }

  std::map<int, arma::vec> beta_map_D;
  {
    const int K_D = static_cast<int>(beta_k_D.size());
    for (int k = 1; k <= K_D; k++) {
      beta_map_D[k] = Rcpp::as<arma::vec>(beta_k_D[k - 1]);
    }
  }
  
  // -------------------------------------------------------------------------
  // 3. Working copy of xi_M as std::vector<int> (1-based throughout)
  // -------------------------------------------------------------------------
  std::vector<int> xi_vec(n);
  for (int i = 0; i < n; i++) xi_vec[i] = xi_M[i];

  // -------------------------------------------------------------------------
  // 4. Running cluster-count map over all n observations
  //    Decrement/re-add observation i each iteration -> reflects xi_M_{-i}
  //    Replaces O(n) table(xi_M[-i]) calls with O(1) updates
  // -------------------------------------------------------------------------
  std::map<int, int> cnt;
  for (int i = 0; i < n; i++) cnt[xi_vec[i]]++;

  // =========================================================================
  // 5. Main sweep (sequential, i = 0 .. n-1, xi_M updated in place)
  // =========================================================================
  for (int i = 0; i < n; i++) {
  
    // --- 5a. Pre-computed values for observation i ---
    const int    mun0        = obs_to_mun[i] - 1;    // 0-based municipality index
    const double kappa2_i    = mun_kappa[mun0];
    const double sqrt_kap_i  = sqrt_kappa[mun0];
    const double kappa2_inv  = 1.0 / kappa2_i;
    const arma::vec z_i      = Z_full.row(i).t();    // (p_M x 1) column vector
    const double zTz_i       = zTz_vec[i];
  
    // --- 5b. Compute E-level and D-level contributions for observation i ---
    // These are fixed while we sample xi_M[i]
    // contrib_E = x_i^T beta_k_E  (individual contribution)
    const arma::vec x_i      = X_E_full.row(i).t();  // (p_E x 1)
    const int    ke_i        = xi_E[i];               // E-level cluster of obs i (1-based)
    const double contrib_E   = arma::dot(x_i, beta_map_E.at(ke_i));
  
    // contrib_D = w_q^T beta_k_D  (department contribution)
    const arma::vec w_i      = W_full.row(i).t();    // (p_D x 1)
    const int    kd_i        = xi_D[i];               // D-level cluster of obs i (1-based)
    const double contrib_D   = arma::dot(w_i, beta_map_D.at(kd_i));
  
    // --- 5c. Residual r_i^(M) = y_i - (beta + contrib_E + contrib_D) ---
    // Used for: (1) new cluster marginal likelihood
    //           (2) new beta_k_M posterior sampling
    const double r_i_M = y[i] - (beta_int + contrib_E + contrib_D);
  
    // --- 5d. Temporarily remove obs i from count map -> cnt = xi_M_{-i} ---
    const int old_cluster  = xi_vec[i];
    cnt[old_cluster]--;
    const bool old_empty   = (cnt[old_cluster] == 0);
    if (old_empty) cnt.erase(old_cluster);
  
    // --- 5e. Unique clusters and counts from xi_M_{-i} ---
    const int K_minus = static_cast<int>(cnt.size());
    std::vector<int> uniq_cl;  uniq_cl.reserve(K_minus);
    std::vector<int> uniq_cnt; uniq_cnt.reserve(K_minus);
    for (const auto& kv : cnt) {
      uniq_cl.push_back(kv.first);   // cluster ID
      uniq_cnt.push_back(kv.second); // count n_{-i, k}
    }
  
    // --- 5f. Log-probability vector: slots 0..K_minus-1 = existing, K_minus = new ---
    arma::vec log_probs(K_minus + 1);
  
    // - Probabilities for EXISTING clusters -
    // p(xi_i^M = k | -) propto n_{-i,k} * N(y_i | beta + contrib_E + z^T beta_k_M + contrib_D, kappa2)
    for (int k = 0; k < K_minus; k++) {
      const arma::vec& bk  = beta_map_M.at(uniq_cl[k]);
      // Mean: beta + contrib_E + z_i^T beta_k_M + contrib_D
      const double mean_k  = beta_int + contrib_E + arma::dot(z_i, bk) + contrib_D;
      // Log probability: log(n_{-i,k}) + log N(y_i | mean_k, kappa2_i)
      log_probs[k] = std::log(static_cast<double>(uniq_cnt[k]))
        + R::dnorm(y[i], mean_k, sqrt_kap_i, /*log=*/1);
    }
  
    // - Probability for NEW cluster -
    // p(xi_i^M = K+1 | -) propto alpha_M * N(r_i^M | z_i^T mu_M, kappa2 + sigma2_M * z_i^T z_i)
    const double var_new  = kappa2_i + (sigma2_M_base * zTz_i); // kappa2 + sigma2_M * z_i^T z_i
    const double mean_new = arma::dot(z_i, mu_M);               // z_i^T mu_M (evaluated on r_i^M)
    log_probs[K_minus] = std::log(alpha_M)
      + R::dnorm(r_i_M, mean_new, std::sqrt(var_new), /*log=*/1);
  
    // --- 5g. Normalize with log-sum-exp trick and sample ---
    const double max_lp = log_probs.max();
    arma::vec probs     = arma::exp(log_probs - max_lp);
    probs              /= arma::sum(probs);
  
    // Inverse-CDF sampling using R's RNG (equivalent to R's sample.int)
    const double u = R::unif_rand();
    int sampled    = K_minus;   // default: new cluster
    {
      double cum = 0.0;
      for (int k = 0; k <= K_minus; k++) {
        cum += probs[k];
        if (u <= cum) { sampled = k; break; }
      }
    }
  
    // --- 5h. Update assignment ---
    int new_cluster;
  
    if (sampled < K_minus) {
      // Assign to existing cluster (use original cluster ID)
      new_cluster = uniq_cl[sampled];
    
    } else {
      // ---- New cluster: sample beta_{K+1}^M from its posterior ----
      
      // Find next available cluster ID
      // Replicates R's: used_ids <- unique(xi_M)  (full xi_M, xi_M[i] = old_cluster)
      std::set<int> used_ids;
      for (const auto& kv : cnt) used_ids.insert(kv.first); // clusters in xi_M_{-i}
      used_ids.insert(old_cluster);  // add back old_cluster (xi_M[i] still at old value)
      
      const int max_id = *used_ids.rbegin(); // largest used ID
      int new_id       = max_id + 1;         // fallback: one past current max
      // Find first gap or use max+1 (replicates R's setdiff logic)
      for (int id = 1; id <= max_id; id++) {
        if (used_ids.find(id) == used_ids.end()) {
          new_id = id;
          break;
        }
      }
      
      // ---- Sample beta_new from posterior ----
      //
      // Precision: A = (1/kappa^2) z_i z_i^T + (1/sigma^2_M) I_{p_M}
      // Location:  b = (1/kappa^2) r_i^M z_i  + (1/sigma^2_M) mu_M
      // Posterior mean: M = A^{-1} b
      // Sample:    beta_new = M + R^{-1} z_rand,  z_rand ~ N(0, I_{p_M})
      //            where A = R^T R  (upper Cholesky)
      //
      // Note: r_i^M is used instead of t_i = y[i] - beta_int because the
      // residual already absorbs the E and D contributions
      
      // Precision matrix: A = (1/kappa^2) z_i z_i^T + (1/sigma^2_M) I_{p_M}
      const arma::mat A      = kappa2_inv * (z_i * z_i.t()) + sigma2_diag;
      
      // Upper Cholesky factor: A = R^T R  (matches R's chol())
      const arma::mat R_chol = arma::chol(A, "upper");
      
      // b = (1/kappa^2) r_i^M z_i + (1/sigma^2_M) mu_M
      // Note: only observation i is in this new cluster, so no sum needed
      const arma::vec b = kappa2_inv * r_i_M * z_i + sigma2_inv * mu_M;
      
      // M = A^{-1} b via two triangular solves
      // R's: backsolve(R, backsolve(R, b, transpose=TRUE))
      // Step 1: solve R^T w = b  (forward substitution, lower triangular)
      const arma::vec w      = arma::solve(arma::trimatl(R_chol.t()), b);
      // Step 2: solve R M = w    (back substitution, upper triangular)
      const arma::vec M_post = arma::solve(arma::trimatu(R_chol), w);
      
      // z_rand ~ N(0, I_{p_M}) via R's RNG
      arma::vec z_rand(p_M);
      for (int j = 0; j < p_M; j++) z_rand[j] = R::rnorm(0.0, 1.0);
      
      // beta_new = M + R^{-1} z_rand  (R's: M_post + backsolve(R, z_rand))
      const arma::vec beta_new = M_post + arma::solve(arma::trimatu(R_chol), z_rand);
      
      // Store new beta at new cluster ID
      beta_map_M[new_id] = beta_new;
      new_cluster        = new_id;
    }
    
    // --- 5i. Commit: update count map and xi_vec ---
    cnt[new_cluster]++;
    xi_vec[i] = new_cluster;
    
  } // end main sweep
  
  // =========================================================================
  // 6. Final cleanup: relabel clusters consecutively 1, 2, ..., K_M_final
  //    Done ONCE at the end, not inside the loop
  //    Mirrors R's as.integer(factor(xi_M, levels=used_clusters, ...))
  // =========================================================================
  
  // Collect used labels (sorted, unique)
  const std::set<int>    used_set(xi_vec.begin(), xi_vec.end()); // R: sort(unique(xi_M))
  const std::vector<int> used_vec(used_set.begin(), used_set.end());
  const int K_final = static_cast<int>(used_vec.size());
  
  // Remap old cluster IDs to new consecutive 1-based IDs
  std::map<int, int> remap;
  for (int idx = 0; idx < K_final; idx++) {
    remap[used_vec[idx]] = idx + 1; // 1-based output (R convention)
  }
  
  // beta_k_M_clean: Rcpp::List of length K_final (consecutive, 1-based in R)
  Rcpp::List beta_k_M_clean(K_final);
  for (int idx = 0; idx < K_final; idx++) {
    const arma::vec& bk = beta_map_M.at(used_vec[idx]);
    beta_k_M_clean[idx] = Rcpp::NumericVector(bk.begin(), bk.end());
  }
  
  // xi_M_final: relabelled, 1-based
  // R equivalent: as.integer(factor(xi_M, levels=used_clusters, labels=seq_along(used_clusters)))
  Rcpp::IntegerVector xi_M_final(n);
  for (int i = 0; i < n; i++) {
    xi_M_final[i] = remap.at(xi_vec[i]);
  }
  
  // R equivalent: return(list(xi_M = xi_M_final, beta_k_M = beta_k_M_clean))
  return Rcpp::List::create(
    Rcpp::Named("xi_M")     = xi_M_final,
    Rcpp::Named("beta_k_M") = beta_k_M_clean
  );
}



// =============================================================================
// sample_xi_D_cpp.cpp
//
// C++ / RcppArmadillo translation of sample_xi_D() (R function).
// D-level cluster-assignment updating for a three-DP Dirichlet Process model.
//
//   Outputs:
//   list with:
//     - xi_D    : updated D-level cluster assignments (relabelled 1,...,K_D)
//     - beta_k_D: updated list of D-level cluster coefficients
// =============================================================================


// [[Rcpp::export]]
Rcpp::List sample_xi_D_cpp(
    const arma::vec&  y,            // (n)        response variable
    const arma::mat&  W_full,       // (n x p_D)  department-level covariate matrix
    arma::ivec        xi_D,         // (n)        current D-level cluster assignments (1-based)
    Rcpp::List        beta_k_D,     // (K_D)      D-level cluster coefficients (consecutive 1-based)
    double            sigma2_D_base,//            base measure variance sigma^2_D
    double            alpha_D,      //            D-level DP concentration parameter
    const arma::vec&  mu_D,         // (p_D)      base measure mean
    const arma::mat&  X_E_full,     // (n x p_E)  individual-level covariate matrix
    const arma::mat&  Z_full,       // (n x p_M)  municipality-level covariate matrix
    const arma::ivec& xi_E,         // (n)        current E-level cluster assignments (1-based)
    const arma::ivec& xi_M,         // (n)        current M-level cluster assignments (1-based)
    Rcpp::List        beta_k_E,     // (K_E)      E-level cluster coefficients (consecutive 1-based)
    Rcpp::List        beta_k_M,     // (K_M)      M-level cluster coefficients (consecutive 1-based)
    const arma::vec&  mun_kappa,    // (M_mun)    kappa^2_{j,q} per municipality, 0-based in C++
    const arma::ivec& obs_to_mun,   // (n)        municipality index per obs, 1-based (R convention)
    const arma::vec&  wTw_vec,      // (n)        pre-computed w_i^T w_i per observation
    const arma::vec&  sqrt_kappa,   // (M_mun)    pre-computed sqrt(kappa^2_{j,q}), 0-based in C++
    double            beta_int      //            global intercept beta
) {
  // -------------------------------------------------------------------------
  // 0. Dimensions and global constants
  // -------------------------------------------------------------------------
  const int    n           = static_cast<int>(y.n_elem);
  const int    p_D         = static_cast<int>(W_full.n_cols);
  
  // Pre-compute constants used in new cluster beta sampling
  const double sigma2_inv  = 1.0 / sigma2_D_base;
  
  // sigma2_inv * I_{p_D} -- used every time a new cluster is born
  const arma::mat sigma2_diag = sigma2_inv * arma::eye<arma::mat>(p_D, p_D);
  
  // -------------------------------------------------------------------------
  // 1. Load beta_k_D into a map indexed by 1-based cluster ID
  //    Input is consecutive [1..K_D_init]; gaps can develop during the sweep
  // -------------------------------------------------------------------------
  std::map<int, arma::vec> beta_map_D;
  {
    const int K_D_init = static_cast<int>(beta_k_D.size());
    for (int k = 1; k <= K_D_init; k++) {
      beta_map_D[k] = Rcpp::as<arma::vec>(beta_k_D[k - 1]); // 0-based R list -> 1-based map
    }
  }
  
  // -------------------------------------------------------------------------
  // 2. Load beta_k_E and beta_k_M into maps for fast contrib lookup
  //    These are fixed during the D-level sweep (not updated here)
  // -------------------------------------------------------------------------
  std::map<int, arma::vec> beta_map_E;
  {
    const int K_E = static_cast<int>(beta_k_E.size());
    for (int k = 1; k <= K_E; k++) {
      beta_map_E[k] = Rcpp::as<arma::vec>(beta_k_E[k - 1]);
    }
  }
  
  std::map<int, arma::vec> beta_map_M;
  {
    const int K_M = static_cast<int>(beta_k_M.size());
    for (int k = 1; k <= K_M; k++) {
      beta_map_M[k] = Rcpp::as<arma::vec>(beta_k_M[k - 1]);
    }
  }
  
  // -------------------------------------------------------------------------
  // 3. Working copy of xi_D as std::vector<int> (1-based throughout)
  // -------------------------------------------------------------------------
  std::vector<int> xi_vec(n);
  for (int i = 0; i < n; i++) xi_vec[i] = xi_D[i];
  
  // -------------------------------------------------------------------------
  // 4. Running cluster-count map over all n observations
  //    Decrement/re-add observation i each iteration -> reflects xi_D_{-i}
  //    Replaces O(n) table(xi_D[-i]) calls with O(1) updates
  // -------------------------------------------------------------------------
  std::map<int, int> cnt;
  for (int i = 0; i < n; i++) cnt[xi_vec[i]]++;
  
  // =========================================================================
  // 5. Main sweep (sequential, i = 0 .. n-1, xi_D updated in place)
  // =========================================================================
  for (int i = 0; i < n; i++) {
    
    // --- 5a. Pre-computed values for observation i ---
    const int    mun0        = obs_to_mun[i] - 1;    // 0-based municipality index
    const double kappa2_i    = mun_kappa[mun0];
    const double sqrt_kap_i  = sqrt_kappa[mun0];
    const double kappa2_inv  = 1.0 / kappa2_i;
    const arma::vec w_i      = W_full.row(i).t();    // (p_D x 1) column vector
    const double wTw_i       = wTw_vec[i];
    
    // --- 5b. Compute E-level and M-level contributions for observation i ---
    // These are fixed while we sample xi_D[i]
    // contrib_E = x_i^T beta_k_E  (individual contribution)
    const arma::vec x_i      = X_E_full.row(i).t();  // (p_E x 1)
    const int    ke_i        = xi_E[i];               // E-level cluster of obs i (1-based)
    const double contrib_E   = arma::dot(x_i, beta_map_E.at(ke_i));
    
    // contrib_M = z_{j,q}^T beta_k_M  (municipality contribution)
    const arma::vec z_i      = Z_full.row(i).t();    // (p_M x 1)
    const int    km_i        = xi_M[i];               // M-level cluster of obs i (1-based)
    const double contrib_M   = arma::dot(z_i, beta_map_M.at(km_i));
    
    // --- 5c. Residual r_i^(D) = y_i - (beta + contrib_E + contrib_M) ---
    // Used for: (1) new cluster marginal likelihood
    //           (2) new beta_k_D posterior sampling
    const double r_i_D = y[i] - (beta_int + contrib_E + contrib_M);
    
    // --- 5d. Temporarily remove obs i from count map -> cnt = xi_D_{-i} ---
    const int old_cluster  = xi_vec[i];
    cnt[old_cluster]--;
    const bool old_empty   = (cnt[old_cluster] == 0);
    if (old_empty) cnt.erase(old_cluster);
    
    // --- 5e. Unique clusters and counts from xi_D_{-i} ---
    const int K_minus = static_cast<int>(cnt.size());
    std::vector<int> uniq_cl;  uniq_cl.reserve(K_minus);
    std::vector<int> uniq_cnt; uniq_cnt.reserve(K_minus);
    for (const auto& kv : cnt) {
      uniq_cl.push_back(kv.first);   // cluster ID
      uniq_cnt.push_back(kv.second); // count n_{-i, k}
    }
    
    // --- 5f. Log-probability vector: slots 0..K_minus-1 = existing, K_minus = new ---
    arma::vec log_probs(K_minus + 1);
    
    // - Probabilities for EXISTING clusters -
    // p(xi_i^D = k | -) propto n_{-i,k} * N(y_i | beta + contrib_E + contrib_M + w^T beta_k_D, kappa2)
    for (int k = 0; k < K_minus; k++) {
      const arma::vec& bk  = beta_map_D.at(uniq_cl[k]);
      // Mean: beta + contrib_E + contrib_M + w_i^T beta_k_D
      const double mean_k  = beta_int + contrib_E + contrib_M + arma::dot(w_i, bk);
      // Log probability: log(n_{-i,k}) + log N(y_i | mean_k, kappa2_i)
      log_probs[k] = std::log(static_cast<double>(uniq_cnt[k]))
        + R::dnorm(y[i], mean_k, sqrt_kap_i, /*log=*/1);
    }
    
    // - Probability for NEW cluster -
    // p(xi_i^D = K+1 | -) propto alpha_D * N(r_i^D | w_i^T mu_D, kappa2 + sigma2_D * w_i^T w_i)
    const double var_new  = kappa2_i + (sigma2_D_base * wTw_i); // kappa2 + sigma2_D * w_i^T w_i
    const double mean_new = arma::dot(w_i, mu_D);               // w_i^T mu_D (evaluated on r_i^D)
    log_probs[K_minus] = std::log(alpha_D)
      + R::dnorm(r_i_D, mean_new, std::sqrt(var_new), /*log=*/1);
    
    // --- 5g. Normalize with log-sum-exp trick and sample ---
    const double max_lp = log_probs.max();
    arma::vec probs     = arma::exp(log_probs - max_lp);
    probs              /= arma::sum(probs);
    
    // Inverse-CDF sampling using R's RNG (equivalent to R's sample.int)
    const double u = R::unif_rand();
    int sampled    = K_minus;   // default: new cluster
    {
      double cum = 0.0;
      for (int k = 0; k <= K_minus; k++) {
        cum += probs[k];
        if (u <= cum) { sampled = k; break; }
      }
    }
    
    // --- 5h. Update assignment ---
    int new_cluster;
    
    if (sampled < K_minus) {
      // Assign to existing cluster (use original cluster ID)
      new_cluster = uniq_cl[sampled];
      
    } else {
      // ---- New cluster: sample beta_{K+1}^D from its posterior ----
      
      // Find next available cluster ID
      // Replicates R's: used_ids <- unique(xi_D)  (full xi_D, xi_D[i] = old_cluster)
      std::set<int> used_ids;
      for (const auto& kv : cnt) used_ids.insert(kv.first); // clusters in xi_D_{-i}
      used_ids.insert(old_cluster);  // add back old_cluster (xi_D[i] still at old value)
      
      const int max_id = *used_ids.rbegin(); // largest used ID
      int new_id       = max_id + 1;         // fallback: one past current max
      // Find first gap or use max+1 (replicates R's setdiff logic)
      for (int id = 1; id <= max_id; id++) {
        if (used_ids.find(id) == used_ids.end()) {
          new_id = id;
          break;
        }
      }
      
      // ---- Sample beta_new from posterior ----
      //
      // Precision: A = (1/kappa^2) w_i w_i^T + (1/sigma^2_D) I_{p_D}
      // Location:  b = (1/kappa^2) r_i^D w_i  + (1/sigma^2_D) mu_D
      // Posterior mean: M = A^{-1} b
      // Sample:    beta_new = M + R^{-1} z_rand,  z_rand ~ N(0, I_{p_D})
      //            where A = R^T R  (upper Cholesky)
      //
      // Note: r_i^D is used instead of t_i = y[i] - beta_int because the
      // residual already absorbs the E and M contributions
      
      // Precision matrix: A = (1/kappa^2) w_i w_i^T + (1/sigma^2_D) I_{p_D}
      const arma::mat A      = kappa2_inv * (w_i * w_i.t()) + sigma2_diag;
      
      // Upper Cholesky factor: A = R^T R  (matches R's chol())
      const arma::mat R_chol = arma::chol(A, "upper");
      
      // b = (1/kappa^2) r_i^D w_i + (1/sigma^2_D) mu_D
      // Note: only observation i is in this new cluster, so no sum needed
      const arma::vec b = kappa2_inv * r_i_D * w_i + sigma2_inv * mu_D;
      
      // M = A^{-1} b via two triangular solves
      // R's: backsolve(R, backsolve(R, b, transpose=TRUE))
      // Step 1: solve R^T w = b  (forward substitution, lower triangular)
      const arma::vec ww     = arma::solve(arma::trimatl(R_chol.t()), b);
      // Step 2: solve R M = w    (back substitution, upper triangular)
      const arma::vec M_post = arma::solve(arma::trimatu(R_chol), ww);
      
      // z_rand ~ N(0, I_{p_D}) via R's RNG
      arma::vec z_rand(p_D);
      for (int j = 0; j < p_D; j++) z_rand[j] = R::rnorm(0.0, 1.0);
      
      // beta_new = M + R^{-1} z_rand  (R's: M_post + backsolve(R, z_rand))
      const arma::vec beta_new = M_post + arma::solve(arma::trimatu(R_chol), z_rand);
      
      // Store new beta at new cluster ID
      beta_map_D[new_id] = beta_new;
      new_cluster        = new_id;
    }
    
    // --- 5i. Commit: update count map and xi_vec ---
    cnt[new_cluster]++;
    xi_vec[i] = new_cluster;
    
  } // end main sweep
  
  // =========================================================================
  // 6. Final cleanup: relabel clusters consecutively 1, 2, ..., K_D_final
  //    Done ONCE at the end, not inside the loop
  //    Mirrors R's as.integer(factor(xi_D, levels=used_clusters, ...))
  // =========================================================================
  
  // Collect used labels (sorted, unique)
  const std::set<int>    used_set(xi_vec.begin(), xi_vec.end()); // R: sort(unique(xi_D))
  const std::vector<int> used_vec(used_set.begin(), used_set.end());
  const int K_final = static_cast<int>(used_vec.size());
  
  // Remap old cluster IDs to new consecutive 1-based IDs
  std::map<int, int> remap;
  for (int idx = 0; idx < K_final; idx++) {
    remap[used_vec[idx]] = idx + 1; // 1-based output (R convention)
  }
  
  // beta_k_D_clean: Rcpp::List of length K_final (consecutive, 1-based in R)
  Rcpp::List beta_k_D_clean(K_final);
  for (int idx = 0; idx < K_final; idx++) {
    const arma::vec& bk = beta_map_D.at(used_vec[idx]);
    beta_k_D_clean[idx] = Rcpp::NumericVector(bk.begin(), bk.end());
  }
  
  // xi_D_final: relabelled, 1-based
  // R equivalent: as.integer(factor(xi_D, levels=used_clusters, labels=seq_along(used_clusters)))
  Rcpp::IntegerVector xi_D_final(n);
  for (int i = 0; i < n; i++) {
    xi_D_final[i] = remap.at(xi_vec[i]);
  }
  
  // R equivalent: return(list(xi_D = xi_D_final, beta_k_D = beta_k_D_clean))
  return Rcpp::List::create(
    Rcpp::Named("xi_D")     = xi_D_final,
    Rcpp::Named("beta_k_D") = beta_k_D_clean
  );
}





// =============================================================================

//                  SOLVE LABEL SWITCHING PROBLEM

// =============================================================================

// =============================================================================
// solve_label_switching_3DP_cpp.cpp   (cost-matrix optimized version)
//
// C++ / RcppArmadillo translation of solve_label_switching_3DP() (R function).
// Sequential greedy label-switching correction for the three-level DP model.
// Each level (E, M, D) is corrected independently in turn, holding the other
// two fixed -- with E first, then M (using corrected E), then D (using both).
//
// =============================================================================
// NAMING NOTE
// =============================================================================
// The variables holding the per-cluster prediction matrices are named
// `Mpred_E`, `Mpred_M`, `Mpred_D` (NOT `M_E`, `M_M`, `M_D`).
//
// The reason: `M_E` is a STANDARD MACRO from <cmath> / <math.h> for Euler's
// number (~2.71828). Since RcppArmadillo pulls in cmath, the preprocessor
// textually replaces every `M_E` in source code with the numeric constant.
// Trying to declare `arma::mat M_E = ...` then becomes
// `arma::mat 2.71828... = ...` which the compiler rightly rejects.
//
// Same precaution applies to M_PI, M_LN2, M_SQRT2, etc.  As a rule, avoid
// `M_*` identifiers in numerical C++ code.
//
// =============================================================================
// OPTIMIZATION NOTES (vs. the previous solve_label_switching_3DP_cpp)
// =============================================================================
// The previous version called compute_contrib() inside each per-permutation
// loop, which built a fresh (n x p) `beta_per_obs` matrix and a length-n
// row-sum for every permutation. For k_mode = 7 this is 5040 per level x 3
// levels = 15,120 of these per iteration, each touching ~n*p doubles.
//
// This version applies the same cost-matrix decomposition that we already
// used in solve_label_switching_cpp, INDEPENDENTLY to each of the three steps.
//
// For STEP 1 (E-level), with M and D held at their original order, define:
//
//      r_i^E      = y_i - beta_int_b - contrib_M_orig[i] - contrib_D_orig[i]
//      Mpred_E    = X_E_full * beta_E_matrix^T              (n x k_E_mode)
//
// where Mpred_E(i, j) = x_i^T * beta_E^{(j)}. Then for any permutation sigma_E:
//
//      sum_i (y_i - y_hat_i(sigma_E))^2
//          = sum_i (r_i^E - Mpred_E(i, sigma_E[xi_E_i]))^2
//          = sum_c sum_{i: xi_E_i = c} (r_i^E - Mpred_E(i, sigma_E[c]))^2
//          = sum_c C_E(c, sigma_E[c])
//
// where C_E is a k_E_mode x k_E_mode cost matrix built ONCE per iteration.
// Scoring any permutation then costs only k_E_mode additions.
//
// Identical argument for STEPS 2 and 3 (with Mpred_M, r_M, C_M and
// Mpred_D, r_D, C_D).
//
// A nice side-benefit: contrib_M_orig, contrib_D_orig, contrib_E_corrected and
// contrib_M_corrected are all just LOOKUPS into the Mpred_* matrices we
// already computed for cost-matrix construction -- so compute_contrib() is
// dropped entirely. The only matrix-matrix products per iteration are the
// three BLAS GEMM calls that build Mpred_E, Mpred_M, Mpred_D once each.
//
// Per-iteration FLOP cost drops by roughly 700x for k=7, p~10, mirroring the
// single-DP optimization, but now applied three times.
//
// =============================================================================
// STATISTICAL EQUIVALENCE
// =============================================================================
//  * Permutation generation is unchanged (same generate_permutations_cpp,
//    same std::random_device seeding for k > 7). Permutation index s has the
//    same meaning as before for each level.
//  * The objective minimised at each step is identical to the original MSE up
//    to the positive constant 1/n, which does not affect argmin.
//  * Floating-point summation order differs slightly, so individual objective
//    values may differ from the original MSE * n by ~1e-13 relative. The
//    selected permutation (best_s_E, best_s_M, best_s_D) is the same in all
//    non-pathological cases. Tie-breaking matches arma::index_min() because
//    we use strict `<`.
//  * The sequential-greedy semantics are preserved exactly:
//        Step 1: corrects E with M and D at original order
//        Step 2: corrects M with E corrected, D at original order
//        Step 3: corrects D with E and M both corrected
//  * The beta_k reorder and xi relabel blocks are identical to the previous
//    version. In particular, the xi relabel uses
//        xi_*_correct_order(b, i) = best_sigma_*[xi_*_b[i] - 1] + 1
//    -- this is the same direct-indexing convention as the previous C++
//    version (which differs from the match()-based inversion in the R code).
//    Preserved here as-is; do NOT silently switch conventions.
//
// Compile from R with:
//   Rcpp::sourceCpp("solve_label_switching_3DP_cpp.cpp")
// =============================================================================

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <algorithm>   // std::next_permutation, std::shuffle
#include <limits>      // std::numeric_limits
#include <numeric>     // std::iota
#include <random>      // std::mt19937, std::random_device
#include <set>         // std::set for unique permutation deduplication
#include <vector>

// -----------------------------------------------------------------------------
// Helper: generate_permutations_cpp   (UNCHANGED from previous version)
// Generates 0-indexed permutations for a given k.
// k <= 7 : exact enumeration via std::next_permutation (all k! permutations)
// k >  7 : 5040 unique random permutations via std::shuffle
//          (uses std::random_device, NOT seed-reproducible from R's set.seed())
// -----------------------------------------------------------------------------
static std::vector<std::vector<int>> generate_permutations_cpp(int k) {
  
  std::vector<std::vector<int>> perms;
  std::vector<int> perm(k);
  std::iota(perm.begin(), perm.end(), 0);  // fill 0, 1, ..., k-1
   
  if (k <= 7) {
    // Enumerate all k! permutations in lexicographic order
    do {
      perms.push_back(perm);
    } while (std::next_permutation(perm.begin(), perm.end()));
    
  } else { 
    // Sample exactly 5040 unique random permutations
    int n_permu = 5040;                          // 7! cap, same as R
    std::set<std::vector<int>> perms_set;
    std::mt19937 gen(std::random_device{}());
     
    while ((int)perms_set.size() < n_permu) {
      std::shuffle(perm.begin(), perm.end(), gen);
      perms_set.insert(perm);
    } 
    perms.assign(perms_set.begin(), perms_set.end());
  }
   
  return perms;
} 

// -----------------------------------------------------------------------------
// Helper: argmin_over_permutations
// Given a square cost matrix C (k x k) and a list of 0-indexed permutations,
// returns the index s* of the permutation minimising sum_c C(c, sigma[c]).
//
// Uses strict `<` so ties are broken by FIRST-MIN (matches arma::index_min
// and therefore the original mse_*.index_min() semantics).
// -----------------------------------------------------------------------------
static int argmin_over_permutations(const arma::mat& C,
                                    const std::vector<std::vector<int>>& permu,
                                    int k_mode) {
  const int n_permu = static_cast<int>(permu.size());
  double best_total = std::numeric_limits<double>::infinity();
  int    best_s     = 0;
   
  for (int s = 0; s < n_permu; s++) {
    const std::vector<int>& sigma = permu[s];
    double total = 0.0;
    for (int c = 0; c < k_mode; c++) {
      total += C(c, sigma[c]);
    } 
    if (total < best_total) {
      best_total = total;
      best_s     = s;
    } 
  }
  return best_s;
} 

// -----------------------------------------------------------------------------
// Helper: build_cost_matrix
// Builds the k x k cost matrix
//
//      C(c, j) = sum over i with xi[i] = c of (r[i] - Mpred(i, j))^2
//
// xi_0based : 0-indexed cluster assignments (vector of int, length n)
// r         : (n) residual vector
// Mpred     : (n x k) precomputed contribution matrix
// -----------------------------------------------------------------------------
static arma::mat build_cost_matrix(const std::vector<int>& xi_0based,
                                   const arma::vec&        r,
                                   const arma::mat&        Mpred,
                                   int                     k_mode) {
  const int n = static_cast<int>(r.n_elem);
  arma::mat C(k_mode, k_mode, arma::fill::zeros);
   
  for (int i = 0; i < n; i++) {
    const int    c  = xi_0based[i];
    const double ri = r[i];
    // Innermost loop is over k_mode (small, e.g. 7). Compiler should vectorise.
    for (int j = 0; j < k_mode; j++) {
      const double d = ri - Mpred(i, j);
      C(c, j) += d * d;
    }
  } 
  return C;
} 

// =============================================================================
// Main function: solve_label_switching_3DP_cpp
// =============================================================================

// [[Rcpp::export]]
Rcpp::List solve_label_switching_3DP_cpp(
    Rcpp::List  cadena_filtered,    // filtered MCMC chain (R list)
    arma::vec   y,                  // response vector (n)
    arma::mat   X_E_full,           // individual-level covariates (n x p_E)
    arma::mat   Z_full,             // municipal-level covariates  (n x p_M)
    arma::mat   W_full,             // department-level covariates (n x p_D)
    int         k_E_mode,           // E-level modal number of clusters
    int         k_M_mode,           // M-level modal number of clusters
    int         k_D_mode) {         // D-level modal number of clusters
  
  // Clone cadena_filtered so we do not modify the original R object
  Rcpp::List result = Rcpp::clone(cadena_filtered);
  
  // ---------------------------------------------------------------------------
  // 1. Extract components from result
  // ---------------------------------------------------------------------------
  Rcpp::IntegerMatrix xi_E_rmat = result["xi_E"];     // n_iter x n  (1-based)
  Rcpp::IntegerMatrix xi_M_rmat = result["xi_M"];
  Rcpp::IntegerMatrix xi_D_rmat = result["xi_D"];
  
  const int n_iter = xi_E_rmat.nrow();
  const int n      = xi_E_rmat.ncol();
  
  Rcpp::NumericVector beta_int_vec = result["beta_int"];   // length n_iter
  
  Rcpp::List beta_k_E_all = result["beta_k_E"];
  Rcpp::List beta_k_M_all = result["beta_k_M"];
  Rcpp::List beta_k_D_all = result["beta_k_D"];
  
  // ---------------------------------------------------------------------------
  // 2. Generate permutations once per level (outside the iteration loop)
  // ---------------------------------------------------------------------------
  std::vector<std::vector<int>> permu_E = generate_permutations_cpp(k_E_mode);
  std::vector<std::vector<int>> permu_M = generate_permutations_cpp(k_M_mode);
  std::vector<std::vector<int>> permu_D = generate_permutations_cpp(k_D_mode);
  
  Rcpp::Rcout << "Solving label switching problem (3 DP, cost-matrix version)\n";
  Rcpp::Rcout << "  Filtered iterations: " << n_iter << "\n";
  Rcpp::Rcout << "  k_E_mode: " << k_E_mode
              << "  | Permutations: " << permu_E.size() << "\n";
  Rcpp::Rcout << "  k_M_mode: " << k_M_mode
              << "  | Permutations: " << permu_M.size() << "\n";
  Rcpp::Rcout << "  k_D_mode: " << k_D_mode
              << "  | Permutations: " << permu_D.size() << "\n\n";
  
  // ---------------------------------------------------------------------------
  // 3. Storage for corrected outputs
  // ---------------------------------------------------------------------------
  Rcpp::List beta_k_E_correct_order(n_iter);
  Rcpp::List beta_k_M_correct_order(n_iter);
  Rcpp::List beta_k_D_correct_order(n_iter);
  
  Rcpp::IntegerMatrix xi_E_correct_order(n_iter, n);
  Rcpp::IntegerMatrix xi_M_correct_order(n_iter, n);
  Rcpp::IntegerMatrix xi_D_correct_order(n_iter, n);
  
  // =========================================================================
  // 4. Main loop over filtered iterations
  // =========================================================================
  for (int b = 0; b < n_iter; b++) {
    
    if ((b + 1) % 100 == 0) {
      Rcpp::Rcout << "  Processing iteration " << (b + 1)
                  << " / " << n_iter << "\n";
      R_FlushConsole();
    }
    
    // -------------------------------------------------------------------------
    // 4a. Extract iteration-b quantities
    //     xi_*_b_0       : 0-based cluster assignments (length n)
    //     beta_int_b     : scalar intercept
    //     beta_k_*_b     : list of k_*_mode numeric vectors
    // -------------------------------------------------------------------------
    std::vector<int> xi_E_b_0(n), xi_M_b_0(n), xi_D_b_0(n);
    for (int i = 0; i < n; i++) {
      xi_E_b_0[i] = xi_E_rmat(b, i) - 1;     // 1-based -> 0-based (once per iteration)
      xi_M_b_0[i] = xi_M_rmat(b, i) - 1;
      xi_D_b_0[i] = xi_D_rmat(b, i) - 1;
    }
    
    const double beta_int_b = beta_int_vec[b];
    
    Rcpp::List beta_k_E_b = Rcpp::as<Rcpp::List>(beta_k_E_all[b]);
    Rcpp::List beta_k_M_b = Rcpp::as<Rcpp::List>(beta_k_M_all[b]);
    Rcpp::List beta_k_D_b = Rcpp::as<Rcpp::List>(beta_k_D_all[b]);
    
    // -------------------------------------------------------------------------
    // 4b. Build the three beta matrices (k_mode x p), row k = beta_k_b[[k]]
    //     Identical to the original do.call(rbind, beta_k_b)
    // -------------------------------------------------------------------------
    const int p_E = Rcpp::as<arma::rowvec>(beta_k_E_b[0]).n_elem;
    const int p_M = Rcpp::as<arma::rowvec>(beta_k_M_b[0]).n_elem;
    const int p_D = Rcpp::as<arma::rowvec>(beta_k_D_b[0]).n_elem;
    
    arma::mat beta_E_matrix(k_E_mode, p_E);
    arma::mat beta_M_matrix(k_M_mode, p_M);
    arma::mat beta_D_matrix(k_D_mode, p_D);
    for (int k = 0; k < k_E_mode; k++)
      beta_E_matrix.row(k) = Rcpp::as<arma::rowvec>(beta_k_E_b[k]);
    for (int k = 0; k < k_M_mode; k++)
      beta_M_matrix.row(k) = Rcpp::as<arma::rowvec>(beta_k_M_b[k]);
    for (int k = 0; k < k_D_mode; k++)
      beta_D_matrix.row(k) = Rcpp::as<arma::rowvec>(beta_k_D_b[k]);
     
    // -------------------------------------------------------------------------
    // 4c. (NEW) Precompute the three per-cluster contribution matrices ONCE.
    //
    //     Mpred_E(i, j) = x_i^T beta_E^{(j)}    (n x k_E_mode)
    //     Mpred_M(i, j) = z_i^T beta_M^{(j)}    (n x k_M_mode)
    //     Mpred_D(i, j) = w_i^T beta_D^{(j)}    (n x k_D_mode)
    //
    //     Each is a single BLAS-backed matrix multiplication. From here on,
    //     every cluster contribution we need -- original or permuted, for
    //     scoring or for the final corrected fit -- is just a LOOKUP.
    //
    //     IMPORTANT: do NOT rename these to M_E / M_M / M_D. `M_E` is a
    //     <cmath> macro for Euler's number; the preprocessor would replace
    //     it with the numeric constant 2.71828... and break compilation.
    // -------------------------------------------------------------------------
    const arma::mat Mpred_E = X_E_full * beta_E_matrix.t();   // n x k_E_mode
    const arma::mat Mpred_M = Z_full   * beta_M_matrix.t();   // n x k_M_mode
    const arma::mat Mpred_D = W_full   * beta_D_matrix.t();   // n x k_D_mode
    
    // -------------------------------------------------------------------------
    // 4d. (NEW) Original M and D contributions via lookup
    //     contrib_M_orig[i] = Mpred_M(i, xi_M_i)   (identity permutation on M)
    //     contrib_D_orig[i] = Mpred_D(i, xi_D_i)   (identity permutation on D)
    //
    //     Replaces the two compute_contrib(... id_M / id_D) calls in the
    //     previous version, at zero cost beyond the Mpred_M / Mpred_D GEMMs.
    // -------------------------------------------------------------------------
    arma::vec contrib_M_orig(n), contrib_D_orig(n);
    for (int i = 0; i < n; i++) {
      contrib_M_orig[i] = Mpred_M(i, xi_M_b_0[i]);
      contrib_D_orig[i] = Mpred_D(i, xi_D_b_0[i]);
    }
    
    // =========================================================================
    // STEP 1: Correct E-level (M and D held at original order)
    //
    // Score for sigma_E:
    //   sum_i (y_i - beta_int - Mpred_E(i, sigma_E[xi_E_i])
    //                         - Mpred_M_orig - Mpred_D_orig)^2
    //   = sum_c C_E(c, sigma_E[c])
    // with r_E_i = y_i - beta_int - contrib_M_orig[i] - contrib_D_orig[i].
    // =========================================================================
    
    // (NEW) Residual vector for E-level scoring
    const arma::vec r_E = y - beta_int_b - contrib_M_orig - contrib_D_orig;
    
    // (NEW) Build cost matrix and find argmin over permutations
    const arma::mat C_E = build_cost_matrix(xi_E_b_0, r_E, Mpred_E, k_E_mode);
    const int best_s_E  = argmin_over_permutations(C_E, permu_E, k_E_mode);
    const std::vector<int>& best_sigma_E = permu_E[best_s_E];
    
    // ---- Reorder beta_k_E (identical to previous version) -------------------
    Rcpp::List beta_k_E_corrected(k_E_mode);
    for (int k = 0; k < k_E_mode; k++)
      beta_k_E_corrected[k] = beta_k_E_b[best_sigma_E[k]];
    beta_k_E_correct_order[b] = beta_k_E_corrected;
    
    // ---- Relabel xi_E (identical to previous version) -----------------------
    //
    // Convention preserved EXACTLY as in the previous C++:
    //     xi_E_correct_order(b, i) = best_sigma_E[xi_E_b[i] - 1] + 1
    //
    // NOTE: the R version uses match(xi_E_b, best_sigma_E) which would give
    // a different (inverse-permutation) result. Preserved as-is on purpose
    // -- do not change without explicit user confirmation.
    for (int i = 0; i < n; i++) {
      xi_E_correct_order(b, i) = best_sigma_E[xi_E_b_0[i]] + 1;
    }
    
    // ---- Compute corrected E contribution (used in Steps 2 and 3) -----------
    //     contrib_E_corrected[i] = Mpred_E(i, best_sigma_E[xi_E_i])
    arma::vec contrib_E_corrected(n);
    for (int i = 0; i < n; i++) {
      contrib_E_corrected[i] = Mpred_E(i, best_sigma_E[xi_E_b_0[i]]);
    }
    
    // =========================================================================
    // STEP 2: Correct M-level (E corrected, D at original order)
    // =========================================================================
    
    const arma::vec r_M = y - beta_int_b - contrib_E_corrected - contrib_D_orig;
    
    const arma::mat C_M = build_cost_matrix(xi_M_b_0, r_M, Mpred_M, k_M_mode);
    const int best_s_M  = argmin_over_permutations(C_M, permu_M, k_M_mode);
    const std::vector<int>& best_sigma_M = permu_M[best_s_M];
    
    Rcpp::List beta_k_M_corrected(k_M_mode);
    for (int k = 0; k < k_M_mode; k++)
      beta_k_M_corrected[k] = beta_k_M_b[best_sigma_M[k]];
    beta_k_M_correct_order[b] = beta_k_M_corrected;
    
    for (int i = 0; i < n; i++) {
      xi_M_correct_order(b, i) = best_sigma_M[xi_M_b_0[i]] + 1;
    }
    
    // Corrected M contribution for Step 3
    arma::vec contrib_M_corrected(n);
    for (int i = 0; i < n; i++) {
      contrib_M_corrected[i] = Mpred_M(i, best_sigma_M[xi_M_b_0[i]]);
    }
    
    // =========================================================================
    // STEP 3: Correct D-level (E and M both corrected)
    // =========================================================================
    
    const arma::vec r_D = y - beta_int_b - contrib_E_corrected - contrib_M_corrected;
    
    const arma::mat C_D = build_cost_matrix(xi_D_b_0, r_D, Mpred_D, k_D_mode);
    const int best_s_D  = argmin_over_permutations(C_D, permu_D, k_D_mode);
    const std::vector<int>& best_sigma_D = permu_D[best_s_D];
    
    Rcpp::List beta_k_D_corrected(k_D_mode);
    for (int k = 0; k < k_D_mode; k++)
      beta_k_D_corrected[k] = beta_k_D_b[best_sigma_D[k]];
    beta_k_D_correct_order[b] = beta_k_D_corrected;
    
    for (int i = 0; i < n; i++) {
      xi_D_correct_order(b, i) = best_sigma_D[xi_D_b_0[i]] + 1;
    }
    
  } // end main iteration loop
  
  // ---------------------------------------------------------------------------
  // 5. Add corrected outputs to result  (unchanged)
  // ---------------------------------------------------------------------------
  result["beta_k_E_correct_order"] = beta_k_E_correct_order;
  result["beta_k_M_correct_order"] = beta_k_M_correct_order;
  result["beta_k_D_correct_order"] = beta_k_D_correct_order;
  
  result["xi_E_correct_order"] = xi_E_correct_order;
  result["xi_M_correct_order"] = xi_M_correct_order;
  result["xi_D_correct_order"] = xi_D_correct_order;
  
  Rcpp::Rcout << "\n cadena_filtered$beta_k_E_correct_order successfully added.\n";
  Rcpp::Rcout << " cadena_filtered$beta_k_M_correct_order successfully added.\n";
  Rcpp::Rcout << " cadena_filtered$beta_k_D_correct_order successfully added.\n";
  Rcpp::Rcout << " cadena_filtered$xi_E_correct_order successfully added.\n";
  Rcpp::Rcout << " cadena_filtered$xi_M_correct_order successfully added.\n";
  Rcpp::Rcout << " cadena_filtered$xi_D_correct_order successfully added.\n";
  
  return result;
}

// =============================================================================

// sample_kappa2_jq_cpp.cpp


// =============================================================================


#include <map>
#include <vector>

// [[Rcpp::export]]
arma::vec sample_kappa2_jq_cpp(
    const arma::vec&  y,            // (n)        response variable
    const arma::mat&  X_E_full,     // (n x p_E)  individual-level covariate matrix
    const arma::mat&  Z_full,       // (n x p_M)  municipality-level covariate matrix
    const arma::mat&  W_full,       // (n x p_D)  department-level covariate matrix
    const arma::ivec& xi_E,         // (n)        E-level cluster assignments (1-based)
    const arma::ivec& xi_M,         // (n)        M-level cluster assignments (1-based)
    const arma::ivec& xi_D,         // (n)        D-level cluster assignments (1-based)
    Rcpp::List        beta_k_E,     // (K_E)      E-level cluster coefficients (1-based)
    Rcpp::List        beta_k_M,     // (K_M)      M-level cluster coefficients (1-based)
    Rcpp::List        beta_k_D,     // (K_D)      D-level cluster coefficients (1-based)
    double            beta_int,     //            global intercept beta
    Rcpp::List        mun_map,      //            municipality mapping list (q,j,start,end,n)
    const arma::vec&  kappa2_q,     // (m)        department-level hyperparameters kappa^2_q
    double            nu_kappa      //            hyperparameter nu_kappa
) {
  // -------------------------------------------------------------------------
  // 0. Dimensions
  // -------------------------------------------------------------------------
  const int n_mun_total = static_cast<int>(mun_map.size());
  
  // -------------------------------------------------------------------------
  // 1. Load beta_k_E, beta_k_M, beta_k_D into maps indexed by 1-based cluster ID
  //    Mirrors the approach in sample_xi_*_cpp for consistent lookup
  // -------------------------------------------------------------------------
  std::map<int, arma::vec> beta_map_E;
  {
    const int K_E = static_cast<int>(beta_k_E.size());
    for (int k = 1; k <= K_E; k++)
      beta_map_E[k] = Rcpp::as<arma::vec>(beta_k_E[k - 1]);
  }

  std::map<int, arma::vec> beta_map_M;
  {
    const int K_M = static_cast<int>(beta_k_M.size());
    for (int k = 1; k <= K_M; k++)
      beta_map_M[k] = Rcpp::as<arma::vec>(beta_k_M[k - 1]);
  }

  std::map<int, arma::vec> beta_map_D;
  {
    const int K_D = static_cast<int>(beta_k_D.size());
    for (int k = 1; k <= K_D; k++)
      beta_map_D[k] = Rcpp::as<arma::vec>(beta_k_D[k - 1]);
  }

  // -------------------------------------------------------------------------
  // 2. Output vector: kappa2_jq[mun_idx] = updated kappa^2_{j,q}
  //    Equivalent to R's: kappa2_jq <- numeric(n_mun_total)
  // -------------------------------------------------------------------------
  arma::vec kappa2_jq_out(n_mun_total);

  // =========================================================================
  // 3. Main loop over municipalities (mirrors R's for (mun_idx in seq_len(n_mun_total)))
  // =========================================================================
  for (int mun_idx = 0; mun_idx < n_mun_total; mun_idx++) {
  
    // --- Extract municipality info from mun_map ---
    // mun_map is a list of lists with named elements: q, j, start, end, n
    Rcpp::List mun_info = Rcpp::as<Rcpp::List>(mun_map[mun_idx]);
  
    const int q   = Rcpp::as<int>(mun_info["q"]);      // department index (1-based)
    // j is stored but not used in the computation (only q is needed for kappa2_q[q])
    const int n_jq  = Rcpp::as<int>(mun_info["n"]);    // number of observations in mun
    const int start = Rcpp::as<int>(mun_info["start"]); // first obs index (1-based, R convention)
    const int end   = Rcpp::as<int>(mun_info["end"]);   // last  obs index (1-based, R convention)
  
    // Convert to 0-based C++ indices
    const int start0 = start - 1;
    const int end0   = end   - 1;
    
    // --- Extract data for municipality (j, q) ---
    // y_jq: response vector for observations in this municipality (length n_jq)
    // Equivalent to R's: y_jq <- y[idx_jq]
    const arma::vec y_jq = y.subvec(start0, end0);     // length n_jq
  
    // --- E-level contribution: x_i^T beta_k_E for each i in municipality ---
    // Equivalent to R's:
    //   beta_E_mat <- do.call(rbind, beta_k_E[xi_E[idx_jq]])  # n_jq x p_E
    //   contrib_E  <- rowSums(X_E_full[idx_jq, ] * beta_E_mat)
    arma::vec contrib_E(n_jq);
    for (int r = 0; r < n_jq; r++) {
      const int obs_r  = start0 + r;                          // 0-based obs index
      const int k_E_r  = xi_E[obs_r];                         // 1-based E-cluster of obs r
      const arma::vec& b_E = beta_map_E.at(k_E_r);            // p_E-vector
      contrib_E[r] = arma::dot(X_E_full.row(obs_r).t(), b_E); // x_i^T beta_k_E
    }
  
    // --- M-level contribution: z_{j,q}^T beta_k_M for each i in municipality ---
    // Equivalent to R's:
    //   beta_M_mat <- do.call(rbind, beta_k_M[xi_M[idx_jq]])  # n_jq x p_M
    //   contrib_M  <- rowSums(Z_full[idx_jq, ] * beta_M_mat)
    arma::vec contrib_M(n_jq);
    for (int r = 0; r < n_jq; r++) {
      const int obs_r  = start0 + r;
      const int k_M_r  = xi_M[obs_r];                         // 1-based M-cluster of obs r
      const arma::vec& b_M = beta_map_M.at(k_M_r);            // p_M-vector
      contrib_M[r] = arma::dot(Z_full.row(obs_r).t(), b_M);   // z_i^T beta_k_M
    }
  
    // --- D-level contribution: w_q^T beta_k_D for each i in municipality ---
    // Equivalent to R's:
    //   beta_D_mat <- do.call(rbind, beta_k_D[xi_D[idx_jq]])  # n_jq x p_D
    //   contrib_D  <- rowSums(W_full[idx_jq, ] * beta_D_mat)
    arma::vec contrib_D(n_jq);
    for (int r = 0; r < n_jq; r++) {
      const int obs_r  = start0 + r;
      const int k_D_r  = xi_D[obs_r];                         // 1-based D-cluster of obs r
      const arma::vec& b_D = beta_map_D.at(k_D_r);            // p_D-vector
      contrib_D[r] = arma::dot(W_full.row(obs_r).t(), b_D);   // w_i^T beta_k_D
    }
  
    // --- Full fitted values: vartheta_{i,j,q} ---
    // Equivalent to R's: fitted_jq <- beta_int + contrib_E + contrib_M + contrib_D
    const arma::vec fitted_jq = beta_int + contrib_E + contrib_M + contrib_D; // length n_jq
  
    // --- Residuals: (y_{j,q} - vartheta_{j,q}) ---
    // Equivalent to R's: residuals_jq <- y_jq - fitted_jq
    const arma::vec residuals_jq = y_jq - fitted_jq;           // length n_jq
  
    // --- Sum of squared residuals ---
    // Equivalent to R's: sum_sq_resid <- sum(residuals_jq^2)
    const double sum_sq_resid = arma::dot(residuals_jq, residuals_jq);
  
    // -----------------------------------------------------------------------
    // Posterior parameters for Inverse Gamma:
    // shape: (nu_kappa + n_{j,q}) / 2
    // rate : (nu_kappa * kappa2_q[q] + sum_sq_resid) / 2
    // -----------------------------------------------------------------------
    // kappa2_q is 1-based from R; q is 1-based -> use q-1 for 0-based C++ index
    const double shape_post = (nu_kappa + static_cast<double>(n_jq)) / 2.0;
    const double rate_post  = (nu_kappa * kappa2_q[q - 1] + sum_sq_resid)  / 2.0;
  
    // --- Sample from Inverse Gamma: IG(a, b) = 1 / Gamma(a, rate=b) ---
    // R's rgamma(shape, scale) where scale = 1/rate
    // Equivalent to R's: 1 / rgamma(1, shape = shape_post, rate = rate_post)
    kappa2_jq_out[mun_idx] = 1.0 / R::rgamma(shape_post, 1.0 / rate_post);
    
  } // end municipality loop

  // Return numeric vector of length n_mun_total
  // Equivalent to R's: return(kappa2_jq)
  return kappa2_jq_out;
}






// =============================================================================

// sample_kappa2_q_cpp.cpp
//
// =============================================================================


// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::vec sample_kappa2_q_cpp(
    const arma::vec&  mun_kappa,    // (n_mun_total) current municipality variances kappa^2_{j,q}
    Rcpp::List        dept_to_mun,  // (m)           dept_to_mun[[q]] = flat municipality indices in dept q (1-based)
    double            nu_kappa,     //               hyperparameter nu_kappa
    double            alpha_kappa,  //               prior shape alpha_kappa
    double            beta_kappa,   //               prior rate beta_kappa
    int               m             //               number of departments
) {
  // -------------------------------------------------------------------------
  // Output vector: kappa2_q[q] = updated kappa^2_q for department q
  // Equivalent to R's: kappa2_q <- numeric(m)
  // -------------------------------------------------------------------------
  arma::vec kappa2_q_out(m);

  // =========================================================================
  // Main loop over departments (mirrors R's for (q in seq_len(m)))
  // =========================================================================
  for (int q = 0; q < m; q++) {
  
    // --- Get pre-computed municipality indices for department q ---
    // dept_to_mun[[q]] in R (1-based) -> dept_to_mun[q] in C++ (0-based list index)
    // Values inside are 1-based municipality indices (R convention)
    // Equivalent to R's: mun_indices_q <- dept_to_mun[[q]]
    const Rcpp::IntegerVector mun_indices_q =
      Rcpp::as<Rcpp::IntegerVector>(dept_to_mun[q]);
  
    // --- Number of municipalities n_q in department q ---
    // Equivalent to R's: n_q <- length(mun_indices_q)
    const int n_q = mun_indices_q.size();
  
    // --- Sum of 1/kappa^2_{j,q} for all municipalities j in department q ---
    // Equivalent to R's: sum_inv_kappa <- sum(1 / mun_kappa[mun_indices_q])
    // mun_indices_q values are 1-based; subtract 1 for 0-based C++ indexing
    double sum_inv_kappa = 0.0;
    for (int idx = 0; idx < n_q; idx++) {
      const int mun_idx_0based = mun_indices_q[idx] - 1;   // convert 1-based to 0-based
      sum_inv_kappa += 1.0 / mun_kappa[mun_idx_0based];
    }
  
    // -----------------------------------------------------------------------
    // Posterior parameters for Gamma distribution:
    // shape: (nu_kappa * n_q + alpha_kappa) / 2
    // rate : beta_kappa/2 + (nu_kappa/2) * sum_{j=1}^{n_q} 1/kappa^2_{j,q}
    // -----------------------------------------------------------------------
    // Equivalent to R's:
    //   shape_post <- (nu_kappa * n_q + alpha_kappa) / 2
    //   rate_post  <- (beta_kappa / 2) + (nu_kappa / 2) * sum_inv_kappa
    const double shape_post = (nu_kappa * static_cast<double>(n_q) + alpha_kappa) / 2.0;
    const double rate_post  = (beta_kappa / 2.0) + (nu_kappa / 2.0) * sum_inv_kappa;
  
    // --- Sample from Gamma(shape, rate) ---
    // R's rgamma(1, shape = shape_post, rate = rate_post)
    // R::rgamma takes (shape, scale) where scale = 1/rate
    // Equivalent to R's: kappa2_q[q] <- rgamma(1, shape = shape_post, rate = rate_post)
    kappa2_q_out[q] = R::rgamma(shape_post, 1.0 / rate_post);
    
  } // end department loop

  // Return numeric vector of length m
  // Equivalent to R's: return(kappa2_q)
  return kappa2_q_out;
} 



