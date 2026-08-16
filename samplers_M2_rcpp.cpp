
// =============================================================================
// sample_xi_optimized.cpp
//
// C++ / RcppArmadillo translation of sample_xi_optimized() (R function).
// Cluster-assignment updating for a Dirichlet Process model.
//
// Compile from R with:
//   Rcpp::sourceCpp("sample_xi_optimized.cpp")
//
// Key design decisions vs. the R original:
//   1. A running std::map<int,int> cluster-count map replaces per-iteration
//      table(xi[-i]) calls, reducing the sweep from O(n^2) to O(n * K).
//   2. beta_k is kept in a std::map<int, arma::vec> keyed by 1-based cluster
//      IDs; gaps can exist during the sweep (same as R) and are cleaned up at
//      the end, exactly once.
//   3. All random draws go through R's RNG interface (R::unif_rand,
//      R::rnorm, R::dnorm) so results are reproducible under set.seed().
//   4. Cholesky-based back-substitution mirrors R's chol() + backsolve().
//   5. obs_to_mun is 1-based (R convention); -1 is applied before indexing
//      0-based C++ vectors.
//   6. When creating a new cluster the "used IDs" set is built from
//      cnt-map keys PLUS {old_cluster}, exactly replicating R's
//      unique(xi) call on the full xi (xi[i] is still at its old value).
// =============================================================================

// RcppArmadillo for interface between R and cpp, and lineal algebra (vectors and matrices)
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <map>
#include <set>
#include <vector>
#include <cmath>

// [[Rcpp::export]]
Rcpp::List sample_xi_optimized_cpp( // Result is a list
    const arma::vec&  y,            // (n)      response variable, numeric vector
    const arma::mat&  X_full,       // (n x p)  covariate matrix
    arma::ivec        xi,           // (n)      current cluster assignments, integer vector
    Rcpp::List        beta_k,       // (K)      cluster coefficients, consecutive 1-based, list
    double            sigma2_base,  //          variance of base distribution (sigma^2), double
    double            alpha,        //          DP concentration parameter, double
    const arma::vec&  mu_vec,       // (p)      mean of base distribution, numeric vector
    const arma::vec&  mun_kappa,    // (M)      kappa^2 per municipality, 0-based in C++, , numeric vector
    const arma::ivec& obs_to_mun,   // (n)      municipality index per obs, 1-based, integer vector
    const arma::vec&  xTx_vec,      // (n)      pre-computed x_i^T x_i, numeric vector
    const arma::vec&  sqrt_kappa,   // (M)      pre-computed sqrt(kappa^2), 0-based in C++, numeric vector
    double            beta_int      //          global intercept, double
) {
  // -------------------------------------------------------------------------
  // 0.  Dimensions and global constants
  // -------------------------------------------------------------------------
  const int    n          = static_cast<int>(y.n_elem);         // .n_elem: number of elements in Armadillo vector y
  const int    p          = static_cast<int>(X_full.n_cols);    // n_cols num of cols. static_cast<int> converts the size type to an integer.
  
  // Pre-compute constants used in new cluster sampling
  const double sigma2_inv = 1.0 / sigma2_base;
  
  // sigma2_inv * I_p  -- used every time a new cluster is born
  const arma::mat sigma2_diag = sigma2_inv * arma::eye<arma::mat>(p, p); //eye<arma::mat>(p, p) creates I_p
  
  // -------------------------------------------------------------------------
  // 1.  Load beta_k list into a map indexed by 1-based cluster ID
  //     input is consecutive [1..K_init], but during the sweep IDs can gap
  //     as clusters may disappear.
  // -------------------------------------------------------------------------
  
  // Some clusters may disappear during the algorithm, maping let us store coefficients indexed by their
  //explicit cluster IDs. Introduces a mapping for efficient cluster ID handling.
  std::map<int, arma::vec> beta_map; //Declares a C++ map
  {
    const int K_init = static_cast<int>(beta_k.size()); // number (integer) of initial cluster coefficient vectors (initial k)
    for (int k = 1; k <= K_init; k++) {                 // From k to initial K. Initialised in 1 because R code initialise 1-based
      beta_map[k] = Rcpp::as<arma::vec>(beta_k[k - 1]); // Rcpp::as<arma::vec> Convert the list into Armadillo vector
    }                                                   //beta_k[k - 1] transform to 0-based C++ indexing
  }
  
  // -------------------------------------------------------------------------
  // 2.  Working copy of xi as std::vector<int> (1-based throughout)
  // -------------------------------------------------------------------------
  std::vector<int> xi_vec(n);                      // declares a vector of length n (integer array)
  for (int i = 0; i < n; i++) xi_vec[i] = xi[i];   // Copies the current labels from xi into xi_vec
  
  // -------------------------------------------------------------------------
  // 3.  Running cluster-count map over all n observations.
  //     We temporarily decrement/re-add observation i each iteration,
  //     so at any point inside the loop the map reflects xi_{-i}.
  //     This replaces O(n) table(xi[-i]) calls -> O(1) updates.
  // -------------------------------------------------------------------------
  
  // Build a map from cluster ID to number of observation currently assigned to that cluster
  std::map<int, int> cnt;                       // stores ordered key–value pairs. cluster_id (1-based) ; cnt = count
  for (int i = 0; i < n; i++) cnt[xi_vec[i]]++; // increments the count for the cluster label of observation i.
  
  // =========================================================================
  // 4.  Main sweep  (sequential, i = 0 .. n-1, xi updated in place)
  // =========================================================================
  for (int i = 0; i < n; i++) {         // Iterate over i = 1,...,n
    
    // --- 4a. Pre-computed values for observation i ---
    const int    mun0      = obs_to_mun[i] - 1;   // 0-based municipality index (-1 to 1-based to 0-based)
    const double kappa2_i  = mun_kappa[mun0]; 
    const double sqrt_k_i  = sqrt_kappa[mun0];
    const arma::vec x_i    = X_full.row(i).t();   // extracts covariates i row and transforms to column vector (px1) (Armadillo is column oriented)
    const double xTx_i     = xTx_vec[i];          // Precomputed xTx for i
    
    // --- 4b. Get cluster assignments excluding observation i cnt  ->  cnt now = xi_{-i} ---
    const int old_cluster  = xi_vec[i];               // Retrieves current cluster of observation i 
    cnt[old_cluster]--;                               // Decrements the count in the mapping. -- subtracts 1
    const bool old_empty   = (cnt[old_cluster] == 0); // boolean that checks if the count is 0
    if (old_empty) cnt.erase(old_cluster);            // if count is zero (empty cluster), ereases it from cnt
    
    // --- 4c. Get unique clusters and counts from xi_{-i} ---
    // std::map iterates in sorted index order (equivalent to R's table())
    const int K_minus = static_cast<int>(cnt.size());      // Counts the number of unique clusters
    std::vector<int> uniq_cl;   uniq_cl.reserve(K_minus);  // Creates an array uniq_cl with k_minus elements
    std::vector<int> uniq_cnt;  uniq_cnt.reserve(K_minus); // Creates an array uniq_cnt with k_minus elements
    for (const auto& kv : cnt) {                           // iterates over every entry in cnt map
      uniq_cl.push_back(kv.first);                         // push_back adds elements to the end of the vector: i cluster ID 
      uniq_cnt.push_back(kv.second);                       // push_back adds elements to the end of the vector: i cluster count
    }
    
    // --- 4d. Log probabilities ---
    //Allocate log probability vector
    // Creates a vector with one slot per existing cluster plus one for new cluster
    arma::vec log_probs(K_minus + 1);   // last slot = new cluster
    
    // - Probabilities for Existing clusters -
    for (int k = 0; k < K_minus; k++) {                         // Iterates over clusters
      //Get beta for this cluster
      const arma::vec& bk = beta_map.at(uniq_cl[k]);            // Looks-up the coefficient vector for cluster uniq_cl[k] in beta_k(mapped)
      //Mean: beta_int + x_i^T beta_k
      const double mean_k = beta_int + arma::dot(x_i, bk);      // performs dot product between x_i and bk (covariates and betas)
      //Log probability: log(n_{-i,k}) + log N(y_i | mean_k, kappa^2)
      log_probs[k] = std::log(static_cast<double>(uniq_cnt[k])) // Prior probability 
        + R::dnorm(y[i], mean_k, sqrt_k_i, /*log=*/1);          // likelihood in log
    }
    
    // - Probability for NEW cluster -
    // # Variance: kappa^2_{j,q} + sigma^2 * x_i^T x_i
    const double var_new  = kappa2_i + (sigma2_base * xTx_i);
    // # Mean: beta_int + x_i^T mu
    const double mean_new = beta_int + arma::dot(x_i, mu_vec);
    // # Log probability: log(alpha) + log N(y_i | mean_new, var_new)
    log_probs[K_minus] = std::log(alpha)                        // Prior probability 
      + R::dnorm(y[i], mean_new, std::sqrt(var_new), 1);        // likelihood in log
    
    // --- 4e. Sample new cluster assignment ---
    //# Normalize probabilities (log-sum-exp trick for numerical stability)
    const double max_lp = log_probs.max();             // Finds the max value in log_probs
    arma::vec probs = arma::exp(log_probs - max_lp);   // Subtracts max value
    probs /= arma::sum(probs);                        // arma::sum(probs) and /= divides by the total sum
    
    // --- 4f. Sample cluster index (0-based) using R's RNG ---
    //        Index 0..K_minus-1 = existing clusters
    //        Index K_minus      = new cluster
    //        Equivalent to R: sample.int(K_new+1, ...) shifted by -1.
    const double u = R::unif_rand();  // random uniform betwen 0 and 1
    int sampled = K_minus;            // default: new cluster
    {  // manually implement inverse-CDF sampling from a uniform draw.
      double cum = 0.0;                       // Initialises a cummulative sum
      for (int k = 0; k <= K_minus; k++) {    // Iterates over k (including k_new)
        cum += probs[k];                      // Accumulates probability
        if (u <= cum) { sampled = k; break; } // Once where the number landed is found, stops
      }
    }
    
    // --- 4g. Update cluster assignment ---
    int new_cluster;                    // Declares the variable with chosen cluster ID
    
    if (sampled < K_minus) {
      // # Assign to existing cluster (use original cluster ID)
      new_cluster = uniq_cl[sampled];  // If sampler choses a already existing cluster, look up for its real ID in uniq_cl
      
    } else {
      // ---- Create a new cluster ----
      //#Find next available cluster ID
      // New-ID search replicates R's:
      //   used_ids <- unique(xi)   # FULL xi, so xi[i] = old_cluster
      //
      // cnt currently holds xi_{-i}.  Adding {old_cluster} back gives
      // the same set as unique(full xi), because xi[i] = old_cluster.
      //
      
      // Replicate used_ids <- unique(xi)
      std::set<int> used_ids;           // keeps numbers sorted and unique. We want to create a list of take IDs
      for (const auto& kv : cnt) used_ids.insert(kv.first);  // Collects IDs in cnt map and 
      used_ids.insert(old_cluster);   // xi[i] is still at its old value. Adds the xi_{-i} ID back (The one taken out)
      
      const int max_id = *used_ids.rbegin();        // finds max_id. rbegin ponts to the largest number in the set (since its sored)
      //# Find first gap or use max+1
      int new_id = max_id + 1;                      // Prepares to use cluster label max_id + 1 (new). fallback: one past current max
      for (int id = 1; id <= max_id; id++) {    
        if (used_ids.find(id) == used_ids.end()) {  // If gap exists,
          new_id = id;                              //uses the smallest missing label, (if a number isn't in the set, takes this as new ID)
          break;                                    //otherwise use max_id + 1
        }
      }
      
      // ---- Sample beta_{new} from its posterior ----
      //
      // Posterior derivation (same as R):
      //   Prior:      beta ~ N(mu, sigma^2 I)
      //   Likelihood: y_i | beta ~ N(beta_int + x_i^T beta, kappa^2_i)
      //
      //   Precision:  A = (1/kappa^2) x_i x_i^T + (1/sigma^2)I_p
      //   Location:   b = (1/kappa^2) t_i x_i  +  (1/sigma^2) mu
      //   Posterior mean: M = A^{-1} b
      //   Sample:     beta_new = M + A^{-1/2} z,  z ~ N(0,I)
      //             = M + R^{-1} z,  where A = R^T R  (Cholesky, upper)
      //

      const double t_i       = y[i] - beta_int;
      const double kappa2_inv = 1.0 / kappa2_i;
      
      // Precision matrix (positive definite)
      //# A = (1/kappa_j,q^2) x_i x_i^T + (1/sigma^2)I_p
      const arma::mat A = kappa2_inv * (x_i * x_i.t()) + sigma2_diag;   // Matrix
      
      // Use cholesky decomposition instead of solve() for A
      // Upper Cholesky factor R  (A = R^T R, matches R's chol())
      const arma::mat R_chol = arma::chol(A, "upper");
      
      //  Mean vector: M = A^{-1} b
      //# b = (1/kappa^2) t_i x_i + (1/sigma^2) mu
      // (Note: notice that the sum of kappa2_inv does not appear because we are sampling
      // based on observation i, the only observation in that cluster)
      const arma::vec b = kappa2_inv * t_i * x_i + sigma2_inv * mu_vec;
      
      // M = A^{-1} b  via two triangular solves
      //   R's: backsolve(R, backsolve(R, b, transpose=TRUE))
      //   Step 1: solve R^T w = b  (forward substitution, lower triangular)
      const arma::vec w = arma::solve(arma::trimatl(R_chol.t()), b);
      //   Step 2: solve R M = w    (back substitution, upper triangular)
      const arma::vec M = arma::solve(arma::trimatu(R_chol), w);
      
      // z ~ N(0, I_p)  Draw a p-vector of iid standar normal: via R's RNG
      arma::vec z(p);                                          // Creates the px1 vector
      for (int j = 0; j < p; j++) z[j] = R::rnorm(0.0, 1.0);   // Draws p numbers from N(0,1)
      
      // beta_new = M + R^{-1} z  (R's: M + backsolve(R, z))
      // To get a sample from N(M, A^-1), we take the mean M and add the white noise z scaled by
      // the square root of the variance
      const arma::vec beta_new = M + arma::solve(arma::trimatu(R_chol), z);  // Gives a draw from N(M,A^-1)
      
      //  Store new beta at new cluster ID
      beta_map[new_id] = beta_new;  // Saves the new coefficient vector under the new cluster ID
      new_cluster      = new_id;    // Assign observation to new cluster
    }
    
    // --- 4h. Commit: update count map and xi ---
    cnt[new_cluster]++;         // Increments the count for the cluster that observation i is now assigned to
    xi_vec[i] = new_cluster;    // Updates the label vector xi
    
  }   // end main sweep
  
  // =========================================================================
  // 5.  Final cleanup: relabel clusters consecutively  1, 2, ..., K_final
  //     Mirrors R's as.integer(factor(xi, levels=used_clusters, ...))
  //     Done ONCE at the end, not inside the loop.
  // =========================================================================
  // Collect used labels Creates a set that removes duplicates and sorts ascending
  const std::set<int>    used_set(xi_vec.begin(), xi_vec.end());     // R equivalent sort(unique(xi))
  const std::vector<int> used_vec(used_set.begin(), used_set.end()); // Copies from the set to a vector, already sorted
  // Compute the final number of clusters
  const int K_final = static_cast<int>(used_vec.size());             // R equivalent length(used_clusters)
  
  // Remapping from old labels to new labels.
  // Old cluster ID  ->  new consecutive 1-based ID
  // i.e used_vec[0] -> new label 1  
  // mapping old cluster D to new consecutive cluster ID; # Relabel xi to be 1, 2, 3, ..., K
  std::map<int, int> remap;                  // R equivalent: factor(xi, levels = used_clusters, labels = seq_along(used_clusters))   
  for (int idx = 0; idx < K_final; idx++) {  // Iterates over indexes 1,...,k_final (K)
    remap[used_vec[idx]] = idx + 1;          // Remap to R indexing 1-based (adds 1)
  }
  
  // beta_k_clean: Rcpp::List of length K_final (consecutive, 1-based in R)
  Rcpp::List beta_k_clean(K_final);                                // creates a list of length k_final (K) (R: vector("list", K_final))
  for (int idx = 0; idx < K_final; idx++) {                        // Iterates over k
    const arma::vec& bk = beta_map.at(used_vec[idx]);              // Retrieves the coefficient vector associated to k (R: beta_k[[used_clusters[idx]]])
    beta_k_clean[idx] = Rcpp::NumericVector(bk.begin(), bk.end()); // Converts an Armadillo vector to a numeric vector
  }
  
  // xi_final: relabelled, 1-based
  // # Relabel xi to be 1, 2, 3, ..., K
  // Equivalent to R: xi_final <- as.integer(factor(xi, levels = used_clusters, labels = seq_along(used_clusters)))
  Rcpp::IntegerVector xi_final(n);        // Creates an integer vector of length n
  for (int i = 0; i < n; i++) {           // Iterates over obervations
    xi_final[i] = remap.at(xi_vec[i]);    // For each observation replaces old cluster label with the new consecutive label
  }
  
  // R equivalent: return(list(xi = xi_final, beta_k = beta_k_clean))
  return Rcpp::List::create(              // Returns a list
    Rcpp::Named("xi")     = xi_final,     // xi vector 
    Rcpp::Named("beta_k") = beta_k_clean  // beta_k list
  );
}





// =============================================================================
// solve_label_switching_cpp.cpp  

// C++ / RcppArmadillo translation of solve_label_switching() (R function).
// Solves the label switching problem for a model with one Dirichlet Process
// by finding, for each posterior iteration, the permutation of cluster labels
// that minimises the MSE between fitted and observed values.
//
// =============================================================================
// OPTIMIZATION NOTES (vs. the previous version)
// =============================================================================
// The previous version recomputed, for every permutation sigma:
//      permuted_beta_matrix (k x p)
//      beta_per_obs         (n x p)         <-- expensive memory traffic
//      y_hat                (n)
//      MSE                  (scalar)
// at total cost O(n_permu * n * p) per iteration.
//
// This version exploits the additive decomposition of the MSE across clusters.
// Let:
//      r_i = y_i - beta_int_b                (residuals against intercept)
//      M   = X_full * beta_matrix^T          (n x k_mode matrix)
//
// where M(i, j) is the dot product x_i^T beta_j (i.e. the contribution to y_hat
// that observation i would receive if it were assigned to cluster j).
//
// Then for any permutation sigma:
//
//     sum_i (y_i - y_hat_i(sigma))^2
//         = sum_i (r_i - M(i, sigma[xi_i]))^2
//         = sum_c sum_{i : xi_i = c} (r_i - M(i, sigma[c]))^2
//         = sum_c C(c, sigma[c])
//
// where C(c, j) = sum_{i : xi_i = c} (r_i - M(i, j))^2 is a k x k cost matrix
// computed ONCE per iteration. The argmin over permutations of
// sum_c C(c, sigma[c]) equals the argmin of MSE(sigma) exactly (dividing by n
// does not affect argmin). Per-iteration cost drops to:
//
//      O(n * p * k)        for X_full * beta^T
//    + O(n * k)            to build C
//    + O(n_permu * k)      to score all permutations on C   <-- 7 adds/perm
//
// For n=40,000, p~10, k=7, n_permu=5040 this is roughly 700x fewer FLOPs and
// far less memory traffic than the original.
//
// =============================================================================
// STATISTICAL EQUIVALENCE
// =============================================================================
//  * The set and order of permutations enumerated is unchanged (same generator,
//    same seed). The permutation index `s` has the same meaning as before.
//  * The objective being minimised is mathematically identical to the original
//    MSE up to the constant factor 1/n, which does not affect argmin.
//    (We drop the 1/n; this is the only deviation from the original arithmetic.)
//  * Floating-point summation order differs slightly, so individual objective
//    values may differ from the original MSE * n by ~1e-13 relative. The
//    selected permutation `best_s` is the same in all non-pathological cases
//    (exact-tie behaviour matches the original: index_min picks the first one;
//     we use strict `<` for the same effect).
//  * xi_correct_order and beta_k_correct_order are constructed by exactly the
//    same blocks (4f / 4g) as before, so given the same `best_s` the outputs
//    are bit-identical.
//
// Compile from R with:
//   Rcpp::sourceCpp("solve_label_switching_cpp.cpp")
// =============================================================================

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <algorithm>   // std::next_permutation
#include <limits>      // std::numeric_limits
#include <numeric>     // std::iota
#include <random>      // std::mt19937, std::shuffle
#include <set>         // std::set for unique permutation deduplication
#include <vector>

// -----------------------------------------------------------------------------
// Helper: generate_permutations    (UNCHANGED from previous version)
// Generates 0-indexed permutations of {0, 1, ..., k-1}.
// k <= 7 : all k! permutations (exact, same set as gtools::permutations())
// k >  7 : exactly 5040 unique random permutations (same cap as R version)
// -----------------------------------------------------------------------------
static std::vector<std::vector<int>> generate_permutations(int k,
                                                           unsigned int seed = 777) {
  std::vector<std::vector<int>> perms;
  std::vector<int> perm(k);
  std::iota(perm.begin(), perm.end(), 0);  // fill 0, 1, ..., k-1
  
  if (k <= 7) {
    // Enumerate all k! permutations in lexicographic order
    // Equivalent to gtools::permutations(n=k, r=k) -- same set, same order
    do {
      perms.push_back(perm);
    } while (std::next_permutation(perm.begin(), perm.end()));
    
  } else { 
    // Sample exactly 5040 unique random permutations
    const int n_permu = 5040;            // 7! cap, same as R
    std::set<std::vector<int>> perms_set;
    std::mt19937 gen(seed);              // fixed seed for reproducibility
     
    while (static_cast<int>(perms_set.size()) < n_permu) {
      std::shuffle(perm.begin(), perm.end(), gen);
      perms_set.insert(perm);            // set guarantees uniqueness
    } 
    perms.assign(perms_set.begin(), perms_set.end());
  }
   
  return perms;
} 

// [[Rcpp::export]]
Rcpp::List solve_label_switching_cpp(
    Rcpp::List  cadena_filtered,  // filtered MCMC chain (R list)
    arma::vec   y,                // response vector (n)
    arma::mat   X_full,           // covariate matrix (n x p)
    int         k_mode            // modal number of clusters
) {
   
  // Clone cadena_filtered so we do not modify the original R object
  Rcpp::List result = Rcpp::clone(cadena_filtered);
   
  // ---------------------------------------------------------------------------
  // 1. Extract components from result
  //    xi      : IntegerMatrix  (n_iter x n), values 1-based (R convention)
  //    beta_k  : List[n_iter], each element is List[k_mode] of NumericVector
  //    beta_int: NumericVector of length n_iter
  // ---------------------------------------------------------------------------
  Rcpp::IntegerMatrix  xi_mat       = result["xi"];        // n_iter x n
  Rcpp::List           beta_k_all   = result["beta_k"];    // list of n_iter lists
  Rcpp::NumericVector  beta_int_vec = result["beta_int"];  // length n_iter
   
  const int n_iter = xi_mat.nrow();
  const int n      = xi_mat.ncol();
  const int p      = static_cast<int>(X_full.n_cols);
   
  // ---------------------------------------------------------------------------
  // 2. Generate permutations (0-indexed)
  // ---------------------------------------------------------------------------
  std::vector<std::vector<int>> permu = generate_permutations(k_mode, 777u);
  const int n_permu = static_cast<int>(permu.size());
   
  Rcpp::Rcout << "Solving label switching problem (cost-matrix version)\n";
  Rcpp::Rcout << "  Filtered iterations: " << n_iter  << "\n";
  Rcpp::Rcout << "  k_mode:              " << k_mode  << "\n";
  Rcpp::Rcout << "  Permutations:        " << n_permu << "\n\n";
   
  // ---------------------------------------------------------------------------
  // 3. Storage for corrected outputs
  // ---------------------------------------------------------------------------
  Rcpp::List          beta_k_correct_order(n_iter);
  Rcpp::IntegerMatrix xi_correct_order(n_iter, n);
   
  // =========================================================================
  // 4. Main loop over filtered iterations
  // =========================================================================
  for (int b = 0; b < n_iter; b++) {
     
    if ((b + 1) % 10 == 0) {
      Rcpp::Rcout << "  Processing iteration " << (b + 1)
                  << " / " << n_iter << "\n"; 
      R_FlushConsole();
    } 
    
    // -------------------------------------------------------------------------
    // 4a. Extract iteration-b quantities
    //     xi_b      : 0-based cluster assignments (vector of ints, length n)
    //     beta_int_b: scalar intercept
    //     beta_k_b  : list of k_mode numeric vectors, each length p
    // -------------------------------------------------------------------------
    std::vector<int> xi_b(n);
    for (int i = 0; i < n; i++) {
      xi_b[i] = xi_mat(b, i) - 1;  // convert to 0-based ONCE per iteration
    } 
    
    const double beta_int_b = beta_int_vec[b];
     
    Rcpp::List beta_k_b = Rcpp::as<Rcpp::List>(beta_k_all[b]);
     
    // -------------------------------------------------------------------------
    // 4b. Build beta matrix (k_mode x p): row k = beta_k_b[[k]]
    //     Identical to the original (R: beta_matrix <- do.call(rbind, beta_k_b))
    // -------------------------------------------------------------------------
    arma::mat beta_matrix(k_mode, p);
    for (int k = 0; k < k_mode; k++) {
      arma::rowvec bk = Rcpp::as<arma::rowvec>(beta_k_b[k]);
      beta_matrix.row(k) = bk;
    }
     
    // -------------------------------------------------------------------------
    // 4c. (NEW) Precompute the per-cluster contribution matrix and residuals.
    //
    //     M(i, j) = x_i^T beta_j         (size n x k_mode)
    //     r(i)    = y_i - beta_int_b     (size n)
    //
    //     Equivalence note:
    //       Original y_hat_i(sigma) = beta_int_b + x_i^T beta_{sigma[xi_i]}
    //                              = beta_int_b + M(i, sigma[xi_i])
    //       so (y_i - y_hat_i)^2 = (r_i - M(i, sigma[xi_i]))^2.
    //
    //     M is computed via a single BLAS-backed matrix multiplication, which
    //     is much faster than the original per-permutation expansion into a
    //     fresh (n x p) `beta_per_obs` matrix.
    // -------------------------------------------------------------------------
    const arma::mat M = X_full * beta_matrix.t();   // n x k_mode  (BLAS GEMM)
    const arma::vec r = y - beta_int_b;             // n
    
    // -------------------------------------------------------------------------
    // 4d. (NEW) Build the k x k cost matrix
    //
    //     C(c, j) = sum over i with xi_i = c of (r_i - M(i, j))^2
    //
    //     This is the crux of the optimisation: once C is built, scoring any
    //     permutation sigma takes only k additions:
    //         score(sigma) = sum_c C(c, sigma[c])
    // -------------------------------------------------------------------------
    arma::mat C(k_mode, k_mode, arma::fill::zeros);
    
    for (int i = 0; i < n; i++) {
      const int    c  = xi_b[i];     // already 0-based
      const double ri = r[i];
      // Innermost loop over k_mode (small, e.g. 7); compiler should vectorise.
      for (int j = 0; j < k_mode; j++) {
        const double d = ri - M(i, j);
        C(c, j) += d * d;
      }
    }
    
    // -------------------------------------------------------------------------
    // 4e. (NEW) Enumerate permutations on the cost matrix and pick the argmin.
    //
    //     This is the SAME minimisation as the original mse_permu.index_min(),
    //     because:
    //         n * MSE(sigma) = sum_i (y_i - y_hat_i(sigma))^2
    //                        = sum_c C(c, sigma[c])
    //     and argmin is invariant to the positive multiplicative constant n.
    //
    //     Tie-breaking: strict `<` keeps the FIRST permutation achieving the
    //     minimum -- matching arma::index_min()'s behaviour and therefore
    //     matching the original implementation exactly on ties.
    // -------------------------------------------------------------------------
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
    
    const std::vector<int>& best_sigma = permu[best_s];  // 0-indexed
    
    // -------------------------------------------------------------------------
    // 4f. Reorder beta_k according to best_sigma
    //     IDENTICAL to the original (block 4e in the previous version):
    //         R: beta_k_correct_order[[b]] <- beta_k_b[best_sigma]
    // -------------------------------------------------------------------------
    Rcpp::List beta_k_corrected(k_mode);
    for (int k = 0; k < k_mode; k++) {
      beta_k_corrected[k] = beta_k_b[best_sigma[k]];   // best_sigma[k] is 0-based
    }
    beta_k_correct_order[b] = beta_k_corrected;
    
    // -------------------------------------------------------------------------
    // 4g. Relabel xi according to best_sigma
    //     IDENTICAL to the original (block 4f in the previous version):
    //         R: xi_correct_order_list[[b]] <- best_sigma[xi_b]
    //
    //     xi_b is now 0-based (we did the -1 in 4a), and best_sigma is 0-based.
    //     We add 1 at the end to return to R's 1-based convention.
    // -------------------------------------------------------------------------
    for (int i = 0; i < n; i++) {
      xi_correct_order(b, i) = best_sigma[xi_b[i]] + 1;
    }
    
  } // end main iteration loop
   
  // ---------------------------------------------------------------------------
  // 5. Add corrected outputs to result (unchanged)
  // ---------------------------------------------------------------------------
  result["beta_k_correct_order"] = beta_k_correct_order;
  result["xi_correct_order"]     = xi_correct_order;
  
  Rcpp::Rcout << "\n cadena_filtered$beta_k_correct_order successfully added.\n";
  Rcpp::Rcout << " cadena_filtered$xi_correct_order successfully added.\n";
  
  return result;
}








