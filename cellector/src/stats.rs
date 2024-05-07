extern crate statrs;


pub struct LogBetaBinomialExpectation{
    pub expectation: f64,
    pub expectation_variance: f64,
}


pub fn expected_log_beta_binomial_pmf(total_count: usize, alpha: f64, beta: f64, precomputed_log_binomial_coefficients: &Vec<Vec<f64>>) -> LogBetaBinomialExpectation {
    let mut log_likelihoods: Vec<f64> = Vec::new();
    for k in 0..(total_count + 1) {
        let log_binomial_coefficient;
        if total_count < precomputed_log_binomial_coefficients.len() {
            log_binomial_coefficient = precomputed_log_binomial_coefficients[total_count][k];
        } else {
            log_binomial_coefficient = statrs::function::factorial::ln_binomial(total_count as u64, k as u64) as f64;
        }
        log_likelihoods.push(log_beta_binomial_pmf(k as f64, (total_count - k) as f64, alpha, beta, log_binomial_coefficient));
    }
    let mut expectation = 2.0*log_likelihoods[0]; // 2.0* in log space to square it
    for k in 1..(total_count+1) {
        expectation = logsumexp(expectation, 2.0*log_likelihoods[k]); // again 2.0* in log to square
    }
    let mut variance: f64 = 0.0; // um how do we compute in log space? or even should we?
    // I guess this doesnt need to be in log
    // so this should be sum_k=0 to n p(k,n|a,b)*(log(p(k,n|a,b) - E(log(D|a,b))))^2
    for log_likelihood in log_likelihoods { 
        variance += log_likelihood.exp() * (log_likelihood - expectation).powi(2);
    }
    return LogBetaBinomialExpectation{
        expectation: expectation,
        expectation_variance: variance,
    };
}


pub fn logsumexp(val_1: f64, val_2: f64) -> f64 {
    let max = val_1.max(val_2);
    let sum = (val_1 - max).exp() + (val_2 - max).exp();
    max + sum.ln()
}

pub fn log_beta_binomial_pmf(alt_count: f64, ref_count: f64, alpha: f64, beta: f64, ln_coefficient: f64) -> f64 {
    let ln_numerator = log_beta_calc(alt_count + alpha, ref_count + beta);
    let ln_denominator = log_beta_calc(alpha, beta);
    let log_likelihood = ln_coefficient + ln_numerator - ln_denominator;
    return log_likelihood;
}

pub fn log_beta_calc(alpha: f64, beta: f64) -> f64 {
    let log_gamma_alpha = statrs::function::gamma::ln_gamma(alpha);
    let log_gamma_beta = statrs::function::gamma::ln_gamma(beta);
    let log_gamma_alpha_beta = statrs::function::gamma::ln_gamma(alpha + beta);
    return log_gamma_alpha + log_gamma_beta - log_gamma_alpha_beta;
}

pub fn precompute_log_binomial_coefficients(max_n: usize) -> Vec<Vec<f64>> {
    let mut log_binomial_coefficients: Vec<Vec<f64>> = Vec::new();
    for n in 0..(max_n+1) {
        let mut k_log_binomial_coefficients: Vec<f64> = Vec::new();
        for k in 0..(n+1) {
            k_log_binomial_coefficients.push(statrs::function::factorial::ln_binomial(n as u64, k as u64) as f64);
        }
        log_binomial_coefficients.push(k_log_binomial_coefficients);
    }
    return log_binomial_coefficients;
}


