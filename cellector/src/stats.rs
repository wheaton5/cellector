extern crate statrs;
extern crate logaddexp;

use self::logaddexp::LogSumExp;

pub fn expected_log_beta_binomial_pmf(total_count: usize, alpha: f64, beta: f64, ln_coefficient: f64) -> f64 {
    return (0..(total_count + 1)).map(|x| log_beta_binomial_pmf(x as f64, (total_count - x) as f64, alpha, beta, ln_coefficient)).ln_sum_exp(); // do I know how this work? no. does it work? no idea.
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
        for k in 0..n {
            k_log_binomial_coefficients.push(statrs::function::factorial::ln_binomial(n as u64, k as u64) as f64);
        }
        log_binomial_coefficients.push(k_log_binomial_coefficients);
    }
    return log_binomial_coefficients;
}


