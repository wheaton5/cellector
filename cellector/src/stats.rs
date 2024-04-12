extern crate statrs;

pub fn expected_log_beta_binomial_pmf(total_count: usize, alpha: f64, beta: f64, precomputed_log_binomial_coefficients: &Vec<Vec<f64>>) -> f64 {
    //let x = (0..(total_count + 1)).map(|k| log_beta_binomial_pmf(k as f64, (total_count - k) as f64, alpha, beta, precomputed_log_binomial_coefficients[total_count][k])).ln_sum_exp();
    let mut log_likelihoods_squared: Vec<f64> = Vec::new();
    for k in 0..(total_count + 1) {
        let log_binomial_coefficient;
        if total_count < precomputed_log_binomial_coefficients.len() {
            log_binomial_coefficient = precomputed_log_binomial_coefficients[total_count][k];
        } else {
            log_binomial_coefficient = statrs::function::factorial::ln_binomial(total_count as u64, k as u64) as f64;
        }
        log_likelihoods_squared.push(2.0*log_beta_binomial_pmf(k as f64, (total_count - k) as f64, alpha, beta, log_binomial_coefficient));
    }
    let mut to_return = log_likelihoods_squared[0];
    for k in 1..(total_count+1) {
        to_return = logsumexp(to_return, log_likelihoods_squared[k]);
    }
    return to_return;
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


