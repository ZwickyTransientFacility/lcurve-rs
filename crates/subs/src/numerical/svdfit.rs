use nalgebra::DMatrix;

/// Data point for SVD fitting.
pub struct FitDatum {
    pub x: f64,
    pub y: f64,
    /// Error (weight). Points with z <= 0 are ignored.
    pub z: f64,
}

/// Solve a linear least-squares problem using SVD.
///
/// Fits the model: y_i = sum_j a[j] * vect[i][j], weighted by 1/z_i.
/// Points with `z <= 0` are skipped.
///
/// * `data` - observations with x, y, z (error)
/// * `vect` - design matrix: `vect[i][j]` is value of basis function j at data point i
///
/// Returns `(coefficients, chisq)`.
pub fn svdfit(data: &[FitDatum], vect: &[Vec<f64>]) -> (Vec<f64>, f64) {
    let nc = if vect.is_empty() { 0 } else { vect[0].len() };
    let ndata = data.len();

    // Count valid data
    let valid: Vec<usize> = (0..ndata).filter(|&i| data[i].z > 0.0).collect();
    let ndat = valid.len();

    if ndat == 0 || nc == 0 {
        return (vec![0.0; nc], 0.0);
    }

    // Build weighted design matrix and RHS
    let mut u = DMatrix::zeros(ndat, nc);
    let mut b = vec![0.0; ndat];

    for (k, &i) in valid.iter().enumerate() {
        let tmp = 1.0 / data[i].z;
        for j in 0..nc {
            u[(k, j)] = tmp * vect[i][j];
        }
        b[k] = tmp * data[i].y;
    }

    // SVD
    let svd = u.svd(true, true);
    let u_mat = svd.u.unwrap();
    let v_mat = svd.v_t.unwrap().transpose();
    let s = svd.singular_values;

    // Edit singular values (threshold)
    let tol = 1.0e-5;
    let wmax = s.iter().copied().fold(0.0f64, f64::max);
    let thresh = tol * wmax;

    // Back-substitution: a = V * diag(1/w) * U^T * b
    let mut a = vec![0.0; nc];
    for j in 0..nc {
        if s[j] > thresh {
            let mut sum = 0.0;
            for k in 0..ndat {
                sum += u_mat[(k, j)] * b[k];
            }
            sum /= s[j];
            for i in 0..nc {
                a[i] += v_mat[(i, j)] * sum;
            }
        }
    }

    // Compute chi-squared on original data
    let mut chisq = 0.0;
    for &i in &valid {
        let mut sum = 0.0;
        for j in 0..nc {
            sum += a[j] * vect[i][j];
        }
        let diff = (data[i].y - sum) / data[i].z;
        chisq += diff * diff;
    }

    (a, chisq)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_linear_fit() {
        // Fit y = a + b*x to exact data: y = 2 + 3*x
        let data: Vec<FitDatum> = (0..10)
            .map(|i| {
                let x = i as f64;
                FitDatum {
                    x,
                    y: 2.0 + 3.0 * x,
                    z: 1.0,
                }
            })
            .collect();

        let vect: Vec<Vec<f64>> = data.iter().map(|d| vec![1.0, d.x]).collect();

        let (a, chisq) = svdfit(&data, &vect);

        assert!((a[0] - 2.0).abs() < 1e-10, "a[0]={}", a[0]);
        assert!((a[1] - 3.0).abs() < 1e-10, "a[1]={}", a[1]);
        assert!(chisq < 1e-20, "chisq={}", chisq);
    }

    #[test]
    fn test_negative_errors_skipped() {
        let data = vec![
            FitDatum { x: 0.0, y: 2.0, z: 1.0 },
            FitDatum { x: 1.0, y: 5.0, z: 1.0 },
            FitDatum { x: 2.0, y: 999.0, z: -1.0 }, // should be skipped
            FitDatum { x: 3.0, y: 11.0, z: 1.0 },
        ];

        let vect: Vec<Vec<f64>> = data.iter().map(|d| vec![1.0, d.x]).collect();

        let (a, _) = svdfit(&data, &vect);

        assert!((a[0] - 2.0).abs() < 1e-10);
        assert!((a[1] - 3.0).abs() < 1e-10);
    }
}
