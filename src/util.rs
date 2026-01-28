pub mod arithmetic;

use num_complex::Complex32;

pub fn approx_eq_f32_slice(a: &[f32], b: &[f32], eps: f32) -> bool {
    if a.len() != b.len() {
        return false;
    }

    a.iter()
        .zip(b.iter())
        .all(|(&x, &y)| (x - y).abs() <= eps)
}

pub fn approx_eq_complex32_slice(
    a: &[Complex32],
    b: &[Complex32],
    eps: f32,
) -> bool {
    if a.len() != b.len() {
        return false;
    }

    a.iter()
        .zip(b.iter())
        .all(|(x, y)| (*x - *y).norm() <= eps)
}
