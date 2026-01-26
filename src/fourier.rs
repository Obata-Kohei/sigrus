use num_complex::Complex32;
use crate::prelude::*;

pub fn dft(input: &[Complex32], output: &mut [Complex32]) {
    let n = input.len();
    assert_eq!(n, output.len());

    for k in 0..n {
        let mut sum = Complex32::ZERO;
        for i in 0..n {
            sum += input[i] * Complex32::from_polar(1.0, -TWO_PI * k as f32 * i as f32 / n as f32);
        }
        output[k] = sum
    }
}

pub fn idft(input: &[Complex32], output: &mut [Complex32]) {
    let n = input.len();
    assert_eq!(n, output.len());

    for i in 0..n {
        let mut sum = Complex32::ZERO;
        for k in 0..n {
            sum += input[k] * Complex32::from_polar(1.0, TWO_PI * k as f32 * i as f32 / n as f32);
        }
        output[i] = sum / n as f32;
    }
}

pub fn fft_cooley_tukey(input: &[Complex32], output: &mut [Complex32]) {
    let n = input.len();
    assert_eq!(n, output.len());
    assert!(n.is_power_of_two());

    // ベースケース
    if n == 1 {
        output[0] = input[0];
        return;
    }

    // 偶数・奇数に分解
    let mut even = Vec::with_capacity(n / 2);
    let mut odd  = Vec::with_capacity(n / 2);

    for i in 0..n / 2 {
        even.push(input[2 * i]);
        odd.push(input[2 * i + 1]);
    }

    // 再帰 FFT
    let mut even_fft = vec![Complex32::ZERO; n / 2];
    let mut odd_fft  = vec![Complex32::ZERO; n / 2];

    fft_cooley_tukey(&even, &mut even_fft);
    fft_cooley_tukey(&odd,  &mut odd_fft);

    // マージ（バタフライ）
    for k in 0..n / 2 {
        let theta = -2.0 * PI * k as f32 / n as f32;
        let w = Complex32::from_polar(1.0, theta);

        let t = w * odd_fft[k];
        output[k] = even_fft[k] + t;
        output[k + n / 2] = even_fft[k] - t;
    }
}

fn bit_reverse_permute(x: &mut [Complex32]) {
    let n = x.len();
    let bits = n.trailing_zeros();

    for i in 0..n {
        let j = i.reverse_bits() >> (usize::BITS - bits);
        if i < j {
            x.swap(i, j);
        }
    }
}

pub fn fft_inplace(x: &mut [Complex32]) {
    let n = x.len();
    assert!(n.is_power_of_two());

    // 1. ビット反転
    bit_reverse_permute(x);

    // 2. 段階的 FFT
    let mut m = 2;
    while m <= n {
        let half = m / 2;
        let theta = -2.0 * PI / m as f32;
        let wm = Complex32::from_polar(1.0, theta);

        for k in (0..n).step_by(m) {
            let mut w = Complex32::new(1.0, 0.0);
            for j in 0..half {
                let t = w * x[k + j + half];
                let u = x[k + j];
                x[k + j] = u + t;
                x[k + j + half] = u - t;
                w *= wm;
            }
        }
        m *= 2;
    }
}

fn make_full_twiddle_table(n: usize) -> Vec<Complex32> {
    assert!(n.is_power_of_two());

    let mut table = Vec::with_capacity(n / 2);
    for k in 0..n / 2 {
        let theta = -2.0 * PI * k as f32 / n as f32;
        table.push(Complex32::from_polar(1.0, theta));
    }
    table
}

pub fn fft_inplace_fulltable(x: &mut [Complex32], twiddle: &[Complex32]) {
    let n = x.len();
    assert!(n.is_power_of_two());
    assert_eq!(twiddle.len(), n / 2);

    bit_reverse_permute(x);

    let mut m = 2;
    while m <= n {
        let half = m / 2;
        let step = n / m;

        for k in (0..n).step_by(m) {
            for j in 0..half {
                let w = twiddle[j * step];

                let t = w * x[k + j + half];
                let u = x[k + j];

                x[k + j] = u + t;
                x[k + j + half] = u - t;
            }
        }
        m *= 2;
    }
}

fn make_twiddle_step(m: usize) -> Complex32 {
    let theta = -2.0 * PI / m as f32;
    Complex32::from_polar(1.0, theta)
}

pub fn fft_inplace_steptwiddle(buf: &mut [Complex32]) {
    let n = buf.len();
    assert!(n.is_power_of_two());

    // 並べ替え
    bit_reverse_permute(buf);

    let mut m = 2;
    while m <= n {
        let half = m / 2;
        let wm = make_twiddle_step(m);

        for k in (0..n).step_by(m) {
            let mut w = Complex32::new(1.0, 0.0);

            for j in 0..half {
                let t = w * buf[k + j + half];
                let u = buf[k + j];

                buf[k + j] = u + t;
                buf[k + j + half] = u - t;

                w *= wm;
            }
        }

        m <<= 1;
    }
}
