use num_complex::Complex32;
use crate::prelude::*;

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

pub struct FftPlan {
    pub n: usize,
    pub table: Vec<(usize, Vec<Complex32>)>,  // Vec< (段数, その段での回転因子) >
}

impl FftPlan {
    pub fn new(n: usize) -> Self {
        assert!(n.is_power_of_two());

        let mut table = Vec::new();
        let mut m = 2;

        while m <= n {
            let mut wm = Vec::with_capacity(m / 2);
            for k in 0..m / 2 {
                let theta = -2.0 * PI * k as f32 / m as f32;
                wm.push(Complex32::from_polar(1.0, theta));
            }
            table.push((m, wm));
            m <<= 1;
        }

        Self { n, table }
    }
}

pub fn fft(x: &mut [Complex32], plan: &FftPlan) {
    assert_eq!(x.len(), plan.n);

    bit_reverse_permute(x);

    for (m, wm) in &plan.table {
        let half = m / 2;

        for k in (0..plan.n).step_by(*m) {
            for j in 0..half {
                let t = wm[j] * x[k + j + half];
                let u = x[k + j];
                x[k + j] = u + t;
                x[k + j + half] = u - t;
            }
        }
    }
}


pub struct RealFftPlan {
    pub n: usize,
    pub half_plan: FftPlan,
    pub wk_table: Vec<Complex32>,
}

impl RealFftPlan {
    pub fn new(n: usize) -> Self {
        assert!(n.is_power_of_two());
        assert!(n >= 2);
        assert!(n % 2 == 0);

        let nh = n / 2;

        let mut wk_table = Vec::with_capacity(nh);
        for k in 0..nh {
            let theta = -2.0 * PI * k as f32 / n as f32;
            wk_table.push(Complex32::from_polar(1.0, theta));
        }

        Self {
            n,
            half_plan: FftPlan::new(n / 2),
            wk_table,
        }
    }
}

pub fn fft_real(
    input: &[f32],  // 長さ n
    output: &mut [Complex32], // 長さ n/2 + 1
    buf: &mut [Complex32],    // 長さ n/2
    real_fft_plan: &RealFftPlan,      // n/2 点 FFT 用
) {
    let n = real_fft_plan.n;
    let nh = real_fft_plan.half_plan.n;

    assert_eq!(input.len(), n);
    assert_eq!(output.len(), nh + 1);
    assert_eq!(buf.len(), nh);

    // 実信号を複素数にパック
    for k in 0..nh {
        buf[k] = Complex32::new(input[2 * k], input[2 * k + 1]);
    }

    // N/2 点 FFT
    fft(buf, &real_fft_plan.half_plan);

    // DC / Nyquist
    let a0 = buf[0];
    output[0] = Complex32::new(a0.re + a0.im, 0.0);
    output[nh] = Complex32::new(a0.re - a0.im, 0.0);

    // 正の周波数のみ生成
    for k in 1..nh {
        let a = buf[k];
        let b = buf[nh - k].conj();

        let t1 = a + b;
        let t2 = (a - b) * Complex32::new(0.0, -1.0);

        output[k] = (t1 + real_fft_plan.wk_table[k] * t2) * 0.5;
    }
}



/* sample

* complex fft *
let n = 8;
let plan = FftPlan::new(n);

let mut x = vec![Complex32::ZERO; n];
for i in 0..n {
    let theta = 2.0 * PI * i as f32 / n as f32;
    x[i] = Complex32::from_polar(1.0, theta);
}

println!("complex input:");
for v in &x {
    println!("{:?}", v);
}

fft(&mut x, &plan);

println!("\ncomplex FFT output:");
for (k, v) in x.iter().enumerate() {
    println!("k = {}: {:?}", k, v);
}

* real fft *
let n = 8;
    let real_plan = RealFftPlan::new(n);

    let mut input = vec![0.0f32; n];
    for i in 0..n {
        input[i] = (2.0 * PI * i as f32 / n as f32).sin();
    }

    let mut output = vec![Complex32::ZERO; n / 2 + 1];
    let mut buf = vec![Complex32::ZERO; n / 2];

    println!("\nreal input:");
    for v in &input {
        println!("{:.6}", v);
    }

    fft_real(&input, &mut output, &mut buf, &real_plan);

    println!("\nreal FFT (half spectrum):");
    for (k, v) in output.iter().enumerate() {
        println!("k = {}: {:?}", k, v);
    }
*/





/*
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

*/