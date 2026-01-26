use rand::Rng;
use crate::prelude::*;

pub trait Oscillator {
    fn oscillate();
}

pub fn sine(
    out: &mut [f32],
    amp: f32,
    freq: f32,
    phase0: f32,
    sample_rate: f32,
) -> f32 {
    let mut p = phase0;
    let dp = TWO_PI * freq / sample_rate;

    for x in out.iter_mut() {
        *x = amp * p.sin();
        p += dp;
        if p >= TWO_PI {
            p -= TWO_PI;
        }
    }

    p
}

pub fn saw(
    out: &mut [f32],
    amp: f32,
    freq: f32,
    phase0: f32,
    sample_rate: f32,
) -> f32 {
    let mut p = phase0;
    let dp = TWO_PI * freq / sample_rate;

    for x in out.iter_mut() {
        let u = p / TWO_PI; // [0,1)
        *x = amp * (2.0 * u - 1.0);

        p += dp;
        if p >= TWO_PI {
            p -= TWO_PI;
        }
    }
    p
}

pub fn triangle(
    out: &mut [f32],
    amp: f32,
    freq: f32,
    phase0: f32,
    sample_rate: f32,
) -> f32 {
    let mut p = phase0;
    let dp = TWO_PI * freq / sample_rate;

    for x in out.iter_mut() {
        let u = p / TWO_PI;
        *x = amp * (1.0 - 4.0 * (u - 0.5).abs());

        p += dp;
        if p >= TWO_PI {
            p -= TWO_PI;
        }
    }
    p
}

pub fn square(
    out: &mut [f32],
    amp: f32,
    freq: f32,
    phase0: f32,
    sample_rate: f32,
) -> f32 {
    let mut p = phase0;
    let dp = TWO_PI * freq / sample_rate;

    for x in out.iter_mut() {
        let u = p / TWO_PI;
        *x = if u < 0.5 { amp } else { -amp };

        p += dp;
        if p >= TWO_PI {
            p -= TWO_PI;
        }
    }
    p
}

pub fn white_noise(output: &mut [f32]) {
    let mut rng = rand::rng();

    for x in output.iter_mut() {
        *x = rng.random_range(-1.0..1.0);
    }
}

pub fn pink_noise(output: &mut [f32], rows: usize) {
    let n = output.len();
    let mut rng = rand::rng();

    // 作業用バッファ
    let mut acc = vec![0.0f32; n];

    for _ in 0..rows {
        let mut sum = 0.0f32;
        for i in 0..n {
            let w: f32 = rng.random(); // 一様分布
            sum += w;
            acc[i] += sum;
        }
    }

    // 正規化（max abs）
    let max = acc
        .iter()
        .map(|x| x.abs())
        .fold(0.0f32, f32::max);

    if max > 0.0 {
        for (o, x) in output.iter_mut().zip(acc.iter()) {
            *o = *x / max;
        }
    }
}

