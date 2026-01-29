use crate::constants::*;

pub fn next_phase_of(phase: f32, freq: f32, sample_rate: f32) -> f32 {
    let mut next_phase = phase + (TWO_PI * freq / sample_rate);
    if next_phase >= TWO_PI {
        next_phase -= TWO_PI;
    }

    next_phase
}

pub fn sin(amp: f32, phase: f32) -> f32 {
    amp * phase.sin()
}

pub fn saw(amp: f32, phase: f32) -> f32 {
    amp * (2.0 * (phase / TWO_PI) - 1.0)
}

pub fn tri(amp: f32, phase: f32) -> f32 {
    amp * (1.0 - 4.0 * ((phase / TWO_PI) - 0.5).abs())
}

pub fn sqr(amp: f32, phase: f32) -> f32 {
    if phase / TWO_PI < 0.5 { amp } else { -amp }
}
