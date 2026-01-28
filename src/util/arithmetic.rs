use crate::prelude::*;

pub fn signal_add(a: &[f32], b: &[f32], output: &mut [f32]) {
    let n = a.len();
    assert_eq!(n, b.len());
    assert_eq!(b.len(), output.len());

    for i in 0..a.len() {
        output[i] = a[i] + b[i];
    }
}

pub fn signal_scale(buf: &[f32], scale: f32, output: &mut [f32]) {
    for i in 0..buf.len() {
        output[i] = buf[i] * scale;
    }
}

