use crate::prelude::*;

pub fn add(lhs: &[f32], rhs: &[f32], output: &mut [f32]) {
    assert_eq!(lhs.len(), rhs.len());
    assert_eq!(rhs.len(), rhs.len());

    for i in 0..lhs.len() {
        output[i] = lhs[i] + rhs[i];
    }
}

