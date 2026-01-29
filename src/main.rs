use num_complex::Complex32;
use sigrus::prelude::*;

fn main() {
    const FS: f32 = 16.;
    const SEC: f32 = 1.;

    let mut s = [0.0f32; (FS * SEC) as usize];
    osc::sin(&mut s, 1., 1., 0., FS);

    let mut c = [0.0f32; (FS * SEC) as usize];
    osc::sin(&mut c, 2.0, 2., -HALF_PI, FS);

    let mut sig = [Complex32::ZERO; (FS * SEC) as usize];
    for i in 0..sig.len() {
        sig[i] = Complex32::new(s[i], c[i]);
    }
    println!("orig sig\n{:?}\n", &sig);
}
