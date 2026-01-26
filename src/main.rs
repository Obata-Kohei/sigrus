use num_complex::Complex32;
use sigrus::prelude::*;

fn main() {
    const FS: f32 = 16.;
    const SEC: f32 = 1.;
    let eps: f32 = 1e-5;


    let mut s = [0.0f32; (FS * SEC) as usize];
    sine(&mut s, 1., 1., 0., FS);

    let mut c = [0.0f32; (FS * SEC) as usize];
    sine(&mut c, 2.0, 2., -HALF_PI, FS);

    let mut sig = [Complex32::ZERO; (FS * SEC) as usize];
    for i in 0..sig.len() {
        sig[i] = Complex32::new(s[i], c[i]);
    }
    println!("orig sig\n{:?}\n", &sig);


    let mut out = [Complex32::ZERO; (FS * SEC) as usize];
    dft(&sig, &mut out);
    println!("DFT\n{:?}\n", &out);


    let mut out2 = [Complex32::ZERO; (FS * SEC) as usize];
    idft(&out, &mut out2);
    println!("iDFT == sig: {}\n", approx_eq_complex32_slice(&sig, &out2, eps));


    let mut out3 = [Complex32::ZERO; (FS * SEC) as usize];
    fft_cooley_tukey(&sig, &mut out3);
    println!("fft_c&t == dft: {}\n", approx_eq_complex32_slice(&out, &out3, eps));

}
