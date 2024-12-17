#![allow(non_snake_case)]

use crate::complex_nums::Cplx;

mod complex_nums;
use lin_alg::f64::Vec3;

const ħ: f64 = 1.;
const K_C: f64 = 1.;
const M_ELEC: f64 = 1.;
pub const KE_COEFF: f64 = -(ħ * ħ) / (2. * M_ELEC);
pub const KE_COEFF_INV: f64 = 1. / KE_COEFF;

// We use this to prevent numerical anomolies and divide-by-0 errors in coulomb interactions, where
// the positions are very close to each other.
const SOFTENING_FACTOR: f64 = 0.000000000000001;

#[derive(Debug)]
struct Nucleus {
    pub mass: f64,
    pub charge: f64, // todo: Integer?
    pub posit: Vec3,
}

#[derive(Clone, Debug)]
struct ElectronTrace {
    pub posit: Vec3,
    /// todo: A/R. Direction the curve opens; amt of curvature. (By what metric?)
    /// todo: Or is this not directly necessary ?
    pub curve: (Vec3, f64),
}


/// Single-point Coulomb potential, from a single point charge.
pub(crate) fn V_coulomb(posit_charge: Vec3, posit_sample: Vec3, charge: f64) -> f64 {
    let diff = posit_sample - posit_charge;
    let r = diff.magnitude();

    K_C * charge / (r + SOFTENING_FACTOR)
}

/// Single-point Coulomb electric field, form a single point charge.
/// todo: Return the result as a vector?
pub(crate) fn E_coulomb(posit_charge: Vec3, posit_sample: Vec3, charge: f64) -> f64 {
    let diff = posit_sample - posit_charge;
    let r = diff.magnitude();

    K_C * charge / (r.powi(2) + SOFTENING_FACTOR)
}

/// Calculate the V that must be acting on a given psi, and its (known to be accurate, eg numerical
/// differentiation) derivative.
pub fn calc_V_on_psi(psi: Cplx, psi_pp: Cplx, E: f64) -> f64 {
    // psi''/psi is always real, due to being an eigenvalue of a Hermitian operator.
    KE_COEFF * (psi_pp / psi).real - E
}


/// A mirror of `calc_V_on_psi`.
/// note: This is identical to calc_V_on_psi.
pub fn calc_E_on_psi(psi: Cplx, psi_pp: Cplx, V: f64) -> f64 {
    calc_V_on_psi(psi, psi_pp, V)
}

/// Alternative API, taking advantage of analytic psi''/psi
/// psipp_div_psi is always real; Hermitian.
pub fn _calc_V_on_psi2(psi_pp_div_psi: f64, E: f64) -> f64 {
    KE_COEFF * psi_pp_div_psi - E
}

/// Alternative API, taking advantage of analytic psi''/psi
/// psipp_div_psi is always real; Hermitian.
pub fn _calc_E_on_psi2(psi_pp_div_psi: f64, V: f64) -> f64 {
    _calc_V_on_psi2(psi_pp_div_psi, V)
}


fn main() {
    let nucs = vec![
        Nucleus{ mass: 1., charge: 1., posit: Vec3::new_zero() }
    ];

    let traces = vec![
        ElectronTrace {posit: Vec3::new_zero(), curve: (Vec3::new_zero(), 0.)}
    ];



}
