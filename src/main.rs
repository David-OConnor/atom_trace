#![allow(non_snake_case)]

use std::ops::Add;
use lin_alg::complex_nums::Cplx;

mod integrate;
mod paraboloid;

use lin_alg::f64::Vec3;
use lin_alg::complex_nums::IM;
use crate::integrate::integrate_rk4;

use barnes_hut::BodyModel;

const ħ: f64 = 1.;
const K_C: f64 = 1.;
const M_ELEC: f64 = 1.;
const Q_ELEC: f64 = -1.;
pub const KE_COEFF: f64 = -(ħ * ħ) / (2. * M_ELEC);
pub const KE_COEFF_INV: f64 = 1. / KE_COEFF;

// We use this to prevent numerical anomolies and divide-by-0 errors in coulomb interactions, where
// the positions are very close to each other.
const SOFTENING_FACTOR: f64 = 1e-6;

#[derive(Debug)]
struct Nucleus {
    pub mass: f64,
    pub charge: f64, // todo: Integer?
    pub posit: Vec3,
}

#[derive(Clone, Debug)]
struct ElectronTrace {
    pub posit: Vec3,
    pub vel: Vec3,
    pub accel: Vec3,
    /// todo: A/R. Direction the curve opens; amt of curvature. (By what metric?)
    /// todo: Or is this not directly necessary ?
    pub curve: (Vec3, f64),
}

impl BodyModel for ElectronTrace {
    fn posit(&self) -> Vec3 {
        self.posit
    }

    fn mass(&self) -> f64 {
        Q_ELEC
    }
}

#[derive(Clone, Default, Debug)]
/// A set of derivatives at a single point
/// Organization: index, da
pub struct DerivativesSingle {
    pub dx: Cplx,
    pub dy: Cplx,
    pub dz: Cplx,
    pub d2x: Cplx,
    pub d2y: Cplx,
    pub d2z: Cplx,
    pub d2_sum: Cplx,
}

impl Add<&Self> for DerivativesSingle {
    type Output = Self;

    fn add(self, rhs: &Self) -> Self::Output {
        let mut result = self;

        result.dx = result.dx + rhs.dx;
        result.dy = result.dy + rhs.dy;
        result.dz = result.dz + rhs.dz;
        result.d2x = result.d2x + rhs.d2x;
        result.d2y = result.d2y + rhs.d2y;
        result.d2z = result.d2z + rhs.d2z;
        result.d2_sum = result.d2_sum + rhs.d2_sum;

        result
    }
}

/// Calculate L^2, given derivatives. Used in one of the two momentum eigenfunctions. See Onenote:
/// Exploring the WF, part 9
/// todo: This is coming out incorrectly. Note that the `l_z` function above appears to work correctly.
pub(crate) fn calc_L_sq(posit: Vec3, d: &DerivativesSingle) -> Cplx {
    let x = posit.x;
    let y = posit.y;
    let z = posit.z;

    // todo: This is fine too; just a different representation.
    // let part0 =
    //     (d.d2y + d.d2z) * x.powi(2) + (d.d2x + d.d2z) * y.powi(2) + (d.d2x + d.d2y) * z.powi(2);

    // todo: This should be equivalent, and it seems to be so.
    let x_sq = x.powi(2);
    let y_sq = y.powi(2);
    let z_sq = z.powi(2);

    let part0 = d.d2z * (x_sq + y_sq) + d.d2x * (y_sq + z_sq) + d.d2y * (z_sq + x_sq);

    // // todo: Experimenting...
    // let l_sq_sq = part0.abs_sq();
    // return Cplx::from_real(l_sq_sq);

    // todo: This is wrong (Multiplying dpsi/dx * dpsi/dy. But shouldn't this come out to 0, so it doesn't matter??
    // let part1 = (d.dx * d.dy * x * y + d.dx * d.dz * x * z + d.dy * d.dz * y * z) * 2.;

    // todo: Guess
    let d_dydz = d.dy * d.dz;
    let d_dzdx = d.dz * d.dx;
    let d_dxdy = d.dx * d.dy;

    // Todo: My best guess: We need the cross terms. Ref ChatGpt. YOu need d^2/(dx dz) and d^2/(dzdx)
    let part1 = (d_dydz * y * z + d_dzdx * z * x + d_dxdy * x * y) * 2.;
    // println!("Part 1:  {:?}", part1); // todo: part1 0 from cancelling cross terms? Appears to be.

    // todo: Maybe you need to take the derivatives in series, instead of multiplying? Yea.
    // todo: dx * dx isn't (dpsi/dx)^2... it's the second derivative.

    // TS note: Part 1 seems to be 0 for several test cases. (Cross-deriv terms cancelling, likely)

    -part0
    // -part0 + part1
}

/// We choose z by convention; we need one of any dimension.
/// Uses the momentum eigenfunction p_a = -i ħ d/da. L_z = -iħ (x d/dy - y d/dx)
pub fn calc_L_z(posit: Vec3, d: &DerivativesSingle) -> Cplx {
    -IM * (d.dy * posit.x - d.dx * posit.y)
}

/// See `calc_L_z`.
pub fn calc_L_x(posit: Vec3, d: &DerivativesSingle) -> Cplx {
    -IM * (d.dz * posit.y - d.dy * posit.z)
}

/// See `calc_L_z`.
pub fn calc_L_y(posit: Vec3, d: &DerivativesSingle) -> Cplx {
    -IM * (d.dx * posit.z - d.dz * posit.x)
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

/// Calcualte psi'', calculated from psi, and E. Note that the V term used must include both
/// electron-electron interactions, and electron-proton interactions.
/// V here is a potential field from the nucleus, so multiply it by the electron's
/// charge to find potential energy
/// At a given i, j, k.
///
/// This solves, analytically, the eigenvalue equation for the Hamiltonian operator.
///
/// Hψ = Eψ. -ħ^2/2m * ψ'' + Vψ = Eψ. ψ'' = [(E - V) / (-ħ^2/2m)] ψ
pub fn calc_ψ_pp(ψ: Cplx, V: f64, E: f64) -> Cplx {
    // Note that V input is potential field; we get potential energy by multiplying it
    // by the charge being acted on (the electron)
    // todo: QC -V sign.
    ψ * (E - V) * KE_COEFF_INV
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

    let mut traces = vec![
        ElectronTrace {posit: Vec3::new_zero(), curve: (Vec3::new_zero(), 0.)}
    ];
    let dt = 0.001;

    for i in 0..10_000 {
        integrate_rk4(&mut traces, dt)
    }


}
