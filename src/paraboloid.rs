//! Codes related to the paraboloids. We use these to represent the
//! local behavior of a sample point that satisfies the second-order differential equation.

// todo: We may need complex Vecs.

use lin_alg::f64::{Vec2, Vec3};
use lin_alg::complex_nums::Cplx;

/// A 1D function. Assumes we've already calculated psi and psi''. `x_ref` is the point where
/// psi and psi'' were chosen. A function on a line.
// todo: Real for now, for simplicity.
// fn trace_1d(x: f64, x_ref: f64, psi: Cplx, psi_pp: Cplx) -> Cplx {
fn trace_1d(x: f64, x_ref: f64, psi: f64, psi_pp: f64) -> f64 {
    // todo: Confirm we don't need a linear term. I believe that is undefined,
    // todo due to our diffeq not specifying psi' in the relation.
    let x_diff = x - x_ref;
    let lin_term = 0.0_f64; // todo?

    let x_sq = 0.5 * x_diff.powi(2);

    psi_pp * x_sq + psi
}

/// A 2D function. Assumes we've already calculated psi and psi''. `x_ref` is the point where
/// psi and psi'' were chosen. A function on a surface.
// todo: Real for now, for simplicity.
// fn trace_2d(x: Vec2, x_ref: Vec2, psi: Cplx, psi_pp: Cplx) -> Cplx {
fn trace_2d(x: Vec2, x_ref: Vec2, psi: f64, psi_pp: f64) -> f64 {
    // todo: This is trickier. You need all combinations. Which is unbounded, I believe, given
    // todo you can balance positives and negatives.
    // let x_sq = 0.5 * (x - x_ref).powi(2);
    // let x_diff = x - x_ref;
    let lin_term = 0.0_f64; // todo?

    // todo: Experimenting
    let x_sq_x = 0.5 * (x.x - x_ref.x).powi(2);
    let x_sq_y = 0.5 * (x.y - x_ref.y).powi(2);

    (psi_pp * x_sq_x + psi) + (psi_pp * x_sq_y + psi)


}


/// A 3D function. Assumes we've already calculated psi and psi''. `x_ref` is the point where
/// psi and psi'' were chosen. A function on a volume.
// todo: Real for now, for simplicity.
// fn trace_3d(x: Vec3, x_ref: Vec3, psi: Cplx, psi_pp: Cplx) -> Cplx {
fn trace_3d(x: Vec3, x_ref: Vec3, psi: f64, psi_pp: f64) -> f64 {
    // let x_sq = 0.5 * (x - x_ref).powi(2);

    let x_diff = x - x_ref;
    let lin_term = 0.0_f64; // todo?

    // todo: Experimenting
    let x_sq_x = 0.5 * (x.x - x_ref.x).powi(2);
    let x_sq_y = 0.5 * (x.y - x_ref.y).powi(2);
    let x_sq_z = 0.5 * (x.z - x_ref.z).powi(2);

    (psi_pp * x_sq_x + psi) + (psi_pp * x_sq_y + psi) + (psi_pp * x_sq_z + psi)
}