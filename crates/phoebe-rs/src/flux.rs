//! Flux computation for EB surface elements.
//!
//! Computes visible flux at a given orbital phase, accounting for:
//! - Projection (cos between normal and line of sight)
//! - Limb darkening (linear law)
//! - Eclipse/occultation via proper Roche geometry (fblink)
//! - Gravity darkening (von Zeipel law for radiative, Lucy for convective)
//! - Simple reflection effect (irradiation from companion)

use lcurve_subs::Vec3;
use lcurve_roche::Star;
use crate::geometry::{StarMesh, SurfaceElement};
use crate::passband::SBTable;

/// Compute the line-of-sight unit vector for a given orbital phase and inclination.
///
/// Phase convention: phase=0.0 is primary eclipse (secondary in front of primary),
/// phase=0.5 is secondary eclipse.
pub fn line_of_sight(phase: f64, inclination_rad: f64) -> Vec3 {
    let sin_i = inclination_rad.sin();
    let cos_i = inclination_rad.cos();
    // Shift by half period to match standard convention (phase=0 = primary eclipse)
    let phi = std::f64::consts::TAU * (phase + 0.5);
    Vec3::new(
        -sin_i * phi.cos(),
        -sin_i * phi.sin(),
        cos_i,
    )
}

/// Compute flux from one star at a given orbital phase.
///
/// Uses proper Roche geometry eclipse checking via `fblink`.
///
/// * `mesh` — surface mesh of the emitting star
/// * `los` — line-of-sight unit vector toward observer
/// * `q` — mass ratio M2/M1
/// * `eclipser` — which star is doing the eclipsing
/// * `eclipser_fillout` — filling factor of the eclipsing star
/// * `t_eff` — effective temperature of the emitting star
/// * `t_eff_companion` — temperature of the companion (for reflection)
/// * `ld_coeff` — limb darkening coefficient (linear law)
/// * `gravity_darkening` — gravity darkening exponent
/// * `reflection_albedo` — reflection/irradiation albedo (0 = none, 1 = full)
pub fn star_flux(
    mesh: &StarMesh,
    los: &Vec3,
    q: f64,
    eclipser: Star,
    eclipser_fillout: f64,
    t_eff: f64,
    t_eff_companion: f64,
    ld_coeff: f64,
    gravity_darkening: f64,
    reflection_albedo: f64,
    sb_table: &SBTable,
) -> f64 {
    let mut flux = 0.0;
    let mean_grav: f64 = mesh.elements.iter().map(|e| e.grav).sum::<f64>()
        / mesh.elements.len() as f64;

    let companion_centre = match eclipser {
        Star::Primary => Vec3::new(0.0, 0.0, 0.0),
        Star::Secondary => Vec3::new(1.0, 0.0, 0.0),
    };

    let sb_companion = sb_table.eval(t_eff_companion);
    let sb_self = sb_table.eval(t_eff);
    let l_ratio = if sb_self > 0.0 { sb_companion / sb_self } else { 0.0 };

    let eclipse_acc = 0.05;

    for elem in &mesh.elements {
        let cos_gamma = elem.normal.x * los.x
            + elem.normal.y * los.y
            + elem.normal.z * los.z;

        if cos_gamma <= 0.0 {
            continue;
        }

        let eclipsed = lcurve_roche::eclipse::fblink(
            q, eclipser, 1.0, eclipser_fillout, eclipse_acc,
            los, &elem.pos,
        ).unwrap_or(false);

        if eclipsed {
            continue;
        }

        let ld = 1.0 - ld_coeff * (1.0 - cos_gamma);

        // Gravity darkening: T_local = T_eff * (g/g_mean)^beta
        let gd = if mean_grav > 0.0 {
            (elem.grav / mean_grav).powf(gravity_darkening)
        } else {
            1.0
        };
        let t_local = t_eff * gd;

        // Passband-dependent surface brightness
        let sb = sb_table.eval(t_local) * ld;

        // Reflection effect
        let reflection_boost = if reflection_albedo > 0.0 {
            let dx = companion_centre.x - elem.pos.x;
            let dy = companion_centre.y - elem.pos.y;
            let dz = companion_centre.z - elem.pos.z;
            let dist_sq = dx * dx + dy * dy + dz * dz;
            if dist_sq > 1e-10 {
                let dist = dist_sq.sqrt();
                let cos_irr = (elem.normal.x * dx + elem.normal.y * dy + elem.normal.z * dz) / dist;
                if cos_irr > 0.0 {
                    reflection_albedo * l_ratio * cos_irr * elem.area / (4.0 * std::f64::consts::PI * dist_sq)
                } else {
                    0.0
                }
            } else {
                0.0
            }
        } else {
            0.0
        };

        flux += sb * elem.area * cos_gamma + reflection_boost * sb_self;
    }

    flux
}
