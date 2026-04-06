//! Roche geometry and surface mesh generation for EB stars.
//!
//! Supports both detached (separate Roche lobes) and overcontact
//! (common envelope) configurations.

use lcurve_subs::{Vec3, PI, TWOPI};
use lcurve_roche::Star;
use crate::PhoebeError;

/// A surface element on a star.
#[derive(Debug, Clone)]
pub struct SurfaceElement {
    /// Position in the corotating frame (units of separation)
    pub pos: Vec3,
    /// Outward surface normal
    pub normal: Vec3,
    /// Projected area (in units of separation^2)
    pub area: f64,
    /// Local effective gravity (for gravity darkening)
    pub grav: f64,
    /// Which star this element belongs to (for temperature assignment)
    pub owner: Star,
}

/// Surface mesh for one or both stars.
#[derive(Debug)]
pub struct StarMesh {
    pub elements: Vec<SurfaceElement>,
    pub star: Star,
}

/// Build surface mesh for a Roche-distorted star (detached case).
pub fn build_mesh(
    q: f64,
    star: Star,
    fillout: f64,
    n_lat: usize,
) -> Result<StarMesh, PhoebeError> {
    let spin = 1.0;
    let (rref, pot_ref) = lcurve_roche::surface::ref_sphere(q, star, spin, fillout)?;

    let mut elements = Vec::new();
    let d_theta = PI / n_lat as f64;

    let x_off = match star {
        Star::Primary => 0.0,
        Star::Secondary => 1.0,
    };

    for i in 0..n_lat {
        let theta = PI * (i as f64 + 0.5) / n_lat as f64;
        let sin_t = theta.sin();
        let cos_t = theta.cos();
        let n_lng = (16.0f64).max(TWOPI * sin_t / d_theta) as usize;
        let d_phi = TWOPI / n_lng as f64;

        for j in 0..n_lng {
            let phi = TWOPI * (j as f64 + 0.5) / n_lng as f64;
            let dir = Vec3::new(
                sin_t * phi.cos(),
                sin_t * phi.sin(),
                cos_t,
            );

            let mut r = rref;
            for _ in 0..20 {
                let pos = Vec3::new(x_off + r * dir.x, r * dir.y, r * dir.z);
                let pot = match star {
                    Star::Primary => lcurve_roche::potential::rpot1(q, spin, &pos)?,
                    Star::Secondary => lcurve_roche::potential::rpot2(q, spin, &pos)?,
                };
                let grad = match star {
                    Star::Primary => lcurve_roche::potential::drpot1(q, spin, &pos)?,
                    Star::Secondary => lcurve_roche::potential::drpot2(q, spin, &pos)?,
                };
                let dpot_dr = grad.x * dir.x + grad.y * dir.y + grad.z * dir.z;
                if dpot_dr.abs() < 1e-15 { break; }
                let dr = (pot - pot_ref) / dpot_dr;
                r -= dr;
                if dr.abs() < 1e-10 * r.abs().max(1e-10) { break; }
            }

            let pos = Vec3::new(x_off + r * dir.x, r * dir.y, r * dir.z);
            let grad = match star {
                Star::Primary => lcurve_roche::potential::drpot1(q, spin, &pos)?,
                Star::Secondary => lcurve_roche::potential::drpot2(q, spin, &pos)?,
            };
            let grad_mag = (grad.x * grad.x + grad.y * grad.y + grad.z * grad.z).sqrt();
            let normal = if grad_mag > 0.0 {
                Vec3::new(grad.x / grad_mag, grad.y / grad_mag, grad.z / grad_mag)
            } else {
                dir
            };

            // Correct area for surface tilt: dA = r² sinθ dθ dφ / |n·r̂|
            // where n is the surface normal and r̂ is the radial direction.
            // For a sphere n=r̂ and |n·r̂|=1. For Roche distortion |n·r̂|<1
            // so the true area is larger than the spherical estimate.
            let cos_tilt = (normal.x * dir.x + normal.y * dir.y + normal.z * dir.z).abs();
            let area = if cos_tilt > 0.01 {
                r * r * sin_t * d_theta * d_phi / cos_tilt
            } else {
                r * r * sin_t * d_theta * d_phi / 0.01 // cap at 100x
            };

            elements.push(SurfaceElement {
                pos, normal, area, grav: grad_mag, owner: star,
            });
        }
    }

    Ok(StarMesh { elements, star })
}


/// Roche potential adapter for marching triangulation.
struct RocheSurface {
    q: f64,
}

impl crate::marching::PotentialSurface for RocheSurface {
    fn potential(&self, p: &Vec3) -> f64 {
        lcurve_roche::potential::rpot(self.q, p).unwrap_or(0.0)
    }
    fn gradient(&self, p: &Vec3) -> Vec3 {
        lcurve_roche::potential::drpot(self.q, p).unwrap_or(Vec3::new(0.0, 0.0, 0.0))
    }
}

/// Build overcontact mesh using marching triangulation from BOTH poles.
///
/// Seeds the marching from each star's pole, then merges.
/// This ensures full coverage including the neck region.
pub fn build_overcontact_mesh_marching(
    q: f64,
    fillout: f64,
    max_triangles: usize,
) -> Result<(StarMesh, StarMesh), PhoebeError> {
    let spin = 1.0;

    let xl1 = lcurve_roche::lagrange::xl11(q, spin)?;
    let l1_pos = Vec3::new(xl1, 0.0, 0.0);
    let pot_l1 = lcurve_roche::potential::rpot(q, &l1_pos)?;

    let xl2_raw = lcurve_roche::lagrange::xl12(q, spin)?;
    let l2_pos = Vec3::new(1.0 + (1.0 - xl2_raw), 0.0, 0.0);
    let pot_l2 = lcurve_roche::potential::rpot(q, &l2_pos)?;

    let pot_surface = pot_l1 + fillout * (pot_l2 - pot_l1);
    let surface = RocheSurface { q };

    let half_tri = max_triangles / 2;
    let r_eff = xl1 * 0.7;
    let delta = (4.0 * PI * r_eff * r_eff / max_triangles as f64).sqrt();

    // Mesh from primary pole
    let seed1 = Vec3::new(0.0, 0.0, xl1 * 0.5);
    let mesh1 = crate::marching::marching_triangulate(
        &surface, pot_surface, &seed1, delta, half_tri)?;

    // Mesh from secondary pole
    let seed2 = Vec3::new(1.0, 0.0, (1.0 - xl1) * 0.5);
    let mesh2 = crate::marching::marching_triangulate(
        &surface, pot_surface, &seed2, delta, half_tri)?;

    // Convert triangle meshes to surface elements, split by L1
    let mut elems1 = Vec::new();
    let mut elems2 = Vec::new();

    for mesh in [&mesh1, &mesh2] {
        let props = mesh.triangle_properties();
        for (i, tri) in mesh.triangles.iter().enumerate() {
            let (area, normal, grav) = &props[i];
            let v0 = &mesh.vertices[tri.v0];
            let v1 = &mesh.vertices[tri.v1];
            let v2 = &mesh.vertices[tri.v2];
            let centroid = Vec3::new(
                (v0.pos.x + v1.pos.x + v2.pos.x) / 3.0,
                (v0.pos.y + v1.pos.y + v2.pos.y) / 3.0,
                (v0.pos.z + v1.pos.z + v2.pos.z) / 3.0,
            );

            let elem = SurfaceElement {
                pos: centroid,
                normal: *normal,
                area: *area,
                grav: *grav,
                owner: if centroid.x < xl1 { Star::Primary } else { Star::Secondary },
            };

            if centroid.x < xl1 {
                elems1.push(elem);
            } else {
                elems2.push(elem);
            }
        }
    }

    Ok((
        StarMesh { elements: elems1, star: Star::Primary },
        StarMesh { elements: elems2, star: Star::Secondary },
    ))
}

/// Build surface mesh for an overcontact (common envelope) binary.
///
/// Each star is meshed from its own centre using the shared overcontact
/// potential. The surface potential interpolates between L1 (fillout=0)
/// and L2 (fillout=1). Elements beyond L1 toward the companion capture
/// the neck region.
pub fn build_overcontact_mesh(
    q: f64,
    fillout: f64,
    n_lat: usize,
) -> Result<(StarMesh, StarMesh), PhoebeError> {
    let spin = 1.0;

    // L1 point and critical potential
    let xl1 = lcurve_roche::lagrange::xl11(q, spin)?;
    let l1_pos = Vec3::new(xl1, 0.0, 0.0);
    let pot_l1 = lcurve_roche::potential::rpot(q, &l1_pos)?;

    // L2 point potential (outer critical surface)
    let xl2_raw = lcurve_roche::lagrange::xl12(q, spin)?;
    // xl12 returns distance from secondary, so L2 position is 1 + (1-xl2)
    let l2_pos = Vec3::new(1.0 + (1.0 - xl2_raw), 0.0, 0.0);
    let pot_l2 = lcurve_roche::potential::rpot(q, &l2_pos)?;

    // Surface potential: fillout interpolates L1→L2
    let pot_surface = pot_l1 + fillout * (pot_l2 - pot_l1);

    // Mesh each star from its own centre using the synchronous Roche potential
    let mesh1 = build_overcontact_half(q, Star::Primary, pot_surface, xl1, n_lat)?;
    let mesh2 = build_overcontact_half(q, Star::Secondary, pot_surface, xl1, n_lat)?;

    Ok((mesh1, mesh2))
}

/// Build one half of the overcontact surface from a star's centre.
fn build_overcontact_half(
    q: f64,
    star: Star,
    pot_surface: f64,
    xl1: f64,
    n_lat: usize,
) -> Result<StarMesh, PhoebeError> {
    let x_off = match star {
        Star::Primary => 0.0,
        Star::Secondary => 1.0,
    };

    // Initial radius guess: half the distance to L1
    let r_guess = match star {
        Star::Primary => xl1 * 0.5,
        Star::Secondary => (1.0 - xl1) * 0.5,
    };

    let mut elements = Vec::new();
    let d_theta = PI / n_lat as f64;

    for i in 0..n_lat {
        let theta = PI * (i as f64 + 0.5) / n_lat as f64;
        let sin_t = theta.sin();
        let cos_t = theta.cos();
        let n_lng = (16.0f64).max(TWOPI * sin_t / d_theta) as usize;
        let d_phi = TWOPI / n_lng as f64;

        for j in 0..n_lng {
            let phi = TWOPI * (j as f64 + 0.5) / n_lng as f64;
            let dir = Vec3::new(
                sin_t * phi.cos(),
                sin_t * phi.sin(),
                cos_t,
            );

            // Newton-Raphson: find r where rpot(x_off + r*dir) = pot_surface
            let mut r = r_guess;
            let mut converged = false;

            for _ in 0..30 {
                let pos = Vec3::new(x_off + r * dir.x, r * dir.y, r * dir.z);
                let pot = lcurve_roche::potential::rpot(q, &pos)?;
                let grad = lcurve_roche::potential::drpot(q, &pos)?;

                let dpot_dr = grad.x * dir.x + grad.y * dir.y + grad.z * dir.z;
                if dpot_dr.abs() < 1e-15 { break; }
                let dr = (pot - pot_surface) / dpot_dr;
                r -= dr;
                r = r.max(0.005).min(1.5);
                if dr.abs() < 1e-9 * r {
                    converged = true;
                    break;
                }
            }

            if !converged { continue; }

            let pos = Vec3::new(x_off + r * dir.x, r * dir.y, r * dir.z);

            // Skip elements that belong to the other star's side
            // (past L1 for primary, before L1 for secondary)
            match star {
                Star::Primary => { if pos.x > xl1 { continue; } },
                Star::Secondary => { if pos.x < xl1 { continue; } },
            }

            let grad = lcurve_roche::potential::drpot(q, &pos)?;
            let grad_mag = (grad.x * grad.x + grad.y * grad.y + grad.z * grad.z).sqrt();
            let normal = if grad_mag > 0.0 {
                Vec3::new(grad.x / grad_mag, grad.y / grad_mag, grad.z / grad_mag)
            } else {
                dir
            };

            // Correct area for surface tilt
            let cos_tilt = (normal.x * dir.x + normal.y * dir.y + normal.z * dir.z).abs();
            let area = if cos_tilt > 0.01 {
                r * r * sin_t * d_theta * d_phi / cos_tilt
            } else {
                r * r * sin_t * d_theta * d_phi / 0.01
            };

            elements.push(SurfaceElement {
                pos, normal, area, grav: grad_mag, owner: star,
            });
        }
    }

    Ok(StarMesh { elements, star })
}
