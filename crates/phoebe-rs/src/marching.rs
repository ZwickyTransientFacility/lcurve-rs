//! Marching triangulation for Roche equipotential surfaces.
//!
//! Based on: E. Hartmann, "A marching method for the triangulation of surfaces",
//! The Visual Computer (1998) 14: 95-108.
//!
//! Adapted from PHOEBE2's triang_marching.h (Martin Horvat, 2016).

use lcurve_subs::Vec3;
use crate::PhoebeError;

/// A triangle mesh vertex.
#[derive(Debug, Clone)]
pub struct Vertex {
    /// Position on the surface
    pub pos: Vec3,
    /// Surface normal (unit vector)
    pub normal: Vec3,
    /// Tangent basis vector 1
    pub e1: Vec3,
    /// Tangent basis vector 2
    pub e2: Vec3,
    /// Gradient magnitude (local gravity)
    pub grav: f64,
}

/// A triangle defined by vertex indices.
#[derive(Debug, Clone, Copy)]
pub struct Triangle {
    pub v0: usize,
    pub v1: usize,
    pub v2: usize,
}

/// Result of marching triangulation.
#[derive(Debug)]
pub struct TriMesh {
    pub vertices: Vec<Vertex>,
    pub triangles: Vec<Triangle>,
}

impl TriMesh {
    /// Compute the area-weighted properties for each triangle.
    /// Returns (area, centroid_normal, centroid_grav) per triangle.
    pub fn triangle_properties(&self) -> Vec<(f64, Vec3, f64)> {
        self.triangles.iter().map(|tri| {
            let v0 = &self.vertices[tri.v0];
            let v1 = &self.vertices[tri.v1];
            let v2 = &self.vertices[tri.v2];

            // Area via cross product
            let e1 = Vec3::new(v1.pos.x - v0.pos.x, v1.pos.y - v0.pos.y, v1.pos.z - v0.pos.z);
            let e2 = Vec3::new(v2.pos.x - v0.pos.x, v2.pos.y - v0.pos.y, v2.pos.z - v0.pos.z);
            let cross = Vec3::new(
                e1.y * e2.z - e1.z * e2.y,
                e1.z * e2.x - e1.x * e2.z,
                e1.x * e2.y - e1.y * e2.x,
            );
            let area = 0.5 * (cross.x * cross.x + cross.y * cross.y + cross.z * cross.z).sqrt();

            // Average normal and gravity
            let normal = Vec3::new(
                (v0.normal.x + v1.normal.x + v2.normal.x) / 3.0,
                (v0.normal.y + v1.normal.y + v2.normal.y) / 3.0,
                (v0.normal.z + v1.normal.z + v2.normal.z) / 3.0,
            );
            let grav = (v0.grav + v1.grav + v2.grav) / 3.0;

            (area, normal, grav)
        }).collect()
    }
}

/// Potential function trait for marching triangulation.
pub trait PotentialSurface {
    /// Evaluate potential at position p.
    fn potential(&self, p: &Vec3) -> f64;
    /// Evaluate gradient of potential at position p.
    fn gradient(&self, p: &Vec3) -> Vec3;
}

/// Project a point onto the equipotential surface via Newton-Raphson.
fn project_onto_surface<P: PotentialSurface>(
    surface: &P,
    target_pot: f64,
    guess: &Vec3,
    hint_normal: &Vec3,
    max_iter: usize,
) -> Option<Vertex> {
    let mut pos = *guess;

    for _ in 0..max_iter {
        let pot = surface.potential(&pos);
        let grad = surface.gradient(&pos);
        let grad_sq = grad.x * grad.x + grad.y * grad.y + grad.z * grad.z;

        if grad_sq < 1e-30 {
            return None;
        }

        let dpot = pot - target_pot;
        if dpot.abs() < 1e-12 * target_pot.abs().max(1.0) {
            // Converged — build vertex
            let grad_mag = grad_sq.sqrt();
            let normal = Vec3::new(grad.x / grad_mag, grad.y / grad_mag, grad.z / grad_mag);

            // Build tangent basis via Gram-Schmidt
            // Pick a vector not parallel to normal
            let seed = if normal.x.abs() < 0.9 {
                Vec3::new(1.0, 0.0, 0.0)
            } else {
                Vec3::new(0.0, 1.0, 0.0)
            };
            let e1_raw = Vec3::new(
                seed.x - (seed.x * normal.x + seed.y * normal.y + seed.z * normal.z) * normal.x,
                seed.y - (seed.x * normal.x + seed.y * normal.y + seed.z * normal.z) * normal.y,
                seed.z - (seed.x * normal.x + seed.y * normal.y + seed.z * normal.z) * normal.z,
            );
            let e1_mag = (e1_raw.x * e1_raw.x + e1_raw.y * e1_raw.y + e1_raw.z * e1_raw.z).sqrt();
            let e1 = Vec3::new(e1_raw.x / e1_mag, e1_raw.y / e1_mag, e1_raw.z / e1_mag);
            let e2 = Vec3::new(
                normal.y * e1.z - normal.z * e1.y,
                normal.z * e1.x - normal.x * e1.z,
                normal.x * e1.y - normal.y * e1.x,
            );

            return Some(Vertex { pos, normal, e1, e2, grav: grad_mag });
        }

        // Newton step along gradient direction
        pos.x -= dpot * grad.x / grad_sq;
        pos.y -= dpot * grad.y / grad_sq;
        pos.z -= dpot * grad.z / grad_sq;
    }

    None
}

/// Perform marching triangulation on an equipotential surface.
///
/// Starts from `seed_pos` (must be close to the surface) and grows
/// triangles until the surface is covered or `max_triangles` is reached.
///
/// `delta` controls the edge length of triangles.
pub fn marching_triangulate<P: PotentialSurface>(
    surface: &P,
    target_pot: f64,
    seed_pos: &Vec3,
    delta: f64,
    max_triangles: usize,
) -> Result<TriMesh, PhoebeError> {
    let max_iter = 50;

    // Project seed onto surface
    let seed_normal = Vec3::new(0.0, 0.0, 1.0);
    let v0 = project_onto_surface(surface, target_pot, seed_pos, &seed_normal, max_iter)
        .ok_or(PhoebeError::InvalidParams("seed projection failed".to_string()))?;

    let mut vertices: Vec<Vertex> = vec![v0.clone()];
    let mut triangles: Vec<Triangle> = Vec::new();

    // Create initial hexagon around seed
    let mut front: Vec<usize> = Vec::new(); // indices of front vertices (circular)
    let pi3 = std::f64::consts::FRAC_PI_3;

    for k in 0..6 {
        let angle = k as f64 * pi3;
        let (sa, ca) = angle.sin_cos();

        // Point in tangent plane at distance delta
        let guess = Vec3::new(
            v0.pos.x + delta * (ca * v0.e1.x + sa * v0.e2.x),
            v0.pos.y + delta * (ca * v0.e1.y + sa * v0.e2.y),
            v0.pos.z + delta * (ca * v0.e1.z + sa * v0.e2.z),
        );

        if let Some(vk) = project_onto_surface(surface, target_pot, &guess, &v0.normal, max_iter) {
            let idx = vertices.len();
            vertices.push(vk);
            front.push(idx);
        } else {
            return Err(PhoebeError::InvalidParams("initial hexagon projection failed".to_string()));
        }
    }

    // Create initial 6 triangles
    for k in 0..6 {
        let k_next = (k + 1) % 6;
        triangles.push(Triangle {
            v0: 0,
            v1: front[k],
            v2: front[k_next],
        });
    }

    // March outward: expand the front polygon
    let mut iteration = 0;
    let max_iterations = max_triangles * 3;

    while !front.is_empty() && triangles.len() < max_triangles && iteration < max_iterations {
        iteration += 1;

        let n_front = front.len();
        if n_front < 3 {
            break;
        }

        // Find the front edge with the largest "opening" angle
        // (simplified: just process the first vertex)
        let fi = iteration % n_front;
        let vi = front[fi];
        let vi_prev = front[(fi + n_front - 1) % n_front];
        let vi_next = front[(fi + 1) % n_front];

        let v_curr = vertices[vi].clone();
        let v_prev = vertices[vi_prev].clone();
        let v_next = vertices[vi_next].clone();

        // Compute angle at this front vertex
        let d_prev = Vec3::new(
            v_prev.pos.x - v_curr.pos.x,
            v_prev.pos.y - v_curr.pos.y,
            v_prev.pos.z - v_curr.pos.z,
        );
        let d_next = Vec3::new(
            v_next.pos.x - v_curr.pos.x,
            v_next.pos.y - v_curr.pos.y,
            v_next.pos.z - v_curr.pos.z,
        );

        let dot = d_prev.x * d_next.x + d_prev.y * d_next.y + d_prev.z * d_next.z;
        let len_prev = (d_prev.x * d_prev.x + d_prev.y * d_prev.y + d_prev.z * d_prev.z).sqrt();
        let len_next = (d_next.x * d_next.x + d_next.y * d_next.y + d_next.z * d_next.z).sqrt();

        if len_prev < 1e-15 || len_next < 1e-15 {
            front.remove(fi);
            continue;
        }

        let cos_angle = (dot / (len_prev * len_next)).clamp(-1.0, 1.0);
        let angle = cos_angle.acos();

        if angle < std::f64::consts::FRAC_PI_3 * 1.5 {
            // Small angle: close with one triangle (ear clipping)
            triangles.push(Triangle { v0: vi_prev, v1: vi, v2: vi_next });
            front.remove(fi);
        } else if angle < std::f64::consts::FRAC_PI_3 * 2.5 {
            // Medium angle: add one new vertex
            let mid_dir = Vec3::new(
                (d_prev.x / len_prev + d_next.x / len_next) * 0.5,
                (d_prev.y / len_prev + d_next.y / len_next) * 0.5,
                (d_prev.z / len_prev + d_next.z / len_next) * 0.5,
            );
            let mid_len = (mid_dir.x * mid_dir.x + mid_dir.y * mid_dir.y + mid_dir.z * mid_dir.z).sqrt();

            let guess = if mid_len > 1e-15 {
                Vec3::new(
                    v_curr.pos.x + delta * mid_dir.x / mid_len,
                    v_curr.pos.y + delta * mid_dir.y / mid_len,
                    v_curr.pos.z + delta * mid_dir.z / mid_len,
                )
            } else {
                // Fallback: use normal direction
                Vec3::new(
                    v_curr.pos.x + delta * v_curr.e1.x,
                    v_curr.pos.y + delta * v_curr.e1.y,
                    v_curr.pos.z + delta * v_curr.e1.z,
                )
            };

            if let Some(v_new) = project_onto_surface(surface, target_pot, &guess, &v_curr.normal, max_iter) {
                let new_idx = vertices.len();
                vertices.push(v_new);
                triangles.push(Triangle { v0: vi_prev, v1: vi, v2: new_idx });
                triangles.push(Triangle { v0: vi, v1: vi_next, v2: new_idx });
                front[fi] = new_idx;
            } else {
                // Projection failed, just close
                triangles.push(Triangle { v0: vi_prev, v1: vi, v2: vi_next });
                front.remove(fi);
            }
        } else {
            // Large angle: add two new vertices
            let angle_step = angle / 3.0;
            let base_angle = d_prev.y.atan2(d_prev.x); // simplified

            let mut new_indices = Vec::new();
            for j in 1..=2 {
                let a = angle_step * j as f64;
                let dir = Vec3::new(
                    (d_prev.x / len_prev * (1.0 - a / angle) + d_next.x / len_next * (a / angle)),
                    (d_prev.y / len_prev * (1.0 - a / angle) + d_next.y / len_next * (a / angle)),
                    (d_prev.z / len_prev * (1.0 - a / angle) + d_next.z / len_next * (a / angle)),
                );
                let dir_len = (dir.x * dir.x + dir.y * dir.y + dir.z * dir.z).sqrt();
                if dir_len < 1e-15 { continue; }

                let guess = Vec3::new(
                    v_curr.pos.x + delta * dir.x / dir_len,
                    v_curr.pos.y + delta * dir.y / dir_len,
                    v_curr.pos.z + delta * dir.z / dir_len,
                );

                if let Some(v_new) = project_onto_surface(surface, target_pot, &guess, &v_curr.normal, max_iter) {
                    let idx = vertices.len();
                    vertices.push(v_new);
                    new_indices.push(idx);
                }
            }

            if new_indices.len() == 2 {
                triangles.push(Triangle { v0: vi_prev, v1: vi, v2: new_indices[0] });
                triangles.push(Triangle { v0: vi, v1: new_indices[1], v2: new_indices[0] });
                triangles.push(Triangle { v0: vi, v1: vi_next, v2: new_indices[1] });
                // Replace vi in front with the two new vertices
                front[fi] = new_indices[1];
                front.insert(fi, new_indices[0]);
            } else {
                // Fallback
                triangles.push(Triangle { v0: vi_prev, v1: vi, v2: vi_next });
                front.remove(fi);
            }
        }
    }

    Ok(TriMesh { vertices, triangles })
}
