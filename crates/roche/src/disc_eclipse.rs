use lcurve_subs::{sqr, Vec3, TWOPI};
use std::f64::consts::PI;

/// Outcome of LOSC (line-of-sight cone) intersection with a circle.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Circle {
    /// Line of sight cone starts at or above the circle of interest
    Above,
    /// Line of sight circle is everywhere inside circle of interest
    Inside,
    /// Line of sight circle is everywhere outside circle of interest
    Outside,
    /// Line of sight circle is separated from the circle of interest
    Separate,
    /// Line of sight circle cone intersects the circle of interest
    Crossing,
}

/// Helper: compute the phase at which two circles intersect.
///
/// Given a circle of radius `rcone` whose centre is at distance `rxy` from
/// the centre of a circle of radius `radius`, returns the half-angle (in
/// phase units, 0-0.5) subtended by the intersection arc.
///
/// Preconditions (not checked):
///   rxy + rcone > radius
///   rxy < radius + rcone
///   rcone < radius + rxy
fn cut_phase(rxy: f64, rcone: f64, radius: f64) -> f64 {
    let cos_val = (sqr(rxy) + sqr(rcone) - sqr(radius)) / (2.0 * rcone * rxy);
    cos_val.clamp(-1.0, 1.0).acos() / TWOPI
}

/// Determines the nature of intersection of the LOSC cone to a point with a
/// circle centred on the z axis.
///
/// * `rxy` - distance of the point from the z axis
/// * `z` - z coordinate of the point
/// * `zcirc` - z ordinate of the circle
/// * `radius` - radius of the circle
/// * `tani` - tangent of the inclination
///
/// Returns (outcome, phase) where phase is only meaningful for `Crossing`.
fn circle_eclipse(rxy: f64, z: f64, zcirc: f64, radius: f64, tani: f64) -> (Circle, f64) {
    // Point above circle
    if z >= zcirc {
        return (Circle::Above, 0.0);
    }

    let rcone = tani * (zcirc - z);

    // Line-of-sight always outside the circle
    if rcone >= rxy + radius {
        return (Circle::Outside, 0.0);
    }

    // Line-of-sight circle separate from the circle
    if rxy >= rcone + radius {
        return (Circle::Separate, 0.0);
    }

    // Line-of-sight always inside the circle
    if rxy + rcone <= radius {
        return (Circle::Inside, 0.0);
    }

    // Crossing case
    let phase = cut_phase(rxy, rcone, radius);
    (Circle::Crossing, phase)
}

/// Compute phase ranges during which a flared accretion disc eclipses point `r`.
///
/// The disc is cylindrically symmetric with height h(R) = height * R^beta.
///
/// * `iangle` - orbital inclination (degrees)
/// * `rdisc1` - inner disc radius (binary separation units)
/// * `rdisc2` - outer disc radius
/// * `beta` - flaring exponent
/// * `height` - height at R=1
/// * `r` - point being tested
///
/// Returns vector of (ingress, egress) phase pairs.
pub fn disc_eclipse(
    iangle: f64,
    rdisc1: f64,
    rdisc2: f64,
    beta: f64,
    height: f64,
    r: &Vec3,
) -> Vec<(f64, f64)> {
    let mut temp = Vec::new();

    // Compute height of disc at outer boundary
    let hout = height * rdisc2.powf(beta);

    // Deal with points too high ever to be eclipsed whatever the inclination
    if r.z >= hout {
        return temp;
    }

    // Compute cosine and sine of inclination
    let irad = iangle * PI / 180.0;
    let sini = irad.sin();
    let cosi = irad.cos();

    // Special case of exactly edge-on
    if cosi == 0.0 {
        if r.z.abs() < hout {
            let rxy = (sqr(r.x) + sqr(r.y)).sqrt();
            if rxy <= rdisc2 {
                temp.push((0.0, 1.0));
            } else {
                let subtend = (rdisc2 / rxy).asin() / TWOPI;
                let pcen = r.y.atan2(-r.x) / TWOPI;
                let ingress = pcen - subtend;
                let ingress = ingress - ingress.floor();
                let egress = ingress + 2.0 * subtend;
                temp.push((ingress, egress));
            }
        }
        return temp;
    }

    // Work out distance from axis
    let rxy = (sqr(r.x) + sqr(r.y)).sqrt();

    if rdisc1 < rxy && rxy < rdisc2 && r.z.abs() < height * rxy.powf(beta) {
        // Point is inside disc and so is eclipsed
        temp.push((0.0, 1.1));
        return temp;
    }

    let tani = sini / cosi;

    if rxy < rdisc2 && r.z >= height * rdisc1.max(rxy).powf(beta) {
        // Point is in approximately conical region above the disc. Just need to check whether
        // it is not occulted by the edge of the disc
        let (result, phase) = circle_eclipse(rxy, r.z, hout, rdisc2, tani);

        if result == Circle::Outside {
            // point will be occulted by the disc edge at all phases
            temp.push((0.0, 1.1));
        } else if result == Circle::Crossing {
            // point partially occulted by disc edge; work out phases
            let phi0 = r.y.atan2(-r.x) / TWOPI;
            let ingress = phi0 + phase;
            let ingress = ingress - ingress.floor();
            let egress = ingress + 1.0 - 2.0 * phase;
            temp.push((ingress, egress));
        }
        return temp;
    }

    // Compute the radius of circle formed by LOSC in the plane of
    // the lower outer rim of the disc
    let rcone_lo = (tani * (-hout - r.z)).max(0.0);

    // Circle encloses rim, so no intersection
    if rcone_lo >= rxy + rdisc2 {
        return temp;
    }

    // Compute the radius of circle formed by LOSC in the plane of
    // the upper outer rim of the disc
    let rcone_hi = tani * (hout - r.z);

    // Circle disjoint from rim, so no intersection
    if rxy >= rcone_hi + rdisc2 {
        return temp;
    }

    // For the moment we pretend that the disc has no hole at its centre, so
    // that we are simply interested in the phases over which eclipse occurs.
    // At this point we are guaranteed that this will happen. All events are
    // symmetrically located around a phase defined by x and y only which will
    // be calculated at the end. We therefore just find the half range which
    // is called 'eclipse_phase' below.

    let eclipse_phase;
    if rxy + rcone_lo <= rdisc2 {
        // Cone swept out by line of sight always inside lower face so total eclipse
        eclipse_phase = 0.5;
    } else if rxy <= rdisc2 {
        // Points that project close to the z axis which are only
        // partially obscured by the disc hovering above them.
        // this means they must be below -HOUT
        eclipse_phase = cut_phase(rxy, rcone_lo, rdisc2);
    } else {
        // Points further from the z axis than the outer rim of the disc that will be eclipsed.
        if sqr(rcone_hi) + sqr(rdisc2) >= sqr(rxy)
            && sqr(rcone_lo) + sqr(rdisc2) <= sqr(rxy)
        {
            // In this case it is the curved outer disc rim that sets the limit
            eclipse_phase = (rdisc2 / rxy).asin() / TWOPI;
        } else if sqr(rcone_hi) + sqr(rdisc2) < sqr(rxy) {
            // In this case it is upper outer rim that sets the limit
            eclipse_phase = cut_phase(rxy, rcone_hi, rdisc2);
        } else {
            // In this case it is lower outer rim that sets the limit
            eclipse_phase = cut_phase(rxy, rcone_lo, rdisc2);
        }
    }

    // At this point we have covered all cases for the eclipse, whilst ignoring the
    // possibility of seeing the point through the hole in the middle of the disc.
    // Now let's calculate the 'appear_phase' if any.

    // First compute height of disc at inner boundary
    let hin = height * rdisc1.powf(beta);

    let mut appear_phase: f64 = -1.0;

    if r.z < -hout {
        // In this case the LOSC has to run through 4 circles which are the upper and
        // lower outer and inner rims.

        // First, the lower outer rim
        let (result, phase) = circle_eclipse(rxy, r.z, -hout, rdisc2, tani);
        if result == Circle::Inside {
            appear_phase = 0.5;
        } else if result == Circle::Crossing {
            appear_phase = phase;
        }

        // Second, the lower inner rim
        if appear_phase > 0.0 {
            let (result, phase) = circle_eclipse(rxy, r.z, -hin, rdisc1, tani);
            if result == Circle::Crossing {
                appear_phase = appear_phase.min(phase);
            } else if result != Circle::Inside {
                appear_phase = -1.0;
            }
        }

        // Third, the upper inner rim
        if appear_phase > 0.0 {
            let (result, phase) = circle_eclipse(rxy, r.z, hin, rdisc1, tani);
            if result == Circle::Crossing {
                appear_phase = appear_phase.min(phase);
            } else if result != Circle::Inside {
                appear_phase = -1.0;
            }
        }

        // Fourth, the upper outer rim
        if appear_phase > 0.0 {
            let (result, phase) = circle_eclipse(rxy, r.z, hout, rdisc2, tani);
            if result == Circle::Crossing {
                appear_phase = appear_phase.min(phase);
            } else if result != Circle::Inside {
                appear_phase = -1.0;
            }
        }
    } else if rxy < rdisc1 {
        if r.z < -hin {
            // Points hovering around underside of disc. Have to consider just three circles

            // First, the lower inner rim
            let (result, phase) = circle_eclipse(rxy, r.z, -hin, rdisc1, tani);
            if result == Circle::Inside {
                appear_phase = 0.5;
            } else if result == Circle::Crossing {
                appear_phase = phase;
            }

            // Second, the upper inner rim
            if appear_phase > 0.0 {
                let (result, phase) = circle_eclipse(rxy, r.z, hin, rdisc1, tani);
                if result == Circle::Crossing {
                    appear_phase = appear_phase.min(phase);
                } else if result != Circle::Inside {
                    appear_phase = -1.0;
                }
            }

            // Third, the upper outer rim
            if appear_phase > 0.0 {
                let (result, phase) = circle_eclipse(rxy, r.z, hout, rdisc2, tani);
                if result == Circle::Crossing {
                    appear_phase = appear_phase.min(phase);
                } else if result != Circle::Inside {
                    appear_phase = -1.0;
                }
            }
        } else if r.z < hin {
            // Points inside hole in middle of disc. Have to consider just two circles

            // First, the upper inner rim
            let (result, phase) = circle_eclipse(rxy, r.z, hin, rdisc1, tani);
            if result == Circle::Inside {
                appear_phase = 0.5;
            } else if result == Circle::Crossing {
                appear_phase = phase;
            }

            // Second, the upper outer rim
            if appear_phase > 0.0 {
                let (result, phase) = circle_eclipse(rxy, r.z, hout, rdisc2, tani);
                if result == Circle::Crossing {
                    appear_phase = appear_phase.min(phase);
                } else if result != Circle::Inside {
                    appear_phase = -1.0;
                }
            }
        }
    }

    // Here is the central phase
    let phi0 = r.y.atan2(-r.x) / TWOPI;

    if appear_phase <= 0.0 {
        let ingress = phi0 - eclipse_phase;
        let ingress = ingress - ingress.floor();
        let egress = ingress + 2.0 * eclipse_phase;
        temp.push((ingress, egress));
    } else if appear_phase < eclipse_phase {
        let ingress = phi0 - eclipse_phase;
        let ingress = ingress - ingress.floor();
        let egress = ingress + (eclipse_phase - appear_phase);
        temp.push((ingress, egress));

        let ingress = phi0 + appear_phase;
        let ingress = ingress - ingress.floor();
        let egress = ingress + (eclipse_phase - appear_phase);
        temp.push((ingress, egress));
    }

    temp
}

/// Simple disc eclipse test for a single line-of-sight.
///
/// Tests if the line from point `r` toward `earth` intersects the disc.
pub fn disc_eclipse_los(
    _iangle: f64,
    rdisc1: f64,
    rdisc2: f64,
    beta: f64,
    height: f64,
    r: &Vec3,
    earth: &Vec3,
) -> bool {
    let nsteps = 100;
    let lam_max = 3.0;

    for i in 1..=nsteps {
        let lam = lam_max * i as f64 / nsteps as f64;
        let point = *r + lam * *earth;
        let rxy = (sqr(point.x) + sqr(point.y)).sqrt();
        if rxy >= rdisc1 && rxy <= rdisc2 {
            let h = height * rxy.powf(beta);
            if point.z.abs() <= h {
                return true;
            }
        }
    }
    false
}
