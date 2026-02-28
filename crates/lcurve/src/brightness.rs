use lcurve_subs::{sqr, Vec3, PI, planck, dlpdlt};
use rayon::prelude::*;
use crate::model::Model;
use crate::types::Point;
use crate::LcurveError;

// EFAC = sqrt(8 ln 2) for FWHM→sigma conversion
const EFAC: f64 = 2.3548200450309493;

/// Compute star continuum brightness*area for both stars.
///
/// Accounts for gravity darkening, irradiation, star spots, and mirror reflection.
pub fn set_star_continuum(
    mdl: &Model,
    star1: &mut [Point],
    star2: &mut [Point],
) -> Result<(), LcurveError> {
    let (mut r1, mut r2) = mdl.get_r1r2();
    let rl1 = lcurve_roche::lagrange::xl1(mdl.q.value)?;
    if r1 < 0.0 { r1 = rl1; }
    let rl2 = 1.0 - rl1;
    if r2 < 0.0 { r2 = rl2; }

    let cofm2 = Vec3::new(1.0, 0.0, 0.0);

    // Gravity darkening exponent for star 1
    let gdcbol1 = if mdl.gdark_bolom1 {
        mdl.gravity_dark1.value
    } else {
        mdl.gravity_dark1.value / dlpdlt(mdl.wavelength, mdl.t1.value)
    };

    // Compute star spot directions for star 1
    let spot11 = spot_info(mdl.stsp11_long, mdl.stsp11_lat, mdl.stsp11_fwhm, mdl.stsp11_tcen);
    let spot12 = spot_info(mdl.stsp12_long, mdl.stsp12_lat, mdl.stsp12_fwhm, mdl.stsp12_tcen);
    let spot13 = spot_info(mdl.stsp13_long, mdl.stsp13_lat, mdl.stsp13_fwhm, mdl.stsp13_tcen);

    // Uniform equatorial spot
    let is_uespot = mdl.uesp_long1.defined && mdl.uesp_long2.defined
        && mdl.uesp_lathw.defined && mdl.uesp_taper.defined && mdl.uesp_temp.defined;
    let (uespot_dir, longhw) = if is_uespot {
        let lhw = (mdl.uesp_long2.value - mdl.uesp_long1.value) / 2.0;
        let lcen = ((mdl.uesp_long1.value + mdl.uesp_long2.value) / 2.0).to_radians();
        (Vec3::new(lcen.cos(), lcen.sin(), 0.0), lhw)
    } else {
        (Vec3::ZERO, 0.0)
    };

    // Star 1 elements
    star1.par_iter_mut().for_each(|pt| {
        let vec = cofm2 - pt.posn;
        let r = vec.length();
        let mu = Vec3::dot(&pt.dirn, &vec) / r;

        // Compute temperature with spots
        let mut t1 = mdl.t1.value;
        if let Some((dir, fwhm, tcen)) = &spot11 {
            let dist = angular_dist(&pt.posn, dir, pt.posn.length());
            t1 += (tcen - mdl.t1.value) * (-sqr(dist / (fwhm / EFAC)) / 2.0).exp();
        }
        if let Some((dir, fwhm, tcen)) = &spot12 {
            let dist = angular_dist(&pt.posn, dir, pt.posn.length());
            t1 += (tcen - mdl.t1.value) * (-sqr(dist / (fwhm / EFAC)) / 2.0).exp();
        }
        if let Some((dir, fwhm, tcen)) = &spot13 {
            let dist = angular_dist(&pt.posn, dir, pt.posn.length());
            t1 += (tcen - mdl.t1.value) * (-sqr(dist / (fwhm / EFAC)) / 2.0).exp();
        }

        if is_uespot {
            let evec = Vec3::new(pt.posn.x, pt.posn.y, 0.0);
            let evecl = evec.length();
            if evecl > 0.0 {
                let rd = pt.posn.length();
                let clat = (Vec3::dot(&pt.posn, &evec) / evecl / rd).clamp(-1.0, 1.0);
                let dlat = clat.acos().to_degrees();
                let clong = (Vec3::dot(&uespot_dir, &evec) / evecl).clamp(-1.0, 1.0);
                let dlong = clong.acos().to_degrees();

                if dlat <= mdl.uesp_lathw.value && dlong <= longhw {
                    t1 = mdl.uesp_temp.value;
                } else if dlat <= mdl.uesp_lathw.value {
                    t1 += (mdl.uesp_temp.value - mdl.t1.value)
                        * (-(dlong - longhw) / mdl.uesp_taper.value).exp();
                } else if dlong <= longhw {
                    t1 += (mdl.uesp_temp.value - mdl.t1.value)
                        * (-(dlat - mdl.uesp_lathw.value) / mdl.uesp_taper.value).exp();
                } else {
                    t1 += (mdl.uesp_temp.value - mdl.t1.value)
                        * (-(dlat - mdl.uesp_lathw.value) / mdl.uesp_taper.value).exp()
                        * (-(dlong - longhw) / mdl.uesp_taper.value).exp();
                }
            }
        }

        let (geom, temp) = irradiation(t1, pt.gravity as f64, gdcbol1,
                                        mdl.absorb.value, mdl.t2.value, mu, r, r2);

        pt.flux = (pt.area as f64 * planck(mdl.wavelength, temp)) as f32;

        if mdl.mirror {
            pt.flux += (pt.area as f64 * geom * planck(mdl.wavelength, mdl.t2.value.abs())) as f32;
        }
    });

    // Star 2 elements
    let cofm1 = Vec3::ZERO;
    let t2_abs = mdl.t2.value.abs();

    let gdcbol2 = if mdl.gdark_bolom2 {
        mdl.gravity_dark2.value
    } else {
        mdl.gravity_dark2.value / dlpdlt(mdl.wavelength, t2_abs)
    };

    let spot21 = spot_info(mdl.stsp21_long, mdl.stsp21_lat, mdl.stsp21_fwhm, mdl.stsp21_tcen);
    let spot22 = spot_info(mdl.stsp22_long, mdl.stsp22_lat, mdl.stsp22_fwhm, mdl.stsp22_tcen);

    star2.par_iter_mut().for_each(|pt| {
        let vec = cofm1 - pt.posn;
        let r = vec.length();
        let mu = Vec3::dot(&pt.dirn, &vec) / r;

        let mut t2 = t2_abs;
        if let Some((dir, fwhm, tcen)) = &spot21 {
            let off = pt.posn - cofm2;
            let dist = angular_dist(&off, dir, off.length());
            t2 += (tcen - t2) * (-sqr(dist / (fwhm / EFAC)) / 2.0).exp();
        }
        if let Some((dir, fwhm, tcen)) = &spot22 {
            let off = pt.posn - cofm2;
            let dist = angular_dist(&off, dir, off.length());
            t2 += (tcen - t2) * (-sqr(dist / (fwhm / EFAC)) / 2.0).exp();
        }

        let (geom, temp) = irradiation(t2, pt.gravity as f64, gdcbol2,
                                        mdl.absorb.value, mdl.t1.value, mu, r, r1);

        pt.flux = (pt.area as f64 * planck(mdl.wavelength, temp)) as f32;

        if mdl.mirror {
            pt.flux += (pt.area as f64 * geom * planck(mdl.wavelength, mdl.t1.value)) as f32;
        }
    });

    Ok(())
}

/// Compute irradiation temperature and geometric factor.
fn irradiation(
    t_base: f64, gravity: f64, gdcbol: f64,
    absorb: f64, t_irr: f64, mu: f64, r: f64, r_irr: f64,
) -> (f64, f64) {
    let t_grav = t_base * (gravity as f64).powf(gdcbol);

    if mu >= r_irr {
        // Full irradiation
        let geom = sqr(r_irr / r) * mu;
        let temp = (t_grav.powi(4) + absorb * t_irr.powi(4) * geom).powf(0.25);
        (geom, temp)
    } else if mu > -r_irr {
        // Sunset case
        let x0 = -mu / r_irr;
        let geom = sqr(r_irr / r) * r_irr
            * ((1.0 - x0 * x0).sqrt() * (2.0 + x0 * x0) / 3.0
                - x0 * (PI / 2.0 - x0.asin()))
            / PI;
        let temp = (t_grav.powi(4) + absorb * t_irr.powi(4) * geom).powf(0.25);
        (geom, temp)
    } else {
        // No irradiation
        (0.0, t_grav)
    }
}

/// Helper: compute angular distance in degrees between a position vector and a spot direction.
fn angular_dist(posn: &Vec3, spot_dir: &Vec3, posn_len: f64) -> f64 {
    if posn_len <= 0.0 {
        return 180.0;
    }
    let cos_ang = (Vec3::dot(posn, spot_dir) / posn_len).clamp(-1.0, 1.0);
    cos_ang.acos().to_degrees()
}

/// Compute spot direction and parameters.
/// Returns Some((direction, fwhm, tcen)) if the spot is defined.
pub fn spot_info(
    long_p: crate::model::Pparam,
    lat_p: crate::model::Pparam,
    fwhm_p: crate::model::Pparam,
    tcen_p: crate::model::Pparam,
) -> Option<(Vec3, f64, f64)> {
    if !long_p.defined || !lat_p.defined || !fwhm_p.defined || !tcen_p.defined {
        return None;
    }
    let clong = long_p.value.to_radians().cos();
    let slong = long_p.value.to_radians().sin();
    let clat = lat_p.value.to_radians().cos();
    let slat = lat_p.value.to_radians().sin();
    Some((Vec3::new(clat * clong, clat * slong, slat), fwhm_p.value, tcen_p.value))
}

/// Set disc surface brightness (power law with radius).
pub fn set_disc_continuum(
    rdisc: f64, tdisc: f64, texp: f64, wave: f64, disc: &mut [Point],
) {
    let bright = planck(wave, tdisc);
    disc.par_iter_mut().for_each(|pt| {
        let r = pt.posn.length();
        pt.flux = (bright * (r / rdisc).powf(texp) * pt.area as f64) as f32;
    });
}

/// Set disc edge surface brightness with irradiation from secondary.
pub fn set_edge_continuum(
    tedge: f64, r2: f64, t2: f64, absorb: f64, wave: f64, edge: &mut [Point],
) {
    let cofm2 = Vec3::new(1.0, 0.0, 0.0);
    edge.par_iter_mut().for_each(|pt| {
        let vec = cofm2 - pt.posn;
        let r = vec.length();
        let mu = Vec3::dot(&pt.dirn, &vec) / r;

        let temp = if mu >= r2 {
            let geom = sqr(r2 / r) * mu;
            (tedge.powi(4) + absorb * t2.powi(4) * geom).powf(0.25)
        } else if mu > -r2 {
            let x0 = -mu / r2;
            let geom = sqr(r2 / r) * r2
                * ((1.0 - x0 * x0).sqrt() * (2.0 + x0 * x0) / 3.0
                    - x0 * (PI / 2.0 - x0.asin()))
                / PI;
            (tedge.powi(4) + absorb * t2.powi(4) * geom).powf(0.25)
        } else {
            tedge
        };
        pt.flux = (pt.area as f64 * planck(wave, temp)) as f32;
    });
}
