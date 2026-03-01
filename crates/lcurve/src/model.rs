use std::collections::HashMap;
use crate::types::{LDC, LDCType};
use crate::LcurveError;

/// Physical parameter: value + fitting metadata.
#[derive(Debug, Clone, Copy)]
pub struct Pparam {
    pub value: f64,
    pub range: f64,
    pub dstep: f64,
    pub vary: bool,
    pub defined: bool,
}

impl Default for Pparam {
    fn default() -> Self {
        Pparam { value: 0.0, range: 0.0, dstep: 0.0, vary: false, defined: false }
    }
}

impl Pparam {
    /// Parse from a string like "0.12405 0.001 0.0001 0 1"
    pub fn from_str(s: &str) -> Result<Self, LcurveError> {
        let parts: Vec<&str> = s.split_whitespace().collect();
        if parts.len() < 4 {
            return Err(LcurveError::Parse(format!("Pparam: too few fields in '{}'", s)));
        }
        let value: f64 = parts[0].parse().map_err(|_| LcurveError::Parse(format!("bad value: {}", parts[0])))?;
        let range: f64 = parts[1].parse().map_err(|_| LcurveError::Parse(format!("bad range: {}", parts[1])))?;
        let dstep: f64 = parts[2].parse().map_err(|_| LcurveError::Parse(format!("bad dstep: {}", parts[2])))?;
        let vary = parse_bool(parts[3])?;
        let defined = if parts.len() > 4 {
            parse_bool(parts[4])?
        } else {
            true
        };
        Ok(Pparam { value, range, dstep, vary, defined })
    }
}

/// Implicit conversion to f64.
impl From<Pparam> for f64 {
    fn from(p: Pparam) -> f64 {
        p.value
    }
}

fn parse_bool(s: &str) -> Result<bool, LcurveError> {
    match s.trim() {
        "1" | "true" | "yes" => Ok(true),
        "0" | "false" | "no" => Ok(false),
        _ => Err(LcurveError::Parse(format!("bad bool: '{}'", s))),
    }
}

fn parse_double(s: &str) -> Result<f64, LcurveError> {
    s.trim().parse::<f64>().map_err(|_| LcurveError::Parse(format!("bad double: '{}'", s)))
}

fn parse_int(s: &str) -> Result<i32, LcurveError> {
    s.trim().parse::<i32>().map_err(|_| LcurveError::Parse(format!("bad int: '{}'", s)))
}

/// The complete binary light curve model.
#[derive(Debug, Clone)]
pub struct Model {
    // Physical parameters
    pub q: Pparam,
    pub iangle: Pparam,
    pub r1: Pparam,
    pub r2: Pparam,
    pub cphi3: Pparam,
    pub cphi4: Pparam,
    pub spin1: Pparam,
    pub spin2: Pparam,
    pub t1: Pparam,
    pub t2: Pparam,
    pub ldc1_1: Pparam,
    pub ldc1_2: Pparam,
    pub ldc1_3: Pparam,
    pub ldc1_4: Pparam,
    pub ldc2_1: Pparam,
    pub ldc2_2: Pparam,
    pub ldc2_3: Pparam,
    pub ldc2_4: Pparam,
    pub velocity_scale: Pparam,
    pub beam_factor1: Pparam,
    pub beam_factor2: Pparam,
    pub t0: Pparam,
    pub period: Pparam,
    pub pdot: Pparam,
    pub deltat: Pparam,
    pub gravity_dark1: Pparam,
    pub gravity_dark2: Pparam,
    pub absorb: Pparam,
    pub slope: Pparam,
    pub quad: Pparam,
    pub cube: Pparam,
    pub third: Pparam,
    pub rdisc1: Pparam,
    pub rdisc2: Pparam,
    pub height_disc: Pparam,
    pub beta_disc: Pparam,
    pub temp_disc: Pparam,
    pub texp_disc: Pparam,
    pub lin_limb_disc: Pparam,
    pub quad_limb_disc: Pparam,
    pub temp_edge: Pparam,
    pub absorb_edge: Pparam,
    pub radius_spot: Pparam,
    pub length_spot: Pparam,
    pub height_spot: Pparam,
    pub expon_spot: Pparam,
    pub epow_spot: Pparam,
    pub angle_spot: Pparam,
    pub yaw_spot: Pparam,
    pub temp_spot: Pparam,
    pub tilt_spot: Pparam,
    pub cfrac_spot: Pparam,
    // Star spots
    pub stsp11_long: Pparam,
    pub stsp11_lat: Pparam,
    pub stsp11_fwhm: Pparam,
    pub stsp11_tcen: Pparam,
    pub stsp12_long: Pparam,
    pub stsp12_lat: Pparam,
    pub stsp12_fwhm: Pparam,
    pub stsp12_tcen: Pparam,
    pub stsp13_long: Pparam,
    pub stsp13_lat: Pparam,
    pub stsp13_fwhm: Pparam,
    pub stsp13_tcen: Pparam,
    pub stsp21_long: Pparam,
    pub stsp21_lat: Pparam,
    pub stsp21_fwhm: Pparam,
    pub stsp21_tcen: Pparam,
    pub stsp22_long: Pparam,
    pub stsp22_lat: Pparam,
    pub stsp22_fwhm: Pparam,
    pub stsp22_tcen: Pparam,
    pub uesp_long1: Pparam,
    pub uesp_long2: Pparam,
    pub uesp_lathw: Pparam,
    pub uesp_taper: Pparam,
    pub uesp_temp: Pparam,

    // Computational parameters
    pub delta_phase: f64,
    pub nlat1f: i32,
    pub nlat2f: i32,
    pub nlat1c: i32,
    pub nlat2c: i32,
    pub npole: bool,
    pub nlatfill: i32,
    pub nlngfill: i32,
    pub lfudge: f64,
    pub llo: f64,
    pub lhi: f64,
    pub phase1: f64,
    pub phase2: f64,
    pub wavelength: f64,
    pub roche1: bool,
    pub roche2: bool,
    pub eclipse1: bool,
    pub eclipse2: bool,
    pub glens1: bool,
    pub use_radii: bool,
    pub tperiod: f64,
    pub gdark_bolom1: bool,
    pub gdark_bolom2: bool,
    pub mucrit1: f64,
    pub mucrit2: f64,
    pub limb1: LDCType,
    pub limb2: LDCType,
    pub mirror: bool,
    pub add_disc: bool,
    pub nrad: i32,
    pub opaque: bool,
    pub add_spot: bool,
    pub nspot: i32,
    pub iscale: bool,
}

impl Model {
    /// Load model from a parameter file.
    pub fn from_file(path: &str) -> Result<Self, LcurveError> {
        let contents = std::fs::read_to_string(path)
            .map_err(|e| LcurveError::Io(e))?;

        let mut pv: HashMap<String, String> = HashMap::new();
        let mut found: HashMap<String, bool> = HashMap::new();

        // Register all expected parameter names
        let phys_names = [
            "q", "iangle", "r1", "r2", "cphi3", "cphi4", "spin1", "spin2",
            "t1", "t2", "ldc1_1", "ldc1_2", "ldc1_3", "ldc1_4",
            "ldc2_1", "ldc2_2", "ldc2_3", "ldc2_4",
            "velocity_scale", "beam_factor1", "beam_factor2",
            "t0", "period", "pdot", "deltat",
            "gravity_dark1", "gravity_dark2", "absorb",
            "slope", "quad", "cube", "third",
            "rdisc1", "rdisc2", "height_disc", "beta_disc",
            "temp_disc", "texp_disc", "lin_limb_disc", "quad_limb_disc",
            "temp_edge", "absorb_edge",
            "radius_spot", "length_spot", "height_spot",
            "expon_spot", "epow_spot", "angle_spot", "yaw_spot",
            "temp_spot", "tilt_spot", "cfrac_spot",
            "stsp11_long", "stsp11_lat", "stsp11_fwhm", "stsp11_tcen",
            "stsp12_long", "stsp12_lat", "stsp12_fwhm", "stsp12_tcen",
            "stsp13_long", "stsp13_lat", "stsp13_fwhm", "stsp13_tcen",
            "stsp21_long", "stsp21_lat", "stsp21_fwhm", "stsp21_tcen",
            "stsp22_long", "stsp22_lat", "stsp22_fwhm", "stsp22_tcen",
            "uesp_long1", "uesp_long2", "uesp_lathw", "uesp_taper", "uesp_temp",
        ];
        let comp_names = [
            "delta_phase", "nlat1f", "nlat2f", "nlat1c", "nlat2c",
            "npole", "nlatfill", "nlngfill", "lfudge", "llo", "lhi",
            "phase1", "phase2", "wavelength",
            "roche1", "roche2", "eclipse1", "eclipse2", "glens1",
            "use_radii", "tperiod", "gdark_bolom1", "gdark_bolom2",
            "mucrit1", "mucrit2", "limb1", "limb2",
            "mirror", "add_disc", "nrad", "opaque", "add_spot", "nspot", "iscale",
        ];
        for name in phys_names.iter().chain(comp_names.iter()) {
            found.insert(name.to_string(), false);
        }

        // Parse file
        for (line_num, line) in contents.lines().enumerate() {
            let trimmed = line.trim();
            if trimmed.is_empty() || trimmed.starts_with('#') || trimmed.starts_with(' ') || trimmed.starts_with('\t') {
                continue;
            }
            if let Some(eq_pos) = trimmed.find('=') {
                let param = trimmed[..eq_pos].trim().to_string();
                let mut value = trimmed[eq_pos + 1..].to_string();

                // Strip comments (# not preceded by \)
                let mut nhash = 0usize;
                while let Some(pos) = value[nhash..].find('#') {
                    let abs_pos = nhash + pos;
                    if abs_pos > 0 && value.as_bytes()[abs_pos - 1] == b'\\' {
                        nhash = abs_pos + 1;
                    } else {
                        value.truncate(abs_pos);
                        break;
                    }
                }
                let value = value.trim().to_string();

                if !found.contains_key(&param) {
                    return Err(LcurveError::Parse(format!(
                        "Parameter '{}' on line {} was not recognised", param, line_num + 1)));
                }
                found.insert(param.clone(), true);
                pv.insert(param, value);
            }
        }

        // Check required parameters, set defaults for optional ones
        let optional_with_defaults = ["pdot", "third", "temp_edge", "absorb_edge"];
        for (name, was_found) in &found {
            if !was_found {
                if optional_with_defaults.contains(&name.as_str()) {
                    pv.insert(name.clone(), "0 1.e-10 1.e-10 0 0".to_string());
                } else if !name.starts_with("stsp") && !name.starts_with("uesp") {
                    return Err(LcurveError::Parse(format!(
                        "Parameter '{}' was not defined", name)));
                }
            }
        }

        // Helper to get Pparam or default
        let get_pparam = |name: &str| -> Result<Pparam, LcurveError> {
            if let Some(val) = pv.get(name) {
                Pparam::from_str(val)
            } else {
                Ok(Pparam::default())
            }
        };

        // Parse star spot groups: if any in group defined, all must be
        let spot_groups = [
            &["stsp11_long", "stsp11_lat", "stsp11_fwhm", "stsp11_tcen"][..],
            &["stsp12_long", "stsp12_lat", "stsp12_fwhm", "stsp12_tcen"],
            &["stsp13_long", "stsp13_lat", "stsp13_fwhm", "stsp13_tcen"],
            &["stsp21_long", "stsp21_lat", "stsp21_fwhm", "stsp21_tcen"],
            &["stsp22_long", "stsp22_lat", "stsp22_fwhm", "stsp22_tcen"],
        ];
        for group in &spot_groups {
            let any = group.iter().any(|n| *found.get(*n).unwrap_or(&false));
            let all = group.iter().all(|n| *found.get(*n).unwrap_or(&false));
            if any && !all {
                return Err(LcurveError::Parse(format!(
                    "One or more of {} parameters were not initialised",
                    group[0].split('_').next().unwrap())));
            }
        }

        let uesp_group = &["uesp_long1", "uesp_long2", "uesp_lathw", "uesp_taper", "uesp_temp"];
        let uesp_any = uesp_group.iter().any(|n| *found.get(*n).unwrap_or(&false));
        let uesp_all = uesp_group.iter().all(|n| *found.get(*n).unwrap_or(&false));
        if uesp_any && !uesp_all {
            return Err(LcurveError::Parse(
                "One or more uniform equatorial spot parameters were not initialised".to_string()));
        }

        let limb1 = match pv.get("limb1").map(|s| s.as_str()) {
            Some("Poly") => LDCType::Poly,
            Some("Claret") => LDCType::Claret,
            Some(v) => return Err(LcurveError::Parse(format!("bad limb1: '{}'", v))),
            None => return Err(LcurveError::Parse("limb1 not defined".to_string())),
        };
        let limb2 = match pv.get("limb2").map(|s| s.as_str()) {
            Some("Poly") => LDCType::Poly,
            Some("Claret") => LDCType::Claret,
            Some(v) => return Err(LcurveError::Parse(format!("bad limb2: '{}'", v))),
            None => return Err(LcurveError::Parse("limb2 not defined".to_string())),
        };

        let roche1_val = parse_bool(pv.get("roche1").unwrap())?;
        let glens1_val = parse_bool(pv.get("glens1").unwrap())?;
        if glens1_val && roche1_val {
            return Err(LcurveError::Parse(
                "glens1=1 and roche1=1 are not simultaneously allowed".to_string()));
        }

        Ok(Model {
            q: Pparam::from_str(pv.get("q").unwrap())?,
            iangle: Pparam::from_str(pv.get("iangle").unwrap())?,
            r1: Pparam::from_str(pv.get("r1").unwrap())?,
            r2: Pparam::from_str(pv.get("r2").unwrap())?,
            cphi3: Pparam::from_str(pv.get("cphi3").unwrap())?,
            cphi4: Pparam::from_str(pv.get("cphi4").unwrap())?,
            spin1: Pparam::from_str(pv.get("spin1").unwrap())?,
            spin2: Pparam::from_str(pv.get("spin2").unwrap())?,
            t1: Pparam::from_str(pv.get("t1").unwrap())?,
            t2: Pparam::from_str(pv.get("t2").unwrap())?,
            ldc1_1: Pparam::from_str(pv.get("ldc1_1").unwrap())?,
            ldc1_2: Pparam::from_str(pv.get("ldc1_2").unwrap())?,
            ldc1_3: Pparam::from_str(pv.get("ldc1_3").unwrap())?,
            ldc1_4: Pparam::from_str(pv.get("ldc1_4").unwrap())?,
            ldc2_1: Pparam::from_str(pv.get("ldc2_1").unwrap())?,
            ldc2_2: Pparam::from_str(pv.get("ldc2_2").unwrap())?,
            ldc2_3: Pparam::from_str(pv.get("ldc2_3").unwrap())?,
            ldc2_4: Pparam::from_str(pv.get("ldc2_4").unwrap())?,
            velocity_scale: Pparam::from_str(pv.get("velocity_scale").unwrap())?,
            beam_factor1: Pparam::from_str(pv.get("beam_factor1").unwrap())?,
            beam_factor2: Pparam::from_str(pv.get("beam_factor2").unwrap())?,
            t0: Pparam::from_str(pv.get("t0").unwrap())?,
            period: Pparam::from_str(pv.get("period").unwrap())?,
            pdot: Pparam::from_str(pv.get("pdot").unwrap())?,
            deltat: Pparam::from_str(pv.get("deltat").unwrap())?,
            gravity_dark1: Pparam::from_str(pv.get("gravity_dark1").unwrap())?,
            gravity_dark2: Pparam::from_str(pv.get("gravity_dark2").unwrap())?,
            absorb: Pparam::from_str(pv.get("absorb").unwrap())?,
            slope: Pparam::from_str(pv.get("slope").unwrap())?,
            quad: Pparam::from_str(pv.get("quad").unwrap())?,
            cube: Pparam::from_str(pv.get("cube").unwrap())?,
            third: Pparam::from_str(pv.get("third").unwrap())?,
            rdisc1: Pparam::from_str(pv.get("rdisc1").unwrap())?,
            rdisc2: Pparam::from_str(pv.get("rdisc2").unwrap())?,
            height_disc: Pparam::from_str(pv.get("height_disc").unwrap())?,
            beta_disc: Pparam::from_str(pv.get("beta_disc").unwrap())?,
            temp_disc: Pparam::from_str(pv.get("temp_disc").unwrap())?,
            texp_disc: Pparam::from_str(pv.get("texp_disc").unwrap())?,
            lin_limb_disc: Pparam::from_str(pv.get("lin_limb_disc").unwrap())?,
            quad_limb_disc: Pparam::from_str(pv.get("quad_limb_disc").unwrap())?,
            temp_edge: Pparam::from_str(pv.get("temp_edge").unwrap())?,
            absorb_edge: Pparam::from_str(pv.get("absorb_edge").unwrap())?,
            radius_spot: Pparam::from_str(pv.get("radius_spot").unwrap())?,
            length_spot: Pparam::from_str(pv.get("length_spot").unwrap())?,
            height_spot: Pparam::from_str(pv.get("height_spot").unwrap())?,
            expon_spot: Pparam::from_str(pv.get("expon_spot").unwrap())?,
            epow_spot: Pparam::from_str(pv.get("epow_spot").unwrap())?,
            angle_spot: Pparam::from_str(pv.get("angle_spot").unwrap())?,
            yaw_spot: Pparam::from_str(pv.get("yaw_spot").unwrap())?,
            temp_spot: Pparam::from_str(pv.get("temp_spot").unwrap())?,
            tilt_spot: Pparam::from_str(pv.get("tilt_spot").unwrap())?,
            cfrac_spot: Pparam::from_str(pv.get("cfrac_spot").unwrap())?,
            stsp11_long: get_pparam("stsp11_long")?,
            stsp11_lat: get_pparam("stsp11_lat")?,
            stsp11_fwhm: get_pparam("stsp11_fwhm")?,
            stsp11_tcen: get_pparam("stsp11_tcen")?,
            stsp12_long: get_pparam("stsp12_long")?,
            stsp12_lat: get_pparam("stsp12_lat")?,
            stsp12_fwhm: get_pparam("stsp12_fwhm")?,
            stsp12_tcen: get_pparam("stsp12_tcen")?,
            stsp13_long: get_pparam("stsp13_long")?,
            stsp13_lat: get_pparam("stsp13_lat")?,
            stsp13_fwhm: get_pparam("stsp13_fwhm")?,
            stsp13_tcen: get_pparam("stsp13_tcen")?,
            stsp21_long: get_pparam("stsp21_long")?,
            stsp21_lat: get_pparam("stsp21_lat")?,
            stsp21_fwhm: get_pparam("stsp21_fwhm")?,
            stsp21_tcen: get_pparam("stsp21_tcen")?,
            stsp22_long: get_pparam("stsp22_long")?,
            stsp22_lat: get_pparam("stsp22_lat")?,
            stsp22_fwhm: get_pparam("stsp22_fwhm")?,
            stsp22_tcen: get_pparam("stsp22_tcen")?,
            uesp_long1: get_pparam("uesp_long1")?,
            uesp_long2: get_pparam("uesp_long2")?,
            uesp_lathw: get_pparam("uesp_lathw")?,
            uesp_taper: get_pparam("uesp_taper")?,
            uesp_temp: get_pparam("uesp_temp")?,

            delta_phase: parse_double(pv.get("delta_phase").unwrap())?,
            nlat1f: parse_int(pv.get("nlat1f").unwrap())?,
            nlat2f: parse_int(pv.get("nlat2f").unwrap())?,
            nlat1c: parse_int(pv.get("nlat1c").unwrap())?,
            nlat2c: parse_int(pv.get("nlat2c").unwrap())?,
            npole: parse_bool(pv.get("npole").unwrap())?,
            nlatfill: parse_int(pv.get("nlatfill").unwrap())?,
            nlngfill: parse_int(pv.get("nlngfill").unwrap())?,
            lfudge: parse_double(pv.get("lfudge").unwrap())?,
            llo: parse_double(pv.get("llo").unwrap())?,
            lhi: parse_double(pv.get("lhi").unwrap())?,
            phase1: parse_double(pv.get("phase1").unwrap())?,
            phase2: parse_double(pv.get("phase2").unwrap())?,
            wavelength: parse_double(pv.get("wavelength").unwrap())?,
            roche1: roche1_val,
            roche2: parse_bool(pv.get("roche2").unwrap())?,
            eclipse1: parse_bool(pv.get("eclipse1").unwrap())?,
            eclipse2: parse_bool(pv.get("eclipse2").unwrap())?,
            glens1: glens1_val,
            use_radii: parse_bool(pv.get("use_radii").unwrap())?,
            tperiod: parse_double(pv.get("tperiod").unwrap())?,
            gdark_bolom1: parse_bool(pv.get("gdark_bolom1").unwrap())?,
            gdark_bolom2: parse_bool(pv.get("gdark_bolom2").unwrap())?,
            mucrit1: parse_double(pv.get("mucrit1").unwrap())?,
            mucrit2: parse_double(pv.get("mucrit2").unwrap())?,
            limb1,
            limb2,
            mirror: parse_bool(pv.get("mirror").unwrap())?,
            add_disc: parse_bool(pv.get("add_disc").unwrap())?,
            nrad: parse_int(pv.get("nrad").unwrap())?,
            opaque: parse_bool(pv.get("opaque").unwrap())?,
            add_spot: parse_bool(pv.get("add_spot").unwrap())?,
            nspot: parse_int(pv.get("nspot").unwrap())?,
            iscale: parse_bool(pv.get("iscale").unwrap())?,
        })
    }

    /// Get r1, r2 accounting for use_radii vs contact phases.
    pub fn get_r1r2(&self) -> (f64, f64) {
        if self.use_radii {
            (self.r1.value, self.r2.value)
        } else {
            let sini = (self.iangle.value.to_radians()).sin();
            let r2pr1 = (1.0 - (sini * (2.0 * std::f64::consts::PI * self.cphi4.value).cos()).powi(2)).sqrt();
            let r2mr1 = (1.0 - (sini * (2.0 * std::f64::consts::PI * self.cphi3.value).cos()).powi(2)).sqrt();
            ((r2pr1 - r2mr1) / 2.0, (r2pr1 + r2mr1) / 2.0)
        }
    }

    /// Set a physical parameter's value by name.
    ///
    /// Returns `true` if the name was recognised, `false` otherwise.
    pub fn set_param_value(&mut self, name: &str, value: f64) -> bool {
        macro_rules! do_set {
            ($( $key:literal => $field:ident ),+ $(,)?) => {
                match name {
                    $( $key => { self.$field.value = value; true } )+
                    _ => false,
                }
            };
        }
        do_set! {
            "q" => q, "iangle" => iangle, "r1" => r1, "r2" => r2,
            "cphi3" => cphi3, "cphi4" => cphi4, "spin1" => spin1, "spin2" => spin2,
            "t1" => t1, "t2" => t2,
            "ldc1_1" => ldc1_1, "ldc1_2" => ldc1_2, "ldc1_3" => ldc1_3, "ldc1_4" => ldc1_4,
            "ldc2_1" => ldc2_1, "ldc2_2" => ldc2_2, "ldc2_3" => ldc2_3, "ldc2_4" => ldc2_4,
            "velocity_scale" => velocity_scale, "beam_factor1" => beam_factor1, "beam_factor2" => beam_factor2,
            "t0" => t0, "period" => period, "pdot" => pdot, "deltat" => deltat,
            "gravity_dark1" => gravity_dark1, "gravity_dark2" => gravity_dark2, "absorb" => absorb,
            "slope" => slope, "quad" => quad, "cube" => cube, "third" => third,
            "rdisc1" => rdisc1, "rdisc2" => rdisc2, "height_disc" => height_disc, "beta_disc" => beta_disc,
            "temp_disc" => temp_disc, "texp_disc" => texp_disc,
            "lin_limb_disc" => lin_limb_disc, "quad_limb_disc" => quad_limb_disc,
            "temp_edge" => temp_edge, "absorb_edge" => absorb_edge,
            "radius_spot" => radius_spot, "length_spot" => length_spot, "height_spot" => height_spot,
            "expon_spot" => expon_spot, "epow_spot" => epow_spot, "angle_spot" => angle_spot,
            "yaw_spot" => yaw_spot, "temp_spot" => temp_spot, "tilt_spot" => tilt_spot, "cfrac_spot" => cfrac_spot,
            "stsp11_long" => stsp11_long, "stsp11_lat" => stsp11_lat, "stsp11_fwhm" => stsp11_fwhm, "stsp11_tcen" => stsp11_tcen,
            "stsp12_long" => stsp12_long, "stsp12_lat" => stsp12_lat, "stsp12_fwhm" => stsp12_fwhm, "stsp12_tcen" => stsp12_tcen,
            "stsp13_long" => stsp13_long, "stsp13_lat" => stsp13_lat, "stsp13_fwhm" => stsp13_fwhm, "stsp13_tcen" => stsp13_tcen,
            "stsp21_long" => stsp21_long, "stsp21_lat" => stsp21_lat, "stsp21_fwhm" => stsp21_fwhm, "stsp21_tcen" => stsp21_tcen,
            "stsp22_long" => stsp22_long, "stsp22_lat" => stsp22_lat, "stsp22_fwhm" => stsp22_fwhm, "stsp22_tcen" => stsp22_tcen,
            "uesp_long1" => uesp_long1, "uesp_long2" => uesp_long2, "uesp_lathw" => uesp_lathw,
            "uesp_taper" => uesp_taper, "uesp_temp" => uesp_temp,
        }
    }

    /// Get a physical parameter's value by name.
    ///
    /// Returns `None` if the name was not recognised.
    pub fn get_param_value(&self, name: &str) -> Option<f64> {
        macro_rules! do_get {
            ($( $key:literal => $field:ident ),+ $(,)?) => {
                match name {
                    $( $key => Some(self.$field.value), )+
                    _ => None,
                }
            };
        }
        do_get! {
            "q" => q, "iangle" => iangle, "r1" => r1, "r2" => r2,
            "cphi3" => cphi3, "cphi4" => cphi4, "spin1" => spin1, "spin2" => spin2,
            "t1" => t1, "t2" => t2,
            "ldc1_1" => ldc1_1, "ldc1_2" => ldc1_2, "ldc1_3" => ldc1_3, "ldc1_4" => ldc1_4,
            "ldc2_1" => ldc2_1, "ldc2_2" => ldc2_2, "ldc2_3" => ldc2_3, "ldc2_4" => ldc2_4,
            "velocity_scale" => velocity_scale, "beam_factor1" => beam_factor1, "beam_factor2" => beam_factor2,
            "t0" => t0, "period" => period, "pdot" => pdot, "deltat" => deltat,
            "gravity_dark1" => gravity_dark1, "gravity_dark2" => gravity_dark2, "absorb" => absorb,
            "slope" => slope, "quad" => quad, "cube" => cube, "third" => third,
            "rdisc1" => rdisc1, "rdisc2" => rdisc2, "height_disc" => height_disc, "beta_disc" => beta_disc,
            "temp_disc" => temp_disc, "texp_disc" => texp_disc,
            "lin_limb_disc" => lin_limb_disc, "quad_limb_disc" => quad_limb_disc,
            "temp_edge" => temp_edge, "absorb_edge" => absorb_edge,
            "radius_spot" => radius_spot, "length_spot" => length_spot, "height_spot" => height_spot,
            "expon_spot" => expon_spot, "epow_spot" => epow_spot, "angle_spot" => angle_spot,
            "yaw_spot" => yaw_spot, "temp_spot" => temp_spot, "tilt_spot" => tilt_spot, "cfrac_spot" => cfrac_spot,
            "stsp11_long" => stsp11_long, "stsp11_lat" => stsp11_lat, "stsp11_fwhm" => stsp11_fwhm, "stsp11_tcen" => stsp11_tcen,
            "stsp12_long" => stsp12_long, "stsp12_lat" => stsp12_lat, "stsp12_fwhm" => stsp12_fwhm, "stsp12_tcen" => stsp12_tcen,
            "stsp13_long" => stsp13_long, "stsp13_lat" => stsp13_lat, "stsp13_fwhm" => stsp13_fwhm, "stsp13_tcen" => stsp13_tcen,
            "stsp21_long" => stsp21_long, "stsp21_lat" => stsp21_lat, "stsp21_fwhm" => stsp21_fwhm, "stsp21_tcen" => stsp21_tcen,
            "stsp22_long" => stsp22_long, "stsp22_lat" => stsp22_lat, "stsp22_fwhm" => stsp22_fwhm, "stsp22_tcen" => stsp22_tcen,
            "uesp_long1" => uesp_long1, "uesp_long2" => uesp_long2, "uesp_lathw" => uesp_lathw,
            "uesp_taper" => uesp_taper, "uesp_temp" => uesp_temp,
        }
    }

    /// Get LDC struct for star 1.
    pub fn get_ldc1(&self) -> LDC {
        LDC::new(
            self.ldc1_1.value, self.ldc1_2.value, self.ldc1_3.value, self.ldc1_4.value,
            self.mucrit1, self.limb1,
        )
    }

    /// Get LDC struct for star 2.
    pub fn get_ldc2(&self) -> LDC {
        LDC::new(
            self.ldc2_1.value, self.ldc2_2.value, self.ldc2_3.value, self.ldc2_4.value,
            self.mucrit2, self.limb2,
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::LDCType;

    /// Build a minimal Model without file I/O.
    fn make_model() -> Model {
        let p = |v: f64| Pparam { value: v, range: 0.1, dstep: 0.001, vary: false, defined: true };
        Model {
            q: p(0.5), iangle: p(82.0), r1: p(0.015), r2: p(-1.0),
            cphi3: p(0.0), cphi4: p(0.0), spin1: p(1.0), spin2: p(1.0),
            t1: p(15000.0), t2: p(3500.0),
            ldc1_1: p(0.0), ldc1_2: p(0.0), ldc1_3: p(0.0), ldc1_4: p(0.0),
            ldc2_1: p(0.0), ldc2_2: p(0.0), ldc2_3: p(0.0), ldc2_4: p(0.0),
            velocity_scale: p(300.0), beam_factor1: p(1.0), beam_factor2: p(1.0),
            t0: p(0.0), period: p(0.1), pdot: p(0.0), deltat: p(0.0),
            gravity_dark1: p(0.0), gravity_dark2: p(0.0), absorb: p(1.0),
            slope: p(0.0), quad: p(0.0), cube: p(0.0), third: p(0.0),
            rdisc1: p(0.0), rdisc2: p(0.0), height_disc: p(0.0), beta_disc: p(0.0),
            temp_disc: p(0.0), texp_disc: p(0.0), lin_limb_disc: p(0.0), quad_limb_disc: p(0.0),
            temp_edge: p(0.0), absorb_edge: p(0.0),
            radius_spot: p(0.0), length_spot: p(0.0), height_spot: p(0.0),
            expon_spot: p(0.0), epow_spot: p(0.0), angle_spot: p(0.0),
            yaw_spot: p(0.0), temp_spot: p(0.0), tilt_spot: p(0.0), cfrac_spot: p(0.0),
            stsp11_long: Pparam::default(), stsp11_lat: Pparam::default(),
            stsp11_fwhm: Pparam::default(), stsp11_tcen: Pparam::default(),
            stsp12_long: Pparam::default(), stsp12_lat: Pparam::default(),
            stsp12_fwhm: Pparam::default(), stsp12_tcen: Pparam::default(),
            stsp13_long: Pparam::default(), stsp13_lat: Pparam::default(),
            stsp13_fwhm: Pparam::default(), stsp13_tcen: Pparam::default(),
            stsp21_long: Pparam::default(), stsp21_lat: Pparam::default(),
            stsp21_fwhm: Pparam::default(), stsp21_tcen: Pparam::default(),
            stsp22_long: Pparam::default(), stsp22_lat: Pparam::default(),
            stsp22_fwhm: Pparam::default(), stsp22_tcen: Pparam::default(),
            uesp_long1: Pparam::default(), uesp_long2: Pparam::default(),
            uesp_lathw: Pparam::default(), uesp_taper: Pparam::default(),
            uesp_temp: Pparam::default(),
            delta_phase: 0.01, nlat1f: 100, nlat2f: 100, nlat1c: 50, nlat2c: 50,
            npole: false, nlatfill: 0, nlngfill: 0, lfudge: 0.0, llo: 0.0, lhi: 0.0,
            phase1: 0.05, phase2: 0.05, wavelength: 4700.0,
            roche1: false, roche2: true, eclipse1: true, eclipse2: true,
            glens1: false, use_radii: true, tperiod: 0.1,
            gdark_bolom1: false, gdark_bolom2: false,
            mucrit1: 0.0, mucrit2: 0.0,
            limb1: LDCType::Poly, limb2: LDCType::Poly,
            mirror: false, add_disc: false, nrad: 50, opaque: false,
            add_spot: false, nspot: 50, iscale: false,
        }
    }

    #[test]
    fn get_param_value_returns_correct_values() {
        let mdl = make_model();
        assert_eq!(mdl.get_param_value("q"), Some(0.5));
        assert_eq!(mdl.get_param_value("iangle"), Some(82.0));
        assert_eq!(mdl.get_param_value("t1"), Some(15000.0));
        assert_eq!(mdl.get_param_value("period"), Some(0.1));
    }

    #[test]
    fn get_param_value_unknown_returns_none() {
        let mdl = make_model();
        assert_eq!(mdl.get_param_value("nonexistent"), None);
        assert_eq!(mdl.get_param_value(""), None);
    }

    #[test]
    fn set_param_value_modifies_value() {
        let mut mdl = make_model();
        assert!(mdl.set_param_value("q", 0.42));
        assert_eq!(mdl.q.value, 0.42);
        // Other Pparam fields unchanged
        assert_eq!(mdl.q.range, 0.1);
        assert_eq!(mdl.q.dstep, 0.001);
    }

    #[test]
    fn set_param_value_unknown_returns_false() {
        let mut mdl = make_model();
        assert!(!mdl.set_param_value("nonexistent", 1.0));
    }

    #[test]
    fn set_then_get_round_trips() {
        let mut mdl = make_model();
        mdl.set_param_value("iangle", 85.5);
        assert_eq!(mdl.get_param_value("iangle"), Some(85.5));
    }
}
