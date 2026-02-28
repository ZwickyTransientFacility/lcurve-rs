use clap::Parser;
use lcurve::model::Model;
use lcurve::types::{Data, Datum, read_data};
use lcurve::orchestration::light_curve_comp;

#[derive(Parser)]
#[command(name = "lroche")]
#[command(about = "Compute light curve of a white dwarf / Roche-distorted star binary")]
struct Cli {
    /// Model parameter file
    model: String,

    /// Data file (use 'none' to generate regularly-spaced times)
    data: String,

    /// First time to compute (when data = 'none')
    #[arg(long, default_value_t = -0.2)]
    time1: f64,

    /// Last time to compute (when data = 'none')
    #[arg(long, default_value_t = 0.8)]
    time2: f64,

    /// Number of times to compute (when data = 'none')
    #[arg(long, default_value_t = 500)]
    ntime: usize,

    /// Exposure time (when data = 'none')
    #[arg(long, default_value_t = 0.001)]
    expose: f64,

    /// Number of exposure subdivisions (when data = 'none')
    #[arg(long, default_value_t = 1)]
    ndivide: i32,

    /// Output file for computed light curve
    #[arg(short, long)]
    output: Option<String>,

    /// Autoscale to minimise chi-squared (only with real data)
    #[arg(long)]
    scale: bool,
}

fn main() {
    if let Err(e) = run() {
        eprintln!("lroche: {}", e);
        std::process::exit(1);
    }
}

fn run() -> Result<(), Box<dyn std::error::Error>> {
    let cli = Cli::parse();

    let mdl = Model::from_file(&cli.model)?;

    let no_file = cli.data.to_lowercase() == "none";

    let data: Data = if no_file {
        let mut d = Vec::with_capacity(cli.ntime);
        for i in 0..cli.ntime {
            let time = cli.time1 + (cli.time2 - cli.time1) * i as f64 / (cli.ntime - 1) as f64;
            d.push(Datum {
                time,
                expose: cli.expose,
                flux: 0.0,
                ferr: 0.0,
                weight: 1.0,
                ndiv: cli.ndivide,
            });
        }
        d
    } else {
        read_data(&cli.data)?
    };

    if data.is_empty() {
        return Err("No data points".into());
    }

    let scale = cli.scale && !no_file;
    let rdata = !no_file;

    let result = light_curve_comp(&mdl, &data, scale, rdata)?;

    // Print results
    if rdata {
        println!("Weighted chi**2 = {:.12e}, wnok = {:.12e}", result.chisq, result.wnok);
        if mdl.iscale {
            println!("Scale factors = {:.12e}, {:.12e}, {:.12e}, {:.12e}, {:.12e}",
                     result.sfac[0], result.sfac[1], result.sfac[2],
                     result.sfac[3], result.sfac[4]);
        } else {
            println!("Scale factor = {:.12e}", result.sfac[0]);
        }
    }
    println!("White dwarf's contribution = {:.12e}", result.wdwarf);
    println!("log10(g1 [cgs]) = {:.12e}", result.logg1);
    println!("log10(g2 [cgs]) = {:.12e}", result.logg2);
    println!("Vol-averaged r1 = {:.12e}", result.rv1);
    println!("Vol-averaged r2 = {:.12e}", result.rv2);

    // Write output
    if let Some(outpath) = &cli.output {
        use std::io::Write;
        let mut f = std::fs::File::create(outpath)?;
        for (i, d) in data.iter().enumerate() {
            writeln!(f, "{:20.15e} {} {:17.10e} {} {} {}",
                     d.time, d.expose, result.calc[i], d.ferr, d.weight, d.ndiv)?;
        }
        println!("Written data to {}", outpath);
    } else {
        // Write to stdout
        for (i, d) in data.iter().enumerate() {
            println!("{:20.15e} {} {:17.10e} {} {} {}",
                     d.time, d.expose, result.calc[i], d.ferr, d.weight, d.ndiv);
        }
    }

    Ok(())
}
