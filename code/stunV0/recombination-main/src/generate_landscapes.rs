pub mod lib;

use clap::{value_t, App, Arg};
use lib::RoughMountFuji;
use std::path::Path;
use glob::glob;

fn main() {
    let matches = App::new("Generate Rough Mount landscapes")
        .version("0.1.0")
        .author("André Amado <andreamado1@gmail.com>")
        .about("Generate Rough Mount Fuji fitness landscapes.")
        .group(clap::ArgGroup::with_name("model"))
        .arg(
            Arg::with_name("L")
                .short("L")
                .value_name("value")
                .help("Number of alleles")
                .default_value("5"),
        )
        .arg(
            Arg::with_name("mu_a")
                .short("m")
                .long("mua")
                .value_name("value")
                .help("mean of the additive component (μa)")
                .default_value("0.1"),
        )
        .arg(
            Arg::with_name("sigma_a")
                .short("S")
                .long("sigmaa")
                .value_name("value")
                .help("Standard variation of the additive component (σa)")
                .default_value("0."),
        )
        .arg(
            Arg::with_name("sigma_b")
                .short("s")
                .long("sigmab")
                .value_name("value")
                .help("Standard variation of the epistatic component (σb)")
                .default_value("0.5"),
        )
        .arg(
            Arg::with_name("nlandscapes")
                .long("landscapes")
                .value_name("value")
                .help("Number of landscapes to sample")
                .default_value("1"),
        )
        .arg(
            Arg::with_name("identifier")
                .short("i")
                .long("identifier")
                .value_name("value")
                .help("Run identifier")
                .default_value(""),
        )
        .get_matches();

    let l = value_t!(matches.value_of("L"), usize).unwrap();

    let mu_a = value_t!(matches.value_of("mu_a"), f64).unwrap();
    let sigma_a = value_t!(matches.value_of("sigma_a"), f64).unwrap();
    let sigma_b = value_t!(matches.value_of("sigma_b"), f64).unwrap();

    let sample_landscapes = value_t!(matches.value_of("nlandscapes"), usize).unwrap();

    let identifier = matches.value_of("identifier").unwrap();

    let generic_filename = format!(
        "landscapes/rmf_landscape_L{}_mua{}_sa{}_sb{}_*{}.bin.fl",
        l, mu_a, sigma_a, sigma_b, identifier
    );
    let existing_landscapes: usize = glob(&generic_filename).unwrap().count();

    for i in 0..sample_landscapes {
        let filename = format!(
            "landscapes/rmf_landscape_L{}_mua{}_sa{}_sb{}_{}{}.dat",
            l, mu_a, sigma_a, sigma_b, i + existing_landscapes, identifier
        );
        if Path::new(&filename).exists() {
            println!("Landscape '{}' already exists!", filename);
        } else {
            let fitness_landscape = RoughMountFuji::new(l, mu_a, sigma_a, sigma_b);
            fitness_landscape.save_to_file(&filename);

            let filename = format!(
                "landscapes/rmf_landscape_L{}_mua{}_sa{}_sb{}_{}{}.bin.fl",
                l, mu_a, sigma_a, sigma_b, i + existing_landscapes, identifier
            );
            fitness_landscape.save(&filename);
        }
    }
}
