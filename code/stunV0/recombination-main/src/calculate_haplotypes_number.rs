pub mod lib;

use clap::{value_t, values_t, App, Arg, ArgMatches};
use lib::{InitialPopulation, Population};

fn main() {
    ///////////////////////////////////////////////////////////////////////////////////////////////
    // configure and read command line parameters
    let matches = get_command_line_matches();
    let p = Parameters::initialize(&matches);
    ///////////////////////////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////////////////////////
    for l in p.l {
        let mut population = Population::new(l, p.size, None);

        let mut pop_size = 0;
        let mut entropy = 0.;
        for _conf in 0..p.nconfs {
            population.generate_initial_population(&p.initial_population);

            pop_size += population.active_genotypes().len();
            entropy  += population.shannon_entropy();
        } //End of conf

        println!("{} {} {} {} {}", p.size, l, (pop_size as f64)/(p.nconfs as f64), (pop_size as f64)/(p.nconfs as f64 * 2f64.powf(l as f64)), (entropy as f64)/(p.nconfs as f64));
    }
}


struct Parameters {
    l: Vec<usize>,
    size: usize,
    initial_population: InitialPopulation,
    nconfs: usize
}

fn get_command_line_matches() -> ArgMatches<'static> {
    App::new("Haploid populations with standing variation")
        .version("0.1.0")
        .author("André Amado <andreamado1@gmail.com>")
        .about("Evolution in haploid populations with standing variation.")
        .group(clap::ArgGroup::with_name("model"))
        .arg(
            Arg::with_name("L")
                .short("L")
                .value_name("value")
                .help("List of L's (separated by comma)")
                .require_delimiter(true)
                .default_value("5"),
        )
        .arg(
            Arg::with_name("size")
                .short("N")
                .long("size")
                .value_name("value")
                .help("Population size")
                .default_value("2000"),
        )
        .arg(
            Arg::with_name("derived_allele_probability")
                .short("d")
                .long("derived_allele_probability")
                .value_name("value")
                .help("Derived allele probability in initial state")
                .default_value("0.5"),
        )
        .arg(
            Arg::with_name("neutralsfs")
            .long("neutralsfs")
        )
        .arg(
            Arg::with_name("nconfs")
                .short("c")
                .long("nconfs")
                .value_name("value")
                .help("Number of configurations")
                .default_value("100"),
        )
        .get_matches()
}


impl Parameters {
    fn initialize(matches: &ArgMatches) -> Parameters {
        // let l = value_t!(matches.value_of("L"), usize).unwrap();
        let l: Vec<usize> = values_t!(matches.values_of("L"), usize).unwrap();

        let size = value_t!(matches.value_of("size"), usize).unwrap();

        let initial_population = if matches.is_present("neutralsfs") {
            InitialPopulation::NeutralSFS
        } else {
            let derived_allele_probability =
                value_t!(matches.value_of("derived_allele_probability"), f64).unwrap();
            InitialPopulation::UniformSymmetrized(derived_allele_probability)
        };

        let nconfs = value_t!(matches.value_of("nconfs"), usize).unwrap();

        Parameters {
            l,
            size,
            initial_population,
            nconfs
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
