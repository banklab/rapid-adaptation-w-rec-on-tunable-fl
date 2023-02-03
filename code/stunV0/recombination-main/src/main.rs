pub mod lib;

use clap::{value_t, values_t, App, Arg, ArgMatches};
use lib::{InitialPopulation, Population, RoughMountFuji};
use std::io::Write;

type BufferedFile = std::io::BufWriter<std::fs::File>;

macro_rules! build_filename {
    ($filename_structure: expr, $p: expr) => {
        format!($filename_structure, $p.l, $p.size, $p.mu_a, $p.sigma_a, $p.sigma_b, $p.identifier)
    }
}


fn main() {
    ///////////////////////////////////////////////////////////////////////////////////////////////
    // configure and read command line parameters
    let matches = get_command_line_matches();
    let p = Parameters::initialize(&matches);
    ///////////////////////////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////////////////////////
    // open files to save data
    let mut writer = open_file(build_filename!("data/data/rmf_L{}_N{}_mua{}_sa{}_sb{}_{}.dat", p));
    write_to_file(&mut writer, "#landscape_id\trecombination_rate\tconf\ttfixation\tfound_maximum\tgenetic_load\tfinal_fitness\n");

    let mut writer_fixations = open_file(build_filename!("data/fixations/rmf_L{}_N{}_mua{}_sa{}_sb{}_{}.dat", p));
    write_to_file(&mut writer_fixations, "#landscape_id\trecombination_rate\tconf\tfixation_time allele\n");

    let mut writer_details = open_file(build_filename!("data/population_details/rmf_L{}_N{}_mua{}_sa{}_sb{}_{}.dat", p));
    write_to_file(&mut writer_details, "#landscape_id\trecombination_rate\tconf\tt\tmean_fitness\tvariance_fitness\tmax_fitness\tpopulation_entropy\tn_fixations\tn_genotypes\n");

    let mut writer_initial_population = open_file(build_filename!("data/initial_population/rmf_L{}_N{}_mua{}_sa{}_sb{}_{}.dat", p));
    write_to_file(&mut writer_initial_population, "#landscape_id\trecombination_rate\tconf\tgenotype\tcount\t...\n");
    ///////////////////////////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////////////////////////
    // Main loop
    // repeat for the desired number of fitness landscapes
    for i in 0..p.sample_landscapes {
        // load the fitness landscape from a file
        // These fitness landscapes are generated with 'generate_landscapes'
        let landscape_filename = format!(
            "landscapes/rmf_landscape_L{}_mua{}_sa{}_sb{}_{}.bin.fl",
            p.l, p.mu_a, p.sigma_a, p.sigma_b, i
        );
        let fitness_landscape = RoughMountFuji::load(&landscape_filename);

        // run the populations for all the recombination rates listed by the user
        for &recombination_rate in p.recombination_rate_list.iter() {
            let recombination_map = vec![recombination_rate; p.l];

            // repeat for nconfs runs
            let mut population = Population::new(p.l, p.size, Some(&recombination_map));
            for conf in 0..p.nconfs {
                population.generate_initial_population(&p.initial_population);
                population.save(&mut writer_initial_population, i, recombination_rate, conf);

                let mut t = 0;
                let final_genotype = loop {
                    ///////////////////////////////////////////////////////////////////////////////
                    // update the population to the new generation
                    // Biology is here!
                    population.recombination();
                    population.wright_fisher(&fitness_landscape);
                    ///////////////////////////////////////////////////////////////////////////////

                    // check for genotype fixations in the population
                    population.check_fixation(t);
                    t += 1;

                    // stop if only one genotype remains
                    let active_genotypes = population.active_genotypes();
                    if active_genotypes.len() == 1 {
                        break active_genotypes[0];
                    }

                    ///////////////////////////////////////////////////////////////////////////////
                    // write results to a file
                    write_details_to_file(&mut writer_details, i, recombination_rate, conf, t, &population, &fitness_landscape);
                    ///////////////////////////////////////////////////////////////////////////////
                };

                ///////////////////////////////////////////////////////////////////////////////////
                // write some more results to data files
                write_final_data_to_file(&mut writer, i, conf, recombination_rate, t, final_genotype, &fitness_landscape);
                write_fixations_to_file(&mut writer_fixations, i, recombination_rate, &population, conf);
                ///////////////////////////////////////////////////////////////////////////////////
            } //End of conf
        } //End of recombination rate
    }
}




///////////////////////////////////////////////////////////////////////////////////////////////////
// auxiliary functions and data structures to write to files and read command line parameters
// unimportant to the Biology
fn open_file(filename: String) -> BufferedFile {
    let path = std::path::Path::new(&filename);
    let file = match std::fs::File::create(&path) {
        Err(why) => panic!("couldn't create {}: {}", path.display(), why),
        Ok(file) => file,
    };
    std::io::BufWriter::new(file)
}

fn write_to_file(file: &mut BufferedFile, contents: &str) {
    match file.write(contents.as_bytes()) {
        Ok(_) => {}
        Err(why) => { panic!("couldn't write data to file: {}", why) }
    };
}

fn write_details_to_file(file: &mut BufferedFile, i: usize, recombination_rate: f64, c: usize, t: usize, population: &Population, fitness_landscape: &RoughMountFuji) {
    let (_, maxf) = population.maximum_fitness(&fitness_landscape);
    let mean_fitness = population.mean_fitness(&fitness_landscape);
    let data = format!(
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
        i,
        recombination_rate,
        c,
        t,
        mean_fitness,
        population.variance_fitness(&fitness_landscape, mean_fitness),
        maxf,
        population.shannon_entropy(),
        population.count_fixations(),
        population.active_genotypes().len()
    );
    write_to_file(file, &data);
}

fn write_fixations_to_file(file: &mut BufferedFile, i: usize, recombination_rate: f64, population: &Population, conf: usize) {
    let mut fixations_line = format!("{}\t{}\t{}", i, recombination_rate, conf);
    for j in 0..population.get_l() {
        let allele_info = match population.get_fixation(j) {
            lib::Allele::NotFixed => "-1 -1".to_string(),
            lib::Allele::Fixed { time, allele } => {
                format!("{} {}", time, if allele { 1 } else { 0 })
            }
        };
        fixations_line = format!("{}\t{}", fixations_line, allele_info);
    }
    fixations_line = format!("{}\n", fixations_line);
    write_to_file(file, &fixations_line);
}

fn write_final_data_to_file(file: &mut BufferedFile, i: usize, conf: usize, recombination_rate: f64, t: usize, final_genotype: usize, fitness_landscape: &RoughMountFuji) {
    let max_sequence = fitness_landscape.get_indices().to_sequence(final_genotype);
    let (max, maxf) = fitness_landscape.maximum();

    let found_maximum = if final_genotype == max { 1 } else { 0 };
    let final_fitness = fitness_landscape.fitness(max_sequence);

    let genetic_load = (maxf - final_fitness) / maxf;

    let data_line = format!(
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
        i, recombination_rate, conf, t, found_maximum, genetic_load, final_fitness
    );
    write_to_file(file, &data_line);
}

struct Parameters<'a> {
    l: usize,
    size: usize,
    initial_population: InitialPopulation,
    recombination_rate_list: Vec<f64>,
    mu_a: f64,
    sigma_a: f64,
    sigma_b: f64,
    nconfs: usize,
    sample_landscapes: usize,
    identifier: &'a str
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
                .help("Number of alleles")
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
                .conflicts_with("neutralsfs"),
        )
        .arg(Arg::with_name("neutralsfs")
                .long("neutralsfs")
                .value_name("drift_threshold")
                .min_values(0)
                .max_values(1)
                .help("Neutral SFS distribution with a genetic drift threshold")
                .conflicts_with("derived_allele_probability"),
            )
        .arg(
            Arg::with_name("recombination_rates")
                .short("r")
                .long("recombination_rates")
                .value_name("value")
                .help("List of recombination rates (separated by comma)")
                .require_delimiter(true)
                .default_value("0.01"),
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
            Arg::with_name("nconfs")
                .short("c")
                .long("nconfs")
                .value_name("value")
                .help("Number of configurations")
                .default_value("100"),
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
        .get_matches()
}


impl<'a> Parameters<'a> {
    fn initialize(matches: &'a ArgMatches) -> Parameters<'a> {
        let l = value_t!(matches.value_of("L"), usize).unwrap();
        let size = value_t!(matches.value_of("size"), usize).unwrap();

        let initial_population = if matches.is_present("neutralsfs") {
            let values: Vec<usize> = values_t!(matches.values_of("neutralsfs"), usize).unwrap();
            let k = if values.len() == 1 {
                values[0]
            } else {
                1
            };
            InitialPopulation::NeutralSFSDriftThreshold(k)
        } else {
            let derived_allele_probability =
                value_t!(matches.value_of("derived_allele_probability"), f64).unwrap_or(0.5);
            InitialPopulation::Uniform(derived_allele_probability)
        };

        let recombination_rate_list: Vec<f64> =
            values_t!(matches.values_of("recombination_rates"), f64).unwrap();

        let mu_a = value_t!(matches.value_of("mu_a"), f64).unwrap();
        let sigma_a = value_t!(matches.value_of("sigma_a"), f64).unwrap();
        let sigma_b = value_t!(matches.value_of("sigma_b"), f64).unwrap();

        let nconfs = value_t!(matches.value_of("nconfs"), usize).unwrap();
        let sample_landscapes = value_t!(matches.value_of("nlandscapes"), usize).unwrap();

        let identifier = matches.value_of("identifier").unwrap();

        Parameters {
            l,
            size,
            initial_population,
            recombination_rate_list,
            mu_a,
            sigma_a,
            sigma_b,
            nconfs,
            sample_landscapes,
            identifier
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
