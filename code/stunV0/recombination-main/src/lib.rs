use rand::distributions::WeightedIndex;
#[allow(dead_code)]
use rand::prelude::*;
use rand::seq::SliceRandom;

use serde::{Deserialize, Serialize};

#[derive(Serialize, Deserialize)]
pub struct Matrix<T> {
    data: Vec<T>,
    rows: usize,
    cols: usize,
}


impl<T: Default + Copy> Matrix<T> {
    pub fn new(rows: usize, cols: usize) -> Self {
        let mut data = Vec::<T>::new();
        data.resize_with(rows * cols, Default::default);
        Self { data, rows, cols }
    }
    #[inline]
    pub fn get_row(&self, i: usize) -> &[T] {
        &self.data[self.row(i)]
    }
    #[inline]
    pub fn get(&self, i: usize, j: usize) -> T {
        self.data[self.position(i, j)]
    }
    #[inline]
    pub fn set(&mut self, i: usize, j: usize, element: T) {
        let position = self.position(i, j);
        self.data[position] = element;
    }
    pub fn add_row(&mut self, row: &[T]) {
        if row.len() != self.cols {
            panic!(
                "row length ({}) is not equal to the number of columns ({})!",
                row.len(),
                self.cols
            );
        }
        self.data.extend_from_slice(row);
        self.rows += 1;
    }
    #[inline]
    pub fn extract_row(&mut self, i: usize) -> Vec<T> {
        self.data.drain(self.row(i)).collect()
    }
    #[inline]
    pub fn swap_remove_row(&mut self, i: usize) {
        for j in self.row(i).rev() {
            self.data.swap_remove(j);
        }
    }
    #[inline]
    pub fn row_iter(&self) -> std::slice::Chunks<T> {
        self.data.chunks(self.cols)
    }
    #[inline]
    fn row(&self, i: usize) -> std::ops::Range<usize> {
        (i * self.cols)..((i + 1) * self.cols)
    }
    #[inline]
    fn position(&self, i: usize, j: usize) -> usize {
        i * self.cols + j
    }
}

#[derive(Serialize, Deserialize)]
pub struct Indices {
    l: usize,
    seqs: Matrix<bool>,
    powers: Vec<usize>,
}

impl Indices {
    pub fn new(l: usize) -> Indices {
        let size = 2usize.pow(l as u32);
        let mut seqs = Matrix::<bool>::new(size, l);
        let mut powers = vec![0; l];

        for i in 0..l {
            powers[i] = 2usize.pow(i as u32);
        }

        for i in 0..size {
            let mut k = i;
            for j in 0..l {
                seqs.set(i, j, (k % 2) == 1);
                k = k / 2;
            }
        }

        let seqs = seqs;
        let powers = powers;

        Indices { l, seqs, powers }
    }
    pub fn to_index(&self, seq: &[bool]) -> usize {
        let mut idx: usize = 0;
        for i in 0..self.l {
            idx += seq[i] as usize * self.powers[i];
        }
        idx
    }
    pub fn to_sequence(&self, idx: usize) -> &[bool] {
        if idx > self.seqs.rows {
            panic!("index out of bounds {}. Maximum {}", idx, self.seqs.rows);
        }
        &self.seqs.get_row(idx)
    }
    pub fn get(&self, i: usize, j: usize) -> bool {
        self.seqs.get(i, j)
    }
}

use rand_distr::{Distribution, Normal};
use std::{
    io::{BufWriter, Write},
    path::Path,
};
use std::fs::File;
type BufferedFile = std::io::BufWriter<std::fs::File>;

// Rough Mount Fuji landscape
#[derive(Serialize, Deserialize)]
pub struct RoughMountFuji {
    l: usize,
    size: usize,
    mu_a: f64,
    sigma_a: f64,
    sigma_b: f64,
    additive_component: Vec<f64>,
    epistatic_component: Vec<f64>,
    landscape: Vec<f64>,
    indices: Indices,
    _maximum: Option<(usize, f64)>,
    _minimum: Option<(usize, f64)>,
}

impl RoughMountFuji {
    pub fn new(l: usize, mu_a: f64, sigma_a: f64, sigma_b: f64) -> RoughMountFuji {
        let size = 2usize.pow(l as u32);
        let indices = Indices::new(l);

        let mut rmf_landscape = RoughMountFuji {
            l,
            size,
            mu_a,
            sigma_a,
            sigma_b,
            additive_component: vec![0.; l],
            epistatic_component: vec![0.; size],
            landscape: vec![0.; size],
            indices,
            _maximum: None,
            _minimum: None,
        };
        rmf_landscape.generate();
        rmf_landscape
    }

    pub fn generate(&mut self) {
        let mut rng = thread_rng();
        let l = self.l;

        let normal_additive = Normal::new(self.mu_a, self.sigma_a).unwrap();
        let normal_epistatic = Normal::new(0., self.sigma_b).unwrap();
        for i in 0..l {
            self.additive_component[i] = normal_additive.sample(&mut rng);
        }
        for i in 0..self.size {
            let epistatic_component = normal_epistatic.sample(&mut rng);
            let mut additive_component = 0.;
            let seq = self.indices.to_sequence(i);
            for (i, &j) in seq.iter().enumerate() {
                additive_component += f64::from(u8::from(j)) * self.additive_component[i];
            }

            self.epistatic_component[i] = epistatic_component;
            self.landscape[i] = (epistatic_component + additive_component).exp();
        }
        self._maximum = Some(self.find_maximum());
        self._minimum = Some(self.find_minimum());
    }
    pub fn fitness(&self, seq: &[bool]) -> f64 {
        self.landscape[self.indices.to_index(seq)]
    }
    pub fn get_indices(&self) -> &Indices {
        &self.indices
    }
    pub fn maximum(&self) -> (usize, f64) {
        match self._maximum {
            Some(max) => max,
            None => self.find_maximum(),
        }
    }
    pub fn minimum(&self) -> (usize, f64) {
        match self._minimum {
            Some(min) => min,
            None => self.find_minimum(),
        }
    }
    fn find_maximum(&self) -> (usize, f64) {
        let mut max = 0;
        let mut maxf = f64::NEG_INFINITY;
        for (i, &f) in self.landscape.iter().enumerate() {
            if f > maxf {
                max = i;
                maxf = f;
            }
        }
        (max, maxf)
    }
    fn find_minimum(&self) -> (usize, f64) {
        let mut min = 0;
        let mut minf = f64::INFINITY;
        for (i, &f) in self.landscape.iter().enumerate() {
            if f < minf {
                min = i;
                minf = f;
            }
        }
        (min, minf)
    }
    pub fn save_to_file(&self, filename: &str) {
        let path = Path::new(filename);
        let display = path.display();
        let file = match File::create(&path) {
            Err(why) => panic!("couldn't create {}: {}", display, why),
            Ok(file) => file,
        };
        let mut writer = BufWriter::new(file);

        for seq in self.indices.seqs.row_iter() {
            let mut line = String::new();
            for &allele in seq.iter() {
                line = format!("{} {}", line, if allele { 1 } else { 0 });
            }
            line = format!("{} {}\n", line, self.fitness(seq));
            match writer.write(line.as_bytes()) {
                Ok(_) => {}
                Err(why) => {
                    panic!("couldn't write data to file: {}", why)
                }
            };
        }
    }
    pub fn save(&self, filename: &str) {
        let file = File::create(filename).unwrap();
        // serde_json::to_writer(&file, &self).unwrap();
        bincode::serialize_into(&file, &self).unwrap();
    }
    pub fn load(filename: &str) -> Self {
        let file = match File::open(filename) {
            Ok(f)  => f,
            Err(_) => panic!("File {} does not exist. Try to create it using 'generate_landscapes' utility.", filename)
        };
        // serde_json::from_reader(file).unwrap()
        bincode::deserialize_from(file).unwrap()
    }
}

#[derive(Clone, Copy)]
pub enum Allele {
    NotFixed,
    Fixed { time: usize, allele: bool },
}

pub struct Population<'a> {
    indices: Indices,
    pop: Vec<usize>,
    occupied_genotypes: Vec<usize>,
    recombination_map: Option<&'a [f64]>,
    l: usize,
    size: usize,
    fixation: Vec<Allele>,
}

pub enum InitialPopulation {
    NeutralSFS,
    NeutralSFSDriftThreshold(usize),
    Uniform(f64), // The value is the probability of derived allele
    UniformSymmetrized(f64),
}

impl<'a> Population<'a> {
    pub fn new(l: usize, size: usize, recombination_map: Option<&'a [f64]>) -> Population<'a> {
        Population {
            l,
            size,
            indices: Indices::new(l),
            pop: vec![0usize; 2usize.pow(l as u32)],
            occupied_genotypes: Vec::with_capacity(size),
            recombination_map,
            fixation: vec![Allele::NotFixed; l],
        }
    }
    pub fn without_recombination(l: usize, size: usize) -> Population<'a> {
        Population::new(l, size, None)
    }

    pub fn add_element(&mut self, i: usize) {
        if self.pop[i] == 0 {
            self.occupied_genotypes.push(i);
        }
        self.pop[i] += 1;
    }
    pub fn remove_element(&mut self, i: usize) {
        if self.pop[i] == 0 {
            panic!("trying to remove an inexistant element ({})!", i);
        }

        self.pop[i] -= 1;
        if self.pop[i] == 0 {
            let (j, _) = self
                .occupied_genotypes
                .iter()
                .enumerate()
                .find(|x| {
                    let (_, &v) = x;
                    v == i
                })
                .unwrap();
            self.occupied_genotypes.swap_remove(j);
        }
    }

    pub fn get_l(&self) -> usize {
        self.l
    }

    pub fn recombination(&mut self) {
        // Check if a recombination map exists
        let recombination_map = match &self.recombination_map {
            Some(map) => *map,
            None => {
                panic!("recombination map does not exist!")
            }
        };

        let mut rng = thread_rng();
        let mut new_population = vec![0usize; self.pop.len()];
        let mut new_occupied_genotypes = Vec::with_capacity(self.size);

        // Generates a list with all the individuals in a random order
        let mut list_of_individuals = Vec::<usize>::with_capacity(self.size);
        for &i in self.occupied_genotypes.iter() {
            for _ in 0..self.pop[i] {
                list_of_individuals.push(i);
            }
        }
        list_of_individuals.shuffle(&mut rng);

        for indices in list_of_individuals.chunks(2) {
            let mut idx1 = indices[0];
            let mut idx2 = indices[1];

            // if individuals are identical recombination skip trivial recombination
            if idx1 == idx2 {
                new_population[idx1] += 2;
                new_occupied_genotypes.push(idx1);
            } else {
                let recombine = (0..self.l).map(|i| -> bool {
                    let sample: f64 = (&mut rng).sample(rand::distributions::Open01);
                    sample < recombination_map[i]
                });

                let mut recombined = false;
                let mut genotype_1 = self.indices.to_sequence(idx1).to_vec();
                let mut genotype_2 = self.indices.to_sequence(idx2).to_vec();

                for (i, rec) in recombine.enumerate() {
                    if rec {
                        recombined = true;
                        let section1 = (&genotype_1[i..self.l]).to_vec();
                        let section2 = (&genotype_2[i..self.l]).to_vec();
                        (&mut genotype_1[i..self.l]).copy_from_slice(&section2[..]);
                        (&mut genotype_2[i..self.l]).copy_from_slice(&section1[..]);
                    }
                }

                if recombined {
                    idx1 = self.indices.to_index(&genotype_1);
                    idx2 = self.indices.to_index(&genotype_2);
                }
                new_population[idx1] += 1;
                new_population[idx2] += 1;
                new_occupied_genotypes.push(idx1);
                new_occupied_genotypes.push(idx2);
            }
        }

        //remove the repeated genotypes in the occupied_genotypes list
        new_occupied_genotypes.sort_unstable();
        new_occupied_genotypes.dedup();

        self.occupied_genotypes = new_occupied_genotypes;
        self.pop = new_population;
    }

    pub fn wright_fisher(&mut self, landscape: &RoughMountFuji) {
        let mut rng = thread_rng();

        let fitness: Vec<f64> = self
            .occupied_genotypes
            .iter()
            .map(|&i| -> f64 { landscape.landscape[i] * (self.pop[i] as f64) })
            .collect();

        // Zero the old population
        for &i in self.occupied_genotypes.iter() {
            self.pop[i] = 0;
        }

        // Create the new population
        let new_indices = rand::distributions::WeightedIndex::new(&fitness).unwrap();
        for _ in 0..self.size {
            let genotype = self.occupied_genotypes[new_indices.sample(&mut rng)];
            self.pop[genotype] += 1;
        }

        // Filter out the extinct genotypes
        self.occupied_genotypes = self
            .occupied_genotypes
            .iter()
            .filter(|&&i| self.pop[i] > 0)
            .map(|&i| i)
            .collect();
    }

    pub fn generate_initial_population(&mut self, initial_population: &InitialPopulation) {
        let mut rng = thread_rng();

        {   // clear previous population
            self.occupied_genotypes.clear();
            let size = self.pop.len();
            self.pop.clear();
            self.pop.resize(size, 0);
        }

        match initial_population {
            InitialPopulation::NeutralSFS => {
                panic!("InitialPopulation::NeutralSFS deprecated! Use InitialPopulation::NeutralSFSDriftThreshold(1) instead.")
            },
            InitialPopulation::NeutralSFSDriftThreshold(k) => {
                // Generates an initial population where the alleles are distributing according
                // to the site frequency spectrum of a neutrally evolving population
                let mut initial_population = Matrix::<bool>::new(self.size, self.l);

                let harmonic_number: f64 = (*k..=(self.size-*k)).map(|j| 1. / (j as f64)).sum();
                let probability: Vec<f64> = (*k..=(self.size-*k))
                    .map(|j| 1. / (j as f64 * harmonic_number))
                    .collect();
                let dist = WeightedIndex::new(&probability).unwrap();

                for j in 0..self.l {
                    let order = rand::seq::index::sample(&mut rng, self.size, self.size);
                    let count = dist.sample(&mut rng) + *k;

                    // let minor_allele = rng.gen_bool(0.5);
                    let minor_allele = if count < self.size/2 { true } else { false };

                    for (c, i) in order.into_iter().enumerate() {
                        if c < count {
                            initial_population.set(i, j, minor_allele);
                        } else {
                            initial_population.set(i, j, !minor_allele);
                        }
                    }
                }

                for genotype in initial_population.row_iter() {
                    let idx = self.indices.to_index(&genotype);

                    self.pop[idx] += 1;
                    self.occupied_genotypes.push(idx);
                }
            },
            InitialPopulation::Uniform(derived_allele_probability) => {
                // Generates an initial population where each individual has a probability equal to
                // derived_allele_probability of carrying the derived allele form for each allele
                // Minor alleles are encoded have value true
                for _ in 0..self.size {
                    let genotype: Vec<bool> = (0..self.l)
                        .map(|_| (&mut rng).gen::<f64>() < *derived_allele_probability)
                        .collect();

                    let idx = self.indices.to_index(&genotype);

                    self.pop[idx] += 1;
                    self.occupied_genotypes.push(idx);
                }

                // Orders the occupied genotypes and removes duplicates
                self.occupied_genotypes.sort_unstable();
                self.occupied_genotypes.dedup();

                // Recodes the alleles such that the minor allele always has value true
                // First checks if an allele should be recoded
                let mut recode = vec![false; self.l];
                for allele in 0..self.l {
                    let mut count = 0;
                    for &g_idx in &self.occupied_genotypes {
                        let g = self.indices.to_sequence(g_idx);
                        count += if g[allele] { self.pop[g_idx] } else { 0 };
                    }
                    recode[allele] = count > self.size/2;
                }

                // Then recodes all relevant alleles sequentially
                for (allele, &recode_allele) in recode.iter().enumerate() {
                    if recode_allele {
                        let mut new_pop = vec![0usize; 2usize.pow(self.l as u32)];
                        let mut new_occupied_genotypes = Vec::<usize>::new();

                        for &g_idx in &self.occupied_genotypes {
                            let mut g = Vec::<bool>::from(self.indices.to_sequence(g_idx));
                            g[allele] = !g[allele];

                            let new_idx = self.indices.to_index(&g);
                            new_pop[new_idx] = self.pop[g_idx];
                            new_occupied_genotypes.push(new_idx)
                        }

                        self.pop = new_pop;
                        self.occupied_genotypes = new_occupied_genotypes;
                    }
                }
            },
            InitialPopulation::UniformSymmetrized(derived_allele_probability) => {
                panic!("InitialPopulation::UniformSymmetrized deprecated!")
            },
            // InitialPopulation::Uniform(derived_allele_probability) => {
            //     // Generates an initial population where each individual has a probability equal to
            //     // derived_allele_probability of carrying the derived allele form for each allele
            //     for _ in 0..self.size {
            //         let genotype: Vec<bool> = (0..self.l)
            //             .map(|_| (&mut rng).gen::<f64>() < *derived_allele_probability)
            //             .collect();
            //         let idx = self.indices.to_index(&genotype);
            //
            //         self.pop[idx] += 1;
            //         self.occupied_genotypes.push(idx);
            //     }
            // },
            // InitialPopulation::NeutralSFS => {
            //     // Generates an initial population where the alleles are distributing according
            //     // to the site frequency spectrum of a neutrally evolving population
            //     let mut initial_population = Matrix::<bool>::new(self.size, self.l);
            //
            //     let harmonic_number: f64 = (1..self.size).map(|j| 1. / (j as f64)).sum();
            //     let probability: Vec<f64> = (1..self.size)
            //         .map(|j| 1. / (j as f64 * harmonic_number))
            //         .collect();
            //     let dist = WeightedIndex::new(&probability).unwrap();
            //
            //     for j in 0..self.l {
            //         let order = rand::seq::index::sample(&mut rng, self.size, self.size);
            //
            //         let count = dist.sample(&mut rng) + 1;
            //
            //         // let minor_allele = rng.gen_bool(0.5);
            //         let minor_allele = if count < self.size/2 { true } else { false };
            //
            //         for (c, i) in order.into_iter().enumerate() {
            //             if c < count {
            //                 initial_population.set(i, j, minor_allele);
            //             } else {
            //                 initial_population.set(i, j, !minor_allele);
            //             }
            //         }
            //     }
            //
            //     for genotype in initial_population.row_iter() {
            //         let idx = self.indices.to_index(&genotype);
            //
            //         self.pop[idx] += 1;
            //         self.occupied_genotypes.push(idx);
            //     }
            // },
            // InitialPopulation::UniformSymmetrized(derived_allele_probability) => {
            //     // Generates an initial population where each individual has a probability equal to
            //     // derived_allele_probability of carrying the derived allele form for each allele
            //     // This is a symmetrized version where for each locus the minor allele can be
            //     // either
            //     let minor_allele: Vec<bool> =
            //         (0..self.l).map(|_| (&mut rng).gen::<bool>()).collect();
            //
            //     for _ in 0..self.size {
            //         let genotype: Vec<bool> = (0..self.l)
            //             .map(|i| {
            //                 if minor_allele[i] {
            //                     (&mut rng).gen::<f64>() < *derived_allele_probability
            //                 } else {
            //                     (&mut rng).gen::<f64>() < 1. - *derived_allele_probability
            //                 }
            //             })
            //             .collect();
            //         let idx = self.indices.to_index(&genotype);
            //
            //         self.pop[idx] += 1;
            //         self.occupied_genotypes.push(idx);
            //     }
            // }
        }
        // Remove repeated genotypes
        self.occupied_genotypes.sort_unstable();
        self.occupied_genotypes.dedup();

        // Check the alleles that may be fixed by chance in the initial population
        self.fixation = vec![Allele::NotFixed; self.l];
        self.check_fixation(0);
    }

    fn write_to_file(file: &mut BufferedFile, contents: &str) {
        match file.write(contents.as_bytes()) {
            Ok(_) => {}
            Err(why) => { panic!("couldn't write data to file: {}", why) }
        };
    }

    pub fn save(&self, file: &mut BufferedFile, i: usize, recombination_rate: f64, c: usize) -> std::io::Result<()> {
        write!(file, "{}\t{}\t{}", i, recombination_rate, c)?;
        for &g in self.occupied_genotypes.iter() {
            write!(file, "\t")?;
            for &a in self.indices.to_sequence(g) {
                write!(file, "{}", if a { 1 } else { 0 })?;
            }
            write!(file, "\t{}", self.pop[g])?;
        }
        write!(file, "\n")?;
        Ok(())
    }

    pub fn check_fixation(&mut self, t: usize) {
        for i in 0..self.l {
            self.check_fixation_allele(t, i);
        }
    }

    pub fn check_fixation_allele(&mut self, t: usize, i: usize) {
        let initial_allele = self.indices.get(self.occupied_genotypes[0], i);

        let mut update = true;
        match self.fixation[i] {
            Allele::NotFixed => {
                // If marked as not fixed the allele needs to be checked
                for j in 1..self.occupied_genotypes.len() {
                    if self.indices.get(self.occupied_genotypes[j], i) != initial_allele {
                        // If a different version is found then the allele is not fixed and does
                        // not need to be updates
                        update = false;
                        break;
                    }
                } // If no element has a different version of the allele, then it needs to be
                  // marked as fixed
            }
            _ => {
                // Otherwise (if already marked as fixed) the allele does not need to be checked or
                // updated
                update = false;
            }
        }

        // fixation time is registered as well as the fixed version of the allele
        if update {
            self.fixation[i] = Allele::Fixed {
                time: t,
                allele: initial_allele,
            };
        }
    }

    pub fn get_fixation(&self, i: usize) -> Allele {
        self.fixation[i]
    }

    pub fn count_fixations(&self) -> usize {
        self.fixation
            .iter()
            .filter(|&x| match x {
                Allele::Fixed { time: _, allele: _ } => true,
                _ => false,
            })
            .count()
    }

    // pub fn count_fixations_filter(&self, filter: fn(&Allele)->bool) -> usize {
    //     self.fixation.iter().filter(filter).count()
    // }

    pub fn active_genotypes(&self) -> &[usize] {
        &self.occupied_genotypes
    }

    pub fn mean_fitness(&self, landscape: &RoughMountFuji) -> f64 {
        self.occupied_genotypes
            .iter()
            .map(|&i| -> f64 {
                let n = self.pop[i];
                (n as f64) * landscape.landscape[i]
            })
            .sum::<f64>()
            / (self.size as f64)
        // self.occupied_genotypes.iter().enumerate().map(|x| -> f64 {
        //     let (i, &n) = x;
        //     (n as f64)*landscape.landscape[i]
        // }).sum::<f64>()/(self.size as f64)
    }

    pub fn variance_fitness(&self, landscape: &RoughMountFuji, f: f64) -> f64 {
        let f2 = self.occupied_genotypes
            .iter()
            .map(|&i| -> f64 {
                let n = self.pop[i];
                (n as f64) * landscape.landscape[i] * landscape.landscape[i]
            })
            .sum::<f64>()
            / (self.size as f64);

        f2 - f * f
    }

    pub fn shannon_entropy(&self) -> f64 {
        self.occupied_genotypes
            .iter()
            .map(|&i| -> f64 {
                let n = self.pop[i];
                (n as f64) / (self.size as f64)
            })
            .map(|f| -> f64 {
                if f > 0. {
                    -f * f.ln()
                } else {
                    0.
                }
            })
            .sum()
    }

    pub fn maximum_fitness(&self, landscape: &RoughMountFuji) -> (usize, f64) {
        let (mut max, mut fmax) = (0usize, std::f64::NEG_INFINITY);
        for &g in self.occupied_genotypes.iter() {
            let fitness = landscape.landscape[g];
            if fitness > fmax {
                fmax = fitness;
                max = g;
            }
        }
        (max, fmax)
    }

    pub fn sequence(&self, idx: usize) -> &[bool] {
        self.indices.to_sequence(idx)
    }
}

impl std::ops::Index<usize> for Population<'_> {
    type Output = usize;

    fn index(&self, idx: usize) -> &Self::Output {
        &self.pop[idx]
    }
}
