use vcf::population::Population;
use vcf::predictor::{Predictor, PredictorFactory, PredictorJoaoFactory};
use vcf::stream::{MetadataReader, VCFData};

use utils::{mean, process_variant, split, stdev};

// use itertools::Itertools;
use std::collections::HashSet;

fn compute_capacity<S, PF>(
    stream: S,
    population: &Population,
    repeats: usize,
    ratio: f32,
    factory: PF,
) where
    S: MetadataReader,
    PF: PredictorFactory,
    <PF as PredictorFactory>::Predictor: Predictor<Prediction = f32>,
{
    let unfolded = stream.unfold().unwrap();

    let individuals = unfolded.individuals;
    let labels = population
        .labels(individuals.iter())
        .map(|s| s.to_string())
        .collect::<Vec<_>>();

    let split_index = (individuals.len() as f32 * ratio) as usize;

    // Split once, since the call to `split` takes a lot of time to execute for
    // every polymorphism and every repeat. We call for the number of repeats
    // and apply the corresponding one, copying the same approach for all
    // polymorphisms.
    let train_target_indices: Vec<(HashSet<usize>, HashSet<usize>)>;

    if repeats > 1 {
        let train_target_indices_list = vec![split(0..individuals.len(), split_index); repeats];
        train_target_indices = train_target_indices_list
            .iter()
            .map(|(i, j)| (i.iter().cloned().collect(), j.iter().cloned().collect()))
            .collect();
    } else {
        train_target_indices = vec![(
            (0..individuals.len()).collect(),
            (0..individuals.len()).collect(),
        )];
    }

    let worker = Worker {
        repeats,
        individuals: &individuals,
        labels: &labels,
        train_target_indices: &train_target_indices,
        factory: &factory,
    };
    
    for polymorphism in unfolded.stream {
        let identifier = polymorphism.identifier().to_string();
        let (m, s) = worker.process_polymorphism(polymorphism);
        if repeats == 1 {
            println!("{}\t{}", identifier, m);
        } else {
            println!("{}\t{}\t{}", identifier, m, s)
        }
    }
}

struct Worker<'a, PF>
where
    PF: PredictorFactory + 'a,
    <PF as PredictorFactory>::Predictor: Predictor<Prediction = f32>,
{
    repeats: usize,
    individuals: &'a [String],
    labels: &'a [String],
    train_target_indices: &'a [(HashSet<usize>, HashSet<usize>)],
    factory: &'a PF,
}

unsafe impl<'a, PF> Send for Worker<'a, PF>
where
    PF: PredictorFactory,
    <PF as PredictorFactory>::Predictor: Predictor<Prediction = f32>,
{}

impl<'a, PF> Worker<'a, PF>
where
    PF: PredictorFactory,
    <PF as PredictorFactory>::Predictor: Predictor<Prediction = f32>,
{
    fn process_polymorphism(&self, polymorphism: VCFData) -> (f32, f32) {
        let mut accuracies = Vec::with_capacity(self.repeats);

        for repeat_idx in 0..self.repeats {
            let genotypes = polymorphism
                .genotypes()
                .map(|s| process_variant(s).unwrap_or(0.0))
                .collect::<Vec<_>>();
            let to_split = izip!(self.individuals, genotypes, self.labels);

            let mut train = Vec::new();
            let mut target = Vec::new();

            for (i, val) in to_split.enumerate() {
                if self.train_target_indices[repeat_idx].0.contains(&i) {
                    train.push(val);
                }
                if self.train_target_indices[repeat_idx].1.contains(&i) {
                    target.push(val);
                }
            }

            let mut genotypes_and_labels = Vec::new();

            for (_individual, genotype, label) in train {
                genotypes_and_labels.push((genotype, label.to_string()));
            }

            let predictor = self.factory.build(genotypes_and_labels.iter().cloned());
            let mut count = 0.0;
            let total = target.len() as f32;

            for (_individual, genotype, label) in target {
                let predictions = predictor.predict(genotype);
                if let Some(prob) = predictions.get::<str>(&label) {
                    count += prob;
                }
            }

            accuracies.push(count / total);
        }

        let m = mean(&accuracies);
        if self.repeats == 1 {
            // println!("{}\t{:.5}", polymorphism.identifier(), m);
            (m, 0.0)
        } else {
            let s = stdev(&accuracies);
            // println!("{}\t{:.5}\t{:.5}", polymorphism.identifier(), m, s);
            (m, s)
        }
    }
}

pub fn main() -> ::std::io::Result<()> {
    use std::fs::File;
    use std::io::{stdin, BufReader};
    use vcf::stream::from_text_stream;

    let predictor_factory = PredictorJoaoFactory { threshold: 0.0 };

    let population = {
        let file = File::open("../populations/superpopulations.txt")?;
        let reader = BufReader::new(file);
        Population::parse(reader)?
    };

    let stdin = stdin();
    let stream = {
        let handle = stdin.lock();
        from_text_stream(handle)
    };

    compute_capacity(stream, &population, 20, 0.1, predictor_factory);

    Ok(())
}
