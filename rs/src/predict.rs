use vcf::population::Population;
use vcf::predictor::{Predictor, PredictorFactory, PredictorJoaoFactory};
use vcf::stream::MetadataReader;

use utils::{process_variant, split};

use itertools::Itertools;
use std::collections::HashMap;

fn predict<S, PF>(mut stream: S, population: &Population, ratio: f32, factory: PF, min_dist: f32)
where
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

    let to_split = 0..individuals.len();
    let split_index = (to_split.len() as f32 * ratio) as usize;
    let (train_idxs, target_idxs) = split(to_split, split_index);

    let mut cumulative_labels: Vec<HashMap<String, f32>> = vec![HashMap::new(); individuals.len()];

    for polymorphism in unfolded.stream {
        let genotypes = polymorphism
            .genotypes()
            .map(|s| process_variant(s).unwrap_or(0.0))
            .collect::<Vec<_>>();
        let genotypes_and_labels = genotypes.iter().cloned().zip(labels.iter().cloned());
        let predictor = factory.build(genotypes_and_labels);

        for &target_idx in target_idxs.iter() {
            let genotype = genotypes[target_idx];
            let predictions = predictor.predict(genotype);
            for (prediction, confidence) in predictions {
                *cumulative_labels[target_idx]
                    .entry(prediction.to_string())
                    .or_default() += confidence;
            }
        }
    }

    let relative = 0.0 < min_dist && min_dist < 1.0;

    for target_idx in target_idxs {
        let predictions = cumulative_labels[target_idx]
            .iter()
            .sorted_by(|(_, f1), (_, f2)| f2.partial_cmp(f1).unwrap());

        let mut diff = (predictions[0].1 - predictions[1].1) as f32;
        if relative {
            diff /= predictions[0].1;
        }

        let final_prediction = if diff >= min_dist {
            predictions[0].0
        } else {
            "---"
        };

        println!(
            "{}\t{}\t{}",
            individuals[target_idx], labels[target_idx], final_prediction
        )
    }
}

pub fn main() -> ::std::io::Result<()> {
    use std::fs::File;
    use std::io::{stdin, BufReader};
    use vcf::stream::from_text_stream;

    let predictor_factory = PredictorJoaoFactory { threshold: 0.0 };

    let file = File::open("../populations/superpopulations.txt")?;
    let reader = BufReader::new(file);
    let population = Population::parse(reader)?;

    let stdin = stdin();
    let handle = stdin.lock();
    let stream = from_text_stream(handle);

    predict(stream, &population, 0.1, predictor_factory, 0.0);

    Ok(())
}
