use std::collections::HashMap;

use itertools::Itertools;

use utils::compute_frequencies;

pub trait PredictorFactory
{
    type Predictor: Predictor;
    
    fn build<I>(&self, genotypes_and_labels: I) -> Self::Predictor
    where
        I: Iterator<Item = (<Self::Predictor as Predictor>::Prediction, String)>;
}

pub trait Predictor {
    type Prediction;
    
    fn predict(&self, genotype: Self::Prediction) -> HashMap<&str, f32>;
}

pub struct PredictorJoao {
    threshold: f32,
    frequencies: HashMap<String, f32>,
}

impl Predictor for PredictorJoao {
    type Prediction = f32;
    
    fn predict(&self, genotype: f32) -> HashMap<&str, f32> {
        let distances = self
            .frequencies
            .iter()
            .map(|(label, freq)| (&**label, f32::abs(genotype - freq)))
            .sorted_by(|(_, f1), (_, f2)| f1.partial_cmp(f2).unwrap());

        if distances[1].1 - distances[0].1 >= self.threshold {
            vec![(distances[0].0, 1.0)].iter().cloned().collect()
        }
        else {
            vec![].iter().cloned().collect()
        }
    }
}

pub struct PredictorJoaoFactory {
    pub threshold: f32,
}

impl PredictorFactory for PredictorJoaoFactory {
    type Predictor = PredictorJoao;
    
    fn build<I>(&self, genotypes_and_labels: I) -> PredictorJoao
    where
        I: Iterator<Item = (f32, String)>,
    {
        PredictorJoao {
            threshold: self.threshold,
            frequencies: compute_frequencies(genotypes_and_labels),
        }
    }
}
