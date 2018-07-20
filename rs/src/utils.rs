use std::collections::hash_map::Entry;
use std::collections::HashMap;
use std::hash::Hash;
use std::sync::Mutex;

use rand::{thread_rng, Rng};

/// Count the number of occurrences of each value in an iterator
pub fn histogram<K, I>(mut iter: I) -> HashMap<K, u32>
where
    K: Ord + Hash,
    I: Iterator<Item = K>,
{
    let mut counter = HashMap::new();

    for key in iter {
        *counter.entry(key).or_insert(0) += 1;
    }

    counter
}

pub fn compute_frequencies<L, I>(genotypes_and_labels: I) -> HashMap<L, f32>
where
    L: Eq + Hash + Clone,
    I: Iterator<Item = (f32, L)>,
{
    let mut counts: HashMap<L, f32> = HashMap::new();
    let mut totals: HashMap<L, u32> = HashMap::new();

    for (genotype, label) in genotypes_and_labels {
        *counts.entry(label.clone()).or_default() += genotype;
        *totals.entry(label.clone()).or_default() += 1u32;
    }

    counts
        .iter()
        .map(|(label, count)| (label.clone(), count / totals[label] as f32))
        .collect()
}

lazy_static! {
    static ref PROCESS_VARIANT_CACHE: Mutex<HashMap<String, Option<f32>>> =
        Mutex::new(HashMap::new());
}

pub fn process_variant(variant: &str) -> Option<f32> {
    let mut lock = PROCESS_VARIANT_CACHE.lock().unwrap();

    let entry = match lock.entry(variant.to_string()) {
        Entry::Occupied(entry) => return *entry.get(),
        Entry::Vacant(entry) => entry,
    };

    let genotype = variant.split(':').nth(0)?;

    let alleles = if genotype.find('|').is_some() {
        genotype.split('|')
    } else {
        genotype.split('/')
    }.collect::<Vec<_>>();

    let result = Some(alleles.iter().filter(|i| **i != "0").count() as f32 / alleles.len() as f32);
    entry.insert(result);
    result
}

pub fn split<I, T>(iter: I, split_index: usize) -> (Vec<T>, Vec<T>)
where
    I: Iterator<Item = T>,
    T: Clone,
{
    let mut vec: Vec<_> = iter.collect();
    let slice = &mut vec;
    thread_rng().shuffle(slice);
    (
        slice[split_index..].iter().cloned().collect(),
        slice[..split_index].iter().cloned().collect(),
    )
}

pub fn mean(xs: &[f32]) -> f32 {
    let total = xs.iter().map(|_| 1.0).sum::<f32>();
    xs.iter().sum::<f32>() / total
}

pub fn stdev(xs: &[f32]) -> f32 {
    let m = mean(xs);
    let num = xs.iter().map(|x| (*x - m).powf(2.0)).sum::<f32>();
    let den = xs.len() as f32;
    (num / den).powf(0.5)
}
