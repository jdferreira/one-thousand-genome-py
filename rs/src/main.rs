#[macro_use]
extern crate lazy_static;
#[macro_use]
extern crate itertools;
extern crate rand;

mod capacity;
mod predict;
mod utils;
mod vcf;

use vcf::filter::{IndividualsFilter, Pipe, PolymorphismFilter};
use vcf::stream::MetadataReader;

use std::io;

fn _main_1() -> Result<(), &'static str> {
    let stdin = io::stdin();
    let handle = stdin.lock();
    let stream = vcf::stream::from_text_stream(handle);

    let f1 = IndividualsFilter::new(vec!["HG00096", "HG00097", "HG00146"]);
    let f2 = PolymorphismFilter::new(vec![
        "rs370363256",
        "rs61869635",
        "rs371609562",
        "rs558785434",
        "rs575614128",
        "rs566764841",
        "rs200634578",
        "rs560955407",
        "rs537964411",
    ]);

    let stream = Pipe::new(stream).pipe(f1).pipe(f2).stream();
    let unfolded = stream.unfold()?;

    for polymorphism in unfolded.stream {
        polymorphism.debug();
        // println!(
        //     "{}\t{}\t{}\t{}",
        //     data.identifier(),
        //     data.genotype(0),
        //     data.genotype(1),
        //     data.genotype(2)
        // );
    }

    Ok(())
}

fn main() -> ::std::io::Result<()> {
    ::capacity::main()
}
