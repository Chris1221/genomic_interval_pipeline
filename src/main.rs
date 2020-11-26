use std::collections::HashMap;
use std::fs::File;
use std::io::BufReader;
use std::io::BufRead;
use log::{info, debug, warn};
use std::cmp;
use env_logger::Env;
use std::vec::Vec;
use ndarray::prelude::*;
use std::path::PathBuf;
use structopt::StructOpt;
use rust_htslib::faidx;
use pbr::ProgressBar;
use flate2::read::GzDecoder;

extern crate ndarray;
extern crate theban_interval_tree;
extern crate structopt;
extern crate flate2;
extern crate hdf5;


use memrange::Range;

/// Genomic interval preprocessing for keras.
#[derive(StructOpt, Debug, Clone)]
#[structopt(name = "scJaccard")]
struct Opt {

    /// Newline seperated list of bed files to process.
    #[structopt(short, long, default_value = "data/metadata.txt")]
    input: PathBuf,

    /// Reference sequence 
    #[structopt(short, long, parse(from_os_str), default_value = "data/blah.fq")]
    fastq: PathBuf,
   
    /// Size of the regions to report.
    #[structopt(long, default_value = "600")]
    length: u64,

    /// Log level. Defaults to Info (useful information and statistics). 
    #[structopt(long, default_value = "info")]
    loglevel: String,
}

fn main() -> std::io::Result<()> {
    // 1. Read in the meta data file. 
    //
    //      This file is simple a one file per line description of all
    //      bed files that are to be included in the database. Create a hash map of
    //      the file names, with each given a numeric index.
    //
    let opt = Opt::from_args();
    
    let env = Env::default()
        .filter_or("MY_LOG_LEVEL", opt.loglevel)
        .write_style_or("MY_LOG_STYLE", "always");
    env_logger::init_from_env(env);

    // Set up labels
    let mut training = Vec::new();
    for i in 1..19 { training.push(i) }
    let mut test = Vec::new();
    for i in 19..21 { test.push(i) } 
    let mut validation = Vec::new();
    for i in 21..23 { validation.push(i) }

    let mut valid_chromosomes = Vec::new();
    for i in 1..23 { valid_chromosomes.push(format!("{}{}", "chr", i)) }

    // Hashmap to hold the bed indices
    let mut metadata = HashMap::new();

    
    // Read in the file paths to the bed files.
    info!("Reading metadata file.");

    let mut i = 1;
    {
        let metadata_file = File::open(opt.input)?;
        let metadata_file_reader = BufReader::new(metadata_file);
        
        for line in metadata_file_reader.lines(){
            //debug!("Adding a bed file named {}", line.unwrap());
            let line_string = line.unwrap();
            debug!("{}", line_string);
            metadata.insert(line_string, i);
            i += 1;
            debug!("{}", i);
        }
    }

    let number_of_labels = i-1;

    // 2. Read in the files and create an interval tree, per chromosome. 
    //
    // For each bed file, 
    //      For each interval in the bed file
    //          If the interval is not in the tree
    //              a.  Add it to the tree
    //              b.  Add a record to the key hashtable, indexed by hashing the
    //                  range. This will contain a vector of the numeric indices 
    //                  of the bed files from which they derive.
    //          If the interval is already in the tree.
    //              a.  Retrieve the existing range
    //              b.  Delete the range from the tree, and from the the 
    //                  hash table, grabbing the vector.
    //              c.  Merge the ranges together.
    //              d.  Add the merged range to the interval tree
    //                  and the hash table with the current index 
    //                  appended.
    

    // Hash map of interval trees by chromosome. 
    
    let mut interval_trees = HashMap::new();
    let mut regions_by_chr = HashMap::new();

    // Make an entry for all the predefined chromsomes. We don't
    // want all of the alternate assemblies.
    for n in 1..23 {
        let idx = format!("{}{}", "chr", n);
        interval_trees.insert(idx.clone(), theban_interval_tree::IntervalTree::<i32>::new());
        let regions: HashMap<Range, Vec<u64>> = HashMap::new();
        regions_by_chr.insert(idx.clone(), regions);
    }

    //let mut regions: HashMap<Range, Vec<u64>>


    for (bed_file, bed_idx) in metadata {
        info!("\tAdding regions from {}...", bed_file);
        let bed = GzDecoder::new(File::open(bed_file)?);
        //let bed = File::open(bed_file)?;
        let bed_reader = BufReader::new(bed);


        for line in bed_reader.lines(){
            let _line = line?;
            let vec = _line.split("\t").collect::<Vec<&str>>();
            let chr = vec[0].to_string();
            let idx_vec = vec![bed_idx];
            let range = Range::new(vec[1].parse::<u64>().unwrap(), vec[2].parse::<u64>().unwrap());

            // If not a numeric chromosome, skip for now.
            if !valid_chromosomes.iter().any(|k| k == &chr) {
                continue;
            }

            if interval_trees[&chr].contains(range){
                let mut to_delete: Vec<Range> = Vec::new();
                for (_i, stored) in interval_trees[&chr]
                                            .range(range.min, range.max)
                                                .enumerate(){
                    debug!("Region intersection, resolving...");
                    if (range.min == stored.0.min) & (range.max == stored.0.max) {
                        // If the range is the same, then just 
                        // add the index.
                        regions_by_chr
                            .get_mut(&chr)
                            .unwrap()
                            .get_mut(&stored.0)
                            .unwrap()
                            .push(bed_idx);
                    } else {
                        // If not, then resolve the intersection.
                        //
                        // Get the current value for that region
                        let mut v = regions_by_chr
                            .get_mut(&chr)
                            .unwrap()
                            .get_mut(&stored.0)
                            .unwrap()
                            .clone();
                        regions_by_chr
                            .get_mut(&chr)
                            .unwrap()
                            .remove(&stored.0);
                        
                        // Create a new region based on both
                        let new = Range::new(
                            cmp::min(range.min, stored.0.min), 
                            cmp::max(range.max, stored.0.max));

                        debug!("\tCreated a new region {:?}", new);

                        // Add back in the new vector of bed files
                        // for the range.
                        v.push(bed_idx);
                        regions_by_chr
                            .get_mut(&chr)
                            .unwrap()
                            .insert(new, v.to_vec());

                        // Add this to the list of regions to be deleted
                        to_delete.push(stored.0) 
                    }
                }
                for del in to_delete.iter(){
                    interval_trees.get_mut(&chr).unwrap().delete(*del);
                    debug!("Removed region {:?}", del)
                }
            } else {
                // Bed_idx isn't used here, just a placeholder.
                // The real indexing is done in the regions hashmap.
                debug!("Added region {:?}", range);
                interval_trees.get_mut(&chr).unwrap().insert(range, 1);
                regions_by_chr
                    .get_mut(&chr)
                    .unwrap()
                    .insert(range, idx_vec);
            }


        }
    

    }

    info!("\tDone.");

    debug!("Regions now looks like: {:?}", regions_by_chr);

    // Need to post-process the regions now
    //  1.  Center them to n bp long (600)
    //  2.  Merge any that overlap (or not... can probably skip this for now 
    //      with minimal impact.)
    //
    //  So basically just loop over all the regions and apply a center 
    //  function to them, doing the same thing as above and replacing 
    //  the regions with the new ones.
    //
    //  Alternatievly could create a second set.
    //
    let mut qc_interval_trees = HashMap::new();
    let mut qc_regions_by_chr = HashMap::new();

    for n in 1..23 {
        let idx = format!("{}{}", "chr", n);
        qc_interval_trees.insert(
            idx.clone(), 
            theban_interval_tree::IntervalTree::<i32>::new());
        let regions: HashMap<Range, Vec<u64>> = HashMap::new();
        qc_regions_by_chr.insert(idx.clone(), regions);
    }

    info!("Correcting regions to the correct length.");
    let mut pb = ProgressBar::new(22);

    let mut counter = Array::zeros((22,2));

    for n in 1..23 {
        let mut i = 0;
        let idx = format!("{}{}", "chr", n);
        for (_i, stored) in interval_trees.get_mut(&idx.clone()).unwrap().iter().enumerate() {
            let corrected_range = center_region(&stored.0, opt.length);
            debug!("Corrected region {:?} to {:?}.", stored.0, corrected_range);
            // Adding range
            qc_interval_trees
                .get_mut(&idx.clone())
                .unwrap()
                .insert(corrected_range, 1);
            // Adding annotation
            qc_regions_by_chr
                .get_mut(&idx.clone())
                .unwrap()
                .insert(
                    corrected_range, 
                    regions_by_chr[&idx.clone()][&stored.0].clone());
            i += 1;
        }
        counter[[n-1, 0]] = n;
        counter[[n-1, 1]] = i;
        pb.inc();
    }
    info!("\tDone.");
    debug!("{:?}", qc_regions_by_chr);

    // Retrieve the sequence of the ranges. One hot encode them.
    
    let file = hdf5::File::create("test.h5").unwrap();
    
    //let mut seqs = HashMap::new();
    //let mut labels = HashMap::new();
    let mut index = HashMap::new();
    let mut j = HashMap::new();
    let mut i = HashMap::new();
    let mut seqs_db = HashMap::new();
    let mut labels_db = HashMap::new();
     
    j.insert("training", get_total(&counter, &training));
    j.insert("test", get_total(&counter, &test));
    j.insert("validation", get_total(&counter, &validation));

    for label in ["training", "test", "validation"].iter() {
        //seqs.insert(label, Array::zeros(( j[label] as usize, 4 as usize, opt.length as usize)));
        //labels.insert(label,  Array::zeros(( number_of_labels as usize, j[label] as usize)));
        index.insert(label, Vec::new());
        i.insert(label,  1);


        seqs_db.insert(label, file.new_dataset::<u64>()
                       .resizable(true)
                       .create(&format!("{}_{}", label, "seqs"), (1, 4, 600) )
                       .unwrap());
        labels_db.insert(label, file.new_dataset::<u64>()
                       .resizable(true)
                       .create(&format!("{}_{}", label, "labels"), (number_of_labels as usize, 1) )
                       .unwrap());

    }
    info!("One hot encoding sequences and writing to HDF5.");

    let mut pb = ProgressBar::new(22);

    let fa = faidx::Reader::from_path(opt.fastq).unwrap();
    for n in 1..23 {

        let dataset: &str = which_dataset(n, &training, &test, &validation);

        let chr_string = format!("{}{}", "chr", n);

        for (region, annotation) in qc_regions_by_chr[&chr_string].iter(){
            let region_string = format!("{}:{}-{}", &chr_string, region.min, region.max);
            let this_i: usize = *i.get(&dataset).unwrap();

            // Get a reference to the writer function of the 
            // database.
            let writer = seqs_db.get(&dataset).unwrap().as_writer();
            let label_writer = labels_db.get(&dataset).unwrap().as_writer();

            // Expand the dimensionality of the database to ensure there's enoug
            // room for the incoming data.
            seqs_db
                .get(&dataset)
                .unwrap()
                .resize((this_i as usize, 4, opt.length as usize));
            labels_db
                .get(&dataset)
                .unwrap()
                .resize((number_of_labels as usize, this_i as usize));        
            
            index.get_mut(&dataset)
                .unwrap()
                .push(region_string.clone());

            // Fetch sequence from fasta
            let seq = fa.fetch_seq_string(&chr_string, region.min as usize, (region.max -1) as usize).unwrap();
            // One hot encode
            let out = one_hot_encode_seq(&seq, 600);
            // One hot encode label vector 
            let label = one_hot_encode_labels(annotation, number_of_labels);
           
            // Write to resized dataset 
            writer.write_slice(&out, s![this_i-1, .., ..]);
            label_writer.write_slice(&label, s![.., this_i-1]);

            
            // One hot encode the labels and enter them into the array
            //labels.get_mut(&dataset)
             //   .unwrap()
              //  .slice_mut(s![.., this_i]).assign(&label);

            *i.get_mut(&dataset).unwrap() += 1;
            
        }
        pb.inc();
    }

    Ok(())
}


fn center_region(region: &memrange::Range, length: u64 ) -> Range {
    let midpoint = region.min + (region.min + region.max) / 2;
    if (length / 2) > midpoint {
        return Range::new(0, length);
    } else  {   
        return Range::new(midpoint - length / 2, midpoint + length / 2);
    }
}

fn one_hot_encode_seq(seq: &String, length: u64) -> ndarray::ArrayBase<ndarray::OwnedRepr<u64>, ndarray::Dim<[usize; 2]>> {
    let mut out = Array::zeros((4, length as usize));
    for (i, c) in seq.chars().enumerate(){
        let mut code = 10;
        if (c == 'a') || (c == 'A'){
            code = 0;
        } else if (c == 'c') || (c == 'C'){
            code = 1;
        } else if (c == 't') || (c == 'T'){
            code = 2;
        } else if (c == 'g') || (c == 'G'){
            code = 3;
        } else if (c == 'n') || (c == 'N'){
            code = 4
        } else {
            warn!("Bad base encountered. Skipping.")
        }

        if code <= 3 {
            out[[code, i]] = 1
        }
    }
    return out;
}

fn one_hot_encode_labels(labels: &Vec::<u64>, length: u64) -> ndarray::ArrayBase<ndarray::OwnedRepr<u64>, ndarray::Dim<[usize; 1]>> {
    let mut out = Array::zeros(length as usize);
    for l in labels.iter(){
        out[[(l-1) as usize]] = 1;
    }
    return out
}

fn which_true(array:ndarray::ArrayBase<ndarray::OwnedRepr<bool>, ndarray::Dim<[usize; 1]>>   ) -> Vec::<usize> {
    let mut out = Vec::new();
    for (i, elem) in array.iter().enumerate(){
        if *elem { out.push(i) } 
    }
    return out;
}


fn get_total(counter: &ndarray::ArrayBase<ndarray::OwnedRepr<usize>, ndarray::Dim<[usize; 2]>>, elements: &Vec::<usize>) -> u64 {
    let mut sum = 0;
    let index = which_true(counter.slice(s![.., 0]).mapv(|i| elements.iter().any(|&j| j == i)));
    for idx in index.iter() {
        sum += counter[[*idx, 1]] 
    }
    return sum as u64
}
 
fn which_dataset(n: usize, training: &Vec::<usize>, test: &Vec::<usize>, _validation: &Vec::<usize>) -> &'static str {
    if training.iter().any(|&i| i == n)  {
        return "training";
    } else if test.iter().any(|&i| i == n) { 
        return "test";
    } else { 
        return "validation";
    }
}
