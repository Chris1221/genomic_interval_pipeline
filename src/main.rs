//! # A Pipeline for Building Genomic Annotation Datasets for Deep Learning
//! This is a pipeline for creating [HDF5](https://www.hdfgroup.org/solutions/hdf5) input to [Keras](https://keras.io/) from genomic regions and annotations in [Rust](https://www.rust-lang.org). It is a (somewhat) drop-in replacement for [Basset's preprocessing pipeline](https://github.com/davek44/Basset/blob/master/docs/preprocess.md), intended to transform a list of BED files into annotated one-hot encoded sequences for use in a deep learning model. The input and output of both pipelines should be similar, with this one being *substantially* faster for larger datasets.
//!
//!## Usage
//!
//!You must provide:
//!
//!1. A newline seperated list of gzipped BED files.
//!2. Path to the reference genome. This must be compressed with `bgzip` and indexed by `samtools faidx`. 
//!
//!An example of the first can be found in `data/metadata.txt`:
//!
//!```
//!data/file1.bed.gz
//!data/file2.bed.gz
//!```
//!
//!To create the relevant reference genome is straightforward. 
//!
//!1. Download your reference genome of choice, for example hg19. Here, I just used the UCSC genome browser. 
//!
//!```sh
//!wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
//!```
//!
//!2. Compress your reference genome with `bgzip` (if this is already true, skip this step.)
//!
//!```sh
//!gunzip hg19.fa.gz
//!bgzip hg19.fa
//!```
//!
//!3. Index with `samtools`.
//!
//!```sh
//!samtools faidx hg19.fa.gz
//!```
//!
//!### Running the pipeline
//!
//!Invoke the binary with the paths to the metadata of BED files and the reference genome (you don't have to specify where the index is). 
//!
//!```
//!genomic_interval_pipeline -i data/metadata.txt -f hg19.fa.gz -o small_dataset
//!```
//!
//!This will create your dataset at `small_dataset.h5`.
//!
//! ### Arguments
//!
//!
//!| Short | Long | Value | Description | 
//!| :-:   | :-:  | :-: | :-- |
//!| `-i` | `--input` | String | Path to a newline seperated list of bed files to process. | 
//!| `-f` | `--fastq` | String | Path to faidx indexed, `bgzip` compressed reference FASTQ file |
//!| `-o` | `--output` | String | Path to the output `.h5` file. |
//!|     |  `--length` | Number | Standardised length of regions (default: `600`) |
//!|     | `--test_chr` | String| Comma seperated list of chromosomes to use in the test set (default: `chr19,chr20`) | 
//!|     | `--valid_chr` | String| Comma seperated list of chromosomes to use in the validation set (default: `chr21,chr22`) | 
//!|     | `--loglevel` | String | Level of logging (default: `info`) | 
//!
//!## Dataset Format
//!
//!HDF5 files are essentially directories of data. There are six tables within the dataset corresponding to the training, test, and validation sequences and their labels. 
//!
//!Sequences are 3D arrays with dimenions `(batch, length, 4)` where length is optionally specified when building th dataset and refers to the standardized length of the segments. 
//!
//!Labels are 2D arrays with dimensions `(batch, number_of_labels)` where number of labels is the length of the metadata file. You can easily recode this dataset inside the `HDF5` file for more bespoke training outputs.
//!
//!## Using the dataset in `Keras`
//!
//!You can use this data in your own neural network with the [TensorFlow I/O](https://www.tensorflow.org/io) API. Here is an example in Python.
//!
//!```py
//!import tensorflow     as tf
//!import tensorflow_io  as tfio
//!
//!dataset = "small_dataset.h5"
//!
//!train_x = tfio.Dataset.from_hdf5(dataset, "training_sequences")
//!train_y = tfio.Dataset.from_hdf5(dataset, "training_labels")
//!train_dataset = tf.data.Dataset.zip((train_x, train_y))
//!```
//!
//!Pass this `train_dataset` (and similarly test and validation) to `model.fit`. 


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
#[structopt(name = "genomic_interval_pipeline")]
struct Opt {

    /// Newline seperated list of bed files to process.
    #[structopt(short = "i", long, default_value = "data/metadata.txt")]
    input: PathBuf,

    /// Reference sequence 
    #[structopt(short = "f", long, parse(from_os_str), default_value = "data/blah.fq")]
    fastq: PathBuf,

    /// Output filename
    #[structopt(short = "o", long, parse(from_os_str), default_value = "dataset.h5")]
    output: PathBuf,
   
    /// Size of the regions to report.
    #[structopt(short = "l", long, default_value = "600")]
    length: u64,

    /// Test chromosomes, seperated by commas.
    #[structopt(long, default_value = "chr19,chr20")]
    test_chr: String,
    
    /// Validation chromosomes, seperated by commas.
    #[structopt(long, default_value = "chr21,chr22")]
    valid_chr: String,

    /// If creating a multi-label dataset, exclude cases where both are present.
    #[structopt(short = "e", long)]
    exclusive: bool,

    /// Minimum overlapping sequence to merge. 
    #[structopt(short = "m", long, default_value = "200")]
    min_overlap: u64, 

    /// Log level. Defaults to Info (useful information and statistics). 
    #[structopt(long, default_value = "info")]
    loglevel: String,
}

/// Main function.
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

    let test = to_numeric_vector(opt.test_chr);
    let validation = to_numeric_vector(opt.valid_chr);
    
    let mut training = Vec::new();
    for i in 1..23 { 
        if !test.contains(&i) {
            if !validation.contains(&i) {
                training.push(i) 
            }
        }
    }

    info!(target: "Input", "Test: {:?}", test);
    info!("Validation: {:?}", validation);
    info!("Training: {:?}", training); 

    let mut valid_chromosomes = Vec::new();
    for i in 1..23 { valid_chromosomes.push(format!("{}{}", "chr", i)) }

    // Hashmap to hold the bed indices
    let mut metadata = HashMap::new();
 
    // Read in the file paths to the bed files.
    info!("Reading metadata file.");

    let mut i = 1;
    let mut values = Vec::new();
    let mut custom = false;
    {
        let metadata_file = File::open(opt.input)?;
        let metadata_file_reader = BufReader::new(metadata_file);
        
        for line in metadata_file_reader.lines(){
            //debug!("Adding a bed file named {}", line.unwrap());
            let line_string = line.unwrap();
            let vec = line_string.split(" ").collect::<Vec<&str>>();
            println!("{:?}", vec);
            if vec.len() == 1 {
                debug!("{}", line_string);
                metadata.insert(line_string, i);
                i += 1;
                debug!("{}", i);
            } else {
                metadata.insert(vec[0].to_string(), vec[1].parse::<u64>().unwrap());
                values.push(vec[1].parse::<u64>().unwrap());
                custom = true;
            }

        }
    }
    let number_of_labels: u64;
    
    if custom {
        values.sort();
        values.dedup();

        println!("Unique: {}", values.len());
        number_of_labels = values.len() as u64;
    } else {
        number_of_labels = i-1;
    }

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
            let mut idx_vec = vec![bed_idx];
            let mut range = Range::new(vec[1].parse::<u64>().unwrap(), vec[2].parse::<u64>().unwrap());

            // If not a numeric chromosome, skip for now.
            if !valid_chromosomes.iter().any(|k| k == &chr) {
                continue;
            }

            let mut to_delete = Vec::new();

            // Loop over all intersections, redefining range and idx_vec as I go
            for (_i, stored) in interval_trees[&chr].range(range.min, range.max).enumerate(){
                // Get the current value for that region
                let overlap = get_overlap(range, stored.0) ;

                if overlap >= opt.min_overlap {
                    let v = regions_by_chr
                        .get_mut(&chr)
                        .unwrap()
                        .get_mut(&stored.0)
                        .unwrap()
                        .clone();
                    regions_by_chr
                        .get_mut(&chr)
                        .unwrap()
                        .remove(&stored.0);
                    
                    // Extend the vector of the original
                    idx_vec.extend(v.iter().cloned());
                    idx_vec.sort();
                    idx_vec.dedup();
                    
                    if !((range.min == stored.0.min) & (range.max == stored.0.max)) {
                        // Create a new region based on both
                        debug!("Merging regions {}-{} and {}-{}", range.min, range.max, stored.0.min, stored.0.max);
                        let new = Range::new(
                            cmp::min(range.min, stored.0.min), 
                            cmp::max(range.max, stored.0.max));
                        debug!("    Now {}-{}", new.min, new.max);
                         
                        // Replace the range with this new one
                        range = new;
                        
                        // Add this to the list of regions to be deleted
                        debug!("Deleting region {}-{}", stored.0.min, stored.0.max);
                        to_delete.push(stored.0);
                    }

                }
            }
            for del in to_delete.iter(){
                interval_trees.get_mut(&chr).unwrap().delete(*del);
                debug!("Removed region {:?}", del)
            } 

            // Add back in the new vector of bed files
            // for the range.
            //if !v.contains(&bed_idx){
            //   v.push(bed_idx);
            //
            //}
            debug!("Inserting range {}-{}", range.min, range.max);
            interval_trees.get_mut(&chr).unwrap().insert(range, 1);
            regions_by_chr
                .get_mut(&chr)
                .unwrap()
                .insert(range, idx_vec);

            
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

            // Check if this overlaps with anything already in the interval trees

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
    
    let file = hdf5::File::create(opt.output).unwrap();
    
    //let mut seqs = HashMap::new();
    //let mut labels = HashMap::new();
    let mut index = HashMap::new();
    let mut j = HashMap::new();
    let mut i = HashMap::new();
    let mut seqs_db = HashMap::new();
    let mut labels_db = HashMap::new();
    let mut coordinates_db = HashMap::new();
     
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
                       .create(&format!("{}_{}", label, "seqs"), (1, 600, 4) )
                       .unwrap());
        labels_db.insert(label, file.new_dataset::<u64>()
                       .resizable(true)
                       .create(&format!("{}_{}", label, "labels"), (1, number_of_labels as usize) )
                       .unwrap());
        coordinates_db.insert(label, file.new_dataset::<u64>()
                       .resizable(true)
                       .create(&format!("{}_{}", label, "coords"), (1,3) )
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
            let coords_writer = coordinates_db.get(&dataset).unwrap().as_writer();

            // Expand the dimensionality of the database to ensure there's enoug
            // room for the incoming data.
            seqs_db
                .get(&dataset)
                .unwrap()
                .resize((this_i as usize, opt.length as usize, 4))
                .unwrap();
            labels_db
                .get(&dataset)
                .unwrap()
                .resize((this_i as usize, number_of_labels as usize))
                .unwrap();
            coordinates_db
                .get(&dataset)
                .unwrap()
                .resize((this_i as usize, 3))
                .unwrap();
            
            index.get_mut(&dataset)
                .unwrap()
                .push(region_string.clone());

            // Fetch sequence from fasta
            let seq = fa.fetch_seq_string(&chr_string, region.min as usize, (region.max -1) as usize).unwrap();
            // One hot encode
            let out = one_hot_encode_seq(&seq, 600);
            // One hot encode label vector 
            let label = one_hot_encode_labels(annotation, number_of_labels);

            let coord = array![n, region.min as usize, region.max as usize];

            // If not exclusive, write everything
            if !opt.exclusive { 
                // Write to resized dataset 
                writer.write_slice(&out, s![this_i-1, .., ..]).unwrap();
                label_writer.write_slice(&label, s![this_i-1, ..]).unwrap();
                coords_writer.write_slice(&coord, s![this_i-1, ..]).unwrap();
                *i.get_mut(&dataset).unwrap() += 1;
            } else {
                if label.iter().sum::<u64>() == 1{
                    println!("Not writing region because of exclusive option.");
                    writer.write_slice(&out, s![this_i-1, .., ..]).unwrap();
                    label_writer.write_slice(&label, s![this_i-1, ..]).unwrap();
                    coords_writer.write_slice(&coord, s![this_i-1, ..]).unwrap();
                    *i.get_mut(&dataset).unwrap() += 1;
                }
            }

            
            // One hot encode the labels and enter them into the array
            //labels.get_mut(&dataset)
             //   .unwrap()
              //  .slice_mut(s![.., this_i]).assign(&label);

            
        }
        pb.inc();
    }

    Ok(())
}


/// Standardises regions to a certain length.
///
/// If the regions are less than the input length (given as a parameter),
/// half of the length is added to each end of the center of the region. 
///
/// Example:
///
/// ```
/// let region = memrange::Region(5000,6000);
/// let centered: memrange::Region = center_region(region, 600);
/// ```
///
/// The resulting region will be inserted into the interval tree.
fn center_region(region: &memrange::Range, length: u64 ) -> Range {

    debug!("Correcting range: {}-{}, Length {}", region.min, region.max, length);
    let midpoint = (region.min + region.max) / 2;
    if (length / 2) > midpoint {
        return Range::new(0, length);
    } else  {   
        return Range::new(midpoint - length / 2, midpoint + length / 2);
    }
}


/// One hot encode a DNA sequence.
///
/// Takes a string of ACTG (any case) as input and returns an (length, 4) array. N (unknown
/// nucleotide) will be treated as 0s across all four categories while any other character
/// will throw an error.
///
/// Example:
///
/// ```
/// let encoding = one_hot_encode_seq(seq = 'ACTGN', length = 5);
/// println!("{:?}", encoding);
/// ```
fn one_hot_encode_seq(seq: &String, length: u64) -> ndarray::ArrayBase<ndarray::OwnedRepr<u64>, ndarray::Dim<[usize; 2]>> {
    let mut out = Array::zeros((length as usize, 4));
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
            out[[i, code]] = 1
        }
    }
    return out;
}

/// One hot encode a label vector
///
/// Takes a vector of labels and the total number of labels and creates a one hot encoding.
///
/// ```
/// one_hot_encode_labels(vec![1,2,3], 3);
/// ```
fn one_hot_encode_labels(labels: &Vec::<u64>, length: u64) -> ndarray::ArrayBase<ndarray::OwnedRepr<u64>, ndarray::Dim<[usize; 1]>> {
    let mut out = Array::zeros(length as usize);
    for l in labels.iter(){
        out[[(l-1) as usize]] = 1;
    }
    return out
}

/// Returns which elements of a vector are true (for indexing).

fn which_true(array:ndarray::ArrayBase<ndarray::OwnedRepr<bool>, ndarray::Dim<[usize; 1]>>   ) -> Vec::<usize> {
    let mut out = Vec::new();
    for (i, elem) in array.iter().enumerate(){
        if *elem { out.push(i) } 
    }
    return out;
}

/// Returns the sum of specific entries in a 2D array.
fn get_total(counter: &ndarray::ArrayBase<ndarray::OwnedRepr<usize>, ndarray::Dim<[usize; 2]>>, elements: &Vec::<usize>) -> u64 {
    let mut sum = 0;
    let index = which_true(counter.slice(s![.., 0]).mapv(|i| elements.iter().any(|&j| j == i)));
    for idx in index.iter() {
        sum += counter[[*idx, 1]] 
    }
    return sum as u64
}

/// Takes a chromosome and returns whether it is currently in the training, test, or validation
/// dataset.
fn which_dataset(n: usize, training: &Vec::<usize>, test: &Vec::<usize>, _validation: &Vec::<usize>) -> &'static str {
    if training.iter().any(|&i| i == n)  {
        return "training";
    } else if test.iter().any(|&i| i == n) { 
        return "test";
    } else { 
        return "validation";
    }
}

/// Parses command line options for test and validation chromosomes.
fn to_numeric_vector(input: String) -> Vec::<usize> {
    let split = input.split(",").collect::<Vec<&str>>();
    let mut out = Vec::new();
    for i in &split {
        out.push(i[3..].parse::<usize>().unwrap());
    }
    return out;
}

#[derive(hdf5::H5Type, Clone, PartialEq, Debug)]
#[repr(C)]
pub struct Coord {
    chr: u64,
    start: u64,
    end: u64
}

fn get_overlap(interval1: Range, 
                    interval2: Range) -> u64 {

    if interval2.min > interval1.max || interval1.min > interval2.max {
        return 0; // does not contribute
    }

    let start = std::cmp::max(interval1.min, interval2.min);
    let stop = std::cmp::min(interval1.max, interval2.max);

    return stop - start;
}
