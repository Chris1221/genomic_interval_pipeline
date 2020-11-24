use std::collections::HashMap;
use std::fs::File;
use std::io::BufReader;
use std::io::BufRead;
use log::{info, debug};
use std::cmp;
use env_logger::Env;
use std::vec::Vec;

extern crate theban_interval_tree;

use memrange::Range;


fn main() -> std::io::Result<()> {
    // 1. Read in the meta data file. 
    //
    //      This file is simple a one file per line description of all
    //      bed files that are to be included in the database. Create a hash map of
    //      the file names, with each given a numeric index.
    
    let env = Env::default()
        .filter_or("MY_LOG_LEVEL", "info")
        .write_style_or("MY_LOG_STYLE", "always");
    env_logger::init_from_env(env);

    // Hashmap to hold the bed indices
    let mut metadata = HashMap::new();

    // Read in the file paths to the bed files.
    info!("Reading metadata file.");
    {
        let metadata_file = File::open("data/metadata.txt")?;
        let metadata_file_reader = BufReader::new(metadata_file);
        
        let mut i = 1;
        for line in metadata_file_reader.lines(){
            //debug!("Adding a bed file named {}", line.unwrap());
            let line_string = line.unwrap();
            debug!("{}", line_string);
            metadata.insert(line_string, i);
            i += 1;
            debug!("{}", i);
        }
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
    for n in 1..22 {
        let idx = format!("{}{}", "chr", n);
        interval_trees.insert(idx.clone(), theban_interval_tree::IntervalTree::<i32>::new());
        let regions: HashMap<Range, Vec<u64>> = HashMap::new();
        regions_by_chr.insert(idx.clone(), regions);
    }

    //let mut regions: HashMap<Range, Vec<u64>>


    for (bed_file, bed_idx) in metadata {
        let bed = File::open(bed_file)?;
        let bed_reader = BufReader::new(bed);


        for line in bed_reader.lines(){
            let _line = line?;
            let vec = _line.split("\t").collect::<Vec<&str>>();
            let chr = vec[0].to_string();
            let idx_vec = vec![bed_idx];
            let range = Range::new(vec[1].parse::<u64>().unwrap(), vec[2].parse::<u64>().unwrap());

            if interval_trees[&chr].contains(range){
                let mut to_delete: Vec<Range> = Vec::new();
                for (_i, stored) in interval_trees[&chr]
                                            .range(range.min, range.max)
                                                .enumerate(){
                    debug!("Region intersection, resolving...");
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

    println!("{:?}", regions_by_chr);

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

    for n in 1..22 {
        let idx = format!("{}{}", "chr", n);
        qc_interval_trees.insert(
            idx.clone(), 
            theban_interval_tree::IntervalTree::<i32>::new());
    }


    // Let's just make it easy for ourselves and split by chromosome.
    // Reserve two for test, two for validation.
    


    // Later I can deal with the rest, which should be easy.
    //      1. Splitting into training, testing, and validation.
    //      2. Retrieving the sequence from a fastq file 
    //      3. Writing to h5py
    //
    Ok(())
}


fn center_region(region: Range, length: u64 ) -> Range {
    midpoint = region.min + (region.min + region.max) / 2;
}
