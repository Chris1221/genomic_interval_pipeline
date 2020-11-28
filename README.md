### A Pipeline for Building Genomic Annotation Datasets for Deep Learning

This is a pipeline for creating [HDF5](https://www.hdfgroup.org/solutions/hdf5) input to [Keras](https://keras.io/) from genomic regions and annotations in [Rust](https://www.rust-lang.org). It is a (somewhat) drop-in replacement for [Basset's preprocessing pipeline](https://github.com/davek44/Basset/blob/master/docs/preprocess.md), intended to transform a list of BED files into annotated one-hot encoded sequences for use in a deep learning model. The input and output of both pipelines should be similar, with this one being *substantially* faster for larger datasets.

## Building from source

Ensure that you have installed [`cargo` and Rust](https://doc.rust-lang.org/cargo/getting-started/installation.html) on your system, then clone this repository. 

```sh
git clone git@github.com:Chris1221/genomic_interval_pipeline.rs.git
cd genomic_interval_pipeline.rs
```

Use `cargo` to build the executable. It should figure out all the dependencies for you.

```sh
cargo build --release
```

The binary will be in `target/release/genomic_interval_pipeline`.

## Usage

You must provide:

1. A newline seperated list of gzipped BED files.
2. Path to the reference genome. This must be compressed with `bgzip` and indexed by `samtools faidx`. 

An example of the first can be found in `data/metadata.txt`:

```
data/file1.bed.gz
data/file2.bed.gz
```

To create the relevant reference genome is straightforward. 

1. Download your reference genome of choice, for example hg19. Here, I just used the UCSC genome browser. 

```sh
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
```

2. Compress your reference genome with `bgzip` (if this is already true, skip this step.)

```sh
gunzip hg19.fa.gz
bgzip hg19.fa
```

3. Index with `samtools`.

```sh
samtools faidx hg19.fa.gz
```

### Running the pipeline

Invoke the binary with the paths to the metadata of BED files and the reference genome (you don't have to specify where the index is). 

```
genomic_interval_pipeline -i data/metadata.txt -f hg19.fa.gz -o small_dataset
```

This will create your dataset at `small_dataset.h5`.

### Arguments


| Short | Long | Value | Description | 
| :-:   | :-:  | :-: | :-- |
| `-i` | `--input` | String | Path to a newline seperated list of bed files to process. | 
| `-f` | `--fastq` | String | Path to faidx indexed, `bgzip` compressed reference FASTQ file |
| `-o` | `--output` | String | Path to the output `.h5` file. |
|     |  `--length` | Number | Standardised length of regions (default: `600`) |
|     | `--test_chr` | String| Comma seperated list of chromosomes to use in the test set (default: `chr19,chr20`) | 
|     | `--valid_chr` | String| Comma seperated list of chromosomes to use in the validation set (default: `chr21,chr22`) | 
|     | `--loglevel` | String | Level of logging (default: `info`) | 

## Dataset Format

HDF5 files are essentially directories of data. There are six tables within the dataset corresponding to the training, test, and validation sequences and their labels. 

Sequences are 3D arrays with dimenions `(batch, length, 4)` where length is optionally specified when building th dataset and refers to the standardized length of the segments. 

Labels are 2D arrays with dimensions `(batch, number_of_labels)` where number of labels is the length of the metadata file. You can easily recode this dataset inside the `HDF5` file for more bespoke training outputs.

## Using the dataset in `Keras`

You can use this data in your own neural network with the [TensorFlow I/O](https://www.tensorflow.org/io) API. Here is an example in Python.

```py
import tensorflow     as tf
import tensorflow_io  as tfio

dataset = "small_dataset.h5"

train_x = tfio.Dataset.from_hdf5(dataset, "training_sequences")
train_y = tfio.Dataset.from_hdf5(dataset, "training_labels")
train_dataset = tf.data.Dataset.zip((train_x, train_y))
```

Pass this `train_dataset` (and similarly test and validation) to `model.fit`. 
