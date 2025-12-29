# pan_ACR_pred: Pangenomic Chromatin Accessibility Prediction

A pipeline for predicting chromatin accessibility across genomes. We give a similarity score between a known accessible chromatin region (ACR) and a region with unknown accessibility (potentially in a different genome of the same species) via the co-linear chaining of motifs which are candidate cis-regulatory elements (CREs). This creates a feature vector for the unknown region, allowing a machine learning algorithm to classify the region as accessible or inaccessible.

The tool currently supports regions of known accessibility from a single genome and regions with unknown accessibility from a single (potentially different) genome.


## Dependencies

All code is written in Python. The following packages must be installed:
- Bio >= 1.8.1
- biopython >= 1.85
- numpy >= 2.4.0
- pandas >= 2.3.3
- scipy >= 1.16.3
- tqdm >= 4.67.1

For machine learning, the following additional package must be installed:
- scikit_learn >= 1.4.2


## Using the Tool

To run the pipeline, simply navigate to this directory and then run ``$ python ./pan_ACR_pred`` followed by your inputs, as detailed below. 

```
usage: pan_ACR_pred [-h] -o OUTPUT -k KNOWN -u UNKNOWN -kf KNOWN_FIMO [-uf UNKNOWN_FIMO]
                    [-s SCORING] [-r]

Pangenomic Chromatin Accessibility Prediction

options:
  -h, --help            show this help message and exit
  -o OUTPUT, --output OUTPUT
                        Output Directory
  -k KNOWN, --known KNOWN
                        List of known ACRs formatted as ChrX_<start>to<end>
  -u UNKNOWN, --unknown UNKNOWN
                        List of unknown regions formatted as ChrX_<start>to<end>. _pos or _neg
                        can also be appended to the end of the region name for supervised
                        machine learning.
  -kf KNOWN_FIMO, --known_fimo KNOWN_FIMO
                        A directory containing a fimo_out_*/ folder with a fimo.tsv file
                        detailing all candidate CRE motifs in the genome from which the known
                        ACRs originate
  -uf UNKNOWN_FIMO, --unknown_fimo UNKNOWN_FIMO
                        A directory containing a fimo_out_*/ folder with a fimo.tsv file
                        detailing all candidate CRE motifs in the genome from which the unknown   
                        regions originate (if different from fimo_known)
  -s SCORING, --scoring SCORING
                        'None' for unweighted scoring, 'Default' for scoring weighted by motif    
                        relevance, filepath for custom scoring
  -r, --reset           Clear all current temp files when this flag is present (rerun all
                        processes).
```

### Output
The output is a csv file called ``acr_pred.csv``. Each row corresponds to an unknown region, and each column corresponds to a known region. Row and column headers with the names of the regions are present. Each entry in the table is the co-linear chaining score between the unknown and known region corresponding to that row and column. If unknown regions end in _pos or _neg (see Input File Formats), then the last column of the table is "acr_label", where 1 indicates an accessible region and 0 indicates an inaccessible region.

Example:
|                         | Chr3_23445to32145 | Chr4_109884to110998 | ... |
|-------------------------|-------------------|---------------------|-----|
| **Chr1_1198to1259**     | 10.23             | 5.24                |     |
| **Chr3_30192to30498**   | 0                 | 23.3                |     |
| **Chr4_134526to135143** | 1.4               | 18.2                |     |
| ...                     |                   |                     |     |

## Input File Formats

### KNOWN and UNKNOWN Region Lists
KNOWN and UNKNOWN regions should be in two separate files. Each file should contain a single region on each line, formatted as ChrX_\<start>to\<end>, where \<start> is the position of the first base pair in the region, and \<end> is the position of the last base pair in the region.

Example:
```
Chr1_5892to6745
Chr2_90232to10234
Chr1_1235to2123
...
```

UNKNOWN regions optionally may contain "_pos" or "_neg" as a suffix to indicate the accessibility of the region. If all regions contain a suffix, the final csv output will contain 
a column "acr_label", where a 1 indicates an accessible (_pos) region and a 0 indicates an inaccessible (_neg) region. This can then be used for supervised machine learning.

Example:
```
Chr1_5892to6745_pos
Chr2_90232to10234_neg
Chr1_1235to2123_pos
...
```

### Score File
You may pass in a custom score file. The file should contain a header or empty first line. All remaining lines should contain the motif name followed by the score for that motif, tab-separated. The file should give a score for every motif in the fimo.out file.

Example:
```
Motif \t Score
motif_1 \t 8.56
motif_2 \t 5.23
...
```

### FIMO Files

The FIMO files should be generated with the MEME Suite. These files can be generated given a motif file and a genome. More information can be found here: https://meme-suite.org/meme/doc/fimo.html

FIMO should be run on the entire KNOWN and UNKNOWN genomes.

#### Reducing Motif Database Redundancies

Depending on the database, you may want fewer, less redundant motifs. To achieve this, you may run TOMTOM (also from the MEME suite) before FIMO between the motif database file and itself (i.e. use the database file as both the query motifs and target motif database). This calculates similarities between the motifs in database (https://web.mit.edu/meme_v4.11.4/share/doc/tomtom.html). Then, you may run ``./motif_clustering.py`` to cluster the motifs in the database and output a new .meme file which can be used for FIMO. 

## Machine Learning

After the csv file with feature vectors for each unknown region has been generated, binary classification may be used to predict chromatin accessibility. ``./machine_learning/binary_class_tune`` trains and tunes multiple models from Sklearn. Then, the final model can be evaluated with ``./machine_learning/binary_class_test``. 