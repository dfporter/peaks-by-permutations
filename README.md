Find peaks by permutations
=======

The old method of finding CLIP peaks by randomizing read positions.

Uses a config.ini file to set up all the paths.
Can attempt to generate an ini file with:

```bash
python generate_config_ini_file.py INPUT_DIR
```

Needs .bed files of every replicate, plus a gtf file.
The gtf file, such as:

Caenorhabditis_elegans.WBcel235.78.gtf

Needs to have its header removed, and be linked to in the config.ini file as gtf_raw.
The default file path for this is:

gtf_raw: %(lib)s/Caenorhabditis_elegans.WBcel235.78.noheader.gtf

A second gtf-like file then needs to be made by running:

```bash
python create_gtf_table_from_raw_gtf.py -c auto.ini
```

The output of that command will be put in the file under 'gtf' in the config.ini.
The defauly value is:

gtf: %(lib)s/gtf_with_names_column.txt

With those files made, the peak caller can be run as:

Example:

```bash
python find_peaks_by_permutations.py -c auto.ini
```

