Find peaks by permutations
=======

The old method of finding CLIP peaks by randomizing read positions.

Uses a config.ini file to set up all the paths.

Needs .bed files of every replicate, plus a couple gtf-based files.

Example: ::

	$ python find_peaks_by_permutations.py -c <directory with config.py>

