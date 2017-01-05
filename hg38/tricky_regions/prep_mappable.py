#!/usr/bin/env python
"""Prepare mappable regions for hg38 using umap k100
https://www.pmgenomics.ca/hoffmanlab/proj/bismap/

Requires bcbio-nextgen for BED preparation.
"""
import os
import sys
import shutil

import requests

from bcbio.bam.callable import get_ref_bedtool
from bcbio.variation import vcfutils
from bcbio.heterogeneity import chromhacks

OUT_FILE = "umap_k100_mappability.bed"
URL = "https://www.pmgenomics.ca/hoffmanlab/proj/bismap/raw/hg38/k100.umap.bed.gz"

def main(ref_file):
    ref_bedtool = get_ref_bedtool(ref_file, {})

    mappable_file = os.path.basename(URL)
    r = requests.get(URL, stream=True)
    with open(mappable_file, "wb") as f:
        shutil.copyfileobj(r.raw, f)

    ref_bedtool.subtract(mappable_file, nonamecheck=True).saveas(OUT_FILE + ".tmp")
    with open(OUT_FILE + ".tmp") as in_handle:
        with open(OUT_FILE, "w") as out_handle:
            for line in in_handle:
                if chromhacks.is_nonalt(line.split()[0]):
                    out_handle.write("%s\tumap_k100_mappability\n" % line.strip())
    os.remove(OUT_FILE + ".tmp")
    vcfutils.bgzip_and_index(OUT_FILE)
    os.remove(mappable_file)

if __name__ == "__main__":
    main(sys.argv[1])
