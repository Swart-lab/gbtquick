gbtquick - Quick blobplots from assembly files
----------------------------------------------

Parse assembly headers for coverage info per contig and generate blobplots
without mapping. This is the lazy person's [BlobToolKit](https://blobtoolkit.genomehubs.org/blobtools2/)

Supported assemblers:
 * SPAdes
 * Flye
 * Megahit

For anything else, use [blobtools2](https://github.com/blobtoolkit/blobtools2)

Optionally calculate coding density with Prodigal to discriminate prokaryotic
from eukaryotic sequences.
