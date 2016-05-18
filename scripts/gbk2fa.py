from Bio import SeqIO
from os import listdir, path, makedirs
from os.path import isfile, join

mydir = path.expanduser("~/github/GenomeDorm/data")

# ignore hidden files
gbkDir = mydir + '/gbk/Bacillus_subtilis_subsp_subtilis_str_168/'
onlyfiles = [f for f in listdir(gbkDir) if isfile(join(gbkDir, f)) and not f.startswith('.')]

if not path.exists(mydir + "/fa/Bacillus_subtilis_subsp_subtilis_str_168"):
    makedirs(mydir + "/fa/Bacillus_subtilis_subsp_subtilis_str_168")



for x in onlyfiles:
    IN = mydir + '/gbk/Bacillus_subtilis_subsp_subtilis_str_168/'+  x

    OUT = mydir + '/fa/Bacillus_subtilis_subsp_subtilis_str_168/' + x.split('.')[0] + str('.fa')
    input_handle  = open(IN, "r")
    output_handle = open(OUT, "w")
    count = SeqIO.convert(input_handle, "gb", output_handle, "fasta")
    #for seq_record in SeqIO.parse(input_handle, "genbank") :
    #    print "Dealing with GenBank record %s" % seq_record.id
    #    output_handle.write(">%s %s\n%s\n" % (
    #        seq_record.id,
    #        seq_record.description,
    #        seq_record.seq))

    output_handle.close()
    input_handle.close()
