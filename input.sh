ssh.exe -m hmac-sha2-512 fdn@bioinformatics.bmb.stud-srv.sdu.dk

# First time only - clone the repo:
# git clone https://github.com/ElFlemmo/tumour_read_mapping.git /data/Flemming/git/tumour_read_mapping

cd /data/Flemming/git/tumour_read_mapping
git pull

conda activate mapping

cd /data/Flemming/tumour

# Remove old index built from wrong fasta (contig names didn't match GFF3)
rm -f combined_revised_phages.*.bt2

nohup bash /data/Flemming/git/tumour_read_mapping/run_all.sh /data/Flemming/tumour combined_revised_phages.fna combined_revised_phages.gff3 > pipeline.log 2>&1 &
tail -f pipeline.log
