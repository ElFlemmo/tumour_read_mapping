ssh.exe -m hmac-sha2-512 fdn@bioinformatics.bmb.stud-srv.sdu.dk

cd /data/Flemming/tumour

ls

conda activate mapping

bash run_all.sh /data/Flemming/tumour combined_revised_phages.fasta combined_revised_phages.gff3
