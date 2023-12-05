
sudo apt update
sudo apt install sra-toolkit
  
######################################################################################## 1h pour telechargement fastq
fasterq-dump --threads 8 --progress SRR10379721
fasterq-dump --threads 8 --progress SRR10379722
fasterq-dump --threads 8 --progress SRR10379723
fasterq-dump --threads 8 --progress SRR10379724
fasterq-dump --threads 8 --progress SRR10379725
fasterq-dump --threads 8 --progress SRR10379726
wget -q -O reference.fasta "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=CP000253.1&rettype=fasta"
wget -O reference.gff "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&report=gff3&id=CP000253.1"
########################################################################################


nextflow hackathon.nf
























