
			#à réaliser avant de lancer ce script !
#1)
#mettre dans le meme repertoire que les "6 echantillons fastq", de "reference.gff" et "reference.fasta"

#2)
#creation de l'environnement hackathon.sh, 
#executer les lignes directement en ligne de commande de "environnement.sh" si le script ne marche pas, 
#pour activer l'environnement : conda activate hackathon_env

#3)
#lancer ce script : bash RUN.sh 


			#Etape du projet Hackathon :
nextflow hackathon.nf
#1)#Trimmomatic		#2)indexation		#3)mapping		#4)count
#produit le dossier de sortie output_dir avec les 6 counts.txt

Rscript deseq2.r
#5)#analyse DESEQ2 

#
#le script creer 2 datas frames cts et coldata à partir des 6 counts.txt 

