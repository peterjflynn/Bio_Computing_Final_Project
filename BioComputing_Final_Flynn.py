
# coding: utf-8

# # Final project - Peter Flynn
# # ECEV 32000 - Introduction to Scientific Computing
# # Professor- Dr. Stefano Allesina, PhD.
# # Winter, 2017

# My goal for this project is to create a script where the researcher will be able to blast any organism genome against all known viral genomes. The organism genome as well as the most up-to-date virus repository need to be downloaded from NCBI. I was trying to mimic the results of a paper where they blasted the pillbug *Armadillidium vulgare* genome against all complete viral genomes (Theze *et al.* 2014: doi:10.1093/gbe/evu163). 
# 
# This project presented some difficulties because I had to learn how to use the FTP server, work with large genome files, and how to use the blast software. 
# 
# Also, this script can only be performed using an os system since several of the commands seemed to be easier to complete in bash than to carry out exclusively in python.
# 

# The first step is to create a directory where all your files will be stored and where you place the script and files. This directory can be made in your home directory for ease of use.

# # Part 1: Virus Genomes

# Since the virus genomes need to be downloaded from NCBI, I first thought that I could use the ENTREZ module to download all (or any) of the virus genomes. It turns out (after days of trying) that even though ENTREZ is a great and easy way to download NCBI nucleotide and protein databases, it does not download complete genomes using the NCBI genome database. Therefore, I needed to figure out a different way to download all up-to-date complete virus genomes on NCBI. In order to do this I used the File Transfer Protocol (FTP) server. The python module"ftplib" is useful for connecting to an FTP server to retrieve files and process them locally. 
# 
# To get onto the FTP server for the NCBI database you have to connect to the host at "ftp.ncbi.nlm.nih.gov". From there you have to get to the virus genome directory within the genome database and find the file "all.fna.tar.gz". This file is refreshed whenever there is a new complete virus genome added to NCBI. The code below will download this file to the directory where this script is located.  

# In[1]:

#download your virus database
import ftplib
from ftplib import FTP
host = "ftp.ncbi.nlm.nih.gov"
ftp= FTP(host)
ftp.login()
ftp.cwd('genomes/Viruses') #find the genome virus directory within database
#ftp.retrlines('LIST') #this lists the directory contents
ftp.retrbinary('RETR all.fna.tar.gz', open('all.fna.tar.gz','wb').write) 
ftp.close() 
print("Downloaded virus .fna file")


# Next I created a subdirectory where the virus genomes can be stored once downloaded from the NCBI database. I titled this file "all_virus_files" for future reference. 

# In[2]:

#create virus file in main directory
import os
import sys
path = r'./all_virus_files' 
os.mkdir(path)


# I wanted to move the virus genome file into a subdirectory since there are more than 4500 seperate files in the all.fna.tar.gz once you untar it.

# In[3]:

#move virus genome file to its own folder
import shutil
src = "./all.fna.tar.gz"
dst = "./all_virus_files/"
shutil.move(src, dst)


# The command below will "untar" the all.fna.tar.gz file. Since the tar.gz file is compressed, one needs to decompress this file. I used the module tarfile, which makes it possible to read and extract the files from a tar.gz compressed file. I wrote a function that would allow me to extract the files or "untar" the files from all.fna.tar.gz

# In[4]:

#untar virus file 
import tarfile,sys 
def untar(file):
    if (file.endswith("tar.gz")):
        tar = tarfile.open(file)
        tar.extractall(path="./all_virus_files/")
        tar.close()
        print("Extracted in Current Directory")

untar("./all_virus_files/all.fna.tar.gz")


# Now that we have the all.fna.tar.gz file decompressed for the complete virus genome database, every complete genome for every virus is inside its own subdirectory folder within all_virus_files. Therefore, we want to pull out the .fna file from each subdirectory.  I wrote a function that moves all the viruses from their own folder into the main all_virus_files folder.

# In[5]:

# pull all virus fna files out of their seperate subfolders
import shutil
RootDir1 = r'./all_virus_files/'
TargetFolder = r'./all_virus_files/'
for root, dirs, files in os.walk((os.path.normpath(RootDir1)), topdown=False):
        for name in files:
            if name.endswith('fna'):
                SourceFolder = os.path.join(root,name)
                shutil.copy2(SourceFolder, TargetFolder) #copies fna to new folder


# Now that all the seperate virus .fna files are no longer in subfolders, you can easily concatenate all of the files into one .fna file to prepare for the BLAST. The easiest way I found to do this was to us the "cat" function in BASH. Once the virus genomes are all merged into one file it can be used as the query in the BLAST.

# In[6]:

#concatenate the virus files 
import os 
os.system("cd ./all_virus_files/ ;cat *.fna > merged_viruses.fna")


# This step moves the merged virus file back into the home directory so it can be easily found when performing the BLAST.

# In[7]:

#move merged_viruses.fna into the home directory
import shutil 
src = "./all_virus_files/merged_viruses.fna"
dst = "./"
shutil.move(src, dst)


# # Part 2: Organism Genome

# For this part of the project I chose to download the *Armadillidium vulgare* genome as my organism genome, however, any complete organism genome on the NCBI database could be input here. I first found where the genome was located on the FTP and then downloaded the specific file which contains the complete and most up-to-date genome. This script only works if you know where the organism's genome is located in the FTP. You can choose any organism genome you want, just replace the ftp.cwd with the path to that genome and wherever it states: GCA_001887335.1_A_vulgare_v1_genomic.fna.gz with the genome file you want to download.

# In[8]:

#importing your organism genome 
import ftplib
from ftplib import FTP
host = "ftp.ncbi.nlm.nih.gov"
ftp= FTP(host)
ftp.login()
ftp.cwd('genomes/genbank/invertebrate/Armadillidium_vulgare/latest_assembly_versions/GCA_001887335.1_A_vulgare_v1')
ftp.retrlines('LIST')
ftp.retrbinary('RETR GCA_001887335.1_A_vulgare_v1_genomic.fna.gz', open('GCA_001887335.1_A_vulgare_v1_genomic.fna.gz','wb').write)
ftp.close()


# This is not part of the script, but if you hypothetically did not know where the organism genome was located, you could download the entire *Armadillidium vulgare* genome folder on the NCBI FTP server and then pick out the .fna files to see which might be the genome that you are looking for. To download the entire genome you could use the **wget** command in BASH and download all the files in the *Armadillidium vulgare* directory. 

# In[ ]:

#potentially another way to download genomes from NCBI
#import os 
#os.system("wget -r -nH -np -R index.html ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/invertebrate/Armadillidium_vulgare/*")


# Once the organism genome is downloaded from the NCBI database you have to decompress the file. Since this is a .gz file not a .tar file you can import the gzip module and decompress the genome file by reading it into the gzip. I also renamed the file genome.fna.

# In[9]:

# ungzip genome file this time use gzip because not a .tar.gz file but just a .gz file
import gzip
inF = gzip.open('GCA_001887335.1_A_vulgare_v1_genomic.fna.gz', 'rb')
outF = open('genome.fna', 'wb')
outF.write( inF.read() )
inF.close()
outF.close()


# # Part 3: BLAST

# One can perform a BLAST search online through a module in Biopython called Blast, however, since the viral genomes file is so large, we cannot use BLAST for this script. Therefore, it is easier to run BLAST locally on your computer. This piece of the script only works for an OS system, however on the FTP server there are files for linux and windows systems.

# In[10]:

#download blast
import ftplib
from ftplib import FTP
host = "ftp.ncbi.nlm.nih.gov"
ftp= FTP(host)
ftp.login()
ftp.cwd('/blast/executables/blast+/LATEST/')
ftp.retrlines('LIST')
ftp.retrbinary('RETR ncbi-blast-2.6.0+-x64-macosx.tar.gz', open('ncbi-blast-2.6.0+-x64-macosx.tar.gz','wb').write) 
ftp.close()
print('downloaded blast')


# Since the BLAST install file is a tar.gz file, it needs to be decompressed, which is again done through the tarfile module, in the same function we used for decompressing the viral genomes. 

# In[11]:

#untar blast install file 
import tarfile,sys 
def untar(file):
    if (file.endswith("tar.gz")):
        tar = tarfile.open(file)
        tar.extractall(path="./")
        tar.close()
        print("Extracted in Current Directory")

untar("./ncbi-blast-2.6.0+-x64-macosx.tar.gz")


# Once the BLAST files are decompressed I put these files into the home directory so there is an easy place to find the BLAST application for future use.

# In[12]:

#copy program files into home directory
import shutil 
src = "./ncbi-blast-2.6.0+"
dst = "../"
shutil.move(src, dst)


# Once the BLAST file is installed, it is time to prepare your files for the BLAST.
# 
# Either the organism genome file or the virus genome file needs to be made into the database. According to the Theze et al. 2014, they used the pillbug genome as their database and the virus genomes file as their query. Therefore, I am making the organism genome the database for this blast. The database is created and called: "genome_db" which consists of three different files in your directory. 
# makeblastdb is the command in BLAST that creates a database
# -in is the input file
# -dbtype specifies what kind of data the input file (nucleotide, protein, etc.)
# -out specifies what you want the output to be called

# In[13]:

#make organism genome the database
import os 
os.system("makeblastdb -in genome.fna -dbtype nucl -out genome_db")


# Since tblastx takes a long time to run, I have shortened the merged viruses file to only 4000 lines (merged_viruses_short.fna)

# In[19]:

#cut file to be small for tblastx
import os 
os.system('head -4000 merged_viruses.fna > merged_viruses_short.fna')


# It is now time to run your BLAST. We will be performing a tblastx, which compares the six-frame translations of a nucleotide query sequence against the six-frame translations of a nucleotide sequence database. This is the most computationally expensive BLAST technique. It is a powerful gene-predicition tool.
# 
# 
# -db is designating which file is the database
# 
# -query is designating which files is the query 
# 
# -evalue describes the number of hits one can "expect" to see by chance when searching a database of a particular size. An evalue of 1 was set by Theze et al 2014.
# 
# -culling_limit discards hits that have a certain number of reads that are better than it (in this script it only saves the best hit)
# 
# -outfmt is the format you want the output, 5 = xml file
# 
# -out specifies what you want the output to be called
# 
# Since tBlastx takes several hours to run for this job, I provided the end result of the tblastx (tblastx.xml) as well as the version with only 4000 lines from the shortened merged viruses file (tblastx_hits_short.xml).

# In[18]:

#run tBlastx
import os 
cmd = "tblastx -query merged_viruses_short.fna -db genome_db -evalue 1 -culling_limit 1 -outfmt 5 -out tblastx_hits_short.xml"
os.system(cmd)


# I have provided the tblastx.xml hits file that would have been produced from the full tblastx.
# 
# This xml file gives us the best protein hit for each virus. If you put the .xml file into Geneious or another software that can read .xml files, it gives you the viral protein fragment hit (hsp_hseq). To perform a subsequent blastp analysis on this file you need to select the genome protein fragment hits instead (hsp_qseq). Cutting out the pillbug genome protein fragment hits can be performed using a regular expression since before every genome protein fragment has a distinctive Hsp_qseq before and after it. I then can join all of these genome protein sequences onto their own line in a new textfile "tblastx.txt. 

# In[14]:

#take out the Hsp_qseq with regex
import Bio
import re
with open("tblastx_hits.xml", "r") as f:  
        strings = (re.findall(r'(?<=<Hsp_qseq>).*(?=\<)', f.read()))
        g = ("\n".join(strings)) #can join each fragment into a newline
        with open("tblastx_hits.txt", "w") as outp:
            outp.write(g)


# Now that we have a text file with every pillbug genome fragment that was a hit on the tblastx, we can add a symbol (>) and unique number identifier for each protein sequence and convert this textfile to a fasta file in order to be able to continue working with it. 

# In[15]:

#add >number and create FASTA file
openFile= open('./tblastx_hits.txt', "r")
newFile = open('tblastx_hits.fasta', 'w')
count = 0
for lines in openFile:
    newFile.write(('>'+str(count)+'\n'+ lines))
    count +=1


# The next step is to blastp the tblastx_hits.fasta file against the non-redundant protein database. The goal now is to screen the selected subset of *A. vulgare* genome sequences as queries to screen for homologous coding sequences in the whole set of nonredundant protein sequences of the National Center for Biotechnology Information (NCBI) database. The non-redundant database is over 100GB of storage, so I do not have that downloaded on my computer. It is best to use a server that has this database already downloaded to run a blastp.
# 
# Since our file is still too large to use online, we will cut a small portion to use run a blastp on.

# In[ ]:

#run Blastp
#import os 
#cmd = "blastp -query tblastx_hits.fasta -db nr -evalue 0.01 -outfmt 5 -out tblastx_hits.xml"
#os.system(cmd)


# This is as far as I have gotten since I have not had a chance to use the Field Museum server to perform a blastp (which has the non-redundant protein database downloaded). 
