# THIS IS A WORKING DRAFT - NOT FOR USE YET


# Tech Cookbook / Cheat-Sheet:

## General Hackathon Info
    
### Got Git?
        
* Get git here:  [https://git-scm.com/downloads](https://git-scm.com/downloads)
    
### Public Keys
        
For access to GitHub and Hackathon Servers, you'll need an _ssh key_.

* [Generating a new SSH Key for Github](https://help.github.com/articles/generating-a-new-ssh-key-and-adding-it-to-the-ssh-agent/)
    
### Markdown, for writing and formatting ReadMe and other documents on GitHub (like this document!)

* [Markdown Help](https://commonmark.org/help/)
        
* [Handy Markdown cheatsheet](https://github.com/adam-p/markdown-here/wiki/Markdown-Cheatsheet)

## Where’s the Data?!

|Data               |Location     |
|-------------------|--------------|
|SRA alignments| `gs://ncbi_sra_rnaseq/*.bam`|
| Counts | `gs://ncbi_sra_rnaseq/genecounts/*.genecounts` |
|Contig FASTAs| `gs://ncbi_sra_rnaseq/*.contigs.fa/`|
|BigQuery Tables|project: `strides-sra-hackathon-data` <br/>table: `rnaseq.genescounts`, `rnaseq.runinfo` <br>project: `strides-sra-hackathon-data`<br/>table: various, see 4a TODO: [below](#)|
|Server information| Refer to the pinned post Slack (`#help-desk` channel)|
|Additional Tools|`gs://ncbi_hackathon_aux_tools/`|

## How the data was generated
We start with SRA run data pulled from NCBI and align it with hisat2 onto gs://ncbi_sra_rnaseq/refs/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.*.ht2
Unmapped reads then are assembled using skesa into denovo contigs and aligned with hisat2 onto the contigs. The merged BAM file is stored on the bucket. 
BAM files are sorted and indexed, for each BAM there is a flagstat and counts files. The counts are also stored in the BigQuery table rnaseq.genescounts.

## General GCP Advice
    
### BigQuery (a.k.a bq)
        
* [https://cloud.google.com/bigquery/docs/bq-command-line-tool](https://cloud.google.com/bigquery/docs/bq-command-line-tool)
        
> For reference use only -- _should not_ be included in an external docker
        
* Search data pre-indexed by NCBI
            
  * *Query format*: `bq --project_id strides-sra-hackathon-data query --nouse_legacy_sql`
                
    > In some case you may need to use `strides-sra-hackathon-data`, more on this below
            
* *Useful options*: 
    ```
    bq --max_rows=100 
       --format=pretty 
       --project_id strides-sra-hackathon-data query     
       --nouse_legacy_sql "<your_standard_sql_query_here"
       ```

> The format flag also takes “json” and “csv” as valid arguments
            

* *Available data*:
    `bq show --schema --format=prettyjson strides-sra-hackathon-data:rnaseq.runinfo`

> Format also accepts “json”

* *Example query*: to get a list of library selection strategies used by included SRRs:
    ```   
    bq --format=csv --project_id strides-sra-hackathon-data query
       --nouse_legacy_sql 
         "select distinct LibrarySelection from rnaseq.runinfo"  
    ```

* > Expected Output:**
    ```
    Waiting on bqjob_r2257ca5f589030fa_000001681a2c36ac_1 ... (0s)
    Current status: DONE  
LibrarySelection
RANDOM
ChIP
PolyA
other
DNase
cDNA
size fractionation
CAGE
Reduced Representation
unspecified
Hybrid Selection
    ```
    
* If you have a complex query you are interested in, especially if you think it might involve tables in `strides-sra-hackathon-data`, please let us know if we can help in the TODO: `#help-desk` Slack Channel.


### Some other useful commands

####  Working with Google Storage (aka `gs`) buckets:

* `gsutil` is a collection of command line tools to access and modify data stored in google storage buckets. [`gsutil` Documentation Link](https://cloud.google.com/storage/docs/gsutil)
        
* `gsutil` Examples:
    * List files:

      `gsutil ls -l`

      `gs://ncbi_sra_rnaseq/DRR016694*.bam`
        
    * Copy a file from the bucket to the current directory: 
    
      `gsutil cp gs://ncbi_sra_rnaseq/DRR016694.bam .`
        
    * Stream a file:

      `gsutil cat gs://ncbi_sra_rnaseq/DRR016694.flagstat | less`
        
    * Copy multiple files from the bucket:

      `gsutil -m cp gs://ncbi_sra_rnaseq/DRR01669\*.bam .`
    

> Tip: To copy data between servers, 
> * use gsutil to copy from the source server to google storage first, then copy to the destination.
> * or, add your ssh keys to the source server and use `scp localfile username@REMOTE-IP:/foo/bar/`


#### Google compute cloud tools (gcloud):
        
* Check your service account and project: `gcloud info`
        
* Pub/sub queue:
  * Create a topic: `gcloud pubsub topics create test1`
            
  * Create a subscription to the topic: `gcloud pubsub subscriptions create --topic test1 sub1`
            
  * Publish a message to the topic: `gcloud pubsub topics publish test1 --message qwe`
            
  * Pull a message from the subscription: `gcloud pubsub subscriptions pull sub1`

# Server Access

For server access, create an SSH key using `ssh-keygen -t rsa -b 4096 -C "your-email@domain.edu"`, then post the public key in the`#public_keys' slack channel.  

Server IP information is listed in the `#help-desk` channel in slack (pinned message).

For most servers, access them using `ssh your-email@1.2.3.4` where `your-email` is the name part of the email address you used when generating your SSH key, and `1.2.3.4` is the IP address of the server

* A quick way to verify access is to use the command `ssh your-email@1.2.3.4 whoami` - this will echo back your username if successful.

* If you have trouble connecting, check the username, and try first with a simple command-line `ssh` client.

## Premade servers

Several solr and other database servers are available.

* For Solr: URL Access: [http://IP-ADDRESS:7790/solr/\#/](http://IP-ADDRESS:7790/solr/#/) (requires browser authentication).

> See the pinned post in the `#help-desk` channel in slack for server IP addresses. Contact a friendly admin for username or password information.

* Solr is installed as a docker image `solr01`.
  * Copy files to the image using `docker cp localfile solr01:/path/`
  * Restart solr using `sudo docker container restart solr01`
  
* Example solr API Query using `curl`

  `$ curl -u USERNAME:PASSWORD "IP-ADDRESS:7790/solr/tstcol01/select?q=\*:\*"`  

  > Response:
    ```
     {  
     "responseHeader":{  
     "Status":0,  
     "QTime":0,  
     "Params":{  
     "q":"\*:\*"}},  
     "response":{"numFound":0,"start":0,"docs":\[\]  
     }}  
    ```

  * Data Import (all on one line)
    
    ``` 
    curl -u USERNAME:PASSWORD 
         "IP-ADDRESS:7790/solr/tstcol01/update/json/docs?commit=true" 
         -X POST 
         -H 'Content-Type:application/json' 
         --data-binary "FULL_PATH_TO_DATA.json"  
    ```

  > For example:  

   `curl -u USERNAME:PASSWORD "IP-ADDRESS:7790/solr/tstcol01/update/json/docs?commit=true" -X POST -H 'Content-Type:application/json' --data-binary test_known4.json`  

  > Where test_known4.json contains:  
 
    ```

    [{"hit_id" : 1,
     "Contig" : "SRR123.contig.1",
     "method" : "mmseq2",
     "parameters" : "some parameters",
     "accession": "AC12345.1",
     "simil":0.8},**

    {"hit_id" : 2,
     "Contig" : "SRR123.contig.1",
     "method" : "rpstblastn",
     "parameters" : "some parameters",
     "accession": "AC12345.1",
     "simil":0.8},

    {"hit_id" : 3,
     "Contig" : "SRR123.contig.2",
     "method" : "mmseq2",
     "parameters" : "some parameters",
     "accession": "AC12356.2",
     "simil":0.9},

    {"hit_id" : 4,
     "Contig" : "SRR123.contig.2", 
     "method" : "rpstblastn",
     "parameters" : "some parameters",
     "accession": "AC12345.1",
     "simil":0.4}, git

    {"hit_id" : 5,
    "Contig" : "SRR123.contig.1",
    "method" : "AC12356.2",
    "parameters" : "this is some parameters",
    "accession": "AC12345.1",
    "simil":0.5}
     ]
    ```

## Working with NCBI Data/Tools in the cloud

### SRA Realign Object Data and Contigs
        
Realigned metagenomic SRA runs are located in the `ncbi_sra_realign` bucket. To access them directly please refer to gsutil.
        
The bucket with realigned runs should be mounted on the VMs during the main hackathon event under the directory `/data/realign/`. Realign files are named based on their corresponding SRA run accession number, for example `/data/realign/SRR1158703.realign`

Realign objects are effectively the same reads (bases) as
the original SRA runs but aligned onto host (human) and
either viral contigs or references and onto denovo contigs
assembled with skesa. 

* For example, `SRR649927.realign` contains 964 reads mapped onto human, 9043281 reads mapped onto denovo contigs and 1392735 unmapped reads:
  ```
  bq --maxrows=10000 --projectid strides-sra-hackathon-data 
  query 'select * from ncbisrarealign.summary where accession="SRR649927"'
  ```

* To check taxonomy of the contigs:  
  ```
  bq --projectid strides-sra-hackathon-data 
  query 'select * from ncbisrarealign.taxonomy where accession="SRR649927"'`
  ```

* To extract contigs from a realign object:  
  ```
  dump-ref-fasta --localref SRR649927.realign > SRR649927.contigs.fa
  ```
* Human references GRCh38.p12 was used to align reads before assembling the rest. Please find the full list of sequences here:
  [https://www.ncbi.nlm.nih.gov/assembly?term=GRCh38&cmd=DetailsSearch](https://www.ncbi.nlm.nih.gov/assembly?term=GRCh38&cmd=DetailsSearch)
        
* There are 2 types of contigs in realign objects: guided and denovo. Guided contigs have been built with guided  assembler using a predefined viral reference set. These are almost 300 sequences including multiple Influenza serotypes, Herpes, Ebola etc.

The names of guided contigscan be filtered based on regex:

`grep -P '^(?\!Contig)\[A-Z0-9.\]+_\\d+$'`

Contig names are based on reference sequence accession
used as a guide. For example contig named `KF021598.1_1`
has been built using `KF021598.1` as a guide reference.
        
Denovo contigs have been assembled with skesa after
filtering out human and viral reads. Their names start
with Contig prefix.
        
Contig names are unique only within a realign object, not
across realign objects.
    
### SRA Toolkit
        
(Already installed on pre-built hackathon VM instances) 

Download from [https://www.ncbi.nlm.nih.gov/sra/docs/toolkitsoft/](https://www.ncbi.nlm.nih.gov/sra/docs/toolkitsoft/)  

Or, copy the pre-release version from google storage:  

  `gsutil cp gs://ncbi_hackathon_aux_tools/sratoolkit.2.9.4.pre.tar.gz .`
        
* Examples:
  * Lookup by accession:  `srapath <accession>`
  * Download and cache run locally: `prefetch <accession>`
  * Get statistics: `sra-stat --quick --xml <accession>`
  * Dump reads into separate fastq files: `fastq-dump --split-files --split-spot <accession>`
  * Extract contigs (local sequences): `dump-ref-fasta -l <accession>`
  * Convert an SRA object to BAM:  `sam-dump --unaligned <accession> | samtools view -Sb -<accession>.bam`
    
### BLAST suite
        
Instructions on running the dockerized BLAST are at [https://github.com/ncbi/docker/blob/master/blast/README.md](https://github.com/ncbi/docker/blob/master/blast/README.md) 

> Make sure the BLAST_DIR env variable is set before you try the commands.

For a complete list of command-line blast arguments see: [https://www.ncbi.nlm.nih.gov/books/NBK279684/#_appendices_Options_for_the_commandline_a_](https://www.ncbi.nlm.nih.gov/books/NBK279684/#_appendices_Options_for_the_commandline_a_)
        
BLAST has a lot of flexibility in how it presents results:

* Use the `-outfmt` option to specify how the results are presented.
  * Default value for `-outfmt` is 0 (zero), which produces the standard BLAST report.
  * `-outfmt 6` produces a tabular report with a standard set of fields. 
  * `-outfmt 7` produces a tabular report with comment fields
  * `-outfmt 10` produces CSV.
  > To customize the fields in the tabular/CSV output, see examples in slides 6-8 of
    [https://ftp.ncbi.nlm.nih.gov/pub/education/public_webinars/2018/10Oct03_Using_BLAST/Using_BLAST_Well2.pdf](https://ftp.ncbi.nlm.nih.gov/pub/education/public_webinars/2018/10Oct03_Using_BLAST/Using_BLAST_Well2.pdf)
            
  * Use `-outfmt 11` to save your results as a BLAST archive that you can then reformat with the
    `blast_formatter` application. See slide 9 in [https://ftp.ncbi.nlm.nih.gov/pub/education/public_webinars/2018/10Oct03_Using_BLAST/Using_BLAST_Well2.pdf](https://ftp.ncbi.nlm.nih.gov/pub/education/public_webinars/2018/10Oct03_Using_BLAST/Using_BLAST_Well2.pdf)
        
BLAST Databases.
            
* To see BLAST databases on a prepared GCP instance, use the comand
   `docker run --rm ncbi/blast update_blastdb.pl --showall pretty --source gcp`
* Use `update_blastdb.pl` to download the databases you need.  
* Example (note that GCP is specified as source, `BLASTDB_DIR` is defined and points at a directory that
    exists): 
    ```
    docker run --rm -v $BLASTDB_DIR:/blast/blastdb:rw -w /blast/blastdb ncbi/blast update_blastdb.pl --source gcp swissprot_v5
    ```
            
* Use `blastdbcmd` to interrogate the databases.
  * Summary of database: `blastdbcmd -db swissprot_v5 -info`
    > full command with docker is `docker run --rm -v $BLASTDB_DIR:/blast/blastdb:ro -w /blast/blastdb ncbi/blast blastdbcmd -db swissprot_v5 -info`

* Print accession, scientific name, and title for all entries: 
  `blastdbcmd -db swissprot_v5 -entry all -outfmt "%a %S %t"`
            
* Print accession and title for all human (taxid 9606) entries:
  `blastdbcmd -db swissprot_v5 -taxids 9606 -outfmt "%a %t"`
            
* Print accession, scientific name, and title for `P02185`:
  `blastdbcmd -db swissprot_v5 -entry p02185 -outfmt "%a %S %t"`

BLAST programs (this is a brief summary)
    
* `blastn`: DNA-DNA comparisons. Add `-task blastn` to make it more sensitive (and much slower). Default is megablast.
    
* `blastp`: protein-protein comparisons. Add `-task blastp-fast` to make it faster (longer words for initial matches) and only marginally less sensitive.

* `blastx`: (translated) DNA query-protein comparison. Add `-task blastx-fast` to make it faster.

* `tblastn`: protein query-(translated) DNA database comparison. Add `-task tblastn-fast` to make it faster.
    
* `rpsblast`: protein-PSSM comparison (PSSM is position-specific-scoring matrix). Great for identifying domains in your proteins.

* `rpstblastn`: (translated) DNA-PSSM comparison. Great for identifying domains in your translated DNA queries.
    
* `magicblast`: mapper for DNA-DNA comparisons. Can quickly map reads (even long PacBio ones) and identify splice sites.
  > Documentation at [https://ncbi.github.io/magicblast/](https://ncbi.github.io/magicblast/)

* You can extract FASTAs from a BLAST DB:
  `blastcmd -db <db_name> -entry all -outfmt %f -out <file_name>`

* You can blast against a subset of a db using a:
  * gi list: `-gilist <file>`
  * Negative gi list: `-negative_gilist <file>`

### Other tools pre-installed on your VMs (including non-NCBI tools)
    
#### Dockerized
        
Some VMs will have these tools installed:
            
* [https://github.com/NCBI-Hackathons/NCBI_PowerTools_Docker/blob/master/Dockerfile](https://github.com/NCBI-Hackathons/NCBI_PowerTools_Docker/blob/master/Dockerfile)
    
#### Non-Dockerized

| hisat2         | skesa         | guidedassembler_graph | compute-coverage |
| -------------- | ------------- | ---------------------- | ---------------- |
| python         | pip           | python3                | pip3             |
| C++            | R             | Anaconda2              | conda            |
| bedtools       | GATK          | picard                 | BWA              |
| MiniMap 2      | BowTie 2      | EDirect                | HMMer            |
| samtools       | bcftools      | HTS JDK                | HTS lib          |
| STAR           | abyss         | plink-ng               | cufflinks        |
| cytoscape-impl | velvet        | tophat                 | FastQC           |
| HTSeq          | mcl           | muscle                 | MrBayes          |
| GARLI          | Clustal omega | Dedupe                 | trintyrnaseq     |

## Jupyter notebook server setup


Based on [this cheatsheet](https://www.digitalocean.com/community/tutorials/how-to-install-anaconda-on-ubuntu-18-04-quickstart).

> make sure that `wget` and `bzip2` are installed.

1. Get the latest version from https://www.anaconda.com/download/#linux (note your version or filename may be different).

> `$ wget https://repo.continuum.io/archive/Anaconda3-2018.12-Linux-x86_64.sh`

2. Verify your download 

> `$ sha256sum Anaconda3-5.2.0-Linux-x86_64.sh`

3. Run the script

> `$ bash Anaconda3-5.2.0-Linux-x86_64.sh`

4. Create a jupyter config file

> `$ jupyter notebook --generate-config`

5. Create a self-signed SSL certificate

> `$ mkdir certs``
> `$ cd certs/`
> `$ sudo openssl req -x509 -nodes -days 365 -newkey rsa:1024 -keyout mycert.pem -out mycert.pem`

6. Generate a password hash (you will use the same password to access your jupyter notebook remotely)

> `$ ipython` 
> `In [1]: from IPython.lib import passwd`
> `In [2]: passwd()                                                                             
>  `   Enter password: 
>  `   Verify password: 
> `Out[3]: 'sha1:0bebc03bff39:25b6ebdkrisbf0cc3f7e6c8e0f82fcbe9e66801a'
> Copy the `sha` string for later use.

7. Edit `.jupyter/jupyter_notebook_config.py` to add the following
> Note the path to your mycert.pem file from step 5.
> Note the `sha` hash from step 6

```
######################
c.NotebookApp.allow_remote_access = True  
# Kernel config
c.IPKernelApp.pylab = 'inline'  # if you want plotting support always in your notebook
# Notebook config
c.NotebookApp.certfile = u'/home/username/certs/mycert.pem' #location of your certificate file
c.NotebookApp.ip = '*'
c.NotebookApp.open_browser = False  #so that the ipython notebook does not opens up a browser by default
c.NotebookApp.password = u'sha1:ff35ddacee39:e9faaaaa1c6d5008cb50160ad1a00527c9dec09'  #edit this with the SHA hash that you generated after typing in Step 9
# This is the port we opened in Step 3.
c.NotebookApp.port = 8080
 ###############
```
 
 8. Start a screen instance so you can close your terminal without killing jupyter
 
 > `$ screen -DR notebook
 
 9. Run the notebook as root
 
 > `$ sudo /home/username/anaconda3/bin/jupyter notebook --config /home/username/.jupyter/jupyter_notebook_config.py --allow-root`
 
 10. Press `ctrl+a` then `d` to detach from the screen session.
 
 11. Browse to https://SERVER:8080/ and enter the password you supplied in step 6, above.
