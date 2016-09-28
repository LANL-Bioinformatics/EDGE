#!/usr/bin/env python 
# -*- coding: utf-8 -*-
import csv
import mysql.connector as my
import sys, getopt
import codecs

#read from EDGE output summary.txt file that has the format:
#DATASET TOOL    LEVEL   TOP1    TOP2    TOP3    TOP4    TOP5
#allReads        gottcha-genDB-b genus   Comamonas       Delftia Klebsiella      Alicycliphilus  Erythrobacter
#allReads        gottcha-genDB-b species N/A     N/A     N/A     N/A     N/A
#allReads        gottcha-genDB-b strain  N/A     N/A     N/A     N/A     N/A
#allReads        gottcha-genDB-v genus   Unassigned genus - Bacillus phage phBC6A52      N/A     N/A     N/A     N/A
#allReads        gottcha-genDB-v species N/A     N/A     N/A     N/A     N/A
#allReads        gottcha-genDB-v strain  N/A     N/A     N/A     N/A     N/A
#allReads        gottcha-speDB-b genus   Comamonas       Delftia Klebsiella      Alicycliphilus  Erythrobacter
#allReads        gottcha-speDB-b species Comamonas testosteroni  Delftia acidovorans     Delftia sp. Cs1-4       Klebsiella pneumoniae   Alicycliphilus denitrificans
#allReads        gottcha-speDB-b strain  N/A     N/A     N/A     N/A     N/A
#allReads        gottcha-speDB-v genus   Unassigned genus - Bacillus phage phBC6A52      N/A     N/A     N/A     N/A
#allReads        gottcha-speDB-v species Bacillus phage phBC6A52 N/A     N/A     N/A     N/A
#allReads        gottcha-speDB-v strain  N/A     N/A     N/A     N/A     N/A
#allReads        gottcha-strDB-b genus   Comamonas       Delftia Alicycliphilus  Erythrobacter   Ralstonia
#allReads        gottcha-strDB-b species Comamonas testosteroni  Delftia acidovorans     Delftia sp. Cs1-4       Alicycliphilus denitrificans    Erythrobacter litoralis
#allReads        gottcha-strDB-b strain  Comamonas testosteroni CNB-2    Delftia acidovorans SPH-1       Acinetobacter baumannii MDR-ZJ06        Delftia sp. Cs1-4 chromosome    Alicycliphilus denitrificans K601
#allReads        gottcha-strDB-v genus   Unassigned genus - Bacillus phage phBC6A52      N/A     N/A     N/A     N/A
#allReads        gottcha-strDB-v species Bacillus phage phBC6A52 N/A     N/A     N/A     N/A
#allReads        gottcha-strDB-v strain  Bacillus prophage phBC6A52      N/A     N/A     N/A     N/A
#allReads        metaphlan       genus   Comamonas       Delftia Pseudomonas     Acinetobacter   Flavobacterium
#allReads        metaphlan       species Comamonas testosteroni  Delftia acidovorans     Pseudomonas putida      Acinetobacter junii     Sphingopyxis alaskensis
#allReads        metaphlan       strain  N/A     N/A     N/A     N/A     N/A
#allReads        bwa     genus   Comamonas       Delftia Pseudomonas     Acidovorax      Bacillus
#allReads        bwa     species Comamonas testosteroni  Pseudomonas putida      Delftia acidovorans     Delftia sp. Cs1-4       Acidovorax sp. KKS102
#allReads        bwa     strain  Comamonas testosteroni CNB-2    Delftia acidovorans SPH-1       Pseudomonas putida H8234        Bacillus toyonensis BCT-7112    Pseudomonas putida BIRD-1
#allReads        kraken_mini     genus   Comamonas       Delftia Pseudomonas     Acidovorax      Bacillus
#allReads        kraken_mini     species Comamonas testosteroni  Pseudomonas putida      Delftia acidovorans     Acidovorax sp. KKS102   Delftia sp. Cs1-4
#allReads        kraken_mini     strain  Comamonas testosteroni CNB-1    Comamonas testosteroni CNB-2    Delftia acidovorans SPH-1       Pseudomonas putida H8234        Pseudomonas putida NBRC 14164

#we want these lines:
#allReads        gottcha-strDB-b strain  Comamonas testosteroni CNB-2    Delftia acidovorans SPH-1       Acinetobacter baumannii MDR-ZJ06        Delftia sp. Cs1-4 chromosome    Alicycliphilus denitrificans K601
#allReads        gottcha-strDB-v strain  Bacillus prophage phBC6A52      N/A     N/A     N/A     N/A

#eventually will need to grab sample name and metadata so output file can be named appropriately or will need to put output file in an appropriate dir
#get input file name as command line arg
def main(argv):
     infile_name = ''
     outfile_name = ''
     host = ''
     db_name = ''
     user = ''
     passwd = ''
     try:
          opts, args = getopt.getopt(sys.argv[1:],"h:i:o:d:n:u:p:",["ifile=", "ofile=", "dbhost=", "dbname=", "dbuser=", "passwd="])
          #print "opts = ", opts
     except getopt.GetoptError:
          print "error getting args: usage: "
          print 'identify_pathogens.py -i <input file> -o <output file> -d <db host> -n <db name> -u <db user> -p <db passwd>'
          sys.exit(2)
     if len(sys.argv) == 1:
          print "usage: "
          print 'identify_pathogens.py -i <input file> -o <output file> -d <db host> -n <db name> -u <db user> -p <db passwd>'
          sys.exit()
     for opt, arg in opts:
          if opt in ("-h", "--help") :
               print 'identify_pathogens.py -i <inputfile> -o <output file> -d <db host> -n <db name> -u <db user> -p <db passwd>'
               sys.exit()
          elif opt in ("-i", "--ifile"):
               infile_name = arg
          elif opt in ("-o", "--ofile"):
               outfile_name = arg
          elif opt in ("-d", "--dbhost"):
               host = arg
          elif opt in ("-n", "--dbname"):
               db_name = arg
          elif opt in ("-u", "--dbuser"):
               user = arg
          elif opt in ("-p", "--passwd"):
               passwd = arg
          else:
               print "Bad argument: ", opt, arg
               sys.exit()
     #print 'input file is ', infile_name
     #print 'output file is ', outfile_name
     #print 'host is ', host
     #print 'db_name is ', db_name
     #print 'user is ', user
     #print 'passwd is ', passwd
     #open file, grab strain names in list
     with open(infile_name) as tabfile: 
          taxreader = csv.reader(tabfile, delimiter="\t")
          btaxa = ''
          vtaxa = ''
          for row in taxreader:
               #print "row = ", row
               if row[0] == 'DATASET':
                    continue
               elif row[1] == 'gottcha-strDB-b' and row[2] == 'strain':
                    btaxa = (row[3], row[4], row[5], row[6], row[7]) #need to make this dynamic in terms of number of elements
               elif row[1] == 'gottcha-speDB-b' and row[2] == 'species':
                    btaxa = (row[3], row[4], row[5], row[6], row[7]) #need to make this dynamic in terms of number of elements
               elif row[1] == 'gottcha-genDB-b' and row[2] == 'genus':
                    btaxa = (row[3], row[4], row[5], row[6], row[7]) #need to make this dynamic in terms of number of elements
               elif row[1] == 'gottcha-strDB-v' and row[2] == 'strain':
                    vtaxa = (row[3], row[4], row[5], row[6], row[7]) #need to make this dynamic in terms of number of elements
                    #print "strain vtaxa = ", vtaxa
               elif row[1] == 'gottcha-speDB-v' and row[2] == 'species':
                    vtaxa = (row[3], row[4], row[5], row[6], row[7]) #need to make this dynamic in terms of number of elements
                    #print "species vtaxa = ", vtaxa
               elif row[1] == 'gottcha-genDB-v' and row[2] == 'genus':
                    vtaxa = (row[3], row[4], row[5], row[6], row[7]) #need to make this dynamic in terms of number of elements
                    #print "genus vtaxa = ", vtaxa
               else:
                    print "no rows match 'gottcha-xxxDB-x' and 'strain, species or genus': ", row[1], row[2], row[3], row[4], row[5], row[6], row[7]
          #combine bacterial and viral taxa
          taxa = btaxa+vtaxa
          #print "line 112 all taxa = ", taxa
          #open output file for writing
          outfile = codecs.open(outfile_name, 'w',encoding='utf8')
          outfile.write("pathogen"+"\t"+"host(s)"+"\t"+"disease(s)"+"\n")
          #open database connection
          try: 
               db = my.connect(host=host, user=user, password=passwd, database=db_name,buffered=True)

               #for each strain in list, query pathogen database and if found and if host is human, write pathogen name, host and disease info to tab delimited file
               for i in taxa:
                    #prepare cursor object
                    cursor = db.cursor()

                    #print "line 124 taxa line = ", i
                    #if the line contains weird characters, such as '/', remove those characters
                    if '/' in i:
                         i = i.replace('/', '')
                    #split pathogen name into genus and species in case strain is not in pathogenDB
                    gs_array = i.split(" ")
                    #print "line 128 gs_array = ", gs_array
                    genus = gs_array[0]
                    if len(gs_array) > 1:
                         species = gs_array[1]
                    else:
                         species = "none"
                    #print "genus, species:", genus, species
                    query=("SELECT id FROM pathogen WHERE genome_name = %s LIMIT 0, 1")
                    cursor.execute(query, (i,))
	            row = cursor.fetchone()
                    if row is None:
                         print "no entry in pathogen table with exact pathogen name: ", i
                         #try to find genus, species in database
                         query=("SELECT id FROM pathogen " "WHERE genus = %s AND species = %s LIMIT 0, 1")
                         cursor.execute(query, (genus,species,))
			 row = cursor.fetchone()
                         if row is None:
                              print "no entry in pathogen table with genus and species: ", genus, species
                              #lastly try to find genus only
   	                      cursor.execute("""SELECT id FROM pathogen WHERE genus = %s LIMIT 0, 1""", (genus,))
                              row = cursor.fetchone()
                              if row is None:
                                   print "no entry in pathogen table with genus: ", genus
                              else: #found genus so grab pathogen name, host and disease, check if host includes human 
                                   #and if so write to file
                                   #print "line 147 genus = ", genus
                                   #cursor.execute("""SELECT id FROM pathogen WHERE genus = %s""", (genus,))
                                   path_PK = row[0]
                                   #print "pathogen id = ", path_PK
                                   cursor.execute(""" select p.genome_name, p.id, ph.host_pathogen_id, ph.host_id, h.name, pd.disease_pathogen_id, d.disease_name from pathogen p, pathogen_host ph, pathogen_disease pd,host h, disease d where p.id = %s and p.id = ph.host_pathogen_id and ph.host_id = h.id and p.id = pd.disease_pathogen_id and pd.disease_id = d.id LIMIT 0, 1""", (path_PK,))
                                   #grab data from cursor
                                   pdata = cursor.fetchone()
                                   pathogen_name =  pdata[0]
                                   host =  pdata[4]
                                   disease =  pdata[6]
				   cursor.close()
                                   #write to file
                                   if 'Homo sapiens' in host and 'None' not in disease and 'Unknown' not in disease:
                                        outfile.write(pathogen_name+"\t"+host+"\t"+disease+"\n")

                         else: #found genus and species so grab pathogen name, host and disease, check if host includes 
                         #human and if so write to file
                              #print "line 160 genus and species are ", genus, species
                              #cursor.execute("""SELECT id FROM pathogen WHERE genus = %s AND species = %s""", (genus,species,))
                              #path_PK = cursor.fetchone()[0]
                              path_PK = row[0]
                              #print "line 169 pathogen id = ", path_PK
                              cursor.execute(""" select p.genome_name, p.id, ph.host_pathogen_id, ph.host_id, h.name, pd.disease_pathogen_id, d.disease_name from pathogen p, pathogen_host ph, pathogen_disease pd,host h, disease d where p.id = %s and p.id = ph.host_pathogen_id and ph.host_id = h.id and p.id = pd.disease_pathogen_id and pd.disease_id = d.id LIMIT 0, 1""", (path_PK,))
                              pdata = cursor.fetchone()
                              #print "line 172 pdata = ", pdata
                              pathogen_name =  pdata[0]
                              #print "line 170 pathogen name = ", pathogen_name
                              host =  pdata[4]
                              disease =  pdata[6]
			      cursor.close()
                              #write to file
                              if 'Homo sapiens' in host and 'None' not in disease and 'Unknown' not in disease:
                                   outfile.write(pathogen_name+"\t"+host+"\t"+disease+"\n")
                    else: #found entire pathogen name so grab pathogen name, host and disease, check if host includes human and if so write to file
                         #print "line 177 pathogen name = ", pathogen_name
                         #cursor.execute("""SELECT id FROM pathogen WHERE genome_name = %s""", (i,))
                         #path_PK = cursor.fetchone()[0]
                         path_PK = row[0]
                         cursor.execute(""" select p.genome_name, p.id, ph.host_pathogen_id, ph.host_id, h.name, pd.disease_pathogen_id, d.disease_name from pathogen p, pathogen_host ph, pathogen_disease pd,host h, disease d where p.id = %s and p.id = ph.host_pathogen_id and ph.host_id = h.id and p.id = pd.disease_pathogen_id and pd.disease_id = d.id LIMIT 0, 1""", (path_PK,))
                         pdata = cursor.fetchone()
                         pathogen_name =  pdata[0]
                         host =  pdata[4]
                         disease =  pdata[6]
			 cursor.close()
                         #write to file
                         if 'Homo sapiens' in host and 'None' not in disease and 'Unknown' not in disease:
                              outfile.write(pathogen_name+"\t"+host+"\t"+disease+"\n")
          except my.Error as e:
               print >> sys.stderr, "MySQL error: ", e
          except:
               print >> sys.stderr, "Unknown error occurred" 
          outfile.close()
          db.close()

if __name__ == "__main__":
   main(sys.argv[1:])
