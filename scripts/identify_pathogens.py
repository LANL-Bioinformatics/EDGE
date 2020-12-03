# -*- coding: utf-8 -*-
#!/usr/bin/python
import csv
import mysql.connector
import sys, getopt
from mysql.connector import errorcode

#read from EDGE GOTTCHA2 output, for example the file allReads-gottcha-speDB-v.list.txt, which has the format:
#LEVEL   TAXA    ROLLUP  ASSIGNED        LINEAR_LENGTH   TOTAL_BP_MAPPED HIT_COUNT       HIT_COUNT_PLASMID       READ_COUNT      LINEAR_DOC      NORM_COV
#species Enterobacteria phage P1 0.3151          5711    202680  5737    0       2986    35.489406408685 0.059995118372865
#species Enterobacteria phage HK022      0.2218          1092    149016  4222    0       2102    136.461538461538        0.230689295252336
#species Escherichia phage TL-2011b      0.0593          3613    39311   1115    0       562     10.8804317741489        0.0183934547881928
#species Enterobacteria phage mEp237     0.0574          3246    32276   935     0       544     9.94331484904498        0.0168092513162218
#species Stx2-converting phage 1717      0.0480          2139    27469   800     0       455     12.8419822346891        0.0217094711430241
#species Enterobacteria phage HK225      0.0477          2263    32688   928     0       452     14.4445426425099        0.0244186120133924
#species Escherichia phage phiV10        0.0459          5689    27409   778     0       435     4.81789418175426        0.0081446876967643
#species Salmonella phage RE-2010        0.0341          2272    16394   474     0       323     7.21566901408451        0.0121981447549225
#species Enterobacteria phage YYZ-2008   0.0250          1848    17327   489     0       237     9.37608225108225        0.0158503401846061
#species Escherichia phage HK639 0.0210          106     6186    201     0       199     58.3584905660377        0.0986554835336536
#species Klebsiella phage phiKO2 0.0193          2214    7820    235     0       183     3.53206865401987        0.00597098961190713
#species Enterobacteria phage N15        0.0180          1551    8647    245     0       171     5.57511283043198        0.00942477172912059

#get input file name as commandline arg
def main(argv):
     bact_infile_name = ''
     virus_infile_name = ''
     outfile_name = ''
     host = ''
     db_name = ''
     user = ''
     passwd = ''
     btaxa = [] #dynamic list to add bacterial species names to
     vtaxa = [] #dynamic list to add viral species names to
     taxa = [] #for combined list of bacterial and viral species names 
     hosts_list = []
     diseases_list = []
     not_in_db = []
     not_in_db_filename = 'species_not_in_pathogen_DB.txt'
     print "running identify_pathogens.py...."
     try:
          opts, args = getopt.getopt(sys.argv[1:],"h:b:v:o:d:n:u:p:",["bifile=", "vifile=", "ofile=", "dbhost=", "dbname=", "dbuser=", "passwd="])
          #print "opts = ", opts
     except getopt.GetoptError:
          print "error getting args: usage: "
          print 'identify_pathogens2.py -b <input bact file> -v <input virus file>  -o <output file> -d <db host> -n <db name> -u <db user> -p <db passwd>'
          sys.exit()
     if len(sys.argv) == 1:
          print "usage: "
          print 'identify_pathogens2.py -b <input bact file> -v <input virus file> -o <output file> -d <db host> -n <db name> -u <db user> -p <db passwd>'
          sys.exit()
     for opt, arg in opts:
          if opt in ("-h", "--help") :
               print 'identify_pathogens2.py -b <input bact file> -v <input virus file> -o <output file> -d <db host> -n <db name> -u <db user> -p <db passwd>'
               sys.exit()
          elif opt in ("-b", "--bifile"):
               bact_infile_name = arg
          elif opt in ("-v", "--vifile"):
               virus_infile_name = arg
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
     #print 'bact input file is ', bact_infile_name
     #print 'virus input file is ', virus_infile_name
     #print 'output file is ', outfile_name
     #print 'host is ', host
     #print 'db_name is ', db_name
     #print 'user is ', user
     #print 'passwd is ', passwd

     #open bacterial input file, grab species (strain) names in list
     with open(bact_infile_name) as tabfile: 
          taxreader = csv.reader(tabfile, delimiter="\t")
          for row in taxreader:
               #print "row = ", row
               if row[0] == 'LEVEL':
                    continue
               elif row[0] == 'species':
                    #add species name to btaxa
                    btaxa.append(row[1]) #dynamic in terms of number of elements
               #else:
                    #print "no rows match 'species': ",row[0], row[1], row[2], row[3], row[4], row[5], row[6], row[7], row[8], row[9]
     tabfile.close()
     #open viral input file, grab species (strain) names in list
     with open(virus_infile_name) as tabfile: 
          taxreader = csv.reader(tabfile, delimiter="\t")
          for row in taxreader:
               #print "row = ", row
               if row[0] == 'LEVEL':
                    continue
               elif row[0] == 'species':
                    #add species name to vtaxa
                    vtaxa.append(row[1]) #dynamic in terms of number of elements
               #else:
                    #print "no rows match 'species' in column 1: ",row[0]
     tabfile.close()
     #combine bacterial and viral taxa
     taxa = btaxa+vtaxa
     #print "all taxa = ", taxa
     #open output file for writing
     outfile = open(outfile_name, 'w')
     outfile.write("pathogen"+"\t"+"host(s)"+"\t"+"disease(s)"+"\n") #file header
     nidb = open(not_in_db_filename, 'w')
     nidb.write("species not in pathogen database\n")
     #open database connection
     try: 
          #use mysql connector
          db = mysql.connector.connect(host = host, user = user, password = passwd, database = db_name)
          #prepare cursor object
          cursor = db.cursor(buffered=True)

          #for each species in list, query pathogen database and if found and if host is human, write pathogen name, host and disease info to tab delimited file
          for i in taxa:
               hosts_list = []
               diseases_list = []
               #print "line 125 taxa line = ", i
               #if the line contains weird characters, such as '/', remove those characters
               #if '/' in i:
               #     i = i.replace('/', '')
               #check for entire pathogen name ('genome_name') in database first
               cursor.execute("""SELECT id FROM pathogen WHERE genome_name = %s""", (i,))
               #print "at line 131, after execute"
               if cursor.fetchone() is None:
                    #print "no entry in pathogen table with genome_name: ", i
                    #if taxa doesn't match name exactly, split into genus, species, strain and search for those
                    #split pathogen name into genus and species in case strain is not in pathogenDB
                    gs_array = i.split(" ")
                    #print "line 137 gs_array = ", gs_array
                    genus = gs_array[0]
                    if len(gs_array) > 1:
                         species = gs_array[1]
                         #print "genus, species:", genus, species
                         cursor.execute("""SELECT id FROM pathogen WHERE genus = %s AND species = %s""", (genus,species,))
                         rows = cursor.fetchall()
                         if not rows:
                              #print "no entry in pathogen table with genus and species: ", genus, species
                              nidb.write(i+"\n")
                              #the below code should not be necessary if the database is comprehensive
                              #the below code will find genus and write the genus to the list if any species in that genus is pathogen
                              #lastly try to find genus only
                              #cursor.execute("""SELECT id FROM pathogen WHERE genus = %s""", (genus,))
                              #rows = cursor.fetchall()
                              #if not rows:
                              #     print "no entry in pathogen table with genus: ", genus
                              #else: #found genus so grab pathogen name, host and disease, check if host includes human 
                                   #and if so write to file
                              #     print "line 156 genus = ", genus
                              #     cursor.execute("""SELECT id FROM pathogen WHERE genus = %s""", (genus,))
                                   #there may be many rows with genus if there are  many strains in the database
                                   #need to fetchall
                              #     pks = cursor.fetchall()
                              #     print "line 161 pathogen ids = ", pks
                                   #convert pks from list of tuples to list of pk values
                              #     pk_list = list(sum(pks, ()))
                              #     print "first element = ", pk_list[0]
                              #     for pk in pk_list:
                              #           path_PK = pk
                              #           print "line 167 pathogen id = ", path_PK
                              #           cursor.execute("""SELECT p.genome_name, p.id, ph.host_pathogen_id, ph.host_id, h.name, pd.disease_pathogen_id, d.disease_name from pathogen p, pathogen_host ph, pathogen_disease pd,host h, disease d where p.id = %s and p.id = ph.host_pathogen_id and ph.host_id = h.id and p.id = pd.disease_pathogen_id and pd.disease_id = d.id""", (path_PK,))
                                         #grab data from cursor
                              #           pdata = cursor.fetchone()
                                         #pathogen_name =  pdata[0]
                              #           host =  pdata[4]
                              #           host =  host.encode('ascii',errors='ignore')
                              #           disease =  pdata[6]
                              #           disease = disease.encode('ascii',errors='ignore')
                                         #write to file and break out of loop
                              #           if 'Homo sapiens' in host and 'None' not in disease and 'Unknown' not in disease:
                              #                outfile.write(genus+"\t"+host+"\t"+disease+"\n")
                              #                break
                         else: #found genus and species so grab host and disease, check if host includes 
                              #human and if so write to file
                              cursor.execute("""SELECT id FROM pathogen WHERE genus = %s AND species = %s""", (genus,species,))
                              pks = cursor.fetchall()
                              #convert pks from list of tuples to list of pk values
                              pk_list = list(sum(pks, ()))
                              #print "line 186 first element = ", pk_list[0]

                              #need to fix this loop - it writes redundant entries
                              #store pathogen, host, disease then write to file later?
                              for pk in pk_list:
                                   path_PK = pk
                                   #print "line 192 pathogen id = ", path_PK
                                   cursor.execute("""SELECT p.genome_name, p.id, ph.host_id, h.name, d.disease_name, d.id from pathogen p, pathogen_host ph, pathogen_disease pd,host h, disease d where p.id = %s and p.id = ph.host_pathogen_id and ph.host_id = h.id and p.id = pd.disease_pathogen_id and pd.disease_id = d.id""", (path_PK,))
                                   #grab data from cursor - this will return a list of tuples, which are the database records retrieved
                                   pdata = cursor.fetchall()
                                   #print "line 196 pdata = ", pdata
                                   #print "line 197 pdata[0] = ", pdata[0]
                                   #have to loop through tuples in pdata - there is 1 tuple for each record retrieved
                                   #there will be multiple records if more than 1 host, disease per pathogen
                                   for j in pdata:
                                        host =  j[3]
                                        host =  host.encode('ascii',errors='ignore')
                                        #print "line 203 host: ", host
                                        #add host to hosts_list
                                        if host not in hosts_list:
                                             hosts_list.append(host)
                                        disease =  j[4]
                                        disease = disease.encode('ascii',errors='ignore')
                                        #print "line 209 disease: ", disease
                                        #add disease to diseases_list
                                        if disease not in diseases_list:
                                             diseases_list.append(disease)
                              #for debugging
                              #print "hosts list: ", hosts_list
                              #print "diseases list: ", diseases_list
                              #print "\n"
                              #write to file
                              if 'Homo sapiens' in hosts_list:
                                   #count number of 'none', 'None', 'Unknown', and 'unknown' in diseases_list
                                   #if number of the above is < total number of diseases, write to file
                                   nones = diseases_list.count('none') + diseases_list.count('None')
                                   unks = diseases_list.count('unknown') + diseases_list.count('Unknown')
                                   #print "# nones = ", nones
                                   #print "#unks = ", unks
                                   if nones+unks < len(diseases_list):
                                        #outfile.write(genus+"\t"+species+"\t"+", ".join(str(i) for i in hosts_list)+"\t"+", ".join(str(i) for i in diseases_list)+"\n")
                                        #uncomment the above line if multiple hosts are to be displayed
                                        outfile.write(genus+" "+species+"\t"+"Homo sapiens"+"\t"+", ".join(str(i) for i in diseases_list)+"\n")
                                        #break
               else: #found entire pathogen name so grab pathogen name, host and disease, check if host includes human and if so write to file
                    cursor.execute("""SELECT id FROM pathogen WHERE genome_name = %s""", (i,))
                    #print "line 232 pathogen name = ", genome_name
                    path_PK = cursor.fetchone()[0]
                    #print "line 234  pathogen id = ", path_PK
                    cursor.execute("""select p.genome_name, p.id, ph.host_id, h.name, d.disease_name, d.id from pathogen p, pathogen_host ph, pathogen_disease pd,host h, disease d where p.id = %s and p.id = ph.host_pathogen_id and ph.host_id = h.id and p.id = pd.disease_pathogen_id and pd.disease_id = d.id""", (path_PK,))
                    pdata = cursor.fetchall()
                    #print "line 237 pdata = ", pdata
                    #have to loop through tuples in pdata - there is 1 tuple for each record retrieved
                    #there will be multiple records if more than 1 host, disease per pathogen
                    for j in pdata:
                         pathogen_name =  j[0]
                         pathogen_name = pathogen_name.encode('ascii',errors='ignore')
                         #print "line 243 pathogen_name = ", pathogen_name      
                         host =  j[3]
                         host =  host.encode('ascii',errors='ignore')
                         #print "line 246 host: ", host
                         #add host to hosts_list
                         if host not in hosts_list:
                              hosts_list.append(host)
                         disease =  j[4]
                         disease = disease.encode('ascii',errors='ignore')
                         #print "line 252 disease: ", disease
                         #add disease to diseases_list
                         if disease not in diseases_list:
                              diseases_list.append(disease)
                    #write to file
                    if 'Homo sapiens' in hosts_list:
                         #count number of 'none', 'None', 'Unknown', and 'unknown' in diseases_list
                         #if number of the above is < total number of diseases, write to file
                         nones = diseases_list.count('none') + diseases_list.count('None')
                         unks = diseases_list.count('unknown') + diseases_list.count('Unknown')
                         #print "# nones = ", nones
                         #print "#unks = ", unks
                         if nones+unks < len(diseases_list):
                              #outfile.write(pathogen_name+"\t"+", ".join(str(i) for i in hosts_list)+"\t"+", ".join(str(i) for i in diseases_list)+"\n")
                              #uncomment the above line if multiple hosts are to be displayed
                              outfile.write(pathogen_name+"\t"+"Homo sapiens"+"\t"+", ".join(str(i) for i in diseases_list)+"\n")                       
                              #break
     except mysql.connector.Error as err:
          if err.errno == errorcode.ER_ACCESS_DENIED_ERROR:
               print("Something is wrong with your user name or password")
          elif err.errno == errorcode.ER_BAD_DB_ERROR:
               print("Database does not exist")
          else:
               print(err)
               db.close()
     nidb.close()
     outfile.close()
     print "all done!\n"
if __name__ == "__main__":
   main(sys.argv[1:])
