#!/usr/bin/env python

# Po-E (Paul) Li
# B-11, Los Alamos National Lab
# Date: 05/15/2016
# Last update: 05/15/2016

import sys
import io
import os.path
import json
import gzip
import subprocess
import fileinput

####################
# Global variables #
####################

libPath = os.path.dirname(os.path.realpath(__file__))
taxonomyDir = libPath + "/database"
DEBUG=0

taxDepths      = {}
taxParents     = {}
taxRanks       = {}
taxNames       = {}
taxMerged      = {}
taxNumChilds   = {}
accTid         = {}
tidLineage     = {}
tidLineageDict = {}

major_level = {
	'superkingdom' : 'k',
	'phylum'       : 'p',
	'class'        : 'c',
	'order'        : 'o',
	'family'       : 'f',
	'genus'        : 'g',
	'species'      : 's',
	'k'            : 'superkingdom',
	'p'            : 'phylum',
	'c'            : 'class',
	'o'            : 'order',
	'f'            : 'family',
	'g'            : 'genus',
	's'            : 'species'
}

####################
#      Methods     #
####################

def taxidStatus( taxID ):
	if taxID in taxMerged:
		return "merged"

	if taxID in taxNames:
		if '.' in taxID:
			return "valid custom"
		return "valid"
	else:
		return "invalid"
		
# Replace the query taxID with the merged taxID if it exists.
def taxid2mergedTid( taxID ):
	_checkTaxonomy()
	if taxID and taxID in taxMerged.keys():
		return taxMerged[taxID]
	else:
		return taxID

def acc2taxid( acc ):
	_checkTaxonomy()
	accession2taxid_file=taxonomyDir+"/accession2taxid.tsv"
	#remove version number#
	acc = acc.split('.')[0]
	
	if DEBUG: sys.stderr.write( "[INFO] acc2taxid from file: %s\n" % accession2taxid_file )
		

	if not acc in accTid:
		with open( accession2taxid_file ) as f:
			f.seek(0, 2)
			start = 0
			end = f.tell()
			accCur = ""
			
			if DEBUG: sys.stderr.write( "[INFO] acc2taxid from file: %s\n" % accession2taxid_file )
			
			while( acc != accCur and start < end ):
				
				posNew = (end+start)/2
				
				f.seek( posNew )
		
				if posNew != start: f.readline()

				line = f.readline()	
				
				if DEBUG: sys.stderr.write( "[INFO] start: %15d, posNew: %15d, end: %15d, line: %s" % (start, posNew, end, line) )
				if line :
					(accNew, tid) = line.split('\t')
				else:
					break

				if acc > accNew and accCur != accNew and accNew:
					if accNew: posNew = f.tell()
					start = posNew
					if start >= end: end = start+1
				else:
					end = posNew
				
				accCur = accNew

			f.close()

			if accCur == acc:
				accTid[acc] = tid.strip()
			else:
				accTid[acc] = ""

	tid = taxid2mergedTid(accTid[acc])

	return tid

def taxid2rank( taxID, guess_strain=True ):
	_checkTaxonomy()
	taxID = taxid2mergedTid(taxID)
	if not taxID in taxRanks:
		return "unknown"

	if taxID == '1':
		return "root"

	if taxRanks[taxID] == "no rank" and guess_strain:
		# a leaf taxonomy is a strain
		if taxidIsLeaf(taxID):
			return "strain"
		# if not
		else:
			nmtid = taxid2nearestMajorTaxid(taxID)
			nmrank = _getTaxRank(nmtid)
			if nmrank == "species":
				return "species - others"
			else:
				return "others"
	
	return taxRanks[taxID]

def taxid2name( taxID ):
	_checkTaxonomy()
	taxID = taxid2mergedTid(taxID)
	return _getTaxName(taxID)

def taxid2depth( taxID ):
	_checkTaxonomy()
	taxID = taxid2mergedTid(taxID)
	return _getTaxDepth(taxID)

def taxid2type( taxID ):
	_checkTaxonomy()
	taxID = taxid2mergedTid(taxID)
	origID = taxID
	lastID = taxID
	
	taxID = taxid2mergedTid(taxID)
	taxID = taxParents[taxID]

	while taxID != '1' and taxRanks[taxID] != 'species':
		lastID = taxID
		taxID = taxParents[taxID]

	if taxRanks[taxID] != 'species':
		taxID = 0
	else:
		taxID = lastID
		if taxID == origID: taxID = 0

	return taxID

def taxid2parent( taxID ):
	_checkTaxonomy()
	taxID = taxid2mergedTid(taxID)
	taxID = taxParents[taxID]
	while taxID != '1' and taxRanks[taxID] == 'no rank':
		taxID = taxParents[taxID]

	return taxID

def taxid2nameOnRank( taxID, r ):
	_checkTaxonomy()
	taxID = taxid2mergedTid(taxID)
	if taxID == 1: return "root"
	if r == "root": return "root"

	rank = _getTaxRank(taxID)
	name = _getTaxName(taxID)

	if r == "strain" and taxidIsLeaf(taxID):
		return name

	while taxID:
		if rank.upper() == r.upper(): return name
		if name == 'root': break
		taxID = _getTaxParent(taxID)
		rank = _getTaxRank(taxID)
		name = _getTaxName(taxID)

	return ""

def taxid2taxidOnRank( taxID, r ):
	_checkTaxonomy()
	taxID = taxid2mergedTid(taxID)
	rank = _getTaxRank(taxID)
	name = _getTaxName(taxID)

	if r == rank or ( r == 'strain' and rank == 'no rank'): return taxID
	if r == "root": return 1

	while taxID:
		if rank.upper() == r.upper(): return taxID
		if name == 'root': break

		taxID = _getTaxParent(taxID)
		rank = _getTaxRank(taxID)
		name = _getTaxName(taxID)

	return ""

def taxidIsLeaf( taxID ):
	taxID = taxid2mergedTid(taxID)
	if not taxID in taxNumChilds:
		return True
	else:
		return False

def taxid2fullLineage( taxID ):
	_checkTaxonomy()
	taxID = taxid2mergedTid(taxID)
	fullLineage = ""

	while taxID != '1':
		rank = _getTaxRank(taxID)
		name = _getTaxName(taxID)
		if not name: break
		fullLineage += "%s|%s|%s|"%(rank,taxID,name)
		taxID = taxParents[taxID]

	return fullLineage

def taxid2fullLinkDict( taxID ):
	_checkTaxonomy()
	taxID = taxid2mergedTid(taxID)
	fullLineage = ""
	link = {}

	while taxID != '1':
		rank = _getTaxRank(taxID)
		name = _getTaxName(taxID)
		if not name: break

		parID = taxParents[taxID]
		link[parID] = taxID
		taxID = parID

	return link

def taxid2nearestMajorTaxid( taxID ):
	_checkTaxonomy()
	taxID = taxid2mergedTid(taxID)
	ptid = _getTaxParent( taxID )
	while ptid != '1':
		tmp = taxid2rank( ptid )
		if tmp in major_level:
			return ptid
		else:
			ptid = _getTaxParent( ptid )

	return "1"

def taxid2lineage( tid, print_all_rank=1, print_strain=0, replace_space2underscore=1, output_type="auto"):
	return _taxid2lineage( tid, print_all_rank, print_strain, replace_space2underscore, output_type)

def taxid2lineageDICT( tid, print_all_rank=1, print_strain=0, replace_space2underscore=0, output_type="DICT" ):
	return _taxid2lineage( tid, print_all_rank, print_strain, replace_space2underscore, output_type )

def _taxid2lineage(tid, print_all_rank, print_strain, replace_space2underscore, output_type):
	_checkTaxonomy()

	if output_type == "DICT":
		if tid in tidLineageDict: return tidLineageDict[tid]
	else:
		if tid in tidLineage: return tidLineage[tid]

	info = _autoVivification()
	lineage = []
	taxID = tid

	taxID = taxid2mergedTid(taxID)

	level = {
		'k' : '',
		'p' : '',
		'c' : '',
		'o' : '',
		'f' : '',
		'g' : '',
		's' : ''
	}

	rank = taxid2rank(taxID)
	orig_rank = rank
	name = _getTaxName(taxID)
	str_name = name
	if replace_space2underscore: str_name.replace(" ", "_")

	while taxID:
		if rank in major_level:
			if replace_space2underscore: name.replace(" ", "_")
			level[major_level[rank]] = name

			#for output JSON
			info[rank]["name"] = name
			info[rank]["taxid"] = taxID

		taxID = _getTaxParent(taxID)
		rank = _getTaxRank(taxID)
		name = _getTaxName(taxID)

		if name == 'root': break

	# try to get the closest "no_rank" taxa to "type" representing subtype/group (mainly for virus)
	typeTID = taxid2type(tid)
	if typeTID:
		info["type"]["name"]  = _getTaxName(typeTID)
		info["type"]["taxid"] = typeTID

	last = str_name

	ranks = ['s','g','f','o','c','p','k']
	idx = 0
	
	# input taxid is a major rank
	if orig_rank in major_level:
		idx = ranks.index( major_level[orig_rank] )
	# if not, find the next major rank
	else:
		nmtid = taxid2nearestMajorTaxid( tid )
		nmrank = taxid2rank( nmtid )
		if nmrank == "root":
			idx = 7
		else:
			idx = ranks.index( major_level[nmrank] )

	for lvl in ranks[idx:]:
		if print_all_rank == 0:
			if not level[lvl]: continue

		if not level[lvl]:
			level[lvl] = "%s - no_%s_rank"%(last,lvl)
			info[major_level[lvl]]["name"]  = "%s - no_%s_rank"%(last,lvl)
			info[major_level[lvl]]["taxid"] = 0

		last=level[lvl]
		#lineage.append( "%s__%s"%(lvl, level[lvl]) )
		lineage.append( level[lvl] )

	lineage.reverse()

	if print_strain:
		if orig_rank == "strain":
			#lineage.append( "n__%s"%(str_name) )
			lineage.append( "%s"%(str_name) )
			info["strain"]["name"]  = str_name
			info["strain"]["taxid"] = tid

	if output_type == "DICT":
		tidLineageDict[tid] = info
		return info
	else:
		tidLineage[tid] = "|".join(lineage)
		return "|".join(lineage)

def _getTaxDepth( taxID ):
	taxID = taxid2mergedTid(taxID)
	return taxDepths[taxID]

def _getTaxName( taxID ):	
	taxID = taxid2mergedTid(taxID)
	return taxNames[taxID]

def _getTaxParent( taxID ):
	taxID = taxid2mergedTid(taxID)
	return taxParents[taxID]

def _getTaxRank( taxID ):
	taxID = taxid2mergedTid(taxID)
	return taxRanks[taxID]

#def loadStrainName( custom_taxonomy_file ):
#	try:
#		with open(custom_taxonomy_file, 'r') as f:
#			for line in f:
#				temp = line.rstrip('\r\n').split('\t')
#				tid = temp[4]
#				if "." in tid:
#					parent, sid = tid.split('.')
#					if not parent in taxNames: continue
#					taxParents[tid] = parent
#					taxDepths[tid] = str(int(taxDepths[parent]) + 1)
#					taxRanks[tid] = "no rank"
#					taxNames[tid] = temp[0]
#					taxNumChilds[tid] = 1
#					if parent in taxNumChilds: del taxNumChilds[parent]
#		f.close()
#	except IOError:
#		_die( "Failed to open custom taxonomy file: %s.\n" % custom_taxonomy_file )

def loadRefSeqCatelog( refseq_catelog_file, seq_type="nc" ):
	try:
		if refseq_catelog_file.endswith(".gz"):
			#p = subprocess.Popen(["zcat", refseq_catelog_file], stdout = subprocess.PIPE)
			#f = io.StringIO(p.communicate()[0])
			#assert p.returncode == 0
			f = gzip.open( refseq_catelog_file, 'r' )
		else:
			f = open( refseq_catelog_file, 'r' )

		for line in f:
			temp = line.rstrip('\r\n').split('\t')
			acc = temp[2]
			if seq_type == "nc" and ( acc[1] == "P" or acc.startswith("NM_") or acc.startswith("NR_") or acc.startswith("XM_") or acc.startswith("XR_") ):
				continue
			else:
				accTid[acc] = temp[0]

		f.close()
	except IOError:
		_die( "Failed to open custom RefSeq catelog file: %s.\n" % refseq_catelog_file )

def loadTaxonomy( dbpath=taxonomyDir ):
	global taxonomyDir

	if dbpath:
		taxonomyDir = dbpath

	if DEBUG: sys.stderr.write( "[INFO] Open taxonomy files from: %s\n"% taxonomyDir )

	#parsed taxonomy tsv file
	taxonomy_file = taxonomyDir+"/taxonomy.tsv"
	cus_taxonomy_file = taxonomyDir+"/taxonomy.custom.tsv"
	merged_taxonomy_file = taxonomyDir+"/taxonomy.merged.tsv"

	#raw taxonomy dmp files from NCBI
	names_dmp_file = taxonomyDir+"/names.dmp"
	nodes_dmp_file = taxonomyDir+"/nodes.dmp"
	merged_dmp_file = taxonomyDir+"/merged.dmp"

	# try to load taxonomy from taxonomy.tsv
	if os.path.isfile( taxonomy_file ):
		if DEBUG: sys.stderr.write( "[INFO] Open taxonomy file: %s\n"% taxonomy_file )
		try:
			with open(taxonomy_file) as f:
				for line in f:
					tid, depth, parent, rank, name = line.rstrip('\r\n').split('\t')
					taxParents[tid] = parent
					taxDepths[tid] = depth
					taxRanks[tid] = rank
					taxNames[tid] = name
					if parent in taxNumChilds:
						taxNumChilds[parent] += 1
					else:
						taxNumChilds[parent] = 1
				f.close()
		except IOError:
			_die( "Failed to open taxonomy file: %s.\n" % taxonomy_file )
	else:
		try:
			# read name from names.dmp
			if DEBUG: sys.stderr.write( "[INFO] Open taxonomy name file: %s\n"% names_dmp_file )
			with open(names_dmp_file) as f:
				for line in f:
					tid, name, tmp, nametype = line.rstrip('\r\n').split('\t|\t')
					if not nametype.startswith("scientific name"):
						continue 
					taxNames[tid] = name
				f.close()
			
			# read taxonomy info from nodes.dmp
			if DEBUG: sys.stderr.write( "[INFO] Open taxonomy node file: %s\n"% nodes_dmp_file )
			with open(nodes_dmp_file) as f:
				for line in f:
					fields = line.rstrip('\r\n').split('\t|\t')
					tid = fields[0]
					parent = fields[2]
					taxParents[tid] = fields[1]
					taxDepths[tid] = taxDepths[parent]+1 if parent in taxDepths else 0 # could have potiential bug if child node is parsed before parent node.
					taxRanks[tid] = parent
					if parent in taxNumChilds:
						taxNumChilds[tid] += 1
					else:
						taxNumChilds[tid] = 1
				f.close()
		except IOError:
			_die( "Failed to open taxonomy files (taxonomy.tsv, nodes.dmp and names.dmp).\n" )

	# try to load custom taxonomy from taxonomy.custom.tsv
	if os.path.isfile( cus_taxonomy_file ):
		if DEBUG: sys.stderr.write( "[INFO] Open custom taxonomy node file: %s\n"% cus_taxonomy_file)
		try:
			with open(cus_taxonomy_file) as f:
				for line in f:
					tid, depth, parent, rank, name = line.rstrip('\r\n').split('\t')
					taxParents[tid] = parent
					taxDepths[tid] = depth
					taxRanks[tid] = rank
					taxNames[tid] = name
					if parent in taxNumChilds:
						taxNumChilds[parent] += 1
					else:
						taxNumChilds[parent] = 1
				f.close()
		except IOError:
			_die( "Failed to open custom taxonomy file: %s.\n" % cus_taxonomy_file )

	#try to load merged taxids
	if os.path.isfile( merged_taxonomy_file ):
		if DEBUG: sys.stderr.write( "[INFO] Open merged taxonomy node file: %s\n"% merged_taxonomy_file )
		with open(merged_taxonomy_file) as f:
			for line in f:
				mtid, tid = line.rstrip('\r\n').split('\t')
				taxMerged[mtid] = tid
			f.close()
	elif os.path.isfile( merged_dmp_file ):
		if DEBUG: sys.stderr.write( "[INFO] Open merged taxonomy node file: %s\n"% merged_dmp_file )
		with open(merged_dmp_file) as f:
			for line in f:
				fields = line.rstrip('\r\n').split('\t')
				mtid = fields[0]
				tid = fields[2]
				taxMerged[mtid] = tid
			f.close()

	if DEBUG: sys.stderr.write( "[INFO] Done parsing taxonomy.tab (%d taxons loaded)\n" % len(taxParents) )

	if taxParents["2"] == "1":
		_die( "Local taxonomy database is out of date." )

##########################
##  Internal functions  ##
##########################

class _autoVivification(dict):
	"""Implementation of perl's autovivification feature."""
	def __getitem__(self, item):
		try:
			return dict.__getitem__(self, item)
		except KeyError:
			value = self[item] = type(self)()
			return value

def _die( msg ):
	sys.exit(msg)

def _checkTaxonomy():
	if not len(taxParents):
		_die("Taxonomy not loaded. \"loadTaxonomy()\" must be called first.\n")

if __name__ == '__main__':
	#loading taxonomy
	loadTaxonomy( sys.argv[1] if len(sys.argv) else "" )

	print("Enter acc/taxid:")

	for inid in sys.stdin:
		inid = inid.rstrip('\r\n')

		if inid[0] in "1234567890":
			taxid = inid
		else:
			taxid = acc2taxid( inid )
			print( "acc2taxid( %s ) => %s"   % (inid, taxid) )

		if taxid:
			print( "taxid2name( %s )                 => %s" % (taxid, taxid2name(taxid)) )
			print( "taxid2rank( %s )                 => %s" % (taxid, taxid2rank(taxid)) )
			print( "taxid2type( %s )                 => %s" % (taxid, taxid2type(taxid)) )
			print( "taxid2depth( %s )                => %s" % (taxid, taxid2depth(taxid)) )
			print( "taxid2parent( %s )               => %s" % (taxid, taxid2parent(taxid)) )
			print( "taxidIsLeaf( %s )                => %s" % (taxid, taxidIsLeaf(taxid)) )
			print( "taxid2nearestMajorTaxida( %s )   => %s" % (taxid, taxid2nearestMajorTaxid(taxid)) )
			print( "taxid2nameOnRank( %s, 'genus')   => %s" % (taxid, taxid2nameOnRank(taxid, "genus")) )
			print( "taxid2taxidOnRank( %s, 'genus')  => %s" % (taxid, taxid2taxidOnRank(taxid, "genus")) )
			print( "taxid2nameOnRank( %s, 'phylum')  => %s" % (taxid, taxid2nameOnRank(taxid, "phylum")) )
			print( "taxid2taxidOnRank( %s, 'phylum') => %s" % (taxid, taxid2taxidOnRank(taxid, "phylum")) )
			print( "taxid2lineage( %s )              => %s" % (taxid, taxid2lineage(taxid)) )
			print( "taxid2lineageDICT( %s )          => %s" % (taxid, taxid2lineageDICT(taxid)) )
			print( "taxid2fullLineage( %s )          => %s" % (taxid, taxid2fullLineage(taxid)) )
			print( "taxid2fullLinkDict( %s )         => %s" % (taxid, taxid2fullLinkDict(taxid)) )
		else:
			print( "No taxid found.\n" )

		print("Enter acc/taxid:")
