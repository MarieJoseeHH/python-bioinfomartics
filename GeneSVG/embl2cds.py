#M.-J. H. Halsør (2014).
"""Script for converting data from EMBL input file
to a CDS table.

Usage: embl2cds.py ifile.embl base_start base_stop
Output: embl2cds.txt
"""

import re
import sys
args = sys.argv
def parse_embl(ifile):
    """Return list of CDS for a EMBL file.

    Each CDS is returned as a list [location,product]
    where location and  product are strings.

    Assumes input file is opened.
    """
    organism=''
    info=[]
    line=ifile.readline()
    while line != '':
        #getting organism name
        if re.search('origid',line)!= None:
            organism=line.split('"')[1]
        #getting CDS data
        if line.find('  CDS  ') != -1:
            location=line.split()[2]
            while "product" not in line:
                line=ifile.readline()
            product=line.split('"')[1]
            line=ifile.readline()
            #Multiline product name
            if line.find('/')==-1:
                product=product+line.split()[1]
                line=ifile.readline()
            product=product.replace('\n',' ').strip('"')
            line=ifile.readline()
            info.append([location,product])
        line=ifile.readline()
    return(info,organism)

def set_strand(cds):
    """Marks reverse strand cds.

    Expects a list of lists [[location,product],[location,product]],
    where location is of the from 'complement(2606..4186)',
    and product of the form 'Peptide chain release factor 3'

    Returns a modified version of list with additional 'reverse' element."""

    for i,p in enumerate(cds):
        #print(p[0])
        if p[0].find('complement') !=-1:
            p.append('reverse')
    return(cds)

def base_positions(cds):
    """Extracts base positions and strand for cds.

    Expects a list of lists [[location,product],[location,product,'reverse']],
    where location is of the from 'complement(2606..4186)',
    and product of the form 'Peptide chain release factor 3'

    Returns a modified version of list."""
    for i,p in enumerate(cds):
        positions=re.findall('([0-9]+)',p[0])
        p.insert(0,positions[0])
        p.insert(1,positions[-1])
        del p[2]
    return(cds)

def base_range(cds,base_start,base_stop):
    """Limits cds to base range."""
    portion=[]
    for i,p in enumerate(cds):
        if (int(p[0])>=base_start and int(p[1])<=base_stop):
            portion.append(p)
            cds=portion
    return(cds)

#calling functions

#getting information from embl file
ifile=open(args[1],'rU')
info,organism=parse_embl(ifile)
ifile.close()
#extracting cds strand
cds=set_strand(info)

#extracting cds base positions
cds=base_positions(cds)

#limiting to the user range
base_start=int(args[2])
base_stop= int(args[3])
cds=base_range(cds,base_start,base_stop)

#writing all as text in output file.
#first line is organism
text_output=organism
for i,p in enumerate(cds):
    text_output=text_output+'\n'+';'.join(p)
ofile=open('embl2cds.txt','w')
ofile.write(text_output)
ofile.close()
#print(text_output)