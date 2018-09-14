#M.-J. H. Halsør (2014).
"""Script for converting data from EMBL input file
to a svg picture. It is the combination of two scripts:
embl2cds delivers text_output, a table with CDS data
cds2svg delivers the svg file according to the data from text_output

Usage: embl2svg.py ifile.embl base_start base_stop
Output: embl2svg.svg
"""

#from embl2cds.py
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

#functions from cds2svg.py
def set_width(cds_table):
    """Extracts start and stop bases from input file.

    Calculates the length of sequence to be displayed, thus allowing horizontal size"""

    start=cds_table[1].split(';')[0]
    stop=cds_table[-1].split(';')[1]
    svg_width=(int(stop)-int(start))//5
    return(start,stop,svg_width)

def set_coordinates(cds_table,organism,start):
    """Extracts required coordinates for svg paths"""

    arrow_data=[]
    x0=len(organism)*19
    xshift=int(start)//6-x0
    for i,p in enumerate(cds_table):
        cds_start=p.split(';')[0]
        cds_stop=p.split(';')[1]
        name=p.split(';')[2]
        path_length=(int(cds_stop)-int(cds_start))//6
        x1=int(cds_start)//6-xshift
        #check for strand oritentation
        if len(p.split(';'))==4:#reverse
            y1=900
            x2,y2=25,-50
            x3,y3=path_length-25,0
            x4,y4=0,100
            x5,y5=-(path_length-25),0
        else:
            y1=850
            x2,y2=path_length-25,0
            x3,y3=25,50
            x4,y4=-25,50
            x5,y5=-(path_length-25),0
        arrow_data.append([name,x1,y1,x2,y2,x3,y3,x4,y4,x5,y5])
    return(arrow_data)

def draw_arrow (arrow_data):
    """Draws arrows from coordinates and labels according to CDS name

    Example of arrow_data:
    [[...],['name',x1,y1,x2,y2,x3,y3,x4,y4,x5,y5],[...]]"""

    svg_arrow=''
    path_style='   fill="none" stroke="black" stroke-width="3"/>\n'
    for i,p in enumerate(arrow_data):
        arrow_path='\n<path d="M'+str(p[1])+','+str(p[2])+' l'+str(p[3])+','+str(p[4])+' l'+str(p[5])+','+str(p[6])+' l'+str(p[7])+','+str(p[8])+' l'+str(p[9])+','+str(p[10])+'z" '
        arrow_text='<text x="'+str(p[1]+25)+' " y="840" transform="rotate(-45 '+str(p[1]+25)+' 840)" font-size="40"> '+p[0]+'</text>\n'
        svg_arrow=svg_arrow+arrow_path+path_style+arrow_text
    return(svg_arrow)


def code_svg (svg_width,organism,svg_arrow):
    """Writes the code for the svg output file in a text string."""

    svg_code="""<?xml version="1.0" standalone="no"?>
<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.1//EN" 
    "http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd">\n"""+'<svg width="'+str(svg_width)+'" height="1000"\n'+'    xmlns="http://www.w3.org/2000/svg" version="1.1">\n'+'<title>CDS organisation-'+organism+'('+start+'..'+stop+')</title>\n'+'\n<text x="10" y="910" font-size="40" font-style="italic">'+organism+'</text>\n'+svg_arrow+'</svg>'
    return(svg_code)

#from generate_plot.py
def main():
    parser = argparse.ArgumentParser(description='Create example SVG')
    parser.add_argument('-o', '--output', required=True, help='Output file')
    parser.add_argument('-i', '--input', required=True, help='Input file')
    args = parser.parse_args()

    with open(args.output, 'w') as f:
        f.write(svg_data)

#calling functions from generate_plot.py
if __name__ == "__main__":
    main()


#calling functions-from embl2cds.py

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

#calling functions from cds2svg.py
ifile=open('embl2cds.txt','r')
cds_table=ifile.readlines()
ifile.close()

start,stop,svg_width=set_width(cds_table)
organism=cds_table[0].replace('\n','')
del cds_table[0]

arrow_data=set_coordinates(cds_table,organism,start)
svg_arrow=draw_arrow(arrow_data)
svg_code=code_svg(svg_width,organism,svg_arrow)

ofile=open('embl2cds.svg','w')
ofile.write(svg_code)
ofile.close()