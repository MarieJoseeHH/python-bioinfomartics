#M.-J. H. Halsoer (2014).
"""Script for converting the cds in base range
to a picture in vector graphics.

#Usage: python cds2svg.py filename
Output:filename.svg
"""

"""Example of input file:
Aliivibrio_R8-70_concat
101249;102085;Putative uncharacterized protein;reverse
102173;102589;Putative uncharacterized protein;reverse
102742;103962;Arginine deiminase - Arginine dihydrolase
104058;105062;Ornithine carbamoyltransferase chain I

"""
import sys

def set_width(cds_table,organsim):
    """Extracts start and stop bases from input file.

    Calculates the length of sequence to be displayed, thus allowing horizontal size"""

    start=cds_table[1].split(';')[0]
    stop=cds_table[-1].split(';')[1]
    svg_width=abs(int(stop)-int(start))//4
    return(start,stop,svg_width)

def set_coordinates(cds_table,organism,start):
    """Extracts required coordinates for svg paths"""

    arrow_data=[]
    x0=len(organism)*20
    xshift=int(start)//6-x0
    for i,p in enumerate(cds_table):
        cds_start=p.split(';')[0]
        cds_stop=p.split(';')[1]
        name=p.split(';')[2]
        path_length=(int(cds_stop)-int(cds_start))//6
        x1=int(cds_start)//6-xshift
        #check for strand oritentation
        if len(p.split(';'))==4:#reverse
            y1=200
            x2,y2=25,-50
            x3,y3=path_length-25,0
            x4,y4=0,100
            x5,y5=-(path_length-25),0
        else:
            y1=150
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
        arrow_text='<text x="'+str(p[1]+25)+' " y="140" transform="rotate(-45 '+str(p[1]+25)+' 140)" font-size="40"> '+p[0]+'</text>\n'
        svg_arrow=svg_arrow+arrow_path+path_style+arrow_text
    return(svg_arrow)


def code_svg (svg_width,organism,svg_arrow):
    """Writes the code for the svg output file in a text string."""

    svg_code="""<?xml version="1.0" standalone="no"?>
<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.1//EN" 
    "http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd">\n"""+'<svg width="'+str(svg_width)+'" height="400"\n'+'    xmlns="http://www.w3.org/2000/svg" version="1.1">\n'+'<title>CDS organisation-'+organism+'('+start+'..'+stop+')</title>\n'+'\n<text x="10" y="210" font-size="40" font-style="italic">'+organism+'</text>\n'+svg_arrow+'</svg>'
    return(svg_code)

inputfile=sys.argv[1]
ifile=open(inputfile,'r')
cds_table=ifile.readlines()
ifile.close()

organism=cds_table[0].replace('\n','')
start,stop,svg_width=set_width(cds_table,organism)
del cds_table[0]

arrow_data=set_coordinates(cds_table,organism,start)
svg_arrow=draw_arrow(arrow_data)
svg_code=code_svg(svg_width,organism,svg_arrow)

outputfile=inputfile.split('.')[0]+'.svg'
ofile=open(outputfile,'w')
ofile.write(svg_code)
ofile.close()
