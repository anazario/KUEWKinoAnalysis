import os
import fileinput
import sys
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-p", "--path", dest="directory",
                  help="Specify input directory containing the .root files.", metavar="PATH")
parser.add_option("-t", "--tag", dest="tag",
                  help="Specify name of sample to make list and output .txt file.", metavar="TAG")
parser.add_option("-y", "--year", dest="year",
                  help="Specify year of produced sample.", metavar="YEAR")

(options, args) = parser.parse_args()

directory = options.directory
tag = options.tag
yeartag = options.year

if not directory:
    sys.exit("You need to specify the directory! (See help).")
if not tag:
    sys.exit("You need to specify sample name.")

for folder in os.listdir(directory):
    if folder.startswith(tag) and not folder.endswith('.root'):
        for rootfile in os.listdir(directory+folder):
            if not yeartag:
                 with open('samples/Efficiency/'+tag+'.txt', 'a') as ofile:
                     ofile.write(directory+folder+'/'+rootfile + os.linesep)
            else:
                with open('samples/Efficiency/'+tag+'_'+yeartag+'.txt', 'a') as ofile:
                    ofile.write(directory+folder+'/'+rootfile + os.linesep)
