import os
import fileinput
import sys
from optparse import OptionParser

#options
parser = OptionParser()
parser.add_option("-p", "--path", dest="directory", 
                  help="Specify input directory containing the .txt files.", metavar="PATH")

(options, args) = parser.parse_args()

directory = options.directory

if not directory:
    sys.exit("You need to specify the directory! (See help).")

txtfiles=[]
isSMS=False

#loop over input files and add xrootd to root file path
for filename in os.listdir(directory):
    if 'SMS' in filename:
        isSMS=True
    if filename.endswith(".txt"):
        txtfiles.append(os.path.join('samples/NANO/{d}/'.format(d=directory.replace('/', '').replace('.','')),filename))
        for line in fileinput.input([os.path.join(directory, filename)], inplace=True):
            if 'root://cmsxrootd.fnal.gov/' in line:
                sys.stdout.write('{l}'.format(l=line))
            else:
                sys.stdout.write('root://cmsxrootd.fnal.gov/{l}'.format(l=line))
    else:
        continue

#make root text file list
txtfiles.sort()
#print txtfiles
if isSMS:
    with open('{d}.list'.format(d=directory.replace('/', '').replace('.','')), 'w') as filehandle:
        for listitem in txtfiles:
            filehandle.write(listitem+'\n')
else:
    with open('{d}_bkg.list'.format(d=directory.replace('/', '').replace('.','')), 'w') as filehandle:
        for listitem in txtfiles:
            filehandle.write(listitem+'\n')
