import os
import sys
from optparse import OptionParser

#options
parser = OptionParser()
parser.add_option("-i", "--ipath", dest="directory",
                  help="Specify input directory containing the dataset .txt files.", metavar="IPATH")

parser.add_option("-o", "--opath", dest="output", default=".",
                  help="Specify output directory where .txt and .list files are going to be written.", metavar="OPATH")

parser.add_option("-v", "--ver", dest="version", default='v5',
                  help="Specify version of files in the dasgo client.", metavar="VER")


(options, args) = parser.parse_args()

directory = options.directory
output = options.output
version = options.version

if not directory:
    sys.exit("You need to specify the directory! (See help).")

for filename in os.listdir(directory):

    if filename.endswith(".txt"):

        if "2016" in filename:
            yeartag = 'Summer16'
        elif "2017" in filename:
            yeartag = 'Fall17'
        elif "2018" in filename:
            yeartag = 'Autumn18'
        else:
            continue

        print("Producing {year} lists for ".format(year=yeartag)+os.path.join(directory,filename))
        with open(os.path.join(directory,filename), 'r') as handle:
            for line in handle:
                if 'SMS' in line and yeartag=='Autumn18':
                    version='v4'
                #print("    Searching DAS for dataset: "+line)
                os.system('dasgoclient -query="dataset=/{ds_name}/*{year}*{ver}*/*NANO*" >> {ds_name}_dataset.txt'.format(ds_name=line.rstrip('\n'),year=yeartag,ver=version))
                with open('{ds_name}_dataset.txt'.format(ds_name=line.rstrip('\n')), 'r') as datasetlist:
                    for dataset in datasetlist:
                        #print("        Producing ROOT file list for dataset: "+dataset)
                        os.system('dasgoclient -query="file dataset={ds_loc}" >> {ds_name}.txt'.format(ds_loc=dataset.rstrip('\n'),ds_name=line.rstrip('\n')))

        #print("Making lists for {year}".format(year=yeartag))
        os.system('rm *dataset.txt')
        if 'SMS' in filename:
            outpath='{opath}/{year}_102X_SMS'.format(opath=output, year=yeartag)
        else:
            outpath='{opath}/{year}_102X'.format(opath=output, year=yeartag)
        if not os.path.exists(outpath):                           
            os.system('mkdir -p {opath}'.format(opath=outpath))
        os.system('mv *.txt {opath}'.format(opath=outpath))
        os.system('python scripts/addPath.py -p {opath}'.format(opath=outpath))
