__author__ = 'Doyel'

import argparse
from xml.etree import ElementTree
from pprint import pprint;
import cPickle as pickle;
import mimetypes
import os
import glob
import urllib2
import urlparse
import shutil
import tarfile, sys
from tempfile import mkstemp
import subprocess
import fileinput
from Utilities import read_config_file

# python XML_Articles_tar_gz.py Eric_SRI_4TestPMIDs.txt ./xmlFTPDetails/ ./downloaded_articles

def main_body():
    global parent_location; global txt_location; global json_location
    parent_location, txt_location, json_location = read_config_file("./config.cfg")

    parser = argparse.ArgumentParser(prog='XML Parser', usage='XML_Articles_tar_gz.py <file_list> <dir> <download_dir>', description='Script to parse XML files and download tar.gz files from ftp')
    parser.add_argument('PMIDfilelist', help='File listing the PubMedIDs that need to be downloaded')
    parser.add_argument('XMLfiledir', help='Dir where the XML files are located')
    parser.add_argument('download_dir', help='Dir where the articles will be downloaded from the internet')

    args = parser.parse_args()

    pmid_pmcid = pickle.load(open('pmid_pmcid_Extras.p', 'rb'))
    # txt_location = "./"     # Temporary: Please comment out after Eric SRI deliverables are sent

    pmcid_openaccess = []; pmid_openaccess = []
    with open(os.path.join(txt_location, args.PMIDfilelist), 'rt') as f1:
        for pubmedid in f1:
            pmcid = pmid_pmcid[pubmedid.strip()][3:]
            fname_string = '%s.xml' % (os.path.join(args.XMLfiledir,'File_Details_PMC'+pmcid.strip()))
            fname = glob.glob(fname_string)
            if len(fname) == 0:
                print "No xml path file found for PMCID: ", pmcid
                continue
            my_file = fname[0].strip()
            #print my_file
            with open(my_file, 'rt') as f:
                xml_tree = ElementTree.parse(f)
                root = xml_tree.getroot() 
                if root.find("./error") != None: 
                    err = root.find("./error")
                    err_c = err.get("code")
                    if err_c == "idIsNotOpenAccess":
                        continue   
                    else:
                        print "There is some error in the xml file for PMCID: ", pmcid, " Please check!"
                        exit()
                node = root.find("./request")
                current_pmcid = node.get('id'); #print current_pmcid
                pmcid_openaccess.append(current_pmcid)
                pmid_openaccess.append(pubmedid.strip())
                node = root.find("./records/record")
                tar_exists = False;
                for each_link in node:
                    if each_link.get('format') == "tgz": 
                        tar_exists = True;
                        tar_location = each_link.get('href')
                        tarfilename = filename_from_url(tar_location.strip())
                        #print tarfilename; 
                        break
                if not tar_exists:
                    print "No tgz file found for PMC id: ", current_pmcid                    
                    continue;
                print tar_location
                download_file_locally(tar_location.strip(), args.download_dir)
                untar(os.path.join(args.download_dir, tarfilename), args.download_dir)
                os.rename(os.path.join(args.download_dir,tarfilename[:-7]), os.path.join(args.download_dir, current_pmcid))

                #Check whether a nxml file is located in the newly untarred directory
                cwd = os.getcwd()
                os.chdir(os.path.join(args.download_dir, current_pmcid))
                fname_string1 = '*.nxml'
                fname1 = glob.glob(fname_string1)
                #print fname1; 
                if len(fname1) == 0:
                    print "No nxml path file found under: ", os.getcwd()
                    os.chdir(cwd)
                    continue
                elif len(fname1) > 1:
                    print "More than one nxml file found under: ", os.getcwd()
                    #print "Using the first one"
                #os.chdir(cwd)

                ###### Converting the NXML file to HTML #######
                found_str1 = False
                my_file1 = fname1[0].strip() 

                found_str = replace(os.path.join("./", my_file1), "JATS-archivearticle1.dtd", os.path.join(parent_location, "jats-archiving-dtd-1.0", "JATS-archivearticle1.dtd"))
                #found_str = find_string(my_file1, "JATS-archivearticle1.dtd")
                if not found_str:          
                    found_str1 = replace(os.path.join("./", my_file1), "archivearticle.dtd", os.path.join(parent_location, "archive-interchange-dtd", "archivearticle.dtd"))
                    #found_str1 = find_string(my_file1, "archivearticle.dtd")                    
                    if not found_str1:   
                        print "DTD file was NOT recognized for ", current_pmcid
                        os.chdir(cwd)
                        exit() #continue
                    else:
                        shutil.copytree(os.path.join(parent_location, "archive-interchange-dtd", "iso9573-13"), "./iso9573-13")
                        shutil.copytree(os.path.join(parent_location, "archive-interchange-dtd", "iso8879"), "./iso8879")

                shutil.copyfile(os.path.join(parent_location, "JATSPreviewStylesheets-master", "xslt", "main", "jats-html-NOSPAN.xsl"), "./jats-html.xsl")
                shutil.copyfile(os.path.join(parent_location, "JATSPreviewStylesheets-master", "jats-preview.css"), "./jats-preview.css")

                # If an article does NOT have a single jpg file, then use /var/www/py/Code/SRI_55000/JATSPreviewStylesheets-master/xslt/main/jats-html_ORIG.xsl for the conversion
                jpg_string = '*.jpg'     
                jpg_fname = glob.glob(jpg_string)            
                if len(jpg_fname) == 0:
                    os.remove("./jats-html.xsl")
                    shutil.copyfile(os.path.join(parent_location, "JATSPreviewStylesheets-master", "xslt", "main", "jats-html_ORIG.xsl"), "./jats-html.xsl")

                cmd = "java -jar $XSLT_HOME/saxon9he.jar -xsl:jats-html.xsl -s:"+my_file1+" > " + current_pmcid +".html"
                #print cmd
                p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                for line in p.stderr.readlines():   # Do NOT comment this for statement
                    #print line.strip()
                    pass
                os.remove("./jats-html.xsl")

                if found_str1:
                    shutil.rmtree("./iso9573-13")
                    shutil.rmtree("./iso8879")               
                if os.stat(current_pmcid +".html").st_size == 0:
                    print "The output html file for ", current_pmcid, " was of 0 bytes"
                else:
                    print "Processed and converted to HTML: ", current_pmcid
                os.chdir(cwd)
            
    pmid_fn1 = open("PMCIDS_with_Extras_OpenAccess.txt", 'w');
    for pmc in pmcid_openaccess:    
        pmid_fn1.write(pmc[3:]+'\n')
    pmid_fn1.close()

    pmid_fn1 = open("PubMedIDS_with_Extras_OpenAccess.txt", 'w');
    for pm in pmid_openaccess:
        pmid_fn1.write(pm + '\n')
    pmid_fn1.close()
            


def filename_from_url(url):
    return os.path.basename(urlparse.urlsplit(url)[2])


def download_file(url):
    """Create an urllib2 request and return the request plus some useful info"""
    name = filename_from_url(url)
    r = urllib2.urlopen(urllib2.Request(url))
    info = r.info()
    if 'Content-Disposition' in info:
        # If the response has Content-Disposition, we take filename from it
        name = info['Content-Disposition'].split('filename=')[1]
        if name[0] == '"' or name[0] == "'":
            name = name[1:-1]
    elif r.geturl() != url:
        # if we were redirected, take the filename from the final url
        name = filename_from_url(r.geturl())
    content_type = None
    if 'Content-Type' in info:
        content_type = info['Content-Type'].split(';')[0]
    # Try to guess missing info
    if not name and not content_type:
        name = 'unknown'
    elif not name:
        name = 'unknown' + mimetypes.guess_extension(content_type) or ''
    elif not content_type:
        content_type = mimetypes.guess_type(name)[0]
    return r, name, content_type


def download_file_locally(url, dest):
    req, filename, content_type = download_file(url)
    dest = os.path.join(dest, filename)
    with open(dest, 'wb') as f:
        shutil.copyfileobj(req, f)
    req.close()


def untar(fname, tardir):
    if (fname.endswith("tar.gz")):
        tar = tarfile.open(fname)
        tar.extractall(path=tardir)
        tar.close()
        #print "Extracted in Current Directory"
    else:
        print "Not a tar.gz file: '%s '" % fname
    
    
def replace(myfile,searchExp,replaceExp):
    found_str = False
    for line in fileinput.input(myfile, inplace=1):
        if searchExp in line:
            line = line.replace(searchExp,replaceExp)
            found_str = True
        print line,
    return found_str


if __name__ == "__main__":
    main_body()

