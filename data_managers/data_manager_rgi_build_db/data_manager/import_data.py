import argparse
import datetime
import json
import os
import shutil
import sys
import tarfile
import urllib.request, urllib.error, urllib.parse
import zipfile
import rgi.database.write_fasta_from_json

path = os.path.join(os.getcwd(), 'rgi-database') 
data_path = path

def url_download(url, workdir):
    file_path = os.path.join(workdir, 'download.dat')
    if not os.path.exists(workdir):
        os.makedirs(workdir)
    src = None
    dst = None
    try:
        req = urllib.request.Request(url)
        src = urllib.request.urlopen(req)
        dst = open(file_path, 'wb')
        while True:
            chunk = src.read(2**10)
            if chunk:
                dst.write(chunk)
            else:
                break
    except Exception as e:
        print(str(e), file=sys.stderr)
    finally:
        if src:
            src.close()
        if dst:
            dst.close()
    if tarfile.is_tarfile(file_path):
        fh = tarfile.open(file_path, 'r:*')
    elif zipfile.is_zipfile(file_path):
        fh = zipfile.ZipFile(file_path, 'r')
    else:
        return
    # extract only one file : card.json
    for member in fh.getmembers():
        if member.isreg():  # skip if the TarInfo is not files
            member.name = os.path.basename(member.name) # remove the path by reset it
            if member.name == 'card.json':
                print('[import_data] extracting file: {}'.format(str(member.name)))
                fh.extract(member.name,workdir)
    os.remove(file_path)

def checkKeyExisted(key, my_dict):
    try:
        nonNone = my_dict[key] is not None
    except KeyError:
        nonNone = False
    return nonNone

def data_version():
    data_version = ''
    with open(os.path.join(data_path, 'card.json')) as json_file:
        json_data = json.load(json_file)
        for item in list(json_data.keys()):
            if item == '_version':
                data_version = json_data[item]
    json_file.close()
    return data_version

def makeBlastDB():
    if os.path.isfile(os.path.join(path, 'proteindb.fsa')) == True:
        print('[import_data] create blast DB.')
        os.system('makeblastdb -in {}/proteindb.fsa -dbtype prot -out {}/protein.db > /dev/null 2>&1'.format(path, path))

def makeDiamondDB():
    if os.path.isfile(os.path.join(path, 'proteindb.fsa')) == True:
        print('[import_data] create diamond DB.')
        os.system('diamond makedb --quiet --in {}/proteindb.fsa --db {}/protein.db'.format(path, path))
        
def _main(args):
    if not os.path.exists(path):
        print('[import_data] mkdir: {}'.format(path))
        os.makedirs(path)
    print('[import_data] path: {}'.format(path))
    print(args)

    if args.url == None:
        url = 'https://card.mcmaster.ca/latest/data'
    else:
        url = args.url
    print('[import_data] url: {}'.format(url))
    workdir = os.path.join(os.getcwd(), 'rgi-database')
    print('[import_data] working directory: {}'.format(workdir))
    url_download(url, workdir)
    rgi.database.write_fasta_from_json()
    makeBlastDB()
    makeDiamondDB()
    version = data_version()
    print('[import_data] data version: {}'.format(version))
    return version

def run():
    parser = argparse.ArgumentParser(description='Create data manager json.')
    parser.add_argument('--url', dest='url', action='store', help='Url for CARD data')
    args = parser.parse_args()
    _main(args)

if __name__ == '__main__':
    run()
