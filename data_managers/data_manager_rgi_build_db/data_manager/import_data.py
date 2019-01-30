import argparse
import datetime
import json
import os
import shutil
import sys
import tarfile
import urllib.request, urllib.error, urllib.parse
import zipfile
import logging

path = os.path.join(os.getcwd(), 'rgi-database') 
data_path = path

level = logging.WARNING
logger = logging.getLogger(__name__)
logger.setLevel(level)

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

def write_fasta_from_json():
		'''Creates a fasta file from card.json file.'''
		if os.path.isfile(os.path.join(path, 'proteindb.fsa')):
			return
		else:
			try:
				with open(os.path.join(data_path, 'card.json'), 'r') as jfile:
					j = json.load(jfile)
			except Exception as e:
				logger.error(e)
				exit()

			with open(os.path.join(path, 'proteindb.fsa'), 'w') as fout:
				for i in j:
					if i.isdigit():
		            	# model_type: protein homolog model
						if j[i]['model_type_id'] == '40292':
							try:
								pass_bit_score = j[i]['model_param']['blastp_bit_score']['param_value']
							except KeyError:
								logger.warning('No bitscore for model (%s, %s). RGI will omit this model and keep running.' \
									% (j[i]['model_id'], j[i]['model_name']))
								logger.info('Please let the CARD Admins know! Email: card@mcmaster.ca')
							else:
								try:
									for seq in j[i]['model_sequences']['sequence']:
										fout.write('>%s_%s | model_type_id: 40292 | pass_bitscore: %s | %s\n' % (i, seq, pass_bit_score, j[i]['ARO_name']))
										fout.write('%s\n' %(j[i]['model_sequences']['sequence'][seq]['protein_sequence']['sequence']))
								except Exception as e:
									logger.warning('No model sequences for model (%s, %s). RGI will omit this model and keep running.' \
										% (j[i]['model_id'], j[i]['model_name']))
									logger.info('Please let the CARD Admins know! Email: card@mcmaster.ca')


		            	# model_type: protein variant model
						elif j[i]['model_type_id'] == '40293':
							try:
								pass_bit_score = j[i]['model_param']['blastp_bit_score']['param_value']
							except KeyError:
								logger.warning('No bitscore for model (%s, %s). RGI will omit this model and keep running.' \
									% (j[i]['model_id'], j[i]['model_name']))
								logger.info('Please let the CARD Admins know! Email: card@mcmaster.ca')
							else:
								try:
									snpList = [j[i]['model_param']['snp']['param_value'][k] for k in j[i]['model_param']['snp']['param_value']]
								except Exception as e:
									logger.warning('No snp for model (%s, %s). RGI will omit this model and keep running.' \
										% (j[i]['model_id'], j[i]['model_name']))
									logger.info('Please let the CARD Admins know! Email: card@mcmaster.ca')

								try:
									for seq in j[i]['model_sequences']['sequence']:
										fout.write('>%s_%s | model_type_id: 40293 | pass_bit_score: %s | SNP: %s | %s\n' \
											% (i, seq, pass_bit_score, ','.join(snpList), j[i]['ARO_name']))
										fout.write('%s\n' % (j[i]['model_sequences']['sequence'][seq]['protein_sequence']['sequence']))
								except Exception as e:
									logger.warning('No model sequences for model (%s, %s). RGI will omit this model and keep running.' \
										% (j[i]['model_id'], j[i]['model_name']))
									logger.info('Please let the CARD Admins know! Email: card@mcmaster.ca')

		            	# model_type: protein overexpression model
						elif j[i]['model_type_id'] == '41091':
							try:
								pass_bit_score = j[i]['model_param']['blastp_bit_score']['param_value']
							except KeyError:
								logger.warning('No bitscore for model (%s, %s). RGI will omit this model and keep running.' \
									% (j[i]['model_id'], j[i]['model_name']))
								logger.info('Please let the CARD Admins know! Email: card@mcmaster.ca')
							else:
								try:
									snpList = [j[i]['model_param']['snp']['param_value'][k] for k in j[i]['model_param']['snp']['param_value']]
								except Exception as e:
									logger.warning('No snp for model (%s, %s). RGI will omit this model and keep running.' \
										% (j[i]['model_id'], j[i]['model_name']))
									logger.info('Please let the CARD Admins know! Email: card@mcmaster.ca')

								try:
									for seq in j[i]['model_sequences']['sequence']:
										fout.write('>%s_%s | model_type_id: 41091 | pass_bit_score: %s | SNP: %s | %s\n' \
											% (i, seq, pass_bit_score, ','.join(snpList), j[i]['ARO_name']))
										fout.write('%s\n' % (j[i]['model_sequences']['sequence'][seq]['protein_sequence']['sequence']))
								except Exception as e:
									logger.warning('No model sequences for model (%s, %s). RGI will omit this model and keep running.' \
										% (j[i]['model_id'], j[i]['model_name']))
									logger.info('Please let the CARD Admins know! Email: card@mcmaster.ca')

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
    write_fasta_from_json()
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
