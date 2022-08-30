import os
import docker
import sys
import json
import logging
import argparse
import psutil
import shutil
import pandas as pd

from post_status import copy_file
from mkdir import mkdir
from post_status import post_url
from post_status import post_pid
from post_status import write_status
from post_status import write_status
from parse_fastqc_html import get_fqc
from parse_quast_html import parse_quast_html

def shell_db_path(db_path_file):
	with open(db_path_file) as h:
		db_path = json.load(h)
	tmp_str = []
	for i, path in db_path.items():
		tmp_str.append(f'echo "{i}={path}" >>/root/metaWRAP/bin/config-metawrap')

	shell_str = '\n'.join(tmp_str)
	return shell_str

def shell_readqc(data1, data2, name, out, threads):
	outdir = os.path.join(out,'READ_QC', name)
	mkdir(outdir)
	shell_str = f'metawrap read_qc -t {threads} -1 {data1} -2 {data2} -o {outdir}'
	tmp1 = os.path.join(outdir, 'final_pure_reads_1.fastq')
	tmp2 = os.path.join(outdir, 'final_pure_reads_2.fastq')
	return shell_str, [tmp1, tmp2]

def shell_clean(data_list, name_list, out):
	outdir = os.path.join(out, 'CLEAN_READS')
	mkdir(outdir)
	totaldir = os.path.join(out, 'CLEAN_READS_T')
	mkdir(totaldir)
	shell_str = ''
	clean_data_list = []
	for d, n in zip(data_list, name_list):
		tmp_1 = os.path.join(outdir, n+"_1.fastq")
		tmp_2 = os.path.join(outdir, n+"_2.fastq")
		shell_str += f'mv {d[0]} {tmp_1}\n'
		shell_str += f'mv {d[1]} {tmp_2}\n'
		clean_data_list.append([tmp_1, tmp_2])
	tmp_1 = os.path.join(totaldir, "ALL_READS_1.fastq")
	tmp_2 = os.path.join(totaldir, "ALL_READS_2.fastq")
	shell_str += f'cat {os.path.join(outdir, "*_1.fastq")} > {tmp_1}\n'
	shell_str += f'cat {os.path.join(outdir, "*_2.fastq")} > {tmp_2}\n'
	return shell_str, clean_data_list, [tmp_1, tmp_2]

def shell_assembly(data1, data2, out, threads, memory_G):
	outdir = os.path.join(out,'ASSEMBLY')
	shell_str = f'metawrap assembly -m {memory_G} -t {threads} -1 {data1} -2 {data2} -o {outdir}'
	contig_fa = os.path.join(outdir, 'final_assembly.fasta')
	return shell_str, contig_fa

def shell_kraken2(data_list, out, threads):	
	outdir = os.path.join(out, 'KRAKEN')
	shell_str = f'metawrap kraken2 -t {threads} -o {outdir} {" ".join(data_list)}'
	return shell_str

def shell_binning(data_list, data_contig, out, threads, memory_G):
	outdir = os.path.join(out, 'INITIAL_BINNING')
	shell_str = f'metawrap binning -m {memory_G} -t {threads} --metabat2 --maxbin2 --concoct -o {outdir} -a {data_contig} {" ".join(data_list)}'
	return shell_str, outdir

def shell_binRefinement(bin_dir, out, threads):
	outdir = os.path.join(out, 'BIN_REFINEMENT')
	metabat2_bins = os.path.join(bin_dir, 'metabat2_bins')
	maxbin2_bins = os.path.join(bin_dir, 'maxbin2_bins')
	concoct_bins = os.path.join(bin_dir, 'concoct_bins')
	shell_str = f'metawrap bin_refinement -c 50 -x 10 -t {threads} -o {outdir} -A {metabat2_bins} -B {maxbin2_bins} -C {concoct_bins} '
	return shell_str, outdir


def shell_blobology(data_list, data_contig, binRedir, out, threads):
	outdir = os.path.join(out, 'BLOBOLOGY')
	meta_bin = os.path.join(binRedir, 'metawrap_50_10_bins')
	shell_str = f'metawrap blobology -t {threads} -o {outdir} --bins {meta_bin} -a {data_contig} {" ".join(data_list)}'
	return shell_str

def shell_quant(data_list, data_contig, binRedir, out):
	outdir = os.path.join(out, 'QUANT_BINS')
	meta_bin = os.path.join(binRedir, 'metawrap_50_10_bins')
	shell_str = f'metawrap quant_bins -o {outdir} -b {meta_bin} -a {data_contig} {" ".join(data_list)}'
	return shell_str

def shell_reassemble(data1, data2, binRedir, out, threads, memory_G):
	outdir = os.path.join(out, 'BIN_REASSEMBLY')
	meta_bin = os.path.join(binRedir, 'metawrap_50_10_bins')
	shell_str = f'metawrap reassemble_bins -t {threads} -m {memory_G} -c 50 -x 10 -o {outdir} -1 {data1} -2 {data2} -b {meta_bin}'
	return shell_str, outdir

def shell_annotate(binRadir, out, threads):
	outdir = os.path.join(out, 'FUNCT_ANNOT')
	reassemble = os.path.join(binRadir, 'reassembled_bins')
	shell_str = f'metaWRAP annotate_bins -t {threads} -o {outdir} -b {reassemble}'
	return shell_str

def shell_classify(binRadir, out, threads):
	outdir = os.path.join(out, 'BIN_CLASSIFICATION')
	reassemble = os.path.join(binRadir, 'reassembled_bins')
	shell_str = f'metaWRAP classify_bins -t {threads} -o {outdir} -b {reassemble}'
	return shell_str


def get_kraken_shell(db_config, data_list, name_list, outdir, threads, memory_G, kraken=False):
	lst_total_shell =  []
	db_path_shell = shell_db_path(db_config)
	lst_total_shell.append(db_path_shell)
	lst_readqc_shell = []
	lst_readqc_file = []
	lst_readqc_name = []
	for fqs, name, in zip(data_list, name_list):
		if len(fqs)<2:
			logging.error(f'data should be paired fastq: {fqs}')
		else:
			tmp_fq = []
			tmp_shell = []
			if os.path.splitext(fqs[0])[1]=='.gz':
				tmp_fq.append(os.path.join(outdir, f'{name}_1.fastq'))
				tmp_shell = f'gunzip -c {fqs[0]} >{tmp_fq[0]}'
				lst_readqc_shell.append(tmp_shell)
			# if os.path.splitext(fqs[1])[1]=='.gz':
				tmp_fq.append(os.path.join(outdir, f'{name}_2.fastq'))
				tmp_shell = f'gunzip -c {fqs[1]} >{tmp_fq[1]}'
				lst_readqc_shell.append(tmp_shell)
			if not tmp_fq:
				tmp_fq = fqs

			tmp_shell, clean_data_list = shell_readqc(tmp_fq[0], tmp_fq[1], name, outdir, threads)
			lst_readqc_shell.append(tmp_shell)
			lst_readqc_file.append(clean_data_list)
			lst_readqc_name.append(name)

	if lst_readqc_shell:
		readqc_shell = '\n'.join(lst_readqc_shell)
		clean_shell, clean_data_list, total_clean_data_list = shell_clean(lst_readqc_file, lst_readqc_name, outdir)
		lst_total_shell.append(readqc_shell)
		lst_total_shell.append(clean_shell)
	if kraken:
		tmp_data_list = [i for item in clean_data_list+[total_clean_data_list] for i in item]
		kraken_shell = shell_kraken2(tmp_data_list, outdir, threads)
		lst_total_shell.append(kraken_shell)
	shell_str = '\n'.join(lst_total_shell)
	return shell_str, clean_data_list, total_clean_data_list

def get_assemble_shell(db_config, data_list, name_list, outdir, threads, memory_G, kraken=False):
	lst_total_shell = []

	clean_shell, clean_data_list, total_clean_data_list = get_kraken_shell(db_config, data_list, name_list, outdir, threads， memory_G)
	lst_total_shell.append(clean_shell)

	assemble_shell, assemble_contig = shell_assembly(total_clean_data_list[0], total_clean_data_list[1], outdir, threads, memory_G)
	lst_total_shell.append(assemble_shell)

	tmp_data_list = [i for item in clean_data_list+[[assemble_contig]] for i in item]
	kraken_shell = shell_kraken2(tmp_data_list, outdir, threads)
	lst_total_shell.append(kraken_shell)

	shell_str = '\n'.join(lst_total_shell)
	return shell_str, assemble_contig, clean_data_list, total_clean_data_list

def get_binning_shell(db_config, data_list, name_list, outdir, threads, memory_G, kraken=False):
	assemble_shell, assemble_contig, clean_data_list, total_clean_data_list = get_assemble_shell(db_config, data_list, name_list, outdir, threads, memory_G)
	clean_data_list = [i for item in clean_data_list for i in item]

	binning_shell, binning_dir = shell_binning(clean_data_list, assemble_contig, outdir, threads, memory_G)
	binRe_shell, binRe_dir = shell_binRefinement(binning_dir, outdir, threads)
	bolobo_shell = shell_blobology(clean_data_list, assemble_contig, binRe_dir, outdir, threads)
	quant_shell = shell_quant(clean_data_list, assemble_contig, binRe_dir, outdir)
	reass_shell,  binRa_dir = shell_reassemble(total_clean_data_list[0], total_clean_data_list[1], binRe_dir, outdir, threads, memory_G)
	class_shell = shell_classify(binRa_dir, outdir, threads)
	annot_shell = shell_annotate(binRa_dir, outdir, threads)

	lst_total_shell = [assemble_shell, binning_shell, binRe_shell, bolobo_shell, quant_shell, reass_shell, class_shell, annot_shell]

	shell_str = '\n'.join(lst_total_shell)
	return shell_str

def run_docker(image, volumes, cmd):
    client = docker.from_env()
    logging.info(cmd)
    path = ("/root/metaWRAP/bin:/root/miniconda3/envs/metawrap-env/bin:"
    	"/root/miniconda3/condabin:/root/anaconda3/bin:/root/anaconda3/condabin:"
    	"/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/root/bin:/root/metaWRAP/bin:/root/miniconda3/envs/metawrap-env/bin")
    return client.containers.run(image, cmd, volumes=volumes, environment=[f"PATH=PATH:{path}"],
    	remove=True, stdout=True, stderr=True)


def run_sh(sh_file, path_list):
    # image_name='mg_metawrap:220507'
    image_name='wnse/metawrap'
    cmd = f'sh -c "/usr/bin/sh {sh_file} >{sh_file}.log"'
    vols = set(path_list)
    vols = [f'{p}:{p}' for p in vols]    
    try:
        docker_out = run_docker(image_name, vols, cmd)
        logging.info(f'docker out: {docker_out}')
    except Exception as e:
        logging.error(e)


def get_file_path(file):
    outdict = {}
    if os.path.isfile(file):
    	file_key = os.path.split(file)[1]
    	file_key = str(file_key).replace('.','_')
    	logging.info(f'get_file_path:{file}')
    	return {file_key:file}
    elif os.path.isdir(file):
        try:
        	for f in os.listdir(file):
        		tmp_f = os.path.join(file, f)
        		outdict.update(get_file_path(tmp_f))
        except Exception as e:
            logging.error(e)
    else:
    	logging.info(f'get_file_path error')

    return outdict

def copy_file_batch(name_list, tmp_dir, outdir, step=2):
	outdict = {}
	res_file0 = {'READ_QC':[],'KRAKEN':['kronagram.html']}
	res_file1 = {'READ_QC':[],
		'KRAKEN':['kronagram.html'],
		'ASSEMBLY':['final_assembly.fasta','QUAST_out'],}
	res_file2 = {
		'READ_QC':[],
		'KRAKEN':['kronagram.html'],
		'ASSEMBLY':['final_assembly.fasta','QUAST_out'],
		'BLOBOLOGY':['blobplot_figures','blobplot_figures_only_binned_contigs','final_assembly.binned.blobplot','final_assembly.blobplot'],
		'QUANT_BINS':['bin_abundance_table.tab','quant_files'],
		'BIN_REASSEMBLY':['original_bins.stats','reassembled_bins','reassembled_bins.png','reassembled_bins.stats','reassembly_results.eps','reassembly_results.png'],
		'FUNCT_ANNOT':['bin_funct_annotations', 'bin_translated_genes', 'bin_untranslated_genes', 'prokka_out'],
		'BIN_CLASSIFICATION':['bin_taxonomy.tab','contig_taxonomy.tab']
		}

	if step == 0:
		res_file = res_file0
	elif step == 1:
		res_file = res_file1
	else:
		res_file = res_file2

	for sample in name_list:
		res_file['READ_QC'].append(sample)

		html_list = [os.path.join(tmp_dir, 'READ_QC', sample, 'pre-QC_report', i) for i in [f'{sample}_1_fastqc.html', f'{sample}_2_fastqc.html'] ]
		outdict['preQc'] = {}
		for html in html_list:
			sample_name = os.path.splitext(os.path.split(html)[1])[0]
			outdict['preQc'][sample_name] = get_fqc(html)

		html_list = [os.path.join(tmp_dir, 'READ_QC', sample, 'post-QC_report', i) for i in [f'final_pure_reads_1_fastqc.html', f'final_pure_reads_2_fastqc.html'] ]
		outdict['postQc'] = {}
		for html in html_list:
			sample_name = os.path.splitext(os.path.split(html)[1])[0]
			outdict['postQc'][sample_name] = get_fqc(html)

	get_file_path_names = ['kronagram.html','blobplot_figures','blobplot_figures_only_binned_contigs','reassembled_bins.png','reassembly_results.png']
	df_abu = pd.DataFrame()
	df_tax = pd.DataFrame()
	for d, fs in res_file.items():
		target_dir = os.path.join(outdir, d)
		mkdir(target_dir)
		tmp_list = [os.path.join(tmp_dir, d, i) for i in fs]
		try:
			copy_file(tmp_list, target_dir)
		except Exception as e:
			logging.error(f'copy_file {e}')

		for i in fs:
			tmp_file = os.path.join(target_dir, i)
			if os.path.isfile(tmp_file) or os.path.isdir(tmp_file):

				if i == 'QUAST_out':
					tmp_html = os.path.join(target_dir, i, 'report.html')
					outdict['ASSEMBLY_report'] = parse_quast_html(tmp_html)

				if i == 'reassembled_bins.stats':
					outdict['bin_sta'] = [v for i, v in pd.read_csv(tmp_file, sep='\t').to_dict(orient='index').items()]

				if i == 'bin_abundance_table.tab':
					df_abu = pd.read_csv(tmp_file,sep='\t',index_col=0)
				if i == 'bin_taxonomy.tab':
					df_tax = pd.read_csv(tmp_file, sep='\t', header=None)
					df_tax.index = df_tax[0].str.split('.').str[:2].str.join('.')

				if i in get_file_path_names:
					# if i == 'kronagram.html':
					# 	new_html_file = os.path.join(outdir, d,'kronagram_new.html')
					# 	tmp_file = parse_krona_html(tmp_file, new_html_file)
					outdict.update(get_file_path(tmp_file))
	try:
		if (not df_abu.empty) and (not df_tax.empty):
			df_merge = pd.merge(df_tax, df_abu, left_index=True, right_index=True, how='outer').reset_index()
			df_merge.columns = ['bin','bin_fa','tax','reads']
			outdict['bin_tax'] = [v for i, v in df_merge.to_dict(orient='index').items()]
	except Exception as e:
		logging.error(f'merge_tax_abu {e}')

	for sample in name_list:
		try:
			os.remove(os.path.join(outdir, 'READ_QC', sample, "host_reads_1.fastq"))
			os.remove(os.path.join(outdir, 'READ_QC', sample, "host_reads_2.fastq"))
		except Exception as e:
			logging.error(f'rm_file {e}')
	return outdict

if __name__ == '__main__':
	bin_dir = os.path.split(os.path.realpath(__file__))[0]
	memory_total = str(int(psutil.virtual_memory().total/1024/1024/1024))
	cpu_num = psutil.cpu_count()
	if cpu_num > 4:
		cpu_num = int(cpu_num * 0.9)
	parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-i', '--input', nargs='+', action='append', required=True, help='fastq file for metawrap')
	parser.add_argument('-n', '--name', nargs='+', required=True, help='sample name for fastq')
	parser.add_argument('-o', '--outdir', default='./', help='out dir for shell')
	parser.add_argument('-t', '--threads', default=cpu_num, help='threads for metawrap')
	parser.add_argument('-m', '--memory', default=memory_total, help='memory(G) for fastqc')
	parser.add_argument('-config', '--config', default=os.path.join(bin_dir, 'DB_config.json'), help='database path json')
	parser.add_argument('-step', '--step', default='mgTax', choices=['mgTax','mgAss', 'mgBin'], help='chose step for analysis')
	parser.add_argument('-tID', '--taskID', default='', help='task ID for report status')
	parser.add_argument('-debug', '--debug', action='store_true')

	args = parser.parse_args()
	step = 0
	step_ana = get_kraken_shell
	if args.step == 'mgAss':
		step = 1 
		step_ana = get_assemble_shell
	elif args.step == 'mgBin':
		step = 2
		step_ana = get_binning_shell

	outdir = args.outdir
	mkdir(outdir)
	logfile = os.path.join(outdir, 'log')
	logging.basicConfig(level=logging.INFO, filename=logfile, format='%(asctime)s %(levelname)s %(message)s',datefmt='%Y-%m-%d %H:%M:%S')
	# logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s %(message)s',datefmt='%Y-%m-%d %H:%M:%S')
	
	taskID = args.taskID
	post_pid(taskID)
	status_report = os.path.join(outdir, 'status_report.txt')

	db_config = args.config
	data_list = args.input
	name_list = args.name 
	outdir = args.outdir
	threads = args.threads
	memory = args.memory
	tmp_dir = os.path.join(outdir, 'tmp_dir')
	mkdir(tmp_dir)

	try:
		s = f'{args.step}\tR\t'

		shell_str = step_ana(db_config, data_list, name_list, tmp_dir, threads, memory, kraken=True)
		sh_file = os.path.join(outdir, 'metawrap.sh')
		with open(sh_file, 'w') as h:
			print(shell_str, file=h)

		with open(db_config) as h:
			db_path = json.load(h)
		# db_path = '/mnt/data/metawrap_db/'
		data_list_tmp = [i for data in data_list for i in data]
		path_list = list(db_path.values()) + [outdir] + data_list_tmp
		logging.info(path_list)
		# sys.exit()
		run_sh(sh_file, path_list)

		outdict = copy_file_batch(name_list, tmp_dir, outdir, step=step)
		with open(os.path.join(outdir, f'{args.step}.json'), 'w') as H:
			json.dump(outdict, H, indent=2)

	except Exception as e:
		logging.error(e)
		s = f'{args.step}\tE\t'

	try:
		write_status(status_report, args.step)
		post_url(taskID, args.step)
	except Exception as e:
		logging.error(f'{args.step} status {e}')
	if not args.debug:
		try:
			shutil.rmtree(tmp_dir)
		except Exception as e:
			logging.error(e)

'''
	res_file = {
		'READ_QC':[],
		'KRAKEN':['kronagram.html'],
		'ASSEMBLY':['final_assembly.fasta','QUAST_out'],
		'BLOBOLOGY':['blobplot_figures','blobplot_figures_only_binned_contigs','final_assembly.binned.blobplot','final_assembly.blobplot'],
		'QUANT_BINS':['bin_abundance_table.tab','quant_files'],
		'BIN_REASSEMBLY':['original_bins.stats','reassembled_bins','reassembled_bins.png','reassembled_bins.stats','reassembly_results.eps','reassembly_results.png'],
		'FUNCT_ANNOT':['bin_funct_annotations', 'bin_translated_genes', 'bin_untranslated_genes', 'prokka_out'],
		'BIN_CLASSIFICATION':['bin_taxonomy.tab','contig_taxonomy.tab']
		}

	for sample in name_list:
		res_file['READ_QC'].append(os.path.join(tmp_dir, 'READ_QC', sample))

	for d, fs in res_file.items():
		target_dir = os.path.join(outdir, d)
		mkdir(target_dir)
		tmp_list = [os.path.join(tmp_dir, d, i) for i in fs]
		try:
			copy_file(tmp_list, target_dir)
		except Exception as e:
			logging.error(f'copy_file {e}')

	for sample in name_list:
		try:
			os.remove(os.path.join(outdir, 'READ_QC', sample, "host_reads_1.fastq"))
			os.remove(os.path.join(outdir, 'READ_QC', sample, "host_reads_2.fastq"))
		except Exception as e:
			logging.error(f'rm_file {e}')

	try:
		shutil.rmtree(tmp_dir)
	except Exception as e:
		logging.error(f'rm_file {e}')

'''

