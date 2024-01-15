#!/ldfssz1/ST_OCEAN/USER/xudayou/program/miniconda3/bin/python3


import os,sys
from subprocess import getoutput
from glob import glob

if len(sys.argv)!=3:
    print('usage: '+sys.argv[0]+' input_sra output_folder')
    print('sra2fq-->SOAPnuk-->megahit_assemble-->metabat2_binning-->checkm')
    print('resource usage: 50G mem, 10CPU')
    exit(1)

#initialize
inp=sys.argv[1]
out=sys.argv[2]
# name=sys.argv[1].split('/')[-1]
# name=name.split('.sra')[0]
name = out.split('/')[-1]
#os.environ['PATH']='/ldfssz1/ST_OCEAN/USER/xudayou/python/bin:'+os.environ['PATH']
current_path = os.path.dirname(os.path.abspath(__file__))

##sra2fq
#if not os.path.exists(inp):
#    print('Error:input '+inp+' not found')
#    exit(1)

#print('---START AT---'+getoutput('date'))
#if int(os.path.getsize(inp)) > 25000000000:
#    print('Error: input is too large')
#    exit(1)

#os.system('mkdir -p '+sys.argv[2])

#if 'bad version while' in getoutput(current_path+'/bin/sratoolkit-2.10.8/bin/fastq-dump --split-e -F --gzip '+inp+' -O '+out):
#    print('There is problem in extract sra, exit.')
#    exit(1)

#try:
#    fq=glob(out+'/'+name+'*.fastq.gz')
#except:
#    raise Exception('cannot find fastq.gz')
#    exit(1)

##SOAPnuke
#os.system('mkdir -p '+out+'/filter/')
##if SE
#if len(fq)==1:
#    while True:
#        if os.system(current_path+'/bin/SOAPnuke filter -c 0 -l 15 -q 0.2 -n 0.05 -Q 2 -d -o '+out+'/filter/'+' -1 '+fq[0]+' -C '+out+'/filter/'+name+'.clean.fq.gz')==0:
#            os.system('perl '+current_path+'/bin/soapnuke_stat.pl '+out+'/filter/Basic_Statistics_of_Sequencing_Quality.txt '+out+'/filter/Statistics_of_Filtered_Reads.txt >  '+out+'/filter/'+name+'.filter.stat.xls')
#        print('SOAPNUKE_DONE')
#        break
##if PE
#else:
#    for i in fq:
#        if '_1.fastq.gz' in i:
#            fq1=i
#        elif '_2.fastq.gz' in i:
#            fq2=i
#    while True:
#        if os.system(current_path+'/bin/SOAPnuke filter -c 0 -l 15 -q 0.2 -n 0.05 -Q 2 -d -o '+out+'/filter/'+' -1 '+fq1+' -2 '+fq2+' -C '+out+'/filter/'+name+'_1.clean.fq.gz -D '+out+'/filter/'+name+'_2.clean.fq.gz')==0:
#            os.system('perl '+current_path+'/bin/soapnuke_stat.pl '+out+'/filter/Basic_Statistics_of_Sequencing_Quality.txt '+out+'/filter/Statistics_of_Filtered_Reads.txt > '+out+'/filter/'+name+'.filter.stat.xls')
#        print('SOAPNUKE_DONE')
#        break

#megahit_assemble & metabat2_binning
try:
    clean_fq=glob(out+'/filter/'+name+'*.clean.fq.gz')
except:
    raise Exception('cannot find clean.fq.gz')
    exit(1)

#if SE
if len(clean_fq)==1:
    while True:
        if os.system(current_path+'/bin/MEGAHIT-1.2.9/bin/megahit -r '+clean_fq[0]+' -o '+out+'/megahit --out-prefix '+name+' --min-count 2 --k-min 33 --k-max 63 --k-step 10 --no-mercy --memory 0.4 -t 15')!=0:
            name=out.split('/')[-1]
            os.system('cp {0}/megahit/{1}.log {0}/'.format(out, name))
            os.system('rm -rf {0}/megahit'.format(out))
            continue
        print('MEGAHIT_DONE')
        os.system('mkdir '+out+'/metabat')
        os.system('gunzip -dc '+clean_fq[0]+' > '+out+'/metabat/'+name+'.clean.fastq')
        if os.stat(out+'/megahit/'+name+'.contigs.fa').st_size == 0:
            print('assembly size is 0, exit.')
            os.system('rm -r '+sys.argv[2])
            os.system('mkdir -p '+sys.argv[2])
            exit(1)
        if os.system(current_path+'/bin/metawrap/metawrap binning --single-end -l 1000 -o '+out+'/metabat -t 15 -a '+out+'/megahit/'+name+'.contigs.fa --metabat2 '+out+'/metabat/'+name+'.clean.fastq') == 0:
            print('MEGABAT_DONE')
            break
#if PE
#else:
#    for j in clean_fq:
#        if '_1.clean.fq.gz' in j:
#            c_fq1=j
#        elif '_2.clean.fq.gz' in j:
#            c_fq2=j
#    while True:
#        if os.system(current_path+'/bin/MEGAHIT-1.2.9/bin/megahit -1 '+c_fq1+'  -2 '+c_fq2+' -o '+out+'/megahit --out-prefix '+name+' --min-count 2 --k-min 33 --k-max 63 --k-step 10 --no-mercy --memory 0.4 -t 15') != 0:
#            name=out.split('/')[-1]
#            os.system('cp {0}/megahit/{1}.log {0}/'.format(out, name))
#            os.system('rm -rf {0}/megahit'.format(out))
#            continue
#        print('MEGAHIT_DONE')
#        os.system('mkdir '+out+'/metabat')
#        os.system('gunzip -dc '+c_fq1+' > '+out+'/metabat/'+name+'_clean_1.fastq')
#        os.system('gunzip -dc '+c_fq2+' > '+out+'/metabat/'+name+'_clean_2.fastq')
#        if os.stat(out+'/megahit/'+name+'.contigs.fa').st_size == 0:
#            print('assembly size is 0, exit.')
#            os.system('rm -r '+sys.argv[2])
#            os.system('mkdir -p '+sys.argv[2])
#            exit(1)
#        if os.system(current_path+'/bin/metawrap/metawrap binning -l 1000 -o '+out+'/metabat -t 10 -a '+out+'/megahit/'+name+'.contigs.fa --metabat2 '+out+'/metabat/'+name+'_clean_1.fastq'+' '+out+'/metabat/'+name+'_clean_2.fastq') == 0:
#            print('MEGABAT_DONE')
#            break
#checkm
if 'OSError: [Errno 2] No such file or directory:' in getoutput(current_path+'/bin/metawrap/checkm lineage_wf -f '+out+'/metabat/CheckM.txt -t 10 -x fa '+out+'/metabat/metabat2_bins '+out+'/metabat/checkm '):
    os.system('rm -r '+out+'/metabat/checkm')
    os.system(current_path+'/bin/metawrap/checkm lineage_wf -f '+out+'/metabat/CheckM.txt -t 10 -x fa '+out+'/metabat/metabat2_bins '+out+'/metabat/checkm ')

##save_results
#os.system('mkdir '+out+'/good_bin')
#os.system('mkdir '+out+'/bad_bin')
#os.system('mv '+out+'/megahit/'+name+'.contigs.fa '+out)
#if os.path.exists(out+"/metabat/CheckM.txt"):
#    os.system("for i in `grep 'bin' "+out+"/metabat/CheckM.txt | awk '{if($13>50&&$14<10)print$1}' `;do mv "+out+"/metabat/metabat2_bins/$i.fa "+out+"/good_bin;done")
#    os.system('mv '+out+'/metabat/metabat2_bins/* '+out+'/bad_bin ')
#    os.system('mv '+out+'/metabat/CheckM.txt '+out)

##remove tmp
#os.system('rm -rf '+out+'/megahit '+out+'/metabat')

print(inp+' analyze_done')
print('---END AT---'+getoutput('date'))
