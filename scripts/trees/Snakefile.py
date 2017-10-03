include: "Snakefile_utils.py"

##### Load config
MY_HOME=                      config["home"]
RESOURCES=                    config["resources"]
SEEDFILE=                     config['seedfile']  # seedfile corresponds to 1 directory per line
SCRATCH=                      config['scratch']  # LOCAL_SATA or LOCAL_SSD for /local10G, LOCAL_SCRATCH for sherlock
PARTITION=                    config['partition'] # partition to run jobs on

##### Paths

workdir:                      MY_HOME+'/log'

ANACONDA=		      RESOURCES+'/anaconda2'
MUSCLE=                       RESOURCES+'/muscle3.8.31/muscle'
FASTTREE=                     RESOURCES+'/FastTree/FastTreeMP'

LOOKUP_SEQ_STR_UIDS=          MY_HOME+'/scripts/trees/lookup_seq_str_uids.py'
ALIGN_UNGAPPED=               MY_HOME+'/scripts/trees/align_ungapped.py'
LOOKUP_GERMLINE=              MY_HOME+'/scripts/trees/lookup_germline.py'
REPAIR=                       MY_HOME+'/scripts/trees/repair_align.py'

##### Parameters
REFS=                         MY_HOME+'/resources'
V_GERMLINE_SEQS=              REFS+'/Vsegments_20150201.fasta.cleaned'
J_GERMLINE_SEQS=              REFS+'/Jsegments_20150201.fasta.cleaned'

DB_SEQS=                      MY_HOME+'/trees/Bcell_flu_high_res.sequences.csv'
DB_SEQ_STRS=                  MY_HOME+'/trees/Bcell_flu_high_res.sequence_strings.csv'

##### Load samples
SEEDS = []
with open(SEEDFILE) as f:
    for line in f:
        SEEDS.append(line.strip())

##### Rules

rule all:
  input: expand("{dir}/done", dir=SEEDS)
  params: name='all', partition='general', mem='1024'

rule lookup_seqids:
  """ Get sequence string uids for each lineage """
  output: '{dir}/sequence_string_uids.txt',
          '{dir}/sequence_string_uids_V6_Full.txt'
  params: name='lookup_seqids', partition=PARTITION, mem="5300"
  run:
      lineage_uid = os.path.dirname(str(output[0])).split("/")[-1]
      shell("source {ANACONDA}/bin/activate {ANACONDA} && "
            "python {LOOKUP_SEQ_STR_UIDS} "
            "{DB_SEQS} {lineage_uid} > {output[0]} && "
            'grep -e "^.99" < {output[0]} > {output[1]}')

rule align_ungapped:
  input: rules.lookup_seqids.output
  output: '{dir}/alignment_ungapped.fasta'
  params: name='align_ungapped', partition=PARTITION, mem="5300"
  run:
      lineage_uid = os.path.dirname(str(output[0])).split("/")[-1]
      shell("source {ANACONDA}/bin/activate {ANACONDA} && "
            "python {ALIGN_UNGAPPED} "
            "{DB_SEQ_STRS} {input[1]} > {output[0]} ")

rule lookup_germline:
  """ Get sequence string uids for each lineage """
  output: '{dir}/germline.fasta'
  params: name='lookup_germline', partition=PARTITION, mem="5300"
  run:
      lineage_uid = os.path.dirname(str(output[0])).split("/")[-1]
      shell("source {ANACONDA}/bin/activate {ANACONDA} && "
            "python {LOOKUP_GERMLINE} "
            "{DB_SEQS} {lineage_uid} "
            "{V_GERMLINE_SEQS} {J_GERMLINE_SEQS} "
            "> {output[0]} ")

rule muscle_refine:
  input: rules.align_ungapped.output
  output: '{dir}/alignment_ungapped_refined.fasta'
  params: name='muscle_refine', partition=PARTITION, mem="16000"
  run:
    shell('{MUSCLE} '
          '-seqtype dna '
          '-refine '
          '-maxiters 1 '
          '-diags '
          '-gapopen -5000.0 '
          '-in {input[0]} '
          '-out {output[0]} ')

rule muscle_pf:
  input: rules.muscle_refine.output,
         rules.lookup_germline.output
  output: '{dir}/alignment_ungapped_refined_germline.fasta'
  params: name='muscle_pf', partition=PARTITION, mem="16000"
  run:
    shell('{MUSCLE} '
          '-seqtype dna '
          '-profile '
          '-maxiters 1 '
          '-diags '
          '-in1 {input[0]} '
          '-in2 {input[1]} '             
          '-out {output[0]} ')

rule fasttree:
  input: rules.muscle_pf.output
  output: '{dir}/fasttree.nwk'
  params: name='fasttree', partition=PARTITION, mem='5300'
  threads: 3
  run:
    shell('{FASTTREE} '
          '-nt -gtr '
          '< {input[0]} '
          '> {output[0]} ')

rule repair:
  input: rules.muscle_pf.output,
         rules.fasttree.output
  output: '{dir}/alignment_ungapped_refined_germline.rep.fasta',
          '{dir}/fasttree.rep.nwk'
  params: name='repair', partition=PARTITION, mem='5300'
  threads: 3
  run:
    shell("source {ANACONDA}/bin/activate {ANACONDA} && "
          "python {REPAIR} "
          "{input[0]} {input[1]} "
          "{output[0]} {output[1]} "
          "0.5 ")
    
rule clean:
  input: rules.repair.output  
  output: "{dir}/done"
  params: name='clean', partition=PARTITION, mem='1024'
  shell: 'touch {output[0]}'

