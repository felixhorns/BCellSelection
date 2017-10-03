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

NAME_INTERNAL_NODES=               MY_HOME+'/scripts/fitnessMutations/name_internal_nodes.py'
RECONSTRUCT_ANCESTRAL_SEQUENCES=   MY_HOME+'/scripts/fitnessMutations/reconstruct_ancestral_sequences.py'
INFER_FITNESS=                     MY_HOME+'/scripts/fitnessMutations/infer_fitness.py'
ANNOTATE_MUTATIONS=                MY_HOME+'/scripts/fitnessMutations/annotate_mutations.py'
CLEAN_ALN=                         MY_HOME+'/scripts/fitnessMutations/clean_aln.py'

##### Parameters
REFS=                         MY_HOME+'/resources'

##### Load samples
SEEDS = []
with open(SEEDFILE) as f:
    for line in f:
        SEEDS.append(line.strip())

##### Rules

rule all:
  input: expand("{dir}/fitnessMutations.clean.done", dir=SEEDS)
  params: name='all', partition='normal', mem='1024'

rule name_internal_nodes:
  """ Name internal nodes of tree """
  input: '{dir}/fasttree.rep.nwk'
  output: '{dir}/fasttree.rep.named.nwk'
  params: name='name_internal_nodes', partition='normal', mem='5300'
  run:
      shell("source {ANACONDA}/bin/activate {ANACONDA} && "
            "python {NAME_INTERNAL_NODES} {input[0]} {output[0]} ")

rule reconstruct_ancestral_sequences:
  """ Reconstruct ancestral sequences """
  input: rules.name_internal_nodes.output,
         '{dir}/alignment_ungapped_refined_germline.rep.fasta'
  output: '{dir}/alignment_ungapped_refined_germline.rep.asr.fasta'
  params: name='reconstruct_ancestral_sequences', partition='normal', mem='5300'
  run:
      shell("source {ANACONDA}/bin/activate {ANACONDA} && "
            "python {RECONSTRUCT_ANCESTRAL_SEQUENCES} "
            "{input[0]} {input[1]} "
            "{output[0]} ")

rule clean_aln:
  """ Clean alignment by removing positions where all observed sequences have a gap """
  input: rules.reconstruct_ancestral_sequences.output
  output: '{dir}/alignment_ungapped_refined_germline.rep.asr.clean.fasta'
  params: name='clean_aln', partition='normal', mem='5300'
  run:
      shell("source {ANACONDA}/bin/activate {ANACONDA} && "
            "python {CLEAN_ALN} "
            "{input[0]} "
            "{output[0]} ")
  
      
rule infer_fitness:
  """ Infer fitness based on tree shape """
  input: rules.name_internal_nodes.output[0],
         rules.reconstruct_ancestral_sequences.output[0]
  output: '{dir}/fitness.nwk',
          '{dir}/df_fitness.csv',
          '{dir}/fitness.pdf',
          '{dir}/fitness.labeled.pdf' 
  params: name='infer_fitness', partition='normal', mem='5300'
  run:
      shell("source {ANACONDA}/bin/activate {ANACONDA} && "
            "python {INFER_FITNESS} "
            "{input[0]} {input[1]} "
            "{output[0]} {output[1]} {output[2]} {output[3]} ")

rule annotate_mutations:
  """ Get and annotate mutations on each branch """
  input: rules.clean_aln.output[0],
         rules.infer_fitness.output[0]
  output: '{dir}/df_mutations.clean.csv'
  params: name='annotate_mutations', partition='unrestricted', mem='60000'
  run:
      shell("source {ANACONDA}/bin/activate {ANACONDA} && "
            "python {ANNOTATE_MUTATIONS} "
            "{input[0]} {input[1]} "
            "{output[0]} ")

rule clean:
  input: rules.annotate_mutations.output  
  output: "{dir}/fitnessMutations.clean.done"
  params: name='clean', partition='normal', mem='1024'
  shell: 'touch {output[0]}'
