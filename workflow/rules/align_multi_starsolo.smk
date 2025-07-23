#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" align_multi_starsolo

Performs multimapping alignment using STARsolo, allowing for a high number of
multimappers for downstream reassignment using stellarscope

"""

def get_input_files(wildcards, read_type):
    _sruns = samples.loc[wildcards.sampid]['runid']
    _layouts = runs.loc[_sruns]['layout']

    if not all(_ == 'PAIRED' for _ in _layouts):
        msg = f'All runs for sample "{wildcards.sampid}" should be PAIRED for this rule'
        raise WorkflowError(msg)

    _sruns = sorted(_sruns)  # Sort run IDs to ensure order

    if read_type == 'R1':
        return expand("results/fastq/{{sampid}}/{sra_run}_1.fastq.gz", sra_run=_sruns)
    elif read_type == 'R2':
        return expand("results/fastq/{{sampid}}/{sra_run}_2.fastq.gz", sra_run=_sruns)
    else:
        raise ValueError(f'Invalid read_type: {read_type}. Expected "R1" or "R2".')


def align_multi_starsolo_args(wildcards, input):
    argd = {}
    _chem = meta_table.loc[wildcards.sampid]['chemistry']
    # readFilesIn
    if config['barcode_mate'][_chem] == 'R1':
        readBC = ','.join(input['R1'])
        readTX = ','.join(input['R2'])
    elif config['barcode_mate'][_chem] == 'R2':
        readBC = ','.join(input['R2'])
        readTX = ','.join(input['R1'])
    argd['readFilesIn'] = f'{readTX}  {readBC}'
    
    # readFilesCommand
    tf = input['R1'][0] if isinstance(input['R1'], list) else input['R1']
    assert isinstance(tf, str), f'test file {tf} is type {type(tf)}'
    if tf.endswith('.gz'):
        argd['readFilesCommand'] = '"gunzip -c"'

    # default args
    for k,v in config['align_multi_starsolo']['_default'].items():
        argd[k] = v
    
    # chemistry-specific args
    for k,v in config['align_multi_starsolo']['chemistry'][_chem].items():
        argd[k] = v
    
    ret = ' '.join(f'--{k} {v}' for k,v in argd.items())
    return ret


rule align_multi_starsolo:
    output:
        bam = "results/align_multi_starsolo/{sampid}/Aligned.sortedByCoord.out.bam",
        gene_raw_mtx = "results/align_multi_starsolo/{sampid}/Solo.out/Gene/raw/matrix.mtx",
        gene_raw_bcode = "results/align_multi_starsolo/{sampid}/Solo.out/Gene/raw/barcodes.tsv",
        gene_raw_feats = "results/align_multi_starsolo/{sampid}/Solo.out/Gene/raw/features.tsv",
        gene_filt_mtx = "results/align_multi_starsolo/{sampid}/Solo.out/Gene/filtered/matrix.mtx",
        gene_filt_bcode = "results/align_multi_starsolo/{sampid}/Solo.out/Gene/filtered/barcodes.tsv",
        gene_filt_feats = "results/align_multi_starsolo/{sampid}/Solo.out/Gene/filtered/features.tsv",
        sj_raw_mtx = "results/align_multi_starsolo/{sampid}/Solo.out/SJ/raw/matrix.mtx",
        sj_raw_bcode = "results/align_multi_starsolo/{sampid}/Solo.out/SJ/raw/barcodes.tsv",
        sj_raw_feats = "results/align_multi_starsolo/{sampid}/Solo.out/SJ/raw/features.tsv",
        velo_raw_spl = "results/align_multi_starsolo/{sampid}/Solo.out/Velocyto/raw/spliced.mtx",
        velo_raw_unspl = "results/align_multi_starsolo/{sampid}/Solo.out/Velocyto/raw/unspliced.mtx",
        velo_raw_ambig = "results/align_multi_starsolo/{sampid}/Solo.out/Velocyto/raw/ambiguous.mtx",
        velo_raw_bcode = "results/align_multi_starsolo/{sampid}/Solo.out/Velocyto/raw/barcodes.tsv",
        velo_raw_feats = "results/align_multi_starsolo/{sampid}/Solo.out/Velocyto/raw/features.tsv",
        velo_filt_spl = "results/align_multi_starsolo/{sampid}/Solo.out/Velocyto/filtered/spliced.mtx",
        velo_filt_unspl = "results/align_multi_starsolo/{sampid}/Solo.out/Velocyto/filtered/unspliced.mtx",
        velo_filt_ambig = "results/align_multi_starsolo/{sampid}/Solo.out/Velocyto/filtered/ambiguous.mtx",
        velo_filt_bcode = "results/align_multi_starsolo/{sampid}/Solo.out/Velocyto/filtered/barcodes.tsv",
        velo_filt_feats = "results/align_multi_starsolo/{sampid}/Solo.out/Velocyto/filtered/features.tsv",
    input:
        R1=lambda wildcards: get_input_files(wildcards, 'R1'),
        R2=lambda wildcards: get_input_files(wildcards, 'R2')
    params:
        cli_args = align_multi_starsolo_args
    threads: snakemake.utils.available_cpu_count() 
    conda: "../envs/star.yaml"
    shell:
        '''
        tdir=$(mktemp -d {config[tmpdir]}/{rule}.XXXXXX)

        ulimit -n 4096

        STAR\
            --runThreadN {threads}\
            --genomeDir {config[star_index]}\
            {params.cli_args}\
            --outFileNamePrefix $tdir/

        mkdir -p $(dirname {output[0]})
        rsync -a --exclude '_STARtmp' $tdir/ $(dirname {output[0]})
        
        # make absolute symlink to relative
        if [[ -L {output.sj_raw_feats} ]]; then
            rm -f {output.sj_raw_feats}
            ln -s ../../../SJ.out.tab {output.sj_raw_feats}
        fi
        '''

#rule index_bam:
#    output:
#        '{prefix}.sortedByCoord.out.bam.bai'    
#    input:
#        '{prefix}.sortedByCoord.out.bam'
#    conda: "../envs/stellarscope_1.4.yaml"
#    shell:
#        '''
#        samtools index {input[0]}
#        '''

