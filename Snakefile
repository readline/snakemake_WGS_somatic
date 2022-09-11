from os.path import join
import os
import pandas as pd
from scripts.Load import somsamplesheet

snakedir = os.getcwd()
print(snakedir)
configfile: 'config.yaml'
print(config)
somdic, ssdic = somsamplesheet(config['samplesheet'])
workdir: config['workdir']
#sbcmd="sbatch --cpus-per-task={threads} --mem={cluster.mem}"
#sbcmd+=" --time={cluster.time} --partition={cluster.partition}"
#sbcmd+=" --out={cluster.out} {cluster.extra}"
print(somdic)
itv4 = ['%.4d'%(itv) for itv in range(1, 50+1)]

rule all:
    input:
        lambda wildcards: ["12.SV.manta/{}/results/variants/candidateSmallIndels.vcf.gz".format(som) for som in somdic],
        lambda wildcards: ["11.Somatic.Strelka/{}/results/variants/somatic.snvs.vcf.gz".format(som) for som in somdic],
        lambda wildcards: ["11.Somatic.Mutect/{}/{}.mutect2.flt.vcf".format(som,som) for som in somdic],
        lambda wildcards: ["11.Somatic.Lofreq/{}/{}.lofreqsomatic_final.snvs.vcf.gz".format(som,som) for som in somdic],
        lambda wildcards: ["11.Somatic.MuSE/{}/{}.MuSE.vcf.gz".format(som,som) for som in somdic],
        lambda wildcards: ["11.Somatic.Varscan2/{}/{}.snp.LOH.vcf.gz".format(som,som) for som in somdic],
        
rule manta_som:
    input:
        bam1=lambda wildcards: "02.Alignment/Level3/{}/{}.BQSR.bam".format(somdic[wildcards.som][0], somdic[wildcards.som][0]),
        bam0=lambda wildcards: "02.Alignment/Level3/{}/{}.BQSR.bam".format(somdic[wildcards.som][1], somdic[wildcards.som][1]),
    output:
        folder=directory("12.SV.manta/{som}"),
        result="12.SV.manta/{som}/results/variants/candidateSmallIndels.vcf.gz",
    log:
        out = snakedir+"/logs/D01.som_manta/{som}.o",
        err = snakedir+"/logs/D01.som_manta/{som}.e",
    threads:  16
    resources:
        mem  = '32g',
        extra = ' --gres=lscratch:10 ',
    shell:
        '''module load {config[modules][manta]}
        configManta.py \
          --normalBam {input.bam0} \
          --tumorBam {input.bam1} \
          --referenceFasta {config[references][fasta]} \
          --runDir {output.folder} >{log.out} 2>{log.err}
        {output.folder}/runWorkflow.py -m local -j {threads}  >>{log.out} 2>>{log.err} '''

rule strelka_som:
    input:
        bam1=lambda wildcards: "02.Alignment/Level3/{}/{}.BQSR.bam".format(somdic[wildcards.som][0], somdic[wildcards.som][0]),
        bam0=lambda wildcards: "02.Alignment/Level3/{}/{}.BQSR.bam".format(somdic[wildcards.som][1], somdic[wildcards.som][1]),
        manta="12.SV.manta/{som}/results/variants/candidateSmallIndels.vcf.gz",
    output:
        folder=directory("11.Somatic.Strelka/{som}"),
        result="11.Somatic.Strelka/{som}/results/variants/somatic.snvs.vcf.gz",
    log:
        out = snakedir+"/logs/D02.som_strelka/{som}.o",
        err = snakedir+"/logs/D02.som_strelka/{som}.e",
    threads:  16
    resources:
        mem  = '32g',
        extra = ' --gres=lscratch:10 ',
    shell:
        '''module load {config[modules][strelka]}
        configureStrelkaSomaticWorkflow.py \
          --tumorBam {input.bam1} \
          --normalBam {input.bam0} \
          --referenceFasta {config[references][fasta]} \
          --runDir {output.folder} \
          --indelCandidates {input.manta} \
          --outputCallableRegions >{log.out} 2>{log.err}
        {output.folder}/runWorkflow.py -m local -j {threads}  >>{log.out} 2>>{log.err}'''
        
        
rule mutect_som:
    input:
        bam1=lambda wildcards: "02.Alignment/Level3/{}/{}.BQSR.bam".format(somdic[wildcards.som][0], somdic[wildcards.som][0]),
        bam0=lambda wildcards: "02.Alignment/Level3/{}/{}.BQSR.bam".format(somdic[wildcards.som][1], somdic[wildcards.som][1]),
    params:
        normal=lambda wildcards: somdic[wildcards.som][1],
        bed=config['references']['wgsscatter']+'/temp_{itv}_of_50.bed',
    output:
        vcf1="11.Somatic.Mutect/{som}/{som}.{itv}.mutect2.vcf",
        f1r2=temp("11.Somatic.Mutect/{som}/{som}.{itv}.f1r2.tar.gz"),
        model=temp("11.Somatic.Mutect/{som}/{som}.{itv}.read-orientation-model.tar.gz"),
        pile=temp("11.Somatic.Mutect/{som}/{som}.{itv}.getpileupsummaries.table"),
        seg=temp("11.Somatic.Mutect/{som}/{som}.{itv}.segments.table"),
        cont=temp("11.Somatic.Mutect/{som}/{som}.{itv}.calculatecontamination.table"),
        vcf2="11.Somatic.Mutect/{som}/{som}.{itv}.mutect2.flt.vcf",
    log:
        out = snakedir+"/logs/D03.som_mutect/{som}.{itv}.o",
        err = snakedir+"/logs/D03.som_mutect/{som}.{itv}.e",
    threads:  8
    resources:
        mem  = '16g',
        extra = ' --gres=lscratch:10 ',
    shell:
        '''module load {config[modules][gatk]}
        gatk --java-options "-Xmx24g -Xms24g -Djava.io.tmpdir=/lscratch/$SLURM_JOB_ID" \
          Mutect2 \
          -R {config[references][fasta]} \
          -L {params.bed} \
          -I {input.bam1} \
          -I {input.bam0} \
          -normal {params.normal} \
          --panel-of-normals {config[references][mutectpon]}\
          -O {output.vcf1} \
          -germline-resource {config[references][afonlygnomad]} \
          --f1r2-tar-gz {output.f1r2} >{log.out} 2>{log.err}
        gatk --java-options "-Xmx24g -Xms24g -Djava.io.tmpdir=/lscratch/$SLURM_JOB_ID" \
          LearnReadOrientationModel \
          -O {output.model} \
          -I {output.f1r2} >>{log.out} 2>>{log.err}
        gatk --java-options "-Xmx24g -Xms24g -Djava.io.tmpdir=/lscratch/$SLURM_JOB_ID" \
          GetPileupSummaries \
          -I {input.bam1} \
          -V {config[references][exaccommon]} \
          -L {params.bed} \
          -O {output.pile} >>{log.out} 2>>{log.err}
        gatk --java-options "-Xmx24g -Xms24g -Djava.io.tmpdir=/lscratch/$SLURM_JOB_ID" \
          CalculateContamination \
          -I {output.pile} \
          -tumor-segmentation {output.seg} \
          -O {output.cont} >>{log.out} 2>>{log.err}
        gatk --java-options "-Xmx24g -Xms24g -Djava.io.tmpdir=/lscratch/$SLURM_JOB_ID" \
          FilterMutectCalls \
          --reference {config[references][fasta]} \
          -V {output.vcf1} \
          --tumor-segmentation {output.seg} \
          --ob-priors {output.model} \
          -O {output.vcf2} >>{log.out} 2>>{log.err} '''

rule merge_mutect:
    input:
        vcf1=expand("11.Somatic.Mutect/{{som}}/{{som}}.{itv}.mutect2.vcf", itv=itv4),
        vcf2=expand("11.Somatic.Mutect/{{som}}/{{som}}.{itv}.mutect2.flt.vcf", itv=itv4),
    output:
        vcf1="11.Somatic.Mutect/{som}/{som}.mutect2.vcf",
        vcf2="11.Somatic.Mutect/{som}/{som}.mutect2.flt.vcf",
    log:
        out = snakedir+"/logs/D03b.merge_mutect/{som}.o",
        err = snakedir+"/logs/D03b.merge_mutect/{som}.e",
    threads:  2
    resources:
        mem  = '16g',
        extra = ' --gres=lscratch:10 ',
    run:
        inputvcfs1 = ' '.join(" {}".format(i) for i in input.vcf1)
        inputvcfs2 = ' '.join(" {}".format(i) for i in input.vcf2)
        shell('''
        ls {inputvcfs1} > {output.vcf1}
        ls {inputvcfs2} > {output.vcf2}
        ''')
        
rule lofreq_som:
    input:
        bam1=lambda wildcards: "02.Alignment/Level3/{}/{}.BQSR.bam".format(somdic[wildcards.som][0], somdic[wildcards.som][0]),
        bam0=lambda wildcards: "02.Alignment/Level3/{}/{}.BQSR.bam".format(somdic[wildcards.som][1], somdic[wildcards.som][1]),
    params:
        bam1='/lscratch/$SLURM_JOB_ID/{som}.{itv}.1.bam',
        bam0='/lscratch/$SLURM_JOB_ID/{som}.{itv}.0.bam',
        bed=config['references']['wgsscatter']+'/temp_{itv}_of_50.bed',
        prefix="11.Somatic.Lofreq/{som}/itv/{itv}.lofreq",
    output:
        vcf1="11.Somatic.Lofreq/{som}/itv/{itv}.lofreqsomatic_final.snvs.vcf.gz",
    log:
        out = snakedir+"/logs/D04.som_lofreq/{som}.{itv}.o",
        err = snakedir+"/logs/D04.som_lofreq/{som}.{itv}.e",
    threads:  16
    resources:
        mem  = '64g',
        extra = ' --gres=lscratch:200 ',
    shell:
        '''module load {config[modules][lofreq]} {config[modules][sambamba]}
        sambamba view \
          -f bam \
          -L {params.bed} \
          -t 4 \
          {input.bam0} |\
        lofreq indelqual \
          -f {config[references][fasta]} \
          --dindel \
          -o {params.bam0} - >> {log.out} 2>> {log.err}
        sambamba view \
          -f bam \
          -L {params.bed} \
          -t 4 \
          {input.bam1} |\
        lofreq indelqual \
          -f {config[references][fasta]} \
          --dindel \
          -o {params.bam1} - >> {log.out} 2>> {log.err}
        sambamba index -t {threads} {params.bam0} >> {log.out} 2>> {log.err}
        sambamba index -t {threads} {params.bam1} >> {log.out} 2>> {log.err}
        lofreq somatic \
          -n {params.bam0} \
          -t {params.bam1} \
          -o {params.prefix} \
          -l {params.bed} \
          -f {config[references][fasta]} \
          --threads {threads} \
          --call-indels \
          -d {config[references][snp138]} >>{log.out} 2>>{log.err} '''

rule merge_lofreq:
    input:
        vcf1=expand("11.Somatic.Lofreq/{{som}}/itv/{itv}.lofreqsomatic_final.snvs.vcf.gz", itv=itv4),
    output:
        vcf1="11.Somatic.Lofreq/{som}/{som}.lofreqsomatic_final.snvs.vcf.gz",
    log:
        out = snakedir+"/logs/D04.merge_lofreq/{som}.o",
        err = snakedir+"/logs/D04.merge_lofreq/{som}.e",
    threads:  2
    resources:
        mem  = '16g',
        extra = ' --gres=lscratch:10 ',
    run:
        inputvcfs1 = ' '.join(" {}".format(i) for i in input.vcf1)
        shell('''
        ls {inputvcfs1} > {output.vcf1}
        ''')
        
rule muse_som:
    input:
        bam1=lambda wildcards: "02.Alignment/Level3/{}/{}.BQSR.bam".format(somdic[wildcards.som][0], somdic[wildcards.som][0]),
        bam0=lambda wildcards: "02.Alignment/Level3/{}/{}.BQSR.bam".format(somdic[wildcards.som][1], somdic[wildcards.som][1]),
    params:
        bed=config['references']['wgsscatter']+'/temp_{itv}_of_50.bed',
        prefix="11.Somatic.MuSE/{som}/itv/{itv}",
    output:
        vcf="11.Somatic.MuSE/{som}/itv/{itv}.MuSE.vcf.gz",
    log:
        out = snakedir+"/logs/D05.som_muse/{som}.{itv}.o",
        err = snakedir+"/logs/D05.som_muse/{som}.{itv}.e",
    threads:  4
    resources:
        mem  = '16g',
        extra = ' --gres=lscratch:20 ',
    shell:
        '''module load {config[modules][muse]} {config[modules][samtools]}
        MuSE call \
          -f {config[references][fasta]} \
          {input.bam1} \
          {input.bam0} \
          -l {params.bed} \
          -O {params.prefix} > {log.out} 2> {log.err}
        MuSE sump \
          -I {params.prefix}.MuSE.txt \
          -G \
          -D {config[references][snp138]} \
          -O {params.prefix}.MuSE.vcf >> {log.out} 2>> {log.err}
        bgzip {params.prefix}.MuSE.vcf
        tabix -p vcf {output.vcf}'''
        
rule merge_muse:
    input:
        vcf1=expand("11.Somatic.MuSE/{{som}}/itv/{itv}.MuSE.vcf.gz", itv=itv4),
    output:
        vcf1="11.Somatic.MuSE/{som}/{som}.MuSE.vcf.gz",
    log:
        out = snakedir+"/logs/D05.merge_muse/{som}.o",
        err = snakedir+"/logs/D05.merge_muse/{som}.e",
    threads:  2
    resources:
        mem  = '16g',
        extra = ' --gres=lscratch:10 ',
    run:
        inputvcfs1 = ' '.join(" {}".format(i) for i in input.vcf1)
        shell('''
        ls {inputvcfs1} > {output.vcf1}
        ''')
        
rule varscan_som:
    input:
        bam1=lambda wildcards: "02.Alignment/Level3/{}/{}.BQSR.bam".format(somdic[wildcards.som][0], somdic[wildcards.som][0]),
        bam0=lambda wildcards: "02.Alignment/Level3/{}/{}.BQSR.bam".format(somdic[wildcards.som][1], somdic[wildcards.som][1]),
    params:
        bed=config['references']['wgsscatter']+'/temp_{itv}_of_50.bed',
        tmp='/lscratch/$SLURM_JOB_ID/{som}.{itv}',
        prefix="11.Somatic.Varscan2/{som}/itv/{itv}",
    output:
        loh="11.Somatic.Varscan2/{som}/itv/{itv}.snp.LOH.vcf",
    log:
        out = snakedir+"/logs/D06.som_varscan/{som}.{itv}.o",
        err = snakedir+"/logs/D06.som_varscan/{som}.{itv}.e",
    threads:  4
    resources:
        mem  = '16g',
        extra = ' --gres=lscratch:200 ',
    shell:
        '''module load {config[modules][varscan]} {config[modules][samtools]}
        samtools mpileup \
          -q 1 \
          -f {config[references][fasta]} \
          -l {params.bed} \
          {input.bam0} \
          {input.bam1} 2>{log.err} | gzip -c - > {params.tmp}.gz 2>>{log.err}
        zcat {params.tmp}.gz|cut -f 1,2,3,4,5,6 > {params.tmp}.0.mpileup 
        zcat {params.tmp}.gz|cut -f 1,2,3,7,8,9 > {params.tmp}.1.mpileup
        rm {params.tmp}.gz > {log.out} 2>>{log.err}
        varscan \
          somatic \
          {params.tmp}.0.mpileup \
          {params.tmp}.1.mpileup \
          {params.prefix} \
          --min-coverage 8 \
          --min-coverage-normal 8 \
          --min-coverage-tumor 6 \
          --min-var-freq 0.10 \
          --min-freq-for-hom 0.75 \
          --output-vcf 1 >> {log.out} 2>> {log.err}
        varscan \
          processSomatic \
          {params.prefix}.snp.vcf \
          --min-tumor-freq 0.1 \
          --max-normal-freq 0.05 \
          --p-value 0.07 >> {log.out} 2>> {log.err}
        varscan \
          processSomatic \
          {params.prefix}.indel.vcf \
          --min-tumor-freq 0.1 \
          --max-normal-freq 0.05 \
          --p-value 0.07 >> {log.out} 2>> {log.err}
        '''

rule merge_varscan:
    input:
        vcf1=expand("11.Somatic.Varscan2/{{som}}/itv/{itv}.snp.LOH.vcf", itv=itv4),
    output:
        vcf1="11.Somatic.Varscan2/{som}/{som}.snp.LOH.vcf.gz",
    log:
        out = snakedir+"/logs/D06.merge_varscan/{som}.o",
        err = snakedir+"/logs/D06.merge_varscan/{som}.e",
    threads:  2
    resources:
        mem  = '16g',
        extra = ' --gres=lscratch:10 ',
    run:
        inputvcfs1 = ' '.join(" {}".format(i) for i in input.vcf1)
        shell('''
        ls {inputvcfs1} > {output.vcf1}
        ''')



