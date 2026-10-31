"""Registered real routes: OTU97 QIIME2/vsearch, shotgun fastp/Bowtie2/Kraken2/Bracken, matrices."""
from pathlib import Path
import json,csv,shutil
import numpy as np,pandas as pd
from workflow.common import ROOT,frame,taxonomy,tsv,dump,validate_file
from workflow.scanner import qc_pair

def amplicon(r,samples,cfg,stop=None):
    out=r.out/'upstream';p=cfg['amplicon'];threads=str(cfg['threads']);mode=cfg['layout']
    if mode not in ['SE','PE']:raise ValueError('layout must explicitly be SE or PE; best-mate selection requires explicit manifest')
    manifest=out/'qiime_manifest.tsv'
    def prepare():
        rows=[];quality=[]
        for s in samples:
            q=qc_pair(s['r1'],s.get('r2') if mode=='PE' else None)
            if q['median_R1']<p['min_median_length'] or (mode=='PE' and q['median_R2']<p['min_median_length']):raise ValueError(f"{s['sample_id']}: median length below {p['min_median_length']}; explicitly select biological mate as SE in sample manifest, never auto reassign")
            rows.append({'sample-id':s['sample_id'],**({'forward-absolute-filepath':str(Path(s['r1']).resolve()),'reverse-absolute-filepath':str(Path(s['r2']).resolve())} if mode=='PE' else {'absolute-filepath':str(Path(s['r1']).resolve())})});quality.append({'sample_id':s['sample_id'],**q})
        tsv(pd.DataFrame(rows),manifest,False);tsv(pd.DataFrame(quality),out/'fastq_qc.tsv',False)
    fastqs=[s[k] for s in samples for k in (['r1','r2'] if mode=='PE' else ['r1'])]
    r.step('prepare',fastqs,[manifest,out/'fastq_qc.tsv'],{'layout':mode,'min_median_length':p['min_median_length'],'sample_mapping':[{k:s.get(k,'') for k in ['sample_id','run_id','layout','r1','r2']} for s in samples]},prepare,plan=['internal: full FASTQ records and mate ID validation; write QIIME manifest'])
    raw=out/'raw.qza';typ='SampleData[PairedEndSequencesWithQuality]' if mode=='PE' else 'SampleData[SequencesWithQuality]';fmt='PairedEndFastqManifestPhred33V2' if mode=='PE' else 'SingleEndFastqManifestPhred33V2'
    def command(tag,inputs,outputs,cmd,params=None,deps=()):
        r.step(tag,inputs,outputs,params or {},lambda:r.cmd(cmd),deps=deps,plan=[list(map(str,cmd))])
    command('import',[manifest,*fastqs],[raw],['qiime','tools','import','--type',typ,'--input-path',manifest,'--output-path',raw,'--input-format',fmt],{'format':fmt},['prepare'])
    demux=raw
    if mode=='PE':
        demux=out/'merged.qza';unmerged=out/'unmerged.qza'
        args=['qiime','vsearch','merge-pairs','--i-demultiplexed-seqs',raw,'--p-threads',threads,'--p-minovlen',str(p['merge_min_overlap']),'--p-maxdiffs',str(p['merge_max_differences']),'--p-minlen',str(p['merge_min_length']),'--o-merged-sequences',demux,'--o-unmerged-sequences',unmerged]
        command('merge',[raw],[demux,unmerged],args,p,['import'])
    dt=out/'derep-table.qza';ds=out/'derep-seqs.qza'
    command('dereplicate',[demux],[dt,ds],['qiime','vsearch','dereplicate-sequences','--i-sequences',demux,'--p-min-seq-length',str(p['min_seq_length']),'--p-min-unique-size',str(p['min_unique_size']),'--o-dereplicated-table',dt,'--o-dereplicated-sequences',ds],p,['merge' if mode=='PE' else 'import'])
    ct=out/'clustered-table.qza';cs=out/'clustered-seqs.qza'
    command('cluster',[dt,ds],[ct,cs],['qiime','vsearch','cluster-features-de-novo','--i-table',dt,'--i-sequences',ds,'--p-perc-identity',str(p['identity']),'--p-threads',threads,'--o-clustered-table',ct,'--o-clustered-sequences',cs],p,['dereplicate'])
    ch=out/'chimeras.qza';non=out/'nonchimeras.qza';st=out/'uchime-stats.qza'
    command('chimera',[ct,cs],[ch,non,st],['qiime','vsearch','uchime-denovo','--i-table',ct,'--i-sequences',cs,'--p-dn',str(p['uchime_dn']),'--p-mindiffs',str(p['uchime_mindiffs']),'--p-mindiv',str(p['uchime_mindiv']),'--p-minh',str(p['uchime_minh']),'--p-xn',str(p['uchime_xn']),'--o-chimeras',ch,'--o-nonchimeras',non,'--o-stats',st],p,['cluster'])
    ft=out/'table.qza';fs=out/'rep-seqs.qza'
    command('filter_table',[ct,non],[ft],['qiime','feature-table','filter-features','--i-table',ct,'--m-metadata-file',non,'--o-filtered-table',ft],deps=['chimera'])
    command('filter_seqs',[cs,non],[fs],['qiime','feature-table','filter-seqs','--i-data',cs,'--m-metadata-file',non,'--o-filtered-data',fs],deps=['chimera'])
    tb=out/'features'/'feature-table.tsv';fasta=out/'features'/'dna-sequences.fasta';biom=out/'features'/'feature-table.biom'
    cmds=[['qiime','tools','export','--input-path',ft,'--output-path',tb.parent],['biom','convert','-i',biom,'-o',tb,'--to-tsv'],['qiime','tools','export','--input-path',fs,'--output-path',fasta.parent]]
    def features():
        for c in cmds:r.cmd(c)
    def valid_features():
        x=frame(tb)
        if not len(x) or (x.sum()<=0).any():raise ValueError('empty/zero final OTU feature table')
        ids={l[1:].split()[0] for l in fasta.open() if l.startswith('>')}
        if ids!=set(x.index):raise ValueError('feature vs representative sequence ID mismatch')
        if not set(x.columns)<=set(s['sample_id'] for s in samples):raise ValueError('unexpected sample IDs')
    r.step('features',[ft,fs],[tb,biom,fasta],{},features,valid_features,['filter_table','filter_seqs'],cmds)
    taxq=out/'taxonomy.qza';tax=out/'taxonomy'/'taxonomy.tsv';classifier=cfg['classifier']
    command('classify',[fs,classifier],[taxq],['qiime','feature-classifier','classify-sklearn','--i-classifier',classifier,'--i-reads',fs,'--p-confidence',str(p['classification_confidence']),'--p-read-orientation',p['read_orientation'],'--p-n-jobs',threads,'--o-classification',taxq],{'confidence':p['classification_confidence'],'orientation':p['read_orientation']},['filter_seqs'])
    def tax_valid():
        if set(taxonomy(tax).index)!=set(frame(tb).index):raise ValueError('taxonomy/feature IDs differ')
    c=['qiime','tools','export','--input-path',taxq,'--output-path',tax.parent]
    r.step('taxonomy',[taxq,tb],[tax],{},lambda:r.cmd(c),tax_valid,['classify'],[list(map(str,c))])
    return [(tb,tax,None,'counts','taxonomy')]

def shotgun(r,samples,cfg):
    """PE routes preserve original host depletion command including SAM flags.
    Bracken P/F/G/S estimates are independent ranks, never summed from species.
    """
    p=cfg['shotgun'];out=r.out/'upstream';threads=str(cfg['threads']);ranks=p['ranks']
    from scripts.abundance.historical_bracken_merge import parse_kraken_report
    byrank={rank:[] for rank in ranks}
    for s in samples:
        sid=s['sample_id'];d=out/sid;d.mkdir(parents=True,exist_ok=True) if r.mode not in ['dry-run','validate'] else None
        paired=s['layout']=='PE'
        if not paired and p['variant']!='original':raise ValueError('v2 shotgun source supports PE; use variant=original for historical SE')
        raw=[Path(s['r1'])]+([Path(s['r2'])] if paired else []);clean=[d/'clean_1.fq.gz']+([d/'clean_2.fq.gz'] if paired else []);report=d/'fastp.json';html=d/'fastp.html'
        args=['fastp','--thread',str(p['fastp_threads']),'--in1',raw[0],'--out1',clean[0],'--length_required',str(p['min_length']),'--compression',str(p['compression']),'--html',html,'--json',report]
        if paired:args+=['--in2',raw[1],'--out2',clean[1],'--detect_adapter_for_pe']
        if p['correction'] and paired:args+=['--correction']
        if p['disable_duplication_eval']:args+=['--dont_eval_duplication']
        def fastp(args=args,raw=raw):qc_pair(*raw);r.cmd(args)
        tag=sid+'.fastp';r.step(tag,raw,[*clean,report,html],p,fastp,lambda clean=clean:qc_pair(*clean),plan=[list(map(str,args))])
        prof=clean;dep=tag
        if not p['skip_host_removal']:
            prof=[d/'nonhost_1.fq']+([d/'nonhost_2.fq'] if paired else []);bam=d/'nonhost.name.bam'
            # Pipe is contained in a packaged script with pipefail; no external source dependency.
            host=[str(ROOT/'scripts/preprocess/host_filter.sh'),p['host_index'],str(cfg['threads']),p['sort_memory'],str(clean[0]),str(clean[1]) if paired else '',str(bam),str(prof[0]),str(prof[1]) if paired else '']
            host=['bash',*host];tag=sid+'.host'
            r.step(tag,[*clean,*sorted(Path(p['host_index']).parent.glob(Path(p['host_index']).name+'.*.bt2*'))],[bam,*prof],p,lambda host=host:r.cmd(host),lambda prof=prof:qc_pair(*prof),[dep],[host]);dep=tag
        kr=d/'kraken.report.tsv';ko=d/'kraken.output.tsv';args=['kraken2','--db',p['kraken_db'],'--threads',threads,'--use-names','--report',kr,'--output',ko]
        if paired:args+=['--paired']
        if p['confidence'] is not None:args+=['--confidence',str(p['confidence'])]
        if str(prof[0]).endswith('.gz'):args+=['--gzip-compressed']
        args+=prof;tag=sid+'.kraken'
        dbfiles=[Path(p['kraken_db'])/n for n in ['hash.k2d','opts.k2d','taxo.k2d']]
        r.step(tag,[*prof,*dbfiles],[kr,ko],p,lambda args=args:r.cmd(args),deps=[dep],plan=[list(map(str,args))]);dep=tag
        for rank in ranks:
            br=d/(rank+'.bracken.tsv');args=['bracken','-d',p['kraken_db'],'-i',kr,'-o',br,'-r',str(p['read_length']),'-l',rank,'-t',str(p['bracken_threshold'])];tag=sid+'.bracken.'+rank
            def brvalid(br=br):
                b=pd.read_csv(br,sep='\t');assert len(b)>0 and {'taxonomy_id','new_est_reads','fraction_total_reads'}<=set(b.columns),'invalid/empty Bracken'
                assert (b.new_est_reads>=0).all(),'negative estimate'
            r.step(tag,[kr,Path(p['kraken_db'])/f"database{p['read_length']}mers.kmer_distrib"],[br],p,lambda args=args:r.cmd(args),brvalid,[dep],[list(map(str,args))]);byrank[rank].append((sid,br,kr,tag))
    results=[];names={'S':'Species','G':'Genus','F':'Family','P':'Phylum'}
    for rank,items in byrank.items():
        tb=out/rank/'counts.tsv';tx=out/rank/'taxonomy.tsv';rel=out/rank/'bracken_native_relative.tsv'
        def merge(items=items,tb=tb,tx=tx,rel=rel):
            counts={};rels={};lineages={}
            for sid,br,kr,tag in items:
                b=pd.read_csv(br,sep='\t',dtype={'taxonomy_id':str});lm=parse_kraken_report(kr)
                if b.taxonomy_id.duplicated().any():raise ValueError(f'{br}: duplicate taxids')
                b.index='taxid_'+b.taxonomy_id
                counts[sid]=b.new_est_reads;rels[sid]=b.fraction_total_reads
                for key,row in b.iterrows():
                    lin=lm.get(row.taxonomy_id,'')
                    if key in lineages and lineages[key]!=lin:raise ValueError(f'taxid lineage conflict {key}')
                    lineages[key]=lin
            # Sparse observed taxid absence represents a true taxonomic zero, not a missing sample/file.
            tsv(pd.DataFrame(counts).fillna(0),tb);tsv(pd.DataFrame(rels).fillna(0),rel)
            tsv(pd.DataFrame({'Feature ID':list(lineages),'Taxon':list(lineages.values())}),tx,False)
        r.step('shotgun_merge_'+rank,[p for _,br,kr,_ in items for p in [br,kr]],[tb,tx,rel],{'rank':rank},merge,lambda tb=tb:frame(tb),[tag for *_,tag in items])
        results.append((tb,tx,names[rank],'counts','shotgun_merge_'+rank))
    return results
