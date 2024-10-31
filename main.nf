$HOSTNAME = ""
params.outdir = 'results'  

evaluate(new File("${params.projectDir}/nextflow_header.config"))
params.metadata.metadata = "${params.projectDir}/tools.json"
if (!params.mate){params.mate = ""} 
if (!params.reads){params.reads = ""} 
if (!params.mate2){params.mate2 = ""} 

Channel.value(params.mate).into{g_8_mate_g_1;g_8_mate_g_2;g_8_mate_g_15;g_8_mate_g_0;g_8_mate_g_16}
if (params.reads){
Channel
	.fromFilePairs( params.reads , size: params.mate == "single" ? 1 : params.mate == "pair" ? 2 : params.mate == "triple" ? 3 : params.mate == "quadruple" ? 4 : -1 )
	.ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
	.set{g_9_reads_g_15}
 } else {  
	g_9_reads_g_15 = Channel.empty()
 }

Channel.value(params.mate2).into{g_10_mate_g_6;g_10_mate_g_7;g_10_mate_g_5;g_10_mate_g_4;g_10_mate_g_3}


process unizp {

input:
 set val(name),file(reads) from g_9_reads_g_15
 val mate from g_8_mate_g_15

output:
 set val(name),file("*.fastq")  into g_15_reads0_g_0

script:

if(mate=="pair"){
	readArray = reads.toString().split(' ')	
	R1 = readArray[0]
	R2 = readArray[1]
	
	"""
	case "$R1" in
	*.gz | *.tgz ) 
	        gunzip -c $R1 > R1.fastq
	        ;;
	*)
	        cp $R1 ./R1.fastq
	        echo "$R1 not gzipped"
	        ;;
	esac
	
	case "$R2" in
	*.gz | *.tgz ) 
	        gunzip -c $R2 > R2.fastq
	        ;;
	*)
	        cp $R2 ./R2.fastq
	        echo "$R2 not gzipped"
	        ;;
	esac
	"""
}else{
	"""
	case "$reads" in
	*.gz | *.tgz ) 
	        gunzip -c $reads > R1.fastq
	        ;;
	*)
	        cp $reads ./R1.fastq
	        echo "$reads not gzipped"
	        ;;
	esac
	"""
}
}


process MaskPrimers {

input:
 val mate from g_8_mate_g_0
 set val(name),file(reads) from g_15_reads0_g_0

output:
 set val(name), file("*_primers-pass.fast*") optional true  into g_0_reads0_g_17
 set val(name), file("*_primers-fail.fast*") optional true  into g_0_reads_failed1_g_16
 set val(name), file("MP_*")  into g_0_logFile22
 set val(name),file("out*")  into g_0_logFile33

script:
method = params.MaskPrimers.method
barcode_field = params.MaskPrimers.barcode_field
primer_field = params.MaskPrimers.primer_field
barcode = params.MaskPrimers.barcode
revpr = params.MaskPrimers.revpr
mode = params.MaskPrimers.mode
failed = params.MaskPrimers.failed
fasta = params.MaskPrimers.fasta
nproc = params.MaskPrimers.nproc
maxerror = params.MaskPrimers.maxerror
umi_length = params.MaskPrimers.umi_length
start = params.MaskPrimers.start
extract_length = params.MaskPrimers.extract_length
maxlen = params.MaskPrimers.maxlen
skiprc = params.MaskPrimers.skiprc
R1_primers = params.MaskPrimers.R1_primers
R2_primers = params.MaskPrimers.R2_primers
//* @style @condition:{method="score",umi_length,start,maxerror}{method="extract",umi_length,start},{method="align",maxerror,maxlen,skiprc}, {method="extract",start,extract_length} @array:{method,barcode_field,primer_field,barcode,revpr,mode,maxerror,umi_length,start,extract_length,maxlen,skiprc} @multicolumn:{method,barcode_field,primer_field,barcode,revpr,mode,failed,nproc,maxerror,umi_length,start,extract_length,maxlen,skiprc}

method = (method.collect().size==2) ? method : [method[0],method[0]]
barcode_field = (barcode_field.collect().size==2) ? barcode_field : [barcode_field[0],barcode_field[0]]
primer_field = (primer_field.collect().size==2) ? primer_field : [primer_field[0],primer_field[0]]
barcode = (barcode.collect().size==2) ? barcode : [barcode[0],barcode[0]]
revpr = (revpr.collect().size==2) ? revpr : [revpr[0],revpr[0]]
mode = (mode.collect().size==2) ? mode : [mode[0],mode[0]]
maxerror = (maxerror.collect().size==2) ? maxerror : [maxerror[0],maxerror[0]]
umi_length = (umi_length.collect().size==2) ? umi_length : [umi_length[0],umi_length[0]]
start = (start.collect().size==2) ? start : [start[0],start[0]]
extract_length = (extract_length.collect().size==2) ? extract_length : [extract_length[0],extract_length[0]]
maxlen = (maxlen.collect().size==2) ? maxlen : [maxlen[0],maxlen[0]]
skiprc = (skiprc.collect().size==2) ? skiprc : [skiprc[0],skiprc[0]]
failed = (failed=="true") ? "--failed" : ""
fasta = (fasta=="true") ? "--fasta" : ""
def args_values = [];
[method,barcode_field,primer_field,barcode,revpr,mode,maxerror,umi_length,start,extract_length,maxlen,skiprc].transpose().each { m,bf,pf,bc,rp,md,mr,ul,s,el,ml,sk -> {
    
    if(m=="align"){
        s = ""
    }else{
        if(bc=="false"){
            s = "--start ${s}"
        }else{
            s = s + ul
            s = "--start ${s}"
        }
    }
    
    el = (m=="extract") ? "--len ${el}" : ""
    mr = (m=="extract") ? "" : "--maxerror ${mr}" 
    ml = (m=="align") ? "--maxlen ${ml}" : "" 
    sk = (m=="align" && sk=="true") ? "--skiprc" : "" 
    
    PRIMER_FIELD = "${pf}"
    
    // all
    bf = (bf=="") ? "" : "--bf ${bf}"
    pf = (pf=="") ? "" : "--pf ${pf}"
    bc = (bc=="false") ? "" : "--barcode"
    rp = (rp=="false") ? "" : "--revpr"
    args_values.add("${m} --mode ${md} ${bc} ${rp} ${mr} ${s} ${el} ${ml} ${sk} ${pf} ${bf}")
    
    
}}

readArray = reads.toString().split(' ')
if(mate=="pair"){
	args_1 = args_values[0]
	args_2 = args_values[1]
	
  


	R1 = readArray[0]
	R2 = readArray[1]
	
	R1_primers = (method[0]=="extract") ? "" : "-p ${R1_primers}"
	R2_primers = (method[1]=="extract") ? "" : "-p ${R2_primers}"
	
	
	"""
	
	MaskPrimers.py ${args_1} -s ${R1} ${R1_primers} --log MP_R1_${name}.log  --nproc ${nproc} ${failed} ${fasta} 2>&1 | tee -a out_${R1}_MP.log & \
	MaskPrimers.py ${args_2} -s ${R2} ${R2_primers} --log MP_R2_${name}.log  --nproc ${nproc} ${failed} ${fasta} 2>&1 | tee -a out_${R1}_MP.log & \
	wait
	"""
}else{
	args_1 = args_values[0]
	
	R1_primers = (method[0]=="extract") ? "" : "-p ${R1_primers}"
	
	R1 = readArray[0]

	"""
	echo -e "Assuming inputs for R1\n"
	
	MaskPrimers.py ${args_1} -s ${reads} ${R1_primers} --log MP_${name}.log  --nproc ${nproc} ${failed} ${fasta} 2>&1 | tee -a out_${R1}_MP.log
	"""
}

}


process MaskPrimers_reverse {

input:
 val mate from g_8_mate_g_16
 set val(name),file(reads) from g_0_reads_failed1_g_16

output:
 set val(name), file("*_primers-pass.fast*") optional true  into g_16_reads0_g_17
 set val(name), file("*_primers-fail.fast*") optional true  into g_16_reads_failed11
 set val(name), file("MP_*")  into g_16_logFile22
 set val(name),file("out*")  into g_16_logFile33

script:
method = params.MaskPrimers_reverse.method
barcode_field = params.MaskPrimers_reverse.barcode_field
primer_field = params.MaskPrimers_reverse.primer_field
barcode = params.MaskPrimers_reverse.barcode
revpr = params.MaskPrimers_reverse.revpr
mode = params.MaskPrimers_reverse.mode
failed = params.MaskPrimers_reverse.failed
fasta = params.MaskPrimers_reverse.fasta
nproc = params.MaskPrimers_reverse.nproc
maxerror = params.MaskPrimers_reverse.maxerror
umi_length = params.MaskPrimers_reverse.umi_length
start = params.MaskPrimers_reverse.start
extract_length = params.MaskPrimers_reverse.extract_length
maxlen = params.MaskPrimers_reverse.maxlen
skiprc = params.MaskPrimers_reverse.skiprc
R1_primers = params.MaskPrimers_reverse.R1_primers
R2_primers = params.MaskPrimers_reverse.R2_primers
//* @style @condition:{method="score",umi_length,start,maxerror}{method="extract",umi_length,start},{method="align",maxerror,maxlen,skiprc}, {method="extract",start,extract_length} @array:{method,barcode_field,primer_field,barcode,revpr,mode,maxerror,umi_length,start,extract_length,maxlen,skiprc} @multicolumn:{method,barcode_field,primer_field,barcode,revpr,mode,failed,nproc,maxerror,umi_length,start,extract_length,maxlen,skiprc}

method = (method.collect().size==2) ? method : [method[0],method[0]]
barcode_field = (barcode_field.collect().size==2) ? barcode_field : [barcode_field[0],barcode_field[0]]
primer_field = (primer_field.collect().size==2) ? primer_field : [primer_field[0],primer_field[0]]
barcode = (barcode.collect().size==2) ? barcode : [barcode[0],barcode[0]]
revpr = (revpr.collect().size==2) ? revpr : [revpr[0],revpr[0]]
mode = (mode.collect().size==2) ? mode : [mode[0],mode[0]]
maxerror = (maxerror.collect().size==2) ? maxerror : [maxerror[0],maxerror[0]]
umi_length = (umi_length.collect().size==2) ? umi_length : [umi_length[0],umi_length[0]]
start = (start.collect().size==2) ? start : [start[0],start[0]]
extract_length = (extract_length.collect().size==2) ? extract_length : [extract_length[0],extract_length[0]]
maxlen = (maxlen.collect().size==2) ? maxlen : [maxlen[0],maxlen[0]]
skiprc = (skiprc.collect().size==2) ? skiprc : [skiprc[0],skiprc[0]]
failed = (failed=="true") ? "--failed" : ""
fasta = (fasta=="true") ? "--fasta" : ""
def args_values = [];
[method,barcode_field,primer_field,barcode,revpr,mode,maxerror,umi_length,start,extract_length,maxlen,skiprc].transpose().each { m,bf,pf,bc,rp,md,mr,ul,s,el,ml,sk -> {
    
    if(m=="align"){
        s = ""
    }else{
        if(bc=="false"){
            s = "--start ${s}"
        }else{
            s = s + ul
            s = "--start ${s}"
        }
    }
    
    el = (m=="extract") ? "--len ${el}" : ""
    mr = (m=="extract") ? "" : "--maxerror ${mr}" 
    ml = (m=="align") ? "--maxlen ${ml}" : "" 
    sk = (m=="align" && sk=="true") ? "--skiprc" : "" 
    
    PRIMER_FIELD = "${pf}"
    
    // all
    bf = (bf=="") ? "" : "--bf ${bf}"
    pf = (pf=="") ? "" : "--pf ${pf}"
    bc = (bc=="false") ? "" : "--barcode"
    rp = (rp=="false") ? "" : "--revpr"
    args_values.add("${m} --mode ${md} ${bc} ${rp} ${mr} ${s} ${el} ${ml} ${sk} ${pf} ${bf}")
    
    
}}

readArray = reads.toString().split(' ')
if(mate=="pair"){
	args_1 = args_values[0]
	args_2 = args_values[1]
	
  


	R1 = readArray[0]
	R2 = readArray[1]
	
	R1_primers = (method[0]=="extract") ? "" : "-p ${R1_primers}"
	R2_primers = (method[1]=="extract") ? "" : "-p ${R2_primers}"
	
	
	"""
	
	MaskPrimers.py ${args_1} -s ${R1} ${R1_primers} --log MP_R1_${name}.log  --nproc ${nproc} ${failed} ${fasta} 2>&1 | tee -a out_${R1}_MP.log & \
	MaskPrimers.py ${args_2} -s ${R2} ${R2_primers} --log MP_R2_${name}.log  --nproc ${nproc} ${failed} ${fasta} 2>&1 | tee -a out_${R1}_MP.log & \
	wait
	"""
}else{
	args_1 = args_values[0]
	
	R1_primers = (method[0]=="extract") ? "" : "-p ${R1_primers}"
	
	R1 = readArray[0]

	"""
	echo -e "Assuming inputs for R1\n"
	
	MaskPrimers.py ${args_1} -s ${reads} ${R1_primers} --log MP_${name}.log  --nproc ${nproc} ${failed} ${fasta} 2>&1 | tee -a out_${R1}_MP.log
	"""
}

}


process combined_mask_primer {

input:
 set val(name),file(reads1) from g_0_reads0_g_17
 set val(name),file(reads2) from g_16_reads0_g_17

output:
 set val(name),file("*_output.fastq")  into g_17_reads0_g_1

script:
// Assign variables, defaulting to null if files are missing
R1 = reads1 ? reads1[0] : null
R2 = reads1 ? reads1[1] : null
R1_rev = reads2 ? reads2[1] : null
R2_rev = reads2 ? reads2[0] : null

// Determine which files to return based on availability
if (R1 && !R1_rev) {
    // R1 is available but R1_rev is missing, return R1 and R2
    """
    cp ${R1} R1_output.fastq
    cp ${R2} R2_output.fastq
    """
} else if (R1_rev && !R1) {
    // R1_rev is available but R1 is missing, return R1_rev and R2_rev
    """
    cp ${R1_rev} R1_output.fastq
    cp ${R2_rev} R2_output.fastq
    """
} else if (R1 && R1_rev) {
    // Both R1 and R1_rev are available, compare sizes
    R1_size = R1.size()
    R1_rev_size = R1_rev.size()
    
    if (R1_size > R1_rev_size) {
        // Return R1 and R2
        """
        cp ${R1} R1_output.fastq
        cp ${R2} R2_output.fastq
        """
    } else if (R1_rev_size > R1_size) {
        // Return R1_rev and R2_rev
        """
        cp ${R1_rev} R1_output.fastq
        cp ${R2_rev} R2_output.fastq
        """
    } else {
        // Sizes are relatively equal; concatenate and remove duplicates
        """
        # Concatenate files
        cat ${R1} ${R1_rev} > R1_concat.fastq
        cat ${R2} ${R2_rev} > R2_concat.fastq
    
        # Remove duplicate sequences by sequence ID
        awk 'NR % 4 == 1 {seen[\$1]++; if(seen[\$1] == 1) keep=1; else keep=0} {if(keep) print}' R1_concat.fastq > R1_output.fastq
        awk 'NR % 4 == 1 {seen[\$1]++; if(seen[\$1] == 1) keep=1; else keep=0} {if(keep) print}' R2_concat.fastq > R2_output.fastq
        """
    }
}
}


process pair_seq {

input:
 set val(name),file(reads) from g_17_reads0_g_1
 val mate from g_8_mate_g_1

output:
 set val(name),file("*_pair-pass.fastq")  into g_1_reads0_g_2
 set val(name),file("out*")  into g_1_logFile11

script:
coord = params.pair_seq.coord
act = params.pair_seq.act
copy_fields_1 = params.pair_seq.copy_fields_1
copy_fields_2 = params.pair_seq.copy_fields_2
failed = params.pair_seq.failed
nproc = params.pair_seq.nproc

if(mate=="pair"){
	
	act = (act=="none") ? "" : "--act ${act}"
	failed = (failed=="true") ? "--failed" : "" 
	copy_fields_1 = (copy_fields_1=="") ? "" : "--1f ${copy_fields_1}" 
	copy_fields_2 = (copy_fields_2=="") ? "" : "--2f ${copy_fields_2}"
	
	readArray = reads.toString().split(' ')	
	R1 = readArray[0]
	R2 = readArray[1]
	"""
	PairSeq.py -1 ${R1} -2 ${R2} ${copy_fields_1} ${copy_fields_2} --coord ${coord} ${act} ${failed} >> out_${R1}_PS.log
	"""
}else{
	
	"""
	echo -e 'PairSeq works only on pair-end reads.'
	"""
}


}


process assemble_pairs {

input:
 set val(name),file(reads) from g_1_reads0_g_2
 val mate from g_8_mate_g_2

output:
 set val(name),file("*_assemble-pass.f*")  into g_2_reads0_g_3
 set val(name),file("AP_*")  into g_2_logFile11
 set val(name),file("*_assemble-fail.f*") optional true  into g_2_reads_failed22
 set val(name),file("out*")  into g_2_logFile33

script:
method = params.assemble_pairs.method
coord = params.assemble_pairs.coord
rc = params.assemble_pairs.rc
head_fields_R1 = params.assemble_pairs.head_fields_R1
head_fields_R2 = params.assemble_pairs.head_fields_R2
failed = params.assemble_pairs.failed
fasta = params.assemble_pairs.fasta
nproc = params.assemble_pairs.nproc
alpha = params.assemble_pairs.alpha
maxerror = params.assemble_pairs.maxerror
minlen = params.assemble_pairs.minlen
maxlen = params.assemble_pairs.maxlen
scanrev = params.assemble_pairs.scanrev
minident = params.assemble_pairs.minident
evalue = params.assemble_pairs.evalue
maxhits = params.assemble_pairs.maxhits
fill = params.assemble_pairs.fill
aligner = params.assemble_pairs.aligner
// align_exec = params.assemble_pairs.// align_exec
// dbexec = params.assemble_pairs.// dbexec
gap = params.assemble_pairs.gap
usearch_version = params.assemble_pairs.usearch_version
assemble_reference = params.assemble_pairs.assemble_reference
head_seqeunce_file = params.assemble_pairs.head_seqeunce_file
//* @style @condition:{method="align",alpha,maxerror,minlen,maxlen,scanrev}, {method="sequential",alpha,maxerror,minlen,maxlen,scanrev,ref_file,minident,evalue,maxhits,fill,aligner,align_exec,dbexec} {method="reference",ref_file,minident,evalue,maxhits,fill,aligner,align_exec,dbexec} {method="join",gap} @multicolumn:{method,coord,rc,head_fields_R1,head_fields_R2,failed,nrpoc,usearch_version},{alpha,maxerror,minlen,maxlen,scanrev}, {ref_file,minident,evalue,maxhits,fill,aligner,align_exec,dbexec}, {gap} 

// args
coord = "--coord ${coord}"
rc = "--rc ${rc}"
head_fields_R1 = (head_fields_R1!="") ? "--1f ${head_fields_R1}" : ""
head_fields_R2 = (head_fields_R2!="") ? "--2f ${head_fields_R2}" : ""
failed = (failed=="false") ? "" : "--failed"
fasta = (fasta=="false") ? "" : "--fasta"
nproc = "--nproc ${nproc}"

scanrev = (scanrev=="false") ? "" : "--scanrev"
fill = (fill=="false") ? "" : "--fill"

// align_exec = (align_exec!="") ? "--exec ${align_exec}" : ""
// dbexec = (dbexec!="") ? "--dbexec ${dbexec}" : ""


ref_file = (assemble_reference!='') ? "-r ${assemble_reference}" : ""



args = ""

if(method=="align"){
	args = "--alpha ${alpha} --maxerror ${maxerror} --minlen ${minlen} --maxlen ${maxlen} ${scanrev}"
}else{
	if(method=="sequential"){
		args = "--alpha ${alpha} --maxerror ${maxerror} --minlen ${minlen} --maxlen ${maxlen} ${scanrev} ${ref_file} --minident ${minident} --evalue ${evalue} --maxhits ${maxhits} ${fill} --aligner ${aligner}"
	}else{
		if(method=="reference"){
			args = "${ref_file} --minident ${minident} --evalue ${evalue} --maxhits ${maxhits} ${fill} --aligner ${aligner}"
		}else{
			args = "--gap ${gap}"
		}
	}
}


readArray = reads.toString().split(' ')	


if(mate=="pair"){
	R1 = readArray[0]
	R2 = readArray[1]
	
	if(R1.contains(""+head_seqeunce_file)){
		R1 = readArray[0]
		R2 = readArray[1]
	}else{
		R2 = readArray[0]
		R1 = readArray[1]
	}
	
	"""
	if [ "${method}" != "align" ]; then
		if  [ "${aligner}" == "usearch" ]; then
			wget -q --show-progress --no-check-certificate https://drive5.com/downloads/usearch${usearch_version}_i86linux32.gz
			gunzip usearch${usearch_version}_i86linux32.gz
			chmod +x usearch${usearch_version}_i86linux32
			mv usearch${usearch_version}_i86linux32 /usr/local/bin/usearch2
			align_exec="--exec /usr/local/bin/usearch2"
			dbexec="--dbexec /usr/local/bin/usearch2"
		else
			align_exec="--exec /usr/local/bin/blastn"
			dbexec="--dbexec /usr/local/bin/makeblastdb"
		fi
	else
		align_exec=""
		dbexec=""
	fi

	AssemblePairs.py ${method} -1 ${R1} -2 ${R2} ${coord} ${rc} ${head_fields_R1} ${head_fields_R2} ${args} \$align_exec \$dbexec ${fasta} ${failed} --log AP_${name}.log ${nproc}  2>&1 | tee out_${R1}_AP.log
	"""

}else{
	
	"""
	echo -e 'AssemblePairs works only on pair-end reads.'
	"""
}

}


process filter_seq_quality {

input:
 set val(name),file(reads) from g_2_reads0_g_3
 val mate from g_10_mate_g_3

output:
 set val(name), file("*_${method}-pass.fast*")  into g_3_reads0_g_4
 set val(name), file("FS_*")  into g_3_logFile11
 set val(name), file("*_${method}-fail.fast*") optional true  into g_3_reads22
 set val(name),file("out*") optional true  into g_3_logFile33

script:
method = params.filter_seq_quality.method
nproc = params.filter_seq_quality.nproc
q = params.filter_seq_quality.q
n_length = params.filter_seq_quality.n_length
n_missing = params.filter_seq_quality.n_missing
window = params.filter_seq_quality.window
fasta = params.filter_seq_quality.fasta
//* @style @condition:{method="quality",q}, {method="length",n_length}, {method="missing",n_missing} @multicolumn:{method,nproc}

if(method=="missing"){
	q = ""
	n_length = ""
	window = ""
	n_missing = "-n ${n_missing}"
}else{
	if(method=="length"){
		q = ""
		n_length = "-n ${n_length}"
		n_missing = ""
		window = ""
	}else{
		if(method=="length"){
			q = "-q ${q}"
			window = "--win ${window}"
			n_length = ""
			n_missing = ""
		}else{
			q = "-q ${q}"
			n_length = ""
			n_missing = ""
			window = ""
		}
	}
}

readArray = reads.toString().split(' ')	

fasta = (fasta=="true") ? "--fasta" : ""

if(mate=="pair"){
	R1 = readArray[0]
	R2 = readArray[1]
	"""
	FilterSeq.py ${method} -s $R1 ${q} ${n_length} ${n_missing} ${window} --nproc ${nproc} --log FS_R1_${name}.log --failed ${fasta} 2>&1 | tee -a out_${R1}_FS.log
	FilterSeq.py ${method} -s $R2 ${q} ${n_length} ${n_missing} ${window} --nproc ${nproc} --log FS_R2_${name}.log --failed ${fasta} 2>&1 | tee -a out_${R1}_FS.log
	"""
}else{
	R1 = readArray[0]
	"""
	FilterSeq.py ${method} -s $R1 ${q} ${n_length} ${n_missing} ${window} --nproc ${nproc} --log FS_${name}.log --failed ${fasta} 2>&1 | tee -a out_${R1}_FS.log
	"""
}


}


process filter_seq_trimqual {

input:
 set val(name),file(reads) from g_3_reads0_g_4
 val mate from g_10_mate_g_4

output:
 set val(name), file("*_${method}-pass.fast*")  into g_4_reads0_g_5
 set val(name), file("FS_*")  into g_4_logFile11
 set val(name), file("*_${method}-fail.fast*") optional true  into g_4_reads22
 set val(name),file("out*") optional true  into g_4_logFile33

script:
method = params.filter_seq_trimqual.method
nproc = params.filter_seq_trimqual.nproc
q = params.filter_seq_trimqual.q
n_length = params.filter_seq_trimqual.n_length
n_missing = params.filter_seq_trimqual.n_missing
window = params.filter_seq_trimqual.window
fasta = params.filter_seq_trimqual.fasta
//* @style @condition:{method="quality",q}, {method="length",n_length}, {method="missing",n_missing} @multicolumn:{method,nproc}

if(method=="missing"){
	q = ""
	n_length = ""
	window = ""
	n_missing = "-n ${n_missing}"
}else{
	if(method=="length"){
		q = ""
		n_length = "-n ${n_length}"
		n_missing = ""
		window = ""
	}else{
		if(method=="length"){
			q = "-q ${q}"
			window = "--win ${window}"
			n_length = ""
			n_missing = ""
		}else{
			q = "-q ${q}"
			n_length = ""
			n_missing = ""
			window = ""
		}
	}
}

readArray = reads.toString().split(' ')	

fasta = (fasta=="true") ? "--fasta" : ""

if(mate=="pair"){
	R1 = readArray[0]
	R2 = readArray[1]
	"""
	FilterSeq.py ${method} -s $R1 ${q} ${n_length} ${n_missing} ${window} --nproc ${nproc} --log FS_R1_${name}.log --failed ${fasta} 2>&1 | tee -a out_${R1}_FS.log
	FilterSeq.py ${method} -s $R2 ${q} ${n_length} ${n_missing} ${window} --nproc ${nproc} --log FS_R2_${name}.log --failed ${fasta} 2>&1 | tee -a out_${R1}_FS.log
	"""
}else{
	R1 = readArray[0]
	"""
	FilterSeq.py ${method} -s $R1 ${q} ${n_length} ${n_missing} ${window} --nproc ${nproc} --log FS_${name}.log --failed ${fasta} 2>&1 | tee -a out_${R1}_FS.log
	"""
}


}


process filter_seq_length {

input:
 set val(name),file(reads) from g_4_reads0_g_5
 val mate from g_10_mate_g_5

output:
 set val(name), file("*_${method}-pass.fast*")  into g_5_reads0_g_6
 set val(name), file("FS_*")  into g_5_logFile11
 set val(name), file("*_${method}-fail.fast*") optional true  into g_5_reads22
 set val(name),file("out*") optional true  into g_5_logFile33

script:
method = params.filter_seq_length.method
nproc = params.filter_seq_length.nproc
q = params.filter_seq_length.q
n_length = params.filter_seq_length.n_length
n_missing = params.filter_seq_length.n_missing
window = params.filter_seq_length.window
fasta = params.filter_seq_length.fasta
//* @style @condition:{method="quality",q}, {method="length",n_length}, {method="missing",n_missing} @multicolumn:{method,nproc}

if(method=="missing"){
	q = ""
	n_length = ""
	window = ""
	n_missing = "-n ${n_missing}"
}else{
	if(method=="length"){
		q = ""
		n_length = "-n ${n_length}"
		n_missing = ""
		window = ""
	}else{
		if(method=="length"){
			q = "-q ${q}"
			window = "--win ${window}"
			n_length = ""
			n_missing = ""
		}else{
			q = "-q ${q}"
			n_length = ""
			n_missing = ""
			window = ""
		}
	}
}

readArray = reads.toString().split(' ')	

fasta = (fasta=="true") ? "--fasta" : ""

if(mate=="pair"){
	R1 = readArray[0]
	R2 = readArray[1]
	"""
	FilterSeq.py ${method} -s $R1 ${q} ${n_length} ${n_missing} ${window} --nproc ${nproc} --log FS_R1_${name}.log --failed ${fasta} 2>&1 | tee -a out_${R1}_FS.log
	FilterSeq.py ${method} -s $R2 ${q} ${n_length} ${n_missing} ${window} --nproc ${nproc} --log FS_R2_${name}.log --failed ${fasta} 2>&1 | tee -a out_${R1}_FS.log
	"""
}else{
	R1 = readArray[0]
	"""
	FilterSeq.py ${method} -s $R1 ${q} ${n_length} ${n_missing} ${window} --nproc ${nproc} --log FS_${name}.log --failed ${fasta} 2>&1 | tee -a out_${R1}_FS.log
	"""
}


}


process filter_seq_maskqual {

input:
 set val(name),file(reads) from g_5_reads0_g_6
 val mate from g_10_mate_g_6

output:
 set val(name), file("*_${method}-pass.fast*")  into g_6_reads0_g_7
 set val(name), file("FS_*")  into g_6_logFile11
 set val(name), file("*_${method}-fail.fast*") optional true  into g_6_reads22
 set val(name),file("out*") optional true  into g_6_logFile33

script:
method = params.filter_seq_maskqual.method
nproc = params.filter_seq_maskqual.nproc
q = params.filter_seq_maskqual.q
n_length = params.filter_seq_maskqual.n_length
n_missing = params.filter_seq_maskqual.n_missing
window = params.filter_seq_maskqual.window
fasta = params.filter_seq_maskqual.fasta
//* @style @condition:{method="quality",q}, {method="length",n_length}, {method="missing",n_missing} @multicolumn:{method,nproc}

if(method=="missing"){
	q = ""
	n_length = ""
	window = ""
	n_missing = "-n ${n_missing}"
}else{
	if(method=="length"){
		q = ""
		n_length = "-n ${n_length}"
		n_missing = ""
		window = ""
	}else{
		if(method=="length"){
			q = "-q ${q}"
			window = "--win ${window}"
			n_length = ""
			n_missing = ""
		}else{
			q = "-q ${q}"
			n_length = ""
			n_missing = ""
			window = ""
		}
	}
}

readArray = reads.toString().split(' ')	

fasta = (fasta=="true") ? "--fasta" : ""

if(mate=="pair"){
	R1 = readArray[0]
	R2 = readArray[1]
	"""
	FilterSeq.py ${method} -s $R1 ${q} ${n_length} ${n_missing} ${window} --nproc ${nproc} --log FS_R1_${name}.log --failed ${fasta} 2>&1 | tee -a out_${R1}_FS.log
	FilterSeq.py ${method} -s $R2 ${q} ${n_length} ${n_missing} ${window} --nproc ${nproc} --log FS_R2_${name}.log --failed ${fasta} 2>&1 | tee -a out_${R1}_FS.log
	"""
}else{
	R1 = readArray[0]
	"""
	FilterSeq.py ${method} -s $R1 ${q} ${n_length} ${n_missing} ${window} --nproc ${nproc} --log FS_${name}.log --failed ${fasta} 2>&1 | tee -a out_${R1}_FS.log
	"""
}


}


process filter_seq_missing {

input:
 set val(name),file(reads) from g_6_reads0_g_7
 val mate from g_10_mate_g_7

output:
 set val(name), file("*_${method}-pass.fast*")  into g_7_reads0_g_14
 set val(name), file("FS_*")  into g_7_logFile11
 set val(name), file("*_${method}-fail.fast*") optional true  into g_7_reads22
 set val(name),file("out*") optional true  into g_7_logFile33

script:
method = params.filter_seq_missing.method
nproc = params.filter_seq_missing.nproc
q = params.filter_seq_missing.q
n_length = params.filter_seq_missing.n_length
n_missing = params.filter_seq_missing.n_missing
window = params.filter_seq_missing.window
fasta = params.filter_seq_missing.fasta
//* @style @condition:{method="quality",q}, {method="length",n_length}, {method="missing",n_missing} @multicolumn:{method,nproc}

if(method=="missing"){
	q = ""
	n_length = ""
	window = ""
	n_missing = "-n ${n_missing}"
}else{
	if(method=="length"){
		q = ""
		n_length = "-n ${n_length}"
		n_missing = ""
		window = ""
	}else{
		if(method=="length"){
			q = "-q ${q}"
			window = "--win ${window}"
			n_length = ""
			n_missing = ""
		}else{
			q = "-q ${q}"
			n_length = ""
			n_missing = ""
			window = ""
		}
	}
}

readArray = reads.toString().split(' ')	

fasta = (fasta=="true") ? "--fasta" : ""

if(mate=="pair"){
	R1 = readArray[0]
	R2 = readArray[1]
	"""
	FilterSeq.py ${method} -s $R1 ${q} ${n_length} ${n_missing} ${window} --nproc ${nproc} --log FS_R1_${name}.log --failed ${fasta} 2>&1 | tee -a out_${R1}_FS.log
	FilterSeq.py ${method} -s $R2 ${q} ${n_length} ${n_missing} ${window} --nproc ${nproc} --log FS_R2_${name}.log --failed ${fasta} 2>&1 | tee -a out_${R1}_FS.log
	"""
}else{
	R1 = readArray[0]
	"""
	FilterSeq.py ${method} -s $R1 ${q} ${n_length} ${n_missing} ${window} --nproc ${nproc} --log FS_${name}.log --failed ${fasta} 2>&1 | tee -a out_${R1}_FS.log
	"""
}


}


process vdjbase_input {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /${chain}$/) "reads/$filename"}
input:
 set val(name),file(reads) from g_7_reads0_g_14

output:
 file "${chain}"  into g_14_germlineDb00

script:
chain = params.vdjbase_input.chain

"""
mkdir ${chain}
mv ${reads} ${chain}/${name}.fasta
"""

}


process metadata {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*.json$/) "metadata/$filename"}

output:
 file "*.json"  into g_12_jsonFile00

script:
metadata = params.metadata.metadata
"""
#!/usr/bin/env Rscript

if (!requireNamespace("jsonlite", quietly = TRUE)) {
  install.packages("jsonlite")
}
library(jsonlite)

data <- read_json("${metadata}") 

versions <- lapply(1:length(data), function(i){
	
	docker <- data[i]
	tool <- names(data)[i]
	
	if(grepl("Custom", docker)){
		ver <- "0.0"
	}else{
		ver <- system(paste0(tool," --version"), intern = TRUE)
		ver <- gsub(paste0(tool,": "), "", ver)
	}
	ver
	
})

names(versions) <- names(data)

json_data <- list(
  sample = list(
    data_processing = list(
      preprocessing = list(
        software_versions = versions
	   )
	 )
  )
)

# Convert to JSON string without enclosing scalar values in arrays
json_string <- toJSON(json_data, pretty = TRUE, auto_unbox = TRUE)
print(json_string)
# Write the JSON string to a file
writeLines(json_string, "pre_processed_metadata.json")
"""

}


workflow.onComplete {
println "##Pipeline execution summary##"
println "---------------------------"
println "##Completed at: $workflow.complete"
println "##Duration: ${workflow.duration}"
println "##Success: ${workflow.success ? 'OK' : 'failed' }"
println "##Exit status: ${workflow.exitStatus}"
}
